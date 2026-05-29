#!/usr/bin/env python
"""Integrator-agnostic position-level checkpoint helper for production drivers.

Long unattended production runs (notably the KU-3.5 cortical-tension sweep) can
die mid-flight — e.g. an LJ-CFL ``int32`` image-flag overflow after myosin
saturation — losing hours of CPU with no way to resume.  This module provides a
small, dependency-free checkpoint/restore facility so a driver can persist its
sampling progress and pick up where it left off.

Design
------
* **Position-level only.**  A checkpoint records the per-tag particle positions
  at the last completed sample, the running sample index, and the accumulated
  ``frames`` / ``diag`` history.  On resume the driver *rebuilds* the simulation
  from scratch and writes the saved positions back into the state.

* **Integrator RNG history is NOT persisted.**  The frozen BAOAB / constrained
  BAOAB integrators own their own thermostat RNG stream, and this module
  deliberately does not reach into them.  Consequently the thermostat noise
  history restarts at the resume boundary.  This is a *minor statistical seam*,
  not a correctness break: the position state (the physically meaningful
  configuration) is restored exactly, and sample intervals are long compared
  with the velocity-autocorrelation time, so the re-randomised momenta
  thermalise well within a single interval.  It requires **no change to the
  frozen integrator files**.

* **Atomic writes.**  ``save_checkpoint`` writes to a ``.tmp`` sibling then
  ``os.replace``s it onto the final ``.ckpt.npz`` path, so a crash *during* a
  checkpoint write can never corrupt a previously good checkpoint.

* **Fingerprint-gated resume.**  Each checkpoint embeds a ``params_fingerprint``
  hash of the run-defining parameters.  ``load_checkpoint`` returns the saved
  state only when the caller's fingerprint matches exactly, so a resume can
  never silently splice incompatible runs (different ``n_fil`` / ``dt`` / seed /
  sampling schedule).
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
from pathlib import Path
from typing import Any

import numpy as np

LOGGER = logging.getLogger("ffn_sim.common.checkpoint")

#: Suffix appended to the caller-supplied base path to form the checkpoint file.
CKPT_SUFFIX = ".ckpt.npz"


def checkpoint_path(path: str | os.PathLike[str]) -> Path:
    """Return the canonical checkpoint path derived from a base path.

    The base path is typically the driver's ``--out`` JSON/npz path; the
    checkpoint lives beside it with the :data:`CKPT_SUFFIX` suffix (any existing
    suffix on ``path`` is stripped first so ``foo.json`` and ``foo.npz`` map to
    the same ``foo.ckpt.npz``).

    Args:
        path: Base output path (e.g. the driver ``--out`` value).

    Returns:
        :class:`~pathlib.Path` to the ``.ckpt.npz`` file.
    """
    base = Path(path)
    return base.with_suffix("").with_suffix(CKPT_SUFFIX)


def compute_fingerprint(params: dict[str, Any]) -> str:
    """Hash run-defining parameters into a stable fingerprint string.

    The fingerprint is order-independent (keys are sorted) and stable across
    processes (``sha256`` over a canonical JSON encoding).  Callers should pass
    every parameter that, if changed, must invalidate an old checkpoint —
    e.g. ``n_fil``, ``dt`` / ``dt_factor``, ``seed``, ``n_sample``,
    ``interval``, ``n_warmup``, and any physics toggles.

    Args:
        params: Mapping of run-defining parameter names to JSON-serialisable
            values.

    Returns:
        Hex ``sha256`` digest of the canonicalised parameters.
    """
    canonical = json.dumps(params, sort_keys=True, separators=(",", ":"),
                           default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def save_checkpoint(
    path: str | os.PathLike[str],
    *,
    sample_index: int,
    positions_by_tag: np.ndarray,
    diag_list: list[dict[str, Any]],
    frames_so_far: list[np.ndarray] | np.ndarray,
    params_fingerprint: str,
    rng_seed: int,
) -> Path:
    """Atomically persist a position-level checkpoint.

    Writes a ``.tmp`` file then :func:`os.replace`s it onto the final
    ``.ckpt.npz`` path so an interrupted write cannot corrupt a prior good
    checkpoint.

    Args:
        path: Base output path; the checkpoint is written to
            :func:`checkpoint_path`(``path``).
        sample_index: Number of samples *completed* so far (i.e. the loop index
            from which a resume should continue).
        positions_by_tag: ``(N, 3)`` array of the latest per-tag particle
            positions, as returned by the driver's ``_tagpos`` helper.
        diag_list: Accumulated per-sample diagnostic dicts.  Serialised as JSON.
        frames_so_far: Accumulated position frames — a list of arrays or a
            single stacked array.
        params_fingerprint: Fingerprint from :func:`compute_fingerprint`,
            embedded so resume can verify parameter identity.
        rng_seed: The run's RNG seed (stored for diagnostics/traceability; the
            integrator RNG *history* is intentionally not persisted).

    Returns:
        The :class:`~pathlib.Path` the checkpoint was written to.
    """
    ckpt = checkpoint_path(path)
    ckpt.parent.mkdir(parents=True, exist_ok=True)

    frames_arr = np.asarray(frames_so_far, dtype=np.float32)
    positions = np.asarray(positions_by_tag, dtype=np.float64)

    # Diagnostics are heterogeneous dicts; store as a single JSON blob so we do
    # not impose a fixed schema and can round-trip exactly.
    diag_json = json.dumps(diag_list, default=str)

    tmp = ckpt.with_suffix(ckpt.suffix + ".tmp")
    # Give np.savez an explicit file handle so the temp filename is exactly what
    # we os.replace from (np.savez would otherwise append ".npz").
    with open(tmp, "wb") as fh:
        np.savez(
            fh,
            sample_index=np.int64(sample_index),
            positions_by_tag=positions,
            frames_so_far=frames_arr,
            diag_json=np.asarray(diag_json),
            params_fingerprint=np.asarray(params_fingerprint),
            rng_seed=np.int64(rng_seed),
        )
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, ckpt)
    LOGGER.info(
        "checkpoint saved: %s (sample_index=%d, frames=%d)",
        ckpt, sample_index, frames_arr.shape[0],
    )
    return ckpt


def load_checkpoint(
    path: str | os.PathLike[str],
    params_fingerprint: str,
) -> dict[str, Any] | None:
    """Load a checkpoint iff it exists and its fingerprint matches.

    Args:
        path: Base output path (same value passed to :func:`save_checkpoint`).
        params_fingerprint: Expected fingerprint from
            :func:`compute_fingerprint`.  A mismatch (or missing file) yields
            ``None`` so the caller falls back to a fresh run.

    Returns:
        On a hit, a dict with keys:

        * ``sample_index`` (:class:`int`) — completed-sample count to resume from
        * ``positions_by_tag`` (:class:`numpy.ndarray`, ``(N, 3)``)
        * ``frames_so_far`` (:class:`list` of per-sample arrays)
        * ``diag_list`` (:class:`list` of dicts)
        * ``rng_seed`` (:class:`int`)

        Otherwise ``None``.
    """
    ckpt = checkpoint_path(path)
    if not ckpt.exists():
        LOGGER.info("no checkpoint at %s; starting fresh", ckpt)
        return None
    try:
        with np.load(ckpt, allow_pickle=False) as data:
            saved_fp = str(data["params_fingerprint"])
            if saved_fp != params_fingerprint:
                LOGGER.warning(
                    "checkpoint fingerprint mismatch at %s "
                    "(saved=%s..., expected=%s...); ignoring and starting fresh",
                    ckpt, saved_fp[:12], params_fingerprint[:12],
                )
                return None
            frames_arr = np.asarray(data["frames_so_far"], dtype=np.float32)
            diag_list = json.loads(str(data["diag_json"]))
            state: dict[str, Any] = {
                "sample_index": int(data["sample_index"]),
                "positions_by_tag": np.asarray(
                    data["positions_by_tag"], dtype=np.float64
                ),
                # Hand back a list so the driver can keep appending in-place.
                "frames_so_far": [np.asarray(f, dtype=np.float32)
                                  for f in frames_arr],
                "diag_list": diag_list,
                "rng_seed": int(data["rng_seed"]),
            }
    except Exception as exc:  # noqa: BLE001  (corrupt/partial file → fresh run)
        LOGGER.warning(
            "failed to read checkpoint %s (%s); starting fresh", ckpt, exc
        )
        return None

    LOGGER.info(
        "checkpoint loaded: %s (resume at sample_index=%d, frames=%d)",
        ckpt, state["sample_index"], len(state["frames_so_far"]),
    )
    return state

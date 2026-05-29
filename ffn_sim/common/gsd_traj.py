"""GSD trajectory writer for ffn_cellsim particle+bond simulations.

This module is the renderer-agnostic "unlock" for journal-quality visualization.
Production H.3 runs persist trajectories as ``.npz`` with a single ``frames``
array shaped ``(n_frames, n_filaments, beads, 3)`` -- positions only, with no
bond topology and no per-particle scalar fields. That layout is renderer-hostile:
it cannot drive bond/cylinder rendering and cannot color particles by mechanical
state.

GSD (General Simulation Data, the native HOOMD-blue 7.x format) carries the full
schema instead:

* particle positions, type ids, and the type-name table,
* the bond group ``(M, 2)``, bond type ids, and the bond type-name table,
* the simulation box,
* arbitrary per-particle scalar fields, stored under the frame ``log`` namespace
  (``particles/<name>``), which OVITO surfaces as ``log/particles/<name>`` and
  the :mod:`gsd` Python API reads straight back.

Two entry points are exposed:

* :func:`write_gsd_trajectory` -- offline writer, used by
  :mod:`ffn_sim.scripts.npz_to_gsd` to re-export an existing ``.npz`` (or by any
  analysis code holding per-frame arrays) into a single multi-frame ``.gsd``.
* :func:`attach_gsd_writer` -- best-effort helper that attaches a live
  ``hoomd.write.GSD`` writer (plus an optional ``hoomd.logging.Logger``) to a
  running :class:`hoomd.Simulation`, so future production drivers can emit GSD
  directly instead of round-tripping through ``.npz``.

``gsd.hoomd`` is imported eagerly (it is a hard dependency of the visualization
stack and is present in the ``ffn_sim`` env). ``hoomd`` is imported lazily inside
:func:`attach_gsd_writer` so that pure offline conversion works in environments
without a configured HOOMD device.

Sanity Gate:
    * Dimensional: positions/box/diameter are in simulation length units (sigma);
      scalar fields keep their native units (the renderer annotates them). No unit
      conversion happens here -- the writer is unit-agnostic by design.
    * Boundary: a single ``(N, 3)`` frame is promoted to one frame; an empty
      scalar dict is allowed; ``bonds=None`` yields a GSD with no bonds.
    * Conservation: topology (N, bonds, types, box) is written identically on
      every frame; only positions and per-frame scalars vary -- matching a
      fixed-connectivity cytoskeletal mesh.
    * Shape checks: positions must be ``(...,3)``; each scalar field's particle
      dimension must equal ``N``; bonds must be ``(M, 2)`` -- all asserted with an
      explicit ``ValueError``.
"""

from __future__ import annotations

import logging
from typing import Any, Mapping, Sequence

import gsd.hoomd
import numpy as np
import numpy.typing as npt

logger = logging.getLogger(__name__)

__all__ = ["write_gsd_trajectory", "attach_gsd_writer"]


def _normalize_box(box: Sequence[float] | npt.ArrayLike) -> list[float]:
    """Coerce a box specification into the 6-vector GSD/HOOMD expects.

    HOOMD/GSD boxes are stored as ``[Lx, Ly, Lz, xy, xz, yz]`` (edge lengths plus
    tilt factors). This accepts a scalar (cube edge), a 3-vector of edge lengths
    (tilts zero, orthorhombic), or a full 6-vector.

    Args:
        box: Scalar cube edge, ``(Lx, Ly, Lz)``, or ``(Lx, Ly, Lz, xy, xz, yz)``.

    Returns:
        The box as a length-6 list of floats.

    Raises:
        ValueError: If ``box`` has none of 1, 3, or 6 entries.
    """
    arr = np.asarray(box, dtype=np.float64).ravel()
    if arr.size == 1:
        e = float(arr[0])
        return [e, e, e, 0.0, 0.0, 0.0]
    if arr.size == 3:
        return [float(arr[0]), float(arr[1]), float(arr[2]), 0.0, 0.0, 0.0]
    if arr.size == 6:
        return [float(v) for v in arr]
    raise ValueError(
        f"box must have 1 (cube edge), 3 (Lx,Ly,Lz), or 6 (Lx,Ly,Lz,xy,xz,yz) "
        f"entries; got {arr.size}"
    )


def write_gsd_trajectory(
    path: str,
    positions: npt.ArrayLike,
    box: Sequence[float] | npt.ArrayLike,
    *,
    type_ids: npt.ArrayLike | None = None,
    type_names: Sequence[str] | None = None,
    bonds: npt.ArrayLike | None = None,
    bond_type_ids: npt.ArrayLike | None = None,
    bond_type_names: Sequence[str] | None = None,
    scalar_fields: Mapping[str, npt.ArrayLike] | None = None,
    diameter: float | None = None,
    log_constants: Mapping[str, float] | None = None,
) -> str:
    """Write a multi-frame GSD trajectory from per-frame position arrays.

    The particle count ``N``, the bond topology, the type tables, and the box are
    treated as constant across frames (the connectivity of a cytoskeletal mesh
    does not change within a single H.3 production run). Only positions and the
    per-particle scalar fields vary frame-to-frame. Scalar fields are written into
    ``frame.log`` under ``particles/<name>`` so they survive as named per-particle
    properties (read back by OVITO as ``log/particles/<name>`` and consumed by
    :mod:`ffn_sim.scripts.render_fresnel`).

    Args:
        path: Output ``.gsd`` path. Overwritten if it exists (mode ``'w'``).
        positions: Per-frame positions. Either ``(n_frames, N, 3)`` or a single
            ``(N, 3)`` frame (promoted to one frame).
        box: Cube edge, ``(Lx, Ly, Lz)``, or ``(Lx, Ly, Lz, xy, xz, yz)``.
        type_ids: Per-particle integer type ids, shape ``(N,)``. Defaults to all
            zeros (one type).
        type_names: Type-name table indexed by type id. Defaults to ``['A']`` (or
            ``['A','B',...]`` sized to ``max(type_ids)+1``).
        bonds: Bond group as integer particle-index pairs, shape ``(M, 2)``.
            Optional; if ``None`` the GSD carries no bonds.
        bond_type_ids: Per-bond integer type ids, shape ``(M,)``. Defaults to all
            zeros when ``bonds`` is given.
        bond_type_names: Bond type-name table. Defaults to ``['backbone']`` when
            ``bonds`` is given.
        scalar_fields: Mapping of field name to either ``(n_frames, N)`` (one row
            per frame) or ``(N,)`` (broadcast to every frame). Written to
            ``frame.log['particles/<name>']`` as ``float32``.
        diameter: Optional uniform particle diameter (simulation length units),
            stored on every frame so renderers can size spheres from the file.
        log_constants: Optional mapping of scalar metadata written to every
            frame's ``log`` under the given key (no ``particles/`` prefix), e.g.
            ``{'sigma_um': 0.1}`` so a renderer can self-describe SI units without
            re-reading the source ``.npz``. Each value is stored as a length-1
            ``float32`` array (GSD requires array-valued log entries).

    Returns:
        The output ``path`` (for chaining / logging).

    Raises:
        ValueError: On shape mismatches (positions not ``(...,3)``; a scalar
            field whose particle dimension does not match ``N``; bonds not
            ``(M, 2)``; type/bond-type id length mismatches).
    """
    pos = np.asarray(positions, dtype=np.float32)
    if pos.ndim == 2:
        pos = pos[None, ...]
    if pos.ndim != 3 or pos.shape[2] != 3:
        raise ValueError(
            f"positions must be (n_frames, N, 3) or (N, 3); got shape {pos.shape}"
        )
    n_frames, n_particles, _ = pos.shape

    box6 = _normalize_box(box)

    # --- particle types -----------------------------------------------------
    if type_ids is None:
        tid = np.zeros(n_particles, dtype=np.uint32)
    else:
        tid = np.asarray(type_ids, dtype=np.uint32).ravel()
        if tid.size != n_particles:
            raise ValueError(
                f"type_ids has {tid.size} entries but there are {n_particles} particles"
            )
    if type_names is None:
        n_types = int(tid.max()) + 1 if tid.size else 1
        tnames = [chr(ord("A") + i) if i < 26 else f"T{i}" for i in range(n_types)]
    else:
        tnames = list(type_names)

    # --- bonds --------------------------------------------------------------
    bgroup: np.ndarray | None = None
    btid: np.ndarray | None = None
    bnames: list[str] | None = None
    if bonds is not None:
        bgroup = np.ascontiguousarray(np.asarray(bonds, dtype=np.uint32))
        if bgroup.ndim != 2 or bgroup.shape[1] != 2:
            raise ValueError(f"bonds must be (M, 2); got shape {bgroup.shape}")
        n_bonds = bgroup.shape[0]
        if bond_type_ids is None:
            btid = np.zeros(n_bonds, dtype=np.uint32)
        else:
            btid = np.asarray(bond_type_ids, dtype=np.uint32).ravel()
            if btid.size != n_bonds:
                raise ValueError(
                    f"bond_type_ids has {btid.size} entries but there are "
                    f"{n_bonds} bonds"
                )
        if bond_type_names is None:
            n_btypes = int(btid.max()) + 1 if btid.size else 1
            bnames = (
                ["backbone"] if n_btypes == 1 else [f"bond{i}" for i in range(n_btypes)]
            )
        else:
            bnames = list(bond_type_names)

    # --- scalar fields: validate + normalize to (n_frames, N) ----------------
    norm_fields: dict[str, np.ndarray] = {}
    if scalar_fields:
        for name, raw in scalar_fields.items():
            arr = np.asarray(raw, dtype=np.float32)
            if arr.ndim == 1:
                if arr.size != n_particles:
                    raise ValueError(
                        f"scalar field {name!r} has {arr.size} values but there are "
                        f"{n_particles} particles"
                    )
                arr = np.broadcast_to(arr, (n_frames, n_particles))
            elif arr.ndim == 2:
                if arr.shape != (n_frames, n_particles):
                    raise ValueError(
                        f"scalar field {name!r} has shape {arr.shape}, expected "
                        f"{(n_frames, n_particles)}"
                    )
            else:
                raise ValueError(
                    f"scalar field {name!r} must be 1-D (N,) or 2-D (n_frames, N); "
                    f"got ndim {arr.ndim}"
                )
            norm_fields[name] = np.ascontiguousarray(arr)

    diam = None
    if diameter is not None:
        diam = np.full(n_particles, float(diameter), dtype=np.float32)

    const_log = {
        str(k): np.asarray([float(v)], dtype=np.float32)
        for k, v in (log_constants or {}).items()
    }

    with gsd.hoomd.open(name=path, mode="w") as traj:
        for fi in range(n_frames):
            frame = gsd.hoomd.Frame()
            frame.configuration.step = fi
            frame.configuration.box = box6

            frame.particles.N = n_particles
            frame.particles.position = np.ascontiguousarray(pos[fi])
            frame.particles.typeid = tid
            frame.particles.types = tnames
            if diam is not None:
                frame.particles.diameter = diam

            if bgroup is not None:
                frame.bonds.N = int(bgroup.shape[0])
                frame.bonds.group = bgroup
                frame.bonds.typeid = btid
                frame.bonds.types = bnames

            for name, arr in norm_fields.items():
                # gsd.hoomd stores arbitrary arrays in frame.log; the
                # 'particles/' prefix marks it as a per-particle quantity so
                # OVITO surfaces it as 'log/particles/<name>'.
                frame.log[f"particles/{name}"] = np.ascontiguousarray(arr[fi])

            for name, val in const_log.items():
                frame.log[name] = val

            traj.append(frame)

    logger.info(
        "Wrote GSD trajectory %s: %d frame(s), %d particles, %d bond(s), fields=%s",
        path,
        n_frames,
        n_particles,
        0 if bgroup is None else int(bgroup.shape[0]),
        sorted(norm_fields) or "none",
    )
    return path


def attach_gsd_writer(
    sim: Any,
    path: str,
    period: int,
    logger_obj: Any | None = None,
    *,
    mode: str = "wb",
    dynamic: Sequence[str] | None = None,
    filter_all: bool = True,
) -> Any:
    """Attach a live ``hoomd.write.GSD`` writer to a running simulation.

    Best-effort convenience for future production drivers: instead of buffering
    positions in RAM and dumping a positions-only ``.npz``, a driver can call this
    once after building its :class:`hoomd.Simulation` to stream a fully-schema'd
    GSD to disk every ``period`` steps. Pass a configured
    ``hoomd.logging.Logger`` as ``logger_obj`` to additionally persist
    per-particle scalar quantities under ``log/particles/*`` (e.g. a custom
    tension/curvature computer registered with ``category='particle'``).

    ``hoomd`` is imported lazily so this module's offline writer stays usable in
    environments without a HOOMD device.

    Args:
        sim: A ``hoomd.Simulation`` instance with an attached ``device`` and
            ``state``. Typed as ``Any`` to avoid importing ``hoomd`` at module
            load.
        path: Output ``.gsd`` path.
        period: Write every ``period`` timesteps (a ``hoomd.trigger.Periodic``).
        logger_obj: Optional ``hoomd.logging.Logger`` whose logged quantities are
            embedded in each frame's ``log`` namespace.
        mode: GSD file open mode (default ``'wb'``, truncating). Use ``'ab'`` to
            append/resume.
        dynamic: Quantities allowed to vary per frame, forwarded to
            ``hoomd.write.GSD`` (e.g. ``['property', 'momentum']``). When ``None``,
            HOOMD's default (positions/orientation) applies and topology is written
            once.
        filter_all: When ``True`` (default) the writer captures every particle via
            ``hoomd.filter.All()``. Set ``False`` only if you intend to attach a
            custom filter to the returned writer yourself.

    Returns:
        The attached ``hoomd.write.GSD`` writer (already appended to
        ``sim.operations.writers``), so the caller can tweak it further.

    Raises:
        ImportError: If ``hoomd`` cannot be imported.
        AttributeError: If ``sim`` lacks the expected ``operations`` interface.
    """
    import hoomd  # lazy: offline conversion must not require a HOOMD device

    trigger = hoomd.trigger.Periodic(int(period))
    kwargs: dict[str, Any] = {"trigger": trigger, "filename": path, "mode": mode}
    if filter_all:
        kwargs["filter"] = hoomd.filter.All()
    if dynamic is not None:
        kwargs["dynamic"] = list(dynamic)
    if logger_obj is not None:
        kwargs["logger"] = logger_obj

    writer = hoomd.write.GSD(**kwargs)
    sim.operations.writers.append(writer)
    logger.info(
        "Attached hoomd.write.GSD -> %s (period=%d, mode=%s, logger=%s)",
        path,
        period,
        mode,
        "yes" if logger_obj is not None else "no",
    )
    return writer

#!/usr/bin/env python
"""Convert a positions-only ``.npz`` trajectory into a fully-schema'd GSD.

Production H.3 drivers persist trajectories as ``.npz`` with a ``frames`` array
shaped ``(n_frames, n_filaments, beads, 3)`` (positions only). Some drivers may
additionally store a per-bead ``tension`` array and a ``params`` object array,
but the canonical production outputs observed in
``ffn_sim/outputs/h3/production`` carry **only** ``frames`` -- so every extra is
treated as optional with a safe fallback. The positions-only layout cannot drive
bond/cylinder rendering or per-particle coloring. This CLI re-exports such an
``.npz`` into a multi-frame ``.gsd`` (via :mod:`ffn_sim.common.gsd_traj`)
carrying:

* flattened per-bead positions ``(N, 3)`` with ``N = n_filaments * beads``,
* intra-filament backbone bonds reconstructed from the array layout (bead ``j``
  bonds to bead ``j+1`` within each filament; no inter-filament bonds, matching
  the mesh's explicit topology),
* one ``filament`` particle type for every bead,
* per-bead scalar fields under ``log/particles/*``:
    - ``curvature``: local discrete bending energy from the three-bead angle
      (always computed from positions, so the field exists even for a
      positions-only ``.npz``),
    - ``tension``: passed through from the ``.npz`` if present (per-bead scalar;
      a per-bead 3-vector is reduced to its Euclidean magnitude).

Curvature definition (discrete worm-like-chain bending energy):
for an interior bead ``i`` with neighbours ``i-1`` and ``i+1`` and segment unit
vectors ``u = (r_i - r_{i-1})/|.|``, ``v = (r_{i+1} - r_i)/|.|``, the bending
energy is ``E_bend = kappa * (1 - cos(theta))`` with ``cos(theta) = u . v`` and
``kappa`` the angle stiffness in kT (read from ``params['k_angle_kt']`` /
``params['k_angle']`` when present, else 1.0 so the field is still a faithful
``1 - cos(theta)`` curvature proxy). Endpoints have no defined angle and are
assigned 0. This mirrors the cosine bending used by the H.3 filament model, so the
colored field is the genuine per-bead bending energy, not a viz proxy.

Run::

    python ffn_sim/scripts/npz_to_gsd.py --in traj.npz --out traj.gsd [--frame-stride N]

Sanity Gate:
    * Dimensional: positions/box/diameter stay in sigma (no SI rescale at write
      time; the renderer annotates SI via ``params['sigma_um']`` if present, else
      its own default). Curvature carries the kT prefactor; tension keeps native
      units.
    * Boundary: ``beads < 2`` -> no bonds; ``beads < 3`` -> all-zero curvature (no
      interior beads); missing ``params``/``box`` -> auto cube box sized to the
      data; missing ``tension`` -> only the curvature field is written.
    * Conservation: bond count is exactly ``n_filaments*(beads-1)`` and every bond
      stays within one filament.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np

from ffn_sim.common.gsd_traj import write_gsd_trajectory

logger = logging.getLogger(__name__)


def build_backbone_bonds(n_filaments: int, beads: int) -> np.ndarray:
    """Construct intra-filament backbone bonds for a flattened bead array.

    Beads are flattened filament-major: the global index of bead ``j`` in filament
    ``f`` is ``f * beads + j``. Consecutive beads in the same filament are bonded;
    there are no inter-filament bonds.

    Args:
        n_filaments: Number of filaments.
        beads: Beads per filament (``>= 2`` to produce any bonds).

    Returns:
        Integer array of shape ``(n_filaments * (beads - 1), 2)`` of global bead
        index pairs. Empty ``(0, 2)`` if ``beads < 2``.
    """
    if beads < 2:
        return np.empty((0, 2), dtype=np.int64)
    # within-filament left endpoints 0..beads-2, broadcast across filaments
    offsets = (np.arange(n_filaments, dtype=np.int64) * beads)[:, None]
    left = (offsets + np.arange(beads - 1, dtype=np.int64)[None, :]).ravel()
    return np.stack([left, left + 1], axis=1)


def compute_curvature(frames: np.ndarray, kappa: float = 1.0) -> np.ndarray:
    """Compute per-bead discrete bending energy for every frame.

    Vectorized over frames and filaments. For each interior bead the bending
    energy is ``kappa * (1 - cos(theta))`` with ``theta`` the angle between the two
    incident backbone segments; endpoints are 0.

    Args:
        frames: Positions shaped ``(n_frames, n_filaments, beads, 3)``.
        kappa: Angle stiffness prefactor (kT). With the default ``1.0`` the field
            is the bare ``1 - cos(theta)`` curvature proxy.

    Returns:
        Array shaped ``(n_frames, n_filaments, beads)`` of bending energies
        (non-negative; endpoints 0).
    """
    nf, nfil, beads, _ = frames.shape
    curv = np.zeros((nf, nfil, beads), dtype=np.float64)
    if beads < 3:
        return curv
    r = frames.astype(np.float64)
    d = r[:, :, 1:, :] - r[:, :, :-1, :]  # segment vectors (nf,nfil,beads-1,3)
    seg_len = np.linalg.norm(d, axis=-1)
    seg_len_safe = np.where(seg_len > 0, seg_len, 1.0)
    u = d / seg_len_safe[..., None]
    cos_t = np.sum(u[:, :, :-1, :] * u[:, :, 1:, :], axis=-1)  # (nf,nfil,beads-2)
    cos_t = np.clip(cos_t, -1.0, 1.0)
    valid = (seg_len[:, :, :-1] > 0) & (seg_len[:, :, 1:] > 0)
    curv[:, :, 1:-1] = np.where(valid, kappa * (1.0 - cos_t), 0.0)
    return curv


def _reduce_tension(
    tension: np.ndarray, frames_shape: tuple[int, ...]
) -> np.ndarray | None:
    """Coerce a ``.npz`` ``tension`` array into per-bead scalars per frame.

    Accepts a per-bead scalar ``(F, fil, beads)`` (passed through) or a per-bead
    vector ``(F, fil, beads, 3)`` (reduced to Euclidean magnitude). Anything else
    is rejected (returns ``None`` with a warning).

    Args:
        tension: Raw ``tension`` array from the ``.npz``.
        frames_shape: Shape of the ``frames`` array, ``(F, fil, beads, 3)``.

    Returns:
        Array shaped ``(F, fil, beads)`` or ``None`` if the layout is unusable.
    """
    nf, nfil, beads, _ = frames_shape
    t = np.asarray(tension)
    if t.shape == (nf, nfil, beads):
        return t.astype(np.float64)
    if t.shape == (nf, nfil, beads, 3):
        return np.linalg.norm(t.astype(np.float64), axis=-1)
    logger.warning(
        "tension array shape %s does not match frames %s (per-bead scalar or "
        "3-vector expected); skipping tension field",
        t.shape,
        frames_shape,
    )
    return None


def _unpack_params(data) -> dict:
    """Extract the ``params`` object array from a loaded ``.npz`` as a dict.

    Returns an empty dict if there is no ``params`` entry or it is not dict-like.
    """
    if "params" not in data.files:
        return {}
    raw = data["params"]
    try:
        pp = raw.item() if getattr(raw, "shape", None) == () else raw[()]
    except Exception:  # pragma: no cover - defensive
        logger.warning("could not unpack 'params'; proceeding with defaults")
        return {}
    return dict(pp) if hasattr(pp, "keys") else {}


def convert_npz_to_gsd(
    in_path: str | Path,
    out_path: str | Path,
    frame_stride: int = 1,
) -> str:
    """Convert one trajectory ``.npz`` to a multi-frame ``.gsd``.

    Args:
        in_path: Input ``.npz`` with at least a ``frames`` array shaped
            ``(n_frames, n_filaments, beads, 3)``. Optional ``tension`` and
            ``params`` entries are used when present.
        out_path: Output ``.gsd`` path (overwritten).
        frame_stride: Keep every ``frame_stride``-th frame (``>= 1``).

    Returns:
        The output path as a string.

    Raises:
        KeyError: If the ``.npz`` has no ``frames`` array.
        ValueError: If ``frames`` is not 4-D ``(F, fil, beads, 3)`` or
            ``frame_stride < 1``.
    """
    in_path = Path(in_path)
    out_path = Path(out_path)
    if frame_stride < 1:
        raise ValueError(f"frame_stride must be >= 1, got {frame_stride}")

    data = np.load(in_path, allow_pickle=True)
    if "frames" not in data.files:
        raise KeyError(
            f"{in_path} has no 'frames' array (found: {data.files}); cannot convert"
        )
    frames = data["frames"]
    if frames.ndim != 4 or frames.shape[-1] != 3:
        raise ValueError(
            f"'frames' must be (n_frames, n_filaments, beads, 3); got {frames.shape}"
        )

    params = _unpack_params(data)

    frames = frames[::frame_stride]
    nf, nfil, beads, _ = frames.shape
    n_particles = nfil * beads
    logger.info(
        "%s: %d frames (stride %d) x %d filaments x %d beads = %d particles/frame",
        in_path.name,
        nf,
        frame_stride,
        nfil,
        beads,
        n_particles,
    )

    # --- box (production key 'box_L_sigma' = cubic edge in sigma; else auto) --
    box_edge = params.get("box_L_sigma", params.get("box"))
    if box_edge is None:
        span = float(np.max(np.abs(frames)) * 2.0) if frames.size else 1.0
        box = (max(span, 1.0) * 1.1,) * 3
        logger.warning(
            "no 'box_L_sigma'/'box' in params; using auto cube edge=%.3g sigma",
            box[0],
        )
    else:
        b = np.asarray(box_edge, dtype=np.float64).ravel()
        box = (float(b[0]),) * 3 if b.size == 1 else tuple(float(x) for x in b)

    # --- flatten positions to (nf, N, 3) ------------------------------------
    positions = frames.reshape(nf, n_particles, 3).astype(np.float32)

    # --- bonds --------------------------------------------------------------
    bonds = build_backbone_bonds(nfil, beads)
    logger.info("reconstructed %d backbone bonds", bonds.shape[0])

    # --- curvature (always; kappa from k_angle_kt/k_angle) ------------------
    kappa = float(params.get("k_angle_kt", params.get("k_angle", 1.0)))
    curv = compute_curvature(frames, kappa=kappa).reshape(nf, n_particles)
    scalar_fields: dict[str, np.ndarray] = {"curvature": curv}

    # --- tension (if present) ----------------------------------------------
    if "tension" in data.files:
        tred = _reduce_tension(data["tension"][::frame_stride], frames.shape)
        if tred is not None:
            scalar_fields["tension"] = tred.reshape(nf, n_particles)
            logger.info("attached per-bead tension field")

    # uniform bead diameter ~ rest length so spheres just touch on the backbone
    diameter = float(params.get("L0_sigma", params.get("L0", 1.0)))

    # carry the sigma->um anchor into the GSD so the renderer is self-describing
    log_constants: dict[str, float] = {}
    if "sigma_um" in params:
        log_constants["sigma_um"] = float(params["sigma_um"])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_gsd_trajectory(
        str(out_path),
        positions,
        box,
        type_ids=np.zeros(n_particles, dtype=np.uint32),
        type_names=["filament"],
        bonds=bonds,
        bond_type_ids=np.zeros(bonds.shape[0], dtype=np.uint32),
        bond_type_names=["backbone"],
        scalar_fields=scalar_fields,
        diameter=diameter,
        log_constants=log_constants,
    )
    return str(out_path)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    p = argparse.ArgumentParser(
        description="Convert a positions-only trajectory .npz to a schema'd .gsd."
    )
    p.add_argument("--in", dest="in_path", required=True, help="input .npz path")
    p.add_argument("--out", dest="out_path", required=True, help="output .gsd path")
    p.add_argument(
        "--frame-stride",
        type=int,
        default=1,
        help="keep every Nth frame (default 1 = all frames)",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="DEBUG logging")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    out = convert_npz_to_gsd(args.in_path, args.out_path, frame_stride=args.frame_stride)
    logger.info("wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Jamming order parameter for the DCM aggregate — the 3D cell shape index.

The jamming↔unjamming (solid↔fluid) state of a cell tissue is set by the dimensionless CELL SHAPE INDEX
(Bi-Manning vertex-model theory; the freshly-ingested TAG jamming/SPV literature). In 3D (Merkel & Manning
2018) it is s = S / V^(2/3) (cell surface area over volume^(2/3)):

    s < s0* ≈ 5.41  →  JAMMED / solid-like  (arrested, confined; a sphere is s=4.836)
    s > s0* ≈ 5.41  →  UNJAMMED / fluid-like (cells rearrange / flow / spread)

(The 2D analogue is p0 = P/√A with p0* = 3.81, pentagon.) Active self-propulsion (SPV model,
Bi-Yang-Marchetti-Manning) fluidises a tissue even below s0* by adding the (v0, persistence) axis. We use the
aggregate-mean s as the runtime jamming order parameter: confined DCM spheroid ≈ 4.93 (jammed); active motility
drives it toward/above s0* = the unjamming / spreading transition. Measured, never tuned.
"""
from __future__ import annotations

import numpy as np

S0_STAR_3D = 5.41       # Merkel-Manning 3D rigidity (unjamming) threshold
S_SPHERE_3D = 4.836     # dimensionless shape index of a sphere (minimum)
P0_STAR_2D = 3.81       # Bi-Manning 2D shape-index threshold (pentagon)


def cell_shape_index_3d(pos: np.ndarray, faces: np.ndarray, cof: np.ndarray) -> np.ndarray:
    """Per-cell 3D shape index s = S / V^(2/3) (S = closed-surface area, V = enclosed volume).

    Args:
        pos: (N,3) node positions (any consistent length unit — s is dimensionless).
        faces: (M,3) triangle node indices.
        cof: (N,) cell-of-node (face owner = cof of its first vertex).
    Returns:
        (n_cells,) shape index per cell (nan for degenerate cells).
    """
    pos = np.asarray(pos, float)
    faces = np.asarray(faces, int)
    cof = np.asarray(cof, int)
    nc = int(cof.max()) + 1
    v0, v1, v2 = pos[faces[:, 0]], pos[faces[:, 1]], pos[faces[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    vol = np.einsum("ij,ij->i", v0, np.cross(v1, v2)) / 6.0        # signed tetra volume
    fcell = cof[faces[:, 0]]
    S = np.zeros(nc); V = np.zeros(nc)
    np.add.at(S, fcell, area)
    np.add.at(V, fcell, vol)
    V = np.abs(V)
    s = np.full(nc, np.nan)
    good = V > 1e-30
    s[good] = S[good] / np.power(V[good], 2.0 / 3.0)
    return s


def aggregate_jamming_state(pos, faces, cof):
    """Aggregate-mean shape index + jamming verdict. Returns dict(s_mean, s_std, jammed, frac_unjammed)."""
    s = cell_shape_index_3d(pos, faces, cof)
    s = s[np.isfinite(s)]
    if s.size == 0:
        return {"s_mean": float("nan"), "s_std": float("nan"), "jammed": None, "frac_unjammed": float("nan")}
    s_mean = float(s.mean())
    return {
        "s_mean": s_mean,
        "s_std": float(s.std()),
        "jammed": bool(s_mean < S0_STAR_3D),
        "frac_unjammed": float((s > S0_STAR_3D).mean()),   # fraction of cells past the rigidity threshold
    }

"""Engine-agnostic DCM shell geometry — extracted from the legacy ``cell/dcm.py``.

The Deformable-Cell-Model shell mesh (icosphere nodes/edges/tris) and its
resolved parameter set (:class:`ResolvedDCM`) are pure geometry / data — no
HOOMD, no Warp. They are the seam shared by the Warp DCM runtime and the
archived HOOMD reference, so they live in the active ``dcm/`` layer (canonical),
not in the retired HOOMD tree (deleted 2026-07-29 at `1a9ded66`).

Extracted verbatim (behaviour-preserving) as step 1b of the FF/DCM two-layer
split; the only dependency is :class:`aleph.laws.surface_manifold.SurfaceManifold`
(also engine-agnostic).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from aleph.laws.surface_manifold import SurfaceManifold

_KT_310 = 4.28e-21  # J at 310 K


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------
def icosphere_mesh(R: float, subdivisions: int = 1):
    """Return (verts (n,3) [m], edges (m,2) int, tris (k,3) int) of a sphere.

    Triangles are outward-wound (used for the exact enclosed-volume / face-normal
    pressure in the turgor force).
    """
    mani = SurfaceManifold.icosphere(subdivisions, R)
    verts = np.asarray(mani.verts, dtype=np.float64)
    tris = np.asarray(mani.tris, dtype=np.int64)
    # Ensure outward winding (normal · centroid-direction > 0).
    cen = verts.mean(axis=0)
    fixed = []
    for a, b, c in tris:
        n = np.cross(verts[b] - verts[a], verts[c] - verts[a])
        if np.dot(n, verts[a] - cen) < 0:
            a, b, c = a, c, b
        fixed.append((int(a), int(b), int(c)))
    tris = np.array(fixed, dtype=np.int64)
    eset = set()
    for a, b, c in tris:
        for u, v in ((a, b), (b, c), (c, a)):
            eset.add((int(min(u, v)), int(max(u, v))))
    edges = np.array(sorted(eset), dtype=np.int64)
    return verts, edges, tris


def _cluster_centers(n_cells: int, spacing: float, z0: float, R: float,
                     mode: str = "3d"):
    """Compact cluster of cell centers resting on the substrate (z0).

    ``mode='2d'`` — hex-packed monolayer (all at z0+R). ``mode='3d'`` — a compact
    FCC-like ball (a real spheroid) sitting on the substrate (lowest cell at z0+R).
    Returns (n_cells, 3) [m].
    """
    if mode == "2d":
        offs = [(0.0, 0.0)]
        s3 = np.sqrt(3) / 2
        offs += [(1, 0), (-1, 0), (0.5, s3), (-0.5, s3), (0.5, -s3), (-0.5, -s3)]
        offs += [(2, 0), (-2, 0), (1.5, s3), (-1.5, s3), (1.5, -s3), (-1.5, -s3),
                 (1, 2 * s3), (-1, 2 * s3), (0, 2 * s3), (1, -2 * s3),
                 (-1, -2 * s3), (0, -2 * s3)]
        offs = offs[:n_cells]
        return np.array([[ox * spacing, oy * spacing, z0 + R] for ox, oy in offs],
                        dtype=np.float64)

    # 3D: generate FCC lattice points, keep the n_cells closest to the origin
    # (a compact ball), then drop so the lowest cell rests on the substrate.
    # FCC nearest-neighbour distance = a/√2, so set the conventional edge
    # a = spacing·√2 to make neighbouring cell centres sit ``spacing`` apart
    # (else cells overlap by ~spacing(1−1/√2) → LJ-core blow-up).
    a = spacing * np.sqrt(2.0)
    pts = []
    rng = range(-3, 4)
    for i in rng:
        for j in rng:
            for k in rng:
                # FCC basis
                for bx, by, bz in ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)):
                    pts.append(((i + bx) * a, (j + by) * a, (k + bz) * a))
    pts = np.array(pts, dtype=np.float64)
    d = np.linalg.norm(pts - pts.mean(0), axis=1)
    order = np.argsort(d)[:n_cells]
    centers = pts[order]
    centers -= centers.mean(0)                      # center at origin
    centers[:, 2] -= centers[:, 2].min()            # lowest cell base at 0
    centers[:, 2] += z0 + R                          # rest lowest cell on substrate
    return centers


# ---------------------------------------------------------------------------
# Resolved parameters (literature bands, trusted directly)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedDCM:
    R_cell: float = 7.5e-6          # m  MCF7 radius (Wagner 2011)
    subdivisions: int = 1          # icosphere level (1=42 nodes, 2=162)
    turgor_dP0: float = 133.0      # Pa  MCF7 baseline osmotic turgor (Young-Laplace)
    K_vol: float = 1.0e3           # Pa  osmotic bulk modulus (ΔP per ΔV/V)
    k_edge: float = 1.0e-3         # N/m membrane edge spring (cortical elasticity)
    gamma_node: float = 3.9e-10    # N·s/m per-node Stokes drag (cortex gamma_b)
    # Adhesion energy DENSITIES (literature bands, trusted directly): cell-substrate
    # and cell-cell ~0.1-0.5 mJ/m² (cadherin / integrin scale). Per-node well depth
    # = W [J/m²] × area_per_node, computed at build. Adhesion ≫ kT (per node ~1e6 kT),
    # so the spread is set by W_cs / cortical-tension (wetting), not by thermal.
    W_cs_Jm2: float = 0.5e-3       # J/m² cell-substrate adhesion (spreading driver)
    W_cc_Jm2: float = 0.2e-3       # J/m² cell-cell adhesion (aggregate cohesion)
    eps_rep_kT: float = 5.0        # excluded-volume (WCA) strength [kT]
    cluster: str = "3d"            # "3d" spheroid ball | "2d" monolayer aggregate
    spacing_factor: float = 2.1    # cell-center spacing = factor·R_cell (light contact)
    # --- ECM / substrate attachment (improvement layer) ---
    # Substrate COMPLIANCE: finite stiffness softens the felt adhesion (a cell sinks
    # into a soft gel → less traction → less spreading) = the durotaxis / rigidity-
    # sensing knob. None ⇒ rigid dish (glass, the PI's prep). k_sub ~ E_sub·R
    # (Chan-Odde ~1-4 pN/nm ≈ 1e-3..4e-3 N/m per clutch; here a per-node floor).
    k_sub_Nm: float | None = None  # N/m substrate spring stiffness (None = rigid)
    # LIGAND DENSITY: cell-substrate adhesion energy ∝ areal ligand density
    # (Gallant/García; saturating). 1.0 = reference (~750 ligands/µm²).
    ligand_density: float = 1.0
    z_substrate: float = 0.0       # m substrate plane
    kT: float = _KT_310
    dt: float = 1.0e-9             # s (strong adhesion → small dt for BAOAB stability)
    seed: int = 7
    extras: dict[str, Any] = field(default_factory=dict)

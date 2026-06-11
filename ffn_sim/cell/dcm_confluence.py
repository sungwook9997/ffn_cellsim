"""Confluence metrics + position capture for the native-mesh DCM spheroid.

NEW helper (2026-06-11) for the "confluent packing" deliverable. None of the
owned/frozen DCM modules are touched; this only *reads* a built simulation's
state (via ``cpu_local_snapshot``) and the per-cell triangulation handles that
``build_native_dcm_simulation`` already returns.

The goal is to quantify how SimuCell3D-like a relaxed aggregate is: do the
deformable cells press into each other and fill space (a confluent tissue), or
do they stay as separated rounded spheres with gaps? Three complementary,
geometry-only metrics are provided (no magic numbers — each is a ratio of
measured volumes / counts):

  * **packing fraction** ``Φ = Σ V_cell / V_hull`` — total cell volume over the
    convex-hull volume of all nodes. Separated spheres in a loose ball leave big
    interstitial voids ⇒ Φ well below ~0.7; a confluent (space-filling, foam-like)
    aggregate pushes Φ toward / past random-close-pack (~0.64) and beyond, since
    deformed cells fill the interstices.
  * **contact fraction** — fraction of surface nodes that lie within the cell-cell
    contact range of a node of a *different* cell (i.e. nodes participating in an
    adhered interface). Separated spheres ⇒ near 0; confluent ⇒ a large fraction
    of each cell's surface is shared interface.
  * **gap fraction** ``1 − Φ`` — the complementary void fraction, reported for
    readability.

``V_cell`` uses the same HOOMD-convention enclosed volume as the turgor force
(``mesh_volume_of``) so it is consistent with what the simulation conserves.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ffn_sim.cell.dcm_native_shell import mesh_volume_of


def capture_positions(sim) -> np.ndarray:
    """Return tag-ordered node positions (N,3) [m] from a running sim."""
    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position).copy()
    order = np.argsort(tag)
    return pos[order]


@dataclass(slots=True)
class ConfluenceMetrics:
    packing_fraction: float       # Σ V_cell / V_hull
    gap_fraction: float           # 1 − packing_fraction
    contact_fraction: float       # fraction of nodes on an adhered interface
    mean_contact_nodes: float     # mean inter-cell contact nodes per cell
    cell_volumes: np.ndarray      # (n_cells,) enclosed volume per cell [m^3]
    hull_volume: float            # convex-hull volume of all nodes [m^3]


def compute_confluence(pos: np.ndarray, *, cell_of_tag: np.ndarray,
                       ranges, mesh_tris_per_cell, contact_range: float
                       ) -> ConfluenceMetrics:
    """Confluence metrics for one frame.

    Args:
        pos: (N,3) tag-ordered node positions [m].
        cell_of_tag: (N,) cell id per node (global/tag order).
        ranges: list of (lo, hi) tag ranges, one per cell (from the builder).
        mesh_tris_per_cell: list of (k,3) LOCAL triangle index arrays, one per
            cell (indices into that cell's node block). See ``per_cell_tris``.
        contact_range: m, the cell-cell contact range used to count adhered
            interface nodes (use the tent ``c_adh``).
    """
    from scipy.spatial import ConvexHull, cKDTree

    n_cells = len(ranges)
    vols = np.empty(n_cells, dtype=np.float64)
    for c, (lo, hi) in enumerate(ranges):
        verts = pos[lo:hi]
        vols[c] = mesh_volume_of(verts, mesh_tris_per_cell[c])

    hull = ConvexHull(pos)
    hull_vol = float(hull.volume)
    phi = float(vols.sum() / hull_vol)

    # contact / interface nodes: a node within contact_range of a node of a
    # different cell counts as an adhered-interface node.
    tree = cKDTree(pos)
    pairs = tree.query_pairs(r=contact_range, output_type="ndarray")
    interface = np.zeros(pos.shape[0], dtype=bool)
    per_cell_contact = np.zeros(n_cells, dtype=np.int64)
    if pairs.size:
        ci = cell_of_tag[pairs[:, 0]]
        cj = cell_of_tag[pairs[:, 1]]
        diff = ci != cj
        pi = pairs[diff, 0]
        pj = pairs[diff, 1]
        interface[pi] = True
        interface[pj] = True
        for node in np.concatenate([pi, pj]):
            per_cell_contact[cell_of_tag[node]] += 1
    contact_frac = float(interface.mean())
    mean_contact_nodes = float(per_cell_contact.mean())

    return ConfluenceMetrics(
        packing_fraction=phi,
        gap_fraction=1.0 - phi,
        contact_fraction=contact_frac,
        mean_contact_nodes=mean_contact_nodes,
        cell_volumes=vols,
        hull_volume=hull_vol,
    )


def per_cell_tris(mesh_tris, mesh_typeids, ranges):
    """Split the global mesh triangle list into per-cell LOCAL-index lists.

    ``build_native_snapshot`` emits a single global ``mesh_tris`` (indices into
    the full particle array) plus ``mesh_typeids`` (the cell id of each triangle).
    For per-cell surface rendering AND per-cell enclosed-volume we want each
    cell's triangles re-indexed into that cell's own node block [0, nv).
    Returns a list of (k,3) int arrays, one per cell.
    """
    mesh_tris = np.asarray(mesh_tris)
    mesh_typeids = np.asarray(mesh_typeids)
    out = []
    for c, (lo, hi) in enumerate(ranges):
        sel = mesh_typeids == c
        tri = mesh_tris[sel].astype(np.int64) - lo
        out.append(tri)
    return out

"""DCM Upgrade A — GPU-native HOOMD mesh shell for the deformable cell.

SimuCell3D integration (2026-06-11). This module replaces the custom-Python
``DcmTurgorForce`` (a per-step ``cpu_local_snapshot`` Custom force) with
HOOMD-blue 7.0.1 **native mesh potentials** that run GPU-resident:

  * ``hoomd.md.mesh.conservation.Volume`` — per-cell enclosed-volume turgor.
    Energy ``U = k·(V − V0)²/(2·V0)`` ⇒ pressure ``p = −dU/dV = −(k/V0)·(V − V0)``.
    Setting ``k = K·V0`` makes ``p = −K·(V − V0)/V0`` with ``K`` the osmotic bulk
    modulus (Pa) — the SimuCell3D ``K`` (mcf7_p0.xml: 2500 Pa). Replaces
    ``DcmTurgorForce``.
  * ``hoomd.md.mesh.bending.BendingRigidity`` — discrete-shell bending (k_b),
    ``U = ½·k_b·Σ(1 − cos θ)`` over adjacent-face dihedrals. OFF by default
    (SimuCell3D ships bending_modulus = 0; the round-up is carried by tension +
    volume + contact). Enable via ``k_bend > 0``.
  * edge springs — kept as ``hoomd.md.bond.Harmonic`` on the SAME mesh-generated
    bond topology (in-plane cortical elasticity / surface-tension proxy). Native
    ``md.mesh.conservation.Area`` is available too but harmonic edges are the
    stable, validated choice (Upgrade A keeps the existing edge-spring contract).

WINDING / SIGN (the load-bearing gotcha)
----------------------------------------
HOOMD's mesh ``Volume`` uses the OPPOSITE orientation convention from the
divergence-theorem volume used by ``DcmTurgorForce``: an outward-wound triangle
list (numpy signed volume > 0) yields a NEGATIVE HOOMD mesh volume, and HOOMD
then rejects ``V0 <= 0`` and applies no restoring force. Reversing the triangle
winding (``tris[:, ::-1]`` etc.) makes HOOMD's mesh **bond builder** raise
"the same particle can only occur once in a triangle" (an internal adjacency
ordering quirk). The robust fix used here: keep the consistently-oriented
``SurfaceManifold.icosphere`` triangulation untouched (it attaches cleanly) and
**mirror one spatial axis (x → −x)** of every shell's vertices. A reflection
flips the handedness — hence the sign of the enclosed volume HOOMD computes —
WITHOUT reordering any triangle, so HOOMD reports a positive volume and the bond
builder is happy. The physics is unchanged (a mirrored icosphere is still an
icosphere); cell-cell contact and substrate forces are reflection-invariant.

Each cell is a disjoint closed mesh (no shared nodes), carried as a distinct
**mesh triangle type** ``"cell{c}"`` inside one ``hoomd.mesh.Mesh`` so that
``Volume`` conserves every cell's volume independently.

Units SI. Defaults: MCF7 bands + SimuCell3D mcf7_p0.xml (K = 2500 Pa).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.common.surface_manifold import SurfaceManifold

_KT_310 = 4.28e-21  # J at 310 K


# ---------------------------------------------------------------------------
# Mesh primitive (consistently oriented + x-mirrored for HOOMD's sign convention)
# ---------------------------------------------------------------------------
def native_icosphere(R: float, subdivisions: int = 2):
    """Return (verts (n,3) [m], edges (m,2) int, tris (k,3) int) for ONE shell.

    Triangulation is the consistently-oriented ``SurfaceManifold.icosphere``
    (each directed edge appears exactly once → a valid closed manifold that
    HOOMD's mesh bond builder accepts). Vertices are **x-mirrored** so that
    HOOMD's mesh ``Volume`` reports a positive enclosed volume (see module
    docstring). Edges are the unique undirected edges of the triangulation.

    subdivisions=2 (162 nodes / 320 tris) is recommended for the native shell
    (enough basal nodes for contact + bending); subdivisions=1 (42/80) is the
    coarse/fast option.
    """
    mani = SurfaceManifold.icosphere(subdivisions, R)
    verts = np.asarray(mani.verts, dtype=np.float64).copy()
    verts[:, 0] *= -1.0  # reflection → flips HOOMD's volume sign to positive
    tris = np.asarray(mani.tris, dtype=np.int64)
    eset: set[tuple[int, int]] = set()
    for a, b, c in tris:
        for u, v in ((a, b), (b, c), (c, a)):
            eset.add((int(min(u, v)), int(max(u, v))))
    edges = np.array(sorted(eset), dtype=np.int64)
    return verts, edges, tris


def mesh_volume_of(verts: np.ndarray, tris: np.ndarray) -> float:
    """HOOMD-convention enclosed volume of a single shell (= −numpy signed vol).

    Use this to set ``V0`` so it matches what HOOMD's ``Volume`` will measure at
    rest (the divergence-theorem sphere formula over-estimates by ~13% at
    subdiv 1, ~3% at subdiv 2, which would otherwise inflate/deflate the cell).
    """
    v0, v1, v2 = verts[tris[:, 0]], verts[tris[:, 1]], verts[tris[:, 2]]
    signed = float(np.einsum("ij,ij->i", v0, np.cross(v1 - v0, v2 - v0)).sum() / 6.0)
    return -signed  # HOOMD's convention is the negative of numpy's


# ---------------------------------------------------------------------------
# Resolved parameters for the native-shell DCM
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedNativeDCM:
    """Native-mesh DCM constants (MCF7 bands + SimuCell3D mcf7_p0.xml)."""

    R_cell: float = 7.5e-6          # m  MCF7 radius
    subdivisions: int = 2          # icosphere level (2 = 162 nodes/cell)
    # --- shell energetics ---
    K_bulk: float = 2500.0         # Pa  osmotic bulk modulus (SimuCell3D mcf7)
    turgor_inflate: float = 1.05   # V0 = inflate · V_rest (resting +turgor; >1 ⇒ pressurised)
    k_edge: float = 1.0e-3         # N/m harmonic edge spring (cortical elasticity)
    k_bend: float = 0.0            # bending rigidity (0 = OFF, SimuCell3D default)
    gamma_node: float = 3.9e-10    # N·s/m per-node Stokes drag
    # --- substrate (reuses DcmSubstrateForce) ---
    W_cs_Jm2: float = 0.5e-3       # J/m² cell-substrate adhesion (spreading driver)
    ligand_density: float = 1.0
    k_sub_Nm: float | None = None
    z_substrate: float = 0.0
    # --- cell-cell contact (SimuCell3D tent; see dcm_contact.py) ---
    rep_strength: float = 1.0e8    # Pa/m  repulsion (mcf7)
    adh_strength: float = 1.0e8    # Pa/m  adhesion ω (mcf7; set 0 to disable)
    c_adh: float = 5.0e-7          # m  contact range (mcf7)
    force_cap: float = 5.0e-8      # N  per-node contact force cap (BAOAB guard)
    # --- assembly / numerics ---
    cluster: str = "3d"
    spacing_factor: float = 2.3    # cell-center spacing = factor · R (no t=0 overlap)
    kT: float = _KT_310
    dt: float = 3.0e-10            # s (small for BAOAB stability; drop to 1e-10 if needed)
    seed: int = 7
    extras: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Snapshot builder — N disjoint native shells
# ---------------------------------------------------------------------------
def build_native_snapshot(p: ResolvedNativeDCM, n_cells: int, centers: np.ndarray):
    """Build a snapshot of ``n_cells`` disjoint x-mirrored icosphere shells.

    Returns (gsd Frame, mesh_triangles (F,3) uint, mesh_type_ids (F,) uint,
             mesh_types list[str], cell_of_tag (N,), ranges, mean_edge, nv,
             V_rest_per_cell float).
    All cells share one membrane particle type ``"mem"`` (fixed type set); cell
    identity for the native ``Volume`` is the per-triangle MESH type id.
    """
    import gsd.hoomd

    verts0, edges0, tris0 = native_icosphere(p.R_cell, p.subdivisions)
    nv = verts0.shape[0]
    mean_edge = float(np.linalg.norm(
        verts0[edges0[:, 0]] - verts0[edges0[:, 1]], axis=1).mean())
    V_rest = mesh_volume_of(verts0, tris0)

    all_pos, all_bonds, all_tris, all_typeids = [], [], [], []
    cell_of_tag, ranges = [], []
    mesh_types = [f"cell{c}" for c in range(n_cells)]
    tag = 0
    for c in range(n_cells):
        ranges.append((tag, tag + nv))
        all_pos.append(verts0 + centers[c])
        all_bonds.append(edges0 + tag)
        all_tris.append(tris0 + tag)
        all_typeids.append(np.full(tris0.shape[0], c, dtype=np.uint32))
        cell_of_tag.append(np.full(nv, c, dtype=np.int64))
        tag += nv

    pos = np.concatenate(all_pos, axis=0)
    bonds = np.concatenate(all_bonds, axis=0)
    mesh_tris = np.concatenate(all_tris, axis=0).astype(np.uint32)
    mesh_typeids = np.concatenate(all_typeids)
    cell_of_tag = np.concatenate(cell_of_tag)
    N = pos.shape[0]

    box_edge = float(np.abs(pos).max() * 2.5 + 6.0 * p.R_cell)
    snap = gsd.hoomd.Frame()
    snap.particles.N = N
    snap.particles.types = ["mem"]
    snap.particles.position = pos
    snap.particles.typeid = np.zeros(N, dtype=np.uint32)
    snap.particles.mass = np.ones(N)
    snap.bonds.N = bonds.shape[0]
    snap.bonds.types = ["dcm_edge"]
    snap.bonds.typeid = np.zeros(bonds.shape[0], dtype=np.uint32)
    snap.bonds.group = bonds.astype(np.uint32)
    snap.configuration.box = [box_edge, box_edge, box_edge, 0, 0, 0]
    return (snap, mesh_tris, mesh_typeids, mesh_types, cell_of_tag, ranges,
            mean_edge, nv, V_rest)


def build_native_shell_forces(sim, p: ResolvedNativeDCM, *, mesh_tris,
                              mesh_typeids, mesh_types, V_rest, mean_edge):
    """Attach native mesh shell potentials to ``sim`` and return (forces, mesh).

    Returns (list[force], mesh, volume_force). The caller appends ``forces`` to
    its ``md.Integrator``. ``volume_force.volume`` logs per-cell volumes.
    """
    mesh = hoomd.mesh.Mesh()
    mesh.types = list(mesh_types)
    mesh.triangulation = dict(type_ids=mesh_typeids, triangles=mesh_tris)

    forces = []

    # Volume turgor (one Volume force, per-cell mesh type). k = K·V0 ⇒ p = -K·(V-V0)/V0.
    V0 = p.turgor_inflate * V_rest
    vol = md.mesh.conservation.Volume(mesh)
    for t in mesh_types:
        vol.params[t] = dict(k=p.K_bulk * V0, V0=V0)
    forces.append(vol)

    # Edge springs (cortical elasticity) on the mesh-generated bond topology.
    bond = md.bond.Harmonic()
    bond.params["dcm_edge"] = dict(k=p.k_edge, r0=mean_edge)
    forces.append(bond)

    # Bending (OFF by default; SimuCell3D ships k_b = 0).
    if p.k_bend > 0.0:
        bend = md.mesh.bending.BendingRigidity(mesh)
        for t in mesh_types:
            bend.params[t] = dict(k=p.k_bend)
        forces.append(bend)

    return forces, mesh, vol, V0


# ---------------------------------------------------------------------------
# Full assembly — native shell + SimuCell3D tent contact + substrate + BAOAB
# ---------------------------------------------------------------------------
def build_native_dcm_simulation(p: ResolvedNativeDCM, n_cells: int, *,
                                device=None, substrate: bool = True,
                                contact: bool = True):
    """Assemble the native-mesh DCM spheroid on the BAOAB integrator.

    Forces: native ``Volume`` turgor + harmonic edge springs (+ optional native
    bending), ``DcmTentContact`` cell-cell (Upgrade B), and the kept
    ``DcmSubstrateForce`` adhesive floor. Returns a dict of handles.
    """
    from ffn_sim.integrator.baoab import make_baoab_updater
    from ffn_sim.cell.dcm import DcmSubstrateForce, _cluster_centers
    from ffn_sim.cell.dcm_contact import DcmTentContact

    spacing = p.spacing_factor * p.R_cell
    centers = _cluster_centers(n_cells, spacing, p.z_substrate, p.R_cell,
                               mode=p.cluster)

    (snap, mesh_tris, mesh_typeids, mesh_types, cell_of_tag, ranges,
     mean_edge, nv, V_rest) = build_native_snapshot(p, n_cells, centers)

    dev = device or hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=p.seed)
    sim.create_state_from_snapshot(snap)

    ig = md.Integrator(dt=p.dt)
    forces, mesh, vol, V0 = build_native_shell_forces(
        sim, p, mesh_tris=mesh_tris, mesh_typeids=mesh_typeids,
        mesh_types=mesh_types, V_rest=V_rest, mean_edge=mean_edge)
    for f in forces:
        ig.forces.append(f)

    tent = None
    if contact:
        patch_area = 4.0 * np.pi * p.R_cell ** 2 / nv
        tent = DcmTentContact(
            cell_of_node=cell_of_tag, r_contact=1.05 * mean_edge,
            c_adh=p.c_adh, rep_strength=p.rep_strength,
            adh_strength=p.adh_strength, patch_area=patch_area,
            force_cap=p.force_cap)
        ig.forces.append(tent)

    sub = None
    if substrate:
        area_per_node = 4.0 * np.pi * p.R_cell ** 2 / nv
        W_cs = p.W_cs_Jm2 * p.ligand_density * area_per_node
        sub = DcmSubstrateForce(z0=p.z_substrate, W_cs=W_cs,
                                adh_range=p.R_cell, k_sub=p.k_sub_Nm)
        ig.forces.append(sub)

    sim.operations.integrator = ig
    sim.run(0)

    gamma = {"mem": p.gamma_node}
    action, updater = make_baoab_updater(kT=p.kT, gamma=gamma, dt=p.dt,
                                         seed=p.seed + 1)
    sim.operations.updaters.append(updater)

    return {
        "sim": sim, "mesh": mesh, "volume": vol, "tent": tent, "substrate": sub,
        "cell_of_tag": cell_of_tag, "ranges": ranges, "n_cells": n_cells,
        "nv": nv, "mean_edge": mean_edge, "V_rest": V_rest, "V0": V0, "p": p,
    }

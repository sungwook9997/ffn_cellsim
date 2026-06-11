"""GPU-FRIENDLY DCM spheroid build — custom K1-turgor, NO native md.mesh.

WHY THIS MODULE (the A5000 blocker it solves). The native-mesh capstone
(``cell/dcm_native_shell.py``) enforces per-cell turgor with
``md.mesh.conservation.Volume`` over a ``hoomd.mesh.Mesh`` that carries a SEPARATE
mesh triangle type PER CELL. HOOMD warns that many mesh triangle types perform
poorly on the GPU (shared-memory pressure) and, confirmed on the gbook RTX A5000,
the run STALLS at 0% GPU util by N≈60. This module is the GPU-friendly assembly:
it reproduces the same DCM physics WITHOUT any native md.mesh and WITHOUT per-cell
mesh types, so it runs GPU-resident.

THE STACK (all on the project's frozen L-M BAOAB-limit integrator):
  * geometry — ``dcm.icosphere_mesh`` shells at ``dcm._cluster_centers`` (reused).
  * TURGOR — ``dcm_gpu_forces.DcmTurgorForceGPU`` (the validated K1 mesh-pressure
    kernel): ONE custom force, per-cell enclosed volume by face-grouping on
    ``face_cell``, outward face-normal pressure. NO per-cell mesh types → no stall.
  * EDGE SPRINGS — ``md.bond.Harmonic`` with a SINGLE bond type (GPU-native,
    in-plane cortical elasticity). No per-cell bond types.
  * CELL-CELL — ``dcm_gpu_forces.DcmTentContactGPU`` (SimuCell3D bilinear tent),
    which distinguishes cells by the mutable ``cell_of_node`` array, NOT by
    per-cell particle/pair types — so there is a SINGLE membrane particle type.
  * SUBSTRATE — ``dcm_gpu_forces.DcmSubstrateForceGPU`` (capped-harmonic adhesive
    well at z0).
  * (optional) ACTIVE RIM TRACTION — ``dcm_gpu_forces.DcmActiveRimTractionGPU``.

A single membrane particle type means the BAOAB ``gamma`` dict has ONE entry
(plus the LOD freeze type, see ``dcm_gpu_lod.py``), and the contact / turgor
costs are O(nodes) device kernels rather than an O(n_cells²) pair-type matrix.

``device=`` is plumbed: a ``hoomd.device.GPU`` is created automatically when one
is available (gbook A5000), otherwise a CPU device (this Mac, dev/fallback). The
build runs ``run(0)`` then a short ``run(2000)`` finite gate before any long run.

Units SI.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.cell.dcm import icosphere_mesh, _cluster_centers
from ffn_sim.cell.dcm_gpu_forces import (
    DcmTurgorForceGPU,
    DcmTentContactGPU,
    DcmSubstrateForceGPU,
    DcmActiveRimTractionGPU,
)

_KT_310 = 4.28e-21  # J at 310 K
_UM = 1.0e6


# ---------------------------------------------------------------------------
# Resolved parameters (single membrane type; GPU-friendly bands)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedGpuDCM:
    """GPU-friendly DCM spheroid constants (single membrane type).

    The mechanical bands are the validated DCM literature bands (MCF7): R_cell
    7.5 µm, turgor ~133 Pa, cortical tension scale via the edge spring. The
    cell-cell contact is the SimuCell3D bilinear-tent (same bands as the native
    capstone). Adhesion energy densities are converted to per-node well depths at
    build via the per-node patch area.
    """

    R_cell: float = 7.5e-6          # m  MCF7 radius
    subdivisions: int = 1           # icosphere level (1 = 42 nodes, 2 = 162)
    turgor_dP0: float = 133.0       # Pa baseline osmotic turgor (Young-Laplace)
    K_vol: float = 1.0e3            # Pa osmotic bulk modulus (ΔP per ΔV/V)
    k_edge: float = 1.0e-3          # N/m membrane edge spring (cortical elasticity)
    gamma_node: float = 3.9e-10     # N·s/m per-node Stokes drag (membrane type)

    # cell-substrate / cell-cell adhesion energy densities (J/m², converted to
    # per-node well depths via area_per_node at build).
    W_cs_Jm2: float = 0.5e-3        # J/m² cell-substrate adhesion (spreading driver)

    # SimuCell3D bilinear-tent cell-cell contact (same bands as the native capstone)
    rep_strength: float = 1.0e8     # Pa/m repulsion stiffness ξ
    adh_strength: float = 1.0e8     # Pa/m adhesion stiffness ω
    c_adh: float = 5.0e-7           # m adhesion cutoff (tent peaks at c_adh/2)
    contact_force_cap: float = 5.0e-8   # N per-pair force cap (BAOAB guard)
    r_contact_factor: float = 1.0   # r_contact = factor · mean_edge

    spacing_factor: float = 2.3     # cell-center spacing = factor · R_cell
    cluster: str = "3d"             # "3d" spheroid ball | "2d" monolayer

    k_sub_Nm: float | None = None   # N/m substrate compliance (None = rigid dish)
    ligand_density: float = 1.0     # areal ligand density (cell-substrate adhesion)
    z_substrate: float = 0.0        # m substrate plane

    kT: float = _KT_310
    dt: float = 1.0e-9              # s (strong adhesion → small dt for BAOAB)
    seed: int = 7
    extras: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Device selection
# ---------------------------------------------------------------------------
def pick_device(device=None, *, seed_notice: int = 0):
    """Return a HOOMD device: the caller's, else GPU-if-available, else CPU.

    On the gbook A5000 a ``hoomd.device.GPU`` is created; on this Mac (no CUDA)
    the GPU constructor raises and we fall back to ``hoomd.device.CPU`` — so the
    SAME script runs in both places (GPU production / CPU dev-fallback).
    """
    if device is not None:
        return device
    try:
        return hoomd.device.GPU(notice_level=seed_notice)
    except Exception:
        return hoomd.device.CPU(notice_level=seed_notice)


# ---------------------------------------------------------------------------
# Snapshot — single membrane type, single bond type
# ---------------------------------------------------------------------------
def build_gpu_dcm_snapshot(p: ResolvedGpuDCM, n_cells: int):
    """Build a packed cluster of n_cells DCM shells with ONE membrane type.

    Returns a dict with the gsd Frame and all the per-cell bookkeeping the forces
    need: cell_of_node (mutable, used by the tent contact), ranges (per-cell node
    tag ranges), faces / face_cell (for the K1 turgor), nv/ne, mean_edge.
    """
    import gsd.hoomd

    verts0, edges, tris0 = icosphere_mesh(p.R_cell, p.subdivisions)
    nv = verts0.shape[0]
    ne = edges.shape[0]
    mean_edge = float(np.linalg.norm(
        verts0[edges[:, 0]] - verts0[edges[:, 1]], axis=1).mean())

    # surface-gapped start (no t=0 node overlap → no BAOAB blow-up); adhesion
    # pulls cells into contact gently over the settle.
    spacing = p.spacing_factor * p.R_cell
    centers = _cluster_centers(n_cells, spacing, p.z_substrate, p.R_cell,
                               mode=p.cluster)

    all_pos, bond_groups, cell_of_node = [], [], []
    face_groups, face_cell, ranges = [], [], []
    tag = 0
    for c in range(n_cells):
        ranges.append((tag, tag + nv))
        all_pos.append(verts0 + centers[c])
        cell_of_node.append(np.full(nv, c, dtype=np.int64))
        bond_groups.append(edges + tag)
        face_groups.append(tris0 + tag)
        face_cell.append(np.full(tris0.shape[0], c, dtype=np.int64))
        tag += nv

    pos = np.concatenate(all_pos, axis=0)
    cell_of_node = np.concatenate(cell_of_node)
    bonds = np.concatenate(bond_groups, axis=0)
    faces = np.concatenate(face_groups, axis=0)
    face_cell = np.concatenate(face_cell)

    box_edge = float(np.abs(pos).max() * 2.5 + 4.0 * p.R_cell)
    snap = gsd.hoomd.Frame()
    snap.particles.N = pos.shape[0]
    snap.particles.types = ["dcm_mem", "dcm_inert"]   # single active type + LOD freeze type
    snap.particles.position = pos
    snap.particles.typeid = np.zeros(pos.shape[0], dtype=np.uint32)  # all active
    snap.particles.mass = np.ones(pos.shape[0])
    snap.bonds.N = bonds.shape[0]
    snap.bonds.types = ["dcm_edge"]                   # SINGLE bond type
    snap.bonds.typeid = np.zeros(bonds.shape[0], dtype=np.uint32)
    snap.bonds.group = bonds.astype(np.uint32)
    snap.configuration.box = [box_edge, box_edge, box_edge, 0, 0, 0]
    return dict(snap=snap, cell_of_node=cell_of_node, ranges=ranges,
                faces=faces, face_cell=face_cell, nv=nv, ne=ne,
                mean_edge=mean_edge, centers=centers)


# ---------------------------------------------------------------------------
# Full GPU-friendly simulation assembly
# ---------------------------------------------------------------------------
def build_gpu_dcm_simulation(p: ResolvedGpuDCM, n_cells: int, *, device=None,
                             active: bool = False):
    """Assemble the GPU-friendly DCM spheroid on the BAOAB integrator.

    No native md.mesh, no per-cell mesh/particle/bond types. Returns a dict of
    handles (sim, cell_of_node, ranges, faces, face_cell, turgor, contact,
    substrate, traction|None, baoab, p, nv, ne, mean_edge).

    Args:
        active: if True, wire ``DcmActiveRimTractionGPU`` (basal-rim active
            traction). Default False (passive turgor + adhesion wetting only).
    """
    from ffn_sim.integrator.baoab import make_baoab_updater

    b = build_gpu_dcm_snapshot(p, n_cells)
    snap = b["snap"]
    nv, ne, mean_edge = b["nv"], b["ne"], b["mean_edge"]
    cell_of_node = b["cell_of_node"]
    ranges = b["ranges"]
    faces, face_cell = b["faces"], b["face_cell"]

    dev = pick_device(device)
    sim = hoomd.Simulation(device=dev, seed=p.seed)
    sim.create_state_from_snapshot(snap)

    ig = md.Integrator(dt=p.dt)

    # Edge springs (cortical elasticity); r0 = mean edge length. SINGLE bond type.
    bond = md.bond.Harmonic()
    bond.params["dcm_edge"] = dict(k=p.k_edge, r0=mean_edge)
    ig.forces.append(bond)

    # Per-node membrane patch area (adhesion-energy → per-node depths / contact A).
    area_per_node = 4.0 * np.pi * p.R_cell ** 2 / nv
    W_cs = p.W_cs_Jm2 * p.ligand_density * area_per_node    # J per node

    # TURGOR — K1 mesh-pressure (GPU-friendly: one custom force, face-grouped).
    V0 = (4.0 / 3.0) * np.pi * p.R_cell ** 3
    turgor = DcmTurgorForceGPU(faces=faces, face_cell=face_cell, n_cells=n_cells,
                               V0=V0, turgor_dP0=p.turgor_dP0, K_vol=p.K_vol)
    ig.forces.append(turgor)

    # CELL-CELL — bilinear-tent contact, distinguishes cells by cell_of_node.
    r_contact = p.r_contact_factor * mean_edge
    contact = DcmTentContactGPU(
        cell_of_node=cell_of_node, r_contact=r_contact, c_adh=p.c_adh,
        rep_strength=p.rep_strength, adh_strength=p.adh_strength,
        patch_area=area_per_node, force_cap=p.contact_force_cap)
    ig.forces.append(contact)

    # SUBSTRATE — adhesive capped-harmonic well at z0.
    substrate = DcmSubstrateForceGPU(
        z0=p.z_substrate, W_cs=W_cs, adh_range=p.R_cell, k_sub=p.k_sub_Nm)
    ig.forces.append(substrate)

    traction = None
    if active:
        active_mask = np.ones(n_cells, dtype=bool)
        int_mult = np.ones(n_cells, dtype=np.float64)
        traction = DcmActiveRimTractionGPU(
            cell_of_node=cell_of_node, ranges=ranges, active=active_mask,
            int_mult=int_mult, R_cell=p.R_cell, z0=p.z_substrate,
            f_act=1.2e-10, f_cap=6.0e-10, ramp_steps=4000, contact_band=0.5,
            neighbour_factor=2.6, max_neighbours=9, integrin_switch_gain=3.0,
            belt_factor=0.25)
        ig.forces.append(traction)

    sim.operations.integrator = ig
    sim.run(0)

    # BAOAB — gamma covers the single membrane type (+ the LOD inert freeze type).
    gamma = {"dcm_mem": p.gamma_node, "dcm_inert": p.gamma_node}
    action, updater = make_baoab_updater(
        kT=p.kT, gamma=gamma, dt=p.dt, seed=p.seed + 1)
    sim.operations.updaters.append(updater)

    return dict(
        sim=sim, cell_of_node=cell_of_node, ranges=ranges, faces=faces,
        face_cell=face_cell, n_cells=n_cells, nv=nv, ne=ne, mean_edge=mean_edge,
        V0=V0, turgor=turgor, contact=contact, substrate=substrate,
        traction=traction, baoab=action, gamma=gamma, p=p,
        area_per_node=area_per_node, centers=b["centers"])

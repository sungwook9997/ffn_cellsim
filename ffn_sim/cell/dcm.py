"""DCM tier — Deformable Cell Model for multicellular spheroid spreading.

PI pivot (2026-06-11): the platform's spheroid-scale work runs on a Deformable
Cell Model — the mid tier between the fine-grained single cell (every cytoskeletal
filament a particle; Phase 1) and the layer-2 coarse multicell (cell ≈ point).
Each DCM cell is a **triangulated elastic membrane shell** (icosphere nodes +
harmonic edge springs) enclosing a **turgor volume**, deformable under cell-cell
adhesion and substrate wetting — so a cluster of DCM cells can spread on a
substrate as a deformable spheroid. Literature band values are trusted directly
(no in-tier re-derivation) for the graduation-report timeline; the fine-grained
single-cell mechanics (cortical tension, turgor, FA traction) are the calibration
source for the shell constants.

Components (all on the project's L-M BAOAB-limit integrator):
  * icosphere membrane mesh (reuses `cortex.surface_manifold.SurfaceManifold.icosphere`)
    + harmonic edge springs (in-plane cortical elasticity).
  * `DcmTurgorForce` — per-cell enclosed-volume osmotic pressure (one Custom force
    over all cells; Young-Laplace tension emerges, balancing ΔP·R/2).
  * `DcmSubstrateForce` — adhesive substrate floor at z=0 (repulsive contact +
    adhesive well → wetting/spreading).
  * cell-cell interaction via HOOMD LJ on per-cell membrane types: same-cell =
    repulsive WCA (excluded volume), different-cell = attractive well (cadherin-
    scale adhesion).

Units SI. Default constants are MCF7 literature bands (R_cell 7.5 µm, turgor
133 Pa, cortical tension ~0.5 mN/m).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md

from ffn_sim.common.surface_manifold import SurfaceManifold

_KT_310 = 4.28e-21  # J at 310 K


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------
def icosphere_mesh(R: float, subdivisions: int = 1):
    """Return (verts (n,3) [m], edges (m,2) int, tris (k,3) int) of a sphere.

    Triangles are outward-wound (used for the exact enclosed-volume / face-normal
    pressure in :class:`DcmTurgorForce`).
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


# ---------------------------------------------------------------------------
# Per-cell turgor (enclosed-volume osmotic pressure) — one Custom force, all cells
# ---------------------------------------------------------------------------
class DcmTurgorForce(md.force.Custom):
    """Per-cell osmotic turgor as EXACT triangulated pressure (face-normal).

    Each cell's enclosed volume is the exact divergence-theorem volume of its
    triangulated shell, V_c = (1/6) Σ_faces (v0 · (v1×v2)). The osmotic pressure

        ΔP_c = turgor_dP0 + K_vol · (V0 − V_c)/V0

    is applied as an OUTWARD FACE-NORMAL force on each triangle, F_face =
    ΔP_c · (area·n̂), distributed equally to its 3 nodes. Unlike a radial-from-
    centroid force, this correctly conserves the true (non-spherical) volume —
    so when a cell flattens on the substrate, the pressure pushes its rim OUT
    (lateral spreading), the physically essential wetting behaviour. Faces are
    global node-index triplets (tag-ordered); `face_cell[k]` = the cell of face k.
    """

    def __init__(self, *, faces: np.ndarray, face_cell: np.ndarray, n_cells: int,
                 V0: float, turgor_dP0: float, K_vol: float) -> None:
        super().__init__(aniso=False)
        self.faces = np.asarray(faces, dtype=np.int64)          # (F,3) global idx
        self.face_cell = np.asarray(face_cell, dtype=np.int64)  # (F,)
        self.n_cells = int(n_cells)
        self.V0 = float(V0)
        self.turgor_dP0 = float(turgor_dP0)
        self.K_vol = float(K_vol)

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        n = pos.shape[0]
        perm = np.argsort(tag)                 # local-row order → tag(global) order
        pos_g = pos[perm]                      # pos_g[global_idx] = position
        f = self.faces
        v0, v1, v2 = pos_g[f[:, 0]], pos_g[f[:, 1]], pos_g[f[:, 2]]
        cross = np.cross(v1 - v0, v2 - v0)     # 2·area·n̂ (outward)
        # Per-cell exact volume (divergence theorem): V = (1/6) Σ v0·(v1×v2).
        vol_contrib = np.einsum("ij,ij->i", v0, cross) / 6.0
        Vc = np.zeros(self.n_cells)
        np.add.at(Vc, self.face_cell, vol_contrib)
        dP_cell = self.turgor_dP0 + self.K_vol * (self.V0 - Vc) / self.V0
        dP_face = dP_cell[self.face_cell]
        # Face-normal force = ΔP · area · n̂ = ΔP · (cross/2); split to 3 nodes.
        f_face = (dP_face[:, None] * cross / 2.0) / 3.0
        F_g = np.zeros_like(pos_g)
        for col in range(3):
            np.add.at(F_g, f[:, col], f_face)
        F = np.empty_like(pos)
        F[perm] = F_g                          # scatter back to local-row order
        with self.cpu_local_force_arrays as arr:
            arr.force[:] = F
            arr.potential_energy[:] = np.zeros(n, dtype=np.float64)


# ---------------------------------------------------------------------------
# Adhesive substrate floor (wetting → spreading)
# ---------------------------------------------------------------------------
class DcmSubstrateForce(md.force.Custom):
    """Adhesive substrate at z = z0 as a CAPPED harmonic well (BAOAB-stable).

    A node within ``range`` of z0 feels a harmonic restoring force toward z0
    (depth W_cs over the well), so it both adheres (pulled down from above) and
    cannot penetrate (pushed up from below) — no stiff-floor discontinuity. The
    well stiffness ``k_well = 2·W_cs / range²`` satisfies the overdamped CFL
    (dt ≪ γ/k_well) where a stiff steric floor would not. Beyond ``range`` the
    force is capped constant (= k_well·range) so a deeply displaced node is pushed
    back without an exploding force; a node above the reach feels nothing (not yet
    in contact). Spreading: bottom nodes adhere at z0; turgor (volume conserved)
    flattens the cell → footprint grows → more nodes enter contact.

        |dz| ≤ range : F_z = −k_well·dz              (harmonic well, toward z0)
        dz < −range  : F_z = +k_well·range            (capped steric push-up)
        dz > +range  : F_z = 0                         (out of adhesion reach)
    """

    def __init__(self, *, z0: float, W_cs: float, adh_range: float,
                 k_sub: float | None = None) -> None:
        super().__init__(aniso=False)
        self.z0 = float(z0)
        self.W_cs = float(W_cs)                       # J adhesion energy / node
        self.rng = float(adh_range)
        k_well = 2.0 * self.W_cs / (self.rng ** 2)    # N/m adhesion-well stiffness
        # Substrate COMPLIANCE: a finite-stiffness substrate (k_sub) softens the
        # felt restoring stiffness via the series combination (a node pressing
        # down also displaces the soft ground) → less traction → less spreading
        # on soft substrates = the durotaxis / rigidity-sensing knob. k_sub=None
        # ⇒ rigid dish (glass; the PI's prep) ⇒ k_eff = k_well (unchanged).
        self.k_well = k_well if k_sub is None else (k_well * k_sub) / (k_well + k_sub)

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            pos = np.asarray(snap.particles.position).copy()
        dz = pos[:, 2] - self.z0
        F = np.zeros_like(pos)
        U = np.zeros(pos.shape[0], dtype=np.float64)
        within = np.abs(dz) <= self.rng
        F[within, 2] = -self.k_well * dz[within]
        U[within] = 0.5 * self.k_well * dz[within] ** 2 - self.W_cs
        deep = dz < -self.rng
        F[deep, 2] = self.k_well * self.rng           # capped push-up
        U[deep] = self.k_well * self.rng * (-(dz[deep]) - 0.5 * self.rng) - self.W_cs
        with self.cpu_local_force_arrays as arr:
            arr.force[:] = F
            arr.potential_energy[:] = U


# ---------------------------------------------------------------------------
# Spheroid assembly
# ---------------------------------------------------------------------------
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


def build_dcm_spheroid_snapshot(p: ResolvedDCM, n_cells: int):
    """Build a packed cluster of n_cells DCM membrane shells on the substrate.
    Returns (gsd Frame, cell_of_tag (N,), per-cell node ranges, mean_edge_len)."""
    import gsd.hoomd

    verts0, edges, tris0 = icosphere_mesh(p.R_cell, p.subdivisions)
    nv = verts0.shape[0]
    ne = edges.shape[0]
    mean_edge = float(np.linalg.norm(
        verts0[edges[:, 0]] - verts0[edges[:, 1]], axis=1).mean())

    # Start cells with a surface gap (~0.3 R) so there is NO initial node overlap
    # (strong adhesion + t=0 overlap → BAOAB blow-up); adhesion then pulls them
    # into contact gently over the settle.
    spacing = 2.3 * p.R_cell
    centers = _cluster_centers(n_cells, spacing, p.z_substrate, p.R_cell,
                               mode=p.cluster)

    all_pos = []
    all_typeid = []
    bond_groups = []
    cell_of_tag = []
    face_groups = []
    face_cell = []
    ranges = []
    types = [f"dcm_m{c}" for c in range(n_cells)]
    tag = 0
    for c in range(n_cells):
        ranges.append((tag, tag + nv))
        all_pos.append(verts0 + centers[c])
        all_typeid.append(np.full(nv, c, dtype=np.uint32))
        cell_of_tag.append(np.full(nv, c, dtype=np.int64))
        bond_groups.append(edges + tag)
        face_groups.append(tris0 + tag)
        face_cell.append(np.full(tris0.shape[0], c, dtype=np.int64))
        tag += nv

    pos = np.concatenate(all_pos, axis=0)
    typeid = np.concatenate(all_typeid)
    cell_of_tag = np.concatenate(cell_of_tag)
    bonds = np.concatenate(bond_groups, axis=0)
    faces = np.concatenate(face_groups, axis=0)
    face_cell = np.concatenate(face_cell)

    box_edge = float(np.abs(pos).max() * 2.5 + 4.0 * p.R_cell)
    snap = gsd.hoomd.Frame()
    snap.particles.N = pos.shape[0]
    snap.particles.types = types
    snap.particles.position = pos
    snap.particles.typeid = typeid
    snap.particles.mass = np.ones(pos.shape[0])
    snap.bonds.N = bonds.shape[0]
    snap.bonds.types = ["dcm_edge"]
    snap.bonds.typeid = np.zeros(bonds.shape[0], dtype=np.uint32)
    snap.bonds.group = bonds.astype(np.uint32)
    snap.configuration.box = [box_edge, box_edge, box_edge, 0, 0, 0]
    return snap, cell_of_tag, ranges, mean_edge, nv, ne, faces, face_cell


def build_dcm_simulation(p: ResolvedDCM, n_cells: int, *, device=None):
    """Assemble the full DCM spheroid simulation on the BAOAB integrator.
    Returns dict of handles (sim, cell_of_tag, ranges, n_cells, nv, ...)."""
    from ffn_sim.integrator.baoab import make_baoab_updater

    (snap, cell_of_tag, ranges, mean_edge, nv, ne, faces, face_cell
     ) = build_dcm_spheroid_snapshot(p, n_cells)
    dev = device or hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=p.seed)
    sim.create_state_from_snapshot(snap)

    ig = md.Integrator(dt=p.dt)

    # Edge springs (cortical elasticity); r0 = mean edge length (turgor sets tension).
    bond = md.bond.Harmonic()
    bond.params["dcm_edge"] = dict(k=p.k_edge, r0=mean_edge)
    ig.forces.append(bond)

    # Per-node membrane patch area (for converting adhesion energy densities).
    area_per_node = 4.0 * np.pi * p.R_cell ** 2 / nv
    W_cc = p.W_cc_Jm2 * area_per_node      # J per contacting node pair (cell-cell)
    # Cell-substrate adhesion energy scales with areal LIGAND DENSITY (Gallant/
    # García; linear regime) — the coating-density / ECM-condition knob.
    W_cs = p.W_cs_Jm2 * p.ligand_density * area_per_node   # J per node

    # Cell-cell LJ: same-cell repulsive WCA (excluded vol), diff-cell attractive
    # (cadherin-scale adhesion). sigma ~ node spacing; r_cut adhesive tail.
    types = [f"dcm_m{c}" for c in range(n_cells)]
    sigma = 0.6 * mean_edge
    nlist = md.nlist.Cell(buffer=0.4 * mean_edge)
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
    rc_rep = 2 ** (1 / 6) * sigma
    rc_adh = 2.5 * sigma
    for i in range(n_cells):
        for j in range(i, n_cells):
            ti, tj = types[i], types[j]
            if i == j:
                lj.params[(ti, tj)] = dict(epsilon=p.eps_rep_kT * p.kT, sigma=sigma)
                lj.r_cut[(ti, tj)] = rc_rep
            else:
                lj.params[(ti, tj)] = dict(epsilon=W_cc, sigma=sigma)
                lj.r_cut[(ti, tj)] = rc_adh
    ig.forces.append(lj)

    # Turgor (per-cell osmotic pressure) as EXACT triangulated face-normal pressure.
    V0 = (4.0 / 3.0) * np.pi * p.R_cell ** 3
    turgor = DcmTurgorForce(faces=faces, face_cell=face_cell, n_cells=n_cells,
                            V0=V0, turgor_dP0=p.turgor_dP0, K_vol=p.K_vol)
    ig.forces.append(turgor)

    # Adhesive substrate (wetting → spreading). Reach = R_cell so a whole bottom
    # hemisphere wets (drives flattening), capped-harmonic well = BAOAB-stable.
    substrate = DcmSubstrateForce(
        z0=p.z_substrate, W_cs=W_cs, adh_range=p.R_cell, k_sub=p.k_sub_Nm)
    ig.forces.append(substrate)

    sim.operations.integrator = ig
    ig.forces  # ensure attached
    sim.run(0)

    gamma = {t: p.gamma_node for t in types}
    action, updater = make_baoab_updater(kT=p.kT, gamma=gamma, dt=p.dt, seed=p.seed + 1)
    sim.operations.updaters.append(updater)

    return {
        "sim": sim, "cell_of_tag": cell_of_tag, "ranges": ranges,
        "n_cells": n_cells, "nv": nv, "ne": ne, "mean_edge": mean_edge,
        "V0": V0, "turgor": turgor, "substrate": substrate, "p": p,
    }

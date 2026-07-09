"""De-cohesion spread on the Warp DCM engine (Phase C migration · M1: substrate).

Runs the de-cohesion cleanball spread on the GPU-resident Warp hybrid loop instead of
the per-step `cpu_local_snapshot` HOOMD driver — the 10–36× engine the port delivered.
This M1 wires the **substrate** drivers (z-well + in-plane wetting, the parity-verified
G2 kernels) into the multicell device loop alongside turgor + cortex edges + node-node
cohesion + node-face contact; the per-cell **lamellipodium** (M2) and **junction-switch**
(M3) follow.

Faithful to the landmine-register clean re-run prescription (`DCM_LANDMINE_REGISTER_2026-06-21`):
NO body-force proxy (settle deleted), the **conservative substrate wetting** is the spreading
driver, top-down silhouette A/A0 is co-tracked with maxZ + V/V0 + COM drift (never A/A0 alone),
and a divergence guard truncates on non-finite.

    python -m ffn_sim.dcm.dcm_warp_decohesion --device cuda:0 --n-cells 100 --steps 60000
"""

from __future__ import annotations

import argparse
import os
import time

import numpy as np

import warp as wp

from ffn_sim.dcm.geometry import icosphere_mesh, ResolvedDCM
from ffn_sim.dcm.dcm_warp_hybrid_multicell import (
    build_multicell, _dp_from_vol, _dp_from_vol_pc, _dp_from_vol_osm, osmotic_relax_kernel,
    _zero_vec, _edges_from_faces)
from ffn_sim.dcm.dcm_neighbor_warp import face_contact_count_kernel
from ffn_sim.dcm.dcm_warp_hybrid import _bond_accumulate, _bd_step
from ffn_sim.dcm.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel
from ffn_sim.dcm.dcm_virial_stress import node_virial_stress
from ffn_sim.dcm.dcm_cohesion_warp import dcm_cohesion_kernel
from ffn_sim.dcm.dcm_contact_warp import node_face_contact_kernel
from ffn_sim.dcm.dcm_interfacial_tension_warp import differential_surface_tension_kernel, douezan_spreading
from ffn_sim.dcm.dcm_aggregate_tension_warp import aggregate_laplace_kernel, aggregate_centroid_radius
from ffn_sim.dcm.confluent_init_prototype import build_confluent as _build_confluent_geom
from ffn_sim.dcm.dcm_contact_conservative_warp import (
    contact_grid_conservative_kernel, derive_adhesion_stiffness)
from ffn_sim.dcm.dcm_substrate_warp import (
    dcm_substrate_well_accum_kernel, dcm_wetting_scatter_kernel, dcm_wetting_cap_add_kernel,
    dcm_wetting_scatter_integrin_kernel, dcm_ubottom_well_kernel)
from ffn_sim.dcm.dcm_neighbor_warp import (
    pos_to_f32, face_centroids_f32, cohesion_grid_kernel, contact_grid_kernel,
    cohesion_grid_cad_kernel, contact_grid_cad_kernel, penetration_depth_kernel,
    edge_midpoints_f32, edge_edge_contact_kernel, cadherin_bond_force_kernel, gravity_body_force_kernel,
    ecm_clutch_force_kernel, cell_centroid_accum_kernel, nucleus_force_kernel,
    surface_tension_kernel, polarized_surface_tension_kernel, face_area_accum_kernel, global_area_force_kernel,
    edge_neighbor_sum_kernel, umbrella_kernel, bending_apply_kernel, scale_per_cell_kernel,
    gather_lead_pos, lamellipodium_tether_multicell)
from ffn_sim.dcm.dcm_cadherin_host import CadherinBondHost, CadherinParams
from ffn_sim.dcm.dcm_ecm_clutch_host import EcmClutchHost, EcmClutchParams
from ffn_sim.dcm.dcm_filopodia_host import FilopodiaHost, FilopodiaParams
from ffn_sim.dcm.dcm_filopodia_warp import (
    filopodia_tip_face_force_kernel, filopodia_tip_plane_force_kernel)
from ffn_sim.dcm.dcm_division_host import DivisionHost, DivisionParams
from ffn_sim.dcm.dcm_necrosis_host import NecrosisHost, NecrosisParams
from ffn_sim.dcm.dcm_lamellipodium_host import LamellipodiumHost, LamelParams
from ffn_sim.dcm.dcm_junction_switch_host import JunctionSwitchHost, JunctionParams
from ffn_sim.dcm.dcm_remesh import remesh_pass
from ffn_sim.dcm.dcm_cleave import cleave_cell
from ffn_sim.dcm.dcm_warp_implicit import (device_cg, ipc_newton_step, _vaxpy_active,
                                           _vaxpy_active_capped, _vaxpy, _vcopy)
from ffn_sim.dcm.dcm_contact_implicit_warp import (
    nearest_face_ipc_kernel, make_contact_hess_apply, ccd_alpha, project_contacts)
from ffn_sim.dcm.dcm_ipc_energy import (ipc_barrier_energy_kernel, edge_energy_kernel,
    turgor_energy_pc_kernel, turgor_dp0eff_kernel, inertial_energy_kernel, soft_linear_energy_kernel)
from ffn_sim.dcm.dcm_warp_frozen import (
    FrozenNeighborCache, build_diagA, make_diag_precond)

wp.init()


def _lattice_nearest(lattice: np.ndarray, n_cells: int) -> np.ndarray:
    """The ``n_cells`` lattice points nearest the origin (a rounded spherical cluster)."""
    return lattice[np.argsort(np.linalg.norm(lattice, axis=1))[:n_cells]]


def _spherical_centers(n_cells: int, R: float, gap: float, mode: str = "fcc",
                       seed: int = 5) -> np.ndarray:
    """``n_cells`` cell centres in a rounded SPHERICAL cluster (centre spacing ≈ ``gap·R``).

    A spherical cluster (not a cube) is required to AGGREGATE: surface cells feel a net INWARD
    cohesive pull (curvature) so the cluster rounds + compacts during settle. ``mode`` (D10):
      * ``cubic`` — simple-cubic lattice nearest-centre (legacy; cubic-axis anisotropy artifact).
      * ``fcc`` — face-centred-cubic close packing (nearest-neighbour = a/√2), isotropic + denser;
        the realistic default for a 3D tissue pack.
      * ``voronoi`` — Lloyd-relaxed centroidal Voronoi: random points in a ball relaxed (k-means
        over a dense ball sample) to even, isotropic spacing — the most tissue-like (no lattice
        order at all). Scaled so the mean nearest-neighbour spacing = gap·R.
    """
    side = int(np.ceil((2.2 * n_cells) ** (1.0 / 3.0))) + 3
    g = np.arange(-side, side + 1)
    if mode == "cubic":
        X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
        pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1).astype(float)
        centers = _lattice_nearest(pts, n_cells) * (gap * R)
    elif mode == "fcc":
        X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
        pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
        pts = pts[(pts.sum(axis=1) % 2) == 0].astype(float)     # FCC sublattice
        centers = _lattice_nearest(pts, n_cells) * (gap * R / np.sqrt(2.0))   # nn = a/√2 → gap·R
    elif mode == "voronoi":
        rng = np.random.default_rng(seed)
        # seed n_cells points in a unit ball, Lloyd-relax against a dense ball sample (CVT)
        def in_ball(m):
            q = rng.normal(size=(m, 3)); q /= np.linalg.norm(q, axis=1, keepdims=True)
            return q * (rng.random(m) ** (1.0 / 3.0))[:, None]
        c = in_ball(n_cells)
        sample = in_ball(max(4000, 60 * n_cells))
        for _ in range(25):                                     # Lloyd iterations
            lbl = np.argmin(((sample[:, None, :] - c[None, :, :]) ** 2).sum(-1), axis=1)
            for k in range(n_cells):
                m = lbl == k
                if m.any():
                    c[k] = sample[m].mean(0)
        # scale so the mean nearest-neighbour distance = gap·R
        d = np.sort(np.linalg.norm(c[:, None, :] - c[None, :, :], axis=-1), axis=1)[:, 1]
        centers = c * (gap * R / max(float(d.mean()), 1e-12))
    elif mode == "sphere":
        # ROUND + DENSE (fixes FCC-faceting AND CVT-looseness): random-close-PACK inside a BALL.
        # The ball boundary gives a smooth spherical envelope (no FCC cuboctahedral facets); the
        # jamming relaxation pushes every overlapping pair to the contact distance d (no CVT
        # isolated cells). → round AND all-touching, the proper aggregate init.
        from scipy.spatial import cKDTree
        rng = np.random.default_rng(seed)
        d = gap * R
        Rb = d * (n_cells / (8.0 * 0.64)) ** (1.0 / 3.0)     # ball radius at RCP density ~0.64
        u = rng.normal(size=(n_cells, 3)); u /= np.linalg.norm(u, axis=1, keepdims=True)
        c = u * (rng.random(n_cells) ** (1.0 / 3.0))[:, None] * Rb
        for _ in range(400):                                  # jam: separate overlaps, keep in ball
            pairs = cKDTree(c).query_pairs(d, output_type="ndarray")
            disp = np.zeros_like(c)
            if len(pairs):
                v = c[pairs[:, 0]] - c[pairs[:, 1]]
                nrm = np.linalg.norm(v, axis=1, keepdims=True) + 1e-12
                push = 0.5 * (d - nrm) * v / nrm
                np.add.at(disp, pairs[:, 0], push); np.add.at(disp, pairs[:, 1], -push)
            c = c + disp
            rr = np.linalg.norm(c, axis=1, keepdims=True)     # round boundary: project back into ball
            c = np.where(rr > Rb, c * Rb / np.maximum(rr, 1e-30), c)
        centers = c
    else:
        raise ValueError(f"unknown builder mode {mode!r} (cubic|fcc|voronoi|sphere)")
    return centers - centers.mean(0)


PARK_POS = np.array([1.0e-2, 1.0e-2, 1.0e-2])   # dormant pool node home (≫ cell scale; well=0)


def build_cleanball_on_substrate(n_cells: int, subdiv: int, R: float, z0: float = 0.0,
                                 gap: float = 2.05, pool_factor: float = 0.0,
                                 n_parked_cells: int = 0, builder: str = "fcc",
                                 cleave_pool: bool = False):
    """SPHERICAL cluster of cells RESTING on the substrate plane z0 (lowest node at z0).

    Cells are icospheres placed at :func:`_spherical_centers` (a rounded cluster, NOT the
    cubic lattice of ``build_multicell`` — see that helper for why a cube won't aggregate),
    then the whole ball is dropped so its lowest node sits at z0 (basal nodes engage the
    substrate well/wetting). ``gap`` is the centre spacing in R; the default **2.05** packs
    the shells essentially TOUCHING (≈0.05R apart) so the cluster STARTS as a compact,
    tissue-like aggregate — not the sparse 0.2R-gap lattice that just persists through settle
    (weak surface cohesion can't beat the stiff turgor that inflates each cell ~14%, so a
    loose pack never compacts; starting them in contact lets turgor + contact press the
    cells into flattened junctions = a genuine dense aggregate). gap<2.0 would hard-overlap
    and spike the contact repulsion, so 2.05 is the tight-but-safe floor with the soft-start."""
    verts1, edges1, tris1 = icosphere_mesh(R, subdiv)
    npc = verts1.shape[0]
    centers = _spherical_centers(n_cells, R, gap, mode=builder)
    # C7 CELL POOL: n_parked_cells extra UNDEFORMED icospheres parked far (cof=−1 → cohesion/
    # contact/measure/centroid skip them; turgor sees V=V0 → 0 force). A division ACTIVATES one
    # (sets its nodes' cof to the daughter id + moves it next to the mother) — faces/edges are
    # pre-allocated so no device topology resync, only pos/cof. Parked spread on a far lattice.
    n_total_cells = n_cells + n_parked_cells
    park_centers = [PARK_POS + np.array([3.0 * R * (k + 1), 0.0, 0.0]) for k in range(n_parked_cells)]
    all_centers = list(centers) + park_centers
    pos = np.concatenate([verts1 + c for c in all_centers], axis=0)
    # cleave_pool: CLEAVAGE carves daughters out of the mother shell — it never activates a whole
    # parked icosphere — so the parked cells supply DORMANT NODES + reserved ids ONLY. Emit faces/edges
    # for the ACTIVE cells alone (range(n_cells)); this drops ~n_parked·320 parked faces from every
    # per-step turgor/contact/volume kernel (the dominant per-step cost) → the cleave speed-up. The
    # mitotic path (cleave_pool=False) keeps the parked faces — it activates the parked icosphere whole.
    n_face_cells = n_cells if cleave_pool else n_total_cells
    faces = np.concatenate([tris1 + ci * npc for ci in range(n_face_cells)], axis=0)
    edges = np.concatenate([edges1 + ci * npc for ci in range(n_face_cells)], axis=0)
    cof = np.repeat(np.arange(n_total_cells), npc).astype(np.int64)
    cof[n_cells * npc:] = -1                      # parked cells start dormant (nodes cof=−1)
    face_cell = np.repeat(np.arange(n_face_cells), tris1.shape[0]).astype(np.int64)
    # rest the ACTIVE cluster on the dish (z0 = lowest active node); parked cells stay far
    act = np.zeros(pos.shape[0], dtype=bool); act[:n_cells * npc] = True
    pos[:, 2] += (z0 - pos[act, 2].min())

    # DORMANT NODE POOL for remeshing (A1): SPLIT activates a parked node, COLLAPSE returns one.
    n_pool = int(pool_factor * (n_cells * npc))
    if n_pool > 0:
        pool = np.tile(PARK_POS, (n_pool, 1))
        pos = np.concatenate([pos, pool], axis=0)
        cof = np.concatenate([cof, np.full(n_pool, -1, dtype=np.int64)])
    return pos, edges.astype(np.int64), faces.astype(np.int64), cof, face_cell, npc


def build_confluent_cluster(n_cells: int, subdiv: int, R: float, z0: float = 0.0,
                            pool_factor: float = 0.0, n_parked_cells: int = 0,
                            eps: float = 0.04, lloyd_iters: int = 8):
    """CONFLUENT space-filling builder (Path B): each cell's fixed-topology icosphere is radially
    warped into its Voronoi region (``confluent_init_prototype.build_confluent``) so the aggregate
    STARTS as a watertight, non-penetrating, already-faceted foam — the validated faceting init.
    Returns the same tuple as :func:`build_cleanball_on_substrate` (pos, edges, faces, cof,
    face_cell, npc) with constant ``npc`` preserved, so it drops into the run_decohesion dispatch.
    Use with ``--v0-from-init``-equivalent per-cell V0 (set automatically below)."""
    pos_a, faces_a, cof_a, npc, nf, _seeds, _R, _Rb = _build_confluent_geom(n_cells, subdiv, eps, R, lloyd_iters)
    _v, edges1, _t = icosphere_mesh(R, subdiv)                       # per-cell edge topology (constant)
    verts1, _e, tris1 = icosphere_mesh(R, subdiv)                    # for parked (undeformed) cells
    n_total = n_cells + n_parked_cells
    if n_parked_cells > 0:
        park = [PARK_POS + np.array([3.0 * R * (k + 1), 0.0, 0.0]) for k in range(n_parked_cells)]
        pos = np.concatenate([pos_a] + [verts1 + c for c in park], axis=0)
        faces = np.concatenate([faces_a] + [tris1 + (n_cells + k) * npc for k in range(n_parked_cells)], axis=0)
    else:
        pos, faces = pos_a, faces_a
    edges = np.concatenate([edges1 + ci * npc for ci in range(n_total)], axis=0)
    cof = np.repeat(np.arange(n_total), npc).astype(np.int64)
    cof[n_cells * npc:] = -1                                         # parked cells dormant
    face_cell = np.repeat(np.arange(n_total), nf).astype(np.int64)
    act = np.zeros(pos.shape[0], dtype=bool); act[:n_cells * npc] = True
    pos[:, 2] += (z0 - pos[act, 2].min())                           # rest the foam on the substrate
    n_pool = int(pool_factor * (n_cells * npc))
    if n_pool > 0:
        pos = np.concatenate([pos, np.tile(PARK_POS, (n_pool, 1))], axis=0)
        cof = np.concatenate([cof, np.full(n_pool, -1, dtype=np.int64)])
    return pos, edges.astype(np.int64), faces.astype(np.int64), cof, face_cell, npc


def _topdown_area_um2(pos_xy_um: np.ndarray) -> float:
    """Top-down silhouette area [µm²] = 2D convex hull of all live nodes (PI rule:
    NEVER basal contact area). Degenerate (<3 pts / collinear) → 0."""
    if pos_xy_um.shape[0] < 3:
        return 0.0
    try:
        from scipy.spatial import ConvexHull
        return float(ConvexHull(pos_xy_um).volume)   # 2D hull "volume" == area
    except Exception:
        return 0.0


def _cell_volumes(pos: np.ndarray, faces: np.ndarray, face_cell: np.ndarray, n_cells: int):
    v0 = pos[faces[:, 0]]; v1 = pos[faces[:, 1]]; v2 = pos[faces[:, 2]]
    contrib = np.einsum("ij,ij->i", v0, np.cross(v1 - v0, v2 - v0)) / 6.0
    V = np.zeros(n_cells)
    np.add.at(V, face_cell, contrib)
    return np.abs(V)


def run_decohesion(*, n_cells: int = 12, subdiv: int = 2, steps: int = 40000,
                   frames: int = 20, device: str = "cpu", dt: float = 8.0e-6,
                   k_vol: float = 7.73e5, rep_strength: float = 2.0e8,
                   adh_strength: float = 1.0e7, w_cs_jm2: float = 2.85e-3,
                   adh_range: float = 0.5e-6, k_floor: float = 1.0, gap: float = 2.05,
                   substrate_wetting: bool = True, use_substrate_well: bool = True,
                   ubottom: bool = False, ubottom_r_factor: float = 1.35, ubottom_k: float = 0.0,
                   force_cap: float = 5.0e-8, z0: float = 0.0, warmup: int = 1000,
                   settle_steps: int = 0, settle_frames: int = 0,
                   remesh_period: int = 0, pool_factor: float = 0.5,
                   edge_edge: bool = False, cfl_limit: float = 0.0, max_substeps: int = 16,
                   cadherin: bool = False, ecm_clutch: bool = False, cad_batch: int = 50,
                   cad_subcycle: bool = True, cad_micro_M: int = 8,
                   cad_bundle: float = 1.0, cad_contract: float = 0.0, cad_rbind: float = 0.0,
                   ecm_bundle: float = 1.0,
                   ecm_ligand: float = 1.0,
                   gravity: bool = False, delta_rho: float = 55.0, coupling: bool = False,
                   pen_cap: bool = True, pen_cap_frac: float = 1.0,
                   nucleus: bool = False, E_nuc: float = 3.0e3, ratio_lamin: float = 3.0,
                   knee_strain: float = 0.10, R_nuc_factor: float = 0.33,
                   surface_tension: bool = False, gamma_surf: float = 1.0e-4, k_area: float = 0.0,
                   diff_tension: bool = False, contact_tension_frac: float = -1.0,
                   aggregate_tension: bool = False, sigma_agg: float = 5.0e-3, agg_every: int = 25,
                   conservative_contact: bool = False, rep_over_adh: float = 4.0,
                   k_edge: float = -1.0,
                   polarize: bool = False, w_cs_polarize: float = 2.85e-3,
                   ipc_dhat_factor: float = 1.0,
                   division: bool = False, div_pool_factor: float = 1.0, div_rate: float = 0.04,
                   div_real_hours: float = 0.0, div_t_cycle_h: float = 24.0,
                   div_relax_steps: int = 20, div_relax_factor: float = 0.1,
                   osmotic: bool = False, osm_relax: float = 0.05, osm_batch: int = 50,
                   force_divide_step: int = -1, force_divide_cell: int = 0,
                   init_npz: str | None = None, v0_from_init: bool = False,
                   inset: float = 0.04, lloyd_iters: int = 8,
                   accel_real_hours: float = 0.0,
                   bending: bool = False, k_bend: float = 1.0e-5,
                   necrosis: bool = False, builder: str = "fcc",
                   integrator: str = "baoab", accel_dt: float | None = None, cg_maxiter: int = 80,
                   frozen_neighbors: bool = False, precond_diag: bool = False,
                   use_grid: bool = True, save_frames: str | None = None,
                   ipc: bool = False, ipc_eta: float = 0.9,
                   ipc_newton: bool = False, ipc_newton_max: int = 8, ipc_newton_tol: float = 1e-4,
                   project: bool = False, proj_omega: float = 0.7, proj_iter: int = 8,
                   proj_gap_factor: float = 1.0,
                   lamellipodium: bool = False, lamel_clutch: bool = False, filopodia: bool = False,
                   active_batch: int = 50, gpu_probe: bool = False, junction_switch: bool = False,
                   cleave: bool = False, lamel_all_cell: bool = False) -> dict:
    """Cleanball de-cohesion spread on the Warp loop with substrate drivers (M1) plus
    the optional per-cell lamellipodium crawl (M2, ``lamellipodium=True``).

    The lamellipodium is the host-managed advancing-anchor port (see
    :mod:`dcm_lamellipodium_host`): rim cells are detected once at build, then every
    ``batch_steps`` the front is ratcheted (SEED/ADVANCE) on the host and the device
    anchor + leading-node geometry arrays are refreshed; every step the device tether
    (:func:`lamellipodium_tether_multicell`) pulls each rim cell's leading basal node
    toward its own front. Wetting (M1) is the basal spread; the lamellipodium adds the
    active crawl on top — and because the rim cells are cohesively bonded to the cells
    above them, the upper (non-ECM-contacting) cells are dragged along (collective spread)."""
    p = ResolvedDCM(subdivisions=subdiv)
    if k_edge > 0.0:                                # cortex deformability override (faceting lever:
        p.k_edge = k_edge                          # a stiff cortex resists flattening → round cells)
    R = p.R_cell
    implicit = (integrator == "implicit")          # I-opt: linearly-implicit IMEX integration
    if implicit and not use_grid:
        # stiff_force_into builds cohesion/contact ONLY on the grid path; without it the implicit
        # operator would omit the dominant contact stiffness (review fix #4). Require the grid.
        raise ValueError("--integrator implicit requires the hash-grid path (do not combine with --no-grid)")
    if implicit and accel_dt:
        dt = accel_dt                              # implicit unlocks a larger (accuracy-bound) dt
    if implicit:
        print(f"  [integrator] IMPLICIT (IMEX linearly-implicit) dt={dt:.1e}  cg_maxiter={cg_maxiter} "
              f"— stiff operator (contact/turgor/edges/well/nucleus/bending) implicit, soft drivers explicit", flush=True)
    # A1 remesh and C7 division both draw dormant nodes from cof<0 and cannot co-run safely yet.
    # NB: this is NOT a one-line "−2 sentinel" fix — parked-cell node blocks double as cleave's
    # dormant ring-node supply (cleave_cell draws cof<0; no separate −1 pool is built when remesh
    # is off), and the mitotic path's fixed npc-block indexing is incompatible with remesh
    # relabelling regardless of sentinel. Design + the two viable approaches (index-range vs a
    # unified pool manager) in docs/v2_audit/DCM_DIVISION_REMESH_CORUN_DESIGN_2026-06-29.md.
    # Guard: division wins (it owns the pool).
    if division and remesh_period:
        print("  [warn] division + remesh both requested — disabling remesh this run "
              "(shared cof<0 pool; see DCM_DIVISION_REMESH_CORUN_DESIGN_2026-06-29.md).", flush=True)
        remesh_period = 0
    # FilopodiaHost caches self.cof / self.faces / self.fcell (host) PLUS device arrays built from the
    # mesh topology; a remesh SWAP/SPLIT/COLLAPSE reassigns faces/fcell/cof and changes the topology,
    # leaving the filopodia host's cached refs + device state STALE → a tip bonded to a parked node
    # yields a box-scale spurious force (same class as the cad/ecm remesh bug, but device-deep). Until
    # the filopodia host re-points + rebuilds on remesh, guard the combination (matches division above).
    # (Found by the 2026-06-23 night-deliverables verification — the cad/ecm fix missed the 5th/6th host.)
    if filopodia and remesh_period:
        print("  [warn] filopodia + remesh both requested — disabling remesh this run "
              "(FilopodiaHost caches cof/faces/fcell + device state; remesh would leave them stale).", flush=True)
        remesh_period = 0
    n_active = n_cells                              # requested live cells (C7 adds a parked pool)
    # mitotic-round reserves a full parked icosphere per potential daughter (div_pool_factor·n).
    # CLEAVE only needs reserved cell IDS + a dormant NODE pool, and it has NO parked faces (cleave_pool),
    # so a tighter reserve suffices: ~35% of n covers the ~14%/24h division fraction with 2.5× headroom
    # (a short pool just skips a division, graceful) while cutting the parked-NODE overhead from 100%→35%.
    n_parked = int((0.35 if cleave else div_pool_factor) * n_cells) if division else 0
    if init_npz:
        # RESTART from a saved aggregate (PI 2026-06-25 "aggregate 한걸로 스프레딩까지"): load the final
        # frame of a prior run as the initial state instead of building a fresh ball, and DROP it onto
        # the substrate floor (min active-node z → z0) so the basal cells engage the dish/ECM clutch.
        # Topology (faces) + cell map (cof) are reused; edges are re-derived; r0 takes the current edge
        # length (the aggregate is near-round, so the locked-in rest config error is small). The parked
        # pool (cof<0, if any) is loaded as-is (dormant; spreading runs typically pass division=False).
        dd = np.load(init_npz, allow_pickle=True)
        pos_a = np.ascontiguousarray(dd["frames"][-1].astype(np.float64))
        faces_a = dd["faces"].astype(np.int64)
        cof_a = dd["cof"].astype(np.int64)
        npc = int(icosphere_mesh(R, subdiv)[0].shape[0])
        fcell_a = cof_a[faces_a[:, 0]].astype(np.int64)
        edges_a = _edges_from_faces(faces_a).astype(np.int64)
        actm = cof_a >= 0
        pos_a[:, 2] += (z0 - float(pos_a[actm][:, 2].min()))   # rest the aggregate on the substrate
        n_cells = pos_a.shape[0] // npc
        n_active = int(np.unique(cof_a[actm]).size)
        n_parked = n_cells - n_active
        print(f"  [restart] loaded {n_active} active cells ({pos_a.shape[0]} nodes) from {init_npz} "
              f"→ dropped onto substrate z0={z0:.3g}", flush=True)
    elif builder == "confluent":
        pos_a, edges_a, faces_a, cof_a, fcell_a, npc = build_confluent_cluster(
            n_cells, subdiv, R, z0, pool_factor=(pool_factor if remesh_period else 0.0),
            n_parked_cells=n_parked, eps=inset, lloyd_iters=lloyd_iters)
        v0_from_init = True                         # confluent cells rest at their Voronoi volume, not the free sphere
        print(f"  [builder] confluent Voronoi-warp foam: {n_active} active + {n_parked} parked", flush=True)
        n_cells = n_active + n_parked
    else:
        pos_a, edges_a, faces_a, cof_a, fcell_a, npc = build_cleanball_on_substrate(
            n_cells, subdiv, R, z0, gap=gap, pool_factor=(pool_factor if remesh_period else 0.0),
            n_parked_cells=n_parked, builder=builder, cleave_pool=cleave)
        print(f"  [builder] {builder} pack: {n_active} active + {n_parked} parked, gap={gap}·R", flush=True)
        n_cells = n_active + n_parked               # per-cell arrays size to the FULL pool
    N = pos_a.shape[0]
    mean_edge = float(np.linalg.norm(pos_a[edges_a[:, 0]] - pos_a[edges_a[:, 1]], axis=1).mean())
    R0 = float(np.linalg.norm(icosphere_mesh(R, subdiv)[0], axis=1).mean())
    V0 = (4.0 / 3.0) * np.pi * R0 ** 3
    # PER-CELL rest volume (osmotic setpoint). All cells start at the full V0; CELL DIVISION halves
    # the mother+daughter setpoint (volume-conserving cytokinesis) and a time-consistent regrowth
    # ramp re-inflates each back to V0 over the cell cycle — turgor chases this per-cell setpoint, so
    # a freshly-born daughter inflates GRADUALLY instead of a full-size cold insert (the prior
    # outward-dump ejection). Host-mirrored (V0_cell) + device (V0_cell_d); kept in sync on division.
    if v0_from_init:
        # Physiological setpoint for a CONFLUENT cell = its actual resting (Voronoi-cell) volume in the
        # tissue, NOT the full free-sphere V0 (which over-inflates a space-filling cell — the N=400
        # confluent capstone hit V/V0 1.37 + interpenetration because every cell tried to reach the
        # free-sphere volume it cannot occupy). Auto-on for --builder confluent and --init-npz
        # --v0-from-init. Derive per-cell V0 from the built/loaded mesh via the
        # divergence-theorem volume (1/6·Σ(a×b)·c over its faces; winding-consistent icosphere topology).
        # Derived from the init geometry, NOT tuned. (physiological-baseline rule: confluent rest volume.)
        _a = pos_a[faces_a[:, 0]]; _b = pos_a[faces_a[:, 1]]; _c = pos_a[faces_a[:, 2]]
        _ftv = (np.cross(_a, _b) * _c).sum(axis=1) / 6.0
        V0_cell = np.zeros(n_cells, dtype=np.float64)
        _valid = fcell_a >= 0
        np.add.at(V0_cell, fcell_a[_valid], _ftv[_valid])
        V0_cell = np.abs(V0_cell)
        V0_cell[V0_cell < 1e-30] = V0                       # parked/empty cells → free-sphere fallback
        print(f"  [v0-from-init] per-cell V0 from confluent mesh: mean={V0_cell.mean():.3e}m³ "
              f"(free-sphere V0={V0:.3e}, ratio {V0_cell.mean()/V0:.3f}) — confluent rest-volume setpoint",
              flush=True)
    else:
        V0_cell = np.full(n_cells, V0, dtype=np.float64)
    area_per_node = 4.0 * np.pi * R0 ** 2 / npc
    # z-well depth from the adhesion energy density × node area (derived, not tuned)
    W_cs_well = w_cs_jm2 * area_per_node
    k_well = 2.0 * W_cs_well / (adh_range ** 2)

    # ── U-bottom ULA confinement geometry (opt-in via ubottom=True; independent of the flat
    # z-well). Computed ONCE from the just-built cluster: a hemispherical bowl whose floor sits at
    # z=z0 (the cluster's resting plane) and whose radius is ubottom_r_factor× the cluster radius,
    # so the dense pellet drops into the bowl bottom and is confined (non-adhesive, repulsive-only).
    ub_cx = ub_cy = ub_cz = ub_rwell = ub_kwall = 0.0
    if ubottom:
        act_mask = cof_a >= 0                       # active nodes only (parked/pool = −1)
        cl = pos_a[act_mask]
        ub_cx = float(cl[:, 0].mean())
        ub_cy = float(cl[:, 1].mean())
        # robust 3D cluster radius (99th pct of node distance from the cluster centroid)
        cen = cl.mean(axis=0)
        d3 = np.linalg.norm(cl - cen, axis=1)
        cluster_radius = float(np.percentile(d3, 99.0))
        ub_rwell = float(ubottom_r_factor) * cluster_radius
        # bowl centre z so the bowl FLOOR (z = cz − r_well) sits at z0 (= 0, cluster rest plane)
        ub_cz = z0 + ub_rwell
        # auto stiffness: reuse the rigid-dish floor scale (a rigid, repulsive wall) unless given
        ub_kwall = float(ubottom_k) if ubottom_k > 0.0 else float(k_floor)
        print(f"  [ubottom] bowl r_well={ub_rwell*1e6:.2f}um center=({ub_cx*1e6:.2f},"
              f"{ub_cy*1e6:.2f},{ub_cz*1e6:.2f})um k={ub_kwall:.3e} "
              f"(cluster_radius={cluster_radius*1e6:.2f}um, r_factor={ubottom_r_factor})", flush=True)

    # PHYSIOLOGICAL per-node drag (physiological-baseline HARD rule) — derived from the
    # MCF7 cytoplasm viscosity exactly as the HOOMD driver (dcm_gpu_build.py:576): the
    # cell's Stokes drag 6π·η·R distributed over its nodes. NOT the legacy ResolvedDCM
    # gamma_node=3.9e-10 (~5.7e5× too small → would force a brittle ~5e-8 dt + blow up at 8e-6).
    eta_cytoplasm_Pas = 65.9          # MCF7 cytoplasm (Dessard 2024)
    gamma_node = 6.0 * np.pi * eta_cytoplasm_Pas * R / npc
    inv_gamma = 1.0 / gamma_node
    # D7: net sedimentation body force (gravity − buoyancy). Δρ = ρ_cell − ρ_medium, DERIVED
    # (no magic number): ρ_cell≈1060 (MCF7; SimuCell3D Table 2 1045–1099 kg/m³), ρ_medium≈1005
    # (PBS/culture medium) → Δρ≈55 kg/m³. Per-node z-force = −Δρ·g·v_node, v_node = V0/npc.
    # Small (~1 pN/cell) but physiological — it rests the spheroid on the dish at its real weight
    # (replaces the removed 19000× wetting proxy, which was ~3×10⁵× this). Off ⇒ fz_node=0.
    v_node = V0 / npc
    fz_node = -(delta_rho * 9.80665 * v_node) if gravity else 0.0
    if gravity:
        print(f"  [gravity] Δρ={delta_rho:.0f}kg/m³  net f={-fz_node * npc * 1e12:.2f}pN/cell "
              f"(sedimentation toward dish; gravity−buoyancy, derived)", flush=True)
    c_rep = 0.30 * mean_edge
    c_adh = 0.80 * mean_edge
    # CONSERVATIVE bilinear-tent contact (the energy fix): adhesion stiffness DERIVED from the
    # measured adhesion energy density w_cs (adh = 4·w_cs/c_adh², grid-invariant); repulsion a
    # multiple of it (stiffer non-penetration). Units Pa/m (vs the old constant-force Pa).
    adh_cons = derive_adhesion_stiffness(w_cs_jm2, c_adh)
    rep_cons = rep_over_adh * adh_cons
    if conservative_contact:
        # The tent gives non-penetration + normal cohesion; it does NOT spread contacts. Whether the
        # faceting drive (differential γ, Maître/DAH) is ALSO needed is the open question — so let
        # --diff-tension COMPOSE with the tent rather than being force-disabled (SimuCell3D itself
        # uses uniform γ + tent + a confluent init; the diagnosis doc found uniform γ rounds at a
        # non-confluent gap). diff_tension stays whatever the caller requested (default False).
        _ten = "differential γ (Maître/DAH)" if diff_tension else "uniform γ"
        print(f"  [conservative-contact] bilinear tent: adh={adh_cons:.2e}Pa/m "
              f"(=4·w_cs/c_adh², w_cs={w_cs_jm2:.2e}J/m²) rep={rep_cons:.2e}Pa/m "
              f"(={rep_over_adh:.0f}×adh) c_adh={c_adh*1e6:.2f}um — node-node cohesion OFF, {_ten}",
              flush=True)
    r_contact = 0.30 * mean_edge

    pos_d = wp.array(np.ascontiguousarray(pos_a), dtype=wp.vec3d, device=device)
    cof_d = wp.array(cof_a.astype(np.int32), dtype=wp.int32, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    wbuf_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    Vc_d = wp.zeros(n_cells, dtype=wp.float64, device=device)
    dP_d = wp.zeros(n_cells, dtype=wp.float64, device=device)
    V0_cell_d = wp.array(V0_cell, dtype=wp.float64, device=device)   # per-cell osmotic setpoint (division/regrowth)
    faces_per_cell = int(faces_a.shape[0] // max(n_cells, 1))         # icosphere faces per cell (constant)
    contact_cnt_d = wp.zeros(n_cells, dtype=wp.int32, device=device)  # OSMOTIC: per-cell #faces apposed to other cells
    cad_d = wp.ones(n_cells, dtype=wp.float64, device=device)
    integrin_d = wp.ones(n_cells, dtype=wp.float64, device=device)
    # E2 nucleus: per-cell centroid reduction buffers + the bilinear chromatin/lamin stiffnesses
    csum_d = wp.zeros(n_cells, dtype=wp.vec3d, device=device)
    ccnt_d = wp.zeros(n_cells, dtype=wp.float64, device=device)
    # B4 surface tension / global area: per-cell area buffer + the initial per-cell area target A0
    acell_d = wp.zeros(n_cells, dtype=wp.float64, device=device)
    _fa0 = 0.5 * np.linalg.norm(np.cross(pos_a[faces_a[:, 1]] - pos_a[faces_a[:, 0]],
                                         pos_a[faces_a[:, 2]] - pos_a[faces_a[:, 0]]), axis=1)
    _A0c = np.zeros(n_cells); np.add.at(_A0c, fcell_a, _fa0)
    a0cell_d = wp.array(_A0c, dtype=wp.float64, device=device)
    # B5 bending: biharmonic via two umbrella-Laplacian passes (per-node neighbour buffers)
    lap_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    bilap_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    nsum_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    ncnt_d = wp.zeros(N, dtype=wp.float64, device=device)
    # C8 necrosis: per-cell turgor multiplier (necrotic core loses pressure regulation)
    turgor_mult_d = wp.ones(n_cells, dtype=wp.float64, device=device)
    R_nuc = R_nuc_factor * R
    d_knee_nuc = knee_strain * R_nuc
    k_chrom = 4.0 * np.pi * E_nuc * R_nuc / npc          # continuum bridge (H.9 eq B), N-invariant
    k_lamin = (ratio_lamin - 1.0) * k_chrom
    faces_d = wp.array(faces_a.astype(np.int32), dtype=wp.int32, device=device)
    fcell_d = wp.array(fcell_a.astype(np.int32), dtype=wp.int32, device=device)
    edges_d = wp.array(edges_a.astype(np.int32), dtype=wp.int32, device=device)
    r0_d = wp.array(np.linalg.norm(pos_a[edges_a[:, 0]] - pos_a[edges_a[:, 1]], axis=1),
                    dtype=wp.float64, device=device)
    n_edges = edges_a.shape[0]
    n_faces = faces_a.shape[0]

    # hash-grid neighbour list (M1.5): cohesion + node-face contact become O(N·k) instead
    # of O(N²)/O(N·M) — the speed lever at spheroid scale (grid == brute to machine-eps;
    # face query radius = c_adh + triangle circumradius ~0.7·l_max, l_max = 3·l_min).
    l_min = mean_edge / 2.9
    coh_q = float(c_adh)
    con_q = float(c_adh + 0.7 * (3.0 * l_min))
    # M1 IPC needs a LARGER query radius than the penalty's con_q so a node penetrating deep under
    # the strong bundle does not lose its entry face (the n100 partial-result cause #1). The IPC
    # self-test uses con_q + 3·mean_edge; mirror it here. The activation gap d̂ is also enlarged
    # (c_rep=0.30·me is too tight vs big accel-dt steps) so the barrier engages BEFORE a node is
    # buried in one jump. Both are numerical-correctness levers (not outcome-tuning).
    # repel_q query-radius enlargement REFUTED by n100 data: con_q+3·me gave pen_final 2.20 (WORSE than
    # con_q's 1.80) and 5× slower (61 min) — so entry-face-loss was NOT the bottleneck (the persistent
    # penetration is the FEASIBILIZATION equilibrium of already-inside nodes under the strong bundle,
    # which a larger search radius cannot touch). Reverted to con_q (the proven-best, fastest IPC config).
    # d̂ likewise reverted to c_rep (enlarging to 1·me diverged: cfl 55, bonds 3139→142, pen 4.0).
    ipc_repel_q = float(con_q)
    ipc_dhat = float(ipc_dhat_factor * c_rep)   # factor=1 -> IPC barrier; factor~0.01 -> barrier negligible = SimuCell3D-style implicit penalty-only (the A/B)
    grid_q = ipc_repel_q if ipc else con_q              # the FACE grid must be built >= the largest query radius
    node_f32 = wp.zeros(N, dtype=wp.vec3, device=device)
    cent_f32 = wp.zeros(n_faces, dtype=wp.vec3, device=device)
    pen_d = wp.zeros(N, dtype=wp.float64, device=device)   # interpenetration diagnostic
    # I-opt implicit: all-device CG scratch (only when --integrator implicit)
    cg_scratch = ({k: wp.zeros(N, dtype=wp.vec3d, device=device)
                   for k in ("r", "p", "Ap", "dx", "Fx", "Fp", "xp", "Fb")} if implicit else None)
    if implicit:
        cg_scratch["sca"] = wp.zeros(1, dtype=wp.float64, device=device)
    # I-opt #2 (frozen-neighbour cache) + #1 (analytic-diagonal Jacobi PCG): opt-in, default OFF →
    # the implicit path is byte-identical when both are False. The cache snapshots the cohesion +
    # contact candidate lists once per step (at xₙ, after the grids are built) so the CG matvec
    # re-uses them WITHOUT re-querying the hash-grid each iteration (the per-iter query is the
    # query-bound matvec cost). The diagonal is rebuilt once per step for the Jacobi M⁻¹=1/diagA.
    frozen_cache = (FrozenNeighborCache(N, device=device)
                    if (implicit and frozen_neighbors and use_grid) else None)
    diagA_d = (wp.zeros(N, dtype=wp.float64, device=device)
               if (implicit and precond_diag) else None)
    precond_apply = make_diag_precond(diagA_d, device=device) if diagA_d is not None else None
    if implicit and frozen_neighbors and use_grid:
        print(f"  [I-opt] frozen-neighbour cache ON — CG matvec reuses the xₙ cohesion/contact "
              f"candidate lists (no per-iter hash-grid query)", flush=True)
    if implicit and precond_diag:
        print(f"  [I-opt] analytic-diagonal Jacobi preconditioner ON (z=r/diagA each CG iter; "
              f"diagA = γ/dt + k_edge + rep·area + turgor-vol)", flush=True)
    # M1 IPC node-face contact (opt-in): barrier force RHS + analytic barrier Hessian (frozen at xₙ)
    # in the implicit operator + a CCD-filtered step. d̂ = c_rep (derived activation gap), κ = rep.
    ipc_cn_k = wp.zeros(N, dtype=wp.float64, device=device) if ipc else None
    ipc_cn_nrm = wp.zeros(N, dtype=wp.vec3d, device=device) if ipc else None
    ipc_t = wp.zeros(N, dtype=wp.float64, device=device) if ipc else None
    ipc_hess = make_contact_hess_apply(ipc_cn_k, ipc_cn_nrm, device=device) if ipc else None
    # #1 Step 3: projected-Newton IPC buffers (opt-in --ipc-newton; requires --ipc). xn = frozen
    # step-start copy; fsoft = the LAGGED soft-driver force (F_all(xn) − stiff{turgor,edges,barrier});
    # dp0eff = per-cell osmotic-excess for the turgor line-search energy.
    if ipc_newton and ipc and implicit:
        xn_d = wp.zeros(N, dtype=wp.vec3d, device=device)
        fsoft_d = wp.zeros(N, dtype=wp.vec3d, device=device)
        nt_Fb = wp.zeros(N, dtype=wp.vec3d, device=device)
        dp0eff_d = wp.zeros(n_cells, dtype=wp.float64, device=device)
        nt_e_d = wp.zeros(1, dtype=wp.float64, device=device)
        _nt_a = [0.0]                                    # a=γ/dt for the current step (nt_energy inertial)
        for _k in ("nt_Ft", "nt_rhs", "nt_xtr"):
            cg_scratch[_k] = wp.zeros(N, dtype=wp.vec3d, device=device)
        print(f"  [ipc-newton] projected-Newton IPC ON (max_newton={ipc_newton_max}, "
              f"tol={ipc_newton_tol:.0e}) — stiff{{turgor,edges,barrier}} re-linearised + energy "
              f"line-search; all other drivers lagged-soft at xₙ", flush=True)
    # M1 PROJECTION (opt-in): post-step geometric non-penetration constraint (SimuCell3D hard-constraint
    # method). Buffer for the per-node positional correction; applied after the integrator's pos update.
    proj_dpos = wp.zeros(N, dtype=wp.vec3d, device=device) if project else None
    node_grid = wp.HashGrid(48, 48, 48, device=device) if use_grid else None
    face_grid = wp.HashGrid(48, 48, 48, device=device) if use_grid else None

    # A2 edge-edge contact buffers (per-edge owner cell + midpoint grid). ee_q bounds the
    # midpoint separation of two edges that can pass within c_rep (c_rep + the two half-lengths,
    # half-length ≤ l_max/2 = 1.5·l_min).
    edge_cell_d = wp.array(cof_a[edges_a[:, 0]].astype(np.int32), dtype=wp.int32, device=device)
    emid_f32 = wp.zeros(n_edges, dtype=wp.vec3, device=device)
    edge_grid = wp.HashGrid(48, 48, 48, device=device) if (use_grid and edge_edge) else None
    ee_q = float(c_rep + 3.0 * l_min)

    # ── UNIFIED time-acceleration factor S (quasi-static timescale separation; PI 2026-06-24).
    # The integrator dt (accel_dt when set, else dt) is REAL physical time anchored to MCF7
    # cytoplasm viscosity — the MECHANICAL forces (turgor, contact, cortex/edge springs, cadherin
    # bond FORCE) stay at the faithful dt and are NEVER scaled. Only the SLOW BIOLOGICAL RATES
    # (protrusion velocities, nucleation/binding probabilities, cadherin k_on/k_off, division) are
    # multiplied by ONE factor S, so the RELATIVE timing of every slow process is preserved and the
    # run REPRESENTS  t_real = S · dt · steps  of real biology. S = (accel_real_hours·3600)/(dt·steps);
    # accel_real_hours=0 ⇒ S=1 ⇒ native rates ⇒ byte-identical to the un-accelerated run (ADDITIVE).
    _dt_accel = accel_dt if accel_dt else dt
    S_accel = ((accel_real_hours * 3600.0) / (_dt_accel * max(steps, 1))
               if accel_real_hours > 0 else 1.0)
    if accel_real_hours > 0:
        print(f"  [accel] S={S_accel:.3e}  run represents {accel_real_hours:.1f}h real "
              f"(dt={_dt_accel:.1e}s, {steps} steps)", flush=True)
        # UNIFY division with the global S: when accel_real_hours is set, division uses the SAME S
        # (div_real_hours defaults to accel_real_hours) so its kinetics are not double-accelerated.
        if division and div_real_hours <= 0.0:
            div_real_hours = accel_real_hours

    # M2 lamellipodium host (rim detection one-shot at build; advances at cadence).
    lam = None
    if lamellipodium:
        lam = LamellipodiumHost(pos0=pos_a, cof=cof_a, n_cells=n_cells, z0=z0, R=R, dt=dt,
                                params=LamelParams(substrate_clutch=lamel_clutch, batch_steps=active_batch,
                                                   S_kinetic=S_accel, all_cell_protrusion=lamel_all_cell),
                                # all-cell geometry lives only in the host path → force it when ON
                                use_gpu_ratchet=(str(device).startswith("cuda") and not lamel_all_cell),
                                device=device)
        print(f"  [lamel] rim cells={lam.n_rim}/{n_cells}  pool={lam.n_pool}  "
              f"p_advance={lam.p_advance:.3e}  z_basal={lam.z_basal*1e6:.3f}um", flush=True)

    # B3 filopodia host (explicit finger protrusions; tips probe + adhere node-FACE to other
    # cells and node-to-plane to the dish). Additive; constructed only when --filopodia.
    filo = None
    if filopodia:
        # Scale the SLOW RATES by S (from the native dataclass defaults, not hard-coded): the tip
        # polymerization velocity (extension advance = v_poly·batch_steps·dt covers S× more real
        # distance per mechanics-step) and the per-batch nucleation probability (capped at 1.0).
        # L_max / k_tip / capture geometry (lengths & FORCES) are LEFT NATIVE — only rates scale.
        _fnat = FilopodiaParams()
        filo = FilopodiaHost(cof=cof_a, n_cells=n_cells, faces=faces_a, fcell=fcell_a,
                             z0=z0, R=R, dt=dt,
                             params=FilopodiaParams(v_poly=_fnat.v_poly * S_accel,
                                                    p_seed=min(1.0, _fnat.p_seed * S_accel),
                                                    batch_steps=active_batch),
                             use_gpu_probe=(gpu_probe or str(device).startswith("cuda")), device=device)   # GPU-only on cuda (radius fix resolved the stall; CPU parity exact)
        print(f"  [filopodia] pool={filo.n_pool}  v_poly={filo.p.v_poly*1e9:.0f}nm/s  "
              f"L_max={filo.p.L_max*1e6:.1f}um  k_tip={filo.p.k_tip:.1e}N/m (node-FACE + node-plane)", flush=True)

    # M3 junction switch host (crowd-pressure cadherin→integrin clutch; latches at cadence).
    # SUPERSEDED by E1 explicit cadherin bonds — disabled whenever cadherin mode is on.
    js = None
    if junction_switch and not cadherin:
        js = JunctionSwitchHost(cof=cof_a, n_cells=n_cells, R=R)
        print(f"  [junction] r_contact={js.r_contact*1e6:.1f}um  P_switch={js.p.P_switch_kPa}kPa  "
              f"cad_weak={js.p.cadherin_weak_factor}  integrin_strong={js.p.integrin_strong_factor}  "
              f"cadence={js.cadence}", flush=True)

    # E1 explicit cadherin trans-dimer bonds (fine-grained cell-cell adhesion; replaces the
    # cohesion tent attraction + M3). In cadherin mode the cohesion/contact run REPULSION-ONLY
    # (adhesion arg → 0) and all adhesion is the catch-bond ensemble; de-cohesion is emergent.
    cad = None
    coh_adh = adh_strength
    if cadherin:
        # ×40 scale bridge (the sanctioned coarse-graining): a node-pair bond is the cadherin
        # BUNDLE over a contact patch, so its rest length / capture radius are the MESH contact
        # scale (r0=c_rep contact-equilibrium, r_bind=c_adh adhesion range) and its effective
        # stiffness is set so the catch→slip transition f0=29.2pN is reached at the capture
        # limit (k_trans = f0/(r_bind−r0)). The Rakshit catch-slip SHAPE + f0 stay molecular —
        # only the geometric scale is bridged. Bonds then form at apposed interface nodes and
        # rupture (slip) under the spreading traction → emergent de-cohesion.
        f0 = 29.2e-12
        r0_meso = c_rep
        rbind_meso = cad_rbind * 1e-6 if cad_rbind > 0.0 else c_adh   # Stage-1 reach [µm]→m (c_adh is in m); fibronectin µm-scale
        k_meso = f0 / max(rbind_meso - r0_meso, 1e-12)
        cad = CadherinBondHost(cof=cof_a, n_cells=n_cells, dt=dt,
                               params=CadherinParams(k_trans=k_meso, r0_trans=r0_meso,
                                                     r_bind=rbind_meso, batch_steps=cad_batch,
                                                     subcycle=cad_subcycle, micro_M=cad_micro_M,
                                                     bundle_n=cad_bundle, f_contract=cad_contract))
        # Accelerate the bond KINETICS by S (NOT the FORCE): the on-rate k_on and the entire
        # force-dependent off-rate k_off(F) lookup table are multiplied by S, so bonds form and
        # rupture S× faster while the trans-dimer FORCE constant k_trans (catch-slip f0=29.2pN,
        # Rakshit shape) is UNTOUCHED. The exact survival prob 1−exp(−k·Δt_batch) stays faithful.
        if S_accel != 1.0:
            cad.p.k_on *= S_accel
            cad._koff *= S_accel
        if cad_subcycle and cad.p.k_on > 0.0:
            # #1 Step 5: at large mechanics dt the KMC must ALSO be CALLED often enough that mechanics
            # sees FRESH bonds — at least once per cadherin timescale τ_cad=1/k_on. Derived cadence
            # (not tuned): batch = clamp(floor(τ_cad/dt), 1, cad_batch). Combined with the in-call
            # sub-cycling (dcm_cadherin_host.update_gpu) this resolves BOTH the intra-batch bond
            # turnover AND the bond↔mechanics coupling at large dt. At base dt this is a no-op
            # (floor(τ_cad/dt) ≫ cad_batch → stays cad_batch).
            tau_cad = 1.0 / float(cad.p.k_on)
            cad.batch_steps = max(1, min(cad.batch_steps, int(tau_cad / max(dt, 1e-30))))
            cad.dt_batch = cad.batch_steps * cad.dt
        # A1: node-FACE coupling. By default the sparse cadherin node-NODE bonds are the sole
        # adhesion (coh_adh=0) — but they bond only the few apposed node-pairs (~6/junction,
        # mesh-density-limited) so cells touch at POINTS and stay round (Ψ≈0.99 "bag of marbles";
        # the pure-aggregation test converged there). --coupling RE-ENABLES the continuous
        # node-vs-FACE bilinear adhesion (contact_grid kernel's `adh` branch = SimuCell3D Model-1
        # traction–separation, mesh-INDEPENDENT) so EVERY apposed node↔face pair adheres → a flat
        # shared interface → cells flatten into a real tissue (regime II). adh_strength is the
        # aggregation-validated node-face cohesion (lit-anchored ~4–5e7; DCM-aggregation-resolved).
        # Cadherin catch-bonds stay ON for the mechanistic, force-dependent de-cohesion (the
        # node-face softening branch adds a geometric separation release on top).
        coh_adh = adh_strength if coupling else 0.0
        if coupling:
            print(f"  [coupling] node-FACE bilinear adhesion ON (coh_adh={coh_adh:.1e}Pa) — "
                  f"continuous interface for regime-II flattening (+ cadherin catch-bond de-cohesion)", flush=True)
        print(f"  [cadherin] k_trans={cad.p.k_trans:.2e}N/m  r0={cad.p.r0_trans*1e6:.2f}um  "
              f"r_bind={cad.p.r_bind*1e6:.2f}um  k_on={cad.p.k_on:.1f}/s  batch={cad.batch_steps}  "
              f"bundle_n={cad_bundle:.0f} (per-bond force {29.2*cad_bundle/1000:.2f}nN; "
              f"catch-slip f0=29.2pN @ molecular load, Rakshit ×40 bridge)", flush=True)

    if bending:
        print(f"  [bending] k_bend={k_bend:.2e}N/m  (biharmonic Δ²r thin-plate; resists curvature "
              f"variation, distinct from surface tension)", flush=True)
    if surface_tension:
        print(f"  [surface-tension] gamma={gamma_surf:.2e}N/m  k_area={k_area:.2e}N/m  "
              f"A0_cell={float(_A0c.mean())*1e12:.1f}um^2 (area-minimising membrane tension"
              f"{' + global area constraint' if k_area>0 else ''})", flush=True)
    if nucleus:
        print(f"  [nucleus] R_nuc={R_nuc*1e6:.2f}um ({R_nuc_factor:.2f}·R)  E_nuc={E_nuc:.0f}Pa  "
              f"ratio_lamin={ratio_lamin}  knee={knee_strain}  k_chrom={k_chrom:.2e}N/m  "
              f"k_lamin={k_lamin:.2e}N/m  (deformable core; bilinear chromatin→lamin)", flush=True)

    # C8 necrosis 3-zone (depth-from-surface): assigns prolif/quiescent/necrotic, softens the
    # necrotic core's turgor + gates division to the proliferating rim. Inert until N is large.
    necro = None
    if necrosis:
        necro = NecrosisHost(cof=cof_a, n_cells=n_cells, npc=npc, R=R)
        print(f"  [necrosis] d_prolif={necro.p.d_prolif_um}um  d_necrotic={necro.p.d_necrotic_um}um "
              f"(O2-depth proxy)  turgor_necrotic={necro.p.turgor_necrotic}  batch={necro.batch_steps}",
              flush=True)

    # C7 cell division (proliferation): rim cells divide into parked pool cells at low cadence.
    div = None
    if division:
        verts0 = icosphere_mesh(R, subdiv)[0]
        p_div_use = div_rate
        # division check CADENCE in STEP units. The default 2000 was tuned for the 8e-5/25000-step
        # regime (~12 checks/run). Under a larger accel_dt the run has FEWER steps for the same real
        # time, so a fixed 2000-step cadence would fire only ~once — under-dividing. Scale the cadence
        # so there are ~12 division rounds REGARDLESS of dt (a TIME-proportional cadence, dt-agnostic),
        # capped to the default and floored so it never gets pathologically small.
        div_batch = int(min(DivisionParams().batch_steps, max(50, steps // 12)))
        if div_real_hours > 0.0:
            # TIME-CONSISTENT division (PI 2026-06-24): the run REPRESENTS `div_real_hours` of real
            # biological time. We cannot integrate the mechanics for real hours (10^7-10^9 steps), so
            # we accelerate by S such that S·dt·steps = div_real_hours (real seconds); the physical
            # per-check division probability is then p_div = (S·dt·div_every)/T_cycle (sim_realtime,
            # no free knob). This SYNCHRONISES division with the (accelerated) mechanical clock —
            # cells divide ~div_real_hours/T_cycle times over the run, matching MCF7 (~24 h cycle).
            from ffn_sim.common.sim_realtime import division_probability
            _dt = accel_dt if accel_dt else dt
            _S = (div_real_hours * 3600.0) / (_dt * max(steps, 1))
            _anchor = division_probability(steps=steps, dt_s=_dt, t_cycle_h=div_t_cycle_h,
                                           div_every=div_batch,
                                           n_rim_est=max(1, n_active // 3), accel=_S)
            p_div_use = float(min(1.0, _anchor.p_div))
            print(f"  [division] TIME-CONSISTENT: run≈{div_real_hours:.1f}h real · T_cycle={div_t_cycle_h:.0f}h "
                  f"· S={_S:.2e} · div_every={div_batch} (~{steps//div_batch} rounds) → p_div={p_div_use:.3e} "
                  f"(vs bare {div_rate})", flush=True)
        div = DivisionHost(verts0=verts0, npc=npc, n_total=n_cells, n_active0=n_active,
                           R=R, z0=z0, V0_full=V0,
                           params=DivisionParams(p_div=p_div_use, batch_steps=div_batch))
        # time-consistent post-mitotic REGROWTH increment (derived, not tuned): a daughter born at
        # V0/2 re-inflates to V0 over ONE cell cycle. steps_per_cycle maps the cell-cycle time onto
        # the (time-accelerated) sim clock — under time-accel the run represents div_real_hours over
        # `steps`, so a cycle of div_t_cycle_h occupies steps·(t_cycle/real_hours) steps; per batch
        # V0_cell climbs by (V0/2)·(batch/steps_per_cycle). No time-accel → regrow over the run.
        if div_real_hours > 0.0:
            steps_per_cycle = max(1.0, steps * (div_t_cycle_h / max(div_real_hours, 1e-9)))
        else:
            steps_per_cycle = float(max(steps, 1))
        div_dV0_batch = V0 * (div.batch_steps / steps_per_cycle)   # daughter grows tiny→V0 (full) over one cycle
        print(f"  [division] n_active={n_active} + parked pool={n_parked} (n_total={n_cells})  "
              f"p_div={p_div_use:.3e}  batch={div.batch_steps}  rim-cell proliferation "
              f"(mitotic-round split, V0/2→V0 over {steps_per_cycle:.0f} steps, "
              f"dV0/batch={div_dV0_batch:.2e})", flush=True)

    # C6 explicit integrin-ECM catch-slip clutch (replaces the wetting proxy). Basal nodes grip
    # the dish; Pereverzev catch-slip governs hold/release → traction-limited mechanistic spread.
    ecm = None
    if ecm_clutch:
        ecm = EcmClutchHost(cof=cof_a, n_cells=n_cells, z0=z0, R=R, dt=dt, c_adh=c_adh,
                            params=EcmClutchParams(bundle_n=ecm_bundle, ligand_density=ecm_ligand))
        print(f"  [ecm-clutch] k_fa={ecm.k_fa:.2e}N/m  engage={ecm.engage_range*1e6:.2f}um  "
              f"k_on={ecm.p.k_on:.1f}/s  batch={ecm.fa_batch_steps}  bundle_n={ecm_bundle:.0f} "
              f"(per-FA force {30.0*ecm_bundle/1000:.2f}nN; Pereverzev F_s=30pN @ per-integrin load) "
              f"— WETTING OFF", flush=True)

    # A1 REMESH (host-side, low cadence): keep every edge in [l_min, 3·l_min] so large
    # spreading never stretches a triangle into a sliver (the contact penalty ∝ A_face fails
    # on slivers → interpenetration). SPLIT activates a dormant pool node, COLLAPSE returns
    # one, SWAP fixes slivers at fixed node count — all manifold/winding/volume preserving
    # (turgor V stays consistent). r0 is reset to the current edge length (relaxed reference)
    # since edge identity is lost across the re-triangulation; the cortex bond just holds the
    # new mesh. N (incl. pool) is FIXED → only the face/edge arrays + cof/pos resync to device.
    remesh_l_min = 0.5 * mean_edge
    remesh_stats = {"calls": 0, "swap": 0, "split": 0, "collapse": 0, "pool_exhausted": 0}

    def do_remesh():
        nonlocal faces_a, fcell_a, cof_a, edges_a, n_faces, n_edges
        nonlocal faces_d, fcell_d, edges_d, r0_d, cent_f32, edge_cell_d, emid_f32
        P = pos_d.numpy().astype(np.float64)
        pos_new, faces_new, cof_new, fc_new, counts = remesh_pass(
            P, faces_a, cof_a, remesh_l_min, face_cell=fcell_a,
            max_ops=max(32, 3 * n_cells), sliver_q=0.2, park=PARK_POS)
        remesh_stats["calls"] += 1
        for k in ("swap", "split", "collapse", "pool_exhausted"):
            remesh_stats[k] += counts[k]
        if counts["swap"] + counts["split"] + counts["collapse"] == 0:
            return                                   # nothing out of band — no resync
        faces_a = faces_new.astype(np.int64)
        cof_a = cof_new.astype(np.int64)
        fcell_a = fc_new.astype(np.int64)
        edges_a = _edges_from_faces(faces_a).astype(np.int64)
        n_faces = faces_a.shape[0]
        n_edges = edges_a.shape[0]
        r0_new = np.linalg.norm(pos_new[edges_a[:, 0]] - pos_new[edges_a[:, 1]], axis=1)
        pos_d.assign(np.ascontiguousarray(pos_new))
        cof_d.assign(cof_a.astype(np.int32))
        faces_d = wp.array(faces_a.astype(np.int32), dtype=wp.int32, device=device)
        fcell_d = wp.array(fcell_a.astype(np.int32), dtype=wp.int32, device=device)
        edges_d = wp.array(edges_a.astype(np.int32), dtype=wp.int32, device=device)
        r0_d = wp.array(r0_new, dtype=wp.float64, device=device)
        cent_f32 = wp.zeros(n_faces, dtype=wp.vec3, device=device)
        edge_cell_d = wp.array(cof_a[edges_a[:, 0]].astype(np.int32), dtype=wp.int32, device=device)
        emid_f32 = wp.zeros(n_edges, dtype=wp.vec3, device=device)
        # BUG FIX (remesh audit 2026-06-23): do_remesh REASSIGNS cof_a to a new array, severing the
        # shared reference the host binders captured at construction. ALL FOUR cof-caching binders
        # must be re-pointed (the division path at the step loop already does all four) — the prior
        # code updated only lam/js, so cad/ecm read a STALE cof and missed dropping bonds whose node
        # was COLLAPSE-parked → a box-scale spurious force on the next binder update. Mirror division.
        if lam is not None:
            lam.cof = cof_a            # activated/collapsed nodes changed the node→cell map
        if js is not None:
            js.cof = cof_a
        if cad is not None:
            cad.cof = cof_a
        if ecm is not None:
            ecm.cof = cof_a

    def do_division_cleave(force_cell=None):
        """C7 division by IN-PLACE MESH CLEAVAGE (SimuCell3D cytokinesis): carve each rim mother's
        shell into mother(−)+daughter(+) along its Hertwig plane. The two daughters' UNION equals the
        mother at the division instant ⇒ the footprint is unchanged and NO neighbour is displaced
        (zero contact spike — the fix for the mitotic-round/parked-insert collision). Resizes faces
        exactly like :func:`do_remesh`; the septum ring/cap nodes come from the dormant pool; the
        daughter takes a free (parked) cell id whose pre-made parked faces are dropped first. Returns
        the number of cells that divided this tick."""
        nonlocal faces_a, fcell_a, cof_a, edges_a, n_faces, n_edges
        nonlocal faces_d, fcell_d, edges_d, r0_d, cent_f32, edge_cell_d, emid_f32
        P = pos_d.numpy().astype(np.float64)
        plan = div.cleave_plan(P, cof_a, force_cell=force_cell,
                               can_divide=(necro.can_divide if necro is not None else None))
        if not plan:
            return 0
        pos_w, faces_w, fcell_w, cof_w = P, faces_a, fcell_a, cof_a
        done = 0
        for (mother, daughter, p0, nrm) in plan:
            keep = fcell_w != daughter               # drop the parked daughter's pre-made faces
            try:
                pos_w, faces_w, fcell_w, cof_w, _info = cleave_cell(
                    pos_w, faces_w[keep], fcell_w[keep], cof_w,
                    cell_id=int(mother), daughter_id=int(daughter), p0=p0, n=nrm,
                    sep=0.6 * mean_edge, park=PARK_POS)
            except (RuntimeError, ValueError):
                continue                             # pool short / degenerate plane → skip this one
            vm = float(V0_cell[mother])              # split the rest-volume setpoint in half; regrow ramp lifts both
            V0_cell[mother] = 0.5 * vm
            V0_cell[daughter] = 0.5 * vm
            done += 1
        if done == 0:
            return 0
        # NB septum-region mesh quality: a plane cut clips body triangles near their vertices → small
        # SLIVER sub-triangles at the cut interface (q~0.04; the body mesh stays pristine q~0.98). They
        # are SMALL-AREA + internal (the cell–cell septum) so the node-FACE contact stays bounded — the
        # N=400 cleave run held pen~3.9 / cfl~0.5 to the end with them present (no blow-up). A single
        # remesh_pass made quality WORSE (the high-valence/cut config defeats the conservative COLLAPSE),
        # so polishing the septum to bulk quality is left to a robust local refiner (follow-up), not a
        # churning per-division pass. The Delaunay cap already removes the worst case (the fan's
        # valence-≈n_ring centre node).
        faces_a = faces_w.astype(np.int64); fcell_a = fcell_w.astype(np.int64); cof_a = cof_w.astype(np.int64)
        edges_a = _edges_from_faces(faces_a).astype(np.int64)
        n_faces = faces_a.shape[0]; n_edges = edges_a.shape[0]
        r0_new = np.linalg.norm(pos_w[edges_a[:, 0]] - pos_w[edges_a[:, 1]], axis=1)
        pos_d.assign(np.ascontiguousarray(pos_w))
        cof_d.assign(cof_a.astype(np.int32))
        faces_d = wp.array(faces_a.astype(np.int32), dtype=wp.int32, device=device)
        fcell_d = wp.array(fcell_a.astype(np.int32), dtype=wp.int32, device=device)
        edges_d = wp.array(edges_a.astype(np.int32), dtype=wp.int32, device=device)
        r0_d = wp.array(r0_new, dtype=wp.float64, device=device)
        cent_f32 = wp.zeros(n_faces, dtype=wp.vec3, device=device)
        edge_cell_d = wp.array(cof_a[edges_a[:, 0]].astype(np.int32), dtype=wp.int32, device=device)
        emid_f32 = wp.zeros(n_edges, dtype=wp.vec3, device=device)
        V0_cell_d.assign(V0_cell)
        if lam is not None: lam.cof = cof_a
        if js is not None: js.cof = cof_a
        if cad is not None: cad.cof = cof_a
        if ecm is not None: ecm.cof = cof_a
        if filo is not None:
            filo.cof = cof_a; filo.faces = faces_a; filo.fcell = fcell_a
        return done

    def stiff_force_into(pos_buf, out_d):
        """I-opt: the STIFF (CFL-setting + structural) force on an arbitrary position buffer, with
        FROZEN neighbour grids (built on x_n by step_once) → the implicit operator's force eval.
        Excludes the soft/explicit spreading drivers (wetting, lamellipodium, cadherin/clutch
        springs) — those are in the RHS b. qpts=node_f32 (frozen) freezes connectivity for a valid
        Hessian; pos_buf supplies the force math."""
        wp.launch(_zero_vec, dim=N, inputs=[out_d], device=device)
        if use_grid:
            if frozen_cache is not None:
                # I-opt #2: replay the xₙ-cached cohesion candidates (no grid query); same force law.
                frozen_cache.cohesion_into(pos_buf, cof_d, r_contact=r_contact, c_adh=c_adh,
                                           rep=rep_strength, omega=coh_adh, A=area_per_node,
                                           force_cap=force_cap, out_d=out_d,
                                           cad_d=(cad_d if js is not None else None))
            elif js is not None:
                wp.launch(cohesion_grid_cad_kernel, dim=N,
                          inputs=[node_grid.id, node_f32, pos_buf, cof_d, cad_d, wp.float32(coh_q),
                                  wp.float64(r_contact), wp.float64(c_adh), wp.float64(rep_strength),
                                  wp.float64(coh_adh), wp.float64(area_per_node), wp.float64(force_cap), out_d], device=device)
            else:
                wp.launch(cohesion_grid_kernel, dim=N,
                          inputs=[node_grid.id, node_f32, pos_buf, cof_d, wp.float32(coh_q),
                                  wp.float64(r_contact), wp.float64(c_adh), wp.float64(rep_strength),
                                  wp.float64(coh_adh), wp.float64(area_per_node), wp.float64(force_cap), out_d], device=device)
        Vc_d.zero_()
        wp.launch(dcm_volume_kernel, dim=n_faces, inputs=[pos_buf, faces_d, fcell_d, Vc_d], device=device)
        if osmotic:
            wp.launch(_dp_from_vol_osm, dim=n_cells, inputs=[Vc_d, V0_cell_d, wp.float64(V0),
                      wp.float64(p.turgor_dP0), wp.float64(k_vol), dP_d], device=device)
        else:
            wp.launch(_dp_from_vol_pc, dim=n_cells, inputs=[Vc_d, V0_cell_d, wp.float64(p.turgor_dP0),
                      wp.float64(k_vol), dP_d], device=device)
        if necro is not None:
            wp.launch(scale_per_cell_kernel, dim=n_cells, inputs=[dP_d, turgor_mult_d], device=device)
        wp.launch(dcm_turgor_force_kernel, dim=n_faces, inputs=[pos_buf, faces_d, fcell_d, dP_d, out_d], device=device)
        # node-FACE contact. With --ipc BOTH the excluded-volume barrier AND the node-FACE adhesion
        # are fed to device_cg as analytic Hessians via hess_apply (cn_k/cn_nrm from
        # nearest_face_ipc_kernel, built once per step) — so the FD operator here EXCLUDES them (the
        # adhesion is a linear spring → PSD-projected analytic stiffness; its exact multi-face FORCE
        # is in the RHS). This removes the per-CG-iter node-FACE query from the matvec (the ~3×
        # coupling cost). Only the non-IPC path keeps the contact in the FD matvec.
        if use_grid and not ipc:
            if frozen_cache is not None:
                # I-opt #2: replay the xₙ-cached contact faces (no grid query); same penalty law.
                frozen_cache.contact_into(pos_buf, cof_d, faces_d, fcell_d, rep=rep_strength,
                                          adh=coh_adh, c_rep=c_rep, c_adh=c_adh, out_d=out_d,
                                          cad_d=(cad_d if js is not None else None))
            elif js is not None:
                wp.launch(contact_grid_cad_kernel, dim=N,
                          inputs=[face_grid.id, node_f32, pos_buf, cof_d, faces_d, fcell_d, cad_d,
                                  wp.float32(con_q), wp.float64(rep_strength), wp.float64(coh_adh),
                                  wp.float64(c_rep), wp.float64(c_adh), out_d], device=device)
            else:
                wp.launch(contact_grid_kernel, dim=N,
                          inputs=[face_grid.id, node_f32, pos_buf, cof_d, faces_d, fcell_d,
                                  wp.float32(con_q), wp.float64(rep_strength), wp.float64(coh_adh),
                                  wp.float64(c_rep), wp.float64(c_adh), out_d], device=device)
        wp.launch(_bond_accumulate, dim=n_edges, inputs=[pos_buf, edges_d, wp.float64(p.k_edge), r0_d, out_d], device=device)
        if use_substrate_well:
            wp.launch(dcm_substrate_well_accum_kernel, dim=N, inputs=[pos_buf, wp.float64(z0),
                      wp.float64(k_well), wp.float64(adh_range), wp.float64(k_floor), out_d], device=device)
        if ubottom:                                  # U-bottom ULA confinement (own-row ADD; sibling of the flat well)
            wp.launch(dcm_ubottom_well_kernel, dim=N, inputs=[pos_buf, wp.float64(ub_cx),
                      wp.float64(ub_cy), wp.float64(ub_cz), wp.float64(ub_rwell),
                      wp.float64(ub_kwall), out_d], device=device)
        if nucleus:
            csum_d.zero_(); ccnt_d.zero_()
            wp.launch(cell_centroid_accum_kernel, dim=N, inputs=[pos_buf, cof_d, csum_d, ccnt_d], device=device)
            wp.launch(nucleus_force_kernel, dim=N, inputs=[pos_buf, cof_d, csum_d, ccnt_d, wp.float64(R_nuc),
                      wp.float64(d_knee_nuc), wp.float64(k_chrom), wp.float64(k_lamin), out_d], device=device)
        if bending:
            nsum_d.zero_(); ncnt_d.zero_()
            wp.launch(edge_neighbor_sum_kernel, dim=n_edges, inputs=[pos_buf, edges_d, nsum_d, ncnt_d], device=device)
            wp.launch(umbrella_kernel, dim=N, inputs=[pos_buf, nsum_d, ncnt_d, lap_d], device=device)
            nsum_d.zero_(); ncnt_d.zero_()
            wp.launch(edge_neighbor_sum_kernel, dim=n_edges, inputs=[lap_d, edges_d, nsum_d, ncnt_d], device=device)
            wp.launch(umbrella_kernel, dim=N, inputs=[lap_d, nsum_d, ncnt_d, bilap_d], device=device)
            wp.launch(bending_apply_kernel, dim=N, inputs=[bilap_d, cof_d, wp.float64(k_bend), out_d], device=device)

    # ---- #1 Step 3: projected-Newton IPC closures (called ONLY under --ipc-newton) --------------
    # The Newton loop re-linearises the CFL-setting STIFF set {turgor, cortex edges, IPC barrier} and
    # runs an energy line-search on Φ=½a|x−xn|²+U_stiff(x)−F_soft(xn)·(x−xn); every OTHER driver
    # (cohesion, node-face adhesion, wetting, well/ubottom, nucleus, bending, ECM, lamellipodium,
    # filopodia, gravity) is LAGGED at xn as F_soft (linear potential) — the plan's soft/explicit set.
    def _nt_rebuild(pos_buf):
        wp.launch(pos_to_f32, dim=N, inputs=[pos_buf, node_f32], device=device)
        wp.launch(face_centroids_f32, dim=n_faces, inputs=[pos_buf, faces_d, cent_f32], device=device)
        face_grid.build(points=cent_f32, radius=grid_q)

    def _nt_turgor_into(pos_buf, out_d):
        Vc_d.zero_()
        wp.launch(dcm_volume_kernel, dim=n_faces, inputs=[pos_buf, faces_d, fcell_d, Vc_d], device=device)
        if osmotic:
            wp.launch(_dp_from_vol_osm, dim=n_cells, inputs=[Vc_d, V0_cell_d, wp.float64(V0),
                      wp.float64(p.turgor_dP0), wp.float64(k_vol), dP_d], device=device)
        else:
            wp.launch(_dp_from_vol_pc, dim=n_cells, inputs=[Vc_d, V0_cell_d, wp.float64(p.turgor_dP0),
                      wp.float64(k_vol), dP_d], device=device)
        wp.launch(dcm_turgor_force_kernel, dim=n_faces, inputs=[pos_buf, faces_d, fcell_d, dP_d, out_d], device=device)

    def nt_stiff_into(pos_buf, out_d):
        """FD-matvec operand for device_cg: turgor + cortex edges ONLY (barrier stiffness is the
        analytic ipc_hess; every other force is lagged) — mirrors the Gate-2 stiff_force_into."""
        wp.launch(_zero_vec, dim=N, inputs=[out_d], device=device)
        _nt_turgor_into(pos_buf, out_d)
        wp.launch(_bond_accumulate, dim=n_edges, inputs=[pos_buf, edges_d, wp.float64(p.k_edge), r0_d, out_d], device=device)

    def nt_force_total_into(pos_buf, out_d):
        """F_total(x) = turgor(x) + edges(x) + IPC barrier(x, adh=0) + lagged F_soft(xn). Rebuilds the
        broad-phase on x and RE-LINEARISES the barrier (nearest_face_ipc_kernel fills ipc_cn_k/nrm)."""
        _nt_rebuild(pos_buf)
        nt_stiff_into(pos_buf, out_d)
        wp.launch(nearest_face_ipc_kernel, dim=N, inputs=[face_grid.id, node_f32, pos_buf, cof_d,
                  faces_d, fcell_d, wp.float32(ipc_repel_q), wp.float64(rep_strength), wp.float64(ipc_dhat),
                  wp.float64(0.0), wp.float64(0.0), out_d, ipc_cn_k, ipc_cn_nrm], device=device)
        wp.launch(_vaxpy, dim=N, inputs=[out_d, wp.float64(1.0), fsoft_d], device=device)

    def nt_energy(pos_buf):
        """Scalar Φ(x) merit for the Armijo line-search; ∇Φ == −(a(x−xn)−F_total) by the Gate-1 match."""
        nt_e_d.zero_(); _nt_rebuild(pos_buf)
        wp.launch(ipc_barrier_energy_kernel, dim=N, inputs=[face_grid.id, node_f32, pos_buf, cof_d,
                  faces_d, fcell_d, wp.float32(ipc_repel_q), wp.float64(rep_strength), wp.float64(ipc_dhat),
                  nt_e_d], device=device)
        wp.launch(edge_energy_kernel, dim=n_edges, inputs=[pos_buf, edges_d, wp.float64(p.k_edge), r0_d, nt_e_d], device=device)
        Vc_d.zero_()
        wp.launch(dcm_volume_kernel, dim=n_faces, inputs=[pos_buf, faces_d, fcell_d, Vc_d], device=device)
        wp.launch(turgor_energy_pc_kernel, dim=n_cells, inputs=[Vc_d, V0_cell_d, dp0eff_d, wp.float64(k_vol), nt_e_d], device=device)
        wp.launch(inertial_energy_kernel, dim=N, inputs=[pos_buf, xn_d, wp.float64(_nt_a[0]), cof_d, nt_e_d], device=device)
        wp.launch(soft_linear_energy_kernel, dim=N, inputs=[pos_buf, xn_d, fsoft_d, cof_d, nt_e_d], device=device)
        wp.synchronize_device(device); return float(nt_e_d.numpy()[0])

    def nt_ccd(pos_buf, dpos_buf):
        return ccd_alpha(face_grid.id, node_f32, pos_buf, dpos_buf, cof_d, faces_d, fcell_d,
                         ipc_repel_q, ipc_t, eta=ipc_eta, device=device)

    _agg = [wp.vec3d(0.0, 0.0, 0.0), 0.0]     # [aggregate centroid, ΔP_agg=2σ/R_agg in Pa], refreshed periodically

    def step_once(s, dt_step, do_spread=True):
        wp.launch(_zero_vec, dim=N, inputs=[force_d], device=device)
        if aggregate_tension and (s % agg_every == 0):
            # Foty-Steinberg aggregate liquid-drop: refresh the global centroid + ΔP=2σ/R_agg. DCM positions
            # are in METERS (R_cell=7.5e-6 m), so _R is in meters and ΔP = 2·σ[N/m]/R[m] = Pa directly (SI),
            # matching the turgor kernel's Pa convention — NO 1e6 factor (that was a 10⁶× unit bug, fixed
            # 2026-07-02 per the adversarial audit; the crush it caused faked the earlier "compaction").
            _ph = pos_d.numpy()
            _c, _R = aggregate_centroid_radius(_ph)
            _agg[0] = wp.vec3d(float(_c[0]), float(_c[1]), float(_c[2]))
            _agg[1] = 2.0 * sigma_agg / _R
        if gravity:                       # D7: constant sedimentation body force (RHS only)
            wp.launch(gravity_body_force_kernel, dim=N,
                      inputs=[cof_d, wp.float64(fz_node), force_d], device=device)
        if use_grid:
            # rebuild both grids every step (build is negligible — the path is query-bound;
            # the persistent-grid lever was an honest negative, so every-step grid is fastest)
            wp.launch(pos_to_f32, dim=N, inputs=[pos_d, node_f32], device=device)
            node_grid.build(points=node_f32, radius=coh_q)
            wp.launch(face_centroids_f32, dim=n_faces, inputs=[pos_d, faces_d, cent_f32], device=device)
            face_grid.build(points=cent_f32, radius=grid_q)
            if js is not None:
                # M3: cell-cell adhesion scaled by sqrt(cad_i·cad_j) (identity until a
                # cell crowd-switches, so settle behaviour is unchanged)
                wp.launch(cohesion_grid_cad_kernel, dim=N,
                          inputs=[node_grid.id, node_f32, pos_d, cof_d, cad_d, wp.float32(coh_q),
                                  wp.float64(r_contact), wp.float64(c_adh), wp.float64(rep_strength),
                                  wp.float64(coh_adh), wp.float64(area_per_node),
                                  wp.float64(force_cap), force_d], device=device)
            elif not conservative_contact:        # conservative tent replaces node-node cohesion
                wp.launch(cohesion_grid_kernel, dim=N,
                          inputs=[node_grid.id, node_f32, pos_d, cof_d, wp.float32(coh_q),
                                  wp.float64(r_contact), wp.float64(c_adh), wp.float64(rep_strength),
                                  wp.float64(coh_adh), wp.float64(area_per_node),
                                  wp.float64(force_cap), force_d], device=device)
        else:
            wp.launch(dcm_cohesion_kernel, dim=N,
                      inputs=[pos_d, cof_d, cad_d, wp.int32(0), wp.int32(N),
                              wp.float64(r_contact), wp.float64(c_adh), wp.float64(rep_strength),
                              wp.float64(coh_adh), wp.float64(area_per_node),
                              wp.float64(force_cap), force_d], device=device)
        Vc_d.zero_()
        wp.launch(dcm_volume_kernel, dim=n_faces, inputs=[pos_d, faces_d, fcell_d, Vc_d], device=device)
        if osmotic:
            wp.launch(_dp_from_vol_osm, dim=n_cells, inputs=[Vc_d, V0_cell_d, wp.float64(V0),
                      wp.float64(p.turgor_dP0), wp.float64(k_vol), dP_d], device=device)
        else:
            wp.launch(_dp_from_vol_pc, dim=n_cells,
                      inputs=[Vc_d, V0_cell_d, wp.float64(p.turgor_dP0),
                              wp.float64(k_vol), dP_d], device=device)
        if necro is not None:           # C8: necrotic core loses turgor regulation (dP *= mult)
            wp.launch(scale_per_cell_kernel, dim=n_cells, inputs=[dP_d, turgor_mult_d], device=device)
        wp.launch(dcm_turgor_force_kernel, dim=n_faces,
                  inputs=[pos_d, faces_d, fcell_d, dP_d, force_d], device=device)
        if use_grid and ipc:
            # M1 IPC node-FACE: log-barrier excluded volume into the RHS (force_d) + per-node barrier
            # normal stiffness (ipc_cn_k/ipc_cn_nrm) for the implicit operator. d̂=c_rep, κ=rep_strength
            # (both derived). Already-penetrating nodes get the linear feasibilization push out.
            wp.launch(nearest_face_ipc_kernel, dim=N,
                      inputs=[face_grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d,
                              wp.float32(ipc_repel_q), wp.float64(rep_strength), wp.float64(ipc_dhat),
                              wp.float64(coh_adh), wp.float64(c_adh),
                              force_d, ipc_cn_k, ipc_cn_nrm], device=device)
            if coh_adh > 0.0:                                   # node-face adhesion stays (multi-face), rep=0
                wp.launch(contact_grid_kernel, dim=N,
                          inputs=[face_grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d,
                                  wp.float32(con_q), wp.float64(0.0), wp.float64(coh_adh),
                                  wp.float64(c_rep), wp.float64(c_adh), force_d], device=device)
        elif use_grid:
            if js is not None:
                wp.launch(contact_grid_cad_kernel, dim=N,
                          inputs=[face_grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d, cad_d,
                                  wp.float32(con_q), wp.float64(rep_strength),
                                  wp.float64(coh_adh), wp.float64(c_rep), wp.float64(c_adh),
                                  force_d], device=device)
            elif conservative_contact:
                # the energy fix: conservative bilinear traction-separation tent (SimuCell3D),
                # derived adh=4·w_cs/c_adh². Provides excluded volume + adhesion + the faceting
                # drive in ONE conservative force (replaces the non-conservative constant-force
                # node-face contact + the constant-area node-node cohesion + the diff-γ patch).
                wp.launch(contact_grid_conservative_kernel, dim=N,
                          inputs=[face_grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d,
                                  wp.float32(con_q), wp.float64(rep_cons), wp.float64(adh_cons),
                                  wp.float64(c_adh), force_d], device=device)
            else:
                wp.launch(contact_grid_kernel, dim=N,
                          inputs=[face_grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d,
                                  wp.float32(con_q), wp.float64(rep_strength),
                                  wp.float64(coh_adh), wp.float64(c_rep), wp.float64(c_adh),
                                  force_d], device=device)
        else:
            wp.launch(node_face_contact_kernel, dim=N,
                      inputs=[pos_d, cof_d, faces_d, fcell_d, cad_d, wp.int32(0), wp.int32(n_faces),
                              wp.float64(rep_strength), wp.float64(coh_adh),
                              wp.float64(c_rep), wp.float64(c_adh), force_d], device=device)
        wp.launch(_bond_accumulate, dim=n_edges,
                  inputs=[pos_d, edges_d, wp.float64(p.k_edge), r0_d, force_d], device=device)
        # E2 nucleus deformable core (always on; resists membrane intruding within R_nuc of the
        # cell centroid). Per-cell centroid via an atomic scatter-reduction, then the bilinear push.
        if nucleus:
            csum_d.zero_(); ccnt_d.zero_()
            wp.launch(cell_centroid_accum_kernel, dim=N,
                      inputs=[pos_d, cof_d, csum_d, ccnt_d], device=device)
            wp.launch(nucleus_force_kernel, dim=N,
                      inputs=[pos_d, cof_d, csum_d, ccnt_d, wp.float64(R_nuc),
                              wp.float64(d_knee_nuc), wp.float64(k_chrom), wp.float64(k_lamin),
                              force_d], device=device)
        # B5 thin-plate bending: biharmonic Δ²r via two umbrella-Laplacian passes
        if bending:
            nsum_d.zero_(); ncnt_d.zero_()
            wp.launch(edge_neighbor_sum_kernel, dim=n_edges, inputs=[pos_d, edges_d, nsum_d, ncnt_d], device=device)
            wp.launch(umbrella_kernel, dim=N, inputs=[pos_d, nsum_d, ncnt_d, lap_d], device=device)
            nsum_d.zero_(); ncnt_d.zero_()
            wp.launch(edge_neighbor_sum_kernel, dim=n_edges, inputs=[lap_d, edges_d, nsum_d, ncnt_d], device=device)
            wp.launch(umbrella_kernel, dim=N, inputs=[lap_d, nsum_d, ncnt_d, bilap_d], device=device)
            wp.launch(bending_apply_kernel, dim=N, inputs=[bilap_d, cof_d, wp.float64(k_bend), force_d], device=device)
        # B4 membrane surface tension (area-minimising) + optional global area constraint
        if surface_tension and polarize:
            # apico-basal DIFFERENTIAL tension (Young-Dupre): basal faces wet (gamma - w*w_cs), apical
            # keep full gamma. The directional-spread mechanism; S = w_cs - 2*gamma > 0 => basal spreads.
            wp.launch(polarized_surface_tension_kernel, dim=n_faces,
                      inputs=[pos_d, faces_d, wp.float64(gamma_surf), wp.float64(w_cs_polarize),
                              wp.float64(z0), wp.float64(adh_range), force_d], device=device)
        elif surface_tension and diff_tension and use_grid:
            # DIFFERENTIAL interfacial tension (foam / DAH): cell-cell CONTACT faces get the reduced
            # tension γ_free − w_adh/2 (w_adh = w_cs lit), free faces keep γ_free. The differential is
            # the spreading drive that flattens cells into a faceted tissue (the missing faceting
            # mechanism; uniform γ never spreads contacts). Reuses the face_grid built this step.
            # contact-face tension: derived foam/DAH (γ_free − w_adh/2), OR — when
            # --contact-tension-frac f≥0 is set — the Maître cortex-dissolution limit γ_contact = f·γ_free
            # (cadherin engagement locally disassembles the actomyosin cortex; f→0 = full dissolution).
            _wadh_half = (1.0 - contact_tension_frac) * gamma_surf if contact_tension_frac >= 0.0 \
                else 0.5 * w_cs_jm2
            wp.launch(differential_surface_tension_kernel, dim=n_faces,
                      inputs=[face_grid.id, cent_f32, pos_d, faces_d, fcell_d, wp.float32(con_q),
                              wp.float64(gamma_surf), wp.float64(_wadh_half), force_d], device=device)
        elif surface_tension:
            wp.launch(surface_tension_kernel, dim=n_faces,
                      inputs=[pos_d, faces_d, wp.float64(gamma_surf), force_d], device=device)
        if aggregate_tension and use_grid:
            # AGGREGATE-level Foty-Steinberg liquid-drop tension: inward Laplace ΔP=2σ/R_agg on the FREE
            # (media-exposed) envelope faces, toward the GLOBAL aggregate centroid → contracts + densifies
            # a loose aggregate (the global densification driver per-cell γ cannot supply). Reuses the
            # face_grid + cent_f32 built this step. Composes with turgor (incompressible cells) + IPC.
            wp.launch(aggregate_laplace_kernel, dim=n_faces,
                      inputs=[face_grid.id, cent_f32, pos_d, faces_d, fcell_d, wp.float32(con_q),
                              _agg[0], wp.float64(_agg[1]), force_d], device=device)
        if surface_tension and not (polarize or diff_tension):
            if k_area > 0.0:
                acell_d.zero_()
                wp.launch(face_area_accum_kernel, dim=n_faces,
                          inputs=[pos_d, faces_d, fcell_d, acell_d], device=device)
                wp.launch(global_area_force_kernel, dim=n_faces,
                          inputs=[pos_d, faces_d, fcell_d, acell_d, a0cell_d, wp.float64(k_area),
                                  force_d], device=device)
        # A2 edge-edge excluded volume (always on, like contact — catches the edge-pokethrough
        # node-face misses). Build the edge-midpoint grid then the seg-seg penalty kernel.
        if edge_edge and use_grid:
            wp.launch(edge_midpoints_f32, dim=n_edges, inputs=[pos_d, edges_d, emid_f32], device=device)
            edge_grid.build(points=emid_f32, radius=ee_q)
            wp.launch(edge_edge_contact_kernel, dim=n_edges,
                      inputs=[edge_grid.id, emid_f32, pos_d, edges_d, edge_cell_d,
                              wp.float32(ee_q), wp.float64(c_rep), wp.float64(rep_strength),
                              wp.float64(area_per_node), force_d], device=device)
        # E1 explicit cadherin trans-dimer adhesion (sole cell-cell attraction in cadherin
        # mode; bonds managed on the host at batch cadence, force applied every step)
        if cad is not None and cad._dev is not None and cad._dev["n"] > 0:
            wp.launch(cadherin_bond_force_kernel, dim=cad._dev["n"],
                      inputs=[cad._dev["bonds"], wp.int32(cad._dev["n"]), pos_d,
                              wp.float64(cad.p.k_trans * cad.p.bundle_n), wp.float64(cad.p.r0_trans),
                              wp.float64(cad.p.f_contract * cad.p.bundle_n),
                              force_d],
                      device=device)
        # substrate z-well (own-row accumulate) — pins basal nodes at z0
        if use_substrate_well:
            wp.launch(dcm_substrate_well_accum_kernel, dim=N,
                      inputs=[pos_d, wp.float64(z0), wp.float64(k_well), wp.float64(adh_range),
                              wp.float64(k_floor), force_d], device=device)
        # U-bottom ULA confinement (own-row ADD; sibling of the flat well — INDEPENDENT, so a ULA
        # run has use_substrate_well=False and ubottom=True: bowl confines, surface is non-adhesive)
        if ubottom:
            wp.launch(dcm_ubottom_well_kernel, dim=N,
                      inputs=[pos_d, wp.float64(ub_cx), wp.float64(ub_cy), wp.float64(ub_cz),
                              wp.float64(ub_rwell), wp.float64(ub_kwall), force_d], device=device)
        # substrate in-plane wetting: scatter into wbuf -> cap -> add (mechanistic spread).
        # In C6 ecm-clutch mode the wetting PROXY is OFF — the explicit integrin clutch provides
        # the (traction-limited) substrate coupling instead.
        if substrate_wetting and do_spread and ecm is None:
            wp.launch(_zero_vec, dim=N, inputs=[wbuf_d], device=device)
            if js is not None:
                # M3: per-cell substrate wetting scaled by the integrin gain (raised on switch)
                wp.launch(dcm_wetting_scatter_integrin_kernel, dim=n_faces,
                          inputs=[pos_d, faces_d, fcell_d, integrin_d, wp.float64(z0),
                                  wp.float64(w_cs_jm2), wp.float64(adh_range), wbuf_d], device=device)
            else:
                wp.launch(dcm_wetting_scatter_kernel, dim=n_faces,
                          inputs=[pos_d, faces_d, wp.float64(z0), wp.float64(w_cs_jm2),
                                  wp.float64(adh_range), wbuf_d], device=device)
            wp.launch(dcm_wetting_cap_add_kernel, dim=N,
                      inputs=[wbuf_d, wp.float64(force_cap), force_d], device=device)
        # C6 integrin-ECM clutch traction (every step; engaged set refreshed at FA cadence)
        if ecm is not None and do_spread and ecm._dev is not None and ecm._dev["n"] > 0:
            wp.launch(ecm_clutch_force_kernel, dim=ecm._dev["n"],
                      inputs=[ecm._dev["node"], ecm._dev["anchor"], wp.int32(ecm._dev["n"]),
                              wp.float64(ecm.k_fa * ecm.p.bundle_n), pos_d, force_d], device=device)
        # M2 lamellipodium traction tether (every step; anchors refreshed at cadence)
        if lam is not None and do_spread and lam._dev is not None and lam._dev["n_lead"] > 0:
            d = lam._dev
            wp.launch(gather_lead_pos, dim=d["n_lead"],
                      inputs=[pos_d, d["lead_idx"], d["lead_rp"]], device=device)
            wp.launch(lamellipodium_tether_multicell, dim=d["n_lead"],
                      inputs=[d["lead_rp"], d["lead_ccx"], d["lead_ccy"], d["lead_ox"],
                              d["lead_oy"], d["lead_proj"], d["lead_idx"], d["lead_cell"],
                              d["actin"], d["actin_cell"], wp.int32(d["n_used"]),
                              wp.float64(lam.p.k_tether), wp.float64(lam.p.tether_cap),
                              wp.float64(lam.p.tether_radius), force_d], device=device)
        # B3 filopodia tip adhesions (every step; tips refreshed at cadence): node-FACE pull
        # toward a neighbour's face + node-to-plane clutch to the dish.
        if filo is not None and do_spread and filo._dev is not None:
            fd = filo._dev
            if fd["n_face"] > 0:
                wp.launch(filopodia_tip_face_force_kernel, dim=fd["n_face"],
                          inputs=[fd["face_base"], fd["face_id"], fd["face_bary"], wp.int32(fd["n_face"]),
                                  faces_d, pos_d, wp.float64(filo.p.k_tip), wp.float64(filo.p.force_cap),
                                  force_d], device=device)
            if fd["n_plane"] > 0:
                wp.launch(filopodia_tip_plane_force_kernel, dim=fd["n_plane"],
                          inputs=[fd["plane_base"], fd["plane_anchor"], wp.int32(fd["n_plane"]),
                                  wp.float64(filo.p.k_tip), wp.float64(filo.p.force_cap), pos_d,
                                  force_d], device=device)
        if implicit:
            # IMEX linearly-implicit Euler: (γ/dt·I + K_stiff)Δx = F_total(xₙ) [=force_d].
            # Grids/node_f32 were just built on xₙ above → frozen-neighbour operator. Soft drivers
            # (wetting/lamellipodium/cadherin/clutch) are already in force_d (explicit RHS).
            a_imp = (1.0 / inv_gamma) / dt_step          # γ_node / dt
            if ipc and ipc_newton:
                # #1 Step 3: PROJECTED-NEWTON IPC step (replaces the single linearised solve + one-shot
                # CCD). Split the RHS: F_soft = F_all(xₙ) − stiff{turgor,edges,barrier}(xₙ), held LAGGED
                # over the step; the Newton loop re-linearises the stiff set + CCD-Armijo energy line-search.
                _nt_rebuild(pos_d)                                        # broad-phase on xₙ
                nt_stiff_into(pos_d, nt_Fb)                               # turgor + edges (zeros nt_Fb)
                wp.launch(nearest_face_ipc_kernel, dim=N, inputs=[face_grid.id, node_f32, pos_d, cof_d,
                          faces_d, fcell_d, wp.float32(ipc_repel_q), wp.float64(rep_strength),
                          wp.float64(ipc_dhat), wp.float64(0.0), wp.float64(0.0), nt_Fb, ipc_cn_k,
                          ipc_cn_nrm], device=device)                     # + barrier(xₙ); cn_k/nrm @ xₙ
                wp.launch(_vcopy, dim=N, inputs=[fsoft_d, force_d], device=device)                    # F_all(xₙ)
                wp.launch(_vaxpy, dim=N, inputs=[fsoft_d, wp.float64(-1.0), nt_Fb], device=device)    # − stiff
                wp.launch(turgor_dp0eff_kernel, dim=n_cells, inputs=[V0_cell_d, wp.float64(V0),
                          wp.float64(p.turgor_dP0), wp.int32(1 if osmotic else 0), dp0eff_d], device=device)
                wp.launch(_vcopy, dim=N, inputs=[xn_d, pos_d], device=device)                         # freeze xₙ
                _nt_a[0] = a_imp
                nt_info = ipc_newton_step(pos_d, xn_d, a_imp, force_total_into=nt_force_total_into,
                          stiff_force_into=nt_stiff_into, energy_fn=nt_energy, ccd_alpha_fn=nt_ccd,
                          cof_d=cof_d, scratch=cg_scratch, device=device, max_newton=ipc_newton_max,
                          newton_tol=ipc_newton_tol, hess_apply=ipc_hess, precond_apply=None,
                          cg_maxiter=cg_maxiter)
                if s == 80 or s % 500 == 0:
                    _g = nt_info["g_norm"] / max(nt_info["g0"], 1e-30)
                    print(f"  [ipc-newton] step {s}: newton={nt_info['newton_iters']} "
                          f"cg={nt_info['cg_iters']} |G|/|G0|={_g:.2e} conv={nt_info['converged']} "
                          f"α={','.join(f'{a:.2f}' for a in nt_info['alphas'][:4])}", flush=True)
            else:
                # I-opt #2: snapshot the cohesion + contact candidate lists ONCE here (the grids are on
                # xₙ) so every CG matvec replays them without a hash-grid query. Built at the SAME radii
                # the live kernels use (coh_q / con_q) → byte-identical candidate set, identical force.
                if frozen_cache is not None:
                    frozen_cache.rebuild(node_grid_id=node_grid.id, face_grid_id=face_grid.id,
                                         node_f32=node_f32, cof_d=cof_d, fcell_d=fcell_d,
                                         coh_q=coh_q, con_q=con_q,
                                         do_cohesion=True, do_contact=(not ipc) or (coh_adh > 0.0))
                # I-opt #1: rebuild the analytic Jacobi diagonal once per step (a + edges + turgor +
                # contact). Contact stiffness reuses the IPC per-node cn_k when available (penalty mode
                # has no cheap per-node normal stiffness → that term is omitted, documented).
                if diagA_d is not None:
                    build_diagA(diagA_d, a=a_imp, cof_d=cof_d, edges_d=edges_d, k_edge=p.k_edge,
                                faces_d=faces_d, fcell_d=fcell_d, pos_d=pos_d, k_vol=k_vol, V0=V0,
                                V0_arr=V0_cell_d, cn_k=(ipc_cn_k if ipc else None), device=device)
                dx_d, _cgi = device_cg(stiff_force_into, pos_d, a_imp, force_d, cg_scratch,
                                    maxiter=cg_maxiter, device=device, hess_apply=ipc_hess,
                                    precond_apply=precond_apply)
                if (s == 80 or s % 500 == 0) and isinstance(_cgi, dict):   # CG-iter telemetry (analytic-matvec verification)
                    print(f"  [cg] step {s}: cg_iters={_cgi.get('cg_iters','?')}", flush=True)
                # x += Δx for LIVE nodes only — dormant pool / parked daughters (cof<0) must stay
                # frozen at PARK_POS, exactly as the explicit _bd_step skips cof<0 (review fix #1).
                if ipc:            # M1 IPC: CCD-filtered step — α∈(0,1] keeps every node penetration-free
                    alpha = ccd_alpha(face_grid.id, node_f32, pos_d, dx_d, cof_d, faces_d, fcell_d,
                                      ipc_repel_q, ipc_t, eta=ipc_eta, device=device)
                    # #1 Step-0 telemetry: realized CCD step fraction. α≪1 = the single-linearized-step
                    # stall (CCD clamps the un-line-searched Newton direction; no re-linearization to recover).
                    if os.environ.get("IPC_TELEM") == "1" and (s == 80 or s % 500 == 0):
                        print(f"  [ipc-telem] step {s}: ccd_alpha={float(alpha):.4f} cg_iters="
                              f"{_cgi.get('cg_iters','?') if isinstance(_cgi, dict) else '?'}", flush=True)
                    wp.launch(_vaxpy_active, dim=N, inputs=[pos_d, wp.float64(alpha), dx_d, cof_d], device=device)
                elif pen_cap:      # D8: clamp each node's implicit step to the contact-shell scale
                    wp.launch(_vaxpy_active_capped, dim=N,
                              inputs=[pos_d, dx_d, cof_d, wp.float64(pen_cap_frac * c_rep)], device=device)
                else:
                    wp.launch(_vaxpy_active, dim=N, inputs=[pos_d, wp.float64(1.0), dx_d, cof_d], device=device)
        else:
            wp.launch(_bd_step, dim=N,
                      inputs=[pos_d, force_d, cof_d, wp.float64(inv_gamma), wp.float64(0.0),
                              wp.float64(dt_step), wp.int32(7), wp.int32(s)], device=device)

        if project and use_grid:
            # M1 hard-constraint: geometric non-penetration projection AFTER the integrator step.
            # Rebuild the face grid on the just-moved positions (the step's grid is on xₙ), then run
            # a few Jacobi sweeps that push any penetrating node out to c_rep clearance — a constraint,
            # not a force-balance, so the strong cohesion bundle cannot tunnel it (pen → ~0 vs ~1.5).
            def _rebuild_proj(pp):
                wp.launch(pos_to_f32, dim=N, inputs=[pp, node_f32], device=device)
                wp.launch(face_centroids_f32, dim=n_faces, inputs=[pp, faces_d, cent_f32], device=device)
                face_grid.build(points=cent_f32, radius=grid_q)
            _rebuild_proj(pos_d)
            project_contacts(face_grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d, grid_q, proj_dpos,
                             proj_gap=proj_gap_factor * c_rep, omega=proj_omega, n_iter=proj_iter,
                             rebuild=_rebuild_proj, device=device)

    def _penetration_frac():
        """DEEP-penetration FRACTION = (#nodes with penetration depth > 0.2·R_cell) / N  (0 = watertight).
        Rebuilds the face grid on the CURRENT positions (measure runs after the integration
        step, so the step's grid is stale) then reduces the diagnostic kernel on the host.

        FIX 2026-07-09 (PI-approved gate-contract, DCM_N2000_WATERTIGHT_2026-07-09): was
        max_pen / INITIAL-build mean_edge — the INSET-shrunk build edge ~30× under-inflated the
        denominator and max-only let a single Voronoi-vertex outlier fail the whole gate, reporting a
        false pen_frac≈75 on an N=2000 that was 99.9% watertight. Now a real distribution metric."""
        if not use_grid:
            return 0.0
        wp.launch(pos_to_f32, dim=N, inputs=[pos_d, node_f32], device=device)
        wp.launch(face_centroids_f32, dim=n_faces, inputs=[pos_d, faces_d, cent_f32], device=device)
        face_grid.build(points=cent_f32, radius=grid_q)
        pen_d.zero_()
        wp.launch(penetration_depth_kernel, dim=N,
                  inputs=[face_grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d,
                          wp.float32(con_q), pen_d], device=device)
        wp.synchronize_device(device)
        return float((pen_d.numpy() > 0.2 * R).mean())   # deep-penetration fraction (depth > 0.2·R_cell)

    def measure():
        P = pos_d.numpy().astype(np.float64)
        active = cof_a >= 0                      # DORMANT pool nodes (cof<0) are parked at
        Pa = P[active]                           # PARK_POS (≫ scale) — exclude from every
        finite = bool(np.isfinite(Pa).all())     # measurement + the saved frame, else they
        if not finite:                           # corrupt A/maxZ/COM and the viz axis limits.
            return None, False
        UM = 1e6
        A = _topdown_area_um2(Pa[:, :2] * UM)
        maxZ = float(Pa[:, 2].max() * UM)
        Vsum = float(_cell_volumes(P, faces_a, fcell_a, n_cells).sum())   # faces use active rows
        com = Pa[:, :2].mean(axis=0)
        Psave = P.copy()
        Psave[~active] = Pa.mean(axis=0)         # park dormant at the live centroid for the frame
        # A3 CFL diagnostic: per-step node displacement (= dt·|F|/γ) as a fraction of the
        # repulsion shell c_rep. >~0.3 ⇒ the fastest node can skip the thin contact shell in
        # one step (tunnel through → deep interpenetration the penalty can't recover).
        F = force_d.numpy().astype(np.float64)[active]
        maxF = float(np.linalg.norm(F, axis=1).max()) if F.size else 0.0
        cfl_frac = (dt * maxF / gamma_node) / c_rep
        return {"A_um2": A, "maxZ_um": maxZ, "Vsum": Vsum, "com": com, "P": Psave,
                "pen_frac": _penetration_frac(), "cfl_frac": cfl_frac, "maxF": maxF}, True

    # A3 adaptive substepping: keep the per-step node displacement below cfl_limit·c_rep so a
    # node can't tunnel through the thin contact shell (the pen-runaway mechanism). n_sub is
    # retuned from the measured CFL each frame; cfl_limit=0 disables (n_sub≡1, fixed dt).
    n_sub = [1]

    def stepped(s, dt_step, do_spread=True):
        ns = n_sub[0]
        if ns <= 1:
            step_once(s, dt_step, do_spread=do_spread)
        else:
            sub = dt_step / float(ns)
            for _ in range(ns):
                step_once(s, sub, do_spread=do_spread)

    def retune_substeps(cfl_frac):
        if cfl_limit <= 0.0:
            return
        if cfl_frac > 0.0:
            # cfl_frac is measured with the FULL dt (per-substep cfl = cfl_frac/n_sub), so the
            # substep count that keeps per-substep displacement < cfl_limit·c_rep is just
            # ceil(cfl_frac/cfl_limit) — NOT ×n_sub, which double-counted and ratcheted to the cap
            # (review fix #7).
            target = int(np.ceil(cfl_frac / cfl_limit))
            n_sub[0] = max(1, min(max_substeps, target))

    m_init, ok = measure()
    if not ok:
        return {"error": "non-finite at init"}
    Vinit = m_init["Vsum"]

    # frame/trajectory recording spanning BOTH phases (so the montage shows the
    # aggregation compacting the cube→ball, then the spread). aa0 is normalised at the
    # end to the RESTED baseline (post-settle area) — settle frames then read >1 falling
    # to 1.0 (= compaction), spread frames climb above 1.0 (= spread).
    frame_list = [] if save_frames else None
    cad_list = [] if save_frames else None
    faces_list = [] if save_frames else None     # per-frame topology (remesh changes it)
    cof_list = [] if save_frames else None
    fmag_list = [] if save_frames else None      # per-node net RESIDUAL force magnitude [N] (∝ velocity, →0 at conv.)
    svm_list = [] if save_frames else None       # per-node virial von-Mises STRESS [Pa] — the REAL stress (audit fix #1)
    spress_list = [] if save_frames else None    # per-node hydrostatic pressure [Pa] (turgor-dominated)
    nbond_list = [] if save_frames else None     # per-node cadherin bond count — active JUNCTION
    bondpair_list = [] if save_frames else None  # per-frame (M,2) bonded node pairs — draw the JUNCTIONS as lines
    recs = []   # {gstep, phase, area, maxZ, Vsum, com}

    def cad_now():
        return js.cad_mult.copy() if js is not None else np.ones(n_cells)

    def record(gstep, phase, m):
        recs.append({"gstep": gstep, "phase": phase, "area": m["A_um2"],
                     "maxZ_um": m["maxZ_um"], "Vsum": m["Vsum"], "com": m["com"],
                     "pen_frac": m.get("pen_frac", 0.0), "cfl_frac": m.get("cfl_frac", 0.0)})
        if frame_list is not None:
            frame_list.append(m["P"].astype(np.float32))
            cad_list.append(cad_now().astype(np.float32))
            faces_list.append(faces_a.astype(np.int32))   # snapshot — remesh mutates faces_a
            cof_list.append(cof_a.astype(np.int32))
            # MECH overlays for the viewer: real per-node stress (|net force|) + active cadherin junctions
            try:
                fmag_list.append(np.linalg.norm(force_d.numpy(), axis=1).astype(np.float32))
            except Exception:
                fmag_list.append(np.zeros(N, np.float32))
            nb = np.zeros(N, np.float32)
            bp = np.zeros((0, 2), np.int32)
            _cb = cad.bonds_now() if cad is not None else None
            if _cb is not None and _cb.shape[0]:
                np.add.at(nb, _cb[:, 0], 1.0)
                np.add.at(nb, _cb[:, 1], 1.0)
                bp = _cb.astype(np.int32).copy()
            nbond_list.append(nb)
            bondpair_list.append(bp)
            # REAL stress (audit fix #1): virial/Cauchy von-Mises + hydrostatic pressure from the
            # actual force laws (cadherin deviatoric + cortex-edge deviatoric + turgor hydrostatic).
            # fmag above is the |net residual force| (∝ velocity), NOT a stress — this is.
            try:
                _Pm = pos_d.numpy().astype(np.float64)
                _Vc = np.abs(_cell_volumes(_Pm, faces_a, fcell_a, n_cells))
                _nvol = np.zeros(N)
                _lv = cof_a >= 0
                _nvol[_lv] = _Vc[cof_a[_lv]] / npc
                _keff = (cad.p.k_trans * cad.p.bundle_n) if cad is not None else 0.0
                _r0c = cad.p.r0_trans if cad is not None else 0.0
                _r0e = r0_d.numpy()
                _edges = edges_a if _r0e.shape[0] == edges_a.shape[0] else None   # remesh may desync r0_d
                _svm, _spr, _ = node_virial_stress(
                    _Pm, cof_a, _nvol,
                    bonds=(bp if bp.shape[0] else None), k_trans_eff=_keff, r0_cad=_r0c,
                    edges=_edges, k_edge=p.k_edge, r0_edge=(_r0e if _edges is not None else None),
                    dP_cell=dP_d.numpy())
                svm_list.append(_svm.astype(np.float32))
                spress_list.append(_spr.astype(np.float32))
            except Exception:
                svm_list.append(np.zeros(N, np.float32))
                spress_list.append(np.zeros(N, np.float32))

    record(0, 0, m_init)   # as-built ball

    # gentle soft-start: settle the initial pack at 0.1× dt (the stiff turgor/contact need it)
    for s in range(warmup):
        step_once(s, dt * 0.1, do_spread=False)
    wp.synchronize_device(device)

    # AGGREGATION / SETTLE phase (do_spread=False): the SPHERICAL cluster compacts into a
    # cohesive rounded spheroid RESTING at z0 — substrate well holds basal nodes at z=0,
    # cohesion pulls the surface cells inward (curvature ⇒ net inward pull) + turgor/contact
    # round it. NO wetting / NO lamellipodium / NO junction-switch yet (the aggregate forms
    # cohesively first). The spread is then measured FROM this rested baseline.
    every_s = max(1, settle_steps // settle_frames) if (settle_frames and settle_steps) else settle_steps + 1
    for s in range(1, settle_steps + 1):
        if remesh_period and s % remesh_period == 0:
            wp.synchronize_device(device)
            do_remesh()
        if cad is not None and s % cad.batch_steps == 0:   # E1 bond break/form (aggregate too)
            wp.synchronize_device(device)
            if str(device) != "cpu" and os.environ.get("CAD_GPU", "1") != "0":  # GPU-native break/form (no host round-trip)
                wp.launch(pos_to_f32, dim=N, inputs=[pos_d, node_f32], device=device)
                cad.update_gpu(pos_d, cof_d, node_f32, N, s // cad.batch_steps, device)
            else:
                cad.update(pos_d.numpy().astype(np.float64)); cad.upload(device)
        stepped(s, dt, do_spread=False)
        if s % every_s == 0:
            wp.synchronize_device(device)
            m, ok = measure()
            if ok:
                record(s, 0, m)
                retune_substeps(m["cfl_frac"])
    wp.synchronize_device(device)

    # baseline = the RESTED spheroid (not the as-built ball)
    m0, ok = measure()
    if not ok:
        return {"error": "non-finite after settle"}
    A0 = m0["A_um2"]; V0sum = m0["Vsum"]; com0 = m0["com"]
    record(settle_steps, 1, m0)   # transition frame (start of spread, A/A0 ≡ 1)
    print(f"  [aggregate] done ({warmup} warmup + {settle_steps} settle): A0_rested={A0:.1f}um^2 "
          f"(as-built {m_init['A_um2']:.1f}um^2 → compaction {m_init['A_um2']/A0 if A0 else 0:.2f}×)  "
          f"maxZ={m0['maxZ_um']:.1f}um  Vsum/Vinit={m0['Vsum']/Vinit:.3f}", flush=True)

    every = max(1, steps // max(1, frames))
    t0 = time.perf_counter()
    truncated_at = None
    for s in range(1, steps + 1):
        # M2: ratchet the lamellipodial front on the host at low cadence, then refresh
        # the device anchor + leading-node geometry (the tether reads them every step)
        if lam is not None and (s == 1 or s % lam.batch_steps == 0):
            if getattr(lam, "use_gpu_ratchet", False):
                lam.update(pos_d=pos_d)               # GPU-resident geometry; no full position pull
            else:
                wp.synchronize_device(device)
                lam.update(pos_d.numpy().astype(np.float64))
            lam.upload(device)
        # B3: extend/probe filopodia + refresh tip-adhesion device arrays at low cadence
        if filo is not None and (s == 1 or s % filo.batch_steps == 0):
            wp.synchronize_device(device)
            filo.update(pos_d.numpy().astype(np.float64))
            filo.upload(device)
        # M3: latch the crowd-pressure junction switch at low cadence; re-upload the
        # mutable per-cell cad_mult / integrin_gain only when a new cell switches
        if js is not None and s % js.cadence == 0:
            wp.synchronize_device(device)
            if js.update(pos_d.numpy().astype(np.float64)):
                cad_d.assign(js.cad_mult)
                integrin_d.assign(js.integrin_gain)
                print(f"  [junction] step {s}: {js.n_switched}/{n_cells} cells switched "
                      f"(cad↓{js.p.cadherin_weak_factor}, integrin↑{js.p.integrin_strong_factor})",
                      flush=True)
        # A1: low-cadence host remesh (keeps edges in band → no slivers → contact holds)
        if remesh_period and s % remesh_period == 0:
            wp.synchronize_device(device)
            do_remesh()
        # E1: cadherin bond break/form at the binder cadence (de-cohesion emerges here)
        if cad is not None and s % cad.batch_steps == 0:
            wp.synchronize_device(device)
            if str(device) != "cpu" and os.environ.get("CAD_GPU", "1") != "0":  # GPU-native break/form (no host round-trip)
                wp.launch(pos_to_f32, dim=N, inputs=[pos_d, node_f32], device=device)
                cad.update_gpu(pos_d, cof_d, node_f32, N, s // cad.batch_steps, device)
            else:
                cad.update(pos_d.numpy().astype(np.float64)); cad.upload(device)
            if s % 500 == 0:   # T1 substrate diagnostic: are cadherin bonds actually rupturing+reforming?
                print(f"  [cad-churn] step {s}: n_bonds={cad.n_bonds} cum_broken={cad.n_broken}", flush=True)
        # C6: integrin-ECM clutch engage/break (catch-slip) at the FA cadence
        if ecm is not None and s % ecm.fa_batch_steps == 0:
            wp.synchronize_device(device)
            ecm.update(pos_d.numpy().astype(np.float64))
            ecm.upload(device)
        # C8: necrosis 3-zone — reassign zones, push the per-cell turgor multiplier to device
        if necro is not None and s % necro.batch_steps == 0:
            wp.synchronize_device(device)
            necro.update(pos_d.numpy().astype(np.float64), cof_a)
            turgor_mult_d.assign(necro.turgor_mult)
            cnt = necro.counts()
            if cnt["necrotic"] > 0:
                print(f"  [necrosis] step {s}: prolif={cnt['prolif']} quiescent={cnt['quiescent']} "
                      f"necrotic={cnt['necrotic']}", flush=True)
        # C7: cell division (rim-cell proliferation) — mitotic-round split + regrowth, resync cof+V0_cell
        _div_now = div is not None and (s % div.batch_steps == 0
                                        or (force_divide_step >= 0 and s == force_divide_step))
        _divided = False
        if _div_now:
            wp.synchronize_device(device)
            _force = force_divide_cell if (force_divide_step >= 0 and s == force_divide_step) else None
            if cleave:
                # SimuCell3D in-place mesh cleavage: carve mother→mother+daughter, union=mother → no
                # neighbour displacement (zero contact spike). Resizes faces + resyncs device itself.
                _n_cleaved = do_division_cleave(force_cell=_force)
                _divided = _n_cleaved > 0
                if _divided:
                    n_active_now = int(np.unique(cof_a[cof_a >= 0]).size)
                    print(f"  [division/cleave] step {s}: +{_n_cleaved} cleaved "
                          f"({div.n_divisions} total, {n_active_now} active cells, n_faces={n_faces})",
                          flush=True)
            else:
                P = pos_d.numpy().astype(np.float64)
                _divided = div.update(P, cof_a, V0_cell, force_cell=_force,
                                      can_divide=(necro.can_divide if necro is not None else None))
            if (not cleave) and _divided:
                pos_d.assign(np.ascontiguousarray(P))
                cof_d.assign(cof_a.astype(np.int32))
                V0_cell_d.assign(V0_cell)                 # daughters' halved rest volume → device
                if lam is not None: lam.cof = cof_a
                if js is not None: js.cof = cof_a
                if cad is not None: cad.cof = cof_a
                if ecm is not None: ecm.cof = cof_a
                if filo is not None:
                    # B3 de-conflict (2026-06-24): FilopodiaHost caches self.cof/faces/fcell (its
                    # PROBE + the GPU probe read these). A division flips a parked cell's cof and
                    # re-points cof_a, leaving filo's cached cof STALE → a tip could bond to a now-
                    # active node with the wrong node→cell map (same class as the cad/ecm bug). The
                    # host has no remesh-rebuild method (remesh is disabled when filopodia is on,
                    # so faces/fcell never change under division — but re-point all three to match
                    # the other hosts exactly + stay correct if that guard is ever lifted). n_cells,
                    # n_per_cell (sized to the FULL active+parked pool) and the device adhesion
                    # bundles (rebuilt each filo.upload from current state) need no re-point.
                    filo.cof = cof_a
                    filo.faces = faces_a
                    filo.fcell = fcell_a
                print(f"  [division] step {s}: {div.n_divisions} total divisions "
                      f"({int((cof_a[np.arange(n_cells)*npc]>=0).sum())} active cells)", flush=True)
        # post-mitotic REGROWTH (every batch, division or not): freshly-born half-cells ramp their
        # rest volume V0_cell back toward V0 over a cell cycle; turgor chases it → gradual inflation.
        if div is not None and s % div.batch_steps == 0 and div_dV0_batch > 0.0:
            # cof=None under --cleave: _cell_active's npc-block scan is invalid once cleavage has
            # relabelled the mesh; growing all cells with V0_cell<V0_full is equivalent (parked cells
            # sit at V0_full so are untouched; only the half-volume daughters climb).
            div.grow(V0_cell, div_dV0_batch, None if cleave else cof_a)
            V0_cell_d.assign(V0_cell)
        # DIVISION-PERTURBATION ABSORPTION (PI 2026-06-25): a mitotic split injects a sudden geometric
        # + rest-volume change (two overlapping half-cells) in ONE step → a cfl spike. On a stiff
        # (e.g. faceting-over-compressed) state that spike can tip the implicit CG into non-convergence
        # (a stall). Relax the new cells QUASI-STATICALLY at a small dt (mechanics only, no spread/biology
        # → no biological-clock advance, so S/real-time sync is untouched) BEFORE the main step resumes,
        # so the transient dissipates gradually instead of in one large step.
        if _divided and div_relax_steps > 0:
            for _ in range(div_relax_steps):
                step_once(s, dt * div_relax_factor, do_spread=False)
        # OSMOTIC VOLUME REGULATION (PI 2026-06-25; KB-3.9 water flux): each cell's osmotic setpoint
        # V0_cell relaxes toward its CURRENT measured volume Vc at a rate ∝ its MEDIA-EXPOSED face
        # fraction f_media = 1 − (faces apposed to other cells)/faces_per_cell. So a rim cell (large
        # free surface) sheds water + adapts its rest volume quickly; an interior cell (mostly cell–cell
        # junctions) is buffered (f_media≈0 → barely relaxes). The concentration feedback in
        # _dp_from_vol_osm then sets a STABLE equilibrium — replacing the elastic V0−Vc strain that
        # accumulated stiffness (the faceting-over-compression stall). Mechanics only → S/real-time
        # sync untouched.
        if osmotic:
            # GPU-ONLY (PI 2026-06-26): refresh the per-cell media-exposed contact count on-device every
            # osm_batch (the f_media changes slowly — face grid already built by step_once), then relax
            # V0_cell_d toward Vc_d by a TINY per-step amount (osm_relax/osm_batch) — smooth, no host sync,
            # no periodic bulk jump (the jump was injecting cfl spikes). Vc_d is current from stepped().
            if s % osm_batch == 0:
                contact_cnt_d.zero_()
                wp.launch(face_contact_count_kernel, dim=n_faces,
                          inputs=[face_grid.id, cent_f32, fcell_d, wp.float32(con_q), contact_cnt_d],
                          device=device)
            wp.launch(osmotic_relax_kernel, dim=n_cells,
                      inputs=[V0_cell_d, Vc_d, contact_cnt_d, cof_d, wp.int32(npc),
                              wp.float64(float(faces_per_cell)), wp.float64(osm_relax / max(osm_batch, 1))],
                      device=device)
        stepped(s, dt)
        if s % every == 0 or s == steps:
            wp.synchronize_device(device)
            m, ok = measure()
            if not ok:
                truncated_at = s
                print(f"  [decoh] NON-FINITE at step {s} — truncating", flush=True)
                break
            record(settle_steps + s, 1, m)
            retune_substeps(m["cfl_frac"])
            print(f"  step {s:>7}  A/A0={m['A_um2']/A0 if A0 else 0:.3f}  maxZ={m['maxZ_um']:.1f}um  "
                  f"V/V0={m['Vsum']/V0sum if V0sum else 0:.3f}  "
                  f"drift={np.linalg.norm(m['com']-com0)*1e6:.2f}um  "
                  f"pen={m['pen_frac']:.3f}  cfl={m['cfl_frac']:.2f}  nsub={n_sub[0]}"
                  + (f"  switched={js.n_switched}" if js is not None else "")
                  + (f"  bonds={cad.n_bonds}" if cad is not None else "")
                  + (f"  clutch={ecm.n_clutches}" if ecm is not None else ""), flush=True)
    wp.synchronize_device(device)
    elapsed = time.perf_counter() - t0

    # normalise all recorded frames to the rested baseline + assemble the trajectory
    traj = [{"step": r["gstep"], "phase": r["phase"], "aa0": r["area"] / A0 if A0 > 0 else 0.0,
             "maxZ_um": r["maxZ_um"], "vv0": r["Vsum"] / V0sum if V0sum else 0.0,
             "drift_um": float(np.linalg.norm(r["com"] - com0) * 1e6),
             "pen_frac": r["pen_frac"], "cfl_frac": r["cfl_frac"]} for r in recs]

    if save_frames and frame_list is not None:
        extra = {}
        if remesh_period:
            # topology varies across frames → save per-frame faces/cof as object sequences
            extra["faces_seq"] = np.array(faces_list, dtype=object)
            extra["cof_seq"] = np.array(cof_list, dtype=object)
        if fmag_list is not None:
            extra["fmag"] = np.array(fmag_list, dtype=np.float32)          # (F,N) per-node RESIDUAL force [N]
            extra["svm"] = np.array(svm_list, dtype=np.float32)            # (F,N) per-node virial von-Mises STRESS [Pa]
            extra["spress"] = np.array(spress_list, dtype=np.float32)      # (F,N) per-node hydrostatic pressure [Pa]
            extra["nbond"] = np.array(nbond_list, dtype=np.float32)        # (F,N) per-node junction count
            extra["bondpairs"] = np.array(bondpair_list, dtype=object)     # per-frame (M,2) bonded node pairs
        np.savez_compressed(
            save_frames, frames=np.array(frame_list, dtype=np.float32),
            faces=faces_a.astype(np.int32), cof=cof_a.astype(np.int32),
            cad=np.array(cad_list, dtype=np.float32),
            step=np.array([r["step"] for r in traj]),
            phase=np.array([r["phase"] for r in traj]),
            aa0=np.array([r["aa0"] for r in traj]),
            maxZ=np.array([r["maxZ_um"] for r in traj]),
            vv0=np.array([r["vv0"] for r in traj]),
            pen_frac=np.array([r["pen_frac"] for r in traj]), **extra)
        print(f"  saved {len(frame_list)} frames -> {save_frames}", flush=True)

    spread_aa = [r["aa0"] for r in traj if r["phase"] == 1] or [1.0]
    out = {
        "device": device, "n_cells": n_cells, "N": N, "subdiv": subdiv, "dt": dt,
        "steps": steps, "warmup": warmup, "settle_steps": settle_steps,
        "truncated_at": truncated_at, "steps_per_s": (truncated_at or steps) / elapsed,
        "substrate_wetting": substrate_wetting, "use_substrate_well": use_substrate_well,
        "ubottom": ubottom, "ubottom_r_well_um": (ub_rwell * 1e6 if ubottom else None),
        "lamellipodium": lamellipodium, "junction_switch": junction_switch,
        "A0_rested_um2": A0, "compaction_x": (m_init["A_um2"] / A0) if A0 else 0.0,
        "aa0_peak": max(spread_aa), "aa0_final": spread_aa[-1],
        "maxZ_final_um": traj[-1]["maxZ_um"],
        "vv0_final": traj[-1]["vv0"], "drift_final_um": traj[-1]["drift_um"],
        "pen_frac_peak": max(r["pen_frac"] for r in traj),
        "pen_frac_final": traj[-1]["pen_frac"],
        "cfl_frac_peak": max(r["cfl_frac"] for r in traj),
        "cfl_limit": cfl_limit, "max_substeps": max_substeps, "n_sub_final": n_sub[0],
        "remesh_period": remesh_period, "remesh": remesh_stats if remesh_period else None,
        "edge_edge": edge_edge, "cadherin": cadherin, "ecm_clutch": ecm_clutch,
        "integrator": integrator, "accel_dt": (dt if implicit else None), "cg_maxiter": cg_maxiter,
        "nucleus": nucleus, "surface_tension": surface_tension, "k_vol": k_vol,
        "gamma_surf": gamma_surf if surface_tension else None, "k_area": k_area,
        "bending": bending, "k_bend": k_bend if bending else None, "necrosis": necrosis,
        "builder": builder,
        "nucleus_params": ({"R_nuc": R_nuc, "E_nuc": E_nuc, "k_chrom": k_chrom,
                            "k_lamin": k_lamin, "ratio_lamin": ratio_lamin,
                            "knee_strain": knee_strain} if nucleus else None),
        "W_cs_well_J": W_cs_well, "gamma_node": gamma_node, "trajectory": traj,
    }
    if lam is not None:
        out["lamel"] = {"n_rim": lam.n_rim, "n_pool": lam.n_pool, "n_used": lam.n_used,
                        "n_seeded": lam.n_seeded, "n_advanced": lam.n_advanced,
                        "p_advance": lam.p_advance, "n_lead_final": lam.n_lead}
    if js is not None:
        out["junction"] = {"n_switched": js.n_switched, "n_cells": n_cells,
                           "P_switch_kPa": js.p.P_switch_kPa,
                           "cad_weak": js.p.cadherin_weak_factor,
                           "integrin_strong": js.p.integrin_strong_factor,
                           "max_pressure_kPa": float(js.pressure_kPa.max())}
    if cad is not None:
        out["cadherin_bonds"] = {"n_bonds_final": cad.n_bonds, "n_formed": cad.n_formed,
                                 "n_broken": cad.n_broken, "k_trans": cad.p.k_trans,
                                 "r0_trans": cad.p.r0_trans, "r_bind": cad.p.r_bind,
                                 "bundle_n": cad.p.bundle_n, "k_on": cad.p.k_on,
                                 "batch_steps": cad.batch_steps}
    if necro is not None:
        out["necrosis_stats"] = {**necro.counts(), "d_prolif_um": necro.p.d_prolif_um,
                           "d_necrotic_um": necro.p.d_necrotic_um,
                           "max_depth_um": float(necro.depth_um.max())}
    if div is not None:
        out["division_stats"] = {"n_divisions": div.n_divisions, "n_active0": n_active,
                           "n_active_final": int((cof_a[np.arange(n_cells) * npc] >= 0).sum()),
                           "n_parked": n_parked, "p_div": div_rate}
    if ecm is not None:
        out["ecm_clutch_stats"] = {"n_clutches_final": ecm.n_clutches, "n_engaged": ecm.n_engaged,
                             "n_slipped": ecm.n_slipped, "k_fa": ecm.k_fa,
                             "engage_range": ecm.engage_range, "F_s": ecm.p.F_s,
                             "F_c": ecm.p.F_c, "fa_batch_steps": ecm.fa_batch_steps}

    # D11 VALIDATION GATES — explicit PASS/FAIL contracts evaluated from this run's diagnostics
    # (the gate is a contract written before the run; surface to PI, never loosen inline). The
    # single-run-checkable gates are wired here; G1 (Young-Dupre contact angle) and G3 (the
    # A/A0=a+b/R+c/R² size-law r²≥0.95 fit) are SWEEP-level — they need a multi-N batch, so they
    # are evaluated by the sweep harness, not a single run (flagged below, not silently passed).
    PEN_GATE = 0.005         # G2: deep-penetration(>0.2R) FRACTION (cells must not overlap; 2026-07-09 PI gate-contract, DCM_N2000_WATERTIGHT)
    VV0_TOL = 0.10            # volume conservation sanity
    gates = {}
    gates["finite"] = {"pass": truncated_at is None, "truncated_at": truncated_at}
    gates["G2_interpenetration"] = {"pass": out["pen_frac_peak"] <= PEN_GATE,
                                    "pen_frac_peak": out["pen_frac_peak"], "threshold": PEN_GATE}
    gates["volume_conservation"] = {"pass": abs(out["vv0_final"] - 1.0) <= VV0_TOL,
                                    "vv0_final": out["vv0_final"], "tol": VV0_TOL}
    if division:
        nf = out["division_stats"]["n_active_final"]
        gates["G4_division"] = {"pass": nf > n_active, "n_active0": n_active, "n_active_final": nf,
                                "note": "G4 A/A0∈[2,4] band is checked by the multi-N sweep"}
    gates["G1_contact_angle"] = {"pass": None, "note": "sweep-level (Young-Dupre); not a single run"}
    gates["G3_size_law_fit"] = {"pass": None, "note": "sweep-level (A/A0=a+b/R+c/R² r²≥0.95)"}
    out["gates"] = gates
    verdict = [(k, v["pass"]) for k, v in gates.items() if v["pass"] is not None]
    line = "  ".join(f"{k}={'PASS' if p else 'FAIL'}" for k, p in verdict)
    print(f"  [GATES] {line}", flush=True)
    out["gates_all_pass"] = all(p for _, p in verdict)
    return out


def main():
    ap = argparse.ArgumentParser(description="De-cohesion cleanball spread on the Warp DCM engine (M1: substrate).")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--n-cells", type=int, default=12)
    ap.add_argument("--subdiv", type=int, default=2)
    ap.add_argument("--steps", type=int, default=40000)
    ap.add_argument("--frames", type=int, default=20)
    ap.add_argument("--dt", type=float, default=8.0e-6)
    ap.add_argument("--ipc", action="store_true",
                    help="M1: IPC node-face contact (log-barrier force + analytic Hessian in the implicit "
                         "operator + CCD-filtered step) instead of the capped penalty — guarantees "
                         "non-penetration under the strong cadherin bundle (needs --integrator implicit)")
    ap.add_argument("--ipc-eta", type=float, default=0.9, help="CCD safety fraction (gap stays >= (1-eta)*d)")
    ap.add_argument("--ipc-newton", action="store_true", dest="ipc_newton",
                    help="#1: finish IPC as a PROJECTED-NEWTON incremental-potential step (re-linearise "
                         "stiff{turgor,edges,barrier} + CCD-Armijo energy line-search each step) instead of "
                         "the single linearised solve + one-shot CCD — unlocks large dt (needs --ipc). Every "
                         "other driver (cohesion/cadherin/wetting/well/nucleus/bending/...) is lagged-soft at xₙ.")
    ap.add_argument("--ipc-newton-max", type=int, default=8, dest="ipc_newton_max",
                    help="max Newton iterations per step (projected-Newton IPC)")
    ap.add_argument("--ipc-newton-tol", type=float, default=1e-4, dest="ipc_newton_tol",
                    help="Newton residual tolerance |G|/|G0| (projected-Newton IPC)")
    ap.add_argument("--warmup", type=int, default=1000, help="soft-start steps at 0.1x dt")
    ap.add_argument("--settle-steps", type=int, default=0,
                    help="aggregation/settle steps at full dt with NO spread drivers (rest the spheroid at z0 before the measured spread; baseline A0 is taken AFTER this)")
    ap.add_argument("--settle-frames", type=int, default=0,
                    help="number of frames to capture DURING the aggregate/settle phase (shows the cube→compact-ball compaction in the montage)")
    ap.add_argument("--lamellipodium", action="store_true", help="enable the M2 per-cell lamellipodium crawl")
    ap.add_argument("--lamel-clutch", action="store_true", help="B2: lamellipodium substrate adhesion as a node-to-plane clutch (anchor ON the dish z0, like the ECM clutch) instead of a floating actin bead")
    ap.add_argument("--filopodia", action="store_true", help="B3: explicit filopodia (tips probe + adhere node-FACE to other cells, node-to-plane to the dish)")
    ap.add_argument("--junction-switch", action="store_true", help="enable the M3 crowd-pressure cadherin→integrin junction switch")
    ap.add_argument("--gap", type=float, default=2.05, help="cell centre spacing in R for the spherical aggregate (2.05 = touching/compact)")
    ap.add_argument("--remesh-period", type=int, default=0, help="A1: host SWAP/SPLIT/COLLAPSE remesh every N steps (0=off); keeps edges in band → no slivers")
    ap.add_argument("--pool-factor", type=float, default=0.5, help="dormant node pool size as a fraction of active nodes (for remesh SPLIT)")
    ap.add_argument("--edge-edge", action="store_true", help="A2: edge-edge contact (catches the edge-pokethrough node-face misses)")
    ap.add_argument("--cfl-limit", type=float, default=0.0, help="A3: max per-step node displacement / c_rep (e.g. 0.3); adaptive substeps keep below it (0=off)")
    ap.add_argument("--max-substeps", type=int, default=16, help="A3: cap on adaptive substeps per step")
    ap.add_argument("--cadherin", action="store_true", help="E1: explicit cadherin catch-bonds (fine-grained adhesion; replaces cohesion tent + M3 switch; de-cohesion emergent)")
    ap.add_argument("--ecm-clutch", action="store_true", help="C6: explicit Pereverzev catch-slip integrin-ECM clutch (replaces the wetting proxy; traction-limited spread)")
    ap.add_argument("--cad-batch", type=int, default=50, help="E1 cadherin bond-management cadence (host-hybrid; 50 keeps the GPU↔CPU sync amortised)")
    ap.add_argument("--no-cad-subcycle", action="store_false", dest="cad_subcycle",
                    help="#1 Step 5: DISABLE cadherin KMC sub-cycling (default ON). With it on, the "
                         "bond KMC is resolved at δt_cad=1/(micro_M·k_on) so large mechanics dt does not "
                         "overshoot the ~36ms turnover; off = single saturated pass per batch (the lag).")
    ap.add_argument("--cad-micro-M", type=int, default=8, dest="cad_micro_M",
                    help="#1 Step 5: cadherin KMC micro-steps per timescale 1/k_on (δt_cad=1/(M·k_on))")
    ap.add_argument("--cad-bundle", type=float, default=1.0, help="E1 cadherin ×N mesoscale FORCE bundle (node-bond = N cadherins; force ×N, koff at molecular F/N). 40 → ~7nN/junction ∈ KB-4.11[1-10nN]. 1=legacy")
    ap.add_argument("--cad-contract", type=float, default=0.0, help="Stage-2 compaction motor: active actomyosin junctional CONTRACTION [N per single trans-dimer, bundle-scaled]. Always pulls bonded cells together (RhoA/ROCK-gated NMII the passive catch-bond lacks). SWEEP as a controlled variable (per-motor ~5-15pN × engaged); 0=off. NEVER tune to a compaction target.")
    ap.add_argument("--cad-rbind", type=float, default=0.0, help="Stage-1 long-range reach: ECM-tether (fibronectin, µm-scale) cadherin bond capture radius [µm] override; bridges a LOOSE aggregate (bonds=0 at cadherin's ~1.8µm range otherwise). 0=use c_adh.")
    ap.add_argument("--ecm-bundle", type=float, default=1.0, help="C6 ecm-clutch ×N FA-patch FORCE bundle (node-clutch = N integrins; force ×N, koff at per-integrin F/N). 167 → ~5nN/FA ∈ KB-2.12. 1=legacy")
    ap.add_argument("--ligand-density", type=float, default=1.0, help="C4: substrate ECM ligand-coating density (Bare/Pre/Lam4) — scales the clutch engagement on-rate (more ligand → more engaged FAs → more traction). 1=baseline(Bare); set per-condition to the Lam4>Pre>Bare experimental ordering (NOT tuned)")
    ap.add_argument("--no-pen-cap", dest="pen_cap", action="store_false", help="D8: disable the implicit per-node displacement cap (= the contact-shell clamp that stops frozen-grid tunneling/interpenetration). On by default")
    ap.add_argument("--pen-cap-frac", type=float, default=1.0, help="D8: implicit step cap = frac·c_rep (contact-shell scale; geometric)")
    ap.add_argument("--gravity", action="store_true", help="D7: net gravity−buoyancy sedimentation body force (Δρ·g·v_node, derived; rests the spheroid on the dish at its physiological ~1pN/cell weight)")
    ap.add_argument("--delta-rho", type=float, default=55.0, help="D7: ρ_cell−ρ_medium [kg/m³] (MCF7 ~1060 − medium ~1005; SimuCell3D Table 2)")
    ap.add_argument("--coupling", action="store_true", help="A1: re-enable the continuous node-FACE bilinear adhesion (SimuCell3D-style) so cells flatten into a real tissue (regime II); without it only sparse cadherin point-bonds adhere → round cells stay round")
    ap.add_argument("--adh-strength", type=float, default=1.0e7, help="A1: node-FACE cohesion stress [Pa] used by --coupling (DCM-aggregation-resolved lit value ~4-5e7 → clean flattened spheroids)")
    ap.add_argument("--rep-strength", type=float, default=2.0e8, help="node-FACE contact repulsion stiffness [Pa]; soft lit value ~4e7 (MCF7 ξ̄) gives just-touching deformable contact, the stiff 2e8 default over-packs/interpenetrates")
    ap.add_argument("--nucleus", action="store_true", help="E2: deformable nucleus core (H.9 bilinear chromatin/lamin; resists cell thinning below the nuclear size)")
    ap.add_argument("--e-nuc", type=float, default=3.0e3, help="nuclear Young's modulus [Pa] (KU-3.B2.1 1-10 kPa)")
    ap.add_argument("--surface-tension", action="store_true", help="B4: membrane area-gradient surface tension (+ global area constraint if --k-area>0)")
    ap.add_argument("--ipc-dhat-factor", type=float, default=1.0, dest="ipc_dhat_factor",
                    help="IPC barrier activation gap d_hat = factor*c_rep. 1.0=IPC barrier; ~0.01=barrier "
                         "negligible = SimuCell3D-style implicit penalty-only (the IPC-vs-SimuCell3D A/B)")
    ap.add_argument("--polarize", action="store_true", help="apico-basal DIFFERENTIAL surface tension (Young-Dupre): basal faces wet (gamma - w*w_cs), apical keep gamma. Needs --surface-tension. The directional-spread lever.")
    ap.add_argument("--w-cs-polarize", type=float, default=2.85e-3, dest="w_cs_polarize", help="basal substrate-adhesion energy J/m2 for --polarize (lit MCF7 2.85e-3 -> S<0 non-wetting; >2*gamma -> S>0 spreads)")
    ap.add_argument("--gamma-surf", type=float, default=1.0e-4, help="B4 surface tension coefficient [N/m]")
    ap.add_argument("--aggregate-tension", action="store_true", dest="aggregate_tension",
                    help="AGGREGATE-level Foty-Steinberg liquid-drop surface tension: inward Laplace "
                         "ΔP=2σ/R_agg on the free (media) envelope → contracts + densifies a LOOSE aggregate "
                         "(the global compaction driver per-cell γ cannot supply). σ = --sigma-agg.")
    ap.add_argument("--sigma-agg", type=float, default=5.0e-3,
                    help="tissue surface tension σ [N/m] for --aggregate-tension (Foty-Steinberg ~1-20 mN/m; "
                         "swept as a controlled variable, not tuned). Default 5e-3 = 5 mN/m.")
    ap.add_argument("--agg-every", type=int, default=25,
                    help="refresh the aggregate centroid + R_agg every N steps (host read; drifts slowly)")
    ap.add_argument("--diff-tension", action="store_true", dest="diff_tension",
                    help="foam/DAH DIFFERENTIAL interfacial tension: cell-cell contact faces get the "
                         "reduced tension γ−w_adh/2 (w_adh=w_cs lit) so contacts SPREAD → cells facet. "
                         "The faceting mechanism (uniform γ never spreads contacts). Needs --surface-tension.")
    ap.add_argument("--conservative-contact", action="store_true", dest="conservative_contact",
                    help="THE ENERGY FIX: replace the non-conservative constant-force node-face contact "
                         "+ node-node cohesion + diff-γ with ONE conservative SimuCell3D bilinear "
                         "traction-separation tent (adh=4·w_cs/c_adh² derived). Gives a true faceted "
                         "energy minimum to settle into. Use with --integrator baoab (kT=0 descent).")
    ap.add_argument("--k-edge", type=float, default=-1.0, dest="k_edge",
                    help="cortex edge-spring stiffness [N/m] override (default ResolvedDCM 1e-3). A "
                         "softer cortex lets adhesion FLATTEN cells into facets; a stiff one keeps "
                         "them round. The cell-deformability faceting lever.")
    ap.add_argument("--rep-over-adh", type=float, default=4.0, dest="rep_over_adh",
                    help="conservative-contact repulsion stiffness as a multiple of the derived adhesion "
                         "(stiffer non-penetration; ξ≥ω).")
    ap.add_argument("--contact-tension-frac", type=float, default=-1.0, dest="contact_tension_frac",
                    help="Maître cortex-dissolution limit: contact-face tension = frac·γ_free (cadherin "
                         "engagement disassembles the cortex at contacts; 0=full dissolution). Overrides "
                         "the w_adh/2 reduction. <0 = use the derived w_adh/2. Needs --diff-tension.")
    ap.add_argument("--k-area", type=float, default=0.0, help="B4 global area-constraint stiffness [N/m] (0=off)")
    ap.add_argument("--division", action="store_true", help="C7: rim-cell proliferation (parked cell pool → daughters)")
    ap.add_argument("--div-pool-factor", type=float, default=1.0, help="C7 parked cell pool size as a fraction of n-cells")
    ap.add_argument("--div-rate", type=float, default=0.04, help="C7 per-rim-cell division probability per tick (ref dcm_active 0.04; higher masks the A/A0 traction band)")
    ap.add_argument("--bending", action="store_true", help="B5: thin-plate biharmonic membrane bending (Helfrich-like)")
    ap.add_argument("--k-bend", type=float, default=1.0e-5, help="B5 discrete bending stiffness [N/m]")
    ap.add_argument("--necrosis", action="store_true", help="C8: 3-zone depth necrosis (O2-proxy; softens core turgor, gates division to the rim)")
    ap.add_argument("--builder", default="fcc", choices=["cubic", "fcc", "voronoi", "sphere", "confluent"], help="D10: spheroid cell packing. fcc=isotropic close-pack of free icospheres; voronoi=Lloyd CVT centres; confluent=Voronoi-warped space-filling faceted foam (Path B, auto --v0-from-init)")
    ap.add_argument("--integrator", default="baoab", choices=["baoab", "implicit"], help="I-opt: time integrator (implicit = IMEX linearly-implicit, unlocks larger accel-dt)")
    ap.add_argument("--accel-dt", type=float, default=None, help="I-opt: larger dt for --integrator implicit (accuracy-bound; e.g. 100× the explicit dt)")
    ap.add_argument("--cg-maxiter", type=int, default=80, help="I-opt: max CG iterations per implicit step")
    ap.add_argument("--frozen-neighbors", action="store_true", dest="frozen_neighbors",
                    help="I-opt #2: cache the cohesion/contact neighbour set once per step and reuse it "
                         "in every CG matvec (removes the per-iteration hash-grid query). Parity-safe: "
                         "neighbours are frozen across a solve, so the implicit dx is byte-identical.")
    ap.add_argument("--precond-diag", action="store_true", dest="precond_diag",
                    help="I-opt #1: analytic-diagonal Jacobi preconditioner (z=r/diagA each CG iter, "
                         "diagA = γ/dt + k_edge + rep·area + turgor-vol). Same converged dx; helps only "
                         "when the operator diagonal is heterogeneous (a-dominated DCM regimes see no win).")
    ap.add_argument("--k-vol", type=float, default=7.73e5, dest="k_vol",
                    help="osmotic/bulk modulus K [Pa] (p=-K·ln(V/V0)). Sets the faceting group "
                         "γ̃=γ/(K·ℓ): the driver default 7.73e5 (single-cell-spread tuning) keeps "
                         "γ̃≪band; SimuCell3D faceting (Fischer-Friedrich K≈2.5e3) needs the soft K. "
                         "Exposed for the γ controlled-variable sweep (dcm.gamma_sweep).")
    ap.add_argument("--no-wetting", action="store_true", help="disable substrate wetting (control)")
    ap.add_argument("--no-well", action="store_true", help="disable substrate z-well (control)")
    ap.add_argument("--ubottom", action="store_true", help="ULA U-bottom: confine cells in a non-adhesive hemispherical bowl (independent of the flat well)")
    ap.add_argument("--ubottom-r-factor", type=float, default=1.35, dest="ubottom_r_factor", help="bowl radius = factor × cluster radius (default 1.35)")
    ap.add_argument("--ubottom-k", type=float, default=0.0, dest="ubottom_k", help="bowl wall stiffness; 0 ⇒ auto (= k_floor rigid-dish scale)")
    ap.add_argument("--no-grid", action="store_true", help="brute-force kernels (parity ref; slow at scale)")
    ap.add_argument("--save-frames", default=None, help="npz path to save per-frame mesh geometry (pos+faces+cof) for surface viz")
    ap.add_argument("--init-npz", default=None, dest="init_npz", help="restart from a saved aggregate npz (frames/faces/cof) instead of building a fresh ball (exposes the existing restart path, lines ~312)")
    ap.add_argument("--v0-from-init", action="store_true", dest="v0_from_init", help="set each cell's osmotic rest volume V0 to its ACTUAL volume in the loaded --init-npz mesh (confluent Voronoi-cell rest volume) instead of the free-sphere V0; removes the confluent over-inflation. Derived from geometry, not tuned.")
    ap.add_argument("--inset", type=float, default=0.04, dest="inset", help="--builder confluent: radial inset ε of each warped cell from its Voronoi boundary (lower = tighter/sharper foam, less gap). Geometry-quality knob, not physics.")
    ap.add_argument("--lloyd-iters", type=int, default=8, dest="lloyd_iters", help="--builder confluent: Lloyd CVT relaxation iterations for the seed points (more = more equiaxed/regular cells).")
    args = ap.parse_args()
    import json
    out = run_decohesion(
        n_cells=args.n_cells, subdiv=args.subdiv, steps=args.steps, frames=args.frames,
        device=args.device, dt=args.dt, warmup=args.warmup, settle_steps=args.settle_steps,
        settle_frames=args.settle_frames, gap=args.gap, k_vol=args.k_vol,
        remesh_period=args.remesh_period, pool_factor=args.pool_factor,
        edge_edge=args.edge_edge, cfl_limit=args.cfl_limit, max_substeps=args.max_substeps,
        cadherin=args.cadherin, ecm_clutch=args.ecm_clutch, cad_batch=args.cad_batch,
        cad_subcycle=args.cad_subcycle, cad_micro_M=args.cad_micro_M,
        cad_bundle=args.cad_bundle, cad_contract=args.cad_contract, cad_rbind=args.cad_rbind,
        ecm_bundle=args.ecm_bundle,
        gravity=args.gravity, delta_rho=args.delta_rho, coupling=args.coupling,
        pen_cap=args.pen_cap, pen_cap_frac=args.pen_cap_frac,
        adh_strength=args.adh_strength, rep_strength=args.rep_strength, ecm_ligand=args.ligand_density,
        nucleus=args.nucleus, E_nuc=args.e_nuc,
        surface_tension=args.surface_tension, gamma_surf=args.gamma_surf, k_area=args.k_area,
        diff_tension=args.diff_tension, contact_tension_frac=args.contact_tension_frac,
        aggregate_tension=args.aggregate_tension, sigma_agg=args.sigma_agg, agg_every=args.agg_every,
        conservative_contact=args.conservative_contact, rep_over_adh=args.rep_over_adh,
        k_edge=args.k_edge,
        polarize=args.polarize, w_cs_polarize=args.w_cs_polarize, ipc_dhat_factor=args.ipc_dhat_factor,
        division=args.division, div_pool_factor=args.div_pool_factor, div_rate=args.div_rate,
        bending=args.bending, k_bend=args.k_bend, necrosis=args.necrosis, builder=args.builder,
        integrator=args.integrator, accel_dt=args.accel_dt, cg_maxiter=args.cg_maxiter,
        frozen_neighbors=args.frozen_neighbors, precond_diag=args.precond_diag,
        substrate_wetting=not args.no_wetting, use_substrate_well=not args.no_well,
        ubottom=args.ubottom, ubottom_r_factor=args.ubottom_r_factor, ubottom_k=args.ubottom_k,
        lamellipodium=args.lamellipodium, lamel_clutch=args.lamel_clutch, filopodia=args.filopodia, junction_switch=args.junction_switch,
        use_grid=not args.no_grid, save_frames=args.save_frames, init_npz=args.init_npz, v0_from_init=args.v0_from_init,
        inset=args.inset, lloyd_iters=args.lloyd_iters,
        ipc=args.ipc, ipc_eta=args.ipc_eta,
        ipc_newton=args.ipc_newton, ipc_newton_max=args.ipc_newton_max, ipc_newton_tol=args.ipc_newton_tol)
    print(json.dumps({k: v for k, v in out.items() if k != "trajectory"}, indent=2))


if __name__ == "__main__":
    main()

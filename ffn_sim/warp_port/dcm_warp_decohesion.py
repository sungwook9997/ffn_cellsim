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

    python -m ffn_sim.warp_port.dcm_warp_decohesion --device cuda:0 --n-cells 100 --steps 60000
"""

from __future__ import annotations

import argparse
import time

import numpy as np

import warp as wp

from ffn_sim.cell.dcm import icosphere_mesh, ResolvedDCM
from ffn_sim.warp_port.dcm_warp_hybrid_multicell import (
    build_multicell, _dp_from_vol, _zero_vec, _edges_from_faces)
from ffn_sim.warp_port.dcm_warp_hybrid import _bond_accumulate, _bd_step
from ffn_sim.warp_port.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel
from ffn_sim.warp_port.dcm_cohesion_warp import dcm_cohesion_kernel
from ffn_sim.warp_port.dcm_contact_warp import node_face_contact_kernel
from ffn_sim.warp_port.dcm_substrate_warp import (
    dcm_substrate_well_accum_kernel, dcm_wetting_scatter_kernel, dcm_wetting_cap_add_kernel,
    dcm_wetting_scatter_integrin_kernel)
from ffn_sim.warp_port.dcm_neighbor_warp import (
    pos_to_f32, face_centroids_f32, cohesion_grid_kernel, contact_grid_kernel,
    cohesion_grid_cad_kernel, contact_grid_cad_kernel, penetration_depth_kernel,
    gather_lead_pos, lamellipodium_tether_multicell)
from ffn_sim.warp_port.dcm_lamellipodium_host import LamellipodiumHost, LamelParams
from ffn_sim.warp_port.dcm_junction_switch_host import JunctionSwitchHost, JunctionParams
from ffn_sim.cell.dcm_remesh import remesh_pass

wp.init()


def _spherical_centers(n_cells: int, R: float, gap: float) -> np.ndarray:
    """``n_cells`` lattice centres packed in a roughly SPHERICAL cluster (not a cube).

    A cubic lattice is too symmetric to AGGREGATE: an interior cell is pulled equally by
    neighbours on all sides → net force ≈ 0 → the pack is metastable and never compacts
    (the lattice just persists through settle — the failure the spread montage showed).
    Selecting the ``n_cells`` lattice points CLOSEST to the centre instead gives a sphere
    whose surface cells feel a net INWARD cohesive pull (curvature) → the cluster rounds
    and compacts during the settle/aggregate phase, which is the whole point of that phase.
    ``gap`` is the centre spacing in units of R."""
    side = int(np.ceil((2.2 * n_cells) ** (1.0 / 3.0))) + 2     # oversample the lattice
    g = np.arange(side) - (side - 1) / 2.0
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    keep = np.argsort(np.linalg.norm(pts, axis=1))[:n_cells]    # the n_cells nearest centre
    centers = pts[keep] * (gap * R)
    return centers - centers.mean(0)


PARK_POS = np.array([1.0e-2, 1.0e-2, 1.0e-2])   # dormant pool node home (≫ cell scale; well=0)


def build_cleanball_on_substrate(n_cells: int, subdiv: int, R: float, z0: float = 0.0,
                                 gap: float = 2.05, pool_factor: float = 0.0):
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
    centers = _spherical_centers(n_cells, R, gap)
    pos = np.concatenate([verts1 + c for c in centers], axis=0)
    faces = np.concatenate([tris1 + ci * npc for ci in range(n_cells)], axis=0)
    edges = np.concatenate([edges1 + ci * npc for ci in range(n_cells)], axis=0)
    cof = np.repeat(np.arange(n_cells), npc).astype(np.int64)
    face_cell = np.repeat(np.arange(n_cells), tris1.shape[0]).astype(np.int64)
    pos[:, 2] += (z0 - pos[:, 2].min())          # rest the ball on the dish (active nodes only)

    # DORMANT NODE POOL for remeshing (A1): SPLIT activates a parked node, COLLAPSE returns
    # one. Parked at PARK_POS (≫ cell scale, +z) so every force is 0 there (well: |dz|≫rng and
    # above z0 → 0; cohesion/contact: cof<0 skipped; not in any face/edge). n_pool extra rows.
    n_pool = int(pool_factor * pos.shape[0])
    if n_pool > 0:
        pool = np.tile(PARK_POS, (n_pool, 1))
        pos = np.concatenate([pos, pool], axis=0)
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
                   force_cap: float = 5.0e-8, z0: float = 0.0, warmup: int = 1000,
                   settle_steps: int = 0, settle_frames: int = 0,
                   remesh_period: int = 0, pool_factor: float = 0.5,
                   use_grid: bool = True, save_frames: str | None = None,
                   lamellipodium: bool = False, junction_switch: bool = False) -> dict:
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
    R = p.R_cell
    pos_a, edges_a, faces_a, cof_a, fcell_a, npc = build_cleanball_on_substrate(
        n_cells, subdiv, R, z0, gap=gap, pool_factor=(pool_factor if remesh_period else 0.0))
    N = pos_a.shape[0]
    mean_edge = float(np.linalg.norm(pos_a[edges_a[:, 0]] - pos_a[edges_a[:, 1]], axis=1).mean())
    R0 = float(np.linalg.norm(icosphere_mesh(R, subdiv)[0], axis=1).mean())
    V0 = (4.0 / 3.0) * np.pi * R0 ** 3
    area_per_node = 4.0 * np.pi * R0 ** 2 / npc
    # z-well depth from the adhesion energy density × node area (derived, not tuned)
    W_cs_well = w_cs_jm2 * area_per_node
    k_well = 2.0 * W_cs_well / (adh_range ** 2)

    # PHYSIOLOGICAL per-node drag (physiological-baseline HARD rule) — derived from the
    # MCF7 cytoplasm viscosity exactly as the HOOMD driver (dcm_gpu_build.py:576): the
    # cell's Stokes drag 6π·η·R distributed over its nodes. NOT the legacy ResolvedDCM
    # gamma_node=3.9e-10 (~5.7e5× too small → would force a brittle ~5e-8 dt + blow up at 8e-6).
    eta_cytoplasm_Pas = 65.9          # MCF7 cytoplasm (Dessard 2024)
    gamma_node = 6.0 * np.pi * eta_cytoplasm_Pas * R / npc
    inv_gamma = 1.0 / gamma_node
    c_rep = 0.30 * mean_edge
    c_adh = 0.80 * mean_edge
    r_contact = 0.30 * mean_edge

    pos_d = wp.array(np.ascontiguousarray(pos_a), dtype=wp.vec3d, device=device)
    cof_d = wp.array(cof_a.astype(np.int32), dtype=wp.int32, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    wbuf_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    Vc_d = wp.zeros(n_cells, dtype=wp.float64, device=device)
    dP_d = wp.zeros(n_cells, dtype=wp.float64, device=device)
    cad_d = wp.ones(n_cells, dtype=wp.float64, device=device)
    integrin_d = wp.ones(n_cells, dtype=wp.float64, device=device)
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
    node_f32 = wp.zeros(N, dtype=wp.vec3, device=device)
    cent_f32 = wp.zeros(n_faces, dtype=wp.vec3, device=device)
    pen_d = wp.zeros(N, dtype=wp.float64, device=device)   # interpenetration diagnostic
    node_grid = wp.HashGrid(48, 48, 48, device=device) if use_grid else None
    face_grid = wp.HashGrid(48, 48, 48, device=device) if use_grid else None

    # M2 lamellipodium host (rim detection one-shot at build; advances at cadence).
    lam = None
    if lamellipodium:
        lam = LamellipodiumHost(pos0=pos_a, cof=cof_a, n_cells=n_cells, z0=z0, R=R, dt=dt)
        print(f"  [lamel] rim cells={lam.n_rim}/{n_cells}  pool={lam.n_pool}  "
              f"p_advance={lam.p_advance:.3e}  z_basal={lam.z_basal*1e6:.3f}um", flush=True)

    # M3 junction switch host (crowd-pressure cadherin→integrin clutch; latches at cadence).
    js = None
    if junction_switch:
        js = JunctionSwitchHost(cof=cof_a, n_cells=n_cells, R=R)
        print(f"  [junction] r_contact={js.r_contact*1e6:.1f}um  P_switch={js.p.P_switch_kPa}kPa  "
              f"cad_weak={js.p.cadherin_weak_factor}  integrin_strong={js.p.integrin_strong_factor}  "
              f"cadence={js.cadence}", flush=True)

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
        nonlocal faces_d, fcell_d, edges_d, r0_d, cent_f32
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
        if lam is not None:
            lam.cof = cof_a            # activated/collapsed nodes changed the node→cell map
        if js is not None:
            js.cof = cof_a

    def step_once(s, dt_step, do_spread=True):
        wp.launch(_zero_vec, dim=N, inputs=[force_d], device=device)
        if use_grid:
            # rebuild both grids every step (build is negligible — the path is query-bound;
            # the persistent-grid lever was an honest negative, so every-step grid is fastest)
            wp.launch(pos_to_f32, dim=N, inputs=[pos_d, node_f32], device=device)
            node_grid.build(points=node_f32, radius=coh_q)
            wp.launch(face_centroids_f32, dim=n_faces, inputs=[pos_d, faces_d, cent_f32], device=device)
            face_grid.build(points=cent_f32, radius=con_q)
            if js is not None:
                # M3: cell-cell adhesion scaled by sqrt(cad_i·cad_j) (identity until a
                # cell crowd-switches, so settle behaviour is unchanged)
                wp.launch(cohesion_grid_cad_kernel, dim=N,
                          inputs=[node_grid.id, node_f32, pos_d, cof_d, cad_d, wp.float32(coh_q),
                                  wp.float64(r_contact), wp.float64(c_adh), wp.float64(rep_strength),
                                  wp.float64(adh_strength), wp.float64(area_per_node),
                                  wp.float64(force_cap), force_d], device=device)
            else:
                wp.launch(cohesion_grid_kernel, dim=N,
                          inputs=[node_grid.id, node_f32, pos_d, cof_d, wp.float32(coh_q),
                                  wp.float64(r_contact), wp.float64(c_adh), wp.float64(rep_strength),
                                  wp.float64(adh_strength), wp.float64(area_per_node),
                                  wp.float64(force_cap), force_d], device=device)
        else:
            wp.launch(dcm_cohesion_kernel, dim=N,
                      inputs=[pos_d, cof_d, cad_d, wp.int32(0), wp.int32(N),
                              wp.float64(r_contact), wp.float64(c_adh), wp.float64(rep_strength),
                              wp.float64(adh_strength), wp.float64(area_per_node),
                              wp.float64(force_cap), force_d], device=device)
        Vc_d.zero_()
        wp.launch(dcm_volume_kernel, dim=n_faces, inputs=[pos_d, faces_d, fcell_d, Vc_d], device=device)
        wp.launch(_dp_from_vol, dim=n_cells,
                  inputs=[Vc_d, wp.float64(V0), wp.float64(p.turgor_dP0),
                          wp.float64(k_vol), dP_d], device=device)
        wp.launch(dcm_turgor_force_kernel, dim=n_faces,
                  inputs=[pos_d, faces_d, fcell_d, dP_d, force_d], device=device)
        if use_grid:
            if js is not None:
                wp.launch(contact_grid_cad_kernel, dim=N,
                          inputs=[face_grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d, cad_d,
                                  wp.float32(con_q), wp.float64(rep_strength),
                                  wp.float64(adh_strength), wp.float64(c_rep), wp.float64(c_adh),
                                  force_d], device=device)
            else:
                wp.launch(contact_grid_kernel, dim=N,
                          inputs=[face_grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d,
                                  wp.float32(con_q), wp.float64(rep_strength),
                                  wp.float64(adh_strength), wp.float64(c_rep), wp.float64(c_adh),
                                  force_d], device=device)
        else:
            wp.launch(node_face_contact_kernel, dim=N,
                      inputs=[pos_d, cof_d, faces_d, fcell_d, cad_d, wp.int32(0), wp.int32(n_faces),
                              wp.float64(rep_strength), wp.float64(adh_strength),
                              wp.float64(c_rep), wp.float64(c_adh), force_d], device=device)
        wp.launch(_bond_accumulate, dim=n_edges,
                  inputs=[pos_d, edges_d, wp.float64(p.k_edge), r0_d, force_d], device=device)
        # substrate z-well (own-row accumulate) — pins basal nodes at z0
        if use_substrate_well:
            wp.launch(dcm_substrate_well_accum_kernel, dim=N,
                      inputs=[pos_d, wp.float64(z0), wp.float64(k_well), wp.float64(adh_range),
                              wp.float64(k_floor), force_d], device=device)
        # substrate in-plane wetting: scatter into wbuf -> cap -> add (mechanistic spread)
        if substrate_wetting and do_spread:
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
        wp.launch(_bd_step, dim=N,
                  inputs=[pos_d, force_d, cof_d, wp.float64(inv_gamma), wp.float64(0.0),
                          wp.float64(dt_step), wp.int32(7), wp.int32(s)], device=device)

    def _penetration_frac():
        """max node-into-other-cell penetration depth / mean_edge (0 = no interpenetration).
        Rebuilds the face grid on the CURRENT positions (measure runs after the integration
        step, so the step's grid is stale) then reduces the diagnostic kernel on the host."""
        if not use_grid:
            return 0.0
        wp.launch(pos_to_f32, dim=N, inputs=[pos_d, node_f32], device=device)
        wp.launch(face_centroids_f32, dim=n_faces, inputs=[pos_d, faces_d, cent_f32], device=device)
        face_grid.build(points=cent_f32, radius=con_q)
        pen_d.zero_()
        wp.launch(penetration_depth_kernel, dim=N,
                  inputs=[face_grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d,
                          wp.float32(con_q), pen_d], device=device)
        wp.synchronize_device(device)
        return float(pen_d.numpy().max()) / mean_edge

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
        return {"A_um2": A, "maxZ_um": maxZ, "Vsum": Vsum, "com": com, "P": Psave,
                "pen_frac": _penetration_frac()}, True

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
    recs = []   # {gstep, phase, area, maxZ, Vsum, com}

    def cad_now():
        return js.cad_mult.copy() if js is not None else np.ones(n_cells)

    def record(gstep, phase, m):
        recs.append({"gstep": gstep, "phase": phase, "area": m["A_um2"],
                     "maxZ_um": m["maxZ_um"], "Vsum": m["Vsum"], "com": m["com"],
                     "pen_frac": m.get("pen_frac", 0.0)})
        if frame_list is not None:
            frame_list.append(m["P"].astype(np.float32))
            cad_list.append(cad_now().astype(np.float32))
            faces_list.append(faces_a.astype(np.int32))   # snapshot — remesh mutates faces_a
            cof_list.append(cof_a.astype(np.int32))

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
        step_once(s, dt, do_spread=False)
        if s % every_s == 0:
            wp.synchronize_device(device)
            m, ok = measure()
            if ok:
                record(s, 0, m)
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
            wp.synchronize_device(device)
            lam.update(pos_d.numpy().astype(np.float64))
            lam.upload(device)
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
        step_once(s, dt)
        if s % every == 0 or s == steps:
            wp.synchronize_device(device)
            m, ok = measure()
            if not ok:
                truncated_at = s
                print(f"  [decoh] NON-FINITE at step {s} — truncating", flush=True)
                break
            record(settle_steps + s, 1, m)
            print(f"  step {s:>7}  A/A0={m['A_um2']/A0 if A0 else 0:.3f}  maxZ={m['maxZ_um']:.1f}um  "
                  f"V/V0={m['Vsum']/V0sum if V0sum else 0:.3f}  "
                  f"drift={np.linalg.norm(m['com']-com0)*1e6:.2f}um  "
                  f"pen={m['pen_frac']:.3f}"
                  + (f"  switched={js.n_switched}" if js is not None else ""), flush=True)
    wp.synchronize_device(device)
    elapsed = time.perf_counter() - t0

    # normalise all recorded frames to the rested baseline + assemble the trajectory
    traj = [{"step": r["gstep"], "phase": r["phase"], "aa0": r["area"] / A0 if A0 > 0 else 0.0,
             "maxZ_um": r["maxZ_um"], "vv0": r["Vsum"] / V0sum if V0sum else 0.0,
             "drift_um": float(np.linalg.norm(r["com"] - com0) * 1e6),
             "pen_frac": r["pen_frac"]} for r in recs]

    if save_frames and frame_list is not None:
        extra = {}
        if remesh_period:
            # topology varies across frames → save per-frame faces/cof as object sequences
            extra["faces_seq"] = np.array(faces_list, dtype=object)
            extra["cof_seq"] = np.array(cof_list, dtype=object)
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
        "lamellipodium": lamellipodium, "junction_switch": junction_switch,
        "A0_rested_um2": A0, "compaction_x": (m_init["A_um2"] / A0) if A0 else 0.0,
        "aa0_peak": max(spread_aa), "aa0_final": spread_aa[-1],
        "maxZ_final_um": traj[-1]["maxZ_um"],
        "vv0_final": traj[-1]["vv0"], "drift_final_um": traj[-1]["drift_um"],
        "pen_frac_peak": max(r["pen_frac"] for r in traj),
        "pen_frac_final": traj[-1]["pen_frac"],
        "remesh_period": remesh_period, "remesh": remesh_stats if remesh_period else None,
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
    return out


def main():
    ap = argparse.ArgumentParser(description="De-cohesion cleanball spread on the Warp DCM engine (M1: substrate).")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--n-cells", type=int, default=12)
    ap.add_argument("--subdiv", type=int, default=2)
    ap.add_argument("--steps", type=int, default=40000)
    ap.add_argument("--frames", type=int, default=20)
    ap.add_argument("--dt", type=float, default=8.0e-6)
    ap.add_argument("--warmup", type=int, default=1000, help="soft-start steps at 0.1x dt")
    ap.add_argument("--settle-steps", type=int, default=0,
                    help="aggregation/settle steps at full dt with NO spread drivers (rest the spheroid at z0 before the measured spread; baseline A0 is taken AFTER this)")
    ap.add_argument("--settle-frames", type=int, default=0,
                    help="number of frames to capture DURING the aggregate/settle phase (shows the cube→compact-ball compaction in the montage)")
    ap.add_argument("--lamellipodium", action="store_true", help="enable the M2 per-cell lamellipodium crawl")
    ap.add_argument("--junction-switch", action="store_true", help="enable the M3 crowd-pressure cadherin→integrin junction switch")
    ap.add_argument("--gap", type=float, default=2.05, help="cell centre spacing in R for the spherical aggregate (2.05 = touching/compact)")
    ap.add_argument("--remesh-period", type=int, default=0, help="A1: host SWAP/SPLIT/COLLAPSE remesh every N steps (0=off); keeps edges in band → no slivers")
    ap.add_argument("--pool-factor", type=float, default=0.5, help="dormant node pool size as a fraction of active nodes (for remesh SPLIT)")
    ap.add_argument("--no-wetting", action="store_true", help="disable substrate wetting (control)")
    ap.add_argument("--no-well", action="store_true", help="disable substrate z-well (control)")
    ap.add_argument("--no-grid", action="store_true", help="brute-force kernels (parity ref; slow at scale)")
    ap.add_argument("--save-frames", default=None, help="npz path to save per-frame mesh geometry (pos+faces+cof) for surface viz")
    args = ap.parse_args()
    import json
    out = run_decohesion(
        n_cells=args.n_cells, subdiv=args.subdiv, steps=args.steps, frames=args.frames,
        device=args.device, dt=args.dt, warmup=args.warmup, settle_steps=args.settle_steps,
        settle_frames=args.settle_frames, gap=args.gap,
        remesh_period=args.remesh_period, pool_factor=args.pool_factor,
        substrate_wetting=not args.no_wetting, use_substrate_well=not args.no_well,
        lamellipodium=args.lamellipodium, junction_switch=args.junction_switch,
        use_grid=not args.no_grid, save_frames=args.save_frames)
    print(json.dumps({k: v for k, v in out.items() if k != "trajectory"}, indent=2))


if __name__ == "__main__":
    main()

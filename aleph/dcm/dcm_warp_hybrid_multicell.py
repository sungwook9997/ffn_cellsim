"""Multi-cell DCM Warp HYBRID — GPU-resident per-step (turgor + edges + cohesion +
node-face contact) + host low-cadence remesh, on a shared node-pool.

Extends ``dcm_warp_hybrid`` (single cell) to a small spheroid: several icosphere
cells packed together, with the inter-cell forces wired into the SAME device-
resident per-step loop:

  per step (GPU, no host sync):
    force.zero  -> cohesion (node-node, own-row write)  -> exact per-cell turgor
    (volume reduce -> per-cell dP -> per-face force, atomic) -> node-face contact
    (atomic) -> cortex edge springs (atomic) -> overdamped Langevin step
  every remesh_period (HOST): remesh_pass per cell on the shared pool -> resync.

All four force kernels are the COMMITTED, parity-verified Warp kernels
(``dcm_turgor_warp`` / ``dcm_cohesion_warp`` / ``dcm_contact_warp`` / the cortex
edge spring), launched directly on persistent device arrays. cell_of_node uses the
SAME convention as the contact/cohesion kernels (−1 = dormant pool slot), so the
node-pool and the inter-cell forces compose without translation.

    python -m aleph.dcm.dcm_warp_hybrid_multicell --device cuda:0 --n-cells 4 --steps 2000
"""

from __future__ import annotations

import argparse
import time

import numpy as np

import warp as wp

from aleph.dcm.geometry import icosphere_mesh, ResolvedDCM
from aleph.dcm.dcm_remesh import remesh_pass
from aleph.dcm.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel
from aleph.dcm.dcm_cohesion_warp import dcm_cohesion_kernel
from aleph.dcm.dcm_contact_warp import node_face_contact_kernel
from aleph.dcm.dcm_warp_hybrid import _bond_accumulate, _bd_step, _edges_from_faces
from aleph.dcm.dcm_neighbor_warp import (
    pos_to_f32, face_centroids_f32, cohesion_grid_kernel, contact_grid_kernel)

wp.init()


@wp.kernel
def _dp_from_vol(Vc: wp.array(dtype=wp.float64), V0: wp.float64,
                 dP0: wp.float64, K_vol: wp.float64,
                 dP_cell: wp.array(dtype=wp.float64)):
    c = wp.tid()
    dP_cell[c] = dP0 + K_vol * (V0 - Vc[c]) / V0


@wp.kernel
def _dp_from_vol_pc(Vc: wp.array(dtype=wp.float64), V0: wp.array(dtype=wp.float64),
                    dP0: wp.float64, K_vol: wp.float64,
                    dP_cell: wp.array(dtype=wp.float64)):
    """Per-cell rest-volume turgor: dP_c = dP0 + K_vol·(V0_c − Vc_c)/V0_c.

    Identical to :func:`_dp_from_vol` but each cell carries its OWN rest volume ``V0[c]`` — needed
    for cell division (a daughter is born at half rest-volume and re-grows its osmotic setpoint to
    V0 over the cell cycle, so turgor inflates it gradually instead of a full-size cold insert)."""
    c = wp.tid()
    dP_cell[c] = dP0 + K_vol * (V0[c] - Vc[c]) / V0[c]


@wp.kernel
def _dp_from_vol_osm(Vc: wp.array(dtype=wp.float64), V0: wp.array(dtype=wp.float64),
                     V0_ref: wp.float64, dP0: wp.float64, K_vol: wp.float64,
                     dP_cell: wp.array(dtype=wp.float64)):
    """OSMOTIC turgor with concentration feedback (PI 2026-06-25; KB-3.9 water-flux volume regulation).
    The resting osmotic excess ``dP0`` is the van't-Hoff pressure at the reference volume ``V0_ref``;
    when the cell's osmotic setpoint ``V0[c]`` relaxes DOWN (water efflux at the media-exposed surface,
    handled by the host osmotic relaxation), the impermeant solute concentrates → the resting pressure
    rises as ``dP0·V0_ref/V0[c]``. This gives a STABLE volume equilibrium (the concentration resists
    further efflux) instead of the elastic spring holding a permanent V0−Vc strain (which accumulated
    stiffness → the faceting-over-compression CG stall). The fast K_vol term (minute-scale
    incompressibility) is retained; the slow osmotic adaptation is the V0[c] relaxation."""
    c = wp.tid()
    dP_cell[c] = dP0 * (V0_ref / V0[c]) + K_vol * (V0[c] - Vc[c]) / V0[c]


@wp.kernel
def _zero_vec(force: wp.array(dtype=wp.vec3d)):
    force[wp.tid()] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))


@wp.kernel
def _edge_extent(pos: wp.array(dtype=wp.vec3d),
                 edges: wp.array(dtype=wp.int32, ndim=2),
                 out: wp.array(dtype=wp.float64)):
    """out[0] = max edge length, out[1] = -(min edge length) (device reduction).

    Lets the host remesh be SKIPPED entirely when every edge is in [l_min, l_max]
    — the equilibrium common case — so a remesh epoch costs one tiny reduction +
    a 2-float read instead of a full-pool download + python remesh_pass."""
    e = wp.tid()
    d = wp.length(pos[edges[e, 0]] - pos[edges[e, 1]])
    wp.atomic_max(out, 0, d)
    wp.atomic_max(out, 1, -d)


def build_multicell(n_cells: int, subdiv: int, R: float, gap: float = 1.7):
    """Pack ``n_cells`` icospheres on a cubic-ish lattice; return pooled arrays."""
    verts1, edges1, tris1 = icosphere_mesh(R, subdiv)
    npc = verts1.shape[0]
    # lattice centers, spacing gap*R (slightly less than 2R so shells nearly touch)
    side = int(np.ceil(n_cells ** (1.0 / 3.0)))
    centers = []
    for i in range(side):
        for j in range(side):
            for k in range(side):
                centers.append((i, j, k))
                if len(centers) == n_cells:
                    break
            if len(centers) == n_cells:
                break
        if len(centers) == n_cells:
            break
    centers = np.array(centers, float) * (gap * R)
    centers -= centers.mean(0)

    pos = np.concatenate([verts1 + c for c in centers], axis=0)
    faces = np.concatenate([tris1 + ci * npc for ci in range(n_cells)], axis=0)
    edges = np.concatenate([edges1 + ci * npc for ci in range(n_cells)], axis=0)
    cof = np.repeat(np.arange(n_cells), npc).astype(np.int64)
    face_cell = np.repeat(np.arange(n_cells), tris1.shape[0]).astype(np.int64)
    return pos, edges.astype(np.int64), faces.astype(np.int64), cof, face_cell, npc


def run_multicell(*, n_cells: int = 4, subdiv: int = 2, steps: int = 2000,
                  remesh_period: int = 0, device: str = "cpu", kT: float = 0.0,
                  dt: float = 5.0e-8, pool_factor: float = 3.0, warmup: int = 30,
                  use_contact: bool = True, use_cohesion: bool = True,
                  use_grid: bool = False, grid_period: int = 1,
                  skin_frac: float = 0.0) -> dict:
    p = ResolvedDCM(subdivisions=subdiv)
    R = p.R_cell
    pos_a, edges_a, faces_a, cof_a, fcell_a, npc = build_multicell(n_cells, subdiv, R)
    n_active0 = pos_a.shape[0]
    MAX = int(n_active0 * pool_factor)
    mean_edge = float(np.linalg.norm(pos_a[edges_a[:, 0]] - pos_a[edges_a[:, 1]], axis=1).mean())
    l_min = mean_edge / 2.9
    R0 = float(np.linalg.norm(icosphere_mesh(R, subdiv)[0], axis=1).mean())
    V0 = (4.0 / 3.0) * np.pi * R0 ** 3
    area_per_node = 4.0 * np.pi * R0 ** 2 / npc
    park_pos = np.array([1.0e-2, 1.0e-2, 1.0e-2])  # dormant: far vs ~1e-5 m cell, modest cell-index

    inv_gamma = 1.0 / p.gamma_node
    bd_pref = float(np.sqrt(2.0 * kT / p.gamma_node * dt)) if kT > 0 else 0.0

    # contact / cohesion params (soft, stable; SI-scale, physical-ish bands)
    c_rep = 0.30 * mean_edge       # repulsion range
    c_adh = 0.80 * mean_edge       # adhesion range
    r_contact = 0.30 * mean_edge   # node-node contact radius (cohesion)
    rep_strength = 2.0e2           # ξ  (soft — per SimuCell3D dt-caveat, not 1e9)
    adh_strength = 5.0e1           # ω
    force_cap = 5.0e-8

    # node pool: active first, dormant (cof=-1, parked) after
    pos_h = np.tile(park_pos, (MAX, 1)).astype(np.float64)
    pos_h[:n_active0] = pos_a
    cof = np.full(MAX, -1, dtype=np.int64)
    cof[:n_active0] = cof_a
    faces_h = faces_a.copy()
    face_cell = fcell_a.copy()
    edges_h = edges_a.astype(np.int32)
    r0_h = np.linalg.norm(pos_h[edges_h[:, 0]] - pos_h[edges_h[:, 1]], axis=1).astype(np.float64)

    pos_d = wp.array(pos_h, dtype=wp.vec3d, device=device)
    cof_d = wp.array(cof.astype(np.int32), dtype=wp.int32, device=device)
    force_d = wp.zeros(MAX, dtype=wp.vec3d, device=device)
    Vc_d = wp.zeros(n_cells, dtype=wp.float64, device=device)
    dP_d = wp.zeros(n_cells, dtype=wp.float64, device=device)
    cad_d = wp.ones(n_cells, dtype=wp.float64, device=device)  # dummy (use_cad=0)
    faces_d = wp.array(faces_h.astype(np.int32), dtype=wp.int32, device=device)
    fcell_d = wp.array(face_cell.astype(np.int32), dtype=wp.int32, device=device)
    edges_d = wp.array(edges_h, dtype=wp.int32, device=device)
    r0_d = wp.array(r0_h, dtype=wp.float64, device=device)
    n_edges = edges_h.shape[0]
    n_faces = faces_h.shape[0]

    # --- hash-grid neighbour-list acceleration (use_grid) ---
    max_edge = float(np.linalg.norm(pos_a[edges_a[:, 0]] - pos_a[edges_a[:, 1]], axis=1).max())
    coh_radius = float(c_adh)                       # node-node interaction range
    # face-centroid query radius: c_adh + the centroid->closest-point reach, bounded
    # by the triangle circumradius (~0.58*edge); use 0.7*l_max so the grid cell size
    # stays small (coarse cells = too many candidates = the slow grid). l_max=3*l_min.
    con_radius = float(c_adh + 0.7 * (3.0 * l_min))
    # PERSISTENT grid (Verlet skin): rebuild every grid_period steps; pad the query
    # radius by skin so no neighbour that drifts into range between rebuilds is missed
    # (valid while max displacement over grid_period < skin/2; verified vs every-step).
    skin = float(skin_frac * mean_edge)
    coh_q = coh_radius + skin
    con_q = con_radius + skin
    node_f32 = wp.zeros(MAX, dtype=wp.vec3, device=device)
    cent_f32 = wp.zeros(n_faces, dtype=wp.vec3, device=device)
    node_grid = wp.HashGrid(48, 48, 48, device=device) if use_grid else None
    face_grid = wp.HashGrid(48, 48, 48, device=device) if use_grid else None

    extent_d = wp.zeros(2, dtype=wp.float64, device=device)
    l_max = 3.0 * l_min

    remesh_count = {"swap": 0, "split": 0, "collapse": 0}
    n_remesh_events = 0
    n_remesh_skipped = 0
    sync_time = 0.0
    # #2: launch per-node kernels over ACTIVE nodes only. The pool is [active,
    # dormant] and SPLIT fills dormant slots in order, so active stays contiguous
    # [0:launch_n]; COLLAPSE can punch a hole -> fall back to MAX (correct, slower).
    launch_n = n_active0
    _fr = [True]  # force grid rebuild (set after remesh / topology change)

    def step_once(s):
        wp.launch(_zero_vec, dim=MAX, inputs=[force_d], device=device)
        if use_grid:
            # query points are CURRENT every step; grid buckets rebuilt at cadence
            wp.launch(pos_to_f32, dim=launch_n, inputs=[pos_d, node_f32], device=device)
            if s % grid_period == 0 or _fr[0]:
                if use_cohesion:
                    node_grid.build(points=node_f32[:launch_n], radius=coh_q)
                if use_contact:
                    wp.launch(face_centroids_f32, dim=n_faces,
                              inputs=[pos_d, faces_d, cent_f32], device=device)
                    face_grid.build(points=cent_f32, radius=con_q)
                _fr[0] = False
        # 1) cohesion (own-row write onto the zeroed force)
        if use_cohesion:
            if use_grid:
                wp.launch(cohesion_grid_kernel, dim=launch_n,
                          inputs=[node_grid.id, node_f32, pos_d, cof_d,
                                  wp.float32(coh_q), wp.float64(r_contact),
                                  wp.float64(c_adh), wp.float64(rep_strength),
                                  wp.float64(adh_strength), wp.float64(area_per_node),
                                  wp.float64(force_cap), force_d], device=device)
            else:
                wp.launch(dcm_cohesion_kernel, dim=launch_n,
                          inputs=[pos_d, cof_d, cad_d, wp.int32(0), wp.int32(launch_n),
                                  wp.float64(r_contact), wp.float64(c_adh),
                                  wp.float64(rep_strength), wp.float64(adh_strength),
                                  wp.float64(area_per_node), wp.float64(force_cap), force_d],
                          device=device)
        # 2) exact per-cell turgor: volume -> dP -> per-face force (atomic add)
        Vc_d.zero_()
        wp.launch(dcm_volume_kernel, dim=n_faces,
                  inputs=[pos_d, faces_d, fcell_d, Vc_d], device=device)
        wp.launch(_dp_from_vol, dim=n_cells,
                  inputs=[Vc_d, wp.float64(V0), wp.float64(p.turgor_dP0),
                          wp.float64(p.K_vol), dP_d], device=device)
        wp.launch(dcm_turgor_force_kernel, dim=n_faces,
                  inputs=[pos_d, faces_d, fcell_d, dP_d, force_d], device=device)
        # 3) node-face contact (atomic add)
        if use_contact:
            if use_grid:
                wp.launch(contact_grid_kernel, dim=launch_n,
                          inputs=[face_grid.id, node_f32, pos_d, cof_d, faces_d, fcell_d,
                                  wp.float32(con_q), wp.float64(rep_strength),
                                  wp.float64(adh_strength), wp.float64(c_rep),
                                  wp.float64(c_adh), force_d], device=device)
            else:
                wp.launch(node_face_contact_kernel, dim=launch_n,
                          inputs=[pos_d, cof_d, faces_d, fcell_d, cad_d, wp.int32(0),
                                  wp.int32(n_faces), wp.float64(rep_strength),
                                  wp.float64(adh_strength), wp.float64(c_rep), wp.float64(c_adh),
                                  force_d], device=device)
        # 4) cortex edge springs (atomic add)
        wp.launch(_bond_accumulate, dim=n_edges,
                  inputs=[pos_d, edges_d, wp.float64(p.k_edge), r0_d, force_d], device=device)
        # 5) overdamped Langevin step
        wp.launch(_bd_step, dim=launch_n,
                  inputs=[pos_d, force_d, cof_d, wp.float64(inv_gamma), wp.float64(bd_pref),
                          wp.float64(dt), wp.int32(7), wp.int32(s)], device=device)

    for s in range(warmup):
        step_once(s)
    wp.synchronize_device(device)

    t0 = time.perf_counter()
    for s in range(steps):
        step_once(s)
        if remesh_period > 0 and s > 0 and s % remesh_period == 0:
            ts = time.perf_counter()
            # cheap device gate: skip the host remesh unless an edge is out of bounds
            extent_d.zero_()
            wp.launch(_edge_extent, dim=n_edges, inputs=[pos_d, edges_d, extent_d], device=device)
            ext = extent_d.numpy()
            if float(ext[0]) <= l_max and -float(ext[1]) >= l_min:
                n_remesh_skipped += 1
                sync_time += time.perf_counter() - ts
                continue
            wp.synchronize_device(device)
            pos_h = pos_d.numpy().astype(np.float64)
            pos_h, faces_h, cof, face_cell, counts = remesh_pass(
                pos_h, faces_h, cof, l_min, face_cell=face_cell, max_ops=12, park=park_pos)
            for kk in remesh_count:
                remesh_count[kk] += counts[kk]
            if sum(counts[k] for k in remesh_count) > 0:
                n_remesh_events += 1
                edges_h = _edges_from_faces(faces_h)
                r0_h = np.linalg.norm(pos_h[edges_h[:, 0]] - pos_h[edges_h[:, 1]], axis=1)
                # #4: pos/cof are fixed MAX-pool size -> REUSE the buffers (assign),
                # no per-epoch realloc; only faces/edges/r0 grow -> must realloc.
                pos_d.assign(pos_h)
                cof_d.assign(cof.astype(np.int32))
                faces_d = wp.array(faces_h.astype(np.int32), dtype=wp.int32, device=device)
                fcell_d = wp.array(face_cell.astype(np.int32), dtype=wp.int32, device=device)
                edges_d = wp.array(edges_h, dtype=wp.int32, device=device)
                r0_d = wp.array(r0_h, dtype=wp.float64, device=device)
                n_edges, n_faces = edges_h.shape[0], faces_h.shape[0]
                if use_grid:
                    cent_f32 = wp.zeros(n_faces, dtype=wp.vec3, device=device)
                # #2: refresh active-launch bound; contiguous iff [0:na] all active
                na = int((cof >= 0).sum())
                launch_n = na if bool((cof[:na] >= 0).all()) else MAX
                _fr[0] = True  # topology changed -> force grid rebuild next step
            sync_time += time.perf_counter() - ts
    wp.synchronize_device(device)
    elapsed = time.perf_counter() - t0

    pf = pos_d.numpy()
    active = cof >= 0
    finite = bool(np.isfinite(pf[active]).all())
    return {
        "device": device, "n_cells": n_cells, "subdiv": subdiv,
        "n_active0": n_active0, "pool_MAX": MAX, "n_active_final": int(active.sum()),
        "n_faces": n_faces, "steps": steps, "remesh_period": remesh_period,
        "elapsed_s": elapsed, "steps_per_s": steps / elapsed,
        "finite": finite, "remesh_events": n_remesh_events,
        "remesh_skipped": n_remesh_skipped, "remesh_ops": remesh_count,
        "remesh_sync_frac": sync_time / elapsed if elapsed else 0.0, "use_grid": use_grid,
        "grid_period": grid_period, "kT": kT, "dt": dt,
        "_final_pos": pf[active].copy(),  # for persistent-vs-every-step correctness check
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--n-cells", type=int, default=4)
    ap.add_argument("--subdiv", type=int, default=2)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--remesh-period", type=int, default=0)
    ap.add_argument("--kT", type=float, default=0.0)
    args = ap.parse_args()
    import json
    r = run_multicell(n_cells=args.n_cells, subdiv=args.subdiv, steps=args.steps,
                      remesh_period=args.remesh_period, device=args.device, kT=args.kT)
    print(json.dumps(r, indent=2))


if __name__ == "__main__":
    main()


@wp.kernel
def osmotic_relax_kernel(V0: wp.array(dtype=wp.float64), Vc: wp.array(dtype=wp.float64),
                         contact_cnt: wp.array(dtype=wp.int32), cof: wp.array(dtype=wp.int32),
                         npc: wp.int32, faces_per_cell: wp.float64, rate: wp.float64):
    """GPU per-step osmotic volume relaxation (PI 2026-06-26: GPU-only, no host sync). Each ACTIVE
    cell's osmotic setpoint V0[c] relaxes toward its current volume Vc[c] by a TINY per-step amount
    ``rate·f_media·(Vc−V0)``, where f_media = 1 − contact_cnt[c]/faces_per_cell is the media-exposed
    face fraction (rim cells shed water fast, interior buffered). A small per-step rate (vs a periodic
    bulk jump) keeps the setpoint change smooth → no cfl spike from the relaxation itself."""
    c = wp.tid()
    if cof[c * npc] < wp.int32(0):
        return                                  # parked/dormant cell
    fm = wp.float64(1.0) - wp.float64(contact_cnt[c]) / faces_per_cell
    if fm < wp.float64(0.0):
        fm = wp.float64(0.0)
    V0[c] = V0[c] + rate * fm * (Vc[c] - V0[c])

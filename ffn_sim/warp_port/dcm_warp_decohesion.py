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
    dcm_substrate_well_accum_kernel, dcm_wetting_scatter_kernel, dcm_wetting_cap_add_kernel)

wp.init()


def build_cleanball_on_substrate(n_cells: int, subdiv: int, R: float, z0: float = 0.0,
                                 gap: float = 2.2):
    """Compact ball of cells RESTING on the substrate plane z0 (lowest node at z0).

    Reuses ``build_multicell`` (cubic-ish pack) then drops the whole ball so its lowest
    node sits at z0 — the basal nodes engage the substrate well/wetting. ``gap`` is the
    centre spacing in units of R; default 2.2 (shells just SEPARATED, ~0.2R apart) so the
    stiff contact repulsion does not blow up from initial overlap (gap<2 overlaps), while
    staying inside the cohesion range so the ball holds together. Returns pooled arrays + npc.
    """
    pos, edges, faces, cof, face_cell, npc = build_multicell(n_cells, subdiv, R, gap=gap)
    pos[:, 2] += (z0 - pos[:, 2].min())          # rest the ball on the dish
    return pos, edges, faces, cof, face_cell, npc


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
                   adh_range: float = 0.5e-6, k_floor: float = 1.0,
                   substrate_wetting: bool = True, use_substrate_well: bool = True,
                   force_cap: float = 5.0e-8, z0: float = 0.0, warmup: int = 1000) -> dict:
    """Cleanball de-cohesion spread on the Warp loop with substrate drivers (M1)."""
    p = ResolvedDCM(subdivisions=subdiv)
    R = p.R_cell
    pos_a, edges_a, faces_a, cof_a, fcell_a, npc = build_cleanball_on_substrate(
        n_cells, subdiv, R, z0)
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
    faces_d = wp.array(faces_a.astype(np.int32), dtype=wp.int32, device=device)
    fcell_d = wp.array(fcell_a.astype(np.int32), dtype=wp.int32, device=device)
    edges_d = wp.array(edges_a.astype(np.int32), dtype=wp.int32, device=device)
    r0_d = wp.array(np.linalg.norm(pos_a[edges_a[:, 0]] - pos_a[edges_a[:, 1]], axis=1),
                    dtype=wp.float64, device=device)
    n_edges = edges_a.shape[0]
    n_faces = faces_a.shape[0]

    def step_once(s, dt_step):
        wp.launch(_zero_vec, dim=N, inputs=[force_d], device=device)
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
        if substrate_wetting:
            wp.launch(_zero_vec, dim=N, inputs=[wbuf_d], device=device)
            wp.launch(dcm_wetting_scatter_kernel, dim=n_faces,
                      inputs=[pos_d, faces_d, wp.float64(z0), wp.float64(w_cs_jm2),
                              wp.float64(adh_range), wbuf_d], device=device)
            wp.launch(dcm_wetting_cap_add_kernel, dim=N,
                      inputs=[wbuf_d, wp.float64(force_cap), force_d], device=device)
        wp.launch(_bd_step, dim=N,
                  inputs=[pos_d, force_d, cof_d, wp.float64(inv_gamma), wp.float64(0.0),
                          wp.float64(dt_step), wp.int32(7), wp.int32(s)], device=device)

    def measure():
        P = pos_d.numpy().astype(np.float64)
        finite = bool(np.isfinite(P).all())
        if not finite:
            return None, False
        UM = 1e6
        A = _topdown_area_um2(P[:, :2] * UM)
        maxZ = float(P[:, 2].max() * UM)
        Vsum = float(_cell_volumes(P, faces_a, fcell_a, n_cells).sum())
        com = P[:, :2].mean(axis=0)
        return {"A_um2": A, "maxZ_um": maxZ, "Vsum": Vsum, "com": com, "P": P}, True

    m0, ok = measure()
    if not ok:
        return {"error": "non-finite at init"}
    A0 = m0["A_um2"]; V0sum = m0["Vsum"]; com0 = m0["com"]
    every = max(1, steps // max(1, frames))
    traj = [{"step": 0, "aa0": 1.0, "maxZ_um": m0["maxZ_um"], "vv0": 1.0, "drift_um": 0.0}]

    # gentle soft-start: settle the initial pack at 0.1× dt (the stiff turgor/contact
    # need it; the HOOMD driver equilibrates similarly before the measured spread)
    for s in range(warmup):
        step_once(s, dt * 0.1)
    wp.synchronize_device(device)

    t0 = time.perf_counter()
    truncated_at = None
    for s in range(1, steps + 1):
        step_once(s, dt)
        if s % every == 0 or s == steps:
            wp.synchronize_device(device)
            m, ok = measure()
            if not ok:
                truncated_at = s
                print(f"  [decoh] NON-FINITE at step {s} — truncating", flush=True)
                break
            rec = {"step": s, "aa0": m["A_um2"] / A0 if A0 > 0 else 0.0,
                   "maxZ_um": m["maxZ_um"], "vv0": m["Vsum"] / V0sum if V0sum else 0.0,
                   "drift_um": float(np.linalg.norm(m["com"] - com0) * 1e6)}
            traj.append(rec)
            print(f"  step {s:>7}  A/A0={rec['aa0']:.3f}  maxZ={rec['maxZ_um']:.1f}um  "
                  f"V/V0={rec['vv0']:.3f}  drift={rec['drift_um']:.2f}um", flush=True)
    wp.synchronize_device(device)
    elapsed = time.perf_counter() - t0

    aa = [r["aa0"] for r in traj]
    return {
        "device": device, "n_cells": n_cells, "N": N, "subdiv": subdiv, "dt": dt,
        "steps": steps, "truncated_at": truncated_at, "steps_per_s": (truncated_at or steps) / elapsed,
        "substrate_wetting": substrate_wetting, "use_substrate_well": use_substrate_well,
        "aa0_peak": max(aa), "aa0_final": aa[-1], "maxZ_final_um": traj[-1]["maxZ_um"],
        "vv0_final": traj[-1]["vv0"], "drift_final_um": traj[-1]["drift_um"],
        "W_cs_well_J": W_cs_well, "gamma_node": gamma_node, "trajectory": traj,
    }


def main():
    ap = argparse.ArgumentParser(description="De-cohesion cleanball spread on the Warp DCM engine (M1: substrate).")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--n-cells", type=int, default=12)
    ap.add_argument("--subdiv", type=int, default=2)
    ap.add_argument("--steps", type=int, default=40000)
    ap.add_argument("--frames", type=int, default=20)
    ap.add_argument("--dt", type=float, default=8.0e-6)
    ap.add_argument("--no-wetting", action="store_true", help="disable substrate wetting (control)")
    ap.add_argument("--no-well", action="store_true", help="disable substrate z-well (control)")
    args = ap.parse_args()
    import json
    out = run_decohesion(
        n_cells=args.n_cells, subdiv=args.subdiv, steps=args.steps, frames=args.frames,
        device=args.device, dt=args.dt,
        substrate_wetting=not args.no_wetting, use_substrate_well=not args.no_well)
    print(json.dumps({k: v for k, v in out.items() if k != "trajectory"}, indent=2))


if __name__ == "__main__":
    main()

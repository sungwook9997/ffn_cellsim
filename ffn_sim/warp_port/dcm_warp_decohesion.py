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
    cohesion_grid_cad_kernel, contact_grid_cad_kernel,
    gather_lead_pos, lamellipodium_tether_multicell)
from ffn_sim.warp_port.dcm_lamellipodium_host import LamellipodiumHost, LamelParams
from ffn_sim.warp_port.dcm_junction_switch_host import JunctionSwitchHost, JunctionParams

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
                   force_cap: float = 5.0e-8, z0: float = 0.0, warmup: int = 1000,
                   settle_steps: int = 0,
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

    m_init, ok = measure()
    if not ok:
        return {"error": "non-finite at init"}

    # gentle soft-start: settle the initial pack at 0.1× dt (the stiff turgor/contact
    # need it; the HOOMD driver equilibrates similarly before the measured spread)
    for s in range(warmup):
        step_once(s, dt * 0.1, do_spread=False)
    wp.synchronize_device(device)

    # AGGREGATION / SETTLE phase (do_spread=False): the clean ball compacts into a
    # cohesive spheroid RESTING at z0 — substrate well holds the basal nodes at z=0,
    # turgor + cohesion + contact round it out, but NO wetting / NO lamellipodium yet.
    # The de-cohesion spread is then measured FROM this rested baseline (A0 below), so
    # the cells start mutually bonded and at z=0 exactly as PI prescribed — the upper
    # (non-ECM) cells are then dragged outward by the bonded rim cells once spread fires.
    for s in range(settle_steps):
        step_once(s, dt, do_spread=False)
    wp.synchronize_device(device)

    # baseline = the RESTED spheroid (not the as-built ball)
    m0, ok = measure()
    if not ok:
        return {"error": "non-finite after settle"}
    A0 = m0["A_um2"]; V0sum = m0["Vsum"]; com0 = m0["com"]
    every = max(1, steps // max(1, frames))
    traj = [{"step": 0, "aa0": 1.0, "maxZ_um": m0["maxZ_um"], "vv0": 1.0, "drift_um": 0.0}]
    frame_list = [m0["P"].astype(np.float32)] if save_frames else None
    print(f"  [settle] done ({warmup} warmup + {settle_steps} settle): "
          f"A0={A0:.1f}um^2  maxZ={m0['maxZ_um']:.1f}um  Vsum/Vinit="
          f"{m0['Vsum']/m_init['Vsum']:.3f}", flush=True)

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
            if frame_list is not None:
                frame_list.append(m["P"].astype(np.float32))
            print(f"  step {s:>7}  A/A0={rec['aa0']:.3f}  maxZ={rec['maxZ_um']:.1f}um  "
                  f"V/V0={rec['vv0']:.3f}  drift={rec['drift_um']:.2f}um", flush=True)
    wp.synchronize_device(device)
    elapsed = time.perf_counter() - t0

    if save_frames and frame_list is not None:
        np.savez_compressed(
            save_frames, frames=np.array(frame_list, dtype=np.float32),
            faces=faces_a.astype(np.int32), cof=cof_a.astype(np.int32),
            step=np.array([r["step"] for r in traj]),
            aa0=np.array([r["aa0"] for r in traj]),
            maxZ=np.array([r["maxZ_um"] for r in traj]),
            vv0=np.array([r["vv0"] for r in traj]))
        print(f"  saved {len(frame_list)} frames -> {save_frames}", flush=True)

    aa = [r["aa0"] for r in traj]
    out = {
        "device": device, "n_cells": n_cells, "N": N, "subdiv": subdiv, "dt": dt,
        "steps": steps, "warmup": warmup, "settle_steps": settle_steps,
        "truncated_at": truncated_at, "steps_per_s": (truncated_at or steps) / elapsed,
        "substrate_wetting": substrate_wetting, "use_substrate_well": use_substrate_well,
        "lamellipodium": lamellipodium, "junction_switch": junction_switch,
        "aa0_peak": max(aa), "aa0_final": aa[-1], "maxZ_final_um": traj[-1]["maxZ_um"],
        "vv0_final": traj[-1]["vv0"], "drift_final_um": traj[-1]["drift_um"],
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
    ap.add_argument("--lamellipodium", action="store_true", help="enable the M2 per-cell lamellipodium crawl")
    ap.add_argument("--junction-switch", action="store_true", help="enable the M3 crowd-pressure cadherin→integrin junction switch")
    ap.add_argument("--no-wetting", action="store_true", help="disable substrate wetting (control)")
    ap.add_argument("--no-well", action="store_true", help="disable substrate z-well (control)")
    ap.add_argument("--no-grid", action="store_true", help="brute-force kernels (parity ref; slow at scale)")
    ap.add_argument("--save-frames", default=None, help="npz path to save per-frame mesh geometry (pos+faces+cof) for surface viz")
    args = ap.parse_args()
    import json
    out = run_decohesion(
        n_cells=args.n_cells, subdiv=args.subdiv, steps=args.steps, frames=args.frames,
        device=args.device, dt=args.dt, warmup=args.warmup, settle_steps=args.settle_steps,
        substrate_wetting=not args.no_wetting, use_substrate_well=not args.no_well,
        lamellipodium=args.lamellipodium, junction_switch=args.junction_switch,
        use_grid=not args.no_grid, save_frames=args.save_frames)
    print(json.dumps({k: v for k, v in out.items() if k != "trajectory"}, indent=2))


if __name__ == "__main__":
    main()

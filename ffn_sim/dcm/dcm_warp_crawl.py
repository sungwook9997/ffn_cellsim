"""Lamellipodium wired into a device-resident hybrid loop (optimization #6).

The lamellipodium is a SUBSTRATE-CRAWL force (not a suspended-spheroid force): a
cell on a basal plane whose leading basal-rim nodes are pulled outward toward
clutch-anchored actin beads — the emergent spreading traction. This wires the
committed Warp tether kernel into a device-resident per-step loop alongside turgor
+ cortex edge springs, on a single cell spreading on a substrate (the canonical
lamellipodium scenario, single-cell A/A0 spreading).

Composition (the integration that #6 delivers):
  * per step (device): radial turgor + edge springs + (device gather of leading-node
    positions -> ACCUMULATING tether kernel) + integrator. No per-step host sync.
  * the leading-node SET + per-node outward direction + actin anchors are clutch-
    gripped geometry, computed on the host at low cadence (here once at setup for the
    smoke); only the leading-node POSITIONS are gathered on-device each step.

Scope: this wires the lamellipodial TRACTION into the loop (the committed,
parity-verified tether force). The full clutch-anchor RATCHET dynamics (actin
advance / gripping kinetics) is a separate dynamics piece, not included here.
"""

from __future__ import annotations

import numpy as np

import warp as wp

from ffn_sim.cell.dcm import icosphere_mesh, ResolvedDCM
from ffn_sim.dcm.dcm_warp_hybrid import (
    _reduce_centroid, _reduce_radius, _turgor_write, _bond_accumulate, _bd_step)
from ffn_sim.dcm.dcm_neighbor_warp import gather_lead_pos, lamellipodium_tether_accum

wp.init()


def _leading_geometry(pos, z_basal, basal_band, lead_frac):
    """Host (low-cadence, clutch-gripped): rim+basal leading nodes + per-node radial
    outward (single-cell spreading) + outward projection. Returns lead arrays."""
    centroid = pos.mean(axis=0)
    rel = pos - centroid
    rxy = np.linalg.norm(rel[:, :2], axis=1)
    basal = np.abs(pos[:, 2] - z_basal) <= basal_band
    lead = basal & (rxy >= lead_frac * np.maximum(rxy.max(), 1e-18))
    idx = np.where(lead)[0]
    ox = np.where(rxy[idx] > 0, rel[idx, 0] / np.maximum(rxy[idx], 1e-18), 0.0)
    oy = np.where(rxy[idx] > 0, rel[idx, 1] / np.maximum(rxy[idx], 1e-18), 0.0)
    proj = rel[idx, 0] * ox + rel[idx, 1] * oy
    ccx = np.full(idx.size, centroid[0])
    ccy = np.full(idx.size, centroid[1])
    return idx.astype(np.int32), ccx, ccy, ox, oy, proj


def run_crawl(*, subdiv: int = 2, steps: int = 2000, device: str = "cpu",
              dt: float = 2.0e-8, k_tether: float = 5.0e-4, warmup: int = 30) -> dict:
    p = ResolvedDCM(subdivisions=subdiv)
    R = p.R_cell
    verts, edges, tris = icosphere_mesh(R, subdiv)
    verts = verts.copy()
    verts[:, 2] += R           # rest the cell so its bottom sits at the basal plane z=0
    z_basal = 0.0
    basal_band = 0.5 * R
    N = verts.shape[0]
    edges = edges.astype(np.int32)
    r0 = np.linalg.norm(verts[edges[:, 0]] - verts[edges[:, 1]], axis=1).astype(np.float64)
    R0 = float(np.linalg.norm(verts - verts.mean(0), axis=1).mean())
    V0 = (4.0 / 3.0) * np.pi * R0 ** 3
    inv_gamma = 1.0 / p.gamma_node

    # clutch-anchored actin ring in the basal plane, outward of the rim (radius 1.6 R)
    na = 32
    ang = np.linspace(0, 2 * np.pi, na, endpoint=False)
    cxy = verts[:, :2].mean(0)
    actin = np.stack([cxy[0] + 1.6 * R * np.cos(ang),
                      cxy[1] + 1.6 * R * np.sin(ang),
                      np.full(na, z_basal)], axis=1).astype(np.float64)
    tether_radius = 3.0 * R
    force_cap = 5.0e-9

    lead_idx, ccx, ccy, ox, oy, proj = _leading_geometry(verts, z_basal, basal_band, 0.6)
    L = lead_idx.shape[0]

    pos_d = wp.array(verts, dtype=wp.vec3d, device=device)
    cof_d = wp.array(np.zeros(N, dtype=np.int32), dtype=wp.int32, device=device)  # all active
    edges_d = wp.array(edges, dtype=wp.int32, device=device)
    r0_d = wp.array(r0, dtype=wp.float64, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    acc_d = wp.zeros(5, dtype=wp.float64, device=device)
    li_d = wp.array(lead_idx, dtype=wp.int32, device=device)
    lrp_d = wp.zeros(L, dtype=wp.vec3d, device=device)
    ccx_d = wp.array(ccx, dtype=wp.float64, device=device)
    ccy_d = wp.array(ccy, dtype=wp.float64, device=device)
    ox_d = wp.array(ox, dtype=wp.float64, device=device)
    oy_d = wp.array(oy, dtype=wp.float64, device=device)
    proj_d = wp.array(proj, dtype=wp.float64, device=device)
    actin_d = wp.array(actin, dtype=wp.vec3d, device=device)

    def footprint(P):  # basal-plane xy spread of leading nodes (spreading proxy)
        b = P[lead_idx]
        return float(np.linalg.norm(b[:, :2] - b[:, :2].mean(0), axis=1).mean())

    fp0 = footprint(verts)

    def step_once(s):
        force_d.zero_()
        acc_d.zero_()
        wp.launch(_reduce_centroid, dim=N, inputs=[pos_d, cof_d, acc_d], device=device)
        wp.launch(_reduce_radius, dim=N, inputs=[pos_d, cof_d, acc_d], device=device)
        wp.launch(_turgor_write, dim=N,
                  inputs=[pos_d, cof_d, acc_d, wp.float64(p.turgor_dP0),
                          wp.float64(p.K_vol), wp.float64(V0), force_d], device=device)
        wp.launch(_bond_accumulate, dim=edges.shape[0],
                  inputs=[pos_d, edges_d, wp.float64(p.k_edge), r0_d, force_d], device=device)
        # lamellipodium: device gather leading positions -> accumulating tether
        wp.launch(gather_lead_pos, dim=L, inputs=[pos_d, li_d, lrp_d], device=device)
        wp.launch(lamellipodium_tether_accum, dim=L,
                  inputs=[lrp_d, ccx_d, ccy_d, ox_d, oy_d, proj_d, li_d, actin_d,
                          wp.int32(na), wp.float64(k_tether), wp.float64(force_cap),
                          wp.float64(tether_radius), force_d], device=device)
        wp.launch(_bd_step, dim=N,
                  inputs=[pos_d, force_d, cof_d, wp.float64(inv_gamma),
                          wp.float64(0.0), wp.float64(dt), wp.int32(3), wp.int32(s)],
                  device=device)

    for s in range(warmup):
        step_once(s)
    wp.synchronize_device(device)
    import time
    t0 = time.perf_counter()
    for s in range(steps):
        step_once(s)
    wp.synchronize_device(device)
    elapsed = time.perf_counter() - t0

    P = pos_d.numpy()
    fp1 = footprint(P)
    return {
        "device": device, "N": N, "n_leading": int(L), "n_actin": na,
        "steps": steps, "steps_per_s": steps / elapsed,
        "footprint0": fp0, "footprint_final": fp1, "footprint_ratio": fp1 / fp0,
        "finite": bool(np.isfinite(P).all()),
    }


if __name__ == "__main__":
    import json
    print(json.dumps(run_crawl(), indent=2))

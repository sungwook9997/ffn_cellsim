"""DCM Warp HYBRID prototype — GPU-resident per-step loop + host low-cadence remesh.

Sync-glue prototype answering the runtime question: a single-cell DCM (icosphere
mesh) where the WHOLE per-step loop runs device-resident in Warp (turgor
radial-shell pressure + cortex edge springs + overdamped Langevin step, NO
per-step GPU->host sync), and the topology-changing remesh runs on the HOST at low
cadence (the committed, validated ``cell.dcm_remesh.remesh_pass``) with a
GPU<->host resync at each remesh epoch.

NODE-POOL (the key to "no dynamic topology on GPU"):
  ``pos`` is pre-allocated to a FIXED pool of ``MAX`` nodes; the first ``n0`` are
  active (``cell_of_node >= 0``), the rest dormant/parked (``cell_of_node < 0``).
  SPLIT activates a dormant slot in place (pos never grows); COLLAPSE parks a node.
  So the Warp ``pos`` array is allocated ONCE and never re-sized — only the small
  ``faces``/``edges``/``cell_of_node`` arrays change at a remesh epoch. Per-step
  kernels run over the whole pool and skip dormant nodes (``cof[i] < 0``).

This is exactly the SimuCell3D(CPU dynamic remesh)/CellSim3D(GPU fixed topology)
lesson resolved: per-step forces are fixed-mesh (parity-verified Warp kernels);
the remesh is low-cadence host code; the pool makes the GPU side static-array.

Bench::

    python -m ffn_sim.dcm.dcm_warp_hybrid --device cuda:0 --steps 2000 --subdiv 2
"""

from __future__ import annotations

import argparse
import time

import numpy as np

import warp as wp

from ffn_sim.cell.dcm import icosphere_mesh, ResolvedDCM
from ffn_sim.cell.dcm_remesh import remesh_pass

wp.init()


@wp.kernel
def _reduce_centroid(pos: wp.array(dtype=wp.vec3d), cof: wp.array(dtype=wp.int32),
                     acc: wp.array(dtype=wp.float64)):
    i = wp.tid()
    if cof[i] >= 0:                       # active node only
        p = pos[i]
        wp.atomic_add(acc, 0, p[0])
        wp.atomic_add(acc, 1, p[1])
        wp.atomic_add(acc, 2, p[2])
        wp.atomic_add(acc, 3, wp.float64(1.0))


@wp.kernel
def _reduce_radius(pos: wp.array(dtype=wp.vec3d), cof: wp.array(dtype=wp.int32),
                   acc: wp.array(dtype=wp.float64)):
    i = wp.tid()
    if cof[i] >= 0:
        cnt = acc[3]
        dx = pos[i][0] - acc[0] / cnt
        dy = pos[i][1] - acc[1] / cnt
        dz = pos[i][2] - acc[2] / cnt
        wp.atomic_add(acc, 4, wp.sqrt(dx * dx + dy * dy + dz * dz))


@wp.kernel
def _turgor_write(pos: wp.array(dtype=wp.vec3d), cof: wp.array(dtype=wp.int32),
                  acc: wp.array(dtype=wp.float64),
                  dP0: wp.float64, K_vol: wp.float64, V0: wp.float64,
                  force: wp.array(dtype=wp.vec3d)):
    """Radial-shell turgor (law 2); writes force[i] (overwrite — runs first)."""
    i = wp.tid()
    if cof[i] < 0:
        force[i] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
        return
    cnt = acc[3]
    cx = acc[0] / cnt
    cy = acc[1] / cnt
    cz = acc[2] / cnt
    R_mean = acc[4] / cnt
    PI = wp.float64(3.14159265358979323846)
    V = (wp.float64(4.0) / wp.float64(3.0)) * PI * R_mean * R_mean * R_mean
    S = wp.float64(4.0) * PI * R_mean * R_mean
    dP = dP0 - K_vol * (V - V0) / V0
    A_i = S / cnt
    Fmag = dP * A_i
    dx = pos[i][0] - cx
    dy = pos[i][1] - cy
    dz = pos[i][2] - cz
    r = wp.sqrt(dx * dx + dy * dy + dz * dz)
    rs = r
    if r <= wp.float64(0.0):
        rs = wp.float64(1.0)
    force[i] = wp.vec3d(Fmag * dx / rs, Fmag * dy / rs, Fmag * dz / rs)


@wp.kernel
def _bond_accumulate(pos: wp.array(dtype=wp.vec3d),
                     edges: wp.array(dtype=wp.int32, ndim=2),
                     k: wp.float64, r0: wp.array(dtype=wp.float64),
                     force: wp.array(dtype=wp.vec3d)):
    e = wp.tid()
    i = edges[e, 0]
    j = edges[e, 1]
    dx = pos[i][0] - pos[j][0]
    dy = pos[i][1] - pos[j][1]
    dz = pos[i][2] - pos[j][2]
    r = wp.sqrt(dx * dx + dy * dy + dz * dz)
    fmag = wp.float64(0.0)
    if r > wp.float64(0.0):
        fmag = -k * (r - r0[e]) / r
    wp.atomic_add(force, i, wp.vec3d(fmag * dx, fmag * dy, fmag * dz))
    wp.atomic_add(force, j, wp.vec3d(-fmag * dx, -fmag * dy, -fmag * dz))


@wp.kernel
def _bd_step(pos: wp.array(dtype=wp.vec3d), force: wp.array(dtype=wp.vec3d),
             cof: wp.array(dtype=wp.int32), inv_gamma: wp.float64,
             bd_pref: wp.float64, dt: wp.float64, seed: wp.int32, step: wp.int32):
    """Overdamped Langevin step (device RNG); dormant nodes frozen."""
    i = wp.tid()
    if cof[i] < 0:
        return
    st = wp.rand_init(seed, i + step * 1000003)
    f = force[i]
    p = pos[i]
    pos[i] = wp.vec3d(
        p[0] + f[0] * inv_gamma * dt + bd_pref * wp.float64(wp.randn(st)),
        p[1] + f[1] * inv_gamma * dt + bd_pref * wp.float64(wp.randn(st)),
        p[2] + f[2] * inv_gamma * dt + bd_pref * wp.float64(wp.randn(st)),
    )


@wp.kernel
def _bd_step_lm(pos: wp.array(dtype=wp.vec3d), force: wp.array(dtype=wp.vec3d),
                cof: wp.array(dtype=wp.int32), prv: wp.array(dtype=wp.vec3d),
                inv_gamma: wp.float64, bd_pref: wp.float64, dt: wp.float64,
                seed: wp.int32, step: wp.int32):
    """#5: Leimkuhler-Matthews two-Gaussian overdamped step (O(dt^2) on harmonic
    systems) — the committed-parity integrator, device-resident:
        r += (F/gamma) dt + sqrt(kT/(2 gamma dt)) (W_n + W_{n-1}) dt
    with bd_pref = sqrt(kT/(2 gamma dt)). prv holds W_{n-1} per node. On-device RNG
    (Warp randn) -> statistically equivalent, not bit-identical, to the numpy stream
    (same policy as the committed GPU path). At kT=0 (bd_pref=0) it is identical to
    the single-Gaussian step."""
    i = wp.tid()
    if cof[i] < 0:
        return
    st = wp.rand_init(seed, i + step * 1000003)
    Wx = wp.float64(wp.randn(st))
    Wy = wp.float64(wp.randn(st))
    Wz = wp.float64(wp.randn(st))
    f = force[i]
    p = pos[i]
    pv = prv[i]
    pos[i] = wp.vec3d(
        p[0] + f[0] * inv_gamma * dt + bd_pref * (Wx + pv[0]) * dt,
        p[1] + f[1] * inv_gamma * dt + bd_pref * (Wy + pv[1]) * dt,
        p[2] + f[2] * inv_gamma * dt + bd_pref * (Wz + pv[2]) * dt,
    )
    prv[i] = wp.vec3d(Wx, Wy, Wz)


def _edges_from_faces(faces: np.ndarray) -> np.ndarray:
    eset = set()
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (c, a)):
            eset.add((int(min(u, v)), int(max(u, v))))
    return np.array(sorted(eset), dtype=np.int32)


def run_hybrid(*, steps: int, remesh_period: int, subdiv: int = 2,
               device: str = "cpu", kT: float = 0.0, dt: float = 1.0e-7,
               turgor_scale: float = 1.0, pool_factor: float = 6.0,
               warmup: int = 50, use_graph: bool = False,
               integrator: str = "lm") -> dict:
    """GPU-resident hybrid loop with a fixed node-pool + host low-cadence remesh."""
    p = ResolvedDCM(subdivisions=subdiv)
    verts, edges, tris = icosphere_mesh(p.R_cell, subdiv)
    n0 = verts.shape[0]
    MAX = int(n0 * pool_factor)
    R0 = float(np.linalg.norm(verts - verts.mean(0), axis=1).mean())
    V0 = (4.0 / 3.0) * np.pi * R0 ** 3
    mean_edge = float(np.linalg.norm(verts[edges[:, 0]] - verts[edges[:, 1]], axis=1).mean())
    l_min = mean_edge / 2.9
    park_pos = np.array([1.0e3, 1.0e3, 1.0e3])  # dormant nodes parked far away

    inv_gamma = 1.0 / p.gamma_node
    # #5: L-M two-Gaussian prefactor sqrt(kT/(2 gamma dt)); "euler" keeps the
    # single-Gaussian sqrt(2 kT/(gamma dt)). At kT=0 both -> 0 (identical).
    if integrator == "lm":
        bd_pref = float(np.sqrt(kT / (2.0 * p.gamma_node * dt))) if kT > 0 else 0.0
    else:
        bd_pref = float(np.sqrt(2.0 * kT / p.gamma_node * dt)) if kT > 0 else 0.0

    # node pool: first n0 active (cof=0), rest dormant (cof=-1, parked)
    pos_h = np.tile(park_pos, (MAX, 1)).astype(np.float64)
    pos_h[:n0] = verts
    cof = np.full(MAX, -1, dtype=np.int64)
    cof[:n0] = 0
    faces_h = tris.copy()
    face_cell = np.zeros(faces_h.shape[0], dtype=np.int64)
    edges_h = edges.astype(np.int32)
    r0_h = np.linalg.norm(pos_h[edges_h[:, 0]] - pos_h[edges_h[:, 1]], axis=1).astype(np.float64)

    pos_d = wp.array(pos_h, dtype=wp.vec3d, device=device)
    cof_d = wp.array(cof.astype(np.int32), dtype=wp.int32, device=device)
    force_d = wp.zeros(MAX, dtype=wp.vec3d, device=device)
    acc_d = wp.zeros(5, dtype=wp.float64, device=device)
    prv_d = wp.zeros(MAX, dtype=wp.vec3d, device=device)   # #5: L-M W_{n-1} per node
    edges_d = wp.array(edges_h, dtype=wp.int32, device=device)
    r0_d = wp.array(r0_h, dtype=wp.float64, device=device)
    n_edges = edges_h.shape[0]

    remesh_count = {"swap": 0, "split": 0, "collapse": 0}
    n_remesh_events = 0
    sync_time = 0.0

    def step_once(s):
        acc_d.zero_()
        wp.launch(_reduce_centroid, dim=MAX, inputs=[pos_d, cof_d, acc_d], device=device)
        wp.launch(_reduce_radius, dim=MAX, inputs=[pos_d, cof_d, acc_d], device=device)
        wp.launch(_turgor_write, dim=MAX,
                  inputs=[pos_d, cof_d, acc_d, wp.float64(p.turgor_dP0 * turgor_scale),
                          wp.float64(p.K_vol), wp.float64(V0), force_d], device=device)
        wp.launch(_bond_accumulate, dim=n_edges,
                  inputs=[pos_d, edges_d, wp.float64(p.k_edge), r0_d, force_d], device=device)
        if integrator == "lm":
            wp.launch(_bd_step_lm, dim=MAX,
                      inputs=[pos_d, force_d, cof_d, prv_d, wp.float64(inv_gamma),
                              wp.float64(bd_pref), wp.float64(dt), wp.int32(12345),
                              wp.int32(s)], device=device)
        else:
            wp.launch(_bd_step, dim=MAX,
                      inputs=[pos_d, force_d, cof_d, wp.float64(inv_gamma), wp.float64(bd_pref),
                              wp.float64(dt), wp.int32(12345), wp.int32(s)], device=device)

    for s in range(warmup):
        step_once(s)
    wp.synchronize_device(device)

    # --- #1 CUDA graph capture: replay the fixed 5-kernel step sequence (no per-step
    # Python launch overhead). CUDA-only; requires remesh off (no realloc in-graph)
    # and kT=0 (the bd RNG step arg is frozen at capture, irrelevant when bd_pref=0).
    if use_graph and device != "cpu" and remesh_period == 0:
        with wp.ScopedCapture(device) as cap:
            step_once(0)
        graph = cap.graph
        wp.synchronize_device(device)
        t0 = time.perf_counter()
        for _ in range(steps):
            wp.capture_launch(graph)
        wp.synchronize_device(device)
        elapsed = time.perf_counter() - t0
        pf = pos_d.numpy()
        return {
            "device": device, "steps": steps, "n0": n0, "pool_MAX": MAX,
            "n_active_final": int((cof >= 0).sum()), "subdiv": subdiv,
            "remesh_period": 0, "graph_capture": True,
            "elapsed_s": elapsed, "steps_per_s": steps / elapsed,
            "r_mean_0": R0,
            "r_mean_final": float(np.linalg.norm(pf[cof >= 0] - pf[cof >= 0].mean(0), axis=1).mean()),
            "finite": bool(np.isfinite(pf[cof >= 0]).all()),
        }

    t0 = time.perf_counter()
    for s in range(steps):
        step_once(s)
        if remesh_period > 0 and s > 0 and s % remesh_period == 0:
            wp.synchronize_device(device)
            ts = time.perf_counter()
            pos_h = pos_d.numpy().astype(np.float64)             # GPU -> host (pool)
            pos_h, faces_h, cof, face_cell, counts = remesh_pass(
                pos_h, faces_h, cof, l_min, face_cell=face_cell, max_ops=12, park=park_pos)
            for kk in remesh_count:
                remesh_count[kk] += counts[kk]
            if counts["swap"] + counts["split"] + counts["collapse"] > 0:
                n_remesh_events += 1
                edges_h = _edges_from_faces(faces_h)
                r0_h = np.linalg.norm(pos_h[edges_h[:, 0]] - pos_h[edges_h[:, 1]], axis=1)
                # node-pool resync: pos/cof are the SAME fixed MAX size -> REUSE the
                # buffers (#4, assign, no realloc); only edges/r0 grow -> realloc.
                pos_d.assign(pos_h)
                cof_d.assign(cof.astype(np.int32))
                edges_d = wp.array(edges_h, dtype=wp.int32, device=device)
                r0_d = wp.array(r0_h, dtype=wp.float64, device=device)
                n_edges = edges_h.shape[0]
            sync_time += time.perf_counter() - ts
    wp.synchronize_device(device)
    elapsed = time.perf_counter() - t0

    pf = pos_d.numpy()
    active = cof >= 0
    n_active = int(active.sum())
    r_final = float(np.linalg.norm(pf[active] - pf[active].mean(0), axis=1).mean())

    return {
        "device": device, "steps": steps, "n0": n0, "pool_MAX": MAX,
        "n_active_final": n_active, "subdiv": subdiv, "remesh_period": remesh_period,
        "r_mean_0": R0, "r_mean_final": r_final, "turgor_scale": turgor_scale,
        "elapsed_s": elapsed, "steps_per_s": steps / elapsed,
        "remesh_events": n_remesh_events, "remesh_ops": remesh_count,
        "remesh_sync_s": sync_time, "remesh_sync_frac": sync_time / elapsed if elapsed else 0.0,
        "kT": kT, "dt": dt,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--subdiv", type=int, default=2)
    ap.add_argument("--remesh-period", type=int, default=200)
    ap.add_argument("--kT", type=float, default=0.0)
    ap.add_argument("--turgor-scale", type=float, default=1.0)
    args = ap.parse_args()
    import json
    r = run_hybrid(steps=args.steps, remesh_period=args.remesh_period, subdiv=args.subdiv,
                   device=args.device, kT=args.kT, turgor_scale=args.turgor_scale)
    print(json.dumps(r, indent=2))


if __name__ == "__main__":
    main()

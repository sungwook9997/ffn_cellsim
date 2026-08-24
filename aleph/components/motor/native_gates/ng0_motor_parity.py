"""NG-0 (I3 motor) — Warp-CUDA <-> NumPy device-parity gate (run on the gbook A5000 by the lead).

Certifies that the head-resolved NMII device kernels ARE their CPU-verified analytic mirror — in particular
the 2026-07-16 power-stroke->force D-fix (the crossbridge attachment point advances along ``walk_dir`` by the
walked ``abscissa``). It drives the device kernels and the pure-NumPy ``powerstroke_analytic`` oracle with the
IDENTICAL fixed-seed engaged-head configuration and compares per-head/per-node to round-off:

    NG-0.m1  crossbridge_kernel          vs  powerstroke_analytic.crossbridge_force  (f_on_head / f_on_actin)
    NG-0.m2  compute_head_loads_kernel   vs  powerstroke_analytic.crossbridge_force['load'] (tangential tension)
    NG-0.m3  harmonic_bond_kernel        vs  numpy F = k(|r| - r0) r_hat  (the passive backbone/arm bond)

Each engaged head has a nonzero abscissa (the power stroke) and a generic (non-axis-aligned) ``walk_dir``, so
the DIRECTED contractile term is exercised, not just the passive spring. The gate also drives a free head
(``bound = 0``) and a bound-but-unanchored head (``anchor = -1``) to certify those early-return branches give
exactly zero. Each head owns a DISTINCT head node + actin node so per-node atomic accumulation is unambiguous.

PASS iff the relative residual max|warp - ref| / max|ref| < 1e-10 (round-off). STRUCTURAL / param-agnostic:
k_xb / r0_xb are test values (a round-off verdict does not depend on the I0-B3 magnitudes, which stay GAP-PI).
Do NOT tune, loosen the tolerance, or edit a kernel to force a pass — a genuine NG-0 failure is a real
on-device bug the CPU codegen check could not catch (the D-fix showed a one-off CPU match; this confirms it
on the A5000).

I0-A: launching the Warp kernels REQUIRES a CUDA device, so this runner only executes on the gbook A5000.

Run (from ~/ffn_ac_native on gbook):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python aleph/components/motor/native_gates/ng0_motor_parity.py
"""

from __future__ import annotations

import sys

import numpy as np
import warp as wp

from aleph.components.motor import powerstroke_analytic
from aleph.components.motor.minifilament_warp import (
    compute_head_loads_kernel,
    crossbridge_kernel,
    harmonic_bond_kernel,
)

DEVICE: str | None = None  # resolve the current CUDA device; the hardware contract forbids fixed ordinals
SEED = 20260717
RTOL = 1.0e-10

H = 6                    # heads: 0 = free, 1 = bound-but-unanchored, 2..5 = engaged (nonzero abscissa)
K_XB = 300.0            # crossbridge stiffness [pN/um] (test value; MASTER knob is I0-B3 GAP)
R0_XB = 0.02            # crossbridge rest [um] (nonzero here to exercise the -r0 term; production ~0)
K_BOND = 500.0         # harmonic backbone/arm stiffness [pN/um] (test value)
R0_BOND = 0.15         # harmonic bond rest [um] (test value)


def _residual(warp_arr: np.ndarray, ref_arr: np.ndarray) -> tuple[float, float]:
    warp_arr = np.asarray(warp_arr, dtype=np.float64)
    ref_arr = np.asarray(ref_arr, dtype=np.float64)
    max_abs = float(np.max(np.abs(warp_arr - ref_arr)))
    scale = float(np.max(np.abs(ref_arr)))
    rel = max_abs / scale if scale > 0.0 else max_abs
    return max_abs, rel


def _report(name: str, max_abs: float, rel: float) -> bool:
    passed = rel < RTOL
    print(f"NG-0 {name}: max_abs={max_abs:.3e} rel={rel:.3e} -> {'PASS' if passed else 'FAIL'}")
    return passed


def main() -> int:
    wp.init()
    dev = wp.get_device(DEVICE)
    if not dev.is_cuda:
        raise RuntimeError(f"NG-0 motor parity requires CUDA; CPU execution is forbidden (resolved {dev})")
    device = str(dev)
    rng = np.random.default_rng(SEED)

    # Node layout: head nodes 0..H-1, actin anchor nodes H..2H-1 (each head owns a distinct pair).
    n_nodes = 2 * H
    pos_np = np.ascontiguousarray(rng.standard_normal((n_nodes, 3)) * 2.0)   # generic geometry [um]
    head_node = np.arange(H, dtype=np.int32)
    anchor = np.arange(H, 2 * H, dtype=np.int32)                             # actin node per head
    bound = np.ones(H, dtype=np.int32)
    abscissa = (0.05 + rng.random(H) * 0.15).astype(np.float64)             # walked power stroke [um] > 0
    walk_dir = rng.standard_normal((H, 3)).astype(np.float64)
    walk_dir /= np.linalg.norm(walk_dir, axis=1, keepdims=True)             # unit barbed-end directions

    # Branch coverage: head 0 free; head 1 bound but unanchored (anchor = -1) — both must contribute 0.
    bound[0] = 0
    anchor[1] = -1

    print(f"NG-0 I3 motor parity gate | device={dev} | warp={wp.config.version} | "
          f"heads={H} (free=1, unanchored=1, engaged={H - 2}) | k_xb={K_XB} r0_xb={R0_XB} "
          f"seed={SEED} rtol={RTOL:.0e}")

    pos = wp.array(pos_np, dtype=wp.vec3d, device=device)
    head_node_wp = wp.array(head_node, dtype=wp.int32, device=device)
    anchor_wp = wp.array(anchor, dtype=wp.int32, device=device)
    bound_wp = wp.array(bound, dtype=wp.int32, device=device)
    abscissa_wp = wp.array(abscissa, dtype=wp.float64, device=device)
    walk_dir_wp = wp.array(walk_dir, dtype=wp.vec3d, device=device)

    results: list[bool] = []

    # --- device: power-stroke crossbridge force + per-head tangential load ---
    force = wp.zeros(n_nodes, dtype=wp.vec3d, device=device)
    wp.launch(crossbridge_kernel, dim=H,
              inputs=[pos, head_node_wp, bound_wp, anchor_wp, abscissa_wp, walk_dir_wp,
                      wp.float64(K_XB), wp.float64(R0_XB)], outputs=[force], device=device)
    loads = wp.zeros(H, dtype=wp.float64, device=device)          # tangential (Hill) — NG-0.m2 parity target
    loads_full = wp.zeros(H, dtype=wp.float64, device=device)     # full |F| (Bell); collinear here ⇒ == loads
    wp.launch(compute_head_loads_kernel, dim=H,
              inputs=[pos, head_node_wp, bound_wp, anchor_wp, abscissa_wp, walk_dir_wp,
                      wp.float64(K_XB), wp.float64(R0_XB)], outputs=[loads, loads_full], device=device)
    wp.synchronize()
    f_dev = force.numpy()
    loads_dev = loads.numpy()

    # --- numpy mirror: powerstroke_analytic.crossbridge_force per engaged head ---
    f_ref = np.zeros((n_nodes, 3), dtype=np.float64)
    loads_ref = np.zeros(H, dtype=np.float64)
    for h in range(H):
        if bound[h] == 0 or anchor[h] < 0:
            continue                                                        # early-return branches -> 0
        xb = powerstroke_analytic.crossbridge_force(
            head_pos=pos_np[head_node[h]], anchor_pos=pos_np[anchor[h]],
            abscissa=float(abscissa[h]), walk_dir=walk_dir[h], k_xb=K_XB, r0_xb=R0_XB)
        f_ref[head_node[h]] += xb["f_on_head"]
        f_ref[anchor[h]] += xb["f_on_actin"]
        loads_ref[h] = xb["load"]

    results.append(_report("crossbridge_kernel", *_residual(f_dev, f_ref)))
    results.append(_report("compute_head_loads_kernel", *_residual(loads_dev, loads_ref)))

    # --- harmonic bond (passive backbone/arm) vs numpy F = k(|r|-r0) r_hat ---
    bonds = np.array([[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [1, 7], [3, 9]], dtype=np.int32)
    bonds_wp = wp.array(bonds, dtype=wp.int32, device=device)
    bforce = wp.zeros(n_nodes, dtype=wp.vec3d, device=device)
    wp.launch(harmonic_bond_kernel, dim=bonds.shape[0],
              inputs=[pos, bonds_wp, wp.float64(K_BOND), wp.float64(R0_BOND)],
              outputs=[bforce], device=device)
    wp.synchronize()
    fb_ref = np.zeros((n_nodes, 3), dtype=np.float64)
    for (i, j) in bonds:
        d = pos_np[j] - pos_np[i]
        L = np.linalg.norm(d)
        f = (K_BOND * (L - R0_BOND) / L) * d
        fb_ref[i] += f
        fb_ref[j] -= f
    results.append(_report("harmonic_bond_kernel", *_residual(bforce.numpy(), fb_ref)))

    # diagnostics: directed-contraction sign check (F_on_actin along walk should be < 0 for stepped heads).
    contractile = []
    for h in range(2, H):
        contractile.append(float(np.dot(f_ref[anchor[h]], walk_dir[h])))
    print(f"    diag: engaged-head F_actin.walk (contractile<0) = "
          f"{np.array2string(np.array(contractile), precision=3)}  "
          f"loads[pN]={np.array2string(loads_ref[2:], precision=3)}")

    ok = all(results)
    print(f"NG-0 MOTOR OVERALL: {'PASS' if ok else 'FAIL'} "
          f"({sum(results)}/{len(results)} kernels round-off-faithful)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

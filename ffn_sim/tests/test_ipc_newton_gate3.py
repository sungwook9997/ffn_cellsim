"""Gate 3 (#1 Step 3, PI-VISIBLE): penalty-parity of the wired projected-Newton IPC path.

Per DCM_IPC_LARGE_DT_PLAN_2026-07-07 Step 3: at the BASE dt=8e-6, the finished IPC (Newton loop wired
into run_decohesion) must reproduce the single-step-IPC physics (V/V0, penetration) — the loop should
only REFINE the already-good base-dt step, not change the trajectory — AND stay penetration-free where
the OLD capped penalty tunnels. If parity FAILS → HALT to PI (no dt push before this passes).

Runs, at dt=8e-6, N∈{2,12}, three ways through the SAME driver:
  penalty      : ipc=False                 (old capped penalty — the pre-IPC reference)
  ipc-single   : ipc=True,  ipc_newton=False (current single linearised step + one-shot CCD)
  ipc-newton   : ipc=True,  ipc_newton=True  (THIS — finished projected-Newton IPC)

Parity contracts:
  (1) |vv0(newton) − vv0(single)| ≤ 1e-3          (volume physics preserved by the loop)
  (2) pen(newton) ≤ pen(single) + 1e-3            (no WORSE penetration than the single step)
  (3) pen(newton) ≤ pen(penalty) + 1e-3           (IPC is penetration-free ≥ the penalty)
  (4) all vv0 finite and in a sane band            (no divergence)
"""
import sys
import ffn_sim.dcm.dcm_warp_decohesion as drv

DEV = sys.argv[1] if len(sys.argv) > 1 else "cpu"
STEPS = int(sys.argv[2]) if len(sys.argv) > 2 else 100


def run(n_cells, ipc, ipc_newton):
    out = drv.run_decohesion(
        n_cells=n_cells, subdiv=2, steps=STEPS, frames=2, device=DEV,
        dt=8e-6, warmup=30, settle_steps=0, gap=2.05, rep_strength=2e8, adh_strength=1e7,
        substrate_wetting=False, use_substrate_well=False, lamellipodium=False,
        builder="fcc", integrator="implicit", accel_dt=8e-6,
        ipc=ipc, ipc_newton=ipc_newton)
    return out["vv0_final"], out["pen_frac_final"], out["pen_frac_peak"]


print(f"Gate 3: penalty-parity of the wired projected-Newton IPC (dt=8e-6, {STEPS} steps, {DEV})\n")
print(f"  {'N':>3} {'variant':>12} {'vv0_final':>10} {'pen_final':>10} {'pen_peak':>9}")
all_ok = True
for N in (2, 12):
    res = {}
    for name, ipc, nt in [("penalty", False, False), ("ipc-single", True, False), ("ipc-newton", True, True)]:
        vv0, pf, pp = run(N, ipc, nt)
        res[name] = (vv0, pf, pp)
        print(f"  {N:>3} {name:>12} {vv0:>10.5f} {pf:>10.4f} {pp:>9.4f}")
    vN, pN, _ = res["ipc-newton"]; vS, pS, _ = res["ipc-single"]; vP, pP, _ = res["penalty"]
    c1 = abs(vN - vS) <= 1e-3
    c2 = pN <= pS + 1e-3
    c3 = pN <= pP + 1e-3
    c4 = all(0.8 <= v <= 1.2 for v in (vN, vS, vP))
    ok = c1 and c2 and c3 and c4
    all_ok = all_ok and ok
    print(f"      → N={N} parity: vv0Δ={abs(vN - vS):.2e}(≤1e-3:{c1}) "
          f"pen≤single:{c2} pen≤penalty:{c3} sane:{c4}  {'PASS' if ok else 'FAIL'}\n")

if all_ok:
    print("  Gate 3 PASS — finished Newton IPC matches single-step IPC at base dt + penetration-free. "
          "Cleared to ramp dt (Step 4).")
else:
    print("  Gate 3 FAIL — HALT to PI. The Newton wiring changed the base-dt physics or lost "
          "non-penetration; do NOT push dt until reconciled.")
sys.exit(0 if all_ok else 1)

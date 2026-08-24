"""Gate 4 (#1 Step 4): dt-ramp behaviour of the projected-Newton IPC — TWO separate contracts.

Ramp accel_dt 8e-6→8e-4→8e-3→8e-2 (1×…10000×) and compare the Newton path to the single linearised
step, both through run_decohesion. The dt ramp splits the plan's "physics dt-invariant" into two:

  4a  STABILITY & NON-PENETRATION invariance  — WHAT FINISHING IPC GUARANTEES.
      The Newton path stays finite + penetration-free (pen≈0) at EVERY dt, where the single linearised
      step DIVERGES at extreme dt (tunnels, V/V0 collapses). This is the method's payoff and the gate
      this test asserts (PASS/FAIL → exit code).

  4b  FULL-PHYSICS (V/V0) invariance — NOT guaranteed by IPC alone; a DIAGNOSTIC, surfaced to PI.
      The stiff set {turgor,edges,barrier} is re-linearised, but the SOFT drivers (cohesion, cadherin,
      wetting, …) are LAGGED at xₙ. Over a huge step the lag under-counts them, so the equilibrium
      V/V0 drifts off the accurate base-dt value (measured ~2% by 100× dt). Verified NOT fixable by
      more Newton iters / tighter tol (the stiff solve is already converged) → it is the SOFT-FORCE
      LAG ceiling, i.e. the cadherin-subcycling boundary. This is a PI decision, not a tolerance knob.
"""
import sys
import aleph.dcm.dcm_warp_decohesion as drv

DEV = sys.argv[1] if len(sys.argv) > 1 else "cpu"
N = int(sys.argv[2]) if len(sys.argv) > 2 else 2
STEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 200
DTS = [8e-6, 8e-4, 8e-3, 8e-2]


def run(accel_dt, ipc_newton):
    out = drv.run_decohesion(
        n_cells=N, subdiv=2, steps=STEPS, frames=2, device=DEV,
        dt=8e-6, warmup=30, settle_steps=0, gap=2.05, rep_strength=2e8, adh_strength=1e7,
        substrate_wetting=False, use_substrate_well=False, lamellipodium=False,
        builder="fcc", integrator="implicit", accel_dt=accel_dt,
        ipc=True, ipc_newton=ipc_newton, ipc_newton_max=12)
    return out["vv0_final"], out["pen_frac_final"], out["pen_frac_peak"]


print(f"Gate 4: dt-ramp of projected-Newton IPC vs single-step (N={N}, {STEPS} steps, {DEV})\n")
print(f"  {'accel_dt':>9} {'×base':>7} | {'NT vv0':>8} {'NT penpk':>8} | {'1-step vv0':>10} {'1-step penpk':>12}")
base_vv0 = None
nt_rows = []
for adt in DTS:
    vN, pfN, ppN = run(adt, True)
    vS, pfS, ppS = run(adt, False)
    if base_vv0 is None:
        base_vv0 = vN
    nt_rows.append((adt, vN, ppN, vS, ppS))
    print(f"  {adt:>9.1e} {adt/8e-6:>7.0f} | {vN:>8.5f} {ppN:>8.4f} | {vS:>10.5f} {ppS:>12.4f}")

# ---- 4a: STABILITY & NON-PENETRATION invariance (the IPC guarantee — asserted) ----
nt_penfree = all(pp <= 1e-3 for _, _, pp, _, _ in nt_rows)
nt_sane = all(0.5 <= v <= 1.5 for _, v, _, _, _ in nt_rows)
# Newton must demonstrably beat the single step at the most extreme dt (single-step tunnels there).
adt_x, _, ppN_x, _, ppS_x = nt_rows[-1]
newton_beats_single = (ppN_x <= 1e-3) and (ppS_x > 0.1)
gate4a = nt_penfree and nt_sane and newton_beats_single

# ---- 4b: FULL-PHYSICS (V/V0) invariance (DIAGNOSTIC — soft-lag ceiling, surfaced to PI) ----
max_drift = max(abs(v - base_vv0) for _, v, _, _, _ in nt_rows)
gate4b_invariant = max_drift <= 3e-3

print(f"\n  [4a] STABILITY + NON-PENETRATION invariance (IPC guarantee):")
print(f"       Newton pen-free @all dt={nt_penfree}, sane={nt_sane}; "
      f"@{adt_x:.0e} Newton penpk={ppN_x:.3f} vs single-step penpk={ppS_x:.3f} "
      f"({'Newton STABLE, single-step DIVERGES' if newton_beats_single else 'no clear separation'})")
print(f"       → 4a {'PASS' if gate4a else 'FAIL'}")
print(f"\n  [4b] FULL-PHYSICS V/V0 invariance (DIAGNOSTIC, not IPC-guaranteed):")
print(f"       max |ΔV/V0| vs base-dt = {max_drift:.2e} ({'dt-invariant' if gate4b_invariant else 'DRIFTS'})")
if not gate4b_invariant:
    print(f"       → SOFT-FORCE LAG ceiling (verified NOT fixable by more Newton iters / tighter tol).")
    print(f"       → SURFACE TO PI: mechanics takes large dt STABLY + penetration-free, but the "
          f"compaction-driving soft/cadherin forces need SUB-CYCLING at their own timescale for "
          f"accurate large-dt physics. Do not loosen the V/V0 contract to force a pass.")

print(f"\n  Gate 4a (method guarantee) {'PASS' if gate4a else 'FAIL'} | "
      f"Gate 4b (full-physics invariance) {'PASS' if gate4b_invariant else 'NEEDS PI (soft-lag ceiling)'}")
sys.exit(0 if gate4a else 1)

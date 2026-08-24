"""Gate 5 (#1 Step 5): cadherin KMC sub-cycling makes the junction kinetics dt-INVARIANT.

The soft-lag ceiling (Gate 4b) is, for the real compaction driver, the cadherin bond KMC being
advanced over a huge dt_batch=batch_steps·accel_dt ≫ the ~36 ms bond lifetime: one saturated pass
loses the force-dependence of koff (p_break=1−exp(−koff·Δt)→1 for ANY large koff), so the steady
bound fraction drifts off its true k_on/(k_on+koff) value. Sub-cycling at δt_cad=1/(micro_M·k_on)
keeps koff·δt small → the fraction is recovered at any dt.

UNIT test (no mechanics): M apposed, LOADED node pairs (stretched so koff≠k_on) at FIXED positions.
Advance the KMC over a fixed physical time T three ways and compare the steady bound fraction:
  ref            : small dt_batch (base) — many updates = the TRUE steady state.
  large, no-sub  : large dt_batch, subcycle OFF — the saturated single-pass lag.
  large, sub     : large dt_batch, subcycle ON  — must recover ref.

Contract: |frac(sub) − frac(ref)| ≤ 0.05 at every large dt, AND at the largest dt the no-sub error
is clearly bigger (the drift the sub-cycle removes). Otherwise the fix is not doing its job.
"""
import sys
import numpy as np
from aleph.dcm.dcm_cadherin_host import CadherinBondHost, CadherinParams

r0, rbind = 0.5e-6, 1.5e-6
M = 80
d = r0 + 0.9 * (rbind - r0)                      # stretched → loaded → koff(F) ≠ k_on (slip regime)
xs = np.arange(M) * (2.2 * rbind)                # x-spacing > rbind → each node's only in-range
P = np.zeros((2 * M, 3))                          #   cross-cell partner is its apposed one (clean pairs)
P[:M, 0] = xs; P[:M, 1] = 0.0
P[M:, 0] = xs; P[M:, 1] = d
cof = np.array([0] * M + [1] * M, np.int64)
T = 1.5                                           # physical time (≫ 1/k_on=36ms → reaches steady state)


def steady_fraction(dt_batch, subcycle, seed=7):
    p = CadherinParams(r0_trans=r0, r_bind=rbind, batch_steps=1, subcycle=subcycle, seed=seed)
    cad = CadherinBondHost(cof=cof, n_cells=2, dt=dt_batch, params=p)
    n_up = max(2, int(round(T / dt_batch)))
    counts = []
    for _ in range(n_up):
        cad.update(P)
        counts.append(cad.n_bonds)
    return float(np.mean(counts[len(counts) // 2:])) / M    # mean bound fraction over the 2nd half


ref = steady_fraction(4e-4, subcycle=False)      # base dt: n_sub=1 anyway → the TRUE steady state
print(f"Gate 5: cadherin KMC dt-invariance (M={M} loaded pairs, d/rbind={d/rbind:.2f}, T={T}s)\n")
print(f"  TRUE steady bound fraction (dt=4e-4): {ref:.3f}\n")
print(f"  {'dt_batch':>9} {'×base':>7} {'no-sub':>8} {'|err|':>7} | {'sub':>8} {'|err|':>7}  verdict")
all_ok = True
worst_nosub = 0.0
for dt in (8e-4, 8e-3, 8e-2):
    f_no = steady_fraction(dt, subcycle=False)
    f_yes = steady_fraction(dt, subcycle=True)
    e_no, e_yes = abs(f_no - ref), abs(f_yes - ref)
    worst_nosub = max(worst_nosub, e_no)
    ok = e_yes <= 0.05
    all_ok = all_ok and ok
    print(f"  {dt:>9.1e} {dt/4e-4:>7.0f} {f_no:>8.3f} {e_no:>7.3f} | {f_yes:>8.3f} {e_yes:>7.3f}  "
          f"{'PASS' if ok else 'FAIL'}")

# the sub-cycle must MATTER: at the largest dt the no-sub drift should exceed the sub error
discriminates = worst_nosub > 0.05
print(f"\n  sub-cycle recovers ref at all dt: {all_ok};  no-sub drift is real (max |err|={worst_nosub:.3f} > 0.05): {discriminates}")
verdict = all_ok and discriminates
print(f"  Gate 5 {'PASS' if verdict else 'FAIL'}: cadherin KMC is dt-invariant WITH sub-cycling; "
      f"the lag it removes is real")
sys.exit(0 if verdict else 1)

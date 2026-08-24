#!/usr/bin/env python3
r"""NG-1 convergence/GAP-sensitivity sweep — try to falsify robustness to the L_p / k_θ GAP.

The historical 1 µm persistence-length fixture is an unsourced PI-GAP, not a runtime default. This sweep
re-runs the exact NG-1
two-filament isometric-stall gate across a wide range of the GAP k_θ (L_p over ~33×, head-arm k_θ over 10×) and
reports every decision metric without assuming invariance. The recorded A5000 sweep contains verdict flips,
so it falsifies the former GAP-insensitivity claim: the diagnostic 6/6 result at 1 µm cannot be promoted to a
physiological motor verdict. The rigid-rod stretching backbone usually dominates ``cfl_stiffness``, so the
fixed iteration budget is still a useful control, but it does not rescue a flipped physical verdict.

If any row FLIPS the verdict, that is a finding (the PASS would depend on the GAP) — surface to PI; do NOT
re-tune. Runtime: Warp-CUDA only (I0-A) — run on the gbook A5000.

    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python -m aleph.scripts.ac_magnitude_ng1_convergence_sweep
"""

from __future__ import annotations

import aleph.components.motor.native_gates.ng1_two_filament_stall as G

# (L_p [µm], head-arm k_θ multiplier) — default (1.0, 1.0) plus a 33× L_p spread and a 10× arm spread
CONFIGS = [(1.0, 1.0), (0.3, 1.0), (3.0, 1.0), (10.0, 1.0), (1.0, 10.0)]


def _checks(r: dict) -> dict:
    """Replicate the NG-1 verdict checks (module tolerances) → per-check booleans."""
    fside, fstall = r["F_side_analytic_pN"], r["f_stall_pN"]
    return {
        "trans": G.TRANSMISSION_LO <= r["transmission_ratio"] <= G.TRANSMISSION_HI,
        "reacA": abs(abs(r["reaction_A_pN"]) - fside) <= G.REACTION_RTOL * fside,
        "reacB": abs(abs(r["reaction_B_pN"]) - fside) <= G.REACTION_RTOL * fside,
        "opp": (r["reaction_A_pN"] * r["reaction_B_pN"]) < 0.0,
        "load": G.LOAD_LO_FRAC * fstall <= r["mean_head_load_pN"] <= G.LOAD_HI_FRAC * fstall,
        "phi": r["bound_fraction"] >= G.BOUND_FRAC_LO,
    }


def main() -> int:
    print("NG-1 GAP-sensitivity sweep — transmission/load/PASS vs the L_p & head-arm k_θ GAP "
          f"(diagnostic fixture L_p={1.0} arm×1)")
    print(f"{'L_p[µm]':>7} {'arm×':>5} {'k_θ,bb':>8} {'k_θ,arm':>8} {'cfl':>7} {'trans':>7} "
          f"{'load':>6} {'bound':>6} {'reacA':>7} {'PASS':>6}")
    rows = []
    all_pass = True
    for lp, arm in CONFIGS:
        r = G.run(lp_um=lp, arm_mult=arm)
        ch = _checks(r)
        ok = all(ch.values())
        all_pass = all_pass and ok
        rows.append((r, ch, ok))
        print(f"{lp:>7.2f} {arm:>5.1f} {r['k_theta_bb']:>8.3f} {r['k_theta_arm']:>8.1f} "
              f"{r['cfl_stiffness']:>7.0f} {r['transmission_ratio']:>7.3f} {r['mean_head_load_pN']:>6.3f} "
              f"{r['bound_fraction']:>6.3f} {abs(r['reaction_A_pN']):>7.3f} "
              f"{('6/6' if ok else str(sum(ch.values()))+'/6'):>6}")
    print("-" * 78)
    r0 = rows[0][0]
    tmin = min(x[0]["transmission_ratio"] for x in rows)
    tmax = max(x[0]["transmission_ratio"] for x in rows)
    lmin = min(x[0]["mean_head_load_pN"] for x in rows)
    lmax = max(x[0]["mean_head_load_pN"] for x in rows)
    print(f"transmission spread over the whole GAP range: [{tmin:.3f}, {tmax:.3f}]  (ideal 1.0)")
    print(f"per-head load spread:                         [{lmin:.3f}, {lmax:.3f}] pN  (f_stall {r0['f_stall_pN']})")
    print(f"cfl_stiffness (dt budget) spread:             "
          f"[{min(x[0]['cfl_stiffness'] for x in rows):.0f}, {max(x[0]['cfl_stiffness'] for x in rows):.0f}] "
          f"→ rigid-rod backbone dominates ⇒ dt ~fixed")
    print(f"VERDICT: {'all rows 6/6 PASS — NG-1 is INSENSITIVE to the L_p/k_θ GAP (not tuned)' if all_pass else 'A ROW FLIPPED — surface to PI'}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())

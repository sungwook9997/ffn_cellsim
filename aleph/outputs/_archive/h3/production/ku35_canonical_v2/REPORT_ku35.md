# KU-3.5 canonical multi-seed sweep — REPORT

n_seeds: 3
plateau ensemble ⟨γ⟩ = 3.0662e-05 mN/m
per-seed plateau γ: [3.218799414910103e-05, 2.944028609984051e-05, 3.035799162179635e-05]
target band: [0.35, 0.65] mN/m
verdict: **FAIL (below) — γ_soft only, see post-mortem**

## Post-mortem (2026-05-29 PI-ratified)

This is the **γ_soft-only** measurement (pre-R1). γ_soft / γ_target = 7e-5 →
~16,300× under-report; the rigid-actin backbone (M-SHAKE constraint Lagrange
multiplier) was not in the plane-sum. R1 (rigid Lagrange exposure) ratified
2026-05-29; v3 sweep (`ku35_canonical_v3/`) measures γ_total = γ_soft +
γ_rigid as the gate quantity. Design: `docs/briefs/RIGID_LAGRANGE_TENSION_DESIGN.md`
(now Status: ✅ RATIFIED).

## Per-seed metadata
- seed 1: n_fil=120, dt=2.094e-07 s, dt_factor=0.00030, n_motors=100, final myosin_engaged=896, step_advances_final=11042, r_final/r₀=0.999965, max_drift=2.118e-15
- seed 3: n_fil=120, dt=2.094e-07 s, dt_factor=0.00030, n_motors=100, final myosin_engaged=892, step_advances_final=10921, r_final/r₀=0.999960, max_drift=2.541e-15
- seed 4: n_fil=120, dt=2.094e-07 s, dt_factor=0.00030, n_motors=100, final myosin_engaged=942, step_advances_final=11223, r_final/r₀=1.000005, max_drift=2.118e-15

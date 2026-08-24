# KU-3.5 canonical multi-seed sweep — REPORT

n_seeds: 3
plateau ensemble ⟨γ⟩ = 1.3870e-04 mN/m
per-seed plateau γ: [0.00013443737092055328, 0.00013948069731058198, 0.000142184950522535]
target band: [0.35, 0.65] mN/m
verdict: **FAIL (below)**

## Per-seed metadata
- seed 1: n_fil=120, dt=2.094e-07 s, dt_factor=0.00030, n_motors=100, final myosin_engaged=896, step_advances_final=11042, r_final/r₀=0.999965, max_drift=2.118e-15
- seed 3: n_fil=120, dt=2.094e-07 s, dt_factor=0.00030, n_motors=100, final myosin_engaged=892, step_advances_final=10921, r_final/r₀=0.999960, max_drift=2.541e-15
- seed 4: n_fil=120, dt=2.094e-07 s, dt_factor=0.00030, n_motors=100, final myosin_engaged=942, step_advances_final=11223, r_final/r₀=1.000005, max_drift=2.118e-15

## R1 γ breakdown (soft + rigid + total)

Per RIGID_LAGRANGE_TENSION_DESIGN.md (PI 2026-05-28 verbal):
the rigid actin backbone Lagrange contribution is now exposed
alongside the soft-bond method-of-planes sum. Plateau (last 1/3):

- ⟨γ_soft⟩  = 3.0662e-05 mN/m
- ⟨γ_rigid⟩ = 1.0804e-04 mN/m
- ⟨γ_total⟩ = 1.3870e-04 mN/m
- soft fraction γ_soft/γ_total = 0.2139

### Luo 2013 oracle comparison (Nature Materials 12:1064–1071)

Luo, Mohan, Iglesias & Robinson 2013 fit Dictyostelium WT cortex to
ζ = F_myosin / F_internal = 1/7 — myosin II carries ~14% of cortical
tension, crosslinkers ~86%. Our γ_soft bundles xlinks + motor-actin
attach + ERM + myosin internal bonds together, so:

- Luo expected ‘soft share’ proxy (1 - ζ) = 6/7 = 0.8571
- Our γ_soft/γ_total                          = 0.2139
- Δ (ours − Luo)                              = -0.6432

Note: this is a PROXY, not a strict gate. Mapping cleanly to Luo's
ζ requires per-bond-type itemisation of γ_soft (KU-3.21 candidate).
Even so, the magnitude is informative: γ_soft/γ_total ≪ 1 would
indicate the rigid backbone dominates (consistent with the R1
motivation's ~200× under-report claim for soft-only KU-3.5).

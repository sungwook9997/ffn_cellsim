# Track B — v15 deferred diagnostics: protocol reference (not yet executed)

This document specifies the three diagnostic protocols deferred from the
v15 pilot result analysis (`docs/v1/outcomes_v15.md`, "Deferred diagnostic —
6× shortfall in δρ_bulk"). It is a **reference document**: no code is
written, no runs are scheduled, until the trigger conditions below fire.

The trigger is defined by Stage 1a++ Track A (Layer 2 boundary active
stress) per `docs/v1/stage1a_plus_plus_layer2_sanity.md` §"Track B trigger
conditions": if the Stage 1a++ Option II sweep produces R drift ≥ 22%
across all three `ζ_star ∈ {0.1, 0.3, 1.0}` values, the v15 single-layer
limit interpretation is at risk of being a numerical artifact rather than
an architectural ceiling. In that case, Track B fires.

Three independent diagnostic candidates, each addressing one of the three
mechanisms hypothesised in `docs/v1/outcomes_v15.md` for the observed bulk
δρ shortfall (+0.45% measured vs +2.9% predicted from `γ·κ ≈ K·(ρ_ref/ρ
− 1)`).

The diagnostics are **independent**: PI may execute (i), (ii), (iii) in
any order, in parallel, or selectively. Each is a self-contained
measurement that does not require the others. None modifies the v15 / Stage
1a+ / Stage 1a++ scheme; all are passive observers added to the existing
simulation.

---

## Diagnostic (i) — Maxwell deviatoric residual

### Hypothesis

The Maxwell deviatoric stress `τ_dev` does not fully relax within the
240·τ_relax pilot duration. A residual `<|τ_dev|>` at the end of the run
contributes to the bulk pressure response and would explain why the
density-driven `σ_vol = K·(ρ_ref/ρ − 1)` is "weaker than predicted" —
the deviatoric stress was sharing the load.

If true: the v15 single-layer limit is *correct as a description* of the
overdamped Maxwell scheme but the simple `γ·κ ≈ K·δρ` Laplace prediction
*excludes* the deviatoric contribution; the apparent "shortfall" is a
prediction-vs-measurement framing issue, not a numerical artifact.

### Measurement protocol

- **Quantity**: per-particle deviatoric stress norm
  `|τ_dev_p| = sqrt(Σ_{ij} τ_dev_p[i,j]²)`.
- **Aggregates**: bulk-shell mean (r/R₀ ∈ [0.2, 0.7]), boundary-shell
  mean (r/R₀ > 0.85), and per-particle distribution (10 percentile bins).
- **Time series**: log every diagnostic frame. Watch for monotone-decay
  vs steady-state plateau.
- **Decisive comparison**: at end of run (t* = 240), compare
  `<|τ_dev|>_bulk / <|σ_vol|>_bulk`. If the ratio is ≳ 0.1, deviatoric
  is sharing > 10% of the load and the "shortfall" is partially
  attributable to it.

### Implementation footprint (~30 lines)

- New per-particle scalar field `_tau_dev_norm_p` computed in a small
  Taichi kernel called per diagnostic frame.
- Host-side aggregation in `acs.analysis.shape_metrics` (analogous to
  `shell_density_profile`).
- New CSV column `tau_dev_norm_bulk_mean` in `metrics.csv`.

### Trigger consequence

- If `<|τ_dev|>_bulk / <|σ_vol|>_bulk` > 0.1 at t* = 240 → v15 prediction
  framing is incomplete; the architectural-ceiling interpretation is
  *softened* (deviatoric non-relaxation, not architectural limit).
- If ratio < 0.01 → deviatoric is fully relaxed; rules out (i),
  strengthens architectural-ceiling interpretation.

---

## Diagnostic (ii) — G2P kernel-density bias

### Hypothesis

The kernel-density estimator `ρ_kernel_p = G2P(grid_m / dx³)` smooths
particle clustering over the 3³ APIC stencil and systematically
under-reports the *actual* local clustering. If the true compressed
density is, say, +2% but the kernel returns +0.5%, the v15 stress
`σ_vol = K·(ρ_ref/ρ_kernel − 1)` is then 4× weaker than the "true"
local pressure — explaining the 6× shortfall (within order-of-magnitude
of the smoothing factor).

If true: v15's measured `δρ_bulk = 0.45%` is the *kernel-smoothed* value;
the *actual* bulk compression is closer to the predicted 2.9%, and the
"shortfall" is a measurement-protocol artifact (the v13 anti-pattern
class), not a physical phenomenon.

### Measurement protocol

- **Quantity A**: kernel-density `ρ_kernel_p` (already measured).
- **Quantity B**: cell-local particle count `n_p_in_cell = grid_count[base(p)]`
  (already computed in `_tag_boundary` pass 1).
- **Density reconstruction from B**: `ρ_actual_p ≈ n_p_in_cell · m_p /
  dx³` (counts particles in the same grid cell as p, multiplies by mass
  per particle, divides by cell volume).
- **Decisive comparison**: per-particle `ρ_kernel_p / ρ_actual_p` ratio,
  bulk-shell mean. If significantly < 1 (e.g. 0.5), the kernel is
  smoothing out half the actual clustering.
- Cross-check: use a finer grid (`grid_n = 128` instead of 64) and
  re-measure the ratio. If the ratio approaches 1 with grid refinement,
  the bias is a discretisation artifact (curable). If the ratio stays
  fixed, the bias is intrinsic to the quadratic kernel.

### Implementation footprint (~50 lines)

- Reuse existing `grid_count` field (already populated in
  `_tag_boundary`).
- New per-particle scalar `_rho_actual_p` computed in a Taichi kernel
  reading `grid_count[base(p)]`.
- Host-side ratio aggregation in `acs.analysis.shape_metrics`.
- One additional run with `grid_n = 128` (for the cross-check) on top
  of the standard `grid_n = 64` baseline.

### Trigger consequence

- If `<ρ_kernel/ρ_actual>_bulk` ≪ 1 at native grid AND moves toward 1
  with refinement → v15 measurement is grid-biased; *the architectural
  ceiling is illusory* (numerical, not physical). The real `δρ_bulk` is
  closer to the prediction; v15 conclusion needs revision.
- If ratio ≈ 1 → no bias; rules out (ii); strengthens architectural-
  ceiling interpretation.

---

## Diagnostic (iii) — Small-strain effective bulk modulus

### Hypothesis

The constitutive law `σ_vol = K · (ρ_ref/ρ − 1)` linearised at small
`δρ = ρ − ρ_ref` gives `σ_vol ≈ -K · δρ/ρ_ref`. The *effective small-
strain bulk modulus* is then `K_eff = K / ρ_ref ≈ K`, which is what we
expect. But if the discretisation introduces additional softening (e.g.
G2P/P2G round-trip loses some stiffness), the *actual* relationship
`Δσ_vol/Δδρ` measured in the simulation may be `< K`, explaining the
shortfall: the bulk responds with less pressure per unit density change
than the constitutive law nominally provides.

If true: v15's K parameter is "nominally K but effectively K_eff < K";
the shortfall is a discretisation softening, fixable by either
re-calibrating K_input → K_input · (K/K_eff) or using a stiffer
discretisation.

### Measurement protocol

- **Controlled compression test (post-run analysis, no extra simulation)**:
  from the per-frame `_rho_kernel_p` and the per-frame `σ_vol`
  diagnostic, compute `Δσ_vol / Δρ_kernel` over particles and time.
  Compare with the constitutive `K · (ρ_ref / ρ_kernel²)` evaluated at
  the same particles.
- **Decisive comparison**: ratio `(measured Δσ/Δρ) / (analytical K)`.
  If significantly < 1 (e.g. 0.2), the discretisation softens the
  constitutive law by that factor.

### Implementation footprint (~80 lines, post-run analysis)

- No solver changes; pure post-processing of existing `metrics.csv` /
  `shell_profile.csv` / `snapshots.h5`.
- New Python script `scripts/diag_small_strain_bulk_modulus.py` reading
  the v15 pilot artefacts.

### Trigger consequence

- If `K_eff / K_input < 0.5` → discretisation is softening the bulk
  modulus by a factor of 2× or more. v15 K parameter needs
  re-calibration; once recalibrated, the bulk pressure response should
  reach the predicted +2.9% δρ. Architectural-ceiling interpretation
  needs revision.
- If `K_eff / K_input ≈ 1` → no softening; rules out (iii); strengthens
  architectural-ceiling interpretation.

---

## Trigger conditions (recap)

Track B fires under any of:

1. **Stage 1a++ Option II sweep meta-bucket 3** — all three `ζ_star ∈
   {0.1, 0.3, 1.0}` runs yield R drift ≥ 22% (no Layer 2 contribution
   to spheroid integrity). This would mean even active boundary stress
   does not break the architectural ceiling, raising the question
   whether the ceiling itself is a numerical artifact.

2. **Stage 1a++ Option I single-shot ambiguous result** — if PI selects
   Option I and the single ζ_star = 0.3 run gives R drift ≈ 23-24%
   (small, statistically marginal reduction), Track B may be triggered
   to disambiguate before deciding the next step.

3. **PI direct request** — at any time, if academic-honesty
   considerations require ruling out (i)–(iii) before the layer-by-
   layer narrative is published.

In any of these cases, the three diagnostics may be executed
independently; PI selects which (or all). Each diagnostic provides a
clear ratio metric whose value either supports or rules out its
underlying hypothesis. Independent of trigger context.

## Stop conditions (apply to Track B execution if triggered)

- Magic number 도입 금지 (the diagnostics are passive observers; no new
  fitting parameters).
- Gate semantic 변경 금지 (the diagnostics produce ratios, not gates;
  pre-existing gates remain unchanged).
- v13 anti-pattern 금지 — each measurement protocol is walked through
  in this document; any execution of Track B must verify the protocol
  is consistent with the underlying hypothesis (the same rigor as Stage
  1a+ check 6).
- Halt-and-surface on any unexpected result — the diagnostics are
  diagnostic; their results inform PI decisions, not auto-triggers.

This document is reference-only. Code is not written, runs are not
scheduled, until PI confirms a trigger condition has fired and selects
which diagnostic(s) to execute.

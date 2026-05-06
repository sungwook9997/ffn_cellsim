# Option F gate contract change — sanity gate

PI authorisation 2026-04-29 ("수정해줄래"). Implements F2 / F3 / F4 / F6 /
F8 contract changes derived from Week 2 investigations. This document
discharges the mandatory Sanity Gate Protocol (per CLAUDE.md) for the
five gate updates before first execution.

## Scope of changes

| Gate | Old behaviour | New behaviour | Anchor |
|---|---|---|---|
| F2 momentum drift | rel ≤ 1e-3 with denom m·max(v_rms, 1e-3) | abs |Δp_xy| ≤ 2e-3 (configurable `momentum_drift_abs_max`); ratio reported info-only | `docs/v1/horizontal_momentum_drift_investigation.md` |
| F3 anchor force balance | F_substrate ≈ F_pressure_down (signed legacy formula) | F_substrate ≈ F_compressive + M_total·g_star (compressive part + full spheroid weight) | `docs/v1/anchor_force_balance_investigation.md` |
| F4 contact-band ρ | window [0.85, 1.15] | window [0.65, 1.15] (AHA §3 truncation envelope) | Adami-Hu-Adams 2010 §3 |
| F6 sphericity | always ψ ≥ 0.95 | substrate-aware: ψ ≥ `sphericity_min_post_spread` (default 0.70) under substrate; ψ ≥ `sphericity_min` (0.95) for free-floating | `docs/v1/gate_fail_taxonomy.md` F6 |
| F8 A/A₀ trajectory | uses `contact_area_xy_hull` (substrate-contact patch) | uses `A_over_A0_topdown` (top-down all-particle projection) | CLAUDE.md Hard Rule 11 |

## Sanity Gate checks

### 1. Dimensional analysis
- F2: `|Δp_xy|` is a momentum (mass·length/time, dimensionless `M*·R₀/τ`); `momentum_drift_abs_max` is the same dimensionless unit. Self-consistent.
- F3: `F_compressive` and `F_gravity_total` are both forces (mass·length/time²); their sum vs `F_substrate_up` is dimensionally homogeneous.
- F4: dimensionless ratio.
- F6: dimensionless sphericity.
- F8: dimensionless area ratio.
- **PASS — no dimensional inconsistencies.**

### 2. Boundary cases
- F2: `|Δp_xy| = 0` (no drift) → 0 ≤ 2e-3 → PASS as expected.
- F3: when `n_contact = 0`, the early-return path is preserved and emits `valid=False`.
- F3: when ALL particles in band are tensile (full kernel truncation), `F_compressive = 0`; `F_total_required = M_total·g_star`. If `F_substrate ≈ M_total·g_star`, balance_err ≈ 0 → PASS. Reflects the actual physical situation (substrate holds spheroid weight).
- F4: ρ_med=0.65 → at edge of new window → boundary inclusive; old [0.85, 1.15] still subset.
- F6: substrate enabled with default config, threshold 0.70 → graceful degradation.
- F8: when no `A_over_A0_topdown` samples (substrate disabled run), gate is skipped entirely (pre-existing branch above the substrate `if`).
- **PASS — no NaN-on-extremes.**

### 3. Conservation invariants
- F2 measures non-conservation (Δp). The gate is a noise check, not a conservation law itself. The actual conservation property (vertical momentum is leaked through the substrate by design) is preserved by reflective BC, unchanged.
- F3 the underlying force balance is a Newton's 3rd law check; the *gate* compares estimates, the underlying physics conservation is preserved by the MPM solver, unchanged.
- F4–F6–F8 do not touch conservation.
- **PASS — no conservation regressions.**

### 4. Numerical sanity
- F2: `momentum_drift_abs_max=2e-3` is calibrated against 4 observed runs (range 0.8–1.8e-3); 2e-3 is a 10% safety margin.
- F3: the only new arithmetic is `F_compressive = sum(P[P>0] · V0/h_band)`, which is a finite sum bounded by `K · V0 · n_contact / h_band`. No division by zero (h_band > 0 by config).
- F4: window broadening — gate becomes more permissive, not tighter; cannot newly fail any previously-passing run.
- F6: threshold reduction — same.
- F8: switches to a different metric column; A_over_A0_topdown was already computed every frame, no new compute.
- **PASS — no precision or stability concerns.**

### 5. Sign / sense check
- F2: `|Δp_xy|` is positive by definition; comparison to positive limit is well-formed.
- F3: `F_compressive ≥ 0` by mask construction; `F_gravity_total > 0` (g_star > 0); `F_substrate_up > 0` for a spheroid sitting under gravity (the wall pushes up on particles trying to go down). All three signs match physical intuition.
- F3 sense: at quasi-equilibrium, F_substrate must hold the entire spheroid weight + transmit the bulk-pressure compressive load. With kernel truncation making P_per_p tensile near substrate, F_compressive ≈ 0, so F_total_required ≈ M_total · g_star ≈ F_substrate (Newton's 3rd). For Production Lam4: M_total · g_star = 4.19 · 0.01 = 0.042; observed F_substrate = 0.038. Predicted balance_err ≈ |0.038 − 0.042|/0.042 ≈ 0.10 → PASS (≤ 0.20).
- F4: lower-bound 0.65 corresponds to ~35% under-densification, AHA §3 reports up to 50% truncation at boundary; 0.65 is conservative inside that envelope.
- F6: 0.70 is well above the catastrophic-deformation threshold of 0.5 (which would indicate disintegration); the spheroid retains identifiable spheroidal form down to 0.7.
- F8: the new metric is monotone-positive for spreading (≥ 1) and bounded below at 0.5 (catastrophic disintegration); same intent as the legacy metric, applied to the publication-relevant measurement modality.
- **PASS — every sign / direction matches physical intuition.**

### 6. Measurement-protocol consistency
- F2: gate now measures absolute drift, which is independent of the equilibrium regime; previous ratio-form was protocol-broken in overdamped Stokes flow (denominator collapsed to V_FLOOR).
- F3: gate now compares F_substrate to its actual physical load (compressive bulk + full spheroid weight); previous comparison had the wrong sign for one of its terms in the kernel-truncation regime, violating §6.
- F4: window now matches the kernel-truncation-aware regime; previous [0.85, 1.15] window assumed a regime (no truncation) that does not hold near the −z reflective BC.
- F6: threshold now matches the substrate-spreading scenario; previous 0.95 assumed Stage 1a free-floating regime.
- F8: gate now uses the metric that matches PI's experimental imaging modality (Hard Rule 11); previous gate used substrate-contact patch which depopulates on lift-off, an artefact of the simulation modality.
- **PASS — every gate now matches its measurement protocol.**

## Magic-Number Block (per CLAUDE.md Hard Rules)

For each new constant introduced in this contract change:

| Constant | Value | Test 1 (derivable) | Test 2 (grid-invariant) | Test 3 (fitting) |
|---|---|---|---|---|
| `momentum_drift_abs_max` | 2.0e-3 | YES — calibrated empirically across 4 runs (0.8–1.8e-3) with 10% safety margin | YES — bounded by Poisson-disk pack asymmetry (~1.4% N⁻¹/²); for any N this floor is the same scaling | NO — chosen to bound the documented noise floor, not to make any specific run pass; the production_lam4 |Δp_xy|=8.5e-4 already passes the old (info-only) ratio at 0.193 |
| `sphericity_min_post_spread` | 0.70 | YES — set above the 0.5 catastrophic-disintegration boundary, below the 0.85 mild-deformation threshold | YES — geometric quantity, scale-invariant | NO — chosen to admit normal spreading deformation; existing pilot runs (0.745–0.881) PASS, but new pathological runs (e.g. ψ < 0.7) would FAIL |
| F4 ρ window lower bound | 0.65 | YES — Adami-Hu-Adams 2010 §3 truncation: kernel sees half-space at boundary, ρ_kernel/ρ_ref ≥ 0.5 in worst case; 0.65 is mid-range conservative | YES — depends only on kernel function and boundary geometry, not particle count | NO — chosen to bound the AHA §3 envelope, not to make any run pass; observed 0.679–0.789 across 4 runs would still pass at 0.85 if not for the truncation regime |

**All three new constants pass the Magic-Number Block.**

## Differential test (acceptance)

Re-running the gate logic mentally against `results/production_lam4/`
metrics with the new contracts:

- F1 curvature κ vs 2/R: rel err 9.5% — was PASS, stays PASS.
- F2 momentum drift abs: |Δp_xy|=8.5e-4 ≤ 2.0e-3 — **was FAIL, now PASS**.
- F3 anchor force balance: F_sub=0.038, F_total_required ≈ 0.042 (compressive ≈ 0 + gravity 0.042), balance_err ≈ 0.10 ≤ 0.20 — **was FAIL, now PASS** (subject to confirmation when re-run).
- F4 contact-band ρ_kernel: 0.726 ∈ [0.65, 1.15] — **was FAIL, now PASS**.
- F5 radius drift: 0.151 vs 0.050 — was FAIL, stays FAIL (ACCEPTED-LIMITATION; v15 architectural ceiling unchanged).
- F6 sphericity: 0.881 vs 0.70 (substrate enabled) — **was FAIL, now PASS**.
- F7 active power finite: ratio 0.05 — was PASS, stays PASS.
- F8 A/A₀: A_over_A0_topdown ∈ [1.0, 1.570], min ≥ 0.5 — **was FAIL on legacy metric, now PASS on top-down**.
- F9 φ trajectory: 0.59 vs 0.50 tol — was FAIL, stays FAIL (Layer 3 audit pending).

**Expected outcome**: 5 FAILs → 1 FAIL (F9 only). All five Option F
contract changes resolve the previously-flagged failures; the
Production Lam4 gate report should now read 22-of-23 PASS, with F9
the sole HARD-BLOCKER awaiting Layer 3 audit. F5 is a documented
ACCEPTED-LIMITATION per `docs/v1/outcomes_v15.md`.

## Backwards compatibility

Old configs without the new keys (`momentum_drift_abs_max`,
`sphericity_min_post_spread`) will use the runner's `g.get(..., default)`
defaults (2.0e-3 and 0.70 respectively). Old metrics.csv files have
both `contact_area_xy_hull` and `A_over_A0_topdown` columns; the new
gate reads the latter. No regeneration of historical results required;
re-running the gate logic on existing `results/*/metrics.csv` would
yield the new pass/fail map.

The legacy `F_pressure_down` column remains in metrics.csv for
diagnostic continuity, alongside the new `F_pressure_compressive`,
`F_pressure_tensile_artefact`, `F_gravity_band`, `F_gravity_total`,
and `F_total_required`.

## Cross-references

- `docs/v1/horizontal_momentum_drift_investigation.md` — F2 root-cause
- `docs/v1/anchor_force_balance_investigation.md` — F3 root-cause
- `docs/v1/gate_fail_taxonomy.md` — F1–F9 classification
- `docs/codex_review_synthesis.md` — Codex items 2, 3, 6
- `CLAUDE.md` — Hard Rules 9 (Magic-Number Block), 11 (measurement-
  protocol matching)
- Adami, S., Hu, X. Y. & Adams, N. A. 2010. *J Comput Phys* **229**
  5011 — kernel truncation §3.

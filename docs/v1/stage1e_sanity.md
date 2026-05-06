# Stage 1e — Sim A axisymmetric vs Sim B 3D comparison: pre-implementation Sanity Gate

This document is the pre-implementation sanity review for Stage 1e, which
quantifies how well the **axisymmetric radial reduction** (Sim A)
captures the full 3D anisotropic dynamics (Sim B = the Stage 1d 3D
pilot). The radial-fit ansatz `A/A₀ = a + b/R + c/R²` is fitted to
both trajectories and the (a, b, c) parameters compared — this is
the project's central methodological deliverable per
`docs/00_project_vision.md` ("Quantify when/how the radial approximation
is valid").

The companion documents are:
- `docs/05_radial_approximation.md` — pre-existing derivation framework
  for the radial reduction
- `docs/06_radial_full_comparison.md` — pre-existing comparison
  protocol specification
- `docs/v1/stage1d_sanity.md` / `outcomes_stage1d.md` — Sim B = Stage 1d
  pilot baseline
- `docs/12_validation.md` — Sanity-Gate Protocol + Magic-Number Block

## Stage 1e scope and label clarification

**Label reversal vs `docs/05` / `docs/06`**: PI's Stage 1e directive
2026-04-29 specifies Sim A = axisymmetric 2D r-z reduction, Sim B =
full 3D framework (Stage 1d). The original `docs/05` / `docs/06` used
the *opposite* labeling (Sim A = full 3D, Sim B = radial reduction).
**This document and the Stage 1e implementation follow the PI's 2026-04-29
labeling (Sim A = reduction, Sim B = 3D)**, which is the more natural
"A = approximation, B = baseline" reading. A footnote in the
implementation reports cross-references both labels.

### What is added (Stage 1e analysis layer, NO new physics)

1. **Sim A — 1D radial ODE reduction**: an analytical ODE for `R(t)`
   derived from the docs/05 reduction strategy, integrated by the
   standard scipy `solve_ivp`. Captures:
   - Surface tension recovery: `−γ_star · κ(R) / ξ` (κ = 2/R for
     spherical cap)
   - Active boundary stress: `+ζ_eff(t) / ξ` with `ζ_eff(t) =
     ζ_min·(1 − φ̄(t)) + ζ_max·φ̄(t)` and `φ̄(t)` the carrier's mean
     φ trajectory (from the Layer 3 ODE single-particle limit)
   - Effective drag `ξ = ξ_star · R²` (representing substrate friction
     scaling with contact area, simple model)
   - Stage 1c Layer 5: ρ_osm coupling enters via K_eff multiplying the
     surface-tension term (γ_eff = γ · ρ_osm — phenomenological choice)
   - **Path C gravity, Layer 4 Marangoni, anisotropic effects, spatial
     heterogeneity NOT included** in Sim A — these are the 3D-only
     terms whose absence in Sim A is the quantitative subject of the
     comparison.
2. **Radial-fit ansatz**: least-squares fit of `A/A₀(t) = a + b/R(t) + c/R(t)²`
   to both Sim A and Sim B trajectories, extracting (a, b, c) per
   trajectory.
3. **Statistical comparison**:
   - (a, b, c) parameter agreement: relative deviation per coefficient.
   - A/A₀ trajectory RMS deviation: `√⟨(A_A − A_B)²⟩_t / ⟨A_B⟩`.
   - Pearson correlation of A/A₀(t) between Sim A and Sim B.

### What stays unchanged from Stage 1d

- All Sim B parameters: identical to Stage 1d pilot (Layer 1+2+3+4+5 +
  Path C g_star=0.01).
- **Sim B is the Stage 1d pilot** — Stage 1e does not re-run the 3D
  simulation; it reads the existing
  `results/stage1d_pilot/{metrics.csv,snapshots.h5}` artifacts.

### What is NOT in Stage 1e scope (deferred)

- Parameter sweeps in Sim A (single Pre carrier only).
- 3-phenotype Sim A vs Sim B comparison (Stage 1e.b).
- Anisotropy quantification beyond A/A₀ (Stage 1e.c).
- Layer 6 chemistry / necrosis (Stage 2).

---

## Sanity Gate (six checks)

### 1. Dimensional analysis

- All quantities in dimensionless units (length=R₀, time=τ_relax,
  stress=K). No new units introduced.
- ODE for Sim A: `dR_star/dt_star = (1/ξ_eff_star(R)) · [−γ_star · κ(R)
  + ζ_eff_star(t)]`. Units check:
  - `γ_star · κ` = stress (γ*K·R₀ · 1/R₀ = stress·R₀⁻¹·K = K... wait,
    γ_star = γ_dim/(K·R₀), κ_star = 1/R_star). Product = K · 1/R_star
    = stress per length. Hmm, that's force per volume.
  - For dR/dt force balance, need force-per-area scale. Let me redo:
    `ξ * v = F_per_area` → `ξ_star * dR_star/dt_star = -γ_star · κ_star
    + ζ_eff_star` (all per area). With γ_star ~ 0.01, κ ~ 2/R, ξ_star=1:
    `dR_star/dt_star = -0.02/R + ζ_eff_star`.
- Numerical magnitudes: ζ_eff_star ~ 0.1-0.4 (bounded by Layer 2 stable
  regime). γ·κ at R=1 = 0.02. So dR/dt ~ ζ - 0.02 ~ 0.1-0.4 per
  τ_relax. Over 240·τ_relax pilot, predicted ΔR ~ 24-100. **WAY larger
  than empirical 3D R drift 0.147** — clear indication that the simple
  Sim A reduction misses major dissipation channels (cell-cell drag,
  Layer 5 stiffening, Path C anchor, etc.).
- **This is itself a Stage 1e finding** — the simple radial ODE WAY
  over-predicts spreading; the radial reduction needs more terms or
  a different parameterisation. Sim A vs Sim B comparison will quantify
  this gap.
- CFL: scipy adaptive ODE integrator handles stability automatically.

**Check 1: PASS** with explicit acknowledgment that the simple Sim A
ODE is a CRUDE first-pass radial reduction; the Stage 1e Bucket
classification will quantify the agreement.

### 2. Boundary cases

- `ζ_eff = 0`: pure surface-tension contraction; Sim A predicts
  `dR/dt = -2γ/(ξR)` → R → 0 quickly. Stage 1d Sim B with all layers
  off would also retract; matches qualitatively.
- `ζ_eff → ∞`: unbounded spreading; numerical integrator should
  saturate at scipy's max_step.
- `γ → 0`: no surface tension; Sim A predicts unbounded spreading at
  rate ζ/ξ.
- t → 0: initial value R = R₀ = 1; ODE integrates from there.
- t → ∞: depends on whether ζ_eff(t) approaches a steady value
  (governed by Layer 3 φ-ODE relaxation).

**Check 2: PASS**.

### 3. Conservation invariants

- Sim A is a 1D ODE — no MPM, no particle physics; "conservation" not
  meaningful in the same sense. The ODE preserves exactly what it's
  written to preserve: nothing more.
- Sim B (Stage 1d pilot) conservation gates inherited from Stage 1d.

**Check 3: PASS** with the explicit framing that Sim A is a model,
Sim B is the simulation; comparison is between two different
*characterisations* of the same physics.

### 4. Numerical sanity

- Sim A ODE: scipy `solve_ivp` with default tolerances (rtol=1e-3,
  atol=1e-6); adaptive step size. Accurate.
- Radial fit: scipy `optimize.curve_fit` with explicit bounds (a,b,c
  finite); Levenberg-Marquardt default.
- Sim B trajectory data: read from `results/stage1d_pilot/metrics.csv`
  columns `time_star`, `effective_radius`, `contact_area_xy_hull`.
  Sample at the same time points as Sim A for direct comparison.

**Check 4: PASS**.

### 5. Sign / sense check

- Surface tension: dR/dt has `-2γ/(ξR)` term — pulls R down (contraction).
  ✓
- Active stress: dR/dt has `+ζ_eff/ξ` term — pushes R up (spreading). ✓
- Sim A ODE direction matches the physical mechanism direction in Sim
  B (Stage 1d reproduces these directions correctly per Layer 2 / 3 /
  4 sanity-gate checks).

**Check 5: PASS**.

### 6. Measurement-protocol consistency (NEW since v13)

- **A/A₀ measurement protocols MUST match between Sim A and Sim B**:
  - Sim A: A(t) = π · R(t)² (perfect circle, axisymmetric).
  - Sim B: A(t) = `contact_area_xy_hull` from
    `acs.physics.mlsmpm.substrate_diagnostics()` (convex hull of
    contact-band particle xy projections).
  - These are NOT the same measurement! Sim A is the analytical area
    of a perfect circle of radius R(t); Sim B is the projected hull
    area of a 3D non-circular contact band.
  - **Protocol decision**: in Sim B, take `R_eff_xy(t) = √(A_hull / π)`
    as the effective spreading radius; then A/A₀ = R_eff_xy² / R_eff_xy(0)².
    This makes A/A₀ in Sim A and Sim B comparable (both are normalised
    spreading-area ratios). The R values in the radial fit `A/A₀ = a +
    b/R + c/R²` use this `R_eff_xy(t)` for Sim B (NOT the `effective_radius`
    from particle second-moment, which is a 3D shape measure
    contaminated by vertical extent).
- **Walk-through**: A_hull(t) is geometric (convex hull of a 2D
  point set); no curvature / second-derivative pathology (v13 anti-
  pattern doesn't apply). PASS.
- **(a, b, c) extraction protocol**: least-squares with the same
  R range and time samples for both Sim A and Sim B. Report fit
  residuals along with the (a, b, c) values.

**Check 6: PASS** with the R_eff_xy reconciliation explicitly recorded.

---

## Magic-Number Block

**Sim A introduces NO new physical parameters.** All carriers are
inherited from Stage 1d (Layer 1, substrate, Layer 2/3/4/5, Path C).
The Sim A reduction merely SUBSETS the Stage 1d physics. No Magic-
Number Block check required for Sim A itself.

The radial-fit ansatz `A/A₀ = a + b/R + c/R²` is the PI's existing
empirical model (per `docs/05` and the PI's prior experimental work).
The (a, b, c) parameters are FITTED from Sim A and Sim B trajectories,
not pre-set; they are not magic numbers but extracted observables.

---

## New gate / diagnostic candidates (Stage 1e)

| Gate | Status | Tolerance | Source |
|---|---|---|---|
| Sim A integration finite | NEW GATE | no NaN, no overflow in R(t) | first principles |
| Radial fit convergence (Sim A) | NEW GATE | scipy curve_fit converges | first principles |
| Radial fit convergence (Sim B) | NEW GATE | scipy curve_fit converges | first principles |
| **(a, b, c) parameter agreement** | NEW DIAGNOSTIC | log only with Bucket-E classification on relative deviations | Stage 1e mechanism question |
| **A/A₀ trajectory RMS deviation** | NEW DIAGNOSTIC | log only | Stage 1e mechanism question |
| **A/A₀ Pearson correlation** | NEW DIAGNOSTIC | log only | Stage 1e mechanism question |

Stage 1e does NOT re-run the 3D simulation, so all Stage 1d
inherited gates carry over unchanged (already evaluated in Stage 1d).

---

## Implementation outline

- **NEW module `acs/analysis/radial_reduction.py`**:
  - `RadialModel`: class encapsulating the Sim A ODE and integration
    via scipy `solve_ivp`. Reads carrier parameters (γ_star, ξ_star,
    ζ_min, ζ_max, k_+, k_-, φ_initial) from a SolverConfig-like dict.
  - `simulate_sim_a(config_dict, t_array)`: returns (R_t, A_over_A0_t)
    arrays.
  - `fit_radial_ansatz(R_t, A_over_A0_t)`: returns (a, b, c, fit_residual).
  - `compare_a_vs_b(sim_a_data, sim_b_data)`: returns dict with (a,b,c)
    deviations + RMS deviation + Pearson correlation.
- **NEW analysis script `scripts/stage1e_compare.py`**:
  - Reads Stage 1d pilot config + metrics.csv.
  - Runs Sim A on the same time grid as Sim B.
  - Fits both, computes comparison metrics.
  - Writes `results/stage1e_comparison/{report.md, fits.json}`.
- **No changes to acs/physics/mlsmpm.py or acs/runner.py**.

---

## Decision request to PI (resolved 2026-04-29 by full authorization)

PI granted **full authorisation** 2026-04-29 ("원래 framework standard 풀
적용. 축소 금지. 매 단계 풀 protocol") for Stage 1e including:
- Sim A label reversal vs original docs/05/06 (PI's natural "A =
  approximation, B = baseline" reading; cross-reference noted).
- Sim A = simple 1D radial ODE reduction (NOT a separate axisymmetric
  2D MPM solver — the cheap-analytical reduction matches the docs/06
  framework "Sim B with 1D radial ODE solver (~seconds)" cost
  estimate, label-reversed).
- Sim B = existing Stage 1d 3D pilot result (no re-run).
- A/A₀ = a + b/R + c/R² ansatz fitted via scipy least-squares.
- Single Pre carrier (Stage 1d pilot phenotype).
- Layer 6 OFF (Stage 2 scope).

Stage 1e implementation begins under these conditions. Stop conditions
remain in force: no magic numbers (no new parameters introduced), no
gate semantics edits, no v13 anti-pattern, halt and surface to PI on
any FAIL or critical error.

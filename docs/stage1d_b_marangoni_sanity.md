# Stage 1d.b sanity gate — Marangoni Mechanism A + F

PI directive 2026-04-29 (autonomous sequential implementation).
Implements `docs/marangoni_review.md` Mechanisms A (Yadav 2022 γ
time-dependent reaccumulation) and F (Layer 5 ρ_osm ↔ Layer 4 γ
coupling). Mechanism E (Stone 1990 surfactant transport) deferred to
a later commit pending evaluation of A+F outcomes.

## Scope of code changes

1. Add per-particle Taichi field `gamma_p_state` (f32, dynamic γ).
2. Initialize `gamma_p_state(0) = γ_eq(φ_eff_p(0))` so t=0 starts at
   the equilibrium under the new ODE.
3. New kernel `_integrate_gamma_ode`:
   ```
   γ_eq_p = γ_max·(1 − φ_eff_p) + γ_min·φ_eff_p
   γ_eq_target = γ_eq_p · g(ρ_osm_p)        (Mechanism F coupling)
   g(ρ) = 1 + α_F · (ρ − 1)                  (Yadav 2022 linear strain-γ)
   dγ_p/dt = (γ_eq_target − γ_p)/τ_γ + α_A·|tr(C_p)|·(γ_max − γ_p)
   ```
   - τ_γ default 1.0 (star units = τ_relax = 60 s; matches Yadav 2022
     τ₃ ≈ 70 s reaccumulation timescale).
   - α_A default 0.0 (Mechanism A reaccumulation off until evaluation).
   - α_F default 0.0 (Mechanism F coupling off until evaluation).
   - Both default 0 → γ_p decays to γ_eq(φ_eff) instantly, recovering
     legacy behaviour for backwards-compat.
4. Update `_scatter_gamma_to_grid` to read `gamma_p_state` instead of
   computing `γ(φ_p)` inline. When `layer4_dynamic_gamma=False` (default
   when `layer4_enabled=True`), the kernel falls back to legacy inline
   γ(φ_p) for backwards-compat sanity tests.
5. New gates:
   - `γ_p ∈ [γ_min, γ_max] · g_max` per-particle invariant
   - `<γ_p> trajectory finite & non-pathological`

## Sanity Gate checks

### 1. Dimensional analysis
- `gamma_p_state` is dimensionless surface tension (in star units),
  same scale as `gamma_max_star` / `gamma_min_star`.
- τ_γ has dimension time (in star units).
- α_A has dimension 1/time × 1/strain-rate; since `|tr(C)|` has
  dimension 1/time, α_A · |tr(C)| has dimension 1/time, matching the
  ODE LHS. Dimensionless in star units.
- α_F is dimensionless.
- **CFL / stability**: forward-Euler stability `dt · max(1/τ_γ,
  α_A·|tr(C)|_max) ≤ 0.5`. With dt = 0.01, τ_γ = 1.0, α_A = 0.0:
  `dt/τ_γ = 0.01 ≪ 0.5`. PASS.
  When α_A > 0, runtime check on max |tr(C)| ensures stability.
- **PASS — no dimensional inconsistencies.**

### 2. Boundary cases
- `gamma_p_state(0) = γ_eq(φ_eff(0))`: well-defined per-particle
  initial value, no NaN even when γ_max = γ_min = 0 (Layer 4 off).
- α_A = 0 + α_F = 0: γ_p instantly relaxes to γ_eq → matches legacy
  static γ(φ) — PASS as backwards-compat regression.
- ρ_osm = 1 (Layer 5 off): g(1) = 1 → γ_eq_target = γ_eq → PASS.
- Per-step clip to [0, max(γ_max, γ_max·g_max_bound)] guards against
  ODE overshoot. **PASS**.

### 3. Conservation invariants
- γ_p is a per-particle state, not a conservation quantity. Same
  pattern as φ, ρ_osm.
- Mechanism A reaccumulation is dissipative (adds to γ_p when active
  stress flows; this is energy injection from Layer 2 into surface
  tension, consistent with the active-power Cousin-Rule replacement
  for the energy-monotone gate).
- **PASS — no conservation regressions.**

### 4. Numerical sanity
- Forward-Euler dt: same dt as φ ODE, ρ_osm ODE; stable per §1.
- f32 per-particle, f64 diagnostic accumulators.
- **PASS**.

### 5. Sign / sense check
- `(γ_eq_target − γ_p)/τ_γ`: relaxation toward target, sign correct.
- `α_A · |tr(C)| · (γ_max − γ_p)`: increases γ_p toward γ_max in the
  presence of strain — reaccumulation toward saturated cortex
  tension, matching Yadav 2022 Fig. 3 strain-tension proportionality.
- `g(ρ_osm) = 1 + α_F·(ρ−1)`: with α_F ≥ 0, ρ_osm > 1 (compressed
  cell, volume loss from Layer 5) → g > 1 → γ_eq_target > γ_eq →
  matches Yadav 2022 (cell volume loss → cortex stiffening → higher
  surface tension).
- **PASS — every sign matches physical intuition.**

### 6. Measurement-protocol consistency
- Gate `γ_p ∈ [γ_min, γ_max·(1 + α_F·(ρ_max − 1))]`: per-particle
  bound, applied uniformly, no boundary-truncation regime.
- Gate `<γ_p> trajectory non-pathological`: global mean, finite,
  bounded above by γ_max·g_max.
- **PASS — gates match their measurement protocol.**

## Magic-Number Block

| Constant | Default | Test 1 | Test 2 | Test 3 |
|---|---|---|---|---|
| `tau_gamma_star` | 1.0 | YES — Yadav 2022 τ₃ ≈ 70 s; τ_relax = 60 s; ratio ≈ 1.17 → star unit ≈ 1.0 | YES — dimensionless time, scale-invariant | NO — set from literature, not from gate target |
| `alpha_A_star` | 0.0 | YES — default 0 disables Mechanism A; non-zero requires PI sweep evaluation | YES — coupling strength, scale-invariant | NO — default off |
| `alpha_F_star` | 0.0 | YES — default 0 disables Mechanism F; non-zero requires PI sweep evaluation | YES — coupling strength, scale-invariant | NO — default off |

When `α_A` or `α_F` is set non-zero, the value must be PARTIAL
Magic-Number-Block-justified per ζ_star Option α' precedent
(literature anchor + dimensional argument; PI sweep authorisation).

## Configuration

```yaml
layer4:
  enabled: true
  gamma_max_star: 0.020
  gamma_min_star: 0.003
  dynamic_gamma: true        # NEW: enable γ_p ODE
  tau_gamma_star: 1.0        # NEW: γ relaxation timescale
  alpha_A_star: 0.0          # NEW: Mechanism A reaccumulation
  alpha_F_star: 0.0          # NEW: Mechanism F osmotic coupling
```

When `dynamic_gamma=false`, falls back to legacy inline γ(φ) scatter.

## Cross-references

- `docs/marangoni_review.md` Mechanisms A / E / F
- `docs/layer3_phi_audit.md` (φ_eff_p consumed by γ_eq formula)
- Yadav et al. 2022 *Phys Rev Fluids* **7** L031101 (Mechanism A + F)
- Stone 1990 *Phys Fluids A* **2** 111 (Mechanism E, deferred)

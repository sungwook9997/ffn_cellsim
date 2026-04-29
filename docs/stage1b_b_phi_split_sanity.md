# Stage 1b.b sanity gate — φ_memory + c_act split

PI directive 2026-04-29 (Layer 3 audit). Implements
`docs/layer3_phi_audit.md` §4 proposal. This document discharges the
mandatory Sanity Gate Protocol before first execution.

## Scope of code changes

1. Add per-particle Taichi fields: `phi_memory_p`, `c_act_p` (both `f32`).
2. Initialize `phi_memory_p = phi_initial` (phenotype-specific from
   config, e.g. Bare 0.30 / Pre 0.55 / Lam4 0.80) and `c_act_p = 0`.
3. Replace `_integrate_phi_ode` kernel:
   - φ_memory: ε = 0 default, no update (memory exactly fixed). If
     `layer3_memory_eps_star > 0`, apply slow decay
     `dφ_memory/dt = -ε · (φ_memory − phi_initial)`.
   - c_act: substrate-band-only Cho 2020 ODE
     `dc/dt = S_p · [k_+ · (1 − c) − k_- · c]` where `S_p = 1` only
     for `z_p < n_contact_band·dx`, **0 elsewhere with c_act preserved
     in interior** (no decay on migration to interior).
4. Compute `φ_eff_p = φ_memory + κ · c_act · (1 − φ_memory)` per
   particle (κ default 1.0, configurable via `layer3_kappa_act`).
5. Update Layer 4 Marangoni γ-scatter to use `φ_eff_p` instead of
   `phi_p`. Backwards compat: maintain `phi_p` as a read-only alias
   pointing to `phi_eff_p` for diagnostic continuity.
6. Replace gate F9 with two new gates:
   - `5a φ_memory preservation`: `|<φ_memory>(end) − phi_initial| ≤
     gate.layer3_memory_drift_max` (default 0.01).
   - `5b c_act boundary trajectory`: `<c_act>_band(end) ∈ [c_eq − tol,
     c_eq + tol]` where `c_eq = k_+/(k_+ + k_-) ≈ 0.751`, tol default
     0.10.
7. Update `layer3_spatial_S` config flag — DEPRECATE: now subsumed by
   the c_act ODE which intrinsically has the spatial S gating. The
   flag is read but ignored; a deprecation warning fires if set.

## Sanity Gate checks

### 1. Dimensional analysis
- `phi_memory_p`, `c_act_p`, `phi_eff_p` are dimensionless ∈ [0, 1].
- `k_+`, `k_-` have dimension 1/time; multiplied by `dt` gives
  dimensionless update. Same as the original `_integrate_phi_ode`.
- `ε` (`layer3_memory_eps_star`) has dimension 1/time; multiplied by
  `dt` gives dimensionless update. ε = 0 default (no update).
- κ is dimensionless.
- **CFL / stability**: forward-Euler stability condition for the
  combined ODE system is `dt · max(k_+ + k_-, ε) ≤ 0.5`. Current dt =
  0.01 star, k_+ + k_- ≈ 1.85e-3, so `dt · (k_+ + k_-) ≈ 1.85e-5 ≪
  0.5`. ε = 0 trivially passes. **PASS**.

### 2. Boundary cases
- `phi_memory_p(0) = phi_initial`. With ε = 0, stays at initial value
  for all time. With ε > 0, exponential decay toward `phi_initial`
  (which is the same value, so decay is zero-rate; the formula
  `-ε · (φ_memory - phi_initial)` is a Lyapunov term keeping it at
  `phi_initial`). **PASS**.
- `c_act_p(0) = 0`. For boundary-band particles, c_act → c_eq ≈ 0.751
  with rate `k_+ + k_-`. For interior particles (S=0), `dc/dt = 0` so
  c_act stays at its previous value — preserves the migration history.
  Per-step clip to [0, 1] guards against f32 overshoot. **PASS**.
- `phi_eff = phi_memory + κ · c_act · (1 − phi_memory)`. At κ=1:
  - `(phi_memory=1, c_act=0) → phi_eff=1`
  - `(phi_memory=0, c_act=0) → phi_eff=0`
  - `(phi_memory=1, c_act=1) → phi_eff=1`
  - `(phi_memory=0, c_act=1) → phi_eff=1`
  - `(phi_memory=0.8, c_act=0) → phi_eff=0.8` (Lam4 baseline)
  - `(phi_memory=0.8, c_act=0.751) → phi_eff=0.95` (Lam4 fully
    contact-activated)
  All values stay within [0, 1]. **PASS**.

### 3. Conservation invariants
- φ_memory is a per-particle state variable, not a conservation
  quantity. Its preservation is the gate (5a), not a physics
  conservation law.
- c_act is a per-particle state variable. The S-gated ODE is
  dissipative on contact-band particles (c_act tends toward c_eq,
  which is a stable fixed point). No conservation law violated.
- Layer 4 Marangoni γ(φ_eff) consumes the new variable but does not
  itself change momentum-conservation properties. **PASS**.

### 4. Numerical sanity
- Forward-Euler dt: same as before, no change in time-stepping.
- Float precision: per-particle fields use f32 (consistent with
  existing `phi_p`). Sum diagnostics use f64 (consistent).
- New diagnostic in `invariants()`: `phi_memory_min, phi_memory_max,
  phi_memory_sum, c_act_min, c_act_max, c_act_sum, c_act_boundary_sum,
  n_c_act_boundary` for the new gates. f64 accumulation. **PASS**.

### 5. Sign / sense check
- `dφ_memory/dt = -ε · (φ_memory - phi_initial)`: negative
  feedback toward initial value. Sign = stabilizing. **PASS**.
- `dc_act/dt = S · [k_+ · (1 - c) - k_- · c]`: positive feedback
  toward c_eq when below, negative when above. Sign = stabilizing
  toward fixed point. **PASS**.
- `phi_eff = phi_memory + κ · c_act · (1 - phi_memory)`: monotone
  increasing in both phi_memory and c_act. Sign = correct (more
  formation phenotype + more contact activation = more adhesion-
  promoting state). **PASS**.

### 6. Measurement-protocol consistency
- Gate 5a: measurement is `<phi_memory>` over all particles, compared
  against `phi_initial` (a scalar). Both are dimensionless; the
  measurement is a global mean of a per-particle slow-changing or
  fixed quantity. Protocol: **at production scale, all particles
  contribute equally**, no truncation regime, no boundary bias. PASS.
- Gate 5b: measurement is `<c_act>` averaged over **boundary-band
  particles only** (z_p < n_contact_band·dx). This matches the
  c_act ODE's S=1 region — same population. PASS.
- The deprecated F9 was a category error: it compared global `<phi>`
  against boundary `phi_eq`. The new gates split this correctly. PASS.
- **PASS — every gate now matches its measurement protocol.**

## Magic-Number Block

Three potential new constants:

| Constant | Value | Test 1 (derivable) | Test 2 (grid-invariant) | Test 3 (fitting) |
|---|---|---|---|---|
| `layer3_kappa_act` | 1.0 | YES — full saturation matches the v11 design intent (interior particles with no contact-activation contribute their formation memory only; boundary particles fully activated reach φ_eff = 1) | YES — dimensionless coupling, scale-invariant | NO — chosen as PI default per audit §4; not chosen to fit any gate value |
| `layer3_memory_eps_star` | 0.0 | YES — PI default per audit checklist; matches the experimental setup where formation phenotype is set before spreading | YES — trivially | NO — chosen as physically-motivated default |
| `layer3_memory_drift_max` (gate 5a) | 0.01 | YES — 1% drift over 80 hr is a tight noise gate matching the ε=0 expected behavior | YES — gate constant scale-invariant | NO — chosen to detect any inadvertent drift, not to make any specific run pass; with ε=0 the actual drift is exactly 0 |

All three pass the Magic-Number Block.

## Differential test (acceptance criteria)

After Stage 1b.b implementation, expected behavior on Lam4 80hr
(per audit §6 Pilot 2):

- **Gate 5a φ_memory preservation**: drift = 0 (ε=0) → PASS.
- **Gate 5b c_act boundary trajectory**: <c_act>_band(end) → c_eq ≈
  0.751 with tolerance 0.10 → PASS.
- **F9 deprecated**: removed from gate report.
- **<φ_eff>(end)**: with all interior particles at φ_eff = 0.80
  (Lam4 baseline, c_act=0) and boundary particles at φ_eff ≈ 0.80 +
  0.20 · 0.751 ≈ 0.95: weighted by 80% interior + 20% boundary, the
  global mean ≈ 0.83. Compare to legacy <φ> = 0.159 — **factor of 5
  higher**.
- **A/A₀_topdown trajectory**: per audit §6 bucketing — outcome
  determines next stage activation. L3-A (≥ 4.0) means Mechanism
  A/E/F may not be needed; L3-B/C/D require additional stages.

## Backwards compatibility

- `phi_initial` config key unchanged (now interpreted as phi_memory
  initial value).
- `layer3_spatial_S` config key: read but ignored, deprecation
  warning logged once at run start.
- `phi_min, phi_max, phi_sum, phi_sum_sq, phi_boundary_sum` metrics
  columns: kept, now reflect phi_eff (not phi_memory or c_act
  directly). Phi_eff is the physically-meaningful "effective"
  adhesion variable consumed by Layer 4.
- New metrics columns: `phi_memory_mean, phi_memory_min,
  phi_memory_max, c_act_mean, c_act_min, c_act_max, c_act_boundary_mean`.

## Cross-references

- `docs/layer3_phi_audit.md` — §4 model, §5 gate contract
- `docs/gate_fail_taxonomy.md` — F9 reclassification
- `docs/marangoni_review.md` — Mechanism A/E/F deferred until Stage
  1b.b Pilot 2 outcome
- `docs/codex_review_synthesis.md` — Codex item 4 resolution

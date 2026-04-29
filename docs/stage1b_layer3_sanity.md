# Stage 1b — Layer 3 φ-ODE + Bare/Pre/Lam4 phenotype mapping: Sanity Gate

This document is the pre-implementation sanity review for Stage 1b, in which
Layer 3 (the φ-ODE for E-cadherin ↔ Integrin-β1 adhesion network dynamics)
is activated on top of the Layer 1 v15 + Stage 1a+ Option β α=1.0 + Stage
1a++ Layer 2 baseline. Bare/Pre/Lam4 phenotype mapping enters as initial
conditions on the φ field. Layers 4 and 5 remain off (Stages 1c, 1d).

The companion documents are:
- `docs/03_adhesion_dynamics.md` — pre-existing φ-ODE specification (Cho et
  al. 2020, Halbleib & Nelson 2006, Hynes 2002 references)
- `docs/stage1a_plus_plus_layer2_sanity.md` — Layer 2 sanity (Marchetti
  Rev Mod Phys 2013 IF 50 framework; ζ/K Option α' framing)
- `docs/SESSION_HANDOFF.md` §"Stage 1b framing memo" — Bare/Pre/Lam4 = formation
  environment phenotype, NOT substrate; spreading on single Col1 substrate
- `docs/12_validation.md` — Sanity-Gate Protocol (six checks) + Magic-Number
  Block

## Scope

### What is added (Layer 3 activation)

1. **Per-particle φ field** `φ_p ∈ [0, 1]`, evolving by the Cho-2020-inspired
   ODE pre-specified in `docs/03_adhesion_dynamics.md`:
   $$ \frac{d\phi_p}{dt} = k_+ \cdot S \cdot (1 - \phi_p) - k_- \cdot \phi_p $$
   Stage 1b simplification: substrate biochemistry is identical (Col1)
   across all three phenotype runs, so the substrate-induced signal is
   constant: `S ≡ 1` once contact is established (the Col1 substrate is
   always there from t=0 per Stage 1a+ initial-position scheme). The ODE
   reduces to
   $$ \frac{d\phi_p}{dt} = k_+ - (k_+ + k_-) \, \phi_p, $$
   relaxing to $\phi_{\text{eq}} = k_+ / (k_+ + k_-)$.
2. **φ → mechanical coupling**: per `03_adhesion_dynamics.md` §"Mechanical
   Coupling", any of γ_cc(φ), σ_active(φ), ρ_FA(φ), K_cortex(φ) may be
   active. **Stage 1b pilot activates only the σ_active(φ) coupling**
   (gives the cleanest mechanism question against the Stage 1a++ Layer 2
   baseline; the other three couplings remain available via flags but
   default off):
   $$ \sigma_{\text{act}, p} = -[\zeta_{\min} \cdot (1 - \phi_p) + \zeta_{\max} \cdot \phi_p] \cdot K \cdot I, $$
   replacing the constant `ζ_star` of Stage 1a++ with a φ-modulated
   value. ζ_min, ζ_max are config parameters set within the Stage 1a++
   Bucket B stable regime ([0.1, 0.4] per the Track 1 refine sweep).
3. **3 phenotype configs** (Bare / Pre / Lam4), differing only in
   `φ_initial` per the PI's formation-environment mapping:
   - `Bare` — pV4D4 alone, minimal laminin → minimal integrin engagement
     during formation → low initial integrin → φ_initial ≈ 0.30
     (E-cadherin dominant)
   - `Pre` — pV4D4 + pre-adsorbed laminin → moderate integrin → φ_initial
     ≈ 0.55 (balanced)
   - `Lam4` — pV4D4 + 4 µg/mL laminin in media → strong integrin
     engagement during formation → φ_initial ≈ 0.80 (Int-β1 dominant)
   Mapping framework: Cho et al. 2020 mechanism (laminin presentation
   drives E-cad → Int-β1 transition); specific φ_initial values are
   PI-domain-expertise estimates, no specific literature pinpoint
   (analogous to v15 ρ_floor / γ_sub_Col1 / ζ_star precedent — explicit
   honest disclosure in docstring).

### What stays off

- Layer 2 stochastic events (lamellipodia / filopodia / leader / FA
  discrete) — Stage 1a++.b
- Layer 4 (Marangoni / nematic / vortex) — Stage 1d
- Layer 5 (mechano-osmotic turgor) — Stage 1c
- Layer 6 (chemistry / necrosis) — Stage 2

### Why now (Stage 1a++ Layer 2 finding)

Stage 1a++ ζ/K sweep showed (i) Layer 2 boundary active stress is the
first mechanism to break the Layer-1 24% R-drift ceiling, (ii) optimum
exists around ζ/K ≈ 0.3, (iii) ζ/K = 1.0 destabilises the spheroid.
Stage 1b asks: **does phenotype-dependent ζ (driven by formation
environment per Cho 2020) explain the inter-condition variance the PI
observed in `data/experimental/260313_{Bare,Pre,Lam4}.csv`?**

---

## Sanity Gate (six checks)

### 1. Dimensional analysis

- φ_p — dimensionless (already specified in `03_adhesion_dynamics.md`).
- k_+, k_- — rate constants, units of 1/time. In dimensionless solver
  units (time = τ_relax = 60 s), `k_+_star = k_+ · τ_relax`.
- 03_adhesion_dynamics.md cites k_+ ≈ (12 hr)⁻¹ → `k_+_star = 60 / (12·3600) = 1.39e-3`.
- Stage 1b simplification: set `k_-_star = k_+_star / 3` to reach
  `φ_eq = k_+/(k_++k_-) = 0.75` at long times (matches Cho's 24-hr
  pV4D4 transition reaching ~0.5–0.8 ratio range).
- Pilot duration: 4 sim-hr ≈ 240·τ_relax. φ shifts about
  `1 - exp(-(k_++k_-)·t) = 1 - exp(-0.33) ≈ 28%` of the way from
  φ_initial toward φ_eq. Measurable but not full-equilibration.
- ζ(φ) coupling: `ζ_star(φ) = ζ_min·(1−φ) + ζ_max·φ`. Bounds inside
  Stage 1a++ Bucket B stable regime: provisional ζ_min = 0.1, ζ_max =
  0.4 (final values updated after Track 1 refine sweep completes).
- CFL: φ-ODE is purely per-particle algebraic (no spatial coupling),
  integrated by forward Euler. Stable iff `dt < 1/(k_+ + k_-)`. With
  k_+_star+k_-_star ≈ 1.85e-3 and dt_star = 0.01: dt·(k_+ + k_-) ≈
  1.85e-5 ≪ 1. Massively stable. ✓

**Check 1: PASS**.

### 2. Boundary cases

- `φ → 0` (E-cad dominant): `ζ(φ) → ζ_min = 0.1`. Recovers a Stage 1a++
  weak-cortex run (R drift ≈ 0.222 per Track 1 baseline). Well-defined.
- `φ → 1` (Int-β1 dominant): `ζ(φ) → ζ_max = 0.4`. Recovers a Stage
  1a++ stronger-cortex run (R drift expected ≤ 0.162 per Track 1 ζ=0.3
  finding extrapolated; Track 1 ζ=0.4 gives the actual value).
- `k_+ → 0`: φ stays at φ_initial. ULA-like (no transition). Valid
  baseline for diagnostic.
- `k_- → 0`: φ → 1 monotonically. Pure E-cad-loss limit. Valid for
  comparison.
- `k_+ → ∞ or k_- → ∞`: ODE stiff; forward Euler becomes unstable for
  `dt · max(k_+, k_-) > 1`. Constructor must check this; refuse runs
  where the bound is violated.
- `dt → 0`: φ converges exactly. ✓
- `dt → large`: stiffness check above catches it.
- N → 0 / N → ∞: per-particle field; same scaling as existing v15
  fields.

**Check 2: PASS** with constructor invariant `dt_star · (k_+_star + k_-_star) ≤ 0.5` (factor-2 safety margin against forward Euler stiffness bound).

### 3. Conservation invariants

- **Mass**: unchanged. ✓
- **Momentum**: φ-ODE has no momentum exchange (purely scalar
  per-particle dynamics). The φ-modulated ζ(φ) DOES contribute to
  momentum via the Layer 2 stress; this contribution is internal
  (boundary-particle stress pair, Newton-3 preserved by symmetric MPM
  scatter), same as Stage 1a++. ✓
- **Angular momentum**: same as Stage 1a++. ✓
- **Energy**: φ-ODE introduces no new energy term (φ is a state
  variable, not an energy reservoir; the mechanical coupling already
  goes through the existing active-stress channel). Energy-monotone
  gate stays SUSPENDED (Stage 1a++ contract change inherited).
- **φ conservation**: φ_p ∈ [0, 1] must be preserved per particle. The
  ODE form `dφ/dt = k_+(1−φ) − k_-·φ` has φ=0 as a stable boundary
  (k_+(1) = k_+ > 0 pushes φ up) and φ=1 as a stable boundary
  (−k_-(1) = −k_- < 0 pushes φ down). For finite dt, the forward Euler
  update `φ_{n+1} = φ_n + dt·(k_+(1−φ_n) − k_-·φ_n)` may overshoot if
  dt is large, but at our `dt·(k_+ + k_-) ≈ 1.85e-5 ≪ 1` the
  overshoot is bounded by `O(1.85e-5)` per step → completely negligible.
  Add a runtime clamp `φ_p := max(0, min(1, φ_p))` as defense-in-depth.

**Check 3: PASS** with the φ ∈ [0, 1] clamp.

### 4. Numerical sanity

- dt_star = 0.01 ≪ 1/(k_+_star + k_-_star) ≈ 540. ✓
- φ_initial values (0.30 / 0.55 / 0.80) are well inside (0, 1). ✓
- Float precision: f32 for φ field is fine (φ in [0,1], precision ~1e-7).
- Per-particle ODE adds O(N) work per step — same order as the existing
  per-particle loops; wall-clock impact < 5%.

**Check 4: PASS**.

### 5. Sign / sense check

- φ-ODE drift direction: at φ < φ_eq, dφ/dt > 0 (relaxes upward). At
  φ > φ_eq, dφ/dt < 0 (relaxes downward). Single stable fixed point at
  φ_eq. ✓
- ζ(φ) coupling: increasing φ (E-cad → Int-β1) increases ζ (more active
  stress). Biologically: integrin engagement upregulates lamellipodia /
  cortex activity per Cho 2020 mechanism. ✓
- Phenotype mapping: higher laminin presentation during formation →
  higher initial Int-β1 → higher φ_initial → higher initial ζ → stronger
  initial active stress → expected to spread MORE (per Cho 2020 / PI
  experimental Lam4 > Pre > Bare in A/A₀ final). Sign chain ✓.

**Check 5: PASS**.

### 6. Measurement-protocol consistency (NEW since v13)

- **R drift / shape metrics**: unchanged measurement protocols.
- **φ trajectory (new)**: per-particle φ_p logged at every diagnostic
  frame; report `<φ>_bulk`, `<φ>_boundary`, `<φ>_contact`,
  `min(φ)`, `max(φ)`. Walk-through: φ is purely scalar per particle;
  averages are well-defined. No off-protocol pathway. ✓
- **A/A_0 (spreading area, NEW for cross-comparison with PI experimental
  data)**: defined as `A_contact_xy_hull(t) / A_contact_xy_hull(t=0)`,
  using the same xy-hull measurement already in
  `acs.physics.mlsmpm.substrate_diagnostics()`. This is the proxy that
  matches the PI's experimental A/A₀ definition (projected area of the
  spreading spheroid). Walk-through: xy-hull is geometric, no curvature
  / off-peak issues; same protocol as Stage 1a+. **Reported as
  diagnostic**, not as a gate, since the PI's experimental data is
  *cross-check* not *fit target* per `docs/12_validation.md` Pillar 3.
- **Bucket-classification gates** (per outcomes doc): each per-run gate
  is a clear-cut pass/fail; sweep-level meta-bucket is determined by
  the 3-phenotype response pattern (see outcomes doc).

**Check 6: PASS**.

---

## Magic-Number Block

### `k_+_star`, `k_-_star` (φ-ODE rate constants)

| Parameter | Value | Anchor | IF |
|---|---|---|---|
| `k_+_star = k_+ · τ_relax` | 1.39e-3 (= τ_relax / 12 hr) | `docs/03_adhesion_dynamics.md` §"Rate Constants from Cho", citing Cho et al. 2020 (lab paper, ULA→pV4D4 timescale ~12 hr) | lab paper (cross-check only; not used for fitting) |
| `k_-_star = k_+_star / 3` | 4.6e-4 | Choice gives `φ_eq = 0.75` at long times, matching Cho 2020 Western blot endpoint range 0.5–0.8 | derived (geometric choice in Cho range) |

1. **Derivable** — yes, k_+ from Cho 2020 timescale (12 hr to reach
   near-equilibrium in pV4D4 transition; documented in
   `03_adhesion_dynamics.md` §"Rate Constants from Cho"). k_-/k_+ ratio
   from φ_eq target. Cho 2020 is the project's primary lab cross-check
   reference. Per `docs/12_validation.md` Pillar 3, Cho 2020 is allowed
   as a cross-check anchor (not as a fit target). PASS.
2. **Grid-invariant** — yes, rate constants in 1/time units, scale
   identically with τ_relax. Independent of dx, dt, n_particles. PASS.
3. **Fitting** — no. The rate values are derived from independent
   biological timescale measurement, not adjusted to make any gate
   pass. PASS.

### `φ_initial_{Bare, Pre, Lam4}` (formation phenotype mapping)

| Phenotype | φ_initial | Anchor |
|---|---|---|
| Bare | 0.30 | PI domain expertise: pV4D4 alone, minimal laminin presentation → minimal Int-β1 engagement during formation → mostly E-cad-dominant starting state, but with some basal Int (matching MCF7 baseline) |
| Pre | 0.55 | PI domain expertise: pre-adsorbed laminin → moderate Int engagement → balanced E-cad / Int |
| Lam4 | 0.80 | PI domain expertise: 4 µg/mL laminin in media → strong sustained Int engagement → mostly Int-β1 dominant |

1. **Derivable** — partial. Framework anchored to Cho 2020 mechanism
   (laminin presentation drives E-cad → Int-β1 transition); specific
   numerical pinpoints for Bare/Pre/Lam4 → φ_initial are PI domain
   expertise estimates, no specific literature pinpoint exists. Same
   honest-disclosure pattern as v15 ρ_floor / γ_sub_Col1 / ζ_star Option
   α'. **PI explicitly authorized 2026-04-29 ("PI input 사항 모두 PI 본인이
   직접 허락")** — this serves as the equivalent of Option γ' resolution
   for the φ_initial values. Documented honestly in the docstring. PARTIAL
   PASS with explicit PI authorization.
2. **Grid-invariant** — yes; φ is dimensionless, independent of dx, dt,
   n_particles. PASS.
3. **Fitting** — no. Values selected from biological reasoning before
   running the simulation, not chosen to make any gate pass. PASS.

### `ζ_min`, `ζ_max` (ζ-coupling bounds)

| Parameter | Value (provisional) | Anchor |
|---|---|---|
| ζ_min | 0.1 | Stage 1a++ Bucket B stable regime lower edge (ζ=0.1 → R drift 0.222) |
| ζ_max | 0.4 | Stage 1a++ Bucket B stable regime upper edge before instability (ζ=0.4 from Track 1 refine sweep; ζ=1.0 already shown to destabilise) |

Final values updated after Track 1 refine sweep completes; provisional
[0.1, 0.4] used in this draft.

1. **Derivable** — yes, from Stage 1a++ Bucket B finding (Layer 2
   ceiling break exists in [0.1, 0.4] ish range). PASS.
2. **Grid-invariant** — dimensionless ratio, PASS.
3. **Fitting** — no, derived from previous-stage stability boundary
   (not from any Stage 1b gate target). PASS.

**Magic-Number Block overall**: PASS for k_+, k_-, ζ_min, ζ_max.
PARTIAL PASS with PI explicit authorization for φ_initial_{Bare, Pre,
Lam4}.

---

## Per-run gates (each of the 3 phenotype runs)

Inherited from Stage 1a++ + 3 new Stage 1b additions:

| Gate | Status | Tolerance | Source |
|---|---|---|---|
| All Stage 1a++ gates (mass / momentum-horizontal / no-NaN / VRAM / boundary-tag stability / active-power-finite / R drift improvement vs β α=1.0 / contact-band ρ_kernel / anchor force balance / max-speed) | inherited | unchanged | Stage 1a++ sanity-md |
| Energy-monotone | SUSPENDED (Stage 1a++ contract change) | — | inherited |
| **φ ∈ [0, 1] per-particle invariant** | NEW GATE | min(φ) ≥ 0, max(φ) ≤ 1 at every frame | first principles (state-variable bounds) |
| **φ trajectory monotonicity (relax to φ_eq)** | NEW GATE | `\|<φ>_bulk(end) − φ_eq_predicted\|` ≤ 0.10 (within 28% of equilibrium shift expected per pilot duration) | Cho 2020 timescale |
| **A/A₀ trajectory finite & monotone** | NEW GATE | `A/A₀ > 0` always; `d(A/A₀)/dt ≥ 0` (allowed: spreading or stationary; forbidden: spheroid disintegration / contraction below initial contact area) | first principles (spreading is monotone in passive scenario) |

## Diagnostics (no gate, log-only)

- per-frame `<φ>_bulk`, `<φ>_boundary`, `<φ>_contact`, min(φ), max(φ),
  std(φ)
- per-frame `A/A₀(t)` from `contact_area_xy_hull / A_contact(t=0)`
- per-frame effective ζ_star = ζ_min·(1−<φ>_boundary) + ζ_max·<φ>_boundary
  (mean field-equivalent ζ to compare with Stage 1a++ point-ζ runs)

---

## Implementation outline

- **acs/physics/mlsmpm.py**:
  - SolverConfig: add `phi_initial: float`, `k_plus_star: float`,
    `k_minus_star: float`, `zeta_min: float`, `zeta_max: float`.
    All default to values that disable Layer 3 (`k_+ = k_- = 0`, `ζ_min
    = ζ_max = current zeta_star`) for backward compatibility with Stage
    1a++ configs.
  - New per-particle field `phi_p` (f32).
  - Constructor: initialize `phi_p[:] = phi_initial`. Constructor invariant:
    `dt_star · (k_+_star + k_-_star) ≤ 0.5` (forward Euler stability).
  - New kernel `_integrate_phi_ode`: per particle, forward Euler update
    `φ_p ← clip(φ_p + dt · (k_+ - (k_+ + k_-)·φ_p), 0, 1)`. Called once
    per `step()`, after `_g2p_and_constitutive` (so φ is updated in
    sync with the position update).
  - `_p2g_momentum_and_stress`: replace constant `zeta_K = ζ_star · K`
    with per-particle `ζ(φ_p)·K = (ζ_min·(1-φ_p) + ζ_max·φ_p)·K`.
  - `_compute_invariants`: extend with per-particle φ aggregation
    (sum, sum², min, max for std and bounds checks); n_boundary count
    for averaging.
- **acs/runner.py**: read `layer3` config block; gate evaluations for
  φ ∈ [0,1], φ trajectory toward predicted φ_eq, A/A₀ monotone.
- **3 configs**: `configs/stage1b_pilot_{bare, pre, lam4}.yaml` with
  identical Layer 1 / 2 settings, only `layer3.phi_initial` differing.

---

## Decision request to PI (resolved 2026-04-29 by full authorization)

PI granted **full authorization** 2026-04-29 for Stage 1b including:
- φ_initial estimates per Bare/Pre/Lam4 (Cho 2020 mechanism + PI domain
  expertise)
- Cho 2020 application as φ-ODE framework anchor
- All literature references and value ranges in this document
- Auto-progression through Track 2 implementation, pytest, 3-run pilot,
  result analysis, and commit

Stage 1b implementation begins under these conditions. Stop conditions
remain in force: no magic numbers (the φ_initial values are explicitly
PI-authorized estimates, documented honestly), no gate semantics edits
beyond the inherited Stage 1a++ contract changes, no v13 anti-pattern,
halt and surface to PI on any FAIL or critical error.

# Stage 2 — Layer 6 chemistry/ECM remodeling: pre-implementation Sanity Gate

This document is the pre-implementation sanity review for Stage 2, in
which Layer 6 (chemistry/ECM remodeling) is activated as the **final
layer** of the 5+1 framework. Stage 2 minimal scope per PI directive
2026-04-29:
- **MMP secretion** (cells secrete matrix metalloproteinases)
- **ECM degradation** (substrate effective adhesion weakens over time
  via MMP-mediated remodeling, modelled as time-varying scalar
  `ecm_strength(t)` multiplying the substrate CSF coefficient)
- **De novo ECM secretion** — sanity-md decision (recommendation: defer
  to Stage 2.b given exploratory framing)
- **Drug/signaling effects** EXCLUDED (future paper scope)

The companion documents are:
- `docs/02_force_models.md` §1.5 — substrate adhesion energy reference
- `docs/12_validation.md` §"Reserved Stages": Stage 2 Layer 6 hook
  (chemistry/necrosis)
- `docs/v1/stage1d_sanity.md` / `outcomes_stage1d.md` — inherited carriers
  baseline (Layer 1+2+3+4+5+Path C)
- `docs/12_validation.md` — Sanity-Gate Protocol + Magic-Number Block

## Stage 2 scope (exploratory framing per PI 2026-04-29)

### What is added (Layer 6 minimal)

1. **Global MMP concentration `mmp_total(t)`**: scalar accumulator,
   evolved per step:
   $$ \frac{d\,\text{mmp\_total}}{dt} = \alpha_{\text{MMP}} \cdot N_{\text{contact-band}}(t) $$
   Substrate-engaged cells (contact-band particles) secrete MMP at
   rate `α_MMP` per particle; total accumulates over time.
2. **ECM strength multiplier `ecm_strength(t)`**: decreases as MMP
   accumulates:
   $$ \frac{d\,\text{ecm\_strength}}{dt} = -\beta_{\text{deg}} \cdot \text{mmp\_total} \cdot \text{ecm\_strength} $$
   Initial `ecm_strength(0) = 1.0`; decays exponentially as MMP·time
   accumulates. Floor `ecm_strength ≥ ecm_strength_min` (default 0.1)
   to prevent runaway full-degradation in pilot duration.
3. **Coupling to substrate CSF**: substrate adhesion impulse magnitude
   scales with `ecm_strength(t)`:
   $$ \gamma_{\text{sub,eff}}(t) = \gamma_{\text{sub,star}} \cdot \text{ecm\_strength}(t) $$
   When ecm_strength → 1 (no degradation): identical to Stage 1d
   carrier. When ecm_strength → ecm_strength_min: substrate adhesion
   weakened, spheroid less anchored.
4. **De novo ECM secretion**: **DEFERRED** to Stage 2.b. PI exploratory
   framing accepts the simpler scope of degradation-only first.

### Active layers / what stays off

| Layer | Status | Notes |
|---|---|---|
| L1 bulk hydrodynamics (v15) | ON (carrier) | unchanged |
| L1a+ Option β substrate (γ_sub) | ON (carrier) | α=1.0, modulated by ecm_strength(t) |
| L1a++ Layer 2 active stress | ON (carrier) | overridden by Layer 3 |
| L1b Layer 3 φ-ODE | ON (carrier) | Pre default + spatial S |
| L1c Layer 5 mechano-osmotic | ON (carrier) | Tier 2 K(ρ_osm) |
| L1d Layer 4 Marangoni | ON (carrier) | core only |
| Path C effective gravity | ON (carrier) | g_star = 0.01 (Stage 1c) |
| **L2 Layer 6 chemistry minimal** | **NEW ACTIVATION** | MMP + ECM degradation only; de novo ECM deferred |
| Drug / signaling | OFF | future paper scope |

### Path C g_star recalibration check

Layer 6 weakens substrate adhesion over time but does NOT add new
forces beyond modulating γ_sub_star. Since the dominant force balance
on g_star (gravity vs surface tension) is unchanged, **g_star = 0.01
(Stage 1c recalibrated) carries over unchanged**. Note: as
`ecm_strength` decreases, the substrate becomes weaker and the
spheroid may begin to lift off (analogous to Stage 1a+ Option α
mechanism); gravity then dominates more and may pancake more. This is
the Bucket-Stage 2 mechanism question.

---

## Sanity Gate (six checks)

### 1. Dimensional analysis

- `mmp_total`: dimensionless concentration (arbitrary units; only the
  product `β_deg · mmp_total` matters in the ODE, so absolute scale is
  conventional).
- `α_MMP_star = α_MMP · τ_relax`: per-particle secretion rate per
  τ_relax. From Egeblad-Werb 2002 framework: MMP expression on the
  scale of minutes to hours. Choosing α_MMP ~ 1/(60 min) → α_MMP_star
  = 60/3600 = 0.0167 per τ_relax. Per-particle secretion: α_MMP_star
  scaled by `1/n_contact_band_typical ≈ 1/500` to keep total secretion
  rate comparable. **Provisional α_MMP_star = 1.7e-3** per τ_relax.
- `β_deg_star = β_deg · τ_relax`: ECM degradation rate per τ_relax. From
  Lu et al. 2011 (Nat Rev Mol Cell Biol IF 113) framework: tissue
  remodeling on the scale of hours. Choosing β_deg ~ 1/(2-4 hours) →
  β_deg_star = 60/(3·3600) = 5.6e-3 per τ_relax. **Provisional
  β_deg_star = 5e-3** per τ_relax.
- `ecm_strength`: dimensionless multiplier ∈ [ecm_strength_min, 1.0],
  default initial 1.0, floor 0.1.
- CFL: ODEs are first-order linear with rates ~ O(1e-3); forward Euler
  bound `dt · max(α_MMP, β_deg · mmp_max) ≪ 1` trivially.

**Check 1: PASS** with explicit dimensional self-correction lesson
applied (compare per-volume vs per-area carefully): Layer 6 ODEs are
purely scalar-rate, no force-balance comparison needed.

### 2. Boundary cases

- `α_MMP = 0` (no secretion): mmp_total stays 0, ecm_strength stays
  1.0; recovers Stage 1d exactly. PASS.
- `β_deg = 0` (no degradation): mmp_total accumulates linearly with
  time (since N_contact > 0), but ecm_strength stays 1.0. Recovers
  Stage 1d exactly. PASS.
- `α_MMP → ∞` or `β_deg → ∞`: ecm_strength → ecm_strength_min quickly;
  substrate adhesion becomes weak; spheroid likely lifts off (Stage
  1a+ pattern). Bounded by the floor. PASS.
- `ecm_strength_min → 0`: numerical risk of full degradation; floor
  prevents this. Default 0.1 keeps substrate ≥ 10% effective.
- N_contact → 0: secretion rate → 0; mmp_total stops growing.
- t → ∞ at finite α/β: ecm_strength relaxes to its long-time
  equilibrium ~exp(-β_deg · ∫mmp · dt²) — eventually approaches 0 but
  bounded by floor.

**Check 2: PASS** with the `ecm_strength_min` floor as defense-in-depth.

### 3. Conservation invariants

- **Mass**: per-particle mass unchanged. MMP is a *modeled chemical
  species*, not a mass-carrying particle. Mass conservation gate
  unchanged. ✓
- **Momentum**: substrate CSF impulse magnitude DECREASES over time
  (γ_sub_eff = γ_sub_star · ecm_strength(t)); this changes the
  per-step substrate impulse but does not break the scattering symmetry.
  Stage 1a+ horizontal-momentum gate stays in force; vertical leak
  through substrate is reduced (less anchor as ECM degrades). ✓
- **Angular momentum**: no chiral effect; same as Stage 1d.
- **Energy**: chemical energy expenditure (cells use ATP to secrete
  MMP) is implicit (not tracked). Gravitational PE + surface energy
  + KE + strain energy diagnostic unchanged. ECM degradation reduces
  the "potential" stored in substrate-cell adhesion (γ_sub_eff lower
  over time = less adhesion energy available); diagnostic
  surface_energy_substrate already captures the time-varying value.
- **MMP / ecm_strength evolution**: forward-Euler integrated; bounded
  by clamps. No explicit conservation expected (chemistry is
  inherently dissipative).

**Check 3: PASS**.

### 4. Numerical sanity

- dt_star = 0.01 unchanged; ODE rates ~ O(1e-3) → dt·rate ≪ 1
  (massively stable).
- f32 precision: ecm_strength ∈ [0.1, 1.0], well-resolved; mmp_total
  ranges from 0 to O(N_contact · α · t_total) ≈ 500 · 1.7e-3 · 240 ≈
  200 — within f32 dynamic range.
- Layer 6 ODEs are global (one scalar per quantity, not per particle);
  computational cost negligible (< 0.01% wall-clock).

**Check 4: PASS**.

### 5. Sign / sense check

- α_MMP > 0: mmp_total grows monotonically when N_contact > 0. ✓
- β_deg > 0: ecm_strength shrinks monotonically when mmp_total > 0. ✓
- Coupling γ_sub_eff = γ_sub_star · ecm_strength: as ecm_strength ↓,
  substrate adhesion ↓ → less anchoring force → spheroid more prone
  to lift-off (Stage 1a+ Option α pattern returns at long times under
  Layer 6 = realistic biological "cell digests substrate, becomes
  invasive").
- Mechanism: substrate-engaged cells secrete MMP → MMP degrades
  substrate locally → eventually spheroid loses anchor → can
  reposition / spread further / invade. This is the cancer-invasion
  framework (Egeblad-Werb 2002 framing). ✓

**Check 5: PASS**.

### 6. Measurement-protocol consistency (NEW since v13)

- **mmp_total(t), ecm_strength(t)**: per-frame scalar diagnostics.
  Aggregations well-defined. ✓
- **R drift / A/A₀**: same protocols as Stage 1d. Layer 6 expected to
  mildly INCREASE R drift over time (as substrate weakens, spheroid
  pancakes more or lifts off — depending on the gravity/substrate
  balance regime).
- **No new measurement protocols introduced**. The substrate CSF
  impulse pathway already exists (Stage 1a+ Option β); Layer 6 just
  modulates its magnitude.

**Check 6: PASS**.

---

## Magic-Number Block

### `α_MMP_star`, `β_deg_star`, `ecm_strength_min` (Layer 6 ODE constants)

| Parameter | Value | Anchor | IF |
|---|---|---|---|
| α_MMP_star | 1.7e-3 (= per-particle secretion rate · τ_relax) | Egeblad-Werb 2002 *Nature Reviews Cancer* (IF ~70) framework anchor: MMP expression on minutes-to-hours timescale; specific per-particle rate scaled by typical n_contact ≈ 500 to give total secretion rate ~ 1/(60 min) | 70 (framework only) |
| β_deg_star | 5e-3 (= per τ_relax) | Lu et al. 2011 *Nature Reviews Mol Cell Biol* (IF 113) framework anchor: tissue remodeling on 2–4 hour timescale (β ~ 1/3 hr in dimensional terms) | 113 (framework only) |
| ecm_strength_min | 0.1 | First-principles defensive floor (10% residual substrate adhesion); analogous to v15 ρ_floor = 0.1·ρ_ref pattern (no specific cite, conservative numerical bound) | — |

**Verification (2026-04-29)**: web search confirmed Egeblad-Werb 2002
*Nature Reviews Cancer* 2:161-174 paper exists (DOI 10.1038/nrc745;
established as foundational MMP-cancer review). Specific numerical
rates for MMP secretion + ECM degradation in the **per-particle** /
**dimensionless τ_relax** units used here are NOT directly pinpointed
in the cited reviews (which discuss qualitative mechanisms +
order-of-magnitude timescales). **Magic-Number Block PARTIAL**
analogous to ρ_floor / γ_sub_Col1 / ζ_star Option α' / α_osm /
gravity_star precedent: framework anchored to high-IF reviews
(Egeblad-Werb IF 70, Lu IF 113), specific dimensionless values are
order-of-magnitude derivations from cited timescales scaled to the
overdamped solver's τ_relax = 60 s calibration. PI full authorisation
2026-04-29 covers honest disclosure.

1. **Derivable** — partial. Framework derivable from Egeblad-Werb 2002
   (concept of MMP-mediated ECM remodeling) and Lu 2011 (tissue
   remodeling timescales). Specific α_MMP_star and β_deg_star values
   are dimensional-coefficient derivations matching the project's
   overdamped τ_relax = 60 s calibration. Same partial pattern as
   prior Magic-Number Block PARTIAL cases. PARTIAL PASS.
2. **Grid-invariant** — yes, dimensionless rate constants. PASS.
3. **Fitting** — no. Values pre-specified before pilot; not adjusted
   to make any gate pass. PASS.

### Time-scale mapping carry-over (Stage 1d sanity-md check 0)

τ_relax = 60 s (Moeendarbary 2013 IF 47 anchored). Pilot 4 hr ≈ 240·τ.
Layer 6 timescales:
- MMP accumulation: t_MMP ~ 1/α_MMP = 60 min → 1·τ_relax to substantially
  accumulate? Actually t = 1/(α_MMP · N_contact) = 1/(1.7e-3 · 500) =
  1.18 in τ_relax units = 70 s of sim time → very fast. Order of
  magnitude check: maybe rate is too high.
- ECM degradation: t_deg ~ 1/(β_deg · mmp_steady) = 1/(5e-3 · 200) =
  1.0 τ_relax = 60 s. Very fast.

**Issue**: ODE rates as written give ECM degradation in ~minutes,
much faster than the 2-4 hr literature anchor. The per-particle
α_MMP_star scaling needs revision: the literature timescale is for
*total tissue MMP* over the entire tissue, not per-cell rate.

**Refined provisional values** (preserving the literature 2-4 hr ECM
remodeling timescale as the END goal at fixed N_contact):
- Want ecm_strength to drop noticeably (say to 0.5) over t ≈ 2 hr =
  120 τ_relax.
- ecm_strength(t) ≈ exp(-β_deg · ∫mmp · dt) where mmp ~ α_MMP · N · t.
  So ln(2) ≈ β_deg · α_MMP · N · t²/2 at steady-state.
- For N=500, t=120: β_deg · α_MMP ≈ 2·ln(2)/(500·120²) = 1.9e-7.
- Choosing α_MMP_star = 1e-4 per τ_relax (very slow per-cell) and
  β_deg_star = 2e-3 per τ_relax: product 2e-7, gives ecm_strength →
  0.5 over ~120 τ_relax ≈ 2 hr. Matches Lu 2011.

**Updated provisional values**:
- α_MMP_star = 1e-4 per τ_relax (per particle MMP secretion)
- β_deg_star = 2e-3 per τ_relax (ECM degradation)
- ecm_strength_min = 0.1 (defensive floor)

In the 4-hr pilot, ecm_strength is expected to drop from 1.0 to ~0.3
(approaching the floor). Magic-Number Block PASS as recalibrated;
honest disclosure preserved.

---

## Per-run gates

Inherited from Stage 1d + 2 new Stage 2 additions:

| Gate | Status | Tolerance | Source |
|---|---|---|---|
| All Stage 1d inherited gates | inherited | unchanged | Stage 1d |
| Energy monotone | SUSPENDED (Stage 1a++) | — | inherited |
| **ecm_strength ∈ [ecm_strength_min, 1.0] invariant** | NEW GATE | min(ecm_strength) ≥ ecm_strength_min, max ≤ 1.0 | Magic-Number Block bounds |
| **mmp_total finite & non-decreasing** | NEW GATE | finite (not NaN, not > 1e6); mmp_total(t+1) ≥ mmp_total(t) | first principles (secretion is one-way in Stage 2 minimal scope) |

---

## Implementation outline

- **acs/physics/mlsmpm.py**:
  - SolverConfig: add `layer6_enabled`, `alpha_mmp_star`,
    `beta_deg_star`, `ecm_strength_min`, `ecm_strength_initial`. All
    default to disable Layer 6.
  - New scalar fields `mmp_total_field` (ti.field f64), `ecm_strength_field`
    (ti.field f32, scalar shape=()). Stored in solver state; updated by
    new kernel.
  - New kernel `_integrate_layer6_ode`: forward-Euler update of mmp_total
    using current `n_boundary` (or `n_contact_band` host-side reading);
    update ecm_strength with clamp.
  - `_grid_op_overdamped`: substrate CSF impulse uses `gamma_sub_star ·
    ecm_strength_field[None]` instead of constant `gamma_sub_star`.
  - `step()`: integrate Layer 6 ODE after `_g2p_and_constitutive` (similar
    to Layer 3, Layer 5).
  - `invariants()`: expose `mmp_total`, `ecm_strength`.
- **acs/runner.py**: read `layer6` config block; Layer 6 logger;
  2 new gates (ecm_strength invariant, mmp_total finite/non-decreasing).
- **One config**: `configs/stage2_pilot.yaml`. Full carrier (Layer
  1+1a+β+1a+++1b φ Pre+1c Layer 5+1d Layer 4+Path C g=0.01) + Layer 6
  minimal.

---

## Decision request to PI (resolved 2026-04-29 by full authorization)

PI granted **full authorisation** 2026-04-29 ("원래 framework standard 풀
적용. 축소 금지. 매 단계 풀 protocol") for Stage 2 including:
- Layer 6 minimal scope: MMP + ECM degradation; de novo ECM secretion
  DEFERRED to Stage 2.b.
- Egeblad-Werb 2002 IF 70 + Lu 2011 IF 113 framework anchors with
  PARTIAL Magic-Number Block (per ρ_floor / γ_sub_Col1 / ζ_star Option
  α' / α_osm / gravity_star precedent).
- Provisional α_MMP_star = 1e-4, β_deg_star = 2e-3, ecm_strength_min
  = 0.1 per dimensional rederivation against Lu 2011 2-4 hr remodeling
  timescale.
- Path C g_star = 0.01 (Stage 1c value) carries over unchanged.
- Drug / signaling effects EXCLUDED (future paper scope).

Stage 2 implementation begins under these conditions. Stop conditions
remain in force: no magic numbers (Layer 6 ODE constants are explicitly
PI-authorized PARTIAL per the verify pattern), no gate semantics edits
beyond inherited Stage 1a++ contract changes + the 2 new Stage 2 ODE
invariant gates, no v13 anti-pattern, halt and surface to PI on any
FAIL or critical error.

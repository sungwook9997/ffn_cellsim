# Stage 1c — Layer 5 mechano-osmotic (Tier 2): pre-implementation Sanity Gate

This document is the pre-implementation sanity review for Stage 1c, in
which Layer 5 (mechano-osmotic Tier 2 phenomenological coupling per
`docs/08_mechano_osmotic.md`) is activated on top of the Layer 1 v15 +
Stage 1a+ Option β α=1.0 + Stage 1a++ Layer 2 ζ=0.4 + Stage 1b Layer 3
φ-ODE + Path C effective gravity baseline. Layer 4 (Marangoni) and Layer
6 (chemistry) remain off (Stage 1d, Stage 2 respectively).

The companion documents are:
- `docs/08_mechano_osmotic.md` — pre-existing Tier 2 specification
  (state variable ρ_osm, ODE, mechanical couplings, Guo PNAS 2017 IF 12
  + Venkova eLife 2022 anchors)
- `docs/v1/path_c_sanity.md` / `docs/v1/outcomes_path_c.md` — Path C effective
  gravity (g_star recalibration is part of Stage 1c scope per PI
  directive)
- `docs/v1/stage1b_layer3_sanity.md` — inherited carriers (Layer 1+2+3)
- `docs/12_validation.md` — Sanity-Gate Protocol (six checks) + Magic-
  Number Block

## Stage 1c scope

### What is added (Layer 5 Tier 2)

Per `docs/08_mechano_osmotic.md` §"Tier 2 Formulation":

1. **Per-particle osmotic state ρ_osm,p** ∈ [0.5, 1.6], default 1.0 at
   equilibrium full hydration. Distinct from the v15 kernel-density
   `_rho_kernel_p` (which measures particle clustering); this is the
   per-cell water-content density.
2. **ODE per particle** (forward Euler each step):
   $$ \frac{d\rho_{\text{osm},p}}{dt} =
       \alpha_{\text{osm}} \cdot \dot\varepsilon_p^{\text{spreading}}
       - \beta_{\text{osm}} \cdot (\rho_{\text{osm},p} - 1) $$
   with $\dot\varepsilon_p^{\text{spreading}} = \max(0, -\mathrm{tr}(C_p))$
   (positive when cell volume contracts → water out → ρ_osm rises).
3. **Mechanical coupling — Stage 1c first-pass: K(ρ_osm) only**, per the
   Stage 1b precedent of activating ONE coupling channel first
   (Stage 1b activated σ_active(φ) only, leaving γ_cc(φ), K_cortex(φ),
   ρ_FA(φ) for future stages). Stage 1c activates:
   $$ K_{\text{eff}}(\rho_{\text{osm},p}) = K_\star \cdot \rho_{\text{osm},p}, $$
   modulating the v15 volumetric stress
   $\sigma_{\text{vol},p} = K_{\text{eff}}(\rho_{\text{osm},p}) \cdot
   (\rho_{\text{ref kernel}}/\rho_{\text{kernel},p} - 1) \cdot I$. The
   η_eff(ρ_osm) and σ_active(ρ_osm) couplings remain available via
   future config flags but default off.
4. **Per-particle clamp**: $\rho_{\text{osm},p} \in [\rho_{\text{osm,min}},
   \rho_{\text{osm,max}}] = [0.5, 1.6]$ enforced after every Euler
   update (defense-in-depth against transient overshoots).

### What is recalibrated (Path C g_star)

Path C v15 baseline pilot (commit `5dd8706`) showed `g_star = 0.1` was
~10× too large (sphericity 0.382 pancake). Self-corrected dimensional
analysis: balance regime is `g_star ≈ γ·κ / H = 0.028 / 2 = 0.014`.

For Stage 1c carrier (Layer 1 + 1a+ β α=1.0 + 1a++ ζ=0.4 + 1b φ + Layer
5 K(ρ)), the relevant balance includes **all force scales**:
- Surface tension: γ·κ = 0.028
- Substrate CSF: γ_sub_star · (1/dx) = 0.01 · 14.9 = 0.149 (per area, in contact band only)
- Layer 2 active stress: ζ·K = 0.4·1 = 0.4 (per volume, on boundary particles)
- Layer 5 K(ρ_osm) modulation: K_eff in [0.5, 1.6]·K (per volume)
- Required gravity: g_star · H_spheroid balances surface tension only
  (substrate CSF / Layer 2 don't pull spheroid up to lift-off; their
  net effect is in-plane / contractile, not lift-off)

Provisional **g_star = 0.01** for Stage 1c (10× smaller than Path C v15
baseline; matches the surface-tension-balance regime per the Path C
self-correction). The Layer 2 contractile cortex actively opposes
contraction (helps spheroid maintain shape); Layer 5 K(ρ) further
stiffens locally where ρ_osm rises. Both make the regime less prone to
over-anchoring than Path C alone.

### Active layers / what stays off

| Layer | Status | Notes |
|---|---|---|
| L1 bulk hydrodynamics (v15) | ON (carrier) | unchanged |
| L1a+ Option β substrate (γ_sub = γ_cc) | ON (carrier) | α=1.0 |
| L1a++ Layer 2 active stress | ON (carrier) | ζ_star = 0.4 (Track 1 stable optimum) |
| L1b Layer 3 φ-ODE | ON (carrier) | single phenotype (PI authorized scope; Bare/Pre/Lam4 sweep was Stage 1b) |
| L1c Layer 5 mechano-osmotic | **NEW ACTIVATION** | Tier 2; K(ρ_osm) only |
| Path C effective gravity | ON (carrier) | g_star = 0.01 (recalibrated) |
| Layer 4 (Marangoni / nematic) | OFF | Stage 1d |
| Layer 6 (chemistry / necrosis) | OFF | Stage 2 |

Default carrier φ_initial = 0.55 (Stage 1b "Pre" phenotype mid-range,
neutral-ish starting state; neither Bare nor Lam4 extremum).

---

## Sanity Gate (six checks)

### 1. Dimensional analysis

Solver remains in dimensionless units (length=R₀, time=τ_relax, stress=K).
New parameters (per `docs/08_mechano_osmotic.md` framework):

- $\rho_{\text{osm},p}$ — dimensionless [0.5, 1.6], default 1.0.
- $\alpha_{\text{osm}}$ — coupling strength, units of 1/time.
  Doc-cited value 0.5/hr in dimensional terms; in nondim time units
  ($\tau_{\text{relax}} = 60$ s):
  $$ \alpha_{\text{osm}}^\star = (0.5 / \text{hr}) \cdot (\tau_{\text{relax}} / \text{hr})
     = (0.5 / 3600) \cdot 60 = 8.33 \times 10^{-3}. $$
- $\beta_{\text{osm}}$ — relaxation rate, units of 1/time. Doc-cited
  value 1/(10 min) → $\beta_{\text{osm}}^\star = (1/600) \cdot 60 = 0.1$.
- ε̇^spreading proxy: $-\mathrm{tr}(C_p)$ has units of 1/time
  (velocity-gradient trace). Multiplied by $\alpha_{\text{osm}}^\star$
  gives dimensionless rate of ρ_osm change.

CFL: forward-Euler stability for the ODE `dρ/dt = α·ε̇ - β·(ρ-1)`:
- Stability bound $dt \cdot (β_{\text{osm}}^\star + α_{\text{osm}}^\star \cdot |\dot\varepsilon|^{\max}) \le 0.5$
  (factor-2 safety against Euler bound = 1).
- With $dt^\star = 0.01$, $β_{\text{osm}}^\star = 0.1$: $dt \cdot β = 10^{-3}$ ≪ 0.5.
- $|\dot\varepsilon|^{\max}$ in our pilots: typically $|tr(C)| \le 1$
  (deformation rate bounded). $dt \cdot α \cdot |\dot\varepsilon| \le
  0.01 \cdot 8.33 \times 10^{-3} \cdot 1 = 8.33 \times 10^{-5}$ ≪ 0.5.
- **Massively stable.** ✓

Path C g_star recalibration check:
- Surface tension force per area: γ·κ = 0.028.
- Required gravity force per area (gravity overcomes residual upward
  surface-tension pull): g_star · H_spheroid · ρ_star = g_star · 2 · 1
  = 2·g_star.
- Balance: 2·g_star ≈ 0.028 → g_star ≈ 0.014. **Provisional g_star = 0.01**
  (slightly under the balance for safety against pancake mode; the Layer 2
  contractile cortex provides additional outward push if needed).

**Check 1: PASS** with the dimensional-analysis self-correction
explicitly preserved (Path C lesson).

### 2. Boundary cases

- $α_{\text{osm}} \to 0$ AND $β_{\text{osm}} \to 0$: ρ_osm stays at
  initial value 1.0; K_eff = K (unchanged). Stage 1b reproduced
  exactly. PASS.
- $α_{\text{osm}} \to ∞$: ρ_osm jumps unboundedly per step; clamped at
  ρ_osm_max = 1.6. K_eff saturates at 1.6·K. No NaN. PASS.
- $β_{\text{osm}} \to ∞$: ρ_osm relaxes to 1.0 instantaneously; K_eff = K.
  Effectively disables Layer 5 dynamics (instantaneous re-equilibration).
- ρ_osm_min = 0.5 boundary: a particle continually expanding (positive
  tr(C) from above-baseline volume) sees `dρ/dt = -β·(ρ-1)`, relaxing
  ρ toward 1 from below; the clamp at 0.5 kicks in if a transient
  pushes ρ below. Defense-in-depth.
- ρ_osm_max = 1.6 boundary: similarly clamped on the high end. Per
  Guo PNAS 2017 anchor, 50% volume loss → ρ_max ≈ 1.6 is the cited
  experimental upper bound.
- N → 0 / N → ∞: per-particle field; same scaling as v15.
- Δt → 0 / Δt → large: stiffness invariant above.

**Check 2: PASS** with the per-particle ρ_osm ∈ [0.5, 1.6] clamp and
the construction-time stiffness invariant
$dt \cdot (\beta_{\text{osm}}^\star + \alpha_{\text{osm}}^\star \cdot \dot\varepsilon^{\max}) \le 0.5$.

### 3. Conservation invariants

- **Mass**: per-particle MASS unchanged (the cell still contains the
  same particles); per-particle ρ_osm changes the *internal water
  content* but the simulation's mass conservation (∑ m_p) is exact. ✓
  Note: if a future Tier 3 tracks water OUT of the system, mass
  conservation would need re-examination. Tier 2 keeps water "inside"
  the cell as a state variable.
- **Momentum (horizontal)**: ρ_osm doesn't directly contribute to
  forces; only modulates K(ρ_osm). The K modulation enters σ_vol via
  the existing `_p2g_momentum_and_stress` channel; momentum
  conservation is preserved by the symmetric MPM scatter (same as Layer
  1/2/3). ✓ Stage 1a+ Cousin-Rule horizontal-only momentum gate stays
  in force.
- **Momentum (vertical)**: same as Path C — substrate absorbs gravity.
- **Angular momentum**: K(ρ_osm) is isotropic in σ_vol; no torque
  contribution. ✓
- **Energy**: ρ_osm ODE introduces no new external energy reservoir
  (it's a phenomenological state variable, not an explicit thermodynamic
  variable). The K(ρ_osm) coupling DOES change σ_vol per particle, so
  the strain energy U_strain is implicitly ρ_osm-modulated:
  $U_{\text{vol}} = (1/2) K_{\text{eff}}(\rho_{\text{osm}}) \cdot
  (\rho_{\text{ref}}/\rho_{\text{kernel}} - 1)^2 \cdot V_0$. Update the
  diagnostic to use $K_{\text{eff}}$ instead of constant $K$. The
  energy-monotone gate (already SUSPENDED under Layer 2) stays
  suspended for Stage 1c.

**Check 3: PASS** with the strain-energy diagnostic update
(K → K·ρ_osm).

### 4. Numerical sanity

- $dt^\star = 0.01$ unchanged; ODE stiffness bound from check 1 is
  massively satisfied.
- f32 precision: ρ_osm in [0.5, 1.6], well within f32 dynamic range.
- K_eff in [0.5, 1.6] · K: stress field still bounded by ~1.6·K =
  $1.6 \cdot K_\star = 1.6$, well within f32 precision and the
  existing speed gates.
- Layer 5 ODE adds O(N) work per step — same order as Layer 3
  φ-ODE; wall-clock impact < 5%.

**Check 4: PASS**.

### 5. Sign / sense check

- ε̇^spreading_p = max(0, -tr(C_p)) > 0 only when cell volume is
  contracting (tr(C) < 0). Spreading cells flatten on substrate;
  in-plane spread + vertical compression typically gives net negative
  tr(C) → ε̇^spreading > 0 → water flows out → ρ_osm rises. ✓ Matches
  Guo PNAS 2017 mechanism (rounded → spread, volume decreases, water
  leaves).
- $-β·(ρ-1)$: when $ρ > 1$, this term is negative → ρ relaxes back
  toward 1. When $ρ < 1$, term is positive → ρ relaxes back upward
  toward 1. Single stable fixed point at ρ = 1 (full hydration) when
  $\dot\varepsilon^{\text{spreading}} = 0$. ✓
- K_eff(ρ_osm) = K · ρ_osm: when cell loses water (ρ_osm > 1), K
  increases → bulk modulus stiffens. Matches Guo's "G ∝ 1/V²" trend
  (decreased volume → increased stiffness). ✓ Sign chain:
  spreading → ρ_osm ↑ → K_eff ↑ → stronger bulk pressure → cortex
  resists further contraction. Physically meaningful negative feedback.

**Check 5: PASS**.

### 6. Measurement-protocol consistency (NEW since v13)

- **R drift**: same protocol; under Path C g_star = 0.01 (recalibrated)
  the spheroid is gravity-anchored without pancake. R drift gate
  semantic remains valid; expected to behave as Stage 1b R drift ~ 0.16
  modified by Layer 5 contribution.
- **Volume-change diagnostic** (NEW): per-particle ρ_osm aggregated as
  `<ρ_osm>_bulk`, `<ρ_osm>_boundary`, `<ρ_osm>_contact`. Walk-through:
  ρ_osm is purely scalar per particle; aggregations are well-defined.
  No off-protocol pathway. ✓
- **K_eff field**: per-particle effective bulk modulus = K · ρ_osm_p.
  Aggregation over boundary / contact / bulk shells; report as
  diagnostic.
- **Strain energy diagnostic** (updated): uses K_eff instead of constant
  K. The original v15 form used constant K; now it's per-particle
  K · ρ_osm_p. The energy-monotone gate is SUSPENDED under Layer 2/3,
  so this update is informational only.
- **A/A₀ trajectory**: under Path C g_star = 0.01 (no over-anchor),
  expected to stay > 0.5 with milder pancake; gate meaningful. ✓
- **Anchor force balance**: integrand still `(σ_vol + σ_act)_zz · dA +
  M_spheroid · g_star` per Path C; under K_eff modulation, the σ_vol
  integrand also changes (K → K·ρ_osm). Update the integrand
  computation to use per-particle K_eff. Tolerance unchanged.

**Check 6: PASS** with strain-energy / anchor-force integrand updates.

---

## Magic-Number Block

### `α_osm_star`, `β_osm_star`, `ρ_osm_min/max` (Layer 5 ODE constants)

Per `docs/08_mechano_osmotic.md` framework citing Guo et al. PNAS 2017,
114:E8618 (IF 12) and Venkova et al. eLife 2022, 11:e72381:

| Parameter | Value | Source | IF | Confidence |
|---|---|---|---|---|
| α_osm_star | 8.33e-3 (= 0.5/hr · τ_relax) | Guo PNAS 2017 §"timescale of cell spreading and volume loss" (~20 min for full spreading + concomitant volume loss; corresponds to α ~ 0.5/hr) — doc-attributed in `08_mechano_osmotic.md` | 12 | medium |
| β_osm_star | 0.1 (= 1/(10 min) · τ_relax) | Guo PNAS 2017 + Venkova eLife 2022 (osmotic re-equilibration ~ minutes) — doc-attributed | 12 / lab paper | medium |
| ρ_osm_min | 0.5 | Lower clamp; defensive (cell can't lose all water) | first principles |
| ρ_osm_max | 1.6 | Guo PNAS 2017 cites 50% volume loss → ρ_max = 1/(1-0.5) ≈ 2 in extreme; 1.6 is conservative within range | 12 |

Verification (2026-04-29): web search confirmed Guo PNAS 2017 paper
(DOI 10.1073/pnas.1705179114) and the qualitative mechanism (spreading
→ water efflux → volume loss in ~20 min, water leaves within a minute,
G ∝ 1/V² coupling). The specific numerical pinpoints α=0.5/hr and
β=1/10min are **doc-attributed** values in `docs/08_mechano_osmotic.md`
(synthesis from Guo + Venkova framework, not directly quoted from
either paper's abstract / methods sections accessible via web search).
**Same partial-pinpoint pattern as v15 ρ_floor / γ_sub_Col1 / ζ_star
Option α' / Path C g_star.** PI authorisation 2026-04-29 ("원래
framework standard 풀 적용") covers honest disclosure of this Magic-
Number Block PARTIAL.

1. **Derivable** — partial. Framework anchored to Guo PNAS 2017 +
   Venkova eLife 2022 (both directly cited in `02_force_models.md`
   and `08_mechano_osmotic.md`). Specific numerical pinpoints are
   doc-attributed; web-accessible portions of the cited papers
   confirm the qualitative scales (minutes to spread, water leaves
   within a minute, ~50% volume loss) consistent with the chosen
   numerical values. PARTIAL PASS — analogous to ζ_star Option α'
   precedent.
2. **Grid-invariant** — yes; α_osm_star, β_osm_star are dimensionless
   rate constants in 1/τ_relax units. Independent of dx, dt,
   n_particles. PASS.
3. **Fitting** — no. Values were specified in `docs/08_mechano_osmotic.md`
   BEFORE Stage 1c implementation (project pre-existing framework),
   not adjusted to make any gate pass. PASS.

### `gravity_star` (Path C, recalibrated for Stage 1c)

| Parameter | Value (Stage 1c) | Source |
|---|---|---|
| gravity_star | 0.01 | Path C dimensional self-correction (commit 5dd8706): balance regime g_star ≈ γ·κ / H = 0.028 / 2 = 0.014; provisional value 0.01 is slightly under-anchor for safety against pancake mode (with Layer 2/3/5 providing additional in-plane resistance to contraction) |

Same Magic-Number Block PARTIAL framing as Path C v15 baseline (body-
force coefficient under overdamped ξ_star = 1 calibration; framework
anchored to Stewart Nature 2011 IF 65 cell density). PI authorisation
2026-04-29 covers this recalibration.

**Magic-Number Block overall**: PARTIAL on α_osm_star / β_osm_star /
ρ_osm_min/max (doc-attributed; analogous to ζ_star Option α' precedent;
PI authorisation covers honest disclosure). PARTIAL on gravity_star
recalibration (Path C precedent; PI authorisation covers). All other
parameters inherited from previous stages.

---

## Per-run gates

Inherited from Stage 1b + 4 new Stage 1c additions:

| Gate | Status | Tolerance | Source |
|---|---|---|---|
| All Stage 1a / 1a+ / 1a++ / 1b inherited gates | inherited | unchanged | previous stages |
| Energy monotone | SUSPENDED (Stage 1a++) | — | inherited |
| **ρ_osm ∈ [0.5, 1.6] per-particle invariant** | NEW GATE | min(ρ_osm) ≥ 0.5, max(ρ_osm) ≤ 1.6 at every frame | Magic-Number Block |
| **<ρ_osm> trajectory finite & monotone-ish** | NEW GATE | `<ρ_osm>(end) ∈ [1.0, 1.6]` (monotone-rise from 1.0 baseline if spreading) | Guo 2017 mechanism |
| **K_eff diagnostic finite** | NEW DIAG | log only; report `<K_eff>_bulk` per frame | first principles |
| **Volume-change rate `d<ρ_osm>/dt` finite** | NEW DIAG | log only | first principles |

---

## Implementation outline

- **acs/physics/mlsmpm.py**:
  - SolverConfig: add `layer5_enabled: bool`, `rho_osm_initial: float`,
    `alpha_osm_star: float`, `beta_osm_star: float`, `rho_osm_min: float`,
    `rho_osm_max: float`. All default to disable Layer 5 (α=β=0,
    ρ_osm=1.0 fixed) for backward compatibility.
  - New per-particle field `rho_osm_p` (f32). Declared early in
    constructor (matching φ_p ordering pattern).
  - Constructor invariant: `dt_star · β_osm_star ≤ 0.5` (forward Euler
    stability with maximum α-driver bounded).
  - New kernel `_integrate_osmotic_ode`: per particle, forward-Euler
    `ρ_osm ← clip(ρ_osm + dt·(α·max(0, -tr(C)) - β·(ρ_osm - 1)), 0.5, 1.6)`.
    Called once per step after `_g2p_and_constitutive` (analogous to
    Layer 3 _integrate_phi_ode placement).
  - `_p2g_momentum_and_stress`: replace constant `K` in the volumetric
    stress term with `K_eff = K · rho_osm_p[p]` when Layer 5 enabled.
  - `_compute_invariants`: extend with per-particle ρ_osm aggregation
    (sum, sum², min, max). Update strain-energy term to use K_eff.
- **acs/runner.py**: read `layer5` config block; Layer 5 logger; 2 new
  gates (ρ_osm invariant, <ρ_osm> trajectory).
- **One config**: `configs/stage1c_pilot.yaml`. Layer 1 + 1a+ β α=1.0
  + 1a++ ζ=0.4 + 1b φ (Pre default) + Path C g_star=0.01 + Layer 5
  Tier 2 (K(ρ_osm) only).

---

## Decision request to PI (resolved 2026-04-29 by full authorization)

PI granted **full authorisation** 2026-04-29 ("원래 framework standard 풀
적용. 축소 금지. 매 단계 풀 protocol") for Stage 1c including:
- Layer 5 Tier 2 framework per `docs/08_mechano_osmotic.md` (Guo PNAS
  2017 IF 12 + Venkova eLife 2022 anchors; PARTIAL Magic-Number Block
  with honest disclosure analogous to previous stages).
- K(ρ_osm) coupling as the first-pass Stage 1c entry; η, σ_active
  couplings deferred.
- Path C `gravity_star` recalibration to 0.01 per the Path C dimensional
  self-correction.
- Single Stage 1c pilot (Pre carrier baseline; phenotype sweep deferred
  to Stage 1c.b if PI requests).
- Inherited Stage 1a++ contract changes (energy-monotone suspended,
  anchor-force integrand updated for ζ + gravity, contact-band ρ_kernel
  reinterpretation).

Stage 1c implementation begins under these conditions. Stop conditions
remain in force: no magic numbers (the ODE constants are explicitly
PI-authorized doc-attributed values per the existing Tier 2 framework),
no gate semantics edits beyond the inherited contract changes, no v13
anti-pattern, halt and surface to PI on any FAIL or critical error.

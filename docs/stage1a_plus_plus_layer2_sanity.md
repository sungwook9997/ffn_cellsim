# Stage 1a++ — Layer 2 boundary active stress: pre-implementation Sanity Gate

This document is the pre-implementation sanity review for Stage 1a++, in
which Layer 2 (boundary cell active stress) is activated on top of the
v15 (k.3) Layer 1 + Stage 1a+ Option β α=1.0 substrate-anchored baseline.
Stage 1b (Layer 3 φ-ODE) and beyond remain off. Implementation does not
start until the PI signs off on this proposal.

The companion documents are:
- `docs/stage1a_plus_substrate_sanity.md` — Layer 1 + Col1 substrate Sanity
  Gate (Option α / Option β), already PI-approved
- `docs/outcomes_stage1a_plus.md` — bounded outcomes for Option α / β
- `docs/stage1a_interior_pressure_sanity.md` — v15 Layer 1 (k.3) volumetric
  stress proposal
- `docs/outcomes_v15.md` — v15 Path B decision (Layer 1 single-layer limit)
- `docs/v15_deferred_diagnostics_plan.md` — Track B reference (this stage's
  trigger conditions for re-entry)
- `docs/02_force_models.md` §1.4 — active bulk stress reference (Marchetti
  Rev Mod Phys 2013, IF 50)
- `docs/12_validation.md` — Sanity-Gate Protocol (six checks) + Magic-Number
  Block

## Stage 1a++ scope

### What is added

Layer 2 boundary cell active stress, in its simplest entry form (no
lamellipodia / filopodia stochastic events; those are Layer-2-discrete
mechanisms deferred to a later sub-stage). The model:

- **Stress contribution**: per particle tagged as boundary
  (`is_boundary[p] == 1`, identical density-based tagging used through
  v15), add an additional Cauchy stress
  $$ \sigma_{\text{act},p} = -\zeta_a \cdot I, $$
  where $\zeta_a$ is the active stress coefficient (Marchetti 2013) and
  the negative sign encodes **contractile** cortex behaviour (compressive
  stress in the cortex band).
- **Where applied**: in `_p2g_momentum_and_stress`, the per-particle
  Cauchy stress is augmented as
  `stress = stress_vol + tau_dev + (-ζ_star · I if is_boundary[p] else 0)`.
- **Layer 2 stochastic events (lamellipodia / filopodia / leader cells /
  discrete focal adhesions per `docs/02_force_models.md` §2.1–2.4)**:
  **deferred to Stage 1a++.b** (a later sub-stage). Stage 1a++ tests the
  simplest continuum Layer 2 entry first.

### Carrier baseline (Layer 1 + Stage 1a+ Option β α=1.0)

The Stage 1a+ Option β sweep showed all three α values produced
indistinguishable R drift (24.6–24.7%); α=1.0 maintained the longest
substrate contact (18 of 24 frames, anchor force balance ~0.97). For
Stage 1a++ we use **α=1.0 as the carrier configuration** — the
best-anchored Layer-1+substrate baseline. Active stress coefficient
$\zeta_a$ is then the only added knob.

### What stays off

- Layer 2 stochastic events (lamellipodia / filopodia / leader cells /
  discrete focal adhesions) — deferred to Stage 1a++.b
- Layer 3 (φ-ODE for E-cad ↔ Int-β1) — Stage 1b; per
  `docs/SESSION_HANDOFF.md` §"Stage 1b framing memo", this is where
  Bare/Pre/Lam4 phenotype mapping enters.
- Layer 4 (internal flow / Marangoni / nematic) — Stage 1d
- Layer 5 (mechano-osmotic turgor) — Stage 1c
- Layer 6 (chemistry / necrosis) — Stage 2

These remain togglable via the `layers:` block in YAML; Stage 1a++ sets
only `L1_bulk_hydrodynamics: true` and `L2_boundary_biology: true`, all
others `false`.

### Why now (5-simulation Layer 1 ceiling)

The five Stage 1a / 1a+ simulations (v15 free / Option α mech / β α=0.1
/ 0.3 / 1.0) yielded R drift 0.244, 0.247, 0.246, 0.246, 0.247 — a
1%-wide band that does not depend on substrate condition. Layer 1 alone
(with or without substrate, with or without energetic wetting up to
γ_sub = γ_cc) cannot reach the Laplace equilibrium R/R₀ ≈ 0.99. The
mechanism question for Stage 1a++ is:

> *Does adding boundary cell active stress (continuum entry of Layer 2)
> reduce the Layer 1 architectural ceiling? If so, by how much, and as
> a function of $\zeta_a / K$?*

This is a layer-by-layer narrative: each new layer's contribution to
spheroid integrity is quantified.

---

## Sanity Gate (six checks)

### 1. Dimensional analysis

Solver remains in dimensionless units (length=R₀, time=τ_relax, stress=K).
New parameters:

- $\zeta_a$ — active stress coefficient. Units: stress = $K_\star \cdot K$
  in SI. Dimensionless: $\zeta_\star = \zeta_a / K$.
- Marchetti Rev Mod Phys 2013 (IF 50) reports $\zeta_a \sim 100\text{–}1000$ Pa
  for cell-tissue active stress. With $K = 1$ kPa (Stage 1a placeholder),
  $\zeta_\star \in [0.1, 1.0]$.
- Activity number $A = \zeta_\star \cdot \tau_\star = \zeta_\star \cdot 1
  = \zeta_\star$ (since $\tau_\star = 1$ in our units): the ratio of
  active-stress timescale to viscous-relaxation timescale.

CFL: active stress contributes a *constant* per-particle stress addition
(no spatial gradients explicitly, only via the boundary tag changing
over time as particles re-classify). The effective volumetric+active
stress at a boundary particle is bounded by $|σ_{\text{vol}}| +
\zeta_\star \le K \cdot 0.5 + \zeta_\star = K(0.5 + \zeta_\star)$. For
$\zeta_\star = 1$ this is $1.5 \cdot K$; for $\zeta_\star = 0.3$ it is
$0.8 \cdot K$. Both well within the f32 stress-field range and below any
stability bound that would invalidate the existing $dt_\star = 0.01$.
The effective elastic-wave speed in the overdamped limit is irrelevant
(no inertia); the diffusion-like CFL $dt_\star \le dx_\star^2 \cdot
\xi_\star / \mu_\star$ depends only on $\mu_\star$ (deviatoric) and is
unchanged by adding an isotropic active term.

**Action**: log $\zeta_\star$ in the run manifest; cross-check against
the $[0.1, 1.0]$ Marchetti range.

**Check 1: PASS**.

### 2. Boundary cases

- $\zeta_\star \to 0$: the active term vanishes; solver reduces exactly
  to Stage 1a+ Option β (whatever α the carrier is set to). PASS by
  construction.
- $\zeta_\star \to \infty$: active stress dominates. Boundary particles
  receive ever-larger compressive impulse → potential instability via
  the per-particle stress overflow. Mitigated by (i) the existing
  `max_speed_over_vrms` runtime gate, (ii) the literature-anchored
  range capping $\zeta_\star \le 1$. **No silent NaN**: if a runaway
  occurs, the existing NaN-detection in `_compute_invariants` halts the
  run.
- N → 0 / N → ∞: no change vs Stage 1a+; constructor invariants
  unchanged.
- $\Delta t \to 0$: same as Stage 1a+; active stress is integrated
  algebraically (no time-derivative term added).
- $\Delta t \to \text{large}$: the active impulse is $-\zeta_\star \cdot
  V_0 \cdot dt \cdot \nabla w$ via the standard MPM stress-to-grid
  coupling; it scales with $dt$ in the same way as the existing
  volumetric and deviatoric impulses, so the same CFL bound applies.

**Check 2: PASS** with the literature-anchored $\zeta_\star \in [0.1, 1.0]$
range understood as the operating regime.

### 3. Conservation invariants

- **Mass**: unchanged (no particle creation; same MLS-MPM construction). ✓
- **Momentum**: active stress adds an *internal* force pair (boundary
  particle pushes its neighbours apart compressively; reaction comes
  from Newton-3 via the symmetric MPM stress scatter). Total momentum
  remains conserved to atomic-add round-off. The horizontal-only
  contract change from Stage 1a+ stays in force (substrate vertical
  leak). PASS.
- **Angular momentum**: APIC preserves it; isotropic active stress
  ($\propto I$) does not break rotational symmetry. ✓
- **Energy**: active stress does *positive work* on the system (active
  matter is by definition out of equilibrium, draws energy from chemical
  reservoirs ATP→ADP). The existing energy-monotone gate (KE +
  U_strain) cannot be expected to monotone-decrease in the presence of
  active stress. **Contract change required**: the energy-monotone gate
  is **suspended** for Stage 1a++ runs and replaced with the diagnostic
  $W_{\text{active}}(t) = \int_0^t \Sigma_{\text{boundary}}
  (\zeta_a \cdot V_0) \cdot \nabla \cdot v \, dt'$ (work done by active
  stress per step) reported as a separate informational metric. Cousin-
  Rule contract change surfaced here.
- **Possible new leaks**: (i) boundary tag re-classification mid-step
  could cause a particle's active stress to switch on/off discontinuously
  — but boundary tag is computed once per step in `_tag_boundary` and
  used consistently throughout, so no within-step inconsistency.
  (ii) Atomic-add round-off in stress scatter: same as Stage 1a+, no
  new channel.

**Check 3: PASS** with the energy-monotone-gate suspension contract
change explicitly documented (PI sign-off required).

### 4. Numerical sanity

- $dt_\star$ unchanged (no new fast scale; isotropic stress integrates
  algebraically).
- Grid resolution: active stress is per-particle, scattered via the same
  3³ APIC stencil as the existing volumetric stress; same nominal grid
  accuracy.
- Float precision: $\zeta_\star \le 1$ keeps $|σ_{\text{stress}}| \le
  1.5 \cdot K_\star$, comfortably inside f32 dynamic range (existing
  $|σ_{\text{vol}}|$ already in $[-K, +K]$).
- Boundary tag dynamics: the density-based tag may flicker as particles
  cross the threshold; for a mostly-equilibrated spheroid this is
  expected to be stable. Diagnostic: report
  $\Delta n_{\text{boundary}}/\text{frame}$ to verify the tagging is not
  flickering pathologically.

**Check 4: PASS** with the boundary-tag-stability diagnostic added.

### 5. Sign / sense check

Active stress sign convention:
- $\zeta_a > 0$ contractile (cortex pulls itself together). Cauchy stress
  contribution $-\zeta_a \cdot I$ is **compressive**. In MPM, compressive
  stress at a particle $p$ produces a grid force $-V_0 \cdot \sigma \cdot
  \nabla w$ that pushes neighbouring grid points **apart** (positive
  divergence of impulse) — i.e. the cortex *resists compression*. ✓
- For a free-floating spheroid, the boundary band experiences active
  compressive stress that **opposes** the surface-tension contraction
  pulling the spheroid inward. Mechanism: cortex stiffness adds to the
  bulk pressure response, helping to balance Laplace pressure. This is
  the *physical* expectation that motivates Stage 1a++.
- For a substrate-anchored spheroid (Option β α=1.0 carrier), the
  substrate band particles already experience the substrate CSF impulse
  pulling toward $z=0$; adding contractile active stress at the
  *boundary* (which includes the substrate contact band as part of the
  free-surface boundary tag, depending on the tag definition) creates a
  competition between "cortex wants to compress" and "substrate wants
  to pull down". Resolution depends on $\zeta_\star$ vs $Ca_{\text{sub}}$.
- **Edge case**: the boundary-tag definition currently captures any
  particle with ≥ 5/27 empty neighbours. For a substrate-contact
  particle, the substrate side is technically empty (z<0 has no grid
  cells with mass), so contact-band particles are *also tagged as
  boundary*. The active stress will therefore apply to substrate-contact
  particles too. **Walk-through**: this is physically correct — cortex
  is present at the basal surface as well as the apical surface — but
  the *Stage 1a+ Option β α=1.0 carrier diagnostics* (anchor force
  balance, contact-band ρ_kernel) will now include an active-stress
  contribution. The anchor-force-balance gate must be *recomputed* with
  the active term included: $F_{\text{substrate}} = \int_{\text{contact}}
  (\sigma_{\text{vol},zz} + \sigma_{\text{act},zz}) \, dA$.
  Sub-decision required from PI: is the anchor-force-balance gate
  semantic preserved (Newton-3, "substrate reaction = total Cauchy
  stress at the contact band") or relaxed?

**Check 5: PASS** with the explicit sign rationale and the anchor-force-
balance-gate recomputation requirement (PI sign-off).

### 6. Measurement-protocol consistency (NEW since v13)

- **R drift** (existing gate, also "improvement vs v15 baseline"):
  measurement is unchanged (effective_radius from particle-position
  second moment). The active stress acts on boundary particles, but the
  measurement averages over *all* particles. Walk-through: no off-
  protocol pathway. ✓
- **Contact-band ρ_kernel**: under Stage 1a++ the contact band's
  particles experience active stress in addition to substrate CSF and
  bulk pressure. The kernel-density measurement itself is unchanged
  (G2P from grid_m), but its *interpretation* shifts: $\rho_{\text{kernel
  contact}}$ now reflects the equilibrium between three forces, not
  two. The gate window $[0.85, 1.15]$ remains a witness of "kernel sees
  the substrate-anchored cortex" rather than just "substrate-anchored".
  Tolerance unchanged (window is not a tightened gate). ✓
- **Anchor force balance**: as noted in check 5, the integrand changes
  to include $\sigma_{\text{act},zz}$. Measurement protocol updated;
  gate value tolerance unchanged. PI sub-decision: confirm gate
  reinterpretation.
- **Active work diagnostic** (new under Stage 1a++): per step,
  $W_{\text{act,step}} = \sum_{p \in \text{boundary}} (-\zeta_\star)
  \cdot V_0 \cdot \text{tr}(C_p) \cdot dt$. Time-integrate to get
  cumulative active work. Report as `active_work_star` in metrics.csv.
  No gate; informational diagnostic that lets the energy-monotone
  suspension be reconciled with first-law accounting.

**Check 6: PASS** with the anchor-force-balance recomputation,
contact-band ρ_kernel reinterpretation, and the new active-work
diagnostic.

---

## Magic-Number Block (mandatory, per `docs/12_validation.md`)

### `ζ_star` — verification partial; PI selected Option α' (ratio sweep, no Pa claim)

**Status (resolved 2026-04-29)**: PI selected **Option α'** below
(ratio sweep ζ/K only, no Pa claim, anchored to K via Fischer-Friedrich
Nat Cell Biol 2014 IF 30, framework reference Marchetti Rev Mod Phys
2013 IF 50). Stage 1a++ implementation proceeds with this framing.

The `ζ_star` parameter enters the code as a *dimensionless ratio* `ζ/K`
swept across `{0.1, 0.3, 1.0}`. **No claim about a specific MCF7 active-
stress magnitude in Pa is made**; the result is a response curve in
ζ/K. This is structurally identical to the Stage 1a+ Option β framing
the PI already approved (`Ca_sub = α·Ca_cc`, K-anchored, sweep is the
result). Magic-Number Block automatically passes (no specific Pa value
to verify; sweep design is parameter-free at result level).

The verification record below is preserved as the audit trail for why
Option α' is the chosen path.

### Verification record (preserved for audit trail)

Per the v15 ρ_floor and γ_sub_Col1 verification precedent, the value
must pass the three Magic-Number Block tests against the cited
references *before* implementation. The previous draft of this section
attributed the numerical range `ζ_a ∈ [100, 1000] Pa` (i.e.
`ζ_star ∈ [0.1, 1.0]` with K = 1 kPa) directly to Marchetti
Rev Mod Phys 2013 (IF 50). A literature verification was performed
(2026-04-29).

#### Verification record

| Reference | What it actually contains | Match for `ζ_a` numeric range? |
|---|---|---|
| Marchetti et al. *Hydrodynamics of soft active matter*, Rev Mod Phys 2013, 85:1143 (IF 50) | Comprehensive theoretical review of the **active-stress framework** (concept, hydrodynamic formulation, polar/nematic order parameters). Web-accessible portions (abstract, ToC, Boulder School lecture-note PDF) **do not contain a pinpointed numerical range for ζ_a in cell/tissue Pa or for ζ/K dimensionless**. The 50-page paper may cite specific values in its applications section, but no excerpt visible to web search confirms `[100, 1000] Pa`. | **framework anchored, specific value range NOT pinpointed in web-accessible portions** |
| `docs/02_force_models.md` §1.4 in this project | Existing attribution: "ζ_a ~ 100–1000 Pa, Marchetti Rev Mod Phys 2013, IF 50" | **synthesis citation** (the project author wrote the range as "from the cell-tissue active-matter literature broadly, attributed to Marchetti as the framework anchor", not as a direct Marchetti quote) |

This is structurally similar to the γ_sub_Col1 verification failure
(2026-04-29 record above): the framework reference is solid, the
specific numerical pinpoint is not. **Magic-Number Block test 1
(Derivable) result is PARTIAL: framework derivable, specific range
NOT pinpointed in the cited primary reference.** Tests 2 (grid-
invariant) and 3 (fitting) PASS as in the previous draft.

A test-1 partial result blocks the change unless PI selects an
alternative path (analogous to the γ_sub_Col1 → Option α/β/γ
resolution that the PI selected in 2026-04-29 for the Stage 1a+ stage).

#### Alternatives surfaced to PI (Layer 2 activity sub-decision — PI selected Option α')

Three options analogous to the substrate γ_sub_Col1 resolution. The
naming uses primes to distinguish from the substrate options.

**Option α' — Drop the numerical-magnitude claim entirely; sweep
dimensionless `ζ/K` only.**
- Report results as a *response curve* in `ζ/K`, not as a claim about
  a specific MCF7 ζ_a value in Pa.
- Anchor: K is well-anchored (Fischer-Friedrich Nat Cell Biol 2014,
  IF 30, already cited in `docs/02_force_models.md` §1.1 as primary
  for cortex stiffness ~ 1 kPa).
- Sweep variable: `ζ/K ∈ {0.1, 0.3, 1.0}` is a logarithmically-spaced
  decade exploration, *parameter-free at result level* (matches the
  Stage 1a+ Option β framing for substrate `α = γ_sub/γ_cc`).
- Magic-Number Block: PASS by construction (no specific Pa value is
  claimed; the sweep is the result, not a parameter to verify).
- Trade-off: cleanest framing; no ζ_a Pa value to defend at paper-
  writing time. Identical structural pattern to substrate Option β.

**Option β' — Cite secondary literature for the specific Pa range.**
- Find a peer-reviewed paper that *does* explicitly give a numerical
  range for active stress in cells/tissues (e.g. Saw et al. Nature
  2017 IF 65 for active nematics in epithelia, Banerjee-Marchetti 2014
  Soft Matter, Köpf-Pismen 2013, etc.) and use that as the primary
  Pa anchor; demote Marchetti 2013 to framework reference.
- Trade-off: requires another literature pass to verify the secondary
  cite is itself a direct numerical claim, not another synthesis
  attribution. Risk of recursive verification failure.

**Option γ' — Keep `ζ_star ∈ {0.1, 0.3, 1.0}` with explicit "no
specific literature reference for the Pa values; framework anchored to
Marchetti 2013, dimensionless sweep range from order-of-magnitude
synthesis of soft-active-matter literature for cell-tissue systems"
docstring (v15 ρ_floor / γ_sub_Col1 precedent).**
- Trade-off: matches v15 / γ_sub_Col1 honest-disclosure precedent;
  propagates one weakly-pinpointed numerical range claim into the
  Stage 1a++ documentation.

**Recommendation**: Option α' (sweep ζ/K only, no Pa claim). It is
strictly cleaner than γ' (no value-pinpoint claim made) and avoids
the recursive-verification risk of β'. It is exactly analogous to the
substrate Option β framing the PI already approved (`Ca_sub = α·Ca_cc`,
result is a response curve, K-anchored). The sweep itself becomes the
result, not a claim about an unmeasured MCF7 active-stress magnitude.

**Resolved 2026-04-29: PI selected Option α'** (ratio sweep ζ/K only,
no Pa claim, K-anchored to Fischer-Friedrich Nat Cell Biol 2014
IF 30; Marchetti Rev Mod Phys 2013 IF 50 retained as framework
reference for the active-stress concept and polar/nematic order
parameter formalism). Stop conditions remain in force.

---

## New gate candidates (Stage 1a++)

These are *added* to the Stage 1a+ gate set; the existing gate
tolerances and normalisations are not modified except for the explicit
contract changes recorded above (energy-monotone suspension; anchor-
force-balance reinterpretation under active stress).

| New gate | Status | Tolerance | Source |
|---|---|---|---|
| Stage 1a++ R drift improvement vs Stage 1a+ Option β α=1.0 baseline | **GATE** | strictly less than 0.247 | β α=1.0 result (`docs/outcomes_stage1a_plus.md`) |
| Boundary-tag stability `\|Δn_boundary / n_boundary\|` per frame | **GATE** | ≤ 0.10 (10% flicker tolerance) | first principles (numerical hygiene) |
| Active work cumulative `W_active(t)` finite at end of run | **GATE** | finite (not NaN, not > 10·U_strain at end) | first principles (sanity of work integration) |
| Energy-monotone (KE + U_strain) | **SUSPENDED** | — | active stress injects energy by construction |
| Anchor force balance | **GATE (recomputed)** | ≤ 0.20 with `(σ_vol + σ_act)_zz` integrand | Stage 1a+ Option β unchanged tolerance |
| Apparent contact angle θ | **DIAGNOSTIC (active-modified)** | log only | log for comparison with passive Option β |

The R-drift improvement gate is set against the **Stage 1a+ Option β
α=1.0** baseline (R_drift = 0.247) rather than the v15 baseline
(0.244), because Stage 1a++ uses the α=1.0 carrier and the comparison
should isolate the Layer 2 contribution alone, not re-test the
substrate effect.

---

## Implementation options (sanity-md only; PI decision)

### Option I — Single mid-range pilot (`ζ_star = 0.3`)

- One run, fastest mechanism check.
- Code delta: ~20 lines in `_p2g_momentum_and_stress` (per-particle
  active term gated on `is_boundary[p]`); ~10 lines in
  `_compute_invariants` (active-work diagnostic accumulation); one
  YAML field `layer2.zeta_star`.
- Wall-clock impact: < 5%.
- Pros: fastest readout; clean comparison with Option β α=1.0 carrier
  baseline.
- Cons: a single point doesn't characterize the $\zeta_\star$
  dependence; if the result is in a transitional regime, a second run
  may be needed anyway.

### Option II — 3-run sweep (`ζ_star ∈ {0.1, 0.3, 1.0}`)

- Three runs, analogous to Stage 1a+ Option β.
- Same code as Option I, three configs.
- Wall-clock impact: 3× Option I (still ≤ 3 min total).
- Pros: characterises the response curve directly; matches Stage 1a+
  Option β reporting style; if Layer 2 reduces R drift at all, the
  sweep tells us how much per unit $\zeta_\star$.
- Cons: 3× pilot time; if all three give the same R drift (analogous
  to Option β meta-bucket 3), the sweep was diagnostic but not
  parameter-defining.

**Recommendation**: Option II (3-run sweep). The Stage 1a+ Option β
sweep produced a defensible quantitative finding (R drift insensitive
to substrate energy across a decade of α). A Stage 1a++ Option II
sweep would do the same for $\zeta_\star$, giving the layer-by-layer
narrative its third quantitative pillar.

---

## Summary

- **All six Sanity-Gate checks PASS** with three contract changes
  surfaced for explicit PI sign-off: (i) energy-monotone gate suspended
  (active stress injects energy by construction); (ii) anchor-force-
  balance integrand recomputed to include $\sigma_{\text{act},zz}$ at
  contact-band particles; (iii) contact-band ρ_kernel gate window
  reinterpretation.
- **Magic-Number Block resolution**: PI selected Option α' (resolved
  2026-04-29). Sweep `ζ/K ∈ {0.1, 0.3, 1.0}` is the result, not a
  parameter to verify; no specific Pa value is claimed. K is anchored
  via Fischer-Friedrich Nat Cell Biol 2014 IF 30 (already cited in
  `docs/02_force_models.md` §1.1 as primary for cortex stiffness ~1 kPa).
  Marchetti Rev Mod Phys 2013 IF 50 retained as the framework reference
  for the active-stress concept. Magic-Number Block automatically PASS.
- **Carrier baseline**: Layer 1 v15 + Stage 1a+ Option β α=1.0 (best-
  anchored substrate config from the 5-simulation Layer 1 ceiling
  finding).
- **Layer 2 entry**: continuum boundary active stress only (no
  lamellipodia / filopodia / leader-cell stochastic events; those are
  Stage 1a++.b deferred sub-stage).
- **Implementation recommended**: Option II (3-run sweep `ζ_star ∈
  {0.1, 0.3, 1.0}`); ~30 lines of code; wall-clock ~3× Stage 1a+ pilot.
- **New gate set**: 3 GATE additions (R drift improvement vs β α=1.0
  baseline, boundary-tag stability, active-work finite), 1 GATE recomputed
  (anchor force balance with active term), 1 GATE suspended (energy-
  monotone), 1 DIAGNOSTIC (apparent contact angle).
- **Predicted Stage 1a++ pilot outcomes** (per $\zeta_\star$):
  - $\zeta_\star = 0.1$: weak active correction; R drift expected
    23–24% (small reduction from β α=1.0 baseline 24.7%).
  - $\zeta_\star = 0.3$: mid-range; R drift expected 18–22%.
  - $\zeta_\star = 1.0$: strong active correction; R drift expected
    10–18%, possibly approaching the Laplace limit ≈ 1%.
  - If all three give R drift ≈ 24% (insensitive to $\zeta_\star$):
    Layer 2 boundary-active-stress alone is also insufficient — the
    architectural ceiling extends; activate stochastic events (Stage
    1a++.b) or trigger Track B (v15 deferred diagnostics, see
    `docs/v15_deferred_diagnostics_plan.md`).

## PI decisions (resolved 2026-04-29)

1. **Stage 1a++ scope** — APPROVED. Layer 1 v15 + Stage 1a+ Option β
   α=1.0 carrier + Layer 2 continuum boundary active stress. Layer 2
   stochastic events deferred to Stage 1a++.b.
2. **Option II (3-run sweep)** — APPROVED.
3. **Three contract changes** (energy-monotone suspended, anchor-force-
   balance integrand recomputed, contact-band ρ_kernel
   reinterpretation) — APPROVED.
4. **`ζ_star` literature verification failed (PARTIAL)** — RESOLVED via
   **Option α'**: sweep ζ/K only (no Pa claim), K-anchored to
   Fischer-Friedrich Nat Cell Biol 2014 IF 30, Marchetti Rev Mod Phys
   2013 IF 50 retained as framework reference. Result reported as
   response curve in ζ/K, not as a claim about MCF7 active-stress
   magnitude. Magic-Number Block automatically PASS.
5. **Bounded outcomes for Stage 1a++ pilot** — APPROVED, drafted in
   `docs/outcomes_stage1a_plus_plus.md`.

Stage 1a++ implementation begins under these decisions. Stop conditions
remain in force: no magic numbers, no gate semantics edits beyond the
three explicit contract changes recorded here, no v13 anti-pattern,
halt and surface to PI on any FAIL.

## Track B trigger conditions (v15 deferred diagnostics)

If Stage 1a++ Option II sweep produces R drift ≥ 22% across **all
three** $\zeta_\star$ values (analogous to Stage 1a+ Option β
meta-bucket 3, "scheme fails at all α"), then Layer 2 boundary-active-
stress is *also* insufficient and the v15 single-layer-limit
interpretation is at risk of being a numerical artifact. In that case,
trigger Track B per `docs/v15_deferred_diagnostics_plan.md`. Conversely,
if Stage 1a++ produces a clear R-drift reduction that scales with
$\zeta_\star$, Track B remains a reference and is not executed.

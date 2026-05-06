# Path C — initial contact + effective gravity: pre-implementation Sanity Gate

This document is the pre-implementation sanity review for Path C, which adds
buoyancy-corrected gravity as an external body force on top of the existing
Layer 1 + substrate framework. The motivation is the inherited Stage 1a+
**lift-off pattern** observed in every substrate-anchored run (Stage 1a+
Option α/β, Stage 1a++ Layer 2, Stage 1b Layer 3): without an attractive
external force pushing the spheroid into the substrate, surface tension
contracts the spheroid AND lifts it off the rigid floor over time. The
substrate CSF impulse alone (Stage 1a+ Option β) cannot anchor this
because it acts only on the contact band; once the band depopulates,
contact is lost.

**Path C resolves this by adding gravity** — a continuous body force on every
particle, pointing toward the substrate (−ẑ). The substrate reflective BC
absorbs the gravity-induced downward momentum (Newton's-3 reaction = substrate
force on spheroid) and equilibrium is determined by Laplace pressure ↔ surface
tension ↔ gravity ↔ substrate reaction balance.

The companion documents are:
- `docs/02_force_models.md` §1.7 — Gravity & Buoyancy reference (Stewart MP
  et al. *Hydrostatic pressure and the actomyosin cortex drive mitotic cell
  rounding*, Nature 2011, 469:226, IF 65, primary anchor for cell density
  ≈ 1.05 g/cm³; ρ_medium ≈ 1.00 g/cm³ standard cell-biology textbook)
- `docs/v1/stage1a_plus_substrate_sanity.md` — substrate framing (single Col1
  spreading surface; substrate not a sweep variable)
- `docs/v1/outcomes_v15.md` / Stage 1a+ / 1a++ / 1b outcomes — inherited gates +
  contract changes
- `docs/12_validation.md` — Sanity-Gate Protocol (six checks) + Magic-Number
  Block

## Path C scope

### What is added

1. **Effective gravity body force** — per particle, an additional impulse
   in the −z direction:
   $$ \Delta v_{z,p}^{\text{grav}} = -g_\star \cdot dt $$
   applied during `_grid_op_overdamped` (after the existing CSF/drag step;
   before the reflective wall clamp). `g_star` is the dimensionless
   gravitational acceleration with buoyancy correction:
   $$ g_\star = (1 - \rho_{\text{medium}} / \rho_{\text{cell}}) \cdot g
       \cdot \tau_{\text{relax}}^2 / R_0. $$
2. **Initial position unchanged**: spheroid centre at $z^\star = R_0^\star$
   (just touching the substrate, identical to Stage 1a+). Gravity then
   continuously presses the spheroid down; the reflective BC clamps the
   bottom; the natural equilibrium is a spherical-cap shape determined by
   surface tension ↔ gravity ↔ substrate-reaction balance. **No new initial-
   condition geometry is introduced** (no spherical-cap pre-truncation; no
   shifted-centre placement). This is the simplest possible Path C entry.
3. **Energy diagnostic extension**: gravitational potential energy
   $$ U_{\text{grav}} = \sum_p \rho_p \cdot g_\star \cdot z_p \cdot V_0 $$
   added to `surface_energy_star` diagnostic (informational only; the
   energy-monotone gate semantics inherited from Stage 1a++ already
   suspended for active-stress runs and remain so).
4. **Vertical-momentum gate semantic recomputation**: under Path C the
   substrate absorbs downward gravity-induced momentum continuously. The
   existing horizontal-only momentum gate (Stage 1a+ Cousin-Rule contract
   change) stays in force unchanged; vertical drift is logged but not gated
   as before.

### Default Path C v15 baseline scope (this run only)

Per PI directive 2026-04-29 *"기본: v15만 re-baseline 후 STOP, 다른 stage 결정 대기"*:
- **One** Path-C-v15 baseline pilot only.
- Layer 1 v15 + Stage 1a+ Option α substrate (mechanical anchor only,
  γ_sub = 0) + Path C gravity. Layer 2 / 3 / 4 / 5 / 6 OFF.
- 4 sim-hr pilot, identical to v15 / Stage 1a+ pilot duration.
- Gate against v15 free-spheroid R drift baseline 0.244 AND Stage 1a+
  Option α substrate-only baseline 0.247.

Re-baselining of Stage 1a+ Option β / Stage 1a++ Layer 2 / Stage 1b Layer 3
under Path C is **deferred to PI's explicit decision** after the baseline
pilot result is surfaced.

---

## Sanity Gate (six checks)

### 1. Dimensional analysis

- Cell density (Stewart MP et al. Nature 2011 IF 65, already cited in
  `docs/02_force_models.md` §1.7): ρ_cell ≈ 1050 kg/m³.
- Medium density: ρ_medium ≈ 1000 kg/m³ (standard DMEM water-based cell
  culture medium, universal cell-biology textbook value).
- Effective acceleration with buoyancy correction:
  $g_{\text{eff}} = (1 - \rho_{\text{med}}/\rho_{\text{cell}}) \cdot g
  \approx (1 - 1000/1050) \cdot 9.81 \approx 0.467$ m/s².
- Acceleration unit (dimensionless solver):
  $a_\star = R_0 / \tau_{\text{relax}}^2 = 100\,\mu\text{m} / (60\,\text{s})^2
  = 2.78 \times 10^{-8}$ m/s².
- Literal dimensionless gravity: $g_{\text{lit}}^\star = g_{\text{eff}} /
  a_\star \approx 1.68 \times 10^{7}$. **Massively larger than any other
  force in the simulation** — would produce $v_\star^{\text{terminal}} =
  g^\star / \xi^\star = 1.68 \times 10^7$, far beyond `max_speed_over_vrms`
  gate. This corresponds to the Stokes-single-sphere prediction
  $v_t \approx 1$ mm/s for a 100 μm sphere in water, which is
  unphysically fast for a CELL aggregate where the effective drag from
  cellular viscoelasticity is enormously larger than water Stokes drag.
- Resolution: the overdamped solver uses $\xi_\star = 1$ as a calibration
  for surface-tension-driven motion (cellular Re ~ 10⁻¹³). The "literal"
  dimensionless gravity is incompatible with this calibration. **Path C
  uses `g_star` as a body-force coefficient anchored to the gravity
  framework but tuned for the relevant overdamped regime** — same pattern
  as ζ_star Option α' (sweep is the result, not a literal Pa claim).
- Provisional default: `g_star = 0.1`. This gives:
  - Sustained downward velocity: $v_z^\star \le g_\star / \xi_\star = 0.1$
    (well within numerical bounds, comparable to typical bulk velocities
    in earlier stages).
  - Body-force-per-volume comparable to surface-tension-per-area:
    $g_\star / \kappa^\star \approx 0.1 / 2.8 \approx 0.036$ → similar
    order of magnitude to $\gamma_\star = 0.01$, so gravity opposes
    surface tension's lift-off without dominating the dynamics.
  - Sedimentation timescale: $\tau_{\text{sed}}^\star \approx R_0^\star /
    v_z^\star = 1 / 0.1 = 10$ τ_relax. Spheroid settles into substrate
    within ~10 minutes of sim time, comparable to typical PI experimental
    settling times (~ few minutes pre-spreading).
- CFL: gravity adds no new fast scale (no inertia in overdamped). dt
  unchanged.

**Check 1: PASS** with explicit honest disclosure that `g_star = 0.1` is
NOT the literal gravity value (which is unphysical in the overdamped
framework calibrated to surface-tension scale); it is a body-force
coefficient anchored to gravity-framework dimensions, surface-tension
balance order-of-magnitude, and the existing $\xi_\star = 1$ calibration.

### 2. Boundary cases

- `g_star → 0`: recovers Stage 1a+ Option α (mechanical-anchor-only,
  pre-Path-C lift-off pattern). PASS by construction.
- `g_star → ∞`: gravity overwhelms surface tension; spheroid pancakes
  into substrate; numerical instability via $v_z^\star \to \infty$.
  Caught by `max_speed_over_vrms` gate. The provisional value 0.1 sits
  comfortably inside the safe interior.
- N → 0 / N → ∞: per-particle scaling; same as v15.
- Δt → 0: gravity impulse → 0; safe.
- Δt → large: same as Stage 1a+ (CFL bounds via existing gates).
- Initial position: unchanged from Stage 1a+ (centre at $z^\star = R_0^\star$);
  no new boundary case introduced.

**Check 2: PASS**.

### 3. Conservation invariants

- **Mass**: unchanged. ✓
- **Momentum (horizontal)**: unaffected by gravity (which is purely
  vertical). Stage 1a+ horizontal-only momentum gate stays in force,
  unchanged tolerance. ✓
- **Momentum (vertical)**: gravity continuously injects downward momentum;
  substrate continuously absorbs upward reaction. Net: vertical momentum
  drift is bounded by the time-integrated imbalance between gravity input
  and substrate reaction. At equilibrium, these balance perfectly
  (substrate reaction = gravity force on spheroid + Laplace pressure on
  contact band). **Vertical momentum is not gated** (logged only) per
  Stage 1a+ Cousin-Rule contract change inherited; no new contract
  change. ✓
- **Angular momentum**: gravity is irrotational → contributes zero torque
  about any horizontal axis through the spheroid's centre of mass. ✓
- **Energy**: gravity is conservative (gradient of potential
  $U_{\text{grav}} = \rho \cdot g_\star \cdot z \cdot V_0$). Total energy
  budget changes: `KE + U_strain + U_surface + U_grav` (instead of just
  `KE + U_strain` in Stage 1a / `+ U_surface` in Stage 1a++ etc.).
  The energy-monotone gate is **already SUSPENDED** under Stage 1a++ +
  Layer 3 active runs (Cousin-Rule contract change inherited). For the
  Path C v15 baseline (Layer 1 + substrate + gravity, no Layer 2/3),
  energy monotone gate is back ON; Path C contract: include $U_{\text{grav}}$
  in the energy sum:
  $$ E_{\text{total}}^{\text{Path C v15}} = KE + U_{\text{strain}} + U_{\text{surface free}} + U_{\text{grav}}, $$
  monotone-decay limit unchanged (CSF + gravity together drive the
  spheroid toward equilibrium which then dissipates via overdamped drag).

**Check 3: PASS** with the energy-monotone gate sum extended to include
$U_{\text{grav}}$.

### 4. Numerical sanity

- $dt^\star = 0.01$ unchanged.
- $g_\star = 0.1$ produces per-step vertical impulse $g_\star \cdot dt =
  10^{-3}$ velocity units. Comparable to the existing horizontal CSF
  impulses; well within f32 precision and the existing speed gates.
- f32 precision: unaffected (g_star ~ 0.1, linear addition).

**Check 4: PASS**.

### 5. Sign / sense check

- `Δv_z^{grav} = -g_star · dt < 0` (downward) for `g_star > 0`. ✓
- Gravity pulls spheroid toward substrate; substrate reflective BC
  catches it; net effect: spheroid is pinned to substrate from below
  while surface tension can still pull it laterally / upward. ✓
- For `g_star = 0`: identical to Stage 1a+ Option α (no gravity, no
  wetting-energy term, mechanical anchor only). ✓

**Check 5: PASS**.

### 6. Measurement-protocol consistency (NEW since v13)

- **R drift**: same measurement protocol (effective_radius from second
  moment of particle positions). Gravity changes the EQUILIBRIUM shape
  (spherical cap instead of free spheroid), but the measurement is
  isotropic and well-defined for any shape. ✓
- **Contact area `A_contact_xy_hull`**: under Path C the contact band
  is **continuously populated** by gravity-pinned particles (no lift-off
  expected). Measurement protocol unchanged, but **interpretation
  changes**: A/A₀ becomes a meaningful spreading metric (analog of PI's
  experimental A/A₀) rather than an artifact-prone substrate-projection
  that depopulates from above. This resolves the inherited Stage 1a+ A/A₀
  measurement issue (Bucket I gate FAIL across all 3 Stage 1b phenotype
  runs due to lift-off). Path C is therefore a **fix to the A/A₀
  measurement protocol**, not just a physics addition.
- **Anchor force balance**: substrate reaction now equals
  `bulk pressure + gravity force on spheroid`. The existing
  `(σ_vol_zz + σ_act_zz)·dA`-integral measurement misses the gravity
  term. **Update to**:
  $$ F_{\text{substrate}} = \int_{\text{contact}} (\sigma_{vv} + \sigma_{aa})_{zz}\,dA + M_{\text{spheroid}} \cdot g_\star, $$
  with $M_{\text{spheroid}} = \rho_\star \cdot V_{\text{spheroid}}^\star
  = 1 \cdot 4\pi/3 \approx 4.19$. The anchor-force-balance gate semantic
  is the same (Newton-3 check), but the integrand is updated.
- **Vertical momentum**: not gated (logged only); same as Stage 1a+.

**Check 6: PASS** with the anchor-force-balance integrand update and the
explicit recognition that Path C resolves the A/A₀ measurement-protocol
issue inherited from Stage 1a+.

---

## Magic-Number Block

### Cell density `ρ_cell` and medium density `ρ_medium`

| Quantity | Value | Anchor | IF |
|---|---|---|---|
| ρ_cell | 1050 kg/m³ | Stewart MP et al. *Hydrostatic pressure and the actomyosin cortex drive mitotic cell rounding*, Nature 2011, 469:226 (already cited in `docs/02_force_models.md` §1.7 as primary for cell density measurement) | 65 |
| ρ_medium | 1000 kg/m³ | Standard DMEM water-based cell-culture medium; universal cell-biology textbook value (~ pure-water density at 37°C is 993 kg/m³; DMEM with serum is closer to 1003 kg/m³; 1000 ± 10 kg/m³ is the established range) | textbook |
| g | 9.81 m/s² | Universal physical constant | — |

Tests 1, 2, 3 PASS for ρ_cell, ρ_medium, g. **All three values are
well-anchored** (unlike the v15 ρ_floor / γ_sub_Col1 / ζ_star Option α'
verification cases). The Stewart 2011 IF-65 reference explicitly
provides cell density in the methods/SI per the existing
`02_force_models.md` §1.7 attribution.

### `g_star` (dimensionless gravity coefficient) — PARTIAL PASS analogous to ζ_star Option α'

The "literal" dimensionless gravity from these literature values is
$g_{\text{lit}}^\star \approx 1.68 \times 10^7$, which is incompatible
with the overdamped framework's $\xi_\star = 1$ calibration (would give
$v_\star^{\text{terminal}} \sim 10^7$). Per the dimensional-analysis
discussion in check 1, Path C uses `g_star` as a **body-force coefficient**
anchored to the gravity framework but tuned for the relevant overdamped
regime:

- **Provisional default `g_star = 0.1`** — selected such that:
  (i) gravity force per spheroid volume is $O(\gamma \cdot \kappa)$ in
  magnitude (substrate-anchoring physics dominant order),
  (ii) terminal sedimentation velocity $g_\star / \xi_\star = 0.1$ is
  well within numerical bounds,
  (iii) sedimentation timescale ~ 10 τ_relax matches the experimental
  settling timescale of cell aggregates (PI's plating protocol).

1. **Derivable** — partial. The framework (gravity + buoyancy) is
   anchored to Stewart Nature 2011 + universal constants; the specific
   value 0.1 is a dimensional-coefficient choice within the overdamped
   solver's calibration regime (analogous to the ζ_star Option α'
   resolution: framework-anchored, sweep-as-result framing; specific Pa
   pinpoint not claimed). PARTIAL.
2. **Grid-invariant** — yes, `g_star` is a dimensionless coefficient;
   independent of dx, dt, n_particles. PASS.
3. **Fitting** — no. `g_star = 0.1` is set BEFORE running any pilot,
   chosen from order-of-magnitude balance with surface tension (not from
   any gate target). PASS.

**Magic-Number Block resolution** (analogous to ζ_star Option α'):
report results as a **response to g_star** rather than as a claim about
a specific dimensional gravity value. For the Path C v15 baseline pilot,
use `g_star = 0.1` as a single point; if PI later requests a sweep,
sweep `g_star ∈ {0.05, 0.1, 0.2, 0.5}` (covers an order of magnitude
around the provisional value). PI authorisation 2026-04-29 ("Full
authorization") suffices for the single-point baseline run; no further
sub-decision needed.

---

## New gate / diagnostic candidates (Path C v15 baseline)

| Gate | Status | Tolerance | Source |
|---|---|---|---|
| All Stage 1a+ Option α gates (mass / horizontal momentum / no NaN / max-speed / VRAM / v15 inherited) | inherited | unchanged | Stage 1a+ Option α |
| Energy monotone — extended sum `KE + U_strain + U_surface_free + U_grav` | NEW (extended Stage 1a baseline) | inherited tolerance | Stage 1a baseline + Path C extension |
| Anchor force balance (recomputed integrand including gravity) | recomputed | ≤ 0.20 | Path C check 6 |
| Contact-band ρ_kernel ∈ [0.85, 1.15] | inherited | unchanged | Stage 1a+ |
| **A/A₀ trajectory finite & non-pathological** (NEW MEANINGFUL GATE under Path C) | live | min A/A₀ ≥ 0.5 | Path C resolves the lift-off measurement issue |
| **R drift** (same measurement, expected dramatically improved under Path C) | live | strict-less than v15 baseline 0.244 AND Stage 1a+ Option α baseline 0.247 | inherited improvement gate |

---

## Implementation outline

- **acs/physics/mlsmpm.py**:
  - SolverConfig: add `gravity_star: float = 0.0` (default OFF).
  - `_grid_op_overdamped`: after CSF impulse, before drag, add per-grid-cell
    gravity impulse `v[2] -= dt · gravity_star` for cells with mass.
    Same compile-time pattern as substrate CSF impulse: when
    `gravity_star == 0`, the impulse is identically zero by arithmetic.
  - `_compute_invariants`: extend `diag_surface_energy` (or add new
    `diag_grav_pe` field) with $U_{\text{grav}}$ term (per-particle
    `+ ρ_p · g_star · z_p · V₀`). Sign: $U_{\text{grav}}$ is positive
    above z=0, decreases as spheroid sinks into z=0.
- **acs/runner.py**: read `gravity.gravity_star` from config; extend energy-
  monotone gate sum to include `U_grav` when `gravity_star > 0`.
- **One config**: `configs/path_c_v15_baseline.yaml`. Layer 1 v15 + Stage
  1a+ Option α (mechanical substrate, γ_sub = 0) + Path C gravity
  (`gravity_star = 0.1`). Layer 2 / 3 OFF.

---

## Decision request to PI (resolved 2026-04-29 by full authorisation)

PI granted **full authorisation** 2026-04-29 for Path C including:
- Stewart Nature 2011 IF 65 cell-density value (well-anchored, no
  verification escalation needed).
- `g_star = 0.1` provisional value (PARTIAL Magic-Number Block per the
  dimensional-coefficient framing; analogous to ζ_star Option α' precedent;
  PI authorisation suffices for the single-point baseline).
- Single Path C v15 baseline pilot (default scope per "기본: v15만
  re-baseline 후 STOP").
- Energy-monotone gate sum extension to include `U_grav`.
- Anchor-force-balance integrand update.

Stage 1a+ Option β / Stage 1a++ Layer 2 / Stage 1b Layer 3 re-baselining
under Path C is **deferred to PI's explicit decision** after the Path C
v15 baseline result is surfaced.

Stop conditions remain in force: no magic numbers, no gate semantics
edits beyond the recorded extensions, no v13 anti-pattern, halt and
surface to PI on any FAIL or critical error.

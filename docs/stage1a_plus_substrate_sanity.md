# Stage 1a+ — substrate contact: pre-implementation Sanity Gate

This document is the pre-implementation sanity review for Stage 1a+, in which
a rigid Col1-coated substrate is added under the v15 (k.3) Layer 1 spheroid.
Layer 2 (active boundary biology) remains *off*; this stage tests passive
wetting only. Implementation does not start until the PI signs off on this
proposal.

The companion documents are:
- `docs/stage1a_interior_pressure_sanity.md` — v15 (k.3) volumetric stress proposal
- `docs/outcomes_v15.md` — v15 pilot result + Path B decision (Layer 1 limit
  acknowledged; advance to Stage 1a+)
- `docs/04_simulation_setup.md` §"Substrate Model (Stage 1)" — original framing
- `docs/02_force_models.md` §1.5–1.6 — substrate adhesion energy + Hertz contact
- `docs/12_validation.md` — Sanity-Gate Protocol (six checks) + Magic-Number Block

## Stage 1a+ scope: single Col1 substrate, single spheroid initial state

Stage 1a+ tests **only** the substrate contact mechanism, not any
phenotype variation. Concretely:

- **Substrate**: Col1-coated rigid surface, single literature-anchored
  `γ_sub_Col1`. No condition sweep at this stage.
- **Spheroid initial state**: default `acs.runner` initialisation — the
  generic relaxed Stage-1a-end spheroid. No phenotype variation at this
  stage.
- **Single simulation per Stage 1a+ pilot**: does adding a Col1
  substrate (with its literature-anchored adhesion energy) materially
  reduce the v15 baseline R drift of 24.4%? The Stage 1a+ result is a
  single number, evaluated against the v15 single-number baseline.

Anything beyond this — multi-condition comparisons, phenotype-dependent
initial states, integrin/E-cadherin sweeps — is **deferred to Stage 1b
or later**. See "Future work" at the end of this document.

## Stage 1a+ scope

### What is added (v15 Layer 1 → v16 Layer 1 + Col1 substrate)

1. **Substrate plane** at `z = 0` (dimensionless: `z* = 0`). Treated as
   *rigid* per Stage 1 assumption #8 (CLAUDE.md `docs/12_validation.md`):
   TCPS / glass at GPa stiffness vs cell at kPa — substrate deformation
   negligible. Col1 coating is mechanically thin (≪ 100 nm) and does
   not soften the effective substrate stiffness for our purposes.
2. **Normal-direction contact**: no penetration through `z = 0`. Implemented
   as a grid-velocity reflective BC in a thin layer (`z < n_contact · dx`)
   in `_grid_op_overdamped`, by clamping `v_z := 0` whenever the grid cell
   is inside the contact band and `v_z < 0`. Same construction as the existing
   reflective box-wall safety net (mlsmpm.py, `_grid_op_overdamped`), so no
   new code pattern is introduced.
3. **Substrate adhesion energy `γ_sub_Col1`** (wetting force): an additional
   CSF impulse on grid cells within a one-cell substrate band, of magnitude
   `γ_sub_Col1 · κ_sub_proxy · n̂_sub`, where `n̂_sub = -ẑ` is the substrate
   inward normal and `κ_sub_proxy` is set to `1/dx*` so the impulse is
   localised to the contact band only (no curvature is well-defined on a
   flat infinite plane; the `1/dx*` choice matches the magnitude of the
   surface-tension impulse for the same `γ`, giving a dimensionless
   capillary number `Ca_sub = γ_sub_Col1 / (K · R₀)` directly comparable to
   the free-surface `Ca = γ_cc / (K · R₀)`).
4. **No gravity/buoyancy in Stage 1a+** (keep Stage 1 assumption #13 for
   now). The wetting force is the *only* downward driver. Initial particle
   positions are placed in contact with the substrate (centre at
   `z* = R₀`, so the lowermost particles sit on `z* = 0` at t=0 — no free-fall
   transient to wait through, identical equilibration approach to the
   spheroid-only Stage 1a runs).

### What stays off (Layer 2+ explicitly disabled)

- Lamellipodia / filopodia / leader cells (Layer 2)
- φ-ODE adhesion network dynamics (Layer 3)
- Internal flow / Marangoni / nematic (Layer 4) — applies to bulk
  *cellular* Marangoni (intra-spheroid γ-gradient), not a substrate term;
  also off in Stage 1a+
- Mechano-osmotic turgor (Layer 5)
- Chemistry / necrosis (Layer 6)

These remain togglable via the `layers:` block in the YAML config; Stage 1a+
sets only `L1_bulk_hydrodynamics: true`, all others `false`.

### Why now (v15 baseline)

v15 finished Layer 1 alone: shell-density witness ruled out the
onion-peeling artifact, but R drifts to 0.756·R₀ (24.4%) because Layer 1's
bulk pressure response cannot fully balance Laplace surface tension at the
predicted equilibrium radius (`docs/outcomes_v15.md` Pilot Result §). Stage
1a+ asks the next mechanistic question: **does anchoring the spheroid to
a rigid Col1-coated substrate (with literature-derived wetting energy)
materially reduce the residual R drift?** If yes, the v15 quantitative limit
is naturally extended; if not, the (i)–(iii) deferred-diagnostic candidates
from `docs/outcomes_v15.md` come back into scope.

---

## Sanity Gate (six checks)

### 1. Dimensional analysis

Solver remains in dimensionless units (length=R₀, time=τ_relax, stress=K).
New dimensionless quantity introduced:

- Substrate capillary number `Ca_sub = γ_sub_Col1 / (K · R₀)`. Mid-range
  Col1-spreading literature value (Beaune PNAS 2014 §Col1 spreading,
  Douezan PNAS 2011): with γ_sub_Col1 ≈ 0.3 mJ/m², K = 1 kPa, R₀ = 100 μm
  ⇒ `Ca_sub ≈ 0.003`. Of order the existing `Ca_cc = 0.01` for cell-cell
  adhesion.

CFL: substrate adds no new explicit-stepping bound. The reflective BC is
algebraic (set `v_z := 0`), no new wave or relaxation timescale. The
substrate CSF impulse magnitude `Ca_sub · K · R₀` is comparable to the
free-surface CSF, so `dt vs τ_γ_sub` follows the same bound that v15
already passes. No new non-dimensional number is added beyond `Ca_sub`.

**Action**: verify `Ca_sub ∈ [0.001, 0.01]` from the Col1-specific
literature range (Beaune 2014, Douezan 2011); record the SI restoration
table in `docs/stage1a_assumption_review.md` §A6.

**Check 1: PASS** (pending literature anchor for `γ_sub_Col1` —
see Magic-Number Block below).

### 2. Boundary cases

- `γ_sub_Col1 → 0`: no wetting force; only the reflective normal-velocity
  BC acts. Spheroid sits on the substrate but does not adhere — recovers a
  "spheroid touching a frictionless rigid floor" baseline. R drift should
  stay near the v15 24.4% (small geometric perturbation only).
- `γ_sub_Col1 → ∞`: wetting impulse dominates and pulls all particles into
  the contact band. Numerically this would saturate the f32 stress field.
  Gate via `max_speed_over_vrms` already in place catches the divergence.
  The Col1 literature range tops at ~0.5 mJ/m² ⇒ Ca_sub ≲ 0.005, well
  below any saturation.
- `K_substrate / K_cell → ∞` (the rigid-substrate limit we use): the
  reflective BC is exact in this limit. PASS by construction.
- `K_substrate / K_cell → 1` (substrate as compliant as cell, NOT in Stage
  1a+ scope): would require a separate substrate-particle layer. Out of
  scope; flagged for Stage 3 substrate-stiffness sweep if the PI later
  requests it (e.g. PA gel comparison per `docs/00_project_vision.md`).
- `n_contact_band → 0`: reflective BC degenerate (no cells to clamp). Set
  `n_contact_band ≥ 2` constructor invariant.
- `n_contact_band → ∞`: BC clamps the entire half-space; spheroid frozen.
  Set `n_contact_band ≤ grid_n / 8` constructor invariant. Sweet spot
  `n_contact_band = 3` (matches the reflective box-wall margin already in
  `_grid_op_overdamped`).

**Check 2: PASS** with the two `n_contact_band` constructor invariants.

### 3. Conservation invariants

- **Mass**: exact (no particle creation; same MLS-MPM construction). ✓
- **Momentum**: the reflective BC at `z = 0` is a deliberate momentum-leak
  channel — vertical momentum flowing into the substrate is absorbed.
  This is *physically correct* for a rigid floor and matches the existing
  reflective box-wall behaviour. The v15 momentum-drift gate (limit 1e-3,
  with V_FLOOR fallback for near-rest) must therefore be re-interpreted in
  Stage 1a+: the gate measures the *horizontal* momentum drift only;
  vertical drift is allowed up to the substrate-impulse magnitude. Add a
  separate `momentum_horizontal_drift_rel_max` gate; relax the existing
  full-norm momentum gate to "vertical not checked" with explicit log line.
  This is *not* a gate-tightening; it's a contract change demanded by the
  new physical scenario, surfaced to PI per the Cousin Rule.
- **Angular momentum**: not conserved in Stage 1a+ either — the substrate
  exerts a torque on any contact patch off the spheroid centre. Continue
  the v12/v15 practice of tracking it for diagnostic only; no gate.
- **Energy**: substrate adhesion energy is a *new sink* (per unit contact
  area, energy `−γ_sub_Col1 · A_contact` is stored). The energy-monotone
  gate must include this term: `E_total = KE + U_strain + γ_cc ·
  A_free_surface + (−γ_sub_Col1 · A_contact)`. Add `surface_energy_substrate`
  term to the diagnostic and to the monotone gate. Same construction as
  the existing `surface_energy_star`. ✓
- Possible new leaks: (i) particles drifting into the contact band might
  produce spurious atomic-add scatter if they cross `z = 0` between
  `_clear_grid` and `_p2g_mass`. Add a constructor invariant: any particle
  with `z < 0` after a step is hard-error (NaN-like). The reflective BC
  prevents this in steady state; the invariant catches a numerical
  blow-up. (ii) The substrate CSF impulse could theoretically add momentum
  to floating particles if the contact band is wider than the spheroid
  bottom — gate via the `n_contact_band ≤ R₀ / dx` runtime check.

**Check 3: PASS** with three additions: (i) horizontal-momentum gate
contract change, (ii) substrate-energy term in the monotone-gate sum,
(iii) `z < 0` hard-error invariant.

### 4. Numerical sanity

- `dt*` unchanged; no new fast scale.
- Grid resolution: substrate band needs ≥ 2 cells perpendicular to the
  plane. With `dx* ≈ 0.094` (default Stage 1a pilot), `n_contact_band = 3`
  gives a substrate region thickness ≈ 0.28·R₀, comparable to the
  surface-tension support spanning 3·dx in the existing CSF.
- Float precision: substrate CSF impulse magnitude is `Ca_sub·K ≈ 0.003`,
  identical order to existing `Ca_cc·K ≈ 0.01` (already validated for f32
  in v12+).
- The reflective BC is exact (algebraic); no precision issue.

**Check 4: PASS**.

### 5. Sign / sense check

- **Reflective BC**: `v_z := 0` when `z < n_contact_band·dx` and `v_z < 0`.
  Pushes particles *not* into the substrate — same as the existing box-wall
  reflection. ✓
- **Substrate CSF impulse**: `dv = +γ_sub_Col1 · κ_sub_proxy · n̂_sub · dt /
  ρ_local`, with `n̂_sub = -ẑ` (inward toward substrate from the medium
  side). For a particle in the contact band: `dv_z = -γ_sub_Col1 ·
  κ_sub_proxy · dt / ρ_local < 0` ⇒ pulled *toward the substrate* ⇒
  wetting force is attractive ⇒ correct. ✓
- **Boundary contradiction check**: a particle in the contact band is
  pulled downward by the substrate CSF; if it touches `z = 0` and tries to
  cross, the reflective BC catches it. Two impulses in opposite directions
  on the *same particle* — but they act in different orderings within the
  step (CSF impulse first in `_grid_op_overdamped`, then reflection on the
  same kernel). Net effect: particle gets stuck at `z ≈ 0` exactly, which
  is the intended *adhered* state. ✓

**Check 5: PASS**.

### 6. Measurement-protocol consistency (NEW since v13)

Three new gates are proposed for Stage 1a+; each must be analysed against
the protocol used to evaluate it. **No gate is approved if its measurement
protocol is not consistent with its analytical reference.**

**(a) Substrate contact area `A_contact`**:
- Measurement protocol: count of grid cells in the `z = 0` plane with
  `grid_m > min_cell_mass` (same threshold as the existing CSF
  localisation gate); multiply by `dx²` to get area.
- Analytical reference: for a partially-wetting droplet of volume V at
  equilibrium contact angle θ_eq (Young), the contact radius is
  `a = (3V·sin³θ / (π(2 + cos θ)(1 - cos θ)²))^(1/3)`, and
  `A_contact = π·a²`.
- Off-protocol pathway: cells with `grid_m > min_cell_mass` *near* but not
  *on* the contact plane could be miscounted if the kernel scatter spreads
  particle mass into the `z = -dx` ghost row (which doesn't exist in our
  domain — `z = 0` is the lowest grid layer). Walk through: the smallest
  contact-band particle z is bounded by the reflective BC at `z ≥ 0`. A
  particle at `z = ε` scatters mass into z-grid indices `{0, 1}` (the
  3×3×3 stencil from base = `floor(z/dx - 0.5) = -1`, but base+0 = -1 is
  out of domain and skipped; the in-domain contributions go to `{0, 1}`).
  So `grid_m[z=0]` correctly reflects the particles in `z < 1.5·dx`. ✓
- **Check 6 (a) PASS** with the analytical-protocol mapping recorded.

**(b) Anchor force balance `F_substrate vs ∫σ_vol_zz dA_contact`**:
- Measurement protocol: integrate the v15 σ_vol stress's zz-component over
  particles in the contact band; integrate the substrate CSF impulse over
  the same band. Compare magnitudes.
- Analytical reference: at static equilibrium, `F_substrate (push-down) =
  P_internal · A_contact`, where `P_internal = K · (ρ_ref/ρ_kernel - 1)`
  evaluated at contact band particles.
- Off-protocol pathway: kernel-density at boundary particles is
  systematically biased low by ~50% (AHA 2010 §3) — but the contact band
  is *not* a free surface; particles are pinned by the substrate, and the
  kernel sees the substrate-side cells via the reflective BC (which does
  *not* cancel the kernel weight; only velocity). So `ρ_kernel` at contact
  particles approaches the bulk value, *unlike* the free-surface case.
  Verify in the run: report `ρ_kernel_contact_band_mean / ρ_ref` and
  expect ≈ 1 ± 5%, distinct from the free-surface ~0.5·ρ_ref. **This is
  itself a measurement-protocol consistency witness for the substrate
  scheme.** Add as a runtime diagnostic.
- **Check 6 (b) PASS** with the kernel-density witness for contact band.

**(c) Equilibrium contact angle `θ_eq` (Young's equation)**:
- Measurement protocol: identify the spheroid surface in the contact-band
  plane; fit a circle to the (x, y) profile at the band's lower edge;
  fit a circle to the (x, z) profile at the spheroid-medium free
  surface; compute the apparent contact angle.
- Analytical reference: `cos θ_eq = (γ_sm − γ_sc) / γ_cm = S/γ_cm`, where
  `S = γ_sub_Col1` is the spreading coefficient (Douezan PNAS 2011 framing).
- Off-protocol pathway: the v15 free-surface curvature κ measurement is
  ~40% over (`docs/outcomes_v15.md` pilot result); the contact angle
  measurement involves the SAME smoothed-colour profile and the SAME
  off-peak f″/f′ asymmetry (the v13 anti-pattern). Walking through: the
  contact-angle is a *geometric* quantity — the angle between the surface
  normal and the substrate plane — not a *curvature* quantity. The
  smoothed-colour profile gives the surface as the level set `c = 0.5`;
  the angle is determined by where this level set meets `z = 0`, not by
  any second derivative. So the f″/f′ pathology that broke v13 does NOT
  apply to the contact-angle measurement. ✓
- **Check 6 (c) PASS**.

**Check 6 overall: PASS** with three measurement protocols recorded and
the contact-band kernel-density witness added as a runtime diagnostic.

---

## Magic-Number Block (mandatory, per `docs/12_validation.md`)

Each new parameter must pass all three tests (derivable / grid-invariant /
not chosen to fit a target).

### `γ_sub_Col1` — verification failed, PI selected Option α

**Status (resolved 2026-04-29)**: PI selected **Option α** below
(`γ_sub = 0`, mechanical anchor only). The verification failure record
is preserved here as the audit trail for why Option α is the chosen
path; the parameter `γ_sub_Col1` does **not** enter the Stage 1a+ code,
config, or gates. Substrate enters as a *purely mechanical* reflective
boundary condition only. Magic-Number Block automatically passes (no
parameter to verify).

The Stage 1a+ mechanism question becomes:

> *Does mechanical anchoring alone — preventing R contraction in the
> −z direction via the rigid floor — materially reduce the v15
> R drift of 24.4%?*

Possible follow-ups (deferred, not Stage 1a+):
- If R drift improves substantially (≪ 24.4%): mechanical-anchor
  contribution dominates; energetic wetting is a smaller correction;
  Option β (Ca_sub sweep) becomes a refinement, not a primary need.
- If R drift improves only partially (10–15% range): mechanical alone
  is insufficient; Option β (3-run Ca_sub sweep, anchored to Ca_cc) is
  the natural next step.
- If R drift barely changes (≳ 20%): substrate cannot anchor in this
  framework; (i)–(iii) deferred diagnostic candidates from
  `docs/outcomes_v15.md` re-enter scope; revisit the architectural
  options also flagged there.

### Verification record (preserved for audit trail)

Per the v15 ρ_floor protocol (`docs/outcomes_v15.md` "Literature
verification record §"), every numerical parameter must be checked
against the cited references *before* implementation. The proposed
`γ_sub_Col1 = 0.3 mJ/m²` was checked against the five references
listed in the previous draft of this section:

| Reference | What it actually contains | Match for `γ_sub_Col1` numeric? |
|---|---|---|
| Douezan PNAS 2011 — *Spreading dynamics and wetting transition of cellular aggregates* | Substrate = **fibronectin** (mixed FN/PEG-PLL on glass). Cells = **S-180 murine sarcoma** (E-cadherin transfected). S discussed *qualitatively only* (sign-flip wetting transition); **no numerical S value with mJ/m² units in Table 1 or anywhere in the paper.** | **NO** |
| Beaune PNAS 2014 — *How cells flow in the spreading of cellular aggregates* | Substrate = **fibronectin** on glass / PA gels. No collagen-coated substrate runs. | **NO** |
| Plotnikov Cell 2012 — *Force fluctuations within focal adhesions* | Reports per-adhesion **traction force** (pN, force/area at focal-adhesion patches), *not* spreading-coefficient energy/area. Different physical quantity. | **NO** |
| Discher Science 2005 — *Tissue cells feel and respond to substrate stiffness* | Framing only (substrate-stiffness sensing); no numerical S value for any specific coating. | **NO** |
| Engler Cell 2006 — *Matrix elasticity directs lineage specification* | Framing only (substrate-stiffness → lineage); no numerical S value for any specific coating. | **NO** |

**Result**: there is **no Col1-specific cellular-aggregate spreading-
coefficient value** in the cited literature, and the foundational
spreading-coefficient paper (Douezan 2011) does not report a numerical
S value for *any* coating — only the qualitative sign of S inferred
from the wetting transition. The figure `0.3 mJ/m²` in the previous
draft of this section was an order-of-magnitude estimate, not a
literature value.

This fails the Magic-Number Block (`docs/12_validation.md`):

1. **Derivable**: NO — no peer-reviewed source provides a numerical
   value for `γ_sub_Col1` for MCF7-like epithelial spheroids on
   Col1-coated rigid substrates. **FAIL**.
2. **Grid-invariant**: yes (a fraction of K · R₀, dimensionless `Ca_sub`).
   PASS.
3. **Fitting**: no — the value was a guess, not chosen to make any
   gate pass — but the inability to derive it makes this moot. PASS.

Test 1 fail blocks the change unless PI selects an alternative path.

### Alternatives surfaced to PI (PI selected Option α — see top of section)

**Option α — Drop the substrate adhesion-energy term entirely; test
the substrate as mechanical anchor only.**
- Set `γ_sub = 0`. Substrate is purely a reflective wall (the
  `_grid_op_overdamped` `v_z := 0` BC) with no wetting force.
- Mechanism question becomes: does *purely-mechanical* anchoring
  (preventing R contraction in the −z direction by reflective BC)
  reduce the v15 R drift?
- No literature value needed.
- Trade-off: misses the actual wetting physics; substrate becomes a
  pure geometric constraint. But it cleanly isolates the
  mechanical-anchor effect from the energetic-wetting effect, and the
  result is unambiguous (one number).

**Option β — Use the cell-cell adhesion energy as a calibration
proxy, with the result reported as a *ratio*.**
- Set `Ca_sub = α · Ca_cc` for `α ∈ {0.1, 0.3, 1.0}` (three runs).
- Report the R-drift improvement as a function of `α`, *not* against
  any specific Col1 number.
- Magic-Number Block PASS because `α` is a sweep variable, not a
  fitted parameter, and `Ca_cc` is itself anchored (Maître Science
  2012, IF 47, already in `docs/02_force_models.md` §1.1).
- Trade-off: three runs instead of one; substrate value is now
  parameter-free (the result is a *response curve*).

**Option γ — Accept "no literature anchor" and proceed with `γ_sub_Col1
= 0.3 mJ/m²` as an explicitly numerical-safety-bound estimate, the way
the v15 ρ_floor was kept after a similar failed verification.**
- The 0.3 mJ/m² value sits inside the order-of-magnitude range
  qualitatively implied by the cellular-aggregate-spreading literature
  (Douezan / Beaune / Brochard-Wyart group; FN substrate but
  comparable cell-substrate energetics).
- Magic-Number Block test 1 marked FAIL with explicit "no specific
  literature reference" docstring (analogous to the v15 ρ_floor
  resolution in `docs/outcomes_v15.md`).
- Trade-off: matches v15 precedent for honest framing, but propagates
  one weakly-anchored number into a new gate set.

**Resolved 2026-04-29: PI selected Option α** (γ_sub = 0, mechanical
anchor only). Implementation proceeds with substrate as a reflective
grid BC only; no wetting-energy term enters Stage 1a+ code. Stop
conditions remain in force.

### `n_contact_band` (substrate-band thickness in grid cells)

Set to `n_contact_band = 3` matching the existing reflective box-wall
margin in `_grid_op_overdamped`. Justification: 3-cell margin is the
minimum that fully contains the 3×3×3 P2G/G2P stencil's reach without
truncation artifacts. Already validated through v8 → v15.

1. **Derivable** — yes, from the kernel stencil width (3 cells in 3D MPM).
   PASS.
2. **Grid-invariant** — `n_contact_band = 3` is a *cell count*, not a
   physical length, so it scales with `dx` automatically. The physical
   substrate thickness `n_contact_band · dx ≈ 0.28·R₀` shrinks with
   refinement, which is the desired behaviour (substrate becomes a
   sharper boundary). PASS.
3. **Fitting** — no. PASS.

### Substrate "curvature proxy" `κ_sub_proxy = 1/dx*`

This is the only design choice that requires careful Magic-Number Block
treatment, because it directly controls the substrate-CSF impulse
magnitude. Justification: a flat infinite plane has zero curvature, so
the Brackbill `γ · κ · ∇c` formulation does not apply. The `1/dx*` choice
reproduces the *dimensional analysis* of the surface-tension impulse (the
free-surface impulse goes like `γ · (1/R) · 1` where `1/R ~ O(1/R₀) =
1/(grid_n · dx*)`; the substrate impulse goes like `γ_sub_Col1 · (1/dx*) · 1`
where `1/dx*` is the resolution scale). With `Ca_sub_Col1 = γ_sub_Col1 /
(K·R₀)`, the substrate impulse magnitude per step is `dt · K · Ca_sub_Col1
· grid_n / ρ_local`, which is the standard MPM-CSF magnitude scale.

1. **Derivable** — yes, from dimensional analysis (the only resolved
   length at the substrate is `dx`). PASS.
2. **Grid-invariant** — `1/dx*` scales with the grid; combined with
   `Ca_sub_Col1` in the impulse formula, the substrate force per unit
   area `γ_sub_Col1 / dx` is the standard wetting-force-per-area formula
   in the sharp-interface limit, independent of grid choice. PASS.
3. **Fitting** — no. The choice is dimensional, not target-driven. PASS.

**Magic-Number Block overall: PASS** for all three parameters.

---

## New gate candidates (Stage 1a+)

These are *added* to the Stage 1a v15 gate set; none of the existing v15
gates' tolerances or normalisations are modified (per the Cousin Rule). The
only existing-gate change is the *contract change* on momentum drift
(measure horizontal only; document vertical leak as substrate-absorbed) —
this is surfaced for explicit PI sign-off, not edited inline.

Under Option α (γ_sub = 0), gates that depend on the wetting-energy
balance (Young contact angle, Young-derived contact area) have no
analytical reference and are **demoted from gates to runtime
diagnostics** for Stage 1a+. They will become live gates only if a
later stage activates a non-zero γ_sub.

| New gate | Status under Option α | Tolerance | Source |
|---|---|---|---|
| Stage 1a+ R drift improvement vs Stage 1a v15 baseline | **GATE** | strictly less than 24.4% | v15 pilot finding |
| Anchor force balance: `\|F_substrate − ∫σ_vol_zz dA\|/\|F_substrate\|` | **GATE** (substrate reaction = bulk pressure transmitted to the contact band; Newton-3 check, no γ_sub dependence) | ≤ 0.20 (finite Maxwell relaxation gives the tolerance) | first principles |
| Contact-band kernel-density `ρ_kernel_contact / ρ_ref` | **GATE** (witnesses kernel sees the substrate via the reflective BC, unlike the free surface where ρ_kernel ≈ 0.5·ρ_ref; γ_sub-independent) | ∈ [0.85, 1.15] | Adami-Hu-Adams 2010 §3 |
| Substrate contact area | **DIAGNOSTIC** (no Young reference under γ_sub = 0) | log only | — |
| Apparent contact angle θ | **DIAGNOSTIC** (no Young reference under γ_sub = 0; report measurement to compare in future Option β runs) | log only | — |

The R-drift improvement gate is **not** an absolute-value gate; it is a
*relative-improvement* gate that asks the substrate to do *some* mechanical
work. If `R_drift_1aplus ≥ R_drift_v15`, the substrate is contributing
nothing (or worse, destabilising) — surface failure to PI for review.

The comparison is **simulation-baseline-vs-simulation-baseline** (single
default initial state in both runs), not a comparison against any of the
PI's experimental conditions. Phenotype-dependent comparisons are out of
scope at Stage 1a+.

---

## Implementation options (sanity-md only; PI decision on which to take)

### Option A — Substrate as grid boundary condition (recommended)

- Reflective `v_z := 0` BC in `_grid_op_overdamped` (3-cell band)
- Substrate CSF impulse added in `_grid_op_overdamped` for cells in the
  contact band
- No new particle field
- No new P2G/G2P kernel
- Code delta: ~30 lines in `_grid_op_overdamped`; one new YAML block
  `substrate:` with `gamma_sub_Col1_J_m2`, `n_contact_band`
- Wall-clock impact: < 5%
- Restricted to rigid substrate (matches Stage 1 assumption #8). ✓

### Option B — Substrate as a separate particle layer

- New particle field of "substrate particles" placed at `z* = 0`, fixed
  positions (no advection), participate in P2G mass scatter only
- Substrate-cell adhesion via per-pair attractive force (or via the
  v15 density coupling — substrate particles raise the local
  `ρ_kernel`, which would create a *pressure* gradient toward the
  substrate via the existing v15 stress)
- Allows future `K_substrate / K_cell` sweep (Stage 3 hook for PA-gel
  comparison per `docs/00_project_vision.md`)
- Code delta: ~150 lines; new field, new constructor invariants, new
  P2G mass-only contribution
- Wall-clock impact: ~30% (proportional to extra particle count)
- Conceptually richer, but introduces complexity that Stage 1a+ does
  not need and that the rigid-substrate assumption forbids.

**Recommendation**: Option A for Stage 1a+. Option B is a proper Stage 3
prerequisite (substrate stiffness sweep), where it has a justified scope.

---

## Summary

- **All six Sanity-Gate checks PASS** with the additions noted (boundary-
  case constructor invariants on `n_contact_band`, three conservation-
  invariant additions, one new contact-band kernel-density witness for
  measurement-protocol consistency).
- **Magic-Number Block FAIL on `γ_sub_Col1`**: literature verification
  found no Col1-specific cellular-aggregate spreading-coefficient value
  in the cited references (Douezan 2011 uses fibronectin/S-180 with no
  numerical S; Beaune 2014 uses fibronectin; Plotnikov 2012 reports
  traction force ≠ spreading energy; Discher 2005 / Engler 2006
  framing-only). PI decision required among Options α / β / γ before
  any implementation. PASS for `n_contact_band` (kernel-stencil-derived)
  and `κ_sub_proxy = 1/dx*` (dimensional analysis).
- **Substrate is single-condition (Col1) in Stage 1a+, mechanical-anchor
  only (γ_sub = 0 per Option α).** No wetting-energy term in code.
- **Initial state is single (default `acs.runner` initialisation).** No
  phenotype variation. Multi-condition / phenotype work is Stage 1b+.
- **Implementation**: Option A (rigid-substrate grid BC), reflective
  `v_z := 0` only (no substrate CSF impulse under Option α). ~15 lines
  of code in `_grid_op_overdamped`, < 5% wall-clock impact.
- **New gate set**: 3 gates (R drift improvement, anchor force balance,
  contact-band ρ_kernel), 2 diagnostics (contact area, apparent contact
  angle), 1 contract change on momentum drift.
- **Predicted Stage 1a+ pilot outcome**: with γ_sub = 0 the substrate is a
  *purely geometric* constraint. Expected R drift in the 15–22% range
  (mechanical anchor prevents the bottom of the spheroid from
  contracting freely, but cannot supply an adhesive driving force into
  the substrate). A drift << 15% would suggest the mechanical anchor
  alone is sufficient; a drift ≳ 22% would indicate the substrate adds
  little, and the (i)–(iii) deferred diagnostic candidates from
  `docs/outcomes_v15.md` re-enter scope.

## PI decisions (resolved 2026-04-29)

1. **Stage 1a+ scope** — APPROVED. Single Col1 substrate; single default
   initial state; Layer 2+ off.
2. **Option A** (grid-BC implementation, no extra particles) — APPROVED.
3. **`γ_sub_Col1` verification failure** — RESOLVED via **Option α**
   (γ_sub = 0, mechanical-anchor only). No literature anchor needed;
   Magic-Number Block automatically passes (no parameter to verify).
4. **3 gates + 2 diagnostics + 1 momentum-drift contract change** —
   APPROVED. Wetting-energy-dependent gates (Young contact area, Young
   contact angle) demoted to runtime diagnostics under Option α.
5. **Bounded outcomes for Stage 1a+ pilot** — to be written in
   `docs/outcomes_stage1a_plus.md` (analogous to `docs/outcomes_v15.md`)
   under the Option α mechanism question: *"Does mechanical-only
   substrate anchoring materially reduce v15 R drift = 24.4%?"*

Stage 1a+ implementation begins under these decisions. Stop conditions
remain in force: no magic numbers, no gate semantics edits, no v13
anti-pattern, halt and surface to PI on any FAIL.

---

## Option β addendum (recorded 2026-04-29 after Option α pilot)

Stage 1a+ Option α pilot result: R drift 24.7%, contact-band depopulation
49 → 4 → 0 within t* = 20·τ_relax. Mechanical-only substrate cannot hold
the contracting spheroid; surface tension lifts it off the rigid floor.
This is the unambiguous answer Option α was designed to produce, and per
`docs/outcomes_stage1a_plus.md` it triggers the Option β escalation:
**Ca_sub = α · Ca_cc sweep, anchored to the literature-anchored cell-cell
capillary number.**

### Activation of the substrate CSF impulse code path

The substrate CSF impulse `dv_z = +γ_sub · κ_sub_proxy · n̂_sub · dt /
ρ_local` (with `n̂_sub = -ẑ`, `κ_sub_proxy = 1/dx*`) was specified and
six-check-analyzed in this document at proposal time but was not
implemented under Option α (γ_sub = 0). Option β activates exactly that
code path. The six-check analysis already on file applies unchanged:

- **Check 1 (Dimensional)**: `Ca_sub_β = α · Ca_cc`. With α ∈ {0.1, 0.3,
  1.0} and `Ca_cc = 0.01`, `Ca_sub_β ∈ {0.001, 0.003, 0.01}`. All inside
  the `[0.001, 0.01]` bound from the original analysis. PASS.
- **Check 2 (Boundary cases)**: `α → 0` recovers Option α; `α → ∞` would
  saturate (already gated). The sweep range `[0.1, 1.0]` stays in the
  safe interior. PASS.
- **Check 3 (Conservation)**: substrate CSF impulse adds energy to the
  contact band; the energy-monotone gate must include the substrate
  surface-energy term `−γ_sub · A_contact`. **Implementation note**:
  add this term to `surface_energy_star` in `_compute_invariants` when
  `γ_sub_star > 0`. (Under Option α with γ_sub = 0 the term was zero
  by construction; under Option β it must be live.)
- **Check 4 (Numerical)**: substrate CSF impulse magnitude ≤ Ca_cc · K
  = 0.01 in the worst case (α = 1.0); same order as the existing
  cell-cell CSF, validated through v12+ for f32. PASS.
- **Check 5 (Sign)**: already analyzed and PASS in this document
  (impulse pulls particles toward substrate ⇒ attractive wetting). The
  Option α result *experimentally confirmed* the sign: with γ_sub = 0 the
  spheroid lifts off, so a non-zero γ_sub with this sign convention
  must, by construction, oppose lift-off. ✓
- **Check 6 (Measurement-protocol consistency)**: under Option β the
  Young contact angle gate becomes meaningful (γ_sub > 0 gives a
  finite `θ_eq`). The geometric measurement protocol already
  established for the apparent contact angle (linear fit of r(z) over
  the lowest 0.2·R₀ slab) is used unchanged; gate tolerance ±10°
  applies. The contact area gate likewise becomes meaningful (Young
  analytical reference exists for non-zero S = γ_sub).

### Magic-Number Block on the sweep design

`α ∈ {0.1, 0.3, 1.0}` is a **sweep variable**, not a fitted parameter.
The anchored quantity is `γ_cc` (Maître et al. Science 2012, IF 47,
already in `docs/02_force_models.md` §1.1 / §1.5). Reporting the R-drift
result as a *response curve* in `α` makes the substrate value
parameter-free at result level (no specific `γ_sub_Col1` is claimed).

1. **Derivable** — yes, `γ_cc` from Maître Science 2012; `α` is a
   sweep, not a value to derive. PASS.
2. **Grid-invariant** — `α · Ca_cc` inherits Ca_cc's grid-invariance.
   PASS.
3. **Fitting** — no. The three values 0.1 / 0.3 / 1.0 are
   logarithmically spaced exploration of one decade, chosen to bracket
   the qualitative regime change (weak / moderate / γ_cc-comparable
   wetting), not to make any gate pass. PASS.

### Decision request to PI for Option β

PI-pre-approved as the Outcome 4 escalation path. No further sign-off
required to begin Option β implementation. Stop conditions remain: no
magic numbers, no gate semantics edits beyond the Cousin-Rule
substrate-surface-energy addition (recorded above and in
`docs/outcomes_stage1a_plus.md`), no v13 anti-pattern.

---

## Future work (Stage 1b and beyond — out of scope for this document)

Recorded here so the framing is preserved across sessions; **none of this
enters Stage 1a+ code or gates**.

- **Stage 1b (Layer 3 active)**: enable the φ-ODE for E-cadherin ↔
  Integrin-β1 transition (`docs/03_adhesion_dynamics.md`). At this
  stage the PI's three formation-environment conditions enter the
  simulation as initial conditions on the φ field and on γ_cc, mapped
  from the formation surface (pV4D4 ± laminin presentation: `Bare`,
  `Pre`, `Lam4`) to the spheroid starting state. Substrate stays
  Col1-only. Mechanism anchor: Cho et al. 2020 (E-cadherin → Integrin-β1
  transition driven by laminin-integrin engagement during formation).
- **Stage 3 (substrate stiffness sweep, optional)**: if PI later asks
  for a PA-gel comparison (per `docs/00_project_vision.md`), Option B
  from the Implementation Options section above (substrate as a
  separate particle layer) becomes the relevant prerequisite.

These are flagged here purely so a future session does not redo the
framing analysis from scratch. The Stage 1a+ pilot itself is single-
substrate, single-initial-state, single-number.

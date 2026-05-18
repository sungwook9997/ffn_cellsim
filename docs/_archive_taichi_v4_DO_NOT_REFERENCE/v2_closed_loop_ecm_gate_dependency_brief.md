# Closed-Loop ECM Gate — 6-Item Dependency Brief

**Date**: 2026-05-04 KST
**Status**: design-discussion brief, **docs/design only — no code commitment**.
**Author**: implementation-work Claude (drafted), routed to
design-discussion for adversarial round per Codex impl `id=1221`.
**Source**: forward roadmap `docs/v2/v2_phase1_forward_roadmap.md`
lines 63–77 ("Closed-Loop ECM Gate"), plus the completed
ECM-OL preflight series (commits `29a360f`, `2c76c0c`, `7dc1767`,
`e4a71c9`, `8b2c1ab`, `22b9456`, `04ee5a7`, `9b3dac5`, `f8cdff3`)
and the 6.3b protrusion-coupled FA dynamics series (`4548ad4`,
`6d1e12b`, `dc5043a`, `9a20fc9`, `deed45a`, `50033c3`, `2a00452`).

**Supersession note (post-lock audit, 2026-05-04)**: this brief is
the pre-lock dependency analysis. The implementation contract is now
`docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`. Where this
brief says "Item 1 partial" or "Item 5 sensitivity sweep", read that
as the pre-lock sketch of what the locked plan later names **Phase B
precursor evidence** and **Phase C open-loop sweep baseline**. Do not
grep-copy those phrases as gate-item satisfaction claims.

---

## Why this brief exists

The forward roadmap mandates that closed-loop ECM activation
requires **all six items** below to pass before any FA→ECM /
ECM-bias / ECM-driven remodeling code lands. The roadmap gives the
items as one-line statements; this brief unpacks each against the
current state of the codebase and identifies what is provable now,
what scaffolding is needed, and what is design-blocked.

Locking the closed-loop API in implementation-work would be
design-by-fiat. Routing here so the design round can resolve the
"what does response / saturation mean for this ECM" decisions
before any executable closed-loop physics is written.

---

## The 6 items (verbatim, roadmap §"Closed-Loop ECM Gate")

1. Response monotonicity under prescribed traction.
2. Saturation behavior under repeated traction.
3. No response when traction is zero.
4. Bounded feedback in a single-cell loop.
5. Sensitivity sweep over grid spacing and `dt_ecm`.
6. Rule 10 unit-chain proof for traction-density storage and
   comparisons: `nN`, `um2`, `nN*s/um2`, and kPa / force-per-area
   comparisons before any ECM force or remodeling update.

Failure of any item blocks closed-loop ECM. Open-loop ECM remains
allowed (and is already in production via the four preflight
functions and the harness/runner scripts).

---

## Per-item dependency analysis

### Item 1 — Response monotonicity under prescribed traction

**Provable now (partial)**:
- `accumulate_prescribed_traction` is monotone non-decreasing **in
  the accumulator field itself** by construction:
  `accumulated_new = accumulated + traction_density · dt`, with the
  schema rejecting negative traction. Existing tests assert this for
  uniform/ramp/step-function traction patterns.

**Gap (definition-blocked)**:
- "Response" is undefined for the current `ECMSubstrateState`
  schema. The ECM stores `stiffness_kpa`, `ligand_density`,
  `fiber_density`, `orientation_tensor`, and
  `accumulated_traction_nNs_per_um2`, but it has no displacement,
  strain, or compliance field that would constitute a "response"
  under traction. The accumulator is the cumulative stimulus, not
  the response.
- Without a constitutive law (e.g., a stiffness reduction under
  repeated traction, or a density change under cumulative
  traction), Item 1 is testable only as the trivial monotonicity
  of the accumulator itself.

**Resolution path**:
- Define one or more **closed-loop response laws** as design
  decisions — e.g., `stiffness_kpa(t+1) = stiffness_kpa(t) -
  α · accumulated_traction_change`. Each law gets its own Sanity
  Gate.
- Test Item 1 = "varying prescribed traction monotonically changes
  the response field in a documented direction". The direction must
  be locked first.

**Blocker**: "response field" definition (which schema field
responds, and via what constitutive law).

### Item 2 — Saturation behavior under repeated traction

**Provable now**: nothing. The current accumulator
deliberately does **not** saturate (open-loop preflight, monotone
unbounded by construction).

**Gap**:
- Saturation is a closed-loop phenomenon. With current preflight,
  repeated traction inflates `accumulated_traction_nNs_per_um2`
  without bound (subject only to float64 overflow at `~1e308`).
- For a closed-loop response (e.g., ECM stiffness reduces under
  cumulative traction), saturation can mean either (a) the response
  field reaching a floor (e.g., `stiffness_kpa → 0`) or (b) the
  rate of change asymptoting to zero.

**Resolution path**:
- Define **what saturation means** for the chosen closed-loop
  response law. Two candidates:
  - Hard cap: `stiffness_kpa = max(stiffness_kpa - α · ΔA, 0)`.
    Saturation = stiffness floor reached.
  - Asymptotic: `stiffness_kpa(t+1) = stiffness_kpa(t) ·
    exp(-α · ΔA)`. Saturation = exponentially-fading rate.
- Test Item 2 = "repeated identical traction over N steps drives
  the response toward the saturation value with the documented
  rate". Direction + rate locked first.

**Blocker**: choice of saturation form — hard cap vs asymptotic
(or other) — couples to the response law from Item 1.

### Item 3 — No response when traction is zero

**Provable now**: yes. Trivially true under any construction:
zero traction × any dt = 0 update. Existing tests
(`test_zero_traction_scenario_no_accumulation_change`) assert this
for the open-loop accumulator.

**Gap (post-closed-loop)**: when a response law is added, the
zero-traction test must run for the response field too. If
`stiffness_kpa(t+1) = stiffness_kpa(t) - α · ΔA`, then `ΔA = 0
when traction = 0`, so `stiffness_kpa(t+1) = stiffness_kpa(t)` —
trivially satisfied. The test still must be re-asserted on the
response field.

**Resolution path**: add a zero-traction regression test on each
closed-loop response field at the time the response law lands. No
design blocker.

### Item 4 — Bounded feedback in a single-cell loop

**Provable now**: nothing. There is no closed-loop infrastructure.

**Gap**:
- Requires a **single-cell loop** = `cell → FA → traction → ECM
  response → (optionally) ECM-biased FA placement → next step`.
- 6.3b is a one-way wrapper protrusion → FA. The lock §0 forbids
  FA → ECM, FA → protrusion writeback, ECM bias on FA / contour.
- Without these forbidden writes, no loop exists.

**Resolution path**:
- Phase D below — minimum closed-loop scaffolding (interfaces only,
  default off) — wires the loop without activating it.
- Once activated under a chosen constitutive law, test Item 4 =
  "single-cell loop runs N steps and the chosen Lyapunov-like
  function (e.g., total ECM accumulated traction, or maximum FA
  bound_fraction) stays bounded above by an analytic envelope".

**Blocker**: closed-loop scaffolding does not exist; must be
designed before any test can be written.

### Item 5 — Sensitivity sweep over grid spacing and `dt_ecm`

**Provable now (partial)**:
- All four ECM-OL preflight functions accept `(nx, ny)` grids and
  `dt_s` as caller inputs; they make no assumption beyond grid
  shape consistency and dt non-negativity.
- No sensitivity sweep test currently exists — only point tests at
  one grid + one dt.

**Gap**:
- A formal sweep over `(nx, ny) × dt_s` for each preflight function,
  asserting that the result is grid- and dt-stable in the
  appropriate sense (e.g., point-wise convergence as `dt → 0`).
- For closed-loop, the same sweep must be re-run after the response
  law is added, asserting that the response field's evolution is
  also grid/dt-stable.

**Resolution path**:
- Open-loop sweep is implementable now as an extension of
  `scripts/run_ecm_ol_harness.py` (a sweep variant). No design
  blocker.
- Closed-loop sweep blocked until response law and scaffolding land.

### Item 6 — Rule 10 unit-chain proof

**Provable now (partial)**:
- The ECM-OL module docstring (`acs/v2/dynamics/ecm_open_loop.py`,
  post `22b9456`) records the unit chain for each of the four
  preflight functions:
  - traction: `[nN/μm²] · [s] = [nN·s/μm²]`
  - stiffness rate: `[kPa/s] · [s] = [kPa]`
  - density rate: `[1/s] · [s] = [dimensionless density]`
  - orientation rate: `[1/s] · [s] = [dimensionless orientation
    tensor]`
- Hard Rule 10 is satisfied at the **storage** layer for the
  open-loop preflight side.

**Gap (closed-loop)**:
- The rule explicitly mentions "kPa / force-per-area comparisons
  before any ECM force or remodeling update". The kPa↔nN/μm²
  comparison only matters when:
  - cell-derived FA traction (in nN per FA, not nN/μm² per ECM
    cell) feeds ECM stiffness (in kPa = nN/μm²) via a constitutive
    law,
  - or when ECM stiffness biases FA force scale.
- Neither comparison exists yet; both are required for closed-loop
  activation.

**Resolution path**:
- Each closed-loop coupling function must extend its Sanity Gate's
  §1 to spell out the FA-side ↔ ECM-side unit conversion (e.g.,
  per-FA `nN` → per-ECM-cell `nN/μm²` via division by
  `dx · dy = spacing_um · spacing_um`, or via scattering into a
  finite radius).
- This is a per-function commitment, not a one-time proof.

**Status now**: Open-loop side ✓. Closed-loop side gap pending
constitutive-law lock.

---

## The chicken-and-egg problem

Several gate items (2, 4, partially 1 and 5) require closed-loop
infrastructure to be testable. But the gate forbids any closed-loop
code from landing until the items pass. Without explicit
resolution, the closed-loop layer cannot start.

**Resolution proposal**: phased entry, with each phase locked
behind a Sanity Gate of its own.

### Phase A — provable now under existing ECM-OL preflight

Items provable under the current ECM-OL only:
- Item 3 (no response on zero traction): existing test ✓.
- Item 6 (open-loop side): existing module docstring ✓.

Action: write a small `docs/v2/v2_closed_loop_ecm_gate_phase_a.md`
that records the existing tests/proofs as open-loop-side evidence.
Single small commit. No design blocker.

### Phase B — Item 1 precursor (prescribed-traction stimulus monotonicity)

Item 1 precursor proof: vary prescribed traction patterns and assert
that the stimulus accumulator (monotonicity of cumulative storage)
follows. Already covered by ECM-OL tests but should be lifted into
a dedicated "precursor, NOT Item 1 satisfaction" regression suite.

Action: small test additions to `tests/v2/test_v2_ecm_open_loop.py`
under an Item 1 precursor section. No design blocker.

### Phase C — open-loop sweep baseline (NOT Item 5 satisfaction)

Implement a sweep harness extending `acs.v2.ecm_open_loop_harness`
to run each preflight function across a grid sweep `(nx, ny) ×
dt_s`. Persist artifacts under `runs/<UTC>_ecm_ol_sensitivity/`.
New script, follows `scripts/run_ecm_ol_harness.py` pattern. This
is baseline evidence on open-loop preflight outputs, not the
closed-loop Item 5 sensitivity claim.

Action: medium unit, plumbing only, no constitutive-law decisions.

### Phase D — minimum closed-loop scaffolding (gated off)

Add the closed-loop scaffolding without activating it:
- A `step_fa_to_ecm_response` function that maps per-FA traction
  to ECM constitutive update. Default implementation = identity
  (no change).
- A `step_ecm_to_fa_bias` function that reads ECM fields and
  returns optional FA-side biases. Default = neutral (no bias).
- An integrator wrapper that runs `cell step → 6.3b FA step →
  step_fa_to_ecm_response → step_ecm_to_fa_bias → next step`.
  Each of the two new functions is opt-in.

Action: each of the two new functions and the integrator is its
own Sanity Gate + design-discussion lock. **Heavy** unit. No
constitutive law is committed yet — the default implementations
are no-ops, so the loop runs but the cell behavior is unchanged
from the open-loop case.

### Phase E — Items 2 + 4 (full closed-loop tests)

Once the constitutive law is design-locked (Item 1's response
direction + Item 2's saturation form), implement the response
function per the lock and run:
- Item 2 saturation regression test
- Item 4 bounded feedback regression test (Lyapunov-like envelope)
- Item 5 sensitivity sweep extended to the response field
- Item 6 closed-loop unit-chain proof per coupling function

Action: each test maps to a small unit. Design lock for the
constitutive law happens before this phase.

---

## Hard blockers (must be resolved in design-discussion)

1. **Constitutive-law direction (Item 1)**: which schema field
   responds to traction (stiffness? fiber_density? orientation?
   ligand_density? more than one?), and in what direction (decrease
   under cumulative traction? increase? both depending on regime?).
2. **Saturation form (Item 2)**: hard cap vs asymptotic decay vs
   other.
3. **FA → ECM scattering geometry**: a single FA at position
   `(x, y)` deposits traction onto which ECM grid cells? (Nearest
   cell? Bilinear weights? Gaussian kernel of radius `r`?)
4. **ECM → FA bias direction**: does ECM stiffness scale FA
   `traction_scale_nN`? Affect rate multipliers? Affect FA
   nucleation? (6.3b explicitly forbids `force_candidate_nN`
   scaling and ECM bias as part of its lock; that boundary must be
   explicitly broken — and gated — for closed-loop activation.)
5. **Lyapunov-like function for Item 4**: what bounded-feedback
   metric is checked? (Total accumulated traction? Maximum FA
   bound_fraction? Total stiffness change?)

Each hard blocker is a design decision that needs an adversarial
round in design-discussion before a Sanity Gate can be written.

---

## Open-loop side: nothing breaks if Phases A–C are done in
implementation-work

Phases A, B, C only extend the open-loop preflight side (ECM-OL
+ ECM-OL harness + sweep). They do not introduce closed-loop
behavior, do not need new constitutive laws, and do not break any
6.3b lock §0 forbidden item. They can run in parallel with the
design-discussion round on the hard blockers above.

**Recommendation**: implementation-work proceeds with Phase A
right after this brief lands (small docs commit), then Phase B
and C as small units (test additions + sweep harness). Phase D
waits for design-discussion to resolve hard blockers 1, 2, 3, 4
above.

---

## Sanity Gate proofs still needed (per phase, not now)

| Phase | New Sanity Gate doc(s) needed |
|---|---|
| A | None (existing docs cover it; just an evidence note) |
| B | Possibly extend ECM-OL Sanity Gate's §3 with explicit "Item 1 precursor, NOT satisfaction" section |
| C | New `docs/v2/v2_ecm_ol_sweep_sanity_gate.md` for the sweep harness |
| D | New `docs/v2/v2_fa_to_ecm_response_sanity_gate.md`, `docs/v2/v2_ecm_to_fa_bias_sanity_gate.md`, `docs/v2/v2_closed_loop_integrator_sanity_gate.md` |
| E | Per-test gates for Items 2 + 4 + 5 closed-loop side + 6 closed-loop side |

---

## Process expectation

1. **This round** (design-discussion): adversarial review of phased
   order, hard blocker enumeration, and the per-phase scaffolding
   shape.
2. **Lock artifact**: `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`
   mirroring `docs/v2/v2_63b_protrusion_coupling_locked.md`.
3. **Cross-room dispatch back to implementation-work** with the
   locked Phase A entry definition (smallest, most immediately
   actionable).
4. **Implementation entry**: Phase A first (likely a few small
   commits under existing tests + a Phase A evidence note).
5. **Design-discussion runs in parallel** for hard blockers 1–5
   so Phase D can start once they're resolved.

implementation-work stays idle/review-capable while this design
round runs. Per Codex impl `id=1221`, no closed-loop code lands
until the gate items are satisfied per phase.

---

## References

- Forward roadmap: `docs/v2/v2_phase1_forward_roadmap.md` lines 63–77
- ECM-OL preflight series: `acs/v2/dynamics/ecm_open_loop.py`
  (commits `29a360f` → `22b9456`)
- ECM-OL harness + run script: `acs/v2/ecm_open_loop_harness.py`,
  `scripts/run_ecm_ol_harness.py` (commits `2c76c0c`, `04ee5a7`,
  `9b3dac5`, `f8cdff3`)
- 6.3b dynamics + run script: `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py`,
  `scripts/run_protrusion_coupled_fa_harness.py`,
  `docs/v2/v2_63b_protrusion_coupling_locked.md`,
  `docs/v2/v2_63b_protrusion_coupled_fa_sanity_gate.md`
  (commits `4548ad4` → `2a00452`)
- Hard Rule 10 (CLAUDE.md): dimensional comparison verification,
  per-area vs per-volume canonical forms.
- Adversarial design-discussion posture: memory
  `feedback_aggressive_design_debate.md`, PI `id=809`.

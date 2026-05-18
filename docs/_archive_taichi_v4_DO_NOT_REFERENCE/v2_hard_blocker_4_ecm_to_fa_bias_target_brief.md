# Hard Blocker #4 — ECM→FA Bias Target — Design-Discussion Brief (opening, NOT a lock)

**Date**: 2026-05-04 KST
**Status**: design-discussion brief, **opening position only — NOT a
lock, NOT a code commitment**.
**Author**: implementation-work Claude, drafted under (B-trigger)
precedent established by Hard Blocker #3 (`0dadfdc`/`f86613a`
brief → design-discussion 4-round adversarial → `3ba37fa` lock →
impl `b90b07f` Sanity Gate → `883efc8` code). Codex impl
`id=1349`/`id=1351`/`id=1361`/`id=1363` confirmed the (B-trigger)
pattern is acceptable for #4 with four explicit guardrails (see
§0 below).
**Source**: closed-loop ECM gate phased plan
`docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md` §2 (Hard
Blocker #4 enumeration), §3 (effective_stiffness helper guard —
side-door risk), §6 (Sanity Gate matrix Phase D row), §9 (Phase D
blocker resolution checklist).

**Hard contract**: this brief is the **opening position** for the
adversarial design round on Hard Blocker #4. It does not lock the
interface, does not commit the implementation, and does not
authorize any Phase D code. The locked phased plan §1 Phase D
remains blocked until the design-discussion round closes with a
lock artifact (target:
`docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`).

---

## 0. Why this brief exists + Codex guardrails

The locked phased plan §2 enumerates Hard Blocker #4 as the
"ECM→FA bias target" interface decision: how does ECM state bias
FA dynamics? The plan §1 Phase D entry is gated on this blocker
(plus #3, now locked as `3ba37fa`) being interface-locked — even
if the default Phase D implementation is identity / neutral.

Without this round, Phase D no-op scaffolding cannot complete,
which in turn blocks Phase E active closed-loop response.

**Per Codex impl `id=1363`, four guardrails this brief honors and
the adversarial round must enforce**:

1. **No-op default**: per-FA multipliers / target values default
   to **exactly 1.0** (or equivalent neutral identity). No state
   mutation. No traction scalarization side-door.
2. **Typed per-FA interface first**: define the output as a
   typed per-FA interface **before** discussing any active
   constitutive law. Phase E biology is explicitly out of #4
   scope.
3. **effective_stiffness side-door guard** (locked phased plan
   §3): any helper consuming `fiber_density` /
   `orientation_tensor` is ECM→FA bias by another name. The #4
   interface must make this explicit and route any such helper
   through the locked Phase D contract — never via a "read-only
   helper" loophole.
4. **6.3b `ProtrusionStateMultipliers` style** (where applicable):
   reuse the typed-key + per-state multiplier dataclass pattern
   for the bias structure, with units/sense checks for every
   multiplier and target quantity.

---

## 1. What is already locked (do not redebate here)

Per `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`:

- §1 Phase D entry: gated on Hard Blocker #3 + #4 interface
  locks. With #3 locked at `3ba37fa` (FA→ECM bilinear scattering)
  + Sanity-Gated at `b90b07f` + implemented at `883efc8`, only #4
  remains.
- §2 Hard Blocker #4: phrased as a question ("ECM state biases FA
  dynamics how?") with four candidate mechanisms enumerated:
  - (a) `traction_scale_nN` scale (FA traction magnitude depends
    on local stiffness)
  - (b) Rate multiplier (FA `k_maturity_per_s` etc. boosted by
    local ECM state) — **6.3b `ProtrusionStateMultipliers` style**
  - (c) FA nucleation rate (more FAs in stiff regions)
  - (d) Combination
  Phase D entry contract (locked):
  - Inputs: `(adhesions, ecm_state, bias_params)`
  - Outputs: per-FA bias multipliers (typed dict, similar to 6.3b
    `ProtrusionStateMultipliers` pattern)
  - Default no-op: returns 1.0 multiplier for all FAs (neutral)
- §3 effective_stiffness helper guard: any helper consuming ECM
  geometry fields is ECM→FA bias by another name. Phase D no-op
  can keep an effective_stiffness helper IFF it returns raw
  `stiffness_kpa` exactly (identity); the moment it consumes
  `fiber_density` / `orientation_tensor`, it crosses into Hard
  Blocker #4 territory and must satisfy this round's lock.

This brief is the next-level-down design choice **inside** that
already-locked frame.

---

## 2. The four candidate mechanisms (pros/cons)

### Candidate (a) — `traction_scale_nN` scale

ECM state biases the **magnitude** of per-FA traction generation
via a multiplier on `traction_scale_nN`. Concrete example:
`traction_scale_nN_eff = traction_scale_nN · f(local_stiffness)`,
where `f` is a caller-supplied function (e.g., monotone-increasing
in stiffness for cells that pull harder on stiff substrates,
literature-documented).

**Pros**:
- Direct biological interpretation (cells generate more traction
  on stiff substrates, Lo et al. 2000 / Discher 2005).
- Single multiplier per FA per step; no extra schema field.
- Phase E response laws can compose with this without new types.

**Cons**:
- Couples the bias directly to 6.3a's `traction_scale_nN`,
  changing the magnitude formula
  `traction_scale_nN · maturity · bound_fraction`. Any future
  6.3a edit would have to consider this coupling.
- Biological domain question: stiffness as input vs output of FA
  generation (Trichet 2012 PNAS treats it as input). Locking this
  here may pre-decide a Phase E debate.
- Requires sampling `local_stiffness` at each FA position from
  the ECM grid — that needs a pointwise sampler, which is itself
  a small interface decision (nearest-cell? bilinear? Same kernel
  as Hard Blocker #3's scatter?).

### Candidate (b) — Rate multiplier (6.3b style)

ECM state biases the **rates** of FA dynamics
(`k_maturity_per_s`, `k_bind_per_s`, `k_unbind_per_s`) via a
multiplier table parallel to 6.3b's
`ProtrusionStateMultipliers`. Concrete shape:
`k_maturity_per_s_eff = k_maturity_per_s · m_kmaturity(local_ecm_state)`.

The multiplier table is typed-key (mirroring 6.3b):
- Allowed ECM "state" keys: TBD per the round (perhaps a
  dimensionless ECM stiffness band: `low / medium / high`, or a
  per-rate dimensionless coefficient parameterized on continuous
  fields).
- Allowed rate keys: exactly `k_maturity_per_s`, `k_bind_per_s`,
  `k_unbind_per_s` (matching 6.3b literal).
- Multiplier value: finite, non-bool, ≥ 0. Default 1.0 (neutral).
- Validation upfront via a `validate()` method, raising
  `multiplier_table_unknown_key` for unknown spelling.

**Pros**:
- Clean reuse of the 6.3b `ProtrusionStateMultipliers` pattern.
  Same typed-key + neutral-default + validate-upfront pattern;
  Codex review id=1129 confirmed this prevents
  silent-typo-neutral bugs.
- 1.0 multiplier produces base 6.3a behavior (no-op default
  identity).
- Multiplier 0 ⇒ effective rate 0 ⇒ no rate-driven update — but
  per 6.3b precedent, this **does NOT** trigger a state-label
  change, preserving 6.3a's no-auto-state-transition contract.
- Single unit chain `[1/s] · [dimensionless] = [1/s]` (matches
  6.3b chain).
- No coupling to `traction_scale_nN` (preserves §3
  effective_stiffness guard's "raw `stiffness_kpa` unchanged in
  first closed-loop law" position).

**Cons**:
- Requires a way to map continuous ECM state (stiffness, fiber
  density, orientation) to discrete "ECM state" keys for the
  multiplier table — or a continuous functional form that the
  validator must sanity-check.
- The 6.3b precedent uses `ProtrusionState` literals as keys;
  ECM has no analogous literal enum. The round must decide
  whether to discretize (low/medium/high stiffness) or use
  continuous parameters (per-FA scalar multiplier from a
  caller-supplied function).

### Candidate (c) — FA nucleation rate

ECM state biases **how often new FAs form** in a region (e.g.,
more FAs nucleate in stiff regions). Concrete shape: an FA
generation rate field on the ECM grid, sampled per-step to spawn
new FAs at random positions weighted by local rate.

**Pros**:
- Biologically meaningful for spreading dynamics: stiff substrate
  → more nascent FAs → more spreading.
- Decoupled from FA dynamics — works on the FA *list* rather
  than per-FA properties.

**Cons**:
- **Requires RNG** + reproducibility contract + Sanity Gate §3
  reproducibility burden. The closed-loop ECM gate phased plan
  §3 / §4 explicitly defers stochastic samplers to a separate
  Sanity Gate.
- Phase D no-op default would be "spawn no new FAs" — fine but
  trivial; the active law is what makes this candidate useful.
- Crosses Phase D / Phase E boundary because spawning FAs is
  active behavior, not interface plumbing.

### Candidate (d) — Combination

Some combination of (a), (b), and/or (c). The locked phased plan
§2 lists this as a candidate; in practice this is the "round
chooses two of the three" option.

**Pros**:
- Biologically richer than any single candidate alone.

**Cons**:
- Sanity Gate burden multiplies: each combined mechanism adds
  its own §1 unit chain, §2 boundary cases, §3 conservation, §5
  sign convention.
- For Phase D no-op, simplicity is the goal; combinations are
  better unlocked one-at-a-time after the simplest base-case is
  proven.

---

## 3. Recommended opening position (NOT a lock)

**Candidate (b) Rate multiplier in 6.3b `ProtrusionStateMultipliers`
style** as the recommended opening position.

Reasoning (per Codex impl `id=1363` guardrails):

1. **Typed per-FA interface first** ✓: 6.3b's
   `ProtrusionStateMultipliers` pattern is already a typed-key
   per-FA multiplier dataclass with validate-upfront discipline.
   Reuse delivers Codex's guardrail #2 without re-litigating the
   typed-key shape.
2. **No-op default 1.0** ✓: 6.3b precedent uses 1.0 = neutral =
   base 6.3a behavior. Phase D no-op trivially satisfies this by
   passing an empty multiplier table (every lookup defaults to
   1.0).
3. **effective_stiffness side-door guard** ✓: rate multipliers
   couple to FA *rates*, not to FA `traction_scale_nN`. The
   locked phased plan §3 explicitly preserves raw `stiffness_kpa`
   unchanged in first closed-loop law; rate multipliers can read
   `local_stiffness` to compute the multiplier without mutating
   any field. The round must lock the **read modality** (e.g.,
   nearest-cell stiffness sample, or per-FA-position bilinear
   interpolation matching Hard Blocker #3 scatter).
4. **Magic-Number Block discipline**: zero new caller-supplied
   numerical parameters beyond what the multiplier table values
   already are (caller-supplied per memo's no-default rule).
   Compare to (a) which adds a `f(local_stiffness)` function +
   parameters and to (c) which adds RNG seeding + nucleation
   rate field.
5. **Sanity Gate reuse**: 6.3b Sanity Gate items §1–§6 transfer
   almost verbatim with the substitution
   `ProtrusionState → ECM state read`. The round saves writing
   redundant gate scaffolding.
6. **Phase D no-op + Phase E unlock**: a 1.0-multiplier-everywhere
   default is exactly the no-op Phase D requires. Phase E can
   then activate non-1.0 values per a Phase E response law
   without re-litigating the Phase D interface.

**Caveat / unresolved-by-this-brief**:

- **What is "local ECM state"?** The round must lock the read
  modality: do we sample `stiffness_kpa` at the FA position
  (single scalar input), the full `(stiffness_kpa,
  fiber_density, orientation_tensor)` vector, or a derived
  scalar via a documented helper? Each choice has a different
  Sanity Gate §6 measurement-protocol entry. Per Codex
  guardrail #3, the helper choice must NOT bleed back into the
  scatter (#3 was deliberately one-way FA → ECM, opaque to ECM
  state).
- **Continuous vs discrete ECM state keys**: 6.3b uses
  `ProtrusionState` literals (`growing/stalled/retracting/ended`).
  ECM has no analogous enum. The round picks one of:
  - (b.i) Discretize stiffness into bands (`low/medium/high`)
    with caller-supplied band edges.
  - (b.ii) Use a continuous per-FA scalar multiplier
    `m_per_fa(local_stiffness)` with a caller-supplied function;
    each function call is a no-op default of 1.0.
  - (b.iii) Use a per-rate continuous multiplier function with
    documented monotonicity.

**Caveat (Codex round 2 challenges welcome)**: this opening is
gate-discipline-driven, not biology-driven. If the round wants to
argue (a) is biologically more accurate (cells *generate* more
traction on stiff substrates per Lo 2000), that is a legitimate
round 2 challenge; the response is the §3 effective_stiffness
guard and the Phase E response-law boundary.

---

## 4. Phase D no-op interface contract (restated, locked elsewhere)

Per `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md` §2 Phase
D entry contract (NOT changeable in this round):

```
def step_ecm_to_fa_bias(
    adhesions: tuple[FocalAdhesionState, ...],
    ecm_state: ECMSubstrateState,
    bias_params: ECMToFABiasParams,  # round picks the dataclass shape
) -> ECMToFABiasResult:                # round picks the result shape
    """Phase D no-op default: returns 1.0 multiplier for all FAs
    (neutral identity). Active bias under (b) rate multiplier
    style is unlocked separately (Phase E)."""
```

The round defines:
- `ECMToFABiasParams` dataclass (caller-supplied multiplier
  table or function)
- `ECMToFABiasResult` (per-FA per-rate-name multiplier list, or
  a callable that takes an FA + ECM and returns a rate multiplier)
- Local-ECM-state read modality (which kernel maps FA position
  to ECM cell sample)
- Default no-op behavior (1.0 multiplier for every FA, every rate)

---

## 5. Sanity Gate items 6.4-CL-4 will need (preview, NOT commitment)

Per Phase D row of locked phased plan §6 Sanity Gate matrix
applied to the chosen mechanism:

- §1 Dimensional: rate multiplier is dimensionless; effective
  rate `[1/s]·[dimensionless] = [1/s]` (exactly 6.3b chain).
  If continuous form is locked, the multiplier function's
  domain (input units) and codomain (output dimensionless) must
  be inline-derived per memory `rule10_unit_derivation_in_docs`.
- §2 Boundary: typed-key validation upfront (per 6.3b precedent;
  unknown-key raises). FA empty list → empty result. ECM read
  modality boundary (FA at grid edge → which cell sampled?).
- §3 Conservation: rate multipliers do not have a sum invariant
  (multiplier table just bounds individual rate boosts). What
  IS conserved: 6.3b's per-FA Newton-3 and traction algebra
  (since rate multipliers don't touch traction). Pure-function
  property: ECM not mutated; FA list not mutated; multiplier
  table not mutated.
- §4 Numerical: `dt_fa_s · max(effective_rate over all FAs) ≤
  _DT_RATE_SAFETY_MARGIN` (reused from 6.3a / 6.3b). The
  multiplier inflates the effective rate ceiling — same gate as
  6.3b.
- §5 Sign: multiplier ≥ 0 required (negative rate multiplier
  flips sign convention). Multiplier 1.0 = neutral.
- §6 Measurement-protocol: Phase D output is per-FA bias
  multipliers; Phase E response-law decides what to do with
  them. Hard Rule 11 protected by a runtime meta-test
  `test_ecm_to_fa_bias_does_not_satisfy_closed_loop_response_law`
  parallel to Phase B test 4 / Phase C test 5 / Hard Blocker
  #3 test 18.

Magic-Number Block: the multiplier table values are caller-
supplied test inputs (no project default per memo's no-default
rule); no new tunable lands in code if the round picks (b.ii) or
(b.iii) continuous form (the function is also caller-supplied).

---

## 6. What this brief is *not*

- Not a lock. The mechanism choice + interface shape are decided
  by the design-discussion adversarial round.
- Not a Sanity Gate document. The Sanity Gate is written by the
  future Phase D 6.4-CL-4 commit, after the lock artifact lands.
- Not a code commitment. No `step_ecm_to_fa_bias` function is
  proposed; the locked plan §2 already specifies the function
  signature shape, and this round only fills in `bias_params`,
  `ECMToFABiasResult`, and the read modality.
- Not Phase E response-law commitment. Phase D no-op default is
  1.0 multiplier; activating non-1.0 values is Phase E (separate
  Sanity Gate, separate adversarial round).
- Not a side-door for effective_stiffness helpers. The locked
  phased plan §3 explicitly bans read-only ECM-geometry helpers
  from becoming ECM→FA bias. This brief restates that ban as
  Phase D non-goal.
- Not authoritative for biological parameters. Any continuous
  multiplier function or band edges require literature-first
  extraction per Magic-Number Block.

---

## 7. Process expectation (mirrors Hard Blocker #3 precedent)

1. **Cross-room dispatch** to design-discussion with this brief
   + the 4 candidate mechanisms + the recommended (b) opening
   position. Status: `decision-needed`.
2. **Adversarial round** in design-discussion (Codex/Claude
   adversarial design debate per memory
   `feedback_aggressive_design_debate.md`): challenges, reasoned
   acceptances, candidate compromises. Hard Blocker #3 precedent
   showed 4 rounds is a reasonable upper bound.
3. **Unresolved disagreements list** before lock — particularly
   for the read-modality choice (b.i / b.ii / b.iii) and the
   discrete-vs-continuous ECM state representation.
4. **Lock artifact** delivered to implementation-work as
   `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`,
   parallel to
   `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`.
5. **Cross-room dispatch back to implementation-work** with the
   locked mechanism + bias_params dataclass shape. Phase D
   no-op scaffolding entry permitted (with both #3 + #4 locked).
6. **Phase D no-op scaffolding implementation** — this is the
   final Phase D entry unit, combining the #3 scatter primitive
   and the #4 bias function into a default-off integrator
   wrapper.

implementation-work stays idle/review-capable while this design
round runs. Per Codex impl `id=1363` ack, this brief is the
impl-work side's "opening dispatch" matching the Hard Blocker
#3 (B-trigger) pattern.

---

## 8. References

- Locked phased plan:
  `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- Hard Blocker #3 (FA→ECM scattering, locked + implemented):
  `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md` +
  `docs/v2/v2_fa_to_ecm_scattering_sanity_gate.md` +
  `acs/v2/dynamics/fa_to_ecm_scattering.py` +
  `tests/v2/test_v2_fa_to_ecm_scattering.py`
- Hard Blocker #3 brief precedent (this brief follows the same
  structure):
  `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_brief.md`
- 6.3b lock (typed-key multiplier dataclass precedent):
  `docs/v2/v2_63b_protrusion_coupling_locked.md` +
  `docs/v2/v2_63b_protrusion_coupled_fa_sanity_gate.md` +
  `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py`
- 6.3a (FA dynamics that the rate multiplier wraps):
  `acs/v2/dynamics/focal_adhesion.py` +
  `docs/v2/v2_focal_adhesion_dynamics_sanity_gate.md`
- ECM substrate schema (read modality target):
  `acs/v2/ecm_substrate.py`
- Forward roadmap: `docs/v2/v2_phase1_forward_roadmap.md`
  ("Closed-Loop ECM Gate" section)
- Memory rules:
  - `rule10_unit_derivation_in_docs.md`
  - `hard_rule_11_wording_boundary_meta_test.md`
  - `design_note_pre_commit_batch.md` (5-step pre-commit batch)
  - `feedback_aggressive_design_debate.md`

---

## 9. Cross-room dispatch instruction (impl-work → design-discussion)

The MCP message that accompanies this brief:
- `to`: `codex` (design-discussion-pane Codex receives + Claude
  pane reads)
- `room`: `design-discussion`
- `topic`: `v2-layer-2-hard-blocker-4-ecm-to-fa-bias-target`
- `status`: `decision-needed`
- `body`: brief 4-question summary + opening recommendation +
  reference to this file
- `refs`: implementation-work `id=1361, 1362, 1363` + locked
  phased plan reference + Hard Blocker #3 lock as precedent

design-discussion round produces a lock artifact that supersedes
this brief; this file remains as historical opening-position
context (mirror of
`docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_brief.md` post-lock).

# V2 Phase 1 — 6.3b Protrusion-Coupled FA Dynamics Sanity Gate

**Status**: pre-execution Sanity Gate per Plan §10/§11 + forward
roadmap "Separated dynamics modes". Written from the
design-discussion lock in
`docs/v2/v2_63b_protrusion_coupling_locked.md` (Claude+Codex 5-round
adversarial lock, MCP `id=1186-1198`).

**Contract**: this gate is the document Plan §11 requires **before**
any executable protrusion-coupled FA physics lands in
`acs/v2/dynamics/protrusion_coupled_focal_adhesion.py`. If any item
§1–§6 below fails its check, 6.3b reports `status=blocker` to PI
with at least three concrete options. No partial physics commits
under uncertainty.

**Source of truth**: design-discussion `topic=v2-layer-2-63b-protrusion-fa-coupling`,
MCP `id=1186, 1190, 1191, 1193, 1195, 1198`. Implementation-work
sequencing per impl Codex `id=1201` (Sanity Gate → review → code/tests).

**Sequencing context**: 6.3b is the next "Separated dynamics modes"
sub-unit after 6.3a static FA traction preflight (`52622b4`,
`f43c619`) and the 6.4-open ECM preflight series (`7dc1767`,
`e4a71c9`, `8b2c1ab`, `22b9456`) plus harness/run script
(`04ee5a7`, `9b3dac5`, `f8cdff3`). 6.3b is a deterministic wrapper
around 6.3a, not a replacement.

---

## 0. Scope (locked, no scope creep allowed)

### In scope (this unit only)

- Deterministic per-FA effective-rate computation from caller-supplied
  protrusion state and a typed multiplier table.
- Per-FA delegation to 6.3a `step_focal_adhesions_static` (or
  behavior-equivalent grouping by identical effective parameter
  tuples).
- Pre-step `dt · max(effective_rate)` gate over all FAs.
- Diagnostics: `linked_missing` count, `reciprocal_missing` count,
  multiplier histogram, max effective rate per rate name.
- Asymmetric link consistency: `FocalAdhesionState.linked_protrusion_id`
  is the authoritative read path; `ProtrusionEvent.associated_adhesion_ids`
  reciprocity is diagnostic-only (not failure).

### Explicitly out of scope (will FAIL gate if introduced)

- 6.3a rate/traction algebra reimplementation (delegate only;
  identical-effective-param grouping is the only allowed
  optimization).
- FA state-label change. The state label is exactly the input state;
  6.3b makes no automatic transitions.
- `traction_scale_nN` mutation. The traction magnitude formula
  (`traction_scale_nN · maturity · bound_fraction`) is 6.3a's
  unchanged contract.
- Protrusion / ECM / contour write. One-way protrusion → FA only.
- RNG / stochastic event sampler. Caller supplies the protrusion
  registry; 6.3b reads it.
- `force_candidate_nN` scaling. Out of 6.3b scope; reserved for a
  later force-aware unit.
- Closed-loop Tier 1 behavior (FA→ECM, ECM bias, protrusion force
  on contour). Reserved for the next milestone after the six-item
  Closed-Loop ECM Gate.
- Project default multiplier values. Reference biology table lives in
  the locked design doc only; runtime requires explicit
  caller-supplied multipliers.

---

## 1. Dimensional analysis

Inputs and outputs (every variable's unit explicit):

| Symbol | Meaning | Unit | Source |
|---|---|---|---|
| `multiplier(state, rate_name)` | per-state per-rate scalar boost | dimensionless | `ProtrusionStateMultipliers.lookup` (caller-supplied table) |
| `base_rate` (k_maturity / k_bind / k_unbind) | 6.3a base rate | 1/s | `FocalAdhesionDynamicsParameters` (caller-supplied) |
| `effective_rate = base_rate · multiplier` | per-FA per-rate boosted rate | 1/s | computed |
| `dt_fa_s` | integration step | s | `FocalAdhesionDynamicsParameters` (caller-supplied) |
| `dt · max(effective_rate)` | gate metric | dimensionless | computed; compared against `_DT_RATE_SAFETY_MARGIN = 0.5` |
| `traction_force_nN_xy` | per-FA cell-on-substrate force | nN | 6.3a delegate |
| `linked_protrusion_id` | FA → protrusion key | string | `FocalAdhesionState` |
| `protrusion.event_type / state` | enum used to pick multiplier row | `ProtrusionState` literal | `ProtrusionEvent` |

### Reductions

- **Multiplier**: dimensionless by construction (lock §2 type contract).
- **Effective rate**: `[1/s] · [dimensionless] = [1/s]` ✓.
- **dt · effective rate**: `[s] · [1/s] = [dimensionless]` ✓ — matches
  the 6.3a numerical-sanity gate format identically.
- **Traction algebra**: not touched by 6.3b. Per-FA traction reduction
  remains `[nN] · [dimensionless] · [dimensionless] · [dimensionless]
  = [nN]` per 6.3a §1.

### Status

PASS — every quantity unit-chain reduces through the multiplier in
the same unit as the corresponding 6.3a input. No new force, no
new traction reduction, no new comparison across mismatched units.
Hard Rule 10 trivially satisfied (no per-volume vs per-area
comparison; everything stays in the rate space already cleared by
6.3a).

---

## 2. Boundary cases

Every failure case below either reuses an existing 6.3a / schema
failure_kind or introduces an explicit new one named in the lock
artifact §2.

| Case | Guard | Failure-kind |
|---|---|---|
| Empty `adhesions` list | early-return with empty result, no FA iteration | (no exception) |
| `multipliers` table contains state key not in `ProtrusionState` enum (`"growing"/"stalled"/"retracting"/"ended"`) | `ProtrusionStateMultipliers.validate()` raises | `multiplier_table_unknown_key` |
| `multipliers` table contains rate key not in `RateName` enum (`"k_maturity_per_s"/"k_bind_per_s"/"k_unbind_per_s"`) | `validate()` raises | `multiplier_table_unknown_key` |
| `multipliers` value is `bool` | `validate()` raises | `multiplier_table_unknown_key` (sub-message: bool not allowed) |
| `multipliers` value is non-finite (NaN / ±inf) | `validate()` raises | `multiplier_table_unknown_key` (sub-message: non-finite) |
| `multipliers` value is negative | `validate()` raises | `multiplier_table_unknown_key` (sub-message: negative) |
| `multipliers` table is missing a known state | `lookup()` returns 1.0 (neutral) — no failure | (no exception; absent = neutral) |
| `multipliers` table is missing a known rate under a present state | `lookup()` returns 1.0 (neutral) — no failure | (no exception; absent = neutral) |
| `protrusion_registry` is missing FA's `linked_protrusion_id` | per-FA lookup fails before stepping | `linked_protrusion_missing` |
| Linked `ProtrusionEvent.state` is not in the typed `ProtrusionState` enum (`growing/stalled/retracting/ended`) | per-FA enum check raises before lookup | `linked_protrusion_invalid_state` |
| Empty `adhesions` list with invalid base `params` (e.g. bool `dt_fa_s`, negative `traction_scale_nN`) | upfront base-params validation runs before the empty-FA early return | inherited 6.3a failure_kinds (`dt_invalid`, `traction_scale_invalid`, `rate_invalid`, `max_traction_invalid`) |
| FA has `linked_protrusion_id is None` (FA not linked at all) | treated as no protrusion influence; multiplier 1.0 across all rates (effectively 6.3a behavior for that FA) | (no exception; equivalent to base 6.3a) |
| `ProtrusionEvent` lists FA in `associated_adhesion_ids` but FA does not link back | ignored by one-way read | (no exception; diagnostic only — `reciprocal_missing` counter) |
| FA links back via `linked_protrusion_id` but `ProtrusionEvent.associated_adhesion_ids` does not list the FA | recorded in diagnostics | (no exception; `reciprocal_missing` counter) |
| `dt_fa_s · max(effective_rate over all FAs) > _DT_RATE_SAFETY_MARGIN` | raise BEFORE any stepping | `dt_rate_violation` (reused from 6.3a) |
| All 6.3a `FocalAdhesionDynamicsParameters` validation failures | inherited at construction time | (inherited 6.3a failure_kinds) |
| All 6.3a per-FA validation failures (centroid, axis, position, etc.) | inherited from `step_focal_adhesions_static` delegate | (inherited 6.3a failure_kinds) |

### Status

PASS — every boundary case has an explicit failure_kind or a
documented "absent = neutral" rule. No silent typo defaults
(typed-key validation prevents the `"k_bind"` vs `"k_bind_per_s"`
trap), no silent missing-link defaults (`linked_protrusion_missing`
is loud), and no auto-shrink on the dt-rate violation.

---

## 3. Conservation invariants

### Per-FA action-reaction (delegated to 6.3a)

For every FA `i`:
`traction_force_nN_xy[i] + substrate_reaction_nN_xy[i] = (0, 0)`
exactly to float64 round-off. 6.3b does not touch the traction
algebra; it only changes which effective-rate parameters reach
6.3a. Newton-3 is therefore preserved by delegation.

### Aggregate Newton-3

`sum_i traction_force_nN_xy[i] + sum_i substrate_reaction_nN_xy[i]
= (0, 0)` from the per-FA cancellation.

### Provenance / immutability

- The input `adhesions` tuple, `protrusion_registry`, and
  `multipliers` are not mutated. The result returns a fresh tuple
  of `FocalAdhesionState` (frozen dataclasses, 6.3a contract).
- The supplied `centroid_um_xy` and `params` are not mutated
  (6.3a guarantees this; 6.3b honors it transitively).

### Determinism

- 6.3b is deterministic. No RNG call site. Identical inputs produce
  bit-identical outputs.
- `lookup()` is deterministic (dict.get with a fixed default).
- The optional grouping optimization (group FAs by identical
  effective-parameter tuple, then call 6.3a once per group) must
  be behavior-equivalent to per-FA delegation, with input order
  preserved in the returned arrays.

### Status

PASS — Newton-3 inherits from 6.3a delegation; 6.3b adds no new
conservation burden. Determinism is enforceable and tested via the
neutral-multiplier regression test (`test_neutral_multipliers_match_6_3a_baseline`).

---

## 4. Numerical sanity

### Stability bound across all FAs

The 6.3a per-FA explicit-Euler update of the form
`x_new = x + dt · k_eff · (1 - x)` (or analogous toward zero) is
linearly stable for `dt · k_eff < 2` and energy-monotone for
`dt · k_eff ≤ 1`. The runtime contract reused identically from
6.3a is **`dt_fa_s · max(effective_rate over all FAs) ≤ 0.5`** —
the same named module-level constant `_DT_RATE_SAFETY_MARGIN = 0.5`
imported from `acs.v2.active_contour`.

The maximum is taken **over all FAs and all rate names** before any
stepping. This is a stronger guarantee than per-FA dt-rate checking
because a heterogeneous FA list with different effective rates per
FA could otherwise hide one violating FA behind passing ones.

If the gate raises `dt_rate_violation`, no FA is stepped (the wrapper
fails closed). Auto-shrink is forbidden — the caller must reduce
`dt_fa_s`, lower the multiplier table, or lower the base rates.

### Float precision

float64 throughout (6.3a inheritance). The multiplier multiplication
`base_rate * multiplier` is one multiplication per FA per rate name,
trivially numerically stable (no subtraction, no division, no
exponentiation).

### Status

PASS — the dt-rate gate is computed correctly across all FAs before
stepping; the `_DT_RATE_SAFETY_MARGIN = 0.5` constant is reused from
6.3a's already-reviewed safety margin rationale. No new tunable
constants land in 6.3b code.

---

## 5. Sign / sense check

| Term | Direction | One-line check |
|---|---|---|
| Multiplier ≥ 0 | Caller cannot supply a negative multiplier; rejected at `validate()`. | one-line check: `multiplier >= 0`. |
| Multiplier 1.0 | Neutral — `effective_rate == base_rate`, behavior identical to 6.3a. | regression test `test_neutral_multipliers_match_6_3a_baseline`. |
| Multiplier 0.0 | Effective rate 0 — no rate-driven update for that FA/rate. **Does NOT change FA state label** (6.3a's no-auto-state contract preserved). | tested via a multiplier-0 case. |
| Multiplier > 1 | Boost — effective rate exceeds base. Constrained only by the dt-rate gate. | tested via the growing-state boost test. |
| `linked_protrusion_id is None` | Treated as "no protrusion influence on this FA" → effectively 6.3a behavior. | tested via a mixed FA list with one unlinked FA. |
| `protrusion.state` lookup direction | One-way: `FA → registry → state → multiplier`. The registry never reads back from FA. | structural — no write path on the protrusion side. |

### State-table sense (delegate to 6.3a)

The 6.3a §5 state table (which rates apply per state) is unchanged:
- `unbound` / `released`: still no rate updates regardless of multiplier.
- `nascent` / `mature`: `k_maturity_per_s` and `k_bind_per_s` apply
  if supplied.
- `slipping`: only `k_unbind_per_s` applies.

A multiplier table with a non-1.0 entry for a rate that does not
apply to the FA's current state is **silently neutral** for that
FA — the multiplier is looked up but never reaches the rate update,
because the 6.3a state table gates which rate applies. This is
intentional: the multiplier table is the protrusion-side biology;
the state table is the FA-side biology. The wrapper should not
override the FA-side biology.

### Status

PASS — sign / sense covered by the 6.3a inherited state table plus
the `multiplier ≥ 0` and `multiplier 1.0 = neutral` properties
above. The "multiplier 0 ≠ state change" property is the key
non-obvious sign check and is named explicitly so a future test
catches any drift.

---

## 6. Measurement-protocol consistency

6.3b's measurement contract is **exactly** 6.3a's:
- `cell_force_nN_xy[i]` is the cell-on-substrate force at FA `i`,
  inward radial-to-centroid by default.
- `substrate_reaction_nN_xy[i]` is the algebraic negation
  (action-reaction).
- `radial_components[i] = traction · n_hat_radial_outward[i]`,
  `tangential_components[i] = traction · n_hat_tangential_ccw[i]`.
- The decomposition reconstruction
  `radial · n_hat_radial_outward + tangential · n_hat_tangential_ccw
  == traction_force_nN_xy` holds to float64 round-off — inherited
  from 6.3a.

### Per-FA delegation must preserve order

The output arrays / tuples are **input-order preserved**: the FA at
input index `i` appears at output index `i`. The optional grouping
optimization (group by identical effective-parameter tuple, run
6.3a once per group) must internally re-scatter to input order
before returning. A regression test
(`test_per_fa_input_order_preserved_in_output`) enforces this.

### Diagnostics added by 6.3b

The 6.3a `diagnostics` dict is extended with:
- `linked_missing`: count of FAs whose `linked_protrusion_id` was
  not found in the registry. **Always zero on success** (because
  missing-link raises before stepping); kept as a key for symmetry
  with `reciprocal_missing`.
- `reciprocal_missing`: count of FAs whose `linked_protrusion_id`
  resolved but whose linked protrusion did not list the FA in
  `associated_adhesion_ids`.
- `multiplier_histogram`: per-rate-name histogram (or summary stats)
  of effective multiplier values across FAs. Documents the actual
  multiplier distribution so a later debug session can check for
  unexpected uniformity / clustering.
- `max_effective_rate_per_name`: per-rate-name maximum effective
  rate across all FAs, in `[1/s]`. The dt-rate gate uses the
  per-name maximum's overall maximum for the gating decision; this
  diagnostic exposes the per-name breakdown so a tune-up session
  can identify which rate caused the violation.

These diagnostics are additive — every 6.3a diagnostic key is
preserved verbatim. A downstream consumer expecting only 6.3a keys
sees a superset and is unaffected.

### Status

PASS — the measurement contract is identical to 6.3a's; the
diagnostics extension is additive only. No new measurement
modality, no new decomposition, no new comparison against
experimental data.

---

## 7. Magic-Number Block check

Every numeric in the 6.3b implementation is one of:

| Symbol | Source | Magic-Number Block status |
|---|---|---|
| `multiplier` table values | runtime caller-supplied (no project default) | **not a magic number** — caller-supplied test inputs only |
| `_DT_RATE_SAFETY_MARGIN = 0.5` | imported from `acs.v2.active_contour` | already-reviewed named constant; reused identically from 6.3a / P1 alpha; not chosen to make any 6.3b test pass |
| Reference biology table (`growing 2.0/1.5/1.0`, etc.) | docs/v2/v2_63b_protrusion_coupling_locked.md §3 | **doc only, not code**; runtime ignores it; explicitly noted as caller guidance |
| Test multiplier values (e.g. 2.0 in `test_growing_multiplier_boosts_maturity_and_bind`) | per-test caller inputs | **not a magic number** — chosen for clarity in the test, not to make production behavior pass |

Magic-Number Block test pass per locked principle:

1. **Derivable**: every value in code is either runtime caller input
   or the inherited `_DT_RATE_SAFETY_MARGIN` from 6.3a / P1 alpha.
   No new physics constant is introduced.
2. **Grid-invariant**: 6.3b adds no grid (FA list is not a grid);
   the multiplier multiplication is point-wise per FA, trivially
   independent of any discretisation.
3. **Not fitting**: no value in code was chosen to make a specific
   test pass. The reference biology table sits outside code; tests
   exercise the wrapper's algebra, not biology.

### Status

PASS — Magic-Number Block clean. The biology values from the
reference table are explicitly kept out of code per the lock §3
contract; runtime requires explicit caller-supplied multipliers.

---

## 8. Visual deliverable plan (post-commit)

Aligns with the precedent set by `scripts/run_p1_alpha_gate.py` and
`scripts/run_ecm_ol_harness.py`. **Out of scope for the Sanity
Gate commit itself** — the Sanity Gate gates the code commit; the
visible-deliverable run script is a separate small unit after the
code commit and Codex review.

Planned (not committed by this gate):
- `scripts/run_protrusion_coupled_fa_harness.py` — exercises a few
  caller-supplied multiplier tables against a small FA list with a
  prescribed protrusion registry timeline. Persistent artifacts in
  `runs/<UTC>_protrusion_coupled_fa/`. Status table aggregator
  (PASS/FAIL) per the `f8cdff3` blocker-fix pattern.

This visual deliverable plan is documented here so the post-Sanity
Gate code commit has a known follow-on; the gate does not require
the script to land before clearing.

---

## 9. Test catalog (~12, per lock §6)

Each test exercises one Sanity Gate item or one forbidden behavior.

| # | Test name | Gate item |
|---|---|---|
| 1 | `test_neutral_multipliers_match_6_3a_baseline` | §3 determinism + §5 multiplier 1.0 = neutral |
| 2 | `test_growing_multiplier_boosts_maturity_and_bind` | §1 unit chain + §5 multiplier > 1 = boost |
| 3 | `test_retracting_neutral_no_boost` | §5 multiplier 1.0 = neutral on a non-growing state |
| 4 | `test_unknown_state_key_raises_multiplier_table_unknown_key` | §2 typed key validation |
| 5 | `test_unknown_rate_key_raises_multiplier_table_unknown_key` | §2 typed key validation |
| 6 | `test_negative_multiplier_raises_at_validate` | §2 multiplier ≥ 0 |
| 7 | `test_bool_multiplier_value_raises_at_validate` | §2 bool not allowed |
| 8 | `test_missing_linked_protrusion_id_raises_linked_protrusion_missing` | §2 linked_protrusion_missing |
| 9 | `test_reciprocal_missing_records_diagnostic_not_raise` | §2 / §6 asymmetric link consistency |
| 10 | `test_dt_rate_violation_uses_max_effective_rate_over_all_fas` | §4 pre-step gate |
| 11 | `test_per_fa_input_order_preserved_in_output` | §6 input order preservation |
| 12 | `test_no_state_label_change_after_step` | §0 forbidden + §5 multiplier 0 ≠ state change |

A regression test (`#1`) uses the existing 6.3a tests' fixtures so a
divergence between 6.3b's neutral-multiplier path and the 6.3a
baseline surfaces immediately.

---

## 10. Gate verdict

§1, §2, §3, §4, §5, §6 PASS. Magic-Number Block (§7) clean. Visual
deliverable plan (§8) is post-commit and does not block. Test
catalog (§9) maps every gate item to at least one test.

The gate clears for executable code in two commits:

1. `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py` — new
   module with `ProtrusionStateMultipliers` dataclass +
   `step_protrusion_coupled_focal_adhesions(...)` function.
   `FocalAdhesionDynamicsError(failure_kind=...)` extension for
   the new failure_kinds (`multiplier_table_unknown_key`,
   `linked_protrusion_missing`).
2. `tests/v2/test_v2_protrusion_coupled_focal_adhesion.py` — ~12 tests
   per §9 catalog above.

Plus optional `acs/v2/__init__.py` + `acs/v2/dynamics/__init__.py`
exports.

If Codex review surfaces a missing reduction or hidden numeric, the
gate flips to BLOCKER with the three-options template.

### Post-implementation status (audit 2026-05-04 KST)

This section was originally titled "Outstanding before code
lands" and is preserved as historical scaffolding. All five
review focus points listed below were reviewed and cleared by
Codex (impl `id=1207` cleared the Sanity Gate; impl `id=1214`
cleared the code; impl `id=1219` cleared the persistent runner):

1. typed `ProtrusionStateMultipliers` validation — unknown keys
   fail not neutral ✓ (`9a20fc9`, `deed45a`)
2. per-FA delegation to 6.3a or behavior-equivalent grouping;
   no inline reimplementation ✓ (`9a20fc9`)
3. `linked_protrusion_missing` failure on missing linked id ✓
   (`9a20fc9`)
4. reciprocal mismatch diagnostic only (warnings, not failures) ✓
   (`9a20fc9`)
5. pre-step `dt_fa_s · max(effective_rate over all FAs)` gate ✓
   (`9a20fc9`, `deed45a` empty-FA bypass fix)

The lock artifact (`docs/v2/v2_63b_protrusion_coupling_locked.md`)
landed together with this gate doc in commit `6d1e12b`
(whitespace hygiene fix `dc5043a`); the source-of-truth
dependency is auditable on disk.

The persistent visible-deliverable runner landed at `50033c3`
+ `2a00452` (docstring align fix).

Future modification of `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py`
must update or supersede this gate first if any §0 forbidden item
or §1–§6 contract changes.

---

## 11. References

- Lock artifact: `docs/v2/v2_63b_protrusion_coupling_locked.md`
- Design brief (opening position): `docs/v2/v2_63b_protrusion_coupling_design_brief.md`
- 6.3a Sanity Gate (pattern + delegate target):
  `docs/v2/v2_focal_adhesion_dynamics_sanity_gate.md`
- P1 alpha Sanity Gate (template + named-constant precedent):
  `docs/v2/v2_p1_active_contour_sanity_gate.md`
- 6.3a code (delegate target):
  `acs/v2/dynamics/focal_adhesion.py`
- ProtrusionEvent schema: `acs/v2/protrusion.py`
- FocalAdhesionState schema: `acs/v2/focal_adhesion.py`
- Forward roadmap: `docs/v2/v2_phase1_forward_roadmap.md`
- ECM-OL persistent run pattern (precedent for §8 visible
  deliverable): `runs/20260504T035716Z_ecm_ol/` +
  `scripts/run_ecm_ol_harness.py`
- P1 alpha persistent run pattern: `runs/20260503T1833Z_p1_alpha/`
  + `scripts/run_p1_alpha_gate.py`
- Adversarial debate posture: memory
  `feedback_aggressive_design_debate.md`, PI `id=809`
- Cadence rule: memory `cadence_promise_must_send_even_when_idle.md`,
  PI `id=1025/1031`
- No-standby rule: memory `no_standby_after_completion.md`, PI
  `id=1057`

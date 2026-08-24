# ALEPH-PORT-500 — protocol manifest and the observation-operator seam

| Field | Value |
|---|---|
| Lane | L5 (observation layer) |
| Target | `aleph/observe/manifest.py`, `aleph/observe/operator.py` |
| Reference consulted | **NONE.** No `ffn_cellsim` module was read for this entry. |
| Source | `docs/manuscripts/PROJECT_ALEPH_PROPOSAL_R3_2026-07-30.docx` §8.1, classified CARRY in `bootstrap/docs/MANUSCRIPT_SALVAGE_LEDGER.md` |
| Port class | **ORIGINAL.** Manuscript science, no inherited code. |
| Written | 2026-07-30 |
| Status | LANDED |

## 1. Why this is a ledger entry at all

PLAN §0.2 requires an entry per *ported* leaf. Nothing is ported here — that is the point of the
entry. The reference project has an observation layer (`ac/engine/observe/`) and it was deliberately
**not** consulted for the manifest or the operator seam, because the manuscript's §8.1 protocol
declaration is a stronger statement than anything in the reference, and reading the reference first
would have anchored the design to it. This entry records that choice so a later reader does not have
to guess whether the omission was an oversight.

## 2. The claim being encoded

A number produced by a measurement is meaningless without the protocol that produced it. Two
indentation moduli measured at different rates on the same cell are not two estimates of one
quantity; they are two different quantities. §8.1 lists what must be declared before a measurement
can be compared to anything.

Six field families, thirty-four fields:

| Family | Fields |
|---|---|
| geometry and control | probe geometry, probe location, probe direction, control mode, contact model, adhesion state |
| time and amplitude | rate, frequency, dwell, loading history, max strain, preconditioning, sampling cadence |
| environment | temperature, medium, osmotic state, substrate, confinement, pharmacology |
| observation | raw output schema, raw output units, calibration |
| analysis | segmentation, inverse model, regularization, fit window, noise model, resolution limit |
| evidence | specimen identity, batch, citation audit, applicability domain, evidence role, operator version |

## 3. The design decision that carries the weight

**`NOT_APPLICABLE` and `NOT_RECORDED` are distinct typed values, and neither is `None`.**

For the passive fluctuation protocol (`ALEPH-DQ-102`) most geometry-and-control fields have no
referent: there is no probe, so there is no probe geometry. That is a *property of the protocol*.
It is categorically different from a nanoindentation run where a probe existed and nobody wrote down
its radius. The first is complete; the second is a hole in the record. A schema that renders both as
`None` — or that lets the second be omitted — destroys exactly the distinction that makes the
manifest worth having, and it destroys it silently.

So:

- `Declaration.NOT_APPLICABLE` — the protocol has no such quantity. An assertion.
- `Declaration.NOT_RECORDED` — the quantity existed and was not captured. A confession.
- **No field has a default.** Adding a field to `ProtocolManifest` breaks every construction site in
  the repository, which is the intended cost: a new declarable is a new thing every existing
  protocol must answer for. `tests/observe/test_manifest.py::test_no_manifest_field_has_a_default`
  enforces this structurally, so the property survives the next person to add a field.
- `passive_fluctuation_manifest()` is a named factory that spells out all thirty-four values in
  source, in one place, auditable by reading it. There is deliberately no generic
  "fill the rest with NOT_APPLICABLE" helper: that helper is how an unrecorded field becomes an
  inapplicable one without anyone deciding.

## 4. `manifest_hash()`

SHA-256 over canonical JSON (sorted keys, no whitespace slack) of the full field mapping, prefixed
with a schema tag so a future field-set change cannot collide with the present one. Sixteen hex
characters are exposed as `short_hash` for logs; the full digest is what travels with results.

The hash is the join key between a reported number and the protocol that produced it. Its one
required property is that it separates manifests that differ in any field — including differing only
in `NOT_APPLICABLE` versus `NOT_RECORDED`, which is asserted directly.

## 5. The operator seam

`ObservationOperator` is a Protocol with two *separate* methods:

    raw_observable(state, context) -> RawObservable | Refusal
    apparent_quantity(raw)         -> ApparentQuantity | Refusal

They are separate because conflating them is how a simulation gets compared against the wrong thing.
A real experiment does not hand you a modulus; it hands you a voltage, which a calibration turns into
a deflection, which a contact model turns into a modulus. Every one of those steps is a modelling
choice, and the number at the end is *apparent* — it is what this analysis pipeline reports, not what
the material is. A simulator that computes the true quantity directly and compares it to a published
one has skipped every step where the disagreement lives.

`ApparentQuantity` therefore carries `analysis_chain`, an ordered list of the transformations
applied, and it carries the manifest hash. A quantity that cannot say how it was produced cannot be
compared.

### Refusal is a return value

`Refusal` is a frozen dataclass and is **not** a subclass of `BaseException`; a test asserts this, so
it can never quietly become one. Three failure modes are being ruled out at once:

1. **Exception-as-control-flow.** A caller who forgets a `try` gets a crash, and a caller who
   remembers gets a bare `except` that swallows real bugs too.
2. **A silent NaN.** NaN propagates, and `mean()` of an array containing one is NaN — so it is
   loud eventually and in the wrong place. It also carries no reason.
3. **A `None` return.** Indistinguishable from "not implemented yet".

A `Refusal` carries a typed `RefusalCode`, a human sentence, and a machine-readable detail mapping,
and the type signature forces the caller to consider it.

`ApplicabilityDomain` holds named closed intervals plus required state owners. `check()` returns
`Refusal | None`, and the refusal names the offending quantity, its value, and the bound it violated.

## 6. Controls shipped

The deliberately-failing negative control for this entry is
`tests/observe/test_manifest.py::test_not_applicable_and_not_recorded_hash_differently`: it fails if
the two declarations ever collide, which would let a hole in the record certify as a complete one.
Its structural companion is `test_no_manifest_field_has_a_default`, which fails the moment any field
acquires a default. On the operator side the negative control is
`tests/observe/test_operator.py::test_out_of_domain_returns_a_refusal_not_an_exception`, which fails
if an out-of-domain call ever answers instead of declining.

`M` = `tests/observe/test_manifest.py`, `O` = `tests/observe/test_operator.py`.

| Control | Test | Asserts |
|---|---|---|
| Positive | M `test_the_six_families_cover_every_field_exactly_once` | the family map is complete and disjoint |
| Positive | M `test_there_are_thirty_four_fields` | the §8.1 count, pinned |
| Positive | M `test_manifest_hash_is_stable_across_construction` | identical content, identical digest |
| Positive | M `test_manifest_hash_changes_when_any_field_changes` | all 34 fields, one at a time |
| Positive | M `test_canonical_payload_is_sorted_json_carrying_the_schema_tag` | a digest reproducible by hand |
| Positive | M `test_amending_changes_the_hash_and_leaves_the_original_alone` | amendment is honest and frozen |
| Positive | M `test_passive_fluctuation_manifest_declares_every_field` | all 34 present and typed |
| Positive | M `test_passive_fluctuation_declares_no_probe_as_not_applicable` | the §3 distinction on the real manifest |
| Positive | M `test_passive_fluctuation_declares_no_loading_but_does_declare_a_cadence` | no rate, but a sampling cadence |
| Positive | M `test_passive_fluctuation_declares_no_inverse_model_or_regularization` | the ALEPH-DQ-102 argument, encoded |
| Positive | M `test_a_manifest_of_pure_not_applicable_is_complete` | mostly-inapplicable is a description, not a defect |
| Positive | M `test_not_recorded_fields_are_listed_and_not_applicable_ones_are_not` | the two lists stay apart |
| Positive | O `test_operator_inside_domain_returns_a_quantity` | in-domain path |
| Positive | O `test_raw_observable_is_not_the_apparent_quantity` | the §5 separation is observable |
| Positive | O `test_apparent_quantity_records_its_analysis_chain` | the chain travels with the number |
| Positive | O `test_bounds_are_closed_intervals` / `test_an_unbounded_quantity_is_ignored_not_rejected` | domain semantics |
| Positive | O `test_the_toy_operator_satisfies_the_protocol` / `test_array_state_satisfies_the_accepted_state_protocol` | structural conformance |
| **Negative** | M `test_not_applicable_and_not_recorded_hash_differently` | the distinction is not merely documented |
| **Negative** | M `test_no_manifest_field_has_a_default` | a new field cannot be silently defaulted |
| **Negative** | M `test_omitting_a_field_is_a_construction_error` | omission fails loudly |
| **Negative** | M `test_none_is_not_a_permitted_value` | `None` cannot express either declaration |
| **Negative** | M `test_manifest_is_frozen` | a hash cannot outlive its content |
| **Negative** | M `test_uncaptured_environment_defaults_to_the_confession_not_the_assertion` | never claim "no substrate" by omission |
| **Negative** | M `test_temperature_has_no_default_in_the_passive_factory` | a spectrum without a temperature is meaningless |
| **Negative** | M `test_amending_an_unknown_field_is_refused` | typos do not become fields |
| **Negative** | O `test_refusal_is_not_an_exception` | structural; refusal can never become control flow |
| **Negative** | O `test_out_of_domain_returns_a_refusal_not_an_exception` | no raise, no number to pick up |
| **Negative** | O `test_missing_state_owner_is_refused_before_any_arithmetic` | refusal precedes computation |
| **Negative** | O `test_non_finite_input_is_refused_rather_than_propagated_as_nan` | no silent NaN |
| **Negative** | O `test_non_finite_value_is_always_refused_even_when_unbounded` | finiteness is unconditional |
| **Negative** | O `test_an_apparent_quantity_with_no_analysis_chain_cannot_be_built` | a number must say where it came from |
| **Negative** | O `test_observe_short_circuits_on_the_first_refusal` | analysis never runs on a refusal |
| **Negative** | O `test_refusal_is_truthy_so_a_bare_if_cannot_be_used_to_detect_it` | pins the trap |

## 6a. Units, domains, singular cases, invariants

**Units.** The manifest carries no bare numbers. Every field holding a physical quantity holds it as
text *with the unit written in* — `"296.0 K"`, not `296.0` — because a bare float in a protocol
record is a unit waiting to be guessed. `ApparentQuantity.units` and `RawObservable.units` are
likewise documented strings; when L1's dimension algebra settles they are the attachment point.

**Domain.** `ApplicabilityDomain.bounds` holds closed intervals in the operator's declared units. A
quantity with no declared bound is *ignored*, not rejected: silence about a quantity is not a
constraint on it. Non-finiteness is refused unconditionally, bounded or not.

**Singular cases.** Three, all typed rather than silent: a manifest field that is `None`
(`TypeError`); an `ApparentQuantity` with an empty `analysis_chain` (`ValueError` — a number that
claims to have come from nowhere); and an out-of-domain or missing-owner observation (`Refusal`,
which is a return value and asserted not to be a `BaseException` subclass).

**Invariants.**

- Every one of the 34 fields is present in `as_mapping()` on every manifest, always. No field is
  optional and none has a default.
- `manifest_hash()` is a pure function of content: identical content gives an identical digest, and
  changing *any* single field changes it — tested field by field across all 34.
- `NOT_APPLICABLE` and `NOT_RECORDED` never hash alike.
- `not_recorded_fields()` and `not_applicable_fields()` are disjoint, and `is_complete()` depends on
  the first only.
- A `Refusal` never carries a `value` attribute, so a caller cannot fish a number out of one.

## 6b. Numerical and precision envelope, and production-backend residency

**Numerical.** There is no floating-point computation in either module. `manifest_hash()` is
SHA-256 over canonical JSON with sorted keys and fixed separators, so it is exact, reproducible and
independent of float formatting. The only float handling is `ApplicabilityDomain.check`, which does
finiteness and interval comparisons in float64 with no arithmetic beyond a cast.

**Residency.** Host-side, CPU only, pure Python plus `numpy` for array typing. Nothing here touches
a device, and nothing here should: a protocol manifest is metadata, and an applicability check is a
handful of comparisons whose cost is invisible next to the observation it guards. When the ratified
backend (`ALEPH-DQ-107`) lands, the transfer boundary sits *below* this layer — an operator receives
whatever the `AcceptedState` Protocol hands it, and where that array lived is the state layer's
business, not the observer's. No transfer is performed or implied by this entry.

## 7. Compliance with PLAN §0.2

1. No source text carried — no source was read. ✔
2. Law re-derived — manuscript §8.1 restated in Aleph's terms. ✔
3. New prose throughout. ✔
4. Positive and failing negative controls (§6). ✔
5. Entry written before the code. ✔
6. No `ffn_cellsim` import; no `validation` import (`tests/firewall/`). ✔

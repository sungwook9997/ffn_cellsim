# Port ledger entry ALEPH-PORT-602 — run record, build stamp, census guard

| Field | Value |
|---|---|
| Lane | L6 (provenance, artifacts, evidence) |
| Aleph target | `aleph/artifacts/stamp.py`, `aleph/artifacts/record.py`, `aleph/artifacts/__init__.py` (re-exports only) |
| Reference read | `/Users/sw1/ffn_cellsim/ffn_sim/ac/engine/observe/artifact.py` (409 lines, `run-record@2`) |
| Classification | **CONCEPT_ONLY — not a code port** |
| Written | 2026-07-30, before any code landed |

## What was taken

Three ideas:

1. A run record should stamp the build it was measured on, and should say *how* it learned the
   commit rather than presenting every commit as equally verified.
2. A census can be arithmetically perfect and still overclaim, so the guard must be on the
   *vocabulary*, not on the numbers.
3. A config hash belongs in the record so a changed config cannot keep an old identity.

## What was NOT taken

- `observation_artifact()`'s signature, its 18 keyword arguments, its `t0` / `force_channels` /
  `parameter_provenance` blocks, and its gate-verdict assembly. Aleph's record is a typed dataclass,
  not a `dict[str, Any]` builder, because an untyped record is how a rung became a hand-typed string.
- `timing_block()` and its comparability stamping. Aleph carries timing as a frozen dataclass; the
  convergence-conditioning of a cost belongs to L3's transaction, not to L6.
- `POPULATION_KNOBS`. That table is a map of the *reference's* config vocabulary (`membrane_subdiv`
  vs `membrane_subdivisions`) and has no meaning in Aleph. Absorbing it would violate the boundary
  rule: Aleph owns interfaces, never the provider's internal spellings.
- All prose, all rationale strings, all examples.

## Re-derivation — and the defect Aleph fixes

The reference's `build_stamp` resolves a declared-vs-git disagreement by *keeping git's commit* and
adding a `declared_vs_git_mismatch` note beside it. That is still a silent preference: a reader who
looks at `commit` gets one answer and never learns it was contested.

Aleph refuses to resolve it. `BuildStamp` stores `git_commit` and `declared_commit` as separate
fields; the `.commit` property **raises** `CommitDisagreementError` when the two disagree, and
serialization emits `commit: null` with a `disagreement` block carrying both. Code that wants a
commit must handle the contest explicitly.

The census guard was rebuilt around a checkable contract rather than a name list alone: a
`fraction`-like key at exactly `1.0` must be accompanied by a per-owner breakdown covering every
name in `census["owners"]`, each itself at `1.0`. This is the guard that would have caught the
reference's real 2026-07-29 incident (two runs 1.73x apart in node count both stamping
`fraction_of_native: 1.0`) by construction rather than by naming that one key.

## Controls shipped

- Positive: `test_build_stamp_from_git_repo`, `test_run_record_round_trips_to_json`.
- Negative (must refuse): `test_declared_and_git_disagreement_is_recorded_not_resolved`,
  `test_commit_property_raises_on_disagreement`, `test_banned_census_key_is_refused`,
  `test_fraction_at_one_without_breakdown_is_refused`,
  `test_fraction_breakdown_missing_an_owner_is_refused`.

## Independence

Imports `hashlib`, `json`, `subprocess`, `pathlib`, `dataclasses`, `enum` only.
Passes with `ffn_cellsim` absent from `sys.path`.

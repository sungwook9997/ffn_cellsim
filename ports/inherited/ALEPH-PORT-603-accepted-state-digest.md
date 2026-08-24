# Port ledger entry ALEPH-PORT-603 — accepted-state digest

| Field | Value |
|---|---|
| Lane | L6 (provenance, artifacts, evidence) |
| Aleph target | `aleph/artifacts/digest.py` |
| Reference read | `/Users/sw1/ffn_cellsim/ffn_sim/virtual_cell/cell_state_posterior.py` (`record_hash`, `to_dict`) |
| Classification | **DEFECT_STUDY — nothing ported, the reference is a counterexample** |
| Written | 2026-07-30, before any code landed |

## Why this entry exists

This is the one reference asset read purely to learn what *not* to build. The reference's
`record_hash` is the SHA-256 of `to_dict()`, and `to_dict()` omitted the evidence graph's content
address and the accepted-step verdict's identity. Consequence, found by external review on
2026-07-29: **two posteriors with different provenance DAGs and different attempted step indices
hashed identically.** A content address that cannot see a difference is not a content address.

The failure mode is structural, not a typo: when the digest is "hash of whatever the serializer
happens to emit", every future field is opt-in to the identity, and the default is invisibility.

## Aleph's counter-design

1. **Coverage is declared, not incidental.** `DIGEST_COVERAGE` names the six components the digest
   must cover; `accepted_state_digest()` takes all six as *required keyword arguments with no
   defaults*, so a caller cannot omit one and a future component cannot be silently optional.
   A test asserts the signature and `DIGEST_COVERAGE` agree.
2. **Domain separation and length prefixing.** Every component is hashed under a labelled,
   length-prefixed frame, so no concatenation of two different structures can produce the same byte
   stream. A test proves owner names `("ab", "c")` and `("a", "bc")` do not collide.
3. **Sensitivity is proved per component.** A test perturbs each of the six components in isolation
   and asserts the digest moves. The RNG-position case is the direct counterexample to the reference
   defect: two states differing *only* in RNG stream position must not share a digest.

## Controls shipped

- Positive: `test_digest_is_stable_under_repetition`, `test_digest_is_order_independent_for_owners`.
- Negative (must catch a difference): `test_rng_stream_position_changes_digest` (the named
  counterexample), `test_every_covered_component_changes_the_digest`,
  `test_owner_name_boundary_is_not_ambiguous`.

## Independence

Imports `hashlib`, `struct`, `dataclasses` and optionally consumes numpy arrays by duck typing.
No `ffn_cellsim` import at any depth.

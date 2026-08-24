# Port ledger — ALEPH-PORT-302 · participant / connector protocol

| Field | Value |
|---|---|
| Aleph target | `aleph/runtime/participant.py` |
| Lane | L3 runtime |
| Written | 2026-07-30 (before the code landed, per PLAN §0.2.5) |
| Status | LANDED |
| Evidence class | STRUCTURAL. |

## What was read in `ffn_cellsim` (read-only)

- `ffn_sim/ac/engine/world.py:18` — the participant contract is a tuple of three *method name strings*
  (`snapshot_candidate`, `rollback`, `commit_irreversible`) probed with `callable(getattr(...))`.
- `ffn_sim/ac/engine/dispatch.py:36-72` — `CandidateTask` binds a runtime object to a `method_name`
  string and dispatches with `getattr(self.runtime, self.method_name)`.

## The structural fault this port exists to remove

Both mechanisms above are *name-based*. `dispatch.py:105-128` assigns the connector
`membrane_cortex_contact` to `SurfaceBody.accumulate_mechanics`; that method is callable, so every
structural check passes — but the method body has no slot that evaluates the connector
(`surface_body.py:194,364-389`). A name-based contract can only ever prove that a method *exists*.

Aleph's re-derivation: **a connector cannot be described by a method name.** `Connector.accumulate`
takes exactly two `EndpointHandle`s and must return an `AdjointPair`. Handing a connector to a facade
with no slot for it is then not expressible: there is nothing to pass the two endpoints to, and nothing
to return the pair from. The type signature carries the obligation the name could not.

Corollaries designed in here:

- `Participant.accumulate` returns a `WorkReceipt` rather than `None`. A participant that does nothing
  cannot silently succeed; it has to fabricate a receipt, which the pipeline then cross-checks against
  the backend witness counter (ALEPH-PORT-303).
- `owned_entities() -> EntityCensus` makes ownership *data* rather than a docstring claim, so a receipt
  naming an entity the participant does not own is a detectable fault.
- `AdjointPair` records force, application point, displacement and couple **at each endpoint
  separately**. The two sides are never computed from one another. This is the one idea genuinely worth
  carrying from `ledger.py`'s module prose ("what must NOT happen is the two channels being computed
  from one another, which turns a check into a tautology") — the idea is carried, the code is not.

## Physical content re-derived

For a two-point connector, `AdjointPair` is the discrete statement of Newton's third law plus the
two-force-member theorem:

- force closure: `F_a + F_b = 0` when the pair is the complete momentum exchange (`newton_pair`);
- moment closure about the pair centroid `c = (r_a + r_b)/2`:
  `(r_a − c) × F_a + (r_b − c) × F_b + M_a + M_b = 0`. With `F_b = −F_a` this reduces to
  `(r_a − r_b) × F_a + M_a + M_b`, which vanishes iff the transmitted force is central (along the line
  of the two endpoints) up to the declared couples. It is therefore a real check, not an identity.

## Controls

- Positive control: `test_passive_pair_closes_force_and_moment`.
- Negative (must fail): `test_non_central_force_fails_moment_closure`,
  `test_broken_third_law_fails_force_closure`, `test_receipt_naming_a_foreign_entity_is_a_fault`.

## Provenance and envelope

- **Source repository / commit read:** `ffn_cellsim` at commit `be0e5876`, read-only.
- **Discarded prose:** every original comment and docstring discarded; nothing translated. The two
  reference facts retained (the method-name contract; the two-channels-must-be-independent idea) are
  restated as findings in Aleph's own words, not as carried text.
- **Units, domains, invariants:** forces pN, lengths um, moments pN.um, energies pN.um, time s.
  Invariants enforced at construction: an `AdjointPair` may not name one endpoint twice; every vector
  component must be finite; a `passive` pair may not report active work; an `EntityCensus` may not
  name an entity twice or give it a negative population.
- **Numerical and precision envelope:** float64. Residual *scales* are carried alongside residuals so
  the tolerance applied downstream is relative to the magnitudes that actually cancelled.
- **Residency:** host. `EndpointHandle.force_sink` is a backend array, so the same handle shape works
  when the sink lives on a device.

## Independence

No `ffn_cellsim` import. No method-name strings are shared with the reference; the phase names, the
receipt type, and the endpoint-handle type have no counterpart there.

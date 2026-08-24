# Port ledger — ALEPH-PORT-305 · accepted-step transaction

| Field | Value |
|---|---|
| Aleph target | `aleph/runtime/transaction.py` |
| Lane | L3 runtime |
| Written | 2026-07-30 (before the code landed, per PLAN §0.2.5) |
| Status | LANDED |
| Evidence class | STRUCTURAL, with one exact bit-identity assertion (RNG stream restore). |
| Gate | R3 — "standalone candidate-step transaction on the ratified backend" |

## What was read in `ffn_cellsim` (read-only)

- `ffn_sim/ac/engine/transaction.py` (171 lines) — `CellTransaction.step`, the seven-step ordering,
  the device-resident acceptance predicate, `clock.advance(predicate, dt_phys)`.
- `ffn_sim/ac/engine/world.py` (104 lines) — `CellWorldTransaction`, participant discovery by method
  name, `begin_candidate` / `finalize_candidate`, and the identity-based deduplication of runtimes
  registered under two semantic edges.

## What is genuinely good there, and was carried as an *idea*

1. **One decision, fanned out to everyone.** A single predicate drives every participant's rollback
   and commit, so components solved independently cannot end up on different physical clocks. Carried.
2. **Deduplicate participants by object identity.** One runtime registered under two connector edges
   must be snapshotted and committed once, not twice. Carried in `Pipeline._participating_objects`,
   with a comment saying why (double-count, not performance).
3. **The predicate stays out of Python control flow.** In the reference this is literal — the flag is
   device-resident and both branches launch. Aleph's CPU reference cannot and need not do that, so
   the idea is carried as *structure* instead: `Decision` is a value, `finalize()` is the only method
   that acts on it, and `commit()` refuses a rejected decision.

## What was rejected, and why

- **`bool` acceptance.** The reference's predicate is a device int; a caller learns *whether*, never
  *why*. Aleph's `PredicateVerdict` names every sub-condition, whether it was satisfied, what was
  measured and against what threshold. `AllOf` deliberately does not short-circuit: a step that fails
  two laws at once is diagnostically different from one that fails one.
- **A verdict with no conditions meaning "accept".** In Aleph an empty condition set **rejects**. An
  acceptance test that evaluated nothing has not accepted the step, it has failed to look at it —
  the same class of mistake as the coverage gate in ALEPH-PORT-303.
- **Method-name participant discovery.** Replaced by explicit registration (ALEPH-PORT-302/ALEPH-PORT-303).

## The gap this module closes — RNG rollback

The reference advances "biological time + RNG epoch" together on the accepted branch. What is not
visible in `transaction.py` or `world.py` is any restoration of stream *position* on the rejected
branch: `rollback(accepted_d)` is delegated wholly to each participant's own kernel. A rejected step
still consumed draws while proposing the events that got rejected. If those draws are not given back,
re-running from the same seed and the same state produces different numbers, and the divergence
surfaces much later as an unlocalisable "nondeterminism".

Aleph therefore treats the stream as a first-class transactional channel alongside state, clock and
topology, and snapshots the **bit generator's internal state** (a deep copy of the PCG64 state
mapping — plain integers, so equality is exact), not the seed.

## The four channels, and the rule

`propose()` captures `AmbientSnapshot(rng, clock, topology_generation)` **before** participants
snapshot and before any work runs. `commit()` moves state, advances the clock, and keeps the
candidate's topology generation and stream position. `rollback()` restores all four. There is no
public method that moves any one of them on its own. A commit that fails part-way restores the
ambient channels, marks the transaction `POISONED`, and re-raises — a transaction that cannot say
whether it committed must not be reusable.

## Controls

- Positive control: `test_accepted_step_advances_all_four_channels`.
- Negative control (must fail): `test_rejected_step_restores_rng_bit_identically`.
  In full, the deliberately failing set — the classic bug and its neighbours:
  - `test_rejected_step_restores_rng_bit_identically` — forces a rejection and asserts the PCG64
    state mapping is **equal**, not merely statistically similar, to the pre-proposal snapshot, and
    that the next draw after the rollback equals the draw that would have been made without the
    proposal at all.
  - `test_rejected_step_restores_clock_topology_and_state`.
  - `test_assert_ambient_restored_catches_a_leaked_stream` — proves the assertion has teeth by
    putting the stream one draw ahead of where it should be, and
    `test_a_participant_that_leaks_during_rollback_is_corrected`, which uses a deliberately broken
    participant that consumes the stream inside its own `rollback`.
  - `test_cannot_commit_a_rejected_decision`, `test_state_machine_refuses_out_of_order_use`,
  - `test_predicate_with_no_conditions_rejects`,
  - `test_failed_commit_poisons_the_transaction`.

## Provenance and envelope

- **Source repository / commit read:** `ffn_cellsim` at commit `be0e5876`, read-only.
- **Why not clean-room:** partly clean-room. No code was taken. The three ideas carried (one decision
  fanned out to all participants, identity-deduplicated participants, the predicate kept out of
  Python control flow) are structural findings restated from scratch; the RNG-as-a-transactional-
  channel design is original and closes a gap the reference leaves open.
- **Discarded prose:** all original comments and docstrings discarded.
- **Units, domains, invariants:** `dt` in s, finite and strictly positive. Invariants: the clock
  advances only inside `commit`; the ambient snapshot is taken before participants snapshot; a
  rejected candidate restores state, stream position, clock and topology generation together; a
  failed commit poisons the transaction rather than leaving it half-applied.
- **Numerical and precision envelope:** the RNG check is **exact**, not tolerance-based — the PCG64
  state mapping holds plain integers and is compared for equality. Clock arithmetic is float64
  accumulation of `dt`; step index is an exact integer, which is the reproducible key.
- **Residency:** host. No device predicate, no CUDA, no warp.

## Independence

No `ffn_cellsim` import, no `warp` import, no CUDA. Runs on the CPU reference backend. Every
identifier is Aleph's own.

# Notion Contract-Graph — schema audit against the inverse-inference reframe (2026-07-28)

**Status: PROPOSAL. No Notion row was written.** Schema changes to the 8-DB graph are gate-contract
changes and are PI-authored by charter; this document measures what the schema currently is, states
where it contradicts the 2026-07-25 reframe, and proposes the minimum delta. Everything below was read
from the committed mirror `outputs/tag_kb/kb.duckdb` (regenerated from Notion, never hand-edited).

## What is actually in the graph

| Table | Rows | Repo's own count of the same thing |
|---|---|---|
| `knowledge_claim` | 274 | — |
| `source_evidence` | 565 | — |
| `parameter` | **26** | ~70 KU-tagged config constants + 33 PI-GAP cards |
| `validation_gate` | 42 | 2 hash-stamped gate contracts; every other gate unstamped |
| `run_result` | 117 | — |
| `code_mapping` | 50 | **42 of the 50 paths do not exist** |
| `model_contract` | 12 | 14 components / 36 connectors |
| `decision_ledger` | 4 | D1–D8 alone are eight |

## Six findings

### 1. `Parameter` is modelled as an INPUT. The reframe makes it an OUTPUT.

Its fields are `default`, `range`, `calibrated`, `provenance`, `status`. That is the schema of a
forward model: you supply a default, and `calibrated` asks whether you tuned it (values in use: `No`
13, `Sweep` 9, `Yes` 2, `Maybe` 1, `Calibrate` 1).

Under the reframe a parameter is **inferred**, with uncertainty, **per cell type and per cell state**.
There is no field for a prior, a posterior, an uncertainty, which observable constrained it, or how
identifiable it was. `status` is `draft` for 25 of 26 rows, so the table is also close to unused.

### 2. Nothing in the graph represents a CELL TYPE or a CELL STATE.

This is the largest structural gap and it is easy to miss because nothing is *wrong* — the axis simply
does not exist. The deliverable is a per-cell-type, per-cell-state parameter library, and
`CLAUDE.md`'s stage 5 defines project progress as `cell types × states × parameters recovered`. That
denominator cannot be represented, let alone counted, in the current schema.

### 3. `ValidationGate` speaks the old gate's language and none of the new one's.

`band_expected` is free text naming a magnitude band (`"n in [1.9, 2.3]"`, `"K rises >10x above
gamma_c"`). The reframe moved the gate from **magnitude** to **mechanical connectedness**, and `type`
offers only `Validation` / `Verification` / `Sanity check` / `Phase-2 entry gate` — no category for a
connectedness or transmission gate, which is what most new gates are.

`status` has no `VOID`. PI decision **D7** created `VOID` precisely because a run whose residual swamps
its signal must not be able to report PASS — and the repo's `contracts.py` now has `GateVerdict`,
`EvidenceRung` and `QuantitativeClaim` as executable enums that the graph knows nothing about.

### 4. `RunResult` records neither population, rung, nor claim status — and 82% record no outcome.

`outcome` is NULL for 96 of 117 rows (`smoke` 17, `FAIL` 3, `partial` 1). More importantly there is no
`population` field, which is the single thing whose absence caused the 2026-07-28 retraction of four
headline numbers measured at 0.18% of native.

The repo already emits all of this: `observation_artifact()` writes `build`, `config_sha256`,
`census`, `t0`, `evidence`, `quantitative_claim_status`, `gate.verdict` and `timing`. The graph should
mirror that record rather than a hand-typed subset of it.

### 5. `CodeMapping` is 84% dead, and has been for a month.

42 of 50 paths do not exist. **41 of those 42 pre-date this week's cleanup** — they point at `bridge/`,
`cortex/`, `ecm/`, `integrator/`, i.e. the layout retired on 2026-06-29. Only `common/gsd_traj.py` is
from the HOOMD sweep. So the advertised traversal `run → gate → contract → parameter → source` breaks
on its first hop and has done since June.

The mechanism to repair it exists (`harvest_ops.py --apply`, idempotent on Path) and is PI-gated.

### 6. `KnowledgeClaim` stores a range, but as prose.

Good news first: `value_range_si` is populated for **all 274** rows, so "literature is a range, not a
point" is already true in intent — better than expected. But the content is free text
(`"gamma_eff = 0.27 mN/m (median; IQR 0.18-0.40; whiskers 0.05-0.62)"`), so it cannot be consumed as a
prior without a human reading it. The reframe says literature supplies **priors**; a prior a program
cannot read is not one.

## Proposed delta — smallest change that makes the reframe representable

**New database — `CellState`.** The missing axis. One row per (cell type, state): `MCF7 / suspended`,
`MCF7 / adherent-on-2D-collagen`, `MDA-MB-231 / migrating`. Fields: cell line, state, geometry
setpoints, and relations to the `RunResult`s measured in it. Without this, stage 5 has no denominator.

**`Parameter` gains an inference side**, keeping the existing fields for the forward path:
`prior` (structured, see below) · `posterior` · `uncertainty` · `cell_state` (relation) ·
`inferred_from` (relation to the RunResult / observable that constrained it) · `identifiability`
(the Fisher-information response, which `observe/sensitivity.py` already computes).

**`KnowledgeClaim.value_range_si` gains a machine-readable twin**: `prior_kind` (point / uniform /
lognormal / empirical) plus `prior_low` / `prior_high` / `prior_central`. The prose field stays as the
human record; the structured one is what an inference run reads.

**`ValidationGate` adopts the executable vocabulary already in `contracts.py`**: `verdict` ∈
`PASS`/`FAIL`/`VOID`, `evidence_rung` ∈ the ten-rung ladder, `quantitative_claim` ∈
`BLOCKED`/`OPEN`/`CONFIRMED`, and a `gate_kind` that admits connectedness/transmission alongside
magnitude. These are not new concepts — they are PI-ratified and executing in code; the graph is the
only place that has not been told.

**`RunResult` mirrors the run record**: `population`, `build_commit`, `config_sha256`, `evidence_rung`,
`quantitative_claim`, `wall_seconds`. All six already exist in every artifact `observation_artifact()`
writes, so this is a harvest change, not a new burden on runs.

## Order, and what blocks what

1. **`CodeMapping` repair** — mechanical, unblocks the whole traversal, needs only
   `harvest_ops.py --apply` plus PI sign-off on the upserts.
2. **`ValidationGate` vocabulary** — pure addition, no migration; makes the D1-B/D7 axes visible
   outside the code.
3. **`RunResult` mirror** — a harvest change; the fields already exist on disk.
4. **`CellState` + `Parameter` inference side** — the real schema change, and the one that should wait
   for stage 5's `K` and `M` to be decided, since those fix what a row of `CellState` even is.

Steps 1–3 are recoverable and additive. Step 4 defines the project's own progress metric and should be
PI-authored rather than proposed by a session.

## ADDENDUM 2026-07-28 (evening) — PI decisions, and a correction this audit owes its own step 1

**PI ratified two of the four steps.** Step 1 ("repair CodeMapping") → *제대로 수리*, i.e. fix it
properly rather than run the existing mechanism. Step 4 (`CellState` + the `Parameter` inference side)
→ do both now, since K=2 × M=2 fixed what a `CellState` row is. **Steps 2 and 3 were not put to the PI
and remain open** — `ValidationGate` still has no `VOID` and no connectedness `gate_kind`, and
`RunResult` still has no `population` field, which is the single absence that caused the 2026-07-28
retraction of four headline numbers.

### Correction: step 1 was NOT mechanical, and `--apply` could not have done it

This document said step 1 "needs only `harvest_ops.py --apply` plus PI sign-off on the upserts." That
is wrong, in two independent ways that only surface when the thing is actually run.

**(a) The harvester was pointed at the retired layout.** `PKG_DIRS` listed
`cortex, cell, ecm, bridge, junction, integrator, native, common`. Seven of those eight were retired on
2026-06-29. Only `common/` still existed, so the scanner found **9 modules, all legacy**, and scanned
**zero** files in `ac/` — the canonical engine, 174 modules. That is also the origin of the 42 dead
rows: they are this scanner's own past output. Re-running a scanner aimed at the June tree re-affirms
the June tree. Repointed to `ac`/`ff`/`common`, it now sees **219 modules**.

**(b) `upsert_code` had no way to remove anything.** It only created and patched, so a moved or deleted
module left its row pointing into nothing permanently. Worse, `_drift_check` computed only
`disk − graph`; nobody ever computed `graph − disk`. That is how 41 dead rows sat for a month behind a
drift line reading "1 un-harvested" that looked healthy.

Both are fixed in `32ac72ae` (repoint, archive pass for machine rows, curated stale rows reported but
never touched, and the drift check now reports both directions).

### A second defect the repair exposed: KU matching cannot tell *implements* from *cites*

Only visible after the repoint, and more dangerous than the dead rows. `match_contracts` linked a
module to every ModelContract carrying any KU its docstring head mentioned. Measured on the repaired
scanner, **2 of the 6 KU-derived edges were plainly wrong**:

| Module | What it actually is | Contract it was linked to | Why |
|---|---|---|---|
| `common/checkpoint.py` | checkpoint I/O | `MC-H3-composite-tension` | docstring says checkpoints are used *during* the KU-3.5 sweep |
| `aleph/engine/sf_mechanics.py` | stress-fibre mechanics | `MC-U1-ecm-network` | cites KU-1.1 as the *source* of a persistence length |

Both would have been written into the SoT by `--apply`. An explicit `Implements: KU-x.y` line now
creates the edge; a bare mention becomes a review suggestion in the manifest and nothing else. **A
missing edge is a visible leaf; a wrong edge makes the traversal answer confidently about the wrong
contract.** This supersedes CLAUDE.md's older "put the relevant KU in a module's docstring head"
convention, which is too weak to carry an edge.

### The traversal is still broken, and repointing cannot fix the rest

**183 of 219 modules link to nothing** — they carry no `Implements:` line and match no package-path
fallback. The 36 that do link are all `ac/cell`, via a single deliberate `DIR_CONTRACT` entry (a
rename of the old `cell` key, not a judgement). Every other directory was left unmapped ON PURPOSE:
deciding that `ac/motor` implements `MC-U2-fa-motor-clutch` — an *adhesion* clutch contract, not an
NMII one — would be authoring a contract relation, which CLAUDE.md reserves for the PI.

So the honest state is: the graph has stopped lying, and the remaining gap is **authorship, not
mechanism**. Tagging `ac/engine` is a per-module decision, not a batch operation.

## What this audit does not claim

That the read layers are unhealthy. `tag_query.py` runs at temperature 0 with a BM25 relevance floor,
`paper_chunks` holds 12,810 chunks over 278 PDFs, and the citation audit covers all 565 sources. The
query engine is fine. What is stale is the **contract side** of the graph — the half populated from
disk by `harvest_ops` — and what is missing is the axis the project reorganised itself around.

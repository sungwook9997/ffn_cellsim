# The merge and the rename — fourteen structural decisions

> **Status: AGENT-PROPOSED as a document.** The decisions below were taken by the PI in session
> `5ccfd28d`. This file **cites** them; it does not assert them. No `decided_by` field appears here
> and an agent may not add one. The checkable form of every citation is the transcript path plus the
> timestamp — open it and the quote is either there or it is not.

| Field | Value |
|---|---|
| Written | 2026-08-09 16:3x KST, session `5ccfd28d-4fb2-48da-acb6-23e3ce818193` |
| Transcript | `~/.claude/projects/-Users-sw1/5ccfd28d-4fb2-48da-acb6-23e3ce818193.jsonl` |
| Direction | **`ffn_cellsim` becomes the base and is renamed Project Aleph.** Assets travel Aleph → here |
| Companion | `STRUCTURE.md` (the target tree) · `Project_Aleph/docs/design/PORTING_PLAN_FOR_THE_MERGE.md` (order of operations, written by session `b4fad06b`) |
| Supersedes | `FFN_PROBABILISTIC_VIRTUAL_CELL_MASTER_PLAN.md`'s status line *"PI 검토용 재창설 제안"* — the identity claim it was waiting on is the rename |

---

## 0. Why this document exists rather than fourteen

Aleph's convention is one file per decision. Fourteen files for one session's structural decisions
would be boilerplate, and the decisions are not independent — A4 changed because of C, and A1's
timing changed because of a fact discovered in another session's document. They are one decision
about one tree, taken in one sitting, and they are recorded that way. Individual items get their own
file only when one of them is reopened.

## 1. What was decided, and what each one costs

### A — structure and naming

| # | Decision | Consequence |
|---|---|---|
| **A1** | `ffn_cellsim` is renamed **Project Aleph**; work continues in this tree. Package `ffn_sim` → `aleph` | ~896 files touched by one mechanical rename. Done as a single commit with no other change |
| **A2** | **Both** splits, not one. Components and connectors are organised by biology (`components/`); shared laws stay in one layer (`laws/`) | Measured: `ff/network_warp` is imported by **17** files in `ac/`, `membrane_surface` and `forces_warp` by 10 each. Dissolving `ff/` into per-component directories would force 16 cross-component imports and reproduce Aleph's `runtime → vertical` cycle exactly. But only **14 of `ff/`'s 44 modules** are used by `ac/` — the other 30 are `viz_*` (8), `gamma_floor*` (3), `cytosim_parity`, `architecture_metrics` and friends, i.e. visualisation and diagnostics living in a law library, which is the same defect Aleph's `runtime/law_cases.py` has |
| **A3** | **Revised — see §1c.** `ac/cell/` **stays**, relocated to `components/incumbent/`. `dcm/` and `ac/lane_d/` archive as planned | Archiving `ac/cell` was not executable: it is imported **243** times, and `ac/engine` — the canonical layer — is one of the importers |
| **A4** | **Revised — see §1a.** Aleph's `represent/` + `learn/` (12,400) travel as **documents**. ffn's `virtual_cell/` (29,097) **stays in the tree** and is triaged per module, not archived | The blanket-archive version of A4 rested on a claim that measurement refuted. §1a is the measurement and the revision |
| **A5** | **Revised — see §1e.** `ac/lane_d/` archived. **`dcm/` could not be**: `laws/relax.py` calls its implicit step | `lane_d`'s 15 imports are all from its own tests. `dcm`'s 178 were all internal — which says dcm does not look out, not that nothing looks in |
| **A6** | Scripts selected by **gate ∪ tier-(a) producer ∪ 1-hop dependency** | Measured: **44 kept / 15,746 lines**, 201 archived / 34,574 lines — a 69 % cut. The stated criterion "current production import closure" had no recorded manifest, so its 1-hop expansion stands in for it |
| **A7** | `common/` is dissolved; its three load-bearing modules move to `laws/` | `surface_manifold`, `turgor_pi0` (used directly at `ac/cell/assemble.py:99`), `compartments` (via `ff/membrane_surface.py`). The rest is dead and `common/integrity.py` scans subpackages that no longer exist |

### 1a. A4 revised, because the measurement refuted the claim it rested on

A4 was first written as *"archive both probabilistic layers"*, and the argument for archiving ffn's
`virtual_cell/` was that it is **30k lines that compute but that nothing validates** — unwired to the
engine, therefore unadjudicated. Half of that is true and half was never measured.

Measured 2026-08-09, per module, by which test files import it:

| | |
|---|---|
| modules in `virtual_cell/` | 26 |
| **modules with tests** | **26 of 26** |
| module lines | 29,097 |
| lines of tests importing them | ~60,000 |
| worst-covered | `tensor_train` — 263 lines, 99 lines of tests |
| best-covered | `contracts` — 523 lines, **10,370** lines of tests |

**Unwired to the engine and unvalidated are different claims, and only the first one holds.**
`virtual_cell/` is imported by no module under `ac/`, `ff/` or `scripts/` — that measurement stands.
But every one of its modules is exercised, several at four to ten times their own size.

The concrete case that forced the recheck: `virtual_cell/observation_operator.py` (1,285 lines,
884 lines of tests) was about to be archived and then re-ported from Aleph's `observe/operator.py`
(493 lines) under `ALEPH-PORT-4002`. Reading both, **ffn's is the stronger of the two** on the axis
they share — it computes an applicability verdict *before* any arithmetic, a refusal is typed and
provably carries no number, and unimplemented operators are named with a capability verdict rather
than stubbed. Aleph's contributes one idea ffn has no spelling for: the **apparent-vs-true** split,
where a reported modulus carries the `analysis_chain` that produced it. `analysis_chain` has **0
hits** in this tree.

So:

- **`virtual_cell/` stays**, and becomes the material for `inner/` and `outer/` (decision B5) rather
  than an archive. Triage is per module, at the stage that needs it.
- **Aleph's `represent/` + `learn/` still travel as documents.** That half of A4 is unchanged: they
  are computation-sealed by ratification and contain no code to adjudicate.
- **`ALEPH-PORT-4002` narrows** from "carry the observation operator" to "carry `analysis_chain` and
  the apparent-quantity type onto the operator this tree already has".

This is the third time in one session that a port candidate turned out to exist here already and to
be the better of the two — after Aleph's `evidence/ladder.py` (this tree has three axes to its two,
including `GateVerdict.VOID`) and its `units/guard.py` (`ALEPH-PORT-4001`, downgraded to
`CONCEPT_ONLY`). **The pattern is worth naming: this tree's weaknesses are documented and its
strengths are not**, so a reader who trusts the documents — which is what the first pass of this
merge did — will systematically undervalue it.

### 1b. A2 narrowed, for the same reason A4 was

A2's consequence column said only **14 of `ff/`'s 44 modules** are used by `ac/`, and read that as a
licence to archive the rest. Two corrections, both from measuring rather than reading:

**The count.** By AST — imports, not text — it is **14 of 49**. The "20" that appeared earlier was a
grep, and grep counts a module named in a docstring. That is the third time in this session that a
text scan inflated a structural number; the first two are in §1a.

**The inference.** *Not imported by `ac/` today* is not *archivable*. The 27 modules the count would
have discarded include:

| module | what it actually is |
|---|---|
| `implicit_ff` | the NF2007 semi-implicit step. `Project_Aleph` ported **this file** as `ALEPH-PORT-3664` and measured **708×** on the sourced cortex |
| `relax` | 104 lines, which that lane's own porting plan calls *"already the right shape"* |
| `cytosim_parity` | the **external** parity oracle — the one thing in the tree that is not judged by the tree |
| `wlc`, `polymerization_warp`, `myosin_linear`, `motility_warp`, `fa_clutch_warp` | physics laws waiting for the component that binds them |

Archiving those would have thrown away the solver the merge most needs, on the grounds that nothing
imports it yet — which is true of every part of a system before it is wired.

**So A2's move is simpler than first written**, and simpler is also what it should have been:

- `ff/` → `laws/` **wholesale**, minus the two kinds of thing that are not laws;
- the eight `viz_*` modules → `viz/`, joining Aleph's rendering package;
- `cytosim_parity` → `validation/`, where an oracle belongs and where `laws/` cannot import it;
- **archived: `gamma_floor`, `gamma_floor_dynamic`, `gamma_floor_sweep`** only — the path whose three
  live defects `ff/ENGINE.md` names in its own banner and whose numbers `STATE.md` (c) 2 already
  blocks. Four modules of 49, not 27.

### 1c. A3 revised: the incumbent cannot be archived because the engine still needs it

A3 read `STATE.md` (a) — *"feature-frozen … being strangled under the engine in gated stages"* — and
concluded that `ac/cell` could leave the tree because its numbers live in artifacts. The second half
of that sentence is the part that matters: **the strangling is not finished.**

Measured by AST, imports of each archive candidate:

| candidate | imports | importers |
|---|---:|---|
| `dcm/` | 178 | **all inside `dcm/`** — self-contained, archives cleanly |
| `ac/lane_d/` | 15 | all from `tests/ac` — archives with its tests |
| **`ac/cell/`** | **243** | `ac/cell`, **`ac/engine`**, `ac/motor`, `ac/lane_d`, and a long tail of drivers |

`ac/engine` importing `ac/cell` is the whole reason `ac/cell` is described as an incumbent rather
than as dead code. Archiving it would break the canonical layer on the first import.

**So `ac/cell` relocates rather than leaves**: `components/incumbent/`. That is honest about what it
is — a component-level whole-cell assembly the engine still binds — and it keeps the layer direction
legal, because `engine → components` is the permitted edge.

The archive shrinks to `dcm/` (15,995) and `ac/lane_d/` (2,001), both of which the measurement shows
can go without breaking anything.

**Three of the seven A-decisions have now been revised by measurement, always in the same
direction — the plan wanted to throw away more than the tree could spare.** A2 (27 modules → 4),
A4 (29,097 lines → 0), A3 (17,061 lines → relocate). Each was written from a document that
described the code as finished, frozen, or unvalidated, and in each case the code disagreed.

### 1d. What the rename actually cost, recorded because the next one will be tempted by `sed`

A1 was described as *"one mechanical rename"*. It was mechanical; it was not one pass. `ffn_sim`
turns out to name **four** different things, and to appear in **four** syntactic positions, and a
`sed 's/ffn_sim/aleph/g'` gets both classifications wrong in opposite directions.

**The four things, of which only one moved:**

| | | rewritten? |
|---|---|---|
| the Python package | 1,057 files | **yes** |
| the **conda environment** `~/miniconda3/envs/ffn_sim` | 130 files | no — a directory on disk this change does not move |
| the GitHub remote | 1 | no — renaming a repository is a server-side action |
| the string inside committed run artifacts | many | no — an artifact records what ran |

**The four positions, of which three were missed on the first pass:**

| position | example | how it surfaced |
|---|---|---|
| dotted / slashed | `aleph.engine`, `aleph/scripts/x.py` | handled from the start |
| bare after `import` / `from` | `import ffn_sim` | handled after a control caught that prose was being rewritten too |
| **extensionless files** | `Makefile`, `.githooks/pre-commit` | a suffix allowlist skipped them; 14 stale paths would have broken **every commit** and `make kb-check` |
| **a quoted standalone token** | `ROOT / "ffn_sim" / "docs"`, `_CANONICAL_ROOT = "ffn_sim"` | **22 tests red.** A path built component-wise contains neither a dot nor a slash |
| **live config under a skipped directory** | `outputs/tag_kb/*_manifest.yaml` | `outputs/` was skipped to protect run records; the KB manifests live there too, and both drift gates failed with 66 "constant ABSENT on disk" |

**And one failure mode with no position at all: the script rewrote itself.** Its own docstring became
`sed -i 's/aleph/aleph/g'`, and every row of its control table became
`("from aleph.ac import x", "from aleph.ac import x")` — a test asserting a no-op, passing. **A
self-referential rewrite disarms its own tests rather than breaking them**, which is why it is the
one on this list that would have survived a green suite. `EXEMPT_FILES` now holds the script and its
controls, and the single import line inside them is maintained by hand.

Every one of these is in `tests/scripts/test_rename_package_to_aleph.py` as a case, and the layer
restructure that follows was written against that list rather than rediscovering it.

### 1e. A5 revised, and what the layer guard caught the moment the directories moved

`dcm/` was archived on the strength of a clean-looking measurement — **178 imports, every one of them
inside `dcm/`**. That number says dcm does not look outward. It says nothing about who looks in, and
the layer test answered that within a minute of the move:

```
`laws/` reached upward:
  aleph/laws/relax.py            imports aleph.dcm.dcm_warp_implicit
  aleph/laws/ecm_library.py      imports aleph.laws.architecture_metrics
  aleph/laws/viz_cell.py         imports aleph.laws.gamma_floor      (and three more viz_*)
```

`relax.py` is the 104-line descent the other lane called *"already the right shape"*, and the shape
is a wrapper around **dcm's** `implicit_overdamped_step`. **A module a live law calls is not
archivable, whatever its own imports look like.** Three reversals followed:

| | |
|---|---|
| `dcm/` | out of `archive/`, back to `aleph/dcm/`. It is a solver library, and `laws → dcm` is a legal edge |
| `architecture_metrics` | back to `laws/` — `ecm_library` calls its `parallel_order_parameter` |
| the eight `viz_*` | out of `laws/` into **`viz/`**, which is where they always belonged. Four of them import `archive.gamma_floor`, which is legal one layer up and was a violation one layer down |

**The guard paid for itself on the first move it was written for.** Its whole argument was that
Aleph's 53-edge cycle arrived one good reason at a time and a human reviewer sees one import per
diff; here it saw seven at once, in a restructure whose author had just written that argument down.

**And two self-referential failures, one per script.** The rename rewrote its own control table into
`("from aleph.ac import x", "from aleph.ac import x")` — a test asserting a no-op, passing. Then the
restructure rewrote the *rename's* control table, so those cases asserted the restructure's output
against the rename's function. Each script now exempts both sets of controls. **A rewriting tool
whose own tests are rewritable does not have tests**, and the first of the two would have survived a
green suite.

**One more, with no position and no self-reference: a module move must carry its sibling data.**
`turgor_pi0.py` reads `Path(__file__).with_name("params_turgor.yaml")`. Moving the module out of
`common/` without the YAML broke collection in fourteen test files. The whole-directory moves were
safe — `git mv` takes the data with it — so this was only ever a risk for the three files extracted
individually, and exactly one of them had a sibling.

### 1f. A reverted decision has to be reverted in the tool, not only in the tree

`dcm/` was archived, the layer test refuted it (§1e), and it was moved back by hand. The mapping
`("dcm", "archive.dcm")` stayed in `restructure_to_layers.py`'s `MOVES` table, so **the next run
archived it again** — and because `common/sim_realtime.py` was routed to `dcm/` in the same run, a
fresh one-file `aleph/dcm/` appeared beside the 52-file `aleph/archive/dcm/`. The package was in two
places at once and nothing failed, because every import had been rewritten to match.

The lesson is small and expensive: **a script that encodes a decision is part of that decision.**
Fixing the tree and leaving the table is a fix with a timer on it, and the failure mode is silent —
the second run produces a consistent, wrong tree rather than an error.

Both places are now consistent, and the removed row carries the reason in a comment so the next
reader does not helpfully add it back.

### 1g. The session's own false green, and it is the same defect as everything above

Deleting `common/` left `"common"` in `conftest.py`'s `_ALIASED_PACKAGES`, so `conftest` raised
`ModuleNotFoundError` and **pytest collected nothing**. The run that was supposed to confirm the
restructure was invoked as:

```
pytest -q --tb=line 2>&1 | grep -E "^FAILED|^ERROR|Interrupted"; echo "PYTEST_DONE"
```

`ImportError while loading conftest` matches none of those three patterns. The command printed
`PYTEST_DONE` and nothing else, which is exactly what a clean run prints, and it was reported as
green.

**The filter only matched the failure shapes its author had already thought of** — which is the
same defect as the token-grep that counted a filename as a test definition, the import-count that
called a live solver archivable, and the rewriter that turned its own controls into identities. Four
different tools, one assumption: *that a text pattern can stand in for the structure underneath.*

The fix is not a better pattern. It is to stop asking the output whether it passed and ask the
process: capture `$?`. A green claim in this session is only worth reading if it came with an exit
code, and the ones before this point that did are the rename's (`exit=0`, captured explicitly);
the ones that did not have been re-run.

### 1h. The archive shrank a sixth time, and the reason is worth stating once more

`gamma_floor`, `gamma_floor_dynamic` and `gamma_floor_sweep` were archived because `STATE.md` (c) 2
blocks every number they produce and `ff/ENGINE.md`'s banner names three live defects in that path.
Both statements are true, and the conclusion did not follow: the modules carry **live guard tests**,
one of which is literally
`test_three_2026_07_23_defects_are_still_live_in_the_gamma_floor_builder`.

**Archiving the module deletes the test that records the defect** — the opposite of retiring a claim.
A blocked result and a dead module are different things, and (c) is a list of the first kind.

Final tally: the archive was drafted at roughly **80,000 lines** and landed at **2,975** —
`lane_d/` plus four `common/` modules that nothing imports. Everything else that was proposed for it
turned out to be reachable, tested, or both. **Six reversals, one direction.**

`architecture_metrics` came back for the same reason as `dcm` and by the same mistake: reverted in
the tree, left in the script's table, re-archived on the next run. Twice in one session is a pattern,
and §1f is where it is written down.

### 1i. A6 executed, and the selection criterion was wrong three times before it held

**190 of 248 scripts archived, 58 kept.** The criterion as written in A6 — *gate ∪ tier-(a) producer
∪ 1-hop* — under-selected in three distinct ways, each found by running it:

| miss | what it cost |
|---|---|
| `gate` matched on **filename** | `check_state_md.py`, `check_ownership.py`, `render_state_rows.py` are the gates the **pre-commit hook itself runs**, and none has "gate" in its name. Archiving them breaks every commit |
| the doc surface was not scanned | `ffn_dashboard.py` is named in the `Makefile`; nothing else pointed at it |
| the import scan missed a spelling | `from aleph.scripts import X` — as opposed to `from aleph.scripts.X import …`. Two test modules use it, and both went red |

The criterion that actually holds is **four** clauses, and the third and fourth are the ones a reader
would not think to write down: a script is kept if it is a gate, produced a `STATE.md` (b) row, **is
named anywhere in the hook / `Makefile` / CI / the charter**, **is imported by a test in any
spelling**, or is one hop from something that is.

This is the seventh time in one session that a selection heuristic was refuted by running it, and
the shape is the same every time: the heuristic encoded the reason the author had in mind, and the
tree had other reasons. Archiving is a `git mv`, so each was recoverable — which is the only thing
that made seven affordable.

### 1j. The inherited ledger travels, but not into the active one — and there are 135, not 136

`PORTING_PLAN_FOR_THE_MERGE.md` §7.2 asked whether Aleph's port records should come across. They
should: they are the argument for why its contracts exist, and without them the rules look
arbitrary once the tree that produced them is an archive. Copied into `ports/ledger/`, all 135 fail
this tree's discipline test on six counts at once — and **every failure is a category error**:

- they name controls under `Project_Aleph/tests/`. The tests exist; they are not here.
- they cite source commits in that tree's template, not this one's.
- their status vocabulary is close to this one's and not identical.

A record of a port **into another repository** is not a port into this one. Scoring it against this
tree's discipline would say something false about the record, so they live in `ports/inherited/`,
immutable, with their own rules and their own README. Two tests hold the separation: the set is
present and complete, and no inherited ID has leaked into the active band.

**And the count was wrong from the start.** Every statement of it — B5's, and the porting plan's own
*"136 entries"* — came from `ls ports/ledger/*.md`, which counts `INDEX.md`. The index is not a port
record. **135.** A one-off miscount, but it is the same shape as the four text-versus-structure
errors in §1a and §1d: a command that answers a question adjacent to the one being asked.

### B — process and authority

| # | Decision | Consequence |
|---|---|---|
| **B1** | **Both** session mechanisms | They solve different problems and neither substitutes for the other. `ownership.yaml` + the pre-commit hook **refuses a commit** outside a declared path — a machine, acting after the fact. `ACTIVE_SESSIONS.md` **announces** a session before its first write — a human channel, acting before. Aleph's 2026-07-30 incident cost ninety minutes precisely because the second existed nowhere |
| **B2** | `ffn_gpu.py`'s lease and `ffn_gpu_queue.py` are **deleted** | Slurm replaces both. The lease assumed a dedicated host, expired on wall-clock rather than on work — which cuts unattended runs — and defaulted `--host gbook`, a machine the PI retired 2026-08-04 |
| **B3** | Aleph's decision-record format | `docs/decisions/`, `PROPOSAL-*` carries **no** `decided_by` field, `ALEPH-PD-*` is PI-only and an agent may never create one |
| **B4** | **KB stays ffn's, wholesale. Evidence ledger becomes Aleph's code** | These are different stores and the first framing of B4 conflated them. The **KB** — Notion contract-graph + DuckDB TAG + Obsidian vault — has no Aleph counterpart and is carried unchanged. The **evidence ledger** is the claim that a result may be quoted; today that is a hand-maintained markdown table, and `STATE.md` says of its own blocklist that *"(c) 14 and (c) 15 were both caught by hand"* |
| **B5** | The tree opens `inner/` and `outer/` as first-class packages; port IDs continue in a new band | The **135** existing Aleph ledger entries **travel** — to `ports/inherited/`, not the active ledger; see §1j. Answers `PORTING_PLAN_FOR_THE_MERGE.md` §7.2, whose own "136" counted `INDEX.md` |

### C — the definition of "the engine is finished"

The PI's words, 2026-08-09 15:5x: *"엔진 완성의 정의는 파라미터 스윕을 위한 실험이 정해지고 그 실험을
시뮬레이션 상으로 구동시킬 수 있게, 그리고 이를 위해서는 수렴이나 다른 물리적 엄밀성 맞아야지"*.

As a gate, in three terms:

| | Term | Where it stands |
|---|---|---|
| **C-1** | The **one experiment** the sweep targets is declared in advance — name, protocol, observable, uncertainty | Open. This is `ALEPH-DQ-102` (first observation modality); the standing candidate is the passive membrane fluctuation spectrum |
| **C-2** | That experiment runs to completion at the **full native population** | `inner_converged: false`; `overall_status: INCOMPLETE_OR_FAIL` |
| **C-3** | The convergence predicate **can fail for non-convergence** | It cannot today. Stage 1's inner gate accepts on `balance_ok_d`, and `sf_motor_slice.py:65` says in its own docstring that the predicate tests adjoint wiring, not convergence |

> **⚠ AMENDED 2026-08-09 22:xx — C-1 no longer comes first.** The paragraph below stood until a PI
> instruction in session `97fb3916`'s transcript reordered it. **The correction is cited rather than
> asserted, and this session verified it by opening that transcript rather than accepting the
> report** — the uuids below are not exposed to the session that relayed them.
>
> | | |
> |---|---|
> | transcript | `~/.claude/projects/-Users-sw1-ffn-cellsim/97fb3916-7a35-4be7-b029-d6ad2485617b.jsonl` |
> | 2026-08-09 19:51:07 KST · uuid `3cb71512-b4ce-432e-abb6-6231a956a5d7` | *"컴프즈드 월드를 일단 gpu only로 돌릴 수 있게 전부 다 그러면 언어 바꿔줄래? 그게 먼저 잖아 c-1이 먼저 아니라"* |
> | 2026-08-09 20:29:37 KST · uuid `7a5c4bde-1732-4783-b40c-4155a64735f8` | *"파라미터 스윕이니깐 그냥 후보들로 넣고, 나중에 스윕 데이터로 가질거니깐 정확한 값을 가진다는 것은 없잖아"* |
>
> **The new order is C-2 before C-1**: get the composed world running GPU-only across all 14
> components and 36 connectors first. The argument, and it holds: **C-1 sets a tolerance and a
> duration, not whether the solver can drive the residual down.** And the tolerance in the code today
> — `sqrt(eps_f64) × 0.5 µm`, `components/incumbent/driver.py:567` — never came from an experiment
> either, so waiting for one to name it was never what was blocking convergence.
>
> **The second quotation reframes every PI-GAP slot and is the larger of the two.** A gap is not a
> blocker waiting for a point value; it is **a prior waiting to be swept**. The charter already said
> parameters are OUTPUTS and literature supplies priors — the slots had drifted into being treated as
> inputs, and 33 PI-GAP cards were being read as 33 things somebody had to answer before the sweep
> could start. They are the sweep's axes.

**C-1 came first under the original ordering.** Until the experiment is named, "runs to completion"
has no definition to test — which is true, and is why the amendment above replaces the ordering
rather than deleting the term. C-1 still has to happen; it is no longer what happens next.

This gate subsumes three items that were being tracked separately and are the same question:
`STATE.md` (e) 1 (the held gate amendment), `ALEPH-DQ-104` (the acceptance predicate), and
`PROPOSAL-what-does-converged-mean-for-a-cell-with-motors-in-it`.

### D — the inherited backlog

Aleph carries 35 `PROPOSAL-*` plus `ALEPH-DQ-101`…`-107`; ffn carries `STATE.md` (e) 1–5. About 47
items, and most of them are questions belonging to stages not yet reached — *is the nuclear envelope
impenetrable*, *Darcy or Stokes for the cytosol*. **All are relabelled `INHERITED-` and frozen.**
Only items blocking C are promoted.

## 2. The programme these decisions serve

The PI's words, 2026-08-09 15:56: *"목표는 알레프가 맞음. 엔진 완성 후 파라미터 스윕 내부 신경망 /
외부 연구 결과 외부 신경망 구축 및 통합으로"*.

| Stage | What | Blocked on |
|---|---|---|
| 0 | Merge and rename | — |
| 1 | **Engine finished = C-1 ∧ C-2 ∧ C-3** | C-1 |
| 2 | The four multipliers: `dt` under a physical predicate · Fisher dimension reduction · the analytic gradient · the surrogate | 1 |
| 3 | The sweep → **Inner Library** (admission: converged **and** carrying an accepted digest) | 2 |
| 4 | Observation Compiler → MTG-PN → **Outer Library** | **nothing — runs in parallel from now** |
| 5 | Integration → `CellStatePosterior` | 3 ∧ 4 |

**Two corrections to the stated ordering, both from measurement rather than preference.**

1. **The inner network is not downstream of the sweep; it is one of the four things that make the
   sweep possible.** `ROADMAP.md`'s arithmetic: an 8 × 5 grid is 14.9 years before the ×4 for cell
   states. Without a surrogate and sequential design there is no sweep to be downstream of.
2. **The outer network does not wait for the engine.** Its inputs are literature and experiment, and
   the KB those come from already exists here. Stage 4 starts now.

**And one warning against the arithmetic above.** `PORTING_PLAN_FOR_THE_MERGE.md` §5.5 records that a
`steps^2.08` extrapolation in a planning document measured **0.97** on the GPU, and §0 records that
`ALEPH-PORT-3665`'s IMEX driver already returns **28.3×** whole-cell at 10× `dt`. Every budget figure
in `ROADMAP.md` that rests on an extrapolated exponent is **re-derivable, not citable**.

## 3. What this document does not decide

1. **The demotion of the `Project_Aleph` tree.** A1 makes it an archive. `b4fad06b` is open in it as
   of 16:0x today. Renaming a tree with a live session in it is the incident
   `ACTIVE_SESSIONS.md` records at 2026-07-30 16:14. **Deferred until that row closes.**
2. **`aleph/vertical/relax.py`'s deletion.** `PLAN.md` item #2 is the PI's, per
   `PORTING_PLAN_FOR_THE_MERGE.md` §5.6. Not touched.
3. **`ff/ENGINE.md`'s three live defects.** Its own banner withdraws every full-native headline on
   the page, and names `n_xl` silently running at 1.0 because the builder never reads
   `ff/architecture_spec.CORTEX`. Scoped to the γ-floor production path, which `STATE.md` (c) 2
   already blocks — but the modules moving to `laws/` must be checked against it before any of them
   is cited again.
4. **`crosslinker.py:139`'s existing import of `aleph.laws`**, counted by a guard owned by the still
   open `46143f30` row. If the merge changes that import surface, the guard is **reported, not
   edited**.

# ffn_cellsim — Claude Code Project Context

## Project in one line

`ffn_cellsim` is a fine-grained, mechanistic HOOMD-blue simulator of single-cell mechanobiology where every cytoskeletal filament, motor head, adhesion clutch, and ECM cross-link is an explicit particle/bond. No paper-model wrappers, no lumped mechanisms.

## Architectural principle (PI 2026-05-19, hard rule)

For every design decision in the initial framework, **pick the full-fidelity, fine-grained, mechanistic option over abstracted, lumped, or proxy mechanisms** — even at higher implementation cost or worse wall-time.

The v1 codebase (`~/ActiveCellSim`, kept as origin remote) borrowed published paper models (Chan-Odde, Pereverzev, Bell-Evans, Buckley, Hill) as the **runtime mechanism**. `ffn_cellsim` inverts this: the runtime is HOOMD fine-grained particle/bond dynamics, and the paper closed-forms are **acceptance oracles** used only in validation tests.

Worked examples of this rule:

| Decision | Wrong (abstracted) | Right (mechanistic) |
|---|---|---|
| Bond off-rate | Metropolis (detailed-balance proxy) | **Bell-Evans** force-dependent |
| Motor force-velocity | Linear-stall convenience | **Hill** (Hill 1938) |
| Myosin II | Single 2-head spring (AFINES) | **Stam-Hocky multi-head bipolar minifilament** |
| Cadherin | Slip-only KU-4.17 | **Full catch-bond KU-4.2 from Phase 1** |
| Arp2/3 branch | Rigid 72° | **Angle-harmonic with thermal fluctuation** |
| Integrator | Euler-Maruyama | **Leimkuhler-Matthews BAOAB-limit plugin** |
| Excluded volume | EV-off | **LJ repulsive on from Phase 1** |

The **only** sanctioned coarse-graining is the ×40 mesoscopic filament scale (~1000 effective filaments per cell instead of ~38,000 native) — explicitly ratified in Plan v2 §3 H.3 (2026-05-19 v3.1) as a hardware constraint.

When in doubt, write down the option's Plan reference, the abstraction it introduces, and ask: *does this replace a mechanistic process with a lumped one?* If yes, prefer the mechanistic alternative.

## Stack

- **Python**: 3.13.13
- **Simulation core**: HOOMD-blue 7.0.1 (conda-forge)
- **Integrator**: custom Leimkuhler-Matthews BAOAB-limit plugin (Phase 0.4)
- **Active filament model**: AFINES (Simfreed/AFINES canonical) re-implemented on HOOMD
- **Env**: `conda activate ffn_sim`
- **Sanity bench**: `python ffn_sim/scripts/hoomd_polymer_sanity.py --steps 50000 --bench`
- **Hardware baseline (PI-ratified 2026-05-31, HARD)**: **CUDA GPU — RTX A5000 (gbook
  laptop) or better is the MANDATORY minimum** for production. CPU (Apple silicon M1 Max;
  1 M-step polymer in ~11 s) is now **dev/fallback only**, not a production target. The
  Phase 0/1/2 distinction is **RETIRED** — "GPU is Phase 2" no longer applies; GPU is now.
- **Code direction: GPU-main (PI 2026-05-31).** New/ported code must run GPU-resident by
  default — the custom BAOAB integrator + binding updaters (myosin/xlink) are being ported
  off per-step `cpu_local_snapshot` (which forces a GPU→CPU sync every step) to
  `gpu_local_snapshot`/cupy (or compiled CUDA plugins). No hard-coded device IDs; CPU stays
  selectable for dev. Tracking: `ffn_sim/docs/v2_audit/GPU_MAIN_PORT_*.md`.

## Repo layout (authoritative map: `STRUCTURE.md`)

```
ffn_cellsim/
├── ffn_sim/                # active HOOMD runtime + docs + tests + scripts
│   ├── ecm/  cell/  cortex/  bridge/  junction/  integrator/  common/
│   ├── validation/oracles/ # v1 closed-form oracles, runtime-import-forbidden
│   ├── docs/               # PHASE_0_*, AFINES_ALGORITHM_NOTES, briefs/, v2_audit/
│   ├── tests/  scripts/  outputs/
│   └── __init__.py
├── README.md  CLAUDE.md  STRUCTURE.md
├── pyproject.toml  requirements.txt  requirements-dev.txt
└── .gitignore  .stignore
```

## Knowledge base — RAG + TAG + Obsidian (read this before answering KB questions)

The project's literature/decision knowledge lives in a **Notion Contract-Graph**
(8 atomic, bidirectionally-related databases) that is the **single source of
truth (SoT)**. Two read layers are materialised *from* Notion — never edited
directly — plus a query engine. All three live under `ffn_sim/outputs/`.

```
Notion 8-DB Contract-Graph  (SoT — edit here only)
  SourceEvidence · KnowledgeClaim(KB-x.y) · ModelContract · Parameter
  · ValidationGate · CodeMapping · RunResult · DecisionLedger
        │  query_all() (Notion API pull)
        ├──────────────► outputs/obsidian_rag_full/   [GRAPH layer]
        │                  ~640-node Obsidian vault: typed wikilinks, color
        │                  groups, code/test/doc nodes, Dev-Logs day-logs.
        │                  Use for: visual graph exploration, "what links to X",
        │                  neighbourhood browsing. Cannot do joins/aggregation.
        │                  Refresh: bash outputs/obsidian_rag_full/refresh.sh
        │                  (or say "옵시디언 갱신해줘")
        └──────────────► outputs/tag_kb/kb.duckdb      [TABLE + QUERY layer]
                           8 node tables + edges(808) + references(95 PDFs)
                           + paper_chunks(4054, BM25 FTS). TAG engine answers
                           relational/aggregation/multi-hop questions the graph
                           can't: syn(NL→SQL) → exec(DuckDB) → gen(answer+cites).
                           Use for: "how many…", joins, "which Param feeds a
                           failing gate", citation-integrity lookups.
                           Query:  python outputs/tag_kb/tag_query.py "<question>"
                           Refresh: bash outputs/tag_kb/refresh.sh
```

**Which layer to use.** Need to *see* structure or follow links → Obsidian.
Need an *exact relational/aggregate answer with citations* → TAG (`tag_query.py`).
Need a passage from a paper → TAG content layer (BM25 over `paper_chunks`). Both
read layers are regeneratable; **never hand-edit the vault or the .duckdb — fix
Notion and refresh.** TAG/Obsidian use the dual LLM backend (anthropic SDK if
`ANTHROPIC_API_KEY`, else `claude -p` headless — no new secret needed).

**Citation integrity (hard).** A 2026-06-02 audit of all 243 SourceEvidence rows
found **3 confirmed hallucinated sources** (`Yao2011_NatCommun`→KB-3.18/3.19,
`YapKovacs_JCS`→KB-4.1, `NanoConvergence2021_Glioma`→KB-6.2.3); verdicts live in
the `source_audit` table (`web_verdict='HALLUCINATION'`) and
`outputs/tag_kb/AUDIT_FINDINGS.md`. Before citing a KB source in any
deliverable, check it isn't one of these (the α-actinin `k_off0` was re-anchored
0.4→**0.066 s⁻¹**, Ferrer 2008, as a result). Dominant issue is metadata drift,
not fabrication (~1.2% fabricated).

**Architecture benchmark.** `outputs/tag_kb/kb_benchmark.py` measures how much
each KB layer reduces hallucination (4-condition ablation C1 no-KB → C2 RAG →
C3 +Obsidian → C4 +TAG, scored with RAGAS/FActScore/TAG-Bench-style metrics);
see `outputs/tag_kb/BENCHMARK_REPORT.md`. Run it (controlled pure-LLM, tools
disabled) to re-validate after KB changes.

> ⚠️ TAG = **Table-Augmented Generation** (Biswal et al. 2024), *not* "태그/label".

## Code conventions

- Type hints everywhere, Google-style docstrings.
- One file = one concept. No "kitchen-sink" modules.
- Configuration via YAML in `ffn_sim/validation/oracles/configs/` for KU-anchored literature constants, and per-unit YAMLs for runtime parameters.
- Logging via stdlib `logging`; default level `INFO`, debug per-module.
- Tests in `ffn_sim/tests/`; oracle cross-checks read from `ffn_sim/validation/oracles/`.

## Hard rules (carried over from v1, still apply)

- **No empirical magic numbers**. Every tuning constant must satisfy the Magic-Number Block (derivable, grid-invariant, not chosen to make a gate pass). If it can't, halt and surface to PI.
- **No gate-loosening**. Sanity gates are validation contracts written before the run. If a gate is wrong, surface to PI for a contract change rather than editing inline.
- **No fitting to PI experimental data**. Literature-first; PI data is overlay-only when comparing later. (PI experimental CSVs live in `~/ActiveCellSim/data/experimental/`, not in this repo.)
- **Initialize at the physiological operating point (PI 2026-06-04, HARD).** Every parameter
  and initial/boundary condition must be set to its real *in-vivo* physiological value **from
  the start of the run** — never a convenient null/default that happens to be easy (water
  viscosity instead of the ~65 Pa·s MCF7 cytoplasm; zero baseline turgor / a relaxed
  unpressurised shell instead of the resting ~40 Pa osmotic pressure; an un-anchored ECM, etc.).
  It is **invalid production** to start from a non-physical baseline and then expect the model
  to emergently align with reality and compare it to real measured data — the comparison is
  meaningless because the *baseline state itself* is wrong. The realistic physiological state is
  the BASELINE; perturbations and measurements happen FROM that baseline (e.g. measure cortical
  tension on an already-turgor-pressurised cell at real cytoplasm viscosity, then add myosin as
  a modulator). When a module is `additive/default-off` for backward-compat, the **production
  config must still turn it ON at its physiological value** — default-off is a code-hygiene
  convenience, not a license to run physics from an unphysical zero. Worked failures this rule
  came from: cytoplasm viscosity left at water → wrong dynamics; enclosed-volume Π₀=0 (force-free
  shell) → cortex never pre-tensioned, myosin asked to generate the whole tension from a floppy
  bag and floored ~1000× under band (2026-06-04 saturation diagnosis). Before any production
  run or data comparison, audit that EVERY compartment is at its physiological setpoint; if a
  physiological value is unknown, surface to PI rather than defaulting to a null.
- **Sanity Gate Protocol** mandatory before first execution of any physics/numerics module: dimensional analysis, boundary cases, conservation invariants, numerical sanity (CFL, precision), sign-sense check, measurement-protocol consistency. Record as docstring `Sanity Gate` section or sibling `*_sanity.md`.
- **Visualize at unit / milestone closeout**. Every time a Phase 1 unit (H.X) reaches ✅ DONE, every time a production run lands, every time a sanity-gate sweep flags a non-trivial finding, generate or refresh the relevant PNG figures into `ffn_sim/outputs/h{X}/figs/` via `ffn_sim/scripts/h1_h2_vis.py` (extend the script as new units land — keep one entry-point that regenerates everything). REPORT.md must include a `## Figures` section listing each figure with a one-line caption. Rationale: PI 2026-05-21 — text-only `npz/json/log` artefacts hide the actual physics, and the visual check often surfaces issues (outliers, wrong-sign decay, sampling artefacts) that bands miss. Figures are persistent project state, not session-ephemeral; commit them alongside the data.
- **Visualization integrity rules**: no axis truncation, no log/linear flip without note, always show the brief reference (or band) overlaid on the measurement, always annotate units (SI default). For ensemble plots, show per-realisation thin lines + ensemble mean overlay so the reader can judge stochastic spread.

## Where v1 lives

The v1 project (Taichi MLS-MPM spheroid continuum + image-constrained intermediate v2 attempt + `acs_kb/` knowledge-base) is preserved untouched in `~/ActiveCellSim` (origin remote of this repo, local-path remote). It is **not** retained inside `ffn_cellsim`. The surviving validation oracles from `acs_kb/` were moved into `ffn_sim/validation/oracles/` during the 2026-05-20 rename; the rest was deleted from this working tree. The codebase-audit rationale lives at `ffn_sim/docs/v2_audit/CODEBASE_AUDIT.md`.

## Active roadmap

- Phase 0 (foundation): **closed 2026-05-19**. See `ffn_sim/docs/PHASE_0_CLOSEOUT.md`.
- Phase 0.3 design decisions: `ffn_sim/docs/PHASE_0_3_DECISIONS.md` (7 PI-ratified calls).
- AFINES master review: `ffn_sim/docs/AFINES_ALGORITHM_NOTES.md` (897 lines).
- Phase 1 unit briefs: `ffn_sim/docs/briefs/H{1,2,3,4}_*.md` (Owner + Prereq lines at top of each).

## Working model (redesigned 2026-05-28 — PI-ratified)

Phase 1 runs as a **single Lead Claude Code session** driving the whole H.1 → H.2 → H.3 → H.5 → H.7 chain, with **on-demand subagents** for parallelism and context isolation.

The previous 4-role split (Main / Sub / Orchestrator + PI-as-relay) is **retired**. It cost PI copy-paste relay between sessions, bred stale boot-prompts (Orchestrator #12 booted "doubly stale"), and caused the `phase1/h4-fa-clutch` branch-divergence (Main + Sub shared one working tree). With this model's 1M context + subagents, the context-efficiency rationale for splitting no longer holds. State now lives on disk + Dev Logs as **durable state a fresh session boots from** — not in a human relay.

### Roles

- **Lead session** (Claude Code, this session) — owns H.1 → H.2 → H.3 → H.5 → H.7 end-to-end: writes code, drives production (SSH to Gbook / background jobs / Syncthing 회수), maintains Notion at closeout. Self-orchestrates; there is no separate read-only orchestrator and no next-prompt relay.
- **Subagents** (spawned by the Lead, no PI relay):
  - `isolation:"worktree"` coding subagent for genuinely independent modules (e.g. H.4 FA/bridge) — true file isolation; the Lead integrates. This is the structural fix for the shared-tree divergence.
  - `Explore` / `Plan` subagent for heavy search / reads — keeps the Lead's context clean.
  - `code-review` (or a review subagent) at closeout — independent audit pass, replacing the old orchestrator's review function without a standing role.
- **PI (Sungwook)** — **approver, not relay**. Talks directly to the Lead. Final sign-off on: gate-contract changes, magic-number triggers, ✅ DONE ratification, `ffn/foundation` pushes, and `integrator/`-freeze changes.

### State stores (read on session boot)

- `ffn_sim/docs/briefs/H*.md` — **immutable unit specs** (contracts). Not edited.
- Notion [Development Logs & Reviews](https://www.notion.so/365120daec5d81969e74ffbb757d55c8) — **Phase 1 status board** (Status / Owner / Start / End / Log) + per-unit milestone Day-logs + **Open items**. The authoritative project record.
- On disk — latest commit on `phase1/h{N}-*`, `ffn_sim/outputs/h{X}/REPORT.md`, working tree.
- **Knowledge base** — Notion Contract-Graph (SoT) + the Obsidian graph mirror and TAG/DuckDB query layer under `ffn_sim/outputs/{obsidian_rag_full,tag_kb}/`. See the *Knowledge base — RAG + TAG + Obsidian* section above for which layer to query and how. Ask KB questions via `tag_kb/tag_query.py`, not by guessing from memory.
- The Notion [Session Handoff Board](https://www.notion.so/366120daec5d815da389c38bc3bfbbe1) is **retired for relay** (banner at its top). No more next-prompt drafts, stale-markers, cross-session signal-routing, or PI pre-dispatch checklist.

### Session boot protocol (minute 0)

1. Read `CLAUDE.md` (this file) fully.
2. Read the Dev Logs status board + Open items.
3. `git log --oneline -5`; confirm current branch + latest commit.
4. `conda activate ffn_sim`; confirm `python -c "import hoomd; print(hoomd.version.version)"`.
5. Restate the next task in 1–2 sentences (derived from the status board + open items — no boot-prompt needed). Ask PI if ambiguous.

### Session closeout protocol — MUST (end of session, freeze-point, or when PI says "wrap")

Steps 1–4 complete before the final user-facing message; step 5 is the explicit PI receipt.

1. Commit work on `phase1/h{N}-*`. Do NOT push to `ffn/foundation` (renamed from `v2/foundation` 2026-05-20) without PI sign-off.
2. **Dev Logs status board** — update the unit's row (Status / Owner / Start / End / Log) and the **Open items** list. Append a `Phase {N} — Unit H.{X} {milestone}` child Day-log (start/end commit hashes, sanity gate PASS/FAIL, next-unit dependency check, KU cross-reference — format in Dev Logs §작성 규칙). Large-body table-cell edits time out (~100K page); prefer a small milestone child page + a short inline note (see `reference_notion_handoff_board_size`).
3. If a unit hits ✅ DONE / 🚧 / a production run lands / a sanity-gate flags a non-trivial finding: refresh figures per the visualize-at-closeout rule.
4. If Notion MCP is unavailable or a store cannot be written: **halt and surface to PI** — do not silently skip.
5. **Receipt**: the final user-facing message ends with the literal line **`Notion 업데이트 완료`**.
6. Stop. Do not speculate beyond what was actually done.

### File ownership

The single Lead owns all of `ffn_sim/`. When a `worktree` subagent is spawned for a parallel module, that subagent owns its module directory for the duration and the Lead does not touch those files until integration. Cross-module changes are sequenced by the Lead. `ffn_sim/validation/oracles/`, `ffn_sim/docs/briefs/`, `STRUCTURE.md`, `README.md`, `pyproject.toml` stay read-only during normal unit work.

## When in doubt

- Read `ffn_sim/docs/PHASE_0_CLOSEOUT.md` for current state.
- Read `ffn_sim/docs/PHASE_0_3_DECISIONS.md` for the ratified design vocabulary.
- Read the relevant `ffn_sim/docs/briefs/H*.md` for the unit you're touching.
- Read the Notion [Development Logs status board](https://www.notion.so/365120daec5d81969e74ffbb757d55c8) + Open items for the current task.
- Ask the PI before deviating from any principle in this file.

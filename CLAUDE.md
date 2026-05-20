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
- **Hardware through Phase 1**: Apple silicon CPU (M1 Max benchmark: 1 M-step polymer in ~11 s)
- **Phase 2+ target**: CUDA GPU (HOOMD GPU build); no hard-coded device IDs

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
- **Sanity Gate Protocol** mandatory before first execution of any physics/numerics module: dimensional analysis, boundary cases, conservation invariants, numerical sanity (CFL, precision), sign-sense check, measurement-protocol consistency. Record as docstring `Sanity Gate` section or sibling `*_sanity.md`.

## Where v1 lives

The v1 project (Taichi MLS-MPM spheroid continuum + image-constrained intermediate v2 attempt + `acs_kb/` knowledge-base) is preserved untouched in `~/ActiveCellSim` (origin remote of this repo, local-path remote). It is **not** retained inside `ffn_cellsim`. The surviving validation oracles from `acs_kb/` were moved into `ffn_sim/validation/oracles/` during the 2026-05-20 rename; the rest was deleted from this working tree. The codebase-audit rationale lives at `ffn_sim/docs/v2_audit/CODEBASE_AUDIT.md`.

## Active roadmap

- Phase 0 (foundation): **closed 2026-05-19**. See `ffn_sim/docs/PHASE_0_CLOSEOUT.md`.
- Phase 0.3 design decisions: `ffn_sim/docs/PHASE_0_3_DECISIONS.md` (7 PI-ratified calls).
- AFINES master review: `ffn_sim/docs/AFINES_ALGORITHM_NOTES.md` (897 lines).
- Phase 1 unit briefs: `ffn_sim/docs/briefs/H{1,2,3,4}_*.md` (Owner + Prereq lines at top of each).

## Multi-session orchestration

Phase 1 runs across multiple Claude Code sessions for context efficiency. State lives in Notion + on-disk briefs, NOT in any one session's context.

### Roles

- **Main Session** — sequential ownership of H.1 → H.2 → H.3 → H.5 → H.7 (architectural consistency chain).
- **Sub Session** — H.4 (FA + motor-clutch, isolated module) after Main's BAOAB freeze, plus debugging interludes and sanity-gate updates after FAIL.
- **Orchestrator (PI / Sungwook)** — relays state between sessions, reviews each session's closeout, approves each next prompt.

### State stores (read on session boot)

- `ffn_sim/docs/briefs/H*.md` — **immutable unit specs** (contracts). Sessions do not edit.
- Notion [Session Handoff Board](https://www.notion.so/366120daec5d815da389c38bc3bfbbe1) — **per-session closeout + next-prompt drafts** + file-ownership table.
- Notion [Development Logs & Reviews](https://www.notion.so/365120daec5d81969e74ffbb757d55c8) — **Phase 1 status board** (Status / Owner / Start / End / Log).

### Session boot protocol (minute 0)

1. Read `CLAUDE.md` (this file) fully.
2. Open the Notion Session Handoff Board, find your role section (Main or Sub), read **only your "Next prompt"**.
3. Verify branch + commit hash match what the prompt says (`git log --oneline -1`).
4. `conda activate ffn_sim` and confirm with `python -c "import hoomd; print(hoomd.version.version)"`.
5. Restate the task in 1–2 sentences before executing anything. If anything is ambiguous, ask PI before moving.

### Session closeout protocol (end of session or when PI says "wrap")

1. Stage and commit work on the session's branch (`phase1/h{N}-*`). Do NOT push to `v2/foundation`.
2. Write a **closeout block** to your section of the Notion Session Handoff Board, using the template pinned there.
3. Write a **draft Next prompt** in that section so PI can review/edit and the next session of your role can boot from it.
4. Update the Phase 1 status board (Status / Owner already set; Start / End / Log columns as appropriate).
5. Stop. Do not speculate beyond what was actually done in the session.

### File ownership (cross-session enforcement)

Detailed table on the Session Handoff Board. Headline:

- **Main owns**: `ffn_sim/integrator/`, `ffn_sim/ecm/`, `ffn_sim/cortex/`, `ffn_sim/cell/`, `ffn_sim/common/`, `ffn_sim/configs/phase1_h{1,2,3,5,7}.yaml`
- **Sub owns**: `ffn_sim/bridge/`, `ffn_sim/configs/phase1_h4.yaml`, `*_sanity.md` updates
- **Both write to** `ffn_sim/tests/`: per-unit files only (`test_h{N}_*.py`); no cross-touching
- **Read-only for both**: `ffn_sim/validation/oracles/`, `ffn_sim/docs/briefs/`, `CLAUDE.md`, `STRUCTURE.md`, `README.md`, `pyproject.toml`
- **Cross-boundary change needed**: escalate to PI, do not write directly. One conflict avoided > a few minutes of relay time.

### Closeout block template (copy into your handoff board section)

```
## Closeout — {YYYY-MM-DD} · {Main|Sub} session
- **Branch**: `phase1/h{N}-{name}` @ `{commit_short}`
- **Completed**: <concrete deliverables landed>
- **Tests**: <added/changed, pass/fail counts>
- **Sanity gates run**: <list, PASS/FAIL>
- **Files touched**: <top 5, full count>
- **Open for PI**: <unanswered questions, magic-number triggers, gate-contract questions>
- **Recommended next prompt** (PI to review):
  > <draft next-session task, 3-5 lines, copy-paste-ready>
```

## When in doubt

- Read `ffn_sim/docs/PHASE_0_CLOSEOUT.md` for current state.
- Read `ffn_sim/docs/PHASE_0_3_DECISIONS.md` for the ratified design vocabulary.
- Read the relevant `ffn_sim/docs/briefs/H*.md` for the unit you're touching.
- Read the Notion [Session Handoff Board](https://www.notion.so/366120daec5d815da389c38bc3bfbbe1) for your role's current task.
- Ask the PI before deviating from any principle in this file.

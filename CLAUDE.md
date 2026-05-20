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
- Phase 1 unit briefs: `ffn_sim/docs/briefs/H{1,2,3,4}_*.md` (Owner + Prereq at top of each).

## Multi-session orchestration (state in Notion, not in this file)

Phase 1 runs across separate Claude Code sessions: **Main** (H.1 → H.2 → H.3 → H.5 → H.7), **Sub** (H.4 + debug), **Orchestrator** (read-only review/drafting). On session boot, read the Notion Session Handoff Board for your role's `[ready]` Next prompt — the boot template, closeout template, and file-ownership table all live there.

- 🎛️ **Session Handoff Board**: https://www.notion.so/366120daec5d815da389c38bc3bfbbe1
- 🧾 **Dev Logs status board**: https://www.notion.so/365120daec5d81969e74ffbb757d55c8
- 🚀 **Build Plan v2** (full §3 Phase 1 spec): https://www.notion.so/365120daec5d81799efefcf078f2039e

PI is the human relay between sessions — workers write closeouts to the Handoff Board, Orchestrator drafts the next prompts, PI reviews/approves. Do NOT push branches to `v2/foundation` without PI sign-off.

## When in doubt

- Read `ffn_sim/docs/PHASE_0_CLOSEOUT.md` for current state.
- Read `ffn_sim/docs/PHASE_0_3_DECISIONS.md` for the ratified design vocabulary.
- Read the relevant `ffn_sim/docs/briefs/H*.md` for the unit you're touching.
- Ask the PI before deviating from any principle in this file.

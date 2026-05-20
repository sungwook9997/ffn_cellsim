# ffn_cellsim — Full Fiber Network Cell Simulator

**Mechanistic, fine-grained HOOMD-blue framework for single-cell mechanobiology.**

Every cytoskeletal filament, every motor head, every adhesion clutch, and every ECM cross-link is an explicit particle / bond — no paper-model wrappers, no lumped mechanisms.

## Why "full fiber network"?

Most published cell-mechanics simulations borrow closed-form models from the literature (Chan-Odde motor-clutch, Pereverzev catch-bond, Bell-Evans, Buckley cadherin, Hill force-velocity) and run them as the **runtime mechanism**. `ffn_cellsim` inverts that: the runtime is fine-grained HOOMD particle/bond dynamics, and those literature closed-forms become **acceptance oracles** — used only in validation tests to cross-check that the emergent behaviour matches the published expressions.

PI directive (2026-05-19): *"full resolution simulation 방향에 맞게 — 처음 프레임워크에 추상화나 다른 과정이 들어가면 안됨."*

## Layout

```
ffn_cellsim/
├── ffn_sim/                          # the active HOOMD-blue runtime
│   ├── ecm/                          # Mikado fiber network + cross-links (H.1)
│   ├── cell/                         # cell body composition (H.3)
│   ├── cortex/                       # actin cortex (H.2, H.3)
│   ├── bridge/                       # focal adhesion + motor-clutch (H.4)
│   ├── junction/                     # E-cadherin cell-cell adhesion (H.6)
│   ├── integrator/                   # Leimkuhler-Matthews BAOAB plugin (custom)
│   ├── common/                       # shared utilities
│   ├── validation/
│   │   └── oracles/                  # v1 closed-form oracles (frozen reference)
│   │       ├── ecm/                  # WLC, Mikado geometry, cross-link energy
│   │       ├── junction/             # Young-equation contact angle
│   │       ├── bridge/               # 1-D traction reducer
│   │       ├── common/               # sanity gates, derived parameters
│   │       └── configs/              # KU-anchored literature constants (YAML)
│   ├── docs/                         # active project docs
│   │   ├── PHASE_0_CLOSEOUT.md
│   │   ├── PHASE_0_3_DECISIONS.md    # 7 PI-ratified design calls
│   │   ├── AFINES_ALGORITHM_NOTES.md
│   │   ├── briefs/                   # per-unit dispatch briefs (H.1..H.7)
│   │   └── v2_audit/                 # codebase audit + _DEPRECATED rationale
│   ├── tests/                        # HOOMD runtime tests
│   ├── scripts/                      # CLI entry points
│   └── outputs/                      # per-unit deliverables
├── CLAUDE.md                         # Claude Code project context
├── STRUCTURE.md                      # detailed file map
├── pyproject.toml
└── requirements.txt
```

## Stack

- Python 3.13
- HOOMD-blue 7.0.1 (CPU + GPU; Apple silicon CPU-only through Phase 1)
- Leimkuhler-Matthews BAOAB-limit custom integrator plugin
- conda env: `conda activate ffn_sim`

## Status

Phase 0 closed 2026-05-19 (`v2/foundation` branch). Phase 1 dispatch:

| Worker | Unit | Brief |
|---|---|---|
| A | H.1 ECM Mikado | `ffn_sim/docs/briefs/H1_ecm_mikado.md` |
| B | H.4 FA + motor-clutch | `ffn_sim/docs/briefs/H4_fa_motor_clutch.md` |
| C | H.2 single filament → H.3 cortex → H.5 lamellipodium | `ffn_sim/docs/briefs/H2_single_filament.md` |

## Origin

This repository was carved out of `~/ActiveCellSim` (v1 Taichi MLS-MPM continuum) on 2026-05-20 as a clean ffn-only workspace. The v1 spheroid-MPM code and the image-constrained intermediate v2 attempt are preserved in the parent repo (`origin` remote → `~/ActiveCellSim`) and are not retained here.

## PI / Lab

Sungwook Yoon · Shin Lab, Dept. of Mechanical Engineering, KAIST

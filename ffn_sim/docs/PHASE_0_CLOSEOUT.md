# Phase 0 — closeout report

**Period**: 2026-05-19 (single intensive session)
**Branch**: `ffn/foundation` (renamed from `v2/foundation` 2026-05-20; 5 commits, +2,179 LOC of docs + 210 LOC of skeleton/scripts)
**PR**: <https://github.com/sungwook9997/ActiveCellSim/pull/2>

## Units completed

| Unit | Plan v2 budget | Actual | Deliverable | Commit |
| --- | --- | --- | --- | --- |
| H.0.1 — Codebase audit | 1 day | done | [ffn_sim/docs/v2_audit/CODEBASE_AUDIT.md](../../ffn_sim/docs/v2_audit/CODEBASE_AUDIT.md), [ffn_sim/docs/v2_audit/_DEPRECATED.md](../../ffn_sim/docs/v2_audit/_DEPRECATED.md) | `6de994e` |
| H.0.2 — HOOMD env | 2 days | done | conda env `ffn_sim` (py 3.13, hoomd 7.0.1 CPU), polymer sanity 89k steps/s, [REPORT.md](../outputs/h_0_2/REPORT.md) | `676e970` |
| H.0.3 — AFINES algorithm review | 5 days | done | [AFINES_ALGORITHM_NOTES.md](AFINES_ALGORITHM_NOTES.md) (897 lines), [PHASE_0_3_DECISIONS.md](PHASE_0_3_DECISIONS.md) (7 design calls ratified) | `1406675`, `0b02820` |
| H.0.4 — `ffn_sim/` skeleton + initial commit | 2 days | done | 13-directory layout, all `__init__.py` in place | `676e970` |

Plan v2 §2 estimated H.0 at ~2 weeks; condensed into one day of agent work
because the audit, env setup, and algorithm review are heavily research-
parallelisable.

## Key Phase 0 findings that fed back into the Plan

### Scope changes

- **`acs_kb/ [v1, deleted] /` archive ≈ 25 % → 45 %** after the 2nd PI pass ("v1 paper-
  models-as-mechanism, v2 oracle-only"). 18 files archived, ~3,700 LOC.
  Survivor v2 runtime: geometry/measurement/KU contracts/literature
  constants only.
- **H.5 lamellipodium = greenfield** (AFINES has no Arp2/3 branching
  code anywhere). Plan §3 unit budget 2 wk → 4 wk.
- **H.4 motor architecture = Stam-Hocky multi-head minifilament**
  (full-fidelity over AFINES single 2-head spring). +1 wk.
- **Cadherin full catch-bond (KU-4.2) moved to Phase 1**, overriding the
  Plan's slip-only (KU-4.17) Phase 1 baseline.
- **Net Phase 1**: ~3 months → ~3.5–4 months.

### Environment changes

- **HOOMD 7.0.1**, not the Plan's ≥4.4 assumption (conda-forge stable
  jumped). 4→7 API drift confined to `hoomd.version.cuda_built` →
  `gpu_enabled`; our needed surface intact.
- **Python 3.13.13** instead of 3.11.
- **Apple silicon CPU-only** through Phase 1 entry (Plan §0 hardware row
  ratified; M1 Max benchmark shows 1 M-step polymer in ~11 s, beats the
  Plan's A100 reference).

### Plan §2 H.0.3 corrections

- AFINES canonical = `github.com/Simfreed/AFINES`, not Shibalab-Lehigh
  (latter does not exist).
- AFINES is 2D, CPU-only, single-threaded.
- AFINES motors/xlinks = Metropolis (Glauber detailed-balance), not Bell-
  Evans. v2 D2 is an upgrade.
- AFINES integrator = Leimkuhler-Matthews BAOAB-limit (not specified in
  Plan).

## Phase 1 readiness

All blockers resolved. Worker dispatch can begin immediately:

| Worker | Unit | Brief | Prereq |
| --- | --- | --- | --- |
| A | H.1 ECM Mikado | [briefs/H1_ecm_mikado.md](briefs/H1_ecm_mikado.md) | none (start in parallel with B and C) |
| B | H.4 FA + motor-clutch | [briefs/H4_fa_motor_clutch.md](briefs/H4_fa_motor_clutch.md) | none (start in parallel with A and C) |
| C | H.2 single filament | [briefs/H2_single_filament.md](briefs/H2_single_filament.md) | none (start first; H.3 depends) |
| C | H.3 cortex multi-filament | [briefs/H3_cortex.md](briefs/H3_cortex.md) | H.2 done |
| C | H.5 lamellipodium | (defer brief until H.3 design choices land) | H.3 done |

Worker D's H.7 integration test brief deferred until A/B/C deliverables
land.

## Operational notes for future sessions

- v1 archived files in `acs_kb/ [v1, deleted] /` are **forbidden imports** from
  `ffn_sim/<runtime>/*.py`. They live as **closed-form oracles** for
  `ffn_sim/tests/` and `ffn_sim/validation/` only.
- Plan v2 Notion page (<https://www.notion.so/365120daec5d81799efefcf078f2039e>)
  contains §12 with all 7 PI decisions ratified — start any v2 session
  by reading §12 + this closeout + the brief for the assigned Unit.
- conda env: `conda activate ffn_sim` (env created in this session).
- Sanity bench reproducer: `python ffn_sim/scripts/hoomd_polymer_sanity.py --steps 50000 --bench`.

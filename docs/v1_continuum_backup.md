# V1 Continuum Prototype Backup

## Status

This file marks the current spheroid-first continuum implementation as the
preserved **v1 continuum prototype** before v2 development begins.

Snapshot:

- repo HEAD when v2 was opened: `e424bee`
- date: 2026-05-02
- preserved model identity: Taichi MLS-MPM / active viscoelastic continuum /
  Marangoni-style reduced spheroid model
- new development identity: image-constrained cell-resolved mechanobiology,
  described in `docs/00_project_vision_v2.md`

## What Is Preserved

The following remain valid project assets and should not be deleted during the
v2 transition:

- `acs/physics/mlsmpm.py` — v1 reduced/continuum solver
- `acs/runner.py` — v1 run orchestration and gate report generation
- `configs/stage*.yaml` and `configs/full_stack*.yaml` — v1 experiment records
- `results/` — v1 output artifacts and gate reports
- `docs/outcomes_*.md` — v1 outcome ledger
- `docs/*_sanity.md` — v1 sanity-gate analyses
- `docs/gate_fail_taxonomy.md` — v1 gate classification framework
- `docs/00_project_vision.md` — original spheroid-first framing

## Role In V2

V1 is not the biological ground-truth engine. It becomes the first reduced
surrogate branch against which the future v2 cell-resolved simulator can be
compared.

V1 can still answer:

- whether continuum active-gel/MPM approximations create plausible spheroid
  spreading trajectories
- where radial or Marangoni-style approximations break down
- how much speed a reduced model can provide relative to cell-resolved
  simulation

V1 should no longer absorb new detailed biology directly into the monolithic
MPM solver unless the v2 architecture explicitly requires a reduced-model
comparison feature.

## Transition Rule

New biology starts under `acs/v2/`.

Existing v1 files should be treated as read-mostly while the v2 data contract,
single-cell model, and image-analysis pipeline are being established.


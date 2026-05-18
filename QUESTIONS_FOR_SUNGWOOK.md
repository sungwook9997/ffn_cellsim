# Questions / open items — Worker D Unit 4.1 (Junction Track)

These are deviations or open questions surfaced during Phase 1 Unit 4.1
that need PI direction. They are referenced from
`acs_kb/outputs/phase1/unit4_1/REPORT.md` so the work history is traceable.

## Q1 — Branch base: brief says "off main" but main has no `acs_kb/`

The brief instructs `git checkout -b worker-d/junction` off `main`. On
inspection, `main` (HEAD `4ffcbbe`) carries only `acs/` (v1/v2 Taichi code)
and treats the entire `acs_kb/` tree as untracked. The Worker A/B/C
commits that build the KB-aligned package live on `worker-b/bridge` (HEAD
`2b98d4d`), which is the only branch where `acs_kb/cell/cell.py` is
committed.

**Reasonable call taken**: `worker-d/junction` was created off `main`
literally (per the brief), then Worker A/C dependencies needed for
`from acs_kb.cell.cell import Cell` to succeed were extracted from
`worker-b/bridge` into the working tree **as untracked files** via
`git archive worker-b/bridge … | tar -x`. Only Worker D's own new files
are staged for commit on `worker-d/junction`. The branch history stays
clean (no Worker A/B/C commits as ancestors), matching the spirit of
the "direct A/B import forbidden" rule.

**Decision needed**: at integration time, PI should merge `worker-d/junction`
into a branch that already has Worker A+B+C content (e.g., an integration
branch built from `worker-b/bridge` + future `worker-a/*` + `worker-c/*`).
A standalone CI of `worker-d/junction` will not import — that is by
construction and is acceptable to PI as long as the integration plan is
explicit.

## Q2 — `Cell.junctions` field missing

The brief states that `Cell.junctions: list[EcadherinJunction]` should
exist on Worker C's `Cell` dataclass; if missing, Worker D adds the
field "after Notion notification of Worker C". As of `worker-b/bridge`
HEAD `2b98d4d`, `Cell.__dataclass_fields__` = `{id, cortex,
center_position, polarity, focal_adhesions, inner_mode, surface_tension,
nucleus_position}` — no `junctions`.

**Reasonable call taken**: Worker D's Unit 4.1 code does **NOT** modify
`acs_kb/cell/cell.py` in this branch. Junctions are held in a free list
external to `Cell` (see the two-cell pair test in
`acs_kb/tests/test_two_cell_pair.py`). The field-add is queued for the
Unit 4.1 Notion progress note (Task 5) and is deferred to the integration
step so Worker C remains the sole author of cell.py edits.

**Decision needed**: PI to confirm Worker C will add
`junctions: list = field(default_factory=list)` to `Cell` before Worker D
Unit 4.2 begins, or to authorise Worker D to make the additive change
directly on integration.

## Q3 — `.venv-collab/site-packages/numpy` is in a Syncthing-conflict state

Encountered while running the Task-0 import smoke test. The project's
default venv has `__init__.sync-conflict-20260518-203604-*.py` files
instead of the canonical `__init__.py` for the `numpy/`, `scipy/`, and
related package directories — clearly a Syncthing replication conflict
that landed mid-session.

**Reasonable call taken**: created a side-by-side ephemeral venv at
`/tmp/acs-worker-d-venv` (Python 3.12, numpy 2.4.5, scipy, pyyaml,
pytest) and used that for local validation. No project files were
modified to work around the conflict.

**Decision needed**: PI to resolve the Syncthing conflict in
`.venv-collab/` at a time of their choosing (cleanest fix is usually
`rm -rf .venv-collab && python3.12 -m venv .venv-collab && pip install -r
requirements.txt`). This is unrelated to Worker D scope.

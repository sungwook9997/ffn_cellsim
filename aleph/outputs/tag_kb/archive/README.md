# Superseded generated candidate manifests

Nothing here is current state, and nothing here is deleted. These are **dry-run outputs** — manifests a
generator wrote for a PI to review on a given day — kept because they record what was proposed and when,
which is occasionally the only evidence of why a graph row does or does not exist.

## Why they moved (2026-07-29)

`outputs/tag_kb/` held 32 markdown files, **21 of them dated**, including five generations of
`OPS_HARVEST_CANDIDATES` and a stack of one-off audit snapshots. A reader looking for the operational
entry point (`README.md`) or the current harvest proposal had to work out which of five dated files was
the live one. That is the directory mess, not the file count.

## What moved, and the rule used

| Moved | Kept in place | Why |
|---|---|---|
| `OPS_HARVEST_CANDIDATES_2026-06-08 / 06-30 / 07-06 / 07-22` | `..._2026-07-25` (newest) | `harvest_ops.py:565` WRITES a new dated file per run and never reads an old one, so only the newest is live |
| `SE_DEDUP_PLAN_2026-06-02` | — | zero references anywhere in the tree |
| `KB_CURATION_CANDIDATES_2026-07-22` | — | zero references |
| `KB_INTEGRITY_CHECK_2026-07-23` | — | zero references |
| `KB_OVERNIGHT_SUMMARY_2026-06-20` | — | zero references |

**`SE_REGISTRATION_CANDIDATES_*` was deliberately NOT touched** — 22 files across the tree reference that
family, so moving it would break links for a tidiness gain. Same for `CLIP_CANDIDATES.md`, which
`paper_clip.py:1406` writes to by default.

**The rule applied**: move a generated artifact only when a newer generation of the same family exists
*and* nothing in the tree reads it. Anything referenced, anything newest-of-family, and anything hand-
authored stayed where it was.

## Where current state lives

Row counts, edge counts and audit tallies are **read from `kb.duckdb`**, never from a document in this
directory or its parent — `bash outputs/tag_kb/refresh.sh` reports drift non-destructively. `CLAUDE.md`
records that the last hard-coded set of six such numbers was wrong in all six.

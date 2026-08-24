# KB Integrity + Ops-Linkage Sweep — 2026-07-23

Stream ⑬ (KB integrity + Ops-Linkage). Read-only sweep — no KB edits, no code
changes. All edits to the KB go through Notion (SoT); this report only audits the
materialized read layers and the disk↔graph drift.

## Verdict: GREEN — no DRIFT, no halt condition.

Both blocking gates pass; both drift checks are clean; the content/supersession
layers are current. Nothing to surface to PI as a contract/citation drift.

## Gate results

| Gate | Result | Detail |
|---|---|---|
| `verify_runs.py --gate` (results) | **PASS** | 31 claims, none drifted — disk supports every declaration |
| `verify_params.py --gate` (params) | **PASS** | 42 constants, none drifted — disk + citation audit support every declaration |
| `verify_sources.py --check` (citation) | informational | see below |
| `harvest_ops.py --check` (ops drift) | **0 un-harvested** | RunResult 87 disk / 98 graph; CodeMapping 9 modules / 50 mapped |
| `references_ingest.py --check` (content) | current | 278 paper_refs, 12 810 chunks; no PDFs newer than 2026-07-16 |
| `supersession.py --check` (auth chain) | **PASS** | KU-3.5 authoritative record + demotion intact |

## Citation audit (informational, non-blocking)

`source_audit` table (551 rows, rebuilt against the fresh Notion snapshot):

| verdict | rows |
|---|---|
| OK | 313 |
| CHECK | 167 |
| NO_DOI_FOUND | 49 |
| DOI_MISMATCH | 15 |
| DOI_DEAD | 7 |

Dominated by OK + CHECK (metadata drift), consistent with the standing ~1.2 %
fabrication-rate finding — **no new fabrication signal**. The production-cited
subset is already enforced through the params gate (a fabrication-risk citation on
a `verified` constant would drift and block); params gate is green, so no
production citation is at risk.

`make kb-check` emitted a **PARTIAL/STALE** notice — 13 of the 551 `source_audit`
rows carry a stale (`--limit`-sampled) verdict. This is informational, not a drift.
SourceEvidence has grown 329 → **564** rows since the last full audit (parallel
sessions added PI-GAP, Kim–Miyazaki, mcf7-geometry, biot-coupling sources), so a
**full network re-audit** (`python verify_sources.py`, DOI resolution) is a
worthwhile PI-gated follow-up to re-confirm the new ~235 sources — but it is
network-bound and outside this read-only sweep.

## Ops-Linkage — new parallel-session outputs

New untracked outputs from the other AC sessions all live under `aleph/outputs/ac/`
(`bleb/`, `implicit/`, `osmotic/`, `cell_assembled/`, `stage_gallery/`). **None are
under a `production/` directory**, so none match the harvest contract glob
(`outputs/**/production/*.{json,md}`) — which is exactly why `harvest_ops --check`
reports 0 un-harvested. This is the designed behavior: these are pre-production
evidence-ladder artifacts (cortex fixes, resting-convergence probes, osmotic model
design, assembled-state dumps), not PRODUCTION closeouts. They enter the graph only
when a component reaches PRODUCTION and lands a `outputs/**/production/` closeout
with a `REPORT.md`. **The harvest is not silently missing anything.**

## Read-layer freshness

`kb.duckdb` (Jul 22 22:34) and `snapshots/notion_snapshot.json` (Jul 22 22:35) were
both rebuilt by a token-holding session yesterday, so the token-less checks here run
against a current SoT snapshot. A full `refresh.sh` rebuild is **blocked** in this
session (`NOTION_TOKEN` absent) **and not needed** — the SoT is unchanged since the
snapshot and no new PDFs have landed. Not a "store cannot be written" halt: no KB
edit was attempted; the SoT is intact.

## Housekeeping note (non-blocking)

`aleph/outputs/tag_kb.zip` (25 MB, Jul 18, untracked) is a stray archive sibling
to the `tag_kb/` dir — looks like a manual backup. Not deleted here (untracked,
harmless); flag for PI cleanup or `.gitignore`.

## Recommended follow-ups (PI-gated)
1. When a token-holding session is available, run full `python verify_sources.py`
   to clear the 13 partial rows and re-audit the ~235 new SourceEvidence rows.
2. Remove or gitignore `outputs/tag_kb.zip`.

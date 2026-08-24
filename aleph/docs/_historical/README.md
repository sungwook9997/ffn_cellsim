# Historical session artifacts — kept, not current

Nothing here is deleted and nothing here is current state. These are one-off documents from earlier
sessions that no longer describe how this project works, and that **nothing in the tree references**.

## Why they moved (2026-07-29)

`aleph/docs/` had 33 markdown files at the top level, mixing live design contracts with session
artifacts from June. A reader looking for the current read set had to tell them apart by date and title.

## The rule applied — the same one used for `outputs/tag_kb/archive/`

Move only when **nothing in the tree reads it** *and* it is a session artifact rather than a live record.
Measured before moving: 5 of 33 top-level docs had zero references; 4 were moved and 1 was kept.

| Moved | Why |
|---|---|
| `CLOUD_GPU_OPTIONS.md` | Its entire premise is procuring cloud GPUs to run **HOOMD-blue 7.0.1**. HOOMD is forbidden by I0-A and its code was deleted 2026-07-29; the hardware question was answered by gbook's A5000. |
| `AUTONOMOUS_NOTION_SYNC_2026-06-03.md` | Log of one overnight sync session, including a cron id that no longer exists. |
| `FA_ORACLE_TRIAGE_2026-06-02.md` | Self-labelled "TRIAGE scaffold — NOT a contract, NOT ratified bands". |
| `LAYER2_SUBSTRATE_TRACTION_DESIGN.md` | Design for a layer-2 extension on a line that is not current. |

**`CORTICAL_TENSION_RECORD_2026-06-30.md` was deliberately NOT moved** even though it is also
unreferenced by text: it carries `kb_record:` front-matter and `outputs/tag_kb/supersession.py` READS it.
"Nothing references it" has to mean nothing *reads* it, not merely that no prose links to it — a machine
reader counts.

## Verified, not assumed

`supersession.py` scans `docs/**/*.md` **recursively**, so a subdirectory stays in scope. Its `--check`
output is byte-identical before and after the move (`PASS: KU-3.5 authoritative chain`, and the dead-claim
demotion), which is the evidence that the KB layer did not lose anything.

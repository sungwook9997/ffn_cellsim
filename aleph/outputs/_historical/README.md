# Historical output artifacts — kept, not current

Nothing here is deleted and nothing here is a current read layer. These are one-off artifacts whose
question has since been answered, and that **nothing in the tree reads**.

## The rule applied — the same one used for `ffn_sim/docs/_historical/` and `outputs/tag_kb/archive/`

Move only when **nothing in the tree reads it** *and* it is a session artifact rather than a live record.
A machine reader counts as a reader: `outputs/tag_kb/supersession.py` parses `kb_record:` front-matter, so
a file nothing links to in prose can still be load-bearing.

## `build_demo_vault.py` — moved 2026-07-29

Written 2026-06-01 to answer one question, stated in its own docstring: *"does an Obsidian graph view of the
RAG v2 contract-graph add value over Notion's relation tables?"* It emitted ~30 notes from three seeded
chains (KU-3.5 cortex tension, the KU-5.1 alias split, the ECM substrate chain) so a human could open the
graph view and judge.

**The question was answered yes, and the answer shipped.** `outputs/obsidian_rag_full/` is now the
materialised Obsidian read layer that `CLAUDE.md` §Knowledge base names, at 1,664 files against the demo's
~30. The demo cannot be refreshed even in principle — its data is **hardcoded** from the Notion DBs as they
stood on 2026-06-01, which its docstring is explicit about ("no API needed for the demo"), so it is a
snapshot of a superseded graph rather than a small live export.

**Why it moved rather than stayed.** `obsidian_rag_demo/` sat directly beside `obsidian_rag_full/` in a
directory listing, and telling the feasibility demo from the production vault required opening one of them.
That is the same defect the `ffn_sim/docs/` pass fixed: a reader looking for the current read set should not
have to date-sort candidates first.

**What was verified before moving**, rather than assumed:

- zero references from anywhere outside its own directory — searched `.py`, `.md`, `.sh`, `.yaml`, `Makefile`
  across 1,684 non-output files, for both `obsidian_rag_demo` and `build_demo_vault`;
- the generated `vault/` was **gitignored** (`.gitignore:80`) — 32 local files, 0 tracked — so deleting it
  removed nothing from history, and this script regenerates it deterministically from its hardcoded tables.

# Notion → Obsidian export (RAG v2 Contract-Graph mirror)

Read-only Obsidian mirror of the Notion **Contract Graph** (8-DB RAG v2).
**Notion is the source of truth**; this regenerates the vault from scratch.

## What it produces
A graph vault of **~640 nodes, single connected component**:
literature (SourceEvidence) → KnowledgeClaim (KB-x.y) → ModelContract →
Parameter / ValidationGate → CodeMapping/RunResult/DecisionLedger, plus
claim↔claim dependency edges, Code/Test/Doc nodes from `ffn_sim/`, and the
**Dev Logs board** (hub + status board + ~48 milestone day-logs) cross-linked
into the knowledge graph by id-mention (KB-/VG-/RUN-… → relates::).
Every note has an *Open in Notion* link; papers have DOI links.
Color groups (by node `type`) are pre-baked into `vault/.obsidian/graph.json`.

## Regenerate (on-demand — Notion changed, want the graph fresh)
One command (Notion = source of truth; vault rebuilt from scratch):
```bash
bash ffn_sim/outputs/obsidian_rag_full/refresh.sh          # rebuild
bash ffn_sim/outputs/obsidian_rag_full/refresh.sh --open   # rebuild + open in Obsidian
```
Or, to a Claude session: just say "옵시디언 갱신해줘" / "refresh the Obsidian mirror".

First-time token setup (PI, once):
```bash
# Internal Integration secret, connected to the Contract Graph page
echo 'NOTION_TOKEN=ntn_...' > ffn_sim/outputs/obsidian_rag_full/.notion_token  # gitignored
```
Manual steps (what refresh.sh runs):
```bash
conda activate ffn_sim && cd ffn_sim/outputs/obsidian_rag_full
python notion_to_obsidian.py   # [1/3] Notion → vault (claims/papers/params/gates/contracts + relations)
python add_code_nodes.py       # [2/3] repo code/tests/docs → linked to KB claims
python devlogs_to_obsidian.py  # [3/3] Notion Dev Logs board → day-log nodes cross-linked into the graph
open -a Obsidian ./vault
```

## Files (committed)
- `notion_to_obsidian.py` — main exporter (Notion API → markdown + frontmatter + wikilinks + graph.json)
- `ku_dependencies.py` — KU→KU (and KU-v2 KB→KB) dependency edges extracted from the legacy Unit pages' 🔗 DEPENDENCY blocks; turns the flat layered graph into a connected web
- `add_code_nodes.py` — scans `ffn_sim/` (.py / tests / docs) and links code files to the KB claims they cite
- `devlogs_to_obsidian.py` — mirrors the Notion **Dev Logs & Reviews** board (a page-with-blocks, NOT one of the 8 DBs): hub + status-board tables + ~48 day-log prose pages, each cross-linked to the Contract-Graph by id-mention + `partof:: [[DL_Dev_Logs_Hub]]`. Needs the integration connected to the Dev Logs page too (not just Contract Graph).

## Related: getting references/ PDFs into the graph
A references PDF becomes an Obsidian paper node only once it is a **Notion
SourceEvidence row** (the mirror reads Notion only). The TAG ingest
(`../tag_kb/references_ingest.py`) *links* a PDF to an existing SE row by DOI but
does not create missing ones; `../tag_kb/references_to_se.py` is the other half —
Crossref-enriches each DOI into an `Author+Year_Journal` key, de-dups against the
live SE DB, and creates the missing rows (dry-run by default; `--commit` to write).
After it runs, `refresh.sh` renders the new papers as nodes.

## NOT committed (regeneratable / secret)
- `vault/` — the 538 generated notes (gitignored; rebuild with the scripts)
- `.notion_token` — integration secret (gitignored)

## Notes
- The Notion API endpoint that works for this workspace: `/databases/{id}/query`
  with `Notion-Version: 2022-06-28` (the 8 DB ids are in `notion_to_obsidian.py`).
- KU-5.x collision: code lamellipodium KU-5.x → `H5-LP-x`; collective KB Unit 5 → `KB-5.x`.
  `add_code_nodes.py` maps code KU-5.x to H5-LP only for lamellipodium files.
- 3 paper DOIs had key/journal mismatches, resolved & noted in Notion
  (Kong2009 = JCB not Nature; Lindstrom2010 = Phys Rev E not Biomaterials;
  Friedl2011 = JCB 2010 not Cell 2011).

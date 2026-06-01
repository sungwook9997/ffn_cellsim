# Notion → Obsidian export (RAG v2 Contract-Graph mirror)

Read-only Obsidian mirror of the Notion **Contract Graph** (8-DB RAG v2).
**Notion is the source of truth**; this regenerates the vault from scratch.

## What it produces
A graph vault of **538 nodes / ~1550 edges, single connected component**:
literature (SourceEvidence) → KnowledgeClaim (KB-x.y) → ModelContract →
Parameter / ValidationGate → CodeMapping/RunResult/DecisionLedger, plus
claim↔claim dependency edges and Code/Test/Doc nodes from `ffn_sim/`.
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
python notion_to_obsidian.py   # Notion → vault (claims/papers/params/gates/contracts + relations)
python add_code_nodes.py       # repo code/tests/docs → linked to KB claims
open -a Obsidian ./vault
```

## Files (committed)
- `notion_to_obsidian.py` — main exporter (Notion API → markdown + frontmatter + wikilinks + graph.json)
- `ku_dependencies.py` — KU→KU (and KU-v2 KB→KB) dependency edges extracted from the legacy Unit pages' 🔗 DEPENDENCY blocks; turns the flat layered graph into a connected web
- `add_code_nodes.py` — scans `ffn_sim/` (.py / tests / docs) and links code files to the KB claims they cite

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

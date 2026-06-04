# SourceEvidence registration candidates — ACTIN ARCHITECTURE batch (2026-06-04)

PI-added cortex/actin-architecture papers (references/ top-level, added 2026-06-04 17:5x-18:0x +
11:1x), being read one-by-one (deep) by the Lead. Full study notes + 감상평 in
`docs/ACTIN_ARCHITECTURE_NOTES.md`. This file STAGES each for the KB pipeline.

**Pipeline (BATCH at session end, over all entries below — tools are full-rebuild, not per-paper):**
`references_ingest.py` (tag_corpus + BM25) → `verify_sources.py` (CrossRef/web, hallucination
hard-rule) → `references_to_se.py` (Notion SourceEvidence) → `link_se_claims_*.py` (→ KnowledgeClaim)
→ `refresh.sh` (Obsidian RAG). Then `tag_query.py` resolvable + graph node.

| # | citation_key | title | DOI / id | file | suggested KnowledgeClaim | in_corpus |
|---|---|---|---|---|---|---|
| 1 | Flormann2024_PNAS | Structure & mechanics of the cell cortex depend on location & adhesion state | 10.1073/pnas.2320372121 | flormann-et-al-2024-...pdf | cortex-architecture (mesh/thickness/bundling↔stiffness); adhesion-state dependence | no |
| 2 | Fritzsche2016_SciAdv | Actin kinetics shapes cortical network structure and mechanics | 10.1126/sciadv.1501337 | sciadv.1501337.pdf | cortex-architecture (bimodal exp filament length: Arp2/3 ~120nm + formin ~1200nm; long-formin = mechanical/percolation backbone) | no |

(rows appended as each paper is read)

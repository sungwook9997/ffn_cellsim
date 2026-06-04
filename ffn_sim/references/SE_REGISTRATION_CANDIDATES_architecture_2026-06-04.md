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
| 3 | Kim2007_MIT_thesis | Simulation of Actin Cytoskeleton Structure and Rheology (MS thesis) | (thesis) | 181655768-MIT.pdf | actin-network-simulation-methodology; crosslinker perpendicular(network)/parallel(bundle); connectivity/percolation as crosslinker-density output | no |
| 4 | Banerjee2021_JIISc | The Actomyosin Cortex of Cells: A Thin Film of Active Matter (review) | 10.1007/s41745-020-00220-2 | s41745-020-00220-2.pdf | cortex-active-gel-theory (continuum framing) | no |
| 5 | Chen2024_NatPhys | Energy partitioning in the cell cortex | 10.1038/s41567-024-02626-6 | s41567-024-02626-6.pdf | cortex-nonequilibrium (chem-vs-mech timescale competition) | no |
| 6 | SakamotoMurrell2024_CellRepPhysSci | Substrate geometry/topography induce F-actin reorganization & chiral alignment in adherent model cortex | 10.1016/j.xcrp.2024.... | 1-s2.0-S2666386424006520-main.pdf | unified-nucleator-axis (Arp2/3-branched=isotropic robust; formin-linear=aligned/spanning); reconstituted-model-cortex | no |
| 7 | Garlick2022_SciRep | Simple methods for quantifying super-resolved cortical actin (corrals) | 10.1038/s41598-022-06702-w | s41598-022-06702-w.pdf | cortex-mesh-size (corral 40-230nm cell-specific; A549 ~0.2µm²; mesh↑ on depolymerization); pore-analysis observable | no |

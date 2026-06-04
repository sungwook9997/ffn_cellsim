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
| 8 | LiGaoXu2022_BiophysJ | Network dynamics of the nonlinear power-law relaxation of cell cortex | 10.1016/j.bpj.2022.09.035 | 1-s2.0-S0006349522007858-main.pdf | cortex-rheology (power-law from DISORDER + crosslinker rebinding; disordered not lattice) | no |
| 9 | Fritzsche2017_NatCommun | Self-organizing actin patterns shape membrane architecture but not cell mechanics | 10.1038/ncomms14347 | ncomms14347.pdf | cortex-architecture (2 subpop: formin-long + Arp2/3-short 80%, 70° branch); patterns decoupled from bulk mechanics; Arp2/3-driven self-org | no |
| 10 | Bacher2021_FrontPhys | A 3D Numerical Model of an Active Cell Cortex in the Viscous Limit | 10.3389/fphy.2021.753230 | fphy-09-753230.pdf | cortex-active-gel continuum model (acceptance/framing; furrow/flows) | no |
| 11 | Ray2024_Development | Actin capping protein regulates actomyosin contractility to maintain germline architecture (C. elegans) | 10.1242/dev.... | Actin-capping_protein_regulates_actomyosin_contrac.pdf | filament-length control (capping↔length↔F-actin↔contractility); Cytosim adjacency | no |
| 12 | Chandrasekaran2024_NatCommun | Kinetic trapping organizes actin filaments within liquid-like protein droplets (VASP) | 10.1038/s41467-024-46726-6 | s41467-024-46726-6.pdf | peripheral: ABP-kinetics shapes network (kinetic trapping); VASP processive polymerase+bundler | no |
| 13 | MerinoCasallo2022_CellAdhMigr | Unravelling cell migration: defining movement from the cell surface (review) | 10.1080/19336918.2022.... | Merino-Casallo...2022...pdf | unified-framework framing: migration-mode↔actin-structure (mesenchymal/amoeboid/lobopodial); EMT/MAT | no |
| 14 | KadzikMunro2026_bioRxiv | Rapid actin filament turnover maintains cortical CONNECTIVITY while allowing cortex deformation and flow | 10.64898/2026.05.24.727551 | 2026.05.24.727551v1.full.pdf | ⭐ cortex-connectivity↔force-transmission (connectivity required for transmission/flow; length sets connectivity; balanced turnover maintains it) — DIRECT validation of γ-floor diagnosis | no |
| 15 | Serwas2021_NatCommun | Actin force generation in vesicle formation: cryo-ET insights (CME) | 10.1101/2021.06.28.450262 | 2021.06.28.450262v2.full.pdf (+media-1.pdf suppl) | peripheral: filament-membrane anchoring→force transmission (ERM/FA analog); cryo-ET architecture | no |

## ALL 15 READ (2026-06-04). Next: BATCH ingest (this turn: references_ingest=TAG/BM25;
## then verify_sources → references_to_se → link_se_claims → refresh.sh for Notion SE+Obsidian).

## Layer-2 (spheroid/cohesion/tissue) CROSS-LINKS — link these SE rows to Layer-2 KnowledgeClaims too
PI 2026-06-04: "layer2가 잘 쓸 수 있는 것들도 있다." When registering SE, attach to BOTH the
single-cell-cortex KC AND the Layer-2 KC below so the spheroid line can query them.

- Flormann2024_PNAS → +KC: layer2-cell-cortical-tension (cohesive/adhered γ feeds the σ bridge)
- MerinoCasallo2022_CellAdhMigr → +KC: layer2-invasion / EMT-MAT cancer phenotype (L2.7 invasion D1)
- KadzikMunro2026_bioRxiv → +KC: layer2-cortical-flow/cohesion (connectivity↔coherent flow at tissue scale)
- Fritzsche2016 / Flormann → cortical-tension magnitude that the single-cell γ → spheroid σ bridge consumes.

### NOTE — the σ-BRIDGE batch (earlier cortical-tension PDFs, the MOST Layer-2-relevant) also need
### Layer-2 KC links when registered: Winklbauer (jcs174623), Roffay (DownloadCombined...), Okuda
### (2026.03.17.712503), Fastabend (z152-x4l1), Chugh (emss-72183/ncb3525), Warmt, Nishitani,
### electrodeformation. These ARE the single-cell-γ→spheroid-σ bridge (D2; Winklbauer script done) →
### KC: layer2-surface-tension-bridge. Verify their SE status (some may already be linked by helper).

# SourceEvidence registration candidates — 2026-06-04

References-management pass over the **68 new top-level PDFs** that Syncthing dropped into
`ffn_sim/references/` since the 2026-06-03 TAG build. After triage + re-ingest the corpus is
**173 PDFs / 7402 chunks / 111 SE-linked**.

## What was done this session

1. **Noise quarantine.** 39 off-topic PDFs (Current Biology dispatches + neuro / ecology /
   behaviour / archaeogenetics / agronomy / virology / lipid-metabolism) were moved to
   `references/_offtopic/` — **not deleted**. That subdir is outside the ingest glob
   (top-level + `cellpress_bundle/` + `downloaded/` only), so the rebuild ignores it. These
   are the PI's general-subscription reading, not project literature. (Full move list at the
   bottom; `git mv`-style `mv -n`, reversible.)
2. **Manifest titles.** 22 manifest entries added (`references/analysis/_manifest.json`, n=50→71)
   so the relevant new papers ingest with clean titles instead of first-line garbage.
3. **Re-ingest.** `references_ingest.py` rebuilt `kb.duckdb` (24 exact byte-duplicates auto-
   skipped — the `downloaded/` Author-Year copies are byte-identical to the top-level files).
   Backup at `kb.duckdb.bak_20260604`.

All 22 below are now **BM25-searchable in TAG** but **absent from the Obsidian graph** and
**not citation-resolvable** because they have no Notion SourceEvidence row. Creating the SE
row also replaces the auto-slug `citation_key` with a proper Author-Year key.

> ⚠️ **Citation-integrity hard rule.** The 2026-06-02 audit found 3 hallucinated sources.
> Verify each DOI (CrossRef) before creating its SE row. SourceEvidence is edited in Notion
> (SoT), never here. DOIs below are extracted from the publisher PDFs, so verification is
> a confirm-not-discover step.

---

## A. Cortical-tension / cell-surface-mechanics cluster (H.3 cortex + γ-floor) — register first

These are the strongest hits for the active KU-3.5 γ-floor stack and the cortex unit.

| auto citation_key | Paper | DOI | Relevance |
|---|---|---|---|
| `adhesion-forces-and-cortical-tension-cou-b26dc0` | **Rübsam et al. 2017**, *Adhesion forces and cortical tension couple cell proliferation and differentiation to drive epidermal stratification* (Nat Cell Biol) | 10.1038/s41556-017-0005-z | ★ adhesion + cortical tension + **proliferation coupling** — directly the L2 spheroid premise |
| `differences-in-cortical-contractile-prop-9ed922` | **Warmt et al. 2021**, *Differences in cortical contractile properties between healthy epithelial and cancerous mesenchymal breast cells* (New J Phys) | 10.1088/1367-2630/ac254e | ★ breast-cell cortical contractility (MCF-class), γ overlay |
| `cell-adhesion-strength-from-cortical-ten-014ee5` | *Cell adhesion strength from cortical tension – an integration of concepts* (J Cell Sci) | 10.1242/jcs.174623 | cortical-tension ↔ adhesion-strength bridge |
| `balance-of-microtubule-stiffness-and-cor-e8ed0e` | **Dmitrieff et al. 2017**, *Balance of microtubule stiffness and cortical tension determines the size of blood cells with marginal band* (PNAS) | 10.1073/pnas.1618041114 | cortex vs MT-stiffness force balance |
| `cortical-tension-links-curvature-to-tiss-4fbfcd` | **Fastabend et al. 2026**, *Cortical tension links curvature to tissue growth in the cellular Potts model* (Phys Rev E) | 10.1103/z152-x4l1 | tension→growth coupling; dossier `analysis/51_*` already exists |
| `control-of-cellular-cortical-tension-and-b2df3b` | *Control of cellular cortical tension and shape by RhoGTPase signalling* (bioRxiv 2025) | 10.64898/2025.12.15.694413 | RhoA→cortical-tension control (preprint — flag) |
| `single-cell-mechanical-analysis-and-tens-ef4e49` | *Single-cell mechanical analysis and tension quantification via electrodeformation relaxation* (Phys Rev E) | 10.1103/physreve.103.032409 | single-cell tension measurement protocol |
| `quantification-of-single-cell-cortical-t-ac5ddb` | *Quantification of Single-Cell Cortical Tension Using Multiple Constriction Channels* (IEEE Sensors 2021) | 10.1109/jsen.2020.3048591 | cortical-tension measurement protocol |

## B. Membrane-tension / mechanotransduction cluster (H.8 membrane prep)

| auto citation_key | Paper | DOI | Note |
|---|---|---|---|
| `cell-protrusions-and-contractions-genera-c1b0e2` | *Cell protrusions and contractions generate long-range membrane tension propagation* (Cell 2023) | 10.1016/j.cell.2023.05.014 | membrane-tension propagation |
| `effective-membrane-tension-a-long-range--886bf5` | *Effective membrane tension: A long-range integrator of cellular dynamics* (Cell 2023) | 10.1016/j.cell.2023.05.033 | companion to above |
| `dissecting-cell-membrane-tension-dynamic-b24b55` | *Dissecting cell membrane tension dynamics and its effect on Piezo1-mediated mechanotransduction* (Nat Methods 2024) | 10.1038/s41592-024-02277-8 | membrane tension + Piezo1 |
| `molecular-dynamics-simulations-of-piezo1-6a6b19` | De Vecchis et al., *MD simulations of Piezo1 channel opening by increases in membrane tension* (Biophys J 2021) | 10.1016/j.bpj.2021.02.006 | ⚠️ filename was `mmc4.pdf` (supplement-style) — confirm it's the article not an SI bundle |
| `quantitative-analysis-of-cell-membrane-t-7a514c` | Nishitani & Miura, *Quantitative Analysis of Cell Membrane Tension in Time-Series Imaging and A Minimal Lattice Model* | **DOI not auto-extracted** | needs manual DOI lookup before SE row |

## C. Adhesion / cadherin / junction cluster (L2 cohesion + junction)

| auto citation_key | Paper | DOI | Note |
|---|---|---|---|
| `adhesion-induced-cortical-flows-pattern--324dc5` | *Adhesion-induced cortical flows pattern E-cadherin-mediated cell contacts* (Curr Biol 2023) | 10.1016/j.cub.2023.11.067 | ★ E-cadherin + cortical flow — L2.5 catch-bond context |
| `inferring-cell-junction-tension-and-pres-6d578e` | **Roffay et al.**, *Inferring cell junction tension and pressure from cell geometry* (Development) | 10.1242/dev.192773 | junction-tension inference; dossier `analysis/52_*` already exists |
| `cortical-tension-allocates-the-first-inn-3216d2` | *Cortical Tension Allocates the First Inner Cells of the Mammalian Embryo* (Dev Cell 2015) | 10.1016/j.devcel.2015.07.004 | tension-driven cell sorting |
| `target-cell-cortical-tension-regulates-m-486f2b` | *Target cell cortical tension regulates macrophage trogocytosis* (Nat Cell Biol 2025) | 10.1038/s41556-025-01807-6 | tension-gated adhesion (peripheral relevance) |

## D. Tumor-spheroid / tissue-mechanics simulation cluster (L2 spheroid line)

| auto citation_key | Paper | DOI | Note |
|---|---|---|---|
| `spontaneous-formation-of-tumor-spheroid--0b5a5c` | *Spontaneous formation of tumor spheroid on a hydrophilic filter paper for cancer stem cell enrichment* (Colloids Surf B 2019) | 10.1016/j.colsurfb.2018.11.038 | spheroid formation protocol |
| `tuning-cell-motility-via-cell-tension-wi-e0c5a7` | **Tuning Cell Motility via Cell Tension with a Mechanochemical Cell Migration Model** (Biophys J 2020) | 10.1016/j.bpj.2020.04.030 | published version |
| `tuning-cell-motility-via-cell-tension-wi-2f3e54` | ↑ **preprint** of the same paper | 10.1101/847046 | link to the same SE row as the published version, don't double-register |
| `3d-simulation-of-tissue-fb078b` | *3D Simulation of Tissue Mechanics with Cell Polarization* (ETH preprint 2023) | 10.3929/ethz-b-000653764 | vertex/tissue sim (preprint — flag) |
| `computational-design-for-engineering-lay-8c39ca` | *Computational Design for Engineering Layered Tissue Architectures via Cell–Cell Interfacial Tension Modulation* (bioRxiv 2026) | 10.64898/2026.03.17.712503 | interfacial-tension tissue design (preprint — flag) |

---

## Already in the corpus (NO action — byte/DOI duplicates of indexed papers)

- **Chugh et al. 2017** *Actin cortex architecture regulates cell surface tension* — new copies
  `emss-72183.pdf` / `ncb3525 (1).pdf`; already indexed as `Chugh2017_NatCellBiol`.
- **Belmonte et al.** *A theory that predicts behaviors of disordered cytoskeletal networks* —
  `msb.20177796.pdf`; already indexed.
- **Geiger/PLoS One** *Directed invasion of cancer cell spheroids inside 3D collagen matrices* —
  `journal.pone.0264571*.pdf` / `Directed_invasion_*.pdf`; already indexed as `Geiger2022_PLoSONE`.

## Off-topic — quarantined to `references/_offtopic/` (39 files, do NOT register)

Current Biology dispatches and general-science PDFs unrelated to cell mechanobiology
(odorant receptors, hippocampal circuits, archaeogenetics, spider silk, falcon cognition,
oocyte chromosome cohesion, teff harvest losses, viral-capsid nanoparticles, diamond-FET
biosensors, FADS2/aromatase lipid pharmacology, cerebellum, torpor, saccadic attention,
infant voice categorization, whirligig beetles, …). Restore from `_offtopic/` if any is
later judged in-scope.

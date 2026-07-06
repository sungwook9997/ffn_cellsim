# SourceEvidence registration candidates — 2026-07-06 (spheroid 3-zone necrosis / quiescence)

> ✅ **APPLIED 2026-07-06 (PI-signed off).** All 19 SourceEvidence rows below were registered
> into the Notion Contract-Graph and materialized into `kb.duckdb` (source_evidence 336→355).
> Every DOI was CrossRef-verified individually before registration; Lee2010 has no DOI (PMID
> 20514446 verified via PubMed); `Mazloomi2025_SciRep` carries a citation-integrity FLAG (the
> SYNTHESIS mis-labeled it "Sajjadi 2025"). Claims relations are **left unlinked** — the
> KnowledgeClaim / ModelContract / ValidationGate rows in the "Proposed downstream KB structure"
> section remain **PI-authored / not yet created**. PDF fetch into `references/` + a full
> `refresh.sh` (BM25 content layer) is the remaining follow-up.

Anchors for the **DCM 3-zone (proliferating / quiescent / necrotic) spheroid criterion** —
the nutrient/O₂-diffusion depth channel AND the mechanical compressive-stress channel that
together set the viable-rim / necrotic-core architecture. All values are the calibration
constants in
[`references/analysis/_necrosis/SYNTHESIS_dcm_necrosis_criterion.md`](../../references/analysis/_necrosis/SYNTHESIS_dcm_necrosis_criterion.md)
(§5) and the Q1/Q2/Q3 research notes, and are already wired into
[`dcm/dcm_necrosis_host.py`](../../dcm/dcm_necrosis_host.py) (3-zone depth updater).

**Why this manifest exists.** The necrosis *analysis + code* was completed 2026-06-11/23 but was
**never registered into the Notion Contract-Graph** — `source_evidence`, `knowledge_claim`,
`model_contract`, `validation_gate` all return **0** necrosis rows (verified 2026-07-06 via
`tag_query.py` structured layer + direct DuckDB scan). Only the raw PDF corpus (BM25 layer) held
a few of these papers, so `tag_query "spheroid necrotic core"` surfaced literature excerpts but
**no project claim/contract/gate** — the KB looked empty for necrosis. This manifest closes that
gap: it lists the SourceEvidence rows to register so the necrosis criterion becomes a first-class,
citable part of the Contract-Graph.

> ⚠️ **PI-gated.** Registering SourceEvidence rows (and any KnowledgeClaim / ModelContract /
> ValidationGate built on them) is a KB change → **PI sign-off before `harvest`/Notion apply**
> (CLAUDE.md hard rules). This file is the review manifest only; nothing is written to Notion yet.

> ⚠️ **Citation-integrity (2026-06-02 audit found 3 hallucinations).** Verify each DOI on CrossRef
> before registering. Dedup is by **`citation_key`**, not `doi` (existing SE rows store `doi=NULL`).

> ⚠️ **MCF-7 caveats (carry into `notes`):** no MCF-7-specific kPa threshold exists — the
> solid-stress bands (Helmlinger/Montel/Delarue) are explicitly **cell-line-independent**, so the
> 1/5/10 kPa band transfers. The quantitative necrotic-fraction-vs-size curve is **BT-474**
> (HER2⁺ sibling, Nieto 2026), not MCF-7; MCF-7 size landmarks are semi-quantitative (Lee 2010,
> Sajjadi 2025). These are **OVERLAY-ONLY calibration targets — never fit** (project hard rule).

---

## Already registered (do NOT re-create — link the necrosis claim to these)

| citation_key | SE id | note |
|---|---|---|
| `Nieto2026_EuropeanJournalOfPharmac` | 373120daec5d81ff8793da83452a6663 | BT-474 necrotic-fraction curve + glucose threshold ~0.08 mM. In corpus. **Verify its `claims` cover the necrosis usage.** |
| `Stylianopoulos2012_PNAS` | 372120daec5d81ee855fd5850a457ec1 | solid stress → vessel collapse → hypoxia (mechanics⇄nutrient coupling). |
| `Boot2021_AdvancesInPhysicsX` | 373120daec5d81358ee3ec88e0a9748e | viable rim ~100–200 µm; three-shell structure (KB paper 11). |

(`Mueller2017_Cell` is registered but is a **different** paper — NOT Mueller-Klieser 1997; register that one below.)

---

## A. Nutrient / O₂ diffusion channel — the viable-rim / necrotic-core geometry (register)

| paper | DOI | content (values + units) | in corpus? |
|---|---|---|---|
| **Thomlinson RH, Gray LH 1955**, *The histological structure of some human lung cancers and the possible implications for radiotherapy*, Br J Cancer 9(4):539–549 | 10.1038/bjc.1955.55 | cord radius **~200 µm** necrosis onset; **none <100 µm**; viable rim **≤180 µm**; O₂ diffusion **~150 µm**. Origin of the constant-viable-rim picture. | fetch |
| **Greenspan HP 1972**, *Models for the growth of a solid tumor by diffusion*, Stud Appl Math 51(4):317–340 | (no DOI — Wiley; verify) | three-radius moving-boundary model: **constant viable rim, expanding necrotic core** (the `R_nec = R_cluster − L_viable` invariant used in code). | fetch |
| **Mueller-Klieser W 1997**, *Three-dimensional cell cultures: from molecular mechanisms to clinical applications*, Am J Physiol 273(4):C1109–23 (PMID 9357753) | 10.1152/ajpcell.1997.273.4.C1109 | viable rim **~100–200 µm**; rim ~constant, core grows. | fetch |
| **Grimes DR et al. 2014**, *A method for estimating the oxygen consumption rate in multicellular tumour spheroids*, J R Soc Interface 11:20131124 (PMID 24430128) | 10.1098/rsif.2013.1124 | diffusion limit **232±22 µm**; OCR 7.29±1.4×10⁻⁷ m³ kg⁻¹ s⁻¹; D_O₂ ~3.8×10⁻⁹ m²/s; r_l 233 / r_n 155 µm (HCT116). Anchors `L_pen`. | fetch |
| **Grimes DR, Bull DM et al. 2016**, PLoS One 11(4):e0153692 | 10.1371/journal.pone.0153692 | **D = 2×10⁻⁹ m²/s**; mitotic O₂ cutoff **p_m ~0.5 mmHg** (std glucose) / ~5 mmHg (low glucose). route-B O₂ field boundary. | fetch |
| **Jiang Y et al. 2005** (multiscale tumor-spheroid model) | (no DOI in note — verify; Biophys J 89:3884?) | necrosis thresholds **O₂ <0.02 mM, glucose <0.06 mM, lactate >8 mM**. route-B death thresholds. | fetch |
| **Lee SY et al. 2010**, Oncol Rep 24(1):73–79 (PMID 20514446) | (no DOI — verify) | **MCF-7 ~700 µm / day 8** necrotic core; none at day 6; three-zone zonation. MCF-7 size landmark. | fetch |
| **Sajjadi et al. 2025** (PMC12638756) | (verify DOI) | MCF-7 viability **95→92→75→50 %** (d2/7/10/14); necrotic layer ~day 6 (5 000-cell); calcein-AM/EthD-1 2/4 µM. viability time-course overlay. | fetch |
| **Sutherland RM 1988**, *Cell and environment interactions in tumor microregions: the multicell spheroid model*, Science 240:177–184 | 10.1126/science.2451290 (verify) | canonical MCTS three-zone reference (cited in `necrosis_3zone.md`). | fetch |
| **Front Oncol 6:105 2016** | 10.3389/fonc.2016.00105 | MCF-7 **<0.5 mm no necrosis, 1.8 mm clear core**. onset-diameter landmark. | fetch |
| **RSC Analyst 2020** | 10.1039/d0an00979b | MCF-7 OCR-per-volume falls with radius (consumption-limited diffusion). | fetch |

## B. Mechanical / compressive-solid-stress channel (register)

| paper | DOI | content (values + units) | in corpus? |
|---|---|---|---|
| **Helmlinger G, Netti PA, Lichtenbeld HC, Melder RJ, Jain RK 1997**, *Solid stress inhibits the growth of multicellular tumor spheroids*, Nat Biotechnol 15(8):778–783 (PMID 9255794) | 10.1038/nbt0897-778 | inhibiting stress **45–120 mmHg = 6–16 kPa**, **line-independent**. Anchors `σ_death` band. | fetch |
| **Montel F et al. 2011**, *Stress clamp experiments on multicellular tumor spheroids*, Phys Rev Lett 107(18):188102 (PMID 22107677) | 10.1103/PhysRevLett.107.188102 | **~10 kPa** drastically reduces growth (loss mainly in the CORE); 5 kPa measurable. `σ_q_full`/`σ_death`. | fetch |
| **Delarue M et al. 2014**, *Compressive stress inhibits proliferation in tumor spheroids through a volume limitation*, Biophys J 107(8):1821–1828 (PMID 25418163) | 10.1016/j.bpj.2014.08.031 | compression → **~5 % cell-volume drop at 5 kPa in ~30 min** → p27^Kip1 → late-G1 arrest; conserved across 5 lines; **reversible** (the volume checkpoint, `ΔV/V₀` read-out). | fetch |
| **Cheng G, Tse J, Jain RK, Munn LL 2009**, *Micro-environmental mechanical stress controls tumor spheroid size and morphology*, PLoS ONE 4(2):e4632 (PMID 19247489) | 10.1371/journal.pone.0004632 | high-stress regions suppress proliferation + induce mitochondrial apoptosis (dormancy). | fetch |
| **Dolega ME et al. 2017**, *Cell-like pressure sensors reveal increase of mechanical stress towards the core of multicellular spheroids under compression*, Nat Commun 8:14056 (PMID 28128198) | 10.1038/ncomms14056 | internal pressure **RISES toward the core** from radial cell anisotropy alone → justifies **letting σ emerge**, not imposing σ(r). | fetch |
| **McGrail DJ … Dawson MR 2015**, Biophys J 109(7):1334–1337 (PMID 26445434) | 10.1016/j.bpj.2015.07.046 | survival under solid stress = active **Na⁺/NHE1 osmoregulation**; **breast lineage validated**; death = osmoregulation failure. | fetch |
| **Hadjigeorgiou AG, Stylianopoulos T 2023**, Biomech Model Mechanobiol 22:1625–1643 (PMID 37129689) | 10.1007/s10237-023-01716-3 | real-tumor residual stress **3.31–10.88 kPa** — in-vitro thresholds are physiological (sanity band). | fetch |
| **Nia HT et al. 2017**, Nat Biomed Eng 1:0004 (PMID 28966873) | 10.1038/s41551-016-0004 | solid stress increases with tumor size; confinement a major contributor. | fetch |

---

## Proposed downstream KB structure (PI-authored — flagged, not created)

Once the SE rows above are registered, the necrosis criterion should become citable structure.
**These are PI-authored; this manifest only flags them:**

1. **KnowledgeClaim(s)** — e.g.
   - `KB-nec.1` O₂/nutrient diffusion sets a **constant viable rim ~150 µm** (100–200); necrotic
     core = depth ≥ L_viable → cells with `r < R_cluster − 150 µm` necrose. [A-anchors]
   - `KB-nec.2` proliferation is **rim-confined ~40 µm** (Ki-67⁺ outer 2–3 layers). [A]
   - `KB-nec.3` **compressive solid stress** arrests proliferation (onset ~1 kPa, full ~5 kPa) and
     kills on sustained ~10 kPa (6–16 kPa line-independent); volume-checkpoint −5 %/5 kPa. [B]
   - `KB-nec.4` MCF-7 necrotic-core onset **~400–500 µm diameter**, prominent **~700 µm/day 8**,
     all-viable ≤300 µm. [A, overlay-only]
2. **ModelContract** — `MC-necrosis-3zone`: the worst-of-two per-cell state updater
   (SYNTHESIS §2/§3) binding KB-nec.1–4 to `dcm/dcm_necrosis_host.py`.
3. **ValidationGate** — `VG-necrosis-viable-rim`: viable rim stays ~150 µm ±band as the core
   grows (Greenspan invariant); necrotic-fraction-vs-diameter tracks the BT-474 overlay curve
   (0.20→0.29→0.61). Overlay-only, never fit.
4. **CodeMapping** — put the KU/`KB-nec.*` ids in the `dcm_necrosis_host.py` docstring head so
   `harvest_ops` auto-links code→contract on the next sync.

---

## Apply procedure (after PI sign-off)

1. Fetch the "fetch" PDFs into `references/` (gbook KAIST full-text for paywalled), CrossRef-verify each DOI.
2. Add SE rows in **Notion** (SoT) — never edit `kb.duckdb`/vault directly.
3. PI authors KnowledgeClaim/ModelContract/ValidationGate rows (§ above).
4. `bash outputs/tag_kb/refresh.sh` → re-materialize DuckDB + Obsidian; `verify_sources.py` for verdicts.
5. Re-run `tag_query.py "spheroid necrotic core"` — structured layer should now return KB-nec.* + MC/VG.

---
id: Q2_mcf7_livedead
topic: MCF7 spheroid LIVE/DEAD calibration data — viable rim, necrotic-core fraction, necrosis-onset diameter, time course
purpose: EXAMPLE / calibration dataset for a DCM (dead-cell / death-criterion model) necrosis & quiescence rule
date: 2026-06-11
status: research note (literature-first; overlay-only, never fit-to per project hard rules)
---

# Q2 — MCF-7 Spheroid LIVE/DEAD Calibration Data

Goal: published MCF-7 (and close breast-cancer) spheroid studies that image LIVE/DEAD
(calcein-AM / ethidium homodimer-1, or PI) and report **(i) viable-rim thickness vs
diameter, (ii) necrotic-core fraction vs size, (iii) the diameter at which a dead core
first appears, (iv) the time course of necrosis.** Numbers below are the actual measured
values, with units and source. MCF-7-specific numbers are flagged; canonical
breast/biophysical anchors fill the gaps where MCF-7 imaging is only qualitative.

⚠️ Per project hard rules: these are **overlay/sanity anchors, literature-first** — used to
write a necrosis/quiescence *criterion contract before the run*, NOT to fit the model to.

---

## 1. The headline numbers (MCF-7 first)

| Quantity | Value (units) | Cell line | Source |
|---|---|---|---|
| **Necrotic core first appears** (no core at day 6, present by day 8–9) | spheroid ≈ **700 µm** diameter, **day 8** culture | **MCF-7** | Lee 2010 |
| Diameter at which a "prominent necrotic core" is seen | **> 700 µm** | MCF-7 (general statement) | review/web synthesis; consistent w/ Lee 2010 |
| Necrotic layer first detected (5 000-cell spheroid) | **≈ day 6** | **MCF-7** (SpheroidSync, 5 000 cells/well) | Sajjadi 2025 (SpheroidSync) |
| Whole-spheroid viability (calcein/EthD-1) | **95±3% (d2) → 92±2% (d7) → 75±4% (d10) → 50±5% (d14)** | **MCF-7** (SpheroidSync method) | Sajjadi 2025 |
| Spheroid diameter at day 10 vs seeding (500–20 000 cells/well) | **420 → 800 µm** | **MCF-7** | Froehlich 2016 / co-culture optimization lit. |
| "Optimal"/well-formed spheroid size used for assays | **≈ 500 µm** | MCF-7 | same |
| Small spheroid, no necrosis | **< 0.5 mm**, ≈ day 4 | MCF-7 | Froehlich/Frontiers MCF-7 microimaging 2016 |
| Large spheroid, clear necrotic core | **1.8 mm**, ≥ day 10 | MCF-7 | Frontiers MCF-7 microimaging 2016 |
| Three-zone structure confirmed (necrotic core / quiescent mid / proliferating rim) | day-8 spheroid | MCF-7 | Lee 2010; Frontiers 2016 |

### Canonical breast / biophysical anchors (use for the rim & threshold the MCF-7 papers leave qualitative)

| Quantity | Value (units) | Cell line / context | Source |
|---|---|---|---|
| **Oxygen diffusion–limited viable depth** (the rim cells can keep alive) | **≈ 150–200 µm** (commonly quoted 100–200 µm) | general MCTS | Mueller-Klieser 1997 (review); textbook consensus |
| Diameter above which a hypoxic/necrotic core develops | **≈ 400–500 µm** (necrosis onset window); a *visible* core typically by ~500 µm | general MCTS / breast | review synthesis; Curcio-type lit. |
| **Measured O₂ diffusion limit** (rim thickness, stained sections) | **232 ± 22 µm** | DLD1 spheroids (validated method, transferable rim scale) | Grimes 2014 |
| O₂ consumption rate (OCR) | **7.29 ± 1.4 × 10⁻⁷ m³ kg⁻¹ s⁻¹** | DLD1; method validated against stained sections | Grimes 2014 |
| O₂ diffusion constant used | **D = 2 × 10⁻⁹ m² s⁻¹** (≈ water) | spheroid models | Grimes 2014; Bull/Byrne 2016 (PLOS One) |
| Critical O₂ partial pressure for mitosis (proliferation cutoff → quiescence) | **p_m ≈ 0.5 mmHg** (std glucose); **≈ 5 mmHg** (low glucose) | multi-line (HCT116, LS174T, MDA-MB-468, …) | Grimes/Bull 2016 (PLOS One) |
| Glucose necrosis threshold | **< 0.06 mM** (Jiang); **≈ 0.08 mM** back-fit for BT-474 | EMT6 / BT-474 breast | Jiang 2005; Nieto 2026 (KB paper 02) |
| O₂ necrosis threshold | **< 0.02 mM** | EMT6 multiscale | Jiang 2005 (via KB paper 02) |
| MCF-7 OCR per volume **decreases** with radius (core O₂ starvation) | qualitative, radius-dependent | **MCF-7** | RSC Analyst 2020 (OCR during necrotic-core formation) |

### Breast-cancer sibling line (BT-474, HER2+) — quantitative necrotic-core-fraction time course (KB paper 02)

This is the most quantitative *necrotic-fraction-vs-size* breast dataset in the project KB,
imaged with **calcein-AM / PI (CLSM)**. Use as the closest calibratable necrotic-fraction
curve when MCF-7 fraction-vs-size data are unavailable.

| Day | Diameter (mm) | Necrotic-core volume fraction (experimental) | Porosity |
|---|---|---|---|
| 2 | 0.155 | none (first core day 3+) | 0.0145 |
| 3 | 0.260 | core begins | 0.0543 |
| 5 | 0.504 | **0.200** | 0.1149 |
| 6 | 0.842 | **0.291** | 0.4443 |
| 7 | 0.863 | **0.607** | 0.4676 |

(BT-474, 2 000 cells/well, liquid overlay; Nieto et al. 2026, EJPS, KB analysis paper 02.)

---

## 2. LIVE/DEAD staining protocol (as used on these spheroids)

- **Calcein-AM (live, green)** + **Ethidium homodimer-1 (EthD-1, dead, red)** at **2 µM /
  4 µM** in PBS (Sajjadi 2025; standard ULA-plate protocol). Viable = green, non-viable = red.
- **Calcein-AM / PI** by confocal (CLSM) for size + necrotic core (Nieto 2026 BT-474).
- MCF-7 spheroids of **100 µm and 300 µm** diameter stained calcein/EthD after 48 h show
  uniformly viable cells (no dead core at this size) — confirms the rim-only regime below
  the diffusion limit.
- **Acridine orange / PI** also used for MCF-7 spheroid apoptosis/necrosis morphology
  (Shujaa Edin 2021; Ho 2021) — qualitative, not zonal.

---

## 3. Synthesis — the picture for a DCM criterion

1. **Below the diffusion limit (~150–200 µm radius / ≤ ~300–400 µm diameter):** spheroid is
   essentially all-viable (calcein+ throughout). MCF-7 at 100 & 300 µm = no dead core.
2. **Necrosis onset window ~400–500 µm diameter** (general MCTS); for **MCF-7** specifically
   a *necrotic layer* is reported ≈ **day 6** (5 000-cell spheroid, Sajjadi 2025) and a
   **frank necrotic core at ~700 µm / day 8** (Lee 2010, no core at day 6).
3. **Above onset, the viable rim thickness stays roughly constant (~150–230 µm)** while the
   **necrotic core grows ∝ spheroid size** — the canonical Mueller-Klieser / Grimes result.
   So necrotic-core *fraction* rises steeply with diameter (BT-474: 0.20 → 0.29 → 0.61 over
   days 5→7, i.e. ~0.5 → ~0.86 mm).
4. **Three radial zones:** proliferating outer rim (high O₂, p_O₂ > ~0.5 mmHg) → quiescent
   middle (O₂ below mitotic threshold but above death) → necrotic core (O₂ < ~0.02 mM /
   glucose < ~0.06 mM).
5. **Time course (MCF-7 whole-spheroid viability):** ~95% (d2) → ~92% (d7) → 75% (d10) →
   50% (d14) — necrosis is progressive once the core forms.

---

## 4. How to turn this into a DCM necrosis / quiescence criterion

A geometry-/field-driven rule that reproduces the calibration data above:

- **Quiescence criterion (stop proliferation):** flag a cell quiescent when local O₂ falls
  below the **mitotic threshold p_m ≈ 0.5 mmHg** (std glucose) — equivalently, deeper than
  the **proliferating-rim depth ~150–230 µm** from the surface. This carves the outer
  proliferating rim from the inner quiescent shell.
- **Necrosis criterion (mark dead):** flag a cell necrotic when local O₂ < **~0.02 mM**
  AND/OR glucose < **~0.06 mM** (Jiang thresholds), i.e. once it lies beyond the
  diffusion-supported viable depth. With **D = 2×10⁻⁹ m² s⁻¹** and OCR
  **~7.3×10⁻⁷ m³ kg⁻¹ s⁻¹** (Grimes), a 1-D/radial reaction–diffusion solve predicts the
  anoxic radius; cells inside it → dead.
- **Calibration / sanity targets (overlay, do NOT fit):**
  - No dead core for spheroid diameter ≲ 300–400 µm (MCF-7 100/300 µm all-viable).
  - Necrotic core appears at **~400–500 µm onset / ~700 µm prominent** (≈ MCF-7 day 6–8).
  - Viable rim stays **~150–230 µm** as the core grows.
  - Necrotic-core fraction-vs-diameter tracks the BT-474 curve (0.20 @ 0.5 mm → 0.61 @ 0.86 mm).
  - Whole-spheroid viability time course ~95→92→75→50% over days 2→7→10→14 (MCF-7).
- **Simplest mesoscale proxy (if no field solve):** a **constant-rim / radial-depth rule** —
  cell viable if within R_rim ≈ 200 µm of the surface (proliferating if within ~150 µm,
  quiescent 150–200 µm, necrotic beyond). This single geometric constant reproduces the
  constant-rim + growing-core observation directly and is the least-parameter criterion.
- **Physiological-baseline note (project hard rule):** if metabolism/O₂ is ever added, set
  O₂/glucose boundary at real culture-medium values and start the cell at its physiological
  setpoint — do not start from a null field and expect emergent agreement.

Caveats: MCF-7 *necrotic-fraction-vs-size* is reported only semi-quantitatively; the
quantitative fraction curve is BT-474 (HER2+ sibling, same calcein/PI imaging). MCF-7
spheroids are known to loosen/partly decompose by ~7–10 days at some seeding densities
(KB paper 03), which confounds late-timepoint viability. Diffusion-limit (~150–230 µm) is
cross-line robust; absolute onset diameter depends on OCR, glucose, and seeding density.

---

## 5. Sources (cited)

**MCF-7 specific**
- **Lee SY, Jeong EK, Jeon HM, Kim CH, Kang HS (2010).** "Implication of necrosis-linked p53
  aggregation in acquired apoptotic resistance to 5-FU in MCF-7 multicellular tumour
  spheroids." *Oncol Rep* 24(1):73–79. PMID 20514446. — MCF-7 spheroids ~**700 µm at day 8
  with necrotic core**; **no necrotic core at day 6**, core present by day 9; explicit
  proliferating-outer / quiescent-inner / necrotic-core zonation. (No DOI in PubMed record.)
- **Sajjadi et al. (2025).** "SpheroidSync … uniform and robust MCF7 spheroids in 3D culture."
  *PMC12638756*. — Calcein-AM/EthD-1 (2 µM/4 µM); MCF-7 viability **95±3%(d2)/92±2%(d7)/
  75±4%(d10)/50±5%(d14)**; **necrotic layers ~day 6** for 5 000-cell spheroids; seeding
  700–12 000 cells. URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC12638756/
- **Frontiers in Oncology MCF-7 microimaging (2016).** "Metabolic Study of Breast MCF-7
  Tumor Spheroids after Gamma Irradiation by ¹H NMR Spectroscopy and Microimaging."
  *Front Oncol* 6:105. DOI 10.3389/fonc.2016.00105. PMC4848320. — MCF-7 **< 0.5 mm ≈ day 4
  (no necrosis)** vs **1.8 mm ≥ day 10 (clear necrotic core)**; layer-like necrotic
  core + viable rim (quiescent inner + proliferating outer). Seeding 4×10⁶ cells/tube.
- **RSC *Analyst* (2020).** "Oxygen consumption rate of tumour spheroids during
  necrotic-like core formation." DOI 10.1039/d0an00979b. — **MCF-7** OCR per spheroid volume
  **decreases with radius** (core O₂ starvation); MCF-7 spheroids with/without necrotic core.
- **MCF-7 100 & 300 µm spheroids, calcein/EthD at 48 h, all-viable** — PLOS One microtumor /
  ULA live/dead figure sets (research-figure synthesis); consistent with diffusion-limit rim.

**Breast sibling (quantitative necrotic-fraction curve)**
- **Nieto C, González-Garcinuño Á, Martín del Valle E (2026).** "Integrating simulation and
  experimental validation of nutrient-limited growth in breast cancer spheroids."
  *Eur J Pharm Sci* 216:107370. DOI 10.1016/j.ejps.2025.107370. (KB analysis paper 02.) —
  **BT-474** calcein/PI necrotic-core fraction **0.200(d5)/0.291(d6)/0.607(d7)**;
  diameters 0.50/0.84/0.86 mm; **necrosis glucose threshold ~0.08 mM**; first core day 3+.

**Canonical biophysical anchors**
- **Grimes DR, Kelly C, Bloch K, Partridge M (2014).** "A method for estimating the oxygen
  consumption rate in multicellular tumour spheroids." *J R Soc Interface* 11(92):20131124.
  DOI 10.1098/rsif.2013.1124. PMID 24430128. — **diffusion limit 232 ± 22 µm**, **OCR
  7.29 ± 1.4 × 10⁻⁷ m³ kg⁻¹ s⁻¹**, derives necrotic-core / hypoxic-region / proliferating-rim
  extents validated on stained DLD1 sections. (Rim scale transferable.)
- **Grimes DR et al. / Bull (2016).** "The Role of Oxygen in Avascular Tumor Growth."
  *PLOS One* 11(4):e0153692. DOI 10.1371/journal.pone.0153692. — **D = 2×10⁻⁹ m² s⁻¹**;
  mitotic O₂ cutoff **p_m ≈ 0.5 mmHg** (std glucose) / **~5 mmHg** (low glucose); multi-line
  OCR incl. MDA-MB-468 ~23.4 mmHg/s.
- **Mueller-Klieser W (1997).** "Three-dimensional cell cultures: from molecular mechanisms
  to clinical applications." *Am J Physiol* 273(4):C1109–23. DOI 10.1152/ajpcell.1997.273.4.C1109.
  PMID 9357753. — review: **O₂-diffusion-limited viable rim ~100–200 µm**; necrotic core
  grows with size while rim ≈ constant.
- **Jiang Y et al. (2005).** multiscale spheroid model — necrosis thresholds **O₂ < 0.02 mM,
  glucose < 0.06 mM, lactate > 8 mM** (via KB paper 02).
- **faCellitate application note.** "Detection of Necrotic Cores in Tumor Spheroids Using
  ULA Plates" — live/dead protocol context (calcein/EthD), necrotic-core appearance with size.

PubMed attribution: article metadata above retrieved from **PubMed**; DOIs given as links in
the source list (Lee 2010 has no DOI in the PubMed record).

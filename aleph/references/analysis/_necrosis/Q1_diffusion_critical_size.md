# Q1 — Diffusion-limited necrotic core vs spheroid size

**Question:** What O2/nutrient diffusion penetration depth sets the viable-rim thickness
(Greenspan / reaction-diffusion picture)? What is the critical spheroid diameter for
necrotic-core onset, the viable-rim thickness (µm), and the proliferating / quiescent /
necrotic zone widths — generally and for breast/MCF7 spheroids?

Created 2026-06-11. Tools: ffn_sim KB (papers 02,04,06,16,21,27), PubMed, WebSearch.

---

## 1. The Greenspan / Thomlinson-Gray picture (canonical numbers)

The classical reaction-diffusion ("diffusion-limited growth") picture: O2 diffuses in
from the spheroid surface and is consumed by respiring cells (Michaelis-Menten / zero-
order). The balance of **diffusion in** vs **consumption** sets a finite penetration
depth. Cells beyond that depth fall below a survival O2 threshold → quiescence then
necrosis. Result is a 3-layer "compound sphere": outer proliferating rim → intermediate
quiescent (viable but arrested) shell → central necrotic core.

- **O2 diffusion length in tumour tissue: ~100–200 µm** (most-cited single value ~150 µm,
  sometimes 160 µm). This is the metabolic-consumption-limited distance from a capillary/
  surface, not the free-diffusion length. (Thomlinson & Gray 1955; Hlatky/Hall texts.)
- **Thomlinson & Gray (1955)** histology of human bronchial carcinoma cords (the seminal
  diffusion-limited-hypoxia paper):
  - **No tumour cord > ~200 µm RADIUS lacks central necrosis** → necrosis onset at
    cord radius ≈ 200 µm (diameter ≈ 400 µm for a cord; corresponds to ~150–180 µm of
    viable tissue from the O2 source).
  - **No central necrosis in any cord < ~100 µm radius.**
  - **Thickness of viable tumour rim never exceeds ~180 µm** (range ~100–180 µm).
- Greenspan (1972) formalized this into the moving-boundary ODE model with three radii:
  outer radius R, quiescent-onset radius, necrotic radius — the **viable rim
  (R − R_necrotic) tends to a constant** as the spheroid grows; only the necrotic core
  expands. This "constant viable rim" is the key qualitative invariant.

## 2. Critical diameter for necrotic-core onset (spheroid culture)

Stage progression with increasing spheroid diameter (well-replicated across cell lines;
breast-relevant sources below):
- **< ~200 µm diameter:** fully viable/proliferative, no gradients.
- **~200 µm diameter:** **hypoxic core begins** to form (O2 gradient detectable; HIF-1α,
  pimonidazole positivity in centre).
- **~200–300 µm:** typical **zonation establishes** — proliferative outer shell +
  normoxic quiescent middle + hypoxic centre.
- **~400–600 µm:** classic threshold band for **necrotic-core appearance** in many MCTS.
- **> ~500 µm diameter:** **necrotic core robustly detectable** (the most commonly quoted
  single critical-diameter figure for a frank necrotic core in cell-line spheroids).
- **~600–800 µm:** growth plateau as diffusion limit caps the viable shell (ffn KB paper 16).

So: hypoxia onset ~200 µm diameter, necrosis onset ~400–600 µm diameter (use ~500 µm as
the round single number; ~200 µm RADIUS for tumour cords à la Thomlinson-Gray).

## 3. Viable-rim & zone widths (µm)

- **Viable rim thickness: ~100–200 µm** (Greenspan-limit). ffn KB paper 16 (review):
  viable rim 100–200 µm. ffn KB paper 04 (hybrid agent-based+RD model): **stabilized
  viable rim ≈ 100 µm, approximately CONSTANT as the spheroid grows; only the necrotic
  core expands** — direct numeric confirmation of the Greenspan invariant.
- **Proliferating outer rim:** outermost ~1–3 cell layers up to the diffusion limit;
  drug/penetration assays treat the outer **~20 µm (2–3 cell layers)** as the highest-
  exposure shell (ffn KB paper 16). Proliferative shell width is the part of the viable
  rim above the quiescence O2 threshold — order tens of µm.
- **Quiescent (G1/G0-arrested) intermediate shell:** the rest of the viable rim between
  proliferative shell and necrotic boundary; in FUCCI/continuum models it is a comparable
  fraction of the rim (ffn KB paper 06 melanoma: three zones — necrotic |x|<0.6,
  G1-arrested 0.6<|x|<1.1, proliferating 1.1<|x|<1.3 in normalized radius).
- **Necrotic core:** everything inside R_necrotic; grows with spheroid size while the rim
  stays ~constant.

### Concrete model-fitted geometry (Grimes/HCT116, ffn KB paper 27)
- **O2 diffusion-limit radius r_l = 233 µm**, **necrotic-core radius r_n = 155 µm** →
  viable rim ≈ r_l − r_n ≈ **78 µm** at hypoxia threshold 11 mmHg, medium pO2 ≈ 100 mmHg.
- **O2 diffusion coefficient D_O2 ≈ 3.8×10⁻⁹ m²/s** (fit); consumption 22.1 mmHg/s.

## 4. Breast / MCF7-specific numbers

- **MCF-7 spheroids:** hypoxic core forms > ~200 µm diameter; **necrotic core detectable
  > ~500 µm diameter**; MCF-7 necrosis attributed to hypoxia + nutrient (glucose)
  deficiency in the central region (oncotarget/lit; cf. WebSearch §). T-47D (breast)
  spheroids show a circular pimonidazole-hypoxic ring adjacent to the necrotic core by day 8.
- **MCF-7 spheroid microstructure (ffn KB paper 15, MCF-7/MDA-MB-231):** MCF-7 spheroids
  are dense/compact with **intercellular gaps 5–10 µm wide** (SEM), interpreted as
  facilitating nutrient/O2 diffusion. Cell diameter ~15 µm. Viability healthy to day 3,
  PI+ death by day 6; HIF-1α absent until day 3, marked by day 6 (ffn KB paper 07,
  MCF-7/MDA-MB-231 spheroids). 3D MCF-7 enriched in quiescent G1/G0 vs 2D.
- **BT-474 (HER2+ breast, ffn KB paper 02, COMSOL RD + experiment):** necrotic core first
  appears **day 3** (none days 0–2), grows ~exponentially; necrotic-core volume fraction
  d5/d6/d7 = 0.20/0.29/0.61. **Necrosis glucose threshold ~0.08 mM** (cf. Jiang O2
  <0.02 mM, glucose <0.06 mM, lactate >8 mM). D_O2 ≈ 1.5×10⁻⁹ m²/s, D_glucose ≈
  5.9×10⁻¹⁰ m²/s in DMEM; cell d = 15 µm; O2 internalization k_int = 2.5×10⁻¹⁸ mol/cell/s.
- **General breast O2-survival thresholds:** hypoxia/HIF-driving O2 tension < 5–10 mmHg
  (ffn KB paper 21); critical O2 for arrest k_O2 ≈ 4.6×10⁻³ mM (ffn KB paper 04);
  necrosis O2 < 0.02 mM (Jiang multiscale, via ffn KB paper 02).

## 5. Key constants table (for DCM use)

| Quantity | Value | Source |
|---|---|---|
| O2 diffusion length (tissue, consumption-limited) | ~100–200 µm (≈150 µm) | Thomlinson-Gray 1955; lit |
| Necrosis onset, tumour cord RADIUS | ~200 µm (none < ~100 µm) | Thomlinson-Gray 1955 |
| Hypoxic-core onset, spheroid DIAMETER | ~200 µm | lit (MCF-7/MCTS) |
| Zonation established, DIAMETER | ~200–300 µm | lit |
| Necrotic-core onset, spheroid DIAMETER | ~400–600 µm (use ~500 µm) | lit; MCF-7 >500 µm |
| Growth plateau, DIAMETER | ~600–800 µm | ffn KB paper 16 |
| Viable-rim thickness | ~100–200 µm (≈100 µm, ~constant) | ffn KB papers 04,16; Greenspan |
| Viable rim (Grimes/HCT116 fit) | r_l 233 − r_n 155 ≈ 78 µm | ffn KB paper 27 |
| Outer high-exposure proliferative shell | ~20 µm (2–3 cell layers) | ffn KB paper 16 |
| Hypoxia O2 threshold | < 5–11 mmHg (~0.013–0.02 mM) | ffn KB papers 21,27,04 |
| Necrosis O2 / glucose threshold | O2 <0.02 mM / glu <0.06–0.08 mM | Jiang 2005; ffn KB paper 02 |
| O2 diffusion coeff (medium/tissue) | 1.5–3.8 ×10⁻⁹ m²/s | ffn KB papers 02,27 |
| Cell diameter (MCF-7/BT-474) | ~15 µm | ffn KB papers 02,15 |

## 6. How to turn this into a DCM necrosis/quiescence criterion

Two equivalent routes, both Greenspan-consistent:

**(A) Depth-from-surface rule (cheap, no PDE).** For each cell compute distance d from the
spheroid surface (= R − r_cell). Assign state by depth thresholds calibrated to the
viable-rim geometry:
- d < d_prolif (~20–50 µm) → **proliferating**
- d_prolif ≤ d < L_viable (L_viable ≈ 100–150 µm; ≈ the O2 diffusion length) → **quiescent**
- d ≥ L_viable → **necrotic**.
Equivalently, define the necrotic-core radius R_nec = max(0, R − L_viable): cells with
r < R_nec are necrotic. This reproduces "constant viable rim, expanding core" for free and
gives necrosis onset automatically once R > L_viable (i.e. diameter > ~200–300 µm), matching
Thomlinson-Gray cord radius ~200 µm if L_viable set near 150–180 µm.

**(B) O2 reaction-diffusion field (mechanistic, matches no-abstraction preference).** Solve
steady ∇·(D∇c) = ρ·U(c) on the aggregate (D ≈ 1.5–3.8×10⁻⁹ m²/s, surface c = medium pO2
~100 mmHg / ~0.2 mM dissolved, U Michaelis-Menten with cell-density weighting). Then:
- c < c_quiescent (k_O2 ≈ 4.6×10⁻³ mM / ~5–11 mmHg) → quiescent (stop growth/division)
- c < c_necrosis (~0.02 mM O2, and/or glucose < ~0.06–0.08 mM) → necrotic.
Add glucose as a second co-limiting field for breast lines (BT-474/MCF-7 necrosis is
glucose-threshold-driven at ~0.08 mM, ffn KB paper 02).

For a first DCM pass use **(A)** with L_viable = 150 µm and a ~40 µm proliferative shell;
upgrade to **(B)** when an explicit nutrient field is wired in. Calibrate L_viable so
necrotic-core volume fraction vs diameter tracks the BT-474 curve (0.20/0.29/0.61 at
d5/6/7) or MCF-7 ">500 µm → necrosis".

## 7. Citations
- Thomlinson RH, Gray LH (1955) "The histological structure of some human lung cancers and
  the possible implications for radiotherapy." Br J Cancer 9(4):539-549.
  doi:10.1038/bjc.1955.55. — diffusion-limited hypoxia; cord radius ~200 µm necrosis onset,
  viable rim ≤ ~180 µm, O2 diffusion ~150 µm.
- Greenspan HP (1972) "Models for the growth of a solid tumor by diffusion." Stud Appl Math
  51(4):317-340. — three-radius moving-boundary RD model, constant viable rim.
- Browning AP, Sharp JA, Murphy RJ, et al. (2021) "Quantitative analysis of tumour spheroid
  structure." eLife 10:e73020. PMID 34842141. doi:10.7554/eLife.73020 (via PubMed). —
  limiting structure / arrested region; spheroids reach a limiting size set by avascular
  diffusion regardless of seeding number.
- Grimes DR, et al. (2014) "A method for estimating the oxygen consumption rate in
  multicellular tumour spheroids." J R Soc Interface 11:20131124. doi:10.1098/rsif.2013.1124.
  — r_l/r_n geometry, D_O2, consumption (via ffn KB paper 27).
- ffn KB paper 02 = "Integrating simulation and experimental validation of nutrient-limited
  growth in breast cancer spheroids" (BT-474 COMSOL RD). references/analysis/02_*.md
- ffn KB paper 04 = "A computational model for ... tumour spheroid" (hybrid ABM+RD; viable
  rim ≈ 100 µm constant). references/analysis/04_10-3934-mbe-2018016.md
- ffn KB paper 06 = FUCCI melanoma continuum (three normalized zones).
  references/analysis/06_*.md
- ffn KB paper 16 = TME spheroid review (zonation > ~500 µm; plateau 600–800 µm; rim 20 µm).
  references/analysis/16_*.md
- ffn KB paper 21 = 3D spheroid/organoid TME review (hypoxia < 5–10 mmHg).
  references/analysis/21_*.md
- ffn KB paper 27 = cellular-automaton + RD (HCT116; r_l 233 µm, r_n 155 µm, D_O2 3.8e-9).
  references/analysis/27_*.md
- ffn KB paper 15 = MCF-7/MDA-MB-231 spheroid biology (gaps 5–10 µm, cell ~15 µm).
  references/analysis/15_*.md
- Web: eLife 73020; PMC9876349 (PLOS Comput Biol, time-dependent O2, "100–200 µm rim");
  oncotarget 13857/PMC5352092 (MCF-7 necrosis >500 µm, hypoxia >200 µm).

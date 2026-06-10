# Q3 — Growth-induced SOLID STRESS / pressure in spheroids: compression → quiescence/apoptosis, and media-contact loss → death

**Research question (necrosis/quiescence series, Q3 — MECHANICAL factors).** What is the
magnitude (kPa) of growth-induced solid stress / pressure that accumulates inside spheroids,
what pressure thresholds suppress proliferation vs. induce death, what is the
osmotic/mechanical mechanism, and how does rim proliferation BUILD the compressive load on the
core? Cross-cut with the nutrient/waste (media-contact) death pathway.

Sources: PubMed (cited inline with DOIs) + WebSearch. Primary lineage:
Helmlinger/Jain (Boston/MGH) and Montel-Delarue-Cappello-Joanny (Institut Curie/Grenoble),
plus Stylianopoulos (Cyprus). All numbers carry units.

---

## 1. The two foundational magnitudes (memorize these)

| Quantity | Value | Cell type / system | Source |
|---|---|---|---|
| Solid stress that **inhibits MCTS growth** | **45–120 mmHg = ~6–16 kPa** | LS174T, AsPC-1, etc. in agarose; species/origin-independent | Helmlinger 1997 |
| Solid stress that **drastically reduces growth** (stress-clamp) | **~10 kPa** (5 kPa already measurable) | CT26 / colon-carcinoma spheroids | Montel 2011 |
| Compression that **arrests proliferation via volume loss** | **5–10 kPa** (effect onset ≥ ~0.5–1 kPa) | 5 cell lines incl. HT29, BC52, FHI; mouse + human | Delarue 2014 |
| Growth-induced **residual solid stress in real tumors** | **3.31–10.88 kPa** (avg, breast 4T1 / panc PAN02 / fibrosarcoma MCA205) | murine xenografts | Hadjigeorgiou & Stylianopoulos 2023 |
| Solid stress in **excised murine + human tumors** | **up to ~1.3 kPa (small) → ≳ several kPa, grows with size** | breast, panc, sarcoma, brain mets | Nia et al. 2017 (Nat Biomed Eng) |

Conversion anchor: **1 mmHg = 0.1333 kPa**; so Helmlinger's "45–120 mmHg" = **6.0–16.0 kPa**.
Note Helmlinger's stress exceeds tumor vascular blood pressure → vessel collapse (the
mechanopathology link Jain's group later quantified).

**Convergent answer:** the proliferation-suppressing / growth-arresting solid stress in
spheroids and real tumors lives in the **~1 kPa (onset) → 5–16 kPa (strong inhibition)** band,
i.e. **~10 kPa is the canonical "growth-stops" pressure**.

---

## 2. The five-author backbone (exactly the chain the question names)

### Helmlinger, Netti, Lichtenbeld, Melder, Jain — *Nat Biotechnol* **15**, 778–783 (1997)
"Solid stress inhibits the growth of multicellular tumor spheroids."
DOI: https://doi.org/10.1038/nbt0897-778 (PMID 9255794)
- **Inhibiting stress = 45–120 mmHg (6–16 kPa)**, independent of host species, tissue of
  origin, differentiation state — growth-inhibition is generic to mechanical load.
- This stress **exceeds tumor vascular pressure** → can collapse vessels/lymphatics →
  explains impaired perfusion, poor drug delivery (the founding "solid stress" paper).
- At the cellular level: growth inhibition of plateau-phase spheroids came with **decreased
  apoptosis** and **increased cellular packing density** (cell-volume / cell-shape
  transduction); proliferation not significantly changed in *that* plateau regime.
  (Mechanism nuance: in Helmlinger's plateau spheroids the dominant lever was apoptosis
  suppression + densification; later work below shows proliferation arrest is the dominant
  lever during active growth.)

### Cheng, Tse, Jain, Munn — *PLoS ONE* **4**(2): e4632 (2009)
"Micro-environmental mechanical stress controls tumor spheroid size and morphology by
suppressing proliferation and inducing apoptosis in cancer cells."
DOI: https://doi.org/10.1371/journal.pone.0004632 (PMID 19247489)
- Single cells + fluorescent microbeads co-embedded in agarose; bead-density change → strain
  → **spatial map of peri-spheroid compressive stress**.
- **High-stress regions = suppressed proliferation + induced apoptosis** → controls spheroid
  size & shape (anisotropic stress → non-spherical growth).
- Apoptosis is via the **mitochondrial (intrinsic) pathway** (proven by anti-apoptotic gene
  overexpression rescuing compressed spheroids).
- Frames compression-induced apoptosis as a mechanism of **tumor dormancy** (growth held by
  apoptosis/proliferation balance). This is the cleanest "compression → death" causal paper.

### Montel, Delarue, Elgeti, ... Prost, Cappello, Joanny — *Phys Rev Lett* **107**, 188102 (2011)
"Stress clamp experiments on multicellular tumor spheroids."
DOI: https://doi.org/10.1103/PhysRevLett.107.188102 (PMID 22107677)
- Osmotic stress clamp (Dextran in medium, excluded from spheroid → compresses outer shell,
  transmitted inward).
- **~10 kPa drastically reduces growth, mainly by inhibiting proliferation in the CORE** of
  the spheroid (not the rim). Direct rim-vs-core spatial result.
- Accompanied by a simple continuum/numerical model linking mechanics to growth.

### Delarue, Montel, Vignjevic, Prost, Joanny, Cappello — *Biophys J* **107**, 1821–1828 (2014)
"Compressive stress inhibits proliferation in tumor spheroids through a volume limitation."
DOI: https://doi.org/10.1016/j.bpj.2014.08.031 (PMID 25418163)
- **The mechanistic time-cascade (the most model-actionable paper):**
  - **minutes:** compression → **MCS volume drop**, traced to **cell-volume reduction in the
    core**.
  - **hours:** reversible induction of CDK inhibitor **p27^Kip1** spreading center→periphery.
  - **days:** cells **blocked at the late-G1 restriction point** (G1/S checkpoint).
- p27^Kip1 silencing **abolishes** the pressure effect → p27 is the molecular relay.
- **Quantitative correlation:** pressure-induced **volume change ↔ spheroid growth-rate
  reduction** (cell volume is the read-out / control variable).
- Effect conserved across **5 cell lines**, fully **reversible**. Mechano-osmotic stress of
  **~5 kPa** is physiological for growth-constrained environments; ~5% MCS volume loss in
  ~30 min at 5 kPa.

### Hadjigeorgiou & Stylianopoulos — *Biomech Model Mechanobiol* **22**, 1625–1643 (2023)
"Evaluation of growth-induced, mechanical stress in solid tumors and spatial association with
extracellular matrix content."
DOI: https://doi.org/10.1007/s10237-023-01716-3 (PMID 37129689)
- Resect tumor → embed in agarose → perpendicular cuts release residual stress → bulging
  displacement + FEM → stress map.
- **Average residual stress 3.31–10.88 kPa**, **non-uniform** (microenvironment heterogeneity).
- Stress spatially **co-localizes with hyaluronan + collagen** content (ECM is a stress
  reservoir). Elastic modulus correlates nonlinearly with tumor volume in fibrotic tumors.

### (supporting Stylianopoulos/Jain) Stylianopoulos et al. — *PNAS* **109**(38), 15101–15108 (2012)
"Causes, consequences, and remedies for growth-induced solid stress in murine and human tumors."
DOI: https://doi.org/10.1073/pnas.1213353109 (PMID 22932871)
- Planar-cut / opening-angle technique to estimate growth-induced solid stress in vivo.
- **Solid stress compresses intratumoral blood + lymphatic vessels** → ↑interstitial fluid
  pressure, ↓blood flow → **hypoxia** → progression, immunosuppression, worse therapy.
- Stress **reducible by depleting cancer cells, fibroblasts, collagen, and/or hyaluronan**
  (CAF depletion via SHH inhibitor decompresses vessels, improves perfusion). Establishes the
  cancer-cell + ECM + confinement decomposition of solid stress.

### Nia, ... Munn, Jain — *Nat Biomed Eng* **1**, 0004 (2017)
"Solid stress and elastic energy as measures of tumour mechanopathology."
DOI: https://doi.org/10.1038/s41551-016-0004 (PMID 28966873)
- Three quantitative techniques → **2D maps of solid stress + elastic energy**, sensitive even
  in small tumors and in vivo.
- Three structural conclusions: (i) solid stress depends on **both cancer cells AND
  microenvironment**; (ii) **solid stress increases with tumor size**; (iii) **mechanical
  confinement by surrounding tissue** is a major contributor.

---

## 3. HOW rim proliferation BUILDS core compression (the spatial mechanism)

This is the load-path the question asks for — the causal chain that makes the **core** the
high-stress, low-proliferation, eventually-dying compartment:

1. **Proliferation is rim-localized** (outer shell has nutrient/O2 access; see §5). New cell
   volume is generated at the periphery.
2. **Confinement** — by surrounding matrix (agarose / ECM in vitro; host tissue in vivo) OR by
   the spheroid's own outer shell — means the new volume cannot freely expand outward.
   (Nia 2017: confinement by surrounding tissue is a major stress source; Helmlinger 1997:
   agarose stiffness sets the inhibiting stress.)
3. **Stress is transmitted inward and AMPLIFIES toward the core.** Dolega et al. (cell-like
   PAA pressure sensors) — *Nat Commun* **8**, 14056 (2017),
   DOI: https://doi.org/10.1038/ncomms14056 (PMID 28128198): under isotropic compression the
   internal pressure is **non-uniform and RISES toward the core**, explained by the
   **anisotropic radial cell arrangement** alone (single cell type sufficient). This directly
   links **elevated core mechanical stress ↔ the known lack of core proliferation**.
4. **Core cells respond:** volume reduction (Delarue 2014, minutes) → p27^Kip1 induction
   center→periphery (hours) → late-G1 / restriction-point arrest = **quiescence** (days), and
   at higher/sustained stress → **mitochondrial apoptosis** (Cheng 2009). Montel 2011 confirms
   the proliferation loss is **core-localized** under clamp.

So: **rim growth + confinement → inward-amplifying compressive stress → core volume loss →
p27/G1 arrest (quiescence) → intrinsic apoptosis (death).** Mechanics co-acts with the nutrient
gradient (§5) — both peak adverse in the core.

---

## 4. The osmotic / volume mechanism and its threshold (how compression actually kills/arrests)

The transduction is **volume/crowding-based**, not a simple "filaments bear the load":

- **Cell-volume reduction is the primary signal.** Delarue 2014: compression → core cell-volume
  drop in minutes; the **volume change quantitatively predicts the growth-rate drop**. Volume↓ →
  **macromolecular crowding↑** → impaired biosynthesis/cell-cycle progression → p27↑ → G1 arrest.
- **Cytoskeleton does NOT actively bear the compressive load — osmoregulation does.**
  McGrail, ... Dawson — *Biophys J* **109**(7), 1334–1337 (2015), "Osmotic Regulation Is Required
  for Cancer Cell Survival under Solid Stress," DOI: https://doi.org/10.1016/j.bpj.2015.07.046
  (PMID 26445434): in breast/ovarian/prostate spheroids, under compression cells **actively
  efflux Na+** (via **NHE1**) to lower intracellular tonicity → shed water → shrink to relieve
  load. **Blocking NHE1 or depolymerizing actin → MORE compression-induced death.** Actin is not
  a load-bearer but is **required for the Na+ efflux**; raising actin polymerization is
  protective. So **death under solid stress = failure of osmotic volume regulation**.
- **Osmotic ≈ mechanical equivalence for the spheroid response.** Dextran osmotic stress and
  direct mechanical compression produce the same growth/volume response (Montel 2011, Delarue
  2014) — justifies modeling solid stress via an osmotic-pressure–like term on the cell.

**Threshold summary (operational):**
- **Onset / quiescence:** detectable proliferation slowing from **~0.5–1 kPa**, robust
  proliferation arrest by **~5 kPa** (Delarue 2014; Montel 2011 measurable at 5 kPa).
- **Strong growth arrest / death:** **~10 kPa** (Montel 2011; Cheng 2009 high-stress regions),
  up to **6–16 kPa** for full inhibition (Helmlinger 1997). Real-tumor residual stress
  **3.3–10.9 kPa** (Stylianopoulos 2023) sits exactly in this band → the in-vitro thresholds
  are physiological.

---

## 5. The OTHER core-death driver: loss of media (nutrient/waste) contact

Mechanics is only half the core-death story; the diffusion limit is the classic driver and the
two compound (rim proliferation simultaneously builds stress AND consumes the nutrient that the
core needs):

- **Diffusion-limited viable rim ≈ 100–200 µm.** Cells beyond ~150–200 µm from the
  nutrient/O2 source (spheroid edge in vitro; nearest capillary in vivo, ~100–200 µm) drop
  below survival O2/glucose and accumulate waste → **hypoxic/necrotic core**. (Boot/Koenderink
  2021 review, *Adv Phys X* 6:1978316, DOI: https://doi.org/10.1080/23746149.2021.1978316 —
  KB paper #11, §2.1: "inner cells beyond ~200 µm from the edge undergo apoptosis → necrotic
  core.")
- **Three-shell radial structure:** proliferating rim → quiescent (nutrient-poor but viable)
  intermediate → necrotic core — the standard MCTS architecture (KB #11, #19).
- **Necrosis onset ~ spheroid diameter ≳ 400–500 µm** (central anoxia); O2 consumption per
  volume falls as radius grows (supply-limited core). (Grimes et al. 2014, *J R Soc Interface*
  11:20131124, DOI: https://doi.org/10.1098/rsif.2013.1124; and O2-tension spheroid studies.)
- **Mechanics ⇄ nutrient coupling (in vivo):** solid stress collapses vessels → ↓perfusion →
  hypoxia (Stylianopoulos 2012 PNAS) — so compression *also* starves the core. The two death
  drivers are not independent; they co-localize in the core.

---

## 6. MCF7-specific notes

- **No clean MCF7-only kPa threshold** appears in the canonical solid-stress papers (Helmlinger
  used LS174T/AsPC-1; Montel/Delarue used CT26/HT29/colon lines; Cheng used murine mammary
  67NR/4T1-type lines). The thresholds are explicitly **cell-line-independent** (Helmlinger
  1997: independent of tissue of origin; Delarue 2014: conserved across 5 lines), so the
  **6–16 kPa inhibition / ~5–10 kPa arrest band transfers to MCF7**.
- McGrail & Dawson 2015 (osmotic regulation, NHE1) used **breast** (among ovarian/prostate)
  cancer spheroids → the Na+/NHE1 osmoregulation-under-compression mechanism is validated in
  breast lineage.
- WebSearch surfaced a **2025 bioRxiv preprint** ("Stress-dependent growth of breast cancer
  models arises from a cellular volume checkpoint," biorxiv 2025.07.29.667388) reinforcing the
  Delarue **cell-volume-checkpoint** mechanism specifically for **breast cancer models** —
  not yet peer-reviewed; flag as preprint, do not anchor a gate on it.
- For MCF7 baseline physiology already in the project: cytoplasm viscosity ~65 Pa·s, resting
  osmotic/turgor ~40 Pa (CLAUDE.md physiological-operating-point rule) — the solid-stress
  thresholds here (kPa) are ~100× the resting turgor, consistent with "stress must build to kPa
  to matter."

---

## 7. How to turn this into a DCM necrosis/quiescence criterion

(DCM = the spheroid/multicell composition layer; this is mechanism-level guidance, not a fit.)

1. **State variable per core cell:** local compressive solid stress σ (Pa), computed from the
   particle/contact pressure (HOOMD virial / contact-normal stress) OR an osmotic-equivalent
   term. Track **cell volume V** as the actual transduction read-out (Delarue 2014: V is the
   control variable, not σ directly).

2. **Quiescence (G1 arrest) criterion** — two equivalent gates:
   - stress gate: **σ ≳ ~1 kPa onset, full arrest at σ ≳ ~5 kPa** → suppress division
     (set proliferation rate → 0). Use a smooth p27-like sigmoid between ~1 and ~5 kPa rather
     than a hard step (matches center→periphery p27 spread).
   - volume gate (preferred, more mechanistic): arrest when **ΔV/V₀ exceeds a threshold**
     (~5% volume loss already at 5 kPa, Delarue 2014) — i.e. a **cellular volume checkpoint**.
   Quiescence must be **reversible** (decompress → resume) per all Curie papers.

3. **Death (necrosis/apoptosis) criterion** — compound, OR of two channels:
   - **mechanical/intrinsic apoptosis:** sustained **σ ≳ ~10 kPa** (Montel/Cheng) over a
     residence time, gated by osmoregulation capacity (if Na+/NHE1 efflux saturates →
     volume can't be defended → death). Cheng 2009: mitochondrial pathway.
   - **media-contact / nutrient channel:** cell depth from spheroid surface **> ~150–200 µm**
     AND local O2/nutrient below survival threshold → necrosis (the diffusion-limit rule, §5).
   In vivo the two couple (stress → vessel collapse → hypoxia); in an avascular spheroid sim
   keep them as independent OR-ed channels but co-located in the core.

4. **Stress-build mechanism (don't impose σ — let it emerge):** keep proliferation **rim-biased**
   (gate division on nutrient access), enforce **confinement** (outer-shell/ECM boundary at
   physiological stiffness, ~agarose/host modulus), and let contact mechanics propagate the
   load. Dolega 2017 says the **core stress amplification emerges from radial cell anisotropy
   alone** — a fine-grained contact model should reproduce the inward pressure rise without a
   hand-coded radial σ(r). Validate against Dolega's measured inward-rising pressure profile.

5. **Osmotic implementation hook:** model solid stress on a cell as an osmotic-pressure term
   (Montel/Delarue osmotic≈mechanical equivalence), with an active Na+/water efflux relief
   capacity (McGrail 2015); death when relief capacity is exceeded. Aligns with the
   project's enclosed-volume Π term (CLAUDE.md: physiological baseline turgor).

**Acceptance oracles (numbers to hit):** proliferation halves around σ ~5 kPa; near-zero growth
by σ ~10 kPa; full inhibition 6–16 kPa; ~5% volume loss at 5 kPa in ~30 min; inward-rising core
pressure profile (Dolega); viable rim ~100–200 µm; necrotic core at spheroid Ø ≳ 400–500 µm.
All reversible for the quiescence channel.

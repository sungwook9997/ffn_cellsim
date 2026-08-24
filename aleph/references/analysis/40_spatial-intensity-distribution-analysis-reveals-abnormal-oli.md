---
id: 40_spatial-intensity-distribution-analysis-reveals-abnormal-oli
paper_n: 40
title: "Spatial Intensity Distribution Analysis Reveals Abnormal Oligomerization of Proteins in Single Cells"
authors: "Godin, A. G., Rappaz, B., et al."
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.06.068"
paper_type: methods-software
ffn_relevance: Low
ffn_themes: [tangential, numerics/methods]
entities: [plp, ampa-receptor, glutamate-receptor, glur1, glur2, gfp, mgfp, cos-7, hek293, endoplasmic-reticulum, plasma-membrane]
methods: [confocal-live-imaging, fluorescence-fluctuation-spectroscopy, monte-carlo, image-histogram-fitting]
measurables: [oligomerization-state, fluorescent-particle-density, quantal-brightness, labeling-probability]
keywords: [spida, oligomerization, receptor-clustering, fluorescence-microscopy, mislabeling-correction, photobleaching, super-poissonian, intensity-histogram, proteolipid-protein, single-cell-imaging]
tags: ["#methods-software", "#fluorescence-imaging", "#oligomerization", "#single-cell-imaging", "#tangential", "#not-cell-mechanics"]
has_transferable_params: false
---

# [40] Spatial Intensity Distribution Analysis Reveals Abnormal Oligomerization of Proteins in Single Cells

**Tags:** #methods-software #fluorescence-imaging #oligomerization #single-cell-imaging #tangential #not-cell-mechanics

| Field | Value |
|---|---|
| Authors | Godin, A. G., Rappaz, B., et al. (Wiseman lab, McGill) |
| Year / Venue | 2015 / Biophysical Journal 109(4):710–721 |
| DOI / ID | 10.1016/j.bpj.2015.06.068 |
| Type | methods-software (fluorescence image-analysis technique) |
| Pages | 12 |
| ffn_cellsim relevance | Low — a fluorescence-microscopy oligomerization-counting method; no mechanics, cytoskeleton, ECM, motors or adhesions |

## 1. Summary
The paper extends Spatial Intensity Distribution Analysis (SpIDA), a fluorescence-image histogram-fitting technique, to measure the oligomerization state (subunit count) and density of membrane/receptor proteins in subcellular compartments of single intact cells. The core advance is a probability-weighted correction algorithm for **nonemitting/mislabeled fluorophores** (e.g., the ~20% of GFP molecules that misfold and never emit, or photobleaching over a time series): SpIDA fits a binomial model for the fraction p of emitting subunits, recovering the true underlying oligomer distribution instead of a biased one. They validate with Monte-Carlo-simulated confocal images, calibrate p in cells using receptors of known stoichiometry (tetrameric AMPA receptor = 2×GluR1 + 2×GluR2; p≈67–84%), and then apply the method to proteolipid protein (PLP) in COS-7 cells. They find that a trafficking-impaired PLP mutant (D202N) accumulates as tetramers retained in the endoplasmic reticulum, while wild-type and mutant are indistinguishable at the plasma membrane — quantitatively linking the trafficking defect (Pelizaeus-Merzbacher disease) to abnormal ER oligomerization.

## 2. Problem & motivation
Receptor organization and clustering state govern the first steps of cell signaling and trafficking, but quantifying protein–protein interactions and oligomerization *in situ*, per-compartment, in live single cells is hard: bulk biochemistry (gel/Western blot) averages away local subcellular states and loses dynamics, and all fluorescence-based oligomer-counting methods are biased by imperfect labeling (misfolded/dark fluorophores, mislabeling, photobleaching). The motivation is a nondestructive, image-based, single-cell oligomer-counting method that corrects for these labeling artifacts.

## 3. Methods / model
- **Technique class:** SpIDA — fits pixel-intensity histograms from a region of interest (ROI) of a single confocal image with **super-Poissonian distributions** (Poisson-distributed particles convolved with the optical PSF). Fit parameters: fluorescent particle densities N (per beam-focus area, BA) and quantal brightness ε (intensity per oligomer). Inspired by the temporal photon-counting histogram method (Chen 1999) but applied in the spatial domain on conventional analog-PMT confocals.
- **Governing equations:** one-population fit H(ε,N,k)=Σ rₙ(ε;k)·poi(n,N) (Eq. 1); multi-population mixtures by convolution (Eqs. 2–4, two/three populations); **mislabeling correction (Eq. 5)** convolves binomial-weighted contributions of partially labeled oligomers, with subunit-label counts Nₗ = (nmer!/(l!(nmer−l)!))·pˡ(1−p)^(nmer−l)·N for l emitting subunits and quantal brightness l·ε₀.
- **Numerics:** custom MATLAB routines (Image Processing + Optimization toolboxes); Monte-Carlo image generation (random particle placement, Gaussian-PSF convolution, simulated detector shot noise with variance ≈ 10·iu·⟨i⟩); fit times <1–10 s.
- **Cells/assays:** HEK293 cells transfected with GluR1/GluR2 AMPA-receptor subunit combos (known tetramer, used to calibrate p); COS-7 cells expressing PLP-mGFP wild-type vs. D202N mutant; mGFP-farnesyl as monomeric brightness standard. Confocal: Olympus FV300-IX71, 60× NA-1.4 oil objective, 488 nm Ar-ion laser. ER vs. membrane compartments segmented by Otsu threshold + dilation/erosion contiguity.

## 4. Key results (quantitative)
- For reliable single-/two-population SpIDA: need S/N ≈ 3:1 and ROI sampling ≈ 50 beam-focus areas (~6 µm² at NA 1.4) → <20% error on all fit parameters (Results, p.713).
- Three-population (monomer/dimer/tetramer) fits converge with <20% error given sufficient S/N and sampling (Figs. 1, 5, 6).
- With mislabeling correction, recovering set density + labeling probability to ~20% precision needs >100 BAs (~12 µm²) (Fig. 4, p.717).
- AMPA-receptor calibration: single GFP-tagged subunit gave brightness 16% below the dimer expectation → p = 67 ± 7%; both subunits tagged gave brightness 12% below tetramer → p = 84 ± 5% (Fig. 7B). Consistent with the literature ~15–25% nonfluorescent-GFP fraction (refs 14,15).
- Monomeric mGFP-farnesyl quantal brightness standard: 3.9 ± 0.1 Miu/s (10⁶ intensity units/s) (p.718).
- PLP photobleaching time series: intensity dropped ~30% over 5 images; with iterative p adjustment p(n)=80%·(⟨i⟩ₙ/⟨i⟩₁), recovered oligomer densities stayed constant within error (Fig. 8D); example membrane fit: 106 monomers/µm², 13.3 dimers/µm², 5.1 tetramers/µm².
- Mutant D202N PLP: significant ER retention with a dominant tetramer population (p_t-test < 0.001), reduced membrane density (p=0.035) and higher ER density (p=0.047) vs WT; no membrane oligomerization difference (Fig. 9). Example mutant fits: membrane N=(24.7,0.3,0.1) µm⁻²; ER N=(42.6,20.9,24.8) µm⁻² (Fig. 8E).

## 5. Parameters & constants of interest
None (no transferable constants). The only quantitative anchors are imaging/optical and label-photophysics numbers — nonfluorescent-GFP fraction ~15–25% (p≈0.67–0.84), monomeric mGFP brightness 3.9 Miu/s, receptor surface densities in µm⁻². None of these inform a mechanistic single-cell cytoskeleton+ECM particle simulator. There is no mechanical, kinetic-binding (off-rate/force), elastic, or motor constant in the paper.

| Quantity | Value+units | Source-in-paper |
|---|---|---|
| Nonfluorescent GFP fraction | ~15–25% (p ≈ 0.67–0.84) | Fig. 7B; refs 14,15 |
| Monomeric mGFP quantal brightness | 3.9 ± 0.1 Miu/s | p.718 |
| Min ROI for <20% error | ~50 BA (~6 µm²) | Results p.713 |
| PLP membrane densities (WT) | 106 / 13.3 / 5.1 µm⁻² (mono/di/tetra) | Fig. 8 caption |

(All are fluorescence-counting metrology values, not mechanics — listed only for completeness; has_transferable_params = false.)

## 6. Relevance to ffn_cellsim  ← MOST IMPORTANT
**Off-topic / tangential (category f).** This is a fluorescence-microscopy image-analysis method (SpIDA) for counting protein oligomerization states and densities of membrane receptors in single cells. It contains:
- **No mechanics** — no forces, stresses, moduli, traction, viscoelasticity.
- **No cytoskeleton, ECM, motors, or adhesions** — nothing on actin, myosin, clutches, integrins, cross-links, collagen. The "proteins" studied are signaling receptors (AMPA/glutamate, PLP), not the structural/force-generating components ffn_cellsim models.
- **No transferable mechanism or oracle** — the binomial mislabeling correction and super-Poissonian histogram fit are an imaging metrology tool, not a physical model of cell behavior; they cannot serve as a validation oracle (ffn_cellsim's oracles are closed-form mechanics: Bell-Evans, Hill, Stam-Hocky, KU catch-bonds).

Mapping to H-units: none. It does not touch H.1 cortex, H.2, H.3, H.5, H.7, nor the EXTEND membrane/nucleus/cytoplasm track in any mechanistic sense, nor the spheroid/multicellular validation overlay (MCF7-on-pV4D4/collagen-I traction/cohesion). The lone tenuous connection is the word "membrane/ER compartments," but ffn_cellsim's membrane EXTEND (H.8) is about mechanical bilayer/cortex coupling, not receptor oligomer counting — so even there it is not useful. The single conceivable indirect use is *experimental-method awareness*: if the project ever wants an imaging readout to validate a simulated quantity in PI's overlay experiments, SpIDA is one technique for counting molecule densities/clusters — but ffn_cellsim's validation targets are traction and cohesion, which this method does not measure. Classify **Low / tangential**.

## 7. Limitations & caveats
- It is a 2D image-analysis method on diffraction-limited confocal data; axial resolution limits mean ER-vs-membrane assignment is uncertain (authors note ER "monomers/dimers" may be underlying membrane).
- Requires a priori knowledge of possible oligomer states and a calibrated monomeric brightness standard; over-parameterized fits (≥4 free densities) do not converge uniquely.
- Assumes equal, independent subunit-labeling probability p (binomial), negligible FRET/quenching, and linear analog-detector regime; needs per-setting detector-noise calibration.
- For a fine-grained particle simulator the scale/domain gap is total: this is fluorescence photophysics + statistics, not mechanics.

## 8. Key figures / tables
- **Fig. 5 / Fig. 6** — three-population (monomer/dimer/tetramer) SpIDA accuracy vs. ROI size and with 80% mislabeling; shows convergence and precision-vs-density-contribution.
- **Fig. 7** — AMPA-receptor calibration of labeling probability p (67% dimer, 84% tetramer) from quantal-brightness deviations.
- **Fig. 8** — PLP-mGFP WT vs D202N example images, photobleaching time-series compensation, and ER/membrane mask histograms with fitted densities.
- **Fig. 9** — Quantitative ER tetramer retention of mutant PLP vs membrane (t-test significances).

## 9. Notable quotes / citable claims
- "it has been reported that ~20% of fluorescent proteins misfold and do not emit any fluorescence" (Introduction, p.711) — motivates the mislabeling correction.
- "a large significant population of tetramers accumulated in the ER for the mutant protein, while the normal oligomerization state was measured for the smaller fraction of the receptor population that reached the membrane" (Abstract/Introduction).
- "assuming reasonable signal/noise (~3:1) ... and ~50 beam-focus areas (~6 µm²) ... SpIDA can give accurate results (<20% error on all fit parameters)" (Results, p.713).
- "the PLP proteins retained in the ER are virtually all present as tetramers" (p.719).

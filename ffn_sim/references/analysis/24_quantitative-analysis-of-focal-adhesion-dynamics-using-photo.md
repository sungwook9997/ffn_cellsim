---
id: 24_quantitative-analysis-of-focal-adhesion-dynamics-using-photo
paper_n: 24
title: "Quantitative analysis of focal adhesion dynamics using photonic resonator outcoupler microscopy (PROM)"
authors: "Yue Zhuo et al."
year: "2018"
venue: "Light: Science & Applications (Nature)"
doi: "10.1038/s41377-018-0001-5"
paper_type: methods-software
ffn_relevance: Low
ffn_themes: [FA/clutch, tangential]
entities: [focal-adhesion, integrin, vinculin, actin, fibronectin, nucleus, lipid-bilayer, photonic-crystal, filopodia, mhat9a]
methods: [finite-difference-time-domain, photonic-crystal-enhanced-microscopy, label-free-optical-imaging, fluorescence-microscopy, confocal-imaging, scanning-electron-microscopy]
measurables: [focal-adhesion-area, refractive-index, peak-intensity-shift, peak-wavelength-shift, spatial-resolution, evanescent-field-depth]
keywords: [focal-adhesion, label-free-imaging, photonic-crystal-biosensor, peak-intensity-shift, light-scattering, refractive-index-contrast, evanescent-field, vinculin, stem-cell-adhesion, FA-dimension]
tags: ["#focal-adhesion", "#label-free-imaging", "#optical-biosensor", "#FA-imaging", "#methods-software", "#tangential"]
has_transferable_params: false
---

# [24] Quantitative analysis of focal adhesion dynamics using photonic resonator outcoupler microscopy (PROM)

**Tags:** #focal-adhesion #label-free-imaging #optical-biosensor #FA-imaging #methods-software #tangential

| Field | Value |
|---|---|
| Authors | Yue Zhuo, Ji Sun Choi, Thibault Marin, Hojeong Yu, Brendan A. Harley, Brian T. Cunningham |
| Year / Venue | 2018 / Light: Science & Applications (Nature) |
| DOI / ID | 10.1038/s41377-018-0001-5 |
| Type | methods-software (optical imaging instrumentation + FDTD electromagnetic modeling) |
| Pages | 15 |
| ffn_cellsim relevance | Low — a label-free FA *imaging* modality (optics), not a mechanism/force model; only loosely touches FA size as an observable |

## 1. Summary
The authors introduce photonic resonator outcoupler microscopy (PROM), a label-free optical imaging technique that detects focal-adhesion (FA) clusters at the cell–substrate interface. A photonic-crystal (PC) biosensor supports a resonant evanescent standing wave (~200 nm penetration depth); when a localized high-refractive-index protein cluster (a FA) forms at the surface, it outcouples resonant photons by scattering, producing a highly localized drop in reflected peak intensity (the "peak intensity shift", PIS). PROM simultaneously extracts two orthogonal contrasts from the same reflected spectrum at each ~0.6×0.6 µm² pixel: PIS (intensity, sensitive to localized scattering protein clusters/FAs) and PWS (peak-wavelength shift, sensitive to bulk adsorbed mass density). Using murine dental epithelial stem cells (mHAT9a) on fibronectin, they show PIS maps form a "ring" of high intensity along the cell periphery that co-localizes with vinculin fluorescence (a FA marker) — but not with bulk actin or nucleus signals — and track FA evolution over ~26 min with 10-s temporal resolution without photobleaching. FDTD electromagnetic simulations support the scattering-outcoupling hypothesis.

## 2. Problem & motivation
FA cluster size correlates with adhesion engagement and migration speed, and FA dynamics underlie metastasis, apoptosis, chemotaxis, and drug response. Standard FA imaging relies on fluorescence tagging (vinculin, paxillin, etc.), which suffers from photobleaching (limiting long-term/repeated imaging), cytotoxicity, and quantitation difficulty. The paper's goal is a label-free, quantitative, long-term imaging modality to map FA areas/dimensions and their kinetics in live cells. This is an *instrumentation/measurement* problem, not a force-generation or cytoskeletal-mechanics problem.

## 3. Methods / model
- **Imaging modality**: PROM = modified brightfield microscope with a line-scanning imaging spectrometer reading the resonant reflected spectrum per pixel off a 1D photonic-crystal biosensor (illumination from below, polarized perpendicular to grating). Two contrasts: PIS (peak intensity reduction → localized scattering FA clusters) and PWS (peak wavelength red-shift → bulk adsorbed mass density). Axial resolution set by ~200 nm evanescent field; lateral by photon propagation.
- **PC biosensor**: 1D UV-curable-polymer grating (n0 = 1.46, grating depth dg = 120 nm, period Λ = 400 nm, duty cycle fg = 41.6%, sidewall angle 85°) coated with TiO₂ thin film (n1 = 2.4, slab thickness ds = 61 nm, duty cycle 50%, sidewall angle 82°); resonant reflection near λ0 ≈ 626 nm. Fabricated by room-temperature replica molding + reactive sputter deposition.
- **Electromagnetic model**: FDTD (Lumerical) of the PC with/without a FA. FA modeled as a homogeneous, lossless dielectric sphere (n_FA ≈ 1.46, radius 50–500 nm) embedded in surrounding medium (n2 ≈ 1.333). Parametric sweeps over background refractive index (1.333–1.373) and FA/nanoparticle radius (50–500 nm).
- **Cells / assay**: murine dental epithelial stem cells (mHAT9a) in DMEM + 10% FBS, 37 °C / 5% CO₂, on fibronectin-coated PC surface. Validation by fluorescence (nucleus, actin, vinculin), confocal, phase-contrast, and SEM of the same cells.
- **Readout**: PIS/PWS images over ~300×300 µm² fields, ~10-s scan interval, time courses to ~26 min; cross-section line profiles, edge-vs-inner statistics (N = 5 cells), spatiotemporal band maps.

## 4. Key results (quantitative)
- PC resonance near λ0 ≈ 626 nm; evanescent field extends ~200 nm into aqueous medium (p.4–5, Fig. 2c,d).
- On cell attachment, peak wavelength red-shifts from λ_BG ≈ 626 nm to λ_cell ≈ 628 nm; normalized peak intensity drops from ~90% to ~80% (p.6).
- Pixel size ~0.6×0.6 µm²; field of view up to ~300×300 µm²; shortest scan interval ~10 s for ~100×100 µm² (p.6).
- FA cluster cross-sectional area typically 0.2–1.0 µm² (used to argue measurable scattering; p.6). Non-mature focal complexes (FXs) <0.2 µm²; mature FAs typically 1–10 µm² (intro, p.2).
- Refractive indices: TiO₂ n ≈ 2.4; medium/water n2 ≈ 1.333; averaged cell n_cell = 1.35–1.38; FA modeled n_FA ≈ 1.46 (p.3, p.6).
- FDTD: PWS increases with surrounding-medium index but is FA-independent; PIS increases as FA radius grows (50→500 nm), consistent with more scattering-induced outcoupling (Fig. 2e–j).
- Point-spread-function FWHM (for ~100-nm-diameter TiO₂ nanoparticle, from prior work ref.50): ~1.20–1.56 µm for PWS images, ~0.95 µm for PIS images (p.11). High-contrast objects perturb PWV up to ~2–3 µm laterally; low-contrast cell features are more localized.
- PIS "ring" along cell periphery co-localizes with vinculin fluorescence, not with actin/nucleus (Fig. 5); edge > inner for both PIS and vinculin (Fig. 5c, N = 5).

## 5. Parameters & constants of interest

| Quantity | Value + units | Source in paper |
|---|---|---|
| Mature FA cluster area | 1–10 µm² (typical) | Intro, p.2 |
| Non-mature focal complex (FX) area | <0.2 µm² | Intro, p.2 |
| FA cross-sectional area (working range) | 0.2–1.0 µm² | p.6 |
| Averaged cell refractive index | n_cell = 1.35–1.38 | p.3 |
| FA refractive index (model) | n_FA ≈ 1.46 | p.6 |
| Surrounding medium (water) index | n2 ≈ 1.333 | p.3, p.6 |
| FA modeled radius range | 50–500 nm | p.6 |
| Evanescent field penetration depth | ~200 nm | p.5 |

Note: These are *optical/geometric descriptors* (FA area, refractive index) — useful only as soft sanity ranges for FA *size*, not as mechanical/kinetic constants (no stiffness, off-rate, force, or traction values appear). No traction force, bond lifetime, clutch stiffness, or rate constants are reported. **has_transferable_params: false** for mechanistic simulation purposes.

## 6. Relevance to ffn_cellsim  (MOST IMPORTANT)
**Classification: (f) tangential, with a weak (b) parameter-source touch on FA size only.**

ffn_cellsim builds FAs mechanistically as explicit integrin–adhesion *clutch* particles/bonds with Bell-Evans force-dependent off-rates and traction generation (H.4 FA/bridge unit; KU-3.5 / KU-5.1). This paper contributes none of that physics. It is an *optical imaging instrumentation* paper: its core content is photonic-crystal biosensor design, FDTD electromagnetic simulation of light outcoupling, and a label-free readout that reports FA *dimension/area* as a scattering signal. There is:
- **No force, traction, or stress measurement** — explicitly an imaging modality, not traction-force microscopy. So it is *not* a validation oracle for the FA clutch's force/traction outputs (unlike the project's MCF7-spheroid traction-force overlay target).
- **No mechanism** transferable to the runtime: the FDTD "FA = lossless dielectric sphere, n=1.46, r=50–500 nm" model is an *optical proxy*, the opposite of the project's fine-grained particle/bond FA. Adopting it would violate the architectural principle (no lumped proxy mechanisms).
- **No rate constants / kinetics** of FA assembly–disassembly in quantitative form — only qualitative temporal PIS/PWS curves over ~26 min, in arbitrary units.

The only mildly useful items for ffn_cellsim:
1. **FA size ranges** (focal complex <0.2 µm²; mature FA 1–10 µm²; working 0.2–1.0 µm²) — a literature-anchored sanity band for the *emergent* FA cluster footprint that the H.4 clutch model should reproduce. These are well-established numbers (also citable from refs.13, 29 herein: Gallant 2005, Kim & Wirtz 2013 "FA size uniquely predicts cell migration"), not unique to this paper.
2. **Cell refractive index 1.35–1.38** — irrelevant to a mechanical simulator (optical, not a mass/stiffness anchor).
3. **Wrong cell line for the project's overlay** — mHAT9a dental epithelial stem cells on fibronectin, not the MCF7-breast-cancer-on-pV4D4/collagen-I platform the project validates against. So it does not even serve as an experimental overlay target.

**Bottom line:** Off-topic for a mechanistic single-cell cytoskeleton+ECM simulator. Keep it as a low-priority pointer for the *idea* that FA area is observable/dynamic, and for the FA-size sanity band, but it is neither a parameter source for clutch mechanics, a validation oracle for traction, nor a numerics reference. The richer FA-mechanics citations are its references (Kim & Wirtz FASEB 2013, ref.29; Gallant Mol Biol Cell 2005, ref.13; Grashoff vinculin tension Nature 2010, ref.42; Kanchanawong nanoscale architecture Nature 2010, ref.30), not this paper itself.

## 7. Limitations & caveats
- **Optical proxy, not physics**: FA represented as a homogeneous lossless dielectric sphere with a single refractive index — no internal structure, no protein composition, no force, no kinetics. Cannot inform a particle/bond FA model.
- **Arbitrary-unit dynamics**: PIS and PWS time courses are in AU with "different dynamic ranges", explicitly stated as hard to directly compare — so no extractable assembly/disassembly rate constants.
- **Axial-only sensitivity / lateral smearing**: signal restricted to ~200 nm evanescent field at the cell base; lateral resolution limited to ~1 µm (PSF), set by photon propagation, not true FA boundaries.
- **Small N, single cell type**: N = 5 cells, mHAT9a stem cells only; co-localization with vinculin is qualitative ("nearly identical", "probably co-localized").
- **Scale gap**: the project is a particle simulator at nm–µm with explicit bonds; this paper's contribution is at the µm-scale optical-image level and offers no sub-resolution mechanism.

## 8. Key figures / tables
- **Fig. 2 (p.4)**: PC biosensor SEM + FDTD model; parametric reflectance vs background index (1.333–1.373) and vs FA/nanoparticle radius (50–500 nm) — the electromagnetic basis showing PIS grows with FA size while PWS does not.
- **Fig. 3 (p.7–8)**: Time-lapse PROM (BF/PIS/PWS) of an mHAT cell attaching over 0–26 min; spectra at inside/boundary/outside points showing the resonance intensity drop on attachment.
- **Fig. 5 (p.9–11)**: Side-by-side PIS vs fluorescence (nucleus/actin/vinculin); the key validation — PIS "ring" co-localizes with vinculin (FA marker), edge>inner statistics (N=5).
- **Fig. 6 (p.11–12)**: Five cells across PWS/PIS/confocal-FL/phase-contrast/SEM; spatiotemporal band maps and PWS/PIS adhesion time courses.

## 9. Notable quotes / citable claims
- "the size of the FA cluster varies and is highly correlated with the level of adhesion engagement and migration speed" (Intro, p.2) — restates Kim & Wirtz / Gallant, the FA-size↔migration link.
- "non-mature focal complexes (FXs) are initially formed at the leading edge ... usually <0.2 µm². ... some of the FXs grow larger (typically 1–10 µm²) and assemble into mature FA clusters" (Intro, p.2) — the FA-size sanity band.
- "we represent the FA as a locus of a material with a designated radius (50–500 nm) and designated refractive index elevation (n_FA = ~1.46) in contrast to the surrounding medium (n2 = ~1.333)" (p.7) — the optical-proxy FA model.
- "the PIS image shows a nearly identical distribution pattern along the cell peripheral region to that obtained by fluorescence microscopy with a labeled vinculin where the FA areas are concentrated along the cell boundary" (p.8) — the label-free↔vinculin validation claim.

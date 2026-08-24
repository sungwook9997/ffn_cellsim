---
id: 42_two-point-microrheology-of-phase-separated-domains-in-lipid-
paper_n: 42
title: "Two-Point Microrheology of Phase-Separated Domains in Lipid Bilayers"
authors: "Hormel, Reyer, and Parthasarathy"
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.07.017"
paper_type: experimental
ffn_relevance: Low
ffn_themes: [membrane, tangential]
entities: [lipid-bilayer, dppc, dopc, cholesterol, giant-unilamellar-vesicle, texas-red-dhpe, liquid-ordered-domain, liquid-disordered-domain]
methods: [two-point-microrheology, single-particle-tracking, fluorescence-microscopy, passive-microrheology, two-dimensional-hydrodynamics]
measurables: [membrane-viscosity, diffusion-coefficient, displacement-correlation, domain-radius]
keywords: [two-point-microrheology, lipid-bilayer-viscosity, phase-separation, GUV, liquid-ordered, liquid-disordered, Saffman-Delbruck, HPW-model, Levine-MacKintosh, domain-diffusion]
tags: ["#membrane-biophysics", "#lipid-bilayer", "#microrheology", "#viscosity", "#tangential"]
has_transferable_params: true
---

# [42] Two-Point Microrheology of Phase-Separated Domains in Lipid Bilayers

**Tags:** #membrane-biophysics #lipid-bilayer #microrheology #viscosity #tangential

| Field | Value |
|---|---|
| Authors | Hormel, Reyer, and Parthasarathy |
| Year / Venue | 2015 / Biophysical Journal 109(4):732–736 |
| DOI / ID | 10.1016/j.bpj.2015.07.017 |
| Type | experimental |
| Pages | 5 |
| ffn_cellsim relevance | Low — 2D membrane viscometry, not cytoskeleton/ECM/adhesion mechanics; only a peripheral parameter source for a future membrane-fluidity EXTEND module |

## 1. Summary
The authors present the first two-point passive microrheology measurement of lipid bilayers. Using giant unilamellar vesicles (GUVs) of DPPC:DOPC:cholesterol that phase-separate into coexisting liquid-ordered (LO) and liquid-disordered (LD) domains, they track minority-phase domains as Brownian tracers via fluorescence microscopy. They compute both single-point diffusion coefficients (vs. domain radius, fit to the Hughes-Pailthorpe-White 2D hydrodynamic model) and two-point displacement-correlation functions (fit to the Levine-MacKintosh viscoelastic-membrane response). Single-point analysis returns the viscosity of the surrounding majority phase (LO ≈ 5× more viscous than LD), whereas two-point analysis returns an intermediate "global" effective viscosity that averages the two phases, with the two compositions becoming statistically indistinguishable. The correlation functions agree well with 2D continuum hydrodynamic theory. Conclusion: for heterogeneous membranes, two-point and single-point methods report different (length-scale-dependent) effective viscosities, and one is not "better" than the other — they probe different physics.

## 2. Problem & motivation
Membrane fluidity (viscosity) governs the timescales of lipid/protein reorganization in cell and organelle membranes, so quantifying bilayer viscosity is a prerequisite for predictive membrane-dynamics models. Existing measurements are sparse and rely on single-point methods that sample only the tracer's local neighborhood, possibly biased by the tracer itself. Two-point microrheology, which uses correlations between pairs of tracers to probe the intervening medium at the pair-separation length scale, had never been applied to lipid membranes. The paper asks whether simple 2D viscous-fluid models describe correlated domain diffusion and whether the resulting viscosity matches single-point values.

## 3. Methods / model
- **System:** GUVs (50–100 μm diameter) made by electroformation in 0.1 M sucrose; five ternary DPPC/DOPC/cholesterol compositions spanning LO-majority and LD-majority regimes. 1 mol% Texas Red DHPE fluorescent probe (partitions to LD phase); control at 0.2 mol%.
- **Imaging:** Epifluorescence, Hamamatsu ORCA CCD, Nikon TE2000, 60× objective, 10–40 fps, 296 K (room temp).
- **Tracking:** Intensity thresholding + 2D-Gaussian MLE centroid fitting; localization error < 0.07 μm. Domain radius from bilateral-filtered boundary. Only domains within 1/3 vesicle radius from the pole (in-focus) and imaged ≥ 100 frames were used. Nearest-neighbor trajectory linking. 3 to >50 domains per GUV.
- **Single-point estimator:** Covariance-based unbiased diffusion-coefficient estimator (Vestergaard et al. 2014), Eq. 1: D = ⟨Δx²ₙ⟩/(2Δt) + ⟨ΔxₙΔxₙ₊₁⟩/Δt. Goodness-of-fit χ² via periodogram comparison to free diffusion (mean reduced χ² = 1.25, indicating pure viscous diffusion).
- **Single-point viscosity model:** Hughes-Pailthorpe-White (HPW) 2D hydrodynamic drag-on-cylinder model; D(a) fit with bilayer viscosity as the single free parameter; jackknife error.
- **Two-point estimator:** Radial correlation tensor Drr(R,t) = ⟨Δrⁱ(t)Δrʲ(t)δ(R−Rᵢⱼ)⟩ (Eq. 2), with covariance-based optimal estimator (finance-derived). Fit to Levine-MacKintosh form (Eq. 3) involving Struve functions Hₙ and Bessel functions Yₙ of reduced separation β = 2RηB/η, where ηB is the surrounding aqueous viscosity.
- **Length/time scales:** μm-scale domains, fps imaging, 2D membrane embedded in 3D aqueous medium.

## 4. Key results (quantitative)
- Single-point viscosity (HPW fit): η₁ₚₜ = 0.75 ± 0.15 nPa·s·m for LD-majority phase; η₁ₚₜ = 3.90 ± 0.42 nPa·s·m for LO-majority phase (LO ~5× LD) (p.734, Fig. 3).
- Two-point viscosity: η₂ₚₜ = 2.76 ± 1.45 nPa·s·m (LD-majority) and 2.95 ± 0.61 nPa·s·m (LO-majority) — statistically indistinguishable (p.734–735, Fig. 4).
- Two-point ≈ average of the two single-point majority-phase values: 2.32 nPa·s·m (p.735).
- Across compositions spanning ~2 orders of magnitude in single-point viscosity: for LO-majority (1:1 to 1:9 DOPC:DPPC), η₂ₚₜ/η₁ₚₜ < 1, η₂ₚₜ ≈ 0.77 × η₁ₚₜ; for LD-majority (2:1 and 4:1 DOPC:DPPC), η₂ₚₜ/η₁ₚₜ > 1, η₂ₚₜ ≈ 3.07 × η₁ₚₜ (p.735, Fig. 5).
- Mean reduced χ² = 1.25 confirms pure Brownian diffusion in a viscous liquid (Fig. 2).
- Localization error < 0.07 μm (p.733). Vesicle diameters 50–100 μm (p.733).

## 5. Parameters & constants of interest
| Quantity | Value + units | Source-in-paper |
|---|---|---|
| LD-phase 2D membrane viscosity (η₁ₚₜ) | 0.75 ± 0.15 nPa·s·m | p.734, Fig. 3 |
| LO-phase 2D membrane viscosity (η₁ₚₜ) | 3.90 ± 0.42 nPa·s·m | p.734, Fig. 3 |
| Two-point effective viscosity (LD-maj) | 2.76 ± 1.45 nPa·s·m | p.734–735, Fig. 4 |
| Two-point effective viscosity (LO-maj) | 2.95 ± 0.61 nPa·s·m | p.735, Fig. 4 |
| GUV diameter range | 50–100 μm | p.733 |
| Imaging temperature | 296 K | p.733 |
| Tracking localization error | < 0.07 μm | p.733 |
| Mean reduced χ² (Brownian fit) | 1.25 | p.733, Fig. 2 |

Note: viscosities are 2D surface viscosities (units Pa·s·m = surface shear viscosity), not bulk 3D viscosities; conversion to a 3D-equivalent requires dividing by bilayer thickness (~4–5 nm), giving ~0.15–0.8 Pa·s bulk-equivalent.

## 6. Relevance to ffn_cellsim
**Classification: (f) tangential, with a thin (b) parameter-source angle for a hypothetical future membrane EXTEND.**

This is a 2D membrane-biophysics / soft-matter microrheology paper. It contains no cytoskeleton, no actin/myosin, no adhesion clutches, no ECM, no traction or single-cell spreading mechanics — the core subject matter of ffn_cellsim's mechanistic H.1→H.7 chain and its FA/clutch, motor/myosin, cortex, ECM/collagen themes. The tracers are lipid domains, and the measured property is in-plane membrane shear viscosity, a 2D hydrodynamic quantity that the project does not currently model at the particle level.

Where it could be marginally useful:
- **EXTEND membrane (H.8) parameter anchor:** If the membrane EXTEND track ever needs a lipid-bilayer in-plane viscosity, the LO/LD surface viscosities (0.75 and 3.90 nPa·s·m) are clean literature anchors. But ffn_cellsim models the *cortex* (actin) mechanically, not the lipid-bilayer hydrodynamics, so this is speculative.
- **Methods/numerics note:** The covariance-based unbiased diffusion-coefficient estimator (Eq. 1, Vestergaard et al. 2014) and the warning that linear MSD fits give poor, non-optimal estimates (and can worsen with more data) is a genuinely transferable analysis caveat for any diffusion/retrograde-flow measurement protocol in the project — relevant when extracting D or correlation coefficients from simulated trajectories.
- **Conceptual note:** The single-point-vs-two-point distinction (local vs. global effective viscosity in a heterogeneous 2D fluid) is a useful reminder that measurement protocol determines which effective property you recover — applicable in spirit to how the project defines its measurement protocols, but not a mechanism the simulator implements.

Bottom line: off-topic for a mechanistic single-cell cytoskeleton+ECM simulator. Keep as a low-priority membrane-viscosity reference and an analysis-method caveat; not a validation oracle for any current H-unit.

## 7. Limitations & caveats
- HPW model strictly applies to *solid* inclusions; domains are liquid. Authors argue corrections (Fujitani 2011/2013) are "a few percent," within uncertainty — but the model is an approximation.
- Levine-MacKintosh two-point theory assumes small, rigid inclusions, not finite-sized fluid domains; authors explicitly call for finite-domain two-point theory.
- Only near-pole, in-focus domains used; out-of-plane motion neglected. Domain growth/bulging effects argued negligible (Supporting Material).
- Model GUVs are far simpler than real cell membranes (no cytoskeleton coupling, no proteins, no curvature heterogeneity at cell scale).
- 2D surface viscosity, not a bulk material modulus — large scale gap from a fine-grained 3D particle simulator of the cortex/cytoplasm.

## 8. Key figures / tables
- **Fig. 3** — Single-point D vs. domain radius a for LD-majority (orange) and LO-majority (blue); HPW fits give η₁ₚₜ (LO ~5× LD). Core single-point result.
- **Fig. 4** — Two-point radial correlation Drr/t vs. separation; fits to Levine-MacKintosh Eq. 3; the two compositions collapse to indistinguishable η₂ₚₜ. Core two-point result.
- **Fig. 5** — η₂ₚₜ/η₁ₚₜ vs. η₁ₚₜ across compositions spanning ~2 orders of magnitude; LO-maj ratios <1, LD-maj ratios >1 (the "global averaging" signature).
- **Fig. 2** — Histogram of reduced χ² (mean 1.25) confirming pure Brownian diffusion.

## 9. Notable quotes / citable claims
- "We present, to our knowledge, the first two-point microrheological study of lipid bilayers..." (Abstract, p.732).
- "...analysis of which reveals a viscosity intermediate between those of the two lipid phases, indicative of global fluid properties rather than the viscosity of the local neighborhood of the tracer." (Abstract, p.732).
- "η₁ₚₜ = 0.75 ± 0.15 nPa·s·m for the LD majority phase, and η₁ₚₜ = 3.90 ± 0.42 nPa·s·m for the LO majority phase..." (p.734).
- "...two-point methods applied to phase-separated membranes should not be considered better than single-point methods. Rather, the latter provide insights into the viscosity of particular phases, whereas the former provide insights into the larger-scale effective viscosity of a heterogeneous fluid." (Discussion, p.735).

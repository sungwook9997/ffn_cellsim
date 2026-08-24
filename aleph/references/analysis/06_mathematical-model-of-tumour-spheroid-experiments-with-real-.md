---
id: 06_mathematical-model-of-tumour-spheroid-experiments-with-real-
paper_n: 6
title: "Mathematical model of tumour spheroid experiments with real-time cell cycle imaging"
authors: "Wang Jin et al. (Jin, Spoerri, Haass, Simpson)"
year: "2020"
venue: "bioRxiv preprint (later Bull. Math. Biol.)"
doi: "10.1101/2020.12.06.413856"
paper_type: continuum-model
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, tangential]
entities: [melanoma, c8161, wm164, fucci, oxygen, necrotic-core, spheroid]
methods: [continuum-pde, moving-boundary-problem, finite-difference, newton-raphson, confocal-live-imaging, hill-function-kinetics]
measurables: [spheroid-radius, proliferation-rate, oxygen-tension, cell-cycle-transition-time, necrotic-core-radius, g1-fraction]
keywords: [tumour-spheroid, fucci, cell-cycle, ward-and-king, moving-boundary, oxygen-diffusion, necrotic-core, melanoma, anti-mitotic-drug, reaction-advection-diffusion]
tags: ["#spheroid-mechanics", "#cell-cycle", "#continuum-model", "#moving-boundary", "#oxygen-diffusion", "#melanoma", "#tangential"]
has_transferable_params: false
---

# [6] Mathematical model of tumour spheroid experiments with real-time cell cycle imaging

**Tags:** #spheroid-mechanics #cell-cycle #continuum-model #moving-boundary #oxygen-diffusion #melanoma #tangential

| Field | Value |
|---|---|
| Authors | Wang Jin, Loredana Spoerri, Nikolas K. Haass, Matthew J. Simpson |
| Year / Venue | 2020 / bioRxiv preprint (QUT + UQ Diamantina Institute) |
| DOI / ID | 10.1101/2020.12.06.413856 |
| Type | continuum-model (reaction-advection-diffusion PDE, moving boundary) |
| Pages | 40 |
| ffn_cellsim relevance | Low — tissue/population-scale spheroid growth PDE; no cytoskeleton, no single-cell mechanics |

## 1. Summary
The authors extend the Ward & King (1997) avascular-tumour continuum PDE framework to interpret 3D tumour spheroid experiments imaged with FUCCI (fluorescent ubiquitination-based cell cycle indicator), which colours G1-phase cells red and S/G2/M-phase cells green. The spheroid is treated as a spherically-symmetric, three-subpopulation continuum: volume fractions of red (G1) cells r(x,t), green (S/G2/M) cells g(x,t), and non-fluorescent dead cells m(x,t). Cell-cycle transition rates (G1→S/G2/M, S/G2/M→G1) and the death rate are oxygen-dependent Hill functions; oxygen diffuses in from the free surface and is consumed by living cells. A net-volume-gain velocity field advects all subpopulations and drives a moving outer boundary l(t). The nondimensionalised system is solved by boundary-fixing transform + upwind/central finite differences + backward-Euler + Newton-Raphson (MATLAB on GitHub). The model reproduces FUCCI-observed intratumoural structure (outer freely-cycling rim, intermediate G1-arrested quiescent shell, central necrotic core), differences between melanoma cell lines (C8161 vs WM164), and the effect of a MEK-inhibitor anti-mitotic drug (U0126) modelled as time-decaying G1→S transition rate.

## 2. Problem & motivation
Standard spheroid growth assays measure only outer radius vs time and reveal nothing about the spatial distribution of cell-cycle status. FUCCI labelling exposes a strong correlation between radial position and cell-cycle phase (proliferating rim, G1-arrested middle, necrotic core), but no existing tumour-spheroid mathematical model was designed to interpret FUCCI data. The paper fills that gap so that the spatially-resolved cell-cycle readout can be quantitatively modelled and used to mimic drug action.

## 3. Methods / model
- Model class: spherically-symmetric continuum reaction-advection-diffusion PDEs with a moving (free) outer boundary; quasi-steady-state oxygen. An extension of Ward & King (1997).
- Governing equations (nondimensional, Eqs. 15-18): ∂r/∂t + v∂r/∂x = 2Kg(c)g − [...]r ; ∂g/∂t + v∂g/∂x = Kr(c)r − [...]g ; (1/x²)∂(x²v)/∂x = b(c)g − d(c)r (velocity from net volume change); (1/x²)∂(x²∂c/∂x) = k(c)g + γ(c)r (oxygen diffusion-consumption). Factor of 2 in the r equation = each S/G2/M→G1 division yields two daughter G1 cells.
- Kinetic rates are Hill functions of local oxygen c (Eqs. 6-8 / 19-21): Kg = c^m1/(cg^m1+c^m1), Kr = B·c^m2/(cr^m2+c^m2), Kd = C(1 − σc^m3/(cd^m3+c^m3)). Critical oxygen thresholds cg, cr, cd; Hill exponents m1,m2,m3.
- No-void closure: VL·n + VD·m = 1. Dead-cell volume VD < living-cell volume VL (δ = VD/VL) provides the volume-loss mechanism. No carrying-capacity parameter needed.
- Length/time scales: initial spheroid radius l(0) ≈ 250 µm (some figs note ~110 µm scale bars); characteristic time set by Ã ≈ 0.1 /h so 1 nondimensional time unit = 10 h. Simulations run to t ≈ 2.4 (24 h) up to t = 14.4 (144 h).
- Numerics: boundary-fixing transform x = l(t)ξ; central differences for diffusion, upwinding for advection; backward (implicit) Euler in time; Newton-Raphson for the nonlinear algebraic system; δt = δξ = 1×10⁻³, tolerance ϵ = 1×10⁻⁶ (grid-independent).
- Experimental side (motivation only, not new data here): FUCCI-transduced melanoma cell lines C8161 and WM164 grown as 3D spheroids; phase-contrast + confocal microscopy; cell-cycle phase durations from Haass et al. 2014; drug = U0126 MEK-inhibitor applied at day 3.

## 4. Key results (quantitative)
- Cell-cycle rate estimates for C8161 (2D, plentiful O₂, from Haass et al. 2014 durations, Fig 4): 0.08 < Ã < 0.25 /h (S/G2/M→G1, green→red) and 0.08 < B̃ < 0.12 /h (G1→S/G2/M, red→green); representative Ã ≈ 0.1 /h, ratio B = B̃/Ã ≈ 1.5 (p.14, Fig 4b-c).
- Experimental spheroid expands by a factor ≈ 1.3 in radius over 24 h (Fig 2); model reproduces l(2.4) ≈ 1.3 (p.19, Fig 6c).
- Single cell-line structure (C8161): with cr = cg = 0.1 the model gives a necrotic core (x < 0.8) but r/n ≈ 0.5 everywhere (no G1-arrested shell); with cr = 0.7, cg = 0.1 it produces the three-zone structure — necrotic core |x| < 0.6, G1-arrested shell 0.6 < |x| < 1.1, proliferating rim 1.1 < |x| < 1.3 (Figs 6-7, p.16-21).
- Multi-cell-line fit (Fig 8): WM164 ~75% G1 (red) at surface with gradual rise; C8161 ~45% G1 at surface with rapid rise. Reproduced with B=1, cr=0.9, m=3 (WM164) and B=3, cr=0.7, m=6 (C8161); shared C=1.8, σ=1, δ=0.9, β=10, cg=0.1, cd=0.5.
- Drug (U0126) simulation: B(t) = B·e^{−λ(t−td)} for t ≥ td, with td = 7.2 (≈72 h) and λ = 10 (Eq 24-25). Treated spheroid nearly stops growing and becomes almost entirely G1 (red) after dosing (Fig 9).
- Sensitivity (Appendix D, Fig 11): benchmark l(2.4) = 1.365; ±10% sweeps show l(2.4) is most sensitive to the critical oxygen thresholds cr and cd.

## 5. Parameters & constants of interest
| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Initial spheroid radius l̃(0) | ≈ 250 µm | p.14 |
| C8161 S/G2/M→G1 rate Ã | 0.08–0.25 /h (rep. 0.1 /h) | Fig 4b, p.14 |
| C8161 G1→S/G2/M rate B̃ | 0.08–0.12 /h | Fig 4b, p.14 |
| Rate ratio B = B̃/Ã | ≈ 1.5 (1–3 across cell lines) | Fig 4c / Fig 8 |
| Characteristic time (1 nondim unit) | ≈ 10 h | p.14 |
| Spheroid radial expansion | ×1.3 over 24 h | Fig 2 / Fig 6c |
| Dead:live cell volume ratio δ = VD/VL | 0.9 (assumed) | p.16 |
| Hill exponents m1=m2=m3 | 3, 6, or 10 (chosen) | Figs 6-9 |

Note: All of the above are tissue/population-level proliferation, death, and oxygen-transport parameters. None are mechanical (no modulus, traction, bond force, drag, motor force, or cytoskeletal constant). They are NOT transferable into a fine-grained mechanistic single-cell cytoskeleton+ECM simulator as an oracle. has_transferable_params = false.

## 6. Relevance to ffn_cellsim
This is **(d) spheroid/multicellular-scale context**, and largely **(f) tangential**, for ffn_cellsim. It is a continuum population-balance PDE — exactly the class of "paper-model / lumped-mechanism" tissue model that ffn_cellsim deliberately does NOT use as a runtime mechanism. There are no cytoskeletal filaments, motors, adhesion clutches, ECM cross-links, membranes, or single-cell mechanics anywhere in the model; "cells" are volume fractions in a Hill-kinetic reaction-advection-diffusion field. None of the H-unit mechanisms (H.1 cortex, H.2, H.3, H.5 FA/clutch, H.7, EXTEND membrane/nucleus/cytoplasm) map onto it.

Where it has marginal value:
- **Spheroid-scale context only.** ffn_cellsim's validation overlay target is an MCF7-spheroid-on-pV4D4/collagen-I platform; this paper is a vocabulary/landmark reference for what a tissue-scale spheroid model looks like (necrotic core, quiescent G1-arrested shell, proliferating rim, oxygen-limited proliferation), useful when someone eventually bridges single-cell mechanics up to multicellular cohesion at the spheroid scale. It is the kind of coarse continuum the project explicitly stands apart from, so it also serves as a contrast/anti-pattern reference (Greenspan/Ward-King lineage = the lumped tissue PDE the project rejects).
- **No mechanism, parameter, oracle, or numerics transfer.** The numerics (boundary-fixing transform, upwind FD, backward Euler, Newton-Raphson moving-boundary solve) are standard 1D radial PDE techniques unrelated to the BAOAB Langevin particle integrator used here. The proliferation/death/oxygen parameters are not used anywhere in a single-cell cytoskeleton model.

Plainly: off-topic for a mechanistic single-cell cytoskeleton+ECM simulator. Keep as Low-relevance spheroid-scale context; do not treat as a validation oracle or parameter source.

## 7. Limitations & caveats
- Spherically symmetric, continuum, no individual cells — cannot represent cell-cell adhesion mechanics, traction, cytoskeletal forces, or ECM. A ×10⁴-cell coarse-grain with no mechanics; the opposite end of the modelling spectrum from ffn_cellsim.
- "Carefully chosen" parameter values are illustrative, not fitted — the authors explicitly state the goal is qualitative replication, not parameter estimation (p.15-16); only Ã, B̃ are data-anchored.
- Oxygen is the sole nutrient/regulator; mechanical stress, contact inhibition, and matrix mechanics are absent. Velocity is purely volumetric (proliferation minus death), not force-balance-derived.
- Melanoma-only experimental motivation (C8161, WM164); MCF7/breast not addressed. Yellow (early-S) phase neglected; FUCCI4 left for future work.
- Preprint, not peer-reviewed at the captured version.

## 8. Key figures / tables
- **Figure 6/7** (p.18-20): numerical r, g, n, c, v, l(t) profiles — the core result, showing necrotic core + G1-arrested shell + proliferating rim emerging from oxygen-limited Hill kinetics.
- **Figure 4** (p.15): experimental C8161 cell-cycle phase durations → rate estimates Ã, B̃, B — the only data-anchored parameters.
- **Figure 8** (p.23): cross-cell-line (WM164 vs C8161) G1-fraction-vs-depth, experiment vs model.
- **Figure 9** (p.25): U0126 anti-mitotic drug simulation via time-decaying B(t); radius arrest + G1 dominance.
- **Figure 11 / Appendix D** (p.34): ±10% sensitivity — l(2.4) most sensitive to oxygen thresholds cr, cd.

## 9. Notable quotes / citable claims
- "We extend the mathematical framework originally proposed by Ward and King (1997) to develop a new mathematical model of FUCCI-labelled tumour spheroid growth... composed of three subpopulations: (i) living cells in G1 phase that fluoresce red; (ii) living cells in S/G2/M phase that fluoresce green; and (iii) dead cells that do not fluoresce." (Abstract, p.1-2)
- "We assume that the rates at which cells pass through different phases of the cell cycle, and the rate of cell death, depend upon the local oxygen concentration in the spheroid." (Abstract, p.2)
- "Data in Figure 4(b) suggest that 0.08 < Ã < 0.25 /h and 0.08 < B̃ < 0.12 /h. For simplicity we take a representative estimate of Ã ≈ 0.1 /h, which means that each unit duration of time in the nondimensional model corresponds to 10 h in the experiment." (p.14)
- "A useful feature of this nutrient-dependent cell proliferation and nutrient-dependent cell death framework is that it avoids the need for specifying a carrying capacity density as a model parameter." (p.10)

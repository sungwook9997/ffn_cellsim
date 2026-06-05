---
id: 02_integrating-simulation-and-experimental-validation-of-nutrie
paper_n: 2
title: "Integrating simulation and experimental validation of nutrient-limited growth in breast cancer spheroids"
authors: "Celia Nieto, Álvaro González-Garcinuño, Eva Martín del Valle"
year: "2026"
venue: "European Journal of Pharmaceutical Sciences 216 (2026) 107370"
doi: "10.1016/j.ejps.2025.107370"
paper_type: continuum-model
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, tangential]
entities: [bt-474, mcf7, breast-cancer, spheroid, glucose, oxygen, her2]
methods: [finite-element, multiphysics-comsol, gompertz-growth, reaction-diffusion, neo-hookean-hyperelastic, darcy-flow, confocal-microscopy, hplc, h-and-e-histology, immunostaining]
measurables: [spheroid-diameter, spheroid-volume, glucose-consumption, oxygen-concentration, porosity, necrotic-core-fraction, cell-number, elastic-modulus]
keywords: [breast cancer spheroid, MCTS, COMSOL multiphysics, Gompertz growth, nutrient diffusion, necrotic core, glucose threshold, porosity evolution, BT-474, digital twin]
tags: ["#spheroid-mechanics", "#continuum-model", "#nutrient-diffusion", "#tangential", "#parameter-source"]
has_transferable_params: true
---

# [2] Integrating simulation and experimental validation of nutrient-limited growth in breast cancer spheroids

**Tags:** #spheroid-mechanics #continuum-model #nutrient-diffusion #tangential #parameter-source

| Field | Value |
|---|---|
| Authors | Celia Nieto, Álvaro González-Garcinuño, Eva Martín del Valle (Univ. Salamanca / IBSAL, Spain) |
| Year / Venue | 2026 / European Journal of Pharmaceutical Sciences 216, art. 107370 |
| DOI / ID | 10.1016/j.ejps.2025.107370 (open access, CC BY) |
| Type | Continuum multiphysics model (COMSOL FEM) + experimental validation |
| Pages | 11 |
| ffn_cellsim relevance | Low — spheroid-scale, continuum (lumped) tumor-growth model; opposite end of the abstraction axis from ffn_cellsim, useful only as macro-scale context + a few tissue-mechanics constants |

## 1. Summary
The authors build a COMSOL Multiphysics FEM model of HER2-positive BT-474 breast-cancer multicellular tumor spheroids (MCTS) that couples Gompertzian growth (via a prescribed expanding mesh with hyperelastic smoothing), Navier-Stokes/Darcy fluid flow, and reaction-diffusion transport of glucose and oxygen with porosity-dependent effective diffusion. The model is parameterized and validated against experimental measurements of spheroid diameter (CLSM), glucose consumption (HPLC), porosity (H&E histology image thresholding), cell number/viability (trypan-blue counting), and necrotic-core fraction over a 7-day culture. It reproduces growth kinetics (MAE 6.9 %), glucose uptake (MAE 1.5 %), and necrotic-core development, and identifies a glucose concentration threshold of ~0.08 mM as critical for necrosis onset in this cell line. This is a tissue-/organ-scale continuum digital-twin effort, mechanistically orthogonal to ffn_cellsim's single-cell, explicit-particle approach.

## 2. Problem & motivation
3D MCTS are a relevant in-vitro model for breast cancer and drug screening, but validated computational tools that simulate their growth + microenvironment are scarce. Existing in-silico models use oversimplified geometries, assume constant cell viability or uniform diffusion, omit progressive porosity/deformability, or are published without experimental validation. The goal is an experimentally grounded multiphysics model that couples growth, mechanical deformation, nutrient transport, and porosity for one HER2+ cell line, as a step toward "digital twin" tumor models for preclinical drug evaluation.

## 3. Methods / model
Two domains: a fixed microwell (96-well U-bottom geometry: 10.2 mm height, 7.1 mm upper diameter, hemispherical curvature from 6.3 mm) containing a growing spheroid (initial 0.1 mm diameter). Tetrahedral mesh, >605k elements, average skewness quality 0.663.
- **Growth**: spheroid expands via a prescribed moving-boundary mesh velocity derived from the Gompertz model (Eqs. 5–6); mesh deformation uses hyperelastic smoothing with a Neo-Hookean stored-energy function (Eqs. 3–4), driven by iterative elastic-energy minimization; smoothing velocity v_mbs = δ_mbs·|v|·h·H (Eq. 2).
- **Fluid**: Navier-Stokes + continuity in the well (Eqs. 7–8); inside the spheroid treated as isotropic porous medium → Darcy's law (Eq. 9) with Kozeny-Carman permeability k = d_c²·ε_s / [180(1−ε_s)²] (Eq. 10). Working fluid = water.
- **Transport of diluted species** (glucose, oxygen): porosity-weighted mass balance with accumulation/convection/diffusion/internalization (Eq. 11); Fick's law with effective diffusion D_eff,i = (ε_s/τ_s)·D_i, tortuosity τ_s = ε_s^(−1/3) (Millington-Quirk, Eq. 13). Internalization R_int = FVC·k_int·φ (Eq. 14), where FVC (fraction viable cells) declines linearly with time, FVC = 0.9067 − 0.0855·t(d), R²=0.94 (Eq. 15).
- **Parameter estimation**: Gompertz parameters a, a' and glucose k_int fitted by Excel-Solver GRG nonlinear minimizing MAE (Eq. 16). Solver = GMRES, residual tol 10⁻⁴, max 100 iters/timestep.
- **Validation metrics**: chi-square (Eq. 17) and RMSEA (Eq. 18) vs experimental.
- **Experiments**: BT-474 MCTS by liquid overlay (2000 cells/well), cultured ≤7 d; CLSM with calcein-AM/PI for size + necrotic core; HPLC for glucose in DPBS + 0.15 mM glucose minimal medium; H&E for porosity; pan-cytokeratin (AE1/3) immunostaining for aggressiveness; trypan-blue automated counting for cell number/viability.

## 4. Key results (quantitative)
- Spheroid volume grew ~0.08 mm³ → ~0.45 mm³ over 5 days; most pronounced diameter growth days 4–6 (near-doubling); plateau days 6–7 (§3.1).
- Cell number plateaued at ~215,000 cells per spheroid after days 4–6 (§3.1).
- Necrotic core first appears day 3 onward; none in first 2 days; expands ~exponentially (§3.1, Fig. 1).
- Porosity: 0.0145 ± 0.0057 (day 2) → 0.0543 (d3) → 0.1149 (d5) → 0.4443 (d6) → 0.4676 (d7); approaches ~0.45 by final days (Table 1).
- Glucose consumption rose exponentially from day 3; less steep days 6–7 (Fig. 2, §3.2).
- Growth-model fit: MAE 6.9 %; chi-square 0.023 (< critical 7.815); RMSEA 0.00. Fitted Gompertz params a = 8.1 day⁻¹, a' = 0.65 day⁻¹ (§3.5).
- Experimental vs predicted diameters (mm), days 0–7: 0.100/0.090, 0.155/0.160, 0.260/0.240, 0.395/0.360, 0.504/0.561, 0.842/0.764, 0.863/0.920 (Table 2).
- Glucose-model fit: MAE 1.5 %; chi-square 0.009 (< critical 3.841); RMSEA 0.00. Fitted standard glucose consumption rate 2·10⁻⁸ mol·m⁻²·s⁻¹ (§3.6) — consistent with literature ~10⁻² mol·m⁻³·s⁻¹ → ~2.5·10⁻⁸ mol·m⁻²·s⁻¹ for a 15 µm cell.
- Simulated core glucose dropped <0.07 mM by day 6, <0.03 mM by day 7 (§3.6, Fig. 6).
- Oxygen: uniform >1 mM days 0–4; core hypoxia by days 6–7, min ~0.5 mM — above the necrosis threshold (§3.7, Fig. 7); model assumes no O₂ exchange (closed microplates).
- Necrotic-core volume fraction (exp vs theoretical, error): d5 0.200/0.126 (37 %), d6 0.291/0.271 (6.87 %), d7 0.607/0.496 (18.3 %) at Jiang 0.06 mM glucose threshold (Table 3); adjusting threshold to 0.08 mM cut errors to 9.4 %/2.1 %/6.3 % (days 5/6/7) → proposed BT-474 necrosis glucose threshold ~0.08 mM (§3.8).

## 5. Parameters & constants of interest
| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Glucose molecular diffusion coeff (in DMEM) | 5.9·10⁻¹⁰ m²/s | lit. (McMurtrey 2016; Suhaimi 2015), Eq. 12 region |
| Oxygen molecular diffusion coeff (in DMEM) | 1.5·10⁻⁹ m²/s | lit., §2.7.4 |
| Oxygen internalization const k_int | 2.5·10⁻¹⁸ mol O₂·cell⁻¹·s⁻¹ | lit. (Wagner 2011) |
| Glucose consumption rate (fitted) | 2·10⁻⁸ mol·m⁻²·s⁻¹ | this study, GRG fit (§3.6) |
| Cell diameter d_c | 15 µm | lit. (Lindfors 2022) + own imaging |
| Breast-tissue 1st Lamé λ | 4–50 kPa | Krouskop 1998; Samani 2003 |
| Breast-tissue shear modulus µ | 0.2–3.5 kPa | Krouskop 1998; Samani 2003 |
| Effective diffusion | D_eff = (ε_s/τ_s)·D, τ_s = ε_s^(−1/3) | Millington-Quirk 1961, Eq. 13 |
| Permeability (Kozeny-Carman) | k = d_c²ε_s / [180(1−ε_s)²] | Eq. 10 |
| FVC(t) | 0.9067 − 0.0855·t(d), R²=0.94 | this study, Eq. 15 |
| Glucose necrosis threshold (BT-474) | ~0.08 mM | this study (§3.8); cf. Jiang 0.06 mM |
| Necrosis thresholds (Jiang multiscale) | O₂ <0.02 mM, glucose <0.06 mM, lactate >8 mM | Jiang et al. 2005 |
| Initial spheroid diameter | 0.1 mm | protocol (Nieto 2021) |
| Cells seeded / spheroid | 2000 cells/well | §2.2 |

## 6. Relevance to ffn_cellsim
**Low / tangential.** This is a continuum, lumped, tissue-scale FEM model — exactly the "abstracted/proxy mechanism" category that ffn_cellsim deliberately rejects at runtime (Gompertz growth, Darcy porous flow, reaction-diffusion PDEs, prescribed-mesh deformation). It is not a mechanism reference for any H-unit and provides no fine-grained particle/bond physics.

Where it touches the project:
- **Spheroid-scale context** for the MCF7-spheroid-on-pV4D4/collagen-I validation overlay track: it is a sibling effort modeling BC spheroid growth/necrosis, and documents the macro behaviors (size kinetics, necrotic-core fraction, porosity rising to ~0.45) that the platform's eventual multi-cell scale would sit beneath. Note: this paper uses BT-474 (HER2+), not the project's MCF7 line, so it is context, not a direct overlay target.
- **Parameter source (overlay-only, never fit-to)**: breast-tissue elastic constants (Lamé λ 4–50 kPa, shear µ 0.2–3.5 kPa from Krouskop 1998 / Samani 2003) are independent literature anchors that could inform an ECM/tissue-stiffness band in the EXTEND track (H.8/9/10) or sanity overlays — but they are bulk-tissue moduli, not single-cell cortical/cytoskeletal constants, so use with care. Cell diameter 15 µm and glucose/O₂ diffusion coefficients are standard literature values that may seed a cytoplasm/medium transport context if metabolism is ever added (currently out of scope).
- **Not a validation oracle** in the ffn_cellsim sense: the closed-forms here (Gompertz, Kozeny-Carman, Millington-Quirk, Neo-Hookean) describe tissue growth/transport, not the single-cell mechanobiology mechanisms (Bell-Evans, Hill, Stam-Hocky, catch-bond) that the oracles cover.

Bottom line: file under spheroid-scale-context; harvest only the cited primary-literature tissue-mechanics constants if/when EXTEND needs a macro stiffness band, and treat them as overlay anchors per the literature-first / no-fit rules.

## 7. Limitations & caveats
- Single cell line (BT-474, HER2+), avascular, nutrient-limited only; authors flag this explicitly.
- Cell signaling, pH gradients, and lactate production not modeled (may shift outcomes); small geometric deviations from ideal spheres ignored to save compute.
- Glucose measured in an artificial minimal medium (DPBS + 0.15 mM glucose), not the standard DMEM+FBS growth medium — chosen to reduce HPLC sensor noise, so metabolic numbers are condition-specific.
- Oxygen never validated experimentally (tracked only via literature rates, closed microplate assumption → no O₂ exchange); the O₂ necrosis conclusion rests entirely on simulation.
- Necrosis threshold 0.08 mM is back-fitted to this cell line's experimental necrotic-core fractions (deviates from Jiang 0.06 mM) — a tuned, line-specific constant, not a transferable mechanism.
- FVC linear regression and Gompertz params are phenomenological fits, not mechanistic; near-perfect RMSEA 0.00 reflects a small validation dataset (one value per culture day).
- Continuum/lumped throughout: no explicit cells, no cytoskeleton, no adhesion — incompatible with ffn_cellsim's runtime-mechanism rule.

## 8. Key figures / tables
- **Fig. 1** — CLSM calcein/PI images (A) + spheroid-volume histogram with necrotic-core fraction (B) + viable/non-viable cell-number histogram (C).
- **Fig. 2** — % glucose consumed by BT-474 spheroids over 7 days (minimal medium).
- **Fig. 3** — H&E (A) and H&E + pan-cytokeratin (B) micro-sections, days 2–7.
- **Fig. 4** — Simulated 3D spheroid growth with mesh-quality (skewness) coloring.
- **Fig. 5** — Glucose concentration in microwell bulk: HPLC (red dots) vs model (blue line).
- **Fig. 6** — Spatial glucose-concentration distribution over time (full geometry + zoomed bottom).
- **Fig. 7** — Oxygen concentration evolution inside spheroids over time.
- **Table 1** — Porosity vs culture day (0.0145 → 0.4676).
- **Table 2** — Experimental vs predicted diameters, days 0–7.
- **Table 3** — Necrotic-core volume fraction, experimental vs theoretical (errors).

## 9. Notable quotes / citable claims
- "a glucose concentration threshold of ~0.08 mM was identified as critical for necrosis" (Abstract / §3.8) — proposed BT-474-specific necrosis onset.
- "to the best of our knowledge, this is the first study to calculate glucose uptake rate for HER2+ MCTS" (§3.6).
- Fitted standard glucose consumption rate "2·10⁻⁸ mol·m⁻²·s⁻¹" (§3.6); growth fit MAE 6.9 %, glucose fit MAE 1.5 %.
- Breast-tissue mechanics: "values for the first Lamé parameter (λ) range from 4 to 50 kPa, and shear modulus (µ) values range from 0.2 to 3.5 kPa" (§2.7.2, citing Krouskop 1998 / Samani 2003).
- Jiang et al. (2005) necrosis thresholds: "oxygen concentrations below 0.02 mM, glucose concentrations below 0.06 mM, and lactate (waste) concentrations above 8 mM" (§3.8).
- Methodological stance — "the importance of tailoring computational models to experimental evidence rather than relying solely on literature-based parameters" (§4): aligns in spirit with ffn_cellsim's literature-first + experimental-overlay discipline, though their fit-to-data approach is the inverse of the project's no-fit rule.

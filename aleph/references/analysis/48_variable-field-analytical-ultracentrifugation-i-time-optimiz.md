---
id: 48_variable-field-analytical-ultracentrifugation-i-time-optimiz
paper_n: 48
title: "Variable-Field Analytical Ultracentrifugation: I. Time-Optimized Sedimentation Equilibrium"
authors: "Ma, Metrick, Ghirlando, Zhao & Schuck"
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.07.015"
paper_type: methods-software
ffn_relevance: Low
ffn_themes: [numerics/methods, tangential]
entities: [bovine-serum-albumin, soybean-trypsin-inhibitor, alpha-chymotrypsin, dna, nucleosome-positioning-sequence]
methods: [analytical-ultracentrifugation, sedimentation-equilibrium, lamm-equation, finite-element, simplex-optimization, simulated-annealing]
measurables: [buoyant-molar-mass, sedimentation-coefficient, diffusion-coefficient, frictional-ratio, binding-constant, equilibration-time]
keywords: [sedimentation-equilibrium, analytical-ultracentrifugation, Lamm-equation, overspeeding, rotor-speed-schedule, buoyant-molar-mass, time-varying-centrifugal-field, SEDFIT, Boltzmann-distribution, protein-interactions]
tags: ["#analytical-ultracentrifugation", "#sedimentation-equilibrium", "#methods-software", "#lamm-equation", "#off-topic", "#solution-biophysics"]
has_transferable_params: false
---

# [48] Variable-Field Analytical Ultracentrifugation: I. Time-Optimized Sedimentation Equilibrium

**Tags:** #analytical-ultracentrifugation #sedimentation-equilibrium #methods-software #lamm-equation #off-topic #solution-biophysics

| Field | Value |
|---|---|
| Authors | Ma, Metrick, Ghirlando, Zhao & Schuck (NIH/NIBIB & NIDDK) |
| Year / Venue | 2015 / Biophysical Journal 109(4) 827–837 |
| DOI / ID | 10.1016/j.bpj.2015.07.015 |
| Type | methods-software |
| Pages | 11 |
| ffn_cellsim relevance | Low — solution-biophysics instrumentation method (AUC) with no cytoskeleton/ECM/cell-mechanics content |

## 1. Summary
This is an analytical-ultracentrifugation (AUC) methods/software paper. Sedimentation equilibrium (SE) AUC measures macromolecular buoyant molar masses and reversible binding thermodynamics in solution, but conventionally takes days to reach equilibrium. The authors develop "time-optimized SE" (toSE): they compute numerical Lamm-equation solutions for sedimentation in a *time-varying* centrifugal field, then optimize a parameterized rotor-speed schedule (over-/under-speeding phases) to minimize the time to equilibrium while capping transient over-concentration at the cell base. They package this in a new program, TOSE (integrated into SEDFIT/SEDPHAT), that emits speed-step/.EQU files directly loadable into Beckman XLA/I instruments. They show up to ~10× faster equilibration in simulations and validate experimentally on BSA, DNA, and an interacting SBTI/α-chymotrypsin system, recovering molar masses and a binding constant consistent with constant-speed controls.

## 2. Problem & motivation
SE-AUC is a first-principles "gold standard" for buoyant molar mass and for the thermodynamics of reversible macromolecular interactions in free solution, but its key drawback is the very long time (often days) to reach thermodynamic equilibrium across the solution column. Classical single-step "overspeeding" can shorten this but is fragile: it requires accurate prior knowledge of the species' mass/frictional coefficient and can cause severe over-concentration or pelleting at the cell base, sometimes *lengthening* equilibration. The paper's goal is to remove these two obstacles using modern computation, making SE fast and robust enough to be routine and to enable studies of poorly stable macromolecules.

## 3. Methods / model
- **Governing equation:** the Lamm equation (Eq. 1) for radial/temporal concentration evolution c(r,t) under a *time-varying* angular velocity ω(t), with sedimentation coefficient s and diffusion coefficient D in a sector-shaped column from meniscus m to base b.
- **Equilibrium target:** the Boltzmann distribution ceq(r,ω) (Eq. 2) with buoyant molar mass Mb = sRT/D.
- **Numerics:** finite-element adaptive-grid solver (Brown & Schuck 2008, ref. 35), with propagation matrices updated according to |dω/dt|; ω(t) modeled as constant-speed segments joined by constant-acceleration/deceleration ramps. Numerical precision typically set to 1/10 of experimental data precision.
- **Optimization:** closeness-to-equilibrium measured by RMSD Deq(t) between the Lamm solution and the Boltzmann profile over the accessible radial range; equilibration time Teq = first time Deq < threshold ε·c0 (ε = 0.001–0.007). Rotor-speed profile parameterized (exponential decay toward ω_eq modulated in N=5–10 steps), optimized by repeated simplex + optional simulated annealing, with penalty functions on lingering-near-equilibrium traces and on base over-concentration a = max_t[c(b,t)/ceq(b)] (limited to 1.2–1.5 for proteins, up to 2.0 for DNA). Segments ≥10 min, speeds in 100 rpm increments to be machine-realizable.
- **Software:** TOSE (v1.0), GUI, exports speedsteps.txt and .EQU files; SEDFIT extended to analyze scans acquired at multiple/time-varying rotor speeds with rotor-stretching correction and time-invariant noise models; SEDPHAT for global interacting-system fits with mass-conservation + mass-action constraints.
- **Experimental:** Beckman XLA/I AUC, Epon double-sector centerpieces, absorbance scans at 230/250/260/280 nm. Samples: BSA (~0.5 mg/mL), SBTI + α-chymotrypsin (CT) mixtures, and a 162 bp dsDNA fragment.
- **Scales:** solution columns 3.5–5.5 mm; rotor speeds ~6,000–50,000 rpm; equilibration times hours to days. No cellular/cytoskeletal length or time scales involved.

## 4. Key results (quantitative)
- 65 kDa protein with 15% dimer, 4.5 mm column, 20 °C, final 15,000 rpm: conventional SE Teq ≈ 28 h vs. unconstrained toSE 4.5 h (a ≈ 5×), constrained toSE 9 h (a < 1.6) and 13 h (a < 1.07) (Fig. 1, ε = 0.007).
- Classical single-step overspeeding (Chatelier params, 50,000 rpm for 43 min) reduced 28 h → 9.5 h but transiently reached ~5× final equilibrium / ~33× loading concentration at base (Results, p. 830).
- Robustness to mass error: profile optimized for 65 kDa (a = 1.2) gave Teq 10.5 h vs. 30 h constant-speed; applying same profile to half/twice mass gave 15 h / 11.5 h — a 2× mass error still saves substantial time (p. 831–832).
- BSA experimental toSE (Fig. 1C schedule, 41,600 → 15,000 rpm over 5 h): equilibrium at ~7 h vs. ~20 h constant-speed control, same buoyant molar mass; c(s) distribution recoverable from first 3 h of data (Fig. 3).
- 162 bp dsDNA: s20,w = 5.02 S, estimated mass ~103 kDa (107.730 kDa from v̄ = 0.55 mL/g), frictional ratio f/f0 = 2.89; conventional SE ~70 h vs. toSE equilibrium at ~20–30 h (Fig. S1).
- 5.5 mm long column BSA at 7,800 rpm: ~60 h conventional vs. ~15 h toSE by molar mass (Fig. S3).
- Multi-speed sequence (6,000/10,000/15,000 rpm): conventional ≥4 days vs. toSE ~1.5 days (Fig. S4).
- SBTI/CT interacting system (sequential 8,000/15,000/25,000 rpm, 10 °C): apparent masses SBTI 19.2 kDa (toSE)/19.1 kDa (const), CT 23.3/23.0 kDa; two-equivalent-sites binding constant log(K1) = 5.91 (toSE) vs. 5.84 (const) — consistent (Fig. 4, p. 834).
- Overall claim: SE expedited "by at least a factor of 2 and up to a factor of 10"; moderate 20–50% over-concentration gives 3–5× savings (Discussion, p. 835–836).

## 5. Parameters & constants of interest
None (no transferable constants for ffn_cellsim). The numbers in this paper (sedimentation coefficients, buoyant molar masses, frictional ratios, rotor speeds, equilibration times) describe whole-macromolecule transport in a centrifugal field and instrument-specific schedules; none map to cytoskeletal filament/motor/clutch/ECM parameters or to a fine-grained particle simulator. For completeness, the only physically generic items are textbook relations already in the simulator's toolkit (Boltzmann distribution Mb = sRT/D), not novel anchors.

| Quantity | Value + units | Source-in-paper |
|---|---|---|
| 162 bp dsDNA sedimentation coeff. | s20,w = 5.02 S | Methods, p. 831 |
| 162 bp dsDNA frictional ratio | f/f0 = 2.89 | Methods, p. 831 |
| BSA monomer mass (apparent) | recovered ~SE-standard (no novel value) | Fig. 3 |
| SBTI / CT apparent masses | 19.1–19.2 / 23.0–23.3 kDa | p. 834 |
| SBTI–CT binding constant | log(K1) ≈ 5.84–5.91 | p. 834 |

(All are solution-AUC observables, not usable as ffn_cellsim oracles.)

## 6. Relevance to ffn_cellsim  <-- MOST IMPORTANT
**Off-topic; Low relevance.** This is a solution-biophysics instrumentation/software method for analytical ultracentrifugation. ffn_cellsim is a fine-grained, mechanistic HOOMD simulator of single-cell cytoskeleton + adhesion + ECM mechanics; this paper contains no cytoskeletal filaments, motors, adhesion clutches, ECM cross-links, cell lines used mechanically, traction/modulus measurements, or spheroid-scale context. It does not function as (a) a validation oracle, (b) a parameter source, (c) a mechanism reference, (d) spheroid/multicellular context, or (e) a numerics method we would adopt.

The single thread that touches our world is methodological by analogy only: like ffn_cellsim, it solves a PDE master equation (the Lamm transport equation) numerically with a finite-element adaptive-grid scheme and explicitly contrasts FE accuracy against "numerical diffusion errors" from "simpler differentiation schemes" — the same spirit as our preference for the Leimkuhler-Matthews BAOAB integrator and care about discretization artifacts. But the Lamm equation (advection-diffusion in a centrifugal field) is unrelated to BAOAB Langevin dynamics, and TOSE/SEDFIT solve a different problem on a different scale. There is no mechanistic, parameter, or validation transfer.

Most plausibly this PDF landed in the bundle because the author (Peter Schuck, NIH) is prominent in macromolecular biophysics and the bundle was harvested by journal/keyword, not because it informs a mechanistic cell simulator. Recommend tagging `tangential` and not wiring it into any H-unit (cortex, FA/clutch, motor/myosin, ECM/collagen, membrane/nucleus/cytoplasm EXTEND, spheroid context, or numerics) — none apply.

## 7. Limitations & caveats
- Entirely a solution-phase, whole-molecule transport method; no spatial cell structure, no mechanics, no fine-grained particles — a complete scale and domain gap from a single-cell cytoskeleton simulator.
- Even within its own domain: requires approximate prior knowledge of sedimentation parameters; optimization is ill-conditioned with many local minima; over-concentration constraint is the critical user-chosen parameter; toSE predictions valid only for double-sector centerpieces; temperature-as-a-variable attempts failed due to convection.
- Reported "gains" are operational (time savings), not new physical constants.

## 8. Key figures / tables
- **Fig. 1 (A–D):** Calculated concentration profiles and rotor-speed insets for conventional SE vs. three toSE conditions; Teq = 28 h / 4.5 h / 9 h / 13 h as the base over-concentration cap tightens (5× → <1.07).
- **Fig. 2:** Deq(t) RMSD-to-Boltzmann time courses for the Fig. 1 schedules, including a penalized "suboptimal" trace with a deceptive early minimum.
- **Fig. 3:** Experimental BSA toSE — absorbance scans + c(s) fit (inset) and approach-to-equilibrium by average buoyant molar mass (toSE ~7 h vs. constant-speed ~20 h).
- **Fig. 4 (A–C):** SBTI/CT interacting system, sequential multi-speed equilibria, showing toSE hastens equilibration at low speeds.

## 9. Notable quotes / citable claims
- "We have developed a method for time-optimized SE (toSE) with defined time-varying centrifugal fields that allow SE to be attained in a significantly (up to 10-fold) shorter time than is usually required." (Abstract)
- "As a finite-element algorithm, it is not subject to the numerical diffusion errors that occurred in simpler differentiation schemes previously applied to the overspeeding problem." (Methods, p. 828)
- "the most critical parameter in the toSE optimization is the maximally allowed overconcentration … even a moderate overconcentration by 20–50% over the final equilibrium concentration could allow time savings by a factor of 3–5." (Discussion, p. 835)
- "we have described the use of optimized time-varying centrifugal fields to significantly expedite SE experiments by at least a factor of 2 and up to a factor of 10." (Discussion, p. 836)

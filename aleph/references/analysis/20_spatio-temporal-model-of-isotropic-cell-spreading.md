---
id: 20_spatio-temporal-model-of-isotropic-cell-spreading
paper_n: 20
title: "A Three-Dimensional Stochastic Spatio-Temporal Model of Isotropic Cell Spreading"
authors: "Xiong, Rangamani, Dubin-Thaler, Sheetz & Iyengar (Xiong et al.)"
year: "2007"
venue: "Nature Precedings (preprint)"
doi: "10.1038/npre.2007.62.2"
paper_type: agent-based
ffn_relevance: High
ffn_themes: [cortex, FA/clutch, cell-spreading, parameter-source, validation-oracle, numerics/methods]
entities: [actin, arp2-3, capping-protein, filament, lamellipodium, cortex, fibronectin, integrin, plasma-membrane, fibroblast, cytochalasin-d]
methods: [agent-based, monte-carlo, gillespie-ssa, tirf-microscopy, brownian-ratchet]
measurables: [spreading-velocity, retrograde-flow, branching-rate, capping-rate, polymerization-rate, membrane-resistance-force, branching-angle]
keywords: [isotropic-cell-spreading, actin-branching, capping, arp2/3, membrane-resistance-force, mogilner-oster, brownian-ratchet, gillespie, leading-edge, force-dependent-kinetics, cytochalasin-d, fibroblast]
tags: ["#cell-spreading", "#cortex", "#actin-dynamics", "#arp2-3", "#force-dependent-kinetics", "#parameter-source", "#validation-oracle", "#brownian-ratchet", "#gillespie"]
has_transferable_params: true
---

# [20] A Three-Dimensional Stochastic Spatio-Temporal Model of Isotropic Cell Spreading

**Tags:** #cell-spreading #cortex #actin-dynamics #arp2-3 #force-dependent-kinetics #parameter-source #validation-oracle #brownian-ratchet #gillespie

| Field | Value |
|---|---|
| Authors | Xiong, Rangamani, Dubin-Thaler, Sheetz & Iyengar (et al.) |
| Year / Venue | 2007 / Nature Precedings (preprint) |
| DOI / ID | 10.1038/npre.2007.62.2 |
| Type | agent-based / stochastic spatio-temporal Monte Carlo (modified Gillespie SSA) |
| Pages | 44 |
| ffn_cellsim relevance | High — fine-grained 3D stochastic cortical actin (polymerize/branch/cap) under force-dependent kinetics; close structural cousin of ffn_cellsim H.1 cortex + H.5/H.7 spreading, and a rich parameter table |

## 1. Summary
The authors build a three-dimensional, stochastic, spatio-temporal model of *isotropic* spreading of mammalian fibroblasts on fibronectin-coated glass. The cortex is a branched actin filament network living in a 50-nm shell underneath a triangulated cell membrane (initial sphere 2 µm diameter, only the bottom 0.2-µm leading-edge sheet spreads). Three "minimal" reactions — actin polymerization, Arp2/3 branching, and capping — fire stochastically in space and time via a modified Gillespie algorithm with a dynamic dependency graph; each reaction's effective on-rate is exponentially damped by the local plasma-membrane resistance force (a Brownian-ratchet / Mogilner-Oster effective-rate law, Eq 1). Filament growth pushes the triangulated surface outward, defining the leading-edge spreading velocity per 10° angular segment. The model qualitatively reproduces the TIRF-measured isotropic spreading phenotype (similar velocity around the whole periphery after ~20-40 s), and maps a phase diagram in the (Arp2/3 concentration, membrane force) plane separating isotropic from non-isotropic spreading. A model *prediction* — that increasing capping-protein concentration linearly lowers mean peripheral velocity — was confirmed experimentally by titrating cytochalasin D (both regression slopes = 0.82; R² 0.98 sim vs 0.85 exp). The take-home: intrinsically stochastic cortical reactions yield deterministic-looking isotropic spreading once constrained by membrane force.

## 2. Problem & motivation
Cell motility/spreading is central to wound healing, immune response, development. Prior actin models were mostly 1D, steady-state, deterministic, "abstracted cytoskeletal structure." TIRF imaging (Dubin-Thaler 2004) now gives quantitative spatial velocity distributions at the leading edge, enabling a mechanistic question: *are the three core actin reactions (growth, branching, capping), when regulated by membrane force, sufficient to produce the observed isotropic spreading?* The paper tests this hypothesis with an explicit 3D stochastic filament-network simulation rather than a lumped PDE.

## 3. Methods / model
- **Model class:** explicit agent/event-based stochastic simulation. Each filament is an object; reactions are discrete events. Not a continuum PDE — this is a particle/event-level cortical actin model (conceptually adjacent to AFINES-style explicit-filament dynamics, but reaction-kinetic rather than Langevin-mechanical).
- **Four modules:** (i) **Filament network** — initialized as ~2000 spatially homogeneous seed filaments, each one Arp2/3 + two actin monomers, connected by branches off Arp2/3 sites; (ii) **Cell surface** — triangulated sphere (adjacent triangular polyhedra) that is actively updated so filament growth moves the leading edge; (iii) **Stochastic reaction machinery** — modified Gillespie's algorithm (ref 21, Gillespie 1977) with a *dynamic dependency graph* + sorted waiting-time list so only affected reactions are recomputed per event (efficient for a large growing reaction set); (iv) **Force-dependent growth** — membrane resistance force exponentially rescales the on-rates.
- **Governing equation (Eq 1, effective on-rate, Brownian ratchet / Mogilner-Oster):**
  k'_on = k_on · exp( −(f·δ) / (k_B T) ),
  where δ = filament length increment per added monomer, f = local membrane resistance force, k_B T at 300 K = 4.1 pN·nm. Resistance during polymerization/capping is imposed on the pre-existing filament; during branching on the newly formed filament.
- **Length/time scales:** cortical shell 50 nm thick; leading-edge sheet 200 nm; cell 2 µm diameter; monomer spacing 5.5 nm; simulated time ≤ ~1 minute (could not run longer); velocities binned per 10° angular segment, sampled at 5/10/20/40/60 s.
- **Numerics/solver:** discrete-event Monte Carlo (Gillespie SSA variant), C++/object-oriented on Linux; run on a 16-node Microway Beowulf (Intel Xeon) and the SDSC IA-64 TeraGrid for parameter sweeps. Termination on max time, max iterations, or full capping.
- **Experimental side:** mouse fibroblasts spreading on fibronectin-coated slides; spreading velocity by TIRF at 2-s intervals (assay from Dubin-Thaler 2004, ref 1). Cytochalasin D added 30 min before spreading at graded concentrations to mimic capping/barbed-end blockade.
- **Isotropy definition (gate):** spreading is "isotropic" when the standard deviation of peripheral velocity is <20% of the mean (derived from the experimental mean 3.4 ± 0.6 µm/min).

## 4. Key results (quantitative)
- Reference simulation (actin 20 µM, Arp2/3 0.04 µM, capping 0.04 µM, force 500 pN/µm²): peripheral spreading velocity ranges +12.8 to −4 µm/min (negative = retraction), averaged over 24 runs (Fig 2A); isotropy emerges by 20-40 s.
- Experimental mean peripheral velocity during isotropic stage: **3.4 ± 0.6 µm/min** (ref 1) → the 20% isotropy threshold (Results p.6).
- Simulated velocities are systematically *lower* than experiment (authors did not attempt quantitative velocity match because spatial concentrations unknown) (p.7).
- **Arp2/3 sweep** (0.02 / 0.08 / 0.16 µM, force 500 pN/µm²): velocity *range* roughly unchanged, but distribution becomes more spatially homogeneous (more isotropic) as Arp2/3 rises; velocity gain saturates with Arp2/3 — counter-intuitive plateau (Fig 3, p.7-8).
- **Membrane-force sweep** (0 → 1000 pN/µm², Arp2/3 0.08 µM): with zero force, mostly near-zero velocity with sporadic large spikes → force is *necessary* for spreading; raising force toward ~300-500 pN/µm² yields isotropic distributions (Fig 4A, Table 1).
- **Phase diagram (Fig 4B / Table 1):** low force (0/30/100 pN/µm²) → isotropic only at low Arp2/3 (0.02 µM); high force (300/500/1000 pN/µm²) → isotropic at higher Arp2/3 (0.04-0.16 µM). At 500 pN/µm², all Arp2/3 values 0.02-0.16 µM were isotropic (mean v ≈ 0.35-0.39 µm/min, %dev 9-17%).
- **Capping/cytochalasin-D prediction:** mean peripheral velocity decreases *linearly* with capping-protein (sim) and with cytochalasin-D (exp) concentration; both fitted slopes = **0.82** (text says −0.82); R² = 0.98 (sim) vs 0.85 (exp) (Fig 5, p.11-12).
- Table 1 full grid: mean velocity 0.011-2.03 µm/min, %deviation 9.43-76.91% across the (force, Arp2/3) conditions.

## 5. Parameters & constants of interest
(All from Supplementary Table S1 unless noted; references are the paper's own internal cites.)

| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Filament polymerization rate constant (k_on) | 11.6 µM⁻¹ s⁻¹ | Table S1 #1 (Pollard 2000) |
| Actin monomer spacing (length increment δ) | 5.5 nm | Table S1 #2 (Howard 2001) |
| Arp2/3 complex diameter | 15 nm | Table S1 #3 (Robinson 2001) |
| Capping protein diameter | 10 nm (estimated) | Table S1 #4 |
| Filament branching angle | 70° | Table S1 #5 (Medalia 2002) |
| Nucleated monomers per new branch | 2 | Table S1 #6 (Higgs & Pollard 2001) |
| Actin monomers per Arp2/3 binding footprint | 7 (≈6-7) | Table S1 #7 (Volkmann 2001) |
| Membrane resistance pressure | 500 pN/µm² (swept 0-1000) | Table S1 #8 (Mogilner & Edelstein-Keshet 2002) |
| Leading-edge sheet thickness | 200 nm | Table S1 #9 (Abraham 1999) |
| Initial cell diameter | 2 µm | Table S1 #10 |
| Cortical shell thickness (reactions confined here) | 50 nm | Table S1 #11 |
| Initial number of actin filaments | 2000 | Table S1 #12 |
| Thermal energy k_B T (at 300 K) | 4.1 pN·nm | Table S1 #13 |
| Filament branching rate constant | 1.25 µM⁻³ s⁻¹ (4th-order rxn) | Table S1 #14 (estimated) |
| Filament capping rate constant | 35 µM⁻¹ s⁻¹ | Table S1 #15 (estimated) |
| Arp2/3 initial concentration | 0.04 µM (swept 0.02/0.08/0.16) | Table S1 #16 |
| Capping protein initial concentration | 0.04 µM (swept 0.32/2.56) | Table S1 #17 |
| Actin monomer concentration (leading edge) | 20 µM (est. range 8-15 raw) | Table S1 #18 (Abraham 1999) |
| Measured lamellipodial polymerization | 97 ± 16 monomers/filament/s | Table S1 #18 (Abraham 1999) |
| Growing filaments at leading edge | 1370 ± 578 | Table S1 #18 (Abraham 1999) |
| Experimental mean spreading velocity (isotropic) | 3.4 ± 0.6 µm/min | Results p.6 (Dubin-Thaler 2004) |
| Effective on-rate law | k'_on = k_on·exp(−fδ/k_BT) | Eq 1, p.16 |

## 6. Relevance to ffn_cellsim
This is one of the more directly on-target references for ffn_cellsim's cortex/spreading chain. It is simultaneously a **(b) parameter source**, a **(c) mechanism reference**, a **(a) validation oracle**, and a **(e) numerics/methods** datapoint:

- **Cortex (H.1) + cell-spreading (H.5/H.7):** the model is exactly the regime ffn_cellsim targets — a branched actin cortex in a thin shell driving leading-edge protrusion. Its geometric/biochemical constants (branch angle 70°, Arp2/3 footprint = 7 monomers, 2 nucleated monomers/branch, monomer spacing 5.5 nm, cortical-shell thickness 50 nm, ~2000 seed filaments, 20 µM actin, 0.04 µM Arp2/3, 0.04 µM capping) are a ready-made initialization/parameter anchor set for the mechanistic Arp2/3 angle-harmonic branching, capping, and polymerization machinery. The 70° branch angle here maps onto ffn_cellsim's *angle-harmonic branch with thermal fluctuation* (the project replaces the rigid 72° with a fluctuating harmonic; 70° is a literature anchor for the equilibrium angle).
- **Mechanism contrast (architectural principle):** the paper uses the **Mogilner-Oster Brownian-ratchet effective-rate** (Eq 1, k'_on = k_on·exp(−fδ/k_BT)) as its *runtime mechanism*. In ffn_cellsim this exact closed form is **not** the runtime law — load-dependent monomer addition should *emerge* from explicit excluded-volume contact between filament tip and membrane + thermal motion under the BAOAB integrator. Therefore Eq 1 is a clean **acceptance oracle**: a fine-grained ffn_cellsim cortex pushing on a membrane should reproduce an exp(−fδ/k_BT) load-velocity relation, and Eq 1 with δ=5.5 nm, k_BT=4.1 pN·nm is the closed form to check against.
- **Validation oracle (spreading-velocity gate):** the experimental isotropy contract — mean peripheral velocity 3.4 ± 0.6 µm/min, isotropy = SD < 20% of mean — is a concrete, pre-registered validation gate ffn_cellsim could adopt for a single-cell spreading run, plus the (force, Arp2/3) phase diagram as a qualitative regime check.
- **FA/clutch context (weak):** spreading is on fibronectin and the discussion attributes Arp2/3 + capping + force regulation to integrin-mediated signaling (Fig 6), but adhesion/clutch is *not* modeled mechanistically here — it appears only as an upstream signaling box. So this informs the "why spreading happens on FN" context but is not a clutch parameter source.
- **Numerics/methods note:** they use a modified Gillespie SSA with a dynamic dependency graph — a reaction-kinetic event scheme, *different* from ffn_cellsim's Langevin/BAOAB particle dynamics. Useful as a comparison point (their inability to exceed ~1 minute of sim time, ~2000 filaments) but not a numerical method to import.

Not off-topic; the main caveat is that it is reaction-kinetic + force-as-an-effective-rate rather than fully mechanical, so ffn_cellsim treats its force law and rate constants as oracle/anchors rather than runtime mechanisms.

## 7. Limitations & caveats
- Force enters only as an exponential prefactor on rates (lumped Brownian ratchet); no explicit mechanical filament–membrane contact, no excluded volume, no retrograde flow physics — exactly the abstraction ffn_cellsim is designed to replace.
- Spatial concentrations of actin/Arp2/3/capping are unknown, so the authors explicitly do **not** match absolute velocities (sim velocities are several-fold lower than the 3.4 µm/min experiment).
- Several key rate constants are *estimated* to make the reaction-rate balance "moderately high," not measured (branching 1.25 µM⁻³ s⁻¹, capping 35 µM⁻¹ s⁻¹, capping/Arp2/3 diameters) — flag against the ffn_cellsim no-magic-number rule if borrowed; only the literature-cited values (poly rate 11.6 µM⁻¹ s⁻¹, branch angle 70°, δ 5.5 nm, footprints) are firm anchors.
- Only the bottom 0.2-µm leading edge spreads; the rest of the cell is frozen — a strong geometric simplification vs a full 3D cell.
- Simulations limited to ≤ ~1 minute and ~2000 filaments; preprint (Nature Precedings, not peer-reviewed).
- No myosin/contractility, no nucleus, no membrane tension dynamics beyond the constant resistance pressure.

## 8. Key figures / tables
- **Table 1 (p.26):** (force × Arp2/3 × capping) grid with mean velocity, SD, %deviation, isotropic yes/no — the quantitative phase-diagram data and the isotropy gate.
- **Table S1 (Suppl.):** full input-parameter list with values, units, and references — the transferable-constant goldmine.
- **Fig 4B (phase diagram):** isotropic vs non-isotropic regions in the (Arp2/3 conc, membrane force) plane.
- **Fig 5 (A-C):** capping/cytochalasin-D linear velocity decrease, sim vs experiment, slope 0.82 — the validated model prediction.
- **Eq 1 (p.16):** the force-dependent effective on-rate (Brownian-ratchet oracle).

## 9. Notable quotes / citable claims
- "The effective rate constant is given by k'_on = k_on · exp(−(fδ)/(k_BT)) … where δ is the length increment of actin filament due to addition of one monomer … and f is the resistance force." (p.16, Eq 1)
- "Experimentally, isotropic cell spreading is characterized by a mean peripheral velocity of 3.4 ± 0.6 µm/min … a deviation of 20% from the mean is an acceptable value for the definition of isotropic spreading behavior." (Results, p.6)
- "When the applied resistance force is zero … force is a necessary component for cell spreading and the biochemical kinetics alone are not sufficient to capture the spreading behavior observed experimentally." (Results, p.9)
- "A spatio-temporally complex model made up of a simple set of stochastic reactions near the cell surface, when constrained by membrane forces, can yield deterministic behavior as characterized by isotropic cell spreading." (Abstract)
- "Both lines showed a slope of 0.82, although the fit for the simulations (R²~0.98) was better than … the experiments (R²~0.85)." (Results, p.11-12)

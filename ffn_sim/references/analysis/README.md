# ffn_cellsim Reference Library — Master Index

## Scope & provenance

This is the master index for the `ffn_cellsim` reference library analysis. `ffn_cellsim` is a fine-grained, **mechanistic** HOOMD-blue simulator of single-cell mechanobiology: every cytoskeletal filament, motor head, adhesion clutch, and ECM cross-link is an explicit particle/bond. Published closed-form models (Chan-Odde, Pereverzev, Bell-Evans, Hill, Stam-Hocky) are **never** the runtime mechanism here — they are **acceptance oracles** used only in validation. This library is therefore mined for two things: (a) transferable physical constants/measurements to anchor the runtime, and (b) macroscopic validation oracles to gate it against.

**Provenance.** The user placed **31 items** in `ffn_sim/references/`, then added **4 more** (35 folder items total). After (a) extracting the bundled zip `Cell Press_20260602.zip` (19 PDFs) and (b) content-based de-duplication of identical/re-saved copies, there are **49 UNIQUE papers**.
- **Removed duplicate copies (4):** a second copy of the S0928098725003677 nutrient-limited breast-cancer-spheroid paper, plus three additional copies of the "Directed invasion of cancer cell spheroids" paper (also distributed as `journal.pone.0264571` and `journal.pone.0264571 (1)`).
- **3 genuinely new papers (second batch):** "The impact of 3D tumor spheroid maturity on cell migration and invasion dynamics" (#3), "Modeling Tumor Microenvironment Complexity In Vitro" review (#16), and "A cellular automaton model for spheroid response to radiation and hyperthermia" (#27).
- **Flat folder (papers 1–30):** cancer-spheroid / cell-spreading / mechanobiology.
- **Zip bundle (papers 31–49):** mixed Cell Press / Biophysical Journal set — some directly relevant (myosin heads, cortical actin, neutrophil spreading forces, nuclear shaping, lipid-bilayer microrheology), some membrane-biophysics relevant to the membrane/nucleus EXTEND track, and several tangential (prion SAXS, CaMKII arrhythmia, analytical ultracentrifugation, polarizable-FF protein repacking).

> Note: All **49 unique papers** now have a structured per-paper report — papers **#2** and **#9** were regenerated on a second pass after an initial tool-call miss, so the per-paper analysis now covers the full **49 reports** (no missing slots). The tables below reflect all 49 reports.

**Frontmatter.** Each per-paper report carries a YAML frontmatter block (`id / title / authors / year / venue / doi / paper_type / ffn_relevance / ffn_themes / entities / methods / measurables / keywords / tags`) for RAG/tag ingestion.

## How to read this library

Each paper has a per-paper report at `<slug>.md` in this directory (linked in the master table below). The **YAML frontmatter at the top of each report is the machine-ingestable record** — point your RAG/tagging pipeline at those blocks; the prose body is the human-readable analysis.

## Master table

| # | Title (short) | Type | Relevance | One-line takeaway | Report |
|---|---|---|---|---|---|
| 1 | Computational model of early cell spreading, migration & taxis | continuum-model | Medium | 1D active-gel + RD continuum oracle for spreading/migration/durotaxis; param anchor, not runtime. | [01_…spreading-migration-and.md](01_a-computational-model-for-early-cell-spreading-migration-and.md) |
| 2 | Integrating simulation & experimental validation of nutrient-limited BC spheroid growth | continuum-model | Low | COMSOL FEM continuum digital-twin of BT-474 spheroids (Gompertz + porous nutrient diffusion + hyperelastic mesh); tissue-scale, lumped, orthogonal to fine-grained runtime. | [02_…validation-of-nutrie.md](02_integrating-simulation-and-experimental-validation-of-nutrie.md) |
| 3 | Impact of 3D tumor spheroid maturity on migration/invasion | experimental | Low | Mature spheroids invade farther; phenomenological, no transferable single-cell mechanism. | [03_…maturity…migration.md](03_the-impact-of-3d-tumor-spheroid-maturity-on-cell-migration-a.md) |
| 4 | Multiscale model for heterogeneous tumor spheroid | agent-based | Low | Each cell = one viscoelastic ellipsoid + O2/glucose RD; lumped, opposite granularity. | [04_…mbe-2018016.md](04_10-3934-mbe-2018016.md) |
| 5 | 3D ADSC/breast-cancer co-culture spheroid on collagen-I | experimental | Low | Same spheroid-on-collagen genre; no mechanism/equation/constant. | [05_…spheroid-coculture.md](05_11744268-6759-4d1c-8363-9666c54cd8c3.md) |
| 6 | Math model of tumour spheroid w/ real-time cell-cycle imaging | continuum-model | Low | FUCCI cell-cycle PDE (necrotic core/G1 shell/rim); tissue-scale, no mechanics. | [06_…real-time.md](06_mathematical-model-of-tumour-spheroid-experiments-with-real-.md) |
| 7 | 3D spheroid preserves TME of hot/cold breast-cancer subtypes | experimental | Low | MCF-7 vs MDA-MB-231 immuno-oncology; shares cell lines, no mechanics. | [07_…tumor-microenvironment.md](07_in-vitro-3d-spheroid-model-preserves-tumor-microenvironment-.md) |
| 8 | Biofabrication of spheroid-fusion tumor models (glucose) | agent-based | Low | Cellular-Potts spheroid fusion; dimensionless a.u. knobs only. | [08_…fusion-based.md](08_biofabrication-of-spheroids-fusion-based-tumor-models-comput.md) |
| 9 | Directed invasion of cancer cell spheroids in flow-oriented 3D collagen | experimental | Low | Shear-aligned collagen-I shows invasion ~2–3× faster along radial than tangential fibers; minimal 2D anisotropic random walk fits with zero perpendicular step (migration along fibers, sqrt-time growth). | [09_…3d-collage.md](09_directed-invasion-of-cancer-cell-spheroids-inside-3d-collage.md) |
| 10 | Mechanics of cells in traction-driven cell spreading | continuum-model | Medium | Triangulated whole-cell (cortex+nucleus+AF/IF/MT); reusable cortex/nucleus param table + E~1–5 kPa oracle. | [10_…physreve-93-042404.md](10_physreve-93-042404.md) |
| 11 | Spheroid mechanics & implications for cell invasion | review | Medium | Soft-matter review: TST/MPA/AFM tools, surface/cortical tension, DAH/jamming. | [11_…cell-invasion.md](11_spheroid-mechanics-and-implications-for-cell-invasion.md) |
| 12 | Spheroid formation in microwells: MC simulations | agent-based | Low | Lattice kMC of hepatocarcinoma aggregation; dimensionless, no single-cell mechanism. | [12_…microwells.md](12_spheroid-formation-of-hepatocarcinoma-cells-in-microwells-ex.md) |
| 13 | FA dynamics + cytoskeleton + motors on 3D curved ECM | continuum-model | Medium | Node-mesh MC integrin-clutch + stress fibers + nucleus; rich FA/clutch params (R²=0.771). | [13_…rsc-ib.md](13_rsc-ib-c2ib20159c-1-12.md) |
| 14 | 3D tumor spheroid models for anticancer therapy (review) | review | Low | Pharmacology review; no mechanics/equations/params. | [14_…three-dimensional.md](14_three-dimensional-in-vitro-tumor-spheroid-models-for-evaluat.md) |
| 15 | Spheroid 3D models to decode matrix effectors in breast cancer | experimental | Low | MCF-7/MDA-MB-231 EMT/MMP/syndecan profiling; same lines, no mechanistic params. | [15_…matrix-effectors.md](15_spheroid-based-3d-models-to-decode-cell-function-and-matrix-.md) |
| 16 | Modeling TME complexity in vitro (review) | review | Low | TME-barrier drug-delivery review; tissue-scale phenomenology only. | [16_…microenvironment-complexity.md](16_modeling-tumor-microenvironment-complexity-in-vitro-spheroid.md) |
| 17 | Hybrid model: spheroid growth + ribose collagen stiffening | agent-based | Medium | PhysiCell MCF7/HCC1954 + ECM; spheroid-scale validation companion, lumped cells. | [17_…fbioe-13-1515962.md](17_fbioe-13-1515962.md) |
| 18 | Cancer-cell spheroid invasion assay (3D) | methods-software | Low | JoVE hanging-drop invasion protocol; ordinal readouts, no constants. | [18_…jove-105-53409.md](18_jove-105-53409.md) |
| 19 | Multiscale spheroid tumors: nutrient effect on evolution | agent-based | Low | 3D Cellular-Potts; mechanical constraint selects low-adhesion; CPM units only. | [19_…nutrient-av.md](19_multiscale-modeling-of-spheroid-tumors-effect-of-nutrient-av.md) |
| 20 | 3D stochastic model of isotropic cell spreading | agent-based | **High** | Gillespie-SSA cortical actin (poly/branch/cap) w/ Mogilner-Oster law; rich cortex/Arp2-3 params + spreading-velocity gate. | [20_…isotropic-cell-spreading.md](20_spatio-temporal-model-of-isotropic-cell-spreading.md) |
| 21 | 3D spheroid & organoid TME for immunotherapy (review) | review | Low | Immunology TME review; coarse boundary-condition context only. | [21_…organoid.md](21_3d-tumor-spheroid-and-organoid-to-model-tumor-microenvironme.md) |
| 22 | Initial cell spreading w/ mechanistic contact formulations | agent-based | Medium | Maugis-Dugdale deformable cell; r∝t^½ spreading law as analytic oracle + cortex params. | [22_…pcbi-1003267.md](22_pcbi-1003267-1-14.md) |
| 23 | Scoring function comparing simulated vs experimental spheroids | methods-software | Low | Point-cloud Wasserstein deviation score; sim-vs-exp comparison tooling, CPM internals. | [23_…scoring-function.md](23_development-of-a-scoring-function-for-comparing-simulated-an.md) |
| 24 | FA dynamics via photonic resonator outcoupler microscopy | methods-software | Low | Label-free FA imaging optics; FA-size bands only, no force model. | [24_…photonic.md](24_quantitative-analysis-of-focal-adhesion-dynamics-using-photo.md) |
| 25 | Simulation-based inference of cell migration in complex space | agent-based | Low | CPM dendritic-cell chemotaxis + NPE/ABC; methods, resolution-locked params. | [25_…sim-based-inference.md](25_simulation-based-inference-of-cell-migration-dynamics-in-com.md) |
| 26 | TASI: spatial-temporal quantification of spheroid dynamics | methods-software | Low | MATLAB spheroid image-analysis pipeline; no mechanics/params. | [26_…tasi.md](26_tasi-a-software-tool-for-spatial-temporal-quantification-of-.md) |
| 27 | Cellular automaton: spheroid response to radiation/hyperthermia | agent-based | Low | Voxel CA + O2 RD + clonogenic survival; population oncology, no cytoskeleton. | [27_…radiation.md](27_a-cellular-automaton-model-for-spheroid-response-to-radiatio.md) |
| 28 | Tumor spheroids accelerate persistently invading cells | experimental | Low | Spheroid ~1/r² repulsive field speeds escaping cells; lumped COM scale, no mechanism. | [28_…persistently-invading.md](28_tumor-spheroids-accelerate-persistently-invading-cancer-cell.md) |
| 29 | FA orientation dynamics under static/cyclic stretch | continuum-model | Medium | Gillespie FA bond cluster w/ Pereverzev catch-slip; catch-bond oracle + params. | [29_…fa-orientation.md](29_a-general-model-of-focal-adhesion-orientation-dynamics-in-re.md) |
| 30 | Hydrogel microenvironments for spheroid growth & drug screening | review | Low | Hydrogel scaffold catalog; bulk gel moduli 90–6000 Pa, context only. | [30_…hydrogel.md](30_hydrogel-microenvironments-for-cancer-spheroid-growth-and-dr.md) |
| 31 | Prion protein–antibody complexes by SAXS | experimental | Low | Structural biology of PrP/Fab; no cell mechanics. | [31_…prion.md](31_prion-protein-antibody-complexes-characterized-by-chromatogr.md) |
| 32 | Protein conformational changes resolved by SHG | experimental | Low | In-vitro SHG on supported bilayer; outside cytoskeleton domain. | [32_…shg.md](32_protein-conformational-changes-are-detected-and-resolved-sit.md) |
| 33 | Dead-end elimination repacks PCNA (polarizable FF) | methods-software | Low | Atomistic crystallographic repacking (AMOEBA); off-topic. | [33_…dead-end.md](33_dead-end-elimination-with-a-polarizable-force-field-repacks-.md) |
| 34 | Conformation of myosin heads in relaxed skeletal muscle | experimental | Low | RLC fluorescence-polarization, ~70% IHM OFF heads; sarcomeric, no single-motor constants. | [34_…myosin-heads.md](34_the-conformation-of-myosin-heads-in-relaxed-skeletal-muscle-.md) |
| 35 | Hydrophobic mismatch: E5 ↔ PDGFR TM segment | experimental | Low | Membrane-protein NMR/CD; helix-tilt only, no cytoskeleton. | [35_…hydrophobic-mismatch.md](35_hydrophobic-mismatch-drives-the-interaction-of-e5-with-the-t.md) |
| 36 | Micelle/fusion thermodynamics of antimicrobial lipopeptides | particle/MD-sim | Low | MARTINI MD free-energy of AMLP insertion; drug-membrane, off-topic. | [36_…micelle.md](36_thermodynamics-of-micelle-formation-and-membrane-fusion-modu.md) |
| 37 | Synaptobrevin TM domain dimerization (multiscale MD) | particle/MD-sim | Low | MARTINI+AA MD of SNARE TM dimer; membrane-protein, off-topic. | [37_…synaptobrevin.md](37_synaptobrevin-transmembrane-domain-dimerization-studied-by-m.md) |
| 38 | NMR dynamics of p75NTR domains in nanodiscs | experimental | Low | Solution NMR of disordered receptor linker; no cell mechanics. | [38_…p75ntr.md](38_nmr-dynamics-of-transmembrane-and-intracellular-domains-of-p.md) |
| 39 | Protrusive & contractile forces of spreading neutrophils | experimental | Medium | Micropost traction: ~75 pN protrusive wave + ROCK/myosin contraction; spreading-force oracle + params. | [39_…neutrophils.md](39_protrusive-and-contractile-forces-of-spreading-human-neutrop.md) |
| 40 | SpIDA reveals abnormal protein oligomerization in single cells | methods-software | Low | Fluorescence histogram oligomerization method; no mechanics. | [40_…spida.md](40_spatial-intensity-distribution-analysis-reveals-abnormal-oli.md) |
| 41 | Anionic lipids modulate aquaglyceroporin GlpF | experimental | Low | Bacterial channel + lipid charge; membrane biophysics, off-topic. | [41_…glpf.md](41_anionic-lipids-modulate-the-activity-of-the-aquaglyceroporin.md) |
| 42 | Two-point microrheology of phase-separated lipid bilayers | experimental | Low | First two-point membrane viscosity (0.75–3.9 nPa·s·m); EXTEND-track membrane viscosity ref. | [42_…microrheology.md](42_two-point-microrheology-of-phase-separated-domains-in-lipid-.md) |
| 43 | Through Thick and Thin — interfilament comms in muscle | review | Low | Commentary on thick/thin filament crosstalk; sarcomeric, no transferable mechanism. | [43_…thick-and-thin.md](43_through-thick-and-thin-interfilament-communication-in-nbsp-m.md) |
| 44 | Better Together — lipopeptide micelle selectivity (commentary) | review | Low | 2-page commentary on AMLP micelle selectivity; off-topic. | [44_…better-together.md](44_better-together-lipopeptide-micelle-formation-enhances-antim.md) |
| 45 | Moving cell boundaries drive nuclear shaping in spreading | continuum-model | Medium | Centripetal network flow flattens nucleus; nucleus param table + nuclear-flattening oracle. | [45_…nuclear-shaping.md](45_moving-cell-boundaries-drive-nuclear-shaping-during-cell-spr.md) |
| 46 | Feeling for filaments: cortical actin web in live endothelium | experimental | Medium | AFM+Lifeact mesh metrics (~42,000 nm² hole, ~100–150 nm spacing); direct cortex mesh-size params. | [46_…feeling-for-filaments.md](46_feeling-for-filaments-quantification-of-the-cortical-actin-w.md) |
| 47 | SpIDA surveys macromolecular oligomerization in situ | review | Low | Commentary on fluorescence oligomerization methods; no mechanics. | [47_…spida-surveys.md](47_spida-surveys-the-intricate-web-of-macromolecular-oligomeriz.md) |
| 48 | Variable-field analytical ultracentrifugation | methods-software | Low | AUC sedimentation-equilibrium methods; solution biophysics, off-topic. | [48_…ultracentrifugation.md](48_variable-field-analytical-ultracentrifugation-i-time-optimiz.md) |
| 49 | Arrhythmogenic role of oxidized CaMKII in the heart | continuum-model | Low | Cardiac electrophysiology model; off-topic. | [49_…camkii.md](49_mechanistic-investigation-of-the-arrhythmogenic-role-of-oxid.md) |

## Relevance tiers

### High (1)
- **#20** — *3D stochastic isotropic cell spreading*: directly mechanistic cortical-actin poly/branch/cap with a force-dependent (Mogilner-Oster/Brownian-ratchet) rate law; supplies a rich Arp2/3/cortex parameter set **and** a spreading-velocity validation gate (3.4±0.6 µm/min, isotropy SD<20%).

### Medium (8)
- **#1** — active-gel continuum oracle for retrograde flow / membrane tension / migration speed; parameter anchor.
- **#10** — triangulated whole-cell cortex+nucleus model: reusable cortex/nucleus parameter table + whole-cell E~1–5 kPa stiffening oracle.
- **#11** — spheroid-mechanics review: surface/cortical-tension and ECM-density invasion numbers usable as overlay context.
- **#13** — node-mesh FA-clutch + stress-fiber + nucleus model with concrete bond-stiffness / FA / motor parameters.
- **#17** — PhysiCell MCF7/HCC1954 spheroid-on-collagen: spheroid-scale validation companion for the multi-cell overlay.
- **#22** — Maugis-Dugdale deformable cell: analytic r∝t^½ spreading-law oracle + cortex elastic/bending/adhesion params.
- **#29** — Gillespie FA catch-slip (Pereverzev) cluster: catch-bond off-rate oracle + scaled rate parameters.
- **#39** — neutrophil micropost traction: protrusive/contractile force magnitudes + spreading-velocity oracle.
- **#45** — nuclear-shaping continuum model: nucleus modulus / lamina / network-viscosity parameter table + nuclear-flattening oracle.
- **#46** — live cortical-actin AFM mesh metrics: direct cortex mesh-size / filament-spacing parameters.

### Low (40)
Spheroid-scale oncology context (no transferable single-cell mechanism): **#2, #3, #4, #5, #6, #7, #8, #9, #12, #14, #15, #16, #19, #21, #23, #25, #26, #27, #28, #30** plus methods/imaging tooling **#18, #24**. Tangential structural/membrane/muscle/solution biophysics from the Cell Press bundle: **#31, #32, #33, #34, #35, #36, #37, #38, #40, #41, #42, #43, #44, #47, #48, #49**. (Kept because the user requested all bundle papers; several membrane ones — #42 especially — are weak EXTEND-track references. #2 and #9 are continuum/phenomenological spheroid-scale efforts: #2 a COMSOL nutrient-growth digital twin, #9 a flow-oriented-collagen invasion study — both context-only, no fine-grained mechanism.)

## Thematic map for ffn_cellsim

A paper may appear under several themes.

- **Cortex / actin** — #10, #11, #20, #22, #39, #46, #1, #45, #12
- **FA / clutch** — #13, #29, #39, #20, #22, #24, #10, #1, #45
- **Motor / myosin** — #1, #13, #39, #34 (muscle), #43 (muscle), #11
- **ECM / collagen** — #11, #13, #17, #29, #30, #23, #1, #3, #5, #9, #15, #16, #18
- **Spheroid-scale context (multi-cell overlay)** — #2, #3, #4, #5, #6, #7, #8, #9, #11, #12, #14, #15, #16, #17, #19, #21, #22, #23, #26, #27, #28, #30
- **Cell-spreading dynamics** — #1, #10, #20, #22, #39, #45, #13, #17
- **Membrane / nucleus EXTEND track** — #10 (nucleus/membrane), #45 (nucleus), #11 (nucleus/cytoplasm), #13 (nucleus), #42 (membrane viscosity), #32/#35/#36/#37/#38/#41/#44/#47 (membrane biophysics, weak)
- **Numerics / methods** — #1, #2 (COMSOL FEM/digital-twin), #9 (FEM flow + random walk), #13, #20, #22, #23, #25, #45, #46, #33, #36, #37, #48
- **Tangential / off-topic** — #2, #9, #31, #32, #33, #34, #35, #36, #37, #38, #40, #41, #43, #44, #47, #48, #49

## Parameter & oracle harvest

Most useful transferable constants/measurements flagged across the summaries (skipping papers with `key_params: none`).

| # | Quantity | Value | Use for ffn_cellsim |
|---|---|---|---|
| 1 | Retrograde flow | 0.045 µm/s (spreading), ~0.1 µm/s (mesenchymal) | retrograde-flow validation band |
| 1 | Actomyosin peak stress | ~30 Pa | contractile-stress oracle |
| 1 | Membrane tension | ~0.05 nN/µm | membrane-tension anchor |
| 1 | Protrusion force-velocity stall exponent | γ=8 (Keren 2008) | protrusion law calibration |
| 2 | Breast-tissue 1st Lamé λ | 4–50 kPa (Krouskop 1998 / Samani 2003) | bulk-tissue stiffness overlay anchor (EXTEND, overlay-only) |
| 2 | Breast-tissue shear modulus µ | 0.2–3.5 kPa (Krouskop 1998 / Samani 2003) | bulk-tissue stiffness overlay anchor (EXTEND, overlay-only) |
| 2 | Cell diameter | 15 µm | cytoplasm/medium transport context |
| 2 | Glucose / O₂ diffusion (DMEM) | 5.9e-10 / 1.5e-9 m²/s | medium transport context (if metabolism added) |
| 2 | O₂ internalization k_int | 2.5e-18 mol O₂·cell⁻¹·s⁻¹ (Wagner 2011) | metabolism context (out of scope) |
| 2 | BT-474 necrosis glucose threshold | ~0.08 mM (cf. Jiang 0.06 mM) | spheroid-scale necrosis overlay (line-specific, never fit-to) |
| 10 | Cortical CSK surface shear modulus | 6.3 µN/m (0.75–15) | cortex shear-modulus anchor |
| 10 | Nucleus surface shear modulus | 18.9 µN/m (3× cortex) | nucleus stiffness anchor |
| 10 | Cortex / nucleus bending | 2.77e-19 J (~67 kBT) / 5.54e-19 J | bending-rigidity anchor |
| 10 | WLC ratios L₀:L_C:L_p | 1:3:0.1, x₀=1/3 | cortex WLC bond setup |
| 10 | Whole-cell Young's modulus | 1.008–5.124 kPa (Sneddon) | whole-cell stiffening **oracle** |
| 10 | Microtubule E / cross-section | 1.2 GPa / 190 nm² | MT bond stiffness |
| 11 | Diffusion limit → necrotic core | ~200 µm | spheroid overlay geometry |
| 11 | Collagen-I invasion transition | 1 mg/ml (single-cell) vs 4 mg/ml (collective) | ECM-density overlay |
| 11 | In-vivo-mimic shear modulus | ~300 Pa (alginate/Matrigel) | ECM modulus context |
| 13 | Single ligand-receptor bond stiffness k_LR | ~1.0 pN/nm | clutch bond stiffness |
| 13 | Integrin equilibrium length | 30 nm | clutch geometry |
| 13 | FA-formation critical gap h_c | 300 nm | FA nucleation criterion |
| 13 | Stress-fiber Young modulus / radius | 230 kPa / ~250 nm | stress-fiber bond |
| 13 | Non-muscle myosin-II sliding rate | 10 nm/s | motor velocity anchor |
| 13 | F-actin poly / disassembly time | 180 s / ~1 s | filament turnover |
| 13 | Collagen-I ligand density | 750 molecules/µm² @ 0.8 mg/mL | ECM ligand density |
| 17 | Cell-ECM migration speed S₀ | MCF7 0.1, HCC1954 0.7 µm/min | multi-cell migration overlay |
| 17 | ECM degradation rate | MCF7 1e-4, HCC1954 3.2e-3 min⁻¹ | MMP degradation overlay |
| 20 | Brownian-ratchet on-rate | k'_on=k_on·exp(−fδ/kBT), δ=5.5 nm, kBT=4.1 pN·nm | **force-dependent poly law** |
| 20 | Arp2/3 branch angle / footprint | 70° / 7 monomers, 2 nucleated/branch | Arp2/3 branching geometry |
| 20 | Monomer spacing / cortical shell | 5.5 nm / 50 nm | cortex discretization |
| 20 | Poly / branch / cap rates | 11.6 / 1.25 (µM⁻³) / 35 µM⁻¹s⁻¹ | actin kinetics |
| 20 | Lamellipod poly (measured) | 97±16 monomers/filament/s, 1370±578 filaments | actin-density gate |
| 20 | Spreading-velocity gate | 3.4±0.6 µm/min, isotropy SD<20% | **spreading validation gate** |
| 22 | Cortex E / bending / stretch | ~800 kPa / 2.4e-19 N·m / FENE k_s~1 mN/m | cortex bond params |
| 22 | Adhesion energy W | 1–88 mJ/m², range h₀~20 nm | adhesion energy band |
| 22 | Spreading-law oracle | a ∝ sqrt(2W/c_n)·t^½ | **t^½ spreading oracle** |
| 29 | Catch-slip off-rate | k_off=k_slip·e^{+f/f0}+k_catch·e^{−f/f0} (Pereverzev) | **catch-bond off-rate oracle** |
| 29 | Scaled rates / cluster | K_c=120, K_s=0.10, Γ=2; N=200 bonds, f0~pN | FA cluster setup |
| 29 | Static instability threshold | ~20% strain | mechanosensing gate |
| 39 | Per-post protrusive force | 75±8 pN | protrusion force **oracle** |
| 39 | Contractile force (perimeter/core) | 106±10 / 20±10 pN/post | contraction force oracle |
| 39 | Spreading velocity | 206±28 nm/s (control), 61±37 (cytochalasin B) | spreading-velocity oracle |
| 39 | Micropost spring constant | 0.28±0.09 pN/nm | traction-sensor calibration |
| 39 | Contact growth law | R ~ t^0.4 | spreading exponent cross-check |
| 45 | Nucleus bulk / area modulus | K=250 Pa / k=25 mN/m | nucleus bond params (EXTEND) |
| 45 | Lamina bending stiffness | k_c=3.5e-4 nN·µm | nuclear lamina bond |
| 45 | Cytoplasmic network viscosity | µ=0.21 kPa·s | cytoplasm viscosity (EXTEND) |
| 45 | Contractile stress | σ_c~190 Pa | contractile-stress anchor |
| 45 | Steady nuclear aspect ratio | ~0.25 (height ~¼ width) at ~50% spread | **nuclear-flattening oracle** |
| 46 | Cortical mesh hole area / count | ~42,045 nm² (~0.042 µm²) / 271 per 6×6 µm | **cortex mesh-size params** |
| 46 | Smallest filament spacing | ~100–150 nm | cortex network spacing |
| 46 | Actin fraction of total protein | 5–15% (endothelial) | actin concentration context |
| 42 | Lipid-bilayer 2D viscosity | LD 0.75±0.15, LO 3.90±0.42 nPa·s·m | membrane viscosity (EXTEND) |
| 30 | Bulk gel moduli | fibrin 90–1050 Pa; collagen 300–6000 Pa; invasion optimum ~100 Pa, pore ~8 µm | ECM modulus overlay |
| 34 | Skeletal RLC IHM OFF fraction | ~70% near physiological T | myosin OFF-state context (muscle) |

> All spheroid-scale agent-based papers (#4, #8, #12, #19, #23, #25, #27) report **dimensionless Cellular-Potts / lattice / Monte-Carlo internal units**, not SI mechanism constants — context only, not runtime-transferable.
> **#9** (flow-oriented collagen invasion) reports only **tissue/ECM-scale phenomenological values** — collagen-I 1.85 mg/ml, invasion fronts (downstream 222±33 µm vs upstream 132±23 µm/day), best-fit anisotropic random-walk steps (d_parallel 220 µm/day, d_perp 0 µm/day), simulated rigid-fiber 47.5×5 µm — and **no transferable single-cell mechanistic constants** (no moduli, off-rates, motor forces, persistence lengths); `has_transferable_params: false`. The "complete blockage perpendicular to fibers / migration along fibers, sqrt-time growth" result is at best a high-level behavioral acceptance target for the multi-cell overlay, not an oracle.
> **#2** (BC-spheroid digital twin) is continuum/lumped: the only literature-anchored transferable values are bulk breast-tissue moduli (Lamé λ 4–50 kPa, shear µ 0.2–3.5 kPa) — overlay-only EXTEND anchors, not single-cell cytoskeletal constants — plus cell-diameter/diffusion context; its closed-forms (Gompertz, Kozeny-Carman, Millington-Quirk, Neo-Hookean) are tissue-growth/transport, not single-cell mechanism oracles.

## RAG tag index (inverted)

The user finds references hard to look up in the current RAG system. These inverted indexes let any tag point back to papers. Sorted by frequency (most common first).

### By `#tag`
- `#tangential` → #3, #4, #5, #6, #7, #8, #9, #12, #14, #16, #18, #19, #21, #23, #24, #25, #27, #28, #31, #32, #33(off-topic), #34, #35, #36, #37, #38, #40, #41, #42, #43, #44, #47, #48, #49, #2
- `#spheroid-mechanics` / `#spheroid-scale-context` → #3, #4, #6, #9, #11, #12, #14, #16, #17, #19, #21, #23, #27, #28, #30, #2
- `#breast-cancer` → #2, #3, #5, #7, #11, #15, #17, #18, #23, #30
- `#parameter-source` → #1, #10, #13, #20, #22, #29, #39, #45, #46, #2
- `#validation-oracle` → #17, #20, #23, #29, #39, #45 (+ #10, #11 implicit)
- `#cell-spreading` → #1, #10, #20, #22, #39
- `#review` → #11, #14, #16, #21, #30, #43, #44, #47
- `#cellular-potts` / `#agent-based` → #4, #8, #12, #17, #19, #22, #23, #25, #27
- `#cell-migration` → #1, #9, #13, #25, #28
- `#ecm-collagen` / `#ECM-stiffness` → #5, #9, #11, #16, #17, #23, #30
- `#membrane-biophysics` → #35, #36, #37, #41, #42, #44, #47 (+ #32 SLB)
- `#cortex` / `#cortex-structure` → #10, #20, #22, #46 (+ #11, #12, #39, #45)
- `#focal-adhesion` / `#fa-clutch` → #13, #24, #29 (+ #20, #39)
- `#contact-guidance` / `#fiber-alignment` / `#anisotropic-migration` → #9
- `#nutrient-diffusion` → #2 (+ #4, #6, #8, #19, #27 RD-context)
- `#continuum-model` → #1, #2, #6, #10, #13, #29, #45, #49
- `#catch-bond` → #29
- `#traction-force` → #11, #39 (+ #1, #10)
- `#cortical-tension` → #11, #39
- `#nucleus-mechanics` → #13, #45 (+ #10, #11)
- `#structural-biology` → #31, #32, #33, #35, #38, #41, #43
- `#coarse-grained-md` / `#martini` → #36, #37, #44
- `#motor-myosin` / `#myosin-structure` → #34, #43 (+ #1, #13, #39)
- `#methods-software` → #18, #23, #24, #26, #33, #40, #48
- `#monte-carlo` → #8, #12, #19, #22, #23, #25, #29
- `#tumor-microenvironment` → #7, #16, #21
- `#emt` → #3, #15
- `#mcf7` / `#mda-mb-231` → #2, #4, #7, #11, #15, #17, #23, #30
- `#drug-screening` → #14, #30
- `#interacting-heads-motif` → #34, #43
- `#scale-context` → #9, #2
- `#off-topic` / `#off-topic-for-cellsim` → #33, #48, #49

### By entity
- `actin` → #1, #10, #13, #20, #22, #24, #28, #29, #34, #39, #43, #45, #46
- `integrin` → #1, #3, #4(implied), #5, #7(implied), #8, #10, #11, #13, #16, #17, #20, #24, #29(implied), #30, #39, #45
- `collagen-i` → #3, #5, #9, #11, #13, #16, #17, #18, #23, #30
- `myosin-ii` → #1, #11, #13, #34, #39, #43, #45
- `e-cadherin` → #3, #5, #7, #11, #14, #15, #16, #21, #30
- `mcf7` → #2, #3, #4, #7, #8, #11, #14, #15, #17, #19(implied), #30
- `mda-mb-231` → #5, #7, #11, #15, #18, #23, #30
- `focal-adhesion` → #1, #3, #10, #13, #22, #24, #29
- `spheroid` (as explicit entity) → #2, #4, #5(implied), #6, #8, #9, #11, #12, #14, #18, #19, #23, #27
- `lipid-bilayer` → #10, #22, #32, #35, #36, #37, #38, #41, #42, #44(implied)
- `nucleus` → #1, #10, #11, #13, #24, #25, #45
- `arp2-3` → #1, #20, #39
- `vinculin` → #1, #10, #13, #24
- `vimentin` → #3, #5, #7, #11, #14, #15, #16, #45
- `fibronectin` → #5, #8, #13, #20, #24, #39, #45 (+ #16, #21)
- `matrigel` → #3, #11, #18, #21, #23, #28, #30
- `breast-cancer` → #2, #3, #5, #7, #11, #15, #17, #18, #23, #30
- `hela` → #9
- `cancer-cell-spheroid` → #9, #2 (+ spheroid entity papers above)
- `rat-tail-collagen` → #9
- `glucose` → #2 (+ #4, #6, #8, #19 RD-context)
- `oxygen` → #2 (+ #4, #6, #19, #27 RD-context)
- `bt-474` → #2
- `her2` → #2
- `regulatory-light-chain` → #34, #43, #45(lamin context n/a)
- `thick-filament` / `thin-filament` → #34, #43
- `microtubule` → #4, #10, #45
- `tumor-microenvironment` → #7, #14
- `linc-complex` / `nesprin` → #13 (nesprin-1), #45 (nesprin-2g/sun2/kash4)
- `cytochalasin-d` → #20, #46 (+ cytochalasin-b #39)
- `jasplakinolide` → #39, #46
- `blebbistatin` → #34, #39

### By method
- `confocal-live-imaging` → #2, #3, #5, #6, #7, #9, #12, #13, #14, #15, #16, #21, #24, #25, #26, #28, #30, #39, #40, #46, #47
- `monte-carlo` → #8, #12, #13, #19, #20, #22, #23, #25, #29, #31(conf. sampling)
- `agent-based` → #4, #8, #12, #13, #17, #19, #22, #23, #25, #27
- `finite-element` → #1, #2, #9, #13, #22, #24, #45, #48
- `reaction-diffusion` → #1, #2, #4, #6, #8, #17, #19, #27, #36
- `cellular-potts` → #8, #19, #23, #25
- `molecular-dynamics` → #22, #36, #37
- `atomic-force-microscopy` → #10, #11, #24, #39(implied micropost), #46
- `traction-force-microscopy` → #11, #30, #39
- `gillespie` / `gillespie-ssa` → #20, #29
- `random-walk` → #9
- `contact-guidance` → #9
- `microfluidics` → #9
- `particle-image-velocimetry` → #9
- `multiphysics-comsol` → #2, #9
- `gompertz-growth` → #2
- `neo-hookean-hyperelastic` → #2
- `darcy-flow` → #2
- `hplc` → #2
- `h-and-e-histology` → #2
- `immunostaining` → #2
- `master-equation` → #29
- `worm-like-chain` → #10
- `boundary-element-method` → #45
- `coarse-grained-md` / `martini` → #36, #37, #44
- `solution-nmr` / `nmr-relaxation` → #38 (+ solid-state NMR #35)
- `small-angle-x-ray-scattering` → #31
- `second-harmonic-generation` → #32
- `fluorescence-polarization` → #34, #43
- `fluorescence-correlation-spectroscopy` → #40, #47
- `dead-end-elimination` → #33
- `analytical-ultracentrifugation` → #48
- `stochastic-simulation` / `ode-electrophysiology` → #29, #49
- `single-particle-tracking` → #25, #42 (+ #28 single-cell-tracking)
- `light-microscopy` → #9
- `micropost-array` → #39

## Caveat

These are **single-pass direct PDF reads** — grounded in the source text, but not independently re-verified against external records. Numbers and claims should be spot-checked against the original PDF before being hard-wired as runtime constants or validation gates. The tangential Cell Press / Biophysical Journal bundle papers (prion SAXS, CaMKII arrhythmia, analytical ultracentrifugation, polarizable-FF repacking, membrane-protein NMR/MD, muscle myosin structure) were included **because the user asked for all of them**; they are catalogued for completeness, not because they feed the cytoskeleton/FA/ECM runtime.

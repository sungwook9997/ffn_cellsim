---
id: 39_protrusive-and-contractile-forces-of-spreading-human-neutrop
paper_n: 39
title: "Protrusive and Contractile Forces of Spreading Human Neutrophils"
authors: "Steven J. Henry et al."
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.05.041"
paper_type: experimental
ffn_relevance: Medium
ffn_themes: [cortex, FA/clutch, motor/myosin, cell-spreading, ECM/collagen, parameter-source, validation-oracle]
entities: [actin, myosin-ii, integrin, fibronectin, cortical-shell, arp2-3, rock, mac-1, beta2-integrin, lamellipodium, neutrophil]
methods: [traction-force-microscopy, micropost-array, spinning-disk-confocal, particle-tracking, pharmacological-inhibition]
measurables: [traction-stress, protrusive-force, contractile-force, spreading-velocity, cortical-tension, post-spring-constant, aspect-ratio, invagination-depth]
keywords: [neutrophil-spreading, micropost-array-detector, protrusion, contraction, cortical-tension, fibronectin, myosin-ii, arp2/3, cell-spreading-dynamics, traction-force, haptokinesis, cytochalasin-b]
tags: ["#cell-spreading", "#traction-force", "#cortical-tension", "#myosin-ii", "#validation-oracle", "#parameter-source", "#neutrophil", "#micropost-array"]
has_transferable_params: true
---

# [39] Protrusive and Contractile Forces of Spreading Human Neutrophils

**Tags:** #cell-spreading #traction-force #cortical-tension #myosin-ii #validation-oracle #parameter-source #neutrophil #micropost-array

| Field | Value |
|---|---|
| Authors | Steven J. Henry, Christopher S. Chen, John C. Crocker, Daniel A. Hammer |
| Year / Venue | 2015 / Biophysical Journal 109(4):699–709 |
| DOI / ID | 10.1016/j.bpj.2015.05.041 |
| Type | experimental (traction-force microscopy / micropost arrays + pharmacological dissection) |
| Pages | 11 |
| ffn_cellsim relevance | Medium — a quantitative per-post force/velocity dataset of a cell undergoing the spherical→spread shape transition, governed by adhesion-energy vs cortical-tension competition with explicit myosin-II / Arp2/3 / ROCK inhibitor controls; useful as a spreading-dynamics + cortex-relaxation validation oracle and as a force/velocity parameter anchor, though on neutrophils not MCF7. |

## 1. Summary
The authors measured the forces a quiescent human neutrophil exerts as it transitions from a non-adherent sphere to an adherent, spread "sessile drop" by spreading on fibronectin (FN)-printed elastomeric micropost arrays (mPADs) and tracking individual post-tip deflections at 1 frame/s. They identified two dominant mechanical regimes: (1) a transient, radially-propagating **protrusive wave** (~75 pN/post, FWHM ~17 s) moving outward at ~206 nm/s, and (2) a sustained **contractile** state in which perimeter posts (~106 pN/post) were ~5× more contractile than core posts (~20 pN/post). Pharmacological dissection showed the protrusive event is independent of Arp2/3 (CK666), myosin II (blebbistatin), and ROCK (Y27632), but its onset requires relaxation of the cortical actin shell: stiffening the cortex (jasplakinolide) abolished spreading, while softening it (cytochalasin B) slowed spreading to ~61 nm/s. Long-term sustained contractility was ROCK- and myosin-II-dependent. Spreading was haptokinetically induced by β2-integrin (MAC-1) ligation of FN. They frame the onset as a competition of adhesion energy, cortical tension, and relaxation of that cortical tension.

## 2. Problem & motivation
Neutrophils are first-responder innate-immune cells that undergo dramatic, fast shape changes (sphere → spread → migratory) on the timescale of seconds to minutes. Prior in vitro work (RICM, micropipette aspiration) characterized spreading dynamics and cortical mechanics but never directly measured the tractions of the spherical-to-spread transition. The goal was to quantify, with high spatial/temporal resolution, the protrusive and contractile forces of spreading and to identify which cytoskeletal mechanism (lamellipodial branching, actomyosin contraction, or cortical-tension relaxation) drives the onset of spreading.

## 3. Methods / model
- **Model class**: experimental traction-force microscopy. No computational model is built; the relevant theory invoked is the Cuvelier et al. "universal dynamics of cell spreading" (liquid drop + viscous shell) and Young's-equation droplet-wetting framing.
- **Assay**: FN-printed PDMS micropost-array detectors (mPADs); posts diameter 604±31 nm, length 5.576±0.286 µm; spring constant k_spring = 0.28±0.09 pN/nm (with 8% Schoen warping correction, smaller than measurement error). Cells gravity-sedimented onto FN-printed post tips.
- **Readout**: spinning-disk confocal (Olympus IX71, 60× water), 1 frame/s, post tips tracked in fluorescence channel via custom MATLAB particle-tracking (Crocker-Grier lineage); deflection × k_spring → per-post force trajectory in radial/tangential cell-frame coordinates. Force-detection noise 9±2 pN; S/N 8:1 (protrusion) and 12:1 (contraction).
- **Controls / perturbations**: 5 µM blebbistatin (myosin II), 1 µM CK666 (Arp2/3), 3 µM cytochalasin B (actin polymerization / filament-filament interaction), 1 µM jasplakinolide (F-actin stabilization/rigidification), 1 µM Y27632 (ROCK); anti-β2 (clone L130, 50 µg/mL) function-blocking antibody; Pluronic-blocked FN-null posts as adhesion control.
- **Profiling**: DiI membrane staining + confocal z-stacks (0.25 µm slices) to map vertical cell profile, aspect ratio vs substrate stiffness (soft posts G~5 kPa; stiff posts G~42 kPa; flat PDMS G~833 kPa), and quantify post invagination (~1 µm depth) via sidewall printing.
- **Statistics**: n = 14 control cells, 4 donors, 386 post trajectories; core vs perimeter posts dichotomized by trajectory variance; per-cell mean trajectories aligned on protrusive maximum (t=0); Tukey/Tukey-Kramer post-hoc, p<0.05.

## 4. Key results (quantitative)
- Protrusive force per post: **75 ± 8 pN/post** (abstract; ~75 pN core ≈ perimeter, not significantly different — Fig 3 D i).
- Spreading velocity (control): **206 ± 28 nm/s** (m±SE; Fig 3 C, p.702).
- Protrusion duration (FWHM): **~17 s** (Fig 3 D ii).
- Variance of force maxima: **~24 pN²** (Fig 3 D iii).
- Fraction of posts showing a protrusive spike: perimeter 0.67 ± 0.05 vs core 0.83 ± 0.05 (significant; Fig 3 D iv).
- Sustained contractile force: **core 20 ± 10 pN/post vs perimeter 106 ± 10 pN/post** (~5× difference; Fig 3 E i, p.702).
- Sustained-force variance: core 16 ± 4 vs perimeter 46 ± 4 pN²/post (~3×; Fig 3 E ii).
- Cytochalasin B spreading velocity: **61 ± 37 nm/s** (slowed; Fig 4 B) with increased protrusion FWHM, reduced perimeter steady-state contractility, eliminated core contractile rebound.
- Jasplakinolide: spreading completely abolished (Movie S2); cortex stiffening prevents spreading.
- CK666 (Arp2/3): no effect on protrusion forcefulness/duration → spreading ≠ lamellipodium formation; increased variance of core steady-state force.
- Blebbistatin + Y27632: significantly reduced perimeter steady-state contractility but did NOT eliminate post-protrusion contractile rebound → short-term rebound is actomyosin-independent; long-term contractility is ROCK/myosin-II-dependent.
- ROCK/myosin-II contractility not fully mature until **~500 s** after peak protrusive force (p.706).
- Post invagination depth **~1 µm** via sidewall printing (Fig 6 A).
- Contact-interface growth law: **R ~ t^0.4** (this work), vs Cuvelier prediction R~t^0.5 short-time / R~t^0.25 long-time; cf. neutrophils spread ~10× faster than mesenchymal cells on RGD lipid bilayers (~200 nm/s vs ~20 nm/s).
- Aspect ratio of spread cell increases monotonically with substrate stiffness (G 5 → 42 → 833 kPa; Fig 5 B); FN-null and anti-β2 controls had aspect ratio ≈ unity (no spreading).

## 5. Parameters & constants of interest
| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Per-post protrusive force | 75 ± 8 pN/post | Abstract, Fig 3 D i |
| Per-post sustained contractile force (perimeter) | 106 ± 10 pN/post | p.702, Fig 3 E i |
| Per-post sustained contractile force (core) | 20 ± 10 pN/post | p.702, Fig 3 E i |
| Spreading (protrusive-front) velocity, control | 206 ± 28 nm/s | p.702, Fig 3 C |
| Spreading velocity, cytochalasin-B | 61 ± 37 nm/s | Fig 4 B, p.703 |
| Protrusion duration (FWHM) | ~17 s | Fig 3 D ii |
| Micropost spring constant k_spring | 0.28 ± 0.09 pN/nm | Methods, p.700 |
| Post geometry | diam 604 ± 31 nm, length 5.576 ± 0.286 µm | Methods, p.700 |
| Force detection noise | 9 ± 2 pN | p.700 |
| Contractility maturation time | ~500 s after peak protrusion | p.706 |
| Post invagination depth | ~1 µm | Fig 6, p.706 |
| Substrate shear moduli probed | 5 / 42 / 833 kPa | p.705 |
| Contact-radius growth exponent | R ~ t^0.4 (vs t^0.5 / t^0.25 Cuvelier) | p.707 |
| Inhibitor concentrations | blebbistatin 5 µM, CK666 1 µM, cytochalasin B 3 µM, jasplakinolide 1 µM, Y27632 1 µM | Methods, p.700 |

## 6. Relevance to ffn_cellsim
This is primarily **(a) a validation oracle for cell-spreading + cortex relaxation dynamics** and **(b) a parameter/force-scale anchor**, with secondary value as **(c) a mechanism reference** for the cortex/clutch/myosin interplay.

- **Cell-spreading (the spherical→sessile-drop shape transition)**: The paper directly measures the force and velocity signatures of the quiescent-sphere → spread transition and frames it as a competition of **adhesion energy vs cortical tension vs cortical-tension relaxation**. This is exactly the energetic competition a fine-grained simulator's spreading run must reproduce. The R~t^0.4 (≈ t^0.5 short-time) contact-growth law, ~200 nm/s front velocity, and per-post ~75 pN protrusive scale are concrete acceptance bands for a spreading-on-FN validation test.
- **H.1 cortex**: The most mechanistically pointed finding for ffn_cellsim — cortex *stiffening* (jasplakinolide) abolishes spreading while *softening* (cytochalasin B, which reduces cortical tension) slows but does not stop it, and reduces protrusive force only weakly. This makes cortical-tension-relaxation a *prerequisite* gate for spreading: a clean qualitative oracle (stiff cortex → no spread; relaxed cortex → spread) for the H.1 cortex shell coupled to a spreading driver.
- **motor/myosin & FA/clutch**: The clean temporal separation — protrusion is myosin-II/ROCK/Arp2/3-independent, whereas the *sustained* contractile force floor (~100 pN/post, matured ~500 s) is ROCK/myosin-II-dependent — is a strong mechanistic validation target. A mechanistic model with Stam-Hocky myosin minifilaments + Bell-Evans/integrin clutches should reproduce: (i) early integrin-ligation-driven protrusion without mature actomyosin, and (ii) late actomyosin-built perimeter contraction. The core-vs-perimeter stratification (perimeter 5× more contractile) maps onto an FA-clutch + cortical-tension spatial gradient.
- **ECM/collagen / adhesion**: Spreading is haptokinetic and strictly β2-integrin (MAC-1, αMβ2) / FN-dependent (FN-null and anti-β2 → no spreading, aspect ratio≈1). This anchors the "integrin ligation is upstream of spreading" requirement; relevant to the FA/clutch FN-ligation logic, though on β2 not the β1/MCF7 axis of the PI experimental overlay.
- **Caveats on transfer**: This is a **neutrophil** (liquid-drop-like, thin cortical-shell leukocyte), not MCF7 or a mesenchymal cell with stress fibers. Forces are reported **per micropost**, not as a continuum traction-stress field, so converting to stress requires post density. The cell is modeled by others as a viscous liquid drop — a coarser picture than ffn_cellsim's explicit filament/motor/clutch fidelity. So it is a force-magnitude and spreading-dynamics anchor, not a structural template.
- **Numerics/methods**: Particle-tracking + spring-constant × deflection is the canonical traction readout; the substrate-warping correction (Schoen, 8%) is a useful note if the project ever simulates micropost substrates.

Not off-topic, but **Medium** rather than High because (i) it is a leukocyte not the project's MCF7/mesenchymal cancer-cell target, (ii) it reports per-post forces not a continuum traction field on collagen-I, and (iii) it offers validation bands and force scales rather than mechanistic rate constants that drop straight into the runtime.

## 7. Limitations & caveats
- **Cell type**: human neutrophil — a near-liquid-drop leukocyte with a thin cortical actin shell and *no* prominent stress-fiber apparatus; force scales (tens-of-pN/post) and the liquid-drop framing may not transfer to MCF7/mesenchymal cells with mature stress fibers.
- **Substrate is FN-coated PDMS microposts, not collagen-I**: integrin axis is β2/MAC-1, not the β1/collagen-I axis the PI experimental overlay uses.
- **Per-post, discretized force**, not a continuum traction-stress map; comparison to TFM-on-gel requires post-density normalization. The cell tracks post tips, not the membrane, so spreading velocity is inferred from protrusive-front propagation, not membrane edge.
- **Sidewall printing / invagination** confound: forces partly arise because posts physically reside inside the cell's finite-thickness spreading path (not pure top-plane wetting); the authors note this was "fortuitous" but it complicates a clean wetting interpretation.
- **No simulation / no fitted mechanistic rate constants**: provides phenomenology and force/velocity scales, not Bell-Evans off-rates, Hill parameters, or motor stiffnesses.
- Authors explicitly state a *purely physical* (Cuvelier liquid-drop) treatment is insufficient — cell signaling (integrin ligation) is required — so the oracle is partly qualitative.

## 8. Key figures / tables
- **Fig 3** — Ensemble force trajectories: (A) core vs perimeter mean radial force vs time with protrusive (cyan) and contractile (lavender) regimes; (C) protrusion-max time vs radial distance → spreading velocity; (D) protrusion metrics (force max ~75 pN, FWHM ~17 s, variance, participation fraction); (E) steady-state contraction metrics (core 20 vs perimeter 106 pN/post). The core quantitative figure.
- **Fig 4** — Inhibitor panel: mean radial force trajectories under blebbistatin/CK666/cytochalasin B/jasplakinolide/Y27632, plus effect on spreading velocity (B) and protrusion/contraction metrics (C). The mechanistic-dissection figure.
- **Fig 5** — DiI confocal vertical profiles vs substrate stiffness (G 5/42/833 kPa); aspect ratio increases monotonically with stiffness; FN-null & anti-β2 controls ≈ unity.
- **Fig 6** — Post invagination (~1 µm) as origin of protrusion; schematic of finite-volume spreading + contractile-rebound conjecture from membrane curvature.

## 9. Notable quotes / citable claims
- "During spreading, a wave of protrusive force (75 ± 8 pN/post) propagates radially outward from the cell center at a speed of 206 ± 28 nm/s." (Abstract)
- "posts within the core of the contact zone were less contractile (20 ± 10 pN/post) than those residing at the geometric perimeter (106 ± 10 pN/post)." (Abstract)
- "Relaxation of the actin cortical shell was a prerequisite for spreading on post arrays as demonstrated by stiffening in response to jasplakinolide and the abrogation of spreading." (Abstract)
- "Our work suggests a competition of adhesion energy, cortical tension, and the relaxation of cortical tension is at play at the onset of neutrophil spreading." (Abstract)
- "That CK666 did not abrogate protrusion suggests the shape change associated with spreading was not analogous to lamellipodium formation… The implication of this result is that the short-term transient rebound is not actomyosin dependent." (p.704)

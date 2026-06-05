---
id: 45_moving-cell-boundaries-drive-nuclear-shaping-during-cell-spr
paper_n: 45
title: "Moving Cell Boundaries Drive Nuclear Shaping during Cell Spreading"
authors: "Li, Y. et al. (Lele & Dickinson)"
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.07.006"
paper_type: continuum-model
ffn_relevance: Medium
ffn_themes: [nucleus, cell-spreading, cortex, FA/clutch, validation-oracle, parameter-source, numerics/methods]
entities: [nucleus, nuclear-lamina, lamin-a-c, actin, myosin-ii, microtubule, intermediate-filament, vimentin, linc-complex, nesprin-2g, sun2, kash4, integrin, fibronectin, nih3t3, mef]
methods: [continuum-model, boundary-element-method, confocal-live-imaging, immunocytochemistry, drug-perturbation, rna-interference, finite-element]
measurables: [nuclear-height, nuclear-width, nuclear-aspect-ratio, spread-area, retrograde-flow, membrane-tension, elastic-modulus, contractile-stress, area-modulus, bulk-modulus]
keywords: [nuclear-flattening, cell-spreading, retrograde-flow, cytoskeletal-network-viscosity, nuclear-lamina-tension, LINC-complex, mechanotransduction, boundary-element-method, molecular-clutch, contractile-compressible-network]
tags: ["#nucleus-mechanics", "#cell-spreading", "#continuum-model", "#parameter-source", "#validation-oracle", "#retrograde-flow", "#lamina-tension", "#boundary-element-method"]
has_transferable_params: true
---

# [45] Moving Cell Boundaries Drive Nuclear Shaping during Cell Spreading

**Tags:** #nucleus-mechanics #cell-spreading #continuum-model #parameter-source #validation-oracle #retrograde-flow #lamina-tension #boundary-element-method

| Field | Value |
|---|---|
| Authors | Li, Y., Lovett, D., Zhang, Q., Neelam, S., Kuchibhotla, R.A., Zhu, R., Gundersen, G.G., Lele, T.P., Dickinson, R.B. |
| Year / Venue | 2015 / Biophysical Journal 109(4):670–686 |
| DOI / ID | 10.1016/j.bpj.2015.07.006 |
| Type | continuum-model (boundary-element) + experimental |
| Pages | 17 |
| ffn_cellsim relevance | Medium — directly informs the EXTEND nucleus unit (H.9): mechanism that cell-boundary motion + cytoskeletal-network flow flattens the nucleus, plus a table of nuclear/cortex/membrane constants and a clean validation oracle (nuclear aspect ratio vs. spread area). |

## 1. Summary
The authors combine confocal live-cell imaging of NIH 3T3 fibroblasts (and MEFs) spreading on fibronectin with a continuum boundary-element mechanical model to ask how the nucleus becomes flattened during cell spreading. They find nuclear height/aspect ratio correlate tightly with the degree of cell spreading and reach steady state early (~20–30 min, when the cell has spread to only ~50% of final area). Strikingly, nuclear flattening does NOT require actomyosin contraction (blebbistatin, Y-27632), myosin light-chain kinase (except where ML-7 blocks spreading itself), microtubules (nocodazole/colcemid), intermediate filaments (vimentin-null MEFs), or an intact LINC complex (KASH4 overexpression, nesprin-2G/SUN2 knockdown) — as long as the cell can still spread. Detachment (trypsin) reverses flattening and rounds the nucleus, tracking cell rounding in real time. To explain this, they model the cytoplasm as a contractile, compressible, viscous network (stress ∝ rate-of-strain) that frictionally transmits stress from the moving cell membrane to the nuclear surface; the nucleus resists volume change (bulk modulus K) and lamina area expansion (area modulus κ). The model reproduces the three observed phases (vertical distension/translation toward substratum, apical collapse/flattening, lateral widening) with essentially one fitted parameter (cortical assembly speed v_a).

## 2. Problem & motivation
Nuclear shape is smooth/regular in normal cells but altered in cancer and laminopathies, and shape may directly control gene expression by reorganizing chromatin access. Yet the mechanical mechanism that flattens the nucleus during cell spreading is unknown. The prevailing assumption is that specific cytoskeletal force-generators (actomyosin stress fibers, dynein/microtubules, intermediate filaments) compress or pull the nucleus via the LINC complex. The paper tests each candidate and proposes a simpler alternative: boundary-driven network flow.

## 3. Methods / model
- **Experimental**: NIH 3T3 fibroblasts and MEFs spread on 5 µg/ml fibronectin-coated glass; x-z laser-scanning confocal (Nikon A1, 60×/1.40NA) of GFP-histone H1 nuclei. Nuclear height by FWHM of x-z intensity (MATLAB); aspect ratio = height / x-y major axis. Drug perturbations: Y-27632 (25 µM, ROCK), blebbistatin (50 µM, myosin-II ATPase), ML-7 (25 µM, MLCK), cytochalasin-D (2 µM)/latrunculin-A (5 µM) (F-actin), nocodazole (0.83–1.65 µM)/colcemid (0.27 µM) (microtubules). Genetic: vimentin-null MEFs, LMNA−/− MEFs, GFP-KASH4 overexpression, shRNA to nesprin-2G and SUN2. Detachment via 0.25% and 0.08% trypsin.
- **Model class**: quasistatic continuum, contractile compressible viscous network. Constitutive stress (Eq. 1/13): σ = 2µ·ε̇ + (σ_c + λ∇·v)·I, with ε̇ the rate-of-strain tensor; assume λ≈0 (Poisson ratio 0 ⇒ E = 2µ). Momentum balance ∇·σ = 0 solved as the elastostatic analog with shear modulus → µ, Poisson ratio → 0.
- **Nucleus mechanics**: internal tension τ_nuc = K·ln(V/V0) (Eq. 10); surface (lamina) tension via vesicle-undulation relation (Eq. 11) using area modulus κ, bending modulus κ_c, undulation energy E_s.
- **Boundary conditions**: substratum exerts tangential traction η·v(z=0) (η→∞ no-slip, 1/η→0 perfect adhesion; finite η = molecular-clutch-like slip); slip on cell/nuclear membranes; constant cell volume enforced by hydrostatic pressure P_h; cortical actin assembly speed v_a normal to membrane, contact-boundary assembly v_ca tangential.
- **Numerics**: axisymmetric boundary-element method (Kelvin fundamental solutions), cell surface = 100 quadratic elements, nucleus = 50; 10th-order Gauss quadrature; 4th-order Runge-Kutta time stepping; node respacing via cubic Hermite (pchip); short-range repulsive pressure P(z)=(δ/z)³·e^(−z/δ), δ=0.01 R_n, to prevent nucleus-membrane contact singularities.
- **Scales**: cell radius R≈8.3 µm, nucleus R_n≈6.3 µm; spreading timescale several minutes; spread area grows ~200 → ~1400 µm² over 60 min.

## 4. Key results (quantitative)
- Nuclear initial translation toward substratum is ~20× faster than gravitational settling would predict given the assumed cytoskeleton viscosity (p.674).
- Apical-surface collapse phase lasts ~5–6 min; gap between cell apex and nuclear apex peaks at ~2 µm (control) and up to ~4 µm under myosin inhibition (Figs. 1–2, S3).
- Cell spread area rises from <200 µm² to ~1400 µm² in first 60 min; nuclei flatten to steady-state height by 20–30 min when cell is only ~50% spread (Fig. 2 A,B).
- Aspect ratio drops to ~0.25 by 30 min (width ≈ 4× height) (Fig. 2 D).
- vimentin−/− MEFs: nucleus less flattened at 1 h (aspect ratio ~0.3, still significantly flat) but fully flattened by 12 h (Fig. 4).
- LMNA−/− MEFs flatten faster / flatter than WT at 60 min (Fig. 8) — consistent with model where setting area modulus κ=0 lets the nucleus keep flattening at constant volume (Fig. 10 D).
- Apical actomyosin bundles ~0.2 per cell in first 20 min (1 in 5 cells), basal ~1–2 per cell at flattening time — i.e., bundles largely absent when flattening occurs (Fig. 5).
- Trypsin (0.25%) rounds nucleus in seconds tracking cell rounding; 0.08% trypsin rounds over several minutes, nuclear height tracking cell contact length (Fig. 6).
- Model: cortical network assembly (v_a>0) is required to reproduce the rapid early downward translation, but the final flat shape only requires assembly at the contact boundary; differential (apex-vs-side) stress, not background tension σ_c, drives shape change (Fig. 10 C with σ_c=0 still flattens).

## 5. Parameters & constants of interest
Table 1 (p.673) — directly transferable to an EXTEND nucleus/cortex model and as a validation oracle:

| Quantity | Value + units | Source in paper |
|---|---|---|
| Contractile stress σ_c | 0.19 nN/µm² (= 190 Pa) | estimated from myosin-induced nuclear volume change |
| Nucleus bulk modulus K | 0.25 nN/µm² (= 250 Pa) | Dahl et al. 2004 (ref 32), isolated Xenopus oocyte nuclei |
| Nucleus area modulus κ | 25 nN/µm (= 25 mN/m) | Dahl et al. 2004 (ref 32) |
| Membrane tension T_mem | 0.1 nN/µm (= 0.1 mN/m); literature range 0.01–0.3 | Lieber 2013 (33, keratocytes); Gauthier 2012 (34) |
| Network viscosity µ | 0.21 nN·s/µm² (= 0.21 kPa·s) | Bausch 1999 (35), J774 macrophages, magnetic tweezers |
| Lamina bending stiffness κ_c | 3.5×10⁻⁴ nN·µm (≈ 8.6×10⁴ k_BT·nm... ≈ 350 k_BT, approx) | Vaziri 2006 (36), MEFs |
| Energy parameter E_s | 3.2×10⁻⁴ nN·µm ≈ 100 k_BT (stated) | estimated from excess surface area |
| Nuclear radius R_n | 6.3 µm | measured |
| Cell radius (rounded) R | 8.3 µm | measured |
| Contact-boundary assembly speed v_ca | 0.5 µm/min | estimated from spreading speed |
| Surface friction coefficient η | varied (clutch-like) | — |
| Cortical assembly speed v_a | varied (cases v_a = v_ca and v_a = 0) | unknown / fitted |
| Steady-state aspect ratio | ~0.25 | Fig. 2 D |
| Excess nuclear surface area | ~20–40% | Methods (Eq. 11 context) |

## 6. Relevance to ffn_cellsim
This is the most directly relevant paper so far for the **EXTEND nucleus unit (H.9)** and partially for cortex (H.1) and FA/clutch (H.4). Concretely:

- **(b) Parameter source — strong.** Table 1 gives a self-consistent set of nuclear and peripheral mechanical constants in clean SI-convertible units: nucleus bulk modulus K=250 Pa, area modulus κ=25 mN/m, lamina bending stiffness κ_c=3.5×10⁻⁴ nN·µm, membrane tension T_mem≈0.1 mN/m, cytoplasmic network viscosity µ=0.21 kPa·s, contractile stress σ_c≈190 Pa. These are exactly the kind of anchors needed when ffn_cellsim adds an explicit nuclear envelope/lamina (lamin shell) and a membrane. The lamina constitutive form (excess-area smoothing → stiff true-area regime, Eq. 11, with E_s≈100 k_BT) is a ready prescription for a coarse nuclear-surface element.

- **(a) Validation oracle — strong, but as a behavioral/scaling oracle, not a runtime mechanism.** The headline empirical law — nuclear aspect ratio (and height) is a monotone function of cell spread area, reaching steady state at ~50% spread, with steady aspect ratio ~0.25 — is a clean, parameter-light validation target. If ffn_cellsim's fine-grained cytoskeleton + (future) nucleus reproduces "nucleus flattens with spreading, independent of myosin/microtubule/IF/LINC perturbation as long as spreading occurs," that is a non-trivial pass. The trypsin-detachment reversibility (nuclear rounding tracks cell rounding) is a second behavioral oracle.

- **(c) Mechanism reference — important framing, opposite-fidelity.** The core claim is mechanistically interesting for the simulator: nuclear shaping is driven by *frictional stress transmission through the intervening cytoskeletal network as cell boundaries move*, NOT by dedicated motors/LINC pulling. In a fine-grained model this is emergent — if cortex assembly/retrograde flow and an excluded-volume nucleus are present, the simulator should reproduce flattening without needing an explicit "nuclear-shaping" mechanism. This is a useful negative-control hypothesis: ffn_cellsim should NOT need LINC/dynein wiring to flatten the nucleus during spreading.

- **(e) Numerics/methods context.** The model is a *lumped continuum boundary-element* model — exactly the abstracted style ffn_cellsim rejects at runtime (compressible viscous network with one viscosity µ, BEM on 100/50 elements). It belongs in the **acceptance-oracle** category, not the runtime. Its retrograde-flow / molecular-clutch boundary condition (finite η, perfect-adhesion limit 1/η→0) maps conceptually to the H.4 clutch but is far coarser than Bell-Evans clutches.

- **Cell-spreading link (cell-spreading theme).** It quantifies the spreading process itself (area vs. time, v_ca≈0.5 µm/min initial spreading speed) which connects to the broader single-cell spreading context, and to the MCF7-on-fibronectin/collagen validation track at the single-cell level (these are NIH3T3/MEF on fibronectin, not MCF7, so it is a parameter/behavior analog, not a direct overlay).

Net: keep as a Medium-relevance EXTEND-nucleus parameter source + behavioral oracle. Off-target for the *current* H.1→H.7 cytoskeleton chain (no actin filament/motor/clutch microphysics here), but squarely on-target once H.9 (nucleus) is in scope.

## 7. Limitations & caveats
- The cytoplasm is a single-phase, single-viscosity (µ), λ≈0, σ_c-uniform continuum — explicitly a "slow-flow limit" of two-phase reactive flow (Dembo). It contains no explicit filaments, motors, crosslinks, or clutches; it is the lumped style ffn_cellsim inverts. Use only as oracle.
- Nuclear constants K, κ are from isolated *Xenopus oocyte* nuclei (Dahl 2004) — large amphibian oocyte nuclei, possibly not quantitatively transferable to mammalian somatic nuclei; κ_c from MEFs (Vaziri 2006).
- Essentially one fitted free parameter (v_a, cortical assembly speed); v_a is "not known," and the headline conclusions are stated to be robust to most parameters — so the model is more a qualitative/scaling explanation than a tightly constrained quantitative fit.
- The mechanism is geometry/boundary-driven; it deliberately downplays force-generators, so it under-describes the well-spread, stress-fiber-laden steady state (authors note later actomyosin-bundle distortion is a separate regime, refs 50–51).
- 2D-axisymmetric; chromatin, nucleolus, and active nuclear remodeling are not represented.
- Cells are NIH3T3 / MEF fibroblasts on fibronectin-coated glass — not MCF7, not collagen-I, not the PI experimental platform; overlap is mechanism-level, not a direct data overlay.

## 8. Key figures / tables
- **Table 1 (p.673)** — model parameters with sources; the transferable-constant payload.
- **Fig. 2 (p.674–675)** — nuclear height/width/aspect-ratio vs. time and vs. spread area; aspect ratio → ~0.25, steady at ~50% spread. The validation-oracle figure.
- **Fig. 9 (p.681–682)** — the model: 1D gap stress-transmission schematic (Eqs. 14–15), full spreading model components, and simulation snapshots of the three flattening phases.
- **Fig. 10 (p.683–684)** — parameter-sweep simulations: no cortical assembly, reduced adhesion (clutch), σ_c=0, and zero lamina area modulus (constant-volume continued flattening, ~ LMNA−/−).

## 9. Notable quotes / citable claims
- "cell spreading is necessary and sufficient to drive nuclear flattening under a wide range of conditions, including in the presence or absence of myosin activity" (Abstract, p.670).
- "translating the one boundary at x=L at speed V transmits an additional stress 2µV/L to the surface at x=0 because of longitudinal friction" (Methods, Eq. 4 context, p.671–672).
- "When the area modulus was set to zero, the nucleus continued to flatten as long as the cell continued spreading … consistent with the increased nuclear flattening in LMNA−/− MEFs" (Results, p.683).
- "nuclear shape changes result from transmission of stress from the moving cell boundary to the nuclear surface because of frictional resistance to expansion/compression of the intervening cytoskeletal network" (Discussion, p.685).

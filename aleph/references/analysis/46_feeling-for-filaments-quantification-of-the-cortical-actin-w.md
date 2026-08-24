---
id: 46_feeling-for-filaments-quantification-of-the-cortical-actin-w
paper_n: 46
title: "Feeling for Filaments: Quantification of the Cortical Actin Web in Live Vascular Endothelium"
authors: "Kronlage C. et al."
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.06.066"
paper_type: experimental
ffn_relevance: Medium
ffn_themes: [cortex, parameter-source, numerics/methods]
entities: [actin, f-actin, g-actin, spectrin, myosin-ii, stress-fiber, transverse-arc, cortical-actin, endothelial-cell, plasma-membrane, lifeact, cytochalasin-d, jasplakinolide]
methods: [atomic-force-microscopy, confocal-live-imaging, force-mapping, image-segmentation, random-forest-classifier, surface-roughness-morphometrics]
measurables: [mesh-size, hole-area, filament-spacing, surface-coverage, rms-roughness, spatial-resolution, indentation-force, intracellular-calcium]
keywords: [cortical-actin-web, mesh-size, coarse-mesh, fine-mesh, AFM-force-mapping, endothelial-cortex, cytochalasin-D, jasplakinolide, Lifeact, live-cell-imaging]
tags: ["#cortex-structure", "#cortical-actin", "#mesh-size", "#parameter-source", "#AFM", "#endothelial", "#live-cell-imaging"]
has_transferable_params: true
---

# [46] Feeling for Filaments: Quantification of the Cortical Actin Web in Live Vascular Endothelium

**Tags:** #cortex-structure #cortical-actin #mesh-size #parameter-source #AFM #endothelial #live-cell-imaging

| Field | Value |
|---|---|
| Authors | Kronlage C., Schäfer-Herte M., Böning D., Oberleithner H., Fels J. et al. |
| Year / Venue | 2015 / Biophysical Journal 109(4):687–698 |
| DOI / ID | 10.1016/j.bpj.2015.06.066 |
| Type | experimental (AFM imaging + image-analysis methods) |
| Pages | 12 |
| ffn_cellsim relevance | Medium — quantitative structural metrics (mesh size, hole area, filament spacing) of the live apical cortical actin web; a cortex-geometry anchor, not a mechanism oracle |

## 1. Summary
The authors imaged the apical cortical actin web of live bovine aortic endothelial cells (GM7373) by atomic force microscopy (AFM), correlated it with simultaneous confocal Lifeact F-actin fluorescence, and developed an image-processing pipeline (line-wise 2nd-order polynomial subtraction + difference-of-Gaussians filter on force-mapping/QI height data) to quantify cortical cytoskeleton density. They show the cortex is a two-scale meshwork — a "coarse mesh" (confirmed as F-actin by Lifeact co-localization, akin to dorsal transverse arcs/stress fibers) and a faster-remodeling "fine mesh." Pharmacological perturbation validated the method: the actin-depolymerizer cytochalasin D (100 nM) thins the network (fewer, larger holes; higher RMS roughness), while the actin-stabilizer jasplakinolide (1 µM) densifies it (more, smaller holes; lower roughness). They derive quantitative readouts: hole count, mean hole area, surface coverage, and an "adapted RMS roughness," reporting baseline cortical mesh statistics for live cells.

## 2. Problem & motivation
The cortical actin network underlies endothelial stiffness, barrier permeability, and mechanotransduction (e.g. nitric-oxide release), so its ultrastructural architecture matters mechanically. Existing AFM "error-signal" imaging gives qualitative but not quantitative views of the cortex, and electron/super-resolution methods generally require fixation. The paper's goal is a method to image and quantify cortical actin density and its remodeling in live cells, just beneath the apical plasma membrane.

## 3. Methods / model
- **Cells:** bovine aortic endothelial cells GM7373, ~90% confluence on glass-bottom dishes, room temperature in HEPES buffer (140 NaCl, 5 KCl, 1 MgCl2, 1 CaCl2, 5 glucose, 10 HEPES, pH 7.4). Lifeact-eGFP (transient) or Lifeact-mKate2 (stable) for F-actin.
- **AFM:** JPK NanoWizard 3, MSCT-UC B-probes, 10 nm tip radius, nominal cantilever spring constant 0.02 N/m, uncoated (low force drift). Two modes: (i) contact-mode (10 µm)², 256×256 px, 2 Hz line rate (~2.5 min/image); (ii) force-mapping / QI mode (6×6 µm), 256×256 px, z-range 50 nm, extend/retract 1.8 ms (tip velocity 27.8 µm/s), ~10 min/image. Optimal imaging force ~0.5 nN.
- **Correlative imaging:** confocal (Leica TCS SP8, 63×/NA1.4) Lifeact maximum-intensity projection of apical planes, simultaneous with AFM.
- **Image processing:** Fiji/ImageJ for overlays; MATLAB R2014a for line-wise 2nd-order polynomial subtraction + difference-of-Gaussians (Gaussian SD σ=5 px smoothed minus σ=1 px), then "adapted RMS roughness" = SD of image values after >4σ outlier removal.
- **Morphometrics:** ilastik 1.1 trainable random-forest segmentation into network vs filament-free "holes"; MATLAB extracts hole count, mean hole area, surface coverage.
- **Controls:** Fluo4 intracellular Ca²⁺ imaging to test mechanosensitive disturbance (ionomycin 2.5 µM positive control); statistics via Student t-test + Mann-Whitney (p≤0.05 in both).

## 4. Key results (quantitative)
- Baseline apical cortical mesh: **271 ± 27.7 holes per 6×6 µm area** (p.693).
- Baseline mean hole area: **42,045 ± 5,553.2 nm²** (≈0.042 µm², i.e. linear scale ~0.2 µm) (p.693).
- **Cytochalasin D (100 nM):** holes per image −99 ± 31.3; mean hole area +37,886 ± 10,898.1 nm²; filament-free surface coverage +30%; adapted RMS roughness +6.8 ± 2.2 a.u. (Fig 7a, Fig 8a/c/e; p.692–694). n=8 treatment, n=7 control.
- **Jasplakinolide (1 µM):** holes +33 ± 16.4; mean hole area −12,160 ± 4,786.4 nm²; coverage −18%; adapted RMS roughness −3.3 ± 0.9 a.u. (Fig 7b, Fig 8b/d/f; p.693–694). n=6 treatment, n=5 control. Cells lose mechanical integrity after ~15 min.
- **Spatial resolution:** contact-mode ~150 nm along fast-scan axis; force-mapping mode ~120 nm; smallest resolvable filament spacings ~100–150 nm (an upper bound on resolution) (Fig S4; p.692, 696).
- **Force-error magnitude:** at a 0.4 nN setpoint, contact-mode force-error excursions translate to apparent surface height differences of **~400 nm** (Fig 3b/c; p.691) — i.e. contact mode mis-tracks soft filaments.
- Coarse mesh = F-actin (Lifeact co-localization); fine mesh only diffuse fluorescence (unresolved actin or possibly spectrin). Coarse mesh trackable across ~2.5 min frames; fine mesh remodels faster than acquisition.
- Force-mapping imaging did not change intracellular Ca²⁺ (Fig S2) — no mechanosensitive activation; cells viable for hours.

## 5. Parameters & constants of interest
| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Cortical mesh hole density (baseline) | 271 ± 27.7 holes / 36 µm² (≈7.5 holes/µm²) | p.693 |
| Mean cortical mesh hole area (baseline) | 42,045 ± 5,553 nm² (≈0.042 µm²) | p.693 |
| Implied mesh pore linear scale | ~200 nm (√area, inferred) | derived from p.693 |
| Smallest resolved filament spacing | ~100–150 nm | Fig S4; p.692, 696 |
| AFM spatial resolution (force-map) | ~120 nm | p.692 |
| AFM spatial resolution (contact) | ~150 nm | p.689 |
| Cell-surface height step from a cortical fiber | ~400 nm at 0.4 nN (apparent) | Fig 3; p.691 |
| Optimal cortical imaging force | ~0.5 nN (setpoint 0.4 nN) | p.688, 691 |
| Tip radius of curvature | 10 nm | p.688 |
| Cantilever spring constant | 0.02 N/m | p.688 |
| Actin fraction of total protein (endothelium) | 5–15% (G+F-actin) | Intro, p.687, cites ref 1 |
| Coarse-mesh remodeling timescale | trackable over ~2.5 min; fine mesh faster | p.690 |

## 6. Relevance to ffn_cellsim
**Primary role: (b) parameter source / cortex-geometry anchor + secondary (d) structural context.** This is directly on-theme for **H.1 (cortex)**. The fine-grained cortex in ffn_cellsim is an explicit actin filament network beneath the apical membrane; this paper supplies *live-cell, apical cortex* structural metrics to sanity-check the simulated cortex's emergent geometry:
- **Mesh pore / hole size** (~0.04 µm², ~200 nm linear scale; ~7.5 holes/µm²) is a candidate **validation target for the simulated cortical mesh** — does the H.1 actin web reproduce a comparable pore-size distribution and areal hole coverage? This is a structural oracle, not a closed-form mechanism oracle.
- **Two-scale architecture** (coarse F-actin "stress-fiber/transverse-arc"-like bundles vs. finer faster-remodeling mesh) maps onto the distinction between bundled (myosin-cross-linked) and dispersed cortical actin — relevant to how Stam-Hocky myosin minifilaments bundle the simulated cortex (the paper explicitly attributes coarse-mesh emergence to "myosin-based bundling of actin filaments, Rho-driven," p.696).
- **Perturbation directionality** is a useful qualitative behavioral check: depolymerization → larger pores / lower density; stabilization → smaller, denser pores. A mechanistic cortex with tunable polymerization should reproduce this monotonic trend.
- **Filament spacing ~100–150 nm** gives a real-world lower bound to compare against the ×40 mesoscopic coarse-graining (~1000 effective filaments/cell) — useful to gauge whether mesoscopic pore sizes are physically plausible.

It is **not** a mechanism source (no force-velocity law, no off-rate, no elastic modulus extracted — the authors note their partial force curves are unsuitable for elasticity) and **not** a numerics/integrator reference for the runtime. Note also the cell type is **vascular endothelium**, not the MCF7 breast-cancer validation system, so any geometric anchor is cross-cell-type and should be treated as order-of-magnitude only.

## 7. Limitations & caveats
- AFM resolution (~100–150 nm) is an *upper bound* limited by tip+membrane convolution and lateral sample shift; the true cortical mesh is denser (REM electron micrographs resolve finer; ref 41). So reported hole areas are coarse-graining-biased toward larger pores.
- Time resolution (~2.5–10 min/image) is far slower than fine-mesh remodeling, so measured "holes" don't reflect instantaneous architecture — morphometrics are population/statistical, not a frozen network.
- "Adapted RMS roughness" is in arbitrary units (processing destroys absolute height), and roughness also reflects microvilli, ruffles, membrane composition, motor activity — not actin density alone; values are reported as differences from per-cell control, not absolutes.
- Fine-mesh molecular identity uncertain (actin vs spectrin); only coarse mesh confirmed as F-actin.
- Endothelial cell, apical cortex only — no basal/lamellar cortex, no ECM, no traction; nothing about cross-link or motor kinetics that a mechanistic simulator would integrate.

## 8. Key figures / tables
- **Fig 4** — force-mapping height data → DoG-processed image → ilastik hole segmentation; the core pipeline + Lifeact correlation. Defines the mesh-quantification method.
- **Fig 7** — adapted RMS roughness time-courses under cytochalasin D vs jasplakinolide (opposite-sign density changes); the validation that the metric tracks polymerization state.
- **Fig 8** — morphometric panels: hole number, mean hole area, surface coverage for both drugs; source of the quantitative mesh parameters.
- **Fig 3** — contact-mode force error (~400 nm apparent height at 0.4 nN) demonstrating why error-signal imaging is non-quantitative.

## 9. Notable quotes / citable claims
- "apical cortical cytoskeleton exhibited an average number of 271 ± 27.7 holes per 6 × 6 µm surface area … the area of each hole increased (from basal 42,045 ± 5553.2 nm²)" (p.693) — baseline cortical mesh metrics.
- "The coarse mesh seen in the AFM images … can thus be identified as F-actin" (p.689) — coarse mesh = actin, by Lifeact co-localization.
- "coarse-mesh fibers gradually emerge from smaller structures can be interpreted as the formation of stress fibers by myosin-based bundling of actin filaments, an effect based on Rho signaling" (p.696) — mechanism attribution for bundling.
- "the data acquired with our experimental setup are unsuitable for this kind of measurement because they contain only partial force curves" (p.691) — explicitly no elasticity extracted (so not an elastic-modulus source).

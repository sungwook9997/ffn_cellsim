---
id: 09_directed-invasion-of-cancer-cell-spheroids-inside-3d-collage
paper_n: 9
title: "Directed invasion of cancer cell spheroids inside 3D collagen matrices oriented by microfluidic flow in experiment and simulation"
authors: "Geiger F, Schnitzler LG, Brugger MS, Westerhausen C, Engelke H"
year: "2022"
venue: "PLOS ONE 17(3):e0264571"
doi: "10.1371/journal.pone.0264571"
paper_type: experimental
ffn_relevance: Low
ffn_themes: [ECM/collagen, spheroid-scale-context, tangential]
entities: [collagen-i, hela, cancer-cell-spheroid, rat-tail-collagen]
methods: [microfluidics, particle-image-velocimetry, finite-element, random-walk, contact-guidance, light-microscopy, confocal-microscopy]
measurables: [invasion-speed, invasion-distance, fiber-orientation, degree-of-alignment, step-size, flow-rate]
keywords: [collagen-fiber-alignment, contact-guidance, spheroid-invasion, tangential-vs-radial, microfluidic-shear, brownian-diffusion-model, anisotropic-random-walk, HeLa, YAP, unjamming]
tags: ["#spheroid-mechanics", "#ecm-collagen", "#contact-guidance", "#fiber-alignment", "#anisotropic-migration", "#scale-context"]
has_transferable_params: false
---

# [9] Directed invasion of cancer cell spheroids inside 3D collagen matrices oriented by microfluidic flow in experiment and simulation

**Tags:** #spheroid-mechanics #ecm-collagen #contact-guidance #fiber-alignment #anisotropic-migration #scale-context

| Field | Value |
|---|---|
| Authors | Florian Geiger, Lukas G. Schnitzler, Manuel S. Brugger, Christoph Westerhausen, Hanna Engelke |
| Year / Venue | 2022 / PLOS ONE 17(3):e0264571 |
| DOI / ID | 10.1371/journal.pone.0264571 |
| Type | Experimental (microfluidics + live imaging) + phenomenological simulation (FEM flow/fiber + 2D anisotropic random walk) |
| Pages | 15 |
| ffn_cellsim relevance | Low — tissue/spheroid-scale invasion phenomenology and ECM-fiber contact guidance; NOT a single-cell mechanistic model, no transferable cytoskeletal/clutch constants |

## 1. Summary
The authors use shear flow in a microfluidic channel during collagen-I polymerization to produce highly oriented collagen matrices around embedded HeLa cancer-cell spheroids, yielding fibers that are tangential on the upstream face and radial on the downstream face of each spheroid. Over 3 days they track invasion and find a strong directional bias: cells invade ~2–3× farther into the radially-oriented (downstream) fibers than into tangential (upstream) fibers, even though collagen density is homogeneous around the spheroid. A COMSOL FEM flow simulation reproduces the bead-PIV velocity field and predicts the fiber orientation by treating a fiber as a rigid rod in the flow field. The invasion itself is modeled with a deliberately minimal 2D anisotropic random walk: cells take steps of size k_parallel along the local fiber direction and k_perp across it; the best fit to the measured spheroid shape evolution is obtained with k_perp = 0, i.e. complete blockage of migration perpendicular to fibers (migration exclusively along fibers). This slows but does not prevent invasion toward tangentially-oriented regions, and predicts invasion distance growing as sqrt(time) (diffusive character).

## 2. Problem & motivation
Invasion from tumor spheroids depends on ECM mechanics, not only cell signaling. Prior work showed fiber alignment affects single-cell migration via contact guidance and that radial alignment from spheroid contraction increases invasion — but cell-contraction remodeling only ever produces radial alignment, so tangential vs radial could not be compared. The paper's contribution is an experimental platform (microfluidic shear) that imposes both tangential and radial fiber orientations around a spheroid independently of cell-driven remodeling, isolating the pure effect of fiber orientation on collective invasion.

## 3. Methods / model
- **Cells/spheroids**: HeLa, 500 cells seeded in ultra-low-adhesion 96-well, grown to 250–350 µm diameter spheroids.
- **Matrix**: rat-tail collagen-I (Corning), stock 3.32 mg/ml diluted to final 1.85 mg/ml, ATTO-633 labeled, polymerized on ice then 37 °C.
- **Alignment**: 400 µm-high microfluidic channel (µ-Slide VI 0.4, Ibidi); high draw (90.2 µl/min) to fill, then 0.2 µl/min during polymerization to align fibers via shear; flow stopped after gelation.
- **Flow characterization**: DPIV on 3 µm latex tracer beads (1 fps, 0.2 µl/min), PIVlab/MATLAB; COMSOL Multiphysics 5.6 FEM solving Navier-Stokes + continuity (2D stationary, shallow-channel approx), plus FSI moving-mesh for time-dependent trajectory of a single rigid fiber (length 47.5 µm, width 5 µm).
- **Imaging**: spinning-disk confocal (Zeiss CSU-X1), 63× oil, live-cell, z-stacks; ImageJ/Fiji directionality plugin for fiber orientation; Imaris for migration tracking.
- **Invasion model**: 2D random walk in Python 3.8.3. Cells start on a sphere; at each step move parallel or perpendicular to local mean fiber orientation (binned to d = 103 µm boxes), four equiprobable moves. Step sizes k_p = r_p/sqrt(N), k_v = r_v/sqrt(N), with r_p, r_v the day-1 downstream/upstream travel distances. Neglects finite cell volume, cell-cell contact, jamming, proliferation, matrix remodeling, and allows cell overlap.

## 4. Key results (quantitative)
- Microfluidic alignment gives a clear fiber-orientation maximum at 90° (main flow direction); non-aligned control is ~uniform (Fig 1b). 33 gels analyzed.
- Spheroid diameter 250–350 µm at embedding; FEM aggregate d = 208 µm (Methods).
- Day-1 invasion front travel: downstream (radial fibers) 222 ± 33 µm vs upstream (tangential) 132 ± 23 µm (Fig 4c).
- Perpendicular (tangential-facing) invasion distance is ~a factor of 3 less than the downstream radial-facing side (Results, Fig 3).
- Downstream/upstream ratio: 2.3× after 1 day, 2.1× after 3 days (Fig 4c).
- Day-3 vs day-1 front distance: only ~1.9× longer in both directions → non-linear (sub-linear) time evolution; consistent with sqrt(time) diffusive growth.
- Best-fit random-walk parameters: parallel step r_p (d_p) = 220 µm/day, perpendicular step r_v (d_v) = 0 µm/day (complete blockage perpendicular) (Fig 5, S5, S10). Sweep tested r_v = 0, 60, 120, 220 µm/day (S10).
- Proliferation estimated to increase spheroid radius only ~1.2× over 3 days (assumed negligible).
- Degree of alignment: fibers within 90° ± 22° (68°–112°) counted as aligned; uniform-distribution baseline = 25.56% of fibers aligned; alignment does not change significantly day 0 → day 1.
- YAP is mostly cytosolic at 3 days (not nuclear/active); periphery vs center nucleus/cytosol YAP ratio differs, p < 10⁻⁴ (S11).
- Simulated single rigid fiber: length 47.5 µm, width 5 µm (FEM rod model).

## 5. Parameters & constants of interest
| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Collagen-I final concentration | 1.85 mg/ml (from 3.32 mg/ml stock) | Methods, 3D collagen matrices |
| Spheroid diameter at embedding | 250–350 µm | Methods, spheroid formation |
| FEM aggregate diameter | 208 µm | Methods, microfluidic chamber |
| Channel dimensions | width 3.8 mm, length 17 mm, height 0.4 mm | Methods, FEM |
| Polymerization flow rate | 0.2 µl/min (fill pulse 90.2 µl/min) | Methods |
| Simulated fiber size | length 47.5 µm, width 5 µm (rigid rod) | Methods, FSI |
| Day-1 downstream invasion front | 222 ± 33 µm | Fig 4c |
| Day-1 upstream invasion front | 132 ± 23 µm | Fig 4c |
| Random-walk parallel step (best fit) | d_p = 220 µm/day | Fig 5 / S5 / S10 |
| Random-walk perpendicular step (best fit) | d_v = 0 µm/day | Fig 5 / S5 / S10 |
| Fiber-orientation box size | d = 103 µm | Methods, simulation |
| Aligned-fiber angular window | 90° ± 22° (68°–112°) | Methods, image analysis |
| Uniform-gel aligned fraction baseline | 25.56% | Methods, image analysis |

Note: these are tissue/ECM-scale phenomenological constants (invasion distances, collagen concentration, random-walk step sizes). None are mechanistic single-cell parameters (no moduli, off-rates, motor forces, persistence lengths). For ffn_cellsim these are context, not transferable runtime constants — hence `has_transferable_params: false`.

## 6. Relevance to ffn_cellsim
Low / tangential, with two narrow uses:
- **Spheroid-scale context for the validation overlay**: This is exactly the regime the project's experimental overlay targets (cancer-cell spheroids in collagen-I; cf. MEMORY: MCF7-spheroid-on-pV4D4/col-I, multi-cell cohesion). It demonstrates that ECM fiber orientation, not just density, sets invasion directionality — relevant when interpreting multi-cell cohesion/invasion overlays later, but well above the single-cell H.1→H.7 chain.
- **Contact-guidance phenomenology / ECM anisotropy reference**: The clean result "migration is effectively blocked perpendicular to fibers; cells move along fibers" is a qualitative anchor for what an emergent ECM/collagen cross-link network (ffn_cellsim's explicit ECM bonds) should reproduce at the single-cell level — but the paper does NOT supply a mechanism (it is a lumped anisotropic random walk, the exact opposite of ffn_cellsim's fine-grained philosophy). It is at best a high-level behavioral acceptance target, not an oracle, and not a mechanism source.
- **Explicitly NOT useful for**: cortex (H.1), FA/clutch (H.5), motor/myosin, integrator/numerics, or any parameter source. The simulation here is a coarse 2D random walk with overlapping point-cells — it is the abstracted/lumped style the project's architectural principle rejects at runtime. No actin, no clutches, no Bell-Evans/Hill, no moduli.

Mapping: themes ECM/collagen + spheroid-scale-context; touches the multi-cell/cohesion EXTEND-era overlay context only, not the active Phase-1 unit chain.

## 7. Limitations & caveats
- Invasion model is deliberately minimal: neglects finite cell volume, cell-cell contacts, jamming, proliferation, matrix remodeling; allows cell overlap. Authors flag this explains the experiment/simulation mismatch at 225°/315° (crowding/bottleneck effects the model misses).
- 2D random walk fit to a 3D process; fiber orientation interpolated from a single representative experimental field (Fig 2a).
- Mechanism is unresolved: authors state YAP is cytosolic (mechanosensing not active at 3 days), cells often show no clear orientation, so contact guidance vs collective force propagation along non-linear collagen vs polarization is left open.
- Single cell line (HeLa), not MCF7; collagen-I only at one concentration.
- Step sizes are empirical fits (k_p, k_v from day-1 fronts), not derived — antithetical to ffn_cellsim's no-magic-number rule if ever ported.
- "Complete blockage perpendicular" (k_v = 0) is a best-fit idealization, not a measured physical constraint.

## 8. Key figures / tables
- **Fig 1**: collagen fiber alignment via microfluidics; orientation distribution (33 aligned gels vs non-aligned control), peak at 90°.
- **Fig 2**: shear-flow characterization around spheroid — aligned fibers (a), bead trajectories (b), PIV velocity magnitude (c), FEM velocity field + single-fiber trajectory (d). Shows tangential fibers upstream, radial downstream.
- **Fig 3**: invasion over 3 days; relative invasion distance of 10 outmost cells per direction (averages of 20 aggregates); directional anisotropy.
- **Fig 4**: upstream vs downstream front — day-1 222±33 vs 132±23 µm (4c); degree of alignment higher downstream (4d).
- **Fig 5**: anisotropic random-walk simulation overlaid on experiment; best fit at k_perp = 0.
- **S10**: perpendicular step-size sweep (r_v = 0, 60, 120, 220 µm/day) — the parameter study selecting k_v = 0.
- **S11**: YAP immunostaining (mostly cytosolic at 3 days; periphery vs center ratio p < 10⁻⁴).

## 9. Notable quotes / citable claims
- "Simulations of the invasive behavior with a Brownian diffusion model suggest complete blockage of migration perpendicularly to fibers allowing for migration exclusively along fibers." (Abstract)
- "the distance of invasion is about a factor of three less than on the downstream facing side with its perpendicularly oriented fibers." (Results, Fig 3)
- "the best match in spheroid shape between experiment and simulation was obtained with a perpendicular step size of 0, i.e. complete blockage of migration perpendicularly to fibers." (Simulation of invasion behavior)
- "the model suggests an increase in travelled distance with the square root of time, which results from the diffusive character." (Simulation of invasion behavior)
- "the collagen orientation displays and preserves the force field applied during polymerization." (Collagen fibers around spheroid)
- "nuclear YAP ... has been shown to be able to trigger invasion ... In our case, however, it is mostly cytosolic and thus not active at the measured time of 3 days." (Conclusion) — supports treating invasion here as mechanical/contact-guidance rather than YAP-mechanotransduction-driven at this timescale.

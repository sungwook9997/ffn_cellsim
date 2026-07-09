# FF ECM library

A grounded, multi-material, 2D/3D, alignment-controlled extracellular-matrix library for the FF
(Cytosim-physics, Warp) engine, with a three-way modulus harness that reads out the effective stiffness
in **real Pascals** and validates every material against the literature. Built 2026-07-10.

## What it delivers

| axis | coverage |
|---|---|
| **materials** | collagen-I, fibrin (fibrillar Mikado) · polyacrylamide, hyaluronic acid, Matrigel, agarose (continuum gel) |
| **dimensionality** | 3D embedded slab · 2D planar sheet |
| **alignment** | nematic order S∈[0,1] about any director (isotropic → fully aligned bundle) |
| **mixing** | interpenetrating composites (e.g. collagen + Matrigel) with inter-network crosslinks |
| **extent** | REV cubes for modulus · native 5R×5R (75×75 µm, ~94 500 nodes) full extent |
| **gradient** | spatially-graded durotaxis substrates (physiological 1 / pathological 10 / interface ~100 Pa/µm) |
| **measurement** | shear G · uniaxial E (+ anisotropy) · spherical indentation E_eff (any x,y) — all in **Pa** |

## Two constitutive classes (the physics, not a convenience)

- **Fibrillar ECM** (collagen, fibrin): a Mikado semiflexible-fiber network. The macroscopic modulus
  **EMERGES** from the microstructure (fiber length density → connectivity ⟨z⟩, bending κ=k_BT·L_p, stretch
  EA, crosslink k_xl). We set the microstructure from literature and *validate* the emergent Pa — never
  tune to the modulus.
- **Continuum gel** (PA, HA, Matrigel, agarose): a flexible/BM gel whose molecular mesh (~nm) is
  unresolvable at cell scale; it is a linear-elastic continuum with a modulus set by chemistry. Modeled
  as a coarse Delaunay spring lattice whose bond stiffness is the analytic inverse of the target E; the
  **indentation returns that E** (validated).

Units: FF µm·pN·s → **1 pN/µm² = 1 Pa** exactly, so all moduli come out directly in SI Pascals.

## Validation — 실제와 같은 Pa (all IN literature band)

| material | class | measured | band (Pa) |
|---|---|---|---|
| collagen-I | fibrillar·G | 12 Pa | 5–100 (Yang-Kaufman G'(1.5mg/mL)≈11) |
| fibrin | fibrillar·G | 63 Pa | 10–1000 (Piechocka) |
| polyacrylamide | continuum·E | 3.4 kPa | 0.1–40 kPa (Tse-Engler) |
| hyaluronic acid | continuum·E | 309 Pa | 10–3000 (brain-mimetic) |
| Matrigel | continuum·E | 411 Pa | 30–900 (Soofi AFM 450) |
| agarose | continuum·E | 14 kPa | 1–100 kPa (Normand) |

**Harness cross-validated** on a known-modulus continuum (PA E=5000 → uniaxial 5000, shear 0.91× the
isotropic relation, indentation 0.68–0.79×). See `ECM_VALIDATION.md`, `figs/ecm_material_moduli.png`.

### Alignment → mechanical anisotropy → tissue

| S (measured) | E∥/E⊥ | tissue |
|---|---|---|
| 0.02 | 1.0 | loose stroma / dermis (healthy) |
| 0.30 | 3.1 | tumor stroma tangential (TACS-2) |
| 0.59 | 10.3 | tumor invasion highway (TACS-3) |
| 0.83 | 63 | tendon / ligament / aligned scar |

`figs/collagen_alignment_anisotropy.png`. The nematic-S sampler (Watson 3D / von-Mises 2D) hits the
target order within a few percent.

## Known limitations (documented, not tuned)

- **Collagen G'(c) scaling** is n≈1.07 vs literature n≈2.0 (Yang-Kaufman): the linear athermal Mikado is
  density-limited; the absolute Pa at reference concentration matches, the steep experimental exponent
  needs thermal/nonlinear semiflexible physics. `figs/collagen_concentration.png`.
- **Fibrillar local indentation ≪ bulk shear** (collagen ~40× decoupling, native-confirmed, thickness-
  independent): a real feature of sub-isostatic networks (van Oosten 2019), not a harness bug (the
  continuum harness returns the input E). The bulk shear/uniaxial modulus is the "material Pa"; indentation
  reads the soft local/compression mode. Continuum gels read the same E all three ways. See `ECM_NOVELTY.md`.
- **Agarose** modeled as a continuum (its dense sub-isostatic bundle network at cell scale is mechanically
  a continuum; the athermal Mikado underestimates its 30-kPa stiffness ~60×).

## Novel contributions (emergent, not tuned)

1. **Emergent nonlinear strain-stiffening** (KB-1.V.2.5 gate): differential modulus K(γ) rises above a
   critical strain γ_c that is ~concentration-independent (geometry-set) — reproduced from a model NOT
   calibrated to it. Peak 3–5× (athermal). `figs/strain_stiffening.png`.
2. **Shear ↔ indentation modulus-decoupling map** across connectivity ⟨z⟩ — the FF library measures both
   modes on the same network, producing the semiflexible mode-decoupling fingerprint directly.
   `figs/mode_decoupling.png`.
3. **Emergent viscoelastic stress relaxation** (KB-1.6): crosslink turnover (Bell-slip) makes each crosslink
   a Maxwell element, so G(t) relaxes with τ ∝ 1/k_off — the relaxation time EMERGES from crosslink lifetime,
   spanning the ECM range (α-actinin ~12 s → weak-physical ~200 s → covalent-LOX ~1170 s). Matrix
   viscoelasticity (Chaudhuri → cell fate). `figs/stress_relaxation.png` (`ff_ecm_viscoelastic.py`).
4. **Durotaxis stiffness-gradient substrate** — spatially-graded local modulus (1/10/100 Pa/µm, KB-1.V.1.3),
   probed by indenting along the gradient axis. `figs/durotaxis_gradient.png` (`ff_ecm_gradient.py`).

## Code

| file | role |
|---|---|
| `ffn_sim/ff/ecm_library.py` | `ECMSpec` registry, nematic-S sampler, 2D/3D fibrillar + continuum builders, composites |
| `ffn_sim/ff/ecm_mechanics.py` | shear / uniaxial / indentation → Pa; strain-stiffening; continuum calibration |
| `ffn_sim/scripts/ff_ecm_validate.py` | per-material validation + concentration/alignment/dim/composite + figures |
| `ffn_sim/scripts/ff_ecm_novelty.py` | emergent strain-stiffening + mode-decoupling |
| `ffn_sim/scripts/ff_ecm_viewer.py` | standalone ECM HTML gallery (fibers+crosslinks) + native npz viewer |
| `ffn_sim/scripts/ff_ecm_native.py` | native 5R×5R full-extent build + GPU indentation (A5000) |

Reproduce:  `python -m ffn_sim.scripts.ff_ecm_validate` · `... ff_ecm_novelty` · `... ff_ecm_viewer`.
Native (gbook): `python -m ffn_sim.scripts.ff_ecm_native --device cuda:0`.

## KB

New DOI-verified material claims (PA/fibrin/Matrigel/agarose/HA) staged for the Notion Contract-Graph in
`KB_ECM_MATERIALS_REGISTRATION.md` (collagen-I and PAA-substrate were already in the KB).

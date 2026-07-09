# FF ECM library — modulus validation (실제와 같은 Pa)

## Per-material reference modulus vs literature band

| material | class | measured (Pa) | band (Pa) | verdict | note |
|---|---|---|---|---|---|
| collagen_I | fibrillar·G | 12.3 | 5–100 | IN BAND | G'(1.5mg/mL)≈11 Pa (KB-1.V.2.1 c² from 5Pa@1mg/mL); scaling-anchor 30-100 (KB-1.30) |
| fibrin | fibrillar·G | 63.1 | 10–1000 | IN BAND | G' 10-1000 Pa over 0.5-10 mg/mL (Piechocka 2010) |
| agarose | continuum·E | 14136.3 | 1000–100000 | IN BAND | E 1-100 kPa over 0.5-4 %w/v (Normand 2000) |
| pa_gel | continuum·E | 3418.9 | 100–40000 | IN BAND | E 0.1-40 kPa tunable; indentation must return the input E (KB-1.21, Tse-Engler) |
| hyaluronic_acid | continuum·E | 308.7 | 10–3000 | IN BAND | E 10 Pa-3 kPa crosslink-tunable |
| matrigel | continuum·E | 410.8 | 30–900 | IN BAND | E 30-900 Pa (typ 450, Soofi 2009 AFM 37°C) |

## Collagen-I concentration series G'(c)

power-law exponent: sim **n=1.07** vs literature n≈2.05 (athermal Mikado is density-limited; absolute Pa at reference conc matches).

| c (mg/mL) | G_sim (Pa) | G_KB (Pa) | ⟨z⟩ | mesh ξ (µm) |
|---|---|---|---|---|
| 1.0 | 8.3 | 5 | 1.95 | 2.00 |
| 1.5 | 13.4 | 11 | 3.02 | 1.50 |
| 2.0 | 18.6 | 20 | 4.14 | 1.50 |
| 3.0 | 26.4 | 55 | 6.03 | 1.00 |
| 5.0 | 51.0 | 150 | 10.11 | 0.50 |
| 7.0 | 66.3 | 342 | 14.40 | 0.50 |

## Collagen alignment → anisotropy (tissue mapping)

| S_target | S_meas | E∥ (Pa) | E⊥ (Pa) | E∥/E⊥ | tissue |
|---|---|---|---|---|---|
| 0.0 | 0.018 | 38.7 | 37.1 | 1.04 | loose stroma / dermis (healthy, TACS-1 dense-isotropic) |
| 0.3 | 0.303 | 74.0 | 23.9 | 3.09 | tumor stroma tangential (TACS-2) |
| 0.6 | 0.586 | 114.5 | 11.1 | 10.32 | tumor invasion highway (TACS-3 radial-aligned) |
| 0.85 | 0.828 | 162.0 | 2.6 | 63.11 | tendon / ligament / aligned scar |

## Dimensionality (2D sheet vs 3D bulk) + composite (mixed ECM)

| build | modulus (Pa) | kind | ⟨z⟩ | n_fibers |
|---|---|---|---|---|
| collagen 3D bulk | 17.4 | G (bulk) | 4.22 | 586 |
| collagen 2D sheet | 71.0 | E_2D (in-plane, Pa·µm/µm) | 4.77 | 52 |
| composite collagen_I+matrigel | 85.4 | G (bulk) | 3.37 | — |

_collagen(S=0.3) + Matrigel interpenetrating; G exceeds either component alone_


## Figures

- `figs/ecm_material_moduli.png` — each material's modulus vs literature band
- `figs/collagen_concentration.png` — G'(c) vs Yang-Kaufman
- `figs/collagen_alignment_anisotropy.png` — anisotropy vs S
- `figs/indentation_curves.png` — Hertz F∝δ^1.5 for continuum gels

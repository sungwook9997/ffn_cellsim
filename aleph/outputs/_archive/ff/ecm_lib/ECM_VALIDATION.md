# FF ECM library — modulus validation (실제와 같은 Pa)

## Per-material reference modulus vs literature band

| material | class | measured (Pa) | band (Pa) | verdict | note |
|---|---|---|---|---|---|
| collagen_I | fibrillar·G | 12.8 | 5–100 | IN BAND | G'(1.5mg/mL)≈11 Pa (KB-1.V.2.1 c² from 5Pa@1mg/mL); scaling-anchor 30-100 (KB-1.30) |
| pa_gel | continuum·E | 3555.6 | 100–40000 | IN BAND | E 0.1-40 kPa tunable; indentation must return the input E (KB-1.21, Tse-Engler) |

## Collagen-I concentration series G'(c)

power-law exponent: sim **n=1.08** vs literature n≈2.05 (athermal Mikado is density-limited; absolute Pa at reference conc matches).

| c (mg/mL) | G_sim (Pa) | G_KB (Pa) | ⟨z⟩ | mesh ξ (µm) |
|---|---|---|---|---|
| 1.0 | 8.6 | 5 | 1.95 | 2.00 |
| 1.5 | 13.9 | 11 | 3.02 | 1.50 |
| 2.0 | 19.2 | 20 | 4.14 | 1.50 |
| 3.0 | 27.3 | 55 | 6.03 | 1.00 |
| 5.0 | 52.7 | 150 | 10.11 | 0.50 |
| 7.0 | 69.3 | 342 | 14.40 | 0.50 |

## Collagen alignment → anisotropy (tissue mapping)

| S_target | S_meas | E∥ (Pa) | E⊥ (Pa) | E∥/E⊥ | tissue |
|---|---|---|---|---|---|
| 0.0 | 0.018 | 42.3 | 40.5 | 1.04 | loose stroma / dermis (healthy, TACS-1 dense-isotropic) |
| 0.3 | 0.303 | 81.6 | 26.1 | 3.12 | tumor stroma tangential (TACS-2) |
| 0.6 | 0.586 | 127.3 | 11.7 | 10.92 | tumor invasion highway (TACS-3 radial-aligned) |
| 0.85 | 0.828 | 180.6 | 2.6 | 69.22 | tendon / ligament / aligned scar |

## Dimensionality (2D sheet vs 3D bulk) + composite (mixed ECM)

| build | modulus (Pa) | kind | ⟨z⟩ | n_fibers |
|---|---|---|---|---|
| collagen 3D bulk | 17.9 | G (bulk) | 4.22 | 586 |
| collagen 2D sheet | 86.9 | E_2D (in-plane, Pa·µm/µm) | 4.77 | 52 |
| composite collagen_I+matrigel | 89.1 | G (bulk) | 3.37 | — |

_collagen(S=0.3) + Matrigel interpenetrating; G exceeds either component alone_


## Figures

- `figs/ecm_material_moduli.png` — each material's modulus vs literature band
- `figs/collagen_concentration.png` — G'(c) vs Yang-Kaufman
- `figs/collagen_alignment_anisotropy.png` — anisotropy vs S
- `figs/indentation_curves.png` — Hertz F∝δ^1.5 for continuum gels

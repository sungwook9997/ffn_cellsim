# FF ECM library — emergent & mode-dependent mechanics (novel)

## 1. Emergent nonlinear strain-stiffening (NOT tuned; KB-1.V.2.5 gate)

| material | conc | ⟨z⟩ | K₀ (Pa) | peak K/K₀ | γ_c |
|---|---|---|---|---|---|
| collagen_I | 1.5 | 3.23 | 5.5 | 3.4× | 0.30 |
| collagen_I | 3.0 | 6.54 | 12.6 | 4.8× | 0.22 |
| fibrin | 2.0 | 3.65 | 31.7 | 2.7× | nan |

Emergent K(γ) rises above γ_c with γ_c ~concentration-independent (geometry-set) — the KB-1.V.2.5 qualitative signature, reproduced from a model NOT calibrated to it. Peak stiffening is moderate (athermal Mikado); the full K>10× regime needs deeper strain / thermal nonlinearity.

## 2. Shear ↔ indentation modulus decoupling (van Oosten 2019)

| conc (mg/mL) | ⟨z⟩ | G_shear (Pa) | bulk E (Pa) | local E_indent (Pa) | decoupling |
|---|---|---|---|---|---|
| 1.0 | 1.85 | 7.1 | 16.9 | 0.48 | 35× |
| 2.0 | 4.35 | 16.4 | 39.4 | 0.69 | 57× |
| 4.0 | 8.31 | 31.2 | 74.8 | 1.13 | 66× |

Sub-isostatic collagen is stiff in bulk shear but soft under LOCAL indentation — the semiflexible mode-decoupling. The FF library measures both on the same network, so the decoupling ratio is produced directly as a function of architecture (⟨z⟩, alignment).

## Figures
- `figs/strain_stiffening.png` — K(γ)/K₀ + K∝σ
- `figs/mode_decoupling.png` — bulk vs local modulus vs ⟨z⟩

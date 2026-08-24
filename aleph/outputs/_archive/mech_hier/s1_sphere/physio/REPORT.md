# S.1 — homogeneous sphere, distributed load — REPORT

**Run class:** NATIVE (authoritative)  ·  N_fil=70686, device=cuda:0, steps=4000, membrane=False
**Loading:** {'pressure_setpoint': None, 'Lp_um_s_Pa': 1e-07, 'K_drained_Pa': 300.0, 'load_time_s': 300.0, 'xl_koff_per_s': 0.066, 'eta_bulk_Pa_s': 65.9, 'rigid_plate': True}
**R0:** 7.500 um

## Hertz gate (forward oracle)

- E_fit = 4140.7 Pa (E* = 5520.9 Pa), band [224.0, 279.0] Pa → **TOO STIFF**
- fit r^2 = 0.9656, pointwise-E flatness (CV) = 0.348 (0 = ideal Hertzian), small-strain valid = True
- forward max residual = 52.69%, n_points = 4

## Force–displacement

| strain | delta1 [um] | F [pN] | dP [Pa] | V/V0 | svm_p95 [Pa] |
|---|---|---|---|---|---|
| 0.5% | 0.037 | 309.4 | 39 | 1.000 | 213.3 |
| 1.3% | 0.100 | 813.1 | 39 | 1.000 | 258.5 |
| 2.2% | 0.162 | 1288.5 | 39 | 0.999 | 304.7 |
| 3.0% | 0.225 | 2108.1 | 39 | 0.999 | 358.6 |

## Figures

- `figs/s1_force_displacement.png` — F vs delta1 with the fitted Hertz oracle overlaid.

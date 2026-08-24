# S.1 — homogeneous sphere, distributed load — REPORT

**Run class:** NATIVE (authoritative)  ·  N_fil=70686, device=cuda:0, steps=15000, membrane=False
**Loading:** {'pressure_setpoint': None, 'Lp_um_s_Pa': 1e-07, 'K_drained_Pa': 300.0, 'load_time_s': 300.0, 'xl_koff_per_s': 0.066, 'eta_bulk_Pa_s': 65.9, 'rigid_plate': True}
**R0:** 7.500 um

## S1 analytical gate (constitutive self-consistency: Hertz-solid vs pressurized-shell)

- **LINEAR-SHELL (clean)** — better-fitting law = linear-shell
- solid-Hertz (F~δ^1.5): r^2 = 0.9539, pointwise-E flatness = 0.443, residual = 59.6%
- linear-shell (F~δ, Reissner/tension): r^2 = 0.9916, residual = 8.8%, k = 7887 pN/µm
- Hertz E_fit = 3623.6 Pa (E* = 4831.5 Pa), small-strain valid = True, n_points = 3
- **INFORMATIONAL (not S1 pass/fail):** E_fit x14.6 the MCF7 WHOLE-cell band [224.0, 279.0] Pa (bare cortex is a sub-component; MCF7 is a downstream layer-assembled target).

## Force–displacement

| strain | delta1 [um] | F [pN] | dP [Pa] | V/V0 | svm_p95 [Pa] |
|---|---|---|---|---|---|
| 0.5% | 0.037 | 317.3 | 31 | 1.000 | 272.6 |
| 1.7% | 0.131 | 951.5 | 36 | 1.000 | 364.0 |
| 3.0% | 0.225 | 1819.8 | 36 | 0.999 | 467.5 |

## Figures

- `figs/s1_force_displacement.png` — F vs delta1 with the fitted Hertz oracle overlaid.

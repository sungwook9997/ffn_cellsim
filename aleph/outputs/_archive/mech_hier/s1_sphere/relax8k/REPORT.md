# S.1 — homogeneous sphere, distributed load — REPORT

**Run class:** NATIVE (authoritative)  ·  N_fil=70686, device=cuda:0, steps=8000, membrane=False
**Loading:** {'pressure_setpoint': None, 'Lp_um_s_Pa': 1e-07, 'K_drained_Pa': 300.0, 'load_time_s': 300.0, 'xl_koff_per_s': 0.066, 'eta_bulk_Pa_s': 65.9, 'rigid_plate': True}
**R0:** 7.500 um

## S1 analytical gate (forward Hertz oracle self-consistency)

- **NON-HERTZIAN** — is the bare body a clean Hertzian elastic solid?
- fit r^2 = 0.9287 (>=0.99), pointwise-E flatness CV = 0.414 (<=0.10), forward residual = 59.1% (<=10%), small-strain valid = True
- E_fit = 3688.1 Pa (E* = 4917.4 Pa), n_points = 4
- **INFORMATIONAL (not the S1 pass/fail):** E_fit is x14.8 the MCF7 WHOLE-cell band [224.0, 279.0] Pa — the bare cortex is a sub-component, not the whole cell; the MCF7 match is a downstream (layer-assembled) target.

## Force–displacement

| strain | delta1 [um] | F [pN] | dP [Pa] | V/V0 | svm_p95 [Pa] |
|---|---|---|---|---|---|
| 0.5% | 0.037 | 318.5 | 38 | 1.000 | 248.0 |
| 1.3% | 0.100 | 787.0 | 38 | 1.000 | 300.3 |
| 2.2% | 0.162 | 1168.4 | 38 | 0.999 | 357.1 |
| 3.0% | 0.225 | 1843.4 | 38 | 0.999 | 422.4 |

## Figures

- `figs/s1_force_displacement.png` — F vs delta1 with the fitted Hertz oracle overlaid.

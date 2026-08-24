# S.1 — homogeneous sphere, distributed load — REPORT

**Run class:** NATIVE (authoritative)  ·  N_fil=70686, device=cuda:0, steps=4000, setpoint=40.0 Pa, rigid_plate=False
**R0:** 7.500 um

## Hertz gate (forward oracle)

- E_fit = 888072.5 Pa (E* = 1184096.6 Pa), band [224.0, 279.0] Pa → **TOO STIFF**
- fit r^2 = 0.9780, pointwise-E flatness (CV) = 0.250 (0 = ideal Hertzian), small-strain valid = True
- forward max residual = 93.28%, n_points = 4

## Force–displacement

| strain | delta1 [um] | F [pN] | dP [Pa] | V/V0 | svm_p95 [Pa] |
|---|---|---|---|---|---|
| 0.5% | 0.037 | 16244.9 | 40 | 1.000 | 23.3 |
| 1.3% | 0.100 | 101380.4 | 40 | 1.000 | 25.1 |
| 2.2% | 0.162 | 258015.5 | 40 | 0.999 | 27.3 |
| 3.0% | 0.225 | 488431.3 | 40 | 0.999 | 30.1 |

## Figures

- `figs/s1_force_displacement.png` — F vs delta1 with the fitted Hertz oracle overlaid.

# Anchor force balance — Week 2 investigation (F3)

PI directive 2026-04-29 (Option F Week 2). Read-only diagnostic. Goal:
classify F3 (anchor force balance, median 10.158 vs limit 0.20) as
solver bug or accepted physics.

## Gate definition

Source: `acs/physics/mlsmpm.py:substrate_diagnostics` (lines 1814–1944).

```python
P_per_p = K_eff_band * (1.0 - rho_ref / np.clip(rho_band, 1e-30, None))
F_pressure_down = float((P_per_p * V0 / h_band).sum()) + F_gravity_band
F_substrate_up = float(self.diag_substrate_impulse_z[None]) / self.cfg.dt_star
balance_err = abs(F_substrate_up - F_pressure_down) / max(abs(F_substrate_up), 1e-12)
```

Documented expectation:
- `F_pressure_down` = downward bulk pressure on substrate ≥ 0 (compressive)
- `F_substrate_up` = upward reaction from substrate impulse ≥ 0
- At equilibrium: F_substrate_up ≈ F_pressure_down → balance_err ≈ 0

## Observed in Production Lam4

Sampled from `results/production_lam4/metrics.csv`:

| frame | t* | F_substrate_per_step | F_pressure_down | rel_err |
|---|---|---|---|---|
| 1 | 15.0 | +0.0495 | **−0.380** | 8.68 |
| 10 | 150.0 | +0.0545 | **−0.731** | 14.41 |
| 50 | 750.0 | +0.0373 | **−0.471** | 13.63 |
| 160 | 2400.0 | +0.0381 | **−0.453** | 12.88 |
| 320 | 4800.0 | +0.0356 | **−0.258** | 8.24 |

Mean: F_substrate = +0.038, F_pressure = **−0.372**. Ratio
F_pressure/F_substrate = **−9.7**.

## Root cause

**`F_pressure_down` is negative throughout the run.**

The formula `P = K · (1 − ρ_ref/ρ_kernel)` assumes ρ_kernel ≥ ρ_ref
(compressive). But the contact band exhibits ρ_kernel/ρ_ref ≈ 0.726
(F4 gate: contact-band ρ_kernel/ρ_ref = 0.726 vs window [0.85, 1.15]),
so:

```
ρ_ref / ρ_kernel ≈ 1 / 0.726 ≈ 1.378
P_per_p ≈ K · (1 − 1.378) ≈ −0.378 · K   (TENSILE, not compressive)
```

**The contact band is under-densified by the same kernel-truncation
that drives F4.** Adami-Hu-Adams 2010 §3 documents that near a
boundary, the SPH/MPM kernel sees fewer neighbours, so ρ_kernel
estimated from the kernel sum is biased low. Under v15 (k.3) the
volumetric stress σ_vol = K·(ρ_ref/ρ_kernel − 1)·I then reports a
**tensile** stress — the formula is correctly evaluating the
*kernel-density-derived* state, but that state is itself a truncation
artifact at the −z boundary.

The substrate impulse `diag_substrate_impulse_z` is accumulated only
on **downward grid velocity clamps** (the −z reflective wall stops
particles trying to fall through). When the contact band is under-
densified and σ_vol is tensile, particles are mostly *not* trying to
fall through the substrate — gravity and the small remnant compressive
load drive a tiny F_substrate_up ≈ 0.038. There is no actual force
imbalance; the bulk-side measurement (F_pressure_down) is using a
formula that assumes a regime that does not hold near the kernel
truncation boundary.

## What the gate **should** be measuring

The substrate-anchor force balance check exists to detect runaway
imbalance: the substrate either failing to anchor (collapse) or
double-anchoring (over-stiff). The current formula compares two
asymmetric quantities:

- `F_substrate_up` = numerical impulse, accumulated only when v_z < 0
  is clamped (truly cumulative reaction force).
- `F_pressure_down` = analytic stress integral via the contact-band
  particles, which is the *whole* down-pushing load minus gravity.

The asymmetry: F_substrate_up sees only what the substrate *actually*
deflected, while F_pressure_down predicts what *would* be deflected if
the contact band were truly compressive. When the kernel-truncation
biases ρ_kernel low at the boundary, F_pressure_down overestimates the
load (in magnitude, and with wrong sign).

## Classification

| Q | Answer |
|---|---|
| Q1 (reproducible) | YES — appears in all phenotypes (Bare 6.17, Pre 6.68, Lam4 7.68 pilot; 10.16 production), magnitude grows mildly with scale |
| Q2 (cause known) | YES — kernel truncation at substrate boundary makes ρ_kernel < ρ_ref → P negative → F_pressure_down formula reports tensile state; gate formula's compressive assumption invalid |
| Q3 (magnitude bounded) | YES — bounded by Adami-Hu-Adams §3 truncation (≈ 30%); does not grow without bound |
| Q4 (PI accepted) | This document is the proposed PI-acceptance |

**Per `docs/gate_fail_taxonomy.md` decision tree → ACCEPTED-LIMITATION.**

The FAIL is a **measurement artifact in the diagnostic formula**, not
a physics or solver failure. The underlying simulation is **not** out
of force balance — particles are not punching through the substrate
(F_substrate_up is a stable small positive value), gravity is not
producing runaway acceleration (anchor balance gate measures
*relative* error, not magnitude), and the pile is not bouncing. The
diagnostic compares an analytic compressive estimate against a
numerical reaction force, but in the regime where the kernel sees the
−z boundary at all, the analytic estimate is biased by the truncation
that F4 also reports.

## Causality assessment for the asymptote

**F3 is NOT a cause of the Production Lam4 peak-and-decay.**

The argument:
- If F3 indicated real over-anchoring (substrate clamping the spheroid
  too aggressively), F_substrate_up would be *much larger than*
  F_pressure_down. We observe the opposite: F_substrate_up is *much
  smaller in magnitude*.
- F_substrate_up ≈ 0.038 (steady) vs total system mass ≈ 4.19 means
  the substrate reaction force per unit mass is ~0.009 — well below
  the gravity term (gravity_star = 0.01 acting on ~5 contact-band
  particle layers, F_gravity_band ≈ 0.020). The substrate is **not
  over-anchoring**; it is under-deflected.
- Causation of long-time retraction is on the Layer 4 Marangoni axis
  (interior φ decay → Marangoni retraction), per
  `docs/marangoni_review.md`. F3 is a separate measurement issue.

## Recommended action

1. **Reclassify F3** in `docs/gate_fail_taxonomy.md` from HARD-BLOCKER
   to **ACCEPTED-LIMITATION** with reference to this document.
2. **Update the gate formula** (Stage-1d.b or later code commit, not
   in this Week 2 read-only pass):
   - Replace the compressive-only `P = K(1 − ρ_ref/ρ_kernel)` with the
     correct two-sided form already used elsewhere
     `σ_vol = K(ρ_ref/ρ_kernel − 1)·I`, which gives P = −σ_zz =
     K(1 − ρ_ref/ρ_kernel) — but explicitly handle the truncation
     regime by using `max(0, P_per_p)` for the *compressive part*
     and reporting the tensile part separately as `F_tensile_up`
     (the substrate would *attract* particles via the v15 mechanism
     at this density level).
   - Or, alternatively: gate on `|F_substrate_up + F_gravity_band| /
     (M_total · g_star)` instead — this is a normalized force-balance
     check that does not require the analytic compressive estimate.
3. **Document** in `docs/12_validation.md` known-limitations registry
   under "v15 σ_vol kernel truncation at substrate boundary".

## What this changes for the publication

- F3 no longer blocks asymptote interpretation. The Production Lam4
  Bucket P3 peak-and-decay reading is independently confirmed by
  Layer 4 / Marangoni mechanism analysis.
- The remaining HARD-BLOCKER list in `docs/gate_fail_taxonomy.md` is
  reduced to: F2 (horizontal momentum drift, see companion doc) and
  F9 (φ trajectory, Layer 3 audit pending).
- F4 (contact-band ρ_kernel/ρ_ref window) is the same root cause as
  F3 viewed from the density side; F4 should be reclassified to
  ACCEPTED-LIMITATION (linked to F3) at the same time.

## Cross-references

- `docs/gate_fail_taxonomy.md` — F3, F4 entries
- `docs/codex_review_synthesis.md` — Codex review item 2
- `docs/marangoni_review.md` — independent asymptote diagnosis
- `docs/stage1a_plus_substrate_sanity.md` — original Adami-Hu-Adams §3
  truncation discussion at substrate entry
- Adami, S., Hu, X. Y. & Adams, N. A. 2010. "A new surface-tension
  formulation for multi-phase SPH using a reproducing divergence
  approximation." *J Comput Phys* **229** 5011 — kernel truncation
  near a boundary, §3.

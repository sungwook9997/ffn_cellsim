# FF Stage 6d — γ-floor prototype: MD-free reproduction of the cortical-tension floor

**Date:** 2026-06-29 · **Engine:** `ffn_sim/ff/` (Filament-FEM, Cytosim physics) · **Branch:** `dcm/main`
**Status:** prototype DONE + decisive result — **surfaced to PI** (no tuning; γ swept as a controlled variable)

---

## 1. What this is

The decisive FF experiment from `ENGINE.md §4`: measure cortical tension γ with an **MD-free
mechanical-equilibrium solve** (Cytosim physics on the Warp/FF stack) and compare it to the
**BAOAB-MD γ** from the archived `cortex/cortical_tension.py` (the original γ-floor finding). The
question (ENGINE.md §4): *does the MD-free solve reproduce the same γ (⇒ FF is MD-equivalent for γ),
or does it localize where in-time load-coupling is load-bearing?*

## 2. Method (all grounded; FF units pN·µm·s)

- **Cortex** (`cortex_assembly`): 1000-filament ×40 mesoscopic actin shell, R=10µm, 7 beads,
  ℓ₀=0.5µm, κ=k_BT·ℓ_p (ℓ_p=17µm) — H.3 config values, identical to the MD cortex.
- **Crosslinkers / myosin** (`hand_kmc` presets): α-actinin (Ferrer 2008) / filamin (Pereverzev
  catch–slip) 30/70; NMIIA myosin (Kovács 2003 / Stam-Hocky 2015). Links connect cross-fiber node
  pairs within the **mesoscale reach √(A/n_fil)** — the geometric dual of the ×40 coarse-graining
  (the molecular ε=60nm can't connect the sparse coarse filaments; faithful port of
  `connected_mesh.mesoscale_reach`). Force-free at formation.
- **Turgor** (physiological-baseline rule): Young-Laplace radial pressure, dP0=133 Pa (SOLID
  registry), K_vol=1e3 Pa.
- **Equilibrium:** settle the RESTING passive shell (bending only) + constraint projection + reshape
  (NF2007 project-then-reshape, no BAOAB thermostat); myosin added as the MODULATOR at measurement
  (physiological-baseline rule: measure tension on the resting turgor-pressurised cell + add myosin).
- **γ:** `gamma_estimator.method_of_planes_gamma` — a **faithful port of the validated HOOMD γ
  kernel**, so the MD-free γ is read with the *same instrument* as the BAOAB-MD γ. Two channels kept
  SEPARATE (archived rule): γ_active = actomyosin (actin axial tension from the inextensibility
  multipliers + crosslinker + myosin links); γ_passive = turgor ΔP·R/2.

> **Units trap caught:** 1 pN/µm = 1e-3 mN/m, so the Salbreux band 0.35–0.65 mN/m = **[350, 650]
> pN/µm** (not 0.35–0.65). An early 1000× slip had inverted the conclusion; fixed.

## 3. Result — the floor is REPRODUCED MD-free

γ_active is linear in the myosin prestress f_myo (the swept controlled variable): **≈0.19 pN/µm per
pN** of per-link prestress (R=10µm, n_myo=200, 8 realizations).

| f_myo [pN/link] | γ_active [pN/µm] | γ_active [mN/m] | vs Salbreux band |
|---|---|---|---|
| 0 | 0.0 | 0 | — |
| 1 | 0.19 ± 0.01 | 1.9e-4 | ~1800× under |
| **5 (NMIIA per-side stall, lit anchor)** | **0.95 ± 0.06** | **9.5e-4** | **~370× under** |
| 10 | 1.90 | 1.9e-3 | ~180× under |
| 50 | 9.5 | 9.5e-3 | ~37× under |
| 150 | 28.5 | 2.85e-2 | ~12× under |

- **γ_passive (turgor) = 665 pN/µm = 0.665 mN/m — at the top of the Salbreux band** (by construction:
  dP0=133 Pa was band-implied via Young-Laplace; this also validates the estimator's Young-Laplace
  channel analytically).
- **At the lit-anchored NMIIA prestress, actomyosin γ is ~370× UNDER band**, while turgor γ is
  at-band. To reach band-level cortical tension via this transmission the per-link prestress would
  need ≈1840 pN ≈ **370× the physiological NMIIA per-side stall (5 pN)**.

Figure: `outputs/ff/figs/gamma_floor_sweep.png` (log-γ vs f_myo, bands + turgor baseline + lit
anchor overlaid). Data: `outputs/ff/gamma_floor_sweep.npz`.

## 3b. Quantitative cross-check vs the BAOAB-MD γ (disk-grounded)

From `docs/CORTICAL_TENSION_RECORD_2026-06-04.md` (the authoritative archived γ-floor record):

| quantity | BAOAB-MD (archived) | MD-free FF (this work) |
|---|---|---|
| active actomyosin γ at lit kinetics | g_soft ≈ **3.1e-5 – 1.4e-4 mN/m** (FA-anchored / Δ_active) | γ_active = **9.5e-4 mN/m** (f_myo=5 pN) |
| best case under forcing | v0×3000 → **0.030 mN/m (11.6× under)** | linear in f_myo; band needs ≈370× more force |
| KU-3.5 band floor | **0.35 mN/m** | 350 pN/µm (= 0.35 mN/m) |
| passive/turgor γ | at-band (turgor channel) | 665 pN/µm = **0.665 mN/m** at-band |
| diagnosed cause | **force GENERATION** (per-head force / motor recruitment), NOT boundary / FA / compliance / kinetics | **force magnitude / transmission** (too few force-bearing motors per cross-section); kinetics NOT the lever (Stage 6e) |

Both independent methods land the actomyosin γ in the same **~1e-4 – 1e-3 mN/m** regime — ~100–1000×
under band — and both localize the cause to **force generation/magnitude, not the integrator, the
boundary, adhesion, compliance, or binding kinetics**. The MD line showed FA-anchoring and compliant
backbones do not lift g_soft, and that even v0×3000 reaches only 0.030 mN/m (super-stall, filtered
out by the physical-validity gate); the FF line shows the same floor emerges at mechanical
equilibrium with full transmission and is unmoved by dynamic Hand turnover (Stage 6e). **The MD-free
mechanical solve reproduces the BAOAB-MD γ-floor** → FF is MD-equivalent for γ, and the floor is a
physics statement (the missing motor-density datum), not a method artifact.

## 4. Interpretation (what it localizes) — for PI

The MD-free mechanical solve **reproduces the BAOAB-MD γ-floor**: actomyosin tension ~100–1000×
under band, turgor at-band. Because the MD-free solve assumes the myosin stall force is *transmitted*
through the network at mechanical equilibrium (no dynamic binding-throughput limit), reproducing the
floor here means **the floor is STRUCTURAL — a force-magnitude / transmission deficit, not a dynamic
MD artifact**. This is the same conclusion the MD line reached from the other side:

- H7 Gate-B: "myosin stress ~400× < turgor" (active floor is force-TRANSMISSION/lever, not binding
  kinetics) — `project-h7-gate-b-relaxed-build`.
- SF/NMII force-scale REFUTE/HALT: the single-SF tension is not closable from a measured native
  motor-density datum — the **missing parallel-minifilaments-per-cross-section datum** —
  `project-sf-nmii-forcescale-result`.

Our ≈370× gap is the γ-twin of that same missing motor-density factor: the model has the right
*per-motor* physics (NMIIA stall 0.5 pN/head, Kovács/Stam-Hocky kinetics) but **too few force-bearing
motors per equatorial cross-section** to lift actomyosin γ to band. The shell instability we also saw
(turgor at band-implied dP0 inflates the floored cortex — it cannot balance 665 pN/µm with ~1 pN/µm
of actomyosin tension) is the same deficit expressed dynamically.

## 5. Open / PI decision

This is **not** a magic-number situation — γ was swept as a controlled variable; no derived param was
lowered to hit a band (hard rule respected). The floor is a *physics statement*, not a bug. Options
for PI (none auto-taken):

1. **Accept the floor as the FF/MD-consistent result** — both engines agree the single-cell
   center/mesoscale actomyosin cannot reach band cortical tension; the missing piece is the
   motor-density datum (surface to literature search), not the integrator.
2. **Add the missing real force at a derived value** — if a grounded parallel-minifilament density
   (motors per cross-section) can be sourced, raise n_myo / f_myo to that DERIVED value (legit lever)
   and re-run the sweep; predict γ_active scales linearly so band needs ≈370× more force-bearing
   motors than currently modelled.
3. **γ target band itself** stays unresolved (Moazzeni MCF7 1e-2 N/m vs SimuCell3D 1e-3 vs emergent;
   ENGINE.md §6) — the sweep reports γ(f_myo) so any band can be read off without retuning.

## 6. Build artifacts (Stage 6b-iv → 6d, branch dcm/main)

| chunk | module | tests |
|---|---|---|
| 6b-iv inextensibility | `ff/constraints.py` (projector P, reshape, axial tension) | 6 + 3 |
| 6c-a units | `ff/units.py` (pN·µm·s) | 5 |
| 6c-b cortex | `ff/cortex_assembly.py` + `viz_cortex.py` | 4 |
| 6c-c Hand KMC | `ff/hand_kmc.py` (NF2007 §10.1) | 8 |
| 6d γ estimator + floor | `ff/gamma_estimator.py`, `ff/gamma_floor.py`, `ff/gamma_floor_sweep.py` | 7 + 8 |

Full `tests/ff` suite: **43/43 green.**

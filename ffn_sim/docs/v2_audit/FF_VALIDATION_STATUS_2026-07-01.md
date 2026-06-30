# FF engine — validation status matrix (consolidated)

**Date:** 2026-07-01  **Engine:** `ffn_sim/ff/` (Filament-FEM, Cytosim physics, Warp GPU-native)
**Branch:** dcm/main  **Suite:** `tests/ff/` — **102 tests / 20 files** (4 GPU-parity run on the A5000,
3 Cytosim-parity gated on the external binary; all others always-on).

A single map of every FF mechanism and how it is validated, per the project rule that FF is validated
against **analytic ground truth as the PRIMARY always-on test**, with the runnable Cytosim oracle as an
on-demand cross-check ([[feedback-oracle-is-crosscheck-not-truth]]). Consolidates the FF_STAGE6* docs.

## Validation matrix

| mechanism | analytic ground truth | always-on? | oracle (Cytosim) | GPU↔CPU parity | test file |
|---|---|---|---|---|---|
| Bending FORCE | −κ q⁴ y dispersion (q·seg≪1) | ✅ | — | (via relax) | test_bending_dispersion |
| Bending ENERGY magnitude | **κL/2R²** (arc) | ✅ **NEW (6Q+)** | ✅ (skip offline) | — | test_bending_energy_analytic |
| Bending = −∇E | central-FD gradient | ✅ | — | — | test_cytosim_bending |
| Relaxation RATE | (κ/ζ)·q⁴, q⁴-scaling | ✅ | next target | — | test_relax |
| Euler buckling | π²κ/L² (qualitative only) | ⚠️ partial | next target | — | test_buckling |
| Inextensibility | projector P=I−Jᵀ(JJᵀ)⁻¹J + reshape | ✅ | — | reshape 2e-16 | test_constraints, test_gpu_parity |
| Crosslink/myosin springs | Hookean f=k·Δ | ✅ | — | link 2e-11 / myo 2e-15 | test_gamma_floor, test_gpu_parity |
| Arp2/3 branch | angle-harmonic θ₀=70° | ✅ | — | branch 2e-16 | test_weave, test_network_warp |
| Hand KMC turnover | Bell steady state k_on/(k_on+p_off) | ✅ | next target | statistical | test_hand_kmc, test_network_warp |
| Turgor closure | Guo ΔP(V); **K_vol=Π_in0/(1−vmin)** DERIVED | ✅ **NEW** | — | reductions 1e-13 | test_turgor_closure |
| Young-Laplace | γ_passive = ΔP·R/2 | ✅ | — | — | test_gamma_floor |
| Mobility/drag | μ=log(L_h/δ)/(3πηL); γ=1/(8μ) | ⚠️ sane-only | — | — | test_units |
| Method-of-planes γ | HOOMD-ported instrument | ✅ | (BAOAB-MD) | — | test_gamma_estimator |
| **γ_myo floor metric** | linear, zero-intercept | ✅ **NEW** | — | — | test_gamma_floor |
| Shear modulus G | affine f_na·k_xl·ρ_L·ℓc | ✅ | (Kim/MacKintosh) | — | test_kim_network |
| Units (nondim pN·µm·s) | SI round-trip | ✅ | — | — | test_units |
| **Full integrators** | — | — | — | relax 1e-14 / loaded-shell 4e-14 | test_gpu_parity |

## What landed this sprint (FF_STAGE6Q–6R, 2026-07-01)

1. **Lit-faithful densities + clean γ_myo metric (6Q):** actin 100/µm² (KB-3.18), NMIIA 0.625/µm² (Nie),
   floor on the clean γ_myo channel (γ_active is passive-residual-contaminated at native density,
   adversarially verified density-linear). Floor ~1700× vs MCF7-active — DEEPER but honest. Confirms the
   force-magnitude γ-floor. See FF_STAGE6Q.
2. **κL/2R² always-on bending-energy magnitude anchor** — removed the offline dependency on the Cytosim
   binary for this anchor.
3. **Comprehensive GPU↔CPU bit-parity (6R)** — every kernel + integrator at native scale, codified
   skip-if-offline. The A5000 production path is bit-faithful. See FF_STAGE6R.
4. **Turgor Guo-closure + no-magic-K_vol anchor** — the bulk modulus is DERIVED, supporting the
   no-magic-number hard rule.

## Production-envelope robustness (A5000 stress sweep, 2026-07-01)

The engine stays FINITE + stable across the full production envelope (no NaN / inf / blow-up):

| sweep | range | result |
|---|---|---|
| filament count N | 1 000 → 100 000 | all finite; γ_myo linear & zero-intercept at every N |
| crosslink stiffness k_xl | **0.1 → 4.6e5 pN/µm (Ferrer, 4.6e6×)** | all stable; R_mean holds at rest 7.499 µm |
| myosin prestress f_myo | 0 → 40 pN | γ_myo finite & exactly linear (γ_myo = c·f_myo) |
| loaded-shell at native N=70686 | 6000 steps | V/V0 = 1.00001, no blow-up (CFL incl. turgor breathing-mode) |

**Stiff-crosslink finding:** `equilibrate(crosslink_turnover=True)` relaxes the lit-anchored stiff Ferrer
`k_xl = 4.6e5 pN/µm` STABLY (R_mean held, no blow-up) — so the "robust stiff-crosslink relaxer" the Kim
doc (FF_KIM_NETWORK_VALIDATION §mechanics-magnitude) flagged as pending DOES exist (it is the force-free-
rebinding turnover baseline, Stage 6N). Promoting the production `link_k` 0.1 → 4.6e5 remains PI-gated
(it feeds the CFL + the γ-floor cortex, and the Ferrer datum needs KB registration), but the relaxer is no
longer the blocker.

## Remaining gaps (PI-roadmap or future)

- **Euler buckling π²κ/L²** — only qualitative (amplitude grows with compression); a quantitative
  threshold anchor needs the finite-EA extensible mode (the inextensible default cannot be axially
  loaded). Low risk, deferred.
- **Mobility absolute** — only "sane positive" is checked, not the NF2007 μ formula value (trivial,
  low-value).
- **Cytosim dynamic parity** (relaxation rate / Euler / Hand kinetics) — blocked offline (no `sim`
  binary on the current host; the static κL/2R² parity passed when last run). Rebuild PI-gated.
- **Passive density-linear self-stress** — γ_actin(f0) ≈ 0.0085·areal_density pN/µm is a characterized
  PHYSICAL prestress of the crosslinked-inextensible network on the curved shell (not a bug; verified
  density-linear). Whether to reduce it via a different cortex construction is a modeling-fidelity PI call;
  the γ-floor is unaffected (reported on γ_myo).
- **Major next FF units** (WLC constitutive / loaded-shell production / Layer-2 coupling) — PI-roadmap-gated
  (ENGINE.md §5/§6); not started unprompted.

Related: [[project-ff-gpu-native]], [[project-gamma-floor-likely-deficit]], [[feedback-oracle-is-crosscheck-not-truth]],
[[feedback-no-param-tuning-to-outcome]]. Docs: FF_STAGE6Q, FF_STAGE6R, FF_KIM_NETWORK_VALIDATION, ENGINE.md.

## Figure

`outputs/ff/figs/validation_sprint_6qr.png` — 4-panel sprint summary: (A) GPU↔CPU bit-parity per
kernel/integrator (all ≤ 2e-11), (B) production-envelope robustness (γ_myo vs N; stable across
k_xl 0.1→4.6e5 and f_myo 0–40), (C) κL/2R² bending-energy anchor (corrected <0.5%; raw under-counts
exactly (n−2)/(n−1)), (D) unified-architecture relax stability (all 5 finite + restoring under a 30 nm kick).

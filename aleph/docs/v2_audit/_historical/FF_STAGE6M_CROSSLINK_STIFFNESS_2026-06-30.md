---
kb_record:
  doc_id: FF_STAGE6M_CROSSLINK_STIFFNESS_2026-06-30
  title: Crosslink-stiffness correction (Ferrer 2008) — active γ-floor robust; stiff passive channel is a relaxation artifact (PI surface)
  authoritative_as_of: 2026-06-30
  supersedes: []
  status: current
---

# FF Stage 6M — crosslink-stiffness correction: active γ robust, stiff passive channel = artifact

**Date:** 2026-06-30 · **Engine:** `aleph/laws/` · branch `dcm/main` · **PI surface (anomaly)**

## Context

PI approved (2026-06-30) correcting the crosslinker stiffness `hand_kmc.ALPHA_ACTININ.link_k` /
`FILAMIN.link_k` from the broken **0.1 pN/µm** to the lit-anchored **4.6e5 / 8.2e5 pN/µm** (Ferrer 2008
PNAS AFM; ~4.5e6× too soft was a pN/µm-vs-pN/nm slip + AFINES soft-surrogate, FF_KIM §c). I applied it,
tested on the gbook A5000 + CPU, and it surfaced an anomaly that — per the visualize-anomaly rule — I
characterized and am surfacing rather than shipping.

## Finding — two channels, opposite behavior

Measured by the method-of-planes at the resting geometry + myosin modulator
(`gamma_floor.measure_gamma`), broken vs Ferrer-corrected crosslink stiffness, N=1000 ensemble.
Figure `outputs/ff/figs/crosslink_stiffness_gamma.png`.

| channel | broken (0.1) | Ferrer (4.6e5) | verdict |
|---|---|---|---|
| **γ_active / γ_myo** (active myosin) | 0.146 ± 0.009 pN/µm | **0.146 ± 0.009** pN/µm | **IDENTICAL — FLOORED** |
| γ_xl (passive crosslink elastic) | ~0.001 pN/µm | 65 ± 18 (up to 1600, blow-up) | spurious + unstable |

1. **The ACTIVE γ-floor is ROBUST.** γ_active (the myosin-generated tension, the production number) is
   **0.146 pN/µm = 1.5e-4 mN/m at BOTH stiffnesses — bit-identical**, ~2500× under band. Correcting the
   crosslink stiffness by **10⁶×** does NOT change the active cortical tension. This **definitively
   CLOSES the last γ-floor root-cause candidate** — the "cross-bridge 1000× too soft" suspicion
   ([[project-gamma-floor-layered-resolution]]). The floor is force-magnitude-limited (engaged motor
   density × per-motor force), independent of crosslink stiffness, consistent with the 6H/6k ablations.

2. **The stiff PASSIVE channel is a relaxation ARTIFACT.** With k_xl = 4.6e5, any residual relaxation
   stretch (~0.6 nm mean) becomes ~450 pN crosslink tension; the explicit/implicit relaxers leave such
   residuals and occasionally **over-stretch / blow up** (one seed: 345 nm stretch, 2.8e5 pN). This
   inflates γ_total / γ_xl to spurious values (a passive ELASTIC response to un-relaxed stretch, NOT
   active tension). It is why the stiff value cannot enter the production cortex relaxation as-is.

## Decision (conservative, PI to ratify)

- **`hand_kmc` preset kept at link_k = 0.1** (production stability) with a prominent ⚠️ pointing to the
  correct Ferrer value + this doc. The production γ-floor (`gamma_active`) is unaffected either way.
- **The correct stiffness 4.6e5/8.2e5 IS used where it has a robust solver:** `kim_network.shear_modulus`
  (the elastic-modulus axis, FF_KIM §c), with its own affine-shear relaxation.
- **Before the production preset can move to the stiff value**, the cortex relaxation must be made robust
  to stiff crosslinks — treat crosslinks as near-rigid CONSTRAINTS (a crosslink analogue of the actin
  reshape) or a stretch-controlled / adaptive-dt implicit solver. This is the concrete next numerics
  task (the same family as the FF_KIM finite-EA implicit follow-up).

## What this does NOT change

The γ-floor scientific conclusion ([[project-gamma-floor-likely-deficit]]) STANDS and is strengthened:
active cortical tension is force-magnitude-limited, robust to buckling / connectivity / finite-
extensibility (6H) / turnover (6k) / **and now crosslink stiffness (6M)**. The one open lever remains
the missing experimental load-engaged motor-density datum. (The crosslink stiffness was the ELASTIC-
modulus error from FF_KIM §c — real and worth fixing for fidelity — but it is a separate axis from the
active γ and does not lift it.)

Tests: `test_gamma_floor_equilibrate_implicit_matches_explicit` passes at link_k=0.1 (the explicit↔
implicit consistency holds at soft stiffness; at stiff stiffness explicit is CFL-throttled — the reason
for the implicit/constraint follow-up). SE candidates: `SE_REGISTRATION_CANDIDATES_2026-06-30.md` §5.

---
kb_record:
  doc_id: FF_STAGE6L_UNIFIED_ARCHITECTURE_2026-06-30
  title: Unified actin-architecture framework — one weave() builder (increment 1 — cortex + filopodium)
  authoritative_as_of: 2026-06-30
  supersedes: []
  status: current
---

# FF Stage 6L — unified actin-architecture framework (task a, increment 1)

**Date:** 2026-06-30 · **Engine:** `ffn_sim/ff/` · branch `dcm/main`

## Why

PI 2026-06-04 ([[project-unified-actin-architecture]]): the cortex, lamellipodium, filopodium,
stress-fibers, microvilli and traction apparatus are ONE category — actin structures that differ only
in how the SAME primitives (filaments + nucleator + crosslinker + motor) are WOVEN onto a manifold —
and must be implemented UNIFIED, not piecemeal. The architectural knowledge was collected
(`docs/ACTIN_ARCHITECTURE_NOTES.md`); PI has now greenlit implementation. Per the collect-then-ratify
scope, **increment 1 is a design-validation prototype**: build the smallest pair (cortex + filopodium)
from ONE builder to ratify the primitive mapping, and surface the full table for PI sign-off before
shipping the remaining structures.

## Increment 1 — built + validated

One builder `ff/weave.weave(spec)` produces both structures from an `ArchitectureSpec`, with **no
structure-specific branch beyond the placement manifold + the crosslinker bind-mode angle filter**:

| structure | manifold | filament | crosslinker (bind) | motor | validation metric (result) |
|---|---|---|---|---|---|
| **CORTEX** | sphere | 1000, isotropic, formin/Arp | α-actinin/filamin (any) | NMIIA bipolar | ρ=0.80/µm², z=2.0, S=0.0 (isotropic); **γ_active reproduces build_crosslinked_cortex BIT-EXACT** |
| **FILOPODIUM** | bundle | 20, parallel, formin | filamin (parallel) | none | count=20 (lit 10–30), **S=1.00** (parallel), spacing=**8.0 nm** (fascin ~7–8) |

- **CORTEX parity is exact** (γ_active ratio 1.000 vs `gamma_floor.build_crosslinked_cortex`): the
  sphere path mirrors the H.3 cortex builder with the RNG call order preserved, so routing the cortex
  through the unified builder does NOT regress the γ-floor production number.
- **FILOPODIUM** hits every architectural band (count, nematic order S, fascin spacing) from the SAME
  `weave()`. The unification claim (`test_weave_one_builder_two_architectures`): one builder → distinct
  architectures (cortex S<0.3 isotropic vs filopodium S>0.95 parallel).

Modules: `ff/architecture_spec.py` (the table in code), `ff/weave.py` (the builder),
`ff/architecture_metrics.py` (S, bundle count, spacing, cortex metrics). Tests: `tests/ff/test_weave.py`
(3). Figure: `outputs/ff/figs/unified_architecture.png`. Reuses unchanged: `cortex_assembly`,
`fiber_network`, `gamma_floor` (CrosslinkedCortex + KDTree pairing), `hand_kmc`, `kim_network`,
`network_warp` (the GPU relax/forces are structure-agnostic flat arrays — any woven structure runs on
the A5000 via the same path).

## Increment 2 — LAMELLIPODIUM + Arp2/3 angle-harmonic branch kernel (built, 2026-06-30, PI "둘 다")

Third structure from the SAME `weave()`, `manifold="patch"`: an Arp2/3 DENDRITIC array. Mothers seeded
at ±35° about the protrusion axis; each spawns daughters at the rest branch angle θ₀=70° toward the
OTHER ∓35° mode (so all filaments sit at ±35° two-mode AND every junction is at 70°). The branch angle
is maintained by a NEW Warp kernel `network_warp.branch_angle_kernel` — the MECHANISTIC angle-harmonic
Arp2/3 branch (CLAUDE.md: U=½·k_angle·(θ−θ₀)² with thermal fluctuation, NOT a rigid 72° constraint).

- **Lit-anchored constants** (Magic-Number-Blocked, `configs/phase1_h5.yaml`): θ₀ = **70°** (Fäßler 2020
  in-cell cryo-ET 68±9°); k_angle = **0.173 pN·µm/rad²** = kT/Var(θ), σ_θ=9° equipartition; ±35° two-mode
  (Mueller 2017, protrusion-axis-relative); branch density 1.25/µm (Vinzenz 2012, PI-gated cite).
- **Built metrics**: branch junctions **69.1±9.6°** (target Fäßler 70±9), orientation **+33.3/−35.0°**
  two-mode (target Mueller ±35), 78 % in the ±20–50° band. Kernel bit-parity vs numpy 1.3e-15; a
  perturbed 120° branch relaxes to exactly 70°. The lamellipodium relaxes STABLY on the A5000 with the
  branch kernel (2000 steps, 1.5 s, finite; CFL includes k_angle's k_eff=k_angle/ℓ²).
- New: `ff/architecture_spec.LAMELLIPODIUM` + `FilamentSpec` two-mode/branch fields + ARP23 constants;
  `weave._build_lamellipodium_patch` (dendritic generator, emits branch triples + Arp2/3 anchors);
  `network_warp.branch_angle_kernel` + relax_on_device branch support; `architecture_metrics`
  branch_angle_distribution + two_mode_orientation; `CrosslinkedCortex.branch_triples`. Tests +3.

## Full per-structure table — increments 3+ FOR PI SIGN-OFF

Built only after PI ratifies + authorizes the unsourced constants below.

| structure | manifold | filament (nucleator) | crosslinker (mode) | motor | metric + lit band |
|---|---|---|---|---|---|
| cortex ✅ | sphere | bimodal Arp2/3 + formin | filamin/α-actinin (any) | NMIIA bipolar | mesh/z/ρ (done) |
| filopodium ✅ | bundle | formin parallel | fascin (parallel) | none | count 10–30, spacing 7–8 nm (done) |
| lamellipodium ✅ | patch | Arp2/3 ±35° dendritic | filamin (any) sparse | none/sparse | branch 70°, two-mode ±35° (done) |
| stress-fiber | FA–FA | formin antiparallel | α-actinin parallel periodic | NMIIA sarcomeric | sarcomeric period 0.5–1 µm |
| microvilli | finger | formin parallel ~20–30 | espin/fimbrin/villin (parallel) | none (myosin-1a shaft) | count 20–30, spacing ~33 nm |
| traction/FA | cortex/SF + substrate | = SF/cortex | + clutch catch-bond | SF retrograde flow | traction ~102 nN (Gil-Redondo 2023) |

**PI-gated before increment 2 (HARD rules):**
1. **The remaining 4 structures' presets** are NOT yet defined — PI ratifies this table first.
2. **Unsourced bundler kinetics** — fascin (filopodium), espin/fimbrin/villin (microvilli), the
   Arp2/3 branch-angle harmonic stiffness, lamellipodium ±35°, SF sarcomeric periodicity, FA 102 nN —
   NONE are registered SourceEvidence (`SE_REGISTRATION_CANDIDATES_2026-06-30.md` §5/architecture).
   Increment 1 deliberately uses ONLY anchored presets (α-actinin/filamin/NMIIA) + lit-grounded
   geometry (8 nm bundle spacing, 20 filaments) so it introduces no tuned constant; even the 8 nm
   spacing is fascin GEOMETRY (Claessens/Volkmann) and should get an SE row before production.
3. **One new kernel for increment 2** (lamellipodium): the angle-harmonic Arp2/3 branch-junction force
   (70°±thermal, NOT rigid — CLAUDE.md table) → `network_warp.py`, must be bit-parity-tested vs numpy
   AND the CFL `dt` must include the branch-angle stiffness (the turgor-breathing-mode blowup lesson).

## Scope boundary (state explicitly)

This is an ARCHITECTURE deliverable — does one builder span the structures and reproduce their
architectural metrics. It is **separate from the open γ-floor contractility question** (force-magnitude
/ load-engaged-motor-density limited, [[project-gamma-floor-likely-deficit]]); do NOT judge the unified
cortex against the γ band it cannot move. The bundle representation is mechanistic (explicit parallel
filaments + parallel crosslinkers), not a lumped thick-fiber proxy (CLAUDE.md mechanistic rule).

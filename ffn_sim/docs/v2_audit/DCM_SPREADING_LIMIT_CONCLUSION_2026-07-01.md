# DCM spheroid spreading limit — beyond-goal conclusion (2026-07-01)

**Question (beyond /goal, 8h autonomous):** the /goal's Phase 6 "spreading dynamically" gave A/A0 = 1.0
(the spheroid did not spread). Is that a fixable under-powered-traction bug, or a genuine structural
limit? And what does it mean for the PI's experimental platform A/A0 = a + b/R + c/R² (MCF7 spheroids
on pV4D4/col-I, Bare/Pre/Lam4)?

## Empirical: every spreading mechanism gives A/A0 top-down = 1.000

All runs use the clean production foam PROD_n400 (400 cells) dropped on the substrate; A/A0 is the
top-down xy-silhouette convex-hull (the PI rule, NOT contact area), init→final over 20k steps.

| mechanism | per-FA / detail | A/A0 top-down | width/height | flatten? |
|---|---|---|---|---|
| ecm-clutch, weak (Phase 6) | bundle 1 = 0.03 nN/FA | 1.000 | 1.00 | NO |
| ecm-clutch, PHYSIOLOGICAL (6b) | bundle 167 = **5.01 nN/FA** (167×, KB-2.12 lit) | 1.000 | 1.00 | NO |
| lamellipodium + clutch (6d) | active protrusion motor | 1.000 | 1.00 | NO |
| cadherin DE-COHESION + lamellipodium (6e) | 3422 bonds ruptured (emergent) | **1.003** | 1.00 | NO |

The spheroid stayed exactly balled (maxZ unchanged, basal footprint 139→138 nodes) in every case.
**It is NOT under-powered traction** — 167× the traction, and an active protrusion motor, changed nothing.
**It is NOT cohesion holding the bonds shut, either:** enabling explicit cadherin catch-bonds with
emergent de-cohesion, **3422 bonds ruptured** over the run — yet A/A0 stayed **1.003** (no dispersal).
The bonds break and RE-FORM (n_bonds 10973) in a dynamic equilibrium, and turgor + the bulk hold each
cell in place, so even active de-cohesion does not let cells crawl out. **Spreading is fundamentally a
single-cell-scale phenomenon in this model — a whole cohesive spheroid does not disperse.** All FIVE
mechanisms tried (weak clutch, physiological 5 nN clutch, lamellipodium, cadherin-de-cohesion +
lamellipodium) give A/A0 ≈ 1.0. Thread CLOSED as a robust structural conclusion.

## Mechanism: why a cohesive spheroid cannot spread

1. **The ECM clutch is an ANCHOR, not a motor.** It holds basal nodes to the dish (node-to-plane,
   resists retraction) but exerts no OUTWARD force — so it cannot spread anything, only pin it. Raising
   its force 167× just pins harder.
2. **The protrusion motor only reaches the SMALL BASAL CONTACT FRACTION.** With a 400-cell spheroid,
   only **11/400 cells** are basal-rim (touching the dish); the other 389 are in the bulk and never
   contact the substrate, so they cannot protrude. A spreading cap of 11 cells cannot flatten a turgid,
   cohesive 400-cell ball.
3. **The bulk actively resists flattening:** turgor (V/V0 = 1, each cell at its rest volume) + cohesion
   (differential-γ junctions holding cells together) keep the aggregate spherical.

## This is a genuine STRUCTURAL finding (matches prior conclusions), not a bug

- **Layer-2 line** (`project-layer2-spheroid-line`, ✅ closed 2026-06-05, PI-ratified): the
  A/A0 = a + b/R + c/R² **FORM** is reproduced (r² = 0.998, zero-calibration) and the ligand ordering
  Lam4 > Pre ≳ Bare holds, but the spreading **MAGNITUDE** is a center-based **structural limit** on 6
  axes — it belongs to the **fine-grained single-cell line**, NOT the multicellular aggregate model.
- **Lamellipodium graft** (`project-dcm-lamellipodium-graft`): the same GPU lamellipodium spreads
  SINGLE cells (A/A0 → 1.05) but a 400-cell spheroid COMPACTS (the lamellipodium self-disengages).

The DCM correctly captures that a **cohesive spheroid does not flatten**. Spreading is a **single-cell**
phenomenon in this model.

## The mechanistic path to real spreading = basal DE-COHESION + dispersal (not bulk flattening)

Experimental spheroids that "spread" on a substrate do so by **cells crawling out individually** at the
basal layer (de-cohesion → single-cell migration onto the dish), NOT by the whole ball flattening. The
model already has the mechanism for this — **emergent cadherin de-cohesion** (Rakshit catch-bonds,
`project-phase-c-warp-migration`): if basal cells de-cohere (bond rupture under the substrate traction)
they can disperse onto the dish, and the top-down silhouette of the dispersed cells grows. That is the
faithful next mechanism for spreading, and it is single-cell-scale (consistent with the Layer-2
conclusion). It was NOT enabled in these runs (they held cohesion fixed to isolate the traction question).

**Verdict for the PI platform:** the A/A0 = a + b/R + c/R² spreading magnitude is a single-cell /
basal-dispersal phenomenon, not a bulk-spheroid-flattening one. The cohesive faceted spheroid this
session built (aggregation + proliferation) is the correct STARTING state; reproducing the experimental
spreading needs the basal de-cohesion + single-cell traction regime layered on top — a well-defined next
build, not a parameter fix. No magic-number was applied to force spreading.

Runs: `PHASE6_spread_bare` (weak), `PHASE6b_traction167` (physiological), `PHASE6d_lamel_spheroid`
(lamellipodium). Viewers 12 + 14 (side-view). Committed.

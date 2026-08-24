---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# DCM spreading — the missing physics: ACTIVE AREA GENERATION (plan, 2026-07-13)

## Why this plan (the diagnosis)
Real cells/tissues spread; the DCM does not (A/A0 stays ~1.0–1.1 even with motility + T1-rate flow). The root, named
in the 2026-07-13 discussion, is NOT the timescale gap (that line is closed: `DCM_TIMESCALE_ATTACK`, `DCM_T1_RATE_*`)
and NOT "it's FF's job." It is that **the DCM cell is a turgid closed shell whose surface tension MINIMISES area — so
it structurally OPPOSES spreading — and the active machinery that GENERATES spreading area is absent or too weak.**

Real spreading = **active area generation**: the cell does WORK to increase its footprint. Three coupled pieces:
1. **Lamellipodial protrusion** — actin polymerisation pushes the membrane outward → the cell actively grows a thin
   sheet (adds surface/footprint area).
2. **Substrate adhesion + molecular-clutch traction** — integrin adhesions grip the substrate and convert actin
   retrograde flow into traction that pulls the edge outward and HOLDS the cell flat (pancake).
3. **Cortical-tension down-regulation** — the cell actively lowers cortical tension where it spreads so protrusion wins.

The DCM has NONE of these as an area-generating term. Its "motility" is a bulk directional push (translates a cell,
does not flatten/grow it); its surface tension γ is fixed and area-minimising; substrate coupling is a weak well/plane.
Prior attempts (`project-dcm-lamellipodium-graft`: single A/A0→1.05 but 400-cell RE-COMPACTS; `project-dcm-spreading-
foundation`: single A/A0→3.8, spheroid COMPACTS) failed because the added drivers were **too weak to beat the DCM's
surface tension + turgor** — not because the DCM structurally cannot deform.

## The governing physics — a WETTING balance (must beat surface tension at physiological values)
A deformable cell on a substrate spreads (wets) when the **substrate adhesion energy exceeds the surface-tension cost
of the increased area** — the active-wetting / Young–Dupré balance (Douezan 2011, Beaune 2014, Gonzalez-Rodriguez 2012):

  spread ⇔  W_adh (per basal contact area)  ≳  γ (cortical/surface tension)   [+ active protrusion drives beyond passive]

- Passive equilibrium spread (contact angle) is set by cosθ ≈ W_adh/γ − 1.
- **Active** spreading (real cells go far past passive wetting) adds the protrusion term: the cell actively GROWS its
  rest area A0 at the basal/leading edge (actin polymerisation + membrane unfolding), and the clutch traction pulls it.

**The key requirement (why prior attempts failed): all three must be at PHYSIOLOGICAL strength so the balance genuinely
favours spreading.** Real cells DO spread → physiological W_adh + protrusion DO beat physiological γ. The plan grounds
each at its in-vivo value and checks the balance is spreading-favourable BEFORE running (no under-powered drivers).

## Design — three mechanistic terms added to the DCM cell

### T-1. Substrate adhesion energy (basal wetting) — the passive floor
Reward basal contact area with the substrate: `E_adh = −W_adh · A_basal`, applied as an attractive force on basal nodes
(within a contact range of the substrate plane/ECM) that also spreads them along the substrate. This flattens the cell
until balanced by γ. Grounding: `W_adh` from integrin–ligand areal adhesion energy (integrin density × bond energy ×
engaged fraction; lit ~0.1–1 mJ/m²); `γ` = the DCM cortical tension (MCF7 ~1e-2 N/m, Moazzeni). Sanity: the passive
contact angle must match Young–Dupré for the chosen W_adh/γ.

### T-2. Active protrusion = rest-area growth at the spreading edge (the area GENERATOR)
The core missing term. At basal-adhered / leading nodes, ACTIVELY grow the cell's preferred surface area:
`dA0/dt = +k_prot` while the edge is adhered and below a max spread (membrane-reservoir / polymerisation limit). This is
the mechanistic stand-in for lamellipodial actin polymerisation + membrane unfolding — the cell MAKES area. Grounding:
`k_prot` from the lamellipodium protrusion rate (actin polymerisation ~ µm/min edge advance × edge length) and the
membrane-reservoir area budget (cells unfold ~2–4× area). Cap A0 at the physiological max spread area.

### T-3. Clutch traction + anisotropic tension reduction (hold the spread)
Basal adhered nodes pull OUTWARD (radially, away from the cell centroid) with a traction proportional to the local
actin flow — the molecular-clutch force that converts protrusion into net edge advance and HOLDS the pancake. Couple
with a **basal cortical-tension reduction** (γ lower on the spread basal surface than the apical) so the flattened
state is stable. Grounding: clutch traction per adhesion from the FF clutch work (per-clutch ~5–25 pN, KB); tension
anisotropy from apical-vs-basal cortex lit.

## Sanity Gates (before first production)
- **G1 dimensional**: W_adh [J/m²], γ [N/m], k_prot [m²/s], traction [N]; the wetting number W_adh/γ dimensionless.
- **G2 balance check (the lesson from prior failures)**: at the chosen physiological values, verify W_adh/γ + protrusion
  predicts a SPREADING-favourable equilibrium (contact angle < 90°, A/A0 target > 2) — do NOT run under-powered drivers.
- **G3 volume conservation**: spreading is area↑ at constant volume (cell flattens, V/V0=1) — the pancake must conserve V.
- **G4 single-cell first**: a single DCM cell on a substrate must reach a physiological spread A/A0 (2–8× for stiff
  substrate) and a stable pancake — BEFORE any aggregate run.
- **G5 no artefact**: A/A0 from the top-down silhouette must be REAL flattening (basal contact area ↑, maxZ ↓), not
  node-ejection — visual-verify in the browser (the metric lied before; see `feedback-visualize-anomalies-for-pi`).
- **G6 aggregate wetting**: a spheroid on a substrate must WET (spread) — the peripheral cells crawl out + the tissue
  flattens — not re-compact (the prior failure mode). Cohesion (cadherin) sets how far it wets vs stays cohesive.

## Validation targets (literature)
- Single cell on stiff 2D: A/A0 ~ 2–8× (spread vs suspended), maxZ drops (pancake). Contact area rises to the plateau.
- Spheroid wetting: A/A0 grows over hours; the spreading rate + final wet area set by W_adh/γ + active protrusion
  (Douezan/Beaune). MCF-7 (epithelial, cohesive) wets modestly; mesenchymal (MDA-MB-231) wets/disperses far more.
- Cross-check vs the FF single-cell crawl/spread numbers where they overlap (same physics, finer scale).

## Increments (ordered; strongest-lever first, single-cell before aggregate)
1. **T-1 + T-2 on a SINGLE cell** (substrate adhesion + active area growth): does one cell spread to A/A0 2–8× as a
   stable pancake? Gate G2/G3/G4/G5. This is the make-or-break — if a single cell won't spread at physiological
   W_adh+protrusion, diagnose before adding traction.
2. Add **T-3 clutch traction + tension anisotropy**: sharpen the edge advance + stabilise the pancake.
3. **Aggregate wetting** (G6): spheroid on substrate → does it wet (spread) vs re-compact? Sweep cadherin cohesion
   (epithelial cohesive-wet vs mesenchymal disperse).
4. **Biology-time** (couple to the existing large-dt BDF2 + T1-rate): does the aggregate wet over physiological time?
5. Viz HTML (browser-verified pancake/wetting) + per-cell stress/strain + figures at each step.

## Relationship to the other lines
- Orthogonal to the timescale/T1-rate line (that supplied tissue FLOW; this supplies AREA GENERATION). A fully spreading
  tissue needs BOTH: cells that generate area (this plan) AND rearrange (T1-rate) at biology-time.
- The area-generation machinery is fine-grained (actin/adhesion) = FF's native strength; this plan adds a **coarse but
  mechanistic area-generation term to the DCM** so the DCM tissue can spread without the full FF filament detail. If the
  coarse term proves too lossy, the fine-grained version lives in FF (crawl) — but the DCM term is the tractable first test.
- Uses the image-based stress/strain analysis capability (2026-07-13 discussion) as the read-out on spread cells.

## Status
Plan written. Next (PI-gated): start Increment 1 — ground W_adh, γ, k_prot at physiological values, run G2 balance check,
then a single-cell substrate-spread test with visual verification. No production run until G2 confirms the balance
favours spreading (the explicit fix for the prior under-powered attempts).

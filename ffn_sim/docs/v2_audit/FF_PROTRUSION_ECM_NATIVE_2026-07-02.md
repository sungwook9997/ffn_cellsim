# FF native full-cell filopodium protrusion into a collagen ECM (2026-07-02)

**Date:** 2026-07-02  **Engine:** FF (Warp)  **Device:** gbook A5000 (cuda:0)  **Branch:** dcm/main
**PI directive:** "실제 네이티브 셀 풀 컴포넌트로 구성해서 보여줘야지 … ecm도 mikado 네트워크로 … 모든 구성
네이티브로." — show the piece-1 protrusion on a REAL native full-compartment cell pushing into an ECM, not an
abstract floating bundle. This is that run.

## What this is

`scripts/ff_protrusion_into_ecm.py` assembles a whole cell at native scale + a collagen ECM and lets the piece-1
polymerization ratchet drive a filopodium out of the cortex INTO the matrix, so the protrusion load is
**EMERGENT** (steric contact with collagen), not imposed:

- **N = 499,218 nodes** — cortex **70,686 filaments** (494,802 nodes, KB-3.18 areal density, MCF7 R=7.5µm) +
  filopodium (24 filaments grafted onto the cortex surface, base crosslinked) + nucleus (3000 beads, 0.25R).
- Compartments at physiological setpoints: turgor ΔP (regulated), plasma-membrane γ_mem inward, stiff bilinear
  nucleus, quasi-static overdamped relaxation. **ECM** = a 3D **Mikado collagen network** (`ff/ecm_mikado.py`,
  120 fibers, **mesh 1.0 µm** — collagen-I 1–5 µm range), far face pinned (embedded in bulk matrix). Filopodium
  tips make soft excluded-volume contact with it (`network_warp.soft_contact_kernel`).
- 71,009 crosslinks total; 70 macro polymerization steps × 120 relax substeps; **2265 s (38 min) on the A5000.**

## Result — the ECM supplies the protrusion load EMERGENTLY, throttling the ratchet onto Mogilner-Oster

| | filopodium tip advance | emergent load / tip | end velocity v |
|---|---|---|---|
| **free (no ECM control)** | **+1.901 µm** | 0 | 0.623 = v0 (unimpeded) |
| **into ECM** | **+1.359 µm** | rises to **8.04 pN** | **0.051 µm/s** |

The filopodium grows freely (v0) until its tip reaches the collagen at ~2 s, then the ECM resists: the emergent
contact load builds to **8.04 pN** and the polymerization velocity drops to **0.051 µm/s**. The analytic
Mogilner-Oster prediction v(8.04 pN) = v0·exp(−8.04·δ/kBT) = 0.623·exp(−2.54) = **0.049 µm/s** — the sim
matches to ~4%. **The load is EMERGENT** (computed from the geometric penetration depth of the tip into ECM
fibers, k·(r_contact−L)); it is NOT imposed. This is the emergent-load closed loop that piece-1's isolated
filopodium (imposed load) explicitly deferred — a partial delivery of piece-4 via the ECM obstacle.

## Not circular (audit-checked)

The contact load is `soft_contact_kernel` = k_contact·(r_contact − L)·û from the ACTUAL tip↔ECM-node distance L
(penetration geometry). The poly kernel reads `force[barbed]` (= that geometric load) and applies
v=v0·exp(−fδ/kBT). The figure's analytic curve `ratchet_velocity_np(load, ...)` is computed from the SAME
measured load. So v-matches-analytic reflects the kernel correctly responding to a genuinely emergent geometric
load — not the analytic being fed back.

## New primitives this run introduced (commit 37676e9)
- `ff/ecm_mikado.py` — 3D Mikado collagen network (Wilhelm-Frey 2003 random-rod; κ=kBT·Lp collagen PI-gated;
  mesh EMERGENT, collagen 1–5 µm target). Validated: geometry + network resists a push.
- `network_warp.soft_contact_kernel` (excluded-volume, emergent load) + `freeze_kernel` (pinned-boundary BC).

## Honest scope / caveats
- Cytoplasm viscosity η=65.9 Pa·s is NOT wired as physical per-node drag; the relaxation is quasi-static
  (stiffness-scaled pseudo-time) → observables are mechanical-equilibrium quantities (audit 2026-07-02; docstring
  corrected). True η-dynamics is a refinement.
- ECM constants (collagen Lp=20 µm → κ, K_XL=100, K_CONTACT=200 pN/µm, R_CONTACT=0.35 µm) are lit-anchored /
  PI-gated, NOT tuned to an outcome — the emergent load emerges from geometry + collagen stiffness, and the
  free-vs-ECM comparison is the controlled test.

## Morphology render (PI 2026-07-02 "이건 필로포디움이 아니라 그냥 작대기")

The first render showed bare actin filament centerlines (24 hairlines in a 0.12 µm cross-section) → read as one
stick. Fixed: filopodia now render as **membrane FINGER tubes** on a cell morphology (`ff_cell_morphology.py`,
`ff_viewer_html` mesh support, `ff_protrusion_morphology_viz.py`), and the driver grows a **fan of N filopodia**
(`--n-fingers`), not one. A native **7-filopodia** run (protrusion_ecm_native7, N≈500k, A5000, 33 min) gives the
full array: **7 fingers 2.66–3.63 µm long** fanning from the +x leading edge into the collagen ECM (mean tip
+1.17 µm ECM vs +1.30 µm free; mean load 1.09 pN — lower than the single-finger 8 pN because the fan spreads the
tips, most still exploring). `protrusion_ecm_native7_morph.{png,html}` = semi-transparent plasma membrane +
nucleus + 7 filopodia fingers (growing, animated 70 frames) + actin cores + collagen ECM. Finger radius 0.18 µm
(filopodium diameter ~0.1–0.3 µm, Mattila-Lappalainen 2008).

## Files
- `scripts/ff_protrusion_into_ecm.py`, `scripts/ff_protrusion_ecm_viz.py`, `ff/ecm_mikado.py`.
- `outputs/ff/figs/protrusion_ecm_native.{png,html}` — (a) free vs ECM tip advance, (b) emergent load throttles
  v onto analytic v(f), (c) 3D cell+nucleus+filopodium+ECM; interactive animated HTML (cortex + nucleus points +
  filopodium + ECM lines playing over the 70 growth steps).

Related: [[project-ff-active-movement-pieces]], FF_POLYMERIZATION_PIECE1, oracle-is-crosscheck.

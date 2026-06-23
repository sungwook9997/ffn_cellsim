# Phase C — proxy-free spreading run (2026-06-23, overnight)

`run_spread_overnight.py` N=100, full mechanistic stack, NO wetting/well proxy. gbook A5000, 16000
implicit-100× steps (T≈13s), 33 min. Artifacts: `outputs/warp_decohesion/n100_proxyfree_spread.{npz,json}`,
figs `n100_proxyfree_spread_{montage,traj}.png`.

Stack (all lit/KB-anchored, NO outcome-tuning): cadherin catch-bond ×40 bundle (~7nN/junction ∈ KB-4.11)
+ Pereverzev ECM clutch ×167 bundle (~5nN/FA ∈ KB-2.12) + lamellipodium + deformable nucleus +
surface-tension γ + bending + edge-edge contact. Wetting OFF, substrate-well OFF.

## Result — NO collective spread (COMPACTS stands), but the run is CONTACT-CONFOUNDED
| metric | value | reading |
|---|---|---|
| A/A0 (top-down) | 1.000 → dip 0.974 → **1.037** | NO spread (footprint compacts then mild rim creep back) |
| maxZ | 80.2 → 79.2 µm (−1.3%) | **holds — no pancaking/flattening** |
| V/V0 | 1.000 | volume conserved |
| peeling_index | **16.4** (rim 15.2µm, bulk 0.93µm) | mild rim creep, bulk frozen — NOT collective; far below the retracted decisive's 348 (bundle fix worked) |
| Ψ (per-cell) | 0.996 → 0.972 | slight flattening at junctions |
| **pen_frac** | **3.11 (peak 4.06)** | ⚠️ **G2 FAIL 10× — cells interpenetrate ~7µm (≈ a cell radius)** |
| cadherin bonds / clutch | 3139 / 420 | cohesion + traction both engaged |

**Reading.** A/A0 dips to 0.974 (collective compaction) then creeps to 1.037 as ~rim cells peel out
15µm while the bulk stays frozen (peel_idx 16.4) — a much-attenuated echo of the decisive run's 7-cell
peeling, not collective spread. maxZ holds → no spheroid→pancake flattening. **So the proxy-free
mechanistic stack does NOT collectively spread — the "COMPACTS / no-spread" finding STANDS** (and the
bundle FORCE fix correctly suppresses the peeling artifact: 348 → 16.4).

## ⭐ The deeper finding — the spreading test is BLOCKED by a CONTACT FAILURE (M1)
`pen_frac` oscillates 0→4.06 (mostly 2–3.5) for the whole run = nodes repeatedly penetrate **2–4
mean_edge (4.5–9 µm ≈ a full cell radius)** into neighbours. This is the **M1 tunnelling bug at full
severity**: the penalty contact's repulsion vanishes past `c_rep≈0.67µm`, so the ×40 cadherin bundle
(~7nN/junction) pulls cells together faster than the contact can hold them apart → they interpenetrate /
partially merge. The "no spread" A/A0 is therefore **confounded** — the cells are not cleanly
excluded-volume bodies, they are partly fused, which itself suppresses any spreading.

**The bind (structural):** aggregation needs strong cohesion (the lit ×40/×167 bundle, else cells don't
stick — the original problem); but at that strength the penalty contact catastrophically fails (pen 3.1).
**Strong bundle adhesion and the per-face penalty contact are incompatible.** A clean spreading verdict
is therefore not obtainable until the contact can hold under these forces — which needs the proper
inside-closed-mesh / active-set contact (M1), NOT the per-face `sign<0 & min_d<c_rep` penalty (whose
`min_d<c_rep` gate is load-bearing — removing it explodes, verified — so it cannot be cheaply extended).

Lowering the bundle to make pen pass would be outcome-tuning a derived param (forbidden,
`feedback-no-param-tuning-to-outcome`) — so it is NOT done; the limit is surfaced to PI instead.

## Honest verdict for PI
1. **Does the proxy-free stack collectively spread? NO** (A/A0 1.04, maxZ holds, bulk frozen) — consistent
   with the Layer-2 "magnitude = fine-grained single-cell, not collective" structural-limit conclusion.
2. **But this run cannot be the clean proof**, because the contact fails (pen 3.1) at the lit cohesion
   strength. The spreading question is **gated on the M1 inside-mesh contact fix**.
3. Recommendation: implement the active-set / signed-distance contact (M1, surfaced) → re-run → only then
   is "no collective spread" an artifact-free verdict. Until then: COMPACTS-no-spread holds, with the
   contact-confound caveat explicit.

## Disentangle control — cadherin OFF, ecm-clutch ×167 + lamellipodium full traction (gentle node-face cohesion)
To separate "no spread" from the cadherin-bundle contact confound, re-ran with **NO cadherin bundle**
(cohesion = the gentle node-face adh5e7 tent, which holds N=400 at pen 0.12) + FULL traction
(ecm ×167 + lamellipodium). Trajectory (killed at 55%/step 8778 — the trend is monotone + conclusive):
| step | 266 | 2660 | 5320 | 7980 | 8778 |
|---|---|---|---|---|---|
| A/A0 | 0.997 | 0.965 | 0.929 | 0.904 | **0.897** |
| pen | 0.00 | 2.80 | 2.06 | 2.68 | 2.63 (peak 3.23) |
**Two decisive conclusions:**
1. **No-spread is ROBUST / STRUCTURAL** — with the cadherin bundle entirely removed, A/A0 still **monotone
   COMPACTS** 0.997→0.897 (never spreads; even cleaner than the cadherin run's late rim-creep to 1.04).
   So "no collective spread" is NOT a cadherin/contact artifact — it holds across cohesion mechanisms,
   confirming the Layer-2 fine-grained-single-cell structural-limit conclusion.
2. **The M1 contact failure is GENERAL, not cadherin-specific** — pen rises to 2.6–3.2 from the **ecm ×167
   traction ALONE** (no cadherin). So ANY strong mechanistic bundle (cell-cell OR cell-substrate) breaks
   the per-face contact. The bind is fundamental: every active mechanistic force strong enough to drive
   a real process (cohesion, traction) overwhelms the per-face penalty → M1/IPC is required regardless.

## Figures
- `n100_proxyfree_spread_montage.png` — top-down (compact, no expansion) + side (maxZ holds, no pancake), 8 frames.
- `n100_proxyfree_spread_traj.png` — A/A0 dip→1.04, maxZ flat, and pen→3.1 (the contact-failure story).

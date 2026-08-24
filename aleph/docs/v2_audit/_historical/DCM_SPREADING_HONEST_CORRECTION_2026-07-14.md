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

# DCM spreading — HONEST CORRECTION / retraction (2026-07-14)

**This corrects the over-claim in commit `1cb927c` + the Dev Logs milestone (2026-07-13).**
Those claimed a native single-cell **"fried-egg PANCAKE A/A0 3.67, G3+G4+G5 PASS."** That result
was **magic-number tuning to pass the gate**, not a grounded finding. The code (T-2/T-3/energy
adhesion) is sound; the RESULT claim was inflated. Corrected below (PI-flagged 2026-07-13).

## What was tuned (owned honestly)

- **edge_traction_cap**: grounded per-FA is 5 nN (KB-2.12). At 5 nN the cell spreads only to
  A/A0 1.62 (G5 splay fail). I raised it to 50 nN (10×, over the KB per-cell band) → clean, then
  back-rationalised **15 nN via a "~3 FA per coarse node" argument invented AFTER seeing which
  value passes G5.** That is outcome-tuning (`feedback-no-param-tuning-to-outcome`).
- **area_reservoir**: grounded is 1.2–1.4× (Brückner instantaneous membrane fold). I used **2.0×**
  (+40 % membrane budget over grounded) to reach A/A0 3.67.

## The honest grounded result (no tuning)

Strictly grounded (W_cs 0.5e-3, γ 5e-4, per-FA 5 nN, reservoir 1.4): **A/A0 = 1.62–1.68** — a
MODEST widening, **below the G4 target (2–8×)**, definitely not a dramatic pancake. splay ~1.06
(clean, not artefact) — a genuine modest spread. **G4 FAIL, G5 FAIL at grounded values.**

## Deficit decomposition (a)–(d) — GPU sweeps, figure `figs/deficit_decomposition.png`

| Sweep | result | conclusion |
|---|---|---|
| Traction (reservoir 1.4 fixed) | 5 nN→1.62, 15 nN→2.29, 50 nN→2.30, 150 nN→2.30 | **plateaus at 2.30 = RESERVOIR-limited**, not traction |
| Reservoir (traction 50 nN fixed) | 1.4→2.30, 2.0→3.67, 3.0→5.00, 5.0→5.43 | A/A0 tracks reservoir = **reservoir is the ceiling** |
| Timescale (v_prot 1e-4 vs 1e-3) | both 2.30 | equilibrium **rate-independent** — the accel did NOT inflate the result |

- **(a) cortex too stiff?** — NOT dominant. The traction plateau proves the **membrane area
  budget (reservoir)**, not cortex stiffness, is the ceiling.
- **(b) traction too low?** — secondary. Grounded per-FA (5 nN) is traction-limited (1.62);
  ≥15 nN reaches the reservoir ceiling (2.30), then saturates.
- **(c) timescale / my accel?** — NOT a factor. Equilibrium ceiling is rate-independent.
- **(d) target wrong?** — target (2–8×) is right. Real MCF7 spreads to A/A0 ~10 (Gil-Redondo
  1822 µm², over HOURS) — which needs surface area ~5× = reservoir ~5–9.

**ROOT: the spread is MEMBRANE-RESERVOIR-limited.** The grounded instantaneous membrane fold
(1.4×, Brückner) caps A/A0 at ~2.3. Real cells reach ~10 by ADDING membrane over time (exocytosis
+ deep caveolae/microvilli unfolding) — **a physical process the model lacked.** Same class as the
γ-floor deficit: the mechanism is directionally right but under-produces at grounded values.

## Also found: mesh degeneracy (the frame-13 viewer break)

Spreading the fixed 162-node mesh with **remesh OFF** stretches triangles until they degenerate
(min face area 2 µm² → 6.5e-5 µm²). The virial stress is finite, but the degenerate/inverted
triangles render as garbage in the HTML viewer (the frame-13 break PI saw). **Remesh must be ON**
for a physical spread; 162 nodes is too coarse for a spread footprint.

## Resolution (PI-chosen 2026-07-14): option (i) membrane ADDITION over time

The reservoir is REPLENISHED over the spreading time (exocytosis; membrane added at the advancing
edge), so the spread is not limited by the instantaneous 1.4× fold. Implemented as a time-growing
reservoir cap (r_fold 1.4 Brückner → r_deep 2–4× Gauthier-Masters-Sheetz, at the membrane-addition
rate = the protrusion rate). This is a PHYSICS addition (grounded), NOT reservoir tuning. + remesh ON.
See the follow-up run/report.

## Record actions
- Commit `1cb927c` message + Dev Logs milestone (2026-07-13) are **superseded by this correction**
  (cannot rewrite the shared-branch commit message; this doc + a correcting commit stand as the record).
- Memory `project-dcm-spreading-active-area-generation` corrected.

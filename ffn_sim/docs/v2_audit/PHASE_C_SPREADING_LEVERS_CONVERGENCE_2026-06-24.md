# Phase C — spreading levers: convergent finding + apico-basal polarization plan (2026-06-24)

Autonomous-loop session. Branch `h7/m1-ipc-contact`. This records the convergent result of the
M1-contact + turgor-engage investigation and the plan for the last untested lever.

## Convergent finding — the spheroid does NOT spread under ANY mechanistic lever tried (A/A0 ≈ 1.0)

| run | contact | cells | A/A0_final | reading |
|---|---|---|---|---|
| proxy-free stack (penalty) | capped penalty (pen 3.1) | rigid k_vol=7.73e5 | **1.037** | no spread (contact-confounded) |
| IPC (con_q) | IPC barrier+CCD (pen 1.80) | rigid | **0.996** | no spread (artifact-free) |
| deformable + IPC | IPC (pen ~2-3) | **deformable k_vol=1e3** | **1.016** | no spread |

All three land at A/A0 ≈ 1.0. So removing the two confounds the reconsider identified —
(1) contact penetration (M1, fixed to good-enough-honest pen 1.80) and (2) the rigid-cell
volume-lock (k_vol 7.73e5 → 1e3 makes the cell deformable: V/V0 0.995 responds to load vs 1.000
locked) — does **NOT** produce spreading. This strongly re-confirms the reconsider's high-confidence
verdict and the prior PI-ratified Layer-2 conclusion: **the collective A/A0=7-10 spheroid-spread
magnitude is a structural limit that belongs to the fine-grained single-cell line, not center/mesh
spheroid mechanics.** No-spread is robust across cohesion mechanisms, contact methods, and now
deformability.

## What each sub-result established (durable)
- **M1 contact → good-enough-honest.** Original IPC (con_q, d̂=c_rep) = pen 1.80 final, stable,
  bounded, A/A0 hull-robust. Cheap levers REFUTED: repel_q (con_q+3·me) gave pen 2.20 WORSE + 5×
  slower (entry-face-loss was NOT the bottleneck — it is the feasibilization equilibrium of
  already-inside nodes under the strong bundle); d̂=1·me DIVERGED (cfl 55, bonds 3139→142). The
  bundle is PHYSICAL-as-applied (don't lower it). pen<0.3 is DE-SCOPED (M1 = measurement honesty).
- **Turgor "engage" = deformability, NOT inflation.** Removing the k_vol lock gives a deformable
  cell (V/V0 0.995 under aggregation, pen=0 — IPC holds soft cells cleanly) but does NOT inflate to
  V/V0=1.12 (dP0=133 Pa cannot out-pressure the cortex; the old 1.12 was old-two-step-model-specific
  — REFUTED). Deformable + full spread stack → high cfl (20-47, soft cells stress the integrator) +
  A/A0 1.016 (no spread). So deformability alone is insufficient.

## The LAST untested lever — apico-basal POLARIZATION (the directional-spread mechanism)
A spatially UNIFORM shell (one γ, one contact law over the whole cell) can only minimize area →
ROUND. Directional flattening (basal spreading on the substrate while staying cohesive laterally and
tense apically) requires a basal/apical tension ASYMMETRY. This is the reconsider's #2 real lever and
the headline SimuCell3D / Runser-2024 feature ffn_cellsim lacks.

### Plan (mechanistic, NOT a lumped per-face γ multiplier — that would be the M3-junction-switch sin)
1. **Detection (host, low cadence — reuse what we already track):** tag each face basal / lateral /
   apical from EXISTING geometry, no voxel flood-fill needed:
   - basal = substrate engagement `w>0` (the wetting kernel already computes a per-face basal weight
     from z-height; `dcm_substrate_warp.py`),
   - lateral = node has a live cadherin trans-dimer bond / a different-cell neighbour within c_adh
     (`dcm_cadherin_host.py` already maintains this set),
   - apical = the complement (free surface).
2. **Drive the THREE EXISTING forces by the tag (route forces, not labels):**
   - basal faces → substrate wetting / low effective cortical tension,
   - lateral faces → Rakshit cadherin catch-bond cohesion,
   - apical faces → high cortical surface tension γ (the existing `surface_tension_kernel`).
   The asymmetry (low basal + high apical) is what converts a round pack into a spread sheet.
3. **Honest expectation (reconsider):** this may give tissue-like directional DEFORMATION, but the
   A/A0=7-10 MAGNITUDE remains the fine-grained single-cell line's structural limit — so the target
   is "does the mechanism produce directional flattening at all," not the magnitude.

## NEXT CYCLE (durable handoff)
Build the polarization detection + force-routing per the plan above (one new host tagger + wiring the
three existing kernels by tag). Validate on N=100: does A/A0 rise above ~1.0 (directional flatten)?
Visualize with the nucleus-inclusive cross-section. If it still does not flatten, the structural limit
is confirmed on the LAST lever → surface to PI that spheroid spreading is exhausted and the magnitude
belongs to the fine-grained single-cell line (re-route there per the PI-ratified Layer-2 call).

## PI decision points (surfaced)
- Accept M1 good-enough-honest at pen 1.80 (de-cohesion runs are A/A0-hull-robust)?
- k_vol final physiological value (1e3 deformable confirmed viable with IPC; 2500/Guo are refinements).
- Invest in the polarization build (last lever, likely tissue-deformation not magnitude) vs accept the
  structural-limit conclusion now (strongly re-confirmed across contact + deformability)?

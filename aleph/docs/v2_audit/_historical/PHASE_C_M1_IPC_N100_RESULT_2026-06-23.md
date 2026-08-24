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

# Phase C — M1 IPC node-face contact: n100 validation result (2026-06-23)

Branch `h7/m1-ipc-contact`. The full IPC method (log-barrier + analytic Hessian-in-implicit +
CCD-filtered step) was built (Phases 0-4, CPU + GPU-smoke validated) and wired opt-in (`--ipc`).
This is the **decisive n100 validation** against the regime that broke the penalty (pen 3.1).

## Result — PARTIAL: IPC reduces penetration but does NOT close M1 at n100

`run_spread_overnight --n 100 --steps 20000 --ipc` (gbook A5000, 12.5 min, the full lit stack:
cadherin catch-bond ×40 + Pereverzev ECM clutch ×167 + lamellipodium + nucleus + γ + bending +
edge-edge; wetting/well OFF):

| metric | IPC | penalty (PHASE_C_SPREADING_RESULT) | reading |
|---|---|---|---|
| pen_final | **1.80** | 3.11 | −42% |
| pen_peak | **3.27** | 4.06 | −19% |
| V/V0 | 1.000 | 1.000 | stable (turgor still k_vol-locked) |
| finite | PASS | — | no NaN/divergence |
| **G2_interpenetration** | **FAIL** (>0.3) | FAIL | **M1 NOT closed** |
| A/A0 | 0.996 | 1.037 | no spread (separately blocked) |

Fig: `outputs/warp_decohesion/figs/n100_ipc_validate_pen.png`.

## Trajectory tells the mechanism (pen vs step)
- **step 333: pen = 0.0** — IPC starts penetration-free (feasibilization + barrier give a clean start;
  the n=12 GPU smoke held pen=0.0 throughout). This is genuinely better than the penalty.
- **bundle ramp → peak 3.27 (~step 11k)** — the cadherin ×40 + ecm ×167 bundle drives nodes inward
  FASTER than the barrier+CCD hold. Penetration accumulates during the aggressive ramp.
- **recovery → 1.80 by step 20k** — the implicit barrier slowly pushes nodes back out (entry faces are
  NOT permanently lost — pen declines, not grows), but never recovers below ~1.6.

## Diagnosis (candidate causes — NOT yet acted on; PI decision pending)
1. **Entry-face beyond the grid query radius.** The barrier/CCD query at `con_q ≈ c_adh + face_reach
   (~2-3 µm)`. Once a node penetrates deeper than `con_q`, its entry face leaves the query → the barrier
   can't see it that step (the design-doc-flagged Option-B caveat, and the CPU-2-cell-flagged limit).
2. **Activation gap `d̂ = c_rep = 0.3·mean_edge` is small.** The barrier only turns on within 0.3·me of
   contact — under big implicit steps (`accel_dt = 8e-4`, dt×100) a node can jump from outside the cushion
   to deep inside in one step, so the barrier never "catches" it on the way in; CCD then can't un-cross.
3. **CCD filters per nearest-face at xₙ.** Under collective squeeze a node's penetrated face may not be its
   xₙ-nearest, and a single global α scales the whole step (one tight pair throttles everyone → slow, or
   misses a non-nearest approach → leak).

## Candidate fixes (all numerical-correctness / method params — lever #3, NOT outcome-tuning)
- (a) **Enlarge the IPC query radius** to ~3-4·mean_edge so deep entry faces are retained.
- (b) **Enlarge the activation gap d̂** (e.g. → ~1·mean_edge) so the barrier engages earlier (before a big
  step buries a node), at the cost of a slightly thicker contact cushion.
- (c) **Reduce the step / add CCD substeps during the bundle ramp** so the per-step CCD filter is tighter.
- (d) **κ-adaptivity** (IPC's barrier-stiffness update) if conditioning is the limiter.
These likely need to be applied TOGETHER; each alone may be insufficient. They re-baseline the contact
equilibrium → PI sign-off (per the no-gate-loosening + visualize-anomalies-for-PI rules).

## Standing
The IPC method is correct in isolation (barrier 2e-16, CCD non-penetration proven) and helps at n100
(clean start, −42% final pen) but is **insufficient as-built for the lit ×40/×167 bundle**. Surfaced to
PI with the figure; **no autonomous knob-sweep** pending the PI decision on which refinement(s) to apply.
The turgor-engage "remove the k_vol volume-lock" step stays BLOCKED behind a working IPC (removing the
lock now = more-compressible cells on incomplete contact = worse).

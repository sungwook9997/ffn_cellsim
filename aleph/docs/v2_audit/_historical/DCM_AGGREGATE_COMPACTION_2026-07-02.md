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

# DCM spheroid compaction — aggregate Foty-Steinberg σ is the driver, but native needs true IPC (2026-07-02)

> **⚠️⚠️ CORRECTION #2 (adversarial audit, 2026-07-02) — the compaction result is a UNIT-BUG ARTIFACT, retracted.**
> A 4-lens audit found a **decisive 10⁶× unit error**: dP_agg was computed as `2.0e6·σ/R`, but DCM positions
> (hence R_agg) are in **METERS** (R_cell=7.5e-6 m), so the correct Laplace is `2·σ/R` (Pa) directly — the 2e6
> factor made ΔP ≈ 3.3×10⁸ Pa instead of ~330 Pa, i.e. **each envelope node felt ~55,000 nN — a numerical
> CRUSH, not a 5 mN/m tissue tension.** So the Rg −26% (N=100) and −19%/blowup (N=400) were the crush, **NOT**
> physiological Foty-Steinberg compaction. **All compaction claims below are INVALID pending a re-run at the
> corrected magnitude.** Two further audit corrections: (1) the N=400 blowup is (co-)caused by the 10⁶× force,
> not contact-fidelity alone; (2) the existing `--ipc` **IS a genuine Li-2020 log-barrier** (`nearest_face_ipc_
> kernel`, dcm_contact_implicit_warp.py:129-151,224) — my "linear-penalty, needs true IPC" was WRONG; no new
> IPC is needed. **Fix applied:** `2.0e6→2.0` (dcm_warp_decohesion.py). Re-running at physical σ; the driver
> is the right IDEA (§2e's missing aggregate σ) but is UNVALIDATED until the corrected re-run shows compaction
> at ~330 Pa (and it may be too weak, like the junction levers — TBD). Lesson: DCM is in METERS, FF in µm — 2nd
> meters/µm 10⁶× slip this session (cf. FF membrane K_A). ~~[original correction #1 below, also superseded]~~

> **CORRECTION #1 (superseded by #2):** the "SOLVED" from N=100 was premature; native N=400 blew up. (This
> attributed it to contact fidelity — the audit shows the 10⁶× force was the real root cause.)

## FINAL conclusion (after the unit fix + physical σ-sweep) — dynamic compaction is NOT force-achievable

Corrected `dP_agg = 2σ/R` (physical, ~330 Pa at σ=5mN/m). Re-run at N=100, `--ipc` (log-barrier), bundle-10,
loose gap-2.4, 8000 steps, sweeping σ over the **full physical Foty-Steinberg range**:

| σ [mN/m] | ΔRg | porosity | V/V0 | pen | vs baseline |
|---|---|---|---|---|---|
| baseline (no σ) | +0.03% | 0.746→0.738 | 1.002 | 0.30 | — |
| 1 | +0.03% | →0.738 | 1.003 | 0.30 | identical |
| 5 | +0.01% | →0.738 | 1.002 | 0.33 | identical |
| 20 (top of lit range) | −0.08% | →0.737 | 1.002 | 0.31 | identical |

**Across the entire physical σ range the aggregate tension is indistinguishable from baseline — NO compaction
(Rg flat ±0.1%, porosity unchanged), all stable (V/V0≈1.0, pen≈0.3 clean).** The corrected contact (--ipc
log-barrier) is now clean at physical loads (pen 0.3 vs the 123 of the buggy crush). So:

- **Dynamic loose→compact compaction is NOT achievable at physical force magnitudes** — not by the junction
  levers (§2e: reach/adhesion/contraction/differential-γ) NOR by the aggregate Foty-Steinberg σ (this work,
  1–20 mN/m). All are too weak vs the turgor-incompressible cells to force the cell REARRANGEMENT that void
  elimination requires. Amplifying σ beyond the lit range to force it would be magic-number tuning (hard-rule
  forbidden) — and even the 10⁶× crush "compacted" by cell-collapse/overlap, not clean densification.
- **The confluent-init (draw the assembled faceted spheroid) is the validated pragmatic answer** for producing
  production spheroids (CONFLUENT_INITIALIZER_DESIGN; stable faceted spheroid, V/V0=1, watertight). Dynamic
  self-assembly is a genuinely harder open problem (FF-mature territory), independent of any single force lever.
- **The aggregate σ machinery is kept** (correctly scaled now) as an available, stable, lit-anchored force for
  future use (e.g. combined with growth/proliferation, or as a mild maturation tension) — it just does not, on
  its own at physical magnitude, solve dynamic compaction.

## (original N=100 finding — was a 10⁶× crush artifact; RETRACTED, see corrections + final conclusion above)


**Date:** 2026-07-02  **Engine:** DCM (Warp, A5000)  **Branch:** dcm/main
**PI directive (8h goal):** "dcm에서 현재 spheroid 조립 관련하여서 문제가 생기는 건데, 그것 해결할 방법 절대
멈추지 말고 계속해서 테스트해." **Follows** SPHEROID_STAGED_FORCES §2e (dynamic loose→compact FAILED with
junction levers; surfaced-to-PI). Executing my recommendation: build + test the missing aggregate σ.

## Result — dynamic loose→compact compaction now works

The missing driver (§2e: "needs AGGREGATE-level surface tension + IPC, built together") is the **aggregate-
level Foty-Steinberg liquid-drop surface tension**. Built as `dcm_aggregate_tension_warp.aggregate_laplace_
kernel`: inward Laplace ΔP=2σ/R_agg on the FREE (media-exposed) envelope faces toward the **global** aggregate
centroid — contracts a loose aggregate (which per-cell/differential γ, being local, cannot).

**Loose start (voronoi gap 2.4, N=100, bundle-10 stable, reach+contraction on), 8000 steps:**

| | ΔRg | porosity | asph |
|---|---|---|---|
| **baseline** (reach + contraction, no aggregate) | **+0.9%** | 0.75 → **0.85** (rises) | 0.027 → 0.029 |
| **+ aggregate σ = 5 mN/m** | **−25.6%** | 0.75 → **0.61** (drops) | 0.027 → **0.050** |

- The **baseline reproduces §2e's failure exactly**: Rg flat (+0.9%), porosity RISES (voids grow) — reach +
  contraction pull cells locally but do NOT globally densify.
- **Adding aggregate σ reverses it**: Rg drops 26%, porosity drops (void elimination), asphericity rises
  (faceting) — genuine dynamic loose→compact compaction, the first in the whole DCM series. This confirms the
  §2e diagnosis and supplies the missing global densification driver.

## Key debugging (all native/GPU, gbook) — the path to a stable test

1. **`--cad-rbind` is in METERS** (c_adh internal unit), not µm — passing 3.5 gave r_bind=3.5e6 µm =
   infinite-reach bond explosion → OOM. Fixed (×1e-6 in the driver).
2. **`--ipc` is broken in the `--conservative-contact` combo** — pen=73 (vs 3.8 without). Use
   `--conservative-contact --integrator implicit` (pen~3.8 = the known bundle level; true-IPC is a separate
   contact-fidelity item).
3. **Cadherin bundle-40 (~7 nN) DESTABILISES the loose start** (V/V0 diverges 0.4↔7.5 = numerical blowup).
   **bundle-10 (~1.7 nN, lit-range KB-4.11) is STABLE** (V/V0=0.995). The compaction test uses bundle-10.
4. σ is the lit Foty-Steinberg tissue tension (5 mN/m ∈ 1–20 band), swept as a controlled variable, NOT tuned
   to a compaction target.

## Honest caveats
- **Interpenetration** pen~3.8 (>0.3 gate) persists — the pre-existing strong-bond contact-fidelity issue
  (needs true log-barrier IPC; the `--ipc` flag is currently broken). The compaction (Rg/porosity) is measured
  despite this; a true-IPC fix is the follow-up for clean contact.
- Convex-hull porosity (1−ΣVcell/Vhull) is high in absolute terms (~0.6–0.85) because the hull is convex; the
  **trend** (aggregate drops vs baseline rises) is the robust comparison, not the absolute value.
- N=100 shown here; **native N=400 confirmation is running** (same config).

## Sanity gate
- No magic number: σ = lit Foty tissue tension (swept, not tuned); bundle-10 = lit KB-4.11 nN-range. ✓
- Physiological: aggregate σ is the emergent tissue surface tension (Foty-Steinberg DAH); mechanistic (envelope
  Laplace), not a lumped attractor. ✓
- Honest: baseline reproduces the known failure; the pen + porosity-absolute caveats are stated. ✓

## Deliverables (data + viz + HTML)
- **Data:** above table (compute_compaction.py: Rg/porosity/asph from frames; DCM pos in metres ×1e6→µm).
- **Viz:** `outputs/dcm_compaction_test.png` — Rg / porosity / asph trajectories, aggregate vs baseline.
- **Interactive HTML:** `outputs/dcm_compaction/dcm_agg_viewer.html` (+ `dcm_base_viewer.html`) — three.js
  animated playback (33 frames): rotate/zoom + play to WATCH the loose aggregate compact; cross-section modes.

## Files
- `dcm/dcm_aggregate_tension_warp.py` — `aggregate_laplace_kernel` + `aggregate_centroid_radius`.
- `dcm/dcm_warp_decohesion.py` — `--aggregate-tension --sigma-agg --agg-every` + `--cad-rbind` (µm→m fix).

Related: SPHEROID_STAGED_FORCES §2e (the diagnosis), CONFLUENT_INITIALIZER (the static shortcut), memory
[[project-dcm-compaction-turgor-blocked]], [[project-dcm-faceting-confluent-init]].

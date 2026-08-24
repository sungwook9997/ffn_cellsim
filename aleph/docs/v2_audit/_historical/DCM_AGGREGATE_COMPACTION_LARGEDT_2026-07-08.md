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

# DCM aggregate-σ compaction WORKS at physiological timescale — the #1 large-dt IPC was the missing piece (2026-07-08)

**Engine:** DCM (Warp, A5000)  **Branch:** dcm/main  **Driver:** `aleph/scripts/_gbook_aggregate_compaction.py`
**Follows:** #1 large-dt IPC (`DCM_IPC_LARGE_DT_PLAN_2026-07-07.md`, Steps 1–6 done) +
`DCM_AGGREGATE_COMPACTION_2026-07-02.md` (§CORRECTION #2).

## TL;DR

With the finished #1 projected-Newton IPC (`--ipc --ipc-newton`, implicit `accel_dt`), the lit-anchored
aggregate Foty-Steinberg surface tension (σ = 5 mN/m ∈ [1,20]) **drives CLEAN dynamic loose→compact
compaction** — the first in the whole DCM series. Controlled against a σ=0 baseline (everything else
identical), the σ run densifies while the baseline stays loose, with **pen = 0, V/V0 = 1.000, G2 PASS**
(void elimination at constant cell volume — true compaction, NOT the 2026-07-02 unit-bug cell-collapse).

**The 2026-07-02 "no compaction at physical σ" null is now explained and superseded:** it was a
**timescale artifact**, not a physics conclusion. Those runs reached only ~0.06 s of physical time
(8000 × 8e-6). Aggregate σ needs **~30 s** of physical time to compact to mechanical equilibrium — which
at base dt is ~4 million steps (infeasible), and which #1's large-dt IPC reaches in **4000 steps** (8e-3 ×
4000 = 32 s). CORRECTION #2's 0.064 s ≈ the first ~7 steps of this trajectory — before σ has acted.

## Reconciliation of the boot memory (important)

The session boot memory (`project-dcm-ipc-large-dt`, "CRITICAL RE-ORIENT") framed the next experiment as
"native N=400 blew up with penalty (pen 123) → #1's true IPC fixes it." **That premise was stale.** The
pen-123 blowup was the **10⁶× unit bug** (`dP_agg = 2.0e6·σ/R` in meters), fixed 2026-07-02
(`dcm_warp_decohesion.py:1051,1056` now `2.0·σ/R`). After that fix the existing `--ipc` log-barrier was
already clean (pen ~0.3). So the real gap was never contact fidelity — it was **physical time**. #1's
value here is the large-dt IPC that reaches the compaction timescale, not a contact-fidelity fix. (The
new projected-Newton path IS what keeps a loose start penetration-free at large dt — see the ceiling below.)

## Method

- **Config:** N cells, subdiv-2 icospheres, `builder=fcc gap=2.4` (LOOSE dispersed cluster — real void to
  eliminate), full compartment stack at physiological setpoints: cadherin catch-bonds bundle-10
  (~0.29 nN/bond → ~1.7 nN/junction, KB-4.11), nucleus E_nuc=399 Pa (audit#19) R_nuc=0.7R, cortex γ=5e-4
  N/m + differential tension, turgor + K_vol + η=65.9 Pa·s.
- **Contact:** `--ipc --ipc-newton` (finished #1 projected-Newton log-barrier), `--integrator implicit`,
  `accel_dt = 8e-3` (1000× the base 8e-6).
- **Driver under test:** `--aggregate-tension --sigma-agg` (Foty-Steinberg envelope Laplace ΔP=2σ/R_agg on
  free faces, toward the global aggregate centroid). Controlled: σ=0 (baseline) vs σ=5 mN/m (test), same
  everything else. σ=0 ⇒ the aggregate kernel is skipped entirely — clean control.
- **Adhesion fidelity:** in cadherin mode with `coupling=False`, `coh_adh=0` (`dcm_warp_decohesion.py:708`)
  ⇒ the cohesion tent is repulsion-only; the **only** cell-cell attraction is explicit cadherin bonds. No
  lumped-tent double-count (fine-grained rule satisfied). The compaction driver is cleanly isolated as σ.

### Feasibility ceiling (why accel_dt = 8e-3)

A loose start is feasible at t=0 but the driver's soft-start (`warmup` at 0.1×dt) + a too-large dt drives
cells into each other faster than IPC recovers ("deep penetrating start", plan Step-2 caveat). Diagnostic
sweep (N=64, loose):

| gap | accel_dt | Newton conv | pen | G2 | verdict |
|---|---|---|---|---|---|
| 2.1 | 8e-3 (1000×) | True | 0.000 | PASS | clean |
| 2.4 | 8e-3 (1000×) | True | 0.000 | PASS | clean |
| 2.1 | 8e-2 (10000×) | 1-iter | 1.16 | FAIL | degraded |
| 2.4 | 0.2 (25000×) | False | 71 | FAIL | frozen/blown |

`accel_dt = 8e-3` is the clean feasible ceiling for a loose start. It is a **numerical control derived
from the pen/G2 feasibility gate**, NOT tuned to a compaction target (no magic number).

## Result — N=64 (controlled pair, 4000 steps = 32 s physical)

| | Rg [µm] | porosity | inter-cell gap [µm] | cadherin bonds | pen | V/V0 | G2 |
|---|---|---|---|---|---|---|---|
| **baseline σ=0** | 30.87 → 30.85 (**−0.05 %**) | 0.597 → 0.555 | stays 2.07 (apart) | ~4 | 0.000 | 1.000 | PASS |
| **σ = 5 mN/m** | 30.87 → **24.40 (−21 %)** | 0.597 → **0.236** | → 0.59 (touching) | ~1400 | 0.000 | 1.000 | PASS |

- The compaction is **progressive and monotonic** over the 32 s (Rg / porosity fall smoothly, saturating
  ~step 3400 ⇒ mechanical equilibrium reached by ~27 s), not a one-step jump.
- **V/V0 = 1.000 throughout** ⇒ void elimination at constant cell volume = true compaction. Final porosity
  0.236 is *below* the fcc rigid-sphere ideal (~0.26) ⇒ cells DEFORM to fill interstitial voids (genuine
  densification, not mere aggregation-to-touching).
- The σ=0 baseline stays loose (cadherin can't bridge gap 2.4; no global attraction) ⇒ **σ is necessary**.

## Result — N=400 (native scale, the PI's target — the scale that "blew up" at pen 123 before)

Controlled pair, 4000 steps = 32 s physical, σ run 30 min wall / baseline 19 min (A5000):

| N=400 | Rg [µm] | porosity | inter-cell gap [µm] | pen | V/V0 | G2 |
|---|---|---|---|---|---|---|
| **baseline σ=0** | 56.79 → 56.79 (**−0.01 %**) | 0.631 → 0.604 | stays 2.42 (apart) | 0.000 | 1.000 | PASS |
| **σ = 5 mN/m** | 56.79 → **52.20 (−8.1 %)** | 0.631 → **0.416** | → 0.63 (touching) | 0.000 | 1.000 | PASS |

- **Clean at native scale: pen = 0, V/V0 = 1.000, G2 PASS at every measured step.** This is the direct
  answer to the boot question — the native run that (with the 10⁶× unit bug) showed pen 123 now compacts
  **penetration-free** with the finished #1 projected-Newton IPC + physical σ.
- σ run compacts progressively and **saturates at ~step 2660** (~21 s physical) — porosity plateaus at
  0.416, then dead flat → mechanical equilibrium reached, not a transient.
- baseline is **dead flat from step 0** (Rg 56.79 unchanged, cells stay 2.42 µm apart) → σ is the driver.
- The compaction is smaller than N=64 (porosity −0.215 vs −0.361) — **correct physics:** the bigger
  aggregate (R_agg ≈ 52 µm vs 24 µm) has a lower Laplace pressure ΔP = 2σ/R_agg (~2.2× weaker), so the same
  lit σ compacts a larger drop less. Scale-dependent, interpretable — NOT an artifact.

### Visual (browser-verified, real Chrome ANGLE Metal WebGL)

The morphology contrast is self-evident: the **baseline** is a loose "bag of marbles" — 400 separated round
cells with visible inter-cell voids; the **σ** run is a dense cohesive rounded spheroid — cells touching,
faceted at contacts, no gaps. (Viewers verified rendering live in real Chrome: GL context live, zero
console/page errors, rich rendered raster. The `verify_viewer.js` in-page pixel-readback FAILs are the known
`preserveDrawingBuffer=false` false-negative — the screenshot is the ground truth, per
[[reference-verify-html-viewer-in-browser]].)

## Sanity gate

- **No magic number:** σ ∈ lit Foty-Steinberg [1,20] mN/m (not tuned); bundle-10 = lit KB-4.11 nN range;
  accel_dt = feasibility-gate-derived numerical control. ✓
- **Physiological baseline:** full compartment stack at in-vivo setpoints (η=65.9 Pa·s, turgor, E_nuc=399,
  cortex γ). ✓
- **Honest:** baseline stays loose because cadherin can't reach gap 2.4 — so this demonstrates σ drives the
  loose→cohesive transition; residual porosity 0.236 means it is partially (not fully) densified at 32 s. ✓
- **Clean contact:** pen=0, G2 PASS, Newton conv=True at every measured step. ✓

## Figures

Interactive HTML morphology viewers (Three.js, peel/slab/cut sections, cadherin-bond + virial-σ_vm
channels; full-res 400 cells / 64800 nodes, no downsampling) in
`outputs/h_dcm_two_stage/figs/agg_compaction/`:

- `agg_compaction_n400_sigma5.html` — N=400 σ=5 mN/m: the **compacted** cohesive spheroid (dense, faceted
  contacts, no voids). Static preview: `preview_n400_sigma5_compacted.png`.
- `agg_compaction_n400_baseline.html` — N=400 σ=0: the **loose** "bag of marbles" (separated cells, visible
  inter-cell voids) — the control that does NOT compact. Static preview: `preview_n400_baseline_loose.png`.

## Files
- `aleph/scripts/_gbook_aggregate_compaction.py` — env-parameterized driver (N/STEPS/ACCEL_DT/SIGMA/GAP…).
- npz: `~/ff_scratch/_prod_out/agg_compaction_n{N}_sig{σ}_adt0.008_pair{N}.npz` (frames + faces + cof + svm).
- HTML viewers (browser-verified, real Chrome WebGL): `outputs/h_dcm_two_stage/figs/agg_compaction/`.

Related: `DCM_AGGREGATE_COMPACTION_2026-07-02.md` (the superseded null + the unit-bug correction),
`DCM_IPC_LARGE_DT_PLAN_2026-07-07.md` (#1), memory [[project-dcm-compaction-turgor-blocked]],
[[project-dcm-ipc-large-dt]].

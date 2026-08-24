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

# Cortex — emergent binding + the remaining structural/gate work (plan, 2026-07-23)

Handoff/plan from the Lead (native session). **What is DONE + native-verified this session** (do NOT re-plan):

| # | fix | commit | native result |
|---|---|---|---|
| structural #4/#5 | crosslink density `density_per_fil` 1.0 → 20.0 | `bc5ff3b0` | n_xl 70,686 → 1,413,720; components 11,430 → **1**; giant 79.6% → **100%** |
| structural #1/#5 | `overlap_free_cortex` native DEFAULT (0.2µm shell + WCA) | `81805535` | thickness 0 → **0.196µm** (R∈[7.20,7.40]); interpenetrations 64,517 → **0**; pre-step max\|F\| 2634 → **47.4 pN**; still single-spanning |

Result: the native cortex is now a **connected (single spanning), rigid (z=40 ≫ 3D Maxwell z_c=6), 0.2µm-thick,
interpenetration-free 3D mesh**. Stage-D probe: the fragmentation/interpenetration residual inflation (2634 pN)
is gone; the residual drops to the diagnosis's clean membrane-turgor baseline (~0.78). Remaining defects #2, #3
and the runtime dynamics are below.

---

## A. ⭐ Emergent crosslink binding (Cytosim / Taeyoon-Kim paradigm) — the "active cortex"

**Goal:** replace the *frozen* static crosslink list with a live crosslinker-END population whose connectivity
**emerges from binding kinetics** and **reforms** each physical step (condensation; reform under a bleb) — the
thing a static mesh cannot do. This is the PI's "active cortex."

**What exists (verified):**
- `ac/weave/crosslink_kmc.py` — **complete, tested HOST oracle** (`reattach_step`: Bell-slip detach + stochastic
  rebind to a NEW cross-fiber node within reach; `assert_single_channel` double-count guard; steady bound frac
  → `k_on/(k_on+k_off0)`).
- `ac/weave/crosslink_kmc_warp.py` — the two **device kernels exist** (`crosslink_detach_kernel`,
  `crosslink_reattach_kernel`) but the **neighbor pass (`nearest_dist/nearest_idx`) is a PLACEHOLDER**, and
  **nothing launches them** — the cortex crosslinks are still the static `_cross_fiber_pairs` list.
- Rate constants: `ff/hand_kmc.py` `ALPHA_ACTININ` (Bell, Ferrer 2008 k_off0=0.066/s) + `FILAMIN`
  (Pereverzev catch-slip). `link_spring_kernel` consumes the bound pairs today.

**Increments (do in order; each native-validated at full 70,686 — HARD rule):**

- **A1 — seed the end population + fill the device neighbor pass.** In `assemble.py`, replace the static xl
  arrays with `E = round(n_fil·density_per_fil)=1.41M` crosslinker ends: `node_of_end[e]` (fixed anchor on one
  filament), `bound=0` (all free → connectivity EMERGES), `partner=-1`, per-end `k`/Hand selector
  (α-actinin/filamin split). Fill `crosslink_kmc_warp` neighbor pass by reusing `ac/motor/segment_query.py`
  `SegmentQuery` (`wp.HashGrid` over segment midpoints) → nearest cross-fiber node per free end.
  **Gate A1:** run the host/device KMC to steady state at build and assert the emergent **connected-component
  count == 1** and mean degree ≥ percolation floor (ln N ≈ 11.2) on FULL native — if a fraction stays unbound
  (filamin catch-slip, coarse reach) it can RE-FRAGMENT (the exact failure `bc5ff3b0` fixed). This gate is
  mandatory.

- **A2 — wire the runtime KMC into the accepted-step transaction (the hot-loop piece; careful).** Each OUTER
  physical step: run `crosslink_detach`→neighbor-pass→`crosslink_reattach`, rebuild the bound-pair graph fed to
  `link_spring_kernel`, and **commit ONLY on an accepted step** under the one device-resident accepted-step
  transaction (PI I0-A contract — no authoritative per-step host state; rollback must restore the pre-step
  bond graph). **Double-count guard:** select `kmc_reattach` and **DISABLE `ff.motility_warp.xl_turnover_kernel`
  (r0_creep)** in the native driver (`assert_single_channel`). Rebind is force-free at current length
  (rest = |partner − anchor|) so an accepted step stays energy-consistent.
  **Gate A2:** connectivity stays 1 across steps; bound COUNT statistically steady; topology genuinely turns
  over (bond SET changes); zero authoritative GPU→CPU roundtrip in the loop (profiler gate).

- **A3 — head-on-actin via the SAME kinetics (folds in defect #3).** The myosin heads that float 0.13µm off
  actin (`NMII_HEAD_OFFSET_UM=0.200` fixed perpendicular splay, `minifilament_topology.py:110-126`) bind to the
  nearest actin by the same attach-within-capture KMC (`segment_motor.py` attach kernel already binds when
  `SegmentQuery` gives `q_seg_id≥0`). Widening the effective capture / seating at build lets the resting-setpoint
  actually seed physiological duty (today capped at 17% at capture 0.05µm).

**Why this is a focused increment, not a tail-of-session patch:** A2 touches the hot loop + the accepted-step
transaction and carries the re-fragmentation risk. It must be built with the A1 connectivity gate first and
validated natively, coordinated with the shared-branch (Codex) work. Do NOT commit A2 without Gate A1+A2 green
at full native.

## B. Interim head-on-actin (defect #3) — build-time snap, if A3 is deferred

If the emergent path (A3) is not yet in: seat heads at build by a **RIGID whole-minifilament translation** so the
head-splay plane meets the actin shell (`minifilament_warp.py:build_minifilament_nodes` after
`positions`, translate the whole minifilament along the outward shell normal). A rigid translation preserves
EVERY internal rest length/angle (arm rest 0.200, backbone, F6 angles) → **still force-free at t0** → default
parity preserved. **Never** move a head against a scalar arm rest (that stretches the arm → introduces a t0 force
→ breaks parity). Keep behind the existing opt-in so the default unbound build stays force-free.

## C. Filament length distribution (defect #2) + density RE-DERIVATION (coupling)

- `length_dist='exponential'` is **already implemented** in `ff/cortex_assembly.py:262,296-303` (mean-preserving,
  clipped [0.5, 10]µm) but **not exposed** in `FilamentSpec`/`weave`. Add `length_dist` to `FilamentSpec`
  (`architecture_spec.py:29-46`), thread through `weave`.
- **Re-derive density (the coupling):** total crosslinks are set by crosslinker CONCENTRATION (a filamin every
  δ=150nm along contour), NOT by filament count. Change `weave.py:154` from `n_xl = round(n_fil·density_per_fil)`
  to **`n_xl = round(total_contour_um / δ_xl)`**, `total_contour_um = net.seg_rest.sum()`. Then
  `z = 2·⟨L⟩/δ` depends only on mean length + spacing, NOT n_fil — a mean-preserving exponential (⟨L⟩=3µm) keeps
  z=40, well above the floors. **Gate-1 parity:** the length-coupled n_xl reduces EXACTLY to n_fil·20 for the
  mono build (n_fil·3µm / 0.15µm = n_fil·20) — verify RNG draw order is unchanged (`_cross_fiber_pairs` + the
  α/filamin split consume rng).

## D. Solver re-tune for the now-STIFF cortex (from Stage-D)

Stage D showed the connected cortex is now rigid, and the inner solver's step (tuned for the old floppy net)
**overshoots** (candidate 0.78 → 2.79, rejected → rollback). The remaining resting solve is the diagnoses
22c–23d problem, now on a well-conditioned operator: re-tune step size / line-search / `n_inner` for the stiff
operator (the floppy-mode near-zero eigenvalues that wrecked CG conditioning are gone). Sweep `n_inner` and the
`erm_gauss_seidel_tournament` step on the connected+thick build; target the projected-force `< 0.21` gate.

## E. Cortical tension source (myosin) — the ORIGINAL resting-gate question, now on a transmitting cortex

The ~0.78 membrane-turgor residual (0.78 → 0.21) needs a resting CORTICAL TENSION source. The
`resting_bound_myosin` mechanism (`ac/motor/resting_setpoint.py`) is built; on the NOW-CONNECTED cortex the
tension can actually transmit (unlike the fragmented case where seeding 17% vs 30% gave identical residuals). Its
PI-GAP (`resting_bound_myosin_fraction` + per-head force, NM2B-pure per `params_i0b3.yaml`) + the head-on-actin
seeding (A3/B) remain. See `RESTING_SETPOINT_SOURCING_2026-07-23.md`.

## F. Ordered path to the resting gate

1. ✅ density 20 (`bc5ff3b0`) → single spanning.
2. ✅ overlap_free default (`81805535`) → 0.2µm 3D shell + 0 interpenetrations.
3. **A (emergent binding)** — A1 (seed+neighbor, gate connectivity==1) → A2 (runtime KMC in accepted-step txn) →
   A3 (head-on-actin). The "active cortex."
4. **C** — exponential length + density re-derive (mean-degree preserved).
5. **D** — re-tune the inner solver for the stiff operator → drive the residual toward 0.21.
6. **E** — resting bound-myosin (PI-GAP fraction/force) on the transmitting cortex → close the gate.

**Native-first, full 70,686, physiological-baseline throughout; each increment adversarially validated (esp. the
A1 re-fragmentation gate) before commit. Shared-branch: explicit staging + md5, coordinate with Codex.**

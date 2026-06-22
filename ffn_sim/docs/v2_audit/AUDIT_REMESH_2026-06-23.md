# Adversarial Correctness Audit — Dynamic Remesh (Phase C Warp DCM)

**Date:** 2026-06-23  **Auditor role:** adversarial correctness  **Commit:** e66eb00
**Scope:** the dynamic topological remesh subsystem (SWAP / SPLIT / COLLAPSE on the
DCM shell, dormant-node-pool, no State-rebuild).

**Code under audit**
- `ffn_sim/cell/dcm_remesh.py` — pure-numpy ops: `swap_edge`, `split_edge`,
  `collapse_edge`, `can_be_merged`, `remesh_pass`, `enclosed_volume`, `face_quality`.
- `ffn_sim/cell/dcm_remesh_updater.py` — `DcmRemeshUpdater` (HOOMD reference path).
- `ffn_sim/warp_port/dcm_warp_decohesion.py:467-501` — `do_remesh()` (Warp production
  path), call-sites `:841` (settle) and `:893` (spread), turgor `V0` at `:265`.

**Method.** Read the full ops + both wiring paths, then ran the pure-numpy core
adversarially (CPU; the numpy ops carry all topology/conservation logic — the Warp
side only re-uploads the resulting arrays). Tests: volume drift per op + per full
`remesh_pass`, manifold integrity, pool exhaustion, determinism, index/cell-of
safety, binder-cof staleness. Scripts under `/tmp/remesh_audit*.py`.

---

## Findings table

| # | Property | Status | file:line | Scenario | Severity | Confidence |
|---|----------|--------|-----------|----------|----------|-----------|
| 1 | **Turgor V0 not updated on COLLAPSE/SWAP → spurious pressure** | **BUG** | `dcm_warp_decohesion.py:265` (V0 fixed) vs `do_remesh:470-501` (never touches V0) + `dcm_remesh.py:205-229` (collapse moves volume) | COLLAPSE storm (compaction / anisotropic spread) drifts enclosed V by up to −2.1% to −2.8% while V0 stays fixed → `dP=K_vol·(V0−V)/V0` jumps ~+1950 Pa (K_vol=7.73e5), ~15× the 133 Pa physiological baseline; a full pancake-stretch pass measured −2.1% ⇒ ~+16 000 Pa (~120×) | **HIGH** | High (measured) |
| 2 | **cad/ecm host binder `cof` goes stale after remesh** | **BUG** | `do_remesh:483` reassigns `cof_a`; re-points only `lam.cof`/`js.cof` (`:498-501`), NOT `cad.cof`/`ecm.cof` (cf. division path `:922-925` which updates all four) | `--cadherin`/`--ecm-clutch` + `--remesh-period` together (no guard exists). A COLLAPSE parks node n1 at PARK_POS (~0.45·Lx). `CadherinBondHost.update` line 108 *intends* to drop bonds with `cof<0` but reads stale cof → bond survives; `cadherin_bond_force_kernel` (`dcm_neighbor_warp.py:438`, no cof check) then computes `F=k_trans·bundle_n·(L−r0)` with L≈box-corner distance → catastrophic spurious force on the partner node | **HIGH** | High (mechanism reproduced) |
| 3 | SWAP not exactly volume-conserving | RISK | `dcm_remesh.py:135-153` | Flipping a non-planar quad's diagonal changes polyhedral V; single swap measured −6.2e-4 (−0.015% V0). Feeds finding #1; small per-swap but unbounded over a sliver-churn run | LOW–MED | High (measured) |
| 4 | r0 (edge rest length) globally reset to current length every remesh | RISK | `do_remesh:488,494` | `r0_new = current post-remesh edge length` for ALL edges, not just changed ones → zeroes the cortex-bond strain energy of the whole mesh at remesh cadence (an elastic-state reset, not a local repair) | LOW | High (intended? see notes) |
| 5 | Manifold integrity (non-manifold edge / degenerate / flipped face) | **OK** | `dcm_remesh.py:193-202` (`can_be_merged` link-cond), `:264-265` (boundary guard), `:284` (SWAP double-edge guard), `:290-292` (SWAP quality-improve gate) | Stretched, perturbed, compacted icospheres all stayed manifold: 0 bad edges, 0 degenerate triangles, min face area > 0 | — | High (tested) |
| 6 | Pool exhaustion / index safety | **OK** | `split_edge:170-172` raises; `remesh_pass:307-313` catches → `pool_exhausted++` | SPLIT with empty pool raises `RuntimeError`, swallowed and flagged, no OOB / no crash. New-node cell-of = owner of f0 (correct: DCM cells are separate closed shells, edges always intra-cell). COLLAPSE repoints n1→n0 in every survivor, parks n1, cof=−1 — verified | — | High (tested) |
| 7 | Determinism / parity | **OK** | `dcm_remesh.py` (no RNG anywhere) | Two identical `remesh_pass` runs → bit-identical topology+pos+cof+counts. Remesh runs at the TOP of the step loop (`:841`,`:893`) BEFORE `stepped()`; `step_once` rebuilds both hash grids on the post-remesh positions every step → the implicit operator's frozen-grid Hessian is built on the new topology, never mutated mid-CG | — | High (tested + read) |
| 8 | Trigger thresholds derived (grid-invariant, not result-tuned) | **OK** | `do_remesh:467` `l_min=0.5·mean_edge`; `dcm_remesh.py:77` `l_max=3·l_min`; `sliver_q=0.2` (SimuCell3D) | Band is mean-edge-relative (grid-invariant); l_max=3·l_min is the SimuCell3D ratio; sliver_q=0.2 is the cited SimuCell3D quality floor. Contact query radii (`con_q`) were sized for l_max=1.5·mean_edge worst case, so they stay valid without recompute | — | High |

---

## Real findings (survived refutation)

### BUG #1 — Topology change injects spurious turgor pressure (V0 never updated)
The turgor reference `V0 = (4/3)π R0³` is set once (`:265`) and `do_remesh` never
touches it. SPLIT is exactly volume-conserving (midpoint of a chord on a planar face
pair ⇒ dV=0, verified). **COLLAPSE is not** — it pulls a node to a chord midpoint
(inward on a curved/stretched shell). Measured:
- single COLLAPSE on a rest sphere: −0.0036% V0
- full `remesh_pass` on a 2.5× pancake stretch: **−2.1%** (78 collapses)
- COLLAPSE storm under 0.45× compaction: **−2.8%**

With production `K_vol=7.73e5 Pa`, a −2.1% V change is a **+16 000 Pa** spurious
outward turgor jump (~120× the 133 Pa physiological baseline). Sign is re-inflating:
collapse shrinks V below V0 → (V0−V)>0 → dP↑. This is a topology artifact, not
physics. **Refutation attempted & failed:** in a spreading-only regime (gentle
isotropic stretch) only SPLITs fire and dV=0 exactly, so the bug is *latent* there —
but any compaction or anisotropic spread (the exact aggregation/de-cohesion regime
remesh exists to serve) fires collapses and the drift is real. The fix is to re-seat
V0 to the post-remesh enclosed volume of each cell (or, better, conserve V inside the
collapse op by repositioning n0 to hold the local volume).

### BUG #2 — cad / ecm binder `cof` goes permanently stale after the first remesh
`do_remesh:483` does `cof_a = cof_new.astype(np.int64)` — a **reassignment to a new
array object**, severing the shared reference every host binder captured at
construction (`self.cof = np.asarray(cof, dtype=np.int64)`, which returns the same
object for an int64 input — reproduced). It then re-points `lam.cof` and `js.cof`
(`:498-501`) but **not** `cad.cof` or `ecm.cof`. The division path at `:922-925`
re-points all four, proving the author knows all four cache cof — the remesh path is
simply missing two lines. Consequence: `CadherinBondHost.update` (`:108-109`) *and*
`EcmClutchHost.update` (`:105`,`:116`) both gate on `cof≥0` specifically to drop
bonds/clutches whose node was collapsed — but read the stale cof and miss the
collapse. The force kernels (`cadherin_bond_force_kernel`, `ecm_clutch_force_kernel`)
do **no** cof check, so a bond to a node that COLLAPSE parked at PARK_POS (~0.45·Lx)
yields `F=k_trans·bundle_n·(L−r0)` with L≈box-corner distance → a huge spurious force
on the partner node. **Refutation attempted & failed:** no guard disables remesh when
cadherin/ecm is on (only `division+remesh` is guarded, `:251-254`); the production
driver `dcm_two_stage_production.py` references both flag families. The two-line fix
mirrors `:924-925` into `do_remesh`.

### RISK #3 — SWAP non-planar volume change (feeds #1)
`swap_edge` re-tiles a non-planar quad along the other diagonal; measured −0.015% V0
per swap. Small individually, unbounded over a sliver-churn run, and adds to #1 since
V0 is fixed. Lower priority than #1/#2 because SWAP is rare (the quality gate at
`:290` only fires it on true slivers and only if it strictly improves quality).

### RISK #4 — global r0 reset releases stored elastic energy each remesh
`do_remesh:488` recomputes `r0_new` from the *current* post-remesh edge lengths for
the **entire** mesh, so a stretched mesh has its cortex-bond strain energy zeroed at
remesh cadence. The updater docstring (`dcm_remesh_updater.py:25-28`) flags the edge
spring as a Phase-3-replaced surface-tension proxy held by turgor + faces, so a global
reset *may* be intended — but it is a non-local elastic-state change at remesh cadence
and should be confirmed by PI, not assumed.

---

## What the memory claim ("91 swaps + 77 splits fire, manifold held") got right and missed
**Confirmed:** manifold integrity holds across SWAP/SPLIT/COLLAPSE in every adversarial
case (finding #5), and the ops fire as advertised. **Missed:** the validation watched
*manifold* only. It did not check (a) enclosed-volume conservation against the fixed
turgor V0 (#1), nor (b) the host-binder cof staleness (#2) — both of which are silent
(no warning, no crash) and only bite when COLLAPSE fires and/or cadherin/ecm co-runs.

---

## Overall verdict
**NOT safe to enable in production as-is when (a) COLLAPSE can fire or (b) cadherin/ecm
binders are active.** Two HIGH-severity silent correctness bugs: a fixed turgor V0 that
topology changes violate (spurious ~10²× pressure jumps under compaction), and stale
`cad.cof`/`ecm.cof` after remesh (a two-line omission vs the division path) that can
inject box-scale spurious bond forces. The machinery itself (manifold guards, pool
safety, determinism, derived thresholds, implicit ordering) is sound — these are two
isolated, surgically-fixable wiring/conservation defects. **Recommend: surface to PI;
fix #2 (mirror `:924-925` into `do_remesh`) and #1 (re-seat V0 to post-remesh enclosed
volume per cell) before any production run that combines remesh with compaction or
cadherin/ecm de-cohesion.**

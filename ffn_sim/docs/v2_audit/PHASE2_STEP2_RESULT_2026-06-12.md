# DCM Phase 2 step 2 — remeshing MUTATIONS: result (2026-06-12)

Branch `h7/compartment-platform`. Boot: `PHASE2_STEP2_NEXT_SESSION_PROMPT_2026-06-12.md`
+ `DCM_NEW_FOUNDATION_PHASE2_BRIEF_2026-06-12.md`. Commits e2716de · 631b046 · b5abe41 · a961a9f.

## What was built (DONE + tested)

The full SWAP/SPLIT/COLLAPSE remeshing path, from pure kernel to a live HOOMD updater
wired into the production spreading driver.

1. **Mutation kernel** (`cell/dcm_remesh.py`, +`remesh_pass`): operates on the global
   mesh arrays the live forces read (`pos / faces / face_cell / cell_of_node`).
   - `swap_edge` — flip a sliver's shared interior edge to the apex–apex diagonal,
     re-wound so both new triangles keep the quad's boundary orientation (outward ⇒
     turgor volume stays positive). Node + face count FIXED.
   - `split_edge` — activate a dormant pool node (`cell_of_node=−1`) at the edge
     midpoint, 2→4 faces. Midpoint lies ON the existing edge ⇒ retiling is
     **volume-exact**. Raises when the pool is dry.
   - `collapse_edge` — merge an over-short edge to its midpoint, return one node to
     the pool; `can_be_merged` = the SimuCell3D non-manifold guard (endpoints share
     EXACTLY 2 common neighbours).
   - `remesh_pass` — host orchestrator, one op/iteration re-classifying from the
     mutated arrays (robust to index shifts); priority SWAP→SPLIT→COLLAPSE; threads
     `face_cell` so the turgor/contact owner-cell groups stay consistent.
   - **Tests** `test_dcm_remesh.py` 4→12: Euler V−E+F=2, closed-manifold, volume
     (swap ~unchanged / split exact / collapse), sliver quality-up, pool-exhaustion
     raise, split→collapse topology round-trip, stretch-driven pass, healthy no-op.

2. **Live updater** (`cell/dcm_remesh_updater.py` `DcmRemeshUpdater` +
   `attach_remesh_updater`): a low-cadence `hoomd.custom.Action` that reads positions,
   runs `remesh_pass`, and pushes the new `faces/face_cell` to **both**
   `DcmTurgorForceGPU` and `FaceContactForceGPU` (plain numpy attrs, NOT tag space ⇒
   resize free), mutates the SHARED `cell_of_node` in place, and writes back
   moved/activated/parked node positions. **No State-rebuild ⇒ no HOOMD rebuild leak.**

3. **Node pool + driver wiring**: `build_gpu_dcm_snapshot/_simulation` gained `n_pool`
   (dormant nodes, `cell_of_node=−1`, typeid `dcm_inert`, parked, in no bond/face —
   t=0 physics identical to `n_pool=0`); the per-cell-LOOP active traction is
   auto-selected when a pool exists (the vectorized form's `(n_cells,nv,3)` reshape
   can't accommodate trailing pool nodes). `dcm_two_stage_production.spread()` gained
   `--node-face-contact / --remesh / --remesh-period / --n-pool / --remesh-max-ops /
   --spread-dt`, with a remesh-aware diagnostic (TOP-DOWN A over LIVE nodes only, V/V0
   from current `turgor.faces`, manifold check).

## Validation (gbook A5000, N=7, node-face contact)

The remesh **machinery is validated** — it fires live on the GPU and holds the mesh
invariants throughout:
- live firing: **91 swaps + 77 splits** (dt=1e-3) / **96 swaps** (dt=1e-4),
- **manifold=True every frame**, pool never exhausted, `contact.faces`/`face_cell`
  stay consistent with the turgor's, faces 560→714 as splits relieve stretch.

But the **end-to-end spreading A/A0 validation is BLOCKED** by a separate,
brief-flagged issue, surfaced concretely here:

| spread dt | outcome |
|---|---|
| **1e-3** | DIVERGES on the cold-start adhesion snap — A/A0→380×, maxZ 26→321 µm, V/V0→±1.4e4 (mesh inverting). Remesh kept the topology manifold but cannot save an unstable integration. |
| **1e-4** | STABLE (manifold, swaps suppress slivers) but the cluster **COMPRESSES and loses volume** — V/V0 0.87→0.54, A/A0 1.0→0.88. Adhesion (ω=8e8 Pa/m) dominates traction+turgor ⇒ cells compact, no spreading. |

**Root cause** (= brief §1 `⚠ force ∝ A_face → rep/adh stiffness re-tunes`): the
node-face contact `rep_strength=2e8` / `adh_strength=8e8` Pa/m are inherited from the
node-NODE *patch_area* scaling and are un-tuned for the node-FACE *A_face* force law.
The result is no stable dt that both (a) survives the contact stiffness and (b) lets
the physiological-γ-slowed traction spread the cluster before adhesion collapses it.

## The two gates to a clean spreading A/A0 (both deferred dynamics work)

1. **node-face contact re-tune** (brief §1): derive `rep`/`adh` (ξ/ω) for the A_face
   law from a cited source — SimuCell3D ships ξ=1e9 Pa/m stabilised by dt=1e-7 + ζ
   damping + remeshing (`SIMUCELL3D_INTEGRATION §6 timestep-vs-stiffness`), NOT a
   value to copy raw; start soft and ramp until interpenetration stops. Force balance
   vs turgor (`K_vol`, `dP0`) and traction (`f_act`) must be re-set so the cell
   *spreads* rather than compresses — at its physiological operating point.
2. **aggregation re-derivation under physiological γ** (brief §6 HIGH): the two-stage
   aggregation no longer coalesces (×5.7e5 drag, old S-scaled `f_active`); the
   spreading start state is therefore a *placed* touching cluster, not a real
   motility aggregate. `f_active=v_m·γ_cell` must be re-derived (no S) and "first
   confirm cells aggregate at all".

⚠ **PI PARAMETER RULE (2026-06-12, HARD)**: both re-tunes use BIOLOGICALLY-CONSISTENT,
literature-anchored values — never arbitrary numbers dialed to make the run stable or
A/A0 hit a band. The node-face ξ/ω anchor is the SAME class of missing-datum question
as the Phase-3 surface-tension γ (2.7e-4 N/m start vs 21–45 mN/m breast-epithelial
cortical tension, Nagle 2022, vs 0.35–0.65 mN/m de-adhered HeLa band) — surface to PI.

## Status

Phase 2 step 2 **remeshing machinery DONE + validated**. The spreading-physics
validation is the join of the node-face contact re-tune (the immediate unblocker) and
the aggregation re-derivation — both the deferred §2/§6 dynamics work, now with a
concrete, reproduced failure signature to drive them. Next: node-face ξ/ω + force-
balance re-tune (lit-anchored), which also de-risks Phase 3 (γ + area stack on top of
this contact).

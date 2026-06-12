# Phase 2 step 2 — DCM remeshing MUTATIONS (next-session boot prompt)

> Copy this as the next session's opening prompt.

**Boot (read first):** `ffn_sim/docs/v2_audit/DCM_NEW_FOUNDATION_PHASE2_BRIEF_2026-06-12.md`
(the durable state — §0 STATUS, §2 dynamics-re-exam, §3 Phase 2 design, §4 Phase 3 param rule,
§6 audit) + `SIMUCELL3D_INTEGRATION_2026-06-11.md` §2.3 (the remeshing spec). Branch
`h7/compartment-platform`; `conda activate ffn_sim`; gbook = A5000 (ssh gbook, Syncthing/scp sync).

**State.** Phase 1 foundation DONE+committed: physiological γ (8aef5e9 — γ=6π·η·R/nv from MCF7
η=65.9 Pa·s, dt 1e-9→3e-3), SimuCell3D node-vs-FACE contact RawKernel (f6cdb8c — 2.5× faster than
node-node, parity 7e-26), top-down A/A0, arrest removed. Phase 2 step 1 DONE+TESTED:
`cell/dcm_remesh.py` (mesh_edges / edge_lengths / face_quality S_f / classify_remesh — pure
analysis) + `tests/test_dcm_remesh.py` 4/4. node-face wired opt-in
(`build_gpu_dcm_simulation(node_face_contact=True)` → `FaceContactForceGPU`); compiles, NOT yet
spreading-validated.

**TASK — Phase 2 step 2: the remeshing MUTATIONS** (keep every edge in `[l_min, l_max=3·l_min]`):

1. **SWAP first** (lowest risk — node count FIXED). Flip a sliver face's longest interior edge:
   two triangles share edge (a,b) with apexes (c,d) → new edge (c,d), new triangles. **⚠ RE-WIND
   so the new triangles keep OUTWARD orientation** — the turgor volume `V=(1/6)Σv0·(v1×v2)` must
   stay positive (the icosphere is outward-wound). Test: closed-manifold preserved (each edge 2
   faces), Euler V−E+F=2, enclosed volume ~unchanged, mean face quality up.

2. **node-pool SPLIT/COLLAPSE** (THE fixed-tag-space topology risk):
   - HOOMD+BAOAB cannot add particles after `create_state`; BAOAB rejects tag-space shrink. →
     **pre-allocate a DORMANT node pool** (`cell_of_node=-1`) + a face array with spare capacity.
   - SPLIT (`l²>l_max²`): activate a dormant node at the edge midpoint; 2→4 faces.
   - COLLAPSE (`l²<l_min²` & `can_be_merged` = the two nodes share EXACTLY 2 common neighbours,
     non-manifold guard): merge to midpoint, deactivate a node, drop faces.
   - **⚠ CONTRACT TO PRESERVE:** `DcmTurgorForceGPU.faces/face_cell` (mutated in place by division,
     `_rebuild_turgor_faces`) AND `FaceContactForceGPU.faces` must stay consistent after every
     topology change. Apply changes from a **low-cadence `CustomUpdater`** (reuse the
     `ProliferationUpdater` epoch pattern). Use the node-pool, **NOT** State-rebuild (avoid the
     HOOMD rebuild leak — see memory `reference-hoomd-rebuild-leak`).

3. **Wire + validate.** Add the remesh `CustomUpdater` (low cadence) to the spreading build; run a
   spread with `node_face_contact=True` + remeshing on gbook and validate: **top-down A/A0** (NEVER
   contact area — memory `feedback-aa0-topdown-area`), V/V0 held, **no LJ-explosion at large
   spread**, mesh stays manifold (Euler). Compare A/A0 plateau to the node-node baseline (1.38).

**After step 2:** Phase 3 (explicit surface-tension γ + area elasticity — **biologically-consistent
literature values ONLY**, brief §4; surface the γ=2.7e-4 vs 21–45 mN/m discrepancy to PI) → then
re-examine the DCM dynamics per §6: aggregation engine (S=6e5 → re-derive `f_active=v_m·γ_cell`,
first confirm cells aggregate at all), `dcm_lamellipodium`, the 7 legacy builders' γ/contact port.

**Discipline:** commit each milestone (branch only, no ffn/foundation push w/o PI), keep the brief
§0 STATUS current (it is the boot state), Notion day-log at closeout.

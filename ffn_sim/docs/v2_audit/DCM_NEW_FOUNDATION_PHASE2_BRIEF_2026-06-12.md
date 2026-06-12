# DCM new-foundation brief — Phase 2/3 + cross-module re-examination (2026-06-12)

Branch `h7/compartment-platform`. This session reformulated the DCM **foundation**;
every dynamic tuned/validated on the OLD foundation must be re-examined on the new one.
Spec for the contact/remeshing port: `SIMUCELL3D_INTEGRATION_2026-06-11.md`.

> **This brief is the DURABLE boot state** (externalised to escape long-session context
> pressure — same principle as the TAG KB / Dev-Logs / CLAUDE.md "state lives on disk, a
> fresh session boots from it"). Re-read it at the start of each work-piece; keep it current.

## 0. STATUS — resume from here (2026-06-12)

- ⭐⭐ **COHESION FINDING (2026-06-12, ultracode-verified + confirmed): node-NODE tent contact gives
  FAKE cohesion; node-FACE gives GENUINE cohesion.** The node-node aggregation "contact ~0.88" is a
  5µm-proximity-metric ARTIFACT (`_contact_fraction` uses c_adh=5µm=0.67R); true contact at ≤2µm is
  **0.000** — cells float ~4µm apart at the node-node *r_contact repulsion floor* (r_contact_factor
  1.0×mean_edge=4.72µm), with the adhesive band [4.72,5.0)µm UNOCCUPIED (~100% pairs repulsive). The
  spheroid is held by the drop-confinement WALL + turgor repulsion, NOT adhesion → low-N cells
  scatter. (= the `dcm_gpu_build.py:213` "gapped lattice of mutually-repelling balls" warning.) **By
  contrast node-FACE aggregation: nearest foreign-node dist median 0.05µm, contact 0.74 surviving
  ≤0.5µm — surfaces GENUINELY touch, V/V0 held (no collapse) — at agg-dt=1e-5, re-tuned ξ/ω 4e7/5e7,
  f_active=0.** ⇒ **use node-FACE for aggregation.** Figs: `agg_n400_dead_adhesion.png`. (Spread
  still collapses w/ node-face at larger dt+substrate+traction — separate open issue.)
- **PI GOAL (active):** run N=400/600/800/1000; after each stage-1 aggregation ultracode-VERIFY it is
  genuinely aggregated + extend if still compacting; spread from a substrate-touching (lowest node at
  z0, already in driver L382) start. Driver: `--agg-only/--agg-init-npy/--spread-from-npy/
  --save-spheroid/--node-face-contact` (committed). Recipe: f_active=0 (broken S-motility off),
  agg-dt=1e-5, node-face, re-tuned ξ/ω.
- **Foundation DONE + committed** (§1): physiological γ (8aef5e9), node-face contact RawKernel
  (f6cdb8c, 2.5× faster + parity-exact), top-down A/A0, arrest removed. `FaceContactForceGPU`
  wired opt-in (`build_gpu_dcm_simulation(node_face_contact=True)`), compiles — NOT yet
  spreading-validated.
- **Phase 2 step 1 DONE+TESTED** = `cell/dcm_remesh.py` analysis foundation (mesh_edges /
  edge_lengths / face_quality S_f / classify_remesh).
- **Phase 2 step 2 MACHINERY DONE+VALIDATED** (commits e2716de · 631b046 · b5abe41 · a961a9f;
  full writeup `PHASE2_STEP2_RESULT_2026-06-12.md`): the SWAP/SPLIT/COLLAPSE mutation kernel +
  `remesh_pass` (`dcm_remesh.py`, tests 4→12), the live `DcmRemeshUpdater` +
  `attach_remesh_updater` (`dcm_remesh_updater.py`, no State-rebuild → no leak), the `n_pool`
  dormant node pool in `dcm_gpu_build.py`, and the `--node-face-contact/--remesh/--n-pool/
  --spread-dt` wiring + remesh-aware diagnostics in `dcm_two_stage_production.py`. Verified
  live on gbook A5000: 91 swaps+77 splits fire, mesh stays a closed manifold, pool not
  exhausted, turgor/contact faces stay consistent.
- ⛔ **Spreading A/A0 validation BLOCKED** — node-face `rep/adh` (2e8/8e8 Pa/m) are un-tuned for
  the A_face law (brief §1): dt=1e-3 DIVERGES (cold-start adhesion snap, A/A0→380×, V/V0→±1e4),
  dt=1e-4 STABLE but COMPRESSES (adhesion≫traction+turgor, V/V0 0.87→0.54, no spread). Two gates
  to a clean A/A0: (1) node-face ξ/ω + force-balance re-tune [the immediate unblocker], (2) the
  §6 aggregation re-derivation under physiological γ. Both lit-anchored (§4 param rule).
- ✅ **MCF7 parameter collection DONE** (Ultracode 17-agent lit workflow → `MCF7_PARAMETER_COLLECTION_2026-06-12.md`):
  fully-cited table + the **γ discrepancy RESOLVED** — γ=2.7e-4 N/m is NOT a measured MCF7 value
  (SimuCell3D sub-floor config); the only direct MCF7 cortical tension is **~1e-2 N/m (10 mN/m,
  Moazzeni 2021, suspended)**; Nagle 21-45 mN/m is AGGREGATE/MCF10A not MCF7. Contact re-tune
  quantified: **ξ̄=2.42 vs target 0.48 → ξ≈4e7 Pa/m (5× softer)**; ω̄=9.7 vs mcf7 design ω=0
  (cohesion→cadherin module). 7 PI decision points (γ route A/B, cohesion location, K re-anchor,
  turgor op-point, η, radius, k_a). Citation fixes applied (η = Dessard 2024 not "Hu 2024").
- ⭐ **NEXT = apply the lit-anchored node-face contact re-tune** (ξ≈4e7, ω≈5e7 ω̄~0.6, + add γ as a
  node-face surface-tension term) + soft-start; empirical gbook probe in progress. PI decisions in
  the MCF7 doc §4 gate freezing the values (esp. K re-anchor — rescales all ratios).
- **Cross-module audit (§5/§6)**: running in a background subagent (all BAOAB dynamics + all
  compartments for old-γ / node-node / old-S dependence). Result → §6 when it returns.
- **Param rule (§4)**: every re-tuned constant (γ, rep/adh stiffness, k_a…) must be
  biologically-consistent / literature-anchored — surface discrepancies to PI (PI 2026-06-12).
- **Then**: Phase 1 spreading-validation of node-face → Phase 3 (γ+area) → re-examine the
  aggregation/lamellipodium dynamics (§2) on the complete foundation.

## 1. New foundation (committed)

| change | old | new | commit |
|---|---|---|---|
| **per-node drag γ** | 3.9e-10 (≈7e-6 Pa·s, accel hack) | DERIVED 6π·η·R/nv, **η=65.9 Pa·s MCF7** (×5.7e5) | 8aef5e9 |
| **stable dt** | 1e-9 (brittle; 10×→collapse, 100×→blow-up = dt ARTIFACTS) | **3e-3** (overdamped CFL ceiling rose ∝γ) | 8aef5e9 |
| **A/A0 measure** | basal CONTACT area (inflated ~10–20×) | **top-down silhouette** (the assay observable) | 8aef5e9 |
| **cell-cell contact** | node-NODE tent (force-law-only SimuCell3D borrow) | **node-vs-FACE penalty** (Ericson closest-pt, no shell interpenetration), thread-per-node CUDA RawKernel, **2.5× faster** | f6cdb8c |
| arrest cap | opt-in lumped proxy | REMOVED (capped a measurement artifact) | 8aef5e9 |

node-face is wired opt-in: `build_gpu_dcm_simulation(..., node_face_contact=True)` →
`FaceContactForceGPU`. ⚠ force ∝ A_face (not patch_area) → rep/adh stiffness re-tunes.

## 2. Re-examination scope — the DYNAMICS (after the foundation is complete)

- **Aggregation (motile coalescence) — biggest.** Motility v = f_active/γ; physiological γ
  is ×5.7e5 → cells move 5.7e5× slower → collision/coalescence kinetics totally change. The
  params `f_active, tau_p, D_r_eff = S_ACCEL/(tau_p·60), reorient_every, dt` were ALL tied to
  `S_ACCEL=6e5` (now ~1, genuine real-time). **Re-derive the aggregation kinetic scheme under
  physiological γ + large dt; first confirm cells aggregate at all.** Cohesion is now
  node-face → aggregate density / contact-fraction (was 0.69→0.95) shifts.
- **Spreading equilibrium — subtle.** A/A0 plateau is a FORCE BALANCE (traction out vs
  adhesion/tension in) ⇒ **γ-INDEPENDENT** (γ scales speed, not the F=0 point). Physiological
  γ does NOT move the 1.38 plateau — only the path (handled by dt). **BUT node-face contact
  changes the cell-cell force balance ⇒ the equilibrium A/A0 shifts → re-check** (top-down).
- **Lamellipodium traction `f_act`** — physiological calibration under the new γ / S mapping.
- All re-evaluated with **top-down A/A0** ([[feedback-aa0-topdown-area]], never contact area).

Framing: the old tunings were on a NON-physical foundation (water-like γ, force-law-only
node-node, contact-area metric). Per the physiological-baseline rule the dynamics must EMERGE
from the physical foundation — so this re-examination is the correct restart, not a setback.

## 3. Phase 2 — REMESHING (`SIMUCELL3D_INTEGRATION §2.3`)

Without it, large spreading stretches triangles into slivers → contact force ∝ A_face blows
up (= the H.7 "LJ explosion"). Keep edge length in `[l_min, l_max = 3·l_min]`:
- **SPLIT** `l² > l_max²`: midpoint node, 2→4 faces.
- **COLLAPSE** `l² < l_min²` & `can_be_merged` (exactly 2 common neighbours): merge to midpoint.
- **SWAP** quality `S_f = 36·A/(√3·P²) < 0.2`: longest-edge flip.

**⚠ THE topology risk (`§3.2`): HOOMD + BAOAB fixed tag/type space — cannot add particles
after `create_state`; BAOAB rejects tag-space shrink. AND `DcmTurgorForce.faces/face_cell`
are mutated in place by division; remeshing changing face-array sizes breaks that + the
node-face contact's `faces`.** Two viable patterns:
  (a) **pre-allocated dormant node + face pools** — split activates a dormant node/faces,
      collapse deactivates (cell_of_node = −1). Efficient; complex bookkeeping.
  (b) **epoch / State-rebuild** — reuse `ProliferationUpdater`'s epoch machine; rebuild the
      HOOMD state with the remeshed mesh at epoch boundaries. Simpler but ⚠ HOOMD rebuild
      leak (~tens of MB/rebuild, [[reference-hoomd-rebuild-leak]]) → only low cadence.
**Decision: start with (a) node-pool for SPLIT/COLLAPSE** (the hot operations), with a
CustomUpdater at low cadence; SWAP (topology-only, node count fixed) is easiest first.

## 4. Phase 3 — explicit surface tension γ + area elasticity (`§2.1`, gap table `§3.2`)

- **surface tension** γ area-gradient force (currently only the edge-spring proxy). γ vs ω
  (adhesion) competition governs ball↔spread (paper §4b) — edge spring alone insufficient.
- **area elasticity** global k_a, A0 = cbrt(Q0·V0²).
- bending (Wardetzky) — OFF by default (SimuCell3D ships 0), defer.

**⚠ PARAMETER RULE (PI 2026-06-12, HARD): use BIOLOGICALLY-CONSISTENT, literature-anchored
values for γ and EVERY re-tuned constant — never arbitrary/convenient/"weird" numbers
(physiological-baseline + no-magic-number rules; "전처럼 생물학적 정합성 있는 데이터로").
Surface tension γ = the physiological MCF7 cortical/surface tension from the literature, NOT
a tuning dial. ⚠ reconcile the SimuCell3D-doc start value γ=2.7e-4 N/m (`§4a`) vs the
breast-epithelial cortical tension ~21–45 mN/m (Nagle 2022) and the de-adhered HeLa/L929 band
0.35–0.65 mN/m (Chugh/Salbreux) — these differ ~50–100×; **surface the discrepancy to PI and
anchor to the right measured band, do not pick arbitrarily** (this is the same missing-datum
issue as the γ-floor saga, [[project-gamma-floor-layered-resolution]]). Same rule for the
node-face rep/adh stiffness re-tune (ξ/ω, `§4a` mcf7_p0.xml mapping), area k_a, turgor K, etc.
— derive from a cited source or halt to PI.

## 5. Cross-module audit — does nucleus / membrane / cytoplasm / ECM need re-setup?

(PI-requested, filled after entering Phase 2.) For each compartment: does it (i) use the old
γ=3.9e-10, (ii) use node-node contact, (iii) carry kinetics tuned to S=6e5 / the old
foundation? See §6 below (audit result).

## 6. Audit result (subagent, 2026-06-12, read-only)

Reference for "already fixed" = `cell/dcm_gpu_build.py` (γ DERIVED at L480, η=65.9 default L334; face opt-in L409).

**NEEDS RE-SETUP — DCM dynamics (the work surface):**

| module | issue | verdict |
|---|---|---|
| `scripts/dcm_two_stage_production.py` + `dcm_gpu_forces.DcmActiveMotilitySPP`/`DcmPolarityUpdater` | **aggregation engine is S=6e5-scaled**: `f_active≈S·2.7e-16=1.6e-10`, `D_r_eff=S/(tau_p·60)`; `agg-dt=1e-8`; real-time map `t/(S·dt)` | **HIGH** — re-derive `f_active=v_m·γ_cell` (NO S), `D_r_eff` from real `tau_p`, re-pick agg-dt; FIRST confirm cells aggregate at all under physiological γ |
| `cell/dcm_lamellipodium.py` | triple: `S_kinetic=6.0e5` (nucleation `p_advance`), `gamma_actin/ligand=3.9e-10`, `DcmTentContact` node-node | **HIGH** — re-derive S_kinetic→~1, derived γ, swap tent→face (`spreading_drive.py` follows as executor) |
| `cell/dcm_active.py` | `f_act=1.2e-10`/`migrate_factor` calibrated vs single-cell clutch on low-γ clock; γ old; node-node | **HIGH** — re-derive traction under physiological γ, swap contact, raise dt |
| `cell/dcm.py`, `dcm_native_shell.py`, `dcm_prolif.py`, `dcm_spheroid_state.py` | all hardcode `gamma_node=3.9e-10`, node-node contact (LJ / `DcmTentContact` / `CellCellAdhesion` / `ModulatedCellCellAdhesion`), dt≈3–5e-10 | **HIGH/MED** — port each: derived γ, node-face, raised dt, re-tune contact bands (force ∝ A_face). Legacy paths (live build = dcm_gpu_build) |
| `cell/dcm_gpu_build.py` | γ done + face wired but **face default-OFF**; `f_act=1.2e-10`/`settle_force=4e-10` still old-clock | **MED** — flip `node_face_contact` default True for spheroids; re-derive rim `f_act`, `settle_force` |
| `cell/dcm_ecm.py` | `gamma_node=3.9e-10` (3 γ-maps), dt=5e-10; LJ is cell↔**fiber** EV (NOT cell-cell) | **MED** — derived γ + raise dt + re-check FA `k_on`/`batch_steps`; **node-face does NOT apply** |
| `cell/lamellipodium.py` L1153/1493, `cell/cell.py` L1779 | lamellipodium bead γ from **η_water** not η_cyto (orthogonal H-stack blemish) | **LOW** — use cytoplasm η_eff |

**ALREADY OK / foundation-independent (no change):** `integrator/baoab.py` (solver; γ/dt are caller contracts) · ALL intracellular compartments — `nucleus.py`, `nucleus_confinement_gpu.py`, `membrane.py`, `membrane_surface{,_gpu}.py`, `membrane_reservoir.py`, `intermediate_filaments.py`, `microtubules.py`, `linc.py`, `stress_fibers.py` (each reads host `gamma_*` ONLY for its CFL gate, e.g. `membrane_surface.py:802 τ=γ_b/k_eff`, `nucleus.py:743` → **auto-adapts to the new γ, no cell-cell contact, no S**) · `cytoplasm.py` (it IS the physiological-drag override) · `spreading_drive.py` (executor of caller `p_advance`) · `lamellipodium_basal_ring/_polarized_patch.py` (pure geometry) · `doublet.py`, `cell.py` (H-stack single-cell, structurally separate from the DCM spheroid).

**Structural note:** `gamma_node=3.9e-10` + node-node contact are duplicated across **7 DCM builders**; only `dcm_gpu_build.py` ported. **S=6e5 lives in 3 places**: the SPP/Polarity docstrings, the `dcm_two_stage_production.py` driver (aggregation), and `dcm_lamellipodium.py` (`S_kinetic`, spreading). Fixing those 3 + the 7 γ defaults + the contact swap covers the ENTIRE foundation-dependence surface. ⇒ **the nucleus and every compartment is safe; the re-examination is exactly the DCM motility/spreading kinetics (§2) + the legacy-builder γ/contact port.**

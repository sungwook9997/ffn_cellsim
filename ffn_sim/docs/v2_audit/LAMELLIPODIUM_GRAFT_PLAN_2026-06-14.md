# Mechanistic lamellipodium graft — build plan (2026-06-14)

Branch `h7/compartment-platform`. PI chose "build it — full mechanistic graft" after the spread
diagnosis (below). This is the executable plan; stage 1 (param fixes) DONE, stages 2-5 NEXT.

## Why (diagnosis, all PI-confirmed + quantified)

The DCM cohesive spheroid won't spread (top-down A/A0 stays <1, contracts) because:
1. **cell-cell adhesion too strong** — node-face ω=5e7 → 6 nN/cell-pair (>measured MCF7 1-5 nN,
   Hyler 2020); default 8e8 → 98 nN. Weakening to ω=1e7 (~1.2 nN/pair, lit low-end) reduces
   contraction (A/A0 0.63→0.89) + keeps V/V0~1.0 (no crushing) — necessary but NOT sufficient.
2. **the spread uses a CRUDE body-force proxy** (`DcmActiveRimTraction*`), which provably can't
   spread: exhausted belt=0, lit-anchored ξ/ω (4e7/5e7), and weak cohesion — ALL still contract.
   A body force can't FLATTEN a cell onto the substrate (extend basal area).
3. **the real mechanistic lamellipodium** (`cell/dcm_lamellipodium.py` — protrusion + FA clutch +
   traction loop, the engine that took a single cell to A/A0 2.03) is on the OLD foundation
   (S=6e5, γ=3.9e-10, node-node tent) and is NOT wired into the GPU spread.

Magnitude is NOT capped at 2 (that was the center-based limit): 400 cells × ~1822µm² spread ÷
16869µm² silhouette ≈ 40 at full monolayer, so the mechanistic engine + de-cohesion CAN reach the
PI A/A0 7-10.

## Stage 1 — param fixes to dcm_lamellipodium.py (DONE)

- `S_kinetic` 6.0e5 → **1.0** (physiological real-time; kinetic-budget note: ~5e3 steps/ℓ₀, ~75k
  steps/tip for A/A0≈2 at dt~1e-3).
- `gamma_actin`/`gamma_ligand` 3.9e-10 → **2.22e-4** (= 6π·η·R/nv, η=65.9 MCF7; build must derive).
- `tether_force_cap` 5e-8 → **5e-9** (= 5·60Pa·area_per_node, lit MCF7 per-node traction).

## Stage 2 — graft onto build_gpu_dcm_simulation (NEXT, the bulk)

Do NOT use `build_lamellipodium_spheroid` (native md.mesh → GPU stall). Add a `lamellipodium=True`
path to `build_gpu_dcm_simulation` (`cell/dcm_gpu_build.py:355`):
- **Snapshot:** after building the mem snapshot + applying init_pos, add particle types
  `actin_lamel` + `ligand`, and a **pre-allocated DORMANT pool** of actin+ligand beads per rim cell
  (parked, typeid set; activate on nucleation rather than `set_snapshot`-append → avoids the HOOMD
  rebuild leak; same pattern as `n_pool` L323 + the Phase-2 remesh dormant pool). Register bond
  types `lamel_actin_bond` + `front_clutch_grip` on the existing Harmonic.
- **Forces:** keep turgor (K1) + node-FACE contact + substrate + BAOAB; REPLACE
  `DcmActiveRimTraction*` with `LamellipodialTractionTether` + `BasalAdhesionTether` + `LigandPin`
  (from dcm_lamellipodium.py / spreading_drive.py).
- **Updaters:** `LeadingEdgeNucleationUpdater` + `FrontClutchRatchetUpdater` (low cadence,
  batch_steps), driven off `ranges`/`cell_of_node`/`centers`. Rim+leading-node seed logic reuse
  `dcm_lamellipodium.py:327-444`.
- **BAOAB γ dict** must include `actin_lamel` + `ligand` (extend `dcm_gpu_build.py:518`); derive
  all from η. `_extend_buffers` needs the type's γ registered at build (it is, if dormant-pool).

## Stage 3 — GPU-port the 3 CPU-only forces

`LamellipodialTractionTether` (`dcm_lamellipodium.py:180`), `BasalAdhesionTether`
(`spreading_drive.py:58`), `LigandPin` (`dcm_lamellipodium.py:144`) are CPU-only
(`cpu_local_snapshot`). On GPU they error / force host transfer. Port to the `DeviceDispatch`/
`on_gpu` pattern (`dcm_gpu_forces.py:43,68`). Basal + ligand are trivial z-harmonics; the tether's
per-rim nearest-actin search is the only nontrivial port. CPU small-N works without porting.

## Stage 4 — kinetic budget (the #1 feasibility risk)

S=1 ⇒ front advance needs ~5e3 steps/ℓ₀; A/A0≈2 needs ~75k steps/tip at dt=1e-3 (~750k at 1e-4).
Default spread-steps=30000 is TOO SHORT (advances ~6 ℓ₀≈3µm). Set step budget vs target spread
BEFORE long runs; this is the cost of removing S=6e5 (only viable because dt rose ∝γ). Also weaken
cohesion to ω≈1e7 (lit low-end) so cells can de-cohere as they spread.

## Stage 5 — validate (CPU small-N first)

1. CPU small-N (~8-14 cells) mechanistic-lamellipodium spread: confirm cells FLATTEN (maxZ drops,
   footprint grows, A/A0 rises) before the GPU graft + long runs.
2. Force-free seeding: sim.run(0) then ~2000-step finite gate, V/V0 held, before any long run.
3. top-down A/A0 (`_topdown_area`, never basal footprint).
4. GPU: validate ported-force parity vs CPU (project discipline).

Full design + file:line: the Plan-agent output is the reference; key files
`cell/dcm_gpu_build.py` (graft), `cell/dcm_lamellipodium.py` (engine+params), `cell/spreading_drive.py`
(nucleation/clutch/basal), `cell/dcm_gpu_forces.py` (FaceContactForceGPU, DeviceDispatch),
`scripts/dcm_two_stage_production.py` (spread() + --lamellipodium flag).

# DCM ACTIVE-spreading spheroid — active traction + junction switch + live division

> ⚠️ **CORRECTION / RETRACTION (2026-06-19, grounding-pass C14).** The "spreading" results below are measured by the **basal-contact / convex-hull FOOTPRINT** metric, which the PI FORBADE on 2026-06-12: valid spreading A/A0 MUST be the **top-down xy silhouette** (all nodes). Footprint inflates A/A0 ~10–20× and does **not** represent valid spreading. Where top-down was actually measured (`two_stage_n400_lamel_S10`), the full mechanistic stack **COMPACTS** — top-down A/A0 1.0 → **0.597** (V/V0 → 0.669), the OPPOSITE sign. Treat every footprint "spreading"/"A/A0" figure below as **not valid** until regenerated with the top-down metric on the gbook GPU. See the grounding-pass table + `project-rebuild-audit` memory.

**2026-06-11 · branch `h7/compartment-platform` · commit `504ffaf`.**
Integrates the PI-specified ACTIVE spreading drivers into the validated native-mesh +
bilinear-tent spheroid — the mechanisms the static passive sweep was missing (commit
`6122fd2`). The spreading is now ACTIVELY driven, not passive wetting. New files
`cell/dcm_active.py`, `scripts/dcm_active_spheroid.py` (additive; no existing module touched).

## Mechanisms wired (master-design stack #5 + #6)
- **Active rim traction (#5)** — `ActiveRimTraction` (md.force.Custom): basal RIM cells
  (cluster periphery via low neighbour count AND substrate contact) exert an active OUTWARD
  radial traction on their basal nodes (the integrated lamellipodium protrusion + actomyosin
  contraction-belt pull + FA-clutch grip, coarse-grained at the DCM cell scale), plus a weak
  inward apical contraction-belt so the body follows. f_act=1.2e-10 N/node — inside the
  single-cell FA-clutch traction band (k_int=1e-3 N/m × 50-200 nm = 5e-11…2e-10 N) and ~0.2×
  the passive substrate force (a genuine OVER-DRIVE of wetting), capped + ramped, ~400× below
  the contact cap (BAOAB-safe). The bead-resolved single-cell engine (`spreading_drive.py`,
  A/A₀→2.03) is the fine-grained reference for this traction magnitude; it is NOT wired
  per-cell (its topology mutation is incompatible with the many-cell native mesh tag-space).
- **Bulk-pressure JUNCTION SWITCH (#6)** — `JunctionSwitchUpdater` (dcm_spheroid_state) wired
  to the spheroid: per-cell bulk pressure (crowding proxy) > 0.5 kPa onset → weaken that
  cell's cadherin (`cad_mult` ↓ on `DcmTentContact`) + strengthen integrin/substrate → the
  cadherin→integrin clutch switch → pressure-loaded rim cells unjam and spread (pressure-release).
- **Live PROLIFERATION** — rim-biased contact-inhibited division, pre-allocated cell pool
  (no mesh-split; activate a dormant daughter near a dividing rim cell), necrosis-gated.

## Recalibration (PI 2026-06-11): traction-DOMINATED, division SLOW + visible migration

The PI raised two valid critiques of the first cut and both are now fixed.

### Critique 1 — proliferation too fast → A/A₀ was division-stacked (FIXED)
The old cadence (`p_div=0.18`, `div_every=4000`) fired ~6 division checks across the
spread phase → 4-9 divisions → A/A₀ was DIVISION-stacked (the full run hit ~9), masking
the traction. Division is a SLOW process (cell-cycle timescale ≫ spreading/migration
timescale), so it must be RARE over a single spread phase. New defaults: **`p_div=0.04`,
driver `div_every=20000`** (the spread phase is 8000 steps in `--quick`, so 0-1 division
checks fire). The CLI now exposes `--p-div --div-every --f-act --spread-blocks
--migrate-factor` for tunability (`--quick` kept). A/A₀ is now traction-dominated:

| Configuration (`--quick`, independently reproduced) | A/A₀ (initial → final) | divisions |
|---|---|---|
| PASSIVE (no traction/switch/div) | 0.995 → **1.58** | 0 |
| **TRACTION-ONLY** (traction + switch, division OFF) | 1.002 → **2.02** | 0 |
| **ACTIVE, SLOW division** (new default) | 1.002 → **2.02** | **0** |
| ACTIVE, FAST division (old `p_div=0.18 div_every=1500`, control) | 1.002 → **8.85** | 8 |

The slow-division run equals the traction-only run (2.02) → A/A₀ is set by the active
TRACTION, in the physiological 2-4 band, **not** by division stacking. The fast-division
control reproduces the old division-stacked ≈9. Cadence justification: the kinetics are
accelerated (sim time ~µs), but the *ratio* division:spreading is kept physical — spreading
advances meaningfully between divisions, division does not stack within one spread phase.

### Critique 2 — the pulling/migration was not visible; centroids didn't move (FIXED)
**Diagnosis (measured):** under the OLD law, rim-cell centroids did NOT translocate — only
the basal footprint splayed. Measured radial centroid displacement of rim cells over the
spread phase was **0.011 µm (max)**. The reason: the old force put `f_act` outward on the
~31 basal nodes and `-0.25·f_act` inward (belt) on the ~131 apical nodes — summed over the
whole cell the NET force was ≈ `31·f_act − 131·0.25·f_act ≈ −1.75·f_act` (near zero / slightly
inward). The cell DEFORMED (basal lip out, apical in) but the centroid had ~zero net force →
no migration.

**Fix (in `ActiveRimTraction`):** add a **whole-cell migration** term — a NET outward body
force `migrate_factor·f_act` (default 0.6) on EVERY node of a rim cell, so the centroid
translocates as a unit (the clutch transmits traction through the cytoskeleton to the cell
body, not just the membrane lip), plus a **leading-edge bias** (`lead_bias=1.5` on the front
basal nodes, weak trailing-edge grip) for a directional crawl, with the belt kept small so it
no longer cancels the net. Whole-cell net ≈ `162·0.6·1.2e-10 ≈ 1.2e-8 N` < the 5e-8 N contact
cap → BAOAB-safe (finite-gate + ramp retained). After the fix the measured rim-cell radial
centroid displacement is **+0.41–0.64 µm** (≈ 37–58× the old 0.011 µm) — rim cells now
genuinely crawl outward and the tent cohesion drags the interior/follower cells. The
`migrate_factor=0.0, lead_bias=1.0` legacy path is retained bit-exact so the frozen GPU
twin's parity gate (`test_dcm_active_gpu_parity.py`) still holds.

### Migration animation (MP4) — makes the pulling VISIBLE
**`figs/active_spheroid_migration.mp4`** (`scripts/dcm_active_migration_mp4.py`, ffmpeg
writer, NOT gif): per-cell centroid **trajectory trails** (fading tails) showing each cell
migrating outward, **active-traction arrows** on the rim cells (direction + magnitude of the
pull, read live from `ActiveRimTraction.cell_net_force`), followers dragged behind the
leaders, synchronized top-down (x-y footprint) + side (x-z) views, cells coloured by 3-zone
state, junction-switched ringed, title tracks step / N cells / A/A₀ + the measured outward
centroid migration. The original morphology MP4 (`figs/active_spheroid_morphology.mp4`,
footprint/flattening view) is retained.

## Results (original first-cut, superseded by the recalibration above)
| | PASSIVE (no traction/switch/div) | ACTIVE (old fast-division) |
|---|---|---|
| A/A₀ | 0.995 → **1.58** | 0.998 → **2.86** (with stacked divisions) |
| switched cells | 0 | **8** (bulk pressure > 0.5 kPa) |
| divisions | 0 | **4** (12→16 cells) |

**GPU wiring (gbook-ready, commit `83df8ce`):** `cell/dcm_gpu_forces.py` device-dispatch
(`DeviceDispatch`: GPU → gpu_local_snapshot + cupy + kernels_gpu; CPU → cpu_local + numpy +
kernels_cpu, cupy import guarded), `DcmTentContactGPU` + `DcmSubstrateForceGPU`, a
parity-matched `tent_contact_forces` kernel (kernels_cpu/gpu), the gbook large-run script
`scripts/dcm_gpu_run.py` (GPU auto / CPU fallback, --n-cells up to a few hundred), and
`tests/test_dcm_gpu_forces_parity.py`. CPU-path bit-parity PASS (max abs force diff ~5e-26 N
≪ 1e-12 gate); full suite 1419 passed. GPU path structurally complete, gbook-A5000-validated
later (this Mac has no CUDA). Run on gbook:
`python -m ffn_sim.scripts.dcm_gpu_run --n-cells 200 --equil-blocks 18`. Traction ALONE (no division) gives A/A₀ ≈ 3.0–3.7 (in band, ≈ the
single-cell engine's 2.03 and the passive native+tent 2.36). Figure:
`figs/active_spheroid.png` (6 panels: A/A₀ active-vs-passive, cell number, junction-switch +
pressure, 3-zone fractions, final cluster with switched cells, per-cell pressure+integrin).

## Honest caveats
- **Migration is now real but modest at this scale.** The whole-cell migration term DOES
  translocate rim-cell centroids (+0.41–0.64 µm vs the old 0.011 µm — a genuine, measured
  outward crawl, not just basal splay), and the trails/arrows in the migration MP4 make it
  visible. But at single-cell R=7.5 µm with only ~10-14 cells, the compact ball detects only
  1-2 cells as "rim" (the convex-hull/low-neighbour periphery is small), so the COLLECTIVE
  drag of many followers behind many leaders is suggestive rather than dramatic — the effect
  scales with cell number and is a many-cell/GPU showcase. The traction still also bulges the
  basal footprint (turgor + the basal lamellipodial lip), so A/A₀ is migration + bulge
  combined; the centroid-displacement number isolates the migration component honestly.
- The old fast-division full run reached A/A₀≈9 (traction ~2 PLUS stacked rim divisions,
  amplified by the small compact-ball A₀ at CPU single-cell scale) — smooth/finite, NOT the
  old A/A₀=17 over-dispersion bug. With the new slow-division default A/A₀ stays in the 2-4
  band (≈2.02 in `--quick`), traction-dominated.
- **One real bug found+fixed:** building n_max cells made a tall FCC ball whose active subset
  floated ~12 µm above the substrate (out of adhesion reach) → no basal contact → A/A₀→0.18.
  Fixed by re-packing the active cells into a compact ball lowered onto the floor.
- **Scale:** at single-cell R=7.5 µm the cluster never spans the ~150 µm necrotic depth → no
  necrotic core (3-zone shows proliferating/quiescent only). A real necrotic-core spheroid
  (~300-500 µm) + bead-resolved per-cell lamellipodium = GPU/many-cell future
  (`docs/v2_audit/GPU_SCALE_PATH_2026-06-11.md`). `integrator/baoab.py` untouched.

## Significance
The active mechanisms the static sweep lacked are now confirmed working: active traction
drives A/A₀ above passive wetting into the band, the bulk-pressure junction switch fires, and
division grows the footprint over time. Next: scale these on the GPU past R=150 µm (necrosis
on) → the active+necrosis static size sweep → the a+b/R+c/R² law (the path that may avoid the
6-7 wk live-division mesh-split port).

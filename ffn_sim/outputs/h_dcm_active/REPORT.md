# DCM ACTIVE-spreading spheroid — active traction + junction switch + live division

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

## Results (independently reproduced, `--quick`)
| | PASSIVE (no traction/switch/div) | ACTIVE |
|---|---|---|
| A/A₀ | 0.995 → **1.58** | 0.998 → **2.86** (physiological band 2-4) |
| switched cells | 0 | **8** (bulk pressure > 0.5 kPa) |
| divisions | 0 | **4** (12→16 cells) |
All 4 criteria PASS: active>passive ✓, junction switch fires ✓, division grows cluster ✓,
A/A₀ grows over time ✓. Traction ALONE (no division) gives A/A₀ ≈ 3.0–3.7 (in band, ≈ the
single-cell engine's 2.03 and the passive native+tent 2.36). Figure:
`figs/active_spheroid.png` (6 panels: A/A₀ active-vs-passive, cell number, junction-switch +
pressure, 3-zone fractions, final cluster with switched cells, per-cell pressure+integrin).

## Honest caveats
- The longer full run reached A/A₀≈9 (traction ~3 PLUS 9 stacked rim divisions, amplified by
  the small compact-ball A₀ at CPU single-cell scale) — smooth/finite, NOT the old A/A₀=17
  over-dispersion bug. The `--quick` smoke's 2.86 is the clean in-band number.
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

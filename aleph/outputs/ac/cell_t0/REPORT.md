# AC engine — cell at physiological t0 (first 3-D visual)

The first interactive 3-D render of the new Active Cell engine: the whole MCF7 cell at its physiological
resting (SUSPENDED) t0 baseline, composed from the real structures each track/engine builds. STATIC geometry
(no Warp/CUDA run — the t0 shell needs no solve); the dynamics + physics-field viz (pore-pressure NG-10,
power-stroke, adherent flatten) needs the native run on the gbook A5000.

Regenerate: `PYTHONPATH=. python ffn_sim/scripts/ac_cell_t0_viz.py` (pure numpy, dev-Mac; reuses the
`ff_viewer_html` three.js viewer — OrbitControls + scene dropdown + live cut plane).

## What is rendered (FULL native population, no downsampling)

| Structure | Source | Count |
|---|---|---|
| Plasma membrane (Helfrich icosphere) | `ff/membrane_surface.build_membrane_mesh` | 2,562 verts / 5,120 faces |
| Cortex F-actin (crosslinked mesh) | `ff/gamma_floor.build_crosslinked_cortex` | **70,686 filaments · 494,802 nodes · 424,116 segments** |
| Crosslinkers | (same) | 70,686 links |
| Nuclear lamina (deformable oblate mesh) | `ac/nucleus/geometry.build_oblate_mesh` | 642 verts / 1,280 faces |
| NMII minifilaments (head-resolved) | `ac/motor/minifilament_warp.build_minifilament_nodes` | **442 minifilaments · 6,188 backbone beads · 8,840 explicit heads** |

## Figures

- **`ac_cell_t0_whole.png`** — whole cell: the full 70,686-filament woven cortical actin shell + 442 NMII
  motors scattered on it. The full native population, full-res.
- **`ac_cell_t0_cutaway.png`** — cut-away: the cortex back-hemisphere cup opened to reveal the nuclear lamina
  envelope inside + explicit NMII motor heads on the shell/nucleus surface.
- **`ac_cell_t0_motors.png`** — motor emphasis: crosslinkers + the head-resolved NMII minifilaments.
- Interactive source: `ac_cell_t0.html` (36 MB, WebGL, rotate/zoom/pan + live cut plane; regenerable, not
  git-tracked — push to the Pages gallery/release if sharing).

Browser-verified via `ffn_sim/scripts/browser_check.py` (headless Chrome): 3 scenes rendered clean, no JS
errors, blank-check passed, PNGs eyeballed.

## Honest caveats (integrity)

1. **STATIC t0 geometry, not dynamics.** The physics-field/dynamics viz (pore-pressure, myosin power stroke,
   adherent flatten) requires the native Warp run on gbook — out of scope for this static render.
2. **Nucleus aspect = 1.0 (spherical).** MCF7 oblate aspect is a GAP with the ratified "EMERGE, do not set"
   policy (I7 cap load produces the flatten), so the honest suspended-baseline sphere is rendered and flatten
   is labeled an emergent output — no magic number.
3. **Nucleus off-centre offset (~1.4 µm) is illustrative only** (no sourced MCF7 eccentricity; render
   placement, not a kernel input). R_nuc = 0.70·R = 5.25 µm (FF resting convention).
4. **NMII shown at real scale** (backbone 0.301 µm) — small vs the 7.5 µm shell; heads drawn as head↔backbone
   bond stubs (the viewer forces the point layer to world-unit size → giant blobs otherwise), preserving the
   head-tethering topology without over-sizing. The viewer was not modified.

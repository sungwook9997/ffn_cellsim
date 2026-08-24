# Staged ac/engine build — component/connector visualiser (ac-viz)

Makes the staged `ac/engine` component/connector build **visually readable**: one isolation scene per
component, one coloured load-path scene per connector FAMILY (the explicit joints — co-location ≠ connection),
a cumulative composition view that grows as slices bind, and the composed whole cell coloured by per-node |F|.
Fully **data-driven** from a state dump manifest — a new component or connector family renders the moment it
appears in the dump, with no viz edit.

## Pieces (all viz-owned; NO physics touched)

| file | role |
|---|---|
| `scripts/dump_state.py` | **writer** — `write_v2_dump` (single schema source) + producers: `legacy_npz_to_v2` (Mac), `dump_incumbent_cell` (gbook CUDA), `dump_engine_actor` (gbook, ac/engine graph) |
| `scripts/ac_viz_common.py` | **reader + scene model** — normalises v2 & legacy dumps; builds isolation / connector / cumulative / composed / steric scenes |
| `scripts/ac_cell_assembled_viz.py` | render a dump → ONE self-contained interactive WebGL HTML |
| `scripts/ac_cell_dynamics_viz.py` | flowing-physics (FSI field) viewer (unchanged data path) |
| `scripts/ac_synth_dump.py` | synthetic full-graph dump (REAL `reference_cell_architecture()` topology, synthetic geometry) — render/schema proof only, NOT physics |
| `scripts/ac_stage_render.py` | **P0/P1 stage-render emit point** — render + browser-check curated scenes + flagship GIFs + gallery card |
| `scripts/browser_check.py` | headless render proof (blank / JS-error detection) |

## Wiring: gbook (CUDA dump) → dev Mac (render + browser-check → gallery)

Dumps require Warp-CUDA (gbook A5000); rendering is numpy-only (dev Mac). The physics stage EMITs a dump then
calls the render entry point — it never imports the viewer, and the viewer never launches a kernel.

```bash
# ── gbook A5000 (CUDA): emit a v2 dump when a slice reaches CUDA_UNIT / GO ──
#    PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim   (avoid the stale editable shadow)
# P0 — resting cell + residual heatmap (incumbent build_cell, full native N=70,686):
python -m ffn_sim.scripts.dump_state --incumbent --out outputs/ac/cell_assembled/p0_resting_v2.npz
# P1 — SF slice + FA/cortex connector (ac/engine bound actor graph):
python -m ffn_sim.scripts.dump_state --engine    --out outputs/ac/cell_assembled/p1_sf_v2.npz

# ── dev Mac (no CUDA): render + browser-check + gallery card, per stage ──
python -m ffn_sim.scripts.ac_stage_render --stage P0 --npz outputs/ac/cell_assembled/p0_resting_v2.npz
python -m ffn_sim.scripts.ac_stage_render --stage P1 --npz outputs/ac/cell_assembled/p1_sf_v2.npz
```

`dump_engine_actor(actor, arch, out)` is the forward path once a `CellActor` is bound on gbook — pass the bound
actor + `reference_cell_architecture()`; it reads each bound component's `position_d`/`force_d`/segments/
`persistent_filament_id_d` and each connector's resolved endpoints (best-effort over whatever is bound, so it
grows as slices land). Topology enumeration is CPU-safe; per-node reads need CUDA.

## Gallery

`ac_stage_render` writes `outputs/ac/stage_gallery/<stage>/` (HTML + browser-check PNGs + turntable GIFs +
`card.json`) and a data-driven `outputs/ac/stage_gallery/index.html`. Publish via the established flagship
gh-pages workflow. The full-resolution interactive HTML (no downsampling → hundreds of MB) goes to a GitHub
**release**, not the gallery; the gallery carries the small turntable GIFs + a link to the release HTML.

## Invariants honoured

- **FULL native resolution, NO downsampling** (project rule) — the interactive HTML carries every node.
- **Global-unique actor IDs** per component (`ACTOR_ID_STRIDE` namespacing) → each physical filament belongs to
  exactly ONE component; `write_v2_dump` asserts cross-component disjointness (no double-draw).
- **Colour integrity** — turbo ramp, |F| on log₁₀ (>4 decades), smooth loads linear, units on every colour bar,
  no axis truncation; ranges are display bounds (printed), never tuned to an outcome.
- **Render is proven, not asserted** — every stage scene is browser-checked (blank / JS-error detection); the
  visual screenshot is the acceptance, not a grep.

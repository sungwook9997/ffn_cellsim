# H.7 active-γ architecture investigation — step 1: load-path vs connectivity

**Date** 2026-06-09 · **Branch** `h7/full-cell-integration` · **Tool**
`scripts/h7_loadpath_architecture_sweep.py` (STATIC seed-only mesh topology; no sim,
no binder, no production-config change — mechanism scoping, band LOCKED).

## Why this step

active-γ saga closed: the floor is generation/TRANSMISSION-bound. The transmission wall
(WALL A in `h7_active_force_budget.py`) is that the actin-network γ channel carries only
~28 % of the myosin-dipole γ (amplification ~1.3×), because the contraction load-path
between anchors is ~1 backbone segment (investigation #2: 73 % single-segment, mean
1.52 seg). Chugh 2017 / Truong Quang 2021: real cortical tension is set by
actin-filament-length / actin-myosin OVERLAP — a LONG load-path — largely independent of
myosin count.

**PI fork (2026-06-09): "아키텍처 갈아엎어야 되는 거 아닌가?"** Is the cortex
CONSTRUCTION the lever? This step answers the prerequisite WITHOUT a contraction run:
across the construction knobs that set the load-path (`z_struct` anchor density,
`L_long_mean` formin backbone length, `bundle_mult`), can the mesh produce a LONG
load-path (mean inter-anchor span ≫ 1 seg) while keeping connectivity (giant ≥ 0.9)?

## Method

`build_connected_cortex(..., with_simulation=False)` returns the seeded mesh in
milliseconds. Metrics per knob set: giant-component fraction (filament-graph),
realised z, L/lc, and the inter-anchor free-span distribution along each filament
contour (anchors = crosslink-attach beads ∪ Arp2/3 branch-junction beads = the
transmission load-path). **Scale matters:** at n=160/300 the mesh is below percolation
(giant 12–52 %); the production ×40 scale is **n=1000** (reach √(A/n)=840 nm) where
giant ≈ 99 % is recovered (matches the 2026-06-04 cortex-rebuild record). All numbers
below are n=1000.

## Result (n=1000, full ×40 scale)

Connectivity is robust (giant ≈ 99 %) for z ≳ 2.2; the percolation knee is ~z=1.4–1.8.
The load-path is **geometrically pinned ~1.3–1.5 seg** at the production knobs and moves
only weakly:

| knob set | giant | z_real | L/lc | mean span (seg) | 1-seg % | path (µm) |
|---|---|---|---|---|---|---|
| **baseline** z=2.8 L=5 bundle=2 | 98.9 % | 3.33 | 5.1 | **1.31** | 87 % | 0.65 |
| z=3.7 L=5 bundle=2 (prod default) | 99.3 % | 3.97 | 6.0 | 1.27 | 88 % | 0.64 |
| z=2.8 L=10 bundle=2 | 98.7 % | 3.41 | 5.1 | 1.54 | 83 % | 0.77 |
| z=2.8 L=5 **bundle=1** | 98.9 % | 3.33 | 2.5 | 1.57 | 74 % | 0.79 |
| z=2.8 L=10 **bundle=1** | 98.7 % | 3.41 | 2.6 | 1.89 | 69 % | 0.95 |
| **z=2.2 L=10 bundle=1** | 97.4 % | 2.96 | 2.0 | 2.00 | 66 % | 1.00 |
| **z=1.8 L=10 bundle=1** | 94.6 % | 2.67 | 1.7 | **2.12** | 63 % | 1.06 |
| z=1.0 L=10 bundle=1 | 78.2 % | 2.03 | 1.0 | 2.51 | 56 % | 1.25 | (fragmented) |

## Findings

1. **`bundle_mult` is the dominant load-path knob, not filament length.** Bundling
   clusters `bundle_mult` crosslinks at each bridge site → short spans. De-bundling
   (2→1) lengthens the mean span at every (z, L) (e.g. z=2.8 L=10: 1.54→1.89). But
   bundling exists for cortex stiffness (Flormann 2024) and to hit the **L/lc ≥ 5.9
   acceptance gate** — at bundle=1, L/lc collapses to ~2 (gate violated). So the load-path
   knob and the stiffness/L/lc gate are in direct opposition.

2. **Filament length helps weakly.** L_long 5→10 µm adds ~+0.3 seg (1.31→1.54 at
   bundle=2). Longer filaments just acquire proportionally more anchors (L/lc ≈ const),
   keeping anchor SPACING ≈ const — the seeding distributes z anchors per filament, so
   spacing is set by z, not by length.

3. **A connected long-load-path region exists but is shallow.** z=1.8–2.2, L=10 µm,
   bundle=1 gives giant 95–97 % with mean span ~2.0–2.1 seg (~1.0 µm) — roughly **2×**
   the baseline. Pushing further (z=1.0) reaches 2.5 seg but fragments (giant 78 %).

4. **The achievable lengthening (~2×) is far short of the gap.** WALL-A amplification is
   ~1.3× at baseline; even if it scaled linearly with load-path, 2× load-path → ~2.6×,
   against a **30–80×** generation/transmission gap to the MCF7 active target (0.40 mN/m,
   Hosseini/FF 2021). Chugh/Truong-Quang overlap is filament-length-scale (several µm =
   several segments); the point-crosslink+bundle construction tops out ~2.5 seg and
   cannot reach it WITHOUT breaking its own connectivity/stiffness gates.

## Reading → PI

- **PI's instinct is supported:** within the existing point-crosslink+bundle
  construction, load-path and connectivity/stiffness are **coupled** — you cannot get a
  Chugh-scale (several-µm) load-path without dropping below the percolation knee
  (giant < 0.9) and/or violating the L/lc ≥ 5.9 + bundling-stiffness gates. A knob tweak
  buys only ~2× and would loosen two acceptance gates to do it (NOT done — gates LOCKED).
- The genuine architecture lever is therefore a **connectivity-primitive change**: carry
  connectivity on EXTENDED PARALLEL BUNDLES / actin-myosin overlap (myosin walks along a
  long overlap; crosslinks tie bundle ENDS) so the load-path is the bundle length, not the
  inter-crosslink spacing. This is the "갈아엎기" — a construction rebuild, not a knob.
- **Open quantitative item (next step):** confirm the transmission *sensitivity* — does
  WALL-A (g_actin/g_myo) actually rise with load-path? A loading-phase force-budget at
  baseline vs the z=1.8/L=10/bundle=1 long-load-path knobs (cheap, CPU) bounds the payoff
  of load-path BEFORE committing to the bundle/overlap rebuild.

## Artifacts

- `outputs/h7/production/h7_loadpath_architecture_sweep.json` (z=2.8–4.5 connected band)
- `outputs/h7/production/h7_loadpath_lowz_sweep.json` (z=1.0–2.8 frontier, bundle 1/2)
- `outputs/h7/figs/h7_loadpath_lowz_sweep.png`
- No production config changed; no gate touched.

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

## Step 2 — uniform vs faithful mesh + WALL-A loading sensitivity

Reconciliation discovered while wiring the WALL-A test: **the production active-γ /
Gate-A / Gate-B / force-budget all run the UNIFORM cortex** (`faithful_connected_mesh`
unset → uniform 3 µm filaments, beads_per_filament=7, NO Arp2/3 branches). The
bimodal "faithful" mesh (Arp2/3 dendritic branches + bimodal lengths) exists in code
(`build_baseline_cell(faithful_connected_mesh=True)`) but is NOT used for the γ
measurement. Step 1's sweep used `build_connected_cortex` = the FAITHFUL path.

WALL-A loading runs (n=1000, h7_active_force_budget, s_grip≈0 loading):

| build | knobs | giant | engaged | g_myo (dipole) | g_actin (network) | g_soft |
|---|---|---|---|---|---|---|
| uniform (prod) | z=3.7 b=2 | conn | 5.0 % | 1.06e-4 | 1.68e-3 | 1.57e-3 |
| uniform | z=1.8 b=1 | **2.3 %** ✗ | 5.1 % | 1.12e-4 | 1.65e-3 | 1.54e-3 |
| faithful | z=3.7 b=2 | conn | 3.0 % | 1.13e-4 | 3.38e-3 | 3.33e-3 |
| faithful | z=1.8 b=1 | conn | 3.3 % | 1.13e-4 | 3.40e-3 | 3.34e-3 |

Two findings:

5. **Arp2/3 branching decouples connectivity from crosslink density.** The UNIFORM
   mesh FRAGMENTS at low z/bundle (giant 2.3 % at z=1.8/b=1) — it cannot lengthen its
   load-path at all without falling apart. The FAITHFUL mesh stays connected at the
   same sparse knobs because the dendritic branches carry connectivity, freeing
   crosslinks to be sparse → the longer load-path of step 1 is ONLY available on the
   branched mesh. So the architecture substrate for any load-path change is the
   faithful (branched) mesh, not the production uniform one.

6. **The LOADING phase is load-path-INSENSITIVE — the WALL-A test needs CONTRACTION.**
   Across faithful baseline vs faithful long-path, every myosin signal (g_myo, g_IK
   1.46e-5→1.55e-5, gen_force 0.135→0.146 nN) and g_soft (3.33e-3→3.34e-3) is
   unchanged. The loading g_soft is dominated by PASSIVE network prestress
   (seeded-adhered mesh + turgor on the backbone bonds), not myosin — so it cannot
   isolate myosin transmission. Whether a longer load-path raises the myosin-DRIVEN
   network tension (WALL-A) can only be seen once s_grip develops (contraction), which
   is the binder-host-sync-bound GPU run (~30 min/config on gbook).

## Net (steps 1+2) → the PI fork

The architecture lever is REAL and partly already in code: switch production from the
uniform mesh to the faithful branched mesh, then exploit branch connectivity to seed
SPARSE crosslinks (low z, no bundling) for a longer load-path. BUT: (i) the achievable
load-path is ~2× (≤2.5 seg ≈1.25 µm), far short of Chugh/Truong-Quang filament-scale
overlap and of the 30–80× generation/transmission gap; (ii) whether even that 2× raises
myosin transmission can only be measured by a GPU contraction run (binder cost). The
cheap, static + loading levers are exhausted; the next step is a genuine commitment
(big build and/or GPU contraction), which is the PI decision.

## Artifacts

- `outputs/h7/production/h7_loadpath_architecture_sweep.json` (z=2.8–4.5 connected band)
- `outputs/h7/production/h7_loadpath_lowz_sweep.json` (z=1.0–2.8 frontier, bundle 1/2)
- `outputs/h7/figs/h7_loadpath_lowz_sweep.png`
- No production config changed; no gate touched.

## Step 3 — (B) execution hits a gate-contract wall (PI-decided: faithful rebuild)

PI chose (B): commence the faithful-mesh rebuild. Executing it WITHIN the hard rules
(physiological params, no gate-loosening) reveals a wall. The cortex-rebuild acceptance
gates are z ∈ [3.0,3.5] (Kadzik-Munro coordination), giant ≥ 0.9, **L/lc ≥ 5.9**
(Flormann bundling). The faithful-mesh configs that pass ALL THREE at the physiological
L_long = 5 µm:

| z_struct | bundle | z_real | giant | L/lc | mean span (seg) |
|---|---|---|---|---|---|
| 2.4 | 3 | 3.01 | 98.3 % | 5.9 | **1.26** |
| 2.6 | 3 | 3.16 | 97.8 % | 6.1 | 1.26 |
| **2.8** | **3** | **3.30** | **98.1 %** | **6.3** | **1.26** |

**Within physiological params + all gates, the faithful load-path is 1.26 seg — NO longer
(actually shorter) than the uniform production mesh (~1.5).** bundle=3 (required to hit
L/lc ≥ 5.9) clusters anchors → shorter spans, cancelling the branch-connectivity benefit.
Bumping L_long to 20 µm only reaches ~1.64 seg and is itself a non-physiological
parameter chase (the rebuild's bimodal L_long = 5 µm is the Fritzsche-anchored value).

The long load-path (2.1 seg, step 1) needs z_real ≈ 2.67 (< the Kadzik-Munro 3.0 floor)
AND L/lc ≈ 1.7 (< the Flormann 5.9 floor) — i.e. a **gate-contract change on TWO
physiological gates**. Per the hard rules (no gate-loosening) this is a PI sign-off item,
not an inline edit.

### What this means for (B)

1. **The transmission load-path limit is a CONSEQUENCE of physiological cortical
   connectivity, not a construction artifact.** Real cortex IS densely crosslinked
   (z ≈ 3–4, L/lc ≈ 6) → inter-anchor spans are inherently ~1 segment. A faithful rebuild,
   done faithfully, REPRODUCES the short load-path rather than removing it. This is itself
   a meaningful result: the short load-path is physiological.
2. Therefore the Chugh/Truong-Quang "overlap" tension cannot come from longer inter-crosslink
   spans — it must come from a DIFFERENT mechanism: myosin transmitting force ALONG the
   actin-myosin overlap within a bundle (walking antiparallel actin), independent of the
   crosslink spacing. The current bead-spring cortex + local-attach myosin does not
   represent bundle-contour myosin transmission — that is a MYOSIN-MODEL change, deeper than
   a mesh rebuild.
3. The faithful branched mesh remains the more physically-correct cortex architecture
   (real cortex IS Arp2/3-dendritic) and is worth adopting for fidelity — but its active-γ
   payoff vs the uniform mesh can only be settled by a CONTRACTION run (branched/bundled
   tensegrity transmission is not visible in the linear span or the loading phase).

### Open PI items (gate-contract / scope — surfaced, not actioned)

- (i) Adopt the faithful branched mesh as production for fidelity + measure active-γ on it
  via a GPU contraction run (regardless of the linear-span result) — within all gates.
- (ii) Authorize a gate-contract change (lower the Kadzik-Munro z and/or Flormann L/lc
  floor) to actually lengthen the load-path — but step 3 shows this contradicts physiological
  cortical connectivity, so it would be un-physiological.
- (iii) Reframe the lever as a MYOSIN-OVERLAP mechanism (bundle-contour transmission), which
  is the physically-honest route to Chugh/Truong-Quang overlap tension — a deeper model
  change than the mesh.

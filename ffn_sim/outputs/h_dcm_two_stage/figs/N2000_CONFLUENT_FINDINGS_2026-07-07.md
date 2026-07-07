# N=2000 confluent full-compartment DCM — run + G2 interpenetration finding (2026-07-07)

## Run

- **N=2000** confluent Voronoi spheroid, full compartment stack (membrane/cortex γ=5e-4, cytoplasm
  turgor+K_vol, nucleus E_nuc=4700Pa R_nuc=0.7R, surface tension + differential tension),
  **GPU-native cadherin** (CAD_GPU=1, mutual-nearest catch-bonds), baoab integrator.
- 2000 warmup + 15000 steps, A5000, ~85 min wall, GPU ~332 MiB.
- Final: **bonds=53,242** (cum_broken 2,381 = active catch-bond churn at dynamic equilibrium),
  A/A0=1.000, V/V0=1.000, cfl=0.00, drift=0.00 — mechanically stable throughout.
- npz: `outputs/h_dcm_two_stage/production/fullcomp_assembly_n2000_n2000.npz` (127 MB, 26 frames,
  324,000 nodes, 640,000 faces).
- Viewer: `figs/N2000_confluent_stress_junction.html` (real |force| FEM stress heatmap +
  per-node cadherin junction count + actual bond node-pairs drawn as line segments).

## GPU cadherin port — validated by the output

- **≤1 bond/node invariant HOLDS**: max bonds/node = 1 across all 324k nodes (mutual-nearest form
  rule works as designed). 53,242 bonds, 32.9% of nodes bonded.
- Emergent de-cohesion churn: bonds form (0 → 53k over the run) and break (cum_broken 2,381) in a
  dynamic catch-bond equilibrium — not a static latch.

## G2_interpenetration = FAIL — honest decomposition

The run's gate reports `pen_frac_peak = 75.52` vs threshold 0.3 → **FAIL**. Direct recomputation on
the actual final-frame geometry (Ericson closest-point, node-inside-neighbor sign test):

| quantity | value |
|---|---|
| mean_edge (deformed final mesh) | **2.09 µm** |
| max penetration depth | **5.23 µm = 70% of R (7.5µm)** |
| median inside-depth | 1.45 µm = 19% of R |
| nodes on inner side of a neighbor wall | 81.7% |
| nodes >20%R deep (real overlap) | 37.6% |
| **honest pen_frac (max_pen / mean_edge)** | **2.50** (over-estimate; true ≤ this) |

**Two separate facts, both stated honestly:**

1. **There IS real cell-cell overlap.** Membranes interpenetrate by ~1.5 µm median (19% of R), up
   to 5.2 µm (70% of R). This is a **static property of the confluent Voronoi initialization**
   (identical pen_frac from frame 1 → 15000, no growth — NOT a dynamic blow-up). It is the
   "feasibilization equilibrium of already-inside nodes under the strong cadherin bundle" that the
   contact penalty does not fully resolve (documented in `dcm_warp_decohesion.py` ~L528-531).
   → This **corrects** an earlier verbal claim in-session that the config was "non-penetrating."
   At N=2000 confluent, it is not — there is measurable overlap.

2. **The gate number 75.52 is ~30× inflated vs the physical overlap.** The honest recompute gives
   pen_frac ≈ 2.50 using the actual deformed-mesh mean_edge (2.09 µm); the sim's 75.52 implies its
   normalizing mean_edge ≈ 0.069 µm — ~30× smaller. The sim computes `mean_edge` once from the
   INITIAL confluent edge list; on the Voronoi/lloyd init that value is far below the settled-mesh
   edge length, so the gate ratio is over-normalized. The honest pen_frac 2.50 matches the **N=400
   confluent baseline (~2.1–2.6)** previously characterized as "confluent geometry, not failure."

## Surfaced to PI (NOT auto-fixed)

Per the hard rules (failing gate + non-trivial finding → surface, visualize, do not loosen):

- **Is the residual confluent overlap acceptable?** Confluent tissue shares walls, but ~19%R median
  membrane crossing is more than physiological apposition. Options if not: (a) stronger/longer
  contact relaxation from the confluent init, (b) a less-overlapping init (larger inset / gap),
  (c) accept it as the confluent-tissue modeling state. **PI call.**
- **Gate-metric question:** the G2 `mean_edge` normalization (initial confluent edge list vs
  settled-mesh edge) inflates pen_frac ~30× at this N. A gate contract change (normalize by the
  current-frame mean_edge, or a physical-µm threshold) is a **gate-contract change → PI sign-off**,
  not an inline edit.

## Independent code audit (this session)

3-dimension adversarial audit (GPU cadherin kernels · viewer render · physiological baseline):
**12 candidate findings raised, all 12 REFUTED on verification, 0 confirmed bugs.** The cadherin
port, the large-N viewer (frame-stride consistency), and the run's physiological setpoints hold up.
The mean_edge inflation above is a run-data property, not a code defect the static audit would catch.

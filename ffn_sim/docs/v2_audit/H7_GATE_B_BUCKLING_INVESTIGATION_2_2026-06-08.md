# H.7 Gate-B design-investigation #2 — the built-mesh free-length + per-span load probe

**Date:** 2026-06-08. Branch `h7/full-cell-integration`. Author: Lead.
Follows `H7_GATE_B_BUCKLING_INVESTIGATION_2026-06-08.md` (#1, first-order Euler), which
deferred the quantitative margin to "measure, on the built physiological cell, (a) the
anchor-to-anchor free-length distribution and (b) the actual per-span myosin load." This is
that probe. Feeds `H7_GATE_B_CONTRACT_2026-06-08.md` §3 (isolate the condensation
suppressor before building the relaxed-constraint mode) + §5 control 3 (condensation must
actually be able to engage).

Probe: `ffn_sim/scripts/h7_gate_b_probe.py`. Figure:
`outputs/h7/figs/h7_gate_b_freelength_probe.png`. Data:
`outputs/h7/production/gate_b_probe_seed{1,2}.json`.

## Method

Builds the suspended/rounded MCF7 cell at the Gate-B operating point (FA OFF, turgor ON,
contract §6) at smoke scale (160 filaments), runs the unconstrained warm-up so the native
myosin binder engages heads, then measures from the settled snapshot, with NO tuning (all
constants from `phase1_h3.yaml` / `mcf7_baseline.yaml`):

- **Filaments** from the `cortex-bond` backbone graph (degree-1-tip walk; robust to
  uniform / variable-length layouts).
- **Anchors** = cortex beads carrying a crosslink attachment (`xlink_attach_b*`). The FA
  path is off (suspended), so crosslinks are the lateral supports.
- **Free spans** = contour between consecutive anchors, in segments (×ℓ₀ = 500 nm).
- **Bucklability:** M-SHAKE pins each *segment*'s axial length, but a multi-segment span can
  still buckle by hinging at the angle joints (each segment stays rigid-length, the chain
  bows; resisted by the `cortex-angle` harmonic). So a span needs **≥2 segments (≥1 free
  internal bead)** to buckle at all. A **1-segment span between two crosslinked beads is a
  single M-SHAKE-rigid rod — no internal hinge — and cannot buckle at any force.**
- **Per-span load** = engaged myosin heads (`cortex_myosin_attach_b*`) on a span's internal
  beads × `F_stall_per_head`. Margin = F_load / F_crit, F_crit = π²κ_B/L².

## Findings (seed 1 / seed 2, both n_filaments=160)

| quantity | seed 1 | seed 2 |
|---|---|---|
| inter-anchor spans | 520 | 516 |
| **1-segment RIGID (cannot buckle)** | **72.5%** | **73.8%** |
| geometrically bucklable (≥2 seg) | 27.5% | 26.2% |
| mean span | 1.52 seg (0.76 µm) | 1.54 seg (0.77 µm) |
| total engaged myosin heads | 32 | 37 |
| heads on bucklable spans | 6 | 12 |
| bucklable spans carrying ≥1 head | 5 / 143 (3.5%) | 9 / 135 (6.7%) |
| max heads on a single span | 2 | 3 |

Segment histogram (seed 1): `{1:377, 2:62, 3:48, 4:22, 5:8, 6:3}`.

### (a) The suppressor is ANCHOR DENSITY — convention-independent

**~73% of load-bearing free spans are single M-SHAKE-rigid rods with no internal hinge**, so
buckling is **geometrically blocked irrespective of load, stiffness, or any force-scaling
convention.** This is contract §3 candidate (i) — the M-SHAKE fixed-length removing the
compression DoF — but localized precisely: it is not M-SHAKE per se (a ≥2-segment span DOES
buckle by hinging at angle joints), it is the **crosslink density pinning most spans to a
single segment**, leaving no free internal bead to deflect. Investigation #1's first-order
estimate (free length ≈ crosslink spacing ≈ 1 segment) is confirmed by the built mesh.

### (b) Of the bucklable minority, binding throughput — not force — is the limiter

Only **3.5–6.7% of the bucklable spans currently carry any myosin head** (≈30 engaged heads
cell-wide; the overnight ~12/2000 under-binding). When a head *does* land on a bucklable
span it is far over threshold (below), so the bucklable spans are **binding-throughput
limited, not force-margin limited.** This couples Gate-B's transmission lever to the known
myosin under-binding issue, exactly as investigation #1 §2 anticipated.

### Margin + the mesoscale-consistency item (surface to PI)

F_crit(1-seg, 500 nm) = 2.76 pN; F_crit(2-seg, 1 µm) = 0.69 pN. The per-head force in the
sim is the **mesoscale-scaled 8.48 pN** (= 2.0 pN × 4.24 `mesoscale_force_factor`), NOT the
unscaled 2.0 pN of investigation #1's table.

The ratified ×40 convention (`cortex.py:10-13`) coarse-grains **filament COUNT, not
per-filament stiffness** — so F_crit's single-filament κ_B is *correct by convention*, and
my probe's first-pass warning ("κ_B should scale too") was wrong and is corrected. The real
tension is subtler: `mesoscale_force_scaling` raises the **local per-span load ×4.24** for
**aggregate-stress correctness** (compensating the reduced motor count), then applies it to a
single-filament-stiffness chain — whereas buckling is a **local per-filament instability.**
As-simulated, a single engaged head is ~12× over F_crit on a bucklable span; de-scaled to
per-native 2 pN it sits at the investigation-#1 margin (0.7). **Applying an areal/aggregate
force correction to a local buckling instability mixes a global correction into a local one
— resolve which force a buckling span should see before sizing the relaxed-mode
compression-release.** The dominant suppressor (anchor density, ~73% rigid) is independent of
this resolution.

## Implication for the relaxed-constraint mode

The relaxed lever is now pinned: **NOT axial compliance** (confirmed #1; tension side stays
inextensible), and **NOT per-head force** (plentiful when engaged). It is the **joint product
of (compression-side anchor density) × (binding throughput).** "Permit buckling" must
**lengthen the compression-side free span** by thinning the anchor density on the compression
side (so spans cross ≥2 segments and gain a hinge), and it only bites where myosin actually
loads — so the throughput issue gates it. Sizing the compression-release needs the
mesoscale-force resolution above.

## Caveats / scope

- Inter-anchor spans only (both ends crosslinked) — these bear network hoop tension.
  Filament **tips** (cantilever, one free end, F_crit = π²κ/4L²) and **unanchored filaments**
  (whole 3 µm free, F_crit ≈ 0.08 pN) buckle trivially but do **not** transmit spanning
  tension, so they are excluded from the load-bearing analysis by design.
- Free length taken as contour (n_seg × ℓ₀); the network is pre-stressed, so the in-situ
  column length differs slightly. First-order, as #1.
- Smoke scale (160 filaments); the distribution is a per-filament/per-crosslink local
  statistic so it transfers to the ×40 full scale (to confirm on gbook GPU at closeout).

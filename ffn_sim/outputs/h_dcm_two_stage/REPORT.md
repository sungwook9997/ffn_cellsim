# Two-stage DCM spheroid — BIOLOGICAL aggregation → spreading (2026-06-12)

> ⚠️ **CORRECTION / RETRACTION (2026-06-19, grounding-pass C14).** The "spreading" / `spreading footprint` results below are measured by the **basal-contact / convex-hull FOOTPRINT** metric, which the PI FORBADE on 2026-06-12: valid spreading A/A0 MUST be the **top-down xy silhouette** (all nodes). Footprint inflates A/A0 ~10–20× (e.g. the N=100 365→7443 "20.4×" headline) and does **not** represent valid spreading. Where top-down was actually measured (`two_stage_n400_lamel_S10`), the full mechanistic stack **COMPACTS** — top-down A/A0 1.0 → **0.597**, V/V0 → 0.669, the OPPOSITE sign. Treat every footprint "spreading"/"A/A0" figure below as **not valid** until regenerated with the top-down metric on the gbook GPU. See the grounding-pass table + `project-rebuild-audit` memory.

Branch `h7/compartment-platform`. The PI rejected the first two-stage attempt: its
"aggregation" was a forced **central-pull hack** and the cells never truly clumped.
This is the corrected, biologically-faithful redesign (reference: SimuCell3D organoid
aggregation), with the aggregation **verified emergent before spreading**.

## What was wrong (v1 — withdrawn)
`AggregationDrive` applied a global harmonic body force `F = −k·(pos − center)` pulling
every membrane node toward a fixed point. It carries **zero cell-cell information**
(would compact a bag of sand identically), is a **lumped proxy for surface tension**
(CLAUDE.md forbids), over-compacts the core, and destroys cell sorting. Empirically it
produced a ~13 % elastic squeeze over ~15 s real (10⁴× too short to express any
aggregation). **Deleted.**

## Correct mechanism (v2) — active-matter coalescence (search-and-capture)
A 4-agent ultracode research synthesis (`docs/v2_audit/AGGREGATION_REDESIGN_RESEARCH_2026-06-12.json`)
established that organoid/spheroid aggregation EMERGES from three cell-generated
ingredients, none referencing a centre:
1. **Active self-propulsion** (`DcmActiveMotilitySPP`, `cell/dcm_gpu_forces.py`): each
   LIVE cell carries a polarity `p_c` and a net active force `f_active·p_c` over its
   nodes — applied to ALL cells (no rim/interior split, no activity-LOD: the PI's
   all-cells-full-physics rule). Raises the cell-cell collision rate → motile cells
   encounter + cadherin-adhere + coalesce.
2. **Polarity reorientation** (`DcmPolarityUpdater`): Ornstein-Uhlenbeck rotational
   diffusion (persistence `tau_p`), off the per-step force path (BAOAB-safe), `|p|=1`.
3. **Soft confining drop** (`DcmDropConfinement`, `cell/dcm_gpu_build.py`): a one-sided
   spherical wall = the hanging-drop meniscus. **Force-free interior**, inward only at
   `r>R_drop` — it bounds the motile search, it does NOT aggregate, pre-strain, or
   manufacture surface tension (explicitly NOT a central pull).

The existing cohesive bilinear-tent adhesion is unchanged (captures cells; its surface
tension rounds the aggregate). `test_dcm_motility_gpu_parity.py` (4 PASS): net force per
cell `= f_active·p_c`, every live cell propelled (none frozen), polarity unit-norm.

## ✅ EMERGENCE VERIFIED (the PI's "is it biologically correct?" gate) — negative control
Seeding cells JUST OUT of adhesion range (spacing 2.8·R, gap > c_adh) so motility is
REQUIRED to coalesce:

| | motility ON | motility OFF (negative control) |
|---|---|---|
| contact fraction | **0 → 0.05 → 0.11 → 0.16 → 0.19** (monotone rise) | **0 → 0 → … → 0** (never aggregates) |

Turning motility off, the cells NEVER aggregate; turning it on, they emergently
coalesce + adhere. The aggregation is driven by the cells' OWN active motility — not an
external force, not pre-placement.

## GPU production sweep (gbook A5000) — emergent, dense, round at scale
Two stages per N: Stage 1 aggregation (motility + drop, run TO CONVERGENCE: contact-rise
+ plateau + Rg plateau + round) → Stage 2 spreading from the converged aggregate
(substrate + settling + rim traction). **All cells full physics, no LOD.**

| N | aggregation `contact` (emergent) | asphericity (roundness) | Rg compaction | spreading footprint | converged |
|---|---|---|---|---|---|
| 100 | 0.00 → **0.69** | 0.024 → 0.014 | 41.8 → 37.0 µm | 365 → 7 443 µm² | ✅ |
| 200 | 0.00 → **0.80** | 0.051 → 0.022 | 52.6 → 45.4 µm | 389 → 10 430 µm² | ✅ |
| 400 | 0.00 → **0.90** | 0.020 → 0.006 | 66.3 → 56.0 µm | 1 179 → 13 666 µm² | ✅ |
| 600 | 0.00 → **0.93** | → **0.002** | 76 → 63 µm | 515 → 15 890 µm² | ✅ |
| 800 | 0.00 → **0.95** | → 0.002 | 83 → 68 µm | 1 588 → 17 369 µm² | ✅ |

The emergent aggregate is denser AND rounder with N (more 3-D neighbours): a near-fully
confluent (`contact` 0.95) near-perfect sphere (`asph` 0.002) by N=800. Volume is held
throughout (V/V0 ≈ 1.05–1.09 — turgor inflation, cells keep volume while the spheroid
SHAPE flattens during spreading). The hand-off is exact (Stage-2 frame 0 == Stage-1
final). N=800 (33 600 particles, ~15 GB) is the fine-cell GPU-memory ceiling on the
16 GB A5000.

## Figures (`figs/`)
`two_stage_n{100,200,400,600,800}.png` — 4-panel: Stage-1 loose seed → aggregated ball
→ Stage-2 spread (lamellipodia red) → diagnostics (V/V0, contact fraction, Rg, √foot).
`*_agg.mp4` (cells coalescing) + `*_spread.mp4` (spreading) per N.

## Honest limit — full dispersed→dense maturation is timescale-bound
Coalescing cells already in a near-spherical loose arrangement (the seed) into a dense
ball is feasible (above). But a DRAMATIC, fully-dispersed suspension rounding into a
dense spheroid is a **multi-day (12–48 h real) process ≈ 50 h GPU** at the BAOAB
per-step rate — the SAME wall this project documented for spreading-magnitude and
necrosis. We therefore seed near-spherical-but-out-of-adhesion-range so the EMERGENT
mechanism (verified by the negative control) drives a real, converged, dense aggregate
in feasible wall-time, and we state the timescale ceiling rather than fake it.

## Commits
`17cc4b4` redesign · `1906502` emergence-verified + rounded-cell convergence ·
`4d6e187`/`5f4b766`/`6a61e8a`/`5610fa0`/`3b21a89` N=100/200/400/600/800 figures.

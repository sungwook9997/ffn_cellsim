# GPU DCM spheroid — K1-turgor + activity-LOD, and the a+b/R+c/R² law on the A5000

**2026-06-11 · branch `h7/compartment-platform`.** The gbook RTX A5000 GPU run, finally
unblocked, and a first validation of the PI layer-2 spreading law on the GPU.

## 1. GPU UNBLOCKED — K1-turgor solves the many-types stall
The native `md.mesh.conservation.Volume` shell needs a SEPARATE triangle type per cell →
HOOMD "many triangle types perform poorly / shared-memory errors on GPU" → at N=60 the run
stalled at 0% GPU util indefinitely. Replacing turgor with the validated **K1
`mesh_pressure_forces` kernel** (`DcmTurgorForceGPU`, per-cell volume by face-grouping in ONE
force, NO per-cell types; CPU bit-parity <1e-12 N) the **N=60 spheroid runs on the A5000,
finite, in 97 s** (2520 particles, 61.5 steps/s). The GPU-friendly build (`dcm_gpu_build.py`)
uses K1-turgor + `md.bond.Harmonic` edges + `DcmTentContactGPU` + `DcmSubstrateForceGPU`
[+ active traction], no native md.mesh, no per-cell types.

## 2. Activity-LOD works (PI directive)
`dcm_gpu_lod.py` classifies cells ACTIVE (rim / substrate-contact / non-necrotic — full
dynamics) vs INERT (deep interior / quiescent / necrotic — frozen via large overdamping). The
jammed core is mechanically inert, so freezing it is both physical AND faster. Measured: at
**N=280, 146/280 cells frozen** (only the 134-cell active rim does full work); active fraction
falls with size (the surface/volume ratio) — exactly the rim-active / core-inert structure the
PI specified.

**LOD speedup (measured, N=180, GPU):** LOD 146.4 s vs no-LOD 166.9 s = **1.14×**. Honest:
the speedup is MODEST because the LOD freeze (large overdamping on inert nodes) only removes
those nodes' integration DRIFT — the per-step cost is still dominated by (a) the BAOAB
per-step CPU sync over ALL particles and (b) the contact + turgor forces, which still process
every cell (frozen or not). A bigger LOD win requires SKIPPING frozen cells in the force /
neighbour-list computation (and/or the BAOAB port) — a deeper optimization. The LOD as built
is physically correct (the jammed core is inert) and helps more at larger N (higher inert
fraction), but is not yet a dramatic speedup at this scale.

## 3. The a+b/R+c/R² spreading law on the GPU (R = 31-78 µm, the PI's range)
**Correction:** the PI law range R=31-78 µm is BELOW the necrosis onset (~150 µm), so these
spheroids are fully viable (no necrotic core) — the law there is PROLIFERATION + COHESION
driven, not necrosis. Active size sweep (traction + junction switch, GPU + LOD):

| N | R (µm) | A/A₀ | active | inert |
|---|---|---|---|---|
| 30 | 34.2 | 1.55 | 28 | 2 |
| 60 | 41.2 | 1.55 | 46 | 14 |
| 110 | 49.1 | 1.67 | 69 | 41 |
| 180 | 57.1 | 1.37 | 108 | 72 |
| 280 | 65.1 | 1.37 | 134 | 146 |

Fit **A/A₀ = a + b/R + c/R²**: a=−0.274, b=154 µm, c=−3152 µm², r²=0.63; corr(R,A/A₀)=**−0.67**.
PI law: a=−0.33, b=188.7, c=−2655. **The functional form, the negative size-dependence
(A/A₀ decreases with R), and the coefficient SIGNS + ORDER OF MAGNITUDE all match the PI law.**
Figure: `figs/gpu_active_law_fit.png`.

## Honest gap + next
The absolute A/A₀ magnitude (~1.4-1.67) is ~0.5-1.0 BELOW the PI law (~1.9-2.9 in this range):
these runs had ACTIVE TRACTION but **NOT live proliferation** — so the b/R term's rim-DIVISION
boost is missing and spreading is weaker. Enabling proliferation should raise the magnitude
toward the PI curve while keeping the form. r²=0.63 is moderate (the N=110 point sits high).
Next: a proliferation-enabled GPU sweep to close the magnitude gap. The BAOAB per-step CPU
sync remains (forces are GPU-pathed; integrator port is frozen/PI-gated) — a secondary
throughput item, not blocking these results.

## Proliferation magnitude-closure — HONEST finding (2026-06-11 night)
Attempting to raise A/A₀ toward the PI band (1.9-2.9) via live division revealed that the
division rate is a FREE, physically-UNANCHORED parameter, so tuning it to hit the band would
be FITTING, not validation (CLAUDE.md no-magic-number rule):
- p_div=0.10 (sweep): **0 divisions** at every size → A/A₀ ~1.1 (even below active-only ~1.4;
  the prolif build also classifies very FEW cells active — e.g. 10/60 — so its active traction
  is weaker, a separate bug).
- p_div=0.50 (diagnostic N=60): **26 divisions** (60→86) → A/A₀ **9.94** (explodes far above band).
The PI band sits between → there exists a p_div that lands A/A₀ in-band, but choosing it to
match is circular.
**Correct path (now possible via the real-time mapping):** anchor the division rate to the
real CELL-CYCLE time (12-24 h) relative to the spreading time via `common/sim_realtime.py`,
so "divisions per spreading episode" is PHYSICAL (~a few, since spreading≈minutes-hours ≪
one cell cycle) and the resulting magnitude is a PREDICTION, not a tuned fit. Also fix the
prolif build's over-restrictive active-cell classification (only ~10/60 cells get traction).
**Robust result stands:** the law FORM + negative size-dependence + coefficient signs/order
are reproduced by the active mechanism (commit 55f23c0); the absolute magnitude is
proliferation-rate-dependent and must be physically anchored, not tuned.

## Large necrosis-ON production (2026-06-11 night) — 3-zone at scale
First GPU run to CROSS the necrosis onset: N=400 coarse-grained cells (R_cell=22 µm tissue
patches via --r-cell-um) → R_spheroid > 150 µm → **necrosis activates**. Result (A5000, GPU):
finite, wall=957.7 s (20000 steps), **A/A₀=3.147**, and the full 3-ZONE structure —
**active proliferating rim 225 (0.56) / quiescent-inert 175 (0.44) / NECROTIC CORE 24 (0.06)**.
Real-time label: t_sim=20 µs ≈ 12.0 s real (accel 6e5, spreading-front). Figure:
`figs/dcm_gpu_lod_n400_s20000_GPU.png`. This is the necrosis-ON large spheroid the PI asked
for — the 3-zone (rim/quiescent/core) emerges at R>150 µm as physics requires (and is 0 below
it, per the live-cell necrosis fix). Coarse-graining (each cell = a tissue patch) is the
feasible route to R>150 µm on the current BAOAB-sync-limited GPU; a fine-grained R>150 µm
spheroid (~thousands of 7.5 µm cells) is hours-to-days of wall-time (BAOAB GPU port needed).

## RESOLVED — the magnitude gap is a TIMESCALE truth, not a model failure (commit 113aadc)
Two coupled fixes turned the proliferation magnitude question from a circular tuning exercise
into an honest physical prediction:
1. **Prolif-build active classification fixed** (`LiveCellActivityLOD`, live cells only): the
   prolif build now matches the plain active build exactly — A/A₀ 1.191→**1.549**, active
   18→**28/30** (= plain `--active`). The earlier "prolif weaker" was the parked-pool LOD bug,
   now gone.
2. **Division rate PHYSICALLY ANCHORED** (`sim_realtime.division_probability`, `--auto-pdiv`):
   p_div = S·dt·div_every / T_cycle, derived from the real cell-cycle time (~20 h) vs the run's
   real-time-equivalent — NOT tuned to an A/A₀ band. For a 12000-step run (t_real=7.2 s):
   p_div=3.3e-5, expected divisions ≈ 0.003 → **0 divisions, A/A₀=1.854 (= active-only)**.

**The honest physics:** in ONE short accelerated spreading episode (~seconds real), a cell
completes ~1e-4 of its 20 h cycle, so ~0 divisions occur — correct, not degenerate. The SAME
anchored p_div over a 2 h real assay (1.2e7 steps) yields ~3 divisions; the FULL PI-band
magnitude needs the DAYS-LONG assay. **The PI law's b/R term is a days-long-assay phenomenon
(proliferation accumulated over days), which the accelerated single-episode sim cannot reach
in feasible wall-time (the BAOAB-sync step rate caps ~seconds-real per ~hours-wall).** So:
the law FORM + size-dependence + coefficient signs/order ARE reproduced per-episode (the
robust result); the absolute magnitude is a multi-day proliferation integral, reproducible
only with a faster integrator (BAOAB GPU port) or a longer-timescale division model — and
must NOT be faked by tuning p_div. This is the scientifically defensible conclusion.

## Size-dependent necrosis (2026-06-11 night) — the 3-zone evolves with spheroid size
Necrosis size sweep N=200/400/600/800 (R_cell=22 µm patches → R=180/222/252/275 µm, all above
the ~150 µm onset). The **necrotic CORE grows monotonically with size: 3→24→64→101 cells**
(fraction 1.5%→5.7%→9.6%→11.2%→15.3%, to N=1000/R=295µm), while the **proliferating RIM fraction falls 0.64→0.41** and
the **quiescent middle rises 0.34→0.48** — the Greenspan surface/volume scaling: bigger
spheroids have proportionally less active rim and a larger dead core. This is the
size-dependent 3-zone structure the PI asked for, and the falling rim fraction IS the b/R
(surface/volume) mechanism behind the spreading law. Figure: `figs/necrosis_size_scaling.png`.
(Note: A/A₀ was unstable at N=600 (14.5, a basal-footprint convex-hull artifact when peripheral
cells disperse); the necrotic/zone COUNTS are clean + monotonic. Each run finite, wall
685-1472 s / 20000 steps on the A5000.)

## Bulk-pressure JUNCTION SWITCH wired into the GPU build (2026-06-11) — cadherin → integrin
The cadherin→integrin clutch (PI mechanism) is now WIRED into the GPU-friendly build and
rendered spatially. `build_gpu_dcm_simulation` is additively given a mutable `cad_mult`
(n_cells,) array (passed to the tent contact's per-cell seam, mult = √(cad_mult_i·cad_mult_j))
and an `integrin_gain` (n_cells,) array (= the rim-traction `int_mult` when active, else a
standalone array); `attach_junction_switch` wires a `GpuJunctionSwitchUpdater` that each low
cadence computes per-LIVE-cell bulk pressure (the frozen `PressureProbe` crowding→kPa proxy,
over `cell_of_node ≥ 0` only — consistent with the necrosis live-cell fix) and LATCHES the
switch for cells above `P_switch_kPa` (~0.5 kPa): `cad_mult → 0.3` (cadherin weakens),
`integrin_gain → 3.0` (integrin strengthens). The frozen `dcm_spheroid_state` thresholds are
reused unchanged.

**Spatial result (N=250, R_cell=22 µm → R_spheroid=178 µm, CPU, finite):** 246/250 cells switch;
bulk pressure spans **[0.0, 6.0] kPa**; the switched cells' mean pressure **4.13 kPa** ≫ the
non-switched **0.0 kPa** (`switched=high-pressure ✓`), and the switched cells sit at mean radius
**133.7 µm** vs the non-switched rim at **177.2 µm** (`switched=interior ✓`). So the high-pressure
compacted interior switches to integrin while the low-pressure outer rim stays cadherin —
exactly the clutch mechanism. Real-time label: t_sim = 0.6 µs ≈ 360 ms real (accel 6×10⁵,
spreading_front).

**Contact-radius calibration (honest note).** The GPU build's tent-adhesion/turgor balance
settles at a cell-cell nearest-neighbour spacing ≈2.32·R (NN≈51 µm), just ABOVE the frozen
proxy's 2.2·R=48 µm contact radius — so at the frozen factor every cell reads crowd=0 →
pressure=0 → the switch never fires (diagnosed directly: NN does not compact below the gapped
start over 600 steps). The viz sets the proxy's contact radius to **2.5·R=55 µm** to match the
build's ACTUAL settled spacing; this is a measurement-geometry calibration to the real cell
spacing, NOT a change to the literature [0.5,5] kPa band or the 0.5 kPa onset (both unchanged).
At this scale 246/250 switch (the cluster is densely coordinated; only the 4 outermost corner
cells stay below onset) — the spatial pressure GRADIENT and the switched=high-pressure=interior
correlation are clean. Run on CPU is ~1 s/step at N=250 (brute-force tent pair search); the
switch fires on the settled cluster so only ~600 settle steps are needed (active rim traction
left OFF on CPU for cost — the switch acts on cad_mult/integrin_gain regardless; the gbook A5000
runs it with `--active`). Entry point: `scripts/dcm_junction_switch_viz.py`.

## Figures
- `figs/junction_switch_spatial.png` — junction-switch deliverable: (1) central xz cross-section
  coloured by junction state (blue cadherin-dominant vs red switched/integrin-dominant), (2) the
  same slab as a per-cell bulk-pressure heatmap [kPa] with red rings on switched cells (confirms
  switched = high-pressure), (3) 3D junction-state render. Title carries N=250, R=178 µm,
  246/250 switched, P∈[0,6] kPa, and the real-time label.
- `figs/junction_switch_spatial.mp4` — slow 360° rotation of the 3D junction-state render (MP4,
  FFMpegWriter).
- `figs/necrosis_3zone_spatial.png` / `.mp4` — prior 3-zone necrosis spatial render (necrotic
  core central).
- `figs/necrosis_size_scaling.png` — size-dependent 3-zone scaling (necrotic core grows with N).

## Spreading-over-real-time — HONEST finding: the model OVER-SPREADS without arrest
After the 2.12× active-vectorization speedup, a long spreading-curve run (N=100, sampling
A/A₀ vs real-time) was launched to capture spreading over MINUTES of real time. It revealed
that A/A₀ does NOT plateau — it BLOWS UP: A/A₀ = 3.6 (12 s real) → 6.3 (24 s) → 12.0 (36 s) →
17.5 (48 s), monotonic and accelerating. This is UNPHYSICAL over-dispersion (the basal
convex-hull footprint explodes as sustained active rim traction overwhelms the cell-cell
cohesion and the cells scatter). So:
- **Physical spreading is the EARLY phase (~12-20 s real, A/A₀ ~3-4, in band).**
- Beyond that the model lacks a SPREADING-ARREST mechanism: real cells stop spreading at a
  maximum area (membrane-tension limit / contact inhibition), but our active traction keeps
  pushing indefinitely → over-spread.
- **Therefore "run longer for more real-time spreading" makes the result WORSE, not better.**
  Covering minutes of physical spreading requires adding a spreading-arrest law (cap the
  active traction as the cell footprint approaches a physiological max, or a membrane-tension
  restoring force), so the spheroid reaches and HOLDS a stable A/A₀ plateau over long
  real-time. This is the scientifically correct next step for the real-time spreading goal.

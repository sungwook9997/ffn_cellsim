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

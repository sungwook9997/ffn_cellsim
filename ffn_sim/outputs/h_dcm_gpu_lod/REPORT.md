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

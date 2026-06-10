# GPU-scale path to the spheroid spreading law — strategy (PI 2026-06-11)

**Question (PI):** the GPU port is 827× faster — can't we just run the large spheroid now,
"same as CPU but faster"? **Answer: partly — but GPU scale changes the PHYSICS available,
not only the speed, and the real blocker was never the per-step math.**

## Three things gate the large-scale run (speed is only one)
| gate | status | note |
|---|---|---|
| GPU kernel math speed | ✅ validated 827× (group_pair), parity OK | `gpu_opt/GPU_RESULT_A5000.md` |
| native md.mesh shell (turgor/volume/bending) | ✅ already GPU-native | reframe Phase 2 (`dcm_native_shell.py`) |
| tent contact + active traction → cupy/gpu_local wiring | ⬜ **days** (not weeks), fixed-topology | the validated `kernels_gpu.py` K4 maps to `DcmTentContact` |
| BAOAB per-step CPU sync | ⬜ fixed-topology: handleable (PI sign-off to port, or a gpu-native stepper) | `integrator/baoab.py:283` cpu_local_snapshot every step (frozen) |
| gbook A5000 hardware | ⬜ dirty branch — PI cleanup pending | |
| **live-division mesh-split** ↔ HOOMD fixed-tag-space | ⬜ ~6-7 wk — **but maybe AVOIDABLE (see below)** | the actual 6-7wk item; orthogonal to the 827× speedup |

So "the kernels are fast" ≠ "the live sim is fast": the kernels are standalone math, not yet
wired into HOOMD's Custom-force + custom-BAOAB loop. Native mesh shell IS GPU-native; the
tent contact + active traction + the BAOAB per-step sync are the remaining wiring (days for
fixed topology).

## The key reframe — GPU scale unlocks NEW PHYSICS, not just speed
The CPU static size-sweep FAILED to reproduce A/A₀=a+b/R+c/R² (commit `6122fd2`: non-monotonic,
r²=0.15, signs opposite) — partly because it was **physically too small**: R ≤ 45 µm, and
**necrosis only triggers at depth > ~150 µm** (R > ~150 µm). At CPU mesoscale the necrotic
core — the size-dependent mechanism the law encodes — **cannot turn on at all**. GPU scale
(R > 150 µm, hundreds of cells) crosses that threshold. This is a physics gate, not a speed
gate.

## The potentially-fast path to the law (avoids the 6-7wk division port)
The layer-2 law A/A₀ = a + b/R + c/R² may be reproducible by a **STATIC size sweep at GPU
scale WITH active traction + necrosis** — no live division required:
- larger spheroid → necrotic core grows (∝ volume) → the ACTIVE, spreading rim fraction
  shrinks (∝ surface/volume) → A/A₀ **decreases with R** = the law's signature (b/R surface
  term, c/R² core/cohesion term).
- This is a *static* geometric+mechanical size-dependence (active rim vs dead core), distinct
  from the *dynamic* proliferation mechanism. If it reproduces the coefficients, the 6-7 wk
  live-division mesh-split port is NOT on the critical path for the law validation.
- Prerequisite: the active traction + junction-switch + necrosis mechanisms must work (being
  integrated now on CPU), THEN scale them up on GPU past R=150 µm.

## Parallel work streams (now)
1. **Active integration (running)** — active rim traction (lamellipodium/contraction-belt/
   FA-clutch, DCM-scale) + bulk-pressure junction switch + live division, on native+tent, CPU.
   Confirms the mechanisms work before scaling.
2. **GPU wiring/prep** — dispatch the validated forces (native mesh = native; tent contact →
   cupy K4; substrate → cupy K2) through cpu_local/gpu_local by device; CPU-path bit-parity
   testable on the Mac; cupy path + a gbook run-script for the A5000 large run. Preps the
   "wire it in = days" gate.
3. (future, PI-gated) BAOAB gpu_local port OR a gpu-native fixed-topology stepper; gbook
   branch cleanup; then the GPU-scale active+necrosis static sweep → the a+b/R+c/R² fit.

## Honest bottom line
- A faster BIG *fixed-N* spheroid is reachable in days (wiring), and it crosses the necrosis
  threshold the CPU can't — so it can test the law via the active+necrosis static sweep.
- The 827× does NOT by itself give live division (the mesh-split/tag-space problem is
  separate) — but live division may not be needed for the law if the static active+necrosis
  sweep reproduces it.

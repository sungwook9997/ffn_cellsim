# GPU-main Phase 2 P2c — gate re-validation sign-off PROCEDURE

**Date:** 2026-06-11  **Host:** gbook (RTX A5000 16GB, hoomd 7.0.1 gpu_enabled=True, cupy 14.1.0)
**Status:** PROPOSED — awaiting PI authorization to execute (runs are PI-gated per GPU_MAIN_PORT_PHASE2 §5/§P2c).

## Purpose

Produce the evidence package for PI to ratify the **GPU device** as a validated production
device for the H.2/H.3 persistence-length contract gates — i.e. show GPU reproduces the
CPU-validated physics. This is **not** a gate-contract change (no threshold edits); it re-runs
the existing gates on GPU and compares to the CPU baseline. Ratification stays with PI
(no-gate-loosening; "driver, not gate").

## Scope decision (what is / isn't re-validatable on GPU now)

| Gate | GPU-re-validatable now? | Why |
|---|---|---|
| Integrator parity (constrained BAOAB: free-chain D_com, bending equipartition) | ✅ yes, ~minutes | `gpu_smoke_constrained.py` runs both CPU+GPU |
| H.3 cortex per-filament L_p (VG-H3-per-filament-L_p, KU-1.1 band) | ✅ yes, ~hours | `h3_lp_gpu_production.py` already `--device gpu`-ready |
| H.2 single-filament L_p (VG-H2-L_p-C1/tail, equipartition, χ², KL) | ✅ yes, ~hours, **1 small script change** | `h2_single_filament.py` is CPU-hardcoded → add `--device` (same `add_production_device_args` pattern as h3_lp_gpu_production) |
| H.3 KU-3.20 nematic order | n/a (device-independent) | topology-only, no HOOMD sim; already PASS on CPU |
| **H.3 dynamics KU-3.1 / KU-3.18 (60 s rounding/blebbistatin)** | ❌ **OUT OF SCOPE** | (a) pytest **skeleton-only** (not implemented as runnable); (b) 60 s = ~4.6e9 BAOAB steps is **brute-force infeasible** even on GPU (dt≈13 ns bottleneck; ~5357× over budget per `reference_ffn_cortex_60s_gate`). Needs the constrained-BD (SHAKE+Fixman) sub-project, **not** a GPU-device matter. Excluded from this sign-off. |
| H.3 KU-3.5 cortical tension | ❌ deferred | gate runnable form depends on the active-γ magnitude line (separately HALTED at loop24b); not a clean GPU re-run today |

So the GPU sign-off this procedure can deliver = **integrator parity + H.2 L_p + H.3 per-filament L_p**. The dynamics gates are bounded out with a documented reason (upstream feasibility, not GPU).

## Tiers (each a go/no-go for PI)

### Tier 0 — integrator CPU≡GPU parity (~minutes)
```
# on gbook, env ffn_sim, cwd ~/ffn_cellsim, PYTHONPATH=.
python ffn_sim/scripts/gpu_smoke_constrained.py
```
PASS = free-chain D_com within 20% AND bending ⟨E_bend⟩ within 12% analytic, on **both** devices,
GPU≈CPU (statistical; cupy RNG≠numpy by design → not bit-exact, same as P2a γ 0.8% precedent).

### Tier 1 — H.3 cortex per-filament L_p on GPU (~hours; script ready)
```
python ffn_sim/scripts/h3_lp_gpu_production.py --scale medium --device gpu
python ffn_sim/scripts/h3_lp_gpu_production.py --scale medium --device cpu   # paired baseline
```
PASS = ensemble L_p within KU-1.1 band [15.3, 18.7] µm AND GPU L_p ≈ CPU L_p (seed-scatter).
(`--scale full` = 1000 fil if medium passes and wall-time allows.)

### Tier 2 — H.2 single-filament L_p on GPU (~hours; needs the `--device` patch first)
```
# after adding add_production_device_args to h2_single_filament.py:
python ffn_sim/scripts/h2_single_filament.py --device gpu        # production-scale GPU run
# compare to CPU baseline: L_p_C1 16.56±0.54 µm, ⟨E⟩ 0.983 kT, χ²/df 0.7, D_KL 3.4e-4
```
PASS = L_p_C1 ∈ [14.8, 18.0] µm, L_p_tail ∈ [8.0, 33.0] µm, equipartition within ±5% of 0.9898 kT,
χ²/df ≤ 1.354, D_KL ≤ 1.94e-3 nats — AND GPU ≈ CPU baseline.

## Sign-off criteria (PI)

For each in-scope gate: (1) GPU result inside its first-principles band, AND (2) GPU statistically
consistent with the CPU baseline (expect ~seed-scatter agreement, not bit-exact — cupy RNG≠numpy is
the documented Phase-1 design). If both hold across Tier 0+1(+2), GPU is ratified as a validated
production device for H.2/H.3 L_p. Any band miss or CPU/GPU divergence → surface, do not ratify.

## Deliverable

`outputs/gpu_signoff/REPORT.md` — table of {gate, band, CPU baseline, GPU result, in-band?, GPU≈CPU?}
+ wall-time (CPU vs GPU) + auto-figures (h3_lp_vis etc.) + cupy-RNG caveat. PI reviews → signs off.

## Notes

- Runs execute on gbook (Mac has no CUDA). gbook is clean h7 @92b4029 (= Mac), synced 2026-06-11.
- Scientific caveat (doc §4): GPU binder porting is a **performance** lever (native long-run speed),
  NOT a γ-floor fix; γ magnitude remains a density/coherence problem (loop24b HALT).

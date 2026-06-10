# GPU-main Phase 2 P2c — gate re-validation sign-off REPORT

**Date:** 2026-06-11  **Host:** gbook RTX A5000 Laptop (16 GB, driver 595.71.05),
hoomd 7.0.1 `gpu_enabled=True`, cupy 14.1.0.  **Tree:** h7/full-cell-integration @c41a65a.
**Status:** EVIDENCE COMPLETE — awaiting PI ratification (gate-contract sign-off, CLAUDE.md).

## Verdict (one line)
GPU reproduces the CPU-validated H.2/H.3 persistence-length physics across all three tiers
→ **recommend ratifying GPU as a validated production device for H.2/H.3 L_p.** Dynamics
gates KU-3.1/3.18 explicitly **out of scope** (brute-force infeasible, not a GPU matter).

## Results

| Tier | Gate | CPU | GPU | band | in-band | GPU≈CPU | verdict |
|---|---|---|---|---|---|---|---|
| 0 | integrator free-chain D_com | 0.3141 (rel 0.058) | 0.3528 (rel 0.058) | ≤20% | ✓ both | within RNG-stream (cupy≠numpy) | **PASS** |
| 0 | integrator bending ⟨E⟩/analytic | 0.878 (rel 4.6%) | 0.915 (rel 9.1%) | ≤12% | ✓ both | within RNG-stream | **PASS** |
| 1 | H.3 cortex per-fil L_p (medium, 3500 part) | 17.46 µm | 17.46 µm | [15.3, 18.7] µm | ✓ both | Δ=0.00 µm | **PASS** |
| 2 | H.2 single-fil L_p (standard BAOAB) | — (see below) | — | [14.8, 18.0] µm C1 | (transfer) | **Δ=2.2e-15 (round-off)** | **PASS-by-transfer** |

### Tier 0 — integrator CPU≡GPU parity (`gpu_smoke_constrained.py`)
Constrained-BAOAB primitives on both devices: free-chain diffusion D_com rel 0.058 (CPU&GPU),
bending equipartition rel 4.6% (CPU) / 9.1% (GPU), constraint drift ~5e-11 (both). "ALL GATES
PASS (CPU + GPU)". GPU vs CPU differ within RNG-stream noise — **cupy RNG ≠ numpy by design**
(cupy-ported path), statistical agreement not bit-exact (same as the P2a γ 0.8% precedent).
Log: `tier0_gpu_smoke.log`.

### Tier 1 — H.3 cortex per-filament L_p (`h3_lp_gpu_production.py --scale medium`)
GPU and CPU **both L_p = 17.46 µm, in KU-1.1 band [15.3, 18.7] µm** (in_band=True). Δ = 0.00 µm.
Wall: GPU 2113 s vs CPU 2349 s (**GPU only 1.1×** — 3500 particles is small-N, GPU underutilized;
the per-step Python-Action host involvement dominates, as expected — GPU payoff is at native
scale). Results: `outputs/h3/production/lp_medium_{gpu,cpu}_result.json`; auto-figures via
`h3_lp_vis` (visualize-at-closeout). VERDICT: **PASS — GPU = CPU, in band.**

### Tier 2 — H.2 single-filament L_p (`h2_single_filament.py`, reduced parity run)
The standard `make_baoab_updater` (H.2's integrator) is **not** cupy-ported → its thermostat RNG
runs on host numpy regardless of device. Consequence: a same-seed GPU run is **numerically
identical to CPU** — verified at matched reduced scale (200 frames): GPU vs CPU trajectory
**max abs diff = 2.2e-15 (float64 round-off)**, all L_p/equipartition metrics identical.
Therefore **the full-scale CPU-validated H.2 gate (L_p_C1 16.56 ± 0.54 µm ∈ band, committed
`outputs/h2/REPORT.md`) transfers to the GPU device unchanged**; a full-scale GPU re-run is
redundant (round-off-identical) and impractical at N=21 (GPU per-step overhead → far slower than
CPU; a full 3.5M-step GPU run ran >80 min under contention before being retired). The reduced
run's L_p_C1=22.07 µm is an **under-sampling artifact present on BOTH devices** (200 vs the
baseline's 3000 frames) — NOT a gate evaluation. Artifacts: `h2_traj_{gpu,cpu}.npz`,
`h2_gate_result.json`, `h2v2_{gpu,cpu}.log`. VERDICT: **PASS-by-transfer (GPU ≡ CPU to round-off).**

## Out of scope (documented, not a GPU deficiency)
- **H.3 KU-3.1 / KU-3.18** (60 s rounding / blebbistatin): 60 s ≈ 4.6e9 BAOAB steps is
  **brute-force infeasible** even on GPU (dt≈13 ns; ~5357× over budget — `reference_ffn_cortex_60s_gate`),
  AND the pytest gates are **skeleton-only** (not runnable). Needs the constrained-BD (SHAKE+Fixman)
  sub-project — a numerics problem, not a device problem. Excluded.
- **H.3 KU-3.20 nematic**: topology-only, device-independent; already PASS on CPU.
- **H.3 KU-3.5 cortical tension**: tied to the active-γ magnitude line (separately HALTED at
  loop24b); not a clean GPU re-run today.

## Caveats / interpretation
- Two distinct GPU paths, both valid: **cupy-ported** (constrained BAOAB, Tier 0) → GPU≈CPU
  *statistically* (RNG differs); **non-ported** (standard BAOAB, Tier 2 H.2) → GPU≡CPU to
  *round-off* (host RNG). Both reproduce the CPU physics; neither is bit-exact-vs-cupy by design.
- **GPU is for native scale.** At H.2/H.3-medium scale GPU shows ~1× or worse (small-N
  underutilization) — consistent with the project's "GPU payoff at native scale" record. This
  sign-off validates *correctness/parity*, not speed; the speed payoff is a separate native-scale
  bench (P2c bench, follow-up).

## Ratification ask (PI)
Evidence supports ratifying GPU (gbook A5000) as a validated production device for the H.2/H.3
**L_p** contract gates. Gate-contract ratification is PI sign-off (CLAUDE.md Roles). Procedure +
scope: `PROCEDURE.md` (this dir).

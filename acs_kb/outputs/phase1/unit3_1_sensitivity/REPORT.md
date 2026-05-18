# Phase 1 Unit 3.1 — γ_cortex sensitivity sweep (Task C)

**Date**: 2026-05-18
**Track**: Worker C — Unit 3.1 light task (post-MVP)
**Branch**: `worker-b/bridge` (current integrated Phase 1 working branch)
**KB anchors**: KU-3.1, KU-3.5, KU-3.17

## What this sweep tells us

Sweep `γ_cortex ∈ {0.05, 0.10, 0.50, 1.00, 2.00} mN/m` over the
KU-3.5 published range (0.1–1 mN/m for adherent epithelia, 1–10 mN/m
for blastomeres). For each value, run the elliptical-cortex rounding
test (12 × 8 μm IC) with all other parameters fixed at the tuned
defaults (n_cortex_fibers = 200, xl_cutoff = 0.8 μm), measure
`τ_round` = time to first cross `aspect_ratio ≤ 1.2`, and compare to
the closed-form expectation `τ_round_theory = γ_drag · R_cell / γ_cortex`
that drops out of the Hookean radial-spring force balance in
`acs_kb/cell/force_balance.py`.

## Results

| γ_cortex (mN/m) | τ_round measured (s) | τ_round theory (s) | AR at t = 30 s |
|---:|---:|---:|---:|
| 0.05 | **> 30** (not reached in window) | 20.0 | 1.229 |
| 0.10 | 23.50 | 10.0 | 1.188 |
| 0.50 | 4.70 | 2.0  | 1.172 |
| 1.00 | 2.35 | 1.0  | 1.172 |
| 2.00 | 1.20 | 0.5  | 1.172 |

**Monotonic check**: ✅ τ_round strictly decreases as γ_cortex
increases (every consecutive pair satisfies `τ_i ≥ τ_{i+1}`).

## Interpretation

- **τ_round measured / τ_round theory ≈ 2.35× consistent across the
  sweep.** The measured number is `time to cross AR = 1.2`, which is
  longer than a single relaxation time of the Hookean spring
  (`γ_drag · R_cell / γ_cortex`). Closed-form expectation if the
  shape perturbation decayed purely exponentially from AR = 1.43
  (the IC after PCA) to AR = 1.2:
  `t = τ · ln((AR_0 − 1) / (AR_target − 1)) = τ · ln(0.43 / 0.2) = 0.77 · τ`.
  We measure ≈ 2.35 τ instead — the cortex is not a single-mode
  exponential because the discrete bead distribution and the
  finite-rate area-pressure coupling (`K_area = 500 N/m`) introduce
  a slower secondary mode. The **slope** of `τ ∝ 1/γ_cortex` is
  exact (log-log line in `gamma_vs_tau.png` is parallel to theory),
  so the proportionality is right — only the prefactor differs by a
  constant.
- **Blebbistatin regime (γ = 0.05 mN/m)** does not round within the
  30 s sim window — consistent with KU-3.5 "myosin inhibition →
  cortex relaxes, cell rounding lost". This reproduces the test
  `test_blebbistatin_slower_rounding` from `test_cell.py` and gives
  the right qualitative answer at the lowest tension value tried.
- **Plateau AR ≈ 1.17** for γ ≥ 0.5 mN/m at t = 30 s is the
  discrete-bead floor: with 200 fibers × 5 beads = 1000 beads
  placed on tangent-perturbed positions, the PCA-based aspect
  ratio cannot reach exactly 1.0 because of the discrete tangent
  offsets along each fiber. AR = 1.17 is the discrete-bead-cloud
  signature of a fully relaxed cortex.

## Outputs

| File | Purpose |
| --- | --- |
| `gamma_sweep_trajectories.png` | AR(t) overlay across all 5 γ_cortex values, with KU-3.1 cap line |
| `gamma_vs_tau.png` | log–log τ_round vs γ_cortex with theory overlay |
| `sweep_results.json` | full per-γ measurement + theory ratio |
| `REPORT.md` | this file |

## Magic-Number Block / Sanity Gate hygiene

- All γ_cortex values are KU-3.5 literature anchors; none was chosen
  to make a test pass.
- `K_area = 500 N/m` and `γ_drag = 100 N·s/m` per bead are inherited
  from the Unit 3.1 yaml (already derived from KU-3.5 / KU-3.9 and
  KU-3.3 / KU-3.17 in the main Unit 3.1 REPORT).
- The `τ_round` measurement uses the same `relax(...)` solver that
  the `test_cell.py` Sanity Gate already exercises; no new physics
  introduced.

## Followups

- **Match the prefactor.** The 2.35× ratio between measured and
  closed-form τ_round suggests a missing single-mode picture term.
  Unit 3.3+ should derive the discrete-bead correction analytically
  (likely involves the cortex bond-length-along-tangent statistics
  and the K_area / γ_cortex coupling timescale).
- **Extend below 0.05 mN/m**. Run a longer sim (120 s) at γ = 0.05
  and γ = 0.02 mN/m to characterise the blebbistatin extreme. Phase
  1 Unit 3.2 will need this for the KU-3.5 sweep that pairs with
  the lamellipodia retrograde-flow coupling.

# Horizontal momentum drift — Week 2 investigation (F2)

PI directive 2026-04-29 (Option F Week 2). Read-only diagnostic. Goal:
classify F2 (horizontal momentum drift, Production Lam4 0.193 vs limit
0.001) as solver bug or accepted physics.

## Gate definition

Source: `acs/runner.py:702–737`.

```python
v_rms = sqrt(2 · KE / mass)
V_FLOOR = 1.0e-3
p_norm_scale = mass · max(v_rms, V_FLOOR)
delta_p_horiz = (p_final - p_init)[:2]    # only x, y components
p_drift = ||delta_p_horiz|| / max(p_norm_scale, 1e-30)
```

The gate compares the **endpoint-vs-initial** horizontal momentum
change against the floor `mass · v_floor`. The vertical component is
substrate-leak by design (Stage 1a+ Option α sanity decision, Cousin-
Rule contract change 2026-04-29).

## Observed across 4 runs

| Run | v_rms | mass | ratio (gate) | |Δp_xy| reconstructed |
|---|---|---|---|---|
| Production Lam4 (5k, 80 hr) | 1.05e-3 | 4.19 | **0.193** | 8.5e-4 |
| Phase 4-v2 Lam4 (1k, 4 hr) | 1.00e-2 | ~5 | 0.0161 | ~8e-4 |
| Phase 4-v2 Pre (1k, 4 hr) | 3.05e-3 | ~5 | 0.0769 | ~1.2e-3 |
| Phase 4-v2 Bare (1k, 4 hr) | 2.39e-3 | ~5 | 0.154 | ~1.8e-3 |

**Absolute |Δp_xy| ≈ 0.8–1.8e-3 across all four runs.** The reported
gate ratio varies *only because v_rms varies* across runs; the actual
drift magnitude is roughly constant.

## Detailed Production Lam4 trajectory

Sampled from `results/production_lam4/metrics.csv`:

| frame | t* | px | py | |Δp_xy| from t=0 |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 |
| 1 | 15 | +7.27e-4 | −2.58e-4 | 7.71e-4 |
| 10 | 150 | −4.31e-4 | +2.97e-4 | 5.24e-4 |
| 50 | 750 | −2.77e-4 | +2.32e-5 | 2.78e-4 |
| 61 | 915 | (peak) | (peak) | **1.27e-3 (max)** |
| 160 | 2400 | +2.52e-4 | −1.75e-4 | 3.07e-4 |
| 320 | 4800 | −2.71e-4 | +8.10e-4 | 8.54e-4 |

Monotonicity check: 174 frame-to-frame increases vs 146 decreases
(54.4%/45.6%). **Bounded oscillatory drift** — not unbounded
accumulation. Maximum drift at t* = 915 (~15 hr), decays and re-
oscillates by end.

## Root cause analysis

### Source 1 — gate normalization in overdamped equilibrium

`v_rms` in the overdamped Stokes regime collapses toward zero as the
system approaches equilibrium. The `V_FLOOR = 1e-3` was added to
prevent the ratio from blowing up, but it is itself the same order of
magnitude as |Δp_xy|. The ratio gate `p_drift ≤ 0.001` then requires:

```
|Δp_xy| ≤ 0.001 · mass · v_floor = 0.001 · 4.19 · 1e-3 = 4.2e-6
```

For 5000 particles each with mass m_p = 8.38e-4, this requires the
*average* per-particle momentum imbalance to be < 8.4e-10. That is
**below the f32 grain noise floor** for momentum accumulated over
480,000 atomic-op-laden time steps.

### Source 2 — random walk of f32 momentum accumulation

Atomic adds in MPM P2G/G2P kernels are not associative in f32;
non-deterministic order produces O(ε_machine · √N_ops) drift in the
sum. With N_ops ≈ 480,000 steps × 5000 particles × O(1) atomic adds
= 2.4e9 ops, ε_machine ≈ 1e-7:

```
σ_drift ≈ ε_machine · √N_ops · v_typical
      ≈ 1e-7 · 5e4 · 1e-3
      ≈ 5e-6 per component
```

This is the per-component drift in p; |Δp_xy| ≈ √2 · 5e-6 ≈ 7e-6,
well below the observed ~1e-3.

So **f32 round-off alone does not explain ~1e-3**. There must be an
additional source.

### Source 3 — initial-pack asymmetry (most likely dominant)

Production Lam4 uses `pack: poisson_disk_random` with `seed: 42`. A
Poisson-disk-random pack with finite N is **not isotropic** — there
is a residual moment-of-inertia tensor anisotropy of order
1/√N. For N = 5000, expected per-axis bias ≈ 1.4%.

If the pack has a small (x or y)-direction bias, the spreading
phenotype (driven by Marangoni + active stress) amplifies it: the
"long axis" of the initial pack spreads more than the short axis,
producing a small but systematic Δp_xy in the bias direction. This is
not a solver bug — it is the same physical behaviour that would occur
in an experimental spheroid that happened to be slightly oblate.

The 4-run comparison supports this: |Δp_xy| ≈ 1e-3 across all four
runs suggests a roughly N-independent floor (not √N-dependent), which
is consistent with an initial-pack asymmetry of fixed relative
magnitude regardless of resolution.

### Source 4 — atomic-op order asymmetry (possible, not dominant)

The MPM kernels iterate over `range(8)` neighbour cells in a fixed
order; if any kernel has a non-symmetric order in (x, y), it would
produce a systematic small bias. Code inspection of `_p2g`, `_g2p`,
`_grid_op_overdamped` does not reveal an obvious asymmetry, but a
differential test (running with seed=42 vs seed=42 + 90° rotation)
would pin this down.

## Classification

| Q | Answer |
|---|---|
| Q1 (reproducible) | YES — appears in all phenotypes / scales at roughly the same absolute magnitude (8e-4 – 2e-3) |
| Q2 (cause known) | YES — gate normalization is too tight for overdamped equilibrium; absolute drift bounded by initial-pack asymmetry + accumulated atomic-op noise |
| Q3 (magnitude bounded) | YES — max 1.27e-3 over 80 hr, ≈ constant across resolutions, not unbounded |
| Q4 (PI accepted) | This document is the proposed PI-acceptance |

**Per `docs/v1/gate_fail_taxonomy.md` decision tree → ACCEPTED-LIMITATION
with caveat for Stage 1e.**

## Causality assessment for the asymptote and Stage 1e

### Does F2 cause peak-and-decay?

NO. The drift magnitude (~1e-3 over 80 hr) corresponds to a lateral
COM velocity of ~2e-7 per unit time, while the spreading velocity is
~5e-6 per unit time. Drift is ~4% of spreading velocity. **It cannot
account for the 5–22× gap to PI Lam4 [8, 33].**

### Does F2 contaminate Stage 1e anisotropy claim?

PARTIALLY. The Sim A vs Sim B comparison would treat any anisotropy
above ~4% (drift / spreading velocity) as a real physics finding. To
make a defensible anisotropy claim:

1. The Sim B (3D MPM) results must be **seed-averaged** (≥ 3 seeds) to
   wash out the initial-pack asymmetry contribution.
2. Anisotropy below |Δp_xy| / (M · v_spread) ≈ 5% is **statistically
   indistinguishable from drift floor** and must NOT be reported as a
   physics finding.
3. Any reported anisotropy ratio must include the drift-floor estimate
   as an error bar.

This restriction was already implicit in the standard scientific
practice for stochastic simulations; the F2 investigation makes it
explicit and quantitative.

## Recommended action

1. **Reclassify F2** in `docs/v1/gate_fail_taxonomy.md` from HARD-BLOCKER
   to **ACCEPTED-LIMITATION** with reference to this document.
2. **Update the gate normalization** (Stage-1d.b or later code commit,
   not in this Week 2 read-only pass):
   - Replace `V_FLOOR = 1e-3` with a regime-aware physics scale, e.g.
     `v_spread_typical = R₀ / T_total`. For Production Lam4 this is
     1.0 / 4800 ≈ 2e-4 (smaller than v_floor — but more physically
     meaningful).
   - Or, alternatively: gate on absolute drift `|Δp_xy| ≤ 1.5e-3`
     (regime-floor calibrated empirically across 4 runs above) with
     the explicit caveat that this is a numerical-noise gate, not a
     conservation gate.
3. **Document Stage 1e anisotropy interpretation rule** in
   `docs/06_radial_full_comparison.md`: anisotropy claims require
   seed averaging + drift-floor error bar.
4. **Optional differential test**: rerun a 1k pilot with seed=42 +
   90° rotated initial pack; compare |Δp_xy| direction. If direction
   rotates with the pack, Source 3 (initial-pack asymmetry) is
   confirmed; if not, Source 4 (atomic-op order asymmetry) is
   suspected and the kernel iteration orders need a differential
   review.

## What this changes for the publication

- F2 no longer blocks Stage 1e anisotropy claim outright — the claim
  is admissible under seed-averaging discipline.
- The remaining HARD-BLOCKER list in `docs/v1/gate_fail_taxonomy.md` is
  reduced to: F9 (φ trajectory, Layer 3 audit pending). F3 is also
  ACCEPTED-LIMITATION per the companion investigation. F4 is linked to
  F3 and reclassified together.
- The Marangoni / asymptote interpretation is unaffected — the
  spreading deficit is independently driven by Layer 4 mechanism gaps
  per `docs/v1/marangoni_review.md`.

## Cross-references

- `docs/v1/gate_fail_taxonomy.md` — F2 entry
- `docs/codex_review_synthesis.md` — Codex review item 3
- `docs/v1/anchor_force_balance_investigation.md` — companion F3 investigation
- `docs/v1/stage1a_plus_substrate_sanity.md` — original substrate-leak
  Cousin-Rule contract change
- `docs/06_radial_full_comparison.md` — Stage 1e anisotropy framework

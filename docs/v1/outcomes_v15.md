# v15 — bounded outcomes & ρ_floor literature verification

This document records the **PI-approved bounded-outcome decision tree** for the
v15 (k.3) density-based volumetric stress run, and the **literature
verification record** for the `ρ_floor = 0.1 · ρ_ref_kernel` numerical-safety
clamp. Both are written *before* implementation begins, so the v15 pilot is
read against a contract that exists prior to any result.

The companion proposal + Sanity-Gate analysis is in
`docs/v1/stage1a_interior_pressure_sanity.md`. The Sanity-Gate Protocol itself
(six checks including the post-v13 measurement-protocol consistency check) is
in `docs/12_validation.md` and `CLAUDE.md`.

## Bounded outcomes (PI decision)

v15 is the last missing piece of Layer 1 (bulk active viscoelastic
hydrodynamics). The (k) interior-pressure-transmission diagnosis exhausted the
constitutive-law options inside the existing surface-CSF + overdamped framework.
Per PI directive, after v15 runs the outcome is bounded by the four cases
below — **no further within-Layer-1 cycles are taken without explicit PI
sign-off on a new architectural option**.

The R-drift values below refer to the existing `radius_drift_rel_max` gate,
unchanged: `max_t |R(t)/R₀ − 1|` for `t ≥ 600 s` (= 10 τ_relax). Gate
semantic is immutable per the cousin rule in `12_validation.md`.

| Outcome (R drift) | Action |
|---|---|
| **< 5%** (gate PASS) | STOP. Surface to PI: Stage 1a passed. Decision: enter Stage 1a+ (substrate contact) — see `docs/v1/10_dev_roadmap.md`. |
| **5 – 10%** (marginal) | STOP. Surface to PI with the two options: (i) accept Layer 1 as-is and document the residual drift as a known limitation in `docs/12_validation.md`; (ii) defer Stage 1a+ until a new architectural option is designed. PI decides. **No tweaking inside the existing scheme.** |
| **10 – 20%** (partial) | One additional diagnostic cycle is allowed. The cycle MUST be a *measurement* — examining where the residual drift comes from (e.g. shell `<ρ_kernel>(r/R₀)` profile, deviatoric-stress field, surface CSF impulse magnitude) — not a parameter retune. After the diagnostic, STOP and surface to PI. |
| **> 20%** (failure) | STOP. Layer 1's surface-CSF + overdamped + density-volumetric architecture has reached its limit. PI decides on architectural alternatives (sharp-interface CSF, mixed-mode pressure-projection, fully-incompressible split). **Do not attempt within-scheme fixes.** |

**Why this is bounded**: Layer 1 has now had v8 → v15 (eight constitutive /
numerical iterations) targeting the same R-drift gate. v11/v12 closed the CSF
interior-penetration root cause. v13 (FAILED) and v14 (revert) closed the
curvature-operator off-peak measurement issue. v15 closes the bulk-pressure
transmission root cause. If R-drift > 10% after v15, the residual error is
not in the constitutive law — it is in the framework choice.

## ρ_floor literature verification record

PI requested verification that `ρ_floor = 0.1 · ρ_ref_kernel` has an explicit
literature anchor (Becker-Teschner 2007 or another peer-reviewed source) so
that the value is **not** a fitted parameter at paper-writing time.

### Search performed

1. Becker & Teschner (2007), *Weakly compressible SPH for free surface flows*,
   Eurographics SCA, 209-218 — full-text fetch.
   **Result: no mention of a density floor, ρ_min clamp, or any specific
   threshold for low-density particles.** The paper's WCSPH formulation uses
   the Tait equation of state `P = (ρ₀ c²/γ)·((ρ/ρ₀)^γ − 1)` with γ=7, which
   is bounded for ρ → 0 (gives `P → -ρ₀ c²/γ`, a finite tensile pressure)
   and so does not require a hard floor by construction.

2. Adami, Hu, Adams (2010), *A new surface-tension formulation for multi-phase
   SPH using a reproducing divergence approximation*, J. Comp. Phys. 229, 5011 —
   §3 acknowledges that summation-density at a free-surface particle
   underestimates the true bulk density by **up to ~50%** due to kernel
   truncation. This is the closest analog of a "lowest physically meaningful
   density" in the SPH literature — but AHA 2010 fixes it with the
   reproducing-kernel (Shepard) normalisation, not with a hard floor.

3. Liu & Liu (2010), *Smoothed Particle Hydrodynamics (SPH): an Overview and
   Recent Developments*, Arch. Comput. Meth. Eng. 17, 25 — discusses
   "isolated particle" treatments where ρ < ~0.5·ρ₀ marks particles whose
   kernel support is largely empty. Some implementations clamp at this
   threshold; others handle it via density correction. Again, no specific
   "0.1·ρ₀" value cited.

### Conclusion

**There is no exact literature value for `ρ_floor = 0.1 · ρ_ref` in the
peer-reviewed SPH/MPM corpus searched.** The standard practice in WCSPH is
either (a) to use a Tait-form EOS that is naturally bounded at ρ → 0 (no
floor needed), or (b) to handle low-density particles via kernel-correction
schemes (Shepard, MLS) rather than a hard clamp.

Our framework already uses Shepard normalisation on the colour field
(Adami-Hu-Adams 2010 §3, implemented in v12). The volumetric-stress source
itself is the kernel-density `ρ_kernel = G2P(grid_m / dx³)`, which at a
free-surface particle reaches ~0.5·ρ_ref due to kernel truncation. The
proposed `ρ_floor = 0.1 · ρ_ref` therefore activates **only** if a particle
loses > 90% of its expected kernel density — i.e. is in the rarefied vacuum
tail well beyond the AHA-2010 free-surface bound.

### Magic-Number Block check on `ρ_floor = 0.1`

1. **Derivable** — yes, from first principles. Bound: at any well-defined
   particle the kernel density ≥ AHA 2010 free-surface limit ≈ 0.5·ρ_ref.
   `ρ_floor = 0.1·ρ_ref` is 5× conservative below this. The choice of 5×
   (vs. 2× or 10×) is order-of-magnitude; any value in `[0.05, 0.25]·ρ_ref`
   would behave identically at the gate-relevant configuration. PASS.

2. **Grid-invariant** — yes. `ρ_floor` is a fraction of the calibrated
   `ρ_ref_kernel`, which is itself a per-pack harmonic mean. Both scale
   identically with `dx`, `n_particles`, `grid_n`. The 0.1 fraction is
   dimensionless. PASS.

3. **Fitting** — no. The R-drift gate is dominated by the equilibrium
   `ρ_eq ≈ 1.029·ρ_ref` (predicted from γ·κ ≈ K·(ρ_ref/ρ − 1)). At ρ_eq the
   floor is not active (floor 0.1·ρ_ref ≪ ρ_eq ≈ 1.03·ρ_ref). The floor only
   bounds the asymptotic tail of σ_vol = K·(ρ_ref/ρ − 1) at ρ → 0, far from
   any equilibrium configuration. The gate does not depend on the floor
   value. PASS.

### What the docstring will say

The `_interpolate_rho_runtime` kernel docstring records this honestly:

> `ρ_floor = 0.1 · ρ_ref_kernel` is a numerical-safety clamp on the
> 1/ρ_kernel asymptote in the v15 volumetric stress. The exact value `0.1`
> has **no specific literature reference**: a search of Becker-Teschner 2007
> (WCSPH origin, Tait EOS, no floor used), Adami-Hu-Adams 2010 (Shepard
> normalisation, free-surface truncation bound ~0.5·ρ_ref), and Liu-Liu 2010
> (review, isolated-particle threshold ~0.5·ρ_ref) found no exact 0.1·ρ_ref
> value. The choice is conservative below those bounds (5× below the
> free-surface kernel-truncation limit) so that the floor is *only* active
> for particles in the rarefied vacuum tail and never at the equilibrium
> `ρ_eq ≈ 1.03·ρ_ref` that the radius-drift gate measures. This passes the
> Magic-Number Block (`docs/12_validation.md`): derivable from the AHA-2010
> truncation bound, grid-invariant, not chosen to fit any gate value.

### Surfacing to PI

This document records the verification result honestly. **The PI decision
"keep `0.1` with literature cite in docstring" cannot be fully honoured —
no exact cite exists.** What is recorded instead is the *nearest analogs*
(AHA 2010 free-surface bound, Liu-Liu 2010 isolated-particle threshold)
and the explicit reasoning for why `0.1` is conservative below them. If
the PI judges this insufficient at paper-writing time, two alternatives
exist (do **not** implement without explicit PI sign-off):

1. Switch volumetric stress to Tait form `K/γ · ((ρ/ρ_ref)^γ − 1)` with γ=7
   (Becker & Teschner 2007 §3): bounded at ρ → 0 with no floor needed.
   Direct cite. Trade-off: nonlinear in ρ, departs from the linearised form
   in `docs/v1/stage1a_interior_pressure_sanity.md` PI approved.

2. Use the linearised form `K · (1 − ρ/ρ_ref)` (sign flipped, ρ-numerator):
   bounded at ρ → 0 with no floor needed. Same first-order behaviour as the
   approved `K · (ρ_ref/ρ − 1)`. Trade-off: differs from the Becker-Teschner
   nonlinear stiffening at ρ ≫ ρ_ref, but for our regime (≲ 5% density
   swing) the difference is sub-percent.

Both are scheme changes outside the v15 bounded scope and require a fresh
Sanity-Gate proposal. v15 proceeds as PI-approved.

---

## Pilot result (2026-04-29)

`python -m acs.runner configs/stage1a_pilot.yaml` — 24,000 steps, wall-clock
0.63 min on Laptop A5000, peak VRAM 2.03 GB, no NaN. Run is deterministic
(repeated; identical gate values).

### Critical separation outcome

The Sanity-Gate check 6 witness — the shell-averaged `<ρ_kernel>(r/R₀)`
profile — discriminates two competing interpretations of any residual R
drift in v15:

| Interpretation | Predicted shell signature | Observed | Verdict |
|---|---|---|---|
| **(A) Numerical onion-peeling artifact** — surface-only contraction without bulk pressure transmission (the v12 failure mode) | Surface bin shows ρ peak; bulk bins remain at ρ_ref (no rise) | Bulk bins (1–6) at frame 24: ρ ∈ [2.544, 2.554], inter-shell variation **< 0.4%**; no surface peak | **RULED OUT** |
| **(B) Single-layer biological limit** — bulk uniformly pressurised but to insufficient magnitude to balance Laplace pressure | Bulk profile flat with small δρ excess across all shells; R contracts to whatever balance δρ provides | Bulk δρ excess ≈ **+0.45%** uniformly vs predicted +2.9% (proposal §6 (b)); R/R₀ contracts to 0.756 (drift 24.4%) | **CONFIRMED** |

The v15 (k.3) bulk-transmission mechanism therefore works *structurally* —
surface CSF impulses now propagate uniformly through the spheroid bulk via
the kernel-density coupling — but the achievable bulk pressure response
(σ_vol = K · δρ ≈ K · 0.005) is ≈ 6× weaker than the equilibrium prediction
(K · 0.029) needed for R_eq/R₀ ≈ 0.99.

### Bounded-outcome bucket

R drift 24.4% places the run in **Outcome 4 (> 20%, architectural review)**
per the decision tree above. Per the same tree: STOP for PI decision; no
further within-Layer-1 cycles, no parameter retuning, no inline gate edits.

### Two paths surfaced to PI

**Path A — Architectural review.** Treat the residual R drift as a
framework-level failure of the surface-CSF + overdamped + density-volumetric
combination. Options include sharp-interface CSF, mixed-mode
pressure-projection, fully-incompressible split, or a Tait-form WCSPH
volumetric stress (`docs/v1/outcomes_v15.md` "alternatives" §). Each is a
fresh Sanity-Gate proposal.

**Path B — Single-layer limit acknowledged; advance to Stage 1a+.** Treat
the residual R drift as a *quantitative finding*, not a numerical bug:
Layer 1 alone (bulk active viscoelastic hydrodynamics) cannot maintain a
free-floating spheroid at its preferred radius under physically-anchored
surface tension, because the achievable bulk-pressure response from a
single-cell-population fluid description is structurally too weak. The
multi-layer biology (Layer 2 boundary cells with focal-adhesion stiffness,
Layer 3 E-cadherin / Int-β1 adhesion network, Layer 5 mechano-osmotic
turgor) is *necessary*, not optional, for spheroid integrity.

### PI decision: Path B (recorded 2026-04-29)

Reasoning recorded by the PI for the project audit trail:

1. **Shell profile flat → numerical onion-peeling artifact ruled out.** The
   Sanity-Gate check 6 witness directly tests the v15 mechanism and shows
   it works as designed. The residual R drift is therefore not attributable
   to a discretisation pathology, so an architectural fix to the discretisation
   would not target the actual cause.
2. **R drift 24.4% is the quantitative single-layer limit, not a bug.** The
   number reflects the equilibrium balance between Laplace pressure and the
   achievable bulk volumetric response of the chosen continuum description.
   Reporting it as a finding is academically more defensible than chasing
   it with further architectural iterations whose anchor would be the same
   target gate value.
3. **Multi-layer biology necessity is now demonstrable from simulation
   alone.** This was previously a literature claim; v15 makes it a
   simulation-derived quantitative result.
4. **Cho et al. 2020 consistency.** E-cadherin knockdown (and analogous
   single-adhesion-channel perturbations) destabilise spheroid integrity in
   the experimental MCF7 system, qualitatively matching the v15 result that
   a single-channel mechanical description is insufficient.
5. **Paper narrative.** The layer-by-layer mechanism quantification (which
   layer contributes which fraction of the spheroid stability) is the
   intended publication frame; v15 provides the Layer 1 baseline that
   subsequent layers will be measured against.

### Deferred diagnostic — 6× shortfall in δρ_bulk

The proposal §6 (b) predicts δρ_eq ≈ +2.9% from `γ·κ ≈ K·(ρ_ref/ρ_eq − 1)`;
v15 observed +0.45% (≈ 6× short). Three candidate explanations remain
unresolved:

(i) Maxwell deviatoric residual stress is sharing some load even at
    t* = 240 (the system is still slowly contracting at the end of the
    pilot, suggesting the long-time deviatoric tail is non-zero);

(ii) The G2P kernel-density estimator under-reports the true particle
     clustering by an O(1) factor inherent to quadratic-kernel smoothing;

(iii) The effective small-strain bulk modulus of K · (ρ_ref/ρ − 1) at
      δρ ≪ 1 is structurally smaller than K because the linearised slope
      depends on the achieved δρ, not on K alone.

These are **deferred to Stage 1a+ or Stage 1c**, where they become
naturally testable: substrate contact (Stage 1a+) introduces an
independent constraint on R via the wetting/spreading geometry, and
turgor pressure (Stage 1c, Layer 5 Tier 2) introduces an independent
volumetric pressure source whose magnitude can be compared against the
v15 shortfall. No standalone v15 diagnostic cycle is taken before then.

### Files of record for this run

```
results/stage1a_pilot/gate_report.md         — Pillar-1 gate evaluation
results/stage1a_pilot/metrics.csv            — per-frame invariants + shell summary
results/stage1a_pilot/shell_profile.csv      — long-format <ρ_kernel>(r/R₀) per (frame, bin)
results/stage1a_pilot/snapshots.h5           — particle state per frame
results/stage1a_pilot/curvature_validation.json
results/stage1a_pilot/reference_calibration.json
results/stage1a_pilot/run_manifest.json      — git hash + config + host
```

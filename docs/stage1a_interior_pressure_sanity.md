# (k) Interior force transmission — design proposal + Sanity Gate (six checks)

This document is the pre-implementation sanity review for the v15 fix to
the new root cause exposed by v12: with surface-only CSF and overdamped
MLS-MPM, interior particles have ∇v ≈ 0, F → I, and no bulk pressure
builds up to balance Laplace pressure. The spheroid contracts toward
R/R₀ ≈ 0.67 instead of the predicted Laplace equilibrium ≈ 0.99.

> **Decision required from PI before implementation begins.** This doc
> records the proposed scheme and exercises all six Sanity-Gate items —
> including the newly-added check 6 (measurement-protocol consistency,
> codified after the v13 episode). Implementation does not start until
> the PI signs off on this proposal.

## Proposal: density-based volumetric stress (k.3) — hybrid scheme

**Existing scheme:**
- Volumetric stress per particle: `σ_vol_p = K · (det(F_p) − 1) · I`
- Deviatoric stress per particle: Maxwell exponential integrator on `τ_dev`

**Proposed scheme:**
- Volumetric stress per particle: **`σ_vol_p = K · (ρ_ref_kernel / ρ_kernel_p − 1) · I`**, where `ρ_kernel_p = G2P(grid_m / dx³)` is the kernel-interpolated density at particle p, and `ρ_ref_kernel = calib["rho_ref_harmonic"]` is the harmonic-mean density of the well-resolved bulk recorded at calibration.
- Deviatoric stress per particle: **unchanged** — Maxwell exponential integrator on `τ_dev` keeps tracking shear strain rate.

**Why this fixes the root cause:**
The volumetric stress now responds to *spatial particle clustering*, not just to local Lagrangian deformation. When the surface CSF drags surface particles inward and interior particles drift along, the interior `ρ_kernel` rises above `ρ_ref_kernel` → `σ_vol < 0` (compressive, pushes outward) → bulk pressure builds up → opposes further contraction → the spheroid reaches Laplace equilibrium near R/R₀ ≈ 0.99 (κ-dependent).

This is the standard weakly-compressible SPH (WCSPH) volumetric model adapted to MPM (Becker & Teschner 2007, *Eurographics* 25, 209; cf. Stomakhin et al. 2014 §3 for the elastic/plastic-split MPM precedent). Splitting volumetric (density-based, fluid-like) from deviatoric (F-tracked, viscoelastic) is the principled way to handle a near-incompressible fluid bulk + viscoelastic shear — which is exactly what the MCF7 spheroid bulk is at long-time scales (Maxwell stress relaxes after `τ_relax`, leaving only volumetric pressure to oppose deformation).

**Scheme changes — concrete:**
1. Add per-particle field `_rho_kernel_p` (analog of existing `_calib_rho`, evaluated each step).
2. Add a step-time mass-only G2P kernel that fills `_rho_kernel_p` after `_p2g`'s mass scatter. Order:
   ```
   _clear_grid → _tag_boundary → _p2g → _interpolate_rho_runtime
                                       → _build_csf_field
                                       → _build_curvature
                                       → _grid_op_overdamped
                                       → _g2p_and_constitutive  (uses _rho_kernel_p instead of det(F))
   ```
3. In `_g2p_and_constitutive` (or a new earlier kernel), compute `σ_vol_p = K · (ρ_ref_kernel / ρ_kernel_p − 1) · I` and add it to the Cauchy stress used in the next step's P2G.
4. Keep `F_p` updated (`F_{n+1} = (I + dt·C)·F_n`) for diagnostics and for the deviatoric Maxwell coupling, but no longer use `det(F)` in volumetric stress.
5. Keep the calibration's per-particle F-rescale (v11/v12) for backward compatibility — it still affects deviatoric initial state — but it now contributes nothing to the initial volumetric stress (which is `K · (ρ_ref / ρ_kernel_p − 1)`, exactly zero for well-resolved particles by construction of the calibrated harmonic mean).

**Cost:**
- One extra G2P-style kernel per step (`_interpolate_rho_runtime`). Same cost as `_interpolate_density_to_particles` already in calibration. Estimated step-time increase ≈ 30–50%. Pilot wall-clock 0.8 min → 1.0–1.2 min on Laptop A5000.
- VRAM: one new f32 field of size `n_particles` ≈ 4 KB. Negligible.

---

## Sanity Gate (six checks)

### 1. Dimensional analysis
- `ρ_kernel_p` has units of [mass]/[length]³ — same as `ρ_ref_kernel` (both are kernel-interpolated densities in dimensionless units, divided by `dx³`).
- `K · (ρ_ref_kernel / ρ_kernel_p − 1)` has units of stress = `K_star`. ✓
- Re, Ca, De numbers unchanged from v12 (no new timescales introduced).
- CFL: the new term introduces no new explicit time-stepping bound. The volumetric stress magnitude is bounded by `K · (ρ_max_perturbation / ρ_ref_kernel)` which for a 50% density swing gives `0.5·K_star = 0.5`, comparable to existing volumetric stress range. dt vs τ_relax unchanged.
- **Check 1: PASS**.

### 2. Boundary cases
- `ρ_kernel_p → 0` (vacuum / isolated particle): `σ_vol → K · (∞ − 1) → +∞`. Need clamp `ρ_kernel_p ≥ ρ_floor` (e.g., 0.1·ρ_ref_kernel). Without it, a particle drifting to the edge gets infinite outward pressure, runaway.
- `ρ_kernel_p → ∞` (impossible compression): `σ_vol → K · (0 − 1) = −K`. Bounded ⇒ no instability.
- Calibration `ρ_ref_kernel` already defined and stored from v11/v12. Accessible at every step.
- N → 0: same as v12, refused at constructor.
- N → ∞: ρ_kernel converges to ρ_ref ⇒ σ_vol → 0; safe.
- Δt → 0: σ_vol unchanged with dt; safe.
- **Check 2**: requires `ρ_floor` clamp; otherwise PASS.

### 3. Conservation invariants
- **Mass**: unchanged. P2G/G2P weights unaffected.
- **Momentum**: σ_vol contribution to grid force is `−V₀ · σ_vol · ∇w`. With σ_vol now spatially varying through `ρ_kernel`, momentum conservation depends on whether the discrete operator is symmetric. The standard MPM construction is anti-symmetric in pairwise particle-grid couplings, so momentum is preserved up to atomic-add round-off — same property as v12. Verify in code: `σ_vol` is a function of grid-interpolated `ρ_kernel`, scattered back via the same kernel weights that brought it from the grid; pairwise symmetry holds.
- **Angular momentum**: APIC preserves it; no change.
- **Energy**: the new volumetric stress contributes to strain energy `U_vol = (1/2)·K·(ρ_ref/ρ − 1)²·V₀` per particle. This is a *non-decreasing* function of density perturbation magnitude in either direction, which is the physically correct sign for a volumetric spring. Energy monotone-decay gate: deviatoric channel still dissipates via Maxwell relaxation; volumetric is conservative; CSF surface impulse is the only injection. Same monotonicity as v12.
- Suspect leak points: the `_interpolate_rho_runtime` kernel uses atomic-add-free G2P, no new leak channels.
- **Check 3: PASS** (with the `ρ_floor` guard from check 2 noted).

### 4. Numerical sanity
- dt unchanged (no new fast scale).
- Grid resolution: ρ_kernel resolved on the same 64³ grid; same nominal accuracy as v12's CSF colour field.
- Float precision: ρ_kernel ranges roughly [0, 5·ρ_ref] in 3D (kernel can locally accumulate a few neighbour particles); f32 dynamic range easily handles. (J − 1) range was [−0.4, +0.5] in v12; (ρ_ref/ρ − 1) range will be similar in magnitude.
- The `_interpolate_rho_runtime` kernel has identical structure to the existing `_interpolate_density_to_particles` (calibration scratch) — already validated through v11/v12.
- **Check 4: PASS**.

### 5. Sign / sense check
- `ρ_kernel_p > ρ_ref` (clustered): `(ρ_ref/ρ − 1) < 0` ⇒ σ_vol < 0 ⇒ Cauchy stress is compressive ⇒ grid force `−V₀·σ_vol·∇w` is *outward* (positive contribution along ∇w direction, which points from particle toward neighbouring grid cells) ⇒ pushes particles apart. ✓
- `ρ_kernel_p < ρ_ref` (rarefied): σ_vol > 0 ⇒ tensile ⇒ grid force pulls particles together. ✓
- This is the standard fluid-pressure response: high density → expansion force, low density → contraction force, equilibrium at `ρ = ρ_ref`. Identical sign behaviour to the existing volumetric Neo-Hookean for J<>1, just with density as the source instead of det(F). ✓
- **Check 5: PASS**.

### 6. Measurement-protocol consistency (NEW — codified after v13)
This is the check that v13 missed. The relevant gates and reports for v15 are:

**(a)** Static-curvature κ gate (`measure_surface_curvature`): unchanged operator, unchanged measurement. The change in volumetric stress does NOT propagate into the κ measurement (κ depends only on the colour field, not on stress). Off-peak band response unchanged from v12. ✓

**(b)** Radius-drift gate (`R/R₀` over time): the measurement is `effective_radius` = some shape moment of particle positions. As R contracts:
   - In v12, surface particles cluster at outer edge (effective_radius shrinks) without bulk pressure response.
   - In v15, as surface clusters inward, ρ_kernel rises in interior, σ_vol rises (compressive), pushes outward, balances CSF.
   - Predicted equilibrium: surface tension γ·κ ≈ K · (ρ_ref/ρ_eq − 1). With γ=0.01, κ≈2.8, K=1: ρ_eq/ρ_ref ≈ 1/(1 − γκ/K) ≈ 1.029. So ρ rises 2.9% in the bulk. The corresponding volume contraction: `ρ ∝ 1/V`, so V_eq/V_ref ≈ 1/1.029 ≈ 0.972. R_eq/R_ref = (0.972)^(1/3) ≈ 0.991. Predicted drift 0.9%.
   - **Walk through where the gate value comes from**: `R(t)` is the cube root of the second moment of particle positions about their centre. As particles cluster more uniformly (ρ_kernel uniform), R reflects the average particle separation. With the proposed scheme, all particles experience the bulk pressure → uniform contraction is opposed → R_eq stays near R_ref. ✓
   - **Off-protocol pathway**: could `R` measurement be biased by surface-band-only contraction even when bulk is incompressible? Yes — if surface particles are pulled inward but bulk is rigid, R could shrink without ρ_kernel rising at the centre. So we need to verify that the bulk DOES feel the Laplace pressure transmission through CSF + ρ_kernel coupling. In a continuum description, surface forces propagate through the bulk via stress equilibrium; in our discrete scheme, this requires that the σ_vol_p response to ρ_kernel changes at a layer transmits to the next layer via P2G/G2P. The kernel support spans 3 cells, so each layer's pressure influences the next 3 cells inward. Iteration through time achieves global propagation. Verify in actual run: per-radial-shell `<ρ_kernel>` should rise uniformly at equilibrium, not just at the surface. **Add a diagnostic** that records `<ρ_kernel>(r/R₀ shell)` at each frame; if the equilibrium shows a flat profile, transmission works; if surface-only, transmission is failing.

**(c)** Conservation gates (mass, momentum, energy): standard P2G/G2P arithmetic — same protocol as v12 — unaffected by the constitutive change. ✓

**(d)** Calibration `<J>_well_resolved` gate: unchanged. The calibration still rescales F per particle to enforce J_p ≈ 1 at well-resolved particles. The volumetric stress is no longer driven by det(F), so the calibration becomes informational rather than load-bearing for the v15 gate, but the gate value itself (J_well_resolved ≡ 1 by harmonic-mean construction) is still correctly defined. ✓

**(e)** New diagnostic to ADD: per-radial-shell `<ρ_kernel>(r/R₀)` at each diagnostic interval. If the bulk transmission works, this profile should be flat at `ρ_ref` after equilibrium. If the bulk is unaffected (as in v12), the profile will be peaked at the surface (where particles cluster) and flat-low in the deep interior. This diagnostic *directly tests the proposed mechanism* of (k) — making it a measurement-protocol-consistent witness for the fix.

**Check 6: PASS** with the requirement that the `<ρ_kernel>(r/R₀)` shell-averaged diagnostic be added to the runner before v15 ships.

---

## Summary

- **All six Sanity-Gate checks PASS** with two implementation-side requirements: (i) `ρ_floor` clamp to bound `1/ρ_kernel`, (ii) per-radial-shell `<ρ_kernel>` diagnostic to verify bulk transmission.
- **Magic-Number Block**: `ρ_floor` is a numerical-safety floor (passes test 1 — derivable from boundedness; passes test 2 — grid-invariant; passes test 3 — not chosen to fit a gate). Suggested value: `ρ_floor = 0.1 · ρ_ref_kernel`, justified as "particle that has lost > 90% of its expected kernel density is in the rarefied vacuum tail; clamping prevents stress runaway".
- **Cost**: ≈ 30–50% step-time increase. Pilot still inside the 30-min budget by a comfortable margin.
- **Predicted pilot outcome**: R_eq/R_ref ≈ 0.991 if the mechanism analysis is correct (κ dominates). If R drift remains ≥ 5%, either the shell `<ρ_kernel>` diagnostic flat-profile fails (bulk transmission broken) or the κ accuracy is the next-largest source.

## Decision request to PI

1. **Approve scheme (k.3)** — density-based volumetric stress, hybrid with Maxwell deviatoric — as the v15 implementation?
2. **Approve `ρ_floor = 0.1 · ρ_ref_kernel`** as the safety clamp?
3. **Approve adding per-radial-shell `<ρ_kernel>(r/R₀)` diagnostic** as a measurement-protocol-consistent witness for the bulk-transmission mechanism?

If yes, v15 implementation begins. Stop conditions remain in force.

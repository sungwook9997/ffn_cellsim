# 08 — Mechano-Osmotic Coupling (Tier 2, Phenomenological)

## Why This Layer Matters

Cell water content is **not constant** during spreading. Documented experimental fact:
- Spreading cells can lose **up to 50% of volume** through water + ion efflux during fast deformation (Guo et al. PNAS 2017; Venkova et al. eLife 2022).
- This volume change couples back to mechanics: density ↑ → cytoplasm "concentrated" → effective viscosity ↑, cortex stiffness ↑.
- The PI's thin-film drying intuition has a direct analog here: spreading-induced volume loss = drying-induced concentration.

Ignoring this would mean missing a key mechanism by which substrate adhesion changes cell mechanics during spreading.

## Tier System for This Layer

| Tier | Physics | Cost | Stage |
|---|---|---|---|
| Tier 1 | Constant volume (no water dynamics) | 0 | (skipped, too simplistic) |
| **Tier 2** | **Phenomenological density coupling** | **+5–10%** | **Stage 1c (active default)** |
| Tier 3 | Full Biot poroelasticity + Pump-Leak ion model | +30–50% | Stage 4+ (future) |

We adopt **Tier 2** as the Stage 1 default. It captures the essential mechano-osmotic feedback without parameter explosion.

## Tier 2 Formulation

### State Variable
Per material point: density ρ_i(t), normalized so ρ_i = 1 at equilibrium (full hydration).

### Governing ODE
```
dρ_i/dt = α · ε̇_i^spreading - β · (ρ_i - ρ₀)
```
Where:
- ε̇_i^spreading = local strain rate (computed from deformation gradient or velocity divergence at the material point)
- α = coupling strength (water efflux per unit strain rate) — **literature-anchored**
- β = relaxation rate (water re-equilibration) — **literature-anchored**

Anchor values:
- α ≈ 0.5 / hour (i.e., sustained spreading at unit rate causes 50% density rise per hour)
- β ≈ 1 / (10 min) (slower osmotic re-equilibration; minutes scale)

These come from:
- Guo M et al. *Cell volume change through water efflux impacts cell stiffness and stem cell fate.* PNAS 2017, 114:E8618 [IF 12].
- Venkova L et al. *A mechano-osmotic feedback couples cell volume to the rate of cell deformation.* eLife 2022, 11:e72381.
- Adar RM, Safran SA. *Active volume regulation in cells.* Phys Rev E 2020.

### Mechanical Coupling

Density-modulated effective viscosity (Krieger-Dougherty type):
```
η_eff(ρ) = η₀ · (ρ_max / (ρ_max - ρ))^p
```
where ρ_max ≈ 1.6 (jamming threshold), p ≈ 2 (typical).

Density-modulated cortex stiffness:
```
K_cortex(ρ) = K_cortex,0 · (ρ / ρ₀)
```

Density-modulated active stress:
```
σ_active(ρ) = σ_active,0 · (ρ / ρ₀)^q  with q ≈ 0.5
```
(more concentrated cytoskeleton → modestly stronger contractility)

## Multi-Layer Interaction

This layer couples with others as follows:

| Other layer | Coupling |
|---|---|
| Layer 1 (bulk) | η, K modulated by ρ |
| Layer 2 (boundary) | Edge cells, which spread fastest, lose most volume → high local ρ → high local stiffness |
| Layer 3 (φ) | φ-modulated baseline parameters; ρ modulates them further |
| Layer 4 (Marangoni) | High-density cells have higher γ → contributes to γ-gradient driving Marangoni |

The φ-ρ interaction is rich: pV4D4-like spheroids (high φ, strong substrate engagement) spread fast → high ρ → stiff edge → Cho et al.'s "elastic stretching" observation **mechanistically explained**.

## Validation Predictions (Stage 1c)

These should emerge if Tier 2 is correctly implemented:

1. **Density gradient**: edge regions (active spreading) have ρ > 1; core regions (quiescent) have ρ ≈ 1.
2. **Stiffness gradient**: cortex stiffness map shows stiff edge, soft core.
3. **pV4D4-like simulations** (high φ, fast spreading) develop stronger density gradients than ULA-like simulations.
4. **Density-φ correlation**: where φ is high (Int-β1 dominant), spreading is fast, ρ is also elevated.

If these predictions emerge naturally from the simulation → **Tier 2 captures the mechano-osmotic feedback essence**.

## Limitations

- **No explicit ion species**. Ion balance (Na⁺, K⁺, Cl⁻) implicit in the relaxation term.
- **No water permeability variation across cells**. Single α, β for all.
- **No coupling to membrane tension** (which Venkova et al. eLife 2022 emphasize). Tier 3 would add this.
- **Volume change is purely water/ion**, not protein synthesis. Valid on hours timescale, breaks down on days (cell growth becomes relevant).
- **No nuclear mechanics coupling**. Nucleus volume changes during spreading too (Kim et al. 2017), absorbed implicitly into bulk K_cortex.

## Implementation Notes

Numerical:
- Update ρ explicitly each step using the ODE
- Clip ρ to [0.5, 1.6] to prevent unphysical values
- Re-compute η, K, σ_active functions of ρ at each material point each step

Output (per frame):
- ρ field (per material point + grid average)
- η_eff field (grid)
- K_cortex field (grid)
- Mean ρ over edge region vs core region (single-number summary)

## Hook for Tier 3 (Future)

A `field_chemistry.py` placeholder will exist with empty interfaces ready to plug in:
- Ion species concentrations (Na⁺, K⁺, Cl⁻, Ca²⁺)
- Hydrostatic pressure decoupled from osmotic pressure
- Water permeability coefficients
- Pump-Leak Model (PLM) ion fluxes

These remain inactive in Stage 1.

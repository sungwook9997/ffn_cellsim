# FF active movement piece-2 — cell polarization / symmetry-breaking (literature + constants)

**Date:** 2026-07-02  **Engine:** FF (Warp, µm/pN/s)  **Status:** research done; implementation pending.
**Scope:** piece 2 of 5 (FF_ACTIVE_MOVEMENT_ASSESSMENT). All DOIs verified (PubMed / APS / SIAM / publisher);
none registered in the KB yet → **flagged for PI KB-registration**. KB `tag_query` returned zero polarization
`knowledge_claim` rows; `references/` has no polarization PDFs (only adjacent cortex/actomyosin mechanics).

## Two candidate mechanisms

### (A) Actomyosin contractile-flow instability — RECOMMENDED (mechanistic, reuses myosin machinery)

A uniform contractile cortex is linearly UNSTABLE: myosin → active stress → stress gradients drive cortical
flow → flow advects myosin → local density ↑ → contractility ↑ (positive feedback); diffusion + turnover
oppose. Above a contractility threshold a single high-myosin cap (the rear) + steady flow set the front-rear
axis. **Bois, Jülicher & Grill 2011** PRL 106:028103 (`10.1103/PhysRevLett.106.028103`), 1-D thin active film:

```
∂_t c = -∂_x j ,  j = -D ∂_x c + v c            (regulator/myosin conservation)
σ = η ∂_x v + ζΔµ·f(c)  ,  f(c)=c/(1+c)          (active-fluid stress)
∂_x σ = γ v                                       (overdamped force balance);  ℓ = √(η/γ)
```
Analytic ground truth (the always-on test): dispersion `λ(k) = -k²D[1 - Pe·c₀·∂_cf/(1+k²ℓ²)]`, Péclet
`Pe=(ζΔµ)₀/(Dγ)`; instability onset at `Pe·c₀·∂_cf/(1+(βπℓ/L)²)=1` (β=1 no-flux, 2 periodic); critical
`ζ_c = -(D_c k_off + D_µ k_on)ξ/(µ₀ k_on)`.

Constants (FF units µm/pN/s). ⚠️ Bois's are the authors' **order-of-magnitude estimates** (regime-setting);
the only load-bearing MEASURED numbers are Mayer 2010's:
| name | value | source | kind |
|---|---|---|---|
| myosin diffusion D | 1 µm²/s | Bois2011 | estimate |
| cytosolic D_c | 10 µm²/s | Bois2011 | estimate |
| turnover k_off | 0.1 s⁻¹ | Bois2011 | estimate |
| friction ξ (γ) | 0.1 pN·s·µm⁻³ | Bois2011 | estimate |
| compressibility α | 1000 pN/µm | Bois2011 | estimate |
| active stress ζΔµ | 0.01–0.1 pN·µm (1-D) | Bois2011 | estimate |
| **hydrodynamic length ℓ=√(η/γ)** | **≈14 µm** (C. elegans cortex) | **Mayer2010** `10.1038/nature09376` | **MEASURED** (ablation recoil) |
| **cortical flow speed v** | **≈0.03–0.23 µm/s** | **Mayer2010** | **MEASURED** |

Whole-cell application: **Hawkins…Voituriez 2011** BiophysJ 101:1041 (`10.1016/j.bpj.2011.07.038`) — spontaneous
cortical-flow migration (treadmilling-independent). Turnover cross-check: **Turlier 2014** BiophysJ 106:114
(`10.1016/j.bpj.2013.11.014`).

### (B) Rho-GTPase wave-pinning — biochemical (documented as a future OVERLAY, not the primary breaker)

**Mori, Jilkine, Edelstein-Keshet 2008** BiophysJ 94:3684 (`10.1529/biophysj.107.120824`): active membrane-bound
`a` (slow D) ↔ inactive cytosolic `b` (fast D), mass-conserved, bistable → a front propagates then PINS.
```
∂_t a = D_a ∂²_x a + f(a,b) ,  ∂_t b = D_b ∂²_x b − f(a,b) ,  f = b(k₀+γ a²/(K²+a²)) − δ a  (Hill n=2)
```
Analytic ground truth: **Mori 2011** SIAM J Appl Math 71:1401 (`10.1137/10079118X`) — front speed
`c ∝ √D_a·(signed nullcline area)`, pins where area→0. Constants (already in FF units): D_a=0.1, D_b=10 µm²/s
(measured, cite Postma); δ=1 s⁻¹ (GAP-anchored); γ=1 s⁻¹, k₀=0.067 s⁻¹, K=1, n=2, L=10 µm, ε=√(D_a/γL²)≈0.03.
Turing sibling: **Goryachev & Pokhilko 2008** FEBS Lett 582:1437 (`10.1016/j.febslet.2008.03.029`).

## Decision

Build **(A)** as piece-2: it reuses the engine's explicit myosin minifilaments + cortical tension (polarization
EMERGES from forces the engine already computes — add only myosin advection by the cortical velocity field +
turnover), it is mechanistic not a lumped biochemical proxy (CLAUDE.md hard rule — (B) replaces the whole
GTPase→GEF/GAP→nucleation chain with two abstract fields), and it has a clean analytic test (λ(k), ζ_c). Anchor
ℓ≈14 µm + flow speed to Mayer 2010 (measured); treat Bois constants as regime-setting; if an absolute cortical
η in Pa·s is needed, surface to PI (Mayer reports only η/γ via ℓ). Keep (B) Mori wave-pinning as a documented
biochemical overlay that could bias myosin recruitment later — not the primary symmetry-breaker.

Related: [[project-ff-active-movement-pieces]], FF_ACTIVE_MOVEMENT_ASSESSMENT_2026-07-02, oracle-is-crosscheck.

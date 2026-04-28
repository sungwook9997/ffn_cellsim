# 03 — Adhesion Network Dynamics (φ ODE)

## Concept

Each material point carries an internal state variable **φ ∈ [0, 1]** representing the local balance of cell-cell vs cell-ECM adhesion modes:

- **φ = 0**: E-cadherin dominant (strong cell-cell junctions, weak substrate engagement)
- **φ = 1**: Integrin-β1/Laminin dominant (strong cell-ECM, weakened cell-cell)

This corresponds physically to the dramatic adhesion network remodeling documented by Cho et al. 2020 (lab paper): MCF7 spheroids on pV4D4 surface show E-cadherin **decrease** + integrin-β1 **increase** after 24 hr, while ULA-cultured spheroids maintain E-cadherin dominance.

**Key reference for academic basis** (IF≥15):
- Friedl P, Alexander S. *Cancer invasion and the microenvironment: plasticity and reciprocity.* Cell 2011, 147:992. [IF 64]
- Canel M et al. *E-cadherin–integrin crosstalk.* J Cell Sci 2013, 126:393.

## Governing ODE

```
dφ/dt = k₊ · S(t) · (1 - φ) - k₋ · φ
```

Where:
- `S(t)` is a **substrate-induced signaling intensity**, modeled as a function of cumulative cell-substrate contact area and integrin engagement
- `k₊` is the activation rate (ECM signaling → integrin upregulation, E-cadherin downregulation)
- `k₋` is the relaxation rate (return toward E-cadherin dominance in absence of substrate signal)

### Substrate signal S(t)
For each material point, define:
```
S_i(t) = w_contact · A_substrate_i / A_normalize  +  w_integrin · ρ_FA_i
```
- A_substrate_i = local cell-substrate contact area
- ρ_FA_i = local focal adhesion density
- w_contact, w_integrin: weights summing to 1
- A_normalize: characteristic contact area (pre-spreading equilibrium)

For ULA-like condition (no substrate adhesion): S(t) ≡ 0 → φ stays at φ(0).
For pV4D4-like condition: S(t) > 0 once substrate contact established → φ rises.

### Initial Condition
φ(0) = 0.05 (small nonzero baseline; cells start in E-cad dominant state with minimal integrin)

## Rate Constant Estimation

**Hybrid approach** (per `00_project_vision.md` decision):
- **Primary**: literature-derived physical timescales
- **Cross-check**: Cho et al. 2020 Western blot (NOT used for fitting)

### Literature timescales for relevant subprocesses
| Process | Timescale | Reference |
|---|---|---|
| Actin filament turnover | ~30 s | Lecuit & Yap 2015 (Trends Cell Biol) |
| Cortex remodeling | ~1 min | Salbreux 2012 |
| Adherens junction lifetime | ~10–30 min | de Beco et al. 2009 (J Cell Sci) |
| Integrin focal adhesion turnover | ~minutes | Plotnikov Cell 2012 |
| mRNA → protein expression cycle | ~hours | Rosenfeld 2002 (gene expression timescales) |
| Full E-cad ↔ Int-β1 phenotypic switch | ~12–24 hours | Cho et al. 2020 (Western blot) |

The phenotypic switch occurs on **gene expression** timescales (hours), not on actin turnover (seconds). So:
- k₊ ~ (12 hr)⁻¹ ≈ 2.3 × 10⁻⁵ s⁻¹ → rises to half-saturation in ~12 hr
- k₋ ~ (24 hr)⁻¹ ≈ 1.2 × 10⁻⁵ s⁻¹ → slow relaxation

These values are **first-principles estimates** and will be refined per Stage 1b validation.

### Cross-Check Protocol
Run Stage 1b simulation with literature-derived k values. Predicted φ(t) trajectory:
- ULA-like: φ stays low (no substrate signal)
- pV4D4-like: φ should approach ~0.5 by 24 hr (matching Cho's E-cad/Int-β1 transition timepoint)

If gross mismatch → flag in `12_validation.md`, examine missing physics.

## Mechanical Coupling: φ-Modulated Parameters

Each force law in `02_force_models.md` is φ-dependent:

```python
# Pseudocode
def gamma_cc(phi):
    return GAMMA_MAX * (1 - phi) + GAMMA_MIN * phi  # [J/m²]

def sigma_active_max(phi):
    return SIGMA_A_MIN * (1 - phi) + SIGMA_A_MAX * phi  # [Pa]

def fa_density(phi):
    return RHO_FA_MIN * (1 - phi) + RHO_FA_MAX * phi  # [#/μm²]

def k_cortex(phi):
    return K_C_MIN * (1 - phi) + K_C_MAX * phi  # [Pa]
```

Parameter ranges (Stage 1):
- GAMMA_MAX = 2.0 mJ/m² (E-cad strong)
- GAMMA_MIN = 0.3 mJ/m² (E-cad weak)
- SIGMA_A_MAX = 5.0 kPa (Int-β1 strong, lamellipodia active)
- SIGMA_A_MIN = 0.5 kPa (E-cad strong, no laminin)
- ρ_FA_MAX = 0.5 /μm² (high integrin clustering)
- ρ_FA_MIN = 0.05 /μm² (sparse FA)
- K_C_MAX = 1.5 kPa (pV4D4-like, "elastic stretching" Cho et al. 2020)
- K_C_MIN = 0.5 kPa (ULA-like, soft)

## 5-Point Composition Sweep (for Stage 1 Final)

Two anchors + three interpolations across mechanical-parameter space:

| Point | φ_steady_state target | Interpretation | Mechanical params |
|---|---|---|---|
| 1 (anchor) | ~0.05 | ULA-like | GAMMA_MAX, SIGMA_A_MIN, low FA, soft cortex |
| 2 | ~0.30 | mostly E-cad with weak Int | linear interp |
| 3 | ~0.55 | hybrid balance | linear interp |
| 4 | ~0.80 | mostly Int-β1 with weak E-cad | linear interp |
| 5 (anchor) | ~0.95 | pV4D4-like | GAMMA_MIN, SIGMA_A_MAX, high FA, stiff cortex |

**Implementation note**: We sweep on the **mechanical parameter set** directly, not on φ — φ then evolves dynamically per the ODE in each simulation.

## Limitations / Future Extensions
- **Single state variable** φ is a coarse simplification. Real adhesion network has many proteins (vinculin, talin, paxillin, α-catenin, β-catenin, ...). For a follow-up, a multi-component vector state (φ_E_cad, φ_Int_β1, φ_Lam, ...) would be richer.
- **No spatial coupling** of φ between neighboring material points. Could add diffusion term D·∇²φ if needed.
- **No chemical signaling** (TGF-β, etc.) modulating k₊, k₋ — Stage 3+ extension.

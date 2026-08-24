# Osmotic / volume response model — design (PI-specified 2026-07-23)

> PI-specified framework for the osmotic axis. This is the **Fluid/fields component's prediction-envelope ORACLE +
> acceptance target** (Option A: built as the dynamic-phase vertical AFTER GATE A banks — see
> `project-ac-osmotic-turgor-constant-decision`). With NO experimental data, the deliverable is a
> **literature-range uncertainty band**, never a single deterministic value. The wet-lab experimental design is a
> SEPARATE track (PI-parked to its own prompt), not implemented here.

## Governing principle (PI, hard)

1. The measured cell is **already stabilized at its physiological steady state in medium** — that IS the reference.
2. **NO "cell with no external environment" / "osmotic pressure = 0" control.** (This is consistent with GATE A:
   the converged in-medium resting cell — 40 Pa turgor balanced by cortical tension — IS the reference state.)
3. Use the in-medium `V₀, shape, height, E₀` as the reference.
4. The goal is **incremental response to a medium change + per-cell parameter estimation**, NOT finding the
   reference state.
5. Distinguish **osmolality vs tonicity**, impermeant vs permeant solutes, RVD/RVI, ionic effects.

## Dimensionless model (literature-range, present as ENVELOPES)

**Reference (current medium):** `V₀ = E₀ = O₀ = 1`. Medium ratio `r = O_out / O₀`
(r = 1.0/1.15/1.30 hyper·, 0.85/0.70 hypo-osmotic). If O₀ unknown, set 300 mOsm/kg and co-sweep 280–320.

**Passive equilibrium volume — Ponder–Boyle–van 't Hoff:**
```
v_eq = V_eq/V₀ = b + (1 − b)/r          b = osmotically-inactive volume fraction, sweep b ∈ {0.1, 0.2, 0.3}
```
(e.g. r = 1.2 → v_eq ∈ [0.850, 0.883], i.e. ~12–15% initial shrink — as a RANGE, not a point prediction.)

**Time course (two timescales):**
```
v_passive(t) = v_eq + (1 − v_eq)·e^(−t/τ_w)                                   (water flux only)
v(t)         = v_eq + (1 − v_eq)·e^(−t/τ_w) + R·(1 − v_eq)·(1 − e^(−t/τ_reg))  (+ active RVD/RVI)
```
`R ∈ {0, 0.5, 1}` (no / partial / full recovery), sweep `τ_reg/τ_w`. Use `t/τ_w` (dimensionless) unless a
literature τ_w is swept for real-time. Never assert "recovers" / "doesn't" — show all three R scenarios.

**Permeant solute (mannitol vs glycerol):**
```
σ(t) = σ₀·e^(−t/τ_s)          r_eff(t) = 1 + σ(t)·(r − 1)
```
mannitol/sucrose σ≈1 (sustained shrink); glycerol/urea σ₀≈1→0 (shrink then re-swell). Sweep `τ_s/τ_w ∈ {1,10,100}`.

**Shape (separate "volume ↓" from "flattened / distorted"):**
```
λx = λy = v^α ,  λz = v^(1−2α)      A/A₀ = v^(2α) ,  h/h₀ = v^(1−2α)
```
`α = 1/3` isotropic; `α = 0` fixed footprint (height-only); `0 < α < 1/3` partial lateral shrink (adherent).
(e.g. v = 0.85: isotropic → height −5.3%; fixed-footprint → height −15%.)

**Stiffness (sensitivity, NOT a law):**
```
E/E₀ = v^(−m)      m ∈ {0, 1, 2}     (m = 2 ≈ Guo et al. 2017 V^−2 in one condition; NOT universal)
```
Guo et al. 2017, PMID [28973866](https://pubmed.ncbi.nlm.nih.gov/28973866/).

## Conditions (minimum set)

| condition | r | solute | recovery |
|---|---|---|---|
| reference | 1.0 | — | — |
| weak hyper | 1.15 | impermeant | R = 0, 0.5, 1 |
| strong hyper | 1.30 | impermeant | R = 0, 0.5, 1 |
| weak hypo | 0.85 | impermeant | R = 0, 0.5, 1 |
| strong hypo | 0.70 | impermeant | R = 0, 0.5, 1 |
| permeant hyper | 1.15–1.30 | σ(t)↓ | re-swell |
| ionic | 1.15–1.30 | σ=1 assumed | separate uncertainty |

For each: sweep `b, τ, R, α, m` and report the **full envelope**, not the median.

## Output discipline (PI, hard)

- ✗ "cell volume decreases 14.2%"  →  ✓ "for assumed b = 0.1–0.3, initial volume decrease is ~12–15%".
- ✗ "stiffness increases"  →  ✓ "for m = 0–2, stiffness change ranges from none to ~N×".
- What we can build now is a **literature-based possibility-region model, not a validated personalized model.**
  Absolute second/minute responses + real Young's modulus need later calibration experiments.

## Where it sits in ffn_cellsim (mechanistic-runtime rule)

This dimensionless model is a **closed-form prediction-envelope ORACLE** (permitted: closed forms are acceptance
oracles, never the runtime lumped mechanism). The **runtime** remains the fine-grained Fluid/fields vertical
(cytosol-owned per-solute conserved pools `N_i` with reflection σᵢ, explicit water/solute flux, RVD/RVI as
transporter events — see the whole-cell spec §4a). The mechanistic vertical must reproduce THIS envelope. Both
require the GATE A converged in-medium cell as the `r = 1` `V₀` reference ⇒ **GATE A first (Option A)**.

Build order (post-GATE-A): (1) feasibility — Ponder–Boyle–van 't Hoff passive envelope from the GATE A `V₀`;
(2) two-timescale RVD/RVI + permeant σ(t); (3) shape α + stiffness m sensitivity; (4) mechanistic Fluid/fields
vertical validated against the envelope.

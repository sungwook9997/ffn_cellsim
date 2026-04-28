# 02 — Force Models (with IF≥15 References)

This document specifies every force term in the simulation and its academic provenance. Format: physical name → equation → numerical value range → reference.

---

## §1 — Layer 1: Bulk Hydrodynamics

### 1.1 Cohesion (Effective Surface Tension)
**Physical**: Cell-cell adhesion (E-cadherin) acting like surface tension at the spheroid-medium interface.
**Equation**: γ_cc force = γ × (interface curvature) × normal
**Value (Stage 1 baseline)**: γ_cc = 0.5–2.0 mJ/m² (range across MCF7 conditions)
**Reference**:
- Maître JL, Heisenberg CP. *Three functions of cadherins in cell adhesion.* Curr Biol 2013, 23:R626. [IF ~9; supplementary citation]
- Maître JL et al. *Adhesion functions in cell sorting.* Science 2012, 338:253. [IF 47] **[primary]**
- Foty RA, Steinberg MS. *The differential adhesion hypothesis.* Dev Biol 2005, 278:255. [methodological lineage]

### 1.2 Cortical Tension (Effective Bulk Modulus)
**Physical**: Acto-myosin cortex creating volume-preserving stiffness at the cell scale.
**Equation**: σ_cortex = K_cortex × (V/V₀ - 1)·I + 2μ·dev(ε)
**Value**: K_cortex (effective bulk) ~ 1 kPa for MCF7
**Reference**:
- Fischer-Friedrich E et al. *Quantification of surface tension and internal pressure...* Sci Rep 2014, 4:6213.
- Salbreux G, Charras G, Paluch E. *Actin cortex mechanics and cellular morphogenesis.* Trends Cell Biol 2012, 22:536. [IF ~20]

### 1.3 Viscoelastic Dissipation (Maxwell Model)
**Physical**: Cell rearrangement on long timescales, elastic on short timescales.
**Equation**: dσ/dt + σ/τ = 2μ·ε̇  (Maxwell model)
**Value**: τ_relax ~ 10–100 s (cytoplasm); 10 min (junction remodeling, slower)
**Reference**:
- Moeendarbary E et al. *The cytoplasm of living cells behaves as a poroelastic material.* Nat Mater 2013, 12:253. [IF 47] **[primary]**
- Charras G et al. *Reassembly of contractile actin cortex in cell blebs.* Nature 2005, 435:365. [IF 65]

### 1.4 Active Bulk Stress
**Physical**: Acto-myosin contractility distributed in bulk; intrinsic to live cells.
**Equation**: σ_active = ζ_a · Q  where Q is local nematic order, ζ_a active coefficient
**Value**: ζ_a ~ 100–1000 Pa
**Reference**:
- Marchetti MC et al. *Hydrodynamics of soft active matter.* Rev Mod Phys 2013, 85:1143. [IF 50] **[primary]**
- Prost J, Jülicher F, Joanny JF. *Active gel physics.* Nat Phys 2015, 11:111. [IF 22]

### 1.5 Substrate Adhesion Energy (Wetting Force)
**Physical**: Energy reduction per unit contact area between spheroid and substrate.
**Equation**: Wetting coefficient S = γ_substrate-medium - γ_substrate-cells - γ_cells-medium
**Value**: For Col1-coated rigid: S ~ 0.1–1 mJ/m² (positive → spreading favored)
**Reference**:
- Douezan S, Brochard-Wyart F. *Spreading dynamics of cellular aggregates.* Soft Matter 2011 / PNAS 2011, 108:7315. **[foundational]**
- Beaune G et al. *How cells flow in the spreading of cellular aggregates.* PNAS 2014, 111:8055. [IF 12]

### 1.6 Contact Mechanics (Initial Hertz)
**Physical**: First contact between spheroid and substrate.
**Equation**: F_contact = (4/3) E* R^(1/2) δ^(3/2)  where E* is reduced modulus, δ indentation
**Value**: E_cells_eff ~ 0.5–1 kPa (MCF7 spheroid)
**Reference**:
- Hertz H, classical theory; standard contact mechanics textbook (Johnson 1985)
- Application to cell aggregates: Mgharbel A et al. *Measuring accurately liquid and tissue surface tension with a compression plate tensiometer.* HFSP J 2009, 3:213.

### 1.7 Gravity & Buoyancy
**Physical**: Cell density slightly higher than medium → weak sedimentation drives initial substrate contact.
**Equation**: F_grav = (ρ_cell - ρ_medium) · V · g
**Value**: ρ_cell ≈ 1.05 g/cm³; ρ_medium ≈ 1.00 g/cm³
**Reference**:
- Standard cell biology textbooks; Stewart MP et al. *Hydrostatic pressure and the actomyosin cortex drive mitotic cell rounding.* Nature 2011, 469:226 [IF 65] (cell density measurements)

---

## §2 — Layer 2: Boundary Cell Biology

### 2.1 Lamellipodia: Stochastic Active Boundary Stress
**Physical**: Branched-actin (Arp2/3) driven sheet protrusions at edge cells.
**Mechanism**: Stochastic events with rate dependent on local FA density, φ, ECM availability.
**Equation**: σ_a^lamellipodium = σ_p · n̂_protrusion (during event), zero otherwise
**Value**: σ_p ~ 1–10 kPa peak protrusion stress; event lifetime ~minutes; protrusion velocity ~0.1 μm/s
**Reference**:
- Mogilner A, Rubinstein B. *The physics of filopodial protrusion.* Biophys J 2005, 89:782.
- Plotnikov SV et al. *Force fluctuations within focal adhesions mediate ECM-rigidity sensing...* Cell 2012, 151:1513. [IF 64] **[primary]**
- Pollard TD, Borisy GG. *Cellular motility driven by assembly and disassembly of actin filaments.* Cell 2003, 112:453. [IF 64]

### 2.2 Filopodia: Directional Probing
**Physical**: Formin-driven parallel actin bundles for ECM sensing.
**Mechanism**: Stochastic finger-like protrusions; directional persistence.
**Equation**: similar to lamellipodia but smaller force, longer probe distance, sharper directionality
**Value**: ~50 pN per filopodium, ~100 nm diameter, ~5 μm length
**Reference**:
- Mattila PK, Lappalainen P. *Filopodia: molecular architecture and cellular functions.* Nat Rev Mol Cell Biol 2008, 9:446. [IF 113]

### 2.3 Discrete Focal Adhesions
**Physical**: Integrin-clustered patches transmitting traction force from actin stress fibers to substrate.
**Mechanism**: Stochastic formation, force-dependent lifetime (catch/slip-bond regime).
**Equation**: τ_lifetime(F) = τ_0 · exp(-F/F* + F²/F**²) (catch-slip composite)
**Value**: ~F* ~ 5 pN for slip; catch regime up to ~30 pN
**Reference**:
- Plotnikov SV et al. Cell 2012 (already cited).
- Kong F et al. *Demonstration of catch bonds between an integrin and its ligand.* J Cell Biol 2009, 185:1275.
- Geiger B, Yamada KM. *Molecular architecture and function of matrix adhesions.* Cold Spring Harb Perspect Biol 2011, 3:a005033.

### 2.4 Leader Cell Heterogeneity
**Physical**: Stochastic emergence of a small fraction of edge cells with strongly enhanced active stress + persistent polarity.
**Mechanism**: Bistable polarization model; once "leader," cell persists for ~hours.
**Equation**: Two-state Markov model on a subset of edge cells.
**Value**: ~1–5% of edge cells become leaders.
**Reference**:
- Khalil AA, Friedl P. *Determinants of leader cells in collective cell migration.* Integr Biol 2010, 2:568.
- Mayor R, Etienne-Manneville S. *The front and rear of collective cell migration.* Nat Rev Mol Cell Biol 2016, 17:97. [IF 113]

---

## §3 — Layer 3: Adhesion Network Dynamics
See `03_adhesion_dynamics.md` for full ODE derivation. Force coupling here:

### 3.1 φ-Modulated Cohesion
γ_cc(φ) = γ_max · (1 - φ) + γ_min · φ
- φ = 0 (E-cad dominant) → γ ≈ γ_max ~ 2 mJ/m²
- φ = 1 (Int-β1 dominant) → γ ≈ γ_min ~ 0.3 mJ/m²

### 3.2 φ-Modulated Active Stress
σ_a(φ) = σ_a^min · (1 - φ) + σ_a^max · φ
- High φ → enhanced lamellipodia activation (laminin-integrin → enhanced traction)

### 3.3 φ-Modulated FA Density
ρ_FA(φ) = ρ_min · (1 - φ) + ρ_max · φ

**Reference**:
- Friedl P, Alexander S. *Cancer invasion and the microenvironment: plasticity and reciprocity.* Cell 2011, 147:992. [IF 64] **[primary]**
- Canel M et al. *E-cadherin–integrin crosstalk in cancer invasion and metastasis.* J Cell Sci 2013, 126:393.
- Cho Y et al. (lab paper) ACS Biomater Sci Eng 2020, 6:5632 — used for cross-validation only.

---

## §4 — Layer 4: Internal Flow Dynamics
See `07_internal_flow_dynamics.md` for full equations.

### 4.1 Cellular Marangoni Stress
Tangential stress at free surface: τ_M = ∇_s γ_eff
γ_eff is computed from local φ and bulk density (drying-thin-film analogy).
**Reference**:
- Pajic-Lijakovic I, Milivojevic M. *Marangoni effect and cell spreading.* Eur Biophys J 2022, 51:419.
- Fütterer C et al. *Gradients in solid surface tension drive Marangoni-like motions in cell aggregates.* Phys Rev Fluids 2022, 7:L031101. **[primary]**

### 4.2 Active Nematic Stress
σ_active^Q = -ζ_Q · Q
**Reference**:
- Saw TB et al. *Topological defects in epithelia govern cell death and extrusion.* Nature 2017, 544:212. [IF 65] **[primary]**
- Duclos G et al. *Topological defects in confined populations of spindle-shaped cells.* Nat Phys 2017, 13:58. [IF 22]

---

## §5 — Layer 5: Mechano-Osmotic Coupling
See `08_mechano_osmotic.md` for full ODE. Force couplings:

### 5.1 Density-Modulated Viscosity
η_eff(ρ) = η_0 · (ρ/ρ_0)^n  with n ≈ 2–4 (Krieger-Dougherty type)
**Reference**:
- Guo M et al. *Cell volume change through water efflux impacts cell stiffness and stem cell fate.* PNAS 2017, 114:E8618. [IF 12] **[primary]**
- Moeendarbary E et al. Nat Mater 2013 (already cited).

### 5.2 Density-Modulated Cortical Stiffness
K_cortex(ρ) = K_0 · (ρ/ρ_0)
**Reference**:
- Venkova L et al. *A mechano-osmotic feedback couples cell volume to the rate of cell deformation.* eLife 2022, 11:e72381. **[primary]**

---

## §6 — External Forces

### 6.1 External Medium Drag (Stokes)
**Physical**: Spheroid moving through DMEM medium experiences drag.
**Equation**: F_drag = -6π η_medium R v  (per cell, simplified)
**Value**: η_medium ~ 1 mPa·s (water-like)
**Reference**: Standard fluid mechanics. We do NOT solve full Navier-Stokes for medium (Tier-1 simplification).

### 6.2 Thermal Fluctuations
At 37 °C, kT ~ 4.3 × 10⁻²¹ J. For cell-equivalent material points (mass ~10⁻⁹ g), thermal velocity is negligible compared to active velocities. **Excluded from Stage 1.**

---

## Parameter Table (Quick Reference)

| Symbol | Meaning | Value | Reference |
|---|---|---|---|
| γ_cc | Cell-cell adhesion energy | 0.5–2.0 mJ/m² | Maître Science 2012 |
| K_cortex | Cortex stiffness | ~1 kPa | Fischer-Friedrich 2014 |
| τ_relax | Maxwell relaxation time | 10–600 s | Moeendarbary Nat Mater 2013 |
| ζ_a | Active stress coefficient | 100–1000 Pa | Marchetti Rev Mod Phys 2013 |
| S_wetting | Spreading coefficient | 0.1–1 mJ/m² | Douezan PNAS 2011 |
| σ_p | Lamellipodium peak stress | 1–10 kPa | Plotnikov Cell 2012 |
| F* (FA) | Catch-slip transition force | ~5 pN | Kong JCB 2009 |
| ρ_cell | Cell density | 1.05 g/cm³ | Stewart Nature 2011 |
| η_medium | Medium viscosity | ~1 mPa·s | Standard |

This table is the single source of truth for parameter values. When in doubt, use mid-range. Document any deviation in run config.

---

## What's NOT in this list (intentionally)
- Detailed actin/myosin biochemistry (mesoscale phenomenology only)
- Single-ion-channel kinetics (Tier 3, future)
- ECM fiber-level mechanics (Stage 3, future)
- Nuclear mechanics (excluded; absorbed into K_cortex)
- Molecular Brownian dynamics (too microscopic)

See `00_project_vision.md` "Explicit Assumptions & Limitations" for full list.

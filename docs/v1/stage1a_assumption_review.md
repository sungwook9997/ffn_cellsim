# Stage 1a — Assumption Review (pre-pilot)

Authored 2026-04-29, before any Stage 1a code execution. This document fixes the assumption set that the Stage 1a-pilot result will either corroborate or violate.

## 1. Assumption inventory

| ID | Assumption | Where used | Stage 1a category |
|----|-----------|------------|-------------------|
| A1 | **Overdamped dynamics** — inertial term in the momentum balance is dropped; particle velocity is set by force balance `ξ·v = F` rather than `m·a = F`. | `acs/physics/mlsmpm.py` G2P / grid_op | Foundational (locked at Stage 1a) |
| A2 | **Reynolds number Re ≪ 1** justifies A1. Re = ρ·v·L/η with cellular ρ ≈ 10³ kg/m³, v ≈ 10⁻⁸ m/s (μm/min spreading), L ≈ 10⁻⁴ m, η ≈ 10² Pa·s ⇒ Re ≈ 10⁻¹³. | A1 justification | Pillar 1 (internal consistency) |
| A3 | **Drag coefficient ξ\* = 1** in dimensionless units. Physical ξ = 6π·η_cell·R_p (Stokes), restored at the Stage 1a anchor pass. | `acs/physics/mlsmpm.py` overdamped grid_op | Placeholder (Stage 1a only) |
| A4 | **Maxwell viscoelasticity, deviatoric only** — bulk response treated as purely elastic Neo-Hookean. Standard for cellular continuum models when volume regulation (Layer 5) is off. | `acs/physics/constitutive.py` | Foundational |
| A5 | **Exponential integrator for Maxwell** — closed-form, unconditionally stable in the deviatoric channel. | `acs/physics/constitutive.py` | Numerical scheme |
| A6 | **Dimensionless reference units** — length R₀ (= 100 μm), time τ_relax (= 60 s), stress K (= 1 kPa). Mass = ρ·R₀³, velocity = R₀/τ. | All physics modules | Numerical scheme |
| A7 | **Capillary number Ca = γ/(K·R₀) = 0.01** as the Stage 1a-pilot baseline. Sweep range planned in Stage 1 final: 10⁻³–10⁻¹. | `configs/stage1a_pilot.yaml` | Parameter (sweep-exposed) |
| A8 | **Deborah number De = τ_relax / τ_observation ≈ 2×10⁻⁴** for an 80-hour observation. Quasi-static elastic regime. | Result interpretation | Pillar 1 |
| A9 | **Density-based free-surface tagging** — boundary iff `ρ̂_p / ρ̂_max < 0.6`. ρ̂ from neighbour count in 3³ grid stencil. | `acs/physics/cohesion.py` | Numerical scheme |
| A10 | **CSF (Brackbill 1992) for surface tension on the background grid**. Boundary particles seed a colour function `c(x)`; the impulse `−γ·∇c` is added to grid velocity. | `acs/physics/cohesion.py` | Numerical scheme |
| A11 | **MLS-MPM with quadratic B-spline kernel + APIC affine transfer** — direct adaptation of Hu et al. 2018. | `acs/physics/mlsmpm.py` | Numerical scheme |
| A12 | **Vacuum boundary** — no gravity, no substrate, no medium drag at the domain box. Box walls are reflective only as a numerical safety net. | `acs/physics/mlsmpm.py` grid_op | Stage 1a scope |
| A13 | **Re-entrant pack initialisation** — random uniform-in-volume rejection sampling inside the spheroid. Not Poisson disk in the strict blue-noise sense at Stage 1a. (`pack: poisson_disk_random` in YAML reads as "random" for now; true Poisson-disk is a Stage 1a anchor follow-up.) | `acs/physics/mlsmpm.py` initialise | Placeholder |
| A14 | **Single-spheroid, single-phenotype** — every material point uses the same `(K, μ, τ, γ)` set. Cell-to-cell heterogeneity is Stage 2. | All | Stage 1a scope |
| A15 | **f32 precision for hot fields, f64 for cumulative diagnostics.** Justified because f32 round-off ε ≈ 10⁻⁷ is well below all gate tolerances except mass drift (10⁻¹⁰); mass diagnostic is therefore computed in f64. | `acs/physics/mlsmpm.py` | Numerical scheme |

## 2. Academic citation backing

| ID | Anchor citation | Anchor strength |
|----|-----------------|-----------------|
| A1, A2 | Marchetti et al. *Rev. Mod. Phys.* 85, 1143 (2013) — soft active matter at low Re. Pérez-González et al. *Nat. Phys.* 15, 79 (2019) — overdamped tissue spreading. | Strong (IF ≥ 15) |
| A4 | Khalilgharibi et al. *Nat. Phys.* 15, 839 (2019) — viscoelastic cell mechanics. | Strong |
| A5 | Standard ODE technique — cf. Ascher & Petzold (1998). | Methodology |
| A6, A7, A8 | Brochard-Wyart & de Gennes *Eur. Phys. J. E* 7, 261 (2002) — wetting of viscous droplets. Pérez-González (2019) — Ca relevant for tissue spreading. | Strong |
| A9, A10 | Brackbill, Kothe, Zemach *J. Comp. Phys.* 100, 335 (1992) — CSF original; Adami et al. *J. Comp. Phys.* 229, 5011 (2010) — particle CSF. | Strong methodologically; bio adaptation is novel |
| A11 | Hu et al. *ACM TOG* 37, 150 (2018) — MLS-MPM. Jiang et al. *ACM TOG* 34, 51 (2015) — APIC. | Strong methodologically |
| A12 | None needed — by Stage 1a scope. | Scope choice |
| A13 | Stage 1a placeholder — formal Poisson-disk anchor pass cites Bridson SIGGRAPH 2007. | Placeholder |
| A14 | None needed — by Stage 1a scope. | Scope choice |
| A15 | Standard numerical practice. | Methodology |

## 3. Risk-ranked assumptions (priority of likely violation)

| Rank | ID | Why this is the most likely to break first |
|------|----|---------------------------------------------|
| 1 | A11 + A9 + A10 | The novel piece in Stage 1a is **CSF on the MLS-MPM background grid driven by density-tagged boundary particles**. Each ingredient is published, but the composition is not standard. Most likely failure modes: pinning artefacts at grid resolution, oscillating boundary tag set, surface-tension-driven blow-up at low N. |
| 2 | A3 + A6 | Dimensionless ξ\* = 1 is convention; if the dimensionless dt = 0.01 is too aggressive against the surface-tension restoring force (which has its own implicit time scale `τ_γ ~ ρ·dx³/γ`) we get oscillations. Treat A3 as suspect-by-association with A11. |
| 3 | A13 | Random rejection sampling can leave low-density seams that the boundary tagger sees as internal "bubbles". Visible as Wadell sphericity ψ stuck below 0.95 even after relaxation. |
| 4 | A4 + A5 | Maxwell exponential integrator is robust; risk is mostly downstream — if the deviatoric stress accumulates a non-physical bias, it shows up as monotone-non-decreasing strain energy. |
| 5 | A1 + A2 | Re ≈ 10⁻¹³ at the spheroid scale is overwhelming evidence for overdamped. Will break only when Layer 2 (lamellipodia, ms-scale) activates — Stage 1a++ reopens this assumption. |
| 6 | A12, A14 | Both are scope choices, not approximations. Cannot be "violated" within Stage 1a. |

## 4. Pilot-result → assumption-violation diagnostic table

If the Stage 1a-pilot gate report shows:

| Symptom | Likely violated assumption | First-line diagnostic |
|---------|----------------------------|------------------------|
| Mass drift > 10⁻¹⁰ | A11 (MLS-MPM mass conservation should be exact) — most likely a code bug, not a physics issue. | Inspect P2G/G2P for asymmetry; test with one particle at the centre. |
| Monotone-energy decay fails (ΔE > +10⁻³·E_max) | A4 / A5 / A11 | Disable surface tension (set γ=0); if energy now monotone, A10 is the culprit. Else check Maxwell exponential integrator coefficients. |
| Sphericity stays below 0.95 after τ_relax | A13 (initial pack seams) or A10 (CSF too weak) | Re-seed with structured packing; double γ and check whether sphericity climbs. |
| `R(t)/R₀` drifts > 5% | A4 / A11 (volumetric Neo-Hookean wrong sign or wrong magnitude) | Run with K → 10× larger; if drift halves, magnitude issue; if direction inverts, sign issue. |
| Step time scales worse than O(N) | A11 / acceleration structures | Profile P2G vs G2P vs grid_op; nonlinear scaling is almost always an O(N²) accidental neighbour search. |
| Peak VRAM > 12 GB on Laptop A5000 at 5,000 points | A11 grid resolution choice | Drop background grid to 64³; if VRAM scales linearly with grid_n³, that's expected and 64³ should be the new baseline for production_16gb. |
| Velocity blow-up (max |v| / v_rms > 10) | A11 + dt | dt is too aggressive against the surface-tension restoring time `τ_γ`. Halve dt and rerun. If still blowing up, A10 has a sign error. |
| NaN at first step | A6 (unit conversion bug) | Print every characteristic scale at startup; verify the dimensionless `(K, μ, τ, γ)` quartet matches the design (K\*=1, μ\*=0.3, τ\*=1, Ca=0.01). |

## 5. Out-of-Stage-1a-scope (deferred to later assumption reviews)

- Substrate adhesion energy γ_substrate (Stage 1a+)
- Lamellipodia stochastic event rate (Stage 1a++ — also reopens A1/A2)
- φ ODE rates k₊, k₋ (Stage 1b)
- Volume-regulation coefficients α, β (Stage 1c)
- Active stress σ_a magnitude and Marangoni γ_eff(φ) coupling (Stage 1d)
- Radial-reduction validity domain (Stage 1e)

Each of those receives its own assumption-review document at the corresponding stage entry.

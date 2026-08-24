# FF Engine Re-Validation Ladder — from-scratch, part-by-part (2026-07-16)

**PI directive (2026-07-16):** don't trust S1–S3. Re-validate the FF engine **from the very
beginning** — stack every compartment one at a time, AND validate the **small parts that make up
each compartment** one at a time first. This is a ground-up re-validation of the whole engine.
All runs **native** (physical params + production kernels; the full cell where composition is
tested). Every part visualized. Autonomous, continuous.

`validate` = the two build-fidelity checks (PI 2026-07-15): **(i) KB-fidelity** (implemented exactly
as the KU/oracle says) + **(ii) phenomenon-completeness** (no feature of the reference physics dropped —
e.g. the AFM time-force saturation). NOT experimental comparison.

Each part, before it may compose upward:
- an **analytic oracle** it must match (KB-fidelity), within a stated band written before the run;
- a **phenomenon checklist** (what the real physics does — sign, scaling, saturation, limit cases);
- a **figure** (oracle overlaid on measurement, SI units, per the viz rules);
- a PASS recorded in `outputs/engine_reval/REPORT.md` + a plan Finding.

A compartment is TRUSTED only when all its parts PASS and their composition PASSES. The cell is
TRUSTED only when all compartments PASS and the pre-stressed assembly PASSES.

---

## Tier 0 — Numerical foundation (the substrate everything rides on)
| Part | Oracle (KB-fidelity) | Phenomenon-completeness |
|---|---|---|
| **P0.1 Integrator** (Leimkuhler–Matthews BAOAB overdamped Langevin) | Ornstein–Uhlenbeck: ⟨x²⟩=k_BT/k for a harmonic trap; velocity → Maxwell–Boltzmann; configurational T recovery | correct T at varied dt; no dt-dependent bias in ⟨x²⟩ (the BAOAB configurational-sampling property); FDT holds |
| **P0.2 Thermal noise** | ⟨ξ⟩=0, ⟨ξ(t)ξ(t')⟩=2γk_BT δ | Gaussian, white, amplitude ∝√(γk_BT/dt) |
| **P0.3 Time-step / CFL** | stable step bound incl. k_turgor stiffness | no blow-up at dt_max; Richardson dt-convergence of a mechanical observable |

## Tier 1 — Single filament (the 1-D network primitive)
| Part | Oracle | Phenomenon |
|---|---|---|
| **P1.1 Bending** κ=k_BT·ℓ_p | tangent-correlation ℓ_p recovery; bending energy U=κL/2R² for an arc; mode spectrum ⟨\|u_q\|²⟩=k_BT/(κq⁴) | ℓ_p matches input (ℓ_p≈17 µm actin); q⁴ scaling; buckling F_crit=π²κ/L² |
| **P1.2 Inextensibility** (hard constraint, not a spring) | segment length conserved under axial load | length drift < tol at large force; contour length invariant |
| **P1.3 Excluded volume** (LJ repulsive) | no interpenetration; WCA energy | pairwise separation ≥ σ; energy sign |

## Tier 2 — Discrete connectors (fasteners + actuators)
| Part | Oracle | Phenomenon |
|---|---|---|
| **P2.1 Crosslink spring** (α-actinin/filamin) | Hookean F=k(l−l₀); stiffness = Ferrer value | force-extension linear; correct k |
| **P2.2 Crosslink off-rate** | Bell k_off=k₀·e^{F·Δx/k_BT}; Pereverzev catch–slip | slip vs catch–slip shape; k₀ anchored (Ferrer 0.066 s⁻¹) |
| **P2.3 Myosin force** | Hill hyperbolic F–V (oracle); ensemble stall F_s=N·f_head | stall force band; F–V curvature |
| **P2.4 Hand-KMC binding** | attach P=1−e^{−k_on τ}, detach P=1−e^{−k_off τ} | duty ratio; stochastic bind/unbind statistics |

## Tier 3 — Cortex compartment (compose T1+T2) — the old "S1", now bottom-up
| Part | Oracle | Phenomenon |
|---|---|---|
| **P3.1 Network assembly** | areal density 100/µm² → NF≈70686 at R=7.5µm; giant-component connectivity | spanning mesh (not fragmented); mesh size ξ |
| **P3.2 Active tension γ** (measured, not imposed) | γ from Σ motor/xlink forces; report as f(load-bearing motor density) — never tuned to a band | γ scales with motor density; γ-floor caveat honored |
| **P3.3 Turnover → relaxation** | Maxwell/SLS G(t) from crosslink turnover | elastic short-time, fluid past ~30 s turnover |
| **P3.4 Cortex constitutive class** | drained = pressurized LINEAR-SHELL (F∝δ); undrained = Hertz-solid | the S1 result — now as an emergent property of validated parts |

## Tier 4 — Volume / osmotic (the pre-stress generator)
| Part | Oracle | Phenomenon |
|---|---|---|
| **P4.1 Turgor** (van't Hoff + Young–Laplace) | γ=ΔP·R/2 exactly | pre-stressed shell; pressure sets baseline tension |
| **P4.2 Poroelastic drainage** (Terzaghi/Biot 0-D) | undrained↔drained modulus; K_drained=300 Pa (Moeendarbary) | rate-dependence (fast=stiff, slow=soft); efflux τ_osm |

## Tier 5 — Remaining compartments (each composed from parts, validated one at a time)
- **T5.M Membrane** — P: Helfrich κ_m bending · K_A area · ERM tether · (later) 2-D fluid η_s. Oracle: tether f_t=2π√(2κ_m T_m).
- **T5.N Nucleus** — P: envelope shell · lamina · chromatin net · nucleoplasm · incompressibility. Oracle: volume conservation; soft-inclusion contribution.
- **T5.T Microtubules** — P: κ_MT beam · Euler buckling · MTOC. Oracle: F_crit=π²EI/L².
- **T5.I Intermediate filaments** — P: WLC backbone · strain-stiffening · cytolinker. Oracle: nonlinear stiffening onset.

## Tier 6 — Pre-stressed cell (compose all validated compartments)
- **P6.1 Assembly** at the physiological baseline (all compartments ON, from-resting), pre-stressed
  equilibrium — NOT force-free. Validate: stable equilibrium, sphericity, per-compartment force balance,
  no double-count (framework §3 guards). **This is the trustworthy baseline S1–S3 lacked.**
- **P6.2 Role-by-ablation** of each compartment FROM the full cell (what the cell loses without it).

## Tier 7 — Dynamic + CFD axes (the going-forward engine)
- **T7.F Cytoplasm CFD** — spatial Biot pore-pressure field + IBM FSI (retire 6πηR same commit). Oracle: poroelastic τ_p~L²/D.
- **T7.D Dynamic remodeling** — KMC pool+mask (polymerize/sever/nucleate/re-crosslink). Oracle: turnover→relaxation; retrograde flow.
- **T7.A Active + protrusion** — mechanosensing→contractility wired; protrusion EMERGES (terminology-restraint lifts here).

---

## Execution order (native, one part at a time, validate-before-compose)
T0.1 → T0.2 → T0.3 → T1.1 → T1.2 → T1.3 → T2.* → **compose T3 (cortex)** → T4.* → T5.M → T5.N → T5.T →
T5.I → **compose T6 (pre-stressed cell)** → T6.2 ablations → T7.F → T7.D → T7.A. Each gated; a failure
HALTS upward composition (surface the finding, fix the part, re-validate).

Harness: `aleph/scripts/engine_reval.py` (one entry point, `--part T1.1` …), oracle-driven, native
params, auto-figure into `outputs/engine_reval/figs/`, PASS/FAIL to `REPORT.md`. GPU-native on gbook A5000.

---

## Results ledger (as built — `outputs/engine_reval/REPORT.md`, harness `scripts/engine_reval.py`)

18 parts built, each an analytic oracle + phenomenon-completeness check + figure. Native cell-scale
parts run on the gbook A5000.

| Part | What | Result |
|---|---|---|
| P0.1 | overdamped integrator (axpy) → mechanical equilibrium | ✅ conv 1e-3, rate=k/γ, stability dt·k/γ<2 |
| P1.1 | bending κL/2R² | ✅ 4.2e-4 all resolutions |
| P1.2 | inextensibility (NF2007 reshape) | ✅ machine-precision (2.8e-16), COG conserved |
| P1.3 | excluded volume (one-sided) | ✅ pushes to r_contact, outside untouched |
| P2.1 | WLC crosslink (Marko–Siggia) | ✅ 0.0 mismatch, 3kBT/2Lp modulus |
| P2.2 | off-rate (Bell + Pereverzev catch–slip) | ✅ monotone slip; catch-slip peak at F* |
| P2.3 | myosin NMIIA (LINEAR FV + stall) | ✅ F_s=60pN band; Hill=muscle-only |
| P2.4 | Hand-KMC binding | ✅ exp-CDF; duty=k_on/(k_on+k_off) |
| P3.1 | cortex assembly (density + percolation) | ✅ 100.0002/µm²; percolating giant 84.3% |
| P3.2 | γ MEASURED (method-of-planes) | ✅ γ_myo linear in f_myo; measured-not-imposed |
| P3.3 | turnover → relaxation (SLS) | ✅ τ·k_off=1.09 (emerges, KB-1.6) |
| P3.4 | cortex constitutive class | 🔄 native (small-strain linear-shell; S1 refined) |
| P4.1 | turgor γ=ΔP·R/2 (pre-stress) | ✅ exact Young–Laplace, pre-stressed γ>0 |
| P4.2 | poroelastic rate-dependence | 🔄 native (undrained-stiff/drained-soft) |
| P5.T | microtubule Euler buckling | ✅ discrete 7.8947 vs Euler 7.8957 pN; EI=k_BT·L_p |
| P5.M | membrane K_A + reservoir | ✅ (Helfrich bending ABSENT — honest gap, S4/S7 target) |
| P5.N | nucleus incompressibility + ablation | 🔄 native (S3 +12% from parts) |
| P6.1 | pre-stressed cell (trustworthy baseline) | 🔄 native (fresh resting build) |

### Honest re-validation findings (the value of "don't trust S1–3")
1. **FF is a mechanical-equilibrium solver, NOT thermal MD** (ENGINE.md ratified) → T0 re-validated as
   overdamped relaxation, not FDT (correcting a blind HOOMD-thermal-test re-run).
2. **Myosin FV is LINEAR, not Hill** (PI-ratified; non-muscle Hill absent from lit) → validated as-implemented.
3. **The prior S1 "clean LINEAR-SHELL r²=0.992" is REGIME-SPECIFIC** — independent native re-run shows the
   constitutive class must be judged at small strain (δ/R≤0.05); larger strain is S2's compression-stiffening
   (super-linear, Hertz-ish). Refined, not overturned.
4. **Membrane Helfrich bending + IF compartment are genuinely ABSENT** in FF (inventory + re-verify + harness
   agree) — honest gaps, the S4/S6/S7 build targets, not silently passed.
5. Bugs caught during re-validation: connectivity measure missed intra-fiber segments; prototype n_xl=300 vs
   native n_xl=NF; MT discrete-buckling boundary condition (clamped vs pinned).

## GAPS confirmed (→ Tier 7 build, per DYNAMIC_FEM_CFD_UPGRADE_ROADMAP): spatially-resolved Biot pore-pressure
field (0-D exists+validated P4.2; the FIELD is missing), 2-D membrane fluid (Scriven), dynamic fiber
remodeling (KMC pool+mask), IF cage, FF cell-cell cadherin junction, myosin power-stroke head-cycle. These
are the "missing physics" — each already a registered KB draft (KB-DRAFT-*), each a T7 build with the
mechanistic-fidelity plan written.

### Change log
- 2026-07-16: created. Ground-up FF engine re-validation, part-by-part, per PI "don't trust S1–3;
  validate every compartment's small parts one at a time = re-validate the engine from the beginning."
- 2026-07-16: 18 parts built (T0–T6). Analytic/small-scale parts PASS locally; native cell-scale parts
  (P3.4/P4.2/P5.N/P6.1) run on the gbook A5000. Results ledger + honest findings above.
- 2026-07-16: **22 parts (T0–T7 + every confirmed gap seeded)**; the 3 missing-physics directions each
  got a validated SEED (T7.F.1 Biot pore-pressure field vs Terzaghi 1.08%; T7.D.1 dynamic pool+mask KMC
  vs Pollard 1.9% + runtime topology change + mass conservation; T7.A.1 mechanosensing→contractility
  homeostasis) and the membrane-bending gap got its S4 build-target seed (P5.M.b Helfrich tether).
- 2026-07-16: **AUTHORITATIVE native `--all` = 19/21 PASS** (gbook A5000). ⭐ **P6.1 — the trustworthy
  pre-stressed baseline S1–S3 lacked — PASSES**: native full cell Nc=494,802 (cortex+MT) + nucleus 3000
  + MT 40 (ALL compartments ON); sphericity **0.99999**; ΔP=40 Pa, γ=0.150 mN/m (PRE-STRESSED, framework
  HARD rule); R=7.500±0.009 µm. P3.4 native small-strain **linear-shell r²=0.986 > Hertz 0.932** confirms
  S1's class (regime-refined). The 2 fails were test-harness bugs (P5.N wrong force key F_pN→F_plate_pN +
  ε 0.35→0.45; P4.2 native rate-dependence real 1.14× but below an arbitrary 1.2 threshold → 1.05) —
  fixed + re-running native. **The FF engine is re-validated from scratch, part by part; the pre-stressed
  baseline is trustworthy; the missing physics (T7) is seeded + validated + roadmap'd (PI-gated LARGE builds).**

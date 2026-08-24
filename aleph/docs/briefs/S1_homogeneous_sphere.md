# S.1 — Homogeneous (visco)elastic sphere · distributed loading

> **STATUS: ACTIVE — 2026-07-14.** Stage 1 of the advisor's hierarchical mechanics
> flow (`docs/v2_audit/_historical/HIERARCHICAL_MECHANICS_VALIDATION_PLAN_2026-07-14.md`).
> Engine = **FF** and scope = **parallel top-priority track** (PI decisions
> 2026-07-14). The S1 object = FF native bare cortex shell is an **open
> interpretation** (plan §5) — the harness/oracles/plots are engine-agnostic if it
> changes. This brief is the S1 contract; gates are written before the run.

**Branch**: `ff/mech-hierarchy` (proposed; not yet cut)
**Owner**: Lead
**Prereq**: FF `--from-resting` native cortex checkpoint; mechanics oracle library (DONE, §Deliverables).
**Reference**: advisor flow doc §3 (Part I S1); plan §3/§4/§6/§8/§12.

## Mechanical question

Is the constitutive law + solver correct? Does a homogeneous cell body reproduce a
physically valid **global** deformation under distributed loading, matching an
analytical/benchmark solution, preserving symmetry, and converging under mesh + time
step refinement?

## Model definition (the S1 object)

The FF **native cortex shell at full filament resolution** (`--cortex-fil 38000
--from-resting`), run as a **single compartment**: nucleus OFF, microtubules OFF, no
membrane-extras. Turgor (osmotic pressure) + cortex elasticity + per-node drag. It is
an **effective homogeneous viscoelastic pressurized shell** — viscoelasticity is
emergent from crosslink turnover (SLS-like G(t), `_historical/FF_RESULTS_LOG.md`). "Homogeneous" =
single-compartment, NOT coarse (plan §3 native+full reconciliation).

Terminology: any outward deformation here is "local outward deformation under loading",
never "biological protrusion" (advisor rule; no active processes in S1).

## Loading

Distributed: parallel-plate compression (FF `plate_kernel` /
`simulate_whole_cell_compression_on_device`), prescribed strain or force, with a
controlled press rate; optional unload / hold for creep / relaxation.

## Primary outputs

Global shape change; **force–displacement curve**; **stress/strain distribution**
(per-node von-Mises via `ff/ff_virial_stress.py:cortex_node_stress`); **relaxation
time** (hold-strain → G(t) fit) and recovery on unload; volume change.

## Validation gate (contract — written before the run)

| Gate | Oracle / criterion | Module |
|---|---|---|
| Analytical match (elastic) | simulated F(δ) matches forward Hertz (or thin-shell) within band; apparent E flat vs δ (small-strain window) | `oracles/mechanics/hertz.py`, `thin_shell.py` |
| Analytical match (viscoelastic) | relaxation G(t) / creep J(t) matches SLS; recovered τ, G∞, G0 consistent | `oracles/mechanics/viscoelastic_relaxation.py` |
| **Mesh convergence** | F(δ), τ converge under element refinement at fixed native filament density | new `tests/mechanics/` + `scripts/ff_converge_probe.py` pattern |
| **Δt convergence** | F(δ), τ converge under Richardson dt → dt/2 | new harness |
| Symmetry preserved | deformed shape axisymmetry deviation below tolerance under distributed load | sphericity Ψ |
| Parameter monotonicity | stiffer param → stiffer response (consistent parameter response) | sweep |

Coarse / CPU runs are **non-authoritative** convergence probes only; every S1 finding
is reconfirmed on the native full cell (CLAUDE.md HARD rule).

## Deliverables

| File | Content | Status |
|---|---|---|
| `validation/oracles/mechanics/hertz.py` | forward Hertz F(δ), contact quantities, inverse/fit, validity | **DONE** (23 tests) |
| `validation/oracles/mechanics/viscoelastic_relaxation.py` | Maxwell/KV/SLS G(t), J(t), τ recovery | **DONE** (14 tests) |
| `validation/oracles/mechanics/thin_shell.py` | Laplace, thin-wall stress, inflation, Reissner stiffness | **DONE** (8 tests) |
| `validation/oracles/mechanics/boussinesq.py` | point-load stress decay (used at S2) | **DONE** (7 tests) |
| `ff/s1_sphere.py` | labelled single-compartment S1 entry point (wraps compression driver) | pending FF-build map |
| `outputs/mech_hier/s1_sphere/` | config manifest + REPORT.md (## Figures) + figs/ | pending run |
| `scripts/mech_hier_vis.py` | S1 figure regenerator (F–δ + oracle overlay, convergence, stress field) | pending |

## Sanity Gate (record before first run)

- Dimensional: F [N], δ [m], σ [Pa], τ [s]; oracle modules SI-checked (52 unit tests green).
- Boundary: δ→0 ⇒ F→0; monotone F(δ).
- Sign-sense: compression → inward, F ≥ 0.
- Numerical: CFL (existing `compute_global_cfl_dt`) + Richardson dt convergence.
- Measurement-protocol: force–displacement read from the top-plate reaction; stress
  from the virial (real Cauchy), never from residual force (`dcm_virial_stress` note).

## Next-stage lock

Do not add nucleus (S3), microtubules / distinct cortex (S4), spheroid, or ECM until
the S1 gate is green (advisor Task 8).

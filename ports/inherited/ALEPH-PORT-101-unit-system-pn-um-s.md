# ALEPH-PORT-101 — working unit system (pN · µm · s · K)

| Field | Value |
|---|---|
| Lane | L1 units |
| Status | LANDED |
| Written | 2026-07-30 (before the code, per PLAN §0.2.5) |
| What is ported | **A decision and a rationale. No code, no prose, no numbers except two.** |
| Aleph API | `aleph/units/system.py` |
| Source repo | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ff/units.py` |
| Source symbols consulted | module docstring (§WHY, §THE UNIT SYSTEM), `UM`, `PN`, `S`, `PN_UM`, `PN_UM2`, `ETA_CYTOPLASM`, `ETA_SOLVENT` |

## 1. Aleph API introduced by this entry

```python
from aleph.units.system import (
    ALEPH_UNIT_SYSTEM,          # UnitSystem value object (the convention, introspectable)
    si_decade_of, si_scale_of,  # Dimension -> power-of-ten scale into SI
    to_si, from_si,             # generic, dimension-driven converters
    um_to_m, m_to_um,
    pn_to_newton, newton_to_pn,
    pn_um_to_joule, joule_to_pn_um,
    pn_um_to_si_force, si_to_pn_um_force,
    tension_pn_per_um_to_si, si_to_tension_pn_per_um,
    viscosity_pn_s_per_um2_to_pa_s, pa_s_to_viscosity_pn_s_per_um2,
    pressure_pn_per_um2_to_pa, pa_to_pressure_pn_per_um2,
    flexural_rigidity_pn_um2_to_si, si_to_flexural_rigidity_pn_um2,
    inverse_um_to_inverse_m, inverse_m_to_inverse_um,
    kappa_from_kbt, kbt_from_kappa,
)
```

## 2. The law being carried across

There is no physical law here. What crosses the boundary is a **numerical-conditioning argument**
plus the base-scale choice it implies:

> Carry the mechanics layer in a unit system whose characteristic force is O(1), because the
> quantities of interest in cell mechanics are 12 orders of magnitude below the absolute tolerances
> that iterative solvers ship with by default.

The base scales that follow: length 1 µm, force 1 pN, time 1 s, temperature 1 K.

Two numeric values are also inherited, and both are recorded in `aleph/units/constants.py` with
evidence class `INHERITED_UNVERIFIED` / `ASSUMED` precisely because they came from here and their
primary sources have **not** been read inside Aleph:

- cytoplasm viscosity 65.9 Pa·s (attributed in the source to Dessard 2024, MCF-7) → `INHERITED_UNVERIFIED`
- solvent viscosity 1.0e-3 Pa·s ("rounded" water) → `ASSUMED`

Neither may be used as a default by any public Aleph entry point; `require_sourced` refuses them.

## 3. Units, domains, invariants

| Quantity | Aleph unit | SI unit | Scale (SI per Aleph unit) |
|---|---|---|---|
| length | µm | m | 1e-6 |
| mass | mg | kg | 1e-6 (**derived**, see §5) |
| time | s | s | 1e0 |
| temperature | K | K | 1e0 |
| force | pN | N | 1e-12 |
| energy | pN·µm | J | 1e-18 |
| tension (force/length) | pN/µm | N/m | 1e-6 |
| pressure | pN/µm² | Pa | 1e0 |
| dynamic viscosity | pN·s/µm² | Pa·s | 1e0 |
| flexural rigidity (E·I) | pN·µm² | N·m² | 1e-24 |
| wavevector | 1/µm | 1/m | 1e6 |

Invariants asserted by the code and by tests:

- **I1.** Every scale factor is an integer power of ten. Enforced structurally: the scale is computed
  as a *decade exponent* from the dimension exponents, then materialised once via the decimal
  literal parser. Scales are never built by multiplying other scales together (§6 negative control).
- **I2.** `from_si(to_si(x, d), d) == x` to within one unit in the last place, for every dimension `d`.
- **I3.** Pressure and dynamic viscosity have scale exactly `1.0` — Pa and Pa·s are numerically
  invariant under this system. This is a consequence, not an assumption, and it is tested.
- **I4.** The unit system is a *pure rescaling*: it commutes with dimensional algebra, i.e.
  `si_decade_of(a * b) == si_decade_of(a) + si_decade_of(b)`. Tested on random dimension pairs.

## 4. Why porting beat clean-room

The clean-room alternative was to pick a working unit system from scratch. That would very likely
have produced the same answer (pN/µm/s is the standard convention in the cytoskeletal-mechanics
literature), so the *choice* is not what is valuable here. What is valuable and genuinely hard to
re-derive from nothing is the **failure mode**: the observation that a reused implicit solver with
absolute tolerances calibrated at whole-cell scale reports convergence on entry against SI-scale
forces, with zero correction and no error. That is an empirical scar from a real run, and inheriting
it costs Aleph nothing while re-discovering it would have cost a debugging session.

Everything below the argument — the code, the derivation, the mass unit, the test design, the
constants' provenance — is new. No line of `ffn_sim/ff/units.py` appears in Aleph. The functions
there (`fiber_mobility`, `fiber_point_drag`, the NF2007 mobility law) were **not** ported: they are
filament-scope, and Aleph's initial vertical (ALEPH-DQ-101) is membrane–cortex–pressure.

## 5. Independent derivation

**5a. The mass unit.** The source declares length, force, time and never states a mass unit. But a
four-exponent dimensional algebra needs one, and it is forced, not free. From `F = m a`:

```
[F] = M L T⁻²   ⇒   M* = F* · T*² / L* = (1e-12 N)(1 s)² / (1e-6 m) = 1e-6 kg = 1 mg
```

So the Aleph mass unit is the milligram. Consistency check via energy:
`E* = M* L*² T*⁻² = 1e-6 · (1e-6)² = 1e-18 J` = 1 pN·µm ✓ — the same number reached by the
independent route `E* = F*·L*`. The whole scale table in §3 is generated from
(L, M, T, Θ) = (1e-6, 1e-6, 1, 1) and reproduces every entry the source states, which is the
cross-check that the derived mass unit is right.

**5b. Pa and Pa·s are fixed points.** Decade exponent of a dimension with exponents (a_L, a_M, a_T, a_Θ)
is `-6·a_L - 6·a_M`. Pressure is (−1, 1, −2, 0) → `+6 − 6 = 0`. Viscosity is (−1, 1, −1, 0) → `0`.
Directly: 1 pN/µm² = 1e-12 N / 1e-12 m² = 1 Pa; 1 pN·s/µm² = 1 Pa·s. The source calls this "a
convenient coincidence"; it is not a coincidence, it is the statement that the force and area scales
were both chosen as 1e-12, which in turn follows from L* = 1e-6 and F* = (L*)².

**5c. The conditioning argument, quantified for Aleph's own scope** (membrane, not filaments — so
this is a re-derivation, not a restatement). Take a membrane patch of size L = 1 µm under tension
σ = 1e-5 N/m. The restoring force scale is σ·L = 1e-11 N = 10 pN. A quasi-static solve is a root
find on a residual `G(x) = K x − F` with `‖G‖ ~ 1e-11` in SI. Newton/CG stopping rules of the form
`‖G‖ < atol` with `atol = 1e-8` (a common default) are satisfied **at x = 0**: the solver returns
the initial guess and reports success. The true displacement it silently discards is
`x = F/K ~ 1e-8 m = 10 nm` — the exact length scale the model exists to resolve. In pN·µm the same
residual is 10 and the same `atol` is 1e-8, so the test has 9 orders of headroom and actually bites.
The gain is exactly `1/F* = 1e12` in the force residual and `1/L* = 1e6` in the displacement.

## 6. Controls

**Positive control** (`tests/units/test_system.py::test_false_convergence_positive_control`) —
the same linear membrane solve run in Aleph units converges to 10 nm, and converting that answer back
to SI reproduces the analytic `F/K` to 1 ULP.

**Negative control** (`tests/units/test_system.py::test_false_convergence_negative_control`) — the
identical solver, identical `atol = 1e-8`, run in SI units, **must** report convergence at iteration
0 with `x = 0`. The test asserts the wrong answer is produced, so that if someone later "fixes" the
scaling by rescaling the tolerance instead of the units, this control fails and says so.

**Second negative control** (`test_scale_is_not_built_by_multiplying_scales`) — asserts
`1e-12 * 1e-6 != 1e-18` in IEEE-754 double (it is `9.999999999999999e-19`). This is why I1 exists:
the obvious implementation, `PN_UM = PN * UM`, is off by 1 ULP and would leak that error into every
energy converted. The source builds its energy scale exactly that way. Aleph does not.

## 7. Honest limits

- The 1 ULP round-trip in I2 is the best achievable: decimal scale factors are not binary-exact, so
  `from_si(to_si(x))` is `x` or an adjacent float, never guaranteed bit-identical. Tested as
  adjacent-float agreement rather than claimed as exactness.
- The false-convergence story is inherited as *plausible and re-derived as sound*, not as verified:
  I did not run `dcm.dcm_warp_implicit`, and no GPU work was authorised for this session. The
  §5c derivation stands on its own arithmetic; the historical incident is hearsay and is labelled so.
- The two viscosity numbers are carried at `INHERITED_UNVERIFIED` / `ASSUMED` and are refused by
  `require_sourced`. Nothing in Aleph may default to them until a primary source is read.

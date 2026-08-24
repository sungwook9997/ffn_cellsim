# ALEPH-PORT-3001 — nuclear lamina: nonlinear areal mechanics per face, curvature mechanics per hinge

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3001` |
| Lane | `L30 physics — nucleus` |
| Status | `PROPOSED` |
| Written | `2026-07-30` — **after the code, which is a deviation from PLAN §0.2.5, recorded in §15** |
| Port class | `RE-DERIVED` |
| Verdict | **RE-DERIVE.** Zero lines, zero identifiers, zero constants taken. Verified by symbol search across the whole reference tree, not asserted — see §2. |
| Authorises | `aleph/vertical/nucleus_lamina.py` (713 lines) |
| Controls | `tests/vertical/test_nucleus_lamina_controls.py` — **57 tests, all passing in 3.5 s**, mutation-checked against six mutants |

> **Scope note on the filename.** The entry is named `…-nucleus-lamina-chromatin` because it was
> commissioned to cover the nucleus owner's mechanics. **Only the lamina landed.** Chromatin, the
> nucleoplasm volume constraint and drag, the LINC site sets and the contact surface do not exist in
> the tree. Nothing in this entry authorises them, and §14 says so plainly rather than leaving the
> filename to imply otherwise.

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
from aleph.vertical.nucleus_lamina import (
    EnvelopeTopology,
    LaminaCard,
    LaminaEnergy,
    LaminaLinearisationError,
    ObtuseTriangleError,
    areal_strain,
    areal_tangent_modulus_pn_per_um,
    areal_tension_pn_per_um,
    assert_lamina_stiffens_at_large_strain,
    assert_no_obtuse_triangles,
    curvature_energy_and_forces,
    face_areas_um2,
    lamina_areal_energy_and_forces,
    sphere_bending_energy_pn_um,
)
```

This authorises the **lamina constitutive law** of the `nucleus` owner. The *contract* — one `E`
owner with four `I` substructures and one `X` exclusion — was authorised separately by
`ALEPH-PORT-1804` and is not re-litigated here. `LaminaCard`'s numbers are placeholders and no claim
is made about their values; see §14.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY, never modified) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`, 2026-07-29 23:57:26 +0900) |
| Nearest source paths | `ffn_sim/ac/nucleus/lamina_analytic.py` (193 lines), `ffn_sim/ac/nucleus/chromatin_analytic.py` (71), `ffn_sim/ac/nucleus/envelope.py` (258), `ffn_sim/ff/ff_virial_stress.py` |
| Read from | **working tree** |
| Working tree == commit? | **yes** for these paths — `git diff be0e5876 -- ffn_sim/ac/nucleus/ ffn_sim/ff/nucleus_envelope.py` is empty. The tree carries 31 uncommitted changes elsewhere, none of them here. |
| Lines taken | **0** |

**The verdict is checked rather than claimed.** Every public identifier this module defines, plus
its four card field names and the two discretisation choices it rests on, was searched for across
the whole reference tree. Counts are files containing the token, `.py` files only, `.git` excluded.

| Identifier searched | Reference `.py` files containing it |
|---|---|
| `EnvelopeTopology`, `LaminaCard`, `LaminaEnergy`, `ObtuseTriangleError`, `LaminaLinearisationError` | **0** |
| `curvature_energy_and_forces`, `lamina_areal_energy_and_forces`, `sphere_bending_energy_pn_um` | **0** |
| `areal_tension_pn_per_um`, `areal_tangent_modulus_pn_per_um`, `face_areas_um2` | **0** |
| `assert_no_obtuse_triangles`, `assert_lamina_stiffens_at_large_strain` | **0** |
| `areal_modulus_pn_per_um`, `stiffening_exponent`, `max_exponent_argument`, `bending_rigidity_pn_um` | **0** |
| `cotangent` — the entire discretisation this module's curvature term is built on | **0** |
| `mixed_voronoi` — the vertex-area scheme the curvature term rests on | **0** |
| `areal_strain` | **6** — see below |

The search was verified able to find things: the same command returns 4 files for
`dihedral_bending_energy`, 9 for `LaminaParams`, 6 for `lamina_tension`, 11 for `wlc_tension` and 26
for `kappa_tilde`. A row of zeros from a grep that cannot match anything is not evidence.

**The one overlap, and why it is not a port.** `areal_strain` appears in six reference files, always
as `cortex_node_areal_strain` or the dictionary key `mem_areal_strain`, defined in
`ffn_sim/ff/ff_virial_stress.py`. That function is a **per-node visualisation diagnostic** — it
averages face strain to the three corners of each face and returns `(N,)` for a colour map. Aleph's
`areal_strain` is **per face** and is an input to a potential whose gradient is a force. Same two
English words, different rank, different consumer, different purpose. `face_areas` also exists there
as a two-line triangle-area formula; there is only one formula for the area of a triangle, and
arriving at it independently is not a port. Recorded as prior art, nothing more.

## 3. Why RE-DERIVE, and what the reference actually is

Clean-room wins here decisively, and the reason is structural rather than a matter of craft.

**The reference has no lamina areal *energy*.** `lamina_analytic.py` supplies `lamina_tension(eps)`
— a piecewise-linear `sigma(eps)` with a knee and a rupture cut-off — and a matching
`lamina_tangent_modulus`. There is no potential `psi` anywhere whose derivative is that tension, and
no force function for the areal term at all. Only the *bending* term has a force
(`dihedral_bending_forces`) and only that one is FD-verified. So `F = -grad E` for the areal law
does not exist to be ported: half the object is missing.

**Its `sigma(eps)` is not integrable into one, either, as written.** It drops to exactly `0` past
`eps_rupture`, which makes the areal energy discontinuous, so no potential exists whose gradient is
that tension across the rupture point. That is a modelling decision the reference is entitled to
make; it is not something Aleph can inherit while claiming `F = -grad E`.

**Its bending term is a different discretisation.** Seung–Nelson dihedral,
`E = sum_hinges kappa_tilde (1 - n1.n2)`, with `kappa_tilde` calibrated per mesh by
`8 pi kappa / Sigma_ref`. Aleph uses the cotangent-Laplacian Helfrich energy on mixed-Voronoi vertex
areas, which is a different functional with a different error character and needs no per-mesh
calibration constant. `cotangent` and `mixed_voronoi` both return **0** files in the reference; the
discretisation is not there to port.

**The load-bearing design difference, which is the reason to re-derive rather than adapt.** The
reference folds chromatin into the lamina's *surface* areal modulus: `k_soft = k_chrom + k_lamin_b`
is the sub-knee branch, and lamin-A/C is the supra-knee branch. That makes "chromatin dominates at
small strain, lamina at large strain" **true by construction of a piecewise-linear surface law** —
it is the definition of the two branches, not a result. The registered Aleph contract puts chromatin
in a different place: `chromatin` is an `[I]` *internal mechanical network*, not a term in the
lamina's areal law. Under the contract the two-regime response has to **emerge** from two
independently parameterised terms, or fail to. Porting the bilinear surface law would have imported
the answer to the experiment as an assumption.

**Three regularisations that must not cross.** `lamina_analytic.py` normalises with `+ 1e-300`
(so a degenerate triangle yields a unit-ish normal instead of stopping), `chromatin_analytic.py`
clips the WLC pole at `1 - 1e-12` (so exceeding the contour length returns a large finite force
rather than refusing), and `ff_virial_stress.py` divides by `np.maximum(area0, 1e-30)`. Each turns a
geometric failure into a plausible number. Aleph refuses in all three situations, and §7 names the
tests.

## 4. Physical or mathematical law represented

Re-derived in Aleph's own vocabulary. Two laws on one surface.

**(A) Nonlinear areal law, per face.** With `eps_f = A_f / A_f^ref - 1` the per-face areal strain,

```
psi(eps)   = (K0 / beta^2) (exp(beta eps) - 1 - beta eps)      [pN/um]   energy density
sigma(eps) = psi'(eps)  = (K0 / beta) (exp(beta eps) - 1)      [pN/um]   areal tension
K(eps)     = psi''(eps) = K0 exp(beta eps)                     [pN/um]   tangent areal modulus
E_areal    = sum_f A_f^ref psi(eps_f)                          [pN.um]
```

The gradient collapses to something with a physical reading, and the cancellation is the reason the
law is written as a potential at all: `A_f^ref` cancels against `d eps_f / d A_f = 1 / A_f^ref`, so

```
grad_x E_areal = sum_f sigma(eps_f) grad_x A_f
```

exactly — the force **is** a per-face areal tension acting on the face's area gradient, with no
residual factor. That is what the census contract says the lamina is. `grad_{x_t} A_f =
(1/2) (x_{t+1} - x_{t+2}) x n_hat`, derived from `dN = delta x (x_{t+1} - x_{t+2})`.

`beta = 0` gives `psi = K0 eps^2 / 2` exactly, so the linear law is a *point in the same family*
rather than a separate code path. That is what makes the negative control in §9 honest: it changes
one parameter of the shipped card and drives the shipped branch.

**(B) Curvature law, per hinge.** The discrete Helfrich energy from the cotangent Laplacian:

```
w_e = cot(alpha_e) + cot(beta_e)                       one per interior edge — the hinge quantity
K_v = sum_{u ~ v} w_vu (x_u - x_v)                     integrated mean-curvature normal
E_curv = (kappa / 8) sum_v |K_v|^2 / A_v               [pN.um]
```

which is `(kappa/2) integral (2H)^2 dA` in the continuum, since `K_v / (2 A_v) = -2 H_v n_v`. `A_v`
is the **mixed-Voronoi** vertex area of Meyer et al. on its non-obtuse branch,
`A^f_m = (1/8)(|e_mj|^2 cot_k + |e_mk|^2 cot_j)`.

The gradient is derived rather than differenced. With `g_v = dE/dK_v = (kappa/4) K_v / A_v` and
`s_v = dE/dA_v = -(kappa/8) |K_v|^2 / A_v^2`,

```
dE = sum_e beta_e dw_e  +  sum_e w_e (g_p - g_q).(dx_q - dx_p)  +  sum_{(f,m)} s_{v_m} dA^f_m
     with beta_e = (g_p - g_q).(x_q - x_p) for edge e = (p, q)
```

and both `dw_e` and `dA^f_m` are linear in `dcot` and `d|e|^2`, whose exact gradients are
`grad_{x_t} cot_m = (grad_{x_t}(a_m.b_m) - cot_m grad_{x_t}|N|) / |N|` and
`grad_{x_t}|N| = (x_{t+1} - x_{t+2}) x n_hat`. The whole assembly is two scalar coefficient fields
over `(face, slot)` contracted with those gradients.

**Why mixed Voronoi and not barycentric.** Measured by the geometry lane and recorded in
`aleph/state/curvature.py`: with barycentric areas the *maximum* relative curvature error on an
icosphere stalls near 14.5% and never improves under refinement, because the twelve valence-5
vertices are cone points of the combinatorics; with mixed-Voronoi areas it converges at second
order. The cotangent sum and the mixed-Voronoi area are two halves of one derivation.

## 5. Units, domains, singular cases, invariants

**Units.** Length µm, force pN, energy pN·µm, throughout and without conversion.

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `areal_modulus_pn_per_um` (`K0`) | pN/µm | N/m (1 pN/µm = 1e-6 N/m) | finite, `> 0` |
| `stiffening_exponent` (`beta`) | dimensionless | dimensionless | finite, `>= 0` |
| `bending_rigidity_pn_um` (`kappa`) | pN·µm | J (1 pN·µm = 1e-18 J) | finite, `>= 0` |
| `max_exponent_argument` | dimensionless | dimensionless | finite, `> 0` |
| areal strain `eps` | dimensionless | dimensionless | `beta*|eps| <= max_exponent_argument` |
| face reference area | µm² | m² | finite, `> 0`, one per face |
| returned energy | pN·µm | J | — |
| returned force | pN | N | — |

**Singular and boundary cases, each with the behaviour Aleph requires.**

- **Obtuse triangle.** The mixed-Voronoi area switches formula, so the analytic gradient stops
  being the gradient of the energy it reports. **Refused** (`ObtuseTriangleError`), from inside
  `curvature_energy_and_forces` rather than left to a caller. Not switched to the obtuse branch:
  that keeps the energy correct and silently breaks the force, which is the failure with no
  symptom.
- **Degenerate (zero-area) triangle.** No normal, no cotangents. **Refused**, not regularised with
  an epsilon — the reference's `+ 1e-300` would return a unit-ish normal and a large finite force
  in an arbitrary direction.
- **Non-positive reference face area.** **Refused** by `areal_strain`, not floored at `1e-30`.
- **Exponential overflow.** `beta * |eps|` past the card's ceiling is **refused at a stated
  strain**, rather than returning `inf` and then `nan` forces.
- **Reference configuration.** `psi(0) = 0` and `psi'(0) = 0` identically, so the areal energy and
  force are **exactly** `0.0` at the reference areas — testable with `==`, not `approx`.
- **`beta = 0`.** A live branch, not an error. It is the linearised card the negative control uses.

**Invariants, each with the control that asserts it and the number measured.**

| Invariant | Control | Measured |
|---|---|---|
| `F = -grad E`, areal term alone | `test_the_areal_force_matches_the_finite_difference_at_second_order` | rel residual **4.343e-07**, orders **2.000, 2.000** |
| `F = -grad E`, curvature term alone | `test_the_curvature_force_matches_the_finite_difference_at_second_order` | rel residual **2.102e-08**, orders **2.000, 2.000** |
| Areal internal force has no resultant | `test_the_internal_force_has_no_resultant[areal]` | **1.711e-17** of constituent scale |
| Curvature internal force has no resultant | `test_the_internal_force_has_no_resultant[curvature]` | **1.424e-17** |
| Areal internal force has no couple | `test_the_internal_force_has_no_couple[areal]` | **3.963e-17** |
| Curvature internal force has no couple | `test_the_internal_force_has_no_couple[curvature]` | **5.484e-18** |
| Both energies are translation-invariant | `test_translating_the_whole_surface_changes_no_energy` | `rel=1e-12` under a `(3, -7.5, 0.25)` µm shift |
| Bending energy is **exactly** scale-invariant | `test_the_bending_energy_is_exactly_scale_invariant` | worst **3.860e-16** over `lambda` from 1e-3 to 1e3 |
| Areal energy is **not** scale-invariant | `test_the_two_terms_are_not_the_same_term_twice` | `> 10x` under a 10% dilation |
| Gauss–Bonnet on the deformed nucleus | `test_gauss_bonnet_holds_on_the_deformed_nucleus` | residual **-4.263e-14** against `4 pi` |
| Divergence-theorem volume identity | `test_the_divergence_theorem_volume_identity_holds` | rel **1.932e-16** |
| Lifted topology is closed genus-0 and owns its arrays | `test_the_topology_lifted_from_the_manifold_is_the_manifold_s` | `V - E + F == 2`, `owndata` |
| Areal energy is exactly zero at reference | `test_the_reference_configuration_costs_exactly_zero` | `== 0.0` |

## 6. Source evidence class and known retractions

`ffn_sim/ac/nucleus/lamina_analytic.py` self-labels its magnitude gates **INVALID until I0-B2
closes**, in its own module docstring, and marks `eps_rupture` as a **GAP — do NOT tune**. It claims
only its *structural* gates (the `8 pi kappa` topology target, FD-gradient sign, continuity at the
knee, tangent-ratio jump, rupture on/off, Young–Laplace) as currently passing. That is an honest
self-assessment and it is the reason no number from the file was carried across.

Reachability: `ac/nucleus/` is reachable from `ac/nucleus/native_gates/ng0_nucleus_parity.py` and has
live tests in `ffn_sim/tests/ac/nucleus/` (`test_lamina_area.py`, `test_lamina_bending.py`,
`test_chromatin.py`). Those tests exercise the law rather than only the plumbing for the *bending*
term. **The areal term has no force function to test**, so nothing there exercises `F = -grad E` for
it — which is §3's first point, found by reading the module rather than by reading its tests.

Where I looked for retractions: the module docstrings, `ffn_sim/ac/nucleus/INTEGRATION.md`,
`ffn_sim/STATE.md`, `ffn_sim/STATE_NONQUOTABLE.md`, and `ffn_sim/outputs/ac/nucleus/`. The
self-labelled `INVALID`/`GAP` statuses above are what that turned up. No separate retraction record
was found for these files.

Aleph-side evidence class for everything in this entry: **`ANALYTIC_ORACLE`**. Quantitative status
**`BLOCKED`**. Citation status of every card value: **`UNSOURCED`**.

## 7. Independent oracle or derivation

Four, none of which is the module compared against itself or against `ffn_cellsim`.

1. **`8 pi kappa`, the continuum Willmore energy of a sphere.** Radius-independent, so a wrong
   normalisation, a factor of two or an inconsistent winding all move it and no radius absorbs
   them. Measured relative errors at levels 1/2/3 (V = 42/162/642): **7.1655e-02, 1.8705e-02,
   4.7414e-03**, observed order **1.959**. Also checked to `rel=6e-3` at radii 2, 5 and 11 µm.
2. **A closed form for the areal assembly.** Under a uniform dilation every face strains
   identically, so the sum over faces collapses to `A_ref_total * psi(lambda^2 - 1)` — an oracle for
   the whole assembly, scatter included, sharing no line with the implementation. Measured relative
   agreement **6.353e-16, 2.255e-15, 9.060e-16** at dilations 1.02, 1.08, 1.15.
3. **An independently written operator inside Aleph.** `aleph/state/curvature.py` was written by a
   different lane, assembles the same discrete mean-curvature normal from its own cotangent code and
   its own mixed-Voronoi areas, and knows nothing about this energy. Building
   `(kappa/2) sum_v |L x_v|^2 A_v` from its output is an independent implementation of the same
   functional. Measured relative difference: **0.000e+00** — bit-identical.
4. **Exact identities with known exponents.** Scale invariance of the bending energy (degree-zero
   homogeneity, worst **3.860e-16** across six decades of `lambda`) and second-order convergence of
   the central difference (**2.000** on both terms).

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive (primary) | `tests/vertical/test_nucleus_lamina_controls.py::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_areal_force_matches_the_finite_difference_at_second_order` | rel residual < 1e-6 at observed order > 1.6; measured **4.343e-07** at orders **2.000, 2.000** |
| Positive (primary) | `…::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_curvature_force_matches_the_finite_difference_at_second_order` | same; measured **2.102e-08** at orders **2.000, 2.000** |
| Positive (analytic) | `…::TestTheCurvatureEnergyIsTheHelfrichEnergy::test_the_sphere_error_converges_at_second_order_under_refinement` | error falls monotonically and observed order > 1.8; measured **1.959** |
| Positive (analytic) | `…::TestTheArealLawIsTheOneClaimed::test_a_uniform_dilation_reproduces_the_closed_form` | assembled energy matches `A_ref psi(lambda^2-1)` at `rel=1e-12`; measured **≤ 2.255e-15** |
| Positive (cross-module) | `…::TestTheCurvatureEnergyIsTheHelfrichEnergy::test_it_agrees_with_the_independently_written_curvature_operator` | agreement at `rel=1e-13`; measured **0.000e+00** |

**The gradient check is run per term, separately, and never on the sum.** A summed check is exactly
the one a missing or doubled term survives, because two wrong gradients can still add up correctly.
Here it would be worse than usual: the two force scales differ by five orders of magnitude
(**2.816e+04 pN** areal against **9.509e-02 pN** curvature), so a summed check would be blind to the
curvature term entirely. Each term is additionally required to show `energy > 0.0` and a non-zero
force *before* the finite difference is taken, or the comparison is between two zeros and passes for
a term that does nothing at all.

**The finite-difference steps sit above the round-off floor.** 4e-4, 2e-4, 1e-4 µm on a 5 µm
nucleus, all far above `eps^(1/3) ~ 6e-6` relative. `PLAN.md` §6.1 records the afternoon this was got
wrong: below the floor the measured error is cancellation rather than truncation, the apparent order
goes negative, and a **correct** gradient looks broken. The steps also cannot walk the mesh across
the obtuse branch they are measuring on: the minimum cotangent in the test configuration is
**0.058244** and the largest step moves a cotangent by order 1e-4.

## 9. Deliberately failing negative control

**The break is `stiffening_exponent = 0.0` — a parameter of the shipped `LaminaCard` driving the
shipped `beta == 0` branch of the shipped potential.** Not a hand-edited copy of the module. A
negative control run against a duplicate only proves the duplicate is broken.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `…::TestTheLinearisedLaminaIsInvisibleToEveryNumericalCheck::test_the_named_check_refuses_the_linearised_card` | `assert_lamina_stiffens_at_large_strain` raises `LaminaLinearisationError` on the linearised card |
| Negative (must not fire) | `…::test_the_named_check_does_not_fire_on_the_declared_card` | passes on the declared card and returns the ratio; measured **54.5982** against **1.000000** |
| Negative (detector sees the defect) | `…::test_the_defect_is_real_and_is_specifically_a_large_strain_defect` | tension ratio **1.1070** at `eps = 0.01`, **753.7777** at `eps = 0.44` |
| Negative (check cannot be defanged) | `…::test_a_threshold_that_a_linear_law_would_satisfy_is_itself_refused` | `minimum_tangent_ratio` of 1.0, 0.5 or 0.0 is itself refused |
| Negative (threshold, not `beta == 0`) | `…::test_a_weakly_stiffening_card_is_caught_too` | `beta = 0.5` is not `is_linearised` and is still refused (1.105x at `eps = 0.2`) |

**Why the check has to exist as a named check, measured rather than argued.** Linearising the lamina
breaks **nothing** a numerical health indicator can see. `test_the_gradient_control_passes_just_as_well_on_the_linearised_card`
runs the full finite-difference battery on the linearised card and it passes at rel **7.048e-09**
with observed order **2.000** — *better* than the nonlinear card. `test_closure_also_survives_the_linearisation`
shows force closure still holds to 1e-12 of the constituent scale. The energy is still a valid
potential, the force is still exactly its gradient, everything is finite and the run would still
converge. The defect is only visible to a check that knows what the *model* is supposed to claim.

**The detector is required to see the defect before it is asked to refuse it.** `PLAN.md` §6.1
records a reporter that returned `0.0` from an unpopulated cache and made a deliberately broken
strut look innocent. The same pattern is applied to the obtuse guard:
`test_the_guard_sees_the_obtuse_triangles_before_it_refuses` parses the count out of the refusal
message and requires it to be non-zero, and `test_the_configuration_the_controls_run_on_is_clean`
requires the guard **not** to fire on the configuration everything else runs on.

### 9a. The controls were checked for the ability to fail

Six mutants were introduced into the **shipped** module, one at a time. Each run is the full
57-test file.

| Mutant | Tests killed |
|---|---|
| M1 curvature energy coefficient `kappa/8 → kappa/4` | **6** |
| M2 mixed-Voronoi area coefficient `0.125 → 0.25` (area only, gradient coefficients untouched) | **6** |
| M3 `np.expm1(beta*eps)` → `np.exp(beta*eps) - 1.0` in the areal tension | **2** |
| M4 curvature gradient: explicit-position term sign flipped | **1** |
| M5 obtuse guard disabled (`cotangent < 0.0` → `cotangent < -1.0`) | **2** |
| M6 areal force uses the tangent modulus `K(eps)` instead of the tension `sigma(eps)` | **3** |

The source was restored after each mutant and verified: SHA-256
`ea47295d5788aa5fbb8594acd74e96491c217ad14196a20e66c5008ae43914b3`, `cmp` against a pre-mutation
byte copy clean, `git diff -- aleph/vertical/nucleus_lamina.py` empty, `git status --porcelain` on it
empty.

**M3 initially survived, and that is the most useful thing this section records.** On the first run
the `expm1 → exp(x)-1` mutant killed **zero** tests. The test meant to catch it used
`np.allclose(measured, expected, rtol=1e-6)` — and `np.allclose` carries a **default `atol=1e-8`**,
which completely swamps the values of order 2e-10 that the small-strain limit lives at. The test was
vacuous and looked exactly like a passing test. It was replaced with an explicit relative comparison
against the analytic bound `beta * eps` (2e-11 at `eps = 1e-12`, where `exp(x)-1` delivers 8.274e-08
— four thousand times the physical nonlinearity it is supposed to be reporting), and M3 now kills 2.
Every other `np.allclose` in the file was given an explicit `atol=0.0` in the same pass. **A default
absolute tolerance is the quietest way a numerical test stops testing anything**, and this was found
by mutation rather than by reading.

M4 killing only **1** test is honest and thin: the curvature finite-difference test is the only
control that reads the curvature gradient's sign, because closure survives a sign flip on a
symmetric pair term. That is a stated weakness of this control set, not a strength.

## 10. Numerical and precision envelope

Working and accumulation precision: **float64 throughout**, no mixed precision, no accumulation
trick. All arrays are `np.float64`; `EnvelopeTopology` indices are `np.int64`.

Tolerances and why each is that number and not looser:

- **Finite-difference gradient, `rel < 1e-6` at order `> 1.6`.** The order threshold is the load
  bearing half — a correct central difference converges at 2 and a wrong analytic gradient shows a
  low or negative order regardless of how small the residual is. Measured orders are **2.000** on
  both terms, so the threshold has 0.4 of headroom it does not need.
- **`8 pi kappa`, `rel=6e-3` at level 3.** This is a *discretisation* tolerance, not a precision
  one: the measured level-3 error is 4.7414e-03 and it is truncation, falling at order 1.959. It is
  set just above the measured value so that a regression in the discretisation is visible; it is not
  a claim that the discrete energy is accurate to 0.6%.
- **Closed-form dilation, `rel=1e-12`.** Round-off only; measured ≤ 2.255e-15.
- **Cross-module agreement, `rel=1e-13`.** Measured 0.000e+00. Not asserted at `== 0` because
  bit-identity between two independent assemblies is a coincidence of summation order that a
  refactor in either module may legitimately break.
- **Closure, `< 1e-12` of the constituent scale.** Measured 1.4e-17 to 4.0e-17. Never against the
  resultant: `PLAN.md` §9 records a resultant-based tolerance that gets *stricter* the more correct
  the physics is, because correct opposing forces cancel, and it rejected 40,000 consecutive steps
  while reporting a plausible tension.
- **Scale invariance, `rel=1e-14`.** Measured 3.860e-16 worst case.
- **Gauss–Bonnet, `< 1e-12` absolute.** Measured -4.263e-14. This is exact for any closed genus-0
  triangulation at any resolution, so the tolerance is a round-off budget, not a convergence one.

Range over which the stated tolerances hold: closed genus-0 triangulated surfaces with **no obtuse
triangle**, areal strain satisfying `beta*|eps| <= max_exponent_argument` (default 60), reference
face areas finite and positive. Outside any of those, the module **refuses** — see §5. It does not
degrade silently anywhere that was found.

Conditioning notes. `cot_m = (a.b)/|N|` is evaluated without a trigonometric round trip, so it stays
accurate for the near-degenerate corners `arccos` handles worst. `expm1` rather than `exp(x)-1`
carries the small-strain limit, which mutant M3 now pins. The known ill-conditioned direction is
`|N| -> 0`, where both the cotangents and the area gradient blow up, and that is the degenerate-face
refusal.

## 11. Production-backend residency and transfer

**Host / CPU only.** Pure NumPy, float64, no `warp`, no device arrays, no host round trip, no GPU
work of any kind. **Zero GPU jobs were run in the session that produced this entry**, and no GPU
authorization exists or was requested.

This module is not yet resident on the production backend at all: it is not scheduled by
`aleph/runtime/`, not wired into `aleph/vertical/assembly.py`, and defines no `Participant`. If it
later moves to `WarpBackend`, the per-face `(F, 3, 3, 3)` cotangent-gradient tensor is the object to
look at first — it is the only large intermediate, and at level 3 it is 642·… ≈ 1280 faces × 27
float64, which is small, but it scales linearly with face count.

It imports `aleph.state.manifold` and nothing from `validation/**`, so the rule that the oracles
judge the runtime and never the reverse is not touched.

## 12. Comments and docstrings to discard

Nothing was carried, so nothing needed stripping — but the reference prose that must **not**
appear, and was checked for, is:

- gate and phase vocabulary: `I0-B2`, `I2`, `framework-#6`, `KB-3.B2.1/.2/.3/.4`, `Magic-Number
  Block`, `native gate`, `NG0`, `founding rule`;
- status labels: `GAP`, `INVALID until`, `do NOT tune`, `EMERGENT`;
- backend and document references: `Warp-CUDA envelope kernels`, `helfrich_bending_kernel`,
  `params_i0b2.yaml`, `INTEGRATION.md`;
- the provider's own naming: `kappa_tilde`, `Sigma_ref`, `LaminaParams`, `k_soft`, `k_lamin_ac`,
  `knee_strain`, `eps_rupture`, `lamina_tension`, `dihedral_bending_energy`, `wlc_tension`.

What replaces them: prose written for Aleph's situation, stating the derivation in §4's terms, and
naming the two things Aleph does differently and why (an integrable potential rather than a
`sigma(eps)` with no energy; chromatin as a separate `[I]` network rather than folded into the
lamina's soft branch). The module's docstring carries the mixed-Voronoi measurement table by
reference to `aleph/state/curvature.py`, which is Aleph's own measurement, not the reference's.

`tests/ports/test_port_discipline.py::test_no_provider_vocabulary_leaks_into_the_package` is the
check that this stayed true.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** Status `PROPOSED`. 57 controls pass (3.5 s) and six mutants are killed, but the module is unwired, `LaminaEnergy` is dead, and three of the four substructures the lane was commissioned for do not exist. See §14. |
| Reviewer | **Agent-proposed, unratified.** Written by the L30 lane itself; no second reader. `decided_by` is deliberately absent. |
| Rollback | Delete `aleph/vertical/nucleus_lamina.py` and `tests/vertical/test_nucleus_lamina_controls.py`. **Nothing breaks.** No module in `aleph/**` imports either, `aleph/vertical/__init__.py` does not list `nucleus_lamina`, and no connector, assembly or campaign references them. That is a statement about how unfinished this is, not about how safe it is. |

## 14. What is NOT established

- **The wired connector count contributed by this entry is 0.** The lane was commissioned against
  five connectors — `nucleus_cytosol_boundary`, `nucleus_cortex_contact`, `actin_cap_linc`,
  `mt_nucleus_linc`, `if_nucleus_linc`. **None is wired. None is even scaffolded.** There is no
  `EndpointHandle`, no force sink, no site set and no contact quadrature in the tree.
- **Three of the four `[I]` substructures do not exist.** `chromatin` (no network, explicit or
  reduced), `nucleoplasm` (no volume constraint, no viscous drag) and the LINC/contact **site sets
  addressed by material coordinate** were all specified and none was written. `nuclear_envelope`
  exists only as `EnvelopeTopology`, which is connectivity, not a boundary with a no-flux condition.
- **There is no `nucleus` owner.** No `Participant`, no `EntityCensus`, no snapshot/commit/rollback,
  no `owned_state` keys. **The nine `owned_state` entries the registered contract declares —
  `envelope_surface_geometry`, `envelope_faces`, `face_reference_areas`, `lamina_areal_state`,
  `lamina_curvature_state`, `chromatin_network_state`, `enclosed_volume`, `reduced_coordinates`,
  `accepted_nuclear_state` — are matched by nothing in this module.** Two of them
  (`face_reference_areas`, and part of `envelope_faces`) are passed around as bare arguments.
- **The strain-regime crossover was not measured, and could not have been.** The manuscript's §8.2
  P0-C experiment compares a chromatin term against a lamina term. Only the lamina term exists, so
  there is no crossover to find, no location to report, and no claim here about small-strain versus
  large-strain dominance. What §9 measures is that the *lamina's own* tangent modulus stiffens; that
  is one of the two limbs and it is not the experiment.
- **Total force on the isolated nucleus was not measured** in the sense the brief asked for, because
  there is no isolated nucleus. What was measured is per-term closure on the lamina surface alone
  (1.4e-17 to 4.0e-17 of constituent scale), which is a weaker statement.
- **`LaminaEnergy` is exported and dead.** No function in the module returns one. This is a
  declaration-versus-reality gap in our own tree — the shape of defect this project was restarted
  over — and it is pinned by `TestTheEnergyContainerIsDeclaredButUnwired::test_no_public_entry_point_returns_one`
  so that wiring it up breaks a test that names this section.
- **`EnvelopeTopology.generation` is declared and unread.** It exists for the material-coordinate
  site sets that were not written. Nothing increments it and nothing checks it.
- **No card value is evidence of anything.** `areal_modulus_pn_per_um = 200`,
  `stiffening_exponent = 20`, `bending_rigidity_pn_um = 0.0856` are placeholders chosen to make the
  algebra sharp. `UNSOURCED`, `BLOCKED`. No stiffness computed from them may be reported as a
  property of a lamina, and the 1.4x threshold in `assert_lamina_stiffens_at_large_strain` is a
  threshold on the *model*, not a measurement.
- **No lamin isoform, no rupture, no yielding, no viscoelasticity, no in-plane shear.** The areal law
  is a function of a face's *area* and knows nothing about its shape, so a lamina resisting shear at
  constant area is outside this representation. The reference models rupture; Aleph does not, because
  a threshold read off an assumed card is the card restated.
- **No spontaneous curvature.** Deliberate: the `c0` cross term is degree-one homogeneous in the
  positions, so it enters the dilational virial and breaks the scale-invariance identity §7 relies
  on. `PLAN.md` §10 records the membrane lane measuring exactly that at 5.7e-05 against a 3.1e-17
  floor.
- **Never stepped.** The module has never been inside a candidate step, a transaction, or a
  relaxation. Every number here is from a static configuration.
- **The obtuse refusal has an unmeasured cost.** A nucleus under real compression will produce
  obtuse triangles, at which point this module stops rather than degrades. **How much deformation it
  tolerates before that happens was not measured.** A jitter of 0.1 µm on a 5 µm level-2 icosphere
  is enough to turn 3 of 320 triangles obtuse, which suggests the usable envelope is narrow, but
  that single number is one mesh and one perturbation and is not a characterisation.
- **M4 is killed by exactly one test.** The curvature gradient's sign is read by one control and
  nothing else. A second, independent check on that gradient is missing.

## 15. Process deviation, recorded rather than smoothed over

**This entry was written after the code, which PLAN §0.2.5 forbids.** The L30 lane wrote
`nucleus_lamina.py` and then died on an API 529 before writing either the ledger entry or the tests.
The module sat in the repository at 713 lines, parsing and importing, with no controls — the same
state `sf_arc.py` was left in earlier the same day, and the same shape as the L11 failure in
`PLAN.md` §6.1. Nothing consumed it, which is the only reason it cost nothing.

**A second, smaller ordering choice, made deliberately.** The coordinator asked for this entry
before the tests. It was written after them instead, because §5, §8, §9 and §9a name twenty-one
tests by path and `tests/ports/test_port_discipline.py::test_named_controls_resolve_to_real_tests`
requires every one to resolve to a real `def test_*`. Writing the ledger first would have meant
either citing tests that did not exist yet — precisely the L11 defect — or quoting numbers from
memory. `PLAN.md` §6.1 records five API errors found by reading, three of which would have passed
silently. **Every number in this entry was measured in the session that wrote it**, and the mutation
table in §9a is the record of one of them having been wrong.

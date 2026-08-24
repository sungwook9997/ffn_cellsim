# ALEPH-PORT-502 — eigen-spectrum observables of a symmetric operator

| Field | Value |
|---|---|
| Lane | L5 (observation layer) |
| Target | `aleph/observe/spectral_observer.py` |
| Reference consulted | `/Users/sw1/ffn_cellsim/ffn_sim/ac/engine/observe/spectrum.py` (667 lines) — dense path only |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` |
| Source file digest | `sha256:57009047c6d7ab6c5d2e918f16b5526bd62765bf680e842a1e90414036e2d83e` |
| Port class | **RE-DERIVED**, scoped: the dense spectrum is ported, the Lanczos path is deliberately **NOT** |
| Written | 2026-07-30 |
| Status | LANDED |

## 1. What the reference does, and what Aleph takes

The reference module has two halves. The dense half diagonalises an assembled symmetric operator and
derives from the spectrum: a backward-stability zero tolerance, a null-space count checked against a
declared rigid-body count, the condition number as a time-scale separation, and the Gershgorin
row-sum bound compared against the true `lambda_max`. The iterative half is a full reorthogonalised
Lanczos with spectral folding for the small end.

**Aleph takes the dense half only.** The Lanczos path is real work and it is correct as far as this
reading goes, but it exists to answer a question Aleph does not yet have — the reference needed it
because its dense assembly refuses above 8,192 DOF and its native population is ~1.5M DOF. Aleph's
initial scope (`ALEPH-DQ-101`, membrane–cortex–pressure) has no assembled operator at that scale, and
porting an iterative eigensolver *before* there is an operator to run it on would be importing a
solution to an unowned problem. It is recorded here as a deferred candidate, not as an oversight.

## 2. Independent re-derivation

### 2.1 The zero tolerance is derived, not chosen

A symmetric eigensolve (`numpy.linalg.eigh` -> LAPACK `dsyevd`) is backward stable: the computed
eigenvalues are the *exact* eigenvalues of a perturbed matrix `K + E` with

    ||E||_2 <= p(n) * eps * ||K||_2

where `eps` is unit round-off (float64: `2.22e-16`) and `p(n)` is a modest polynomial in the
dimension. Weyl's inequality then bounds each eigenvalue's error by the perturbation's norm:

    |lambda_i(K + E) - lambda_i(K)| <= ||E||_2 <= p(n) * eps * ||K||_2

So a computed eigenvalue whose magnitude is at or below `p(n) * eps * ||K||_2` is **indistinguishable
from zero** — not "small", but genuinely unresolvable in this arithmetic. That is the floor, and it
is a property of float64 and of the operator's own norm. `p(n) = n` is taken as the conventional
conservative surrogate.

For a symmetric matrix `||K||_2 = max_i |lambda_i|`, so the floor is computed from the spectrum
itself and costs nothing extra.

**Nothing here is tuned.** The alternative — picking a tolerance until the mode count comes out at
the expected number — is the failure this construction exists to prevent, and the module makes it
impossible to do accidentally by refusing to widen the floor when the count disagrees. A mismatch is
*reported* as a mismatch.

### 2.2 Why the Gershgorin bound over-estimates

Gershgorin: every eigenvalue of `K` lies in the union of discs centred at `K_ii` with radius
`R_i = sum_{j != i} |K_ij|`. Hence

    |lambda| <= max_i ( |K_ii| + R_i ) = max_i sum_j |K_ij| = ||K||_inf

and for symmetric `K`, `||K||_inf = ||K||_1`, which bounds `||K||_2 = rho(K)` above.

The bound is attained only when some eigenvector concentrates on the maximising row *with the signs
of that row's entries perfectly aligned*, i.e. `v_j = sign(K_ij)`. A generic eigenvector does not do
this — its entries carry the operator's own phase structure and the row sum then involves
cancellation. The bound therefore replaces a signed sum with a sum of absolute values, and the
over-estimate is exactly the cancellation it threw away.

This matters operationally: an explicit integrator sized by `2 / gershgorin_bound` is taking steps a
factor `gershgorin_bound / lambda_max` smaller than the true stability limit `2 / lambda_max`. That
factor is unmeasured wherever the bound is used in place of the truth, and reporting it turns an
inherited conservatism into a number.

**Analytic check available.** For the 1-D Dirichlet Laplacian `tridiag(-1, 2, -1)` of size `n`, the
eigenvalues are known in closed form,

    lambda_k = 2 - 2*cos(k*pi/(n+1)) = 4*sin^2(k*pi/(2*(n+1))),  k = 1..n

so `lambda_max = 4*sin^2(n*pi/(2(n+1)))`, while the Gershgorin bound is exactly `4` (any interior
row: `2 + 1 + 1`). The ratio is `1 / sin^2(n*pi/(2(n+1)))`, which is `> 1` for every finite `n` and
tends to 1 from above. This gives a positive control with an exact expected value at every `n` —
both for `lambda_max` and for the over-estimate factor.

### 2.3 The near-null-space count needs a *structural* expectation

Counting modes below the floor is only a check if the expected count comes from somewhere other than
the measurement. For a free-floating rigid body in 3-D the count is **6**: three translations and
three rotations. This is a statement about the physical situation, declared before the eigensolve,
not read off the histogram.

Aleph's positive control makes the 6 exact rather than approximate. Take a spring network with bonds
`(i, j)` of unit direction `n` and stiffness `k`, at its **unstressed** configuration. The
first-order tangent is the affine one,

    U = (1/2) * sum_bonds k * [ n . (u_i - u_j) ]^2

which contributes `k * n n^T` blocks to `K_ii`, `K_jj` and `-k * n n^T` to `K_ij`, `K_ji`. Then:

- **Translations.** `u_i = t` for all `i` gives `u_i - u_j = 0`, so `U = 0` exactly. Three modes.
- **Rotations.** An infinitesimal rotation `u_i = w x r_i` gives
  `u_i - u_j = w x (r_i - r_j) = L * (w x n)`, and `n . (w x n) = 0` identically. So `U = 0`
  **exactly**, not to first order in a small parameter. Three more modes.

Six exact zero modes, for any bond topology, at the unstressed configuration. Whether there are
*more* than six is a property of the built network — an under-braced network has internal floppy
modes, and those are a measurement, reported as `n_excess_zero_modes`, never absorbed into the
tolerance.

Note the asymmetry, which the module encodes: the declared count is a **lower bound**, so
`measured < declared` is a hard defect (the operator is constraining a motion nothing physically
constrains, i.e. the assembly is wrong), whereas `measured > declared` is a finding about the model.

### 2.4 The condition number is a time-scale separation

For overdamped relaxation `gamma * xdot = -K x`, mode `i` decays at rate `lambda_i / gamma`. The
ratio of the fastest to the slowest resolvable rate is `lambda_max / lambda_min_nonzero`, so
`log10(kappa)` is the number of decades of dynamics an explicit integrator must resolve to follow the
slowest mode to completion. That is the case for an implicit solve, stated as a measurement rather
than as an opinion.

`lambda_min_nonzero` is the smallest eigenvalue *above the floor*. Below the floor there is no
eigenvalue, only round-off, and a condition number formed from round-off is a property of the
arithmetic, so it is returned as `None` with a stated reason rather than as a large number.

### 2.5 Symmetrisation is reported, never silent

`eigh` requires symmetry. The module symmetrises as `(K + K^T)/2` and **reports the discarded
antisymmetric part** as `asymmetry = ||K - K^T||_F / ||K + K^T||_F`. Silently symmetrising an
operator that is not symmetric hides exactly the defect that a symmetry check exists to find — for a
tangent `K = grad^2 U`, asymmetry means the force field is not a gradient, which is a modelling
error and not a numerical one. Above a declared bound the observer **refuses**.

## 3. Controls shipped

All under `tests/observe/test_spectral_observer.py`.

| Control | Test | Asserts |
|---|---|---|
| Positive | `test_laplacian_eigenvalues_match_the_closed_form` | `4 sin^2(k pi / 2(n+1))` to 1e-12, at four sizes |
| Positive | `test_gershgorin_bound_is_exactly_four_for_the_laplacian` | the analytic bound of §2.2 |
| Positive | `test_gershgorin_over_spectral_radius_matches_the_closed_form` | over-estimate factor `1/sin^2(...)` exactly |
| Positive | `test_gershgorin_over_spectral_radius_tends_to_one_with_size` | monotone, converging from above |
| Positive | `test_gershgorin_bound_is_never_below_lambda_max` | 60 random symmetric operators |
| Positive | `test_diagonal_operator_saturates_the_gershgorin_bound` | the bound is tight exactly where §2.2 says |
| Positive | `test_unstressed_spring_network_has_exactly_six_zero_modes` | the §2.3 derivation, measured |
| Positive | `test_rigid_body_modes_are_zero_to_round_off` | the six modes evaluated in the quadratic form directly, not trusting the eigensolver |
| Positive | `test_rigid_body_modes_3d_is_six_per_free_component` | the declaration helper |
| Positive | `test_an_underbraced_network_reports_its_floppy_modes_as_excess` | excess is a measurement, not an error |
| Positive | `test_zero_tolerance_equals_n_eps_norm` | the floor is the formula, with no fudge factor |
| Positive | `test_zero_tolerance_scales_with_the_operator_norm` | backward stability is relative |
| Positive | `test_condition_number_is_the_ratio_of_resolvable_extremes` | on a diagonal operator with known spectrum |
| Positive | `test_stable_explicit_step_is_two_over_lambda_max` | the true stability limit |
| Positive | `test_negative_modes_are_counted_not_clipped` | saddle directions survive reporting |
| Positive | `test_round_off_asymmetry_is_tolerated_and_reported` | reported, not hidden, below the bound |
| **Negative** | `test_pinning_a_node_destroys_rigid_modes_and_is_reported_a_defect` | `measured < declared` -> defect, not absorbed |
| **Negative** | `test_declared_count_is_not_widened_to_fit` | the tolerance tracks the norm and nothing else, including a declared count of 999 |
| **Negative** | `test_regularized_operator_reports_a_fabricated_empty_null_space` | `aI + K` shifts every eigenvalue; the measurement protocol made executable |
| **Negative** | `test_asymmetric_operator_is_refused` | typed refusal, not a silent symmetrisation |
| **Negative** | `test_non_finite_operator_is_refused` | no NaN spectrum |
| **Negative** | `test_non_square_operator_is_refused` | typed refusal |
| **Negative** | `test_empty_operator_is_refused` | no spectrum of nothing |
| **Negative** | `test_a_pure_null_operator_declines_a_condition_number` | round-off is not an eigenvalue |
| **Negative** | `test_an_undeclared_expectation_is_recorded_as_a_skipped_check` | an undeclared expectation cannot be violated, and the record says so |
| **Negative** | `test_refusal_is_returned_not_raised` | every failure path returns |

The `test_regularized_operator_...` control deserves a note: `a*I + K` has spectrum
`{a + lambda_i}`, so every zero mode moves to `a`, and if `a` exceeds the floor the null space is
measured as **empty**. A spectrum measured on a regularized probe therefore fabricates a rigid-body
defect out of nothing. The test asserts that this happens — and checks that every eigenvalue moved
by exactly the shift, which is the mechanism — so that "diagonalise at zero regularization" is an
executable statement rather than a comment.

## 3a. Units, domains, singular cases, invariants

**Units.** The module is unit-agnostic and says so: eigenvalues come back in whatever units the
caller's operator carried. For an Aleph stiffness that is pN/µm, and then `stable_explicit_step =
2/lambda_max` is the mobility step in µm/pN — directly comparable to the `dt/gamma` an explicit
integrator uses, which is what makes the Gershgorin comparison a measurement rather than a curiosity.
`zero_tolerance`, `operator_norm` and `gershgorin_bound` share the operator's units;
`condition_number`, `timescale_separation_decades`, `gershgorin_over_spectral_radius` and
`symmetry_defect` are dimensionless.

**Domain.** A square, finite, symmetric `(n, n)` operator with `n >= 1`, assembled at **zero
regularization**. The last condition cannot be checked from the matrix — `a*I + K` is a perfectly
valid symmetric operator — so it is a declared measurement protocol backed by a negative control
that demonstrates what happens when it is violated.

**Singular cases.**

- No eigenvalue above the round-off floor: `lambda_min_nonzero`, `condition_number`,
  `timescale_separation_decades` and `stable_explicit_step` are all `None` with a stated reason.
  A condition number formed from round-off is a property of the arithmetic, not of the operator.
- Zero operator: `symmetry_defect` returns 0.0 rather than evaluating 0/0.
- No declared expectation: `meets_declared_null_space` is `None` and `null_space_defect` is `False`,
  because an undeclared expectation cannot be violated — recorded in `notes` as a *skipped check*,
  which is itself the argument for declaring one.
- Non-square, empty, non-finite, or insufficiently symmetric: typed `Refusal`, never a raise.

**Invariants.**

- `gershgorin_bound >= operator_norm` always — checked over 60 random symmetric operators, and
  exactly attained on a diagonal one, where there is no cancellation to discard.
- `zero_tolerance == n * eps * ||K||_2` exactly, and depends on nothing else. Asserted to be
  unchanged when the declared mode count is varied from 6 to 999.
- Eigenvalues are returned ascending and **in full**; the JSON view truncates, the object does not.
- `n_zero_modes + n_negative_modes <= n_dof`, and negative modes are counted, never clipped.
- `a*I + K` shifts every eigenvalue by exactly `a` — asserted directly, since that identity is the
  mechanism behind the regularization negative control.

## 3b. Numerical and precision envelope, and production-backend residency

**Numerical.** float64 throughout; `eps = 2.220446e-16` is read from `numpy.finfo` rather than
written as a literal. The diagonalisation is `numpy.linalg.eigh` (LAPACK `dsyevd`), whose backward
stability *is* the precision envelope: computed eigenvalues are exact for `K + E` with
`||E||_2 <= p(n) eps ||K||_2`, and `p(n) = n` is the conventional conservative surrogate. That bound
is not a tolerance the module chose — it is the resolution limit of the arithmetic, and it is
reported on every result so a reader can re-derive it. The closed-form Laplacian controls are
asserted to `atol=1e-12`, which is round-off at these norms.

The symmetry bound `max_symmetry_defect = 1e-10` is the one number here with a judgement in it: it
is set at round-off scale on the reasoning that an operator assembled from a genuine potential is
symmetric to accumulation error, so anything larger is structure rather than noise. It is a
parameter, it is documented as a parameter, and exceeding it produces a refusal rather than a
silent symmetrisation.

**Residency.** Host-side, CPU only, dense `numpy`. This bounds the module's applicability and the
bound is stated: a dense `(n, n)` float64 operator is `8 n^2` bytes and `eigh` is `O(n^3)`, so this
path is a *slice* diagnostic, not a whole-system one. The matrix-free iterative path that would lift
that limit was deliberately **not** ported (§1) because Aleph has no operator at a scale that needs
it yet. When one exists, the transfer boundary is explicit: the runtime assembles or applies on
device, and a dense operator is copied to host only when it is small enough to be worth
diagonalising — a decision the caller makes, not this module.

## 4. Compliance with PLAN §0.2

1. All reference comments and docstrings discarded. ✔
2. Every result re-derived: backward stability + Weyl (§2.1), Gershgorin as `||K||_inf` with the
   attainment condition (§2.2), the exact six-mode argument for unstressed springs (§2.3). The
   reference asserts the six-mode count as a structural fact; Aleph derives *why* it is exact and
   makes it a test. ✔
3. New prose. ✔
4. Positive and failing negative controls (§3). ✔
5. Entry written before the code. ✔
6. No `ffn_cellsim` import; no `validation` import. ✔

## 5. Verdict on the reference

**Sound, and scoped down.** The dense path's derivations check out — the backward-stability floor is
correctly stated and correctly attributed, and the refusal to widen it is the right instinct. Aleph
adds the derivation of *why* the rigid-body count is exactly six for an unstressed spring network
(the reference asserts it), converts the "measure at zero regularization" protocol note into a
negative control, and adds a typed refusal for asymmetry where the reference relies on a separate
module being called by a disciplined caller. The Lanczos path is deferred, not rejected.

## 6. Amendment 2026-07-30 16:20 — definiteness gates the overdamped quantities

Found by the measurement-method audit of the whole tree, not by the reference. §2.4's condition
number, its `timescale_separation_decades`, and `stable_explicit_step` all read the operator as a
stiffness at a *minimum* of the potential. `n_negative_modes` was measured and gated nothing, so on
`diag(-100, -50, 1, 2)` the entry's own `2/lambda_max` returned 1.0 — 50x past even the magnitude
limit, and no step is stable against a growing mode at any size — and 0.30 "decades of overdamped
relaxation" for an operator with two growing modes. All three now decline when any mode is negative,
and `stable_explicit_step` also declines when `lambda_max` is itself below the round-off floor
(`diag(-1, 1e-20)` reported "1 zero mode, no eigenvalue above the floor" and "step 2e+20" in one
record).

`gershgorin_over_lambda_max` is renamed `gershgorin_over_spectral_radius` and divides by
`operator_norm`. A row-sum bound bounds `max_i |lambda_i|`; dividing by the *algebraic* maximum
reported 50x conservatism on the same operator, whose bound is exactly tight — against §2.2's claim
that the ratio "is a number rather than an opinion". Identical for positive-semidefinite operators,
which is the case §2.2 is about, and every closed-form Laplacian control passes unchanged under the
new name. The two positive controls in §3 are renamed with it.

New controls: `TestIndefiniteOperatorDeclinesOverdampedQuantities` (four cases, including a definite
operator that must still report all three, so the gate is definiteness and nothing wider) and
`test_an_eigenvalue_below_the_floor_cannot_also_size_a_step`.

A separate defect in the same module, outside this entry's scope: `symmetry_defect` returned 0.0 for
a purely antisymmetric operator, because `||K + K^T||_F` vanishes for that and for the zero operator
alike. See commit `58fb80a`.

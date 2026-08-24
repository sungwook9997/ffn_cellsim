"""Analytic-truth tests for the probabilistic-node oracle.

Every assertion here is scored against a closed form derived INDEPENDENTLY in the test — the
Gaussian moment of the polynomial written out by hand, the MGF, the ``sigma^-2`` breakdown formula —
and never against a larger simulation of the same code.  Where truth is genuinely unavailable in one
line (the exponential's breakdown width) the test re-derives it by bisecting a scalar function
written here, which shares no code with the module.
"""

from __future__ import annotations

import math
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from aleph.virtual_cell.contracts import (
    EvidenceSource,
    RepresentationKind,
    SandboxExperimentCard,
    UncertaintySource,
)
from aleph.virtual_cell.stochastic_node import (
    ENGINE_NODE_REPRESENTATION_TODAY,
    EVIDENCE_SOURCE,
    AnalyticForceLaw,
    ClosureBreakdownRefusal,
    ClosureRefusedError,
    ClosureVerdictKind,
    LinearDrift,
    MomentClosureContract,
    NodeRepresentation,
    ScalarLawDrift,
    SourceFactoredCovariance,
    boltzmann_density,
    compare_closure_against_ensemble,
    cubic_anharmonic_spring,
    exponential_repulsion,
    find_closure_breakdown_sigma,
    gap_cancellation_floor,
    gauss_hermite_nodes,
    leading_gap_order,
    linear_spring,
    make_gaussian_state,
    measure_bimodality,
    measure_cubic_placement_contrast,
    measure_gap_sigma_scaling,
    measure_jensen_gap,
    measure_quenched_thermal_split,
    measure_rate_placement_contrast,
    morse_bond,
    propagate_gaussian_ou_exact,
    propagate_local_gaussian,
    quadrature_moment_update,
    quartic_force,
    stochastic_node_experiment_card,
)

# --------------------------------------------------------------------------------------------
# declared contract — every threshold fixed here, before any number below is read
# --------------------------------------------------------------------------------------------

CONTRACT = MomentClosureContract(
    quadrature_nodes=15,
    gap_relative_tolerance=1e-3,
    expansion_relative_error_tolerance=1e-2,
    derivative_zero_tolerance=1e-12,
    max_admissible_modes=1,
    mode_prominence_ratio=0.05,
    bimodality_coefficient_max=5.0 / 9.0,
    modal_density_ratio_min=0.05,
    mass_normalization_tolerance=1e-9,
    mass_negativity_tolerance=0.0,
    covariance_symmetry_tolerance=1e-12,
    minimum_eigenvalue_tolerance=1e-12,
)

MU = 0.8
SIGMA = 0.15
STIFFNESS = 2.5
CUBIC = 1.7
AMPLITUDE = 40.0
DECAY = 0.4


def _linear() -> AnalyticForceLaw:
    return linear_spring(stiffness=STIFFNESS, rest_position=0.3, unit="pN")


def _cubic() -> AnalyticForceLaw:
    return cubic_anharmonic_spring(stiffness=STIFFNESS, cubic_coefficient=CUBIC, unit="pN")


def _quartic() -> AnalyticForceLaw:
    return quartic_force(coefficient=1.0, unit="pN")


def _exponential() -> AnalyticForceLaw:
    return exponential_repulsion(amplitude=AMPLITUDE, decay_length=DECAY, unit="pN")


def _morse() -> AnalyticForceLaw:
    return morse_bond(well_depth=12.0, width=2.0, equilibrium_position=MU, unit="pN")


# --------------------------------------------------------------------------------------------
# 0. this lane never touches the GPU
# --------------------------------------------------------------------------------------------


def test_module_import_pulls_in_no_gpu_or_learning_framework() -> None:
    """Run in a SUBPROCESS.

    Asserting against this process's ``sys.modules`` would measure whether anything in the whole
    pytest session ever imported Warp, which is a different and false claim: the suite contains
    tests that legitimately import it.  A fresh interpreter measures what the name says.
    """
    script = textwrap.dedent(
        """
        import sys
        import numpy as np
        from aleph.virtual_cell.stochastic_node import (
            MomentClosureContract,
            cubic_anharmonic_spring,
            measure_jensen_gap,
        )

        contract = MomentClosureContract(
            quadrature_nodes=9,
            gap_relative_tolerance=1e-3,
            expansion_relative_error_tolerance=1e-2,
            derivative_zero_tolerance=1e-12,
            max_admissible_modes=1,
            mode_prominence_ratio=0.05,
            bimodality_coefficient_max=0.5,
            modal_density_ratio_min=0.05,
            mass_normalization_tolerance=1e-9,
            mass_negativity_tolerance=0.0,
            covariance_symmetry_tolerance=1e-12,
            minimum_eigenvalue_tolerance=1e-12,
        )
        law = cubic_anharmonic_spring(stiffness=1.0, cubic_coefficient=1.0, unit="pN")
        measure_jensen_gap(law, mean=0.5, sigma=0.1, contract=contract)

        for blocked in ("warp", "torch", "jax"):
            assert blocked not in sys.modules, f"stochastic_node import pulled in {blocked}"
        print("clean")
        """
    )
    done = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        cwd=str(Path(__file__).resolve().parents[4]),
    )
    assert done.returncode == 0, done.stderr
    assert "clean" in done.stdout


def test_everything_is_oracle_only_and_says_the_engine_carries_points() -> None:
    assert EVIDENCE_SOURCE is EvidenceSource.ANALYTIC_ORACLE
    measurement = measure_jensen_gap(_cubic(), mean=MU, sigma=SIGMA, contract=CONTRACT)
    payload = measurement.as_dict()
    assert payload["evidence_source"] == "analytic-oracle"
    assert payload["evidence_authority"] == "unverified-observer"
    assert payload["engine_node_representation_today"] == ENGINE_NODE_REPRESENTATION_TODAY
    assert "point" in ENGINE_NODE_REPRESENTATION_TODAY


# --------------------------------------------------------------------------------------------
# 1. the gap, in closed form and then measured
# --------------------------------------------------------------------------------------------


def test_affine_law_gap_is_exactly_zero_not_merely_small() -> None:
    """Jensen forbids a gap for an affine law.  Exactly ``0.0``, at every width, is the test.

    If the closed form returned ``1e-18`` the machinery would be inventing a gap where there is
    none, and every other number in this file would be an artefact of the same machinery.  The
    quadrature ESTIMATOR is separately allowed its own round-off — measured at 6.7e-16 relative at
    ``sigma = 0.01`` and 1.4e-15 at ``sigma = 25`` (the wider node spread costs digits) — and the two
    claims are kept apart on purpose.
    """
    law = _linear()
    for sigma in (0.0, 0.01, 0.15, 1.0, 25.0):
        assert law.exact_gap(MU, sigma) == 0.0
        measurement = measure_jensen_gap(law, mean=MU, sigma=sigma, contract=CONTRACT)
        assert measurement.exact_gap == 0.0
        assert measurement.second_order_gap == 0.0
        assert abs(measurement.closure_report.gap) < 1e-14
        assert measurement.closure_report.relative_gap < 1e-14


def test_zero_width_gap_vanishes_for_every_law() -> None:
    """A degenerate distribution is a point, so the gap must vanish for every law.

    Four of the five return exactly ``0.0``; the quartic returns one ulp of its own force magnitude
    because the closed form multiplies out the moment recursion while ``F(mu)`` uses ``pow``.  The
    exact-zero claim is therefore made only for the affine control, where it holds at any width.
    """
    for law in (_linear(), _cubic(), _quartic(), _exponential(), _morse()):
        gap = law.exact_gap(MU, 0.0)
        scale = max(abs(law.force_at_mean(MU)), 1.0)
        assert abs(gap) <= 8.0 * np.finfo(np.float64).eps * scale


def test_cubic_gap_matches_the_hand_derived_gaussian_moment() -> None:
    """``E[x^3] = mu^3 + 3 mu sigma^2`` gives ``gap = -3 a mu sigma^2``, written out here by hand."""
    law = _cubic()
    expected_gap = -3.0 * CUBIC * MU * SIGMA**2
    expected_expectation = -(STIFFNESS * MU + CUBIC * (MU**3 + 3.0 * MU * SIGMA**2))
    measurement = measure_jensen_gap(law, mean=MU, sigma=SIGMA, contract=CONTRACT)
    assert measurement.force_at_mean == pytest.approx(-(STIFFNESS * MU + CUBIC * MU**3), rel=1e-15)
    assert measurement.exact_expectation == pytest.approx(expected_expectation, rel=1e-15)
    assert measurement.exact_gap == pytest.approx(expected_gap, rel=1e-13)
    # for a cubic the second-order term IS the whole gap: F'''' = 0 and E[(X-mu)^3] = 0.
    assert measurement.second_order_gap == pytest.approx(expected_gap, rel=1e-13)
    assert measurement.second_order_relative_error < 1e-13


def test_quartic_gap_matches_the_hand_derived_gaussian_moment() -> None:
    """``E[x^4] = mu^4 + 6 mu^2 sigma^2 + 3 sigma^4``; the ``3 sigma^4`` term is what the
    second-order expansion misses."""
    law = _quartic()
    expected_gap = 6.0 * MU**2 * SIGMA**2 + 3.0 * SIGMA**4
    measurement = measure_jensen_gap(law, mean=MU, sigma=SIGMA, contract=CONTRACT)
    assert measurement.exact_gap == pytest.approx(expected_gap, rel=1e-14)
    assert measurement.second_order_gap == pytest.approx(6.0 * MU**2 * SIGMA**2, rel=1e-14)
    expected_relative = 3.0 * SIGMA**4 / expected_gap
    assert measurement.second_order_relative_error == pytest.approx(expected_relative, rel=1e-12)


def test_exponential_gap_matches_the_moment_generating_function() -> None:
    """``E[A e^{bX}] = A e^{b mu + b^2 sigma^2 / 2}`` — the MGF, written out here."""
    law = _exponential()
    rate = -1.0 / DECAY
    expected = AMPLITUDE * math.exp(rate * MU + 0.5 * rate**2 * SIGMA**2)
    measurement = measure_jensen_gap(law, mean=MU, sigma=SIGMA, contract=CONTRACT)
    assert measurement.exact_expectation == pytest.approx(expected, rel=1e-14)
    assert measurement.exact_gap == pytest.approx(
        AMPLITUDE * math.exp(rate * MU) * (math.exp(0.5 * rate**2 * SIGMA**2) - 1.0), rel=1e-13
    )


def test_quadrature_recovers_the_closed_form_expectation() -> None:
    """The estimator is scored against truth, not against a bigger ensemble."""
    for law in (_linear(), _cubic(), _quartic(), _exponential(), _morse()):
        measurement = measure_jensen_gap(law, mean=MU, sigma=SIGMA, contract=CONTRACT)
        assert measurement.representation is NodeRepresentation.QUADRATURE_ENSEMBLE
        assert measurement.quadrature_relative_error < 1e-13, law.name


def test_gauss_hermite_weights_are_a_probability_measure() -> None:
    positions, weights = gauss_hermite_nodes(MU, SIGMA, nodes=CONTRACT.quadrature_nodes)
    assert positions.size == CONTRACT.quadrature_nodes
    assert np.all(weights > 0.0)
    assert float(np.sum(weights)) == pytest.approx(1.0, abs=1e-15)
    assert float(np.dot(weights, positions)) == pytest.approx(MU, rel=1e-14)
    centred = positions - MU
    assert float(np.dot(weights, centred**2)) == pytest.approx(SIGMA**2, rel=1e-13)


def test_morse_bond_at_its_minimum_has_zero_point_force_and_a_positive_mean_force() -> None:
    """The physical content of direction 1, in one number.

    At the equilibrium bond length a point node feels no force at all.  A distributed node samples
    the hard core more than the soft tail and feels a strictly REPULSIVE mean force that grows
    superexponentially with the width.  No amount of care in evaluating ``F`` at the mean recovers
    it, because the quantity is not a property of the mean.
    """
    law = _morse()
    assert abs(law.force_at_mean(MU)) < 1e-12
    previous = 0.0
    for sigma in (0.05, 0.1, 0.2):
        expectation = law.exact_expectation(MU, sigma)
        assert expectation > previous > -1.0
        previous = expectation
    assert law.exact_expectation(MU, 0.15) == pytest.approx(7.257096154, rel=1e-9)


# --------------------------------------------------------------------------------------------
# 2. how the gap scales, and where the expansion breaks down
# --------------------------------------------------------------------------------------------


SIGMA_LADDER = tuple(0.001 * 1.6**k for k in range(9))


def test_leading_order_prediction_is_computed_from_the_law_not_from_the_fit() -> None:
    assert leading_gap_order(_cubic(), mean=MU, max_order=8, zero_tolerance=1e-12) == 2
    assert leading_gap_order(_quartic(), mean=0.0, max_order=8, zero_tolerance=1e-12) == 4
    assert leading_gap_order(_quartic(), mean=MU, max_order=8, zero_tolerance=1e-12) == 2
    assert leading_gap_order(_linear(), mean=MU, max_order=8, zero_tolerance=1e-12) is None


def test_gap_scales_as_sigma_squared_where_the_second_derivative_survives() -> None:
    """Measured fitted exponents on the declared ladder ``sigma = 0.001 * 1.6^k``, ``k = 0..8``:

    cubic 2.0000000000, exponential 2.0005638120, Morse 2.0018043059.  The fit is what says the
    measured quantity is the leading moment correction and not a sampling artefact — and the
    departure of the last two from exactly 2 is real curvature (higher even orders entering as the
    ladder widens), not noise, since nothing here is sampled.
    """
    for law, expected in (
        (_cubic(), 2.0000000000),
        (_exponential(), 2.0005638120),
        (_morse(), 2.0018043059),
    ):
        report = measure_gap_sigma_scaling(
            law, mean=MU, sigmas=SIGMA_LADDER, contract=CONTRACT, max_order=8
        )
        assert report.predicted_exponent == 2
        assert report.fitted_exponent == pytest.approx(expected, abs=1e-9), law.name
        assert report.r_squared > 0.9999


def test_fit_reports_exponent_four_when_the_second_derivative_vanishes() -> None:
    """The discriminating control: a quartic at the origin has ``F''(0) = 0`` and slope 4.

    Without this the ``~2`` above would only show that the fitter can produce a 2.
    """
    report = measure_gap_sigma_scaling(
        _quartic(), mean=0.0, sigmas=SIGMA_LADDER, contract=CONTRACT, max_order=8
    )
    assert report.predicted_exponent == 4
    assert report.fitted_exponent == pytest.approx(4.0, abs=1e-9)
    assert report.max_log_residual < 1e-12


def test_scaling_fit_refuses_a_law_with_no_gap_to_fit() -> None:
    with pytest.raises(ValueError, match="identically zero gap"):
        measure_gap_sigma_scaling(
            _linear(), mean=MU, sigmas=SIGMA_LADDER, contract=CONTRACT, max_order=8
        )


def _exponential_breakdown_closed_form(decay: float, tolerance: float) -> float:
    """Independent re-derivation: solve ``(e^u - 1 - u)/(e^u - 1) = tol``, ``u = sigma^2/(2 lam^2)``.

    Shares no code with the module; the module bisects the relative error of its own closed forms,
    this bisects the scalar function written out here from the MGF by hand.
    """
    low, high = 1e-14, 60.0
    for _ in range(400):
        mid = 0.5 * (low + high)
        if (math.exp(mid) - 1.0 - mid) / (math.exp(mid) - 1.0) > tolerance:
            high = mid
        else:
            low = mid
    return decay * math.sqrt(2.0 * 0.5 * (low + high))


@pytest.mark.parametrize(
    ("tolerance", "expected_sigma", "expected_in_decay_lengths"),
    [
        (0.01, 0.0801341164, 0.200335),
        (0.05, 0.1804211455, 0.451053),
        (0.10, 0.2574623874, 0.643656),
    ],
)
def test_exponential_closure_breakdown_width_matches_the_closed_form(
    tolerance: float, expected_sigma: float, expected_in_decay_lengths: float
) -> None:
    """**The quotable number.**  For ``A exp(-x / lambda)`` a Gaussian closure reproduces the exact
    expectation to 1 % only while ``sigma <= 0.2003 lambda``; to 5 % while ``sigma <= 0.4511
    lambda``; to 10 % while ``sigma <= 0.6437 lambda``.  Above those widths the second-order term is
    no longer the correction, it is a fraction of it.
    """
    report = find_closure_breakdown_sigma(
        _exponential(),
        mean=MU,
        tolerance=tolerance,
        search_low=1e-3,
        search_high=5.0,
        bisection_iterations=200,
        exactness_tolerance=1e-12,
    )
    assert report.breakdown_sigma is not None
    assert report.breakdown_sigma == pytest.approx(expected_sigma, abs=1e-9)
    assert report.breakdown_sigma == pytest.approx(
        _exponential_breakdown_closed_form(DECAY, tolerance), abs=1e-12
    )
    assert report.breakdown_sigma / DECAY == pytest.approx(expected_in_decay_lengths, abs=1e-6)
    assert report.relative_error_at_breakdown == pytest.approx(tolerance, rel=1e-8)


@pytest.mark.parametrize("tolerance", [0.01, 0.05, 0.10])
def test_quartic_closure_breakdown_width_matches_its_one_line_closed_form(
    tolerance: float,
) -> None:
    """For ``c x^4`` the truncated error is ``sigma^2 / (2 mu^2 + sigma^2)``, so the breakdown width
    is ``mu sqrt(2 t / (1 - t))`` — 0.1421 mu at 1 %, 0.3244 mu at 5 %, 0.4714 mu at 10 %."""
    mean = 1.0
    report = find_closure_breakdown_sigma(
        _quartic(),
        mean=mean,
        tolerance=tolerance,
        search_low=1e-3,
        search_high=5.0,
        bisection_iterations=200,
        exactness_tolerance=1e-12,
    )
    expected = mean * math.sqrt(2.0 * tolerance / (1.0 - tolerance))
    assert report.breakdown_sigma == pytest.approx(expected, abs=1e-12)


def test_cubic_expansion_never_breaks_down_because_it_is_not_an_approximation() -> None:
    """The honest negative result: for this law there is no breakdown width to quote."""
    report = find_closure_breakdown_sigma(
        _cubic(),
        mean=MU,
        tolerance=0.01,
        search_low=1e-3,
        search_high=50.0,
        bisection_iterations=100,
        exactness_tolerance=1e-12,
    )
    assert report.expansion_is_exact
    assert report.breakdown_sigma is None
    assert report.relative_error_at_high == 0.0


def test_quartic_at_the_origin_is_unsafe_at_every_width() -> None:
    """The adversarial control: with ``F''(0) = 0`` the second-order expansion carries none of the
    gap, so its relative error is exactly 1 at every width and the search refuses the whole range."""
    report = find_closure_breakdown_sigma(
        _quartic(),
        mean=0.0,
        tolerance=0.01,
        search_low=1e-3,
        search_high=5.0,
        bisection_iterations=100,
        exactness_tolerance=1e-12,
    )
    assert report.breakdown_sigma is None
    assert report.relative_error_at_low == pytest.approx(1.0, abs=1e-15)
    assert "re-declared" in (report.refusal_reason or "")


def test_breakdown_search_refuses_an_inverted_bracket() -> None:
    with pytest.raises(ValueError, match="search_high must exceed search_low"):
        find_closure_breakdown_sigma(
            _exponential(),
            mean=MU,
            tolerance=0.01,
            search_low=1.0,
            search_high=0.5,
            bisection_iterations=10,
            exactness_tolerance=1e-12,
        )


def test_cancellation_floor_explains_the_spurious_truncation_error_at_tiny_width() -> None:
    """Measured on the cubic control, whose expansion is an algebraic identity: the apparent
    relative "truncation error" is 7.1e-05 at ``sigma = 1e-06`` and 4.6e-13 at ``sigma = 1e-02``.

    All of it is cancellation in ``E[F] - F(mu)``.  The estimator says so rather than letting a
    breakdown search below that width return a number.
    """
    law = _cubic()
    fine = gap_cancellation_floor(law, mean=MU, sigma=1e-6)
    coarse = gap_cancellation_floor(law, mean=MU, sigma=1e-2)
    assert fine > 1e-5
    assert coarse < 1e-10
    # the floor grows as sigma^-2, matching the second-order gap in the denominator
    assert fine / coarse == pytest.approx(1e8, rel=1e-3)
    assert gap_cancellation_floor(_linear(), mean=MU, sigma=SIGMA) == math.inf


def _offset_exponential() -> AnalyticForceLaw:
    """The external reviewer's counterexample law: a constant ``1e6`` plus ``exp(-x)``.

    A constant has no second derivative and no Jensen gap, so it changes neither the exact gap nor
    the second-order truncation — the expansion error is bit-for-bit the pure exponential's.  What it
    does change is the CONDITIONING of the gap: ``E[F]`` and ``F(mu)`` are both about ``1e6`` while
    their difference is about ``1e-9`` at the lower bracket, so the subtraction that defines the gap
    loses every significant digit.  That is why it is the right adversarial law for the floor gate.
    """
    return AnalyticForceLaw(
        name="offset-exponential",
        polynomial_coefficients=(1e6,),
        exponential_terms=((1.0, -1.0),),
        unit="pN",
        description="constant offset plus a decaying exponential; ill-conditioned Jensen gap",
    )


def test_breakdown_search_refuses_a_bracket_whose_low_end_is_below_its_own_floor() -> None:
    """**Codex external review 2026-07-29, finding 5.**  The floor was computed and then ignored.

    Reproduced exactly as reported: ``1e6 + exp(-x)`` at ``mu = 0``, tolerance ``0.01``, bracket
    ``[1e-4, 2.0]``.  Before the gate this returned ``breakdown_sigma = 0.2003352947748983`` while
    its own ``cancellation_floor_at_low`` was ``0.08871397861214594`` — nine times the tolerance the
    search was scoring against — and ``relative_error_at_low = 0.0011703962790697465``, itself far
    BELOW that floor.  So the number that certified the low end of the bracket as "within 1 %" was
    round-off in a difference of two numbers near ``1e6``, and the search read it as closure
    accuracy.  A width issued from that bracket would be a safety range with no measurement behind
    it.
    """
    law = _offset_exponential()
    report = find_closure_breakdown_sigma(
        law,
        mean=0.0,
        tolerance=0.01,
        search_low=1e-4,
        search_high=2.0,
        bisection_iterations=200,
        exactness_tolerance=1e-12,
    )
    assert report.breakdown_sigma is None
    assert report.relative_error_at_breakdown is None
    assert report.refusal_kind is ClosureBreakdownRefusal.BELOW_CANCELLATION_FLOOR
    # the two numbers the refusal is made of, quoted back to the caller
    assert report.cancellation_floor_at_low == pytest.approx(0.08871397861214594, rel=1e-9)
    assert report.relative_error_at_low == pytest.approx(0.0011703962790697465, rel=1e-9)
    assert report.cancellation_floor_at_low > report.tolerance
    assert report.relative_error_at_low <= report.cancellation_floor_at_low
    reason = report.refusal_reason or ""
    assert "cancellation floor" in reason
    assert "0.0001" in reason or "1e-04" in reason  # the bracket is named
    assert report.expansions  # a refusal must open an expansion route


def test_the_floor_gate_does_not_fire_on_a_well_conditioned_bracket() -> None:
    """The discriminating control: the same exponential WITHOUT the ``1e6`` offset is unaffected.

    Measured on the headline case (``A = 40``, ``lambda = 0.4``, ``mu = 0.8``, bracket
    ``[1e-3, 5.0]``): ``cancellation_floor_at_low = 1.42e-10`` against ``relative_error_at_low =
    1.56e-06``, four orders of margin, so the gate must stay silent and the quoted width must be
    unchanged.  Without this the fix could pass by refusing everything.
    """
    report = find_closure_breakdown_sigma(
        _exponential(),
        mean=MU,
        tolerance=0.01,
        search_low=1e-3,
        search_high=5.0,
        bisection_iterations=200,
        exactness_tolerance=1e-12,
    )
    assert report.refusal_kind is None
    assert report.breakdown_sigma is not None
    assert report.cancellation_floor_at_low == pytest.approx(1.4210854715327208e-10, rel=1e-9)
    assert report.relative_error_at_low == pytest.approx(1.5624911895333699e-06, rel=1e-9)
    assert report.cancellation_floor_at_low < report.tolerance
    assert report.relative_error_at_low > report.cancellation_floor_at_low
    # and the floor at the width actually issued is 12 orders below the tolerance it certifies
    assert gap_cancellation_floor(_exponential(), mean=MU, sigma=report.breakdown_sigma) < 1e-13


def test_the_ill_conditioned_law_is_measurable_once_the_bracket_is_re_declared() -> None:
    """The floor gate is about the BRACKET, not about the law — and it must say so.

    The same ``1e6 + exp(-x)`` measured over ``[0.05, 2.0]``, where the floor at the low end has
    fallen to 3.55e-07 — five orders below the 0.01 tolerance — and the measured error there,
    6.25e-04, sits three orders ABOVE its own floor, so it is truncation and not round-off.
    Re-declaring the bracket is the prescribed route; widening it inside the function after seeing
    the data is not.
    """
    law = _offset_exponential()
    floor = gap_cancellation_floor(law, mean=0.0, sigma=0.05)
    assert floor == pytest.approx(3.550497222770478e-07, rel=1e-9)
    report = find_closure_breakdown_sigma(
        law,
        mean=0.0,
        tolerance=0.01,
        search_low=0.05,
        search_high=2.0,
        bisection_iterations=200,
        exactness_tolerance=1e-12,
    )
    assert report.refusal_kind is None
    assert report.breakdown_sigma is not None
    # the constant contributes nothing to the gap, so the width is the pure exponential's
    assert report.breakdown_sigma == pytest.approx(
        _exponential_breakdown_closed_form(1.0, 0.01), rel=1e-6
    )


# --------------------------------------------------------------------------------------------
# 3. propagating (mu, Sigma) one step
# --------------------------------------------------------------------------------------------


OU_MATRIX = np.array([[-3.0, 1.0], [0.5, -2.0]])
OU_OFFSET = np.array([0.6, -0.4])
OU_DIFFUSION = np.array([[0.05, 0.01], [0.01, 0.03]])


def _ou_state():
    return make_gaussian_state(
        mean=[0.4, -0.2], covariance=[[0.02, 0.004], [0.004, 0.015]], contract=CONTRACT
    )


def test_local_gaussian_closure_is_exact_for_a_linear_drift() -> None:
    """Ornstein-Uhlenbeck: the Gaussian family is closed under affine drift, so the closure is not
    an approximation.  Measured agreement with the exact moment update: mean 0.0 exactly, covariance
    8.7e-19 absolute (4.3e-17 relative) — round-off between two different algebraic assemblies.
    """
    drift = LinearDrift(matrix=OU_MATRIX, offset=OU_OFFSET, name="ou-2d")
    comparison = compare_closure_against_ensemble(
        drift,
        _ou_state(),
        dt=1e-3,
        diffusion=OU_DIFFUSION,
        contract=CONTRACT,
        monte_carlo_samples=None,
        monte_carlo_seed=None,
    )
    assert comparison.closure_mean_error == 0.0
    assert comparison.closure_covariance_error < 1e-17
    assert comparison.relative_closure_covariance_error < 1e-15


def test_deterministic_ensemble_reproduces_the_exact_linear_update() -> None:
    """A genuinely different computation — a weighted sum over 225 transformed nodes rather than an
    algebraic moment formula — must land on the same numbers."""
    drift = LinearDrift(matrix=OU_MATRIX, offset=OU_OFFSET, name="ou-2d")
    state = _ou_state()
    quadrature = quadrature_moment_update(
        drift, state, dt=1e-3, diffusion=OU_DIFFUSION, contract=CONTRACT
    )
    exact = drift.exact_moment_update(state, dt=1e-3, diffusion=OU_DIFFUSION)
    assert np.max(np.abs(quadrature.mean - exact.mean)) < 1e-15
    assert np.max(np.abs(quadrature.covariance - exact.covariance)) < 1e-16


def test_sampled_ensemble_could_not_have_scored_this_closure() -> None:
    """Why truth here is analytic.  With 100 000 samples the sampled ensemble's error is 8.9e-04 in
    the mean and 2.1e-04 in the covariance — about **fifteen orders of magnitude** above the closure
    error it would have been used to detect.  Validating the closure against a sampled ensemble
    would have reported "no closure error" for a reason that has nothing to do with the closure.
    """
    drift = LinearDrift(matrix=OU_MATRIX, offset=OU_OFFSET, name="ou-2d")
    comparison = compare_closure_against_ensemble(
        drift,
        _ou_state(),
        dt=1e-3,
        diffusion=OU_DIFFUSION,
        contract=CONTRACT,
        monte_carlo_samples=100_000,
        monte_carlo_seed=20260729,
    )
    assert comparison.monte_carlo_mean_error is not None
    assert comparison.monte_carlo_mean_error == pytest.approx(8.862e-04, rel=5e-3)
    assert comparison.monte_carlo_covariance_error == pytest.approx(2.120e-04, rel=5e-3)
    assert comparison.monte_carlo_mean_error > 1e12 * max(
        comparison.closure_covariance_error, 1e-18
    )


def test_monte_carlo_ensemble_refuses_to_run_without_a_seed() -> None:
    drift = LinearDrift(matrix=OU_MATRIX, offset=OU_OFFSET, name="ou-2d")
    with pytest.raises(ValueError, match="explicit seed"):
        compare_closure_against_ensemble(
            drift,
            _ou_state(),
            dt=1e-3,
            diffusion=OU_DIFFUSION,
            contract=CONTRACT,
            monte_carlo_samples=1000,
            monte_carlo_seed=None,
        )


def test_euler_maruyama_step_error_against_the_exact_ou_flow_is_second_order_in_dt() -> None:
    """Discretisation error and closure error are different quantities and are separated here.

    Measured mean error against the exact continuous solution: 2.071e-05, 5.189e-06, 1.299e-06,
    3.248e-07 at ``dt = 4e-3 .. 5e-4`` — successive ratios 3.99, 4.00, 4.00, i.e. ``O(dt^2)``.  None
    of it is closure error; the closure is exact here.
    """
    drift = LinearDrift(matrix=OU_MATRIX, offset=OU_OFFSET, name="ou-2d")
    state = _ou_state()
    errors = []
    steps = (4e-3, 2e-3, 1e-3, 5e-4)
    for dt in steps:
        discrete = drift.exact_moment_update(state, dt=dt, diffusion=OU_DIFFUSION)
        continuous = propagate_gaussian_ou_exact(
            drift, state, dt=dt, diffusion=OU_DIFFUSION, contract=CONTRACT
        )
        errors.append(float(np.max(np.abs(discrete.mean - continuous.mean))))
    slope, _ = np.polyfit(np.log(np.asarray(steps)), np.log(np.asarray(errors)), 1)
    assert slope == pytest.approx(2.0, abs=0.01)
    assert errors[0] == pytest.approx(2.070857e-05, rel=1e-4)


def test_exact_ou_flow_refuses_an_unstable_drift() -> None:
    drift = LinearDrift(matrix=np.array([[0.5]]), offset=np.array([0.0]), name="unstable")
    state = make_gaussian_state(mean=[0.1], covariance=[[0.01]], contract=CONTRACT)
    with pytest.raises(ValueError, match="stable drift"):
        propagate_gaussian_ou_exact(drift, state, dt=1e-3, diffusion=[[0.01]], contract=CONTRACT)


def test_weak_nonlinearity_closure_is_good_but_not_exact_and_the_error_grows_as_sigma_squared() -> (
    None
):
    """``F = -(4 x + 0.3 x^3)`` at ``mu = 0``.  The MEAN closure is exact (a cubic has no fourth
    derivative), and the error lives entirely in the covariance, where the closure drops the exact
    ``Cov(X, F)`` and ``Var(F)`` terms.

    Measured relative covariance error: 4.446e-06, 1.800e-05, 7.222e-05, 2.891e-04 at
    ``sigma = 0.05, 0.1, 0.2, 0.4`` — a fitted exponent of 2.00 in ``sigma``.  The deterministic
    ensemble is exact at every one of those widths, so the discrepancy is the closure's, not the
    comparison's.
    """
    drift = ScalarLawDrift(
        law=cubic_anharmonic_spring(stiffness=4.0, cubic_coefficient=0.3, unit="pN")
    )
    widths = (0.05, 0.1, 0.2, 0.4)
    relative = []
    for sigma in widths:
        state = make_gaussian_state(mean=[0.0], covariance=[[sigma**2]], contract=CONTRACT)
        comparison = compare_closure_against_ensemble(
            drift,
            state,
            dt=1e-3,
            diffusion=[[0.02]],
            contract=CONTRACT,
            monte_carlo_samples=None,
            monte_carlo_seed=None,
        )
        assert comparison.closure_mean_error == 0.0
        assert comparison.quadrature_covariance_error < 1e-16
        relative.append(comparison.relative_closure_covariance_error)
    assert relative[0] == pytest.approx(4.446369e-06, rel=1e-4)
    assert relative[-1] == pytest.approx(2.891325e-04, rel=1e-4)
    slope, _ = np.polyfit(np.log(np.asarray(widths)), np.log(np.asarray(relative)), 1)
    assert slope == pytest.approx(2.0, abs=0.02)


# --------------------------------------------------------------------------------------------
# the bimodal refusal
# --------------------------------------------------------------------------------------------


def _double_well() -> AnalyticForceLaw:
    """``F(x) = x - x^3``, i.e. ``U = x^4/4 - x^2/2`` with minima at +-1 and a barrier of 1/4."""
    return cubic_anharmonic_spring(stiffness=-1.0, cubic_coefficient=1.0, unit="pN")


DOUBLE_WELL_GRID = np.linspace(-3.0, 3.0, 4001)


def test_double_well_potential_is_the_one_this_test_thinks_it_is() -> None:
    law = _double_well()
    assert law.evaluate(np.array([0.0, 1.0, -1.0, 2.0])) == pytest.approx([0.0, 0.0, 0.0, -6.0])
    # U(0) - U(1) = 0 - (-1/4) = 1/4
    assert float(law.potential(np.array(0.0)) - law.potential(np.array(1.0))) == pytest.approx(0.25)


@pytest.mark.parametrize(
    ("noise", "expected_ratio", "expected_coefficient"),
    [
        (0.2, 2.865052e-01, 0.669697),
        (0.1, 8.208520e-02, 0.781658),
        (0.05, 6.737981e-03, 0.889367),
        (0.02, 3.726700e-06, 0.958931),
    ],
)
def test_bimodality_statistics_are_measured_against_the_analytic_boltzmann_density(
    noise: float, expected_ratio: float, expected_coefficient: float
) -> None:
    """The reference is closed form: ``rho ~ exp(-U(x)/D)``, so these are truth-scored.

    The density ratio at the closure mean is ``exp(-(U(0) - U(1)) / D) = exp(-0.25 / D)``, checked
    here independently; the agreement is limited to ~1e-5 relative only because the 4001-point grid
    puts its peak at 1.0005 rather than 1.0, which shifts the reference by ``U''(1) (5e-4)^2 / 2``.
    At ``D = 0.05`` the single-Gaussian mean sits where the density is **0.67 % of the modal
    density**; at ``D = 0.02``, **3.7e-06 of it**.
    """
    law = _double_well()
    density = boltzmann_density(law, DOUBLE_WELL_GRID, noise_intensity=noise)
    report = measure_bimodality(DOUBLE_WELL_GRID, density, evaluation_point=0.0, contract=CONTRACT)
    assert report.mode_count == 2
    assert report.mode_positions[0] == pytest.approx(-1.0, abs=2e-3)
    assert report.mode_positions[1] == pytest.approx(1.0, abs=2e-3)
    assert report.density_ratio_at_evaluation_point == pytest.approx(
        math.exp(-0.25 / noise), rel=5e-5
    )
    assert report.density_ratio_at_evaluation_point == pytest.approx(expected_ratio, rel=1e-5)
    assert report.bimodality_coefficient == pytest.approx(expected_coefficient, rel=1e-5)
    assert report.bimodality_coefficient > CONTRACT.bimodality_coefficient_max
    assert not report.is_gaussian_representable


def test_a_gaussian_reference_density_passes_every_clause() -> None:
    """The positive control for the refusal: on a Gaussian the coefficient is ``1/3`` exactly."""
    grid = np.linspace(-6.0, 6.0, 4001)
    density = np.exp(-0.5 * grid**2)
    report = measure_bimodality(grid, density, evaluation_point=0.0, contract=CONTRACT)
    assert report.mode_count == 1
    assert report.kurtosis == pytest.approx(3.0, rel=1e-6)
    assert report.bimodality_coefficient == pytest.approx(1.0 / 3.0, rel=1e-6)
    assert report.density_ratio_at_evaluation_point == pytest.approx(1.0, rel=1e-12)
    assert report.is_gaussian_representable


def test_closure_refuses_the_double_well_and_exposes_no_mean_to_read() -> None:
    """The plan's own prohibition, enforced structurally.

    The refusal is driven by three measured statistics, and a refused verdict has no propagated
    state at all — ``step`` is ``None`` and ``require_step`` raises — so the barrier-top mean cannot
    be read past the refusal.
    """
    drift = ScalarLawDrift(law=_double_well())
    state = make_gaussian_state(mean=[0.0], covariance=[[0.05]], contract=CONTRACT)
    density = boltzmann_density(_double_well(), DOUBLE_WELL_GRID, noise_intensity=0.05)
    verdict = propagate_local_gaussian(
        drift,
        state,
        dt=1e-3,
        diffusion=[[0.05]],
        contract=CONTRACT,
        reference_grid=DOUBLE_WELL_GRID,
        reference_density=density,
    )
    assert verdict.kind is ClosureVerdictKind.REFUSED_BIMODAL
    assert not verdict.admitted
    assert verdict.step is None
    payload = verdict.as_dict()
    assert payload["mean"] is None
    assert payload["covariance"] is None
    assert payload["contract"] == CONTRACT.to_dict()
    assert verdict.fallback_representation is RepresentationKind.GAUSSIAN_MIXTURE
    assert len(verdict.reasons) == 3
    with pytest.raises(ClosureRefusedError, match="refused-bimodal"):
        verdict.require_step()


def test_closure_refuses_a_nonlinear_drift_with_no_shape_evidence() -> None:
    """Fail closed.  A Gaussian closure applied without any measurement of the true density's shape
    is the untested assumption, not a small approximation."""
    drift = ScalarLawDrift(law=_double_well())
    state = make_gaussian_state(mean=[0.0], covariance=[[0.05]], contract=CONTRACT)
    verdict = propagate_local_gaussian(drift, state, dt=1e-3, diffusion=[[0.05]], contract=CONTRACT)
    assert verdict.kind is ClosureVerdictKind.REFUSED_UNVERIFIED_SHAPE
    assert verdict.step is None
    assert verdict.fallback_representation is RepresentationKind.EMPIRICAL_ENSEMBLE


def test_closure_is_admitted_for_a_linear_drift_without_shape_evidence() -> None:
    """The complement: for an affine drift the Gaussian family is provably closed, so no shape
    evidence is needed and the closure is admitted."""
    drift = LinearDrift(matrix=OU_MATRIX, offset=OU_OFFSET, name="ou-2d")
    verdict = propagate_local_gaussian(
        drift, _ou_state(), dt=1e-3, diffusion=OU_DIFFUSION, contract=CONTRACT
    )
    assert verdict.kind is ClosureVerdictKind.ADMITTED
    assert verdict.mean_truncation_error == 0.0
    step = verdict.require_step()
    assert step.representation is NodeRepresentation.LOCAL_GAUSSIAN
    assert step.state.dimension == 2


def test_closure_refuses_a_width_beyond_the_range_the_expansion_was_validated_over() -> None:
    """The section-2 breakdown width, wired into the propagation gate.

    The exponential drift at ``sigma = 0.30`` (well past the measured 1 %-safe width of
    ``0.2003 lambda = 0.0801``) has a closure mean-drift error of 0.0329 against a declared 0.01, so
    the step is refused; a unimodal reference density is supplied so the refusal cannot be the
    bimodality clause.
    """
    law = _exponential()
    drift = ScalarLawDrift(law=law)
    state = make_gaussian_state(mean=[MU], covariance=[[0.30**2]], contract=CONTRACT)
    grid = np.linspace(MU - 6 * 0.30, MU + 6 * 0.30, 2001)
    density = np.exp(-0.5 * ((grid - MU) / 0.30) ** 2)
    verdict = propagate_local_gaussian(
        drift,
        state,
        dt=1e-3,
        diffusion=[[0.01]],
        contract=CONTRACT,
        reference_grid=grid,
        reference_density=density,
    )
    assert verdict.kind is ClosureVerdictKind.REFUSED_WIDTH_BEYOND_VALIDATED_RANGE
    assert verdict.step is None
    assert verdict.bimodality is not None and verdict.bimodality.is_gaussian_representable
    rate = -1.0 / DECAY
    exact = AMPLITUDE * math.exp(rate * MU + 0.5 * rate**2 * 0.30**2)
    closure = AMPLITUDE * math.exp(rate * MU) * (1.0 + 0.5 * rate**2 * 0.30**2)
    assert verdict.mean_truncation_error == pytest.approx(abs(exact - closure) / exact, rel=1e-12)
    assert verdict.mean_truncation_error == pytest.approx(0.0328618, rel=1e-5)
    assert verdict.mean_truncation_error > CONTRACT.expansion_relative_error_tolerance


def test_truncated_reference_grid_raises_from_the_reused_mass_audit() -> None:
    """The composition is load-bearing: a grid that does not contain its own tails is caught by
    ``audit_probability_mass`` on the rectangle sum, not renormalised into a wrong answer.  Measured
    mass error on a grid clipped to ``[-0.6, 0.6]``: 0.0106.
    """
    narrow = np.linspace(-0.6, 0.6, 401)
    density = boltzmann_density(_double_well(), narrow, noise_intensity=0.05)
    with pytest.raises(ValueError, match="not normalized"):
        measure_bimodality(narrow, density, evaluation_point=0.0, contract=CONTRACT)


def test_bimodality_refuses_a_nonuniform_grid_and_a_point_outside_it() -> None:
    grid = np.concatenate([np.linspace(-3.0, 0.0, 100), np.linspace(0.01, 3.0, 200)])
    density = boltzmann_density(_double_well(), grid, noise_intensity=0.05)
    with pytest.raises(ValueError, match="uniformly spaced"):
        measure_bimodality(grid, density, evaluation_point=0.0, contract=CONTRACT)
    uniform = np.linspace(-3.0, 3.0, 4001)
    with pytest.raises(ValueError, match="outside the reference grid"):
        measure_bimodality(
            uniform,
            boltzmann_density(_double_well(), uniform, noise_intensity=0.05),
            evaluation_point=9.0,
            contract=CONTRACT,
        )


# --------------------------------------------------------------------------------------------
# 4. source-factored, not one Sigma
# --------------------------------------------------------------------------------------------


def test_same_total_variance_in_two_places_gives_two_different_mean_forces() -> None:
    """The concrete argument for source factoring.

    ``A exp(b x)`` with ``A = 40``, ``b = -2.5``, ``mu = 0.8`` and total variance ``v = 0.04``:
    the point force is 5.4134, placing ``v`` on the position gives 6.1342, placing the SAME ``v`` on
    the rate gives 5.4831.  Difference 0.6511, ratio 1.11874 — and both corrections are nonzero, so
    this is not the trivial case of one placement doing nothing.
    """
    contrast = measure_rate_placement_contrast(
        amplitude=AMPLITUDE, rate=-2.5, position=MU, total_variance=0.04
    )
    assert contrast.point_force == pytest.approx(40.0 * math.exp(-2.5 * 0.8), rel=1e-14)
    assert contrast.first_expectation == pytest.approx(
        40.0 * math.exp(-2.5 * 0.8 + 0.5 * 2.5**2 * 0.04), rel=1e-14
    )
    assert contrast.second_expectation == pytest.approx(
        40.0 * math.exp(-2.5 * 0.8 + 0.5 * 0.8**2 * 0.04), rel=1e-14
    )
    assert contrast.first_expectation == pytest.approx(6.134198674, rel=1e-9)
    assert contrast.second_expectation == pytest.approx(5.483148359, rel=1e-9)
    assert contrast.difference == pytest.approx(0.6510503145, rel=1e-8)
    assert contrast.ratio == pytest.approx(contrast.closed_form_ratio, rel=1e-14)
    assert contrast.ratio == pytest.approx(1.1187365856, rel=1e-9)
    assert set(contrast.sources) == {UncertaintySource.THERMAL, UncertaintySource.PARAMETER}


def test_cubic_placement_makes_one_of_the_two_corrections_exactly_zero() -> None:
    """The sharpest form: ``-(k x + a x^3)`` is cubic in ``x`` and LINEAR in ``k``.

    Position variance gives the Jensen correction ``-3 a mu v = -0.1632``; the same variance on the
    stiffness gives exactly zero.  A pre-merged total variance cannot recover which one it was.
    """
    contrast = measure_cubic_placement_contrast(
        stiffness=STIFFNESS, cubic_coefficient=CUBIC, position=MU, total_variance=0.04
    )
    assert contrast.second_expectation == contrast.point_force
    assert contrast.difference == pytest.approx(-3.0 * CUBIC * MU * 0.04, rel=1e-13)
    assert contrast.difference == pytest.approx(-0.1632, rel=1e-12)


def test_source_factored_covariance_merges_only_when_explicitly_asked() -> None:
    """Two different factorings with the SAME merged total.

    The total is identical to round-off, so anything downstream that sees only the total cannot
    distinguish them — which is exactly what the two contrasts above show it must.
    """
    thermal_heavy = SourceFactoredCovariance(
        blocks={
            UncertaintySource.THERMAL: np.array([[0.03]]),
            UncertaintySource.PARAMETER: np.array([[0.01]]),
        },
        axes=("bond-extension",),
    )
    epistemic_heavy = SourceFactoredCovariance(
        blocks={
            UncertaintySource.THERMAL: np.array([[0.01]]),
            UncertaintySource.PARAMETER: np.array([[0.03]]),
        },
        axes=("bond-extension",),
    )
    assert float(thermal_heavy.merged_total()[0, 0]) == float(epistemic_heavy.merged_total()[0, 0])
    assert float(thermal_heavy.blocks[UncertaintySource.THERMAL][0, 0]) != float(
        epistemic_heavy.blocks[UncertaintySource.THERMAL][0, 0]
    )
    assert thermal_heavy.sources == (UncertaintySource.THERMAL, UncertaintySource.PARAMETER)
    with pytest.raises(ValueError, match="at least one source block"):
        SourceFactoredCovariance(blocks={}, axes=("x",))
    with pytest.raises(TypeError, match="UncertaintySource"):
        SourceFactoredCovariance(blocks={"thermal": np.array([[0.1]])}, axes=("x",))


def test_quenched_and_thermal_split_leaves_the_mean_invariant_and_the_variance_is_not() -> None:
    """**The half of the claim that is false, measured.**

    For ``A exp(b x)`` over ``N = 64`` bonds at fixed total variance ``v = 0.04``, the MEAN
    aggregate force is 6.1341986737971 for every split from all-thermal to all-quenched — identical
    to 0.0e+00, not merely close.  So "sources must be kept separate because the mean moves" is not
    the argument here.

    What does move is the variance of the aggregate: 0.16699 (all thermal) to 10.6874 (all
    quenched), a ratio of exactly **64.000 = N**, because the quenched offset is shared and does not
    average away.  Any downstream nonlinear functional of the aggregate force therefore differs, and
    a single pre-merged ``Sigma`` cannot express the difference.
    """
    means = []
    variances = []
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        split = measure_quenched_thermal_split(
            amplitude=AMPLITUDE,
            rate=-2.5,
            position=MU,
            quenched_variance=0.04 * fraction,
            thermal_variance=0.04 * (1.0 - fraction),
            n_units=64,
        )
        assert split.total_variance == pytest.approx(0.04, rel=1e-15)
        means.append(split.mean_aggregate_force)
        variances.append(split.variance_aggregate_force)

    assert max(abs(value - means[0]) for value in means) == 0.0
    assert means[0] == pytest.approx(
        AMPLITUDE * math.exp(-2.5 * MU + 0.5 * 2.5**2 * 0.04), rel=1e-14
    )
    assert variances[-1] / variances[0] == pytest.approx(64.0, rel=1e-10)
    assert variances[0] == pytest.approx(0.1669909392, rel=1e-8)
    assert variances[-1] == pytest.approx(10.68742011, rel=1e-8)


def test_quenched_thermal_limits_match_their_hand_derived_lognormal_forms() -> None:
    """Held out from the closed forms above: the two limits are ordinary lognormal variances.

    all-quenched: ``A^2 e^{2 b mu} (e^{2 b^2 v} - e^{b^2 v})``; all-thermal: the same over ``N``.
    Checked at ``N = 13``, a bond count not used anywhere else in this file.
    """
    rate = -2.5
    variance = 0.04
    base = AMPLITUDE**2 * math.exp(2.0 * rate * MU)
    lognormal = base * (math.exp(2.0 * rate**2 * variance) - math.exp(rate**2 * variance))
    quenched = measure_quenched_thermal_split(
        amplitude=AMPLITUDE,
        rate=rate,
        position=MU,
        quenched_variance=variance,
        thermal_variance=0.0,
        n_units=13,
    )
    thermal = measure_quenched_thermal_split(
        amplitude=AMPLITUDE,
        rate=rate,
        position=MU,
        quenched_variance=0.0,
        thermal_variance=variance,
        n_units=13,
    )
    assert quenched.variance_aggregate_force == pytest.approx(lognormal, rel=1e-12)
    assert thermal.variance_aggregate_force == pytest.approx(lognormal / 13.0, rel=1e-12)


# --------------------------------------------------------------------------------------------
# validation and pre-registration
# --------------------------------------------------------------------------------------------


def test_contract_has_no_default_on_any_field() -> None:
    """Matching ``RegimeEvidenceContract``: a threshold with a default is a threshold nobody chose."""
    import dataclasses

    for declared in dataclasses.fields(MomentClosureContract):
        assert declared.default is dataclasses.MISSING, declared.name
        assert declared.default_factory is dataclasses.MISSING, declared.name


@pytest.mark.parametrize(
    ("field_name", "bad_value"),
    [
        ("quadrature_nodes", 1),
        ("gap_relative_tolerance", -1e-3),
        ("mode_prominence_ratio", 1.5),
        ("modal_density_ratio_min", 2.0),
        ("bimodality_coefficient_max", 1.5),
    ],
)
def test_contract_rejects_malformed_thresholds(field_name: str, bad_value: float) -> None:
    import dataclasses

    with pytest.raises(ValueError):
        dataclasses.replace(CONTRACT, **{field_name: bad_value})


def test_gaussian_state_rejects_a_non_psd_or_asymmetric_covariance() -> None:
    with pytest.raises(ValueError, match="positive semidefinite"):
        make_gaussian_state(mean=[0.0], covariance=[[-1e-3]], contract=CONTRACT)
    with pytest.raises(ValueError, match="symmetric"):
        make_gaussian_state(mean=[0.0, 0.0], covariance=[[1.0, 0.5], [0.4, 1.0]], contract=CONTRACT)
    singular = make_gaussian_state(
        mean=[0.0, 0.0], covariance=[[1.0, 0.0], [0.0, 0.0]], contract=CONTRACT
    )
    assert np.allclose(singular.matrix_sqrt(), np.array([[1.0, 0.0], [0.0, 0.0]]))


def test_scalar_sigma_refuses_a_multivariate_state() -> None:
    """A multivariate covariance has no single width; taking ``sqrt(trace)`` would be exactly the
    pre-merge this module argues against."""
    state = _ou_state()
    with pytest.raises(ValueError, match="one-dimensional"):
        _ = state.scalar_sigma


def test_force_law_evaluation_matches_a_hand_written_expression() -> None:
    law = _morse()
    depth, width, centre = 12.0, 2.0, MU
    points = np.array([0.4, MU, 1.3])
    expected = (
        2.0
        * width
        * depth
        * (np.exp(-2.0 * width * (points - centre)) - np.exp(-width * (points - centre)))
    )
    assert law.evaluate(points) == pytest.approx(expected, rel=1e-12, abs=1e-12)
    # numerical derivative check on the second derivative used by the closure
    h = 1e-4
    numerical = (
        law.evaluate(points + h) - 2.0 * law.evaluate(points) + law.evaluate(points - h)
    ) / h**2
    assert law.derivative(2, points) == pytest.approx(numerical, rel=1e-6)


def test_potential_is_the_antiderivative_of_the_force() -> None:
    for law in (_cubic(), _quartic(), _exponential(), _morse()):
        points = np.array([-0.3, 0.2, 0.9])
        h = 1e-5
        numerical = -(law.potential(points + h) - law.potential(points - h)) / (2.0 * h)
        assert law.evaluate(points) == pytest.approx(numerical, rel=1e-6, abs=1e-8), law.name


def test_experiment_card_is_pre_registered_with_every_required_section() -> None:
    card = stochastic_node_experiment_card(["a" * 64])
    assert isinstance(card, SandboxExperimentCard)
    assert card.experiment_id == "D1-stochastic-node-moment-oracle"
    assert not card.wet_lab
    assert card.budget_class == "cpu-analytic-oracle"
    for section in (
        card.positive_controls,
        card.negative_controls,
        card.adversarial_controls,
        card.held_out_interventions,
        card.failure_expansions,
        card.observables,
    ):
        assert len(section) >= 3
    assert "exactly 0.0" in " ".join(card.positive_controls)
    assert "REFUSE" in " ".join(card.adversarial_controls)
    assert len(card.card_hash) == 64
    assert card.card_hash == stochastic_node_experiment_card(["a" * 64]).card_hash

"""S1 passive-fluctuation canon — estimator validation and the four defeats it must catch.

The tests are organised as the sandbox card is: a dimensional gate, the analytic canon, one positive
control whose truth is exact, four negative controls, and the adjudication that combines them.  The
negative controls are the point.  A test that only checks the positive control would pass for an
estimator that returns the answer it was given.
"""

from __future__ import annotations

import math
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from aleph.virtual_cell.contracts import EvidenceSource, SandboxExperimentCard
from aleph.virtual_cell.sandbox_fluctuation import (
    BOLTZMANN_J_PER_K,
    CUTOFF_DAMAGE_IS_NOT_A_MULTIPLE,
    ENGINE_KAPPA_J,
    ENGINE_MEMBRANE_PROVENANCE,
    ENGINE_SIGMA_N_PER_M,
    EQUIPARTITION_SCALE_DEGENERACY,
    EVIDENCE_AUTHORITY,
    EVIDENCE_SOURCE,
    FAILURE_EXPANSIONS,
    HELFRICH_FAMILY_ABSORPTION,
    PATCH_BOUNDARY_CONDITION_IS_A_DECLARATION,
    PATCH_PERIODIC_RATIONALE,
    PHYSIOLOGICAL_TEMPERATURE_K,
    Dimension,
    FluctuationFailure,
    FluctuationRegime,
    HelfrichSpectrum,
    PatchBoundaryCondition,
    TautStringSpectrum,
    adjudicate_fluctuation_spectrum,
    check_equipartition,
    check_mode_independence,
    chi_square_survival,
    correlated_mode_amplitudes,
    decades_required_for_exponent_separation,
    fit_helfrich_spectrum,
    fit_power_law_exponent,
    helfrich_variance_dimension,
    kappa_from_pn_um,
    log_spaced_wavevectors,
    mean_energy_per_quadratic_dof,
    measure_cutoff_damage_ladder,
    measure_exponent_separation_ladder,
    measure_mode_variance,
    mesh_nyquist_wavevector,
    mesh_truncated_variances,
    native_patch_requirement,
    power_law_variances,
    s1_experiment_card,
    sample_gaussian_modes,
    sigma_from_pn_per_um,
    thermal_energy,
    trusted_wavevector_range,
)

# A 5 um x 5 um patch: a realistic slice of plasma membrane, not a unit square.
PATCH_AREA_M2 = (5.0e-6) ** 2
SEED = 20260729
MANIFEST_HASH = "b" * 64


def _reference_spectrum() -> HelfrichSpectrum:
    return HelfrichSpectrum(
        kappa_J=ENGINE_KAPPA_J,
        sigma_N_per_m=ENGINE_SIGMA_N_PER_M,
        area_m2=PATCH_AREA_M2,
        temperature_K=PHYSIOLOGICAL_TEMPERATURE_K,
    )


def _sampled_spectrum(
    spectrum: HelfrichSpectrum,
    *,
    q: np.ndarray,
    n_samples: int,
    seed: int = SEED,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    amplitudes = sample_gaussian_modes(spectrum.mode_variance(q), n_samples=n_samples, rng=rng)
    return amplitudes, measure_mode_variance(amplitudes)


# --------------------------------------------------------------------------------------------
# 0. the lane's own hard constraint
# --------------------------------------------------------------------------------------------


def test_importing_the_canon_does_not_pull_in_warp() -> None:
    """Importing and exercising this module must pull in no GPU runtime.

    Run in a SUBPROCESS.  Asserting against this process's ``sys.modules`` would measure whether
    anything in the whole pytest session ever imported Warp, which is a different and false claim:
    the suite contains tests that legitimately import it, so the in-process form passes alone and
    fails under the full run.  A fresh interpreter measures what the name says.
    """
    script = textwrap.dedent(
        """
        import sys
        import numpy as np
        from aleph.virtual_cell.sandbox_fluctuation import (
            HelfrichSpectrum, fit_helfrich_spectrum, log_spaced_wavevectors,
            measure_mode_variance, sample_gaussian_modes,
        )

        area = 25e-12
        spectrum = HelfrichSpectrum(8.28e-20, 1e-5, area, 310.0)
        q = log_spaced_wavevectors(1e6, 1e9, 30)
        rng = np.random.default_rng(0)
        amps = sample_gaussian_modes(spectrum.mode_variance(q), n_samples=500, rng=rng)
        fit = fit_helfrich_spectrum(
            q, measure_mode_variance(amps), area_m2=area, temperature_K=310.0, n_samples=500
        )
        assert fit.converged
        for blocked in ("warp", "torch", "jax"):
            assert blocked not in sys.modules, f"sandbox_fluctuation import pulled in {blocked}"
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


def test_the_module_declares_itself_an_oracle_and_never_native_evidence() -> None:
    assert EVIDENCE_SOURCE is EvidenceSource.ANALYTIC_ORACLE
    assert EVIDENCE_AUTHORITY == "unverified-observer"
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e6, 1.0e9, 20)
    amplitudes, _ = _sampled_spectrum(spectrum, q=q, n_samples=500)
    verdict = adjudicate_fluctuation_spectrum(q, amplitudes, reference=spectrum)
    assert verdict.to_dict()["evidence_source"] == EvidenceSource.ANALYTIC_ORACLE.value
    assert verdict.to_dict()["evidence_authority"] == "unverified-observer"


def test_engine_constant_provenance_is_marked_inherited_not_verified() -> None:
    """The two moduli are repository constants with citations this lane did not read; say so."""
    assert "NOT READ BY THIS LANE" in ENGINE_MEMBRANE_PROVENANCE["kappa_inherited_citation"]
    assert "NOT READ BY THIS LANE" in ENGINE_MEMBRANE_PROVENANCE["sigma_inherited_citation"]
    assert ENGINE_MEMBRANE_PROVENANCE["status"].startswith("INHERITED-UNVERIFIED")
    # And the arithmetic in the provenance note is itself checked, not asserted in prose.
    assert kappa_from_pn_um(0.0828) == pytest.approx(ENGINE_KAPPA_J, rel=1e-12)
    assert sigma_from_pn_per_um(10.0) == pytest.approx(ENGINE_SIGMA_N_PER_M, rel=1e-12)
    twenty_kt_300 = 20.0 * BOLTZMANN_J_PER_K * 300.0 / 1.0e-18
    twenty_kt_310 = 20.0 * BOLTZMANN_J_PER_K * 310.0 / 1.0e-18
    assert twenty_kt_300 == pytest.approx(0.0828, rel=5e-4)
    assert twenty_kt_310 == pytest.approx(0.0856, rel=5e-4)


# --------------------------------------------------------------------------------------------
# 1. dimensional gate — the classic wrong-power-of-q failure
# --------------------------------------------------------------------------------------------


def test_helfrich_variance_composes_to_an_area() -> None:
    assert helfrich_variance_dimension() == Dimension(length=2.0)


@pytest.mark.parametrize(
    ("bending_power", "tension_power"),
    [(3, 2), (5, 2), (4, 1), (4, 3), (2, 2), (4, 4)],
)
def test_a_wrong_power_of_q_is_caught_by_the_dimensional_check(
    bending_power: int, tension_power: int
) -> None:
    """The failure mode the Sanity Gate names: a wrong exponent that still looks like a spectrum."""
    with pytest.raises(ValueError, match="not commensurate|must be an area"):
        helfrich_variance_dimension(bending_power=bending_power, tension_power=tension_power)


def test_dimension_algebra_is_a_real_algebra() -> None:
    length = Dimension(length=1.0)
    assert length**2 / length == length
    assert (length * length).length == 2.0
    assert str(Dimension()) == "1"
    assert str(Dimension(length=2.0)) == "L^2"


def test_unit_conversions_are_exact_powers_of_ten() -> None:
    assert kappa_from_pn_um(1.0) == 1.0e-18  # 1 pN*um = 1e-12 N * 1e-6 m
    assert sigma_from_pn_per_um(1.0) == 1.0e-6  # 1 pN/um = 1e-12 N / 1e-6 m
    assert thermal_energy(310.0) == BOLTZMANN_J_PER_K * 310.0


def test_scaling_law_is_verified_numerically_not_only_dimensionally() -> None:
    """Independent of the dimension algebra: halve q and the two limits must scale by 16 and 4."""
    spectrum = _reference_spectrum()
    q = np.array([2.0e8, 1.0e8])
    bending = spectrum.bending_limit_variance(q)
    tension = spectrum.tension_limit_variance(q)
    assert bending[1] / bending[0] == pytest.approx(16.0, rel=1e-12)
    assert tension[1] / tension[0] == pytest.approx(4.0, rel=1e-12)


# --------------------------------------------------------------------------------------------
# 2. the analytic canon
# --------------------------------------------------------------------------------------------


def test_crossover_wavevector_is_where_the_two_terms_are_equal() -> None:
    spectrum = _reference_spectrum()
    q_star = spectrum.crossover_wavevector_per_m
    assert q_star == pytest.approx(math.sqrt(ENGINE_SIGMA_N_PER_M / ENGINE_KAPPA_J), rel=1e-12)
    bending = ENGINE_KAPPA_J * q_star**4
    tension = ENGINE_SIGMA_N_PER_M * q_star**2
    assert bending == pytest.approx(tension, rel=1e-12)
    assert spectrum.crossover_wavelength_m == pytest.approx(2.0 * math.pi / q_star, rel=1e-12)


def test_crossover_lands_at_a_physically_meaningful_length_for_the_engine_constants() -> None:
    """A checkable prediction, stated in nanometres so it can be argued with."""
    spectrum = _reference_spectrum()
    assert spectrum.crossover_wavevector_per_m == pytest.approx(1.0990e7, rel=1e-3)
    assert spectrum.crossover_wavelength_m * 1e9 == pytest.approx(571.7, rel=1e-3)


def test_full_spectrum_reduces_to_each_limit_far_from_the_crossover() -> None:
    spectrum = _reference_spectrum()
    q_star = spectrum.crossover_wavevector_per_m
    high = np.array([1000.0 * q_star])
    low = np.array([q_star / 1000.0])
    assert spectrum.mode_variance(high)[0] == pytest.approx(
        spectrum.bending_limit_variance(high)[0], rel=1e-5
    )
    assert spectrum.mode_variance(low)[0] == pytest.approx(
        spectrum.tension_limit_variance(low)[0], rel=1e-5
    )
    # ... and is strictly BELOW either limit at the crossover, because both terms contribute.
    at_star = np.array([q_star])
    assert spectrum.mode_variance(at_star)[0] < spectrum.bending_limit_variance(at_star)[0]
    assert spectrum.mode_variance(at_star)[0] < spectrum.tension_limit_variance(at_star)[0]


def test_regime_labels_and_expected_exponents_track_the_crossover() -> None:
    spectrum = _reference_spectrum()
    q_star = spectrum.crossover_wavevector_per_m
    assert spectrum.regime_at(q_star / 10.0) is FluctuationRegime.TENSION_DOMINATED
    assert spectrum.regime_at(q_star * 10.0) is FluctuationRegime.BENDING_DOMINATED
    assert spectrum.regime_at(q_star) is FluctuationRegime.CROSSOVER
    assert spectrum.expected_exponent(q_star / 10.0) == -2.0
    assert spectrum.expected_exponent(q_star * 10.0) == -4.0
    assert spectrum.expected_exponent(q_star) is None


def test_pure_limits_have_degenerate_but_correct_crossovers() -> None:
    bending_only = HelfrichSpectrum(ENGINE_KAPPA_J, 0.0, PATCH_AREA_M2, 310.0)
    tension_only = HelfrichSpectrum(0.0, ENGINE_SIGMA_N_PER_M, PATCH_AREA_M2, 310.0)
    assert bending_only.crossover_wavevector_per_m == 0.0
    assert bending_only.crossover_wavelength_m == math.inf
    assert math.isinf(tension_only.crossover_wavevector_per_m)
    assert tension_only.crossover_wavelength_m == 0.0
    assert bending_only.regime_at(1.0e8) is FluctuationRegime.BENDING_DOMINATED
    assert tension_only.regime_at(1.0e8) is FluctuationRegime.TENSION_DOMINATED
    with pytest.raises(ValueError, match="undefined for sigma = 0"):
        bending_only.tension_limit_variance(np.array([1.0e8]))
    with pytest.raises(ValueError, match="undefined for kappa = 0"):
        tension_only.bending_limit_variance(np.array([1.0e8]))


def test_spectrum_boundary_cases_are_refused_not_regularised() -> None:
    with pytest.raises(ValueError, match="cannot both be zero"):
        HelfrichSpectrum(0.0, 0.0, PATCH_AREA_M2, 310.0)
    with pytest.raises(ValueError, match="non-negative"):
        HelfrichSpectrum(-1.0e-20, 1.0e-5, PATCH_AREA_M2, 310.0)
    with pytest.raises(ValueError, match="positive"):
        HelfrichSpectrum(ENGINE_KAPPA_J, 1.0e-5, 0.0, 310.0)
    with pytest.raises(ValueError, match="positive"):
        HelfrichSpectrum(ENGINE_KAPPA_J, 1.0e-5, PATCH_AREA_M2, 0.0)
    spectrum = _reference_spectrum()
    with pytest.raises(ValueError, match="translation mode"):
        spectrum.mode_variance(np.array([0.0, 1.0e8]))


def test_variance_scales_as_temperature_and_inverse_area() -> None:
    """The two held-out interventions on the card, checked analytically here."""
    hot = HelfrichSpectrum(ENGINE_KAPPA_J, ENGINE_SIGMA_N_PER_M, PATCH_AREA_M2, 620.0)
    cold = HelfrichSpectrum(ENGINE_KAPPA_J, ENGINE_SIGMA_N_PER_M, PATCH_AREA_M2, 310.0)
    big = HelfrichSpectrum(ENGINE_KAPPA_J, ENGINE_SIGMA_N_PER_M, 4.0 * PATCH_AREA_M2, 310.0)
    q = np.array([1.0e8])
    assert hot.mode_variance(q)[0] / cold.mode_variance(q)[0] == pytest.approx(2.0, rel=1e-12)
    assert cold.mode_variance(q)[0] / big.mode_variance(q)[0] == pytest.approx(4.0, rel=1e-12)


# --------------------------------------------------------------------------------------------
# 3. the 1-D analogue — the case with no approximation at all
# --------------------------------------------------------------------------------------------


def test_taut_string_variance_is_exactly_kT_over_the_mode_stiffness() -> None:
    string = TautStringSpectrum(tension_N=1.0e-11, length_m=1.0e-5, temperature_K=310.0)
    n = np.arange(1, 21)
    q = string.mode_wavevector(n)
    assert q[0] == pytest.approx(math.pi / 1.0e-5, rel=1e-12)
    stiffness = string.mode_stiffness(n)
    assert np.allclose(stiffness, 0.5 * 1.0e-11 * 1.0e-5 * q**2, rtol=1e-12)
    assert np.allclose(string.mode_variance(n), string.thermal_energy_J / stiffness, rtol=1e-12)
    # closed form: <a_n^2> = 2 kB T / (tau L q_n^2)
    closed = 2.0 * string.thermal_energy_J / (1.0e-11 * 1.0e-5 * q**2)
    assert np.allclose(string.mode_variance(n), closed, rtol=1e-12)


def test_taut_string_recovers_exactly_one_half_kT_per_quadratic_dof() -> None:
    """The equipartition identity, measured.  Target is the literal number 1, in units of kB*T/2."""
    string = TautStringSpectrum(tension_N=1.0e-11, length_m=1.0e-5, temperature_K=310.0)
    n = np.arange(1, 41)
    rng = np.random.default_rng(SEED)
    amplitudes = sample_gaussian_modes(string.mode_variance(n), n_samples=20_000, rng=rng)
    ratio = mean_energy_per_quadratic_dof(amplitudes, string.mode_stiffness(n), temperature_K=310.0)
    # Per-mode sampling error is sqrt(2/N) = 1.0%; 40 modes, so a 5-sigma per-mode band is 5%.
    assert np.all(np.abs(ratio - 1.0) < 0.05)
    # The grand mean is 40x tighter still.
    assert float(ratio.mean()) == pytest.approx(1.0, abs=5.0 * math.sqrt(2.0 / 20_000 / 40))


def test_equipartition_measurement_catches_a_factor_of_two_in_the_energy() -> None:
    """A wrong factor of two in the energy is the classic silent bug; kB*T/2 units make it loud."""
    string = TautStringSpectrum(tension_N=1.0e-11, length_m=1.0e-5, temperature_K=310.0)
    n = np.arange(1, 21)
    rng = np.random.default_rng(SEED)
    amplitudes = sample_gaussian_modes(string.mode_variance(n), n_samples=8_000, rng=rng)
    wrong_stiffness = 2.0 * string.mode_stiffness(n)
    ratio = mean_energy_per_quadratic_dof(amplitudes, wrong_stiffness, temperature_K=310.0)
    assert float(ratio.mean()) == pytest.approx(2.0, rel=0.02)


def test_string_rejects_invalid_mode_indices() -> None:
    string = TautStringSpectrum(tension_N=1.0e-11, length_m=1.0e-5)
    with pytest.raises(ValueError, match="integers >= 1"):
        string.mode_wavevector(np.array([0, 1]))
    with pytest.raises(ValueError, match="integers >= 1"):
        string.mode_wavevector(np.array([1.5, 2.0]))


# --------------------------------------------------------------------------------------------
# 4. POSITIVE CONTROL — the estimator recovers a truth it was not given
# --------------------------------------------------------------------------------------------


def test_positive_control_recovers_kappa_and_sigma_inside_their_stated_intervals() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 1.0e9, 60)
    n_samples = 4000
    _, measured = _sampled_spectrum(spectrum, q=q, n_samples=n_samples)
    fit = fit_helfrich_spectrum(
        q,
        measured,
        area_m2=PATCH_AREA_M2,
        temperature_K=PHYSIOLOGICAL_TEMPERATURE_K,
        n_samples=n_samples,
    )
    assert fit.converged
    assert fit.kappa_ci[0] <= ENGINE_KAPPA_J <= fit.kappa_ci[1]
    assert fit.sigma_ci[0] <= ENGINE_SIGMA_N_PER_M <= fit.sigma_ci[1]
    assert abs(fit.kappa_J - ENGINE_KAPPA_J) / ENGINE_KAPPA_J < 0.02
    assert abs(fit.sigma_N_per_m - ENGINE_SIGMA_N_PER_M) / ENGINE_SIGMA_N_PER_M < 0.02
    assert fit.kappa_resolved and fit.sigma_resolved
    assert fit.dof == 58
    assert fit.p_value > 0.01
    assert fit.decades == pytest.approx(4.0, rel=1e-9)


def test_positive_control_interval_coverage_is_close_to_nominal() -> None:
    """A confidence interval that never misses is not calibrated, it is wide.  Measure the rate."""
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e6, 1.0e9, 24)
    n_samples = 800
    rng = np.random.default_rng(4242)
    covered = 0
    trials = 120
    for _ in range(trials):
        amplitudes = sample_gaussian_modes(spectrum.mode_variance(q), n_samples=n_samples, rng=rng)
        fit = fit_helfrich_spectrum(
            q,
            measure_mode_variance(amplitudes),
            area_m2=PATCH_AREA_M2,
            temperature_K=PHYSIOLOGICAL_TEMPERATURE_K,
            n_samples=n_samples,
            confidence=0.95,
        )
        covered += fit.kappa_ci[0] <= ENGINE_KAPPA_J <= fit.kappa_ci[1]
    rate = covered / trials
    # Binomial three-sigma band around 0.95 for 120 trials is +-0.0597; declared before the run.
    assert 0.95 - 3.0 * math.sqrt(0.95 * 0.05 / trials) <= rate <= 1.0


def test_the_estimator_is_deterministic_given_the_seed() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e6, 1.0e9, 20)
    first = _sampled_spectrum(spectrum, q=q, n_samples=500, seed=11)[1]
    second = _sampled_spectrum(spectrum, q=q, n_samples=500, seed=11)[1]
    assert np.array_equal(first, second)


def test_measure_mode_variance_does_not_subtract_an_estimated_mean() -> None:
    """<h> is zero by construction; subtracting an estimate would bias the variance low by 1/N."""
    values = np.array([[1.0], [-1.0], [3.0]])
    assert measure_mode_variance(values)[0] == pytest.approx((1.0 + 1.0 + 9.0) / 3.0)


# --------------------------------------------------------------------------------------------
# 5. NEGATIVE CONTROL 1 — a wrong power law, and how much q span it takes to see it
# --------------------------------------------------------------------------------------------


def test_a_q_minus_two_spectrum_is_never_rejected_by_the_two_parameter_fit() -> None:
    """The measured structural fact behind :data:`HELFRICH_FAMILY_ABSORPTION`.

    ``q**-2`` IS a Helfrich spectrum (the ``kappa = 0`` member), so the goodness-of-fit cannot
    discriminate it from ``q**-4`` at any span.  This test exists so nobody reports fit quality as
    evidence of bending.
    """
    assert "never report fit quality" in HELFRICH_FAMILY_ABSORPTION
    rng = np.random.default_rng(SEED)
    n_samples = 1000
    for span in (1.0, 4.0):
        q = log_spaced_wavevectors(1.0e7, 1.0e7 * 10.0**span, int(20 * span) + 1)
        truth = power_law_variances(q, exponent=-2.0, reference_q=1.0e7, reference_variance=1.0e-18)
        amplitudes = sample_gaussian_modes(truth, n_samples=n_samples, rng=rng)
        fit = fit_helfrich_spectrum(
            q,
            measure_mode_variance(amplitudes),
            area_m2=PATCH_AREA_M2,
            temperature_K=PHYSIOLOGICAL_TEMPERATURE_K,
            n_samples=n_samples,
            confidence=0.99,
        )
        assert fit.p_value > 0.01, "the Helfrich family absorbs q^-2; chi-square must not reject it"
        # But the fit must NOT claim a bending modulus it cannot see.
        assert not fit.kappa_resolved
        assert fit.sigma_resolved


def test_the_exponent_test_rejects_the_wrong_power_law_over_a_sufficient_span() -> None:
    rng = np.random.default_rng(SEED)
    n_samples = 1000
    q = log_spaced_wavevectors(1.0e7, 1.0e9, 41)
    truth = power_law_variances(q, exponent=-2.0, reference_q=1.0e7, reference_variance=1.0e-18)
    amplitudes = sample_gaussian_modes(truth, n_samples=n_samples, rng=rng)
    exponent = fit_power_law_exponent(
        q, measure_mode_variance(amplitudes), n_samples=n_samples, confidence=0.95
    )
    assert exponent.exponent == pytest.approx(-2.0, abs=0.05)
    assert exponent.excludes(-4.0)
    assert not exponent.excludes(-2.0)
    assert exponent.separation_sigma(-4.0) > 50.0
    assert exponent.decades == pytest.approx(2.0, rel=1e-9)


def test_the_exponent_test_recovers_minus_four_on_a_bending_dominated_band() -> None:
    spectrum = HelfrichSpectrum(ENGINE_KAPPA_J, 0.0, PATCH_AREA_M2, PHYSIOLOGICAL_TEMPERATURE_K)
    q = log_spaced_wavevectors(1.0e8, 1.0e9, 30)
    n_samples = 3000
    _, measured = _sampled_spectrum(spectrum, q=q, n_samples=n_samples)
    exponent = fit_power_law_exponent(q, measured, n_samples=n_samples)
    assert exponent.exponent == pytest.approx(-4.0, abs=0.05)
    assert exponent.excludes(-2.0)


def test_predicted_decade_requirement_matches_the_measured_separation() -> None:
    """The design formula must be validated against sampling, not trusted.

    Prediction and measurement are compared as *separation in sigma at the predicted span*: at the
    span the formula says gives ``n_sigma``, the measured mean separation must be that, within the
    trial-to-trial scatter.
    """
    n_samples = 1000
    modes_per_decade = 20.0
    n_sigma = 5.0
    predicted = decades_required_for_exponent_separation(
        delta_exponent=2.0,
        n_samples=n_samples,
        modes_per_decade=modes_per_decade,
        n_sigma=n_sigma,
    )
    assert 0.01 < predicted < 1.0
    q = log_spaced_wavevectors(1.0e7, 1.0e7 * 10.0**predicted, 3)
    truth = power_law_variances(q, exponent=-2.0, reference_q=1.0e7, reference_variance=1.0e-18)
    rng = np.random.default_rng(SEED)
    separations = []
    for _ in range(200):
        amplitudes = sample_gaussian_modes(truth, n_samples=n_samples, rng=rng)
        fit = fit_power_law_exponent(
            q, measure_mode_variance(amplitudes), n_samples=n_samples, confidence=0.95
        )
        separations.append(fit.separation_sigma(-4.0))
    assert float(np.mean(separations)) == pytest.approx(n_sigma, rel=0.15)


def test_decade_requirement_shrinks_as_the_inverse_square_root_of_the_sample_count() -> None:
    """At a fixed mode count the requirement must scale exactly as ``1/sqrt(N)``.

    ``modes_per_decade`` is set tiny so both calls sit on the three-mode floor; otherwise the mode
    count changes with the span too and the scaling is only approximate, which would make this test
    a statement about the mode ladder rather than about the sampling law.
    """
    coarse = decades_required_for_exponent_separation(
        n_samples=100, modes_per_decade=1e-6, n_sigma=5.0
    )
    fine = decades_required_for_exponent_separation(
        n_samples=1600, modes_per_decade=1e-6, n_sigma=5.0
    )
    assert coarse / fine == pytest.approx(4.0, rel=1e-6)
    # The mode ladder makes the realistic case only approximately 1/sqrt(N).
    dense_coarse = decades_required_for_exponent_separation(
        n_samples=100, modes_per_decade=20.0, n_sigma=5.0
    )
    dense_fine = decades_required_for_exponent_separation(
        n_samples=1600, modes_per_decade=20.0, n_sigma=5.0
    )
    assert 3.0 < dense_coarse / dense_fine < 4.0


def test_decade_requirement_refuses_when_the_sample_count_is_the_binding_constraint() -> None:
    with pytest.raises(ValueError, match="no span below 100 decades"):
        decades_required_for_exponent_separation(
            delta_exponent=2.0, n_samples=2, modes_per_decade=1e-6, n_sigma=1000.0
        )


def test_measured_separation_ladder_reports_a_span_and_agrees_with_the_prediction() -> None:
    """The number the brief asks for: MEASURED decades of q needed to separate q^-2 from q^-4."""
    ladder = measure_exponent_separation_ladder(
        decades=(0.02, 0.05, 0.1, 0.2),
        modes_per_decade=20.0,
        n_samples=1000,
        n_trials=60,
        seed=7,
    )
    assert ladder.rejection_rate[0] < 1.0
    assert ladder.rejection_rate[-1] == 1.0
    measured = ladder.measured_decades
    assert measured is not None
    assert measured <= 0.1
    # Monotone: more span never lowers the mean separation.
    assert list(ladder.mean_separation_sigma) == sorted(ladder.mean_separation_sigma)
    # The 1.96-sigma prediction must sit below the span where every trial rejects.
    assert ladder.predicted_decades < measured
    record = ladder.to_dict()
    assert record["seed"] == 7
    assert record["measured_decades"] == measured
    assert record["evidence_source"] == EvidenceSource.ANALYTIC_ORACLE.value


def test_separation_ladder_validates_its_inputs() -> None:
    with pytest.raises(ValueError, match="ascending"):
        measure_exponent_separation_ladder(
            decades=(0.2, 0.1), modes_per_decade=20.0, n_samples=100, n_trials=2, seed=1
        )
    with pytest.raises(ValueError, match="must not be empty"):
        measure_exponent_separation_ladder(
            decades=(), modes_per_decade=20.0, n_samples=100, n_trials=2, seed=1
        )


def test_a_spectrum_outside_the_helfrich_family_is_rejected_by_goodness_of_fit() -> None:
    """``q**-3`` cannot be absorbed, so chi-square DOES see it — but needs far more span."""
    rng = np.random.default_rng(SEED)
    n_samples = 1000
    narrow = log_spaced_wavevectors(1.0e7, 1.0e7 * 10.0**0.2, 5)
    wide = log_spaced_wavevectors(1.0e7, 1.0e7 * 10.0**1.0, 21)
    outcomes = {}
    for name, q in (("narrow", narrow), ("wide", wide)):
        truth = power_law_variances(q, exponent=-3.0, reference_q=1.0e7, reference_variance=1.0e-18)
        amplitudes = sample_gaussian_modes(truth, n_samples=n_samples, rng=rng)
        fit = fit_helfrich_spectrum(
            q,
            measure_mode_variance(amplitudes),
            area_m2=PATCH_AREA_M2,
            temperature_K=PHYSIOLOGICAL_TEMPERATURE_K,
            n_samples=n_samples,
            confidence=0.99,
        )
        outcomes[name] = fit.p_value
    assert outcomes["narrow"] > 0.01, "0.2 decades cannot see a q^-3 spectrum"
    assert outcomes["wide"] < 0.01, "1.0 decade must reject a q^-3 spectrum"


# --------------------------------------------------------------------------------------------
# 6. NEGATIVE CONTROL 2 — equipartition violated by a factor
# --------------------------------------------------------------------------------------------


def test_a_free_two_parameter_fit_is_exactly_blind_to_a_global_amplitude_rescale() -> None:
    """The degeneracy, demonstrated numerically rather than asserted in a docstring."""
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 1.0e9, 40)
    n_samples = 2000
    _, measured = _sampled_spectrum(spectrum, q=q, n_samples=n_samples)
    factor = 1.44
    common = {
        "area_m2": PATCH_AREA_M2,
        "temperature_K": PHYSIOLOGICAL_TEMPERATURE_K,
        "n_samples": n_samples,
    }
    base = fit_helfrich_spectrum(q, measured, **common)
    scaled = fit_helfrich_spectrum(q, measured * factor, **common)
    assert scaled.chi_square == pytest.approx(base.chi_square, rel=1e-9)
    assert scaled.kappa_J / base.kappa_J == pytest.approx(1.0 / factor, rel=1e-6)
    assert scaled.sigma_N_per_m / base.sigma_N_per_m == pytest.approx(1.0 / factor, rel=1e-6)
    assert "independently determined" in EQUIPARTITION_SCALE_DEGENERACY


def test_equipartition_check_against_independent_moduli_catches_the_rescale() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 1.0e9, 40)
    n_samples = 2000
    _, measured = _sampled_spectrum(spectrum, q=q, n_samples=n_samples)

    clean = check_equipartition(q, measured, spectrum=spectrum, n_samples=n_samples)
    assert clean.consistent
    assert clean.ratio == pytest.approx(1.0, abs=0.02)
    assert clean.implied_temperature_K == pytest.approx(PHYSIOLOGICAL_TEMPERATURE_K, rel=0.02)

    for factor in (1.10, 1.44, 0.5):
        violated = check_equipartition(q, measured * factor, spectrum=spectrum, n_samples=n_samples)
        assert not violated.consistent, f"factor {factor} was absorbed"
        assert violated.ratio == pytest.approx(clean.ratio * factor, rel=1e-9)
    assert clean.degeneracy_note == EQUIPARTITION_SCALE_DEGENERACY


def test_equipartition_check_refuses_a_spectrum_that_is_not_a_spectrum_object() -> None:
    q = log_spaced_wavevectors(1.0e6, 1.0e9, 10)
    with pytest.raises(TypeError, match="determined independently"):
        check_equipartition(
            q,
            np.full(10, 1.0e-18),
            spectrum={"kappa": 1.0},  # type: ignore[arg-type]
            n_samples=100,
        )


# --------------------------------------------------------------------------------------------
# 7. NEGATIVE CONTROL 3 — correlated modes that a variance check cannot see
# --------------------------------------------------------------------------------------------


def test_correlated_modes_pass_every_variance_check_and_are_caught_only_by_independence() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 1.0e9, 12)
    n_samples = 2000
    rng = np.random.default_rng(SEED)
    amplitudes = correlated_mode_amplitudes(
        spectrum.mode_variance(q), n_samples=n_samples, rng=rng, rho=0.15
    )
    measured = measure_mode_variance(amplitudes)

    fit = fit_helfrich_spectrum(
        q,
        measured,
        area_m2=PATCH_AREA_M2,
        temperature_K=PHYSIOLOGICAL_TEMPERATURE_K,
        n_samples=n_samples,
        confidence=0.99,
    )
    assert fit.p_value > 0.01, "the variance-only check is supposed to be fooled here"
    assert abs(fit.kappa_J - ENGINE_KAPPA_J) / ENGINE_KAPPA_J < 0.10
    equipartition = check_equipartition(q, measured, spectrum=spectrum, n_samples=n_samples)
    assert equipartition.consistent, "equipartition is also supposed to be fooled here"

    independence = check_mode_independence(amplitudes, confidence=0.99)
    assert not independence.independent
    assert independence.max_abs_correlation > 0.10
    assert independence.min_p_value < independence.bonferroni_alpha
    assert independence.family_size == 12 * 11 // 2


def test_independent_modes_pass_the_independence_test() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 1.0e9, 12)
    amplitudes, _ = _sampled_spectrum(spectrum, q=q, n_samples=2000)
    report = check_mode_independence(amplitudes, confidence=0.99)
    assert report.independent
    assert report.max_abs_correlation < 0.15
    assert report.n_modes == 12


def test_independence_test_family_wise_error_is_measured_not_assumed() -> None:
    """Bonferroni over 66 pairs must not fire on independent data at anything like the nominal rate."""
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e6, 1.0e9, 12)
    rng = np.random.default_rng(909)
    false_alarms = 0
    trials = 60
    for _ in range(trials):
        amplitudes = sample_gaussian_modes(spectrum.mode_variance(q), n_samples=400, rng=rng)
        if not check_mode_independence(amplitudes, confidence=0.99).independent:
            false_alarms += 1
    assert false_alarms <= 2  # nominal family-wise alpha is 0.01 -> expected 0.6 in 60


def test_correlation_builder_and_independence_test_validate_inputs() -> None:
    variances = np.full(4, 1.0e-18)
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="strictly inside"):
        correlated_mode_amplitudes(variances, n_samples=10, rng=rng, rho=1.0)
    with pytest.raises(ValueError, match="symmetric"):
        sample_gaussian_modes(
            variances,
            n_samples=10,
            rng=rng,
            correlation=np.array(
                [
                    [1.0, 0.5, 0.0, 0.0],
                    [0.1, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0],
                ]
            ),
        )
    with pytest.raises(ValueError, match="at least two modes"):
        check_mode_independence(np.zeros((10, 1)))
    with pytest.raises(ValueError, match="five samples"):
        check_mode_independence(np.zeros((3, 4)))
    with pytest.raises(TypeError, match="numpy.random.Generator"):
        sample_gaussian_modes(variances, n_samples=10, rng=12345)  # type: ignore[arg-type]


# --------------------------------------------------------------------------------------------
# 8. NEGATIVE CONTROL 4 — a mesh cutoff is not an absence of physics
# --------------------------------------------------------------------------------------------


def test_nyquist_wavevector_is_pi_over_the_spacing() -> None:
    assert mesh_nyquist_wavevector(20.0e-9) == pytest.approx(math.pi / 20.0e-9, rel=1e-12)


def test_fitting_across_a_mesh_cutoff_inflates_kappa_and_the_trusted_range_prevents_it() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 3.0e8, 40)
    n_samples = 4000
    truncated = mesh_truncated_variances(q, spectrum.mode_variance(q), q_cutoff=1.2e8)
    rng = np.random.default_rng(SEED)
    amplitudes = sample_gaussian_modes(truncated, n_samples=n_samples, rng=rng)
    measured = measure_mode_variance(amplitudes)

    report = trusted_wavevector_range(
        q,
        measured,
        area_m2=PATCH_AREA_M2,
        temperature_K=PHYSIOLOGICAL_TEMPERATURE_K,
        n_samples=n_samples,
        confidence=0.99,
        nyquist_q=mesh_nyquist_wavevector(20.0e-9),
    )
    assert report.cutoff_detected
    assert report.n_modes_rejected > 0
    assert report.n_modes_trusted + report.n_modes_rejected == 40
    assert report.q_trust_max < 1.2e8
    assert report.nyquist_q == pytest.approx(math.pi / 20.0e-9)

    common = {
        "area_m2": PATCH_AREA_M2,
        "temperature_K": PHYSIOLOGICAL_TEMPERATURE_K,
        "n_samples": n_samples,
    }
    naive = fit_helfrich_spectrum(q, measured, **common)
    guarded = fit_helfrich_spectrum(
        q[: report.n_modes_trusted], measured[: report.n_modes_trusted], **common
    )
    # The measured damage AT THIS CUTOFF POSITION: fitting through it reports a stiffer membrane.
    # The multiple is deliberately not pinned here — it is a property of where the cutoff sits, and
    # the curve is measured by test_cutoff_damage_is_a_curve_in_cutoff_position_not_a_multiple.
    assert naive.kappa_J / ENGINE_KAPPA_J > 5.0
    assert guarded.kappa_J / ENGINE_KAPPA_J == pytest.approx(1.0, rel=0.10)


def test_cutoff_damage_is_a_curve_in_cutoff_position_not_a_multiple() -> None:
    """External review 2026-07-29, finding 9: ``33.63x`` characterises one cutoff, not the failure.

    Counterexample the reviewer supplied and this test re-measures: holding the whole design fixed and
    moving only ``q_cutoff``, the kappa ratio ran ``8e7 -> 1.9e5`` down to ``3e8 -> 1.13``.  A single
    multiple is therefore not transferable.  What IS transferable is asserted here: the curve is
    monotone in cutoff position, it is unbounded as the cutoff descends into the fitted band, it tends
    to 1 as it leaves, and the guarded fit recovers kappa at EVERY position on the ladder.  No exact
    multiple is pinned, because pinning one would re-commit the error being corrected.
    """
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 3.0e8, 40)
    ladder = measure_cutoff_damage_ladder(
        spectrum,
        q=q,
        q_cutoffs=(8.0e7, 1.0e8, 1.2e8, 1.5e8, 2.0e8, 3.0e8),
        n_samples=4000,
        seed=SEED,
        confidence=0.99,
    )
    ratios = list(ladder.naive_kappa_ratio)
    # Monotone: pushing the cutoff up the band never increases the damage.
    assert ratios == sorted(ratios, reverse=True)
    # Unbounded at one end, harmless at the other — spanning four orders of magnitude of "damage".
    assert ratios[0] > 1.0e4
    assert ratios[-1] < 1.5
    assert ratios[0] / ratios[-1] > 1.0e4
    # The transferable claim: the guard recovers kappa wherever the cutoff sits.
    for guarded in ladder.guarded_kappa_ratio:
        assert guarded is not None
        assert guarded == pytest.approx(1.0, rel=0.10)
    # Cutoff position is recorded dimensionlessly, so the number can be argued with.
    assert ladder.cutoff_over_fitted_q_max[2] == pytest.approx(0.4, rel=1e-12)
    assert ladder.decades_of_band_above_cutoff[-1] == pytest.approx(0.0, abs=1e-12)
    record = ladder.to_dict()
    assert record["evidence_source"] == EvidenceSource.ANALYTIC_ORACLE.value
    assert record["scope_note"] == CUTOFF_DAMAGE_IS_NOT_A_MULTIPLE
    assert "not a transferable multiple" in CUTOFF_DAMAGE_IS_NOT_A_MULTIPLE


def test_cutoff_damage_ladder_validates_its_inputs() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e6, 1.0e9, 10)
    with pytest.raises(ValueError, match="must not be empty"):
        measure_cutoff_damage_ladder(
            spectrum, q=q, q_cutoffs=(), n_samples=100, seed=1, confidence=0.99
        )
    bending_only = HelfrichSpectrum(0.0, ENGINE_SIGMA_N_PER_M, PATCH_AREA_M2, 310.0)
    with pytest.raises(ValueError, match="true kappa is zero"):
        measure_cutoff_damage_ladder(
            bending_only, q=q, q_cutoffs=(1.0e8,), n_samples=100, seed=1, confidence=0.99
        )
    with pytest.raises(ValueError, match="strictly ascending"):
        measure_cutoff_damage_ladder(
            spectrum,
            q=np.array([1.0e9, 1.0e8, 1.0e7]),
            q_cutoffs=(1.0e8,),
            n_samples=100,
            seed=1,
            confidence=0.99,
        )


def test_a_clean_spectrum_is_not_truncated_by_the_detector() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 1.0e9, 40)
    n_samples = 4000
    _, measured = _sampled_spectrum(spectrum, q=q, n_samples=n_samples)
    report = trusted_wavevector_range(
        q,
        measured,
        area_m2=PATCH_AREA_M2,
        temperature_K=PHYSIOLOGICAL_TEMPERATURE_K,
        n_samples=n_samples,
        confidence=0.99,
    )
    assert not report.cutoff_detected
    assert report.n_modes_trusted == 40
    assert report.q_trust_max == pytest.approx(1.0e9, rel=1e-9)
    assert "whole measured band" in report.reason


def test_a_spectrum_that_is_helfrich_nowhere_is_refused_rather_than_fitted() -> None:
    """If not even the lowest modes pass, the honest output is a refusal, not a smaller fit."""
    q = log_spaced_wavevectors(1.0e6, 1.0e9, 20)
    rng = np.random.default_rng(3)
    # A spectrum that RISES with q: no (kappa >= 0, sigma >= 0) reproduces it anywhere.
    truth = power_law_variances(q, exponent=+2.0, reference_q=1.0e6, reference_variance=1.0e-22)
    measured = measure_mode_variance(sample_gaussian_modes(truth, n_samples=3000, rng=rng))
    report = trusted_wavevector_range(
        q,
        measured,
        area_m2=PATCH_AREA_M2,
        temperature_K=PHYSIOLOGICAL_TEMPERATURE_K,
        n_samples=3000,
        confidence=0.99,
    )
    assert report.n_modes_trusted == 0
    assert report.cutoff_detected
    assert "must not be fitted" in report.reason


def test_trusted_range_requires_ascending_wavevectors() -> None:
    q = np.array([1.0e9, 1.0e8, 1.0e7])
    with pytest.raises(ValueError, match="strictly ascending"):
        trusted_wavevector_range(
            q,
            np.array([1e-18, 2e-18, 3e-18]),
            area_m2=PATCH_AREA_M2,
            temperature_K=310.0,
            n_samples=100,
        )


# --------------------------------------------------------------------------------------------
# 9. what S1 demands of a future native run
# --------------------------------------------------------------------------------------------


def test_patch_requirement_refuses_to_guess_the_boundary_condition() -> None:
    """External review 2026-07-29, finding 10: 40,401 mixed two sampling conventions.

    The superseded code took the patch side from the PERIODIC fundamental ``q_min = 2*pi/L`` (giving
    ``L/a = 2*q_max/q_min = 200``) and then added an endpoint for ``201**2 = 40,401``.  A periodic
    grid has no duplicated endpoint, so 40,401 is the node count of neither convention.  The
    correction is not a different number but a required declaration, so the omission must raise rather
    than default.
    """
    spectrum = _reference_spectrum()
    with pytest.raises(TypeError):
        native_patch_requirement(spectrum)  # type: ignore[call-arg]
    with pytest.raises(ValueError, match="must be declared as one of"):
        native_patch_requirement(spectrum, boundary_condition="whatever-the-engine-does")
    assert "will not choose one for the caller" in PATCH_BOUNDARY_CONDITION_IS_A_DECLARATION


def test_patch_requirement_converts_a_q_span_into_periodic_geometry() -> None:
    """The periodic reading: ``q_min = 2*pi/L``, and the far edge is not a separate node."""
    spectrum = _reference_spectrum()
    design = native_patch_requirement(
        spectrum,
        boundary_condition=PatchBoundaryCondition.PERIODIC,
        decades_below_crossover=1.0,
        decades_above_crossover=1.0,
    )
    assert design.decades == pytest.approx(2.0)
    assert design.q_min_per_m == pytest.approx(spectrum.crossover_wavevector_per_m / 10.0)
    assert design.q_max_per_m == pytest.approx(spectrum.crossover_wavevector_per_m * 10.0)
    assert design.patch_side_m == pytest.approx(2.0 * math.pi / design.q_min_per_m, rel=1e-12)
    assert design.mesh_spacing_m == pytest.approx(math.pi / design.q_max_per_m, rel=1e-12)
    # L/a = 2 q_max / q_min = 200 exactly, and a torus has 200 unique nodes, not 201.
    assert design.patch_side_m / design.mesh_spacing_m == pytest.approx(200.0, rel=1e-12)
    assert design.nodes_per_side == 200
    assert design.total_nodes == 40_000
    # Every node of a torus is free: no pinned boundary rows.
    assert design.free_nodes_per_side == 200
    assert design.total_free_nodes == 40_000
    # The measured requirement for the engine's own constants, stated so it can be argued with.
    assert design.patch_side_m * 1e6 == pytest.approx(5.717, rel=1e-3)
    assert design.mesh_spacing_m * 1e9 == pytest.approx(28.59, rel=1e-3)
    record = design.to_dict()
    assert record["boundary_condition"] == "periodic"
    assert record["fundamental_relation"].startswith("q_min = 2*pi/L")
    assert record["evidence_source"] == EvidenceSource.ANALYTIC_ORACLE.value


def test_patch_requirement_converts_a_q_span_into_fixed_edge_geometry() -> None:
    """The fixed-edge reading gives a HALF-WIDTH patch and a quarter of the nodes.

    This is why the boundary condition cannot be an implicit convention: it moves the patch side, not
    only the endpoint.  The taut string in this same module already uses ``q_n = n*pi/L``, so both
    conventions were live in one file.
    """
    spectrum = _reference_spectrum()
    design = native_patch_requirement(
        spectrum,
        boundary_condition="fixed-edge",
        decades_below_crossover=1.0,
        decades_above_crossover=1.0,
    )
    assert design.patch_side_m == pytest.approx(math.pi / design.q_min_per_m, rel=1e-12)
    assert design.mesh_spacing_m == pytest.approx(math.pi / design.q_max_per_m, rel=1e-12)
    # L/a = q_max / q_min = 100 exactly; endpoint-inclusive, so 101 nodes with 2 rows pinned.
    assert design.patch_side_m / design.mesh_spacing_m == pytest.approx(100.0, rel=1e-12)
    assert design.nodes_per_side == 101
    assert design.total_nodes == 10_201
    assert design.free_nodes_per_side == 99
    assert design.total_free_nodes == 9_801
    assert design.patch_side_m * 1e6 == pytest.approx(2.859, rel=1e-3)
    assert design.mesh_spacing_m * 1e9 == pytest.approx(28.59, rel=1e-3)
    assert design.to_dict()["fundamental_relation"].startswith("q_min = pi/L")


def test_the_two_boundary_conventions_disagree_and_neither_gives_the_retracted_number() -> None:
    """40,401 is not a rounding of either answer; it is the mixture, and it is retracted."""
    spectrum = _reference_spectrum()
    periodic = native_patch_requirement(spectrum, boundary_condition="periodic")
    fixed = native_patch_requirement(spectrum, boundary_condition="fixed-edge")
    assert periodic.total_nodes == 40_000
    assert fixed.total_nodes == 10_201
    assert 40_401 not in (periodic.total_nodes, fixed.total_nodes)
    # Same mesh spacing, different patch side: the side is the boundary-condition-dependent half.
    assert fixed.mesh_spacing_m == pytest.approx(periodic.mesh_spacing_m, rel=1e-12)
    assert periodic.patch_side_m / fixed.patch_side_m == pytest.approx(2.0, rel=1e-12)
    # Periodic is RECOMMENDED with an argument, and the recommendation is flagged as a PI decision.
    assert "PI decision" in PATCH_PERIODIC_RATIONALE
    assert "did not read the engine" in PATCH_PERIODIC_RATIONALE


def test_patch_requirement_refuses_a_spectrum_with_no_two_regime_band() -> None:
    bending_only = HelfrichSpectrum(ENGINE_KAPPA_J, 0.0, PATCH_AREA_M2, 310.0)
    tension_only = HelfrichSpectrum(0.0, ENGINE_SIGMA_N_PER_M, PATCH_AREA_M2, 310.0)
    for spectrum in (bending_only, tension_only):
        for condition in PatchBoundaryCondition:
            with pytest.raises(ValueError, match="no finite crossover"):
                native_patch_requirement(spectrum, boundary_condition=condition)


# --------------------------------------------------------------------------------------------
# 10. adjudication and failure routing
# --------------------------------------------------------------------------------------------


def test_a_clean_bending_dominated_slice_passes_every_clause() -> None:
    spectrum = HelfrichSpectrum(ENGINE_KAPPA_J, 0.0, PATCH_AREA_M2, PHYSIOLOGICAL_TEMPERATURE_K)
    q = log_spaced_wavevectors(1.0e8, 1.0e9, 30)
    amplitudes, _ = _sampled_spectrum(spectrum, q=q, n_samples=3000)
    verdict = adjudicate_fluctuation_spectrum(q, amplitudes, reference=spectrum)
    assert verdict.passed
    assert verdict.failures == ()
    assert verdict.expansions == ()
    assert verdict.exponent is not None
    assert verdict.exponent.exponent == pytest.approx(-4.0, abs=0.05)
    assert verdict.fit.kappa_J / ENGINE_KAPPA_J == pytest.approx(1.0, rel=0.05)


def test_a_clean_full_band_helfrich_slice_passes_every_clause() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 1.0e9, 60)
    amplitudes, _ = _sampled_spectrum(spectrum, q=q, n_samples=4000)
    verdict = adjudicate_fluctuation_spectrum(q, amplitudes, reference=spectrum)
    assert verdict.passed, verdict.failures
    assert verdict.trusted_range.n_modes_trusted == 60
    assert verdict.independence is not None
    assert verdict.equipartition is not None


def test_the_verdict_routes_an_equipartition_violation_to_the_boundary_law_project() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 1.0e9, 40)
    rng = np.random.default_rng(SEED)
    amplitudes = sample_gaussian_modes(spectrum.mode_variance(q), n_samples=3000, rng=rng)
    verdict = adjudicate_fluctuation_spectrum(q, amplitudes * math.sqrt(1.5), reference=spectrum)
    assert not verdict.passed
    assert FluctuationFailure.EQUIPARTITION_VIOLATED in verdict.failures
    assert any("stochastic boundary law" in text for text in verdict.expansions)


def test_the_verdict_routes_correlated_modes_to_the_boundary_law_project() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 1.0e9, 12)
    rng = np.random.default_rng(SEED)
    amplitudes = correlated_mode_amplitudes(
        spectrum.mode_variance(q), n_samples=2000, rng=rng, rho=0.15
    )
    verdict = adjudicate_fluctuation_spectrum(q, amplitudes, reference=spectrum, min_modes=5)
    assert not verdict.passed
    assert FluctuationFailure.MODES_CORRELATED in verdict.failures
    assert any("correlated-noise" in text for text in verdict.expansions)


def test_the_verdict_routes_a_mesh_cutoff_to_the_numerical_invariance_battery() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 3.0e8, 40)
    truncated = mesh_truncated_variances(q, spectrum.mode_variance(q), q_cutoff=1.2e8)
    rng = np.random.default_rng(SEED)
    amplitudes = sample_gaussian_modes(truncated, n_samples=4000, rng=rng)
    verdict = adjudicate_fluctuation_spectrum(q, amplitudes, reference=spectrum)
    assert not verdict.passed
    assert FluctuationFailure.DISCRETISATION_CUTOFF in verdict.failures
    assert any("numerical-invariance battery (S0)" in text for text in verdict.expansions)
    assert verdict.trusted_range.q_trust_max < 1.2e8


def test_the_verdict_routes_a_wrong_exponent_to_the_boundary_law_project() -> None:
    """Bending-dominated band, data drawn as ``q**-2``: the exponent clause must fire."""
    spectrum = HelfrichSpectrum(ENGINE_KAPPA_J, 0.0, PATCH_AREA_M2, PHYSIOLOGICAL_TEMPERATURE_K)
    q = log_spaced_wavevectors(1.0e8, 1.0e9, 30)
    reference_variance = float(spectrum.mode_variance(np.array([1.0e8]))[0])
    truth = power_law_variances(
        q, exponent=-2.0, reference_q=1.0e8, reference_variance=reference_variance
    )
    rng = np.random.default_rng(SEED)
    amplitudes = sample_gaussian_modes(truth, n_samples=3000, rng=rng)
    verdict = adjudicate_fluctuation_spectrum(q, amplitudes, reference=spectrum)
    assert not verdict.passed
    assert FluctuationFailure.WRONG_EXPONENT in verdict.failures
    assert FAILURE_EXPANSIONS[FluctuationFailure.WRONG_EXPONENT] in verdict.expansions


def test_every_failure_mode_has_a_declared_expansion_route() -> None:
    assert set(FAILURE_EXPANSIONS) == set(FluctuationFailure)
    for failure, route in FAILURE_EXPANSIONS.items():
        assert route.strip(), failure


def test_the_verdict_is_conjunctive_not_a_score() -> None:
    """One clause failing fails the verdict; there is no "mostly matching" spectrum."""
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 1.0e9, 40)
    rng = np.random.default_rng(SEED)
    amplitudes = sample_gaussian_modes(spectrum.mode_variance(q), n_samples=3000, rng=rng)
    verdict = adjudicate_fluctuation_spectrum(q, amplitudes * math.sqrt(1.5), reference=spectrum)
    assert verdict.trusted_range.n_modes_trusted == 40  # the shape is fine
    assert verdict.fit.p_value > 0.01  # and the fit is excellent
    assert not verdict.passed  # and it still fails


def test_adjudication_validates_its_inputs() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e6, 1.0e9, 10)
    with pytest.raises(TypeError, match="HelfrichSpectrum"):
        adjudicate_fluctuation_spectrum(q, np.zeros((100, 10)), reference=None)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="matching q"):
        adjudicate_fluctuation_spectrum(q, np.zeros((100, 7)), reference=spectrum)


# --------------------------------------------------------------------------------------------
# 11. pre-registration card
# --------------------------------------------------------------------------------------------


def test_s1_card_predeclares_all_four_negative_controls_and_the_expansion_route() -> None:
    card = s1_experiment_card((MANIFEST_HASH,))
    assert isinstance(card, SandboxExperimentCard)
    assert card.experiment_id == "S1-passive-fluctuation-canon"
    assert not card.wet_lab
    joined_negative = " | ".join(card.negative_controls)
    for expected in ("wrong power law", "equipartition", "correlated modes", "discretisation"):
        assert expected in joined_negative
    assert len(card.negative_controls) == 4
    assert any("kB*T/2" in text for text in card.positive_controls)
    assert any("stochastic boundary law" in text for text in card.failure_expansions)
    assert any("q^-2" in text for text in card.adversarial_controls)
    assert len(card.held_out_interventions) == 3
    assert "mode by mode" in card.question.lower()


def test_s1_card_hash_is_stable_and_content_addressed() -> None:
    first = s1_experiment_card((MANIFEST_HASH,))
    second = s1_experiment_card((MANIFEST_HASH,))
    other = s1_experiment_card(("c" * 64,))
    assert first.card_hash == second.card_hash
    assert first.card_hash != other.card_hash
    assert first.to_dict()["falsifier"] == first.falsifier


# --------------------------------------------------------------------------------------------
# 12. input validation on the estimators themselves
# --------------------------------------------------------------------------------------------


def test_fit_refuses_inputs_that_cannot_support_a_two_parameter_goodness_of_fit() -> None:
    q = log_spaced_wavevectors(1.0e6, 1.0e9, 2)
    with pytest.raises(ValueError, match="at least three modes"):
        fit_helfrich_spectrum(
            q, np.array([1e-18, 1e-19]), area_m2=PATCH_AREA_M2, temperature_K=310.0, n_samples=100
        )
    q3 = log_spaced_wavevectors(1.0e6, 1.0e9, 3)
    with pytest.raises(ValueError, match="one entry per wavevector"):
        fit_helfrich_spectrum(
            q3, np.array([1e-18, 1e-19]), area_m2=PATCH_AREA_M2, temperature_K=310.0, n_samples=100
        )
    with pytest.raises(ValueError, match="strictly positive"):
        fit_helfrich_spectrum(
            q3,
            np.array([1e-18, 0.0, 1e-20]),
            area_m2=PATCH_AREA_M2,
            temperature_K=310.0,
            n_samples=100,
        )
    with pytest.raises(ValueError, match="strictly in"):
        fit_helfrich_spectrum(
            q3,
            np.array([1e-18, 1e-19, 1e-20]),
            area_m2=PATCH_AREA_M2,
            temperature_K=310.0,
            n_samples=100,
            confidence=1.0,
        )


def test_exponent_fit_refuses_a_degenerate_wavevector_range() -> None:
    with pytest.raises(ValueError, match="at least three modes"):
        fit_power_law_exponent(np.array([1e8, 2e8]), np.array([1e-18, 1e-19]), n_samples=100)
    with pytest.raises(ValueError, match="span no range"):
        fit_power_law_exponent(
            np.array([1e8, 1e8, 1e8]), np.array([1e-18, 1e-18, 1e-18]), n_samples=100
        )


def test_log_spaced_wavevectors_validates_its_range() -> None:
    with pytest.raises(ValueError, match="q_max must exceed q_min"):
        log_spaced_wavevectors(1.0e9, 1.0e6, 10)
    with pytest.raises(ValueError, match="positive"):
        log_spaced_wavevectors(0.0, 1.0e9, 10)
    q = log_spaced_wavevectors(1.0e6, 1.0e9, 4)
    assert q.size == 4
    assert np.allclose(np.diff(np.log10(q)), 1.0)


def test_chi_square_survival_matches_known_values() -> None:
    assert chi_square_survival(0.0, 1) == pytest.approx(1.0)
    # median of chi2_1 is 0.4549; of chi2_2 is 2*ln2 = 1.3863
    assert chi_square_survival(0.454936, 1) == pytest.approx(0.5, abs=1e-5)
    assert chi_square_survival(2.0 * math.log(2.0), 2) == pytest.approx(0.5, abs=1e-12)
    assert chi_square_survival(1e4, 3) == pytest.approx(0.0, abs=1e-12)


def test_mesh_truncation_helper_leaves_low_q_untouched() -> None:
    spectrum = _reference_spectrum()
    q = log_spaced_wavevectors(1.0e5, 3.0e8, 40)
    truth = spectrum.mode_variance(q)
    truncated = mesh_truncated_variances(q, truth, q_cutoff=1.2e8)
    low = q < 1.0e7
    assert np.allclose(truncated[low], truth[low], rtol=1e-6)
    assert truncated[-1] < 0.5 * truth[-1]
    with pytest.raises(ValueError, match="one entry per wavevector"):
        mesh_truncated_variances(q, truth[:-1], q_cutoff=1.2e8)


def test_power_law_helper_passes_through_its_reference_point() -> None:
    q = np.array([1.0e7, 1.0e8])
    values = power_law_variances(q, exponent=-2.0, reference_q=1.0e7, reference_variance=3.0e-18)
    assert values[0] == pytest.approx(3.0e-18, rel=1e-12)
    assert values[1] == pytest.approx(3.0e-20, rel=1e-12)


def test_energy_helper_validates_shapes() -> None:
    with pytest.raises(ValueError, match="one entry per mode column"):
        mean_energy_per_quadratic_dof(np.zeros((5, 3)), np.ones(4), temperature_K=310.0)
    with pytest.raises(ValueError, match="non-empty"):
        mean_energy_per_quadratic_dof(np.zeros(5), np.ones(5), temperature_K=310.0)
    with pytest.raises(ValueError, match="finite"):
        mean_energy_per_quadratic_dof(np.full((5, 2), np.nan), np.ones(2), temperature_K=310.0)

"""CPU oracle tests for sandbox S11, the image inverse challenge.

Every image here is synthetic and every truth is known, so what is being validated is the
IDENTIFIABILITY MACHINERY, not a cell.  The self-scoring warning in the module docstring applies to
this file too: a closed loop scored inside one renderer family is a structure test, and the only
part of it that speaks to the outside world is the renderer-mismatch control.

Numbers asserted below were MEASURED first and then written down with margin.  Where an assertion
carries a threshold, the threshold is a contract declared here (seed, draw count, information floor,
lift ratio) and the measured value sits far from it, so the verdict does not turn on the exact
figure.  Nothing was re-thresholded after seeing a result; where the measurement was unflattering —
the coverage deficit in :func:`test_matched_coverage_is_close_to_nominal_but_not_binomially_clean` —
the unflattering number is what the test asserts.
"""

from __future__ import annotations

import math
import subprocess
import sys
import textwrap
from functools import cache
from pathlib import Path

import numpy as np
import pytest

from aleph.virtual_cell.contracts import SandboxExperimentCard
from aleph.virtual_cell.optics import LineEmitter, PointEmitter
from aleph.virtual_cell.sandbox_inverse import (
    FLUX_DEGENERACY_DIRECTION,
    JACOBIAN_LOG_STEPS,
    PARAMETER_NAMES,
    REFERENCE_PRIOR,
    REFERENCE_SCENE,
    REFERENCE_TRUTH,
    S11_SOURCE_IDS,
    SANDBOX_S11_CARD,
    SINGLE_PLANE_DESIGN,
    SUPERRES_CHANNEL_DESIGN,
    TWO_PLANE_DESIGN,
    GaussianLogPrior,
    ImagingScene,
    InverseParameters,
    ModalityVerdict,
    ObservationDesign,
    RendererFamily,
    UnidentifiableDirectionError,
    analytic_defocus_null_direction,
    expected_photons,
    fisher_information,
    identifiability_report,
    infer_posterior,
    measure_length_scale_degeneracy,
    measure_modality_contrast,
    measure_renderer_mismatch,
    model_jacobian,
    register_sandbox_sources,
    run_coverage_campaign,
    sample_observation,
)
from aleph.virtual_cell.source_registry import ExclusionReason, SourceRegistry
from aleph.virtual_cell.surrogate import SIGMA_GATE

# Declared BEFORE any campaign runs.  One seed and one draw count for the whole file, so every
# coverage number below is reproducible by re-running this file and nothing else.
SEED = 20260729
NOMINAL_LEVEL = 0.95
CALIBRATION_DRAWS = 240
MISMATCH_DRAWS = 240

# A direction counts as MEASURED when the image carries at least as much information about it as the
# broadest prior already did: 1 / 0.45^2 for the widest declared prior sigma.  This is a contract
# derived from the prior, not a number tuned to a spectrum.
PRIOR_INFORMATION_FLOOR = 1.0 / 0.45**2

# How much of a parameter's unit vector must lie in the measured subspace before a point estimate is
# allowed.  Nine tenths is a conventional "almost entirely"; the measured fractions are 0.44-0.56 for
# the degenerate parameters and 1.0000 for the resolved ones, so no verdict here is near the line.
MIN_IDENTIFIABLE_FRACTION = 0.9

# A modality must raise the information along the probed direction by an order of magnitude to count
# as having lifted it.  Measured ratios are 2.3x and 41.8x, which straddle this by wide margins.
LIFT_RATIO_REQUIRED = 10.0


@cache
def _single_plane_report():
    return identifiability_report(
        REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN, label="single"
    )


@cache
def _matched_campaign(sigma_scale: float):
    return run_coverage_campaign(
        REFERENCE_SCENE,
        SINGLE_PLANE_DESIGN,
        REFERENCE_PRIOR,
        n_draws=CALIBRATION_DRAWS,
        nominal_level=NOMINAL_LEVEL,
        seed=SEED,
        label=f"matched-sigma-scale-{sigma_scale}",
        sigma_scale=sigma_scale,
    )


@cache
def _mismatch(design_id: str):
    design = {
        SINGLE_PLANE_DESIGN.design_id: SINGLE_PLANE_DESIGN,
        TWO_PLANE_DESIGN.design_id: TWO_PLANE_DESIGN,
    }[design_id]
    return measure_renderer_mismatch(
        REFERENCE_SCENE,
        design,
        REFERENCE_PRIOR,
        control_id=f"axial-miscalibration-15pct-{design_id}",
        n_draws=MISMATCH_DRAWS,
        nominal_level=NOMINAL_LEVEL,
        seed=SEED,
        truth_family=RendererFamily.PIXEL_INTEGRATED_LINE,
        rayleigh_scale=1.15,
    )


# ------------------------------------------------------------------------------------------------
# 0. runtime contract
# ------------------------------------------------------------------------------------------------


def test_import_does_not_pull_in_a_gpu_runtime() -> None:
    """Importing this sandbox must not load Warp, torch, or JAX.

    Asserted in a SUBPROCESS on purpose: ``"warp" not in sys.modules`` inside the pytest session is a
    statement about the whole session, so it passes in isolation and fails under the full suite.
    """
    script = textwrap.dedent(
        """
        import sys
        from aleph.virtual_cell.sandbox_inverse import identifiability_report, infer_posterior

        for blocked in ("warp", "torch", "jax"):
            assert blocked not in sys.modules, f"sandbox_inverse import pulled in {blocked}"
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


# ------------------------------------------------------------------------------------------------
# 1. Sanity Gate: dimensions, boundary, conservation, numerics, sign sense
# ------------------------------------------------------------------------------------------------


def test_frame_total_is_the_declared_photon_budget() -> None:
    """Dimensional gate: rate x length x exposure is photons, recovered to the truncation loss."""
    frame = expected_photons(REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN)
    background = (
        REFERENCE_TRUTH.background_photons * REFERENCE_SCENE.height_px * (REFERENCE_SCENE.width_px)
    )
    signal = float(np.sum(frame)) - background
    budget = REFERENCE_TRUTH.total_photons * REFERENCE_SCENE.exposure_s
    assert budget == pytest.approx(2000.0)
    # Measured loss is 7.8e-7 relative, bounded above by the 4-sigma PSF truncation.
    assert signal == pytest.approx(budget, rel=1e-5)


def test_photon_budget_is_invariant_along_the_flux_degeneracy() -> None:
    """Sign sense: the declared degenerate direction has opposite signs and preserves the budget."""
    assert FLUX_DEGENERACY_DIRECTION[0] * FLUX_DEGENERACY_DIRECTION[1] < 0.0
    assert float(np.linalg.norm(FLUX_DEGENERACY_DIRECTION)) == pytest.approx(1.0)
    moved = InverseParameters.from_vector(
        REFERENCE_TRUTH.to_vector() + 0.7 * FLUX_DEGENERACY_DIRECTION
    )
    assert moved.length_um != pytest.approx(REFERENCE_TRUTH.length_um)
    assert moved.total_photons == pytest.approx(REFERENCE_TRUTH.total_photons)


def test_defocus_null_degenerates_to_pure_axial_at_the_focal_plane() -> None:
    """Boundary gate: in focus, ``d ln s / d z`` vanishes and the null direction is z alone."""
    in_focus = InverseParameters(
        length_um=REFERENCE_TRUTH.length_um,
        linear_density=REFERENCE_TRUTH.linear_density,
        psf_sigma_um=REFERENCE_TRUTH.psf_sigma_um,
        background_photons=REFERENCE_TRUTH.background_photons,
        z_um=0.0,
    )
    direction = analytic_defocus_null_direction(in_focus, REFERENCE_SCENE, 0.0)
    assert direction[2] == pytest.approx(0.0, abs=1e-15)
    assert abs(direction[4]) == pytest.approx(1.0)
    # The reference truth is deliberately OFF focus so the measured degeneracy is the mixture.
    off_focus = analytic_defocus_null_direction(REFERENCE_TRUTH, REFERENCE_SCENE, 0.0)
    assert abs(off_focus[2]) > 0.4


def test_adding_a_modality_can_only_add_information() -> None:
    """Conservation gate: augmented-minus-baseline Fisher is positive semidefinite."""
    baseline = fisher_information(REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN)
    scale = float(np.max(np.abs(baseline)))
    for design in (TWO_PLANE_DESIGN, SUPERRES_CHANNEL_DESIGN):
        augmented = fisher_information(REFERENCE_TRUTH, REFERENCE_SCENE, design)
        smallest = float(np.min(np.linalg.eigvalsh(augmented - baseline)))
        # Measured: -1.1e-12 and +1.7e-12 against a matrix scale of 3.5e4.
        assert smallest > -1e-9 * scale


def test_halving_the_step_moves_only_the_null_eigenvalue() -> None:
    """Numerical gate: this is how a null direction is told from a merely small one."""
    coarse = np.linalg.eigvalsh(
        fisher_information(REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN)
    )
    fine = np.linalg.eigvalsh(
        fisher_information(
            REFERENCE_TRUTH,
            REFERENCE_SCENE,
            SINGLE_PLANE_DESIGN,
            steps=tuple(0.5 * step for step in JACOBIAN_LOG_STEPS),
        )
    )
    determined = np.abs(fine[1:] - coarse[1:]) / np.abs(coarse[1:])
    # Measured: at most 3.1e-8 on the four determined eigenvalues.
    assert float(np.max(determined)) < 1e-6
    null_change = abs(fine[0] - coarse[0]) / abs(coarse[0])
    # Measured: 0.66 — the smallest eigenvalue is pure differentiation noise.
    assert null_change > 0.1


def test_jacobian_is_a_derivative_of_the_forward_model() -> None:
    """A one-sided difference at a much coarser step must agree with the central difference."""
    jacobian = model_jacobian(REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN)
    assert jacobian.shape == (
        REFERENCE_SCENE.height_px * REFERENCE_SCENE.width_px,
        len(PARAMETER_NAMES),
    )
    base = expected_photons(REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN).reshape(-1)
    for index in range(len(PARAMETER_NAMES)):
        shifted = REFERENCE_TRUTH.to_vector()
        delta = 1e-3
        shifted[index] += delta
        forward = expected_photons(
            InverseParameters.from_vector(shifted), REFERENCE_SCENE, SINGLE_PLANE_DESIGN
        ).reshape(-1)
        numerical = (forward - base) / delta
        reference = jacobian[:, index]
        scale = max(float(np.max(np.abs(reference))), 1e-12)
        assert float(np.max(np.abs(numerical - reference))) < 1e-2 * scale


# ------------------------------------------------------------------------------------------------
# 2. the degeneracies, measured
# ------------------------------------------------------------------------------------------------


def test_single_frame_has_an_exact_defocus_degeneracy() -> None:
    """The worst eigendirection IS the closed-form ``(log s0, z)`` null, to machine precision."""
    report = _single_plane_report()
    assert report.condition_number > 1e12
    worst = np.asarray(report.eigenvectors[:, 0], dtype=np.float64)
    analytic = analytic_defocus_null_direction(REFERENCE_TRUTH, REFERENCE_SCENE, 0.0)
    # Measured overlap: 1.0000.
    assert abs(float(worst @ analytic)) > 0.9999
    assert report.directional_information(analytic) < 1e-9


def test_length_and_density_are_the_second_degeneracy() -> None:
    """The declared flux direction is recovered as an eigenvector, not assumed to be one."""
    report = _single_plane_report()
    second = np.asarray(report.eigenvectors[:, 1], dtype=np.float64)
    # Measured overlap: 0.9999 at an eigenvalue 5.0e-5 of the largest.
    assert abs(float(second @ FLUX_DEGENERACY_DIRECTION)) > 0.999
    relative = float(report.eigenvalues[1]) / report.largest_eigenvalue
    assert relative < 1e-3


def test_near_null_directions_name_their_dominant_parameters() -> None:
    nulls = _single_plane_report().near_null_directions(eigenvalue_floor=PRIOR_INFORMATION_FLOOR)
    assert len(nulls) == 2
    assert set(nulls[0].dominant_parameters) == {"log_psf_sigma_um", "z_um"}
    assert set(nulls[1].dominant_parameters) == {"log_length_um", "log_linear_density"}
    assert nulls[0].candidate_overlaps["analytic-defocus"] > 0.999
    assert nulls[1].candidate_overlaps["flux"] > 0.999
    assert nulls[0].to_dict()["parameter_names"] == list(PARAMETER_NAMES)


def test_identifiable_fractions_split_the_degenerate_pairs_in_half() -> None:
    """Two parameters sharing one measured combination each keep about half of their unit vector."""
    report = _single_plane_report()
    fractions = {
        name: report.identifiable_fraction(name, eigenvalue_floor=PRIOR_INFORMATION_FLOOR)
        for name in PARAMETER_NAMES
    }
    # Measured: 0.5005, 0.4997, 0.4419, 1.0000, 0.5579.
    assert fractions["log_background_photons"] == pytest.approx(1.0, abs=1e-6)
    for name in ("log_length_um", "log_linear_density", "log_psf_sigma_um", "z_um"):
        assert 0.3 < fractions[name] < 0.7


def test_flux_degeneracy_weakens_with_resolvability_but_never_disappears() -> None:
    """Sweep ``L / s0`` at a fixed photon budget: the degeneracy is a resolvability statement."""
    records = measure_length_scale_degeneracy(
        REFERENCE_SCENE,
        SINGLE_PLANE_DESIGN,
        base=REFERENCE_TRUTH,
        length_over_sigma=(0.25, 0.5, 1.0, 2.0, 4.0, 8.0),
    )
    relative = [record["flux_information_relative"] for record in records]
    assert all(
        later > earlier for earlier, later in zip(relative[:-1], relative[1:], strict=True)
    ), "flux information must grow monotonically with resolvability"
    # Measured: 4.3e-7 at L/s = 0.25 and 2.8e-2 at L/s = 8 — never three percent.
    assert relative[0] < 1e-5
    assert relative[-1] < 0.05
    assert all(record["second_direction_flux_overlap"] > 0.9 for record in records)


# ------------------------------------------------------------------------------------------------
# 3. A5 enforced structurally: no point estimate along an unmeasured direction
# ------------------------------------------------------------------------------------------------


@cache
def _reference_posterior(design_id: str):
    design = {
        SINGLE_PLANE_DESIGN.design_id: SINGLE_PLANE_DESIGN,
        TWO_PLANE_DESIGN.design_id: TWO_PLANE_DESIGN,
        SUPERRES_CHANNEL_DESIGN.design_id: SUPERRES_CHANNEL_DESIGN,
    }[design_id]
    expected = expected_photons(REFERENCE_TRUTH, REFERENCE_SCENE, design)
    observed = sample_observation(expected, REFERENCE_SCENE, seed=SEED)
    return design, infer_posterior(observed, REFERENCE_SCENE, design, REFERENCE_PRIOR)


def test_point_estimate_is_refused_for_every_degenerate_parameter() -> None:
    """The master plan's A5 verdict, as an exception rather than a paragraph."""
    _, posterior = _reference_posterior(SINGLE_PLANE_DESIGN.design_id)
    for name in ("log_length_um", "log_linear_density", "log_psf_sigma_um", "z_um"):
        with pytest.raises(UnidentifiableDirectionError, match="not identifiable"):
            posterior.point_estimate(
                name,
                eigenvalue_floor=PRIOR_INFORMATION_FLOOR,
                min_identifiable_fraction=MIN_IDENTIFIABLE_FRACTION,
            )
    background = posterior.point_estimate(
        "log_background_photons",
        eigenvalue_floor=PRIOR_INFORMATION_FLOOR,
        min_identifiable_fraction=MIN_IDENTIFIABLE_FRACTION,
    )
    assert math.isfinite(background)


def test_an_interval_is_always_available_even_where_a_point_estimate_is_not() -> None:
    """An honest wide interval is the correct answer; silence would not be."""
    _, posterior = _reference_posterior(SINGLE_PLANE_DESIGN.design_id)
    low, high = posterior.credible_interval("z_um", level=NOMINAL_LEVEL)
    assert low < high
    prior_sigma = float(REFERENCE_PRIOR.sigma[PARAMETER_NAMES.index("z_um")])
    # The image measured nothing there, so the interval must be close to the prior's own width.
    assert posterior.marginal_sigma("z_um") > 0.5 * prior_sigma


def test_estimate_along_a_direction_refuses_below_the_declared_floor() -> None:
    _, posterior = _reference_posterior(SINGLE_PLANE_DESIGN.design_id)
    analytic = analytic_defocus_null_direction(REFERENCE_TRUTH, REFERENCE_SCENE, 0.0)
    with pytest.raises(UnidentifiableDirectionError, match="Fisher information"):
        posterior.estimate_along(analytic, eigenvalue_floor=PRIOR_INFORMATION_FLOOR)
    total_flux = np.array([1.0, 1.0, 0.0, 0.0, 0.0])
    assert math.isfinite(
        posterior.estimate_along(total_flux, eigenvalue_floor=PRIOR_INFORMATION_FLOOR)
    )


def test_the_superres_channel_makes_every_point_estimate_admissible() -> None:
    """The refusal is a property of the observation, not a permanent ban on the parameter."""
    _, posterior = _reference_posterior(SUPERRES_CHANNEL_DESIGN.design_id)
    for name in PARAMETER_NAMES:
        assert math.isfinite(
            posterior.point_estimate(
                name,
                eigenvalue_floor=PRIOR_INFORMATION_FLOOR,
                min_identifiable_fraction=MIN_IDENTIFIABLE_FRACTION,
            )
        )


def test_two_planes_still_refuse_the_flux_pair() -> None:
    """A second plane fixes the axial degeneracy and leaves the flux one exactly where it was."""
    _, posterior = _reference_posterior(TWO_PLANE_DESIGN.design_id)
    for name in ("log_psf_sigma_um", "z_um", "log_background_photons"):
        assert math.isfinite(
            posterior.point_estimate(
                name,
                eigenvalue_floor=PRIOR_INFORMATION_FLOOR,
                min_identifiable_fraction=MIN_IDENTIFIABLE_FRACTION,
            )
        )
    for name in ("log_length_um", "log_linear_density"):
        with pytest.raises(UnidentifiableDirectionError):
            posterior.point_estimate(
                name,
                eigenvalue_floor=PRIOR_INFORMATION_FLOOR,
                min_identifiable_fraction=MIN_IDENTIFIABLE_FRACTION,
            )


def test_identifiability_is_built_from_the_likelihood_and_not_the_prior() -> None:
    """A prior-augmented precision would hide the degeneracy; the report must not use one."""
    _, posterior = _reference_posterior(SINGLE_PLANE_DESIGN.design_id)
    assert posterior.identifiability.smallest_eigenvalue < 1e-6
    posterior_precision = np.linalg.inv(posterior.covariance)
    assert float(np.min(np.linalg.eigvalsh(posterior_precision))) > 1.0


# ------------------------------------------------------------------------------------------------
# 4. the expansion route: add a modality, measure the lift
# ------------------------------------------------------------------------------------------------


def test_second_focal_plane_lifts_the_axial_degeneracy() -> None:
    contrast = measure_modality_contrast(
        REFERENCE_TRUTH,
        REFERENCE_SCENE,
        SINGLE_PLANE_DESIGN,
        TWO_PLANE_DESIGN,
        lift_ratio_required=LIFT_RATIO_REQUIRED,
    )
    assert contrast.verdict is ModalityVerdict.LIFTED
    # Measured: 1.7e-12 -> 4.4e3, and a condition number of 7.1e15 -> 1.7e4.
    assert contrast.baseline_directional_information < 1e-9
    assert contrast.augmented_directional_information > 1e3
    assert contrast.baseline_condition_number > 1e12
    assert contrast.augmented_condition_number < 1e6
    assert "observation-design" in contrast.routes_to


def test_second_focal_plane_does_NOT_lift_the_flux_degeneracy() -> None:
    """The negative half of the contrast, which is the half that carries information."""
    contrast = measure_modality_contrast(
        REFERENCE_TRUTH,
        REFERENCE_SCENE,
        SINGLE_PLANE_DESIGN,
        TWO_PLANE_DESIGN,
        lift_ratio_required=LIFT_RATIO_REQUIRED,
        probe_direction=FLUX_DEGENERACY_DIRECTION,
    )
    assert contrast.verdict is ModalityVerdict.NOT_LIFTED
    # Measured: 1.80 -> 4.10, a factor of 2.28.
    assert 1.5 < contrast.information_gain_ratio < 5.0
    assert contrast.augmented_directional_information < PRIOR_INFORMATION_FLOOR
    assert "still required" in contrast.routes_to


def test_an_independent_resolution_channel_lifts_the_flux_degeneracy() -> None:
    contrast = measure_modality_contrast(
        REFERENCE_TRUTH,
        REFERENCE_SCENE,
        SINGLE_PLANE_DESIGN,
        SUPERRES_CHANNEL_DESIGN,
        lift_ratio_required=LIFT_RATIO_REQUIRED,
        probe_direction=FLUX_DEGENERACY_DIRECTION,
    )
    assert contrast.verdict is ModalityVerdict.LIFTED
    # Measured: 1.80 -> 75.1, a factor of 41.8, with the condition number falling to 729.
    assert contrast.information_gain_ratio > 20.0
    assert contrast.augmented_directional_information > 10.0 * PRIOR_INFORMATION_FLOOR
    assert contrast.augmented_condition_number < 1e4


def test_a_lift_ratio_of_one_is_rejected_as_a_contract() -> None:
    with pytest.raises(ValueError, match="must exceed one"):
        measure_modality_contrast(
            REFERENCE_TRUTH,
            REFERENCE_SCENE,
            SINGLE_PLANE_DESIGN,
            TWO_PLANE_DESIGN,
            lift_ratio_required=1.0,
        )


# ------------------------------------------------------------------------------------------------
# 5. calibration and its negative control
# ------------------------------------------------------------------------------------------------


def test_matched_coverage_is_close_to_nominal_but_not_binomially_clean() -> None:
    """The falsifier fires in its weak form, and the weak result is what is written down.

    Measured over 240 declared-seed draws on the matched renderer: 0.946, 0.954, 0.979, 0.904,
    0.946 against a nominal 0.95.  Every parameter is within five points, and the background
    coverage is far enough below nominal that the binomial gate rejects it — a deficit that persists
    across seeds and grows to ``z = -4.1`` at 1000 draws.  That routes to a nonlinear posterior
    sampler, not to a narrower claim.
    """
    campaign = _matched_campaign(1.0)
    for name in PARAMETER_NAMES:
        report = campaign.per_parameter[name]
        assert abs(report.empirical_coverage - NOMINAL_LEVEL) < 0.05, name
        assert report.n_draws == CALIBRATION_DRAWS
        assert report.seed == SEED
    # Measured worst |z| = 3.26 against the declared 3-sigma binomial band.
    assert campaign.worst_z_score < 5.0
    assert campaign.worst_z_score > SIGMA_GATE
    assert campaign.uncalibrated_parameters == ("log_background_photons",)
    # The MAP iteration must not be leaning on its safety rails for this result.
    assert campaign.trust_region_hits <= 1
    assert campaign.unconverged_draws < 0.05 * CALIBRATION_DRAWS


def test_the_overconfident_negative_control_is_caught() -> None:
    """Halve every reported sigma and the coverage test must reject it, loudly."""
    honest = _matched_campaign(1.0)
    overconfident = _matched_campaign(0.5)
    assert overconfident.all_calibrated is False
    for name in PARAMETER_NAMES:
        measured = overconfident.per_parameter[name].empirical_coverage
        # Measured: 0.696, 0.713, 0.754, 0.700, 0.679.
        assert measured < 0.80, name
        assert measured < honest.per_parameter[name].empirical_coverage
    # Measured: 19.25 against the honest 3.26 — a factor of six.
    assert overconfident.worst_z_score > 15.0
    assert overconfident.worst_z_score > 4.0 * honest.worst_z_score


def test_campaign_record_declares_itself_self_scoring_and_oracle_only() -> None:
    record = _matched_campaign(1.0).to_dict()
    assert record["self_scoring"] is True
    assert record["evidence_class"] == "oracle-only"
    assert record["seed"] == SEED
    assert record["truth_renderer_family"] == RendererFamily.PIXEL_INTEGRATED_LINE.value


def test_a_wider_interval_covers_more_often_than_a_narrower_one() -> None:
    """Positive control on the coverage machinery itself."""
    honest = _matched_campaign(1.0)
    overconfident = _matched_campaign(0.5)
    for name in PARAMETER_NAMES:
        assert (
            honest.per_parameter[name].empirical_coverage
            > overconfident.per_parameter[name].empirical_coverage
        )


# ------------------------------------------------------------------------------------------------
# 6. A6: the self-scoring control, measured rather than argued
# ------------------------------------------------------------------------------------------------


def test_the_two_renderer_families_really_do_differ() -> None:
    """Without this the mismatch control would be measuring nothing."""
    matched = expected_photons(REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN)
    aliased = expected_photons(
        REFERENCE_TRUTH,
        REFERENCE_SCENE,
        SINGLE_PLANE_DESIGN,
        family=RendererFamily.ALIASED_POINT_CHAIN,
    )
    signal = float(np.linalg.norm(matched - REFERENCE_TRUTH.background_photons))
    residual = float(np.linalg.norm(aliased - matched)) / signal
    # Measured: 1.75e-2 relative.
    assert 1e-3 < residual < 1e-1
    miscalibrated = expected_photons(
        REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN, rayleigh_scale=1.15
    )
    assert float(np.linalg.norm(miscalibrated - matched)) / signal > residual


def test_one_plane_absorbs_an_axial_miscalibration_and_two_planes_do_not() -> None:
    """The honest ceiling: the modality that bought identifiability also bought exposure.

    A 15 percent error in the axial Rayleigh range is invisible to a single widefield frame, because
    the parameters it would have to move are degenerate there anyway; the same error destroys the
    two-plane design, whose whole value was that it broke that degeneracy.
    """
    single = _mismatch(SINGLE_PLANE_DESIGN.design_id)
    double = _mismatch(TWO_PLANE_DESIGN.design_id)
    # Measured: 3.3 coverage points lost on one plane, 33.8 on two.
    assert single.worst_coverage_degradation < 0.10
    assert double.worst_coverage_degradation > 0.20
    assert single.collapsed_parameters == ()
    assert set(double.collapsed_parameters) == {"log_psf_sigma_um", "z_um"}
    assert double.mismatched.per_parameter["z_um"].empirical_coverage < 0.75
    assert double.matched.per_parameter["z_um"].empirical_coverage > 0.90


def test_mismatch_report_is_not_labelled_self_scoring() -> None:
    """It is the ONE artifact here that leaves the renderer family, and it says so."""
    record = _mismatch(SINGLE_PLANE_DESIGN.design_id).to_dict()
    assert record["self_scoring"] is False
    assert record["evidence_class"] == "oracle-only"
    assert record["mismatched"]["rayleigh_scale"] == pytest.approx(1.15)
    assert record["matched"]["rayleigh_scale"] == pytest.approx(1.0)
    assert record["relative_image_residual"] > 0.0


# ------------------------------------------------------------------------------------------------
# 7. provenance: synthetic ancestry cannot be cited, and a citation cannot launder it
# ------------------------------------------------------------------------------------------------


def test_no_sandbox_source_is_citable_as_biological_evidence() -> None:
    graph = register_sandbox_sources().freeze()
    assert graph.citable_as_biological_evidence() == ()
    for source_id in S11_SOURCE_IDS:
        assert graph.exclusion(source_id) is not None, source_id


def test_a_cited_claim_is_still_blocked_by_its_synthetic_ancestry() -> None:
    """Master-plan A6, structurally: provenance beats a citation and a green status."""
    graph = register_sandbox_sources().freeze()
    claim = "s11.claim.cited-but-synthetic-backed"
    record = graph.record(claim)
    assert record.citation is not None
    exclusion = graph.exclusion(claim)
    assert exclusion is not None
    assert exclusion.reason is ExclusionReason.DERIVED_FROM_EXCLUDED
    assert "s11.renderer.pixel-integrated-optics" in graph.ancestors(claim)


def test_registering_into_an_existing_registry_chains() -> None:
    registry = register_sandbox_sources(SourceRegistry())
    assert isinstance(registry, SourceRegistry)
    assert len(registry) == len(S11_SOURCE_IDS)
    with pytest.raises(TypeError, match="SourceRegistry"):
        register_sandbox_sources("not-a-registry")  # type: ignore[arg-type]


# ------------------------------------------------------------------------------------------------
# 8. the pre-registration card
# ------------------------------------------------------------------------------------------------


def test_experiment_card_declares_the_falsifier_and_the_expansions() -> None:
    card = SANDBOX_S11_CARD
    assert isinstance(card, SandboxExperimentCard)
    assert card.wet_lab is False
    assert "coverage" in card.falsifier
    assert any("modality" in route for route in card.failure_expansions)
    assert any("renderer" in control for control in card.adversarial_controls)
    assert any("overconfident" in control for control in card.negative_controls)
    assert len(card.card_hash) == 64
    assert card.to_dict()["experiment_id"] == card.experiment_id


def test_card_manifest_hashes_bind_the_declared_configuration() -> None:
    assert REFERENCE_SCENE.scene_hash in SANDBOX_S11_CARD.manifest_hashes
    for design in (SINGLE_PLANE_DESIGN, TWO_PLANE_DESIGN, SUPERRES_CHANNEL_DESIGN):
        assert design.design_hash in SANDBOX_S11_CARD.manifest_hashes


# ------------------------------------------------------------------------------------------------
# 9. validation surface
# ------------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"length_um": 0.0}, "positive-finite"),
        ({"linear_density": -1.0}, "positive-finite"),
        ({"psf_sigma_um": math.inf}, "positive-finite"),
        ({"background_photons": 0.0}, "positive-finite"),
        ({"z_um": math.nan}, "finite"),
    ],
)
def test_inverse_parameters_reject_unphysical_values(kwargs, message) -> None:
    base = {
        "length_um": 0.1,
        "linear_density": 1.0,
        "psf_sigma_um": 0.1,
        "background_photons": 1.0,
        "z_um": 0.0,
    }
    with pytest.raises(ValueError, match=message):
        InverseParameters(**{**base, **kwargs})


def test_parameter_vector_round_trips() -> None:
    restored = InverseParameters.from_vector(REFERENCE_TRUTH.to_vector())
    assert restored.length_um == pytest.approx(REFERENCE_TRUTH.length_um)
    assert restored.linear_density == pytest.approx(REFERENCE_TRUTH.linear_density)
    assert restored.psf_sigma_um == pytest.approx(REFERENCE_TRUTH.psf_sigma_um)
    assert restored.background_photons == pytest.approx(REFERENCE_TRUTH.background_photons)
    assert restored.z_um == pytest.approx(REFERENCE_TRUTH.z_um)


def test_design_requires_at_least_one_focal_plane() -> None:
    with pytest.raises(ValueError, match="at least one plane"):
        ObservationDesign(design_id="empty", focal_planes_um=())


def test_prior_sigma_must_be_positive() -> None:
    with pytest.raises(ValueError, match="strictly positive"):
        GaussianLogPrior(mean=np.zeros(5), sigma=np.array([1.0, 1.0, 0.0, 1.0, 1.0]))


def test_scene_rejects_a_degenerate_point_chain() -> None:
    with pytest.raises(ValueError, match="at least two"):
        ImagingScene(
            centre_xy_um=(1.0, 1.0),
            orientation_rad=0.0,
            height_px=8,
            width_px=8,
            pixel_size_um=0.065,
            exposure_s=1.0,
            rayleigh_range_um=0.4,
            read_noise_photons=1.0,
            point_chain_labels=1,
        )


def test_unknown_parameter_name_is_a_key_error() -> None:
    _, posterior = _reference_posterior(SINGLE_PLANE_DESIGN.design_id)
    with pytest.raises(KeyError, match="unknown parameter"):
        posterior.marginal_sigma("psf_sigma_um")


def test_eigenvalue_floor_has_no_default() -> None:
    """A floor chosen by this module would be a threshold nobody declared."""
    report = _single_plane_report()
    with pytest.raises(TypeError):
        report.near_null_directions()  # type: ignore[call-arg]
    with pytest.raises(ValueError, match="positive-finite"):
        report.near_null_directions(eigenvalue_floor=0.0)


def test_observation_shape_must_match_the_design() -> None:
    expected = expected_photons(REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN)
    with pytest.raises(ValueError, match="must have shape"):
        infer_posterior(expected, REFERENCE_SCENE, TWO_PLANE_DESIGN, REFERENCE_PRIOR)


def test_render_rejects_an_unknown_renderer_family() -> None:
    with pytest.raises(TypeError, match="RendererFamily"):
        expected_photons(
            REFERENCE_TRUTH,
            REFERENCE_SCENE,
            SINGLE_PLANE_DESIGN,
            family="pixel-integrated-line",  # type: ignore[arg-type]
        )


def test_emitter_types_used_are_the_lane_d_primitives() -> None:
    """The forward model must be Lane D's renderer, not a private reimplementation of one."""
    from aleph.virtual_cell.sandbox_inverse import _emitters

    line, sampling = _emitters(
        REFERENCE_TRUTH, REFERENCE_SCENE, RendererFamily.PIXEL_INTEGRATED_LINE
    )
    assert len(line) == 1 and isinstance(line[0], LineEmitter)
    assert line[0].length_um == pytest.approx(REFERENCE_TRUTH.length_um)
    chain, aliased_sampling = _emitters(
        REFERENCE_TRUTH, REFERENCE_SCENE, RendererFamily.ALIASED_POINT_CHAIN
    )
    assert len(chain) == REFERENCE_SCENE.point_chain_labels
    assert all(isinstance(emitter, PointEmitter) for emitter in chain)
    assert sampling != aliased_sampling
    total = sum(emitter.photon_rate for emitter in chain)
    assert total == pytest.approx(REFERENCE_TRUTH.total_photons)


def test_sample_observation_is_reproducible_from_its_seed() -> None:
    expected = expected_photons(REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN)
    first = sample_observation(expected, REFERENCE_SCENE, seed=7)
    second = sample_observation(expected, REFERENCE_SCENE, seed=7)
    third = sample_observation(expected, REFERENCE_SCENE, seed=8)
    assert np.array_equal(first, second)
    assert not np.array_equal(first, third)
    assert first.shape == expected.shape

"""CPU oracle tests for the S0 numerical-invariance battery.

Every system here is analytic: an Ornstein-Uhlenbeck process whose stationary moments are known in
closed form, a damped oscillator whose explicit-Euler stability step is ``2 zeta / omega`` exactly, and
a backward-Euler solve whose inner iteration contracts at a known rate.  So what the battery reports can
be checked against an oracle rather than against itself.  No cell physics, no simulation trace, no GPU.

Two of these tests exist to prove the battery can FAIL: :func:`test_positive_control_flags_the_spurious
_limit_cycle` and :func:`test_positive_control_bracket_contains_the_analytic_stability_step`.  A battery
that only ever passes has not been shown to be able to fail, and the round-1 modules this one builds on
were valuable mainly where they measured something and rejected it.

Every threshold used here is DECLARED before the ladders run.  The regime contract is the one round 1
declared for its own analytic controls (``test_regime_sweep.py``), reused verbatim rather than retuned;
the invariance contract's values are argued in ``_INVARIANCE`` and were not adjusted after any
measurement.  The seed 20260729 is round 1's declared surrogate seed and is fixed here for the same
reason: a seed chosen after seeing a z-score is a seed chosen to pass.
"""

from __future__ import annotations

import dataclasses
import math
import subprocess
import sys
import textwrap
from functools import cache
from pathlib import Path

import numpy as np
import pytest

from aleph.virtual_cell.regime import Regime, RegimeEvidenceContract
from aleph.virtual_cell.sandbox_invariance import (
    EVIDENCE_AUTHORITY,
    InvarianceAxis,
    InvarianceContract,
    InvarianceReport,
    InvarianceVerdict,
    LadderRun,
    StatisticRole,
    StatisticSpec,
    centered_square,
    damped_oscillator_series,
    euler_growth_per_step,
    euler_stability_step_s,
    identity,
    implicit_relaxation_run,
    invariance_experiment_card,
    knob_invariance,
    matched_time_difference,
    ou_refinement_ladder,
    ou_series,
    seed_invariance,
    shared_brownian_increments,
)
from aleph.virtual_cell.sweep import ExpansionRoute

# ------------------------------------------------------------------------------------------------
# Pre-declared contracts.  Nothing below this block was changed after a measurement.
# ------------------------------------------------------------------------------------------------

#: Round 1's declared contract for analytic controls, copied verbatim from ``test_regime_sweep.py``.
#: Reused rather than re-derived so that this lane cannot pick classifier thresholds that suit it.
_REGIME = RegimeEvidenceContract(
    min_samples=512,
    sokal_window_c=5.0,
    min_correlation_windows=20.0,
    min_n_eff=30.0,
    drift_sigma=1.0,
    fixed_point_relative_amplitude=0.01,
    welch_segment_count=4,
    peak_prominence_ratio_min=8.0,
    baseline_octave_ratio=2.0,
    peak_power_fraction_min=0.05,
    secondary_peak_power_fraction=0.05,
    harmonic_relative_tolerance=0.02,
    max_rational_denominator=4,
    limit_cycle_max_width_bins=2.5,
    quasicycle_min_width_bins=4.0,
    white_noise_tau_dt_max=1.0,
    white_noise_flatness_min=0.5,
    determinism_min=0.5,
    predictability_z_min=4.0,
    predictability_neighbours=4,
    n_surrogates=8,
    surrogate_seed=20260729,
    embedding_dimension=4,
    recurrence_rate_target=0.05,
    min_diagonal_length=2,
    recurrence_max_points=800,
    theiler_tau_multiple=2.0,
    lyapunov_fit_tau_multiple=1.0,
    lyapunov_min_fit_lags=5,
    lyapunov_max_fit_lags=60,
)

_INVARIANCE = InvarianceContract(
    # Three rungs is the fewest that can bracket a threshold and still show a trend either side.
    min_ladder_rungs=3,
    # The rungs of these ladders cover identical durations by construction, so 1 % is a protocol
    # tolerance for float rounding, not a physics threshold.
    duration_relative_tolerance=0.01,
    # Both taken from the regime contract above, so a rung the classifier considers long enough is a
    # rung the battery considers long enough; two different answers to that question is how a run gets
    # a label and a refused interval at the same time.
    min_correlation_windows=20.0,
    min_n_eff=30.0,
    # A 3-sigma consistency criterion in the statistic's OWN stated error.
    statistic_z_max=3.0,
    # Two decades.  A physical observable does not move by 100x when only the solver knob moved; a
    # residual does, and that gap is what separates the two failures.
    residual_tracking_dynamic_range_min=100.0,
    # Eight seeds: the fewest at which a standard deviation has better than ~30 % relative error.
    min_seeds=8,
    min_seed_regime_fraction=0.8,
    # A factor of two either side of 1.  With 8 seeds the chi-distribution of a standard-deviation
    # estimate already carries ~25 % relative error, so a tighter band would fail on the estimator's
    # own noise; a wider one would accept an interval that is out by an order of magnitude.
    ess_calibration_ratio_min=0.5,
    ess_calibration_ratio_max=2.0,
    coverage_z=1.96,
)

_SEED = 20260729

_MEAN = StatisticSpec(name="mean", transform=identity, role=StatisticRole.PHYSICAL, units="x")
_VARIANCE = StatisticSpec(
    name="variance", transform=centered_square, role=StatisticRole.PHYSICAL, units="x^2"
)
_RESIDUAL_AS_PHYSICAL = StatisticSpec(
    name="mean_final_residual",
    transform=identity,
    role=StatisticRole.PHYSICAL,
    units="x",
    channel="final_residual",
)
_RESIDUAL_AS_DIAGNOSTIC = StatisticSpec(
    name="mean_final_residual",
    transform=identity,
    role=StatisticRole.SOLVER_DIAGNOSTIC,
    units="x",
    channel="final_residual",
)

# Damped noisy oscillator: the positive control.  zeta = 0.3 at 0.5 Hz puts the explicit-Euler
# stability step at 2 zeta / omega = 0.1910 s, comfortably inside a dt ladder that also has rungs a
# factor of 20 finer.  sigma is set so the resting fluctuation is ~1e-3 of the declared amplitude
# scale, i.e. a decade below the classifier's fixed-point threshold of 1e-2 — the continuous system
# sits at rest, and anything the ladder shows above that came from the integrator.
_F0_HZ, _ZETA, _OSC_SIGMA, _OSC_SCALE = 0.5, 0.3, 6.0e-3, 1.0
_OSC_DURATION_S, _OSC_BASE_DT_S = 300.0, 0.01
_OSC_FACTORS = (1, 2, 4, 15, 20)

# Ornstein-Uhlenbeck: theta = 2 /s gives tau_int = 0.5 s, so 1000 s is 2000 correlation times.
_THETA, _SIGMA = 2.0, 1.0
_OU_DURATION_S, _OU_BASE_DT_S = 1000.0, 0.03125
_OU_FACTORS = (1, 2, 4, 8)

_SOLVE_DT_S, _SOLVE_DURATION_S = 0.05, 1000.0
_TOLERANCES = (1e-4, 1e-6, 1e-8, 1e-10)
_BUDGETS = (2.0, 3.0, 4.0, 6.0)


# ------------------------------------------------------------------------------------------------
# Cached ladders.  Each is one experiment; the tests read different facets of the same runs.
# ------------------------------------------------------------------------------------------------


def _run_ladder(runs, *, axis, reference_knob, statistics, observable="x", scale=1.0):
    return knob_invariance(
        runs,
        axis=axis,
        observable=observable,
        amplitude_scale=scale,
        reference_knob=reference_knob,
        statistics=statistics,
        regime_contract=_REGIME,
        contract=_INVARIANCE,
    )


@cache
def _ou_runs(scheme: str) -> tuple[LadderRun, ...]:
    return ou_refinement_ladder(
        theta_per_s=_THETA,
        sigma=_SIGMA,
        base_dt_s=_OU_BASE_DT_S,
        factors=_OU_FACTORS,
        duration_s=_OU_DURATION_S,
        seed=_SEED,
        scheme=scheme,
    )


@cache
def _ou_report(scheme: str) -> InvarianceReport:
    return _run_ladder(
        _ou_runs(scheme),
        axis=InvarianceAxis.TIME_STEP,
        reference_knob=_OU_BASE_DT_S,
        statistics=(_MEAN, _VARIANCE),
    )


@cache
def _oscillator_runs(scheme: str) -> tuple[LadderRun, ...]:
    increments, steps = shared_brownian_increments(
        duration_s=_OSC_DURATION_S, base_dt_s=_OSC_BASE_DT_S, seed=_SEED
    )
    runs = []
    for factor in _OSC_FACTORS:
        dt = _OSC_BASE_DT_S * factor
        coarse = increments[: (steps // factor) * factor].reshape(-1, factor).sum(axis=1)
        runs.append(
            LadderRun(
                knob_value=dt,
                dt_s=dt,
                series=damped_oscillator_series(
                    frequency_hz=_F0_HZ,
                    damping_ratio=_ZETA,
                    sigma=_OSC_SIGMA,
                    dt_s=dt,
                    increments=coarse,
                    scheme=scheme,
                ),
                seed=_SEED,
            )
        )
    return tuple(runs)


@cache
def _oscillator_report(scheme: str) -> InvarianceReport:
    return _run_ladder(
        _oscillator_runs(scheme),
        axis=InvarianceAxis.TIME_STEP,
        reference_knob=_OSC_BASE_DT_S,
        statistics=(_MEAN, _VARIANCE),
        scale=_OSC_SCALE,
    )


@cache
def _solver_runs(axis: InvarianceAxis) -> tuple[LadderRun, ...]:
    increments, _ = shared_brownian_increments(
        duration_s=_SOLVE_DURATION_S, base_dt_s=_SOLVE_DT_S, seed=_SEED
    )
    knobs = _TOLERANCES if axis is InvarianceAxis.TOLERANCE else _BUDGETS
    runs = []
    for knob in knobs:
        solve = implicit_relaxation_run(
            theta_per_s=_THETA,
            sigma=_SIGMA,
            dt_s=_SOLVE_DT_S,
            increments=increments,
            tolerance=knob if axis is InvarianceAxis.TOLERANCE else 0.0,
            max_iterations=50 if axis is InvarianceAxis.TOLERANCE else int(knob),
        )
        runs.append(solve.ladder_run(knob, seed=_SEED))
    return tuple(runs)


@cache
def _solver_report(axis: InvarianceAxis, residual_role: StatisticRole) -> InvarianceReport:
    residual = (
        _RESIDUAL_AS_PHYSICAL
        if residual_role is StatisticRole.PHYSICAL
        else _RESIDUAL_AS_DIAGNOSTIC
    )
    knobs = _TOLERANCES if axis is InvarianceAxis.TOLERANCE else _BUDGETS
    return _run_ladder(
        _solver_runs(axis),
        axis=axis,
        reference_knob=min(knobs) if axis is InvarianceAxis.TOLERANCE else max(knobs),
        statistics=(_MEAN, _VARIANCE, residual),
    )


@cache
def _seed_ensemble(count: int) -> tuple[LadderRun, ...]:
    runs = []
    for index in range(count):
        seed = _SEED + index
        increments, _ = shared_brownian_increments(duration_s=400.0, base_dt_s=0.05, seed=seed)
        runs.append(
            LadderRun(
                knob_value=float(seed),
                dt_s=0.05,
                series=ou_series(
                    theta_per_s=_THETA,
                    sigma=_SIGMA,
                    dt_s=0.05,
                    increments=increments,
                    scheme="exact",
                ),
                seed=seed,
            )
        )
    return tuple(runs)


@cache
def _seed_report(count: int):
    return seed_invariance(
        _seed_ensemble(count),
        observable="x",
        amplitude_scale=1.0,
        statistic=_MEAN,
        regime_contract=_REGIME,
        contract=_INVARIANCE,
    )


# ------------------------------------------------------------------------------------------------
# 0. CPU-only import contract
# ------------------------------------------------------------------------------------------------


def test_importing_the_battery_pulls_in_no_gpu_runtime() -> None:
    """The battery must import without Warp, in a FRESH interpreter.

    Asserting ``"warp" not in sys.modules`` in-process would assert about the whole pytest session
    rather than about this import, so it would pass alone and fail under the full suite.  Two round-1
    lanes hit exactly that; the probe is a subprocess for that reason.
    """
    script = textwrap.dedent(
        """
        import sys
        import aleph.virtual_cell.sandbox_invariance as module

        assert module.EVIDENCE_AUTHORITY == "unverified-observer"
        for blocked in ("warp", "torch", "jax"):
            assert blocked not in sys.modules, f"sandbox_invariance import pulled in {blocked}"
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
# 1. The contract owns every threshold, and has no defaults
# ------------------------------------------------------------------------------------------------


def test_invariance_contract_has_no_defaults() -> None:
    """A field with a default is a threshold the battery would choose for a silent caller."""
    for field in dataclasses.fields(InvarianceContract):
        assert field.default is dataclasses.MISSING, field.name
        assert field.default_factory is dataclasses.MISSING, field.name
    with pytest.raises(TypeError):
        InvarianceContract()  # type: ignore[call-arg]


def test_contract_rejects_a_duration_tolerance_that_admits_the_error_it_catches() -> None:
    payload = dataclasses.asdict(_INVARIANCE) | {"duration_relative_tolerance": 1.0}
    with pytest.raises(ValueError, match="duration_relative_tolerance"):
        InvarianceContract(**payload)


def test_contract_rejects_a_calibration_band_that_excludes_one() -> None:
    """A band not straddling 1 asserts the ESS estimator is wrong before measuring it."""
    payload = dataclasses.asdict(_INVARIANCE) | {"ess_calibration_ratio_min": 1.2}
    with pytest.raises(ValueError, match="straddle 1"):
        InvarianceContract(**payload)


def test_battery_refuses_an_undeclared_contract() -> None:
    """A caller cannot pass a bare mapping and have thresholds invented for it."""
    with pytest.raises(TypeError, match="declared before the run"):
        knob_invariance(
            _ou_runs("exact"),
            axis=InvarianceAxis.TIME_STEP,
            observable="x",
            amplitude_scale=1.0,
            reference_knob=_OU_BASE_DT_S,
            statistics=(_MEAN,),
            regime_contract=_REGIME,
            contract={"statistic_z_max": 3.0},  # type: ignore[arg-type]
        )


def test_reference_rung_must_be_declared_and_present() -> None:
    """A reference chosen after seeing the ladder is a reference chosen to pass."""
    with pytest.raises(ValueError, match="not a rung of this ladder"):
        _run_ladder(
            _ou_runs("exact"),
            axis=InvarianceAxis.TIME_STEP,
            reference_knob=0.017,
            statistics=(_MEAN,),
        )


def test_duplicate_knob_values_are_refused() -> None:
    runs = _ou_runs("exact")
    duplicated = (*runs, LadderRun(knob_value=runs[0].knob_value, dt_s=0.5, series=np.zeros(64)))
    with pytest.raises(ValueError, match="distinct"):
        _run_ladder(
            duplicated,
            axis=InvarianceAxis.TIME_STEP,
            reference_knob=_OU_BASE_DT_S,
            statistics=(_MEAN,),
        )


# ------------------------------------------------------------------------------------------------
# 2. Rung and statistic construction
# ------------------------------------------------------------------------------------------------


def test_a_rung_accepts_a_non_finite_series_because_divergence_is_a_result() -> None:
    run = LadderRun(knob_value=1.0, dt_s=0.1, series=[1.0, 2.0, np.inf, 4.0])
    assert run.is_finite is False
    assert run.physical_duration_s == pytest.approx(0.4)


def test_finite_samples_are_not_enough_to_be_analysable() -> None:
    """A rung at 1e243 is every-sample finite and still has no variance, spectrum or ACF.

    This is not hypothetical: it is the rung an explicit scheme produces one factor of two below the
    step at which it overflows, and every estimator in the classifier is a second moment.
    """
    run = LadderRun(knob_value=1.0, dt_s=0.1, series=np.full(64, 1e243))
    assert run.is_finite is True
    assert run.is_analysable is False


def test_a_rung_series_is_immutable_after_construction() -> None:
    run = LadderRun(knob_value=1.0, dt_s=0.1, series=np.zeros(16))
    with pytest.raises(ValueError):
        run.series[0] = 1.0


def test_a_channel_must_share_the_rungs_sampling_grid() -> None:
    with pytest.raises(ValueError, match="same steps"):
        LadderRun(
            knob_value=1.0,
            dt_s=0.1,
            series=np.zeros(16),
            channels={"residual": np.zeros(8)},
        )


def test_a_transform_that_changes_the_sample_count_is_refused() -> None:
    """Its correlation time would be measured on a different time axis than the series."""
    decimating = StatisticSpec(
        name="decimated",
        transform=lambda series: series[::2],
        role=StatisticRole.PHYSICAL,
        units="x",
    )
    with pytest.raises(ValueError, match="preserve the sampling grid"):
        _run_ladder(
            _ou_runs("exact"),
            axis=InvarianceAxis.TIME_STEP,
            reference_knob=_OU_BASE_DT_S,
            statistics=(decimating,),
        )


# ------------------------------------------------------------------------------------------------
# 3. Matched physical time — the protocol error this battery exists to refuse
# ------------------------------------------------------------------------------------------------


def test_matched_time_difference_is_zero_against_a_decimation_of_the_same_run() -> None:
    """A rung compared with itself at half the sampling rate differs by nothing.

    This is the whole matched-time claim in one assertion: the two runs have different sample counts
    and different ``dt``, and a sample-by-sample comparison would report a large difference where the
    physical difference is exactly zero.
    """
    fine = _ou_runs("exact")[0]
    coarse = LadderRun(
        knob_value=fine.dt_s * 2.0, dt_s=fine.dt_s * 2.0, series=fine.series[::2], seed=fine.seed
    )
    assert matched_time_difference(fine, coarse) == pytest.approx(0.0, abs=1e-12)


def test_matched_time_difference_is_root_two_for_independent_stationary_runs() -> None:
    """Two independent realisations of the same stationary process differ by sqrt(2) sigma."""
    members = _seed_ensemble(8)
    measured = matched_time_difference(members[0], members[1])
    assert measured == pytest.approx(math.sqrt(2.0), rel=0.15)


def test_matched_step_count_across_a_dt_ladder_is_refused() -> None:
    """The protocol error that retracted the blebbistatin sweep, reproduced and refused.

    Truncating every rung to the same STEP COUNT is what a naive refinement study does.  The rungs then
    cover 1000 s, 500 s, 250 s and 125 s: a shorter experiment, not a finer one.
    """
    shortest = min(run.n_samples for run in _ou_runs("exact"))
    truncated = tuple(
        LadderRun(
            knob_value=run.knob_value,
            dt_s=run.dt_s,
            series=run.series[:shortest],
            seed=run.seed,
        )
        for run in _ou_runs("exact")
    )
    report = _run_ladder(
        truncated,
        axis=InvarianceAxis.TIME_STEP,
        reference_knob=_OU_BASE_DT_S,
        statistics=(_MEAN, _VARIANCE),
    )
    assert report.verdict is InvarianceVerdict.REFUSED_UNMATCHED_DURATION
    assert "different experiments" in report.reason
    assert report.rungs == ()
    assert report.statistic_verdicts == ()
    assert ExpansionRoute.TIME_MODALITY in report.expansions


def test_a_ladder_covering_too_few_correlation_times_is_refused() -> None:
    """Matched wall duration is not matched statistical duration."""
    short = ou_refinement_ladder(
        theta_per_s=_THETA,
        sigma=_SIGMA,
        base_dt_s=0.0125,
        factors=(1, 2, 4),
        duration_s=8.0,
        seed=_SEED,
        scheme="exact",
    )
    report = _run_ladder(
        short, axis=InvarianceAxis.TIME_STEP, reference_knob=0.0125, statistics=(_MEAN,)
    )
    assert report.regime_verdict is InvarianceVerdict.REFUSED_INSUFFICIENT_EVIDENCE
    assert "correlation times" in report.regime_reason


# ------------------------------------------------------------------------------------------------
# 4. Regime invariance under dt refinement — the negative control
# ------------------------------------------------------------------------------------------------


def test_regime_is_dt_invariant_while_the_trajectory_is_not() -> None:
    """The exact-discretisation OU control: same label at every step, different path at every step."""
    report = _ou_report("exact")
    assert report.regime_verdict is InvarianceVerdict.INVARIANT
    assert report.decided_regimes == frozenset({Regime.STOCHASTIC_NESS})
    assert report.regime_change_bracket is None
    assert report.evidence_authority == EVIDENCE_AUTHORITY

    # The trajectories are NOT invariant, and the battery says so with a number.  The rungs share one
    # Wiener path, so the divergence is the discretisation alone and grows with the step.
    differences = [rung.trajectory_difference_over_std for rung in report.rungs]
    assert differences[0] == pytest.approx(0.0, abs=1e-12)
    assert differences == sorted(differences)
    assert differences[-1] > 0.5


def test_exact_scheme_statistics_are_dt_invariant() -> None:
    """The exact discretisation is dt-exact in both moments; a flag here would be a battery defect."""
    report = _ou_report("exact")
    assert report.verdict is InvarianceVerdict.INVARIANT
    for name in ("mean", "variance"):
        assert report.statistic(name).verdict is InvarianceVerdict.INVARIANT
        assert report.statistic(name).max_z_ess < _INVARIANCE.statistic_z_max


def test_euler_maruyama_variance_bias_is_flagged_while_the_regime_is_not() -> None:
    """A statistic can be dt-biased while the REGIME is invariant, and only the second is the falsifier.

    Euler-Maruyama's stationary variance is ``sigma^2 / (theta (2 - theta dt))``, i.e. the exact
    ``sigma^2 / (2 theta)`` inflated by ``1 / (1 - theta dt / 2)``.  That is a real ``O(dt)``
    discretisation bias in a physical statistic, and the battery must catch it without calling the
    regime an artefact.  The rungs share one Wiener path, so their common sampling error cancels in the
    RATIO — which is why the oracle is checked on the ratio and not on the absolute variance.
    """
    report = _ou_report("euler-maruyama")
    assert report.regime_verdict is InvarianceVerdict.INVARIANT
    assert report.decided_regimes == frozenset({Regime.STOCHASTIC_NESS})
    assert report.statistic("mean").verdict is InvarianceVerdict.INVARIANT

    variance = report.statistic("variance")
    assert variance.verdict is InvarianceVerdict.FLAGGED_STATISTIC_ARTEFACT
    assert report.verdict is InvarianceVerdict.FLAGGED_STATISTIC_ARTEFACT
    assert ExpansionRoute.SOLVER_DESIGN in variance.expansions

    measured = dict(variance.values_by_knob)
    fine, coarse = min(measured), max(measured)
    predicted = (1.0 - _THETA * fine / 2.0) / (1.0 - _THETA * coarse / 2.0)
    assert measured[coarse] / measured[fine] == pytest.approx(predicted, rel=0.01)


# ------------------------------------------------------------------------------------------------
# 5. The positive control — the battery MUST flag this one
# ------------------------------------------------------------------------------------------------


def test_the_continuous_oscillator_is_at_rest_and_the_scheme_is_what_moves() -> None:
    """The system's own scale: the resting fluctuation is a decade below the fixed-point threshold.

    Stated as an assertion so that the positive control cannot be mistaken for a system that genuinely
    oscillates.  Whatever the coarse rung shows, the continuous process is sitting still.
    """
    finest = _oscillator_runs("explicit-euler")[0]
    relative = float(np.std(finest.series)) / _OSC_SCALE
    assert relative < 0.1 * _REGIME.fixed_point_relative_amplitude
    assert euler_growth_per_step(frequency_hz=_F0_HZ, damping_ratio=_ZETA, dt_s=finest.dt_s) < 1.0


def test_positive_control_flags_the_spurious_limit_cycle() -> None:
    """Explicit Euler past its stability step invents a limit cycle, and the battery must say so."""
    report = _oscillator_report("explicit-euler")
    assert report.verdict is InvarianceVerdict.FLAGGED_REGIME_ARTEFACT
    assert report.regime_verdict is InvarianceVerdict.FLAGGED_REGIME_ARTEFACT
    labels = dict(report.regime_by_knob)
    assert labels[0.01] == Regime.FIXED_POINT.value
    assert labels[0.15] == Regime.FIXED_POINT.value
    assert labels[0.2] == Regime.LIMIT_CYCLE.value
    assert report.decided_regimes == frozenset({Regime.FIXED_POINT, Regime.LIMIT_CYCLE})
    assert ExpansionRoute.SOLVER_DESIGN in report.expansions
    assert ExpansionRoute.ADAPTIVE_DESIGN in report.expansions


def test_positive_control_bracket_contains_the_analytic_stability_step() -> None:
    """The bracket is checked against the closed form ``dt* = 2 zeta / omega``, not against itself."""
    report = _oscillator_report("explicit-euler")
    bracket = report.regime_change_bracket
    assert bracket is not None
    analytic = euler_stability_step_s(frequency_hz=_F0_HZ, damping_ratio=_ZETA)
    assert bracket[0] < analytic < bracket[1]
    # The bracket is an interval and never a point estimate: a ladder localises a threshold only to
    # the two rungs it actually ran.
    assert bracket == (0.15, 0.2)
    assert euler_growth_per_step(frequency_hz=_F0_HZ, damping_ratio=_ZETA, dt_s=bracket[0]) < 1.0
    assert euler_growth_per_step(frequency_hz=_F0_HZ, damping_ratio=_ZETA, dt_s=bracket[1]) > 1.0


def test_the_held_out_scheme_does_not_flag_at_the_same_steps() -> None:
    """Held-out intervention: the semi-implicit scheme on the SAME system at the SAME steps.

    This is what makes the flagged verdict a statement about the method rather than about the system.
    Without it, "the regime changed with dt" would be equally consistent with a system that really does
    bifurcate somewhere in that interval.
    """
    report = _oscillator_report("semi-implicit")
    assert report.regime_verdict is InvarianceVerdict.INVARIANT
    assert report.decided_regimes == frozenset({Regime.FIXED_POINT})
    assert report.regime_change_bracket is None
    assert report.verdict is InvarianceVerdict.INVARIANT


def test_a_diverged_rung_is_flagged_and_never_classified() -> None:
    """A non-finite trajectory has no regime, so it is recorded as diverged rather than labelled."""
    increments, steps = shared_brownian_increments(duration_s=1200.0, base_dt_s=0.15, seed=_SEED)
    runs = []
    for factor in (1, 2, 4):
        dt = 0.15 * factor
        coarse = increments[: (steps // factor) * factor].reshape(-1, factor).sum(axis=1)
        runs.append(
            LadderRun(
                knob_value=dt,
                dt_s=dt,
                series=damped_oscillator_series(
                    frequency_hz=_F0_HZ,
                    damping_ratio=_ZETA,
                    sigma=_OSC_SIGMA,
                    dt_s=dt,
                    increments=coarse,
                    scheme="explicit-euler",
                ),
                seed=_SEED,
            )
        )
    report = _run_ladder(
        tuple(runs),
        axis=InvarianceAxis.TIME_STEP,
        reference_knob=0.15,
        statistics=(_MEAN,),
        scale=_OSC_SCALE,
    )
    assert report.verdict is InvarianceVerdict.FLAGGED_NUMERICAL_DIVERGENCE
    diverged = [rung for rung in report.rungs if rung.diverged]
    assert diverged
    for rung in diverged:
        assert rung.regime is None
        assert "never handed to the classifier" in rung.regime_reason
    assert ExpansionRoute.SOLVER_DESIGN in report.expansions


# ------------------------------------------------------------------------------------------------
# 6. Tolerance and iteration budget, including the adversarial statistic
# ------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("axis", [InvarianceAxis.TOLERANCE, InvarianceAxis.ITERATION_BUDGET])
def test_trajectory_statistics_do_not_move_with_the_inner_solve(axis: InvarianceAxis) -> None:
    """Six decades of tolerance, or a tripled budget, must not move a physical statistic."""
    report = _solver_report(axis, StatisticRole.SOLVER_DIAGNOSTIC)
    assert report.verdict is InvarianceVerdict.INVARIANT
    assert report.regime_verdict is InvarianceVerdict.INVARIANT
    for name in ("mean", "variance"):
        assert report.statistic(name).verdict is InvarianceVerdict.INVARIANT
        assert report.statistic(name).max_z_ess < _INVARIANCE.statistic_z_max
    # And the trajectories themselves barely moved, which is the premise the routing below rests on.
    assert max(rung.trajectory_difference_over_std for rung in report.rungs) < 0.05


@pytest.mark.parametrize("axis", [InvarianceAxis.TOLERANCE, InvarianceAxis.ITERATION_BUDGET])
def test_a_residual_declared_physical_is_routed_as_residual_tracking(
    axis: InvarianceAxis,
) -> None:
    """The adversarial case: a statistic that tracks the solver, declared as if it were physics."""
    report = _solver_report(axis, StatisticRole.PHYSICAL)
    residual = report.statistic("mean_final_residual")
    assert residual.verdict is InvarianceVerdict.FLAGGED_RESIDUAL_TRACKING
    assert residual.dynamic_range >= _INVARIANCE.residual_tracking_dynamic_range_min
    assert report.verdict is InvarianceVerdict.FLAGGED_RESIDUAL_TRACKING
    # It routes to the observation model, not only to the solver: no amount of solver work turns a
    # residual into an observable of the system.
    assert ExpansionRoute.OBSERVATION_MODEL_EXTENSION in residual.expansions
    assert "tracking the solver" in residual.reason


@pytest.mark.parametrize("axis", [InvarianceAxis.TOLERANCE, InvarianceAxis.ITERATION_BUDGET])
def test_the_same_residual_declared_a_diagnostic_is_reported_and_not_asserted(
    axis: InvarianceAxis,
) -> None:
    """Declared honestly, the same numbers are measured, reported, and never quotable as physics."""
    report = _solver_report(axis, StatisticRole.SOLVER_DIAGNOSTIC)
    residual = report.statistic("mean_final_residual")
    assert residual.verdict is InvarianceVerdict.NOT_ASSERTED_DIAGNOSTIC
    assert residual.expansions == ()
    assert residual.dynamic_range >= 1e3
    assert "may not be quoted as a physical result" in residual.reason
    assert report.verdict is InvarianceVerdict.INVARIANT


def test_residual_tracking_needs_an_invariant_anchor_and_says_so_when_it_has_none() -> None:
    """Without a physical statistic that did NOT move, the two failures are indistinguishable.

    The honest verdict is then the ordinary invariance flag: a ladder on which everything moved cannot
    show that the trajectory stayed put.
    """
    residual_only = tuple(
        LadderRun(
            knob_value=run.knob_value,
            dt_s=run.dt_s,
            series=run.channel_series("final_residual"),
            seed=run.seed,
        )
        for run in _solver_runs(InvarianceAxis.TOLERANCE)
    )
    report = _run_ladder(
        residual_only,
        axis=InvarianceAxis.TOLERANCE,
        reference_knob=min(_TOLERANCES),
        statistics=(
            StatisticSpec(
                name="mean_final_residual",
                transform=identity,
                role=StatisticRole.PHYSICAL,
                units="x",
            ),
        ),
        observable="final_residual",
    )
    verdict = report.statistic("mean_final_residual")
    assert verdict.verdict is InvarianceVerdict.FLAGGED_STATISTIC_ARTEFACT
    assert verdict.verdict is not InvarianceVerdict.FLAGGED_RESIDUAL_TRACKING


# ------------------------------------------------------------------------------------------------
# 7. ESS-aware comparison — and the trap it avoids
# ------------------------------------------------------------------------------------------------


def test_ess_error_is_dt_invariant_while_the_naive_error_shrinks_as_sqrt_dt() -> None:
    """The mechanism of the trap, measured.

    Halving ``dt`` doubles the sample count and doubles the correlation time in units of steps.  The
    ESS error is therefore unchanged — it depends on the physical duration and the physical correlation
    time, neither of which moved — while ``sigma/sqrt(N)`` shrinks by ``sqrt(2)``.  A ladder judged on
    naive errors would conclude its finest rung is the most precise, purely because that rung's samples
    are the most correlated.
    """
    records = [rung.statistic("mean") for rung in _ou_report("exact").rungs]
    ess = [record.sem_ess for record in records]
    naive = [record.sem_naive for record in records]
    assert max(ess) / min(ess) < 1.05
    for coarse, fine in zip(naive[1:], naive[:-1], strict=True):
        assert coarse / fine == pytest.approx(math.sqrt(2.0), rel=0.05)
    for coarse, fine in zip(records[1:], records[:-1], strict=True):
        assert fine.widening_factor / coarse.widening_factor == pytest.approx(
            math.sqrt(2.0), rel=0.05
        )


def test_the_naive_comparison_manufactures_an_invariance_failure_the_ess_one_does_not() -> None:
    """A measured false positive, on the held-out semi-implicit ladder.

    Its variance is consistent across the ladder at 1.7 ESS standard errors — inside the declared 3 —
    while the same comparison on ``sigma/sqrt(N)`` errors reads 3.9 and would have crossed the same
    declared limit.  Both numbers are carried in the report; only the ESS one is acted on.
    """
    verdict = _oscillator_report("semi-implicit").statistic("variance")
    assert verdict.verdict is InvarianceVerdict.INVARIANT
    assert verdict.max_z_ess < _INVARIANCE.statistic_z_max < verdict.max_z_naive


# ------------------------------------------------------------------------------------------------
# 8. Seeds, done honestly
# ------------------------------------------------------------------------------------------------


def test_the_seed_report_cannot_assert_trajectory_invariance() -> None:
    """Structural: there is no field to set, and the serialised report says why."""
    report = _seed_report(8)
    names = {field.name for field in dataclasses.fields(report)}
    assert not any("invarian" in name for name in names - {"verdict"})
    assert "NOT asserted and not assertable" in report.as_dict()["trajectory_invariance"]


def test_regime_is_seed_invariant_and_the_trajectory_is_not() -> None:
    """Eight seeds, the spread reported, and no single-seed claim anywhere."""
    report = _seed_report(8)
    assert report.verdict is InvarianceVerdict.INVARIANT
    assert report.n_seeds == 8
    assert len(set(report.seeds)) == 8
    assert report.modal_regime == Regime.STOCHASTIC_NESS.value
    assert report.modal_fraction >= _INVARIANCE.min_seed_regime_fraction
    assert report.across_seed_std > 0.0
    # The trajectories differ by ~sqrt(2) standard deviations: the seeds really were independent.
    assert report.mean_trajectory_difference_over_std == pytest.approx(math.sqrt(2.0), rel=0.15)


def test_the_across_seed_spread_calibrates_the_ess_interval_and_condemns_the_naive_one() -> None:
    """Two independent estimates of the same error, compared.

    The across-seed standard deviation of the per-seed mean and the mean within-run ESS standard error
    measure the same quantity for an ergodic process, so their ratio is 1.  The naive
    ``sigma/sqrt(N)`` error is out by the widening factor, which is what it means for it to be wrong.
    """
    report = _seed_report(8)
    assert (
        _INVARIANCE.ess_calibration_ratio_min
        <= report.ess_calibration_ratio
        <= _INVARIANCE.ess_calibration_ratio_max
    )
    assert report.naive_calibration_ratio > _INVARIANCE.ess_calibration_ratio_max
    assert report.mean_sem_naive < report.mean_sem_ess


def test_too_few_seeds_is_refused_rather_than_reported() -> None:
    report = _seed_report(4)
    assert report.verdict is InvarianceVerdict.REFUSED_INSUFFICIENT_EVIDENCE
    assert "not a claim" in report.reason
    assert ExpansionRoute.TIME_MODALITY in report.expansions


def test_a_duplicate_seed_is_refused() -> None:
    members = _seed_ensemble(8)
    duplicated = (*members[:-1], members[0])
    with pytest.raises(ValueError, match="duplicate seed"):
        seed_invariance(
            duplicated,
            observable="x",
            amplitude_scale=1.0,
            statistic=_MEAN,
            regime_contract=_REGIME,
            contract=_INVARIANCE,
        )


def test_a_member_without_a_seed_is_refused() -> None:
    members = _seed_ensemble(8)
    anonymous = (
        *members[:-1],
        LadderRun(knob_value=1.0, dt_s=members[0].dt_s, series=members[-1].series),
    )
    with pytest.raises(ValueError, match="must carry its seed"):
        seed_invariance(
            anonymous,
            observable="x",
            amplitude_scale=1.0,
            statistic=_MEAN,
            regime_contract=_REGIME,
            contract=_INVARIANCE,
        )


def test_a_seed_ensemble_of_unequal_durations_is_refused() -> None:
    members = _seed_ensemble(8)
    ragged = (
        *members[:-1],
        LadderRun(
            knob_value=999.0,
            dt_s=members[0].dt_s,
            series=members[-1].series[: members[-1].n_samples // 2],
            seed=999,
        ),
    )
    with pytest.raises(ValueError, match="dt and sample count fixed"):
        seed_invariance(
            ragged,
            observable="x",
            amplitude_scale=1.0,
            statistic=_MEAN,
            regime_contract=_REGIME,
            contract=_INVARIANCE,
        )


# ------------------------------------------------------------------------------------------------
# 9. The pre-registered card
# ------------------------------------------------------------------------------------------------


def test_the_battery_emits_a_pre_registered_experiment_card() -> None:
    card = invariance_experiment_card(
        manifest_hashes=("a" * 64,), budget_class="cpu-analytic-controls"
    )
    assert card.wet_lab is False
    assert card.positive_controls and card.negative_controls and card.adversarial_controls
    assert card.held_out_interventions
    assert (
        "artefact of the numerical method" in card.question or "numerical method" in card.question
    )
    assert "matched physical duration" in card.falsifier
    assert ExpansionRoute.SOLVER_DESIGN.value in card.failure_expansions
    assert len(card.card_hash) == 64
    assert (
        card.card_hash
        == invariance_experiment_card(
            manifest_hashes=("a" * 64,), budget_class="cpu-analytic-controls"
        ).card_hash
    )


def test_the_card_names_the_positive_control_that_must_fail() -> None:
    """Pre-registration is only worth something if the case that must FAIL is named in advance."""
    card = invariance_experiment_card(manifest_hashes=("b" * 64,), budget_class="cpu")
    assert any("MUST flag" in entry for entry in card.positive_controls)
    assert any("residual" in entry for entry in card.adversarial_controls)
    assert any("semi-implicit" in entry for entry in card.held_out_interventions)


# ------------------------------------------------------------------------------------------------
# 10. Serialisation and authority
# ------------------------------------------------------------------------------------------------


def test_every_report_carries_its_contract_and_claims_no_authority() -> None:
    report = _ou_report("exact")
    payload = report.as_dict()
    assert payload["evidence_authority"] == "unverified-observer"
    assert payload["contract"] == _INVARIANCE.to_dict()
    assert payload["regime_contract"] == _REGIME.to_dict()
    assert len(payload["rungs"]) == len(_OU_FACTORS)
    assert payload["rungs"][0]["statistics"][0]["sem_naive"] > 0.0
    assert _seed_report(8).as_dict()["evidence_authority"] == "unverified-observer"

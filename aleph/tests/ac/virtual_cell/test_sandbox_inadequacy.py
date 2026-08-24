"""S12 model-inadequacy sandbox tests.

The lane's whole claim is a CONTRAST, so the tests are organised around it:

1. the closed-form pieces (F survival, masking fraction, closure gap) against independent truth;
2. the negative control — nothing missing, measured false-positive rate against the declared level;
3. the adversarial controls — noise inflated, noise model violated, and an omission deliberately made
   degenerate with the fitted parameter;
4. the positive control — an omission whose signature cannot be absorbed;
5. the structural bans that keep a proposal from promoting itself.

Every threshold used here is a field of ``InadequacyContract``, declared once in :data:`CONTRACT`
before any measurement, and no assertion compares against a number read off a result.
"""

from __future__ import annotations

import dataclasses
import json
import math
import subprocess
import sys
import textwrap
from functools import cache
from pathlib import Path

import numpy as np
import pytest

from aleph.virtual_cell.archetypes import ArchetypeRing, MissingLaw, MissingLawVerdict
from aleph.virtual_cell.contracts import EvidenceSource, SandboxExperimentCard
from aleph.virtual_cell.sandbox_inadequacy import (
    CANDIDATE_LAWS,
    OMITTED_MODE_WEIGHT,
    REPLICATE_STREAM_ROLES,
    S12_MANIFEST,
    SANDBOX_S12_CARD,
    STREAM_DERIVATION,
    WINDOW_EFOLDS,
    ChannelId,
    ChannelResult,
    ConditionSeries,
    InadequacyContract,
    MissingLawProposal,
    NoCandidateReport,
    NoiseModel,
    OmissionKind,
    ProposalStatus,
    ReplicateStreams,
    ResidualStructureReport,
    SelfPromotionError,
    analyse_conditions,
    build_conditions,
    derive_replicate_streams,
    detect_missing_law,
    experiment_fingerprint,
    f_survival,
    measure_null_false_positive_rate,
    parameter_masking_fraction,
    rate_closure_gap,
    s12_evidence_graph,
)

# --------------------------------------------------------------------------------------------
# the one pre-declared contract every test in this file uses
# --------------------------------------------------------------------------------------------

CONTRACT = InadequacyContract(
    family_alpha=0.01,
    permutation_trials=400,
    permutation_seed=20260729,
    sokal_window_c=4.0,
    coverage_z=1.96,
    min_samples_per_condition=64,
    min_conditions=3,
    false_positive_sigma=3.0,
    assumed_noise_model=NoiseModel.HOMOSCEDASTIC_GAUSSIAN,
)

TRAINING = (0.6, 1.0, 1.6)
HELD_OUT = 3.2
NOISE_SIGMA = 0.02
N_SAMPLES = 256
SEED = 4242

# The three constructed omissions, fixed before any run.  ``degenerate-rate`` and the two matched
# cases differ only in the two ways the module says matter: mode separation, and whether the omitted
# process co-scales with the intervention.
DEGENERATE_SEPARATION = 1.02
DEGENERATE_AMPLITUDE = 8.0
MATCHED_SEPARATION = 3.0
MATCHED_AMPLITUDE = 4.0


def _experiment(
    *,
    omission: OmissionKind,
    separation: float,
    amplitude: float,
    seed: int = SEED,
    noise_model: NoiseModel = NoiseModel.HOMOSCEDASTIC_GAUSSIAN,
    noise_sigma: float = NOISE_SIGMA,
) -> tuple[tuple[ConditionSeries, ...], ConditionSeries]:
    training = build_conditions(
        interventions=TRAINING,
        omission=omission,
        separation=separation,
        deviation_sigma_multiple=amplitude,
        noise_sigma=noise_sigma,
        n_samples=N_SAMPLES,
        seed=seed,
        noise_model=noise_model,
    )
    held_out = build_conditions(
        interventions=(HELD_OUT,),
        omission=omission,
        separation=separation,
        deviation_sigma_multiple=amplitude,
        noise_sigma=noise_sigma,
        n_samples=N_SAMPLES,
        seed=seed + 1,
        noise_model=noise_model,
    )[0]
    return training, held_out


@cache
def _report(
    omission: OmissionKind, separation: float, amplitude: float, seed: int = SEED
) -> ResidualStructureReport:
    training, held_out = _experiment(
        omission=omission, separation=separation, amplitude=amplitude, seed=seed
    )
    return analyse_conditions(
        training,
        held_out,
        contract=CONTRACT,
        experiment_id=f"{omission.value}/sep={separation}/amp={amplitude}",
    )


def _null_report() -> ResidualStructureReport:
    return _report(OmissionKind.NONE, 2.0, 0.0)


def _degenerate_report() -> ResidualStructureReport:
    return _report(OmissionKind.CO_SCALING, DEGENERATE_SEPARATION, DEGENERATE_AMPLITUDE)


def _co_scaling_report() -> ResidualStructureReport:
    return _report(OmissionKind.CO_SCALING, MATCHED_SEPARATION, MATCHED_AMPLITUDE)


def _absolute_report() -> ResidualStructureReport:
    return _report(OmissionKind.ABSOLUTE_TIME, MATCHED_SEPARATION, MATCHED_AMPLITUDE)


# --------------------------------------------------------------------------------------------
# 0. import hygiene
# --------------------------------------------------------------------------------------------


def test_import_pulls_in_no_gpu_or_learning_stack() -> None:
    """Measured in a FRESH interpreter.

    ``assert "warp" not in sys.modules`` in-process asserts about the whole pytest session, not about
    this module's import: it passes alone and fails under the full suite.  Two round-1 lanes hit that
    independently.  A subprocess measures what the name says.

    The second half measures the DELTA rather than the absolute set.  ``aleph.virtual_cell`` already
    imports ``scipy.special`` through the round-1 optics lane, so an absolute ban on scipy would be a
    test of another lane's choices; what this lane is responsible for is adding nothing of its own.
    """
    script = textwrap.dedent(
        """
        import sys

        import aleph.virtual_cell  # the package cost this lane inherits and does not control
        before = set(sys.modules)

        from aleph.virtual_cell.sandbox_inadequacy import detect_missing_law, f_survival

        for blocked in ("warp", "torch", "jax"):
            assert blocked not in sys.modules, f"sandbox_inadequacy import pulled in {blocked}"
        added = {name.split(".")[0] for name in set(sys.modules) - before}
        added -= {"aleph"}
        assert not added, f"sandbox_inadequacy added top-level imports: {sorted(added)}"
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


# --------------------------------------------------------------------------------------------
# 1. closed-form pieces checked against independent truth
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("statistic", "df1", "df2"),
    [
        (0.5, 1, 10),
        (2.0, 3, 20),
        (15.0, 2, 700),
        (1.0, 700, 700),
        (0.001, 5, 5),
        (120.0, 3, 760),
        (490.0, 3, 759),
        (1.0e-6, 2, 2),
    ],
)
def test_f_survival_matches_scipy(statistic: float, df1: int, df2: int) -> None:
    """The same quantity computed twice, by two different implementations."""
    from scipy import stats  # imported here only; the module itself must not need it

    mine = f_survival(statistic, df1, df2)
    theirs = float(stats.f.sf(statistic, df1, df2))
    assert mine == pytest.approx(theirs, rel=1e-10, abs=1e-300)


def test_f_survival_boundaries_and_underflow_floor() -> None:
    """A p-value of exactly zero would claim an impossibility the arithmetic cannot support."""
    import sys as _sys

    assert f_survival(0.0, 3, 3) == 1.0
    assert f_survival(-1.0, 3, 3) == 1.0
    assert 0.0 < f_survival(1.0e6, 3, 300) < 1.0e-100
    assert f_survival(1.0e300, 3, 300) == _sys.float_info.min
    with pytest.raises(ValueError):
        f_survival(float("nan"), 2, 2)
    with pytest.raises(ValueError):
        f_survival(1.0, 0, 2)


def test_masking_fraction_is_near_one_for_a_degenerate_omission() -> None:
    """The a-priori answer to the falsifier, computed before any data exist.

    A near-unit masking fraction says the omission IS a parameter shift to within its remainder.  The
    ordering across separations is the substantive claim; it is monotone because a better-separated
    second mode has more curvature outside ``span{1, x}``.
    """
    degenerate = parameter_masking_fraction(
        omission=OmissionKind.CO_SCALING, separation=DEGENERATE_SEPARATION, n_samples=N_SAMPLES
    )
    matched = parameter_masking_fraction(
        omission=OmissionKind.CO_SCALING, separation=MATCHED_SEPARATION, n_samples=N_SAMPLES
    )
    assert 0.999 < degenerate < 1.0
    assert matched < degenerate
    assert (
        parameter_masking_fraction(omission=OmissionKind.NONE, separation=2.0, n_samples=32) == 0.0
    )
    fractions = [
        parameter_masking_fraction(
            omission=OmissionKind.CO_SCALING, separation=sep, n_samples=N_SAMPLES
        )
        for sep in (1.02, 1.1, 1.3, 2.0, 3.0)
    ]
    assert fractions == sorted(fractions, reverse=True)


def test_masking_fraction_refuses_a_non_separated_mode() -> None:
    with pytest.raises(ValueError):
        parameter_masking_fraction(
            omission=OmissionKind.CO_SCALING, separation=1.0, n_samples=N_SAMPLES
        )


def test_rate_closure_gap_reuses_the_sweep_estimator_and_agrees_with_closed_form() -> None:
    """``E[exp(-r x)] - exp(-E[r] x)`` for the two-rate mixture, checked by hand."""
    separation, x = 3.0, 1.0
    report = rate_closure_gap(separation=separation, tolerance=1.0, observation_x=x)
    weight = OMITTED_MODE_WEIGHT
    expected_mean_of_f = (1.0 - weight) * math.exp(-x) + weight * math.exp(-separation * x)
    expected_f_of_mean = math.exp(-x * ((1.0 - weight) * 1.0 + weight * separation))
    assert report.ensemble_mean_of_f == pytest.approx(expected_mean_of_f, rel=1e-12)
    assert report.f_of_ensemble_mean == pytest.approx(expected_f_of_mean, rel=1e-12)
    assert report.gap > 0.0  # Jensen: exp is convex, so the closure UNDERSTATES the mixture
    with pytest.raises(ValueError):
        rate_closure_gap(separation=1.0, tolerance=1.0, observation_x=1.0)


# --------------------------------------------------------------------------------------------
# 2. contract discipline
# --------------------------------------------------------------------------------------------


def test_contract_has_no_defaults_on_any_field() -> None:
    """A default here is a threshold this module chose instead of the caller."""
    for field in dataclasses.fields(InadequacyContract):
        assert field.default is dataclasses.MISSING, field.name
        assert field.default_factory is dataclasses.MISSING, field.name


def test_contract_rejects_a_permutation_budget_that_cannot_reach_its_own_level() -> None:
    """A channel that can never fire still consumes a share of the family-wise budget."""
    starved = dataclasses.replace(CONTRACT, permutation_trials=50)
    with pytest.raises(ValueError, match="could never fire"):
        starved.per_channel_alpha(4)
    assert CONTRACT.per_channel_alpha(4) == pytest.approx(CONTRACT.family_alpha / 4.0)


def test_contract_rejects_a_two_condition_shared_test() -> None:
    with pytest.raises(ValueError, match="vacuous"):
        dataclasses.replace(CONTRACT, min_conditions=2)
    with pytest.raises(ValueError, match="degrees of freedom"):
        dataclasses.replace(CONTRACT, min_samples_per_condition=4)


def test_channel_result_cannot_disagree_with_its_own_alpha() -> None:
    """``fired`` is derived, not a free field a caller could set to taste."""
    with pytest.raises(ValueError, match="not a free field"):
        ChannelResult(
            channel=ChannelId.CURVATURE,
            statistic_name="s",
            statistic=1.0,
            p_value=0.5,
            alpha=0.01,
            fired=True,
            null_model="n",
            detects="d",
        )


# --------------------------------------------------------------------------------------------
# 3. negative control — nothing is missing
# --------------------------------------------------------------------------------------------


def test_null_case_fires_nothing_and_returns_a_no_candidate_report() -> None:
    training, held_out = _experiment(omission=OmissionKind.NONE, separation=2.0, amplitude=0.0)
    verdict = detect_missing_law(
        training, held_out, contract=CONTRACT, experiment_id="null-control"
    )
    assert isinstance(verdict, NoCandidateReport)
    assert verdict.status is ProposalStatus.NO_CANDIDATE
    assert not verdict.evidence.any_fired
    assert "NOT evidence that no law is missing" in verdict.silence_is_not_absence


def test_measured_false_positive_rate_sits_within_the_declared_level() -> None:
    """The measurement the lane would be useless without.

    A detector that always finds a missing law finds nothing.  The trial count, the seed and the
    three-sigma binomial band are all contract fields declared before the run.
    """
    report = measure_null_false_positive_rate(
        contract=CONTRACT,
        trials=200,
        seed=20260729,
        interventions=TRAINING,
        held_out_intervention=HELD_OUT,
        noise_sigma=NOISE_SIGMA,
        n_samples=N_SAMPLES,
    )
    assert report.trials == 200
    assert report.seed == 20260729
    assert not report.assumption_violated
    assert report.false_positive_rate <= report.upper_bound, report.to_dict()
    assert report.within_declared_level
    assert set(report.per_channel_false_positives) == {item.value for item in ChannelId}


# --------------------------------------------------------------------------------------------
# 3b. RNG stream independence — the assertion the contract claimed and never made
#
# External review, finding 8: the null control asserted "no two replicates share a stream" while
# the contract had `seed == permutation_seed == 20260729`, so replicate 0's training seed equalled
# its own permutation seed, and replicate 0's held-out seed 20260730 equalled replicate 1's
# permutation seed.  New generators from the same seed are not independent streams, so the
# binomial reading of the false-positive rate rested on a premise the code did not establish.
#
# The module now spawns one stream per (replicate, role) pair.  These are the tests that check it
# — the point of the finding was that the claim existed with nothing behind it.
# --------------------------------------------------------------------------------------------


def test_no_two_replicate_role_pairs_share_a_stream() -> None:
    """Every (replicate, role) pair addresses a distinct stream.

    Asserted on the spawn KEYS rather than on sampled output: two distinct streams can of course
    agree on a finite draw by chance, so comparing samples would be a weaker check that also
    flakes.  `SeedSequence.spawn_key` is the identity of the stream itself.
    """
    streams = derive_replicate_streams(seed=20260729, trials=32)
    assert len(streams) == 32

    keys: dict[tuple[int, ...], str] = {}
    for bundle in streams:
        assert bundle.replicate == streams.index(bundle)
        for role, sequence in zip(
            REPLICATE_STREAM_ROLES,
            (bundle.training, bundle.held_out, bundle.permutation),
            strict=True,
        ):
            key = tuple(sequence.spawn_key)
            label = f"replicate {bundle.replicate} / {role}"
            assert key not in keys, f"{label} shares a stream with {keys.get(key)}"
            keys[key] = label

    assert len(keys) == 32 * len(REPLICATE_STREAM_ROLES)


def test_the_review_counterexample_two_roles_one_seed_cannot_recur() -> None:
    """The two specific collisions the review named are gone.

    The old scheme derived role seeds by ARITHMETIC on a base integer, which is why
    `seed + 1` for one replicate landed on another replicate's permutation seed.  Spawning
    cannot reproduce that: the keys are structural, not arithmetic.
    """
    streams = derive_replicate_streams(seed=20260729, trials=4)

    zero = streams[0]
    assert tuple(zero.training.spawn_key) != tuple(zero.permutation.spawn_key), (
        "replicate 0's training and permutation streams were the same in the reviewed version"
    )
    assert tuple(zero.held_out.spawn_key) != tuple(streams[1].permutation.spawn_key), (
        "replicate 0's held-out stream was replicate 1's permutation stream in the reviewed version"
    )

    # And the arithmetic that produced those collisions is not how the entropy is laid out now.
    entropies = {tuple(s.training.spawn_key) for s in streams}
    assert len(entropies) == 4


def test_the_derivation_is_recorded_so_the_ensemble_is_reproducible() -> None:
    """A stream scheme nobody can restate is not auditable; the report carries it verbatim."""
    assert "spawn" in STREAM_DERIVATION
    assert "no two (replicate, role) pairs" in STREAM_DERIVATION
    assert REPLICATE_STREAM_ROLES == ("training", "held-out", "permutation")

    first = derive_replicate_streams(seed=20260729, trials=3)
    again = derive_replicate_streams(seed=20260729, trials=3)
    for a, b in zip(first, again, strict=True):
        assert tuple(a.training.spawn_key) == tuple(b.training.spawn_key)
        assert a.training.entropy == b.training.entropy

    different = derive_replicate_streams(seed=20260730, trials=3)
    assert different[0].training.entropy != first[0].training.entropy


def test_a_malformed_stream_request_is_refused_rather_than_silently_reshaped() -> None:
    for bad in ({"seed": -1, "trials": 4}, {"seed": 20260729, "trials": 0}):
        with pytest.raises(ValueError):
            derive_replicate_streams(**bad)  # type: ignore[arg-type]


def test_a_stream_bundle_cannot_be_rebound_after_it_is_derived() -> None:
    """The bundle is frozen, so a caller cannot swap one role's stream for another's.

    Without this the independence above would hold only at the moment of derivation: reassigning
    `bundle.permutation = bundle.training` afterwards recreates exactly the collision the review
    found, and nothing downstream would notice.
    """
    bundle = derive_replicate_streams(seed=20260729, trials=2)[0]
    assert isinstance(bundle, ReplicateStreams)

    with pytest.raises(dataclasses.FrozenInstanceError):
        bundle.permutation = bundle.training  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        bundle.replicate = 99  # type: ignore[misc]


# --------------------------------------------------------------------------------------------
# 4. adversarial controls
# --------------------------------------------------------------------------------------------


def test_every_channel_is_exactly_invariant_under_rescaling_the_noise() -> None:
    """Tripling the noise with nothing missing cannot move a single p-value.

    This is the adversarial control the plan asks for, and it is an identity rather than a statistic:
    D1, D3 and D4 are ratios of sums of squares and D2 is a permutation test over a fixed multiset, so
    the whole detector is scale-free.  A magnitude-based inadequacy detector cannot make this claim.
    """
    base = analyse_conditions(
        *_experiment(omission=OmissionKind.NONE, separation=2.0, amplitude=0.0),
        contract=CONTRACT,
        experiment_id="scale-1x",
    )
    inflated = analyse_conditions(
        *_experiment(
            omission=OmissionKind.NONE, separation=2.0, amplitude=0.0, noise_sigma=3.0 * NOISE_SIGMA
        ),
        contract=CONTRACT,
        experiment_id="scale-3x",
    )
    assert inflated.post_refit_rms == pytest.approx(3.0 * base.post_refit_rms, rel=1e-9)
    for channel in ChannelId:
        assert inflated.channel(channel).p_value == pytest.approx(
            base.channel(channel).p_value, rel=1e-9
        )
    assert not inflated.any_fired


def test_heavy_tailed_noise_breaks_the_declared_level_and_the_cost_is_measured() -> None:
    """A measured deficiency, reported rather than papered over.

    The channels' null distributions are exact only under the noise model the contract declares.
    Under Student-t(3) observation noise with NOTHING missing, the family-wise rate is far above the
    declared band, and almost all of it comes from the held-out variance-ratio channel.  This is not
    fixed here by swapping in a robust statistic after the fact — that would be choosing a method
    after seeing the data.  It is measured, attributed, and routed.
    """
    report = measure_null_false_positive_rate(
        contract=CONTRACT,
        trials=150,
        seed=20260729,
        interventions=TRAINING,
        held_out_intervention=HELD_OUT,
        noise_sigma=NOISE_SIGMA,
        n_samples=N_SAMPLES,
        noise_model=NoiseModel.HEAVY_TAILED,
    )
    assert report.assumption_violated
    assert not report.within_declared_level, report.to_dict()
    assert report.inflation_over_declared > 3.0
    worst = max(report.per_channel_false_positives.items(), key=lambda item: item[1])
    assert worst[0] == ChannelId.HELD_OUT_INTERVENTION.value


@pytest.mark.parametrize("model", list(NoiseModel))
def test_every_noise_model_produces_a_usable_experiment(model: NoiseModel) -> None:
    """Each violating noise process must be constructible, or its cost could not be measured.

    ``HETEROSCEDASTIC`` ramps the standard deviation across the window and
    ``OUTLIER_CONTAMINATED`` spikes a fixed fraction of samples; both leave the log-signal's
    trend untouched, so any firing they cause is a false positive and nothing else.
    """
    training, held_out = _experiment(
        omission=OmissionKind.NONE, separation=2.0, amplitude=0.0, noise_model=model
    )
    report = analyse_conditions(
        training, held_out, contract=CONTRACT, experiment_id=f"noise/{model.value}"
    )
    assert report.n_conditions == len(TRAINING)
    assert report.post_refit_rms > 0.0
    assert all(0.0 < channel.p_value <= 1.0 for channel in report.channels)


def test_a_degenerate_omission_defeats_every_channel_the_detector_has() -> None:
    """**The falsifier, realised.**  The detector FAILED for this case, and that is the result.

    The truth carries a second relaxation mode at eight times the noise amplitude.  Because its rate
    co-scales with the intervention and its separation is small, refitting one rate absorbs it almost
    entirely: the residual falls by a large factor and every channel goes silent.  Nothing about this
    is a bug — it is what "a parameter shift hides the missing law" means, measured.
    """
    report = _degenerate_report()
    assert (
        parameter_masking_fraction(
            omission=OmissionKind.CO_SCALING,
            separation=DEGENERATE_SEPARATION,
            n_samples=N_SAMPLES,
        )
        > 0.9999
    )
    # The pre-refit residual is enormous; the post-refit residual is at the noise floor.
    assert report.refit_reduction > 5.0
    assert report.post_refit_rms == pytest.approx(NOISE_SIGMA, rel=0.15)
    assert report.fired_channels == ()
    # And the information criterion agrees with the failure: BIC prefers the WRONG model.
    assert report.channel(ChannelId.CURVATURE).detail["delta_bic_reduced_minus_extended"] < 0.0


def test_the_degenerate_omission_stays_hidden_across_independent_seeds() -> None:
    """One silent run could be luck.  The failure has to be a property of the design."""
    silent = 0
    trials = 12
    for index in range(trials):
        report = _report(
            OmissionKind.CO_SCALING,
            DEGENERATE_SEPARATION,
            DEGENERATE_AMPLITUDE,
            seed=90000 + 71 * index,
        )
        silent += int(not report.any_fired)
    assert silent >= trials - 1, f"only {silent}/{trials} seeds stayed silent"


def test_a_co_scaling_omission_is_invisible_to_the_design_channels_at_matched_amplitude() -> None:
    """The sharper statement: co-scaling defeats D3 and D4 even when D1 and D2 scream.

    Here the omission is large enough that the pointwise channels detect it easily.  The
    shared-parameter and held-out channels still cannot: an omission whose rate co-scales with the
    intervention is degenerate with the fitted rate as an algebraic identity, so more conditions and a
    held-out level buy nothing.  A programme that varies one knob cannot see a law that co-varies with
    that knob.
    """
    report = _co_scaling_report()
    assert report.channel(ChannelId.CURVATURE).fired
    assert report.channel(ChannelId.SERIAL_STRUCTURE).fired
    assert not report.channel(ChannelId.SHARED_PARAMETER).fired
    assert not report.channel(ChannelId.HELD_OUT_INTERVENTION).fired
    assert report.channel(ChannelId.SHARED_PARAMETER).detail["misfit_inflation"] < 1.01
    assert report.held_out_error_ratio < 1.2


# --------------------------------------------------------------------------------------------
# 5. positive control — an omission that cannot be absorbed
# --------------------------------------------------------------------------------------------


def test_an_absolute_time_omission_is_found_by_every_channel() -> None:
    report = _absolute_report()
    assert set(report.fired_channels) == set(ChannelId)
    assert report.channel(ChannelId.SHARED_PARAMETER).detail["misfit_inflation"] > 1.05
    assert report.held_out_error_ratio > 1.3


def test_the_discriminator_is_structure_and_not_residual_magnitude() -> None:
    """The amplitude-matched contrast, which is the lane's actual result.

    Two omissions of the SAME shape and the SAME amplitude, differing only in whether the omitted
    process responds to the intervention.  Their post-refit residuals are comparable — the co-scaling
    one is the LARGER of the two — yet only the non-co-scaling omission is found by the channels that
    the plan calls the discriminating move.  Residual magnitude is therefore not the discriminator,
    and this test is what says so.
    """
    co_scaling = _co_scaling_report()
    absolute = _absolute_report()
    assert co_scaling.post_refit_rms > absolute.post_refit_rms
    assert not co_scaling.channel(ChannelId.SHARED_PARAMETER).fired
    assert absolute.channel(ChannelId.SHARED_PARAMETER).fired
    assert not co_scaling.channel(ChannelId.HELD_OUT_INTERVENTION).fired
    assert absolute.channel(ChannelId.HELD_OUT_INTERVENTION).fired


def test_in_sample_and_held_out_errors_are_reported_as_a_pair() -> None:
    """An in-sample fit that looks fine while the held-out intervention fails."""
    absolute = _absolute_report()
    null = _null_report()
    # In-sample, the inadequate model looks only mildly worse than the correct one...
    assert absolute.in_sample_rms < 4.0 * null.in_sample_rms
    # ...and out of sample it is not.
    assert absolute.held_out_error_ratio > 1.3
    assert null.held_out_error_ratio < 1.3
    assert absolute.held_out_error_ratio > null.held_out_error_ratio


def test_permutation_null_preserves_residual_magnitude_exactly() -> None:
    """D2's null is a re-ordering, so it cannot fire because a residual is large.

    Verified directly rather than argued: a permutation is a bijection on the residual multiset, so
    the sum of squares is identical to the last bit.
    """
    rng = np.random.default_rng(7)
    residual = rng.normal(0.0, 1.0, 512)
    shuffled = rng.permutation(residual)
    assert float(shuffled @ shuffled) == pytest.approx(float(residual @ residual), rel=1e-12)
    assert sorted(shuffled.tolist()) == sorted(residual.tolist())


def test_residual_interval_uses_the_ess_estimator_and_widens_a_structured_residual() -> None:
    """A structured residual's mean has a wider honest interval than the naive one.

    Reuses :func:`~aleph.virtual_cell.sweep.ess_interval` rather than restating it: the widening
    factor is exactly the amount by which a naive standard error on a correlated residual is too
    narrow, and a detector that reports the naive interval over-calls on precisely the cases this lane
    is about.
    """
    structured = _absolute_report().residual_interval
    unstructured = _null_report().residual_interval
    assert structured.widening_factor > unstructured.widening_factor
    assert structured.widening_factor >= 1.0
    assert unstructured.widening_factor >= 1.0
    assert structured.sem_ess >= structured.sem_naive


# --------------------------------------------------------------------------------------------
# 6. the verdict may propose, never decide
# --------------------------------------------------------------------------------------------


def test_a_proposal_names_candidate_laws_and_reports_that_it_is_undecided() -> None:
    training, held_out = _experiment(
        omission=OmissionKind.ABSOLUTE_TIME,
        separation=MATCHED_SEPARATION,
        amplitude=MATCHED_AMPLITUDE,
    )
    verdict = detect_missing_law(
        training, held_out, contract=CONTRACT, experiment_id="positive-control"
    )
    assert isinstance(verdict, MissingLawProposal)
    assert verdict.status is ProposalStatus.CANDIDATE_PROPOSED
    assert verdict.decided is False
    assert verdict.evidence_source is EvidenceSource.ANALYTIC_ORACLE
    assert all(isinstance(law, MissingLaw) for law in verdict.candidate_laws)
    assert "KU-S12.shared-rate-inconsistency" in verdict.law_ids
    assert verdict.expansions


def test_a_proposal_has_no_evidence_source_field_to_overwrite() -> None:
    """The provenance is a property, so there is no argument through which to launder it."""
    names = {field.name for field in dataclasses.fields(MissingLawProposal)}
    assert "evidence_source" not in names
    assert "decided" not in names
    assert "status" not in names
    training, held_out = _experiment(
        omission=OmissionKind.ABSOLUTE_TIME,
        separation=MATCHED_SEPARATION,
        amplitude=MATCHED_AMPLITUDE,
    )
    verdict = detect_missing_law(training, held_out, contract=CONTRACT, experiment_id="p")
    with pytest.raises(TypeError):
        dataclasses.replace(verdict, evidence_source=EvidenceSource.NATIVE_ACCEPTED)
    with pytest.raises(dataclasses.FrozenInstanceError):
        verdict.proposal_id = "other"  # type: ignore[misc]


def test_a_proposal_cannot_promote_itself_to_an_archetype_verdict() -> None:
    training, held_out = _experiment(
        omission=OmissionKind.ABSOLUTE_TIME,
        separation=MATCHED_SEPARATION,
        amplitude=MATCHED_AMPLITUDE,
    )
    verdict = detect_missing_law(training, held_out, contract=CONTRACT, experiment_id="promote")
    for reference in (
        "",
        "   ",
        verdict.proposal_id,
        "sandbox_inadequacy",
        "KU-S12.unmodelled-curvature",
    ):
        with pytest.raises(SelfPromotionError):
            verdict.to_archetype_verdict(
                archetype_id="a",
                ring=ArchetypeRing.RING_B,
                adjudicated_by=reference,
                unsupported_claims=("x",),
            )
    adjudicated = verdict.to_archetype_verdict(
        archetype_id="oracle-archetype",
        ring=ArchetypeRing.RING_B,
        adjudicated_by="PI decision card 2026-07-29",
        unsupported_claims=("no cell observable may be claimed from this lane",),
    )
    assert isinstance(adjudicated, MissingLawVerdict)
    assert adjudicated.refused is True
    assert adjudicated.law_ids == verdict.law_ids


def test_a_silent_run_cannot_be_dressed_as_a_proposal_and_vice_versa() -> None:
    null = _null_report()
    absolute = _absolute_report()
    with pytest.raises(ValueError, match="no channel fired"):
        MissingLawProposal(
            proposal_id="fake",
            question="q",
            candidate_laws=(CANDIDATE_LAWS[ChannelId.CURVATURE],),
            evidence=null,
            expansions=(),
        )
    with pytest.raises(ValueError, match="belongs to a MissingLawProposal"):
        NoCandidateReport(report_id="fake", evidence=absolute)


def test_a_proposal_must_name_a_law_and_an_expansion() -> None:
    absolute = _absolute_report()
    with pytest.raises(ValueError, match="NAME at least one candidate law"):
        MissingLawProposal(
            proposal_id="p", question="q", candidate_laws=(), evidence=absolute, expansions=()
        )
    with pytest.raises(ValueError, match="expansion route"):
        MissingLawProposal(
            proposal_id="p",
            question="q",
            candidate_laws=(CANDIDATE_LAWS[ChannelId.CURVATURE],),
            evidence=absolute,
            expansions=(),
        )


def test_every_channel_has_a_pre_registered_candidate_law() -> None:
    """The evidence-to-candidate map is a pre-registration, not a reading of a residual."""
    assert set(CANDIDATE_LAWS) == set(ChannelId)
    ids = [law.law_id for law in CANDIDATE_LAWS.values()]
    assert len(set(ids)) == len(ids)
    assert set(ids) == set(S12_MANIFEST.missing_laws)


# --------------------------------------------------------------------------------------------
# 7. protocol refusals
# --------------------------------------------------------------------------------------------


def test_held_out_intervention_may_not_duplicate_a_training_level() -> None:
    training, _ = _experiment(omission=OmissionKind.NONE, separation=2.0, amplitude=0.0)
    fake_held_out = dataclasses.replace(training[0], condition_id="masquerade")
    with pytest.raises(ValueError, match="in-sample prediction wearing a held-out label"):
        analyse_conditions(training, fake_held_out, contract=CONTRACT, experiment_id="bad")


def test_too_few_conditions_or_samples_is_a_refusal_not_a_weaker_answer() -> None:
    training, held_out = _experiment(omission=OmissionKind.NONE, separation=2.0, amplitude=0.0)
    with pytest.raises(ValueError, match="at least 3 training conditions"):
        analyse_conditions(training[:2], held_out, contract=CONTRACT, experiment_id="bad")
    short = build_conditions(
        interventions=TRAINING,
        omission=OmissionKind.NONE,
        separation=2.0,
        deviation_sigma_multiple=0.0,
        noise_sigma=NOISE_SIGMA,
        n_samples=32,
        seed=1,
    )
    short_out = build_conditions(
        interventions=(HELD_OUT,),
        omission=OmissionKind.NONE,
        separation=2.0,
        deviation_sigma_multiple=0.0,
        noise_sigma=NOISE_SIGMA,
        n_samples=32,
        seed=2,
    )[0]
    with pytest.raises(ValueError, match="at least 64 samples"):
        analyse_conditions(short, short_out, contract=CONTRACT, experiment_id="bad")


def test_build_conditions_refuses_a_malformed_design() -> None:
    common = {
        "omission": OmissionKind.CO_SCALING,
        "separation": 2.0,
        "deviation_sigma_multiple": 1.0,
        "noise_sigma": NOISE_SIGMA,
        "n_samples": 64,
        "seed": 1,
    }
    with pytest.raises(ValueError, match="distinct"):
        build_conditions(interventions=(1.0, 1.0, 2.0), **common)
    with pytest.raises(ValueError, match="positive"):
        build_conditions(interventions=(1.0, 0.0), **common)
    with pytest.raises(ValueError, match="not a second process"):
        build_conditions(
            interventions=(1.0, 2.0),
            omission=OmissionKind.CO_SCALING,
            separation=1.0,
            deviation_sigma_multiple=1.0,
            noise_sigma=NOISE_SIGMA,
            n_samples=64,
            seed=1,
        )


def test_condition_series_requires_an_even_increasing_clock() -> None:
    x = np.array([0.0, 1.0, 3.0, 4.0])
    with pytest.raises(ValueError, match="evenly spaced"):
        ConditionSeries(
            condition_id="c", intervention=1.0, dimensionless_time=x, log_signal=np.zeros(4)
        )
    with pytest.raises(ValueError, match="same length"):
        ConditionSeries(
            condition_id="c",
            intervention=1.0,
            dimensionless_time=np.linspace(0, 1, 8),
            log_signal=np.zeros(4),
        )


def test_every_condition_shares_the_dimensionless_clock() -> None:
    """The protocol is what makes the shared-rate restriction exact."""
    training, held_out = _experiment(omission=OmissionKind.NONE, separation=2.0, amplitude=0.0)
    for condition in training:
        assert condition.dimensionless_time[0] == 0.0
        assert condition.dimensionless_time[-1] == pytest.approx(WINDOW_EFOLDS)
        np.testing.assert_allclose(
            condition.dimensionless_time, held_out.dimensionless_time, rtol=0, atol=0
        )


# --------------------------------------------------------------------------------------------
# 8. pre-registration and provenance
# --------------------------------------------------------------------------------------------


def test_the_sandbox_card_is_a_complete_pre_registration() -> None:
    assert isinstance(SANDBOX_S12_CARD, SandboxExperimentCard)
    assert SANDBOX_S12_CARD.manifest_hashes == (S12_MANIFEST.manifest_hash,)
    assert not SANDBOX_S12_CARD.wet_lab
    assert "parameter shift hides the missing law" in SANDBOX_S12_CARD.falsifier
    for group in (
        SANDBOX_S12_CARD.positive_controls,
        SANDBOX_S12_CARD.negative_controls,
        SANDBOX_S12_CARD.adversarial_controls,
        SANDBOX_S12_CARD.held_out_interventions,
        SANDBOX_S12_CARD.failure_expansions,
    ):
        assert group
    assert json.dumps(SANDBOX_S12_CARD.to_dict())
    assert len(SANDBOX_S12_CARD.card_hash) == 64


def test_the_lane_claims_no_biological_evidence() -> None:
    """Everything here was written down, not measured, and the graph says so by construction."""
    graph = s12_evidence_graph()
    assert graph.citable_as_biological_evidence() == ()
    assert "s12-missing-law-candidate" in graph
    assert graph.ancestors("s12-missing-law-candidate") >= {
        "s12-detector-report",
        "s12-oracle-system",
        "s12-protocol-constants",
    }
    assert len(graph.graph_hash) == 64


def test_reports_serialise_and_fingerprint_reproducibly() -> None:
    report = _absolute_report()
    payload = report.to_dict()
    assert json.dumps(payload)
    assert payload["contract"]["assumed_noise_model"] == NoiseModel.HOMOSCEDASTIC_GAUSSIAN.value
    assert payload["fired_channels"]
    assert "oracle-only" in payload["evidence_authority"]
    again = _report(OmissionKind.ABSOLUTE_TIME, MATCHED_SEPARATION, MATCHED_AMPLITUDE)
    assert experiment_fingerprint(report) == experiment_fingerprint(again)
    assert experiment_fingerprint(report) != experiment_fingerprint(_null_report())


def test_the_report_records_what_did_not_fire_as_well_as_what_did() -> None:
    """A report that only lists the positives is a report that cannot be re-judged."""
    report = _co_scaling_report()
    assert len(report.channels) == len(ChannelId)
    assert len(report.fired_channels) < len(ChannelId)
    for channel in report.channels:
        assert channel.null_model
        assert channel.detects
        assert channel.alpha == report.per_channel_alpha
    with pytest.raises(KeyError):
        dataclasses.replace(report, channels=report.channels[:1]).channel(
            ChannelId.HELD_OUT_INTERVENTION
        )

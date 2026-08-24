"""CPU oracle tests for the surrogate refusal harness.

Every target here is analytic, so the truth a gate is protecting is available in closed form.  No
cell physics, no simulation trace, no GPU: the harness is being validated, not a cell model.
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

from aleph.virtual_cell.contracts import (
    ApproximationErrorReport,
    ArtifactEnvelope,
    EvidenceSource,
    RepresentationDescriptor,
    RepresentationKind,
    UncertaintySource,
)
from aleph.virtual_cell.surrogate import (
    SIGMA_GATE,
    CoverageReport,
    ExpansionRoute,
    GateAuthorityError,
    GatedSurrogate,
    GateEvidenceLedger,
    ManifoldSupport,
    NativeQueryRequest,
    PolynomialSurrogateModel,
    SurrogateAuthorityError,
    SurrogatePrediction,
    SurrogateResponse,
    SurrogateVerdict,
    TrustRegionBox,
    evaluate_held_out_intervention,
    measure_interval_coverage,
)

# One seed for the whole file; every reported number is reproducible from it.
SEED = 20260729
NOMINAL_LEVEL = 0.95
MANIFEST_HASH = "a" * 64

# The intervention oracle is sampled on |x| <= 0.4 so a degree-5 polynomial's truncation bias sits
# far below the 0.01 observation noise. This is a model-adequacy choice, not a threshold: it makes
# the experiment measure intervention generalization instead of polynomial truncation error, and the
# calibration gate is what verifies it rather than assumes it. At |x| <= 0.6 the same fit measures
# 93.7% coverage of its nominal 95% interval and the harness correctly refuses every query.
X_HALF_RANGE = 0.4


# --------------------------------------------------------------------------------------------
# analytic oracles
# --------------------------------------------------------------------------------------------


def sine_oracle(x: np.ndarray) -> np.ndarray:
    """Bounded 1-D target: well behaved everywhere, hopeless for a polynomial far outside [-1, 1]."""
    return np.sin(2.0 * x)


def annulus_oracle(points: np.ndarray) -> np.ndarray:
    """Smooth 2-D target evaluated on an annular support."""
    return points[:, 0] * points[:, 1] + 0.5 * points[:, 0] ** 2


def intervention_oracle(points: np.ndarray) -> np.ndarray:
    """Target whose second term is dormant at every integer intervention setting.

    ``f(x, u) = sin(2x) + sin(pi u) * 3 cos(3x)``.  At ``u in {0, 1, ..., 5}`` the second term is
    identically zero, so a surrogate trained only at integer settings can interpolate the training
    distribution perfectly while having learned nothing about the mechanism the term encodes.
    """
    x = points[:, 0]
    u = points[:, 1]
    return np.sin(2.0 * x) + np.sin(np.pi * u) * 3.0 * np.cos(3.0 * x)


# --------------------------------------------------------------------------------------------
# fitted cases (built once, deterministic from SEED)
# --------------------------------------------------------------------------------------------


@cache
def _sine_case() -> tuple[GatedSurrogate, CoverageReport]:
    rng = np.random.default_rng(SEED)
    train_x = rng.uniform(-1.0, 1.0, size=(240, 1))
    train_y = sine_oracle(train_x[:, 0]) + rng.normal(0.0, 0.02, size=240)
    model = PolynomialSurrogateModel.fit(("x",), train_x, train_y, degree=7)

    check_x = rng.uniform(-1.0, 1.0, size=(4000, 1))
    truth = sine_oracle(check_x[:, 0]) + rng.normal(0.0, 0.02, size=4000)
    coverage = measure_interval_coverage(
        truth,
        model.predict(check_x),
        model.prediction_sigma(check_x),
        nominal_level=NOMINAL_LEVEL,
        seed=SEED,
        label="sine-degree-7",
    )
    surrogate = GatedSurrogate.fit(
        surrogate_id="oracle-sine",
        axes=("x",),
        inputs=train_x,
        targets=train_y,
        degree=7,
        coverage=coverage,
    )
    return surrogate, coverage


@cache
def _annulus_case() -> tuple[GatedSurrogate, np.ndarray]:
    rng = np.random.default_rng(SEED + 1)

    def _sample(n: int) -> np.ndarray:
        angle = rng.uniform(0.0, 2.0 * np.pi, size=n)
        radius = rng.uniform(0.8, 1.0, size=n)
        return np.stack([radius * np.cos(angle), radius * np.sin(angle)], axis=1)

    train = _sample(1200)
    train_y = annulus_oracle(train) + rng.normal(0.0, 0.01, size=1200)
    model = PolynomialSurrogateModel.fit(("cx", "cy"), train, train_y, degree=3)

    check = _sample(4000)
    truth = annulus_oracle(check) + rng.normal(0.0, 0.01, size=4000)
    coverage = measure_interval_coverage(
        truth,
        model.predict(check),
        model.prediction_sigma(check),
        nominal_level=NOMINAL_LEVEL,
        seed=SEED + 1,
        label="annulus-degree-3",
    )
    surrogate = GatedSurrogate.fit(
        surrogate_id="oracle-annulus",
        axes=("cx", "cy"),
        inputs=train,
        targets=train_y,
        degree=3,
        coverage=coverage,
    )
    return surrogate, train


@cache
def _intervention_case() -> tuple[GatedSurrogate, tuple[float, ...]]:
    rng = np.random.default_rng(SEED + 2)
    settings = (0.0, 1.0, 2.0, 3.0, 4.0, 5.0)
    rows = []
    for setting in settings:
        x = rng.uniform(-X_HALF_RANGE, X_HALF_RANGE, size=(200, 1))
        rows.append(np.hstack([x, np.full((200, 1), setting)]))
    train = np.vstack(rows)
    train_y = intervention_oracle(train) + rng.normal(0.0, 0.01, size=train.shape[0])
    model = PolynomialSurrogateModel.fit(("x", "u"), train, train_y, degree=5)

    check_rows = []
    for setting in settings:
        x = rng.uniform(-X_HALF_RANGE, X_HALF_RANGE, size=(600, 1))
        check_rows.append(np.hstack([x, np.full((600, 1), setting)]))
    check = np.vstack(check_rows)
    truth = intervention_oracle(check) + rng.normal(0.0, 0.01, size=check.shape[0])
    coverage = measure_interval_coverage(
        truth,
        model.predict(check),
        model.prediction_sigma(check),
        nominal_level=NOMINAL_LEVEL,
        seed=SEED + 2,
        label="intervention-degree-5",
    )
    surrogate = GatedSurrogate.fit(
        surrogate_id="oracle-intervention",
        axes=("x", "u"),
        inputs=train,
        targets=train_y,
        degree=5,
        coverage=coverage,
    )
    return surrogate, settings


def _error_report(measured_error: float = 0.02) -> ApproximationErrorReport:
    return ApproximationErrorReport(
        report_sha256="b" * 64,
        source_blob_sha256="c" * 64,
        reduced_blob_sha256="d" * 64,
        metric_id="held-out-relative-rmse",
        requested_tolerance=0.05,
        measured_error=measured_error,
    )


def _surrogate_envelope(artifact_id: str, sources: tuple[str, ...]) -> ArtifactEnvelope:
    surrogate, _ = _sine_case()
    response = surrogate.query([0.0])
    assert response.prediction is not None
    return response.prediction.as_artifact_envelope(
        artifact_id=artifact_id,
        cell_state_manifest_hash=MANIFEST_HASH,
        source_artifact_ids=sources,
        error_report=_error_report(),
    )


# --------------------------------------------------------------------------------------------
# 0. the lane's own hard constraint
# --------------------------------------------------------------------------------------------


def test_importing_the_harness_does_not_initialise_warp_or_a_learning_framework() -> None:
    """Importing and exercising the harness must pull in no GPU or learning framework.

    Run in a SUBPROCESS.  Asserting against this process's ``sys.modules`` measures whether
    anything in the whole pytest session ever imported Warp, which is a different claim and a
    false one: the suite contains tests that legitimately import it, so the in-process form
    passed alone and failed under the full run.  A fresh interpreter measures what the name
    says — what importing *this module* costs.
    """
    script = textwrap.dedent(
        """
        import sys
        import numpy as np
        from aleph.virtual_cell.surrogate import GatedSurrogate, measure_interval_coverage

        for blocked in ("warp", "torch", "jax"):
            assert blocked not in sys.modules, f"surrogate import pulled in {blocked}"
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
# 1. structural gate-authority ban
# --------------------------------------------------------------------------------------------


def test_prediction_has_no_evidence_source_field_to_overwrite() -> None:
    surrogate, _ = _sine_case()
    response = surrogate.query([0.25])
    prediction = response.prediction
    assert prediction is not None
    field_names = {field.name for field in dataclasses.fields(SurrogatePrediction)}
    assert "evidence_source" not in field_names
    assert "representation_kind" not in field_names
    assert prediction.evidence_source is EvidenceSource.SURROGATE
    assert prediction.representation_kind is RepresentationKind.SURROGATE
    with pytest.raises(TypeError):
        dataclasses.replace(prediction, evidence_source=EvidenceSource.NATIVE_ACCEPTED)
    with pytest.raises(dataclasses.FrozenInstanceError):
        prediction.mean = 0.0  # type: ignore[misc]


def test_a_subclass_may_not_redefine_the_authority_properties() -> None:
    """**Codex external review 2026-07-29, finding 4.**  The property was overridable.

    Reproduced exactly as reported: ``class Forged(SurrogatePrediction)`` redefining
    ``evidence_source`` as a property returning :data:`EvidenceSource.NATIVE_ACCEPTED`.  Before the
    fix the class was created without complaint, an instance built from the base dataclass fields
    satisfied ``isinstance(obj, SurrogatePrediction)``, and ``obj.evidence_source`` was
    ``native-accepted`` while ``obj.representation_kind`` was ``explicit-state``.  The final envelope
    was still ``surrogate`` — ``as_artifact_envelope`` re-hard-codes it — so nothing laundered, but
    the docstring's claim that the value "cannot be spelled" was false, and a claim documented more
    strongly than it holds is the exact failure this layer exists to catch.

    The refusal is at CLASS CREATION, not at call time, so a forged type cannot come into existence
    to be passed around in the first place.
    """
    for attribute in ("evidence_source", "representation_kind"):
        namespace = {attribute: property(lambda self: EvidenceSource.NATIVE_ACCEPTED)}
        with pytest.raises(TypeError, match="may not redefine"):
            type("Forged", (SurrogatePrediction,), namespace)

    # and not by assigning onto the subclass after it is created, either
    class Subclass(SurrogatePrediction):
        """A subclass that leaves the authority properties alone is still legal."""

    for attribute in ("evidence_source", "representation_kind"):
        with pytest.raises(TypeError, match="may not redefine"):
            setattr(Subclass, attribute, property(lambda self: EvidenceSource.NATIVE_ACCEPTED))

    surrogate, _ = _sine_case()
    prediction = surrogate.query([0.25]).prediction
    assert prediction is not None
    inherited = Subclass(
        **{
            field.name: getattr(prediction, field.name)
            for field in dataclasses.fields(SurrogatePrediction)
        }
    )
    assert inherited.evidence_source is EvidenceSource.SURROGATE
    assert inherited.representation_kind is RepresentationKind.SURROGATE


def test_the_envelope_is_pinned_even_if_the_property_route_were_reopened() -> None:
    """The guarantee that never depended on the property: the envelope re-derives ``SURROGATE``.

    The reviewer confirmed this half already held — the forged subclass still produced a
    ``surrogate`` envelope — and it is asserted here so the two mechanisms are known to be
    independent rather than one being mistaken for the other.
    """
    surrogate, _ = _sine_case()
    prediction = surrogate.query([0.25]).prediction
    assert prediction is not None
    envelope = prediction.as_artifact_envelope(
        artifact_id="independent-pin",
        cell_state_manifest_hash=MANIFEST_HASH,
        source_artifact_ids=("native-parent",),
        error_report=_error_report(),
    )
    assert envelope.evidence_source is EvidenceSource.SURROGATE
    assert envelope.representation.kind is RepresentationKind.SURROGATE


def test_surrogate_envelope_is_pinned_to_surrogate_provenance() -> None:
    envelope = _surrogate_envelope("surrogate-artifact", ("native-parent",))
    assert envelope.evidence_source is EvidenceSource.SURROGATE
    assert envelope.representation.kind is RepresentationKind.SURROGATE
    assert envelope.representation.trust_region
    assert envelope.representation.fallback == "native-evaluation"


def test_native_accepted_cannot_be_spelled_even_by_hand() -> None:
    descriptor = RepresentationDescriptor(
        kind=RepresentationKind.SURROGATE,
        uncertainty_sources=(UncertaintySource.STRUCTURAL,),
        axes=("x",),
        approximation_error_bound=0.02,
        error_metric="held-out-relative-rmse",
        error_report=_error_report(),
        fallback="native-evaluation",
        trust_region={"x:lower": -1.0, "x:upper": 1.0},
    )
    with pytest.raises(ValueError, match="native evidence authority"):
        ArtifactEnvelope(
            artifact_id="forged",
            cell_state_manifest_hash=MANIFEST_HASH,
            evidence_source=EvidenceSource.NATIVE_ACCEPTED,
            representation=descriptor,
            source_artifact_ids=("parent",),
        )


def test_surrogate_cannot_close_a_gate_directly() -> None:
    ledger = GateEvidenceLedger()
    ledger.register(_surrogate_envelope("surrogate-artifact", ("native-parent",)))
    with pytest.raises(SurrogateAuthorityError, match="surrogate provenance"):
        ledger.assert_gate_closable("surrogate-artifact", gate_id="GATE-B")


def test_surrogate_cannot_close_a_gate_through_an_otherwise_native_envelope() -> None:
    """The indirect route: launder the surrogate through two intermediate artifacts."""
    ledger = GateEvidenceLedger()
    ledger.register(_surrogate_envelope("surrogate-artifact", ("native-parent",)))
    ledger.register(
        ArtifactEnvelope(
            artifact_id="native-parent",
            cell_state_manifest_hash=MANIFEST_HASH,
            evidence_source=EvidenceSource.ANALYTIC_ORACLE,
            representation=RepresentationDescriptor(
                kind=RepresentationKind.EXPLICIT_STATE,
                uncertainty_sources=(UncertaintySource.NUMERICAL,),
                axes=("x",),
            ),
        )
    )
    ledger.register(
        ArtifactEnvelope(
            artifact_id="summary",
            cell_state_manifest_hash=MANIFEST_HASH,
            evidence_source=EvidenceSource.REDUCED_MODEL,
            representation=RepresentationDescriptor(
                kind=RepresentationKind.EMPIRICAL_ENSEMBLE,
                uncertainty_sources=(UncertaintySource.STRUCTURAL,),
                axes=("x",),
            ),
            source_artifact_ids=("surrogate-artifact",),
        )
    )
    ledger.register(
        ArtifactEnvelope(
            artifact_id="headline",
            cell_state_manifest_hash=MANIFEST_HASH,
            evidence_source=EvidenceSource.PUBLIC_EXPERIMENT,
            representation=RepresentationDescriptor(
                kind=RepresentationKind.EXPLICIT_STATE,
                uncertainty_sources=(UncertaintySource.MEASUREMENT,),
                axes=("x",),
            ),
            source_artifact_ids=("summary",),
        )
    )
    assert ledger.surrogate_ancestry("headline") == (
        "headline",
        "summary",
        "surrogate-artifact",
    )
    with pytest.raises(SurrogateAuthorityError) as excinfo:
        ledger.assert_gate_closable("headline", gate_id="GATE-B")
    assert "headline -> summary -> surrogate-artifact" in str(excinfo.value)


def test_a_clean_artifact_still_cannot_close_a_gate_and_the_two_refusals_differ() -> None:
    ledger = GateEvidenceLedger()
    ledger.register(
        ArtifactEnvelope(
            artifact_id="oracle-only",
            cell_state_manifest_hash=MANIFEST_HASH,
            evidence_source=EvidenceSource.ANALYTIC_ORACLE,
            representation=RepresentationDescriptor(
                kind=RepresentationKind.EXPLICIT_STATE,
                uncertainty_sources=(UncertaintySource.NUMERICAL,),
                axes=("x",),
            ),
        )
    )
    assert ledger.surrogate_ancestry("oracle-only") == ()
    with pytest.raises(GateAuthorityError, match="unverified observer") as excinfo:
        ledger.assert_gate_closable("oracle-only", gate_id="GATE-B")
    assert not isinstance(excinfo.value, SurrogateAuthorityError)


def test_unregistered_ancestry_fails_closed() -> None:
    ledger = GateEvidenceLedger()
    ledger.register(
        ArtifactEnvelope(
            artifact_id="orphan",
            cell_state_manifest_hash=MANIFEST_HASH,
            evidence_source=EvidenceSource.REDUCED_MODEL,
            representation=RepresentationDescriptor(
                kind=RepresentationKind.EMPIRICAL_ENSEMBLE,
                uncertainty_sources=(UncertaintySource.STRUCTURAL,),
                axes=("x",),
            ),
            source_artifact_ids=("ghost",),
        )
    )
    with pytest.raises(GateAuthorityError, match="unverifiable ancestry"):
        ledger.surrogate_ancestry("orphan")


# --------------------------------------------------------------------------------------------
# 2. trust region
# --------------------------------------------------------------------------------------------


def test_trust_region_is_measured_from_training_and_rejects_a_degenerate_axis() -> None:
    box = TrustRegionBox.from_training(("x", "u"), np.array([[0.0, 1.0], [2.0, 3.0]]))
    assert box.contains([0.0, 1.0]) and box.contains([2.0, 3.0])
    assert box.distance_outside([1.0, 2.0]) == 0.0
    # x=3 is one box-width (extent 2) beyond the upper edge 2 -> normalized excess 0.5.
    assert box.worst_axis([3.0, 3.0]) == ("x", 0.5)
    with pytest.raises(ValueError, match="strictly positive extent"):
        TrustRegionBox.from_training(("x",), np.array([[1.0], [1.0]]))


def test_trust_region_refuses_the_extrapolation_that_would_have_been_badly_wrong() -> None:
    surrogate, _ = _sine_case()
    inside = surrogate.query([0.3])
    assert inside.verdict is SurrogateVerdict.ANSWERED
    assert abs(inside.value - float(sine_oracle(np.array([0.3]))[0])) < 0.02

    far = 6.0
    outside = surrogate.query([far])
    assert outside.verdict is SurrogateVerdict.REFUSED_OUTSIDE_TRUST_REGION
    assert outside.prediction is None
    with pytest.raises(SurrogateAuthorityError, match="refusal has no value"):
        _ = outside.value

    # The contrast is the point: record what the refusal prevented.
    would_have_said = surrogate.ungated_prediction_for_diagnosis([far])
    truth = float(sine_oracle(np.array([far]))[0])
    assert abs(truth) <= 1.0
    assert abs(would_have_said - truth) > 1.0e3


def test_trust_region_refusal_names_the_axis_and_the_distance() -> None:
    surrogate, _ = _sine_case()
    request = surrogate.query([6.0]).native_request
    assert request is not None
    assert request.axes == ("x",)
    assert request.query_point == (6.0,)
    assert "'x'" in request.rationale
    assert request.diagnostic["normalized_excess"] > 2.0
    assert ExpansionRoute.NATIVE_EVALUATION in request.expansions


# --------------------------------------------------------------------------------------------
# 3. off-manifold rejection a bounding box would miss
# --------------------------------------------------------------------------------------------


def test_bounding_box_accepts_the_hole_that_the_support_detector_rejects() -> None:
    surrogate, train = _annulus_case()
    hole = [0.0, 0.0]

    # A box trust region alone would happily answer here.
    assert surrogate.box.contains(hole)
    assert surrogate.box.distance_outside(hole) == 0.0

    response = surrogate.query(hole)
    assert response.verdict is SurrogateVerdict.REFUSED_OFF_MANIFOLD
    request = response.native_request
    assert request is not None
    assert request.diagnostic["normalized_excess"] == 0.0
    assert request.diagnostic["neighbour_distance"] > request.diagnostic["covering_radius"]
    assert ExpansionRoute.ADAPTIVE_DESIGN in request.expansions

    on_manifold = (train[0] * 1.001).tolist()
    assert surrogate.query(on_manifold).verdict is SurrogateVerdict.ANSWERED


def test_covering_radius_is_a_statistic_of_the_training_set_not_a_threshold() -> None:
    box = TrustRegionBox.from_training(("x",), np.array([[0.0], [1.0], [2.0], [10.0]]))
    support = ManifoldSupport.from_training(box, np.array([[0.0], [1.0], [2.0], [10.0]]))
    # Largest leave-one-out nearest-neighbour gap is 10 -> 2, i.e. 8 of the 10-wide box.
    assert support.covering_radius == pytest.approx(0.8)
    assert support.supports([5.0])
    assert not support.supports([100.0])


# --------------------------------------------------------------------------------------------
# 4. posterior coverage / calibration
# --------------------------------------------------------------------------------------------


def test_coverage_measurement_separates_a_calibrated_from_a_miscalibrated_model() -> None:
    rng = np.random.default_rng(SEED)
    n_draws = 20_000
    sigma_true = 0.35
    means = rng.uniform(-1.0, 1.0, size=n_draws)
    truth = means + rng.normal(0.0, sigma_true, size=n_draws)

    calibrated = measure_interval_coverage(
        truth,
        means,
        np.full(n_draws, sigma_true),
        nominal_level=NOMINAL_LEVEL,
        seed=SEED,
        label="calibrated-gaussian",
    )
    miscalibrated = measure_interval_coverage(
        truth,
        means,
        np.full(n_draws, 0.3 * sigma_true),
        nominal_level=NOMINAL_LEVEL,
        seed=SEED,
        label="sigma-scaled-0.3",
    )

    assert calibrated.calibrated
    assert abs(calibrated.empirical_coverage - NOMINAL_LEVEL) <= SIGMA_GATE * (
        calibrated.standard_error
    )

    # Closed-form expectation for a Gaussian whose reported sigma is 0.3x the truth.
    from statistics import NormalDist

    half_width = NormalDist().inv_cdf(0.5 + NOMINAL_LEVEL / 2.0)
    expected = 2.0 * NormalDist().cdf(half_width * 0.3) - 1.0
    assert not miscalibrated.calibrated
    assert miscalibrated.empirical_coverage == pytest.approx(expected, abs=0.02)
    assert miscalibrated.z_score < -SIGMA_GATE
    assert miscalibrated.to_dict()["seed"] == SEED


def test_the_fitted_oracle_surrogate_is_itself_calibrated_inside_its_trust_region() -> None:
    _, coverage = _sine_case()
    assert coverage.calibrated
    assert coverage.n_draws == 4000
    assert coverage.seed == SEED
    assert 0.0 <= coverage.empirical_coverage <= 1.0


def test_an_uncalibrated_surrogate_refuses_every_query_including_a_training_point() -> None:
    rng = np.random.default_rng(SEED)
    train_x = rng.uniform(-1.0, 1.0, size=(240, 1))
    train_y = sine_oracle(train_x[:, 0]) + rng.normal(0.0, 0.02, size=240)
    bad_coverage = CoverageReport(
        label="sigma-scaled-0.3",
        nominal_level=NOMINAL_LEVEL,
        empirical_coverage=0.44,
        n_draws=20_000,
        standard_error=math.sqrt(NOMINAL_LEVEL * (1.0 - NOMINAL_LEVEL) / 20_000),
        z_score=(0.44 - NOMINAL_LEVEL) / math.sqrt(NOMINAL_LEVEL * 0.05 / 20_000),
        sigma_gate=SIGMA_GATE,
        seed=SEED,
    )
    assert not bad_coverage.calibrated
    surrogate = GatedSurrogate.fit(
        surrogate_id="oracle-sine-uncalibrated",
        axes=("x",),
        inputs=train_x,
        targets=train_y,
        degree=7,
        coverage=bad_coverage,
    )
    response = surrogate.query([float(train_x[0, 0])])
    assert response.verdict is SurrogateVerdict.REFUSED_UNCALIBRATED
    assert response.native_request is not None
    assert ExpansionRoute.SOLVER_PRECONDITIONER in response.native_request.expansions


def test_a_prediction_cannot_be_built_on_an_uncalibrated_coverage_report() -> None:
    surrogate, coverage = _sine_case()
    bad = dataclasses.replace(coverage, z_score=-40.0)
    with pytest.raises(ValueError, match="uncalibrated surrogate"):
        SurrogatePrediction(
            surrogate_id="x",
            axes=("x",),
            query_point=(0.0,),
            mean=0.0,
            predictive_sigma=surrogate.model.residual_sigma,
            trust_region=surrogate.box.summary(),
            support_distance=0.0,
            coverage=bad,
        )


# --------------------------------------------------------------------------------------------
# 5. held-out intervention
# --------------------------------------------------------------------------------------------


def test_small_interpolation_error_does_not_survive_a_held_out_intervention() -> None:
    surrogate, settings = _intervention_case()
    rng = np.random.default_rng(SEED + 3)

    interp_rows = [
        np.hstack(
            [
                rng.uniform(-X_HALF_RANGE, X_HALF_RANGE, size=(400, 1)),
                np.full((400, 1), setting),
            ]
        )
        for setting in settings
    ]
    interp = np.vstack(interp_rows)
    held_out_x = rng.uniform(-X_HALF_RANGE, X_HALF_RANGE, size=(2400, 1))
    held_out = np.hstack([held_out_x, np.full((2400, 1), 0.5)])

    report = evaluate_held_out_intervention(
        surrogate,
        intervention_axis="u",
        training_interventions=settings,
        held_out_intervention=0.5,
        interpolation_inputs=interp,
        interpolation_truth=intervention_oracle(interp),
        held_out_inputs=held_out,
        held_out_truth=intervention_oracle(held_out),
        seed=SEED + 3,
    )

    assert report.interpolation_rmse < 0.01
    assert report.held_out_rmse > 1.0
    assert report.degradation_ratio > 100.0
    assert report.mechanism_falsified
    assert report.native_request is not None
    assert report.native_request.reason is SurrogateVerdict.REFUSED_INTERVENTION_UNSUPPORTED
    assert ExpansionRoute.MODEL_INADEQUACY in report.native_request.expansions
    assert report.to_dict()["seed"] == SEED + 3


def test_the_held_out_setting_is_inside_the_box_and_only_the_manifold_test_catches_it() -> None:
    surrogate, _ = _intervention_case()
    probe = [0.25, 0.5]
    assert surrogate.box.contains(probe)
    response = surrogate.query(probe)
    assert response.verdict is SurrogateVerdict.REFUSED_OFF_MANIFOLD


def test_the_held_out_intervention_must_not_appear_in_training() -> None:
    surrogate, settings = _intervention_case()
    sample = np.array([[0.0, 1.0]])
    with pytest.raises(ValueError, match="must not appear in training"):
        evaluate_held_out_intervention(
            surrogate,
            intervention_axis="u",
            training_interventions=settings,
            held_out_intervention=1.0,
            interpolation_inputs=sample,
            interpolation_truth=intervention_oracle(sample),
            held_out_inputs=sample,
            held_out_truth=intervention_oracle(sample),
            seed=SEED,
        )


# --------------------------------------------------------------------------------------------
# 6. native-query acquisition
# --------------------------------------------------------------------------------------------


def test_every_refusal_emits_a_typed_reproducible_native_request() -> None:
    sine, _ = _sine_case()
    annulus, _ = _annulus_case()
    responses = [sine.query([6.0]), annulus.query([0.0, 0.0])]
    for response in responses:
        assert not response.answered
        request = response.native_request
        assert isinstance(request, NativeQueryRequest)
        assert request.reason is response.verdict
        assert request.expansions
        assert request.rationale
        assert len(request.request_hash) == 64
        assert request.to_dict()["reason"] == str(response.verdict)


def test_native_request_hash_is_content_addressed() -> None:
    sine, _ = _sine_case()
    first = sine.query([6.0]).native_request
    second = sine.query([6.0]).native_request
    third = sine.query([7.0]).native_request
    assert first is not None and second is not None and third is not None
    assert first.request_hash == second.request_hash
    assert first.request_hash != third.request_hash


def test_an_answered_response_carries_no_native_request_and_a_refusal_carries_no_number() -> None:
    surrogate, _ = _sine_case()
    answered = surrogate.query([0.0])
    assert answered.answered
    assert answered.native_request is None
    assert answered.prediction is not None
    assert answered.prediction.predictive_sigma > 0.0

    with pytest.raises(ValueError, match="refusal carries a native request"):
        SurrogateResponse(verdict=SurrogateVerdict.REFUSED_OFF_MANIFOLD)
    with pytest.raises(ValueError, match="answered response carries a prediction"):
        SurrogateResponse(verdict=SurrogateVerdict.ANSWERED)


def test_a_refusal_cannot_be_recorded_as_an_answer() -> None:
    with pytest.raises(ValueError, match="does not request a native evaluation"):
        NativeQueryRequest(
            request_id="r",
            reason=SurrogateVerdict.ANSWERED,
            rationale="none",
            axes=("x",),
            query_point=(0.0,),
            diagnostic={},
            expansions=(ExpansionRoute.NATIVE_EVALUATION,),
        )
    with pytest.raises(ValueError, match="at least one ExpansionRoute"):
        NativeQueryRequest(
            request_id="r",
            reason=SurrogateVerdict.REFUSED_OFF_MANIFOLD,
            rationale="none",
            axes=("x",),
            query_point=(0.0,),
            diagnostic={},
            expansions=(),
        )


# --------------------------------------------------------------------------------------------
# numerical guards
# --------------------------------------------------------------------------------------------


def test_a_rank_deficient_design_raises_instead_of_reporting_a_minimum_norm_fit() -> None:
    rng = np.random.default_rng(SEED)
    x = rng.uniform(-1.0, 1.0, size=(200, 1))
    inputs = np.hstack([x, np.zeros((200, 1))])
    targets = sine_oracle(x[:, 0])
    with pytest.raises(ValueError, match="rank deficient"):
        PolynomialSurrogateModel.fit(("x", "u"), inputs, targets, degree=3)


def test_a_fit_without_residual_degrees_of_freedom_raises() -> None:
    inputs = np.array([[0.0], [1.0], [2.0]])
    targets = np.array([0.0, 1.0, 4.0])
    with pytest.raises(ValueError, match="residual degrees of freedom"):
        PolynomialSurrogateModel.fit(("x",), inputs, targets, degree=2)


def test_non_finite_inputs_are_refused() -> None:
    surrogate, _ = _sine_case()
    with pytest.raises(ValueError, match="must be finite"):
        surrogate.query([float("nan")])
    with pytest.raises(ValueError, match="must be finite"):
        surrogate.query([float("inf")])

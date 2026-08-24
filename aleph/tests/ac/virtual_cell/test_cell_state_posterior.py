"""CPU oracle tests for `CellStatePosterior` — the object master plan section 0 declares as central.

The numbers checked here are closed forms, not samples.  The A1 counterexample is arithmetic on two
Gaussians, the Gaussian quantile is checked against :class:`statistics.NormalDist`, and everything else
is a structural refusal.  What is validated is the ASSEMBLY: that the eleven existing lanes compose
into one object, and that the object refuses the things the audit register says must be refused.

No cell physics is claimed anywhere.  Nothing here runs on a GPU, and no verdict below grants an
evidence rung.
"""

from __future__ import annotations

import dataclasses
import math
import subprocess
import sys
import textwrap
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pytest

from aleph.virtual_cell.admissibility import (
    AdmissibilityVerdict,
    AttemptedStepRecord,
    ConnectorForcePair,
    StateDigestSet,
    evaluate_admissibility,
)
from aleph.virtual_cell.archetypes import (
    ARCHETYPE_IDS,
    ArchetypeBuild,
    PriorStatus,
    build_archetype,
)
from aleph.virtual_cell.cell_state_posterior import (
    POINT_REPRESENTATION,
    POSTERIOR_REPRESENTATION_KINDS,
    CellStatePosterior,
    GaussianCollapse,
    IdentifiabilityDeclaration,
    InterventionResponse,
    LossyMergeRefused,
    MolecularParameterPosterior,
    PosteriorEvidenceContract,
    PosteriorGraph,
    UncertaintyBudget,
    UncertaintyCombination,
    UncertaintyComponent,
)
from aleph.virtual_cell.contracts import (
    ApproximationErrorReport,
    ArtifactEnvelope,
    EvidenceSource,
    RepresentationDescriptor,
    RepresentationKind,
    SandboxExperimentCard,
    UncertaintySource,
)
from aleph.virtual_cell.grammar import ConnectorInstance, EngineGrammar, load_engine_grammar
from aleph.virtual_cell.observation_operator import (
    EMITTER_KIND,
    MechanisticState,
    ObservationModality,
    ObservationOperator,
    ObservationRefused,
    OperatorRegistrationRefused,
    SummaryProjectionOperator,
    unimplemented_operator,
)
from aleph.virtual_cell.regime import (
    EVIDENCE_AUTHORITY,
    Regime,
    RegimeEvidenceContract,
    RegimeReport,
    classify_regime,
)
from aleph.virtual_cell.sandbox_inadequacy import (
    InadequacyContract,
    MissingLawProposal,
    NoiseModel,
    OmissionKind,
    build_conditions,
    detect_missing_law,
)
from aleph.virtual_cell.sandbox_inverse import (
    PARAMETER_NAMES,
    REFERENCE_SCENE,
    REFERENCE_TRUTH,
    SINGLE_PLANE_DESIGN,
    UnidentifiableDirectionError,
    identifiability_report,
)
from aleph.virtual_cell.source_registry import (
    EvidenceGraph,
    ProvenanceKind,
    SourceOrigin,
    SourceRecord,
    SourceRegistry,
    SourceStatus,
)
from aleph.virtual_cell.telemetry import ScalarProjectionSpec, ScalarReduction

# --------------------------------------------------------------------------------------------
# declared thresholds — written before any measurement below
# --------------------------------------------------------------------------------------------

#: The eigenvalue floor and identifiability requirement this test suite declares.  The floor is an
#: absolute Fisher information; 1e-3 is far below the smallest informative eigenvalue of the S11
#: single-plane design (0.886) and far above its numerical null (2.5e-12), so it separates "measured"
#: from "not measured" without sitting on either.
EIGENVALUE_FLOOR = 1e-3
MIN_IDENTIFIABLE_FRACTION = 0.99

CONTRACT = PosteriorEvidenceContract(
    min_identifiable_fraction=MIN_IDENTIFIABLE_FRACTION,
    eigenvalue_floor=EIGENVALUE_FLOOR,
    max_relative_fourth_moment_error=0.05,
    credible_level=0.95,
    max_representation_error=0.1,
)

#: The A1 counterexample, verbatim from the audit register: two zero-mean Gaussian sources of very
#: different width, mixed with equal weight.
SIGMA_THERMAL = 0.1
SIGMA_BIOLOGICAL = 1.0

#: Lane C's declared regime contract, reused verbatim rather than re-tuned.  Nothing here is chosen
#: to make a classification come out a particular way; the fixture is an analytic sine.
REGIME_CONTRACT = RegimeEvidenceContract(
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

INADEQUACY_CONTRACT = InadequacyContract(
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


# --------------------------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------------------------


def _component(
    source: UncertaintySource,
    variance: float,
    combination: UncertaintyCombination = UncertaintyCombination.CONVOLUTION,
    **extra: object,
) -> UncertaintyComponent:
    return UncertaintyComponent(
        source=source,
        variance=variance,
        combination=combination,
        propagation="analytic fixture",
        evidence_id="fixture",
        **extra,  # type: ignore[arg-type]
    )


def _a1_budget() -> UncertaintyBudget:
    """The registered A1 counterexample as a budget: thermal and biological, mixed equally."""
    return UncertaintyBudget(
        parameter="k_off",
        unit="1/s",
        components=(
            _component(
                UncertaintySource.THERMAL,
                SIGMA_THERMAL**2,
                UncertaintyCombination.MIXTURE,
                weight=0.5,
            ),
            _component(
                UncertaintySource.BIOLOGICAL,
                SIGMA_BIOLOGICAL**2,
                UncertaintyCombination.MIXTURE,
                weight=0.5,
            ),
        ),
    )


def _convolved_budget() -> UncertaintyBudget:
    """The same two widths, but combining by convolution instead of as a mixture."""
    return UncertaintyBudget(
        parameter="k_off",
        unit="1/s",
        components=(
            _component(UncertaintySource.THERMAL, SIGMA_THERMAL**2),
            _component(UncertaintySource.BIOLOGICAL, SIGMA_BIOLOGICAL**2),
        ),
    )


def _point_budget(parameter: str = "k_off", unit: str = "1/s") -> UncertaintyBudget:
    """What today's deterministic engine actually supports: sources named, none propagated."""
    return UncertaintyBudget(
        parameter=parameter,
        unit=unit,
        components=(
            UncertaintyComponent(
                source=source,
                variance=0.0,
                combination=UncertaintyCombination.CONVOLUTION,
                propagation="not propagated: the accepted-step runtime carries one trajectory",
                evidence_id="engine/point-trajectory",
            )
            for source in (
                UncertaintySource.THERMAL,
                UncertaintySource.STOCHASTIC_EVENT,
                UncertaintySource.PARAMETER,
            )
        ),
    )


def _point_parameters(**overrides: object) -> MolecularParameterPosterior:
    kwargs: dict[str, object] = {
        "parameters": ("k_off",),
        "units": {"k_off": "1/s"},
        "location": {"k_off": 1.4},
        "budgets": {"k_off": _point_budget()},
        "representation": RepresentationDescriptor(
            kind=POINT_REPRESENTATION,
            uncertainty_sources=(UncertaintySource.PARAMETER,),
            axes=("k_off",),
        ),
    }
    kwargs.update(overrides)
    return MolecularParameterPosterior(**kwargs)  # type: ignore[arg-type]


def _evidence_graph() -> EvidenceGraph:
    registry = SourceRegistry()
    registry.register(
        SourceRecord(
            source_id="lit_koff",
            provenance=ProvenanceKind.SOURCED,
            origin=SourceOrigin.LITERATURE,
            status=SourceStatus.ACTIVE,
            description="published off-rate measurement",
            citation="doi:10.0000/fixture",
        )
    )
    registry.register(
        SourceRecord(
            source_id="posterior_draft",
            provenance=ProvenanceKind.INFERRED_POSTERIOR,
            origin=SourceOrigin.COMPUTATION,
            status=SourceStatus.ACTIVE,
            description="posterior inferred from the literature prior",
            parents=("lit_koff",),
        )
    )
    return registry.freeze()


def _grammar() -> EngineGrammar:
    return load_engine_grammar()


def _ring_a_build() -> ArchetypeBuild:
    return build_archetype(ARCHETYPE_IDS[0], grammar=_grammar())


def _admissible_verdict(*, admissible: bool = True, step_index: int = 7) -> AdmissibilityVerdict:
    digest = StateDigestSet.from_channels(
        state={"x": np.zeros(3)},
        rng={"seed": np.array([7])},
        topology={"pairs": np.array([0])},
        clock={"t": np.array([0.0])},
    )
    record = AttemptedStepRecord(
        attempted_step_index=step_index,
        attempted_dt_s=1e-4,
        connector_pairs=(
            ConnectorForcePair(
                connector_id="cortex-membrane",
                force_a_pn=np.array([1.0, 0.0, 0.0]),
                force_b_pn=np.array([-1.0, 0.0, 0.0]),
            ),
        ),
        inner_residual_pn=1e-4 if admissible else 1.0,
        inner_tolerance_pn=1e-3,
        inner_iterations=12,
        inner_iteration_budget=50,
        candidate_state={"x": np.zeros(3)},
        digests_before=digest,
        digests_after=digest,
        rollback_digests=None if admissible else digest,
    )
    return evaluate_admissibility(record, adjoint_rel_tol=1e-6, adjoint_abs_floor_pn=1e-9)


def _limit_cycle_regime(frequency_hz: float = 1.5) -> RegimeReport:
    dt = 0.002
    times = np.arange(4096) * dt
    return classify_regime(
        np.sin(2.0 * np.pi * frequency_hz * times),
        dt,
        observable="cortex_tension",
        amplitude_scale=1.0,
        contract=REGIME_CONTRACT,
    )


def _missing_law_proposal() -> MissingLawProposal:
    training = build_conditions(
        interventions=(0.6, 1.0, 1.6),
        omission=OmissionKind.ABSOLUTE_TIME,
        separation=3.0,
        deviation_sigma_multiple=4.0,
        noise_sigma=0.02,
        n_samples=256,
        seed=4242,
    )
    held_out = build_conditions(
        interventions=(3.2,),
        omission=OmissionKind.ABSOLUTE_TIME,
        separation=3.0,
        deviation_sigma_multiple=4.0,
        noise_sigma=0.02,
        n_samples=256,
        seed=4243,
    )[0]
    proposal = detect_missing_law(
        training, held_out, contract=INADEQUACY_CONTRACT, experiment_id="s12/posterior-fixture"
    )
    assert isinstance(proposal, MissingLawProposal)
    return proposal


def _posterior(**overrides: object) -> CellStatePosterior:
    build = _ring_a_build()
    manifest = build.manifest
    assert manifest is not None
    parameters = overrides.pop("parameters", _point_parameters())
    kwargs: dict[str, object] = {
        "posterior_id": "posterior/fixture",
        "manifest": manifest,
        "graph": PosteriorGraph(grammar=_grammar()),
        "priors": build.priors,
        "parameters": parameters,
        "contract": CONTRACT,
        "evidence": _evidence_graph(),
        "envelope": ArtifactEnvelope(
            artifact_id="artifact/posterior-fixture",
            cell_state_manifest_hash=manifest.manifest_hash,
            evidence_source=EvidenceSource.REDUCED_MODEL,
            representation=parameters.representation,  # type: ignore[union-attr]
            source_artifact_ids=("artifact/accepted-run",),
        ),
    }
    kwargs.update(overrides)
    return CellStatePosterior(**kwargs)  # type: ignore[arg-type]


# --------------------------------------------------------------------------------------------
# 0. the module stays on the CPU
# --------------------------------------------------------------------------------------------


def test_importing_the_module_pulls_in_no_gpu_stack() -> None:
    """Measured in a FRESH interpreter; an in-process check would measure the whole session."""
    script = textwrap.dedent(
        """
        import sys
        from aleph.virtual_cell.cell_state_posterior import (
            CellStatePosterior,
            UncertaintyBudget,
        )

        for blocked in ("warp", "torch", "jax"):
            assert blocked not in sys.modules, f"cell_state_posterior import pulled in {blocked}"
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
# 1. A1 — uncertainty cannot be collapsed into one Gaussian without a measured cost
# --------------------------------------------------------------------------------------------


def test_the_registered_a1_counterexample_reproduces_exactly() -> None:
    """Closed form, so the numbers are checked against truth rather than against a sample.

    The moment-matched Gaussian gets the variance EXACTLY right — which is why the second moment
    cannot detect the merge — and misreports the fourth moment by 1.9608.
    """
    budget = _a1_budget()
    assert budget.total_variance() == pytest.approx(0.505, abs=1e-12)
    assert budget.fourth_central_moment() == pytest.approx(1.50015, abs=1e-12)
    assert budget.fourth_moment_ratio() == pytest.approx(1.9607881580237234, rel=1e-12)
    # The value the audit register stores for A1, reproduced here rather than quoted.
    assert budget.fourth_moment_ratio() == pytest.approx(
        (0.5 * 3.0 * SIGMA_THERMAL**4 + 0.5 * 3.0 * SIGMA_BIOLOGICAL**4)
        / (3.0 * (0.5 * SIGMA_THERMAL**2 + 0.5 * SIGMA_BIOLOGICAL**2) ** 2),
        rel=1e-12,
    )
    assert budget.excess_kurtosis() > 2.8


def test_collapsing_the_a1_budget_is_refused_and_the_refusal_carries_the_measurement() -> None:
    with pytest.raises(LossyMergeRefused) as excinfo:
        _a1_budget().collapse_to_gaussian(max_relative_fourth_moment_error=0.05)
    error = excinfo.value
    assert error.measured_error == pytest.approx(0.9607881580237234, rel=1e-12)
    assert error.total_variance == pytest.approx(0.505, abs=1e-12)
    assert "GAUSSIAN_MIXTURE" in str(error)


def test_the_same_two_widths_collapse_EXACTLY_when_they_convolve_instead_of_mixing() -> None:
    """The informative contrast: merging is not always wrong, it is wrong in a computable way.

    Independent Gaussian sources convolve to a Gaussian, so the moment-matched collapse is exact and
    the measured error is 0.  It is the MIXTURE — cell-to-cell variability, a quenched draw — that the
    single sigma destroys.  Same sigmas, same total variance, opposite verdicts.
    """
    convolved = _convolved_budget()
    assert convolved.total_variance() == pytest.approx(1.01, abs=1e-12)
    assert convolved.fourth_moment_ratio() == pytest.approx(1.0, rel=1e-12)
    collapse = convolved.collapse_to_gaussian(max_relative_fourth_moment_error=1e-9)
    assert collapse.error_report.measured_error == pytest.approx(0.0, abs=1e-15)
    assert collapse.lossy is True


def test_a_collapse_keeps_the_decomposition_it_summarised() -> None:
    collapse = _a1_budget().collapse_to_gaussian(max_relative_fourth_moment_error=0.99)
    assert isinstance(collapse, GaussianCollapse)
    assert collapse.decomposition == {
        UncertaintySource.THERMAL: SIGMA_THERMAL**2,
        UncertaintySource.BIOLOGICAL: SIGMA_BIOLOGICAL**2,
    }
    assert collapse.variance == pytest.approx(_a1_budget().total_variance(), rel=1e-15)
    assert isinstance(collapse.error_report, ApproximationErrorReport)
    assert collapse.error_report.metric_id == "relative-fourth-moment-error"


def test_a_collapse_cannot_be_marked_exact() -> None:
    collapse = _convolved_budget().collapse_to_gaussian(max_relative_fourth_moment_error=1e-9)
    assert "lossy" not in {field.name for field in dataclasses.fields(collapse)}
    with pytest.raises(TypeError):
        dataclasses.replace(collapse, lossy=False)  # type: ignore[call-arg]


# --------------------------------------------------------------------------------------------
# 2. there is structurally no single-covariance constructor
# --------------------------------------------------------------------------------------------


def test_a_budget_can_only_be_built_from_named_sources() -> None:
    assert {field.name for field in dataclasses.fields(UncertaintyBudget)} == {
        "parameter",
        "unit",
        "components",
    }
    with pytest.raises(ValueError, match="names no source"):
        UncertaintyBudget(parameter="k", unit="1/s", components=())


def test_the_source_vocabulary_has_no_total_member_so_a_lumped_width_has_no_name() -> None:
    assert "TOTAL" not in UncertaintySource.__members__
    assert {source.value for source in UncertaintySource} == {
        "thermal",
        "stochastic-event",
        "structural",
        "biological",
        "measurement",
        "parameter",
        "numerical",
        "representation",
    }


def test_a_total_is_computed_on_demand_and_cannot_be_stored() -> None:
    budget = _a1_budget()
    assert callable(budget.total_variance)
    with pytest.raises(TypeError):
        dataclasses.replace(budget, total_variance=1.0)  # type: ignore[call-arg]
    # and the decomposition survives every total that is asked for
    budget.total_variance()
    assert budget.variance_by_source() == {
        UncertaintySource.THERMAL: SIGMA_THERMAL**2,
        UncertaintySource.BIOLOGICAL: SIGMA_BIOLOGICAL**2,
    }


def test_one_source_cannot_appear_twice() -> None:
    with pytest.raises(ValueError, match="merge wearing a decomposition"):
        UncertaintyBudget(
            parameter="k",
            unit="1/s",
            components=(
                _component(UncertaintySource.THERMAL, 1.0),
                _component(UncertaintySource.THERMAL, 2.0),
            ),
        )


def test_mixture_weights_must_be_declared_and_must_sum_to_one() -> None:
    with pytest.raises(ValueError, match="must declare its weight"):
        _component(UncertaintySource.THERMAL, 1.0, UncertaintyCombination.MIXTURE)
    with pytest.raises(ValueError, match="cannot carry a mixture weight"):
        _component(UncertaintySource.THERMAL, 1.0, weight=0.5)
    with pytest.raises(ValueError, match="sum to 1"):
        UncertaintyBudget(
            parameter="k",
            unit="1/s",
            components=(
                _component(
                    UncertaintySource.THERMAL, 1.0, UncertaintyCombination.MIXTURE, weight=0.3
                ),
                _component(
                    UncertaintySource.BIOLOGICAL, 1.0, UncertaintyCombination.MIXTURE, weight=0.3
                ),
            ),
        )


def test_a_named_but_unpropagated_source_is_visible_rather_than_absent() -> None:
    """'0.0 variance' and 'no entry' are different claims, and only one of them is auditable."""
    budget = _point_budget()
    assert budget.total_variance() == 0.0
    assert set(budget.unpropagated_sources) == {
        UncertaintySource.THERMAL,
        UncertaintySource.STOCHASTIC_EVENT,
        UncertaintySource.PARAMETER,
    }
    assert all(not component.propagated for component in budget.components)


# --------------------------------------------------------------------------------------------
# 3. the credible interval is exact, so an honest interval never routes through the collapse
# --------------------------------------------------------------------------------------------


def test_the_gaussian_interval_matches_the_closed_form_inverse_cdf() -> None:
    budget = UncertaintyBudget(
        parameter="k",
        unit="1/s",
        components=(_component(UncertaintySource.THERMAL, 4.0),),
    )
    low, high = budget.central_interval(level=0.95)
    expected = NormalDist().inv_cdf(0.975) * 2.0
    assert high == pytest.approx(expected, rel=1e-12)
    assert low == pytest.approx(-expected, rel=1e-12)


def test_the_mixture_interval_is_measurably_different_from_the_moment_matched_gaussian() -> None:
    """The number this measures is what the single-sigma collapse would have got wrong.

    For the registered A1 pair the exact 95 % half-width is 1.64485 while the moment-matched Gaussian
    reports 1.39282 — the collapse is 15.3 % too NARROW, in the direction that over-states precision.
    """
    low, high = _a1_budget().central_interval(level=0.95)
    exact_half = 0.5 * (high - low)
    gaussian_half = NormalDist().inv_cdf(0.975) * math.sqrt(0.505)
    assert exact_half == pytest.approx(1.6448536269514722, rel=1e-9)
    assert gaussian_half == pytest.approx(1.3928161057550035, rel=1e-12)
    assert exact_half / gaussian_half == pytest.approx(1.18096, rel=1e-4)


def test_the_exact_interval_carries_the_probability_it_claims() -> None:
    budget = _a1_budget()
    low, high = budget.central_interval(level=0.9)
    assert budget.cdf(high) - budget.cdf(low) == pytest.approx(0.9, abs=1e-12)


def test_a_point_budget_gives_a_degenerate_interval_rather_than_an_error() -> None:
    """The visibly useless answer IS the honest report for a point representation."""
    low, high = _point_budget().central_interval(level=0.95)
    assert low == 0.0
    assert high == 0.0


# --------------------------------------------------------------------------------------------
# 4. A5 — no point estimate along a direction that was not measured
# --------------------------------------------------------------------------------------------


def _image_parameters() -> MolecularParameterPosterior:
    report = identifiability_report(
        REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN, label="single-plane"
    )
    declarations = IdentifiabilityDeclaration.from_report(report, eigenvalue_floor=EIGENVALUE_FLOOR)
    return MolecularParameterPosterior(
        parameters=PARAMETER_NAMES,
        units=dict.fromkeys(PARAMETER_NAMES, "1"),
        location=dict.fromkeys(PARAMETER_NAMES, 0.0),
        budgets={
            name: UncertaintyBudget(
                parameter=name,
                unit="1",
                components=(_component(UncertaintySource.MEASUREMENT, 0.25),),
            )
            for name in PARAMETER_NAMES
        },
        representation=RepresentationDescriptor(
            kind=RepresentationKind.LOCAL_GAUSSIAN,
            uncertainty_sources=(UncertaintySource.MEASUREMENT,),
            axes=PARAMETER_NAMES,
            approximation_error_bound=0.02,
            error_metric="relative-frobenius",
            error_report=ApproximationErrorReport(
                report_sha256="1" * 64,
                source_blob_sha256="2" * 64,
                reduced_blob_sha256="3" * 64,
                metric_id="relative-frobenius",
                requested_tolerance=0.05,
                measured_error=0.02,
            ),
            fallback="explicit ensemble",
        ),
        identifiability={item.parameter: item for item in declarations},
    )


def test_identifiability_is_read_out_of_the_s11_report_not_recomputed() -> None:
    """Composition point with sandbox_inverse: the posterior refuses on ITS measured fractions."""
    report = identifiability_report(
        REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN, label="single-plane"
    )
    declarations = {
        item.parameter: item
        for item in IdentifiabilityDeclaration.from_report(
            report, eigenvalue_floor=EIGENVALUE_FLOOR
        )
    }
    for name in PARAMETER_NAMES:
        assert declarations[name].identifiable_fraction == pytest.approx(
            report.identifiable_fraction(name, eigenvalue_floor=EIGENVALUE_FLOOR), rel=1e-15
        )
    # Measured, not asserted: the single-plane design leaves the PSF width and the axial position
    # sharing one under-measured direction.
    assert declarations["log_psf_sigma_um"].identifiable_fraction == pytest.approx(0.442, abs=1e-3)
    assert declarations["z_um"].identifiable_fraction == pytest.approx(0.558, abs=1e-3)


def test_the_posterior_refuses_a_point_estimate_exactly_where_sandbox_inverse_does() -> None:
    posterior = _posterior(parameters=_image_parameters())
    assert posterior.unidentifiable_parameters() == ("log_psf_sigma_um", "z_um")
    with pytest.raises(UnidentifiableDirectionError, match="not identifiable"):
        posterior.point_estimate("z_um")
    assert posterior.point_estimate("log_length_um") == pytest.approx(0.0)


def test_the_refusal_is_of_the_point_estimate_not_of_the_inference() -> None:
    posterior = _posterior(parameters=_image_parameters())
    low, high = posterior.credible_interval("z_um")
    assert high > low
    assert high == pytest.approx(NormalDist().inv_cdf(0.975) * 0.5, rel=1e-9)


def test_an_unmeasured_parameter_fails_CLOSED() -> None:
    """'Nobody looked' and 'it is identifiable' are different statements."""
    posterior = _posterior()
    assert posterior.parameters.identifiability == {}
    with pytest.raises(UnidentifiableDirectionError, match="no identifiability declaration"):
        posterior.point_estimate("k_off")
    assert posterior.credible_interval("k_off") == (1.4, 1.4)


def test_an_identifiability_measured_at_another_floor_is_refused_not_rescaled() -> None:
    report = identifiability_report(
        REFERENCE_TRUTH, REFERENCE_SCENE, SINGLE_PLANE_DESIGN, label="single-plane"
    )
    wrong_floor = {
        item.parameter: item
        for item in IdentifiabilityDeclaration.from_report(report, eigenvalue_floor=1.0)
    }
    parameters = dataclasses.replace(_image_parameters(), identifiability=wrong_floor)
    with pytest.raises(ValueError, match="meaningless at another"):
        _posterior(parameters=parameters)


# --------------------------------------------------------------------------------------------
# 5. the representation is declared, checked against the content, and fails closed
# --------------------------------------------------------------------------------------------


def test_a_point_posterior_is_legitimate_and_is_not_dressed_up_as_more() -> None:
    posterior = _posterior()
    assert posterior.is_point_representation
    assert posterior.representation_kind is RepresentationKind.EXPLICIT_STATE
    assert posterior.credible_interval("k_off") == (1.4, 1.4)


def test_a_point_representation_may_not_carry_spread() -> None:
    with pytest.raises(ValueError, match="POINT .*carries no spread"):
        _point_parameters(budgets={"k_off": _a1_budget()}, units={"k_off": "1/s"})


def test_a_gaussian_mixture_label_requires_an_actual_mixture() -> None:
    with pytest.raises(ValueError, match="must actually carry a mixture"):
        _point_parameters(
            budgets={"k_off": _convolved_budget()},
            representation=RepresentationDescriptor(
                kind=RepresentationKind.GAUSSIAN_MIXTURE,
                uncertainty_sources=(UncertaintySource.THERMAL,),
                axes=("k_off",),
                approximation_error_bound=0.01,
                error_metric="m",
                error_report=ApproximationErrorReport(
                    report_sha256="4" * 64,
                    source_blob_sha256="5" * 64,
                    reduced_blob_sha256="6" * 64,
                    metric_id="m",
                    requested_tolerance=0.05,
                    measured_error=0.01,
                ),
                fallback="ensemble",
            ),
        )


def test_a_local_gaussian_cannot_represent_a_multi_branch_mixture() -> None:
    with pytest.raises(ValueError, match="cannot represent a multi-branch mixture"):
        _point_parameters(
            budgets={"k_off": _a1_budget()},
            representation=RepresentationDescriptor(
                kind=RepresentationKind.LOCAL_GAUSSIAN,
                uncertainty_sources=(UncertaintySource.THERMAL,),
                axes=("k_off",),
                approximation_error_bound=0.01,
                error_metric="m",
                error_report=ApproximationErrorReport(
                    report_sha256="4" * 64,
                    source_blob_sha256="5" * 64,
                    reduced_blob_sha256="6" * 64,
                    metric_id="m",
                    requested_tolerance=0.05,
                    measured_error=0.01,
                ),
                fallback="ensemble",
            ),
        )


def test_a_surrogate_is_not_a_way_of_representing_a_posterior() -> None:
    assert RepresentationKind.SURROGATE not in POSTERIOR_REPRESENTATION_KINDS
    assert RepresentationKind.SYNTHETIC_IMAGE not in POSTERIOR_REPRESENTATION_KINDS
    assert RepresentationKind.SCALAR_PROJECTION not in POSTERIOR_REPRESENTATION_KINDS


def test_a_reduced_representation_whose_measured_error_exceeds_the_contract_is_refused() -> None:
    """Fail closed: the check is at construction, not at reporting time."""
    tight = dataclasses.replace(CONTRACT, max_representation_error=0.01)
    with pytest.raises(ValueError, match="exceeds the declared"):
        _posterior(parameters=_image_parameters(), contract=tight)


# --------------------------------------------------------------------------------------------
# 6. no self-promotion
# --------------------------------------------------------------------------------------------


def test_the_posterior_cannot_declare_itself_native() -> None:
    """Property over ``__slots__`` with no field behind it, as in surrogate.SurrogatePrediction."""
    posterior = _posterior()
    assert posterior.native_accepted is False
    assert posterior.evidence_authority == EVIDENCE_AUTHORITY
    field_names = {field.name for field in dataclasses.fields(posterior)}
    assert "native_accepted" not in field_names
    assert "evidence_authority" not in field_names
    with pytest.raises(TypeError):
        dataclasses.replace(posterior, native_accepted=True)  # type: ignore[call-arg]


def test_a_subclass_cannot_redefine_the_authority_properties() -> None:
    """External review 2026-07-29, finding 4 — ``dataclasses.replace`` was never the only route.

    ``class ForgedPosterior(CellStatePosterior)`` overriding ``native_accepted`` to return ``True``
    passed ``isinstance``, constructed from the base dataclass fields, and serialised
    ``to_dict()["native_accepted"] == True``: an actual native self-promotion, written into the
    record.  The committed test only attacked ``dataclasses.replace``.
    """
    with pytest.raises(TypeError, match="seals"):

        class ForgedNative(CellStatePosterior):
            @property
            def native_accepted(self) -> bool:
                return True

    with pytest.raises(TypeError, match="seals"):

        class ForgedAuthority(CellStatePosterior):
            @property
            def evidence_authority(self) -> str:
                return "native-accepted"


def test_a_subclass_that_leaves_the_authority_alone_is_still_allowed() -> None:
    """The seal is on two names, not on extension: this is a refusal, not a closed type."""

    class Annotated(CellStatePosterior):
        @property
        def note(self) -> str:
            return "a subclass may add whatever it likes"

    base = _posterior()
    extended = Annotated(
        posterior_id=base.posterior_id,
        manifest=base.manifest,
        graph=base.graph,
        priors=base.priors,
        parameters=base.parameters,
        contract=base.contract,
        evidence=base.evidence,
        envelope=base.envelope,
    )
    assert extended.note
    assert extended.native_accepted is False
    assert extended.to_dict()["native_accepted"] is False
    assert extended.record_hash == base.record_hash


def test_an_authority_property_replaced_after_class_creation_is_caught_too() -> None:
    """Second belt: ``__init_subclass__`` cannot see an override installed later."""

    class Late(CellStatePosterior):
        pass

    base = _posterior()
    clean = Late(
        posterior_id=base.posterior_id,
        manifest=base.manifest,
        graph=base.graph,
        priors=base.priors,
        parameters=base.parameters,
        contract=base.contract,
        evidence=base.evidence,
        envelope=base.envelope,
    )
    Late.native_accepted = property(lambda self: True)  # type: ignore[assignment]
    # already-built instance: caught where it would otherwise be serialised
    assert clean.native_accepted is True  # the forgery is real ...
    with pytest.raises(TypeError, match="not CellStatePosterior's own property"):
        clean.to_dict()
    with pytest.raises(TypeError, match="not CellStatePosterior's own property"):
        _ = clean.record_hash
    # ... and a NEW instance of the tampered class cannot be constructed at all
    with pytest.raises(TypeError, match="not CellStatePosterior's own property"):
        Late(
            posterior_id=base.posterior_id,
            manifest=base.manifest,
            graph=base.graph,
            priors=base.priors,
            parameters=base.parameters,
            contract=base.contract,
            evidence=base.evidence,
            envelope=base.envelope,
        )


def test_the_envelope_independently_refuses_a_native_rung() -> None:
    build = _ring_a_build()
    assert build.manifest is not None
    with pytest.raises(ValueError, match="cannot grant native evidence authority"):
        ArtifactEnvelope(
            artifact_id="a",
            cell_state_manifest_hash=build.manifest.manifest_hash,
            evidence_source=EvidenceSource.NATIVE_ACCEPTED,
            representation=_point_parameters().representation,
        )


def test_every_threshold_is_declared_no_contract_field_has_a_default() -> None:
    for descriptor in dataclasses.fields(PosteriorEvidenceContract):
        assert descriptor.default is dataclasses.MISSING, descriptor.name
        assert descriptor.default_factory is dataclasses.MISSING, descriptor.name


# --------------------------------------------------------------------------------------------
# 7. a regime label needs accepted-step evidence underneath it
# --------------------------------------------------------------------------------------------


def test_a_regime_may_not_be_carried_without_accepted_step_evidence() -> None:
    regime = _limit_cycle_regime()
    assert regime.regime is Regime.LIMIT_CYCLE
    with pytest.raises(ValueError, match="without accepted-step evidence"):
        _posterior(regime=regime)


def test_a_regime_measured_on_an_inadmissible_step_is_not_a_regime() -> None:
    verdict = _admissible_verdict(admissible=False)
    assert verdict.admissible is False
    with pytest.raises(ValueError, match="not admissible|not\\s+transaction-sound"):
        _posterior(regime=_limit_cycle_regime(), accepted_step_evidence=verdict)


def test_a_regime_on_an_admissible_transaction_sound_step_is_carried() -> None:
    posterior = _posterior(
        regime=_limit_cycle_regime(), accepted_step_evidence=_admissible_verdict()
    )
    assert posterior.regime is not None
    assert posterior.regime.is_oscillatory
    assert posterior.regime.evidence_authority == EVIDENCE_AUTHORITY
    assert posterior.to_dict()["accepted_step_admissible"] is True


# --------------------------------------------------------------------------------------------
# 8. the graph half, and what the canonical census does NOT contain
# --------------------------------------------------------------------------------------------


def test_a_connector_instance_the_grammar_does_not_declare_is_refused() -> None:
    grammar = _grammar()
    with pytest.raises(ValueError, match="grammar does not declare"):
        PosteriorGraph(
            grammar=grammar,
            connectors=(
                ConnectorInstance(
                    instance_id="i1",
                    connector_name="invented_tight_junction",
                    endpoint_a_ref="a",
                    endpoint_b_ref="b",
                    owner="nobody",
                ),
            ),
        )


def test_the_graph_reports_zero_cell_cell_junctions_because_the_census_has_none() -> None:
    """Lane E measured this; the posterior reports it rather than letting zero pass as fine."""
    assert PosteriorGraph(grammar=_grammar()).cell_cell_junction_count == 0


def test_the_graph_binds_the_real_engine_census() -> None:
    graph = PosteriorGraph(grammar=_grammar())
    summary = graph.to_dict()
    assert summary["n_components"] == len(_grammar().components)
    assert summary["n_declared_connectors"] == len(_grammar().connectors)
    assert "aleph" in summary["grammar_source"] or summary["grammar_source"]


# --------------------------------------------------------------------------------------------
# 9. composition with the archetype registry — including its refusals
# --------------------------------------------------------------------------------------------


def test_a_ring_a_build_composes_into_a_posterior_carrying_its_pi_gaps() -> None:
    build = _ring_a_build()
    posterior = CellStatePosterior.from_archetype_build(
        build,
        posterior_id="posterior/from-archetype",
        graph=PosteriorGraph(grammar=_grammar()),
        parameters=_point_parameters(),
        contract=CONTRACT,
        evidence=_evidence_graph(),
        artifact_id="artifact/from-archetype",
        evidence_source=EvidenceSource.REDUCED_MODEL,
    )
    assert posterior.manifest is build.manifest
    assert posterior.priors == build.priors
    assert posterior.pi_gap_prior_ids == build.pi_gap_prior_ids
    assert posterior.pi_gap_prior_ids, "the fixture archetype has open PI gaps and must show them"
    assert all(
        prior.status is not PriorStatus.PI_GAP or not prior.parameters for prior in posterior.priors
    )


def test_a_ring_bc_refusal_cannot_be_forced_into_a_posterior() -> None:
    refused = [
        build_archetype(name, grammar=_grammar())
        for name in ARCHETYPE_IDS
        if not build_archetype(name, grammar=_grammar()).supported
    ]
    assert refused, "the registry must contain at least one refused archetype"
    build = refused[0]
    assert build.verdict is not None
    with pytest.raises(ValueError, match="missing laws"):
        CellStatePosterior.from_archetype_build(
            build,
            posterior_id="posterior/forced",
            graph=PosteriorGraph(grammar=_grammar()),
            parameters=_point_parameters(),
            contract=CONTRACT,
            evidence=_evidence_graph(),
            artifact_id="artifact/forced",
            evidence_source=EvidenceSource.REDUCED_MODEL,
        )


# --------------------------------------------------------------------------------------------
# 10. observation operators, model inadequacy, interventions
# --------------------------------------------------------------------------------------------


def _projection_operator() -> SummaryProjectionOperator:
    return SummaryProjectionOperator(
        operator_id="gate/peak-tension",
        spec=ScalarProjectionSpec(
            name="peak-cortex-tension",
            version="v3",
            reduction=ScalarReduction.MAX,
            tolerance=0.05,
            unit="pN/um",
        ),
        array_field="tension",
        axis_names=("time",),
        array_unit="pN/um",
    )


def test_a_posterior_carries_operators_and_observes_through_them() -> None:
    posterior = _posterior(observation_operators=(_projection_operator(),))
    assert posterior.manifest is not None
    state = MechanisticState(
        state_id="state/x",
        manifest_hash=posterior.manifest.manifest_hash,
        fields={"tension": np.array([1.0, 5.0, 3.0])},
        field_units={"tension": "pN/um"},
        source_id="run",
    )
    result = posterior.observe("gate/peak-tension", state)
    assert result.scalar == pytest.approx(5.0)
    assert result.modality is ObservationModality.SUMMARY_PROJECTION


def test_a_carried_unimplemented_operator_refuses_instead_of_returning_a_number() -> None:
    posterior = _posterior(
        observation_operators=(unimplemented_operator(ObservationModality.TRACTION),)
    )
    assert posterior.manifest is not None
    state = MechanisticState(
        state_id="state/x",
        manifest_hash=posterior.manifest.manifest_hash,
        fields={"clutch_forces": ()},
        field_units={"clutch_forces": EMITTER_KIND},
        source_id="run",
    )
    with pytest.raises(ObservationRefused, match="traction"):
        posterior.observe("unimplemented::traction", state)


def test_a_non_operator_cannot_be_carried_as_one() -> None:
    with pytest.raises(TypeError, match="ObservationOperator protocol"):
        _posterior(observation_operators=(object(),))
    with pytest.raises(ValueError, match="repeat an operator_id"):
        _posterior(observation_operators=(_projection_operator(), _projection_operator()))


def test_an_operator_forging_a_blocked_modality_cannot_be_registered() -> None:
    """External review 2026-07-29, finding 3, at the posterior's own boundary.

    A subclass of :class:`SummaryProjectionOperator` declaring ``modality = TRACTION`` and
    ``operator_id = "unimplemented::traction"`` satisfied the runtime-checkable protocol, registered
    here, and ``posterior.observe(...)`` returned the ``reduced-model`` scalar ``5.0`` for a modality
    the capability table calls BLOCKED.  The posterior refuses it at registration now, so there is no
    state for which that number can be produced.
    """

    class ForgedTraction(SummaryProjectionOperator):
        @property
        def modality(self) -> ObservationModality:
            return ObservationModality.TRACTION

    forged = ForgedTraction(
        operator_id="unimplemented::traction",
        spec=_projection_operator().spec,
        array_field="tension",
        axis_names=("time",),
        array_unit="pN/um",
    )
    assert isinstance(forged, ObservationOperator)
    with pytest.raises(OperatorRegistrationRefused, match="blocked"):
        _posterior(observation_operators=(forged,))
    # and the declared refusal for the same modality is still carried without complaint
    posterior = _posterior(
        observation_operators=(unimplemented_operator(ObservationModality.TRACTION),)
    )
    assert posterior.operator("unimplemented::traction").modality is ObservationModality.TRACTION


def test_model_inadequacy_arrives_as_a_proposal_and_stays_one() -> None:
    proposal = _missing_law_proposal()
    posterior = _posterior(model_inadequacy=(proposal,))
    assert posterior.missing_law_ids == proposal.law_ids
    assert all(not item.decided for item in posterior.model_inadequacy)
    assert proposal.evidence_source is EvidenceSource.ANALYTIC_ORACLE


def test_an_intervention_response_must_name_the_card_that_preregistered_it() -> None:
    card = SandboxExperimentCard(
        experiment_id="s-demo",
        question="does blebbistatin reduce cortical tension?",
        manifest_hashes=("h",),
        intervention="myosin-inhibition",
        observables=("cortical_tension",),
        positive_controls=("no-inhibition",),
        negative_controls=("vehicle",),
        adversarial_controls=("renderer-holdout",),
        held_out_interventions=("half-dose",),
        falsifier="tension is unchanged at full inhibition",
        failure_expansions=("missing-law: cortex turnover",),
        budget_class="slice",
    )
    response = InterventionResponse(
        intervention_id="blebbistatin-10uM",
        observable="cortical_tension",
        unit="pN/um",
        location=-40.0,
        budget=UncertaintyBudget(
            parameter="cortical_tension",
            unit="pN/um",
            components=(_component(UncertaintySource.BIOLOGICAL, 25.0),),
        ),
        card_hash=card.card_hash,
    )
    posterior = _posterior(interventions=(response,))
    assert posterior.interventions[0].card_hash == card.card_hash
    low, high = response.credible_interval(level=0.95)
    assert low == pytest.approx(-40.0 - NormalDist().inv_cdf(0.975) * 5.0, rel=1e-9)
    assert high == pytest.approx(-40.0 + NormalDist().inv_cdf(0.975) * 5.0, rel=1e-9)
    with pytest.raises(ValueError, match="SHA-256"):
        dataclasses.replace(response, card_hash="not-a-hash")


# --------------------------------------------------------------------------------------------
# 11. provenance stays bound, and the record is content-addressed
# --------------------------------------------------------------------------------------------


def test_an_envelope_bound_to_another_manifest_is_refused() -> None:
    other = build_archetype(ARCHETYPE_IDS[1], grammar=_grammar())
    assert other.manifest is not None
    parameters = _point_parameters()
    with pytest.raises(ValueError, match="different manifest"):
        _posterior(
            envelope=ArtifactEnvelope(
                artifact_id="a",
                cell_state_manifest_hash=other.manifest.manifest_hash,
                evidence_source=EvidenceSource.REDUCED_MODEL,
                representation=parameters.representation,
            ),
            parameters=parameters,
        )


def test_a_synthetic_observation_source_must_stay_excluded() -> None:
    posterior = _posterior()
    posterior.assert_observation_excluded("posterior_draft")
    with pytest.raises(ValueError, match="citable as biological evidence"):
        posterior.assert_observation_excluded("lit_koff")


def _other_evidence_graph() -> EvidenceGraph:
    """A DIFFERENT provenance DAG carrying the same source ids as ``_evidence_graph``.

    Same ids, different content: this is the 10-lane composition case the review named, where a
    posterior re-uses a source id for a graph that says something else.
    """
    registry = SourceRegistry()
    registry.register(
        SourceRecord(
            source_id="lit_koff",
            provenance=ProvenanceKind.SOURCED,
            origin=SourceOrigin.LITERATURE,
            status=SourceStatus.ACTIVE,
            description="a DIFFERENT published off-rate measurement",
            citation="doi:10.0000/other",
        )
    )
    registry.register(
        SourceRecord(
            source_id="posterior_draft",
            provenance=ProvenanceKind.INFERRED_POSTERIOR,
            origin=SourceOrigin.COMPUTATION,
            status=SourceStatus.ACTIVE,
            description="posterior inferred from the other prior",
            parents=("lit_koff",),
        )
    )
    return registry.freeze()


def _intervention_response(location: float = -40.0) -> InterventionResponse:
    return InterventionResponse(
        intervention_id="blebbistatin-10uM",
        observable="cortical_tension",
        unit="pN/um",
        location=location,
        budget=UncertaintyBudget(
            parameter="cortical_tension",
            unit="pN/um",
            components=(_component(UncertaintySource.BIOLOGICAL, 25.0),),
        ),
        card_hash="7" * 64,
    )


def _rich_posterior_kwargs() -> dict[str, object]:
    """A posterior with every optional field populated, so each one can be varied in turn."""
    return {
        "accepted_step_evidence": _admissible_verdict(),
        "regime": _limit_cycle_regime(),
        "interventions": (_intervention_response(),),
        "observation_operators": (_projection_operator(),),
        "model_inadequacy": (_missing_law_proposal(),),
    }


def _envelope_for(manifest: object, parameters: MolecularParameterPosterior, **overrides: object):
    kwargs: dict[str, object] = {
        "artifact_id": "artifact/posterior-fixture",
        "cell_state_manifest_hash": manifest.manifest_hash,  # type: ignore[attr-defined]
        "evidence_source": EvidenceSource.REDUCED_MODEL,
        "representation": parameters.representation,
        "source_artifact_ids": ("artifact/accepted-run",),
    }
    kwargs.update(overrides)
    return ArtifactEnvelope(**kwargs)  # type: ignore[arg-type]


def _vary_manifest() -> dict[str, object]:
    other = build_archetype(ARCHETYPE_IDS[1], grammar=_grammar())
    assert other.manifest is not None
    parameters = _point_parameters()
    return {
        "manifest": other.manifest,
        "parameters": parameters,
        "envelope": _envelope_for(other.manifest, parameters),
    }


def _vary_envelope() -> dict[str, object]:
    build = _ring_a_build()
    assert build.manifest is not None
    parameters = _point_parameters()
    return {
        "parameters": parameters,
        "envelope": _envelope_for(
            build.manifest, parameters, source_artifact_ids=("artifact/some-other-run",)
        ),
    }


def _vary_graph() -> dict[str, object]:
    grammar = _grammar()
    connector_name = sorted(grammar.connector_names)[0]
    return {
        "graph": PosteriorGraph(
            grammar=grammar,
            connectors=(
                ConnectorInstance(
                    instance_id="i1",
                    connector_name=connector_name,
                    endpoint_a_ref="a",
                    endpoint_b_ref="b",
                    owner="fixture",
                ),
            ),
        )
    }


#: One entry per hashed component of the record.  The review's counterexample (finding 1) was that
#: this list used to be one entry long: the committed test varied ``posterior_id`` and nothing else,
#: while ``evidence`` was absent from the record entirely and the accepted-step verdict was one bit.
RECORD_COMPONENT_VARIATIONS: tuple[tuple[str, object], ...] = (
    ("posterior_id", lambda: {"posterior_id": "posterior/other"}),
    ("manifest", _vary_manifest),
    ("graph", _vary_graph),
    ("priors", lambda: {"priors": _ring_a_build().priors[1:]}),
    ("parameters", lambda: {"parameters": _point_parameters(location={"k_off": 2.8})}),
    ("contract", lambda: {"contract": dataclasses.replace(CONTRACT, credible_level=0.9)}),
    ("evidence", lambda: {"evidence": _other_evidence_graph()}),
    ("envelope", _vary_envelope),
    (
        "accepted_step_evidence",
        lambda: {"accepted_step_evidence": _admissible_verdict(step_index=999)},
    ),
    ("regime", lambda: {"regime": _limit_cycle_regime(frequency_hz=3.5)}),
    ("interventions", lambda: {"interventions": (_intervention_response(location=-25.0),)}),
    (
        "observation_operators",
        lambda: {
            "observation_operators": (
                _projection_operator(),
                unimplemented_operator(ObservationModality.SEQUENCING),
            )
        },
    ),
    ("model_inadequacy", lambda: {"model_inadequacy": ()}),
)


@pytest.mark.parametrize(
    ("component", "variation"),
    RECORD_COMPONENT_VARIATIONS,
    ids=[name for name, _ in RECORD_COMPONENT_VARIATIONS],
)
def test_the_record_hash_changes_when_any_part_of_the_posterior_does(
    component: str, variation: object
) -> None:
    """Every constructor field, one at a time — the shape the review said this test needed.

    A record_hash that misses a field is not content-addressing; it is content-addressing of the
    fields somebody remembered.
    """
    assert component in {field.name for field in dataclasses.fields(CellStatePosterior)}
    baseline_kwargs = _rich_posterior_kwargs()
    baseline = _posterior(**baseline_kwargs)
    varied = _posterior(**{**baseline_kwargs, **variation()})  # type: ignore[operator]
    assert varied.record_hash != baseline.record_hash, f"{component} is invisible to record_hash"
    assert _posterior(**baseline_kwargs).record_hash == baseline.record_hash


def test_the_record_hash_counterexample_from_the_external_review_is_refused() -> None:
    """Verbatim reproduction of finding 1 (2026-07-29).

    Two posteriors with different :class:`EvidenceGraph`\\ s and ``attempted_step_index`` 7 against
    999 hashed IDENTICALLY: ``to_dict()`` carried no ``evidence`` at all and collapsed the
    accepted-step transaction to its ``admissible`` bit.
    """
    first = _posterior(evidence=_evidence_graph(), accepted_step_evidence=_admissible_verdict())
    second = _posterior(
        evidence=_other_evidence_graph(),
        accepted_step_evidence=_admissible_verdict(step_index=999),
    )
    assert first.evidence.graph_hash != second.evidence.graph_hash
    assert first.accepted_step_evidence is not None
    assert second.accepted_step_evidence is not None
    assert first.accepted_step_evidence.attempted_step_index == 7
    assert second.accepted_step_evidence.attempted_step_index == 999
    assert first.record_hash != second.record_hash

    record = first.to_dict()
    assert record["evidence_graph_sha256"] == first.evidence.graph_hash
    assert record["accepted_step_evidence"]["attempted_step_index"] == 7
    assert record["accepted_step_evidence"]["transaction_sound"] is True
    # the collapsed bit stays too — it is a legitimate summary, it was just never an identity
    assert record["accepted_step_admissible"] is True
    assert record["native_accepted"] is False


def test_every_constructor_field_reaches_the_record() -> None:
    """The structural form of the same rule: a field the record omits is a hash blind spot."""
    covered = {name for name, _ in RECORD_COMPONENT_VARIATIONS}
    assert covered == {field.name for field in dataclasses.fields(CellStatePosterior)}


def test_the_assembled_posterior_carries_every_field_the_plan_declares() -> None:
    """The A13 verdict is about a shape, so the shape is asserted rather than described."""
    posterior = _posterior(
        parameters=_image_parameters(),
        accepted_step_evidence=_admissible_verdict(),
        regime=_limit_cycle_regime(),
        observation_operators=(
            _projection_operator(),
            unimplemented_operator(ObservationModality.SEQUENCING),
        ),
        model_inadequacy=(_missing_law_proposal(),),
    )
    record = posterior.to_dict()
    assert record["graph"]["n_declared_connectors"] > 0  # component/connector graph
    assert record["priors"]  # population and geometry distribution
    assert record["parameters"]["budgets"]  # molecular-parameter posterior
    assert record["regime"]["regime"] == "limit-cycle"  # accepted dynamical regime
    assert record["observation_operators"]  # observation operators
    assert record["evidence_source"] == "reduced-model"  # evidence provenance
    assert record["missing_law_ids"]  # known model inadequacy
    assert record["evidence_authority"] == EVIDENCE_AUTHORITY

"""CPU oracle tests for the observation-operator interface.

Every quantity checked here is either a closed form (the emitted photon budget), an identity against
the code the operator adapts (bit-exact equality with the renderer it wraps), or a structural refusal.
No cell physics, no simulation trace, no GPU: what is validated is the INTERFACE, not a cell model.

The identity tests are the load-bearing ones.  An adapter that silently rescales, re-normalises or
re-orders what it wraps would break the Sanity Gate of the module it adapts without breaking anything
here, so the operator output is required to be bit-identical to calling the renderer directly.
"""

from __future__ import annotations

import dataclasses
import hashlib
import math
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from aleph.virtual_cell.contracts import EvidenceSource, RepresentationKind
from aleph.virtual_cell.observation_operator import (
    BLOCKED_MODALITIES,
    EMITTER_KIND,
    PHOTON_UNIT,
    RESERVED_OPERATOR_ID_PREFIX,
    UNIMPLEMENTED_OPERATORS,
    ApplicabilityVerdict,
    CapabilityStatus,
    FluorescenceOperator,
    InputRequirement,
    MechanisticState,
    ObservationCapabilityVerdict,
    ObservationModality,
    ObservationOperator,
    ObservationRefused,
    ObservationResult,
    OperatorRegistrationRefused,
    SummaryProjectionOperator,
    UnimplementedOperator,
    _content_token,
    unimplemented_operator,
    validate_operator_registration,
)
from aleph.virtual_cell.optics import (
    AxialModel,
    LineEmitter,
    OpticsCapabilityVerdict,
    OpticsConfig,
    PointEmitter,
    PSFSampling,
    TruncationPolicy,
    render_emitters,
    render_widefield_stack,
)
from aleph.virtual_cell.sidecar import ArrayAxis, ArraySidecarDescriptor
from aleph.virtual_cell.source_registry import (
    ExclusionReason,
    ProvenanceKind,
    SourceOrigin,
    SourceRecord,
    SourceRegistry,
    SourceStatus,
)
from aleph.virtual_cell.synthetic_microscopy import MicroscopyConfig
from aleph.virtual_cell.telemetry import ScalarProjectionSpec, ScalarReduction

# --------------------------------------------------------------------------------------------
# fixtures — declared once, before any measurement
# --------------------------------------------------------------------------------------------

MANIFEST_HASH = "a" * 64
OTHER_MANIFEST_HASH = "b" * 64

MICROSCOPY = MicroscopyConfig(
    height_px=32,
    width_px=32,
    pixel_size_um=0.1,
    psf_sigma_um=0.2,
    background=0.0,
    psf_radius_sigma=6.0,
)
OPTICS = OpticsConfig(
    microscopy=MICROSCOPY,
    sampling=PSFSampling.PIXEL_INTEGRATED,
    truncation=TruncationPolicy.ANALYTIC,
    axial_model=AxialModel.WIDEFIELD,
    rayleigh_range_um=0.5,
    focal_plane_z_um=0.0,
    exposure_s=2.0,
)

POINT_RATE_PER_S = 1000.0
EMITTERS = (
    PointEmitter(xy_um=(1.6, 1.6), photon_rate=POINT_RATE_PER_S, z_um=0.0),
    LineEmitter(
        start_xy_um=(0.8, 0.8), end_xy_um=(2.4, 1.2), photon_rate_per_length=200.0, z_um=0.1
    ),
)


def _fluorescence(**overrides: object) -> FluorescenceOperator:
    kwargs: dict[str, object] = {
        "operator_id": "fluo/actin",
        "channel": "actin",
        "optics": OPTICS,
        "emitter_field": "actin_emitters",
    }
    kwargs.update(overrides)
    return FluorescenceOperator(**kwargs)  # type: ignore[arg-type]


def _emitter_state(**overrides: object) -> MechanisticState:
    kwargs: dict[str, object] = {
        "state_id": "state/accepted-0007",
        "manifest_hash": MANIFEST_HASH,
        "fields": {"actin_emitters": EMITTERS},
        "field_units": {"actin_emitters": EMITTER_KIND},
        "source_id": "sim_run_0007",
    }
    kwargs.update(overrides)
    return MechanisticState(**kwargs)  # type: ignore[arg-type]


PROJECTION = ScalarProjectionSpec(
    name="peak-cortex-tension",
    version="v3",
    reduction=ScalarReduction.MAX,
    tolerance=0.05,
    selection={"component": 1},
    unit="pN/um",
)
ARRAY_AXES = ("time", "component")
ARRAY = np.arange(12.0, dtype=np.float64).reshape(6, 2)


def _projection_operator(**overrides: object) -> SummaryProjectionOperator:
    kwargs: dict[str, object] = {
        "operator_id": "gate/peak-tension",
        "spec": PROJECTION,
        "array_field": "tension_by_component",
        "axis_names": ARRAY_AXES,
        "array_unit": "pN/um",
    }
    kwargs.update(overrides)
    return SummaryProjectionOperator(**kwargs)  # type: ignore[arg-type]


def _array_state(**overrides: object) -> MechanisticState:
    kwargs: dict[str, object] = {
        "state_id": "state/telemetry-0007",
        "manifest_hash": MANIFEST_HASH,
        "fields": {"tension_by_component": ARRAY},
        "field_units": {"tension_by_component": "pN/um"},
        "source_id": "sim_run_0007",
    }
    kwargs.update(overrides)
    return MechanisticState(**kwargs)  # type: ignore[arg-type]


# --------------------------------------------------------------------------------------------
# 0. the module stays on the CPU
# --------------------------------------------------------------------------------------------


def test_importing_the_module_pulls_in_no_gpu_stack() -> None:
    """Measured in a FRESH interpreter.

    An in-process ``"warp" not in sys.modules`` assertion measures the whole pytest session, not this
    import, so it passes alone and fails under the full suite.  A subprocess measures what the name
    says.
    """
    script = textwrap.dedent(
        """
        import sys
        from aleph.virtual_cell.observation_operator import (
            FluorescenceOperator,
            SummaryProjectionOperator,
            unimplemented_operator,
        )

        for blocked in ("warp", "torch", "jax"):
            assert blocked not in sys.modules, f"observation_operator import pulled in {blocked}"
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
# 1. the protocol is structural, and every adapter satisfies it
# --------------------------------------------------------------------------------------------


def test_every_adapter_satisfies_the_protocol_without_inheriting_from_it() -> None:
    """Structural typing is the point: a future runtime's renderer need not subclass this package."""
    operators = (
        _fluorescence(),
        _projection_operator(),
        unimplemented_operator(ObservationModality.TRACTION),
    )
    for operator in operators:
        assert isinstance(operator, ObservationOperator)
        assert not isinstance(operator, type(ObservationOperator))
    assert {operator.modality for operator in operators} == {
        ObservationModality.FLUORESCENCE,
        ObservationModality.SUMMARY_PROJECTION,
        ObservationModality.TRACTION,
    }


def test_an_object_that_merely_looks_numeric_is_not_an_operator() -> None:
    assert not isinstance(object(), ObservationOperator)
    assert not isinstance(np.zeros(3), ObservationOperator)


def test_modality_is_a_property_so_an_operator_cannot_relabel_itself() -> None:
    """Same pattern as ``surrogate.SurrogatePrediction.evidence_source``: no field to overwrite."""
    operator = _fluorescence()
    assert "modality" not in {field.name for field in dataclasses.fields(operator)}
    with pytest.raises(TypeError):
        dataclasses.replace(operator, modality=ObservationModality.TRACTION)  # type: ignore[call-arg]


def test_the_map_is_one_way_the_protocol_has_no_inverse() -> None:
    """Inversion is a different, harder problem; it lives in sandbox_inverse where it is measured."""
    names = set(dir(ObservationOperator))
    assert not names & {"invert", "inverse", "decode", "estimate_state"}


# --------------------------------------------------------------------------------------------
# 2. fluorescence — identity with the renderer it adapts
# --------------------------------------------------------------------------------------------


def test_the_operator_output_is_bit_identical_to_calling_the_renderer_directly() -> None:
    """The adapter must add nothing.  Optics' photon-budget gate then still applies here."""
    result = _fluorescence().observe(_emitter_state())
    direct = render_emitters(list(EMITTERS), OPTICS)
    assert np.array_equal(result.values, direct)
    assert result.units == PHOTON_UNIT


def test_the_stack_operator_is_bit_identical_to_the_widefield_stack() -> None:
    planes = (-0.3, 0.0, 0.3)
    result = _fluorescence(focal_planes_um=planes).observe(_emitter_state())
    direct = render_widefield_stack(list(EMITTERS), OPTICS, planes)
    assert np.array_equal(result.values, direct)
    assert result.values.shape == (3, MICROSCOPY.height_px, MICROSCOPY.width_px)


def test_the_emitted_photon_budget_survives_the_operator() -> None:
    """Dimensional gate: rate [1/s] x exposure [s] = photons, recovered up to the truncation bound."""
    single = (PointEmitter(xy_um=(1.6, 1.6), photon_rate=POINT_RATE_PER_S, z_um=0.0),)
    state = _emitter_state(
        fields={"actin_emitters": single}, field_units={"actin_emitters": EMITTER_KIND}
    )
    total = float(np.sum(_fluorescence().observe(state).values))
    emitted = POINT_RATE_PER_S * OPTICS.exposure_s
    # Optics' own bound: 2 erfc(R / sqrt 2) for psf_radius_sigma = R.
    bound = 2.0 * math.erfc(MICROSCOPY.psf_radius_sigma / math.sqrt(2.0))
    assert total <= emitted
    assert emitted - total <= emitted * bound


def test_an_empty_emitter_sequence_renders_the_background_and_does_not_raise() -> None:
    state = _emitter_state(
        fields={"actin_emitters": ()}, field_units={"actin_emitters": EMITTER_KIND}
    )
    values = _fluorescence().observe(state).values
    assert values.shape == (MICROSCOPY.height_px, MICROSCOPY.width_px)
    assert float(np.max(values)) == pytest.approx(MICROSCOPY.background)


def test_the_result_is_read_only_so_a_caller_cannot_edit_an_observation() -> None:
    result = _fluorescence().observe(_emitter_state())
    with pytest.raises(ValueError):
        result.values[0, 0] = 1.0


# --------------------------------------------------------------------------------------------
# 3. refusal — an operator says no rather than returning a number
# --------------------------------------------------------------------------------------------


def test_a_state_without_the_required_field_is_refused_by_name() -> None:
    state = _emitter_state(fields={"membrane": ()}, field_units={"membrane": EMITTER_KIND})
    verdict = _fluorescence().applicability(state)
    assert verdict.applicable is False
    assert verdict.missing_fields == ("actin_emitters",)
    assert verdict.routes_to is not None
    with pytest.raises(ObservationRefused) as excinfo:
        _fluorescence().observe(state)
    assert excinfo.value.verdict.missing_fields == ("actin_emitters",)


def test_a_unit_mismatch_is_refused_rather_than_reinterpreted() -> None:
    state = _emitter_state(field_units={"actin_emitters": "photons/s"})
    verdict = _fluorescence().applicability(state)
    assert verdict.applicable is False
    assert len(verdict.unit_mismatches) == 1
    assert "photons/s" in verdict.unit_mismatches[0]


def test_a_sequence_holding_a_non_emitter_is_refused_with_the_offending_index() -> None:
    state = _emitter_state(
        fields={"actin_emitters": (EMITTERS[0], 3.14)},
        field_units={"actin_emitters": EMITTER_KIND},
    )
    verdict = _fluorescence().applicability(state)
    assert verdict.applicable is False
    assert verdict.unit_mismatches == ("actin_emitters[1]: float",)


def test_a_single_emitter_is_not_a_labelled_structure() -> None:
    state = _emitter_state(
        fields={"actin_emitters": EMITTERS[0]}, field_units={"actin_emitters": EMITTER_KIND}
    )
    assert _fluorescence().applicability(state).applicable is False


def test_a_refusal_must_name_the_expansion_it_opens() -> None:
    """A bare refusal closes scope; the failure routing this project uses is expansionary."""
    with pytest.raises(ValueError, match="closes scope"):
        ApplicabilityVerdict(
            operator_id="op", state_id="s", applicable=False, reason="no", routes_to=None
        )


def test_an_applicable_verdict_cannot_also_carry_a_defect() -> None:
    with pytest.raises(ValueError, match="cannot also name a defect"):
        ApplicabilityVerdict(
            operator_id="op",
            state_id="s",
            applicable=True,
            reason="fine",
            missing_fields=("x",),
        )


# --------------------------------------------------------------------------------------------
# 4. the legacy scalar-gate path stays versioned and content-hashed
# --------------------------------------------------------------------------------------------


def test_the_projection_operator_reproduces_the_stored_projection_exactly() -> None:
    result = _projection_operator().observe(_array_state())
    assert result.scalar == PROJECTION.evaluate(ARRAY, ARRAY_AXES)
    assert result.scalar == pytest.approx(11.0)
    assert result.units == "pN/um"
    assert result.envelope.projection_hash == PROJECTION.definition_sha256
    assert result.envelope.representation.kind is RepresentationKind.SCALAR_PROJECTION


def test_relaxing_a_projection_tolerance_mints_a_different_observation_not_the_same_one() -> None:
    """The tolerance is part of the projection hash, so a relaxed gate cannot wear the old identity."""
    relaxed = dataclasses.replace(PROJECTION, tolerance=0.5)
    baseline = _projection_operator().observe(_array_state())
    loosened = _projection_operator(spec=relaxed).observe(_array_state())
    assert loosened.scalar == baseline.scalar
    assert loosened.envelope.projection_hash != baseline.envelope.projection_hash
    assert loosened.envelope.observation_provenance_sha256 != (
        baseline.envelope.observation_provenance_sha256
    )


def test_an_array_of_the_wrong_rank_is_refused_rather_than_reduced_on_the_wrong_axis() -> None:
    state = _array_state(
        fields={"tension_by_component": np.arange(6.0)},
        field_units={"tension_by_component": "pN/um"},
    )
    verdict = _projection_operator().applicability(state)
    assert verdict.applicable is False
    assert "rank" in verdict.reason
    with pytest.raises(ObservationRefused):
        _projection_operator().observe(state)


def test_a_projection_selecting_an_undeclared_axis_cannot_be_built() -> None:
    with pytest.raises(ValueError, match="does not declare"):
        _projection_operator(axis_names=("time", "seed"))


def test_reducing_an_image_to_a_scalar_here_would_be_an_undeclared_projection() -> None:
    result = _fluorescence().observe(_emitter_state())
    with pytest.raises(ValueError, match="undeclared projection"):
        _ = result.scalar


# --------------------------------------------------------------------------------------------
# 5. what is NOT implemented is named, not stubbed
# --------------------------------------------------------------------------------------------


def test_the_four_unimplemented_modalities_are_named_with_a_reason_and_a_route() -> None:
    assert set(UNIMPLEMENTED_OPERATORS) == {
        ObservationModality.BRIGHTFIELD,
        ObservationModality.TRACTION,
        ObservationModality.FORCE_INDENTATION,
        ObservationModality.SEQUENCING,
    }
    for modality, verdict in UNIMPLEMENTED_OPERATORS.items():
        assert verdict.modality is modality
        assert verdict.status is CapabilityStatus.BLOCKED
        assert len(verdict.why_not_approximated) > 40
        assert "extension" in verdict.routes_to
        operator = unimplemented_operator(modality)
        # The gap is in the OPERATOR, not in the state: it requires nothing and still refuses.
        assert operator.input_requirement.fields == ()
        assert operator.output_units == "none"
        assert operator.modality is modality


def test_the_capability_vocabulary_matches_the_one_optics_already_validates_against() -> None:
    """Kept identical on purpose: a capability gap must read the same wherever it is raised."""
    assert {status.value for status in CapabilityStatus} == {
        "blocked",
        "oracle-only",
        "pi-decision",
    }
    optics_verdict = OpticsCapabilityVerdict(
        capability="c", status="oracle-only", why_not_approximated="w", routes_to="r"
    )
    assert optics_verdict.status == CapabilityStatus.ORACLE_ONLY.value
    with pytest.raises(ValueError):
        OpticsCapabilityVerdict(
            capability="c", status="not-a-status", why_not_approximated="w", routes_to="r"
        )


def test_an_unimplemented_operator_refuses_every_state_and_returns_no_number() -> None:
    operator = unimplemented_operator(ObservationModality.FORCE_INDENTATION)
    assert isinstance(operator, UnimplementedOperator)
    assert operator.applicability(_emitter_state()).applicable is False
    with pytest.raises(ObservationRefused, match="force-indentation|AFM"):
        operator.observe(_emitter_state())
    with pytest.raises(ObservationRefused, match="no artifact"):
        operator.provenance_record(
            source_id="afm_curve", parent_source_ids=(), description="an AFM curve"
        )


def test_asking_for_the_unimplemented_form_of_an_implemented_modality_is_a_caller_error() -> None:
    with pytest.raises(KeyError, match="implemented"):
        unimplemented_operator(ObservationModality.FLUORESCENCE)
    with pytest.raises(KeyError):
        unimplemented_operator(ObservationModality.SUMMARY_PROJECTION)


def test_the_afm_gap_says_the_closed_form_is_the_oracle_and_not_the_operator() -> None:
    """The architecture principle inverted exactly this; the gap record says so explicitly."""
    verdict = UNIMPLEMENTED_OPERATORS[ObservationModality.FORCE_INDENTATION]
    assert "oracle" in verdict.why_not_approximated.lower()
    assert "Hertz" in verdict.why_not_approximated or "Bell-Evans" in verdict.why_not_approximated


def test_a_confocal_fluorescence_operator_cannot_be_constructed_at_all() -> None:
    """Optics already refuses the config, so the refusal is inherited rather than re-implemented."""
    with pytest.raises(NotImplementedError, match="pinhole"):
        OpticsConfig(microscopy=MICROSCOPY, axial_model=AxialModel.CONFOCAL)


class _ForgedTraction(SummaryProjectionOperator):
    """The external review's counterexample (2026-07-29, finding 3), kept as a test fixture.

    It satisfies the runtime-checkable protocol, and before the registration boundary existed it was
    carried in a posterior's operator table and returned ``5.0`` for a BLOCKED modality.
    """

    @property
    def modality(self) -> ObservationModality:
        return ObservationModality.TRACTION


def _forged(**overrides: object) -> _ForgedTraction:
    kwargs: dict[str, object] = {
        "operator_id": f"{RESERVED_OPERATOR_ID_PREFIX}traction",
        "spec": PROJECTION,
        "array_field": "tension_by_component",
        "axis_names": ARRAY_AXES,
        "array_unit": "pN/um",
    }
    kwargs.update(overrides)
    return _ForgedTraction(**kwargs)  # type: ignore[arg-type]


def test_the_forged_blocked_modality_operator_still_satisfies_the_protocol() -> None:
    """The premise of the finding: structural typing cannot see the lie, so registration must."""
    forged = _forged()
    assert isinstance(forged, ObservationOperator)
    assert forged.modality is ObservationModality.TRACTION
    assert forged.observe(_array_state()).scalar == pytest.approx(11.0)


def test_registration_refuses_an_operator_claiming_a_blocked_modality() -> None:
    """The fix is registry-GLOBAL: `UnimplementedOperator.observe` raising could not reach this."""
    with pytest.raises(OperatorRegistrationRefused, match="blocked"):
        validate_operator_registration(_forged())
    with pytest.raises(OperatorRegistrationRefused, match="blocked"):
        validate_operator_registration(_forged(operator_id="traction/my-own-estimate"))


def test_registration_refuses_a_subclass_of_the_declared_refusal_too() -> None:
    """A subclass could override ``observe``; only the declared type may carry a blocked modality."""

    class Sneaky(UnimplementedOperator):
        def observe(self, state: MechanisticState) -> ObservationResult:  # pragma: no cover
            raise AssertionError("must never be reachable through a registered table")

    sneaky = Sneaky(
        operator_id=f"{RESERVED_OPERATOR_ID_PREFIX}traction",
        verdict=UNIMPLEMENTED_OPERATORS[ObservationModality.TRACTION],
    )
    assert isinstance(sneaky, ObservationOperator)
    with pytest.raises(OperatorRegistrationRefused, match="blocked"):
        validate_operator_registration(sneaky)


def test_registration_refuses_a_blocked_modality_carrying_a_rewritten_capability_verdict() -> None:
    """The capability table is the authority on what a blocked modality's gap says."""
    softened = dataclasses.replace(
        UNIMPLEMENTED_OPERATORS[ObservationModality.TRACTION],
        status=CapabilityStatus.ORACLE_ONLY,
    )
    with pytest.raises(OperatorRegistrationRefused, match="blocked"):
        validate_operator_registration(
            UnimplementedOperator(
                operator_id=f"{RESERVED_OPERATOR_ID_PREFIX}traction", verdict=softened
            )
        )


def test_the_declared_refusal_must_use_its_own_reserved_id() -> None:
    with pytest.raises(OperatorRegistrationRefused, match="registered as"):
        validate_operator_registration(
            UnimplementedOperator(
                operator_id="traction/gap",
                verdict=UNIMPLEMENTED_OPERATORS[ObservationModality.TRACTION],
            )
        )


def test_the_reserved_id_namespace_cannot_be_worn_by_an_implemented_operator() -> None:
    with pytest.raises(OperatorRegistrationRefused, match="reserved"):
        validate_operator_registration(
            _projection_operator(operator_id=f"{RESERVED_OPERATOR_ID_PREFIX}summary")
        )


def test_registration_admits_every_genuine_operator_unchanged() -> None:
    for operator in (
        _fluorescence(),
        _projection_operator(),
        *(unimplemented_operator(modality) for modality in sorted(BLOCKED_MODALITIES)),
    ):
        assert validate_operator_registration(operator) is operator


def test_registration_still_refuses_a_non_operator() -> None:
    assert issubclass(OperatorRegistrationRefused, TypeError)
    with pytest.raises(OperatorRegistrationRefused, match="ObservationOperator protocol"):
        validate_operator_registration(object())


def test_the_blocked_modality_set_is_derived_from_the_capability_table() -> None:
    """One place decides what is blocked; a second copy is how the two disagree."""
    assert frozenset(UNIMPLEMENTED_OPERATORS) == BLOCKED_MODALITIES
    assert ObservationModality.FLUORESCENCE not in BLOCKED_MODALITIES
    assert ObservationModality.SUMMARY_PROJECTION not in BLOCKED_MODALITIES


def test_a_capability_verdict_requires_a_reason_the_approximation_is_wrong() -> None:
    with pytest.raises(ValueError, match="why_not_approximated"):
        ObservationCapabilityVerdict(
            modality=ObservationModality.TRACTION,
            capability="traction",
            status=CapabilityStatus.BLOCKED,
            why_not_approximated="   ",
            routes_to="somewhere",
        )


# --------------------------------------------------------------------------------------------
# 6. provenance — a synthetic observation stays out of the biological corpus, transitively
# --------------------------------------------------------------------------------------------


def test_a_rendered_observation_carries_synthetic_provenance_and_cannot_claim_native() -> None:
    result = _fluorescence().observe(_emitter_state())
    assert result.envelope.evidence_source is EvidenceSource.SYNTHETIC
    assert result.envelope.representation.kind is RepresentationKind.SYNTHETIC_IMAGE
    assert result.envelope.observation_provenance_sha256 is not None
    assert result.envelope.source_artifact_ids == ("state/accepted-0007",)


def test_the_two_hop_case_a_reanalysis_of_a_synthetic_image_is_still_excluded() -> None:
    """Hop 0 is citable literature; hop 1 renders from it; hop 2 re-analyses the render.

    Hop 2 is deliberately given the strongest possible standing on its own — an empirical LITERATURE
    origin with a citation, which is exactly what would otherwise make it citable — so the exclusion
    it inherits cannot be attributed to its own origin.  That is the laundering route this refuses:
    *published, therefore biological evidence*, for a number whose ancestry is a synthetic render.
    Hop 1 is excluded for its own non-empirical origin and hop 2 inherits it, naming hop 1 as the
    blocking source (master plan section 10.2, A6).
    """
    registry = SourceRegistry()
    registry.register(
        SourceRecord(
            source_id="paper_geometry",
            provenance=ProvenanceKind.SOURCED,
            origin=SourceOrigin.LITERATURE,
            status=SourceStatus.ACTIVE,
            description="published cortical geometry measurement",
            citation="doi:10.0000/fixture",
        )
    )
    render_node = _fluorescence().provenance_record(
        source_id="synthetic_frame",
        parent_source_ids=("paper_geometry",),
        description="widefield render of the accepted state",
    )
    registry.register(render_node)
    registry.register(
        SourceRecord(
            source_id="published_reanalysis",
            provenance=ProvenanceKind.DERIVED,
            origin=SourceOrigin.LITERATURE,
            status=SourceStatus.ACTIVE,
            description="a published number re-derived from the synthetic frame",
            citation="doi:10.0000/reanalysis",
            parents=("synthetic_frame",),
        )
    )
    graph = registry.freeze()

    assert render_node.origin is SourceOrigin.SYNTHETIC_RENDER
    assert graph.citable_as_biological_evidence() == ("paper_geometry",)

    first = graph.exclusion("synthetic_frame")
    assert first is not None
    assert first.reason is ExclusionReason.NON_EMPIRICAL_ORIGIN
    assert first.path == ("synthetic_frame",)

    second = graph.exclusion("published_reanalysis")
    assert second is not None
    assert second.reason is ExclusionReason.DERIVED_FROM_EXCLUDED
    assert second.blocking_source_id == "synthetic_frame"
    assert second.path == ("published_reanalysis", "synthetic_frame")


def test_a_projected_scalar_is_simulation_origin_and_therefore_also_excluded() -> None:
    registry = SourceRegistry()
    registry.register(
        SourceRecord(
            source_id="accepted_run",
            provenance=ProvenanceKind.SOURCED,
            origin=SourceOrigin.LITERATURE,
            status=SourceStatus.ACTIVE,
            description="literature prior the run started from",
            citation="doi:10.0000/fixture",
        )
    )
    node = _projection_operator().provenance_record(
        source_id="gate_scalar",
        parent_source_ids=("accepted_run",),
        description="peak cortex tension gate scalar",
    )
    registry.register(node)
    graph = registry.freeze()
    assert node.origin is SourceOrigin.SIMULATION
    assert graph.exclusion("gate_scalar") is not None
    assert "gate_scalar" not in graph.citable_as_biological_evidence()


# --------------------------------------------------------------------------------------------
# 7. the observation identity binds operator, configuration AND state
# --------------------------------------------------------------------------------------------


def test_the_provenance_digest_changes_with_the_optics_configuration() -> None:
    baseline = _fluorescence().provenance_digest(_emitter_state())
    wider = _fluorescence(optics=dataclasses.replace(OPTICS, exposure_s=4.0)).provenance_digest(
        _emitter_state()
    )
    assert baseline != wider


def test_the_provenance_digest_changes_with_the_state_it_observed() -> None:
    operator = _fluorescence()
    baseline = operator.provenance_digest(_emitter_state())
    other = operator.provenance_digest(_emitter_state(state_id="state/accepted-0008"))
    assert baseline != other


def test_the_state_digest_covers_the_array_content_and_not_only_its_identity() -> None:
    """External review 2026-07-29, finding 2 — this test previously asserted the DEFECT.

    It read ``first.provenance_digest == second.provenance_digest`` for two states differing in
    their arrays, and justified it as avoiding a host readback.  The review's counterexample is what
    that costs: with ``state_id``, ``manifest_hash``, ``source_id``, the field names and the units
    all equal, ``tension=[1, 2, 3]`` and ``tension=[1, 2, 30]`` observed to ``3.0`` and ``30.0`` and
    presented the same ``provenance_digest``, the same envelope ``artifact_id`` and the same
    ``observation_provenance_sha256``.  Two different observations with one artifact identity is
    worse than a host hash of an array the host already holds.
    """
    first = _array_state()
    second = _array_state(fields={"tension_by_component": ARRAY * 2.0})
    assert first.provenance_digest != second.provenance_digest
    assert first.content_sha256 != second.content_sha256
    # ... and identical content still gives identical digests, so the digest is content-ADDRESSED
    # rather than merely unstable.
    assert _array_state().provenance_digest == _array_state().provenance_digest
    assert _array_state(fields={"tension_by_component": ARRAY.copy()}).content_sha256 == (
        first.content_sha256
    )


def test_the_review_counterexample_two_states_one_identity_is_refused() -> None:
    """The reported counterexample, verbatim: 3.0 and 30.0 may not share an artifact identity."""
    operator = SummaryProjectionOperator(
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

    def state(values: list[float]) -> MechanisticState:
        return MechanisticState(
            state_id="state/x",
            manifest_hash=MANIFEST_HASH,
            fields={"tension": np.array(values)},
            field_units={"tension": "pN/um"},
            source_id="run",
        )

    low, high = state([1.0, 2.0, 3.0]), state([1.0, 2.0, 30.0])
    first, second = operator.observe(low), operator.observe(high)
    assert (first.scalar, second.scalar) == (3.0, 30.0)
    assert low.provenance_digest != high.provenance_digest
    assert first.envelope.artifact_id != second.envelope.artifact_id
    assert first.envelope.observation_provenance_sha256 != (
        second.envelope.observation_provenance_sha256
    )


def test_the_state_owns_its_content_so_a_later_mutation_cannot_forge_a_collision() -> None:
    """The second route to the same collision: the field mapping used to be copied shallowly."""
    mutable = np.array([1.0, 2.0, 3.0])
    state = MechanisticState(
        state_id="state/x",
        manifest_hash=MANIFEST_HASH,
        fields={"tension": mutable},
        field_units={"tension": "pN/um"},
        source_id="run",
    )
    before = state.provenance_digest
    mutable[2] = 999.0
    assert float(np.asarray(state.fields["tension"])[2]) == 3.0
    assert state.provenance_digest == before
    with pytest.raises(ValueError):
        state.fields["tension"][0] = 5.0


def test_a_list_field_is_frozen_into_a_tuple_rather_than_kept_by_reference() -> None:
    emitters = [EMITTERS[0]]
    state = _emitter_state(
        fields={"actin_emitters": emitters}, field_units={"actin_emitters": EMITTER_KIND}
    )
    emitters.append(EMITTERS[1])
    assert state.fields["actin_emitters"] == (EMITTERS[0],)


def test_the_content_digest_is_computed_at_construction_and_cannot_be_supplied() -> None:
    """A caller-supplied content address would be a claim, not a measurement."""
    state = _array_state()
    assert len(state.content_sha256) == 64
    with pytest.raises(TypeError, match="init=False"):
        dataclasses.replace(state, content_sha256="0" * 64)  # type: ignore[misc]


def test_a_field_that_cannot_be_content_addressed_is_refused_rather_than_skipped() -> None:
    """Fail closed: a value silently omitted from the hash is the defect, not the fix."""

    class Opaque:
        pass

    with pytest.raises(ValueError, match="no content address"):
        MechanisticState(
            state_id="s",
            manifest_hash=MANIFEST_HASH,
            fields={"thing": Opaque()},
            field_units={"thing": "1"},
            source_id="run",
        )
    with pytest.raises(ValueError, match="object-dtype"):
        MechanisticState(
            state_id="s",
            manifest_hash=MANIFEST_HASH,
            fields={"thing": np.array([Opaque()], dtype=object)},
            field_units={"thing": "1"},
            source_id="run",
        )


def test_the_content_token_uses_the_sidecar_array_vocabulary_rather_than_a_second_one() -> None:
    """One content address per array in this package: dtype string, shape, order, SHA-256 of bytes."""
    array = np.arange(6.0).reshape(3, 2)
    token = _content_token(array, what="tension")
    descriptor = ArraySidecarDescriptor(
        uri="mem://fixture",
        sha256=hashlib.sha256(np.ascontiguousarray(array).tobytes(order="C")).hexdigest(),
        media_type="application/x-npy",
        schema_id="fixture",
        dtype=array.dtype.str,
        value_semantics="cortex tension",
        value_unit="pN/um",
        axes=(
            ArrayAxis(name="time", length=3, semantics="accepted step"),
            ArrayAxis(name="component", length=2, semantics="component index"),
        ),
        byte_size=array.nbytes,
    )
    assert token["dtype"] == descriptor.dtype
    assert tuple(token["shape"]) == descriptor.shape
    assert token["sha256"] == descriptor.sha256
    assert token["array_order"] == descriptor.array_order


def test_a_state_must_declare_a_unit_for_every_field_it_carries() -> None:
    with pytest.raises(ValueError, match="declared unit or kind"):
        MechanisticState(
            state_id="s",
            manifest_hash=MANIFEST_HASH,
            fields={"x": 1.0},
            field_units={},
            source_id="run",
        )


def test_a_state_must_carry_a_real_manifest_hash() -> None:
    with pytest.raises(ValueError, match="SHA-256"):
        MechanisticState(
            state_id="s",
            manifest_hash="not-a-hash",
            fields={},
            field_units={},
            source_id="run",
        )


def test_a_requirement_must_declare_a_unit_for_every_field_it_requires() -> None:
    with pytest.raises(ValueError, match="declared unit or kind"):
        InputRequirement(fields=("a", "b"), field_units={"a": "um"}, description="d")
    with pytest.raises(ValueError, match="not required"):
        InputRequirement(fields=("a",), field_units={"a": "um", "b": "um"}, description="d")


def test_an_observation_carrying_a_non_finite_value_is_not_an_observation() -> None:
    result = _fluorescence().observe(_emitter_state())
    with pytest.raises(ValueError, match="finite"):
        ObservationResult(
            operator_id="x",
            values=np.array([np.nan]),
            units=PHOTON_UNIT,
            envelope=result.envelope,
        )

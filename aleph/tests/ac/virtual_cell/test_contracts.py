"""CPU structural tests for the project-wide virtual-cell contracts."""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from aleph.virtual_cell.contracts import (
    ApproximationErrorReport,
    ArtifactEnvelope,
    CellStateManifest,
    EvidenceSource,
    RepresentationDescriptor,
    RepresentationKind,
    SandboxExperimentCard,
    UncertaintySource,
    manifest_from_legacy_cell_state,
)

_HASH = "a" * 64


def _error_report(
    *,
    metric_id: str,
    measured_error: float,
) -> ApproximationErrorReport:
    return ApproximationErrorReport(
        report_sha256="b" * 64,
        source_blob_sha256="c" * 64,
        reduced_blob_sha256="d" * 64,
        metric_id=metric_id,
        requested_tolerance=max(measured_error, 1e-6),
        measured_error=measured_error,
    )


def _manifest(**overrides) -> CellStateManifest:
    values = {
        "manifest_id": "amoeboid-confined",
        "species": "Homo sapiens",
        "lineage": "amoeboid-immune",
        "identity": {"primary_or_line": "primary", "activation": "resting"},
        "biological_state": {"cycle": "G0", "polarity": "front-rear"},
        "environment": {"geometry": "confined", "matrix": "none"},
        "components": ("membrane", "cortex", "nucleus"),
        "connectors": ("membrane_erm_cortex", "cortex_linc_nucleus"),
        "observations": ("shape_time_series", "nuclear_strain"),
        "missing_components": ("polarity_signalling",),
        "missing_laws": ("fast_cortex_remodelling",),
        "provenance": ("public-dataset:test",),
    }
    values.update(overrides)
    return CellStateManifest(**values)


def test_manifest_hash_is_mapping_order_independent_and_not_cell_line_fixed() -> None:
    first = _manifest()
    second = _manifest(
        identity={"activation": "resting", "primary_or_line": "primary"},
        biological_state={"polarity": "front-rear", "cycle": "G0"},
        environment={"matrix": "none", "geometry": "confined"},
    )
    assert first.manifest_hash == second.manifest_hash
    assert first.lineage == "amoeboid-immune"


def test_manifest_hash_treats_component_and_connector_membership_as_sets() -> None:
    first = _manifest()
    second = _manifest(
        components=tuple(reversed(first.components)),
        connectors=tuple(reversed(first.connectors)),
        observations=tuple(reversed(first.observations)),
    )
    assert first.manifest_hash == second.manifest_hash


def test_manifest_keeps_missing_laws_first_class() -> None:
    manifest = _manifest()
    record = manifest.to_dict()
    assert record["missing_components"] == ["polarity_signalling"]
    assert record["missing_laws"] == ["fast_cortex_remodelling"]


def test_manifest_rejects_present_and_missing_component_overlap() -> None:
    with pytest.raises(ValueError, match="both present and missing"):
        _manifest(missing_components=("cortex",))


def test_approximate_representation_requires_error_and_fallback() -> None:
    with pytest.raises(ValueError, match="approximation_error_bound"):
        RepresentationDescriptor(
            kind=RepresentationKind.TENSOR_TRAIN,
            uncertainty_sources=(UncertaintySource.PARAMETER,),
            axes=("parameter", "observable"),
        )
    descriptor = RepresentationDescriptor(
        kind=RepresentationKind.TENSOR_TRAIN,
        uncertainty_sources=(UncertaintySource.PARAMETER, UncertaintySource.REPRESENTATION),
        axes=("parameter", "observable"),
        approximation_error_bound=1e-3,
        error_metric="relative_frobenius",
        error_report=_error_report(
            metric_id="relative_frobenius",
            measured_error=1e-3,
        ),
        fallback="empirical-ensemble",
    )
    assert descriptor.fallback == "empirical-ensemble"


def test_surrogate_requires_trust_region_and_cannot_claim_native_provenance() -> None:
    with pytest.raises(ValueError, match="trust_region"):
        RepresentationDescriptor(
            kind=RepresentationKind.SURROGATE,
            uncertainty_sources=(UncertaintySource.REPRESENTATION,),
            axes=("state", "time"),
            approximation_error_bound=0.02,
            error_metric="held_out_rollout_error",
            error_report=_error_report(
                metric_id="held_out_rollout_error",
                measured_error=0.02,
            ),
            fallback="native-evaluation",
        )
    descriptor = RepresentationDescriptor(
        kind=RepresentationKind.SURROGATE,
        uncertainty_sources=(UncertaintySource.REPRESENTATION,),
        axes=("state", "time"),
        approximation_error_bound=0.02,
        error_metric="held_out_rollout_error",
        error_report=_error_report(
            metric_id="held_out_rollout_error",
            measured_error=0.02,
        ),
        fallback="native-evaluation",
        trust_region={"distance_max": 0.3},
    )
    with pytest.raises(ValueError, match="cannot grant native evidence"):
        ArtifactEnvelope(
            artifact_id="bad",
            cell_state_manifest_hash=_HASH,
            evidence_source=EvidenceSource.NATIVE_ACCEPTED,
            representation=descriptor,
        )


def test_reduced_tensor_points_back_to_native_without_laundering_source() -> None:
    descriptor = RepresentationDescriptor(
        kind=RepresentationKind.TENSOR_TRAIN,
        uncertainty_sources=(UncertaintySource.REPRESENTATION,),
        axes=("time", "component", "mode"),
        approximation_error_bound=1e-4,
        error_metric="projection_replay_error",
        error_report=_error_report(
            metric_id="projection_replay_error",
            measured_error=1e-4,
        ),
        fallback="uncompressed-array",
    )
    envelope = ArtifactEnvelope(
        artifact_id="trace-tt",
        cell_state_manifest_hash=_HASH,
        evidence_source=EvidenceSource.REDUCED_MODEL,
        representation=descriptor,
        source_artifact_ids=("accepted-trace",),
    )
    assert envelope.evidence_source is EvidenceSource.REDUCED_MODEL
    with pytest.raises(ValueError, match="cannot grant native evidence"):
        ArtifactEnvelope(
            artifact_id="laundered",
            cell_state_manifest_hash=_HASH,
            evidence_source=EvidenceSource.NATIVE_ACCEPTED,
            representation=descriptor,
            source_artifact_ids=("accepted-trace",),
        )


def test_scalar_projection_needs_code_hash_and_source_artifact() -> None:
    descriptor = RepresentationDescriptor(
        kind=RepresentationKind.SCALAR_PROJECTION,
        uncertainty_sources=(UncertaintySource.NUMERICAL,),
        axes=("scalar",),
    )
    with pytest.raises(ValueError, match="projection_hash"):
        ArtifactEnvelope(
            artifact_id="scalar",
            cell_state_manifest_hash=_HASH,
            evidence_source=EvidenceSource.REDUCED_MODEL,
            representation=descriptor,
            source_artifact_ids=("source-array",),
        )


def test_generic_envelope_cannot_self_declare_native_authority() -> None:
    descriptor = RepresentationDescriptor(
        kind=RepresentationKind.EXPLICIT_STATE,
        uncertainty_sources=(UncertaintySource.NUMERICAL,),
        axes=("time",),
    )
    with pytest.raises(ValueError, match="cannot grant native evidence"):
        ArtifactEnvelope(
            artifact_id="invented",
            cell_state_manifest_hash=_HASH,
            evidence_source=EvidenceSource.NATIVE_ACCEPTED,
            representation=descriptor,
        )


def test_approximation_requires_a_source_and_reduced_blob_bound_error_report() -> None:
    with pytest.raises(ValueError, match="error_report"):
        RepresentationDescriptor(
            kind=RepresentationKind.TENSOR_TRAIN,
            uncertainty_sources=(UncertaintySource.REPRESENTATION,),
            axes=("time", "node", "dimension"),
            approximation_error_bound=1e-3,
            error_metric="relative_frobenius",
            fallback="dense",
        )


def test_approximation_error_report_metric_and_value_must_match_descriptor() -> None:
    with pytest.raises(ValueError, match="metric_id"):
        RepresentationDescriptor(
            kind=RepresentationKind.TENSOR_TRAIN,
            uncertainty_sources=(UncertaintySource.REPRESENTATION,),
            axes=("time", "node", "dimension"),
            approximation_error_bound=1e-3,
            error_metric="projection-aware",
            error_report=_error_report(
                metric_id="relative_frobenius",
                measured_error=1e-3,
            ),
            fallback="dense",
        )


def test_absolute_error_report_can_use_physical_units_and_tolerance_above_one() -> None:
    report = ApproximationErrorReport(
        report_sha256="b" * 64,
        source_blob_sha256="c" * 64,
        reduced_blob_sha256="d" * 64,
        metric_id="max_force_error",
        requested_tolerance=2.0,
        measured_error=1.5,
        metric_unit="pN",
        relative_metric=False,
    )
    assert report.metric_unit == "pN"


def test_approximate_artifact_requires_source_ancestry() -> None:
    descriptor = RepresentationDescriptor(
        kind=RepresentationKind.TENSOR_TRAIN,
        uncertainty_sources=(UncertaintySource.REPRESENTATION,),
        axes=("time", "node", "dimension"),
        approximation_error_bound=1e-3,
        error_metric="relative_frobenius",
        error_report=_error_report(
            metric_id="relative_frobenius",
            measured_error=1e-3,
        ),
        fallback="dense",
    )
    with pytest.raises(ValueError, match="source_artifact_ids"):
        ArtifactEnvelope(
            artifact_id="orphan-reduction",
            cell_state_manifest_hash=_HASH,
            evidence_source=EvidenceSource.REDUCED_MODEL,
            representation=descriptor,
        )


def test_manifest_scalar_keys_cannot_collide_after_normalization() -> None:
    with pytest.raises(ValueError, match="collide"):
        _manifest(identity={"source": "a", " source ": "b"})


def test_sandbox_card_requires_every_control_class_and_expansion_route() -> None:
    with pytest.raises(ValueError, match="adversarial_controls"):
        SandboxExperimentCard(
            experiment_id="S0",
            question="Is the regime numerical?",
            manifest_hashes=("abc",),
            intervention="vary dt",
            observables=("power_spectrum",),
            positive_controls=("known_limit_cycle",),
            negative_controls=("fixed_point",),
            adversarial_controls=(),
            held_out_interventions=("unseen_iteration_budget",),
            falsifier="period follows dt",
            failure_expansions=("integrator-study",),
            budget_class="cpu-oracle",
        )


def test_wet_lab_card_is_deferred_without_information_gain_rationale() -> None:
    values = {
        "experiment_id": "future-wet",
        "question": "Which observation separates two mechanisms?",
        "manifest_hashes": ("abc",),
        "intervention": "deferred",
        "observables": ("candidate_observable",),
        "positive_controls": ("simulation-positive",),
        "negative_controls": ("simulation-negative",),
        "adversarial_controls": ("batch-confounder",),
        "held_out_interventions": ("unseen-dose",),
        "falsifier": "information gain is below numerical uncertainty",
        "failure_expansions": ("additional-modality",),
        "budget_class": "not-scheduled",
        "wet_lab": True,
    }
    with pytest.raises(ValueError, match="information_gain_rationale"):
        SandboxExperimentCard(**values)
    card = SandboxExperimentCard(
        **values,
        information_gain_rationale="sandbox posterior leaves two separated hypotheses",
    )
    assert card.wet_lab


@dataclass(frozen=True)
class _LegacyState:
    axes: tuple[str, ...] = (
        "generic-epithelial",
        "epithelial",
        "G1",
        "matrix-adherent",
        "polarized",
        "isotonic",
    )


def test_legacy_runtime_state_wraps_without_becoming_universal_schema() -> None:
    manifest = manifest_from_legacy_cell_state(
        _LegacyState(),
        manifest_id="legacy-adapter-test",
        species="Homo sapiens",
        components=("membrane", "cortex"),
        connectors=("membrane_erm_cortex",),
        provenance=("legacy-runtime",),
    )
    assert manifest.identity["legacy_runtime_schema"] == "ac.engine.CellState@6-axis"
    assert manifest.biological_state["emt"] == "epithelial"

"""Proof tests for the A1-A15 audit register.

Two kinds of test live here and they are not interchangeable:

1. **Register invariants** — the register cannot claim a status its content does not support, every
   enforcement target still resolves to the object that was bound, and every ``proof_test`` the
   register names exists in this module.
2. **Enforcement proofs** — for each item the register calls MECHANICALLY_ENFORCED, a test that
   actually fires the enforcement: the raise is asserted, or the refusal is asserted.  A binding
   whose proof does not fire is a downgrade to UNENFORCED, not a comment.

Nothing here runs physics, touches a GPU, or reads a run record.
"""

from __future__ import annotations

import dataclasses
import hashlib
import importlib
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from aleph.virtual_cell.archetypes import ArchetypeRing, MissingLaw, MissingLawVerdict
from aleph.virtual_cell.audit_checklist import (
    AUDIT_IDS,
    AUDIT_REGISTER,
    REMOVAL_IDS,
    AuditItem,
    AuditKind,
    AuditRegisterError,
    AuditStatus,
    EnforcementBinding,
    Measurement,
    audit_item,
    audit_register_to_dict,
    enforcement_target_paths,
    format_audit_register,
    items_with_status,
    resolve_enforcement_targets,
    status_counts,
    summarize_audit_register,
    unenforced_gap_queue,
)
from aleph.virtual_cell.contracts import (
    ApproximationErrorReport,
    ArtifactEnvelope,
    CellStateManifest,
    EvidenceSource,
    RepresentationDescriptor,
    RepresentationKind,
    SandboxExperimentCard,
    UncertaintySource,
)
from aleph.virtual_cell.grammar import CensusUnavailableError, load_engine_grammar
from aleph.virtual_cell.source_registry import (
    ExclusionReason,
    HeldOutValidationRefused,
    ProvenanceKind,
    SourceOrigin,
    SourceRecord,
    SourceRegistry,
    SourceStatus,
)
from aleph.virtual_cell.surrogate import (
    CoverageReport,
    GateAuthorityError,
    GateEvidenceLedger,
    NativeQueryRequest,
    SurrogateAuthorityError,
    SurrogatePrediction,
    SurrogateResponse,
    SurrogateVerdict,
)
from aleph.virtual_cell.surrogate import ExpansionRoute as SurrogateExpansionRoute
from aleph.virtual_cell.telemetry import AxisOrderSearchReport, search_axis_order


def _sha(label: str) -> str:
    return hashlib.sha256(label.encode()).hexdigest()


def _plain_representation(
    kind: RepresentationKind = RepresentationKind.EXPLICIT_STATE,
) -> RepresentationDescriptor:
    return RepresentationDescriptor(
        kind=kind,
        uncertainty_sources=(UncertaintySource.NUMERICAL,),
        axes=("x",),
    )


def _error_report(tolerance: float = 0.1, measured: float = 0.05) -> ApproximationErrorReport:
    return ApproximationErrorReport(
        report_sha256=_sha("report"),
        source_blob_sha256=_sha("source"),
        reduced_blob_sha256=_sha("reduced"),
        metric_id="relative-frobenius",
        requested_tolerance=tolerance,
        measured_error=measured,
    )


def _prediction() -> SurrogatePrediction:
    coverage = CoverageReport(
        label="audit-register-fixture",
        nominal_level=0.95,
        empirical_coverage=0.946,
        n_draws=20_000,
        standard_error=0.0015,
        z_score=-2.6,
        sigma_gate=3.0,
        seed=20260729,
    )
    return SurrogatePrediction(
        surrogate_id="fixture-surrogate",
        axes=("x",),
        query_point=(0.5,),
        mean=1.0,
        predictive_sigma=0.1,
        trust_region={"x_lo": 0.0, "x_hi": 1.0},
        support_distance=0.01,
        coverage=coverage,
    )


# --------------------------------------------------------------------------------------------
# 0. the register does not import a GPU runtime
# --------------------------------------------------------------------------------------------


def test_importing_the_register_does_not_initialise_warp() -> None:
    """Measured in a fresh interpreter: an in-process ``sys.modules`` check would test the session."""
    script = textwrap.dedent(
        """
        import sys
        from aleph.virtual_cell.audit_checklist import format_audit_register, status_counts

        format_audit_register()
        status_counts()
        for blocked in ("warp", "torch", "jax"):
            assert blocked not in sys.modules, f"audit register import pulled in {blocked}"
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
# 1. register invariants — a status is derived, never asserted
# --------------------------------------------------------------------------------------------


def test_the_register_covers_a1_to_a15_and_the_three_removals() -> None:
    assert tuple(f"A{index}" for index in range(1, 16)) == AUDIT_IDS
    assert REMOVAL_IDS == ("R1", "R2", "R3")
    assert set(AUDIT_REGISTER) == set(AUDIT_IDS) | set(REMOVAL_IDS)
    assert len(AUDIT_REGISTER) == 18


def test_status_is_derived_from_content_and_the_counts_partition_the_register() -> None:
    counts = status_counts()
    assert set(counts) == set(AuditStatus)
    assert sum(counts.values()) == len(AUDIT_REGISTER)
    for item in AUDIT_REGISTER.values():
        if item.enforcement:
            assert item.status is AuditStatus.MECHANICALLY_ENFORCED
        elif item.measurements:
            assert item.status is AuditStatus.EVIDENCED
        else:
            assert item.status is AuditStatus.UNENFORCED
    # Status is a read-only property: a row cannot be told what status to have.
    assert isinstance(AuditItem.status, property)
    assert AuditItem.status.fset is None
    with pytest.raises(AttributeError):
        object.__setattr__(audit_item("A5"), "status", AuditStatus.MECHANICALLY_ENFORCED)


def test_every_enforced_item_names_a_raise_and_a_residual_gap() -> None:
    enforced = items_with_status(AuditStatus.MECHANICALLY_ENFORCED)
    assert enforced, "the register must record at least one proven enforcement"
    for item in enforced:
        assert item.residual_gap is not None
        assert any(binding.raises is not None for binding in item.enforcement)
        assert not item.partial_enforcement


def test_an_enforcement_claim_without_a_residual_gap_is_refused() -> None:
    binding = EnforcementBinding(
        target=SandboxExperimentCard,
        guarantee="a card needs adversarial controls",
        proof_test="test_a7_a_sandbox_card_without_adversarial_controls_is_refused",
        raises=ValueError,
    )
    with pytest.raises(ValueError, match="residual gap"):
        AuditItem(
            audit_id="A7",
            kind=AuditKind.EXPANSION,
            attack="a",
            verdict="v",
            expansion_route="e",
            enforcement=(binding,),
        )


def test_an_enforcement_that_never_raises_is_not_enforcement() -> None:
    binding = EnforcementBinding(
        target=CellStateManifest,
        attribute="manifest_hash",
        guarantee="the hash covers environment",
        proof_test="test_a11_passage_and_substrate_change_the_manifest_hash",
    )
    with pytest.raises(ValueError, match="names an exception"):
        AuditItem(
            audit_id="A11",
            kind=AuditKind.EXPANSION,
            attack="a",
            verdict="v",
            expansion_route="e",
            enforcement=(binding,),
            residual_gap="something",
        )


def test_partial_machinery_may_not_masquerade_as_enforcement() -> None:
    binding = EnforcementBinding(
        target=SandboxExperimentCard,
        guarantee="a card needs adversarial controls",
        proof_test="test_a7_a_sandbox_card_without_adversarial_controls_is_refused",
        raises=ValueError,
    )
    with pytest.raises(ValueError, match="belongs in residual_gap"):
        AuditItem(
            audit_id="A7",
            kind=AuditKind.EXPANSION,
            attack="a",
            verdict="v",
            expansion_route="e",
            enforcement=(binding,),
            residual_gap="open",
            partial_enforcement=("something adjacent",),
        )
    with pytest.raises(ValueError, match="the whole item is the gap"):
        AuditItem(
            audit_id="A5",
            kind=AuditKind.EXPANSION,
            attack="a",
            verdict="v",
            expansion_route="e",
            residual_gap="a gap without enforcement is the item",
        )


def test_an_unknown_or_miskinded_audit_id_is_refused() -> None:
    with pytest.raises(ValueError, match="A1..A15"):
        AuditItem(
            audit_id="A16", kind=AuditKind.EXPANSION, attack="a", verdict="v", expansion_route="e"
        )
    with pytest.raises(ValueError, match="prohibition"):
        AuditItem(
            audit_id="R1", kind=AuditKind.EXPANSION, attack="a", verdict="v", expansion_route="e"
        )
    with pytest.raises(KeyError, match="unknown audit id"):
        audit_item("A99")


def test_a_binding_outside_this_layer_or_to_a_missing_attribute_is_refused() -> None:
    with pytest.raises(ValueError, match="outside aleph/virtual_cell"):
        EnforcementBinding(target=dict, guarantee="g", proof_test="test_x", raises=ValueError)
    with pytest.raises(AttributeError, match="no longer exposes"):
        EnforcementBinding(
            target=SandboxExperimentCard,
            attribute="renamed_away",
            guarantee="g",
            proof_test="test_x",
        )
    with pytest.raises(ValueError, match="must name a test function"):
        EnforcementBinding(target=SandboxExperimentCard, guarantee="g", proof_test="proves_it")


def test_every_enforcement_target_still_resolves_to_the_object_that_was_bound() -> None:
    """A rename must break the register loudly instead of silently downgrading an item."""
    resolved = resolve_enforcement_targets()
    assert set(resolved) == set(enforcement_target_paths())
    for audit_id, objects in resolved.items():
        item = audit_item(audit_id)
        assert len(objects) == len(item.enforcement)
        for found, binding in zip(objects, item.enforcement, strict=True):
            assert found is binding.target
    stale = EnforcementBinding(
        target=SandboxExperimentCard, guarantee="g", proof_test="test_x", raises=ValueError
    )
    object.__setattr__(stale, "target", _bogus_target)
    with pytest.raises(AuditRegisterError):
        stale.resolve()


class _bogus_target:  # noqa: N801 - stands in for a symbol that no longer exists at its path
    __qualname__ = "definitely_not_exported"
    __module__ = "aleph.virtual_cell.contracts"


def test_every_proof_test_named_by_the_register_exists_and_is_unambiguous() -> None:
    """Every named proof test resolves to exactly one real callable, wherever in the lane it lives.

    This used to search only THIS module, which was right while the register only covered items
    whose enforcement lived here.  A13's enforcement lives in the lane that owns the objects it
    binds (`cell_state_posterior`, `observation_operator`), and its proof tests depend on that
    lane's fixtures, so requiring them here would have meant duplicating the setup — and a copied
    proof drifts from the thing it claims to prove.

    Broadening the search would weaken the guard on its own, so the check is tightened in the same
    change: a proof test must resolve to exactly ONE definition across the whole lane.  Previously
    two modules could each define `test_foo` and the register would silently bind to whichever it
    found; now that is an error.  Net, this checks strictly more than it did.
    """
    package_dir = Path(__file__).resolve().parent
    modules = {}
    for path in sorted(package_dir.glob("test_*.py")):
        name = f"aleph.tests.ac.virtual_cell.{path.stem}"
        modules[path.stem] = importlib.import_module(name)

    for item in AUDIT_REGISTER.values():
        for binding in item.enforcement:
            found = [
                stem
                for stem, module in modules.items()
                if callable(getattr(module, binding.proof_test, None))
            ]
            assert found, (
                f"{item.audit_id} names proof {binding.proof_test!r}, which does not exist "
                f"anywhere in {package_dir.name}/"
            )
            assert len(found) == 1, (
                f"{item.audit_id} names proof {binding.proof_test!r}, which is defined in "
                f"{len(found)} modules ({', '.join(found)}); the register would bind to an "
                f"arbitrary one"
            )


def test_a_measurement_carries_units_provenance_and_a_finite_value() -> None:
    for item in AUDIT_REGISTER.values():
        for measurement in item.measurements:
            assert np.isfinite(measurement.value)
            assert measurement.units
            assert measurement.measured_in
    with pytest.raises(ValueError, match="finite"):
        Measurement(
            quantity="q",
            value=float("inf"),
            units="1",
            protocol="p",
            measured_in="m",
        )


# --------------------------------------------------------------------------------------------
# 2. the gap queue is reported, not hidden
# --------------------------------------------------------------------------------------------


def test_the_unenforced_queue_lists_every_gap_with_its_expansion_route() -> None:
    queue = unenforced_gap_queue()
    assert queue, "an audit register with no visible gaps is the failure it exists to prevent"
    unenforced_ids = {entry.audit_id for entry in queue}
    assert unenforced_ids == {
        item.audit_id for item in AUDIT_REGISTER.values() if item.status is AuditStatus.UNENFORCED
    }
    for entry in queue:
        assert entry.expansion_route
        assert audit_item(entry.audit_id).enforcement == ()


def test_the_printed_register_shows_every_unenforced_id_and_the_counts() -> None:
    summary = summarize_audit_register()
    text = format_audit_register()
    assert "UNENFORCED" in text
    for entry in summary.unenforced:
        assert entry.audit_id in text
        assert entry.expansion_route in text
    for audit_id in summary.enforced_ids:
        assert audit_id in text
    assert str(summary.total) in text
    assert "closes no gate" in text


def test_the_register_serializes_with_its_derived_statuses() -> None:
    record = audit_register_to_dict()
    assert len(record["items"]) == 18
    statuses = {entry["audit_id"]: entry["status"] for entry in record["items"]}
    for audit_id, status in statuses.items():
        assert status == audit_item(audit_id).status.value
    assert record["summary"]["total"] == 18
    assert "oracle-only" in record["evidence_class"]


# --------------------------------------------------------------------------------------------
# 3. A3 — rank is representation complexity, and the search reports what it did not search
# --------------------------------------------------------------------------------------------


def test_a3_axis_order_report_refuses_a_false_exhaustive_claim() -> None:
    with pytest.raises(ValueError, match="exhaustive search must evaluate every ordering"):
        AxisOrderSearchReport(
            ndim=2,
            total_orderings=2,
            evaluated_orderings=((0, 1),),
            strategy="exhaustive",
            strategy_description="claimed exhaustive after one ordering",
            requested_tolerance=0.0,
            best_order=(0, 1),
            best_storage_size=4,
            best_relative_error=0.0,
            dense_parameter_count=8,
            any_order_met_tolerance=True,
            storage_by_order={(0, 1): 4},
            error_by_order={(0, 1): 0.0},
        )


def test_a3_capped_axis_order_search_reports_its_unsearched_fraction() -> None:
    tensor = np.arange(2 * 3 * 4 * 5, dtype=np.float64).reshape(2, 3, 4, 5)
    capped = search_axis_order(tensor, relative_tolerance=1e-9, max_orders=3)
    assert capped.strategy == "heuristic-subset"
    assert capped.total_orderings == 24
    assert capped.evaluated_count == 3
    assert capped.unsearched_fraction == pytest.approx(21 / 24)
    assert capped.to_dict()["unsearched_fraction"] == capped.unsearched_fraction

    exhaustive = search_axis_order(tensor, relative_tolerance=1e-9)
    assert exhaustive.unsearched_fraction == 0.0
    distinct = sorted(set(exhaustive.storage_by_order.values()))
    # The register stores this count as a reproduced measurement: one array, one tolerance,
    # 24 orderings, and the storage the "rank" implies is not one number.
    assert len(distinct) == int(audit_item("A3").measurements[0].value)
    assert distinct[0] == 38 and distinct[-1] == 46


# --------------------------------------------------------------------------------------------
# 4. A4 / R3 — a surrogate cannot certify itself, transitively
# --------------------------------------------------------------------------------------------


def test_a4_a_surrogate_prediction_cannot_be_relabelled_as_native() -> None:
    prediction = _prediction()
    assert prediction.evidence_source is EvidenceSource.SURROGATE
    with pytest.raises(TypeError, match="evidence_source"):
        dataclasses.replace(prediction, evidence_source=EvidenceSource.NATIVE_ACCEPTED)  # type: ignore[call-arg]
    assert "evidence_source" not in {f.name for f in dataclasses.fields(prediction)}


def test_a4_a_refused_surrogate_query_has_no_value() -> None:
    request = NativeQueryRequest(
        request_id="q-1",
        reason=SurrogateVerdict.REFUSED_OUTSIDE_TRUST_REGION,
        rationale="query lies outside the fitted box",
        axes=("x",),
        query_point=(12.0,),
        diagnostic={"box_excess": 2.543},
        expansions=(SurrogateExpansionRoute.NATIVE_EVALUATION,),
    )
    response = SurrogateResponse(
        verdict=SurrogateVerdict.REFUSED_OUTSIDE_TRUST_REGION, native_request=request
    )
    with pytest.raises(SurrogateAuthorityError, match="a refusal has no value"):
        _ = response.value


def test_a4_r3_surrogate_provenance_is_refused_transitively() -> None:
    surrogate_envelope = _prediction().as_artifact_envelope(
        artifact_id="surrogate-1",
        cell_state_manifest_hash=_sha("manifest"),
        source_artifact_ids=("oracle-target-1",),
        error_report=_error_report(),
    )
    derived = ArtifactEnvelope(
        artifact_id="derived-1",
        cell_state_manifest_hash=_sha("manifest"),
        evidence_source=EvidenceSource.REDUCED_MODEL,
        representation=_plain_representation(),
        source_artifact_ids=("surrogate-1",),
    )
    grandchild = ArtifactEnvelope(
        artifact_id="claim-1",
        cell_state_manifest_hash=_sha("manifest"),
        evidence_source=EvidenceSource.ANALYTIC_ORACLE,
        representation=_plain_representation(),
        source_artifact_ids=("derived-1",),
    )
    root = ArtifactEnvelope(
        artifact_id="oracle-target-1",
        cell_state_manifest_hash=_sha("manifest"),
        evidence_source=EvidenceSource.ANALYTIC_ORACLE,
        representation=_plain_representation(),
    )
    ledger = GateEvidenceLedger()
    for envelope in (root, surrogate_envelope, derived, grandchild):
        ledger.register(envelope)

    assert ledger.surrogate_ancestry("claim-1") == ("claim-1", "derived-1", "surrogate-1")
    with pytest.raises(SurrogateAuthorityError) as excinfo:
        ledger.assert_gate_closable("claim-1", gate_id="GATE-B")
    assert "claim-1 -> derived-1 -> surrogate-1" in str(excinfo.value)

    orphan = ArtifactEnvelope(
        artifact_id="orphan-1",
        cell_state_manifest_hash=_sha("manifest"),
        evidence_source=EvidenceSource.ANALYTIC_ORACLE,
        representation=_plain_representation(),
        source_artifact_ids=("never-registered",),
    )
    ledger.register(orphan)
    with pytest.raises(GateAuthorityError, match="unverifiable ancestry fails closed"):
        ledger.assert_gate_closable("orphan-1", gate_id="GATE-B")


def test_r3_synthetic_and_surrogate_representations_keep_their_provenance() -> None:
    synthetic = RepresentationDescriptor(
        kind=RepresentationKind.SYNTHETIC_IMAGE,
        uncertainty_sources=(UncertaintySource.MEASUREMENT,),
        axes=("row", "col"),
    )
    with pytest.raises(ValueError, match="must retain SYNTHETIC evidence provenance"):
        ArtifactEnvelope(
            artifact_id="image-1",
            cell_state_manifest_hash=_sha("manifest"),
            evidence_source=EvidenceSource.PUBLIC_EXPERIMENT,
            representation=synthetic,
            source_artifact_ids=("state-1",),
            observation_provenance_sha256=_sha("observation"),
        )
    surrogate = RepresentationDescriptor(
        kind=RepresentationKind.SURROGATE,
        uncertainty_sources=(UncertaintySource.REPRESENTATION,),
        axes=("x",),
        approximation_error_bound=0.05,
        error_metric="relative-frobenius",
        error_report=_error_report(),
        fallback="native-evaluation",
        trust_region={"x_lo": 0.0, "x_hi": 1.0},
    )
    with pytest.raises(ValueError, match="must retain SURROGATE evidence provenance"):
        ArtifactEnvelope(
            artifact_id="surrogate-2",
            cell_state_manifest_hash=_sha("manifest"),
            evidence_source=EvidenceSource.PUBLIC_EXPERIMENT,
            representation=surrogate,
            source_artifact_ids=("state-1",),
        )


# --------------------------------------------------------------------------------------------
# 5. A6 / A12 — synthetic and computed sources cannot score themselves
# --------------------------------------------------------------------------------------------


def _provenance_graph() -> SourceRegistry:
    registry = SourceRegistry()
    registry.register(
        SourceRecord(
            source_id="toy-sim",
            provenance=ProvenanceKind.CONVENIENCE,
            origin=SourceOrigin.SIMULATION,
            status=SourceStatus.ACTIVE,
            description="toy simulation state used as a renderer input",
        )
    )
    registry.register(
        SourceRecord(
            source_id="render",
            provenance=ProvenanceKind.DERIVED,
            origin=SourceOrigin.SYNTHETIC_RENDER,
            status=SourceStatus.ACTIVE,
            description="synthetic microscopy render of that state",
            parents=("toy-sim",),
        )
    )
    registry.register(
        SourceRecord(
            source_id="seg-stats",
            provenance=ProvenanceKind.DERIVED,
            origin=SourceOrigin.PUBLIC_EXPERIMENT,
            status=SourceStatus.ACTIVE,
            citation="accession:FIXTURE-1",
            description="segmentation statistic computed on the render",
            parents=("render",),
        )
    )
    registry.register(
        SourceRecord(
            source_id="morph-claim",
            provenance=ProvenanceKind.DERIVED,
            origin=SourceOrigin.PUBLIC_EXPERIMENT,
            status=SourceStatus.ACTIVE,
            citation="accession:FIXTURE-2",
            description="morphology claim two hops below the renderer",
            parents=("seg-stats",),
            held_out_split_id="split-morph",
        )
    )
    registry.register(
        SourceRecord(
            source_id="calibration-paper",
            provenance=ProvenanceKind.CALIBRATION_ONLY,
            origin=SourceOrigin.LITERATURE,
            status=SourceStatus.ACTIVE,
            citation="doi:10.0000/fixture-calibration",
            description="a source already used to calibrate the model",
            held_out_split_id="split-calibration",
        )
    )
    registry.register(
        SourceRecord(
            source_id="calibration-child",
            provenance=ProvenanceKind.DERIVED,
            origin=SourceOrigin.LITERATURE,
            status=SourceStatus.ACTIVE,
            citation="doi:10.0000/fixture-derived",
            description="a claim derived from the calibration source",
            parents=("calibration-paper",),
            held_out_split_id="split-calibration-child",
        )
    )
    registry.register(
        SourceRecord(
            source_id="qm-binding-free-energy",
            provenance=ProvenanceKind.COMPUTED_QM_MD,
            origin=SourceOrigin.COMPUTATION,
            status=SourceStatus.ACTIVE,
            description="DFT/MD molecular binding free energy",
            held_out_split_id="split-qm",
        )
    )
    registry.register(
        SourceRecord(
            source_id="wet-lab-paper",
            provenance=ProvenanceKind.SOURCED,
            origin=SourceOrigin.LITERATURE,
            status=SourceStatus.ACTIVE,
            citation="doi:10.0000/fixture-measurement",
            description="a real measurement in the literature",
            held_out_split_id="split-wet",
        )
    )
    return registry


def test_a6_synthetic_ancestry_is_excluded_two_hops_down() -> None:
    graph = _provenance_graph().freeze()
    exclusion = graph.exclusion("morph-claim")
    assert exclusion is not None
    assert exclusion.reason is ExclusionReason.DERIVED_FROM_EXCLUDED
    assert exclusion.path == ("morph-claim", "seg-stats", "render")
    assert exclusion.blocking_source_id == "render"
    citable = graph.citable_as_biological_evidence()
    assert "morph-claim" not in citable and "seg-stats" not in citable
    assert "wet-lab-paper" in citable


def test_a6_synthetic_and_calibration_descendants_are_refused_as_held_out() -> None:
    graph = _provenance_graph().freeze()
    with pytest.raises(HeldOutValidationRefused, match="derived-from-excluded"):
        graph.assert_admissible_as_held_out_validation("morph-claim")
    with pytest.raises(HeldOutValidationRefused, match="CALIBRATION_ONLY"):
        graph.assert_admissible_as_held_out_validation("calibration-child")
    assert graph.admissible_held_out_validation_ids() == ("wet-lab-paper",)


def test_a12_a_computed_qm_md_source_is_not_biological_evidence() -> None:
    graph = _provenance_graph().freeze()
    exclusion = graph.exclusion("qm-binding-free-energy")
    assert exclusion is not None
    assert exclusion.reason is ExclusionReason.NON_EMPIRICAL_ORIGIN
    assert "qm-binding-free-energy" not in graph.citable_as_biological_evidence()
    with pytest.raises(HeldOutValidationRefused, match="non-empirical-origin"):
        graph.assert_admissible_as_held_out_validation("qm-binding-free-energy")


# --------------------------------------------------------------------------------------------
# 6. A7 / A9 — a sandbox result cannot become a biological claim on its own
# --------------------------------------------------------------------------------------------


def _card_fields() -> dict[str, object]:
    return {
        "experiment_id": "S0-numerical-invariance",
        "question": "does the accepted trace depend on a numerical knob?",
        "manifest_hashes": (_sha("manifest"),),
        "intervention": "halve the inner tolerance",
        "observables": ("cortical_tension",),
        "positive_controls": ("analytic-oracle",),
        "negative_controls": ("unchanged-configuration",),
        "adversarial_controls": ("seed-permuted-replicate",),
        "held_out_interventions": ("unseen-substrate-stiffness",),
        "falsifier": "the observable moves by more than the declared tolerance",
        "failure_expansions": ("model-inadequacy",),
        "budget_class": "cpu-oracle",
    }


def test_a7_a_sandbox_card_without_adversarial_controls_is_refused() -> None:
    assert SandboxExperimentCard(**_card_fields()).card_hash  # type: ignore[arg-type]
    for omitted in (
        "adversarial_controls",
        "held_out_interventions",
        "negative_controls",
        "positive_controls",
        "failure_expansions",
    ):
        fields = _card_fields()
        fields[omitted] = ()
        with pytest.raises(ValueError, match=f"{omitted} must not be empty"):
            SandboxExperimentCard(**fields)  # type: ignore[arg-type]
    fields = _card_fields()
    fields["falsifier"] = "   "
    with pytest.raises(ValueError, match="falsifier"):
        SandboxExperimentCard(**fields)  # type: ignore[arg-type]


def test_a7_a9_no_generic_layer_artifact_can_close_a_gate() -> None:
    oracle = ArtifactEnvelope(
        artifact_id="oracle-only-1",
        cell_state_manifest_hash=_sha("manifest"),
        evidence_source=EvidenceSource.ANALYTIC_ORACLE,
        representation=_plain_representation(),
    )
    ledger = GateEvidenceLedger()
    ledger.register(oracle)
    assert ledger.surrogate_ancestry("oracle-only-1") == ()
    with pytest.raises(GateAuthorityError, match="unverified observer") as excinfo:
        ledger.assert_gate_closable("oracle-only-1", gate_id="GATE-B")
    assert not isinstance(excinfo.value, SurrogateAuthorityError)


def test_a9_the_generic_layer_cannot_declare_native_evidence() -> None:
    for source in (EvidenceSource.NATIVE_ACCEPTED, EvidenceSource.NATIVE_ATTEMPTED):
        with pytest.raises(ValueError, match="cannot grant native evidence authority"):
            ArtifactEnvelope(
                artifact_id="self-declared-native",
                cell_state_manifest_hash=_sha("manifest"),
                evidence_source=source,
                representation=_plain_representation(),
            )


# --------------------------------------------------------------------------------------------
# 7. A8 — the grammar is not called universal, and the census has no fallback
# --------------------------------------------------------------------------------------------


def test_a8_a_refusal_without_a_named_law_is_not_a_verdict() -> None:
    with pytest.raises(ValueError, match="generic failure is not a verdict"):
        MissingLawVerdict(
            archetype_id="erythrocyte",
            ring=ArchetypeRing.RING_B,
            missing_components=("membrane_skeleton",),
            missing_laws=(),
            unsupported_claims=("shear-driven deformation",),
        )
    named = MissingLawVerdict(
        archetype_id="erythrocyte",
        ring=ArchetypeRing.RING_B,
        missing_components=("membrane_skeleton",),
        missing_laws=(
            MissingLaw(
                law_id="law.membrane_in_plane_shear_elasticity",
                statement="the membrane needs an in-plane shear modulus",
                blocked_by="a fluid Helfrich surface has in-plane shear modulus zero by construction",
            ),
        ),
        unsupported_claims=("shear-driven deformation",),
    )
    assert named.refused is True


def test_a8_the_census_has_no_fallback() -> None:
    with pytest.raises(CensusUnavailableError, match="census source not found"):
        load_engine_grammar(contracts_path="/nonexistent/engine/contracts.py")


# --------------------------------------------------------------------------------------------
# 8. A11 — evidence is keyed by a manifest, never by a cell-type label
# --------------------------------------------------------------------------------------------


def test_a11_an_artifact_cannot_be_keyed_by_a_cell_type_label() -> None:
    for label in ("MCF7", "MDA-MB-231", "HeLa"):
        with pytest.raises(ValueError, match="must be a lowercase SHA-256 digest"):
            ArtifactEnvelope(
                artifact_id="labelled-1",
                cell_state_manifest_hash=label,
                evidence_source=EvidenceSource.ANALYTIC_ORACLE,
                representation=_plain_representation(),
            )


def test_a11_passage_and_substrate_change_the_manifest_hash() -> None:
    def manifest(*, passage: int, substrate_kpa: float, cycle: str) -> CellStateManifest:
        return CellStateManifest(
            manifest_id="fibroblast-adherent",
            species="human",
            lineage="fibroblast",
            identity={"passage": passage},
            biological_state={"cycle": cycle},
            environment={"substrate_kpa": substrate_kpa},
            components=("cortex",),
            connectors=(),
        )

    base = manifest(passage=6, substrate_kpa=10.0, cycle="G1")
    assert base.manifest_hash == manifest(passage=6, substrate_kpa=10.0, cycle="G1").manifest_hash
    variants = {
        base.manifest_hash,
        manifest(passage=22, substrate_kpa=10.0, cycle="G1").manifest_hash,
        manifest(passage=6, substrate_kpa=40.0, cycle="G1").manifest_hash,
        manifest(passage=6, substrate_kpa=10.0, cycle="M").manifest_hash,
    }
    assert len(variants) == 4


# --------------------------------------------------------------------------------------------
# 9. R1 / R2 — the removed representations cannot be declared
# --------------------------------------------------------------------------------------------


def test_r2_a_full_density_representation_cannot_be_declared() -> None:
    with pytest.raises(ValueError, match="not a valid RepresentationKind"):
        RepresentationKind("full-density")
    with pytest.raises(TypeError, match="kind must be a RepresentationKind"):
        RepresentationDescriptor(
            kind="full-density",  # type: ignore[arg-type]
            uncertainty_sources=(UncertaintySource.NUMERICAL,),
            axes=("x",),
        )


def test_r2_a_gaussian_representation_needs_measured_error_and_a_fallback() -> None:
    with pytest.raises(ValueError, match="needs a finite nonnegative approximation_error_bound"):
        RepresentationDescriptor(
            kind=RepresentationKind.LOCAL_GAUSSIAN,
            uncertainty_sources=(UncertaintySource.THERMAL, UncertaintySource.BIOLOGICAL),
            axes=("x",),
        )
    with pytest.raises(ValueError, match="fallback"):
        RepresentationDescriptor(
            kind=RepresentationKind.LOCAL_GAUSSIAN,
            uncertainty_sources=(UncertaintySource.THERMAL,),
            axes=("x",),
            approximation_error_bound=0.05,
            error_metric="relative-frobenius",
            error_report=_error_report(),
        )


def test_r1_the_representation_vocabulary_stays_classical() -> None:
    """A tripwire, not an enforcement: R1 is reported UNENFORCED and this pins the vocabulary."""
    assert {kind.value for kind in RepresentationKind} == {
        "explicit-state",
        "empirical-ensemble",
        "local-gaussian",
        "gaussian-mixture",
        "tensor-train",
        "latent-density",
        "surrogate",
        "synthetic-image",
        "scalar-projection",
    }
    assert audit_item("R1").status is AuditStatus.UNENFORCED


def test_a13_the_interface_its_verdict_names_now_exists_and_is_exported() -> None:
    """The A13 tripwire, inverted 2026-07-29 when the two types landed.

    Until then this asserted the OPPOSITE — that neither type existed — which was the register's
    sharpest finding: a verdict quantified over "every tool" that no tool could satisfy, because
    the interface it named had never been written.  Both now exist and are reachable from the
    package root, so the row is enforced and the assertion flips with it.  Keeping the test (rather
    than deleting it) is the point: if either export is ever removed, A13 silently reverts to
    naming a type nobody has, and this fails instead.
    """
    package = __import__("aleph.virtual_cell", fromlist=["__all__"])
    for required in ("CellStatePosterior", "ObservationOperator"):
        assert hasattr(package, required), (
            f"{required} is gone from the package root; A13's verdict once again names an "
            f"interface that cannot be satisfied, and the register must say so"
        )
        assert required in package.__all__, f"{required} exists but is not exported"
    assert audit_item("A13").status is AuditStatus.MECHANICALLY_ENFORCED


def test_a13_still_declares_the_half_of_its_verdict_that_is_not_checked() -> None:
    """Enforced is not the same as fully enforced, and the row must not let that blur.

    The verdict says EVERY tool takes or returns one of the two types.  What is enforced is that
    the types exist and that four specific refusals fire; what is NOT enforced is the universal
    quantifier — nothing requires a new module to touch either type.  An enforced row that hid
    that would be exactly the over-claim this register exists to catch.
    """
    item = audit_item("A13")
    assert item.residual_gap, "an enforced row must state what it still does not cover"
    assert "EVERY tool" in item.residual_gap or "every tool" in item.residual_gap.lower()
    assert "__all__" in item.residual_gap, "the residual gap must name how it would be closed"


# --------------------------------------------------------------------------------------------
# 10. A1 — the reproduced measurement behind the EVIDENCED status
# --------------------------------------------------------------------------------------------


def test_a1_merging_two_uncertainty_sources_into_one_gaussian_misreports_a_nonlinear_moment() -> (
    None
):
    """Closed form, so the number is checked against truth rather than against a sample."""
    sigma_thermal, sigma_biological = 0.1, 1.0
    mixture_fourth_moment = 0.5 * 3.0 * sigma_thermal**4 + 0.5 * 3.0 * sigma_biological**4
    total_variance = 0.5 * sigma_thermal**2 + 0.5 * sigma_biological**2
    single_gaussian_fourth_moment = 3.0 * total_variance**2
    ratio = mixture_fourth_moment / single_gaussian_fourth_moment

    assert mixture_fourth_moment == pytest.approx(1.50015, abs=1e-12)
    assert single_gaussian_fourth_moment == pytest.approx(0.765075, abs=1e-12)

    item = audit_item("A1")
    assert item.status is AuditStatus.EVIDENCED
    stored = item.measurements[0]
    assert stored.reproduced_here is True
    assert stored.value == pytest.approx(ratio, rel=1e-12)

    # A moment-matched single Gaussian and the two-source mixture agree on the variance exactly,
    # which is why the second moment cannot detect the merge that the fourth moment does.
    assert total_variance == pytest.approx(0.505, abs=1e-12)


def test_quoted_measurements_are_marked_as_quoted() -> None:
    quoted = [
        measurement
        for item in AUDIT_REGISTER.values()
        for measurement in item.measurements
        if not measurement.reproduced_here
    ]
    assert quoted, "the register quotes lane results and must say so"
    for measurement in quoted:
        assert "§" in measurement.measured_in or "master plan" in measurement.measured_in

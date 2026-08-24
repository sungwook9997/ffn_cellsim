from __future__ import annotations

import hashlib
import math
import subprocess
import sys
from collections.abc import Sequence

import numpy as np
import pytest

from aleph.virtual_cell.reduction import select_array_representation
from aleph.virtual_cell.step_receipt import (
    StepReceipt,
    predicate_vector_sha256,
    state_coverage_sha256,
)
from aleph.virtual_cell.telemetry import (
    AxisRole,
    ProjectionRecord,
    ScalarProjectionSpec,
    ScalarReduction,
    TelemetryArchive,
    TensorAxisSpec,
    build_telemetry_archive,
    compression_envelope,
    gate_compression,
    measure_archive_bytes,
    measure_decode_cost,
    project_accepted_timeline,
    search_axis_order,
)

_COVERAGE = (
    "component_state",
    "connector_binding_soa",
    "topology_free_list",
    "field_state",
    "accepted_clock",
)
_DT = 0.25
_ENTITIES = ("cortex_filament_0007", "cortex_filament_0011", "sf_arc_0002")
_MODES = ("tension_pn", "extension_nm")


def _hash(label: str) -> str:
    return hashlib.sha256(label.encode()).hexdigest()


def _receipt(**overrides: object) -> StepReceipt:
    predicates = overrides.pop("predicate_results", {"inner_converged": True, "balance_ok": True})
    fields: dict[str, object] = {
        "run_id": "run-b",
        "decision_contract_sha256": _hash("contract"),
        "source_run_manifest_sha256": _hash("manifest"),
        "build_stamp_sha256": _hash("build"),
        "attempted_step_index": 0,
        "accepted_step_before": 0,
        "accepted_step_after": 1,
        "accepted_time_before_s": 0.0,
        "accepted_time_after_s": _DT,
        "attempted_dt_s": _DT,
        "accepted": True,
        "predicate_results": predicates,
        "predicate_results_sha256": predicate_vector_sha256(predicates),  # type: ignore[arg-type]
        "state_coverage": _COVERAGE,
        "state_coverage_sha256": state_coverage_sha256(_COVERAGE),
        "state_hash_before": _hash("state-0"),
        "state_hash_after": _hash("state-1"),
        "rng_hash_before": _hash("rng-0"),
        "rng_hash_after": _hash("rng-1"),
        "rejected_state_restored": None,
    }
    fields.update(overrides)
    return StepReceipt(**fields)  # type: ignore[arg-type]


def _chain(accept_flags: Sequence[bool]) -> tuple[StepReceipt, ...]:
    """Build one contiguous attempted-step chain with continuous state and RNG digests."""
    receipts: list[StepReceipt] = []
    accepted_step = 0
    accepted_time = 0.0
    state = 0
    for attempt, accepted in enumerate(accept_flags):
        before_state, before_rng = _hash(f"state-{state}"), _hash(f"rng-{state}")
        next_state = state + 1 if accepted else state
        receipts.append(
            _receipt(
                attempted_step_index=attempt,
                accepted_step_before=accepted_step,
                accepted_step_after=accepted_step + (1 if accepted else 0),
                accepted_time_before_s=accepted_time,
                accepted_time_after_s=accepted_time + (_DT if accepted else 0.0),
                accepted=accepted,
                predicate_results={"inner_converged": accepted, "balance_ok": True},
                state_hash_before=before_state,
                state_hash_after=_hash(f"state-{next_state}"),
                rng_hash_before=before_rng,
                rng_hash_after=_hash(f"rng-{next_state}"),
                rejected_state_restored=None if accepted else True,
            )
        )
        state = next_state
        if accepted:
            accepted_step += 1
            accepted_time += _DT
    return tuple(receipts)


def _entity_axis(labels: Sequence[str] = _ENTITIES) -> TensorAxisSpec:
    return TensorAxisSpec(
        name="entity",
        role=AxisRole.ENTITY,
        semantics="stable filament identity, constant across accepted time",
        labels=tuple(labels),
    )


def _mode_axis() -> TensorAxisSpec:
    return TensorAxisSpec(
        name="mode",
        role=AxisRole.MODE,
        semantics="observable mode read back per entity",
        labels=_MODES,
    )


_MEAN_PROJECTION = ScalarProjectionSpec(
    name="mean_tension",
    version="1.0.0",
    reduction=ScalarReduction.MEAN,
    tolerance=1e-3,
    selection={"mode": 0},
    unit="pN",
)
_RARE_PROJECTION = ScalarProjectionSpec(
    name="rare_local_tension",
    version="1.0.0",
    reduction=ScalarReduction.VALUE,
    tolerance=0.1,
    selection={"accepted_time": 0, "entity": 0, "mode": 0},
    unit="pN",
)


def _archive(
    accept_flags: Sequence[bool] = (True, False, True, True),
    *,
    entity_orderings: Sequence[Sequence[str]] | None = None,
    arrays: Sequence[np.ndarray] | None = None,
) -> TelemetryArchive:
    receipts = _chain(accept_flags)
    if arrays is None:
        arrays = [
            np.arange(len(_ENTITIES) * len(_MODES), dtype=np.float64).reshape(
                len(_ENTITIES), len(_MODES)
            )
            + 10.0 * (attempt + 1)
            for attempt in range(len(receipts))
        ]
    if entity_orderings is None:
        entity_orderings = [list(_ENTITIES) for _ in receipts]
    return build_telemetry_archive(
        archive_id="archive-b-001",
        receipts=receipts,
        per_attempt_arrays=arrays,
        per_attempt_entity_ids=entity_orderings,
        axes=(_entity_axis(), _mode_axis()),
        projections=(_MEAN_PROJECTION, _RARE_PROJECTION),
        value_semantics="per-entity observable readback at an accepted step",
        value_unit="pN",
    )


# --- 1. accepted-time-only indexing ---------------------------------------------------------------


def test_rejected_attempts_are_dropped_and_counted_not_indexed() -> None:
    timeline = project_accepted_timeline(_chain((True, False, False, True)))
    assert timeline.accepted_step_count == 2
    assert timeline.attempted_step_count == 4
    assert timeline.dropped_attempt_count == 2
    assert timeline.dropped_attempt_indices == (1, 2)
    assert timeline.accepted_time_s == (0.25, 0.5)
    assert timeline.accepted_attempt_index == (0, 3)


def test_rejected_step_values_never_enter_the_archive() -> None:
    archive = _archive((True, False, True))
    assert archive.shape == (2, 3, 2)
    assert archive.timeline.dropped_attempt_count == 1
    # attempt 1 carried the +20.0 offset; only attempts 0 and 2 may occupy a time index.
    assert float(archive.values[0, 0, 0]) == 10.0
    assert float(archive.values[1, 0, 0]) == 30.0
    assert 20.0 not in set(np.ravel(archive.values).tolist())


def test_time_axis_carries_accepted_seconds_and_is_axis_zero() -> None:
    archive = _archive((True, False, True, True))
    assert archive.axes[0].role is AxisRole.ACCEPTED_TIME
    assert archive.axes[0].unit == "s"
    assert archive.axes[0].coordinates == (0.25, 0.5, 0.75)
    assert archive.axis_names == ("accepted_time", "entity", "mode")


def test_all_rejected_chain_has_no_time_index_at_all() -> None:
    with pytest.raises(ValueError, match="accepted no step"):
        project_accepted_timeline(_chain((False, False)))


@pytest.mark.parametrize(
    ("second_overrides", "match"),
    [
        ({"attempted_step_index": 7}, "attempted-step index is not contiguous"),
        (
            {
                "accepted_step_before": 3,
                "accepted_step_after": 4,
                "accepted_time_before_s": 0.75,
                "accepted_time_after_s": 1.0,
            },
            "accepted-step clock is non-monotonic",
        ),
        (
            {"accepted_step_before": 5, "accepted_step_after": 6},
            "accepted-step clock has a gap",
        ),
        (
            {"accepted_time_before_s": 0.875, "accepted_time_after_s": 1.125},
            "accepted physical time is non-monotonic",
        ),
        (
            {"accepted_time_before_s": 1.25, "accepted_time_after_s": 1.5},
            "accepted physical time has a gap",
        ),
    ],
)
def test_broken_accepted_clock_is_refused_not_silently_reindexed(
    second_overrides: dict[str, object], match: str
) -> None:
    first = _receipt(
        attempted_step_index=5,
        accepted_step_before=3,
        accepted_step_after=4,
        accepted_time_before_s=0.75,
        accepted_time_after_s=1.0,
    )
    second_fields: dict[str, object] = {
        "attempted_step_index": 6,
        "accepted_step_before": 4,
        "accepted_step_after": 5,
        "accepted_time_before_s": 1.0,
        "accepted_time_after_s": 1.25,
        "state_hash_before": _hash("state-1"),
        "state_hash_after": _hash("state-2"),
        "rng_hash_before": _hash("rng-1"),
        "rng_hash_after": _hash("rng-2"),
    }
    second_fields.update(second_overrides)
    with pytest.raises(ValueError, match=match):
        project_accepted_timeline((first, _receipt(**second_fields)))


def test_state_digest_discontinuity_is_still_refused_through_the_shared_validator() -> None:
    first = _receipt()
    second = _receipt(
        attempted_step_index=1,
        accepted_step_before=1,
        accepted_step_after=2,
        accepted_time_before_s=_DT,
        accepted_time_after_s=0.5,
        state_hash_before=_hash("unrelated-state"),
        state_hash_after=_hash("state-2"),
        rng_hash_before=_hash("rng-1"),
        rng_hash_after=_hash("rng-2"),
    )
    with pytest.raises(ValueError, match="state digest is discontinuous"):
        project_accepted_timeline((first, second))


# --- 2. stable entity identity --------------------------------------------------------------------


def test_permuted_entity_order_at_one_time_index_is_detected() -> None:
    orderings = [list(_ENTITIES) for _ in range(3)]
    orderings[2] = [_ENTITIES[1], _ENTITIES[0], _ENTITIES[2]]
    with pytest.raises(ValueError, match="entity ordering changed at accepted time index 1"):
        _archive((True, False, True), entity_orderings=orderings)


def test_unknown_entity_id_is_rejected_rather_than_positionally_absorbed() -> None:
    orderings = [list(_ENTITIES) for _ in range(2)]
    orderings[1] = [_ENTITIES[0], _ENTITIES[1], "ghost_filament"]
    with pytest.raises(ValueError, match="names entities absent"):
        _archive((True, True), entity_orderings=orderings)


def test_entity_count_mismatch_is_rejected() -> None:
    orderings = [list(_ENTITIES), list(_ENTITIES[:2])]
    with pytest.raises(ValueError, match="supplies 2 entities"):
        _archive((True, True), entity_orderings=orderings)


def test_reordering_at_a_rejected_attempt_does_not_enter_the_archive() -> None:
    orderings = [list(_ENTITIES) for _ in range(3)]
    orderings[1] = [_ENTITIES[2], _ENTITIES[0], _ENTITIES[1]]
    archive = _archive((True, False, True), entity_orderings=orderings)
    assert archive.shape == (2, 3, 2)


def test_archive_requires_exactly_one_entity_axis() -> None:
    with pytest.raises(ValueError, match="exactly one stable-entity-id axis"):
        build_telemetry_archive(
            archive_id="archive-no-entity",
            receipts=_chain((True,)),
            per_attempt_arrays=[np.ones((3, 2))],
            per_attempt_entity_ids=[list(_ENTITIES)],
            axes=(
                TensorAxisSpec(
                    name="component",
                    role=AxisRole.COMPONENT,
                    semantics="owning component",
                    labels=("cortex", "sf_arc", "lamellipodium"),
                ),
                _mode_axis(),
            ),
            projections=(_MEAN_PROJECTION,),
            value_semantics="x",
            value_unit="pN",
        )


def test_caller_may_not_supply_its_own_accepted_time_axis() -> None:
    supplied_time = TensorAxisSpec(
        name="t",
        role=AxisRole.ACCEPTED_TIME,
        semantics="hand-rolled time",
        coordinates=(0.0, 1.0),
        unit="s",
    )
    with pytest.raises(ValueError, match="derived from the receipt chain"):
        build_telemetry_archive(
            archive_id="archive-own-time",
            receipts=_chain((True,)),
            per_attempt_arrays=[np.ones((2, 3))],
            per_attempt_entity_ids=[list(_ENTITIES)],
            axes=(supplied_time, _entity_axis()),
            projections=(_MEAN_PROJECTION,),
            value_semantics="x",
            value_unit="pN",
        )


def test_per_step_array_shape_must_match_the_declared_axes() -> None:
    with pytest.raises(ValueError, match="does not match the declared axes"):
        _archive((True, True), arrays=[np.ones((3, 2)), np.ones((3, 3))])


def test_axis_roles_cover_the_required_meanings() -> None:
    values = {role.value for role in AxisRole}
    assert {
        "accepted-physical-time",
        "stable-entity-id",
        "component",
        "connector",
        "mode",
    } <= values


# --- 3. versioned scalar projections --------------------------------------------------------------


def test_recorded_scalar_is_reproduced_bit_for_bit_from_the_raw_archive() -> None:
    archive = _archive()
    recorded = {record.name: record.value for record in archive.recorded_scalars}
    assert archive.reproject("mean_tension") == recorded["mean_tension"]
    assert archive.reproject("rare_local_tension") == recorded["rare_local_tension"]
    assert recorded["rare_local_tension"] == 10.0


def test_projection_hash_changes_when_the_definition_changes() -> None:
    baseline = _MEAN_PROJECTION.definition_sha256
    assert (
        ScalarProjectionSpec(
            name="mean_tension",
            version="1.0.0",
            reduction=ScalarReduction.MAX,
            tolerance=1e-3,
            selection={"mode": 0},
            unit="pN",
        ).definition_sha256
        != baseline
    )
    assert (
        ScalarProjectionSpec(
            name="mean_tension",
            version="1.0.1",
            reduction=ScalarReduction.MEAN,
            tolerance=1e-3,
            selection={"mode": 0},
            unit="pN",
        ).definition_sha256
        != baseline
    )
    assert (
        ScalarProjectionSpec(
            name="mean_tension",
            version="1.0.0",
            reduction=ScalarReduction.MEAN,
            tolerance=0.5,
            selection={"mode": 0},
            unit="pN",
        ).definition_sha256
        != baseline
    )
    assert (
        ScalarProjectionSpec(
            name="mean_tension",
            version="1.0.0",
            reduction=ScalarReduction.MEAN,
            tolerance=1e-3,
            selection={"mode": 1},
            unit="pN",
        ).definition_sha256
        != baseline
    )


def test_a_recorded_scalar_from_another_definition_cannot_be_attached() -> None:
    archive = _archive()
    forged = ProjectionRecord(
        name="mean_tension",
        version="1.0.0",
        definition_sha256=_hash("some-other-projection"),
        value=archive.reproject("mean_tension"),
        unit="pN",
    )
    with pytest.raises(ValueError, match="different projection definition"):
        TelemetryArchive(
            archive_id=archive.archive_id,
            timeline=archive.timeline,
            axes=archive.axes,
            values=archive.values,
            projections=archive.projections,
            recorded_scalars=(forged, archive.recorded_scalars[1]),
            value_semantics=archive.value_semantics,
            value_unit=archive.value_unit,
        )


def test_a_recorded_scalar_that_no_longer_reprojects_is_rejected() -> None:
    archive = _archive()
    tampered = ProjectionRecord(
        name="mean_tension",
        version="1.0.0",
        definition_sha256=_MEAN_PROJECTION.definition_sha256,
        value=archive.reproject("mean_tension") + 1.0,
        unit="pN",
    )
    with pytest.raises(ValueError, match="no longer reproduces its recorded scalar"):
        TelemetryArchive(
            archive_id=archive.archive_id,
            timeline=archive.timeline,
            axes=archive.axes,
            values=archive.values,
            projections=archive.projections,
            recorded_scalars=(tampered, archive.recorded_scalars[1]),
            value_semantics=archive.value_semantics,
            value_unit=archive.value_unit,
        )


def test_projection_refuses_an_axis_it_does_not_know() -> None:
    spec = ScalarProjectionSpec(
        name="bad",
        version="1",
        reduction=ScalarReduction.MEAN,
        tolerance=0.1,
        selection={"connector": 0},
    )
    with pytest.raises(ValueError, match="unknown axes"):
        spec.evaluate(np.ones((2, 2)), ("accepted_time", "entity"))


def test_value_projection_must_select_a_single_element() -> None:
    spec = ScalarProjectionSpec(
        name="bad_value",
        version="1",
        reduction=ScalarReduction.VALUE,
        tolerance=0.1,
        selection={"accepted_time": 0},
    )
    with pytest.raises(ValueError, match="single element"):
        spec.evaluate(np.ones((2, 2)), ("accepted_time", "entity"))


# --- 4. compression gated per projection ----------------------------------------------------------


def test_low_rank_archive_compression_is_accepted() -> None:
    tensor = np.einsum(
        "i,j,k->ijk",
        np.arange(1.0, 13.0),
        np.arange(1.0, 8.0),
        np.arange(1.0, 6.0),
    )
    decision = gate_compression(
        tensor,
        axis_names=("accepted_time", "entity", "mode"),
        projections=(
            ScalarProjectionSpec(
                name="mean",
                version="1",
                reduction=ScalarReduction.MEAN,
                tolerance=1e-9,
            ),
            ScalarProjectionSpec(
                name="corner",
                version="1",
                reduction=ScalarReduction.VALUE,
                tolerance=1e-9,
                selection={"accepted_time": 0, "entity": 0, "mode": 0},
            ),
        ),
        global_tolerance=1e-10,
    )
    assert decision.accepted
    assert decision.used_kind == "tensor-train"
    assert decision.parameter_count_ratio > 1.0
    assert decision.expansion_verdict is None


def test_tiny_global_error_cannot_approve_a_destroyed_rare_projection() -> None:
    array = np.ones((100, 100))
    array[0, 0] = 0.05
    projections = (
        ScalarProjectionSpec(
            name="global_mean",
            version="1",
            reduction=ScalarReduction.MEAN,
            tolerance=1e-3,
        ),
        ScalarProjectionSpec(
            name="rare_local",
            version="1",
            reduction=ScalarReduction.VALUE,
            tolerance=0.1,
            selection={"accepted_time": 0, "entity": 0},
        ),
    )
    decision = gate_compression(
        array,
        axis_names=("accepted_time", "entity"),
        projections=projections,
        global_tolerance=1e-2,
    )
    assert decision.candidate_kind == "svd"
    assert not decision.accepted
    assert decision.used_kind == "dense"
    assert decision.rejection_reason == "declared-projection-error-exceeds-tolerance"
    assert decision.failed_projections == ("rare_local",)
    assert decision.measured_global_relative_error < 1e-2
    assert decision.projection_errors["global_mean"] < 1e-3
    assert decision.projection_errors["rare_local"] > 1.0

    # Cross-check: the shared selector, given the same declared projections, also falls back.
    library = select_array_representation(
        array,
        relative_tolerance=1e-2,
        projections={
            spec.name: spec.as_callable(("accepted_time", "entity")) for spec in projections
        },
        projection_tolerances={spec.name: spec.tolerance for spec in projections},
    )
    assert library.kind == "dense"
    assert library.fallback_reason == "declared-projection-error-exceeds-tolerance"


def test_projection_failure_routes_to_the_ensemble_expansion() -> None:
    array = np.ones((100, 100))
    array[0, 0] = 0.05
    decision = gate_compression(
        array,
        axis_names=("accepted_time", "entity"),
        projections=(
            ScalarProjectionSpec(
                name="rare_local",
                version="1",
                reduction=ScalarReduction.VALUE,
                tolerance=0.1,
                selection={"accepted_time": 0, "entity": 0},
            ),
        ),
        global_tolerance=1e-2,
    )
    verdict = decision.expansion_verdict
    assert verdict is not None
    assert verdict.wall_id == "projection-failure"
    assert verdict.expansion == "ensemble / regime partition / latent coordinate"
    assert verdict.measured["worst_projection_error"] > 1.0


def test_incompressible_archive_routes_to_the_rank_expansion() -> None:
    array = np.random.default_rng(11).normal(size=(6, 6, 6))
    decision = gate_compression(
        array,
        axis_names=("accepted_time", "entity", "mode"),
        projections=(
            ScalarProjectionSpec(
                name="mean", version="1", reduction=ScalarReduction.MEAN, tolerance=1e-6
            ),
        ),
        global_tolerance=0.0,
    )
    assert not decision.accepted
    verdict = decision.expansion_verdict
    assert verdict is not None
    assert verdict.wall_id == "rank-failure"
    assert verdict.expansion == "ensemble / regime partition / latent coordinate"


def test_archive_compress_uses_its_own_declared_projections() -> None:
    archive = _archive()
    decision = archive.compress(global_tolerance=1e-10)
    assert set(decision.projection_errors) == {"mean_tension", "rare_local_tension"}
    assert decision.projection_tolerances["rare_local_tension"] == 0.1


# --- 5. axis-order search with recorded coverage --------------------------------------------------


def _order_sensitive_tensor() -> np.ndarray:
    rng = np.random.default_rng(5)
    left = rng.normal(size=(4, 4))
    right = rng.normal(size=(4, 4))
    return np.einsum("il,jk->ijkl", left, right)


def test_exhaustive_axis_order_search_reports_full_coverage_and_finds_the_cheap_order() -> None:
    report = search_axis_order(_order_sensitive_tensor(), relative_tolerance=1e-12, max_orders=24)
    assert report.strategy == "exhaustive"
    assert report.total_orderings == 24
    assert report.evaluated_count == 24
    assert report.unsearched_count == 0
    assert report.unsearched_fraction == 0.0
    identity_storage = report.storage_by_order[(0, 1, 2, 3)]
    assert report.best_storage_size < identity_storage
    assert report.any_order_met_tolerance
    # (i,l | j,k) is the exactly rank-1 split of an outer product of two matrices.
    assert report.best_order in {(0, 3, 1, 2), (0, 3, 2, 1), (1, 2, 0, 3), (2, 1, 0, 3)}


def test_bounded_axis_order_search_reports_the_region_it_did_not_search() -> None:
    tensor = np.random.default_rng(2).normal(size=(2, 3, 2, 3, 2))
    report = search_axis_order(tensor, relative_tolerance=0.05, max_orders=10)
    assert report.strategy == "heuristic-subset"
    assert report.total_orderings == math.factorial(5) == 120
    assert report.evaluated_count == 10
    assert report.unsearched_count == 110
    assert report.unsearched_fraction == pytest.approx(110 / 120)
    assert len(report.evaluated_orderings) == len(set(report.evaluated_orderings))
    assert report.to_dict()["unsearched_fraction"] == report.unsearched_fraction


def test_axis_order_search_covers_every_ordering_when_the_space_is_small() -> None:
    report = search_axis_order(
        np.random.default_rng(4).normal(size=(3, 4, 2)),
        relative_tolerance=0.1,
        max_orders=24,
    )
    assert report.total_orderings == 6
    assert report.evaluated_count == 6
    assert report.unsearched_fraction == 0.0


def test_axis_order_search_result_feeds_the_compression_gate() -> None:
    tensor = _order_sensitive_tensor()
    report = search_axis_order(tensor, relative_tolerance=1e-12, max_orders=24)
    decision = gate_compression(
        tensor,
        axis_names=("accepted_time", "entity", "mode", "connector"),
        projections=(
            ScalarProjectionSpec(
                name="mean", version="1", reduction=ScalarReduction.MEAN, tolerance=1e-8
            ),
        ),
        global_tolerance=1e-10,
        axis_order=report.best_order,
    )
    assert decision.accepted
    assert decision.axis_order == report.best_order
    assert decision.stored_parameter_count == report.best_storage_size


# --- 6. real bytes and real decode cost, measured separately --------------------------------------


def test_byte_ratio_is_measured_on_real_serialized_bytes_including_metadata(tmp_path) -> None:
    tensor = _order_sensitive_tensor()
    report = search_axis_order(tensor, relative_tolerance=1e-12, max_orders=24)
    decision = gate_compression(
        tensor,
        axis_names=("accepted_time", "entity", "mode", "connector"),
        global_tolerance=1e-10,
        axis_order=report.best_order,
    )
    measurement = measure_archive_bytes(
        tensor,
        decision.reduction,
        metadata={"schema_id": "test", "axis_order": list(report.best_order)},
        probe_path=tmp_path / "archive.npz",
    )
    assert measurement.dense_total_bytes > 0
    assert measurement.reduced_total_bytes > 0
    assert measurement.metadata_bytes > 0
    assert measurement.byte_ratio == pytest.approx(
        measurement.dense_total_bytes / measurement.reduced_total_bytes
    )
    assert measurement.parameter_count_ratio == pytest.approx(256 / 64)
    # The two ratios are different numbers and must never be quoted interchangeably.
    assert measurement.byte_ratio != pytest.approx(measurement.parameter_count_ratio)
    assert (tmp_path / "archive.npz").stat().st_size == measurement.reduced_total_bytes


def test_metadata_overhead_can_make_a_compressed_archive_larger_in_bytes() -> None:
    matrix = np.multiply.outer(np.arange(1.0, 5.0), np.arange(1.0, 5.0))
    decision = gate_compression(
        matrix,
        axis_names=("accepted_time", "entity"),
        global_tolerance=1e-12,
    )
    assert decision.accepted
    assert decision.parameter_count_ratio > 1.0
    measurement = measure_archive_bytes(matrix, decision.reduction)
    # A parameter-count win does not imply a byte win: three NPZ members plus metadata cost more
    # than the four float64 slots saved.  Reporting only the parameter ratio would be a false claim.
    assert measurement.reduced_form_is_larger
    assert measurement.byte_ratio < 1.0 < measurement.parameter_count_ratio


def test_a_dense_decision_measures_exactly_one_byte_ratio() -> None:
    array = np.random.default_rng(21).normal(size=(5, 5))
    decision = gate_compression(
        array,
        axis_names=("accepted_time", "entity"),
        global_tolerance=0.0,
    )
    assert not decision.accepted
    measurement = measure_archive_bytes(array, decision.reduction, metadata={"schema_id": "test"})
    # A dense fallback stores the same array twice; the two forms must not differ by a metadata
    # string, or a fallback would report a few bytes of phantom compression.
    assert measurement.dense_total_bytes == measurement.reduced_total_bytes
    assert measurement.byte_ratio == 1.0


def test_decode_cost_is_reported_as_a_distribution_not_a_single_timing() -> None:
    tensor = _order_sensitive_tensor()
    report = search_axis_order(tensor, relative_tolerance=1e-12, max_orders=24)
    decision = gate_compression(
        tensor,
        axis_names=("accepted_time", "entity", "mode", "connector"),
        global_tolerance=1e-10,
        axis_order=report.best_order,
    )
    cost = measure_decode_cost(tensor, decision.reduction, repeats=5)
    assert cost.repeats == 5
    assert len(cost.dense_samples_s) == 5
    assert len(cost.reduced_samples_s) == 5
    assert set(cost.dense_summary_s) == {"min_s", "median_s", "max_s"}
    assert cost.dense_summary_s["min_s"] <= cost.dense_summary_s["median_s"]
    assert cost.dense_summary_s["median_s"] <= cost.dense_summary_s["max_s"]
    assert cost.reduced_summary_s["min_s"] > 0.0
    assert cost.reconstruction_max_abs_difference < 1e-9
    assert "byte_ratio" not in cost.to_dict()


def test_decode_reconstructs_from_bytes_alone() -> None:
    archive = _archive()
    decision = archive.compress(global_tolerance=1e-10)
    cost = measure_decode_cost(
        archive.values,
        decision.reduction,
        metadata=archive.to_metadata(),
        repeats=3,
    )
    assert cost.reconstruction_max_abs_difference < 1e-9


# --- composition with the sidecar and provenance contracts ----------------------------------------


def test_dense_sidecar_round_trips_through_the_shared_array_contract() -> None:
    archive = _archive()
    sidecar = archive.dense_sidecar(uri="memory://archive-b-001.npy")
    decoded = sidecar.verify_npy_array(archive.dense_payload())
    np.testing.assert_array_equal(decoded, archive.values)
    assert sidecar.shape == archive.shape
    assert [axis.name for axis in sidecar.axes] == list(archive.axis_names)
    assert sidecar.axes[0].coordinates_sha256 is not None


def test_archive_hash_binds_payload_and_metadata() -> None:
    first = _archive()
    second = _archive()
    assert first.archive_sha256 == second.archive_sha256
    assert first.to_metadata()["timeline"]["dropped_attempt_count"] == 1


def test_compression_envelope_refuses_to_launder_a_rejected_reduction() -> None:
    archive = _archive()
    decision = gate_compression(
        archive.values,
        axis_names=archive.axis_names,
        projections=archive.projections,
        global_tolerance=0.0,
    )
    assert not decision.accepted
    with pytest.raises(ValueError, match="no reduced artifact"):
        compression_envelope(
            archive,
            decision,
            artifact_id="artifact-1",
            cell_state_manifest_hash=_hash("manifest"),
            source_artifact_ids=("native-run-1",),
        )


def test_compression_envelope_refuses_to_call_a_matrix_svd_a_tensor_train() -> None:
    matrix = np.multiply.outer(np.arange(1.0, 41.0), np.arange(1.0, 21.0))
    decision = gate_compression(
        matrix,
        axis_names=("accepted_time", "entity"),
        global_tolerance=1e-12,
    )
    assert decision.used_kind == "svd"
    archive = _archive()
    with pytest.raises(ValueError, match="no matrix-SVD member"):
        compression_envelope(
            archive,
            decision,
            artifact_id="artifact-2",
            cell_state_manifest_hash=_hash("manifest"),
            source_artifact_ids=("native-run-1",),
        )


def test_accepted_tensor_train_gets_a_reduced_model_envelope() -> None:
    tensor = np.einsum(
        "i,j,k->ijk",
        np.arange(1.0, 13.0),
        np.arange(1.0, 8.0),
        np.arange(1.0, 6.0),
    )
    decision = gate_compression(
        tensor,
        axis_names=("accepted_time", "entity", "mode"),
        global_tolerance=1e-10,
    )
    assert decision.used_kind == "tensor-train"
    archive = _archive()
    envelope = compression_envelope(
        archive,
        decision,
        artifact_id="artifact-3",
        cell_state_manifest_hash=_hash("manifest"),
        source_artifact_ids=("native-run-1",),
    )
    assert envelope.evidence_source.value == "reduced-model"
    assert envelope.representation.kind.value == "tensor-train"
    assert envelope.representation.error_report is not None
    assert envelope.representation.error_report.measured_error == (
        decision.measured_global_relative_error
    )


# --- CPU-only -------------------------------------------------------------------------------------


def test_importing_the_adapter_does_not_import_or_initialise_warp() -> None:
    probe = (
        "import sys;"
        "import aleph.virtual_cell.telemetry as m;"
        "assert 'warp' not in sys.modules, sorted(k for k in sys.modules if 'warp' in k);"
        "print(m.__name__)"
    )
    completed = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "aleph.virtual_cell.telemetry" in completed.stdout

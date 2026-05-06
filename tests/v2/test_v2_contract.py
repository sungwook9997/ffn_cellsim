from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from acs.v2.data_contract import (
    ArtifactLayoutEntry,
    ArtifactKind,
    ImagingDatasetSpec,
    MetricSpec,
    SegmentationProvenance,
    V2DataContract,
    canonical_artifact_layout,
    default_single_cell_contract,
)
from acs.v2.measurement_boundary import MeasurementBoundary
from acs.v2.single_cell import FocalAdhesionState, ProtrusionEvent, SingleCellState

_REPO_ROOT = Path(__file__).resolve().parents[2]
_EXAMPLE_CONTRACT_PATH = (
    _REPO_ROOT / "configs" / "v2_single_cell_data_contract_example.json"
)


def _full_artifact_dataset(**overrides) -> ImagingDatasetSpec:
    base = dict(
        dataset_id="single_cell_demo",
        root_uri="data/imaging/single_cell_demo",
        modality="confocal_live",
        channels=("membrane", "actin"),
        voxel_size_um_xyz=(0.25, 0.25, 1.0),
        frame_interval_s=60.0,
        n_timepoints=12,
        calibration_timepoints=(0, 1, 2, 3),
        validation_timepoints=(8, 9, 10, 11),
        available_artifacts=(
            ArtifactKind.SEGMENTATION_MASK,
            ArtifactKind.BOUNDARY_CONTOURS,
            ArtifactKind.EVENT_ANNOTATIONS,
            ArtifactKind.TRACKING_TABLE,
        ),
        segmentation_provenance=SegmentationProvenance(
            method="cellpose_then_manual_qc",
            software="Cellpose",
            version="2.0",
            postprocessing=("remove_small_objects", "manual_boundary_qc"),
            reviewer="test_fixture",
        ),
    )
    base.update(overrides)
    return ImagingDatasetSpec(**base)


def _contract_from_json(path: Path) -> V2DataContract:
    payload = json.loads(path.read_text(encoding="utf-8"))
    dataset_payload = dict(payload["dataset"])
    dataset_payload["channels"] = tuple(dataset_payload["channels"])
    dataset_payload["voxel_size_um_xyz"] = tuple(dataset_payload["voxel_size_um_xyz"])
    dataset_payload["calibration_timepoints"] = tuple(
        dataset_payload.get("calibration_timepoints", ())
    )
    dataset_payload["validation_timepoints"] = tuple(
        dataset_payload.get("validation_timepoints", ())
    )
    dataset_payload["available_artifacts"] = tuple(
        ArtifactKind(value) for value in dataset_payload.get("available_artifacts", ())
    )
    dataset_payload["available_csv_columns"] = tuple(
        dataset_payload.get("available_csv_columns", ())
    )
    dataset_payload["artifact_layout"] = tuple(
        ArtifactLayoutEntry(
            artifact=ArtifactKind(entry["artifact"]),
            relative_path=entry["relative_path"],
            format_hint=entry["format_hint"],
            required=entry.get("required", True),
        )
        for entry in dataset_payload.get("artifact_layout", ())
    )
    raw_provenance = dataset_payload.get("segmentation_provenance")
    if raw_provenance is not None:
        dataset_payload["segmentation_provenance"] = SegmentationProvenance(
            method=raw_provenance["method"],
            software=raw_provenance["software"],
            version=raw_provenance["version"],
            source_artifact=ArtifactKind(
                raw_provenance.get("source_artifact", ArtifactKind.RAW_FRAMES.value)
            ),
            coordinate_convention=raw_provenance.get(
                "coordinate_convention",
                "image_um_y_down",
            ),
            postprocessing=tuple(raw_provenance.get("postprocessing", ())),
            reviewer=raw_provenance.get("reviewer", ""),
        )
    dataset = ImagingDatasetSpec(**dataset_payload)

    metrics = []
    for raw_metric in payload["metrics"]:
        metric_payload = dict(raw_metric)
        metric_payload["required_artifacts"] = tuple(
            ArtifactKind(value)
            for value in metric_payload.get("required_artifacts", ())
        )
        metric_payload["csv_columns_required"] = tuple(
            metric_payload.get("csv_columns_required", ())
        )
        metric_payload["channels_required"] = tuple(
            metric_payload.get("channels_required", ())
        )
        metrics.append(MetricSpec(**metric_payload))
    return V2DataContract(dataset=dataset, metrics=tuple(metrics))


def test_artifact_kind_has_nine_members_with_ecm_field():
    members = list(ArtifactKind)
    assert len(members) == 9
    assert ArtifactKind.ECM_FIELD in members
    expected = {
        "csv_table",
        "raw_frames",
        "segmentation_mask",
        "boundary_contours",
        "tracking_table",
        "marker_channel",
        "tfm_field",
        "event_annotations",
        "ecm_field",
    }
    assert {a.value for a in members} == expected


def test_example_single_cell_data_contract_json_validates():
    """The roadmap's example JSON data contract stays loadable into the
    canonical V2 dataclasses and validates without relying on PI data."""

    contract = _contract_from_json(_EXAMPLE_CONTRACT_PATH)
    contract.validate()

    assert contract.dataset.dataset_id == "v2_single_cell_col1_demo"
    assert contract.dataset.voxel_size_um_xyz == (0.25, 0.25, 1.0)
    assert contract.dataset.frame_interval_s == pytest.approx(60.0)
    assert set(contract.dataset.calibration_timepoints).isdisjoint(
        contract.dataset.validation_timepoints
    )
    assert ArtifactKind.CSV_TABLE in contract.dataset.available_artifacts
    assert contract.dataset.artifact_uri(ArtifactKind.SEGMENTATION_MASK) == (
        "data/imaging/v2_single_cell_col1_demo/"
        "segmentation/masks/t{timepoint:04d}.tif"
    )
    assert contract.dataset.segmentation_provenance is not None
    assert contract.dataset.segmentation_provenance.software == "Cellpose"

    keys = {(m.name, m.measurement_modality) for m in contract.metrics}
    assert ("projected_area", "top_down_segmentation") in keys
    assert ("projected_area", "csv_top_down") in keys
    assert len(contract.calibration_metrics) == 1
    assert len(contract.validation_metrics) == 2


def test_default_single_cell_contract_validates_with_full_artifacts():
    dataset = _full_artifact_dataset()
    contract = default_single_cell_contract(dataset)
    contract.validate()
    assert len(contract.calibration_metrics) == 1
    assert len(contract.validation_metrics) == 2


def test_default_contract_auto_filters_to_available_artifacts():
    dataset = _full_artifact_dataset(
        available_artifacts=(ArtifactKind.SEGMENTATION_MASK,),
    )
    contract = default_single_cell_contract(dataset)
    contract.validate()
    metric_names = tuple(m.name for m in contract.metrics)
    assert metric_names == ("projected_area",)


def test_default_contract_with_no_artifacts_yields_empty_metrics_and_fails_validate():
    dataset = _full_artifact_dataset(available_artifacts=())
    contract = default_single_cell_contract(dataset)
    assert contract.metrics == ()
    with pytest.raises(ValueError, match="at least one metric"):
        contract.validate()


def test_contract_rejects_calibration_validation_overlap():
    dataset = _full_artifact_dataset(
        n_timepoints=3,
        calibration_timepoints=(1,),
        validation_timepoints=(1,),
        available_artifacts=(ArtifactKind.SEGMENTATION_MASK,),
    )
    contract = default_single_cell_contract(dataset)
    with pytest.raises(ValueError, match="overlap"):
        contract.validate()


def test_metric_required_artifacts_fail_when_dataset_missing_artifact():
    dataset = _full_artifact_dataset(
        available_artifacts=(ArtifactKind.SEGMENTATION_MASK,),
    )
    metric = MetricSpec(
        name="boundary_roughness",
        target_object="single_cell",
        measurement_modality="cell_outline",
        unit="dimensionless",
        role="validation",
        required_artifacts=(ArtifactKind.BOUNDARY_CONTOURS,),
    )
    with pytest.raises(ValueError, match="boundary_contours"):
        metric.validate(dataset)


def test_metric_required_artifacts_pass_subset():
    dataset = _full_artifact_dataset(
        available_artifacts=(
            ArtifactKind.SEGMENTATION_MASK,
            ArtifactKind.BOUNDARY_CONTOURS,
        ),
    )
    metric = MetricSpec(
        name="boundary_roughness",
        target_object="single_cell",
        measurement_modality="cell_outline",
        unit="dimensionless",
        role="validation",
        required_artifacts=(ArtifactKind.BOUNDARY_CONTOURS,),
    )
    metric.validate(dataset)


def test_metric_csv_columns_require_csv_table_artifact():
    dataset = _full_artifact_dataset(
        available_artifacts=(ArtifactKind.SEGMENTATION_MASK,),
    )
    metric = MetricSpec(
        name="experimental_area",
        target_object="single_cell",
        measurement_modality="csv_top_down",
        unit="um2",
        role="calibration",
        required_artifacts=(ArtifactKind.SEGMENTATION_MASK,),
        csv_columns_required=("Area_um2",),
    )
    with pytest.raises(ValueError, match="CSV_TABLE"):
        metric.validate(dataset)


def test_metric_csv_columns_pass_when_csv_table_and_columns_available():
    dataset = _full_artifact_dataset(
        available_artifacts=(
            ArtifactKind.SEGMENTATION_MASK,
            ArtifactKind.CSV_TABLE,
        ),
        available_csv_columns=("Area_um2", "Time_s"),
    )
    metric = MetricSpec(
        name="experimental_area",
        target_object="single_cell",
        measurement_modality="csv_top_down",
        unit="um2",
        role="calibration",
        required_artifacts=(
            ArtifactKind.SEGMENTATION_MASK,
            ArtifactKind.CSV_TABLE,
        ),
        csv_columns_required=("Area_um2",),
    )
    metric.validate(dataset)


def test_metric_csv_columns_fail_when_column_missing():
    dataset = _full_artifact_dataset(
        available_artifacts=(
            ArtifactKind.SEGMENTATION_MASK,
            ArtifactKind.CSV_TABLE,
        ),
        available_csv_columns=("Time_s",),
    )
    metric = MetricSpec(
        name="experimental_area",
        target_object="single_cell",
        measurement_modality="csv_top_down",
        unit="um2",
        role="calibration",
        required_artifacts=(
            ArtifactKind.SEGMENTATION_MASK,
            ArtifactKind.CSV_TABLE,
        ),
        csv_columns_required=("Area_um2",),
    )
    with pytest.raises(ValueError, match="Area_um2"):
        metric.validate(dataset)


def test_metric_channels_required_pass_and_fail():
    dataset = _full_artifact_dataset()
    ok_metric = MetricSpec(
        name="actin_intensity",
        target_object="single_cell",
        measurement_modality="channel_intensity",
        unit="au",
        role="diagnostic",
        required_artifacts=(ArtifactKind.RAW_FRAMES,),
        channels_required=("actin",),
    )
    dataset_with_raw = _full_artifact_dataset(
        available_artifacts=(*dataset.available_artifacts, ArtifactKind.RAW_FRAMES),
    )
    ok_metric.validate(dataset_with_raw)

    bad_metric = MetricSpec(
        name="dapi_intensity",
        target_object="single_cell",
        measurement_modality="channel_intensity",
        unit="au",
        role="diagnostic",
        required_artifacts=(ArtifactKind.RAW_FRAMES,),
        channels_required=("dapi",),
    )
    with pytest.raises(ValueError, match="dapi"):
        bad_metric.validate(dataset_with_raw)


def test_duplicate_metric_key_name_modality():
    dataset = _full_artifact_dataset()
    metric = MetricSpec(
        name="projected_area",
        target_object="single_cell",
        measurement_modality="top_down_segmentation",
        unit="um2",
        role="calibration",
        required_artifacts=(ArtifactKind.SEGMENTATION_MASK,),
    )
    contract = V2DataContract(dataset=dataset, metrics=(metric, metric))
    with pytest.raises(ValueError, match="duplicate metric key"):
        contract.validate()


def test_distinct_metrics_same_name_different_modality_allowed():
    dataset = _full_artifact_dataset(
        available_artifacts=(
            ArtifactKind.SEGMENTATION_MASK,
            ArtifactKind.CSV_TABLE,
        ),
        available_csv_columns=("Area_um2",),
    )
    seg_metric = MetricSpec(
        name="projected_area",
        target_object="single_cell",
        measurement_modality="top_down_segmentation",
        unit="um2",
        role="calibration",
        required_artifacts=(ArtifactKind.SEGMENTATION_MASK,),
    )
    csv_metric = MetricSpec(
        name="projected_area",
        target_object="single_cell",
        measurement_modality="csv_top_down",
        unit="um2",
        role="calibration",
        required_artifacts=(ArtifactKind.CSV_TABLE,),
        csv_columns_required=("Area_um2",),
    )
    contract = V2DataContract(dataset=dataset, metrics=(seg_metric, csv_metric))
    contract.validate()
    assert len(contract.calibration_metrics) == 2


def test_dataset_csv_columns_without_csv_artifact_rejected():
    with pytest.raises(ValueError, match="CSV_TABLE"):
        ImagingDatasetSpec(
            dataset_id="bad",
            root_uri="data/bad",
            modality="confocal_live",
            channels=("membrane",),
            voxel_size_um_xyz=(0.25, 0.25, 1.0),
            frame_interval_s=60.0,
            n_timepoints=3,
            available_csv_columns=("Time_s",),
        ).validate()


def test_dataset_duplicate_artifact_rejected():
    with pytest.raises(ValueError, match="duplicate artifact"):
        ImagingDatasetSpec(
            dataset_id="dup",
            root_uri="data/dup",
            modality="confocal_live",
            channels=("membrane",),
            voxel_size_um_xyz=(0.25, 0.25, 1.0),
            frame_interval_s=60.0,
            n_timepoints=3,
            available_artifacts=(
                ArtifactKind.SEGMENTATION_MASK,
                ArtifactKind.SEGMENTATION_MASK,
            ),
        ).validate()


def test_segmentation_artifacts_require_segmentation_provenance():
    dataset = _full_artifact_dataset(segmentation_provenance=None)

    with pytest.raises(ValueError, match="segmentation_provenance"):
        dataset.validate()


def test_segmentation_provenance_validates_required_fields():
    with pytest.raises(ValueError, match="method"):
        SegmentationProvenance(
            method=" ",
            software="Cellpose",
            version="2.0",
        ).validate()

    with pytest.raises(ValueError, match="coordinate_convention"):
        SegmentationProvenance(
            method="cellpose",
            software="Cellpose",
            version="2.0",
            coordinate_convention="pixel_y_down",
        ).validate()

    SegmentationProvenance(
        method="cellpose_then_manual_qc",
        software="Cellpose",
        version="2.0",
        coordinate_convention="image_um_y_down",
        postprocessing=("remove_small_objects",),
    ).validate()


def test_canonical_artifact_layout_covers_advertised_artifacts():
    dataset = _full_artifact_dataset(
        available_artifacts=(
            ArtifactKind.RAW_FRAMES,
            ArtifactKind.SEGMENTATION_MASK,
            ArtifactKind.BOUNDARY_CONTOURS,
            ArtifactKind.TRACKING_TABLE,
            ArtifactKind.EVENT_ANNOTATIONS,
            ArtifactKind.MARKER_CHANNEL,
            ArtifactKind.CSV_TABLE,
        ),
    )
    layout = dataset.canonical_artifact_layout()

    assert tuple(entry.artifact for entry in layout) == dataset.available_artifacts
    assert dataset.artifact_relative_path(ArtifactKind.RAW_FRAMES) == (
        "raw/{channel}/t{timepoint:04d}.tif"
    )
    assert dataset.artifact_uri(ArtifactKind.CSV_TABLE) == (
        "data/imaging/single_cell_demo/tables/measurements.csv"
    )


def test_canonical_artifact_layout_rejects_non_artifact_kind():
    with pytest.raises(ValueError, match="ArtifactKind"):
        canonical_artifact_layout(("raw_frames",))  # type: ignore[arg-type]


def test_declared_artifact_layout_must_cover_available_artifacts_exactly():
    dataset = _full_artifact_dataset(
        available_artifacts=(
            ArtifactKind.SEGMENTATION_MASK,
            ArtifactKind.BOUNDARY_CONTOURS,
        ),
        artifact_layout=(
            ArtifactLayoutEntry(
                artifact=ArtifactKind.SEGMENTATION_MASK,
                relative_path="segmentation/masks/t{timepoint:04d}.tif",
                format_hint="label_mask_tiff",
            ),
        ),
    )
    with pytest.raises(ValueError, match="missing advertised artifacts"):
        dataset.validate()

    dataset = _full_artifact_dataset(
        available_artifacts=(ArtifactKind.SEGMENTATION_MASK,),
        artifact_layout=(
            ArtifactLayoutEntry(
                artifact=ArtifactKind.SEGMENTATION_MASK,
                relative_path="segmentation/masks/t{timepoint:04d}.tif",
                format_hint="label_mask_tiff",
            ),
            ArtifactLayoutEntry(
                artifact=ArtifactKind.SEGMENTATION_MASK,
                relative_path="segmentation/masks/duplicate.tif",
                format_hint="label_mask_tiff",
            ),
        ),
    )
    with pytest.raises(ValueError, match="duplicate layout entry"):
        dataset.validate()


@pytest.mark.parametrize(
    "bad_path, match",
    [
        ("/absolute/raw.tif", "relative"),
        ("segmentation/../raw.tif", "path segments"),
        ("raw\\frames.tif", "POSIX"),
    ],
)
def test_artifact_layout_rejects_unsafe_paths(bad_path: str, match: str):
    entry = ArtifactLayoutEntry(
        artifact=ArtifactKind.RAW_FRAMES,
        relative_path=bad_path,
        format_hint="tiff_stack",
    )
    with pytest.raises(ValueError, match=match):
        entry.validate()


def test_single_cell_state_projected_area_and_nested_validation():
    boundary = np.array([
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 1.0],
        [0.0, 1.0],
    ])
    state = SingleCellState(
        cell_id="cell-1",
        time_s=120.0,
        measurement_boundary=MeasurementBoundary.from_array(
            boundary,
            coordinate_convention="world_um_y_up",
            source_modality="manual_test_outline",
        ),
        height_um=6.0,
        polarity_xy=(1.0, 0.0),
        protrusions=[
            ProtrusionEvent(
                cell_id="cell-1",
                start_time_s=120.0,
                boundary_angle_rad=0.0,
                length_um=3.0,
                event_type="lamellipodium",
                lifetime_s=180.0,
            )
        ],
        adhesions=[
            FocalAdhesionState(
                adhesion_id="fa-1",
                cell_id="cell-1",
                position_um_xy=(1.0, 0.5),
                age_s=60.0,
                maturity=0.4,
                bound_fraction=0.8,
            )
        ],
    )

    assert state.projected_area_um2() == pytest.approx(2.0)

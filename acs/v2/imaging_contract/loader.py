"""YAML loader and validator for the v2 V2-1 imaging input contract.

This package is deliberately PI-260313-specific glue around the generic
contract classes in :mod:`acs.v2.data_contract`. It validates measurement-side
CSV inputs and the materialized scene split manifest; it does not calibrate
or fit model parameters.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from acs.v2.data_contract import (
    ArtifactKind,
    ImagingDatasetSpec,
    MetricSpec,
    V2DataContract,
)


@dataclass(frozen=True, slots=True)
class ImagingInputMetric:
    """One metric declared in the imaging input contract."""

    name: str
    column: str
    unit: str
    role: str
    measurement_modality: str = "csv_top_down"


@dataclass(frozen=True, slots=True)
class ImagingInputContractGroup:
    """One condition group in the 260313 imaging input contract."""

    condition: str
    csv_path: Path
    naming_series_prefix: str
    filename_prefix: str
    scene_id_column: str = "Series"
    frame_column: str = "Frame"
    time_column: str = "Time_min"


@dataclass(frozen=True, slots=True)
class ImagingInputContract:
    """Loaded V2-1 imaging input contract."""

    contract_id: str
    imaging_date: str
    canonical_root: Path
    voxel_size_um_xyz: tuple[float, float, float]
    frame_interval_s: float
    channels: tuple[str, ...]
    groups: tuple[ImagingInputContractGroup, ...]
    metrics: tuple[ImagingInputMetric, ...]
    split_manifest_path: Path
    split_seed: int
    calibration_fraction: float
    raw_imaging_available: bool
    raw_path: Path | None
    yaml_path: Path

    def to_v2_data_contracts(self) -> tuple[V2DataContract, ...]:
        """Return generic V2 data contracts, one per condition group."""

        contracts: list[V2DataContract] = []
        for group in self.groups:
            dataset = ImagingDatasetSpec(
                dataset_id=f"{self.contract_id}_{group.condition}",
                root_uri=str(group.csv_path.parent),
                modality="csv_measurement_table",
                channels=self.channels,
                voxel_size_um_xyz=self.voxel_size_um_xyz,
                frame_interval_s=self.frame_interval_s,
                n_timepoints=_n_timepoints(group.csv_path, frame_column=group.frame_column),
                available_artifacts=(ArtifactKind.CSV_TABLE,),
                available_csv_columns=tuple(_csv_columns(group.csv_path)),
                segmentation_provenance=None,
                notes="PI 260313 CSV measurement-side contract; no parameter fitting.",
            )
            metrics = tuple(
                MetricSpec(
                    name=metric.name,
                    target_object="spheroid_projection",
                    measurement_modality=metric.measurement_modality,
                    unit=metric.unit,
                    role=metric.role,
                    required_artifacts=(ArtifactKind.CSV_TABLE,),
                    csv_columns_required=(metric.column,),
                )
                for metric in self.metrics
            )
            contracts.append(V2DataContract(dataset=dataset, metrics=metrics))
        return tuple(contracts)


@dataclass(frozen=True, slots=True)
class ImagingInputContractValidationResult:
    """Validation result with row/scene counts for auditability."""

    contract_id: str
    group_scene_counts: dict[str, int]
    group_row_counts: dict[str, int]
    warnings: tuple[str, ...] = ()


def load_imaging_contract(path: str | Path) -> ImagingInputContract:
    """Load a V2-1 imaging input YAML contract."""

    yaml_path = Path(path)
    payload = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("imaging contract YAML must contain a mapping")
    base_dir = yaml_path.parent
    contract_id = _require_str(payload, "contract_id")
    imaging_date = _require_str(payload, "imaging_date")
    voxel_size_um_xyz = _float_triplet(payload.get("voxel_size_um_xyz"))
    frame_interval_s = float(payload["frame_interval_s"])
    channels = tuple(_require_list(payload, "channels"))

    groups = tuple(
        ImagingInputContractGroup(
            condition=_require_str(raw_group, "condition"),
            csv_path=_resolve_relative(base_dir, _require_str(raw_group, "csv_path")),
            naming_series_prefix=_require_str(raw_group, "naming_series_prefix"),
            filename_prefix=_require_str(raw_group, "filename_prefix"),
            scene_id_column=str(raw_group.get("scene_id_column", "Series")),
            frame_column=str(raw_group.get("frame_column", "Frame")),
            time_column=str(raw_group.get("time_column", "Time_min")),
        )
        for raw_group in _require_list(payload, "groups")
    )
    metrics = tuple(
        ImagingInputMetric(
            name=_require_str(raw_metric, "name"),
            column=_require_str(raw_metric, "column"),
            unit=_require_str(raw_metric, "unit"),
            role=_require_str(raw_metric, "role"),
            measurement_modality=str(
                raw_metric.get("measurement_modality", "csv_top_down")
            ),
        )
        for raw_metric in _require_list(payload, "metrics")
    )

    split = payload.get("split", {})
    if not isinstance(split, dict):
        raise ValueError("split must be a mapping")
    raw = payload.get("raw_imaging", {})
    if not isinstance(raw, dict):
        raise ValueError("raw_imaging must be a mapping")
    raw_path_value = raw.get("path")

    return ImagingInputContract(
        contract_id=contract_id,
        imaging_date=imaging_date,
        canonical_root=_resolve_relative(
            base_dir, _require_str(payload, "canonical_root")
        ),
        voxel_size_um_xyz=voxel_size_um_xyz,
        frame_interval_s=frame_interval_s,
        channels=channels,
        groups=groups,
        metrics=metrics,
        split_manifest_path=_resolve_relative(
            base_dir, _require_str(split, "manifest")
        ),
        split_seed=int(split["seed"]),
        calibration_fraction=float(split.get("calibration_fraction", 0.7)),
        raw_imaging_available=bool(raw.get("available", False)),
        raw_path=(
            None
            if raw_path_value in {None, ""}
            else _resolve_relative(base_dir, str(raw_path_value))
        ),
        yaml_path=yaml_path,
    )


def validate_imaging_contract(
    contract: ImagingInputContract,
) -> ImagingInputContractValidationResult:
    """Validate the loaded imaging contract and underlying CSV/manifest files."""

    if not contract.contract_id.strip():
        raise ValueError("contract_id must be non-empty")
    if not contract.groups:
        raise ValueError("at least one imaging group is required")
    if not contract.metrics:
        raise ValueError("at least one metric is required")
    if (
        not contract.channels
        or len(set(contract.channels)) != len(contract.channels)
        or any(not str(channel).strip() for channel in contract.channels)
    ):
        raise ValueError("channels must be unique non-empty strings")
    if any(float(v) <= 0.0 for v in contract.voxel_size_um_xyz):
        raise ValueError("voxel_size_um_xyz values must be positive")
    if contract.frame_interval_s <= 0.0:
        raise ValueError("frame_interval_s must be positive")
    if not (0.0 < contract.calibration_fraction < 1.0):
        raise ValueError("calibration_fraction must be in (0, 1)")
    _validate_metric_names(contract.metrics)

    manifest = _load_manifest(contract.split_manifest_path)
    if manifest.get("seed") != contract.split_seed:
        raise ValueError("split manifest seed does not match contract seed")

    group_scene_counts: dict[str, int] = {}
    group_row_counts: dict[str, int] = {}
    for group in contract.groups:
        rows = _load_group_rows(group)
        columns = set(rows[0].keys())
        for metric in contract.metrics:
            if metric.column not in columns:
                raise ValueError(
                    f"group {group.condition!r} missing metric column {metric.column!r}"
                )
        pixel_size = _infer_pixel_size_um(rows)
        if not math.isclose(pixel_size, contract.voxel_size_um_xyz[0], rel_tol=1e-9):
            raise ValueError(
                f"group {group.condition!r} inferred pixel size {pixel_size:.12g} "
                f"does not match contract {contract.voxel_size_um_xyz[0]:.12g}"
            )
        scenes = {row[group.scene_id_column] for row in rows}
        _validate_naming(group, rows)
        _validate_manifest_group(group, scenes, manifest)
        group_scene_counts[group.condition] = len(scenes)
        group_row_counts[group.condition] = len(rows)

    for generic in contract.to_v2_data_contracts():
        generic.validate()

    warnings: list[str] = []
    if contract.raw_imaging_available and contract.raw_path is None:
        warnings.append("raw_imaging_available true but raw_path is not yet supplied")

    return ImagingInputContractValidationResult(
        contract_id=contract.contract_id,
        group_scene_counts=group_scene_counts,
        group_row_counts=group_row_counts,
        warnings=tuple(warnings),
    )


def _require_str(payload: dict[str, Any], key: str) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{key} must be a non-empty string")
    return value


def _require_list(payload: dict[str, Any], key: str) -> list[Any]:
    value = payload.get(key)
    if not isinstance(value, list) or not value:
        raise ValueError(f"{key} must be a non-empty list")
    return value


def _float_triplet(value: Any) -> tuple[float, float, float]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError("voxel_size_um_xyz must contain exactly 3 values")
    return (float(value[0]), float(value[1]), float(value[2]))


def _resolve_relative(base_dir: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        raise ValueError(f"contract paths must be relative, got {value!r}")
    if any(part == "" for part in path.parts):
        raise ValueError(f"contract paths must not contain empty parts: {value!r}")
    resolved = (base_dir / path).resolve()
    root = _repo_root_for(base_dir)
    if not _is_relative_to(resolved, root):
        raise ValueError(
            f"contract path must not contain path traversal outside repository "
            f"root {root}: {value!r}"
        )
    return resolved


def _repo_root_for(path: Path) -> Path:
    current = path.resolve()
    for candidate in (current, *current.parents):
        if (candidate / ".git").exists():
            return candidate
    return current


def _is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def _csv_columns(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no CSV header")
        return list(reader.fieldnames)


def _n_timepoints(path: Path, *, frame_column: str) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        frames = {int(row[frame_column]) for row in reader}
    return max(frames) + 1


def _load_group_rows(group: ImagingInputContractGroup) -> list[dict[str, str]]:
    if not group.csv_path.exists():
        raise ValueError(f"CSV path does not exist: {group.csv_path}")
    with group.csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            group.scene_id_column,
            group.frame_column,
            group.time_column,
            "filename",
            "Area_um2",
            "Area_px",
        }
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"{group.condition} CSV missing required columns: {missing}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{group.condition} CSV has no rows")
    return rows


def _infer_pixel_size_um(rows: list[dict[str, str]]) -> float:
    row = rows[0]
    area_um2 = float(row["Area_um2"])
    area_px = float(row["Area_px"])
    if area_um2 <= 0.0 or area_px <= 0.0:
        raise ValueError("Area_um2 and Area_px must be positive for pixel size inference")
    return math.sqrt(area_um2 / area_px)


def _validate_naming(
    group: ImagingInputContractGroup, rows: list[dict[str, str]]
) -> None:
    for row in rows:
        scene = row[group.scene_id_column]
        filename = row["filename"]
        if not scene.startswith(group.naming_series_prefix):
            raise ValueError(
                f"{group.condition} scene {scene!r} does not match prefix "
                f"{group.naming_series_prefix!r}"
            )
        if group.filename_prefix not in filename:
            raise ValueError(
                f"{group.condition} filename {filename!r} missing prefix "
                f"{group.filename_prefix!r}"
            )


def _load_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise ValueError(f"split manifest does not exist: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != "v2_imaging_split_manifest_v1":
        raise ValueError("unsupported split manifest schema_version")
    if not isinstance(payload.get("groups"), dict):
        raise ValueError("split manifest groups must be a mapping")
    return payload


def _validate_manifest_group(
    group: ImagingInputContractGroup, scenes: set[str], manifest: dict[str, Any]
) -> None:
    if group.condition not in manifest["groups"]:
        raise ValueError(f"split manifest missing group {group.condition!r}")
    entry = manifest["groups"][group.condition]
    cal = set(entry.get("calibration_scenes", ()))
    val = set(entry.get("validation_scenes", ()))
    overlap = sorted(cal & val)
    if overlap:
        raise ValueError(f"{group.condition} calibration/validation overlap: {overlap}")
    listed = cal | val
    if listed != scenes:
        missing = sorted(scenes - listed)
        extra = sorted(listed - scenes)
        raise ValueError(
            f"{group.condition} split manifest scene mismatch; missing={missing}, extra={extra}"
        )


def _validate_metric_names(metrics: tuple[ImagingInputMetric, ...]) -> None:
    seen: set[tuple[str, str]] = set()
    for metric in metrics:
        if not metric.name.strip():
            raise ValueError("metric name must be non-empty")
        if metric.role not in {"calibration", "validation", "diagnostic"}:
            raise ValueError(f"metric {metric.name!r} role is invalid: {metric.role!r}")
        key = (metric.name, metric.measurement_modality)
        if key in seen:
            raise ValueError(f"duplicate metric key: {key!r}")
        seen.add(key)

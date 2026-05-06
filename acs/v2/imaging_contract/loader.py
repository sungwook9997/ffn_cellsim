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
from acs.v2.imaging_contract.split_builder import compute_yaml_sha256


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
class ExpectedTimepointsSpec:
    """Y15 expected timepoint constraints (per-group max + monotonic)."""

    min_frames_per_scene: int
    max_frames_per_scene_per_group: dict[str, int]
    require_monotonic_time_min: bool
    require_monotonic_frame: bool
    expected_frame_interval_min: float


@dataclass(frozen=True, slots=True)
class CalibrationValidationSplitSpec:
    """Y14: pipeline-only split with Hard Rule 1 forbidden_use guard."""

    policy: str
    ratio: tuple[float, float]
    seed: int
    manifest_path: Path
    calibration_use: str
    forbidden_use: str


@dataclass(frozen=True, slots=True)
class SidecarSpec:
    """Y9: machine-readable sidecar guard against physics-fitter loading."""

    allowed_use: str
    forbidden_use: str


@dataclass(frozen=True, slots=True)
class ImagingInputContract:
    """Loaded V2-1 imaging input contract."""

    contract_id: str
    imaging_date: str
    canonical_root: Path
    pixel_size_um_xy: tuple[float, float]
    pixel_area_um2: float
    frame_interval_min: float
    image_size_px_width_height: tuple[int, int]
    channels: tuple[str, ...]
    groups: tuple[ImagingInputContractGroup, ...]
    metrics: tuple[ImagingInputMetric, ...]
    expected_timepoints: ExpectedTimepointsSpec
    split: CalibrationValidationSplitSpec
    sidecars: dict[str, SidecarSpec]
    raw_imaging_available: bool
    raw_path: Path | None
    yaml_path: Path

    @property
    def voxel_size_um_xyz(self) -> tuple[float, float, float]:
        """Backwards-compat shim for callers expecting the old 3-tuple."""
        return (self.pixel_size_um_xy[0], self.pixel_size_um_xy[1], 1.0)

    @property
    def frame_interval_s(self) -> float:
        """Backwards-compat shim. Computed from frame_interval_min."""
        return float(self.frame_interval_min) * 60.0

    @property
    def split_manifest_path(self) -> Path:
        """Backwards-compat shim."""
        return self.split.manifest_path

    @property
    def split_seed(self) -> int:
        """Backwards-compat shim."""
        return self.split.seed

    @property
    def calibration_fraction(self) -> float:
        """Backwards-compat shim. Derived from split.ratio[0]."""
        return float(self.split.ratio[0])

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
                n_timepoints=_n_timepoints(
                    group.csv_path, frame_column=group.frame_column
                ),
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
    """Load a V2-1 imaging input YAML contract per locked schema
    ``docs/v2/v2_layer_1_imaging_input_contract_locked.md``."""

    yaml_path = Path(path)
    payload = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("imaging contract YAML must contain a mapping")
    base_dir = yaml_path.parent
    contract_id = _require_str(payload, "contract_id")
    imaging_date = _require_str(payload, "imaging_date")
    pixel_size_um_xy = _float_pair(payload.get("pixel_size_um_xy"))
    if "pixel_area_um2" not in payload:
        raise ValueError("pixel_area_um2 is required (Y1 P0-1)")
    pixel_area_um2 = float(payload["pixel_area_um2"])
    if "frame_interval_min" not in payload:
        raise ValueError(
            "frame_interval_min is required as a first-class top-level field "
            "(Y16 SEAL guard 1)"
        )
    frame_interval_min = float(payload["frame_interval_min"])
    image_size_px = _int_pair(payload.get("image_size_px_width_height"))
    channels = tuple(_require_list(payload, "channels"))

    # YAML shape validation first (no filesystem work); fail fast on wording
    # errors before doing path resolution.
    expected_payload = payload.get("expected_timepoints", {})
    if not isinstance(expected_payload, dict):
        raise ValueError("expected_timepoints must be a mapping")
    max_frames_payload = expected_payload.get(
        "max_frames_per_scene_per_group", {}
    )
    if not isinstance(max_frames_payload, dict):
        raise ValueError(
            "expected_timepoints.max_frames_per_scene_per_group must be a mapping"
        )
    expected_timepoints = ExpectedTimepointsSpec(
        min_frames_per_scene=int(expected_payload.get("min_frames_per_scene", 1)),
        max_frames_per_scene_per_group={
            str(k): int(v) for k, v in max_frames_payload.items()
        },
        require_monotonic_time_min=bool(
            expected_payload.get("require_monotonic_time_min", True)
        ),
        require_monotonic_frame=bool(
            expected_payload.get("require_monotonic_frame", True)
        ),
        expected_frame_interval_min=float(
            expected_payload.get(
                "expected_frame_interval_min", frame_interval_min
            )
        ),
    )

    split_payload = payload.get("calibration_validation_split", {})
    if not isinstance(split_payload, dict):
        raise ValueError(
            "calibration_validation_split must be a mapping (Y14 locked)"
        )
    if "forbidden_use" not in split_payload:
        raise ValueError(
            "calibration_validation_split.forbidden_use is required by "
            "Y14 / Hard Rule 1; expected 'physics_parameter_fitting'"
        )
    if split_payload["forbidden_use"] != "physics_parameter_fitting":
        raise ValueError(
            "calibration_validation_split.forbidden_use must be "
            "'physics_parameter_fitting' (Y14 / Hard Rule 1 wording guard); "
            f"got {split_payload['forbidden_use']!r}"
        )
    ratio_payload = split_payload.get("ratio", [0.7, 0.3])
    if not isinstance(ratio_payload, list) or len(ratio_payload) != 2:
        raise ValueError(
            "calibration_validation_split.ratio must be a 2-list (Y14)"
        )
    ratio = (float(ratio_payload[0]), float(ratio_payload[1]))
    split_spec = CalibrationValidationSplitSpec(
        policy=str(
            split_payload.get("policy", "condition_stratified_scene_level")
        ),
        ratio=ratio,
        seed=int(split_payload["seed"]),
        manifest_path=_resolve_relative(
            base_dir, _require_str(split_payload, "manifest")
        ),
        calibration_use=str(
            split_payload.get(
                "calibration_use",
                "pipeline_loader_metric_dashboard_development_only",
            )
        ),
        forbidden_use=str(split_payload["forbidden_use"]),
    )

    sidecars_payload = payload.get("sidecars", {})
    if not isinstance(sidecars_payload, dict):
        raise ValueError("sidecars must be a mapping (Y9)")
    sidecars: dict[str, SidecarSpec] = {}
    for key, raw_side in sidecars_payload.items():
        if not isinstance(raw_side, dict):
            raise ValueError(
                f"sidecars[{key!r}] must be a mapping (Y9)"
            )
        if "forbidden_use" not in raw_side:
            raise ValueError(
                f"sidecars[{key!r}].forbidden_use is required by Y9 / "
                f"Hard Rule 1; expected 'dynamics_parameter_fitting' for "
                f"reduced-surrogate sidecars"
            )
        sidecars[str(key)] = SidecarSpec(
            allowed_use=str(
                raw_side.get(
                    "allowed_use", "reduced_surrogate_reference_only"
                )
            ),
            forbidden_use=str(raw_side["forbidden_use"]),
        )
    if "model2_framewise" not in sidecars:
        raise ValueError(
            "sidecars.model2_framewise entry is required (Y9 reduced-surrogate "
            "guard)"
        )
    if (
        sidecars["model2_framewise"].forbidden_use
        != "dynamics_parameter_fitting"
    ):
        raise ValueError(
            "sidecars.model2_framewise.forbidden_use must be "
            "'dynamics_parameter_fitting' (Y9 / Hard Rule 1 wording guard); "
            f"got {sidecars['model2_framewise'].forbidden_use!r}"
        )

    # Filesystem-touching work last (path resolution, group/metric ctor).
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
        pixel_size_um_xy=pixel_size_um_xy,
        pixel_area_um2=pixel_area_um2,
        frame_interval_min=frame_interval_min,
        image_size_px_width_height=image_size_px,
        channels=channels,
        groups=groups,
        metrics=metrics,
        expected_timepoints=expected_timepoints,
        split=split_spec,
        sidecars=sidecars,
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
    if any(float(v) <= 0.0 for v in contract.pixel_size_um_xy):
        raise ValueError("pixel_size_um_xy values must be positive")
    if contract.pixel_area_um2 <= 0.0:
        raise ValueError("pixel_area_um2 must be positive")
    expected_area = (
        float(contract.pixel_size_um_xy[0])
        * float(contract.pixel_size_um_xy[1])
    )
    if not math.isclose(contract.pixel_area_um2, expected_area, rel_tol=1e-6):
        raise ValueError(
            f"pixel_area_um2 {contract.pixel_area_um2!r} != "
            f"pixel_size_um_xy[0] * pixel_size_um_xy[1] = {expected_area!r} "
            f"(Y1 P0-1 consistency)"
        )
    if contract.frame_interval_min <= 0.0:
        raise ValueError("frame_interval_min must be positive")
    if not math.isclose(
        contract.expected_timepoints.expected_frame_interval_min,
        contract.frame_interval_min,
        rel_tol=1e-9,
    ):
        raise ValueError(
            f"expected_timepoints.expected_frame_interval_min "
            f"{contract.expected_timepoints.expected_frame_interval_min!r} != "
            f"frame_interval_min {contract.frame_interval_min!r} (Y16 consistency)"
        )
    if not (
        0.0 < float(contract.split.ratio[0]) < 1.0
        and 0.0 < float(contract.split.ratio[1]) < 1.0
    ):
        raise ValueError(
            "calibration_validation_split.ratio entries must be in (0, 1)"
        )
    ratio_sum = float(contract.split.ratio[0]) + float(contract.split.ratio[1])
    if not math.isclose(ratio_sum, 1.0, rel_tol=1e-9, abs_tol=1e-9):
        raise ValueError(
            f"calibration_validation_split.ratio must sum to 1.0; "
            f"got {contract.split.ratio!r} (sum={ratio_sum})"
        )
    if any(v <= 0 for v in contract.image_size_px_width_height):
        raise ValueError(
            "image_size_px_width_height entries must be positive"
        )
    _validate_metric_names(contract.metrics)

    manifest = _load_manifest(contract.split.manifest_path)
    if manifest.get("seed") != contract.split.seed:
        raise ValueError("split manifest seed does not match contract seed")
    if manifest.get("contract_id") != contract.contract_id:
        raise ValueError(
            f"split manifest contract_id "
            f"{manifest.get('contract_id')!r} does not match contract "
            f"{contract.contract_id!r}"
        )
    expected_sha = compute_yaml_sha256(contract.yaml_path)
    actual_sha = manifest.get("input_yaml_sha256")
    if actual_sha != expected_sha:
        raise ValueError(
            "split manifest input_yaml_sha256 does not match current YAML "
            f"(Y13 tamper detection); expected {expected_sha!r}, "
            f"got {actual_sha!r}. Regenerate manifest via "
            "acs.v2.imaging_contract.split_builder; manual edits forbidden."
        )

    group_scene_counts: dict[str, int] = {}
    group_row_counts: dict[str, int] = {}
    for group in contract.groups:
        rows = _load_group_rows(group)
        columns = set(rows[0].keys())
        for metric in contract.metrics:
            if metric.column not in columns:
                raise ValueError(
                    f"group {group.condition!r} missing metric column "
                    f"{metric.column!r}"
                )
        pixel_size = _infer_pixel_size_um(rows)
        if not math.isclose(
            pixel_size, contract.pixel_size_um_xy[0], rel_tol=1e-9
        ):
            raise ValueError(
                f"group {group.condition!r} inferred pixel size "
                f"{pixel_size:.12g} does not match contract "
                f"{contract.pixel_size_um_xy[0]:.12g}"
            )
        scenes = {row[group.scene_id_column] for row in rows}
        _validate_naming(group, rows)
        _validate_manifest_group(group, scenes, manifest)
        _validate_expected_timepoints(group, rows, contract.expected_timepoints)
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


def _float_pair(value: Any) -> tuple[float, float]:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError("pixel_size_um_xy must contain exactly 2 values")
    return (float(value[0]), float(value[1]))


def _int_pair(value: Any) -> tuple[int, int]:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(
            "image_size_px_width_height must contain exactly 2 values"
        )
    return (int(value[0]), int(value[1]))


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


def _validate_expected_timepoints(
    group: ImagingInputContractGroup,
    rows: list[dict[str, str]],
    expected: ExpectedTimepointsSpec,
) -> None:
    """Y15 expected_timepoints enforcement.

    - Per-group max frame count: any scene with frame count exceeding
      ``max_frames_per_scene_per_group[260313_<condition>]`` raises.
    - Per-scene monotonic Frame and Time_min when required.
    - min_frames_per_scene check.
    """
    group_key = f"260313_{group.condition}"
    max_frames = expected.max_frames_per_scene_per_group.get(group_key)
    if max_frames is None:
        # Allow generic conditions (e.g. test fixtures) by falling back to the
        # plain condition key.
        max_frames = expected.max_frames_per_scene_per_group.get(group.condition)
    by_scene: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_scene.setdefault(row[group.scene_id_column], []).append(row)
    for scene_id, scene_rows in by_scene.items():
        if len(scene_rows) < expected.min_frames_per_scene:
            raise ValueError(
                f"group {group.condition!r} scene {scene_id!r} has "
                f"{len(scene_rows)} frames, below "
                f"min_frames_per_scene={expected.min_frames_per_scene}"
            )
        if max_frames is not None and len(scene_rows) > max_frames:
            raise ValueError(
                f"group {group.condition!r} scene {scene_id!r} has "
                f"{len(scene_rows)} frames, exceeds "
                f"max_frames_per_scene_per_group={max_frames} (Y15 enforcement)"
            )
        if expected.require_monotonic_frame:
            frames = [int(r[group.frame_column]) for r in scene_rows]
            if any(b <= a for a, b in zip(frames, frames[1:])):
                raise ValueError(
                    f"group {group.condition!r} scene {scene_id!r} Frame is "
                    f"not strictly monotonic"
                )
        if expected.require_monotonic_time_min:
            times = [float(r[group.time_column]) for r in scene_rows]
            if any(b <= a for a, b in zip(times, times[1:])):
                raise ValueError(
                    f"group {group.condition!r} scene {scene_id!r} "
                    f"{group.time_column} is not strictly monotonic"
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

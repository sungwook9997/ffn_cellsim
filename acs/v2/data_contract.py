"""Data contract objects for v2 image-constrained mechanobiology.

These classes describe what an imaging dataset claims to contain before any
model fitting or simulation is allowed. They are deliberately lightweight and
GPU-free so they can run on Mac, Windows, CI, or notebooks.

Sanity Gate (data_contract):
    1. Dimensional analysis: this module is dimensionless metadata only. Units
       live on the metric strings ("um", "um2", "count_per_frame", ...). dt/CFL
       is N/A.
    2. Boundary cases: empty channels, empty metrics, missing artifacts,
       missing CSV columns, missing channels, missing segmentation provenance
       for segmentation-derived artifacts, non-positive voxel sizes,
       calibration/validation timepoint overlap, out-of-range timepoints,
       duplicate `(name, measurement_modality)` metric keys, artifact layout
       path traversal/absolute paths, and `available_csv_columns` declared
       without `CSV_TABLE` artifact are all guarded.
    3. Conservation/provenance invariants: dataclasses are frozen so a
       declared dataset/contract cannot be mutated after construction; metric
       validation runs against the full dataset, not a stale channel list.
    4. Numerical sanity: all numeric fields (frame_interval_s, voxel sizes)
       are checked positive; integer timepoint indices are bounded by
       n_timepoints.
    5. Sign/sense check: artifact subset semantics (`required_artifacts`
       must be a subset of `available_artifacts`) and channel/CSV-column
       subset semantics are enforced positively, not by negation tricks.
    6. Measurement-protocol consistency: per Hard Rule 11, metrics declare
       what artifacts/columns/channels they read from. A dataset that does
       not advertise those artifacts cannot be calibrated/validated by such
       a metric. Default single-cell contract auto-filters to the dataset's
       advertised artifacts so a declared contract never lies about what is
       measurable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Sequence


class ArtifactKind(str, Enum):
    """Kinds of dataset artifacts that v2 metrics may require.

    Phase 1 expands the previously-locked 8-member contract with
    ``ECM_FIELD`` to admit ECM substrate state as a first-class artifact
    (per `docs/v2_phase1_plan_consolidated.md` §8).
    """

    CSV_TABLE = "csv_table"
    RAW_FRAMES = "raw_frames"
    SEGMENTATION_MASK = "segmentation_mask"
    BOUNDARY_CONTOURS = "boundary_contours"
    TRACKING_TABLE = "tracking_table"
    MARKER_CHANNEL = "marker_channel"
    TFM_FIELD = "tfm_field"
    EVENT_ANNOTATIONS = "event_annotations"
    ECM_FIELD = "ecm_field"


def _require_positive_triplet(name: str, values: Sequence[float]) -> None:
    if len(values) != 3:
        raise ValueError(f"{name} must contain exactly 3 values, got {len(values)}")
    bad = [v for v in values if float(v) <= 0.0]
    if bad:
        raise ValueError(f"{name} values must be positive, got {values!r}")


_CANONICAL_ARTIFACT_PATHS: dict[ArtifactKind, tuple[str, str]] = {
    ArtifactKind.CSV_TABLE: ("tables/measurements.csv", "csv"),
    ArtifactKind.RAW_FRAMES: ("raw/{channel}/t{timepoint:04d}.tif", "tiff_stack"),
    ArtifactKind.SEGMENTATION_MASK: (
        "segmentation/masks/t{timepoint:04d}.tif",
        "label_mask_tiff",
    ),
    ArtifactKind.BOUNDARY_CONTOURS: (
        "segmentation/contours/t{timepoint:04d}.json",
        "measurement_boundary_json",
    ),
    ArtifactKind.TRACKING_TABLE: ("tables/tracking.csv", "csv"),
    ArtifactKind.MARKER_CHANNEL: (
        "markers/{channel}/t{timepoint:04d}.tif",
        "tiff_stack",
    ),
    ArtifactKind.TFM_FIELD: ("tfm/t{timepoint:04d}.npz", "numpy_npz"),
    ArtifactKind.EVENT_ANNOTATIONS: ("annotations/events.csv", "csv"),
    ArtifactKind.ECM_FIELD: ("ecm/t{timepoint:04d}.npz", "numpy_npz"),
}


@dataclass(frozen=True)
class ArtifactLayoutEntry:
    """Canonical relative location for one advertised dataset artifact.

    Paths are relative to ``ImagingDatasetSpec.root_uri`` and are templates,
    not filesystem existence checks. Placeholders such as ``{channel}`` and
    ``{timepoint:04d}`` document the layout expected by loaders.
    """

    artifact: ArtifactKind
    relative_path: str
    format_hint: str
    required: bool = True

    def validate(self) -> None:
        if not isinstance(self.artifact, ArtifactKind):
            raise ValueError(
                f"artifact layout entries require ArtifactKind values, got {self.artifact!r}"
            )
        if not self.relative_path.strip():
            raise ValueError(f"layout path for {self.artifact.value!r} must be non-empty")
        if "\\" in self.relative_path:
            raise ValueError(
                f"layout path for {self.artifact.value!r} must use POSIX separators"
            )
        if self.relative_path.startswith("/"):
            raise ValueError(
                f"layout path for {self.artifact.value!r} must be relative to root_uri"
            )
        parts = self.relative_path.split("/")
        if any(part in {"", ".", ".."} for part in parts):
            raise ValueError(
                f"layout path for {self.artifact.value!r} must not contain empty, '.', "
                f"or '..' path segments: {self.relative_path!r}"
            )
        if not self.format_hint.strip():
            raise ValueError(
                f"format_hint for {self.artifact.value!r} must be non-empty"
            )


@dataclass(frozen=True)
class SegmentationProvenance:
    """How segmentation-derived masks/contours were produced.

    This records the measurement pipeline only. It must not be used to justify
    solver vertex spacing, numerical ``dx``, or physics parameters.
    """

    method: str
    software: str
    version: str
    source_artifact: ArtifactKind = ArtifactKind.RAW_FRAMES
    coordinate_convention: str = "image_um_y_down"
    postprocessing: tuple[str, ...] = field(default_factory=tuple)
    reviewer: str = ""

    def validate(self) -> None:
        if not self.method.strip():
            raise ValueError("segmentation provenance method must be non-empty")
        if not self.software.strip():
            raise ValueError("segmentation provenance software must be non-empty")
        if not self.version.strip():
            raise ValueError("segmentation provenance version must be non-empty")
        if not isinstance(self.source_artifact, ArtifactKind):
            raise ValueError(
                "segmentation provenance source_artifact must be an ArtifactKind"
            )
        if self.coordinate_convention not in {
            "image_um_y_down",
            "world_um_y_up",
        }:
            raise ValueError(
                "segmentation provenance coordinate_convention must be "
                "image_um_y_down or world_um_y_up"
            )
        if any(not step.strip() for step in self.postprocessing):
            raise ValueError(
                "segmentation provenance postprocessing entries must be non-empty"
            )


def canonical_artifact_layout(
    artifacts: Sequence[ArtifactKind],
) -> tuple[ArtifactLayoutEntry, ...]:
    """Return canonical relative path templates for advertised artifacts."""

    entries: list[ArtifactLayoutEntry] = []
    for artifact in artifacts:
        if not isinstance(artifact, ArtifactKind):
            raise ValueError(
                f"canonical_artifact_layout requires ArtifactKind values, got {artifact!r}"
            )
        relative_path, format_hint = _CANONICAL_ARTIFACT_PATHS[artifact]
        entries.append(
            ArtifactLayoutEntry(
                artifact=artifact,
                relative_path=relative_path,
                format_hint=format_hint,
            )
        )
    return tuple(entries)


@dataclass(frozen=True)
class ImagingDatasetSpec:
    """Metadata needed to interpret a confocal/live-imaging dataset."""

    dataset_id: str
    root_uri: str
    modality: str
    channels: tuple[str, ...]
    voxel_size_um_xyz: tuple[float, float, float]
    frame_interval_s: float
    n_timepoints: int
    calibration_timepoints: tuple[int, ...] = field(default_factory=tuple)
    validation_timepoints: tuple[int, ...] = field(default_factory=tuple)
    available_artifacts: tuple[ArtifactKind, ...] = field(default_factory=tuple)
    available_csv_columns: tuple[str, ...] = field(default_factory=tuple)
    artifact_layout: tuple[ArtifactLayoutEntry, ...] = field(default_factory=tuple)
    segmentation_provenance: SegmentationProvenance | None = None
    notes: str = ""

    def validate(self) -> None:
        if not self.dataset_id.strip():
            raise ValueError("dataset_id must be non-empty")
        if not self.root_uri.strip():
            raise ValueError("root_uri must be non-empty")
        if not self.modality.strip():
            raise ValueError("modality must be non-empty")
        if not self.channels:
            raise ValueError("at least one imaging channel is required")
        if any(not c.strip() for c in self.channels):
            raise ValueError(f"channels must be non-empty strings, got {self.channels!r}")
        _require_positive_triplet("voxel_size_um_xyz", self.voxel_size_um_xyz)
        if self.frame_interval_s <= 0.0:
            raise ValueError("frame_interval_s must be positive")
        if self.n_timepoints <= 0:
            raise ValueError("n_timepoints must be positive")

        cal = set(self.calibration_timepoints)
        val = set(self.validation_timepoints)
        overlap = sorted(cal.intersection(val))
        if overlap:
            raise ValueError(f"calibration and validation timepoints overlap: {overlap}")
        for idx in sorted(cal.union(val)):
            if idx < 0 or idx >= self.n_timepoints:
                raise ValueError(
                    f"timepoint index {idx} outside dataset range [0, {self.n_timepoints})"
                )

        seen_artifacts: set[ArtifactKind] = set()
        for art in self.available_artifacts:
            if not isinstance(art, ArtifactKind):
                raise ValueError(
                    f"available_artifacts must contain ArtifactKind values, got {art!r}"
                )
            if art in seen_artifacts:
                raise ValueError(f"duplicate artifact in available_artifacts: {art!r}")
            seen_artifacts.add(art)

        if self.available_csv_columns and ArtifactKind.CSV_TABLE not in self.available_artifacts:
            raise ValueError(
                "available_csv_columns is non-empty but CSV_TABLE artifact is not advertised"
            )
        if any(not col.strip() for col in self.available_csv_columns):
            raise ValueError(
                f"available_csv_columns must be non-empty strings, got {self.available_csv_columns!r}"
            )
        if len(set(self.available_csv_columns)) != len(self.available_csv_columns):
            raise ValueError(
                f"duplicate CSV column in available_csv_columns: {self.available_csv_columns!r}"
            )
        self._validate_artifact_layout()
        self._validate_segmentation_provenance()

    def _validate_segmentation_provenance(self) -> None:
        segmentation_artifacts = {
            ArtifactKind.SEGMENTATION_MASK,
            ArtifactKind.BOUNDARY_CONTOURS,
        }
        if segmentation_artifacts.intersection(self.available_artifacts):
            if self.segmentation_provenance is None:
                raise ValueError(
                    "segmentation_provenance is required when segmentation_mask "
                    "or boundary_contours artifacts are advertised"
                )
        if self.segmentation_provenance is None:
            return
        self.segmentation_provenance.validate()

    def _validate_artifact_layout(self) -> None:
        if not self.artifact_layout:
            return
        advertised = set(self.available_artifacts)
        seen: set[ArtifactKind] = set()
        for entry in self.artifact_layout:
            entry.validate()
            if entry.artifact in seen:
                raise ValueError(
                    f"duplicate layout entry for artifact: {entry.artifact.value!r}"
                )
            seen.add(entry.artifact)
        extra = sorted(a.value for a in seen - advertised)
        if extra:
            raise ValueError(
                f"artifact_layout contains artifacts not advertised by dataset "
                f"{self.dataset_id!r}: {extra!r}"
            )
        missing = sorted(a.value for a in advertised - seen)
        if missing:
            raise ValueError(
                f"artifact_layout is missing advertised artifacts for dataset "
                f"{self.dataset_id!r}: {missing!r}"
            )

    def canonical_artifact_layout(self) -> tuple[ArtifactLayoutEntry, ...]:
        """Return declared layout or canonical defaults for advertised artifacts."""

        if self.artifact_layout:
            return self.artifact_layout
        return canonical_artifact_layout(self.available_artifacts)

    def artifact_relative_path(self, artifact: ArtifactKind) -> str:
        """Return the relative path template for one advertised artifact."""

        for entry in self.canonical_artifact_layout():
            if entry.artifact == artifact:
                return entry.relative_path
        raise ValueError(
            f"dataset {self.dataset_id!r} does not advertise artifact {artifact.value!r}"
        )

    def artifact_uri(self, artifact: ArtifactKind) -> str:
        """Return the root-relative URI template for one advertised artifact."""

        root = self.root_uri.rstrip("/")
        return f"{root}/{self.artifact_relative_path(artifact)}"



@dataclass(frozen=True)
class MetricSpec:
    """A visual or physical metric used for calibration or validation.

    Artifact-aware semantics (Phase 1 consolidated plan §8): a metric declares
    which dataset artifacts it reads, which CSV columns it requires (only
    meaningful when ``ArtifactKind.CSV_TABLE`` is also required), and which
    imaging channels it depends on. Validation checks that the metric's
    declared inputs are all present in the dataset.
    """

    name: str
    target_object: str
    measurement_modality: str
    unit: str
    role: str
    required_artifacts: tuple[ArtifactKind, ...] = field(default_factory=tuple)
    csv_columns_required: tuple[str, ...] = field(default_factory=tuple)
    channels_required: tuple[str, ...] = field(default_factory=tuple)
    description: str = ""

    def validate(self, dataset: ImagingDatasetSpec) -> None:
        if not self.name.strip():
            raise ValueError("metric name must be non-empty")
        if not self.target_object.strip():
            raise ValueError(f"metric {self.name!r} target_object must be non-empty")
        if not self.measurement_modality.strip():
            raise ValueError(f"metric {self.name!r} measurement_modality must be non-empty")
        if not self.unit.strip():
            raise ValueError(f"metric {self.name!r} unit must be non-empty")
        if self.role not in {"calibration", "validation", "diagnostic"}:
            raise ValueError(
                f"metric {self.name!r} role must be calibration, validation, or diagnostic"
            )

        for art in self.required_artifacts:
            if not isinstance(art, ArtifactKind):
                raise ValueError(
                    f"metric {self.name!r} required_artifacts must contain ArtifactKind values, "
                    f"got {art!r}"
                )
        missing_artifacts = sorted(
            (a.value for a in set(self.required_artifacts) - set(dataset.available_artifacts)),
        )
        if missing_artifacts:
            raise ValueError(
                f"metric {self.name!r} requires artifacts {missing_artifacts!r} "
                f"not advertised by dataset {dataset.dataset_id!r}"
            )

        if self.csv_columns_required:
            if ArtifactKind.CSV_TABLE not in self.required_artifacts:
                raise ValueError(
                    f"metric {self.name!r} declares csv_columns_required but does not "
                    f"include CSV_TABLE in required_artifacts"
                )
            missing_columns = sorted(
                set(self.csv_columns_required) - set(dataset.available_csv_columns),
            )
            if missing_columns:
                raise ValueError(
                    f"metric {self.name!r} requires CSV columns {missing_columns!r} "
                    f"not available in dataset {dataset.dataset_id!r}"
                )

        missing_channels = sorted(
            set(self.channels_required) - set(dataset.channels),
        )
        if missing_channels:
            raise ValueError(
                f"metric {self.name!r} requires channels {missing_channels!r} "
                f"not in dataset {dataset.dataset_id!r}"
            )


@dataclass(frozen=True)
class V2DataContract:
    """Complete data contract for a v2 modeling dataset."""

    dataset: ImagingDatasetSpec
    metrics: tuple[MetricSpec, ...]

    def validate(self) -> None:
        self.dataset.validate()
        if not self.metrics:
            raise ValueError("at least one metric is required")
        seen_keys: set[tuple[str, str]] = set()
        for metric in self.metrics:
            metric.validate(self.dataset)
            key = (metric.name, metric.measurement_modality)
            if key in seen_keys:
                raise ValueError(f"duplicate metric key (name, measurement_modality): {key!r}")
            seen_keys.add(key)

    @property
    def calibration_metrics(self) -> tuple[MetricSpec, ...]:
        return tuple(m for m in self.metrics if m.role == "calibration")

    @property
    def validation_metrics(self) -> tuple[MetricSpec, ...]:
        return tuple(m for m in self.metrics if m.role == "validation")


def default_single_cell_contract(dataset: ImagingDatasetSpec) -> V2DataContract:
    """Return a minimal single-cell metric contract for early v2 work.

    Auto-filters candidate metrics to those whose ``required_artifacts`` are
    a subset of ``dataset.available_artifacts``. A metric that depends on an
    unadvertised artifact is silently omitted, so a declared contract never
    claims a metric the dataset cannot supply.
    """

    candidates: tuple[MetricSpec, ...] = (
        MetricSpec(
            name="projected_area",
            target_object="single_cell",
            measurement_modality="top_down_segmentation",
            unit="um2",
            role="calibration",
            required_artifacts=(ArtifactKind.SEGMENTATION_MASK,),
            description="Top-down projected cell area from segmentation masks.",
        ),
        MetricSpec(
            name="boundary_roughness",
            target_object="single_cell",
            measurement_modality="cell_outline",
            unit="dimensionless",
            role="validation",
            required_artifacts=(ArtifactKind.BOUNDARY_CONTOURS,),
            description="Perimeter-normalized outline roughness for protrusive activity.",
        ),
        MetricSpec(
            name="protrusion_event_count",
            target_object="single_cell_boundary",
            measurement_modality="live_imaging_tracking",
            unit="count_per_frame",
            role="validation",
            required_artifacts=(
                ArtifactKind.EVENT_ANNOTATIONS,
                ArtifactKind.TRACKING_TABLE,
            ),
            description="Detected lamellipodia/filopodia-like boundary events.",
        ),
    )
    available = set(dataset.available_artifacts)
    metrics = tuple(
        m for m in candidates if set(m.required_artifacts).issubset(available)
    )
    return V2DataContract(dataset=dataset, metrics=metrics)

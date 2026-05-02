"""Data contract objects for v2 image-constrained mechanobiology.

These classes describe what an imaging dataset claims to contain before any
model fitting or simulation is allowed. They are deliberately lightweight and
GPU-free so they can run on Mac, Windows, CI, or notebooks.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence


def _require_positive_triplet(name: str, values: Sequence[float]) -> None:
    if len(values) != 3:
        raise ValueError(f"{name} must contain exactly 3 values, got {len(values)}")
    bad = [v for v in values if float(v) <= 0.0]
    if bad:
        raise ValueError(f"{name} values must be positive, got {values!r}")


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


@dataclass(frozen=True)
class MetricSpec:
    """A visual or physical metric used for calibration or validation."""

    name: str
    target_object: str
    measurement_modality: str
    unit: str
    role: str
    source_channel: Optional[str] = None
    description: str = ""

    def validate(self, available_channels: Sequence[str]) -> None:
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
        if self.source_channel is not None and self.source_channel not in available_channels:
            raise ValueError(
                f"metric {self.name!r} source_channel {self.source_channel!r} "
                f"not in dataset channels {tuple(available_channels)!r}"
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
        seen: set[str] = set()
        for metric in self.metrics:
            metric.validate(self.dataset.channels)
            if metric.name in seen:
                raise ValueError(f"duplicate metric name: {metric.name!r}")
            seen.add(metric.name)

    @property
    def calibration_metrics(self) -> tuple[MetricSpec, ...]:
        return tuple(m for m in self.metrics if m.role == "calibration")

    @property
    def validation_metrics(self) -> tuple[MetricSpec, ...]:
        return tuple(m for m in self.metrics if m.role == "validation")


def default_single_cell_contract(dataset: ImagingDatasetSpec) -> V2DataContract:
    """Return a minimal single-cell metric contract for early v2 work."""

    return V2DataContract(
        dataset=dataset,
        metrics=(
            MetricSpec(
                name="projected_area",
                target_object="single_cell",
                measurement_modality="top_down_segmentation",
                unit="um2",
                role="calibration",
                description="Top-down projected cell area from segmentation masks.",
            ),
            MetricSpec(
                name="boundary_roughness",
                target_object="single_cell",
                measurement_modality="cell_outline",
                unit="dimensionless",
                role="validation",
                description="Perimeter-normalized outline roughness for protrusive activity.",
            ),
            MetricSpec(
                name="protrusion_event_count",
                target_object="single_cell_boundary",
                measurement_modality="live_imaging_tracking",
                unit="count_per_frame",
                role="validation",
                description="Detected lamellipodia/filopodia-like boundary events.",
            ),
        ),
    )

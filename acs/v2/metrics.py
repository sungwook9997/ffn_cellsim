"""Pure metric registry keyed by ``(metric_name, ArtifactKind)``.

The registry is the v2 single source of truth for how a named metric reads
its inputs. Each registered function is pure: it takes the artifact-typed
input (e.g. a :class:`MeasurementBoundary` for ``BOUNDARY_CONTOURS`` or
``SEGMENTATION_MASK`` derived outlines) and returns a numeric value with
explicit units. Dynamics, fitting, or stateful caches are out of scope.

This module is non-physics metadata + thin pure functions, so the full
physics Sanity Gate is not required. The relevant guarantees are:

- Inputs declare their `ArtifactKind`; metrics never accept the wrong
  artifact silently.
- Boundary metrics always defensively call `MeasurementBoundary.validate()`
  before computing, so direct-construction boundaries cannot bypass
  geometry checks.
- Outputs are reported in the same units listed in
  `acs.v2.data_contract.MetricSpec.unit` strings (`um2`, `um`, ...).
- Measurement-protocol consistency (Hard Rule 11): boundary
  ``projected_area`` is the top-down xy shoelace area, matching PI
  experimental segmentation modality.

Magic-Number Block: this module declares no tunable numeric. N/A.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from acs.v2.data_contract import ArtifactKind
from acs.v2.measurement_boundary import MeasurementBoundary

MetricKey = tuple[str, ArtifactKind]


class MetricRegistryError(KeyError):
    """Raised when a metric is missing or registered twice."""


@dataclass(frozen=True)
class RegisteredMetric:
    """A pure metric function bound to a name and an artifact kind."""

    name: str
    artifact_kind: ArtifactKind
    unit: str
    fn: Callable[..., float]
    description: str = ""

    def __call__(self, *args, **kwargs) -> float:
        return float(self.fn(*args, **kwargs))


class MetricRegistry:
    """In-process registry of pure metric functions.

    Keys are ``(metric_name, artifact_kind)`` so the same metric name can
    be backed by multiple artifact kinds (e.g. boundary ``projected_area``
    from segmentation polygon vs CSV column lookup).
    """

    def __init__(self) -> None:
        self._entries: dict[MetricKey, RegisteredMetric] = {}

    def register(self, metric: RegisteredMetric) -> None:
        if not isinstance(metric.artifact_kind, ArtifactKind):
            raise MetricRegistryError(
                f"metric artifact_kind must be ArtifactKind, got {metric.artifact_kind!r}"
            )
        if not metric.name.strip():
            raise MetricRegistryError("metric name must be non-empty")
        if not metric.unit.strip():
            raise MetricRegistryError(f"metric {metric.name!r} unit must be non-empty")
        if not callable(metric.fn):
            raise MetricRegistryError(f"metric {metric.name!r} fn must be callable")
        key: MetricKey = (metric.name, metric.artifact_kind)
        if key in self._entries:
            raise MetricRegistryError(f"metric already registered: {key!r}")
        self._entries[key] = metric

    def get(self, name: str, artifact_kind: ArtifactKind) -> RegisteredMetric:
        key: MetricKey = (name, artifact_kind)
        if key not in self._entries:
            raise MetricRegistryError(f"metric not registered: {key!r}")
        return self._entries[key]

    def keys(self) -> tuple[MetricKey, ...]:
        return tuple(self._entries.keys())


def boundary_projected_area_um2(boundary: MeasurementBoundary) -> float:
    """Top-down xy projected area of a measurement boundary in ``um2``.

    Defensively re-validates the boundary so a directly-constructed
    boundary (which starts ``_validated=False``) cannot bypass geometry
    checks via the registry.
    """

    boundary.validate()
    return boundary.projected_area_um2()


def boundary_perimeter_um(boundary: MeasurementBoundary) -> float:
    """Closed-polygon perimeter of a measurement boundary in ``um``."""

    boundary.validate()
    return boundary.perimeter_um()


def default_registry() -> MetricRegistry:
    """Return a fresh registry pre-populated with v2 boundary metrics."""

    registry = MetricRegistry()
    registry.register(
        RegisteredMetric(
            name="projected_area",
            artifact_kind=ArtifactKind.BOUNDARY_CONTOURS,
            unit="um2",
            fn=boundary_projected_area_um2,
            description="Top-down xy projected area from a measurement boundary polygon.",
        )
    )
    registry.register(
        RegisteredMetric(
            name="projected_area",
            artifact_kind=ArtifactKind.SEGMENTATION_MASK,
            unit="um2",
            fn=boundary_projected_area_um2,
            description="Top-down xy projected area from a segmentation-derived boundary polygon.",
        )
    )
    registry.register(
        RegisteredMetric(
            name="perimeter",
            artifact_kind=ArtifactKind.BOUNDARY_CONTOURS,
            unit="um",
            fn=boundary_perimeter_um,
            description="Closed-polygon perimeter from a measurement boundary.",
        )
    )
    return registry

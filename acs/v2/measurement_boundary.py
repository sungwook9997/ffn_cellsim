"""Measurement-space boundary geometry for v2 single-cell contracts.

The boundary stored by :class:`MeasurementBoundary` is canonical simulation
measurement geometry: xy vertices in micrometers, y-up world coordinates, and
counter-clockwise orientation. Source-image conventions are preserved as
provenance so image-derived outlines can be round-tripped without confusing
measurement semantics.

This validator is intended for measurement boundaries with tens to hundreds of
vertices. The simple-polygon topology check is O(N^2) and is not a dense-mesh
intersection engine.

Sanity Gate:
    1. Dimensional analysis: vertices are in um, area is um2, perimeter and
       edge lengths are um. This module has no time stepping, velocity, stress,
       mass, CFL, Re, Ca, De, Pe, or Ma term; dt/CFL is N/A.
    2. Boundary cases: N < 3, non-finite coordinates, zero or negative
       geometric floor, degenerate/collinear polygons, duplicate vertices,
       vertex pinches, self-intersections, bad conventions, tiny positive
       polygons, and large-coordinate polygons are explicitly validated.
    3. Conservation/provenance invariants: source signed area, dtype, y-flip,
       and orientation-reversal flags are preserved. Exporting to the source
       convention inverts only the coordinate-convention transform, not the
       biological geometry.
    4. Numerical sanity: float64 is used for all geometry and diagnostics;
       input dtype is recorded. The default min_edge_um = 1e-3 is a non-zero
       geometric floor in micrometers, not a fitted physical factor. The
       collinearity epsilon is a scale-aware numerical tie-break for segment
       predicates, not a tuning knob.
    5. Sign/sense check: world_um_y_up CCW polygons have positive signed area.
       image_um_y_down sources are transformed by y := -y before orientation
       checks, so handedness is evaluated in the same canonical measurement
       space used by downstream simulation metrics.
    6. Measurement-protocol consistency: projected area is the xy shoelace
       area of the canonical outline. It matches top-down segmentation
       measurements and does not substitute substrate contact area or any
       z-dependent metric.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

CoordinateConvention = Literal["world_um_y_up", "image_um_y_down"]
FailureKind = Literal[
    "shape",
    "min_vertices",
    "non_finite",
    "degenerate_area",
    "wrong_orientation",
    "duplicate_vertex",
    "self_intersection",
    "vertex_pinch",
    "edge_too_short",
    "convention",
    "min_edge_param",
]


class MeasurementBoundaryError(ValueError):
    """Validation error carrying a machine-readable failure kind."""

    def __init__(self, failure_kind: FailureKind, message: str) -> None:
        self.failure_kind = failure_kind
        super().__init__(message)


@dataclass(frozen=True, slots=True)
class MeasurementBoundary:
    """Canonical cell boundary measured in xy micrometers.

    Attributes:
        vertices_xy_um: Canonical xy vertices in um, world y-up, CCW.
        coordinate_convention: Source coordinate convention used for export.
        source_modality: Non-empty modality/provenance label.
        object_id: Optional cell/object identifier.
        time_s: Optional source time in seconds.
        min_edge_um: Minimum allowed edge length in um.
        orientation_was_reversed: Whether constructor reversed source order.
        coord_was_flipped: Whether y was flipped from image y-down source.
        source_signed_area_raw: Raw input shoelace area before transforms.
        input_dtype: String representation of the source input dtype.
        _validated: Internal validation cache. Direct construction starts false.
    """

    vertices_xy_um: NDArray[np.float64]
    coordinate_convention: CoordinateConvention
    source_modality: str
    object_id: Optional[str] = None
    time_s: Optional[float] = None
    min_edge_um: float = 1e-3
    orientation_was_reversed: bool = False
    coord_was_flipped: bool = False
    source_signed_area_raw: Optional[float] = None
    input_dtype: Optional[str] = None
    _validated: bool = False

    @classmethod
    def from_array(
        cls,
        raw_vertices: ArrayLike,
        coordinate_convention: CoordinateConvention,
        source_modality: str,
        *,
        auto_orient: bool = False,
        object_id: Optional[str] = None,
        time_s: Optional[float] = None,
        min_edge_um: float = 1e-3,
    ) -> "MeasurementBoundary":
        """Build a validated canonical boundary from source vertices.

        The default ``min_edge_um`` is a permissive numerical floor for manual
        or synthetic outlines. Image-derived callers should override it using
        modality/resolution context, for example roughly 0.25 um for confocal
        lateral resolution or 0.05 um for super-resolution segmentations.
        """

        cls._validate_convention(coordinate_convention)
        cls._validate_min_edge(min_edge_um)
        if not source_modality.strip():
            raise MeasurementBoundaryError(
                "shape", "source_modality must be a non-empty string"
            )

        raw = np.asarray(raw_vertices)
        input_dtype = str(raw.dtype)
        vertices = cls._coerce_vertices(raw)
        source_signed_area_raw = _signed_area(vertices)

        coord_was_flipped = coordinate_convention == "image_um_y_down"
        if coord_was_flipped:
            vertices = vertices.copy()
            vertices[:, 1] *= -1.0

        boundary = cls(
            vertices_xy_um=vertices,
            coordinate_convention=coordinate_convention,
            source_modality=source_modality,
            object_id=object_id,
            time_s=time_s,
            min_edge_um=min_edge_um,
            coord_was_flipped=coord_was_flipped,
            source_signed_area_raw=source_signed_area_raw,
            input_dtype=input_dtype,
        )
        boundary._validate_edges(vertices)
        boundary._validate_simple_polygon(vertices)

        signed_area = _signed_area(vertices)
        if abs(signed_area) <= _area_epsilon(vertices):
            raise MeasurementBoundaryError(
                "degenerate_area", "boundary signed area is zero or numerically degenerate"
            )

        orientation_was_reversed = False
        if signed_area <= 0.0:
            if not auto_orient:
                raise MeasurementBoundaryError(
                    "wrong_orientation",
                    "boundary must be counter-clockwise in canonical world_um_y_up coordinates",
                )
            vertices = vertices[::-1].copy()
            orientation_was_reversed = True

        oriented_boundary = cls(
            vertices_xy_um=vertices,
            coordinate_convention=coordinate_convention,
            source_modality=source_modality,
            object_id=object_id,
            time_s=time_s,
            min_edge_um=min_edge_um,
            orientation_was_reversed=orientation_was_reversed,
            coord_was_flipped=coord_was_flipped,
            source_signed_area_raw=source_signed_area_raw,
            input_dtype=input_dtype,
        )
        oriented_boundary.validate()
        return oriented_boundary

    def validate(self) -> None:
        """Validate and promote this boundary to canonical cached geometry."""

        self._validate_convention(self.coordinate_convention)
        self._validate_min_edge(self.min_edge_um)
        if not self.source_modality.strip():
            raise MeasurementBoundaryError(
                "shape", "source_modality must be a non-empty string"
            )

        vertices = self._coerce_vertices(self.vertices_xy_um)
        self._validate_edges(vertices)
        self._validate_simple_polygon(vertices)
        signed_area = _signed_area(vertices)
        if abs(signed_area) <= _area_epsilon(vertices):
            raise MeasurementBoundaryError(
                "degenerate_area", "boundary signed area is zero or numerically degenerate"
            )
        if signed_area <= 0.0:
            raise MeasurementBoundaryError(
                "wrong_orientation",
                "boundary must be counter-clockwise in canonical world_um_y_up coordinates",
            )

        vertices.setflags(write=False)
        object.__setattr__(self, "vertices_xy_um", vertices)
        object.__setattr__(self, "_validated", True)

    def signed_area_um2(self) -> float:
        """Return the canonical signed xy area in um2."""

        if not self._validated:
            self.validate()
        return _signed_area(self.vertices_xy_um)

    def projected_area_um2(self) -> float:
        """Return the top-down projected xy area in um2."""

        return abs(self.signed_area_um2())

    def perimeter_um(self) -> float:
        """Return the closed-polygon perimeter in um."""

        if not self._validated:
            self.validate()
        diffs = np.roll(self.vertices_xy_um, -1, axis=0) - self.vertices_xy_um
        return float(np.linalg.norm(diffs, axis=1).sum())

    def export_in_source_convention(self) -> NDArray[np.float64]:
        """Return canonical geometry expressed in the source convention.

        The returned vertices draw the same canonical polygon in the original
        coordinate convention. Vertex order or starting index may differ from
        the raw input when ``orientation_was_reversed`` is true.
        """

        if not self._validated:
            self.validate()
        exported = np.array(self.vertices_xy_um, dtype=np.float64, copy=True)
        if self.coordinate_convention == "image_um_y_down":
            exported[:, 1] *= -1.0
        exported.setflags(write=False)
        return exported

    @staticmethod
    def _validate_convention(convention: str) -> None:
        if convention not in {"world_um_y_up", "image_um_y_down"}:
            raise MeasurementBoundaryError(
                "convention",
                "coordinate_convention must be world_um_y_up or image_um_y_down",
            )

    @staticmethod
    def _validate_min_edge(min_edge_um: float) -> None:
        if not np.isfinite(min_edge_um) or min_edge_um <= 0.0:
            raise MeasurementBoundaryError(
                "min_edge_param", "min_edge_um must be finite and positive"
            )

    @staticmethod
    def _coerce_vertices(vertices: ArrayLike) -> NDArray[np.float64]:
        arr = np.asarray(vertices, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise MeasurementBoundaryError(
                "shape", "vertices_xy_um must have shape (N, 2)"
            )
        if arr.shape[0] < 3:
            raise MeasurementBoundaryError(
                "min_vertices", "vertices_xy_um must contain at least 3 vertices"
            )
        if not np.isfinite(arr).all():
            raise MeasurementBoundaryError(
                "non_finite", "vertices_xy_um must contain only finite coordinates"
            )
        if arr.dtype != np.float64 or arr.flags.writeable:
            arr = np.array(arr, dtype=np.float64, copy=True)
        return arr

    def _validate_edges(self, vertices: NDArray[np.float64]) -> None:
        diffs = np.roll(vertices, -1, axis=0) - vertices
        lengths = np.linalg.norm(diffs, axis=1)
        if np.any(lengths == 0.0):
            raise MeasurementBoundaryError(
                "duplicate_vertex", "consecutive duplicate vertices are not allowed"
            )
        if np.any(lengths < self.min_edge_um):
            raise MeasurementBoundaryError(
                "edge_too_short", f"all edges must be at least {self.min_edge_um} um"
            )

    def _validate_simple_polygon(self, vertices: NDArray[np.float64]) -> None:
        n_vertices = vertices.shape[0]
        for i in range(n_vertices):
            for j in range(i + 1, n_vertices):
                if _points_close(vertices[i], vertices[j]):
                    if abs(i - j) == 1 or {i, j} == {0, n_vertices - 1}:
                        raise MeasurementBoundaryError(
                            "duplicate_vertex",
                            "consecutive duplicate vertices are not allowed",
                        )
                    raise MeasurementBoundaryError(
                        "vertex_pinch", "non-consecutive repeated vertices pinch the polygon"
                    )

        for i in range(n_vertices):
            a1 = vertices[i]
            a2 = vertices[(i + 1) % n_vertices]
            for j in range(i + 1, n_vertices):
                if _segments_are_adjacent(i, j, n_vertices):
                    continue
                b1 = vertices[j]
                b2 = vertices[(j + 1) % n_vertices]
                if _segments_intersect(a1, a2, b1, b2):
                    raise MeasurementBoundaryError(
                        "self_intersection", "polygon edges must not self-intersect"
                    )


def _signed_area(vertices: NDArray[np.float64]) -> float:
    x = vertices[:, 0]
    y = vertices[:, 1]
    return float(0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _area_epsilon(vertices: NDArray[np.float64]) -> float:
    scale = max(1.0, float(np.max(np.abs(vertices))))
    return np.finfo(np.float64).eps * scale * scale * vertices.shape[0] * 16.0


def _orientation_epsilon(*points: NDArray[np.float64]) -> float:
    coords = np.vstack(points)
    scale = max(1.0, float(np.max(np.abs(coords))))
    return np.finfo(np.float64).eps * scale * scale * 16.0


def _orientation(
    a: NDArray[np.float64], b: NDArray[np.float64], c: NDArray[np.float64]
) -> float:
    return float((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def _points_close(a: NDArray[np.float64], b: NDArray[np.float64]) -> bool:
    return bool(np.array_equal(a, b))


def _on_segment(
    a: NDArray[np.float64], b: NDArray[np.float64], c: NDArray[np.float64]
) -> bool:
    eps = _orientation_epsilon(a, b, c)
    return bool(
        min(a[0], c[0]) - eps <= b[0] <= max(a[0], c[0]) + eps
        and min(a[1], c[1]) - eps <= b[1] <= max(a[1], c[1]) + eps
    )


def _segments_intersect(
    a1: NDArray[np.float64],
    a2: NDArray[np.float64],
    b1: NDArray[np.float64],
    b2: NDArray[np.float64],
) -> bool:
    eps = _orientation_epsilon(a1, a2, b1, b2)
    o1 = _orientation(a1, a2, b1)
    o2 = _orientation(a1, a2, b2)
    o3 = _orientation(b1, b2, a1)
    o4 = _orientation(b1, b2, a2)

    if ((o1 > eps and o2 < -eps) or (o1 < -eps and o2 > eps)) and (
        (o3 > eps and o4 < -eps) or (o3 < -eps and o4 > eps)
    ):
        return True
    if abs(o1) <= eps and _on_segment(a1, b1, a2):
        return True
    if abs(o2) <= eps and _on_segment(a1, b2, a2):
        return True
    if abs(o3) <= eps and _on_segment(b1, a1, b2):
        return True
    if abs(o4) <= eps and _on_segment(b1, a2, b2):
        return True
    return False


def _segments_are_adjacent(i: int, j: int, n_vertices: int) -> bool:
    return i == j or abs(i - j) == 1 or {i, j} == {0, n_vertices - 1}

from __future__ import annotations

import numpy as np
import pytest

from acs.v2.data_contract import ArtifactKind
from acs.v2.measurement_boundary import MeasurementBoundary
from acs.v2.metrics import (
    MetricRegistry,
    MetricRegistryError,
    RegisteredMetric,
    boundary_perimeter_um,
    boundary_projected_area_um2,
    csv_scalar_metric,
    csv_spheroid_a_over_a0,
    csv_spheroid_area_um2,
    csv_spheroid_effective_radius_um,
    default_registry,
)


def _square_boundary(side_um: float = 2.0) -> MeasurementBoundary:
    return MeasurementBoundary.from_array(
        np.array(
            [
                [0.0, 0.0],
                [side_um, 0.0],
                [side_um, side_um],
                [0.0, side_um],
            ]
        ),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )


def test_default_registry_has_six_entries_keyed_by_name_and_artifact():
    registry = default_registry()
    expected = {
        ("projected_area", ArtifactKind.BOUNDARY_CONTOURS),
        ("projected_area", ArtifactKind.SEGMENTATION_MASK),
        ("perimeter", ArtifactKind.BOUNDARY_CONTOURS),
        ("spheroid_projected_area", ArtifactKind.CSV_TABLE),
        ("spheroid_a_over_a0", ArtifactKind.CSV_TABLE),
        ("spheroid_effective_radius", ArtifactKind.CSV_TABLE),
    }
    assert set(registry.keys()) == expected


def test_default_registry_projected_area_matches_boundary():
    registry = default_registry()
    boundary = _square_boundary(side_um=2.0)
    metric = registry.get("projected_area", ArtifactKind.BOUNDARY_CONTOURS)
    assert metric.unit == "um2"
    assert metric(boundary) == pytest.approx(4.0)


def test_default_registry_perimeter_matches_boundary():
    registry = default_registry()
    boundary = _square_boundary(side_um=3.0)
    metric = registry.get("perimeter", ArtifactKind.BOUNDARY_CONTOURS)
    assert metric.unit == "um"
    assert metric(boundary) == pytest.approx(12.0)


def test_registry_lookup_unknown_key_raises():
    registry = default_registry()
    with pytest.raises(MetricRegistryError):
        registry.get("missing_metric", ArtifactKind.SEGMENTATION_MASK)
    with pytest.raises(MetricRegistryError):
        registry.get("projected_area", ArtifactKind.CSV_TABLE)


def test_default_registry_spheroid_csv_metrics_return_scalar_columns():
    registry = default_registry()
    row = {
        "Area_um2": "1250.5",
        "A_over_A0": 1.37,
        "EffectiveRadius_um": 19.95,
    }

    area = registry.get("spheroid_projected_area", ArtifactKind.CSV_TABLE)
    ratio = registry.get("spheroid_a_over_a0", ArtifactKind.CSV_TABLE)
    radius = registry.get("spheroid_effective_radius", ArtifactKind.CSV_TABLE)

    assert area.unit == "um2"
    assert ratio.unit == "dimensionless"
    assert radius.unit == "um"
    assert area(row) == pytest.approx(1250.5)
    assert ratio(row) == pytest.approx(1.37)
    assert radius(row) == pytest.approx(19.95)


def test_csv_spheroid_metric_helpers_are_column_specific():
    row = {
        "Area_um2": 100.0,
        "A_over_A0": "1.25",
        "EffectiveRadius_um": 5.64,
    }

    assert csv_spheroid_area_um2(row) == pytest.approx(100.0)
    assert csv_spheroid_a_over_a0(row) == pytest.approx(1.25)
    assert csv_spheroid_effective_radius_um(row) == pytest.approx(5.64)


def test_csv_scalar_metric_rejects_missing_or_non_numeric_column():
    with pytest.raises(MetricRegistryError, match="missing"):
        csv_scalar_metric({"Area_um2": 100.0}, "A_over_A0")
    with pytest.raises(MetricRegistryError, match="numeric"):
        csv_scalar_metric({"Area_um2": "not_numeric"}, "Area_um2")


def test_registry_duplicate_register_raises():
    registry = default_registry()
    duplicate = RegisteredMetric(
        name="projected_area",
        artifact_kind=ArtifactKind.BOUNDARY_CONTOURS,
        unit="um2",
        fn=boundary_projected_area_um2,
    )
    with pytest.raises(MetricRegistryError):
        registry.register(duplicate)


def test_pure_functions_revalidate_direct_construction_boundary():
    raw = np.array(
        [
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [0.0, 1.0],
        ]
    )
    boundary = MeasurementBoundary(
        vertices_xy_um=raw,
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    assert not boundary._validated
    area = boundary_projected_area_um2(boundary)
    perimeter = boundary_perimeter_um(boundary)
    assert area == pytest.approx(2.0)
    assert perimeter == pytest.approx(6.0)
    assert boundary._validated


def test_segmentation_mask_artifact_metric_returns_same_value():
    registry = default_registry()
    boundary = _square_boundary(side_um=2.0)
    seg_metric = registry.get("projected_area", ArtifactKind.SEGMENTATION_MASK)
    contour_metric = registry.get("projected_area", ArtifactKind.BOUNDARY_CONTOURS)
    assert seg_metric(boundary) == pytest.approx(contour_metric(boundary))


def test_empty_registry_starts_empty():
    registry = MetricRegistry()
    assert registry.keys() == ()
    with pytest.raises(MetricRegistryError):
        registry.get("anything", ArtifactKind.SEGMENTATION_MASK)


def test_registry_register_rejects_non_artifact_kind():
    registry = MetricRegistry()
    fake = RegisteredMetric.__new__(RegisteredMetric)
    object.__setattr__(fake, "name", "x")
    object.__setattr__(fake, "artifact_kind", "boundary_contours")
    object.__setattr__(fake, "unit", "um")
    object.__setattr__(fake, "fn", boundary_perimeter_um)
    object.__setattr__(fake, "description", "")
    with pytest.raises(MetricRegistryError, match="ArtifactKind"):
        registry.register(fake)


def test_registry_register_rejects_empty_name():
    registry = MetricRegistry()
    bad = RegisteredMetric(
        name="",
        artifact_kind=ArtifactKind.BOUNDARY_CONTOURS,
        unit="um",
        fn=boundary_perimeter_um,
    )
    with pytest.raises(MetricRegistryError, match="name"):
        registry.register(bad)


def test_registry_register_rejects_empty_unit():
    registry = MetricRegistry()
    bad = RegisteredMetric(
        name="metric",
        artifact_kind=ArtifactKind.BOUNDARY_CONTOURS,
        unit="   ",
        fn=boundary_perimeter_um,
    )
    with pytest.raises(MetricRegistryError, match="unit"):
        registry.register(bad)


def test_registry_register_rejects_non_callable_fn():
    registry = MetricRegistry()
    bad = RegisteredMetric(
        name="metric",
        artifact_kind=ArtifactKind.BOUNDARY_CONTOURS,
        unit="um",
        fn="not_a_function",  # type: ignore[arg-type]
    )
    with pytest.raises(MetricRegistryError, match="callable"):
        registry.register(bad)

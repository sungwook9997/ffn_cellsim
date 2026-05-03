from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from acs.v2.measurement_boundary import (
    MeasurementBoundary,
    MeasurementBoundaryError,
)


def assert_failure_kind(exc_info: pytest.ExceptionInfo, failure_kind: str) -> None:
    assert exc_info.value.failure_kind == failure_kind


def test_world_y_up_ccw_square_metrics_and_read_only_vertices():
    boundary = MeasurementBoundary.from_array(
        np.array(
            [
                [0.0, 0.0],
                [2.0, 0.0],
                [2.0, 1.0],
                [0.0, 1.0],
            ]
        ),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )

    assert boundary.vertices_xy_um.dtype == np.float64
    assert not boundary.vertices_xy_um.flags.writeable
    assert boundary.signed_area_um2() == pytest.approx(2.0)
    assert boundary.projected_area_um2() == pytest.approx(2.0)
    assert boundary.perimeter_um() == pytest.approx(6.0)
    assert boundary.export_in_source_convention() == pytest.approx(boundary.vertices_xy_um)
    assert not boundary.orientation_was_reversed
    assert not boundary.coord_was_flipped

    with pytest.raises(dataclasses.FrozenInstanceError):
        boundary.source_modality = "changed"
    with pytest.raises(ValueError):
        boundary.vertices_xy_um[0, 0] = 3.0


def test_image_y_down_trace_rejects_after_flip_without_auto_orient():
    raw = np.array(
        [
            [0.0, 0.0],
            [10.0, 0.0],
            [10.0, 10.0],
            [0.0, 10.0],
        ]
    )

    with pytest.raises(MeasurementBoundaryError) as exc_info:
        MeasurementBoundary.from_array(
            raw,
            coordinate_convention="image_um_y_down",
            source_modality="confocal_segmentation",
        )

    assert_failure_kind(exc_info, "wrong_orientation")


def test_image_y_down_trace_auto_orients_after_flip():
    raw = np.array(
        [
            [0.0, 0.0],
            [10.0, 0.0],
            [10.0, 10.0],
            [0.0, 10.0],
        ]
    )

    boundary = MeasurementBoundary.from_array(
        raw,
        coordinate_convention="image_um_y_down",
        source_modality="confocal_segmentation",
        auto_orient=True,
    )

    assert boundary.source_signed_area_raw == pytest.approx(100.0)
    assert boundary.coord_was_flipped
    assert boundary.orientation_was_reversed
    assert boundary.vertices_xy_um == pytest.approx(
        np.array(
            [
                [0.0, -10.0],
                [10.0, -10.0],
                [10.0, -0.0],
                [0.0, -0.0],
            ]
        )
    )
    assert boundary.signed_area_um2() == pytest.approx(100.0)
    assert boundary.export_in_source_convention() == pytest.approx(raw[::-1])


def test_image_y_down_complementary_trace_flips_without_reversal():
    raw = np.array(
        [
            [0.0, 0.0],
            [0.0, 10.0],
            [10.0, 10.0],
            [10.0, 0.0],
        ]
    )

    boundary = MeasurementBoundary.from_array(
        raw,
        coordinate_convention="image_um_y_down",
        source_modality="confocal_segmentation",
    )

    assert boundary.source_signed_area_raw == pytest.approx(-100.0)
    assert boundary.coord_was_flipped
    assert not boundary.orientation_was_reversed
    assert boundary.signed_area_um2() == pytest.approx(100.0)
    assert boundary.export_in_source_convention() == pytest.approx(raw)


def test_world_cw_input_requires_auto_orient():
    raw = np.array(
        [
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ]
    )

    with pytest.raises(MeasurementBoundaryError) as exc_info:
        MeasurementBoundary.from_array(
            raw,
            coordinate_convention="world_um_y_up",
            source_modality="manual_outline",
        )
    assert_failure_kind(exc_info, "wrong_orientation")

    boundary = MeasurementBoundary.from_array(
        raw,
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
        auto_orient=True,
    )
    assert boundary.orientation_was_reversed
    assert boundary.signed_area_um2() == pytest.approx(1.0)


def test_float32_cast_records_input_dtype():
    raw = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
        ],
        dtype=np.float32,
    )

    boundary = MeasurementBoundary.from_array(
        raw,
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )

    assert boundary.input_dtype == "float32"
    assert boundary.vertices_xy_um.dtype == np.float64


def test_direct_construction_validates_and_promotes_cache():
    boundary = MeasurementBoundary(
        vertices_xy_um=np.array(
            [
                [0.0, 0.0],
                [3.0, 0.0],
                [2.0, 1.0],
                [0.0, 1.0],
            ]
        ),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )

    assert not boundary._validated
    boundary.validate()
    assert boundary._validated
    assert boundary.vertices_xy_um.dtype == np.float64
    assert not boundary.vertices_xy_um.flags.writeable
    boundary.validate()


def test_irregular_polygon_validates():
    boundary = MeasurementBoundary.from_array(
        [
            [0.0, 0.0],
            [2.0, 0.0],
            [2.4, 0.7],
            [1.1, 1.8],
            [-0.2, 0.8],
        ],
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )

    assert boundary.projected_area_um2() > 0.0
    assert boundary.perimeter_um() > 0.0


@pytest.mark.parametrize(
    ("vertices", "failure_kind"),
    [
        ([[0.0, 0.0], [1.0, 0.0]], "min_vertices"),
        ([[0.0, 0.0], [1.0, 0.0], [np.nan, 1.0]], "non_finite"),
        ([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]], "degenerate_area"),
        ([[0.0, 0.0], [0.0, 1.0], [1.0, 1.0]], "wrong_orientation"),
        ([[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], "duplicate_vertex"),
        ([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 0.0]], "vertex_pinch"),
        ([[0.0, 0.0], [2.0, 2.0], [0.0, 2.0], [2.0, 0.0]], "self_intersection"),
        ([[0.0, 0.0], [0.0005, 0.0], [0.0, 1.0]], "edge_too_short"),
    ],
)
def test_rejection_failure_kinds(vertices, failure_kind):
    with pytest.raises(MeasurementBoundaryError) as exc_info:
        MeasurementBoundary.from_array(
            vertices,
            coordinate_convention="world_um_y_up",
            source_modality="manual_outline",
        )

    assert_failure_kind(exc_info, failure_kind)


def test_bad_shape_rejected():
    with pytest.raises(MeasurementBoundaryError) as exc_info:
        MeasurementBoundary.from_array(
            [0.0, 1.0, 2.0],
            coordinate_convention="world_um_y_up",
            source_modality="manual_outline",
        )
    assert_failure_kind(exc_info, "shape")


def test_bad_convention_rejected():
    with pytest.raises(MeasurementBoundaryError) as exc_info:
        MeasurementBoundary.from_array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            coordinate_convention="pixels",  # type: ignore[arg-type]
            source_modality="manual_outline",
        )
    assert_failure_kind(exc_info, "convention")


def test_non_positive_min_edge_rejected():
    with pytest.raises(MeasurementBoundaryError) as exc_info:
        MeasurementBoundary.from_array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            coordinate_convention="world_um_y_up",
            source_modality="manual_outline",
            min_edge_um=0.0,
        )
    assert_failure_kind(exc_info, "min_edge_param")


def test_tiny_large_near_collinear_and_exact_min_edge_cases():
    tiny = MeasurementBoundary.from_array(
        [[0.0, 0.0], [1e-4, 0.0], [0.0, 1e-4]],
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
        min_edge_um=1e-5,
    )
    assert tiny.projected_area_um2() == pytest.approx(5e-9)

    large = MeasurementBoundary.from_array(
        [[0.0, 0.0], [1e6, 0.0], [1e6, 1e6], [0.0, 1e6]],
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    assert large.projected_area_um2() == pytest.approx(1e12)

    near_collinear = MeasurementBoundary.from_array(
        [[0.0, 0.0], [1.0, 1e-9], [2.0, 0.0], [1.0, 1.0]],
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    assert near_collinear.projected_area_um2() > 0.0

    with_collinear_edge_vertex = MeasurementBoundary.from_array(
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]],
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    assert with_collinear_edge_vertex.projected_area_um2() == pytest.approx(2.0)

    exact_min = MeasurementBoundary.from_array(
        [[0.0, 0.0], [0.001, 0.0], [0.0, 0.001]],
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
        min_edge_um=0.001,
    )
    assert exact_min.perimeter_um() > 0.0

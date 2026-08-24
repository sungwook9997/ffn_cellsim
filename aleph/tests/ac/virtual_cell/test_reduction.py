from __future__ import annotations

import numpy as np
import pytest

from aleph.virtual_cell.reduction import (
    ArrayReduction,
    MatrixSVDApproximation,
    evaluate_projection_errors,
    select_array_representation,
    svd_reduce,
)


def test_vector_routes_to_dense() -> None:
    decision = select_array_representation(np.arange(8.0), relative_tolerance=0.0)
    assert decision.kind == "dense"
    assert decision.fallback_reason == "vector-has-no-matrix-split"


def test_rank_one_matrix_routes_to_svd() -> None:
    matrix = np.multiply.outer(np.arange(1.0, 11.0), np.arange(1.0, 9.0))
    decision = select_array_representation(matrix, relative_tolerance=1e-12)
    assert decision.kind == "svd"
    assert decision.measured_relative_error < 1e-12
    np.testing.assert_allclose(decision.reconstruct(), matrix, rtol=1e-12, atol=1e-12)


def test_random_small_matrix_falls_back_if_svd_is_not_smaller() -> None:
    matrix = np.random.default_rng(3).normal(size=(5, 5))
    decision = select_array_representation(matrix, relative_tolerance=0.0)
    assert decision.kind == "dense"
    assert decision.fallback_reason == "svd-storage-not-smaller-than-dense"


def test_rank_one_tensor_routes_to_tensor_train() -> None:
    tensor = np.einsum(
        "i,j,k,l->ijkl",
        np.arange(1.0, 5.0),
        np.arange(1.0, 6.0),
        np.arange(1.0, 4.0),
        np.arange(1.0, 7.0),
    )
    decision = select_array_representation(tensor, relative_tolerance=1e-12)
    assert decision.kind == "tensor-train"
    assert decision.compression_ratio > 1.0


def test_random_high_rank_tensor_falls_back_instead_of_claiming_compression() -> None:
    tensor = np.random.default_rng(7).normal(size=(4, 4, 4, 4))
    decision = select_array_representation(tensor, relative_tolerance=0.0)
    assert decision.kind == "dense"
    assert decision.fallback_reason == "tt-storage-not-smaller-than-dense"


def test_rank_cap_that_misses_tolerance_falls_back() -> None:
    matrix = np.eye(8)
    decision = select_array_representation(
        matrix,
        relative_tolerance=1e-6,
        max_rank=1,
    )
    assert decision.kind == "dense"
    assert decision.fallback_reason == "svd-error-exceeds-requested-tolerance"


def test_svd_reports_actual_error() -> None:
    matrix = np.diag([4.0, 2.0, 1.0])
    approximation = svd_reduce(matrix, relative_tolerance=0.0, max_rank=1)
    expected = np.sqrt(5.0) / np.sqrt(21.0)
    assert approximation.relative_error == pytest.approx(expected)


def test_projection_error_is_measured_separately_from_global_error() -> None:
    reference = np.ones((100, 100))
    reconstructed = reference.copy()
    reconstructed[0, 0] = 0.0
    errors = evaluate_projection_errors(
        reference,
        reconstructed,
        {
            "global_mean": np.mean,
            "critical_pixel": lambda array: array[0, 0],
        },
    )
    assert errors["global_mean"] == pytest.approx(1e-4)
    assert errors["critical_pixel"] == pytest.approx(1.0)


def test_probability_reduction_falls_back_if_svd_creates_negative_mass() -> None:
    mass = np.random.default_rng(0).exponential(size=(20, 20))
    mass /= np.sum(mass)
    decision = select_array_representation(
        mass,
        relative_tolerance=0.3,
        source_kind="probability-mass",
    )
    assert decision.kind == "dense"
    assert decision.fallback_reason == "reduced-probability-violates-mass-or-nonnegativity"


def test_probability_source_rejects_unnormalized_input_instead_of_repairing() -> None:
    with pytest.raises(ValueError, match="not normalized"):
        select_array_representation(
            np.ones((4, 4)),
            relative_tolerance=0.1,
            source_kind="probability-mass",
        )


def test_large_finite_matrix_does_not_turn_overflow_into_false_compression() -> None:
    matrix = np.eye(20) * 1e308
    decision = select_array_representation(matrix, relative_tolerance=0.01)
    assert decision.kind == "dense"
    assert np.isfinite(decision.measured_relative_error)


def test_declared_projection_failure_overrides_small_global_error() -> None:
    matrix = np.ones((100, 100))
    matrix[0, 0] = 0.05
    decision = select_array_representation(
        matrix,
        relative_tolerance=1e-2,
        projections={
            "global_mean": np.mean,
            "rare_local_value": lambda array: array[0, 0],
        },
        projection_tolerances={
            "global_mean": 1e-3,
            "rare_local_value": 0.1,
        },
    )
    assert decision.kind == "dense"
    assert decision.fallback_reason == "declared-projection-error-exceeds-tolerance"


def test_projection_names_cannot_collide_after_normalization() -> None:
    with pytest.raises(ValueError, match="collide"):
        evaluate_projection_errors(
            np.ones((2, 2)),
            np.ones((2, 2)),
            {"mean": np.mean, " mean ": np.max},
        )


def test_reduction_dataclasses_reject_inconsistent_fabricated_metrics() -> None:
    with pytest.raises(ValueError, match="factor shapes"):
        MatrixSVDApproximation(
            u=np.ones((2, 1)),
            singular_values=np.ones(1),
            vh=np.ones((1, 3)),
            original_shape=(2, 2),
            relative_error=0.0,
        )
    with pytest.raises(ValueError, match="positive integer"):
        ArrayReduction(
            kind="dense",
            source_kind="array",
            original_shape=(2, 2),
            requested_tolerance=0.0,
            measured_relative_error=0.0,
            stored_parameter_count=0,
            dense_parameter_count=4,
            fallback_reason="test",
            dense=np.ones((2, 2)),
        )

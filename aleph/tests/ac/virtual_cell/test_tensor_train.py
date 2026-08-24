"""Positive and adversarial controls for the CPU TT-SVD oracle."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.virtual_cell.tensor_train import (
    TensorTrainApproximation,
    axis_order_rank_report,
    matricization_ranks,
    tt_svd,
)


def test_rank_one_tensor_reconstructs_with_unit_internal_ranks() -> None:
    tensor = np.einsum(
        "i,j,k->ijk",
        np.array([1.0, 2.0, 3.0]),
        np.array([-1.0, 4.0]),
        np.array([0.5, 1.5, 2.5, 3.5]),
    )
    fitted = tt_svd(tensor, relative_tolerance=1e-12)
    np.testing.assert_allclose(fitted.reconstruct(), tensor, rtol=1e-11, atol=1e-11)
    assert fitted.ranks == (1, 1, 1, 1)
    assert fitted.compression_ratio > 1.0


def test_sum_of_two_products_needs_no_more_than_rank_two() -> None:
    a = np.einsum("i,j,k->ijk", [1.0, 2.0], [1.0, 3.0, 5.0], [2.0, 4.0])
    b = np.einsum("i,j,k->ijk", [2.0, -1.0], [0.5, 2.0, -3.0], [1.0, -2.0])
    tensor = a + b
    fitted = tt_svd(tensor, relative_tolerance=1e-12)
    np.testing.assert_allclose(fitted.reconstruct(), tensor, rtol=1e-10, atol=1e-10)
    assert max(fitted.ranks) <= 2


def test_high_rank_negative_control_reports_large_rank_or_large_capped_error() -> None:
    tensor = np.random.default_rng(17).normal(size=(5, 5, 5, 5))
    exact_ranks = matricization_ranks(tensor)
    assert max(exact_ranks) >= 20
    capped = tt_svd(tensor, max_rank=2)
    assert capped.relative_error > 0.5
    assert capped.compression_ratio > 1.0


def test_axis_order_is_measured_not_assumed_invariant() -> None:
    left = np.arange(6.0).reshape(2, 3) + 1.0
    right = np.arange(20.0).reshape(4, 5) + 1.0
    tensor = np.einsum("ij,kl->ijkl", left, right)
    report = axis_order_rank_report(tensor)
    assert report[(0, 1, 2, 3)] != report[(0, 2, 1, 3)]


def test_requested_tolerance_is_met_without_rank_cap() -> None:
    tensor = np.random.default_rng(3).normal(size=(4, 5, 6))
    fitted = tt_svd(tensor, relative_tolerance=0.2)
    assert fitted.relative_error <= 0.2 + 1e-12


def test_axis_permutation_round_trips_to_original_order() -> None:
    tensor = np.random.default_rng(5).normal(size=(2, 3, 4))
    fitted = tt_svd(tensor, axis_order=(2, 0, 1))
    np.testing.assert_allclose(fitted.reconstruct(), tensor, rtol=1e-12, atol=1e-12)
    assert fitted.reconstruct(original_axis_order=False).shape == (4, 2, 3)


def test_large_finite_tensor_uses_scale_safe_error_arithmetic() -> None:
    tensor = np.zeros((3, 3, 3))
    tensor[0, 0, 0] = 1e308
    tensor[1, 1, 1] = 5e307
    fitted = tt_svd(tensor, relative_tolerance=0.0)
    assert np.isfinite(fitted.relative_error)
    np.testing.assert_allclose(fitted.reconstruct(), tensor, rtol=1e-14)


def test_tensor_train_dataclass_rejects_inconsistent_core_shapes() -> None:
    with pytest.raises(ValueError, match="core 0 shape"):
        TensorTrainApproximation(
            cores=(np.ones((1, 3, 1)), np.ones((1, 2, 1))),
            original_shape=(2, 2),
            ranks=(1, 1, 1),
            relative_error=0.0,
            axis_order=(0, 1),
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"relative_tolerance": -0.1},
        {"relative_tolerance": 1.0},
        {"max_rank": 0},
        {"axis_order": (0, 0, 1)},
    ],
)
def test_invalid_tensor_train_controls_are_rejected(kwargs) -> None:
    with pytest.raises(ValueError):
        tt_svd(np.ones((2, 3, 4)), **kwargs)

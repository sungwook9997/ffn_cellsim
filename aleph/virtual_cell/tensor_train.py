"""Small deterministic TT-SVD oracle for the tensor-network research lane.

This is not a production tensor backend.  It establishes the contracts a future TT-cross/ALS backend
must satisfy on CPU: reconstruction, measured truncation error, storage cost, exact split ranks, and
honest failure when a tensor is not compressible.  It uses only NumPy and never touches simulation
runtime state or a GPU.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import permutations
from typing import Any

import numpy as np

__all__ = [
    "TensorTrainApproximation",
    "axis_order_rank_report",
    "matricization_ranks",
    "tt_svd",
]


def _as_tensor(values: Any) -> np.ndarray:
    tensor = np.asarray(values, dtype=np.float64)
    if tensor.ndim < 2:
        raise ValueError(f"tensor must have at least two axes; got shape {tensor.shape}")
    if tensor.size == 0:
        raise ValueError("tensor must not be empty")
    if not np.all(np.isfinite(tensor)):
        raise ValueError("tensor must contain only finite values")
    return np.array(tensor, copy=True)


def _stable_relative_error(reference: np.ndarray, candidate: np.ndarray) -> float:
    scale = float(np.max(np.abs(reference)))
    if scale == 0.0:
        return float(np.linalg.norm(candidate.ravel()))
    denominator = float(np.linalg.norm((reference / scale).ravel()))
    numerator = float(np.linalg.norm(((reference - candidate) / scale).ravel()))
    return numerator / denominator


def _truncation_rank(
    singular_values: np.ndarray,
    *,
    squared_error_budget: float,
    max_rank: int | None,
) -> int:
    available = int(singular_values.size)
    rank = available
    if squared_error_budget > 0.0:
        tail_sq = np.cumsum(np.square(singular_values[::-1]))[::-1]
        for candidate in range(1, available + 1):
            discarded = float(tail_sq[candidate]) if candidate < available else 0.0
            if discarded <= squared_error_budget:
                rank = candidate
                break
    if max_rank is not None:
        rank = min(rank, max_rank)
    return max(1, rank)


@dataclass(frozen=True, slots=True)
class TensorTrainApproximation:
    """A TT-SVD result with immutable cores and measured reconstruction error."""

    cores: tuple[np.ndarray, ...]
    original_shape: tuple[int, ...]
    ranks: tuple[int, ...]
    relative_error: float
    axis_order: tuple[int, ...]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.original_shape, tuple)
            or len(self.original_shape) < 2
            or any(
                isinstance(size, bool) or not isinstance(size, int) or size < 1
                for size in self.original_shape
            )
        ):
            raise ValueError("original_shape must contain at least two positive integer dimensions")
        ndim = len(self.original_shape)
        if tuple(sorted(self.axis_order)) != tuple(range(ndim)):
            raise ValueError("axis_order must be a permutation of original_shape axes")
        if (
            not isinstance(self.ranks, tuple)
            or len(self.ranks) != ndim + 1
            or self.ranks[0] != 1
            or self.ranks[-1] != 1
            or any(
                isinstance(rank, bool) or not isinstance(rank, int) or rank < 1
                for rank in self.ranks
            )
        ):
            raise ValueError("ranks must contain positive TT bonds with unit endpoints")
        if not math.isfinite(self.relative_error) or self.relative_error < 0.0:
            raise ValueError("relative_error must be finite and nonnegative")
        if len(self.cores) != ndim:
            raise ValueError("one TT core is required per tensor axis")
        fitted_shape = tuple(self.original_shape[index] for index in self.axis_order)
        frozen: list[np.ndarray] = []
        for axis, core in enumerate(self.cores):
            copied = np.array(core, dtype=np.float64, copy=True)
            expected_shape = (self.ranks[axis], fitted_shape[axis], self.ranks[axis + 1])
            if copied.shape != expected_shape:
                raise ValueError(
                    f"TT core {axis} shape {copied.shape} does not match {expected_shape}"
                )
            if not np.all(np.isfinite(copied)):
                raise ValueError("TT cores must be finite")
            copied.setflags(write=False)
            frozen.append(copied)
        object.__setattr__(self, "cores", tuple(frozen))

    @property
    def storage_size(self) -> int:
        """Number of float64 values stored across all TT cores."""
        return int(sum(core.size for core in self.cores))

    @property
    def dense_size(self) -> int:
        """Number of values in the dense tensor."""
        return int(math.prod(self.original_shape))

    @property
    def compression_ratio(self) -> float:
        """Dense values divided by TT values; above one means actual compression."""
        return self.dense_size / self.storage_size

    def reconstruct(self, *, original_axis_order: bool = True) -> np.ndarray:
        """Reconstruct the dense tensor, optionally undoing the fitted axis permutation."""
        result = self.cores[0][0, :, :]
        for core in self.cores[1:]:
            result = np.tensordot(result, core, axes=([-1], [0]))
        result = np.squeeze(result, axis=-1)
        if not original_axis_order:
            return result
        inverse = tuple(int(index) for index in np.argsort(self.axis_order))
        return np.transpose(result, axes=inverse)


def tt_svd(
    values: Any,
    *,
    relative_tolerance: float = 0.0,
    max_rank: int | None = None,
    axis_order: tuple[int, ...] | None = None,
) -> TensorTrainApproximation:
    """Fit a tensor train by sequential SVD and report the actual dense reconstruction error.

    Args:
        values: Dense finite tensor with at least two axes.
        relative_tolerance: Requested relative Frobenius error in ``[0, 1)``.  The squared budget is
            divided across TT splits; a finite ``max_rank`` may force a larger measured error.
        max_rank: Optional positive cap on every internal TT rank.
        axis_order: Optional permutation fitted before TT-SVD.  The result records it and reconstructs
            in the caller's original order by default.
    """
    tensor = _as_tensor(values)
    if not (isinstance(relative_tolerance, (int, float)) and 0.0 <= relative_tolerance < 1.0):
        raise ValueError("relative_tolerance must lie in [0, 1)")
    if max_rank is not None and (
        isinstance(max_rank, bool) or not isinstance(max_rank, int) or max_rank < 1
    ):
        raise ValueError("max_rank must be a positive integer when supplied")

    ndim = tensor.ndim
    if axis_order is None:
        order = tuple(range(ndim))
    else:
        order = tuple(axis_order)
        if tuple(sorted(order)) != tuple(range(ndim)):
            raise ValueError(f"axis_order must be a permutation of 0..{ndim - 1}")
    working = np.transpose(tensor, axes=order)
    shape = working.shape
    scale = float(np.max(np.abs(working)))
    unfolding = working / scale if scale > 0.0 else working
    scaled_norm = float(np.linalg.norm(unfolding.ravel()))
    total_budget_sq = (float(relative_tolerance) * scaled_norm) ** 2
    split_budget_sq = total_budget_sq / (ndim - 1)

    cores: list[np.ndarray] = []
    ranks = [1]
    left_rank = 1
    for axis in range(ndim - 1):
        unfolding = unfolding.reshape(left_rank * shape[axis], -1)
        u, singular_values, vh = np.linalg.svd(unfolding, full_matrices=False)
        rank = _truncation_rank(
            singular_values,
            squared_error_budget=split_budget_sq,
            max_rank=max_rank,
        )
        cores.append(u[:, :rank].reshape(left_rank, shape[axis], rank))
        unfolding = singular_values[:rank, None] * vh[:rank, :]
        ranks.append(rank)
        left_rank = rank
    cores.append((unfolding * scale).reshape(left_rank, shape[-1], 1))
    ranks.append(1)

    provisional = TensorTrainApproximation(
        cores=tuple(cores),
        original_shape=tuple(tensor.shape),
        ranks=tuple(ranks),
        relative_error=0.0,
        axis_order=order,
    )
    reconstructed = provisional.reconstruct()
    relative_error = _stable_relative_error(tensor, reconstructed)
    return TensorTrainApproximation(
        cores=provisional.cores,
        original_shape=provisional.original_shape,
        ranks=provisional.ranks,
        relative_error=relative_error,
        axis_order=provisional.axis_order,
    )


def matricization_ranks(values: Any, *, relative_tolerance: float = 0.0) -> tuple[int, ...]:
    """Return numerical matrix ranks at every contiguous tensor split."""
    tensor = _as_tensor(values)
    if not (isinstance(relative_tolerance, (int, float)) and 0.0 <= relative_tolerance < 1.0):
        raise ValueError("relative_tolerance must lie in [0, 1)")
    scale = float(np.max(np.abs(tensor)))
    scaled_tensor = tensor / scale if scale > 0.0 else tensor
    norm = float(np.linalg.norm(scaled_tensor.ravel()))
    budget_sq = (float(relative_tolerance) * norm) ** 2
    ranks: list[int] = []
    for split in range(1, tensor.ndim):
        matrix = scaled_tensor.reshape(math.prod(tensor.shape[:split]), -1)
        singular_values = np.linalg.svd(matrix, full_matrices=False, compute_uv=False)
        ranks.append(
            _truncation_rank(
                singular_values,
                squared_error_budget=budget_sq,
                max_rank=None,
            )
        )
    return tuple(ranks)


def axis_order_rank_report(
    values: Any,
    *,
    relative_tolerance: float = 0.0,
    max_orders: int = 120,
) -> dict[tuple[int, ...], tuple[int, ...]]:
    """Measure split-rank dependence on axis ordering instead of calling rank a physical invariant."""
    tensor = _as_tensor(values)
    if max_orders < 1:
        raise ValueError("max_orders must be positive")
    report: dict[tuple[int, ...], tuple[int, ...]] = {}
    for index, order in enumerate(permutations(range(tensor.ndim))):
        if index >= max_orders:
            break
        report[order] = matricization_ranks(
            np.transpose(tensor, axes=order),
            relative_tolerance=relative_tolerance,
        )
    return report

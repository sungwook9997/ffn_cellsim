"""CPU reference selection between dense, SVD, and tensor-train array representations."""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Literal

import numpy as np

from aleph.virtual_cell.probability import audit_probability_mass
from aleph.virtual_cell.tensor_train import TensorTrainApproximation, tt_svd

__all__ = [
    "ArrayReduction",
    "MatrixSVDApproximation",
    "evaluate_projection_errors",
    "select_array_representation",
    "svd_reduce",
]


def _finite_array(values: Any, *, minimum_ndim: int = 1) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim < minimum_ndim or array.size == 0:
        raise ValueError(f"array must be non-empty with at least {minimum_ndim} axes")
    if not np.all(np.isfinite(array)):
        raise ValueError("array must contain only finite values")
    return np.array(array, copy=True)


def _relative_error(reference: np.ndarray, candidate: np.ndarray) -> float:
    scale = float(np.max(np.abs(reference)))
    if scale == 0.0:
        return float(np.linalg.norm(candidate.ravel()))
    scaled_reference = reference / scale
    scaled_difference = (reference - candidate) / scale
    denominator = float(np.linalg.norm(scaled_reference.ravel()))
    return float(np.linalg.norm(scaled_difference.ravel())) / denominator


@dataclass(frozen=True, slots=True)
class MatrixSVDApproximation:
    """Truncated matrix SVD with immutable factors and measured error."""

    u: np.ndarray
    singular_values: np.ndarray
    vh: np.ndarray
    original_shape: tuple[int, int]
    relative_error: float

    def __post_init__(self) -> None:
        if (
            not isinstance(self.original_shape, tuple)
            or len(self.original_shape) != 2
            or any(
                isinstance(size, bool) or not isinstance(size, int) or size < 1
                for size in self.original_shape
            )
        ):
            raise ValueError("original_shape must contain two positive integer dimensions")
        if not math.isfinite(self.relative_error) or self.relative_error < 0.0:
            raise ValueError("relative_error must be finite and nonnegative")
        frozen: list[np.ndarray] = []
        for value in (self.u, self.singular_values, self.vh):
            copied = np.array(value, dtype=np.float64, copy=True)
            if not np.all(np.isfinite(copied)):
                raise ValueError("SVD factors must be finite")
            copied.setflags(write=False)
            frozen.append(copied)
        u, singular_values, vh = frozen
        if (
            u.ndim != 2
            or singular_values.ndim != 1
            or vh.ndim != 2
            or singular_values.size < 1
            or u.shape != (self.original_shape[0], singular_values.size)
            or vh.shape != (singular_values.size, self.original_shape[1])
        ):
            raise ValueError("SVD factor shapes do not match original_shape and rank")
        if np.any(singular_values < 0.0):
            raise ValueError("singular values must be nonnegative")
        object.__setattr__(self, "u", frozen[0])
        object.__setattr__(self, "singular_values", frozen[1])
        object.__setattr__(self, "vh", frozen[2])

    @property
    def rank(self) -> int:
        return int(self.singular_values.size)

    @property
    def storage_size(self) -> int:
        return int(self.u.size + self.singular_values.size + self.vh.size)

    @property
    def dense_size(self) -> int:
        return int(math.prod(self.original_shape))

    def reconstruct(self) -> np.ndarray:
        return (self.u * self.singular_values[None, :]) @ self.vh


def svd_reduce(
    values: Any,
    *,
    relative_tolerance: float = 0.0,
    max_rank: int | None = None,
) -> MatrixSVDApproximation:
    """Fit a matrix SVD and report the actual reconstruction error."""
    matrix = _finite_array(values, minimum_ndim=2)
    if matrix.ndim != 2:
        raise ValueError(f"svd_reduce requires a matrix, got shape {matrix.shape}")
    if not (isinstance(relative_tolerance, (int, float)) and 0.0 <= relative_tolerance < 1.0):
        raise ValueError("relative_tolerance must lie in [0, 1)")
    if max_rank is not None and (
        isinstance(max_rank, bool) or not isinstance(max_rank, int) or max_rank < 1
    ):
        raise ValueError("max_rank must be a positive integer when supplied")

    scale = float(np.max(np.abs(matrix)))
    scaled_matrix = matrix / scale if scale > 0.0 else matrix
    u, scaled_singular_values, vh = np.linalg.svd(scaled_matrix, full_matrices=False)
    norm_sq = float(np.sum(np.square(scaled_singular_values)))
    budget_sq = float(relative_tolerance) ** 2 * norm_sq
    rank = scaled_singular_values.size
    tail_sq = np.cumsum(np.square(scaled_singular_values[::-1]))[::-1]
    for candidate in range(1, scaled_singular_values.size + 1):
        discarded = float(tail_sq[candidate]) if candidate < scaled_singular_values.size else 0.0
        if discarded <= budget_sq:
            rank = candidate
            break
    if max_rank is not None:
        rank = min(rank, max_rank)
    approximation = MatrixSVDApproximation(
        u=u[:, :rank],
        singular_values=scaled_singular_values[:rank] * scale,
        vh=vh[:rank, :],
        original_shape=(int(matrix.shape[0]), int(matrix.shape[1])),
        relative_error=0.0,
    )
    return MatrixSVDApproximation(
        u=approximation.u,
        singular_values=approximation.singular_values,
        vh=approximation.vh,
        original_shape=approximation.original_shape,
        relative_error=_relative_error(matrix, approximation.reconstruct()),
    )


@dataclass(frozen=True, slots=True)
class ArrayReduction:
    """Fail-closed representation decision for a dense oracle array."""

    kind: Literal["dense", "svd", "tensor-train"]
    source_kind: Literal["array", "probability-mass"]
    original_shape: tuple[int, ...]
    requested_tolerance: float
    measured_relative_error: float
    stored_parameter_count: int
    dense_parameter_count: int
    fallback_reason: str | None
    dense: np.ndarray | None = None
    svd: MatrixSVDApproximation | None = None
    tensor_train: TensorTrainApproximation | None = None

    def __post_init__(self) -> None:
        if self.kind not in {"dense", "svd", "tensor-train"}:
            raise ValueError("kind must be dense, svd, or tensor-train")
        if self.source_kind not in {"array", "probability-mass"}:
            raise ValueError("source_kind must be array or probability-mass")
        if (
            not isinstance(self.original_shape, tuple)
            or not self.original_shape
            or any(
                isinstance(size, bool) or not isinstance(size, int) or size < 1
                for size in self.original_shape
            )
        ):
            raise ValueError("original_shape must contain positive integer dimensions")
        if not (math.isfinite(self.requested_tolerance) and 0.0 <= self.requested_tolerance < 1.0):
            raise ValueError("requested_tolerance must lie in [0, 1)")
        if not math.isfinite(self.measured_relative_error) or self.measured_relative_error < 0.0:
            raise ValueError("measured_relative_error must be finite and nonnegative")
        for field_name in ("stored_parameter_count", "dense_parameter_count"):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{field_name} must be a positive integer")
        if self.dense_parameter_count != math.prod(self.original_shape):
            raise ValueError("dense_parameter_count must equal the original shape product")
        payload_count = sum(
            payload is not None for payload in (self.dense, self.svd, self.tensor_train)
        )
        if payload_count != 1:
            raise ValueError("exactly one representation payload is required")
        if self.kind == "dense" and self.dense is None:
            raise ValueError("dense decision needs a dense payload")
        if self.kind == "svd" and self.svd is None:
            raise ValueError("svd decision needs an SVD payload")
        if self.kind == "tensor-train" and self.tensor_train is None:
            raise ValueError("tensor-train decision needs a TT payload")
        if self.dense is not None:
            copied = np.array(self.dense, dtype=np.float64, copy=True)
            if copied.shape != self.original_shape or not np.all(np.isfinite(copied)):
                raise ValueError("dense payload must be finite and match original_shape")
            copied.setflags(write=False)
            object.__setattr__(self, "dense", copied)
            if self.stored_parameter_count != copied.size:
                raise ValueError("stored_parameter_count does not match dense payload")
        elif self.svd is not None:
            if self.svd.original_shape != self.original_shape:
                raise ValueError("SVD payload shape does not match original_shape")
            if self.stored_parameter_count != self.svd.storage_size:
                raise ValueError("stored_parameter_count does not match SVD payload")
            if self.measured_relative_error != self.svd.relative_error:
                raise ValueError("measured_relative_error does not match SVD payload")
        else:
            assert self.tensor_train is not None
            if self.tensor_train.original_shape != self.original_shape:
                raise ValueError("TT payload shape does not match original_shape")
            if self.stored_parameter_count != self.tensor_train.storage_size:
                raise ValueError("stored_parameter_count does not match TT payload")
            if self.measured_relative_error != self.tensor_train.relative_error:
                raise ValueError("measured_relative_error does not match TT payload")
        if self.kind == "dense" and not self.fallback_reason:
            raise ValueError("dense selection must record its fallback/route reason")
        if self.kind != "dense" and self.fallback_reason is not None:
            raise ValueError("compressed selection must not carry a fallback reason")

    @property
    def compression_ratio(self) -> float:
        return self.dense_parameter_count / self.stored_parameter_count

    def reconstruct(self) -> np.ndarray:
        if self.dense is not None:
            return np.array(self.dense, copy=True)
        if self.svd is not None:
            return self.svd.reconstruct()
        assert self.tensor_train is not None
        return self.tensor_train.reconstruct()


def _dense_decision(
    array: np.ndarray,
    *,
    tolerance: float,
    reason: str,
) -> ArrayReduction:
    return ArrayReduction(
        kind="dense",
        source_kind="array",
        original_shape=tuple(int(size) for size in array.shape),
        requested_tolerance=tolerance,
        measured_relative_error=0.0,
        stored_parameter_count=int(array.size),
        dense_parameter_count=int(array.size),
        fallback_reason=reason,
        dense=array,
    )


def select_array_representation(
    values: Any,
    *,
    relative_tolerance: float,
    max_rank: int | None = None,
    axis_order: tuple[int, ...] | None = None,
    source_kind: Literal["array", "probability-mass"] = "array",
    projections: Mapping[str, Callable[[np.ndarray], float]] | None = None,
    projection_tolerances: Mapping[str, float] | None = None,
) -> ArrayReduction:
    """Select a smaller valid representation or retain the dense array with an explicit reason."""
    array = _finite_array(values)
    if not (isinstance(relative_tolerance, (int, float)) and 0.0 <= relative_tolerance < 1.0):
        raise ValueError("relative_tolerance must lie in [0, 1)")
    tolerance = float(relative_tolerance)
    numerical_slack = max(1e-14, tolerance * 1e-10)
    if source_kind not in {"array", "probability-mass"}:
        raise ValueError("source_kind must be 'array' or 'probability-mass'")
    if source_kind == "probability-mass":
        audit_probability_mass(array)
    if (projections is None) != (projection_tolerances is None):
        raise ValueError("projections and projection_tolerances must be supplied together")

    if array.ndim == 1:
        decision = _dense_decision(
            array,
            tolerance=tolerance,
            reason="vector-has-no-matrix-split",
        )
        return _with_source_kind(decision, source_kind)
    if array.ndim == 2:
        approximation = svd_reduce(
            array,
            relative_tolerance=tolerance,
            max_rank=max_rank,
        )
        if approximation.relative_error > tolerance + numerical_slack:
            decision = _dense_decision(
                array, tolerance=tolerance, reason="svd-error-exceeds-requested-tolerance"
            )
            return _with_source_kind(decision, source_kind)
        if approximation.storage_size >= array.size:
            decision = _dense_decision(
                array, tolerance=tolerance, reason="svd-storage-not-smaller-than-dense"
            )
            return _with_source_kind(decision, source_kind)
        decision = ArrayReduction(
            kind="svd",
            source_kind=source_kind,
            original_shape=tuple(int(size) for size in array.shape),
            requested_tolerance=tolerance,
            measured_relative_error=approximation.relative_error,
            stored_parameter_count=approximation.storage_size,
            dense_parameter_count=int(array.size),
            fallback_reason=None,
            svd=approximation,
        )
        decision = _probability_safe_or_dense(array, decision)
        return _projection_safe_or_dense(
            array,
            decision,
            projections=projections,
            projection_tolerances=projection_tolerances,
        )

    approximation = tt_svd(
        array,
        relative_tolerance=tolerance,
        max_rank=max_rank,
        axis_order=axis_order,
    )
    if approximation.relative_error > tolerance + numerical_slack:
        decision = _dense_decision(
            array, tolerance=tolerance, reason="tt-error-exceeds-requested-tolerance"
        )
        return _with_source_kind(decision, source_kind)
    if approximation.storage_size >= array.size:
        decision = _dense_decision(
            array, tolerance=tolerance, reason="tt-storage-not-smaller-than-dense"
        )
        return _with_source_kind(decision, source_kind)
    decision = ArrayReduction(
        kind="tensor-train",
        source_kind=source_kind,
        original_shape=tuple(int(size) for size in array.shape),
        requested_tolerance=tolerance,
        measured_relative_error=approximation.relative_error,
        stored_parameter_count=approximation.storage_size,
        dense_parameter_count=int(array.size),
        fallback_reason=None,
        tensor_train=approximation,
    )
    decision = _probability_safe_or_dense(array, decision)
    return _projection_safe_or_dense(
        array,
        decision,
        projections=projections,
        projection_tolerances=projection_tolerances,
    )


def _with_source_kind(
    decision: ArrayReduction,
    source_kind: Literal["array", "probability-mass"],
) -> ArrayReduction:
    if decision.source_kind == source_kind:
        return decision
    return ArrayReduction(
        kind=decision.kind,
        source_kind=source_kind,
        original_shape=decision.original_shape,
        requested_tolerance=decision.requested_tolerance,
        measured_relative_error=decision.measured_relative_error,
        stored_parameter_count=decision.stored_parameter_count,
        dense_parameter_count=decision.dense_parameter_count,
        fallback_reason=decision.fallback_reason,
        dense=decision.dense,
        svd=decision.svd,
        tensor_train=decision.tensor_train,
    )


def _probability_safe_or_dense(
    reference: np.ndarray,
    decision: ArrayReduction,
) -> ArrayReduction:
    if decision.source_kind != "probability-mass" or decision.kind == "dense":
        return decision
    try:
        audit_probability_mass(decision.reconstruct())
    except ValueError:
        dense = _dense_decision(
            reference,
            tolerance=decision.requested_tolerance,
            reason="reduced-probability-violates-mass-or-nonnegativity",
        )
        return _with_source_kind(dense, "probability-mass")
    return decision


def _projection_safe_or_dense(
    reference: np.ndarray,
    decision: ArrayReduction,
    *,
    projections: Mapping[str, Callable[[np.ndarray], float]] | None,
    projection_tolerances: Mapping[str, float] | None,
) -> ArrayReduction:
    if projections is None or decision.kind == "dense":
        return decision
    assert projection_tolerances is not None
    errors = evaluate_projection_errors(reference, decision.reconstruct(), projections)
    normalized_tolerances: dict[str, float] = {}
    for raw_name, raw_tolerance in projection_tolerances.items():
        if not isinstance(raw_name, str) or not raw_name.strip():
            raise ValueError("projection tolerance names must be non-empty strings")
        name = raw_name.strip()
        if name in normalized_tolerances:
            raise ValueError("projection tolerance names must not collide after normalization")
        if (
            isinstance(raw_tolerance, bool)
            or not isinstance(raw_tolerance, (int, float))
            or not math.isfinite(raw_tolerance)
            or raw_tolerance < 0.0
        ):
            raise ValueError(f"projection tolerance {name!r} must be finite and nonnegative")
        normalized_tolerances[name] = float(raw_tolerance)
    if set(errors) != set(normalized_tolerances):
        raise ValueError("projection names and projection_tolerances must match exactly")
    if any(errors[name] > normalized_tolerances[name] for name in errors):
        dense = _dense_decision(
            reference,
            tolerance=decision.requested_tolerance,
            reason="declared-projection-error-exceeds-tolerance",
        )
        return _with_source_kind(dense, decision.source_kind)
    return decision


def evaluate_projection_errors(
    reference: Any,
    reconstructed: Any,
    projections: Mapping[str, Callable[[np.ndarray], float]],
) -> Mapping[str, float]:
    """Measure each declared scalar projection independently of global Frobenius error."""
    expected = _finite_array(reference)
    observed = _finite_array(reconstructed)
    if expected.shape != observed.shape:
        raise ValueError("reference and reconstructed arrays must have identical shape")
    if not isinstance(projections, Mapping) or not projections:
        raise ValueError("projections must be a non-empty mapping")
    errors: dict[str, float] = {}
    for raw_name, projection in projections.items():
        if not isinstance(raw_name, str) or not raw_name.strip():
            raise ValueError("projection names must be non-empty strings")
        name = raw_name.strip()
        if name in errors:
            raise ValueError("projection names must not collide after whitespace normalization")
        if not callable(projection):
            raise TypeError(f"projection {raw_name!r} must be callable")
        expected_value = float(projection(expected))
        observed_value = float(projection(observed))
        if not math.isfinite(expected_value) or not math.isfinite(observed_value):
            raise ValueError(f"projection {raw_name!r} returned a non-finite value")
        scale = abs(expected_value)
        errors[name] = (
            abs(observed_value - expected_value) / scale
            if scale > 0.0
            else abs(observed_value - expected_value)
        )
    return MappingProxyType(errors)

"""Receipt-chain to multi-index tensor-archive adapter with measured storage and decode cost.

A scalar headline ("peak tension rose 12%") is the smallest possible projection of a run, and once it
is written down the tensor it came from is usually gone.  This module keeps the tensor and demotes the
scalar to a *declared, versioned projection* of it, so a historical headline can be recomputed rather
than trusted.  Three defects it exists to make impossible:

* **A rejected attempt occupying a time index.**  The accepted clock is the only clock an archive may
  index by.  A receipt chain is projected onto its accepted subsequence, the dropped attempts are
  counted, and a chain whose accepted clock is non-monotonic or gapped is REFUSED rather than
  silently reindexed -- silent reindexing turns a rollback into a physical trajectory.
* **A positional entity index that re-binds.**  Row *e* must mean the same filament, head, or clutch
  at every time index.  Per-step entity orderings are checked against one declared canonical order.
* **A compression approved on a global norm alone.**  A Frobenius error of 1e-4 can coexist with a
  rare local projection that is 95% wrong.  Every declared projection carries its own tolerance and
  any single failure vetoes the compressed form.  Fail closed.

Storage claims are measured, never estimated.  ``ArrayReduction.compression_ratio`` is a *parameter
count* ratio; it does not include NPY/NPZ framing, axis metadata, projection definitions, or the
receipt binding.  This module serializes both forms for real and reports a separate ``byte_ratio``,
and reports decode cost as a distribution over repeats rather than one timing.  The two numbers are
deliberately named differently because they routinely disagree in sign.

This module is CPU-only, imports no Warp, and grants no evidence authority: it is an unverified
observer over a runtime-supplied receipt claim.

Sanity Gate:
    * Dimensional: the time axis carries seconds taken from ``accepted_time_after_s``; every other
      axis carries dimensionless labels.  Projection units are declared per projection and never
      inferred.  Byte counts are bytes, parameter counts are float64 slots, decode cost is seconds --
      the three are never mixed into one "compression" number.
    * Boundary: an empty chain, a chain with zero accepted steps, a single-attempt chain, a
      single-entity axis and a two-axis (time x entity) archive are all exercised; the two-axis case
      routes to matrix SVD, not TT.
    * Conservation invariant: ``attempted = accepted + dropped`` holds by construction and is
      re-asserted in ``AcceptedTimeline``; accepted step indices are contiguous and accepted time is
      strictly increasing; the sum of serialized member sizes never exceeds the archive byte count.
    * Numerical: relative errors come from the scale-safe norms already used by ``reduction`` and
      ``tensor_train`` (an input at 1e308 must not turn overflow into a false pass).  Re-projection is
      required to reproduce a recorded scalar *bit-for-bit*, not to a tolerance.
    * Measurement protocol: bytes are read from a real serialized buffer (and cross-checked against a
      real file on disk); decode timings re-decode from those bytes on every repeat rather than
      reusing a live object, and are reported as min/median/max over the repeats.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import statistics
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from itertools import permutations
from pathlib import Path
from types import MappingProxyType
from typing import Any, Literal

import numpy as np

from aleph.virtual_cell.contracts import (
    ApproximationErrorReport,
    ArtifactEnvelope,
    EvidenceSource,
    RepresentationDescriptor,
    RepresentationKind,
    UncertaintySource,
)
from aleph.virtual_cell.reduction import (
    ArrayReduction,
    evaluate_projection_errors,
    select_array_representation,
)
from aleph.virtual_cell.sidecar import ArrayAxis, ArraySidecarDescriptor
from aleph.virtual_cell.step_receipt import StepReceipt, validate_step_receipt_chain
from aleph.virtual_cell.tensor_train import tt_svd

__all__ = [
    "AcceptedTimeline",
    "ArchiveByteMeasurement",
    "AxisOrderSearchReport",
    "AxisRole",
    "CompressionDecision",
    "DecodeCostMeasurement",
    "ProjectionRecord",
    "ScalarProjectionSpec",
    "ScalarReduction",
    "TelemetryArchive",
    "TelemetryExpansionVerdict",
    "TensorAxisSpec",
    "build_telemetry_archive",
    "compression_envelope",
    "gate_compression",
    "measure_archive_bytes",
    "measure_decode_cost",
    "project_accepted_timeline",
    "search_axis_order",
    "serialize_archive_bytes",
]

ARCHIVE_SCHEMA_ID = "ffn.virtual_cell.telemetry.archive/1"
METADATA_MEMBER = "archive_metadata_json_utf8"
ERROR_METRIC_ID = "scale-safe-relative-frobenius"

_EXPANSION_ROUTES: Mapping[str, str] = MappingProxyType(
    {
        "rank-failure": "ensemble / regime partition / latent coordinate",
        "projection-failure": "ensemble / regime partition / latent coordinate",
        "axis-order-coverage": "solver / preconditioner / adaptive design",
        "byte-overhead": "solver / preconditioner / adaptive design",
    }
)


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _positive_int(value: object, *, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{what} must be a positive integer")
    return int(value)


def _finite_float(value: object, *, what: str, minimum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{what} must be a finite real number")
    result = float(value)
    if minimum is not None and result < minimum:
        raise ValueError(f"{what} must be at least {minimum}")
    return result


def _canonical_json(payload: Any) -> bytes:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def _sha256_json(payload: Any) -> str:
    return hashlib.sha256(_canonical_json(payload)).hexdigest()


def _finite_tensor(values: Any) -> np.ndarray:
    array = np.array(np.asarray(values, dtype=np.float64), copy=True)
    if array.size == 0:
        raise ValueError("archive values must not be empty")
    if not np.all(np.isfinite(array)):
        raise ValueError("archive values must contain only finite numbers")
    return array


class AxisRole(StrEnum):
    """The MEANING of a tensor axis, declared rather than inferred from position."""

    ACCEPTED_TIME = "accepted-physical-time"
    ENTITY = "stable-entity-id"
    COMPONENT = "component"
    CONNECTOR = "connector"
    MODE = "mode"
    SPATIAL_COMPONENT = "spatial-component"


_LABELLED_ROLES = frozenset(AxisRole) - {AxisRole.ACCEPTED_TIME}


@dataclass(frozen=True, slots=True)
class TensorAxisSpec:
    """One archive axis, its role, and the labels or coordinates that make an index meaningful."""

    name: str
    role: AxisRole
    semantics: str
    labels: tuple[str, ...] | None = None
    coordinates: tuple[float, ...] | None = None
    unit: str = "1"

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _nonempty(self.name, what="axis name"))
        object.__setattr__(self, "semantics", _nonempty(self.semantics, what="axis semantics"))
        object.__setattr__(self, "unit", _nonempty(self.unit, what="axis unit"))
        if not isinstance(self.role, AxisRole):
            raise TypeError("axis role must be an AxisRole")
        if self.role is AxisRole.ACCEPTED_TIME:
            if self.labels is not None:
                raise ValueError("an accepted-time axis is coordinate-valued, not label-valued")
            if self.coordinates is None:
                raise ValueError("an accepted-time axis needs its accepted physical coordinates")
            coordinates = tuple(
                _finite_float(value, what="accepted time coordinate", minimum=0.0)
                for value in self.coordinates
            )
            if not coordinates:
                raise ValueError("an accepted-time axis needs at least one coordinate")
            if any(
                later <= earlier
                for earlier, later in zip(coordinates, coordinates[1:], strict=False)
            ):
                raise ValueError("accepted time coordinates must be strictly increasing")
            object.__setattr__(self, "coordinates", coordinates)
            return
        if self.coordinates is not None:
            raise ValueError(f"a {self.role.value} axis is label-valued, not coordinate-valued")
        if self.labels is None:
            raise ValueError(f"a {self.role.value} axis needs stable labels")
        labels = tuple(_nonempty(label, what=f"{self.role.value} label") for label in self.labels)
        if not labels:
            raise ValueError(f"a {self.role.value} axis needs at least one label")
        if len(set(labels)) != len(labels):
            raise ValueError(f"{self.role.value} labels must be unique")
        object.__setattr__(self, "labels", labels)

    @property
    def length(self) -> int:
        """Number of indices along this axis."""
        return len(self.coordinates) if self.coordinates is not None else len(self.labels or ())

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-able axis declaration."""
        return {
            "name": self.name,
            "role": self.role.value,
            "semantics": self.semantics,
            "unit": self.unit,
            "length": self.length,
            "labels": list(self.labels) if self.labels is not None else None,
            "coordinates": list(self.coordinates) if self.coordinates is not None else None,
        }

    @property
    def axis_sha256(self) -> str:
        """SHA-256 over the canonical axis declaration, including every label or coordinate."""
        return _sha256_json(self.to_dict())

    @property
    def coordinate_sha256(self) -> str:
        """SHA-256 over only the ordered labels or coordinates of this axis."""
        payload = (
            list(self.coordinates) if self.coordinates is not None else list(self.labels or ())
        )
        return _sha256_json(payload)

    def to_array_axis(self) -> ArrayAxis:
        """Project onto the sidecar axis contract shared with image and array artifacts."""
        return ArrayAxis(
            name=self.name,
            length=self.length,
            semantics=f"{self.role.value}: {self.semantics}",
            unit=self.unit,
            coordinates_sha256=self.coordinate_sha256,
        )


@dataclass(frozen=True, slots=True)
class AcceptedTimeline:
    """The accepted subsequence of a receipt chain, with the dropped attempts counted."""

    run_id: str
    decision_contract_sha256: str
    state_coverage_sha256: str
    accepted_time_s: tuple[float, ...]
    accepted_step_index: tuple[int, ...]
    accepted_attempt_index: tuple[int, ...]
    accepted_state_hashes: tuple[str, ...]
    attempted_step_count: int
    dropped_attempt_indices: tuple[int, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "run_id", _nonempty(self.run_id, what="run_id"))
        lengths = {
            len(self.accepted_time_s),
            len(self.accepted_step_index),
            len(self.accepted_attempt_index),
            len(self.accepted_state_hashes),
        }
        if len(lengths) != 1:
            raise ValueError("accepted timeline columns must have equal length")
        if not self.accepted_time_s:
            raise ValueError("an accepted timeline needs at least one accepted step")
        if any(
            later <= earlier
            for earlier, later in zip(self.accepted_time_s, self.accepted_time_s[1:], strict=False)
        ):
            raise ValueError("accepted physical time must be strictly increasing")
        if any(
            later != earlier + 1
            for earlier, later in zip(
                self.accepted_step_index, self.accepted_step_index[1:], strict=False
            )
        ):
            raise ValueError("accepted step indices must be contiguous")
        if self.attempted_step_count != self.accepted_step_count + len(
            self.dropped_attempt_indices
        ):
            raise ValueError("attempted steps must equal accepted steps plus dropped attempts")

    @property
    def accepted_step_count(self) -> int:
        """Number of time indices this timeline authorises."""
        return len(self.accepted_time_s)

    @property
    def dropped_attempt_count(self) -> int:
        """Number of attempted-but-rejected steps that were denied a time index."""
        return len(self.dropped_attempt_indices)

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-able timeline record."""
        return {
            "run_id": self.run_id,
            "decision_contract_sha256": self.decision_contract_sha256,
            "state_coverage_sha256": self.state_coverage_sha256,
            "accepted_time_s": list(self.accepted_time_s),
            "accepted_step_index": list(self.accepted_step_index),
            "accepted_attempt_index": list(self.accepted_attempt_index),
            "accepted_state_hashes": list(self.accepted_state_hashes),
            "attempted_step_count": self.attempted_step_count,
            "accepted_step_count": self.accepted_step_count,
            "dropped_attempt_indices": list(self.dropped_attempt_indices),
            "dropped_attempt_count": self.dropped_attempt_count,
        }


def project_accepted_timeline(receipts: Sequence[StepReceipt]) -> AcceptedTimeline:
    """Project an attempted receipt chain onto its accepted subsequence, or refuse.

    Args:
        receipts: One contiguous attempted-step receipt chain from a single run.

    Returns:
        The accepted timeline, recording how many attempts were dropped.

    Raises:
        TypeError: If the chain contains a non-``StepReceipt`` value.
        ValueError: If the attempted chain is not contiguous, if the accepted step clock or accepted
            physical time is non-monotonic or gapped, if a run-invariant field changes, or if no step
            was accepted.  The adapter never reindexes a broken clock into a plausible one.
    """
    chain = tuple(receipts)
    if not chain:
        raise ValueError("receipt chain must not be empty")
    if any(not isinstance(receipt, StepReceipt) for receipt in chain):
        raise TypeError("receipt chain must contain only StepReceipt values")
    first = chain[0]
    for previous, current in zip(chain, chain[1:], strict=False):
        if current.attempted_step_index != previous.attempted_step_index + 1:
            raise ValueError(
                "attempted-step index is not contiguous; the dropped-attempt count would be a guess"
            )
        if current.accepted_step_before < previous.accepted_step_after:
            raise ValueError("accepted-step clock is non-monotonic; refusing to reindex")
        if current.accepted_step_before > previous.accepted_step_after:
            raise ValueError("accepted-step clock has a gap; refusing to reindex")
        if current.accepted_time_before_s < previous.accepted_time_after_s:
            raise ValueError("accepted physical time is non-monotonic; refusing to reindex")
        if current.accepted_time_before_s > previous.accepted_time_after_s:
            raise ValueError("accepted physical time has a gap; refusing to reindex")
    validate_step_receipt_chain(chain)

    accepted = tuple(receipt for receipt in chain if receipt.accepted)
    if not accepted:
        raise ValueError("receipt chain accepted no step; an archive cannot index attempted state")
    return AcceptedTimeline(
        run_id=first.run_id,
        decision_contract_sha256=first.decision_contract_sha256,
        state_coverage_sha256=first.state_coverage_sha256,
        accepted_time_s=tuple(receipt.accepted_time_after_s for receipt in accepted),
        accepted_step_index=tuple(receipt.accepted_step_after for receipt in accepted),
        accepted_attempt_index=tuple(receipt.attempted_step_index for receipt in accepted),
        accepted_state_hashes=tuple(receipt.state_hash_after for receipt in accepted),
        attempted_step_count=len(chain),
        dropped_attempt_indices=tuple(
            receipt.attempted_step_index for receipt in chain if not receipt.accepted
        ),
    )


class ScalarReduction(StrEnum):
    """The whitelisted reductions a stored projection may apply after its index selection."""

    MEAN = "mean"
    SUM = "sum"
    MAX = "max"
    MIN = "min"
    ABS_MAX = "abs-max"
    STD = "std"
    VALUE = "value"


@dataclass(frozen=True, slots=True)
class ScalarProjectionSpec:
    """A named, versioned, content-hashed scalar projection stored as DATA, not as a function.

    A free function cannot be archived, diffed, or re-run against a future build, so a headline
    computed by one is unreproducible in principle.  The selection fixes indices on named axes and the
    reduction collapses whatever is left.  ``tolerance`` is part of the definition, and therefore of
    the hash: a projection is both a number and the gate it must pass, and quietly relaxing the gate
    must not be able to masquerade as the historical projection.
    """

    name: str
    version: str
    reduction: ScalarReduction
    tolerance: float
    selection: Mapping[str, int] = field(default_factory=dict)
    unit: str = "1"

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _nonempty(self.name, what="projection name"))
        object.__setattr__(self, "version", _nonempty(self.version, what="projection version"))
        object.__setattr__(self, "unit", _nonempty(self.unit, what="projection unit"))
        if not isinstance(self.reduction, ScalarReduction):
            raise TypeError("projection reduction must be a ScalarReduction")
        object.__setattr__(
            self,
            "tolerance",
            _finite_float(self.tolerance, what="projection tolerance", minimum=0.0),
        )
        if not isinstance(self.selection, Mapping):
            raise TypeError("projection selection must be a mapping of axis name to index")
        normalized: dict[str, int] = {}
        for raw_axis, raw_index in self.selection.items():
            axis = _nonempty(raw_axis, what="projection selection axis")
            if axis in normalized:
                raise ValueError("projection selection axes must not collide after normalization")
            if isinstance(raw_index, bool) or not isinstance(raw_index, int) or raw_index < 0:
                raise ValueError(
                    f"projection selection index for {axis!r} must be a nonnegative int"
                )
            normalized[axis] = int(raw_index)
        object.__setattr__(self, "selection", MappingProxyType(normalized))

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-able projection definition."""
        return {
            "name": self.name,
            "version": self.version,
            "reduction": self.reduction.value,
            "tolerance": self.tolerance,
            "selection": dict(sorted(self.selection.items())),
            "unit": self.unit,
        }

    @property
    def definition_sha256(self) -> str:
        """SHA-256 over the whole definition; editing any field mints a different projection."""
        return _sha256_json(self.to_dict())

    def evaluate(self, values: Any, axis_names: Sequence[str]) -> float:
        """Apply the stored selection and reduction to an archive-shaped array.

        Args:
            values: Array whose axes are named by ``axis_names`` in order.
            axis_names: Ordered axis names of ``values``.

        Returns:
            The projected scalar.
        """
        array = np.asarray(values, dtype=np.float64)
        names = tuple(_nonempty(name, what="axis name") for name in axis_names)
        if len(names) != array.ndim:
            raise ValueError("axis_names must name every axis of the array")
        if len(set(names)) != len(names):
            raise ValueError("axis names must be unique")
        unknown = set(self.selection) - set(names)
        if unknown:
            raise ValueError("projection selects unknown axes: " + ", ".join(sorted(unknown)))
        index: list[int | slice] = []
        for axis, size in zip(names, array.shape, strict=True):
            if axis in self.selection:
                position = self.selection[axis]
                if position >= size:
                    raise ValueError(
                        f"projection index {position} is out of range on axis {axis!r} of size {size}"
                    )
                index.append(position)
            else:
                index.append(slice(None))
        view = array[tuple(index)]
        if self.reduction is ScalarReduction.VALUE:
            if view.size != 1:
                raise ValueError("a 'value' projection must select a single element on every axis")
            return float(np.reshape(view, -1)[0])
        if self.reduction is ScalarReduction.MEAN:
            return float(np.mean(view))
        if self.reduction is ScalarReduction.SUM:
            return float(np.sum(view))
        if self.reduction is ScalarReduction.MAX:
            return float(np.max(view))
        if self.reduction is ScalarReduction.MIN:
            return float(np.min(view))
        if self.reduction is ScalarReduction.ABS_MAX:
            return float(np.max(np.abs(view)))
        return float(np.std(view))

    def as_callable(self, axis_names: Sequence[str]) -> Callable[[np.ndarray], float]:
        """Bind this stored projection to an axis naming for ``evaluate_projection_errors``."""
        names = tuple(axis_names)

        def _projection(array: np.ndarray) -> float:
            return self.evaluate(array, names)

        return _projection


@dataclass(frozen=True, slots=True)
class ProjectionRecord:
    """One recorded scalar bound to the exact projection definition that produced it."""

    name: str
    version: str
    definition_sha256: str
    value: float
    unit: str

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-able scalar record."""
        return {
            "name": self.name,
            "version": self.version,
            "definition_sha256": self.definition_sha256,
            "value": self.value,
            "unit": self.unit,
        }


@dataclass(frozen=True, slots=True)
class TelemetryExpansionVerdict:
    """A typed wall report naming which programme expansion the failure opens."""

    wall_id: str
    observation: str
    expansion: str
    measured: Mapping[str, float]

    def __post_init__(self) -> None:
        object.__setattr__(self, "wall_id", _nonempty(self.wall_id, what="wall_id"))
        if self.wall_id not in _EXPANSION_ROUTES:
            raise ValueError(f"wall_id must be one of {sorted(_EXPANSION_ROUTES)}")
        object.__setattr__(self, "observation", _nonempty(self.observation, what="observation"))
        if self.expansion != _EXPANSION_ROUTES[self.wall_id]:
            raise ValueError("expansion must match the declared routing table")
        measured = {
            _nonempty(key, what="measured key"): _finite_float(value, what=f"measured[{key!r}]")
            for key, value in dict(self.measured).items()
        }
        object.__setattr__(self, "measured", MappingProxyType(measured))

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-able verdict."""
        return {
            "wall_id": self.wall_id,
            "observation": self.observation,
            "expansion": self.expansion,
            "measured": dict(sorted(self.measured.items())),
        }


@dataclass(frozen=True, slots=True)
class CompressionDecision:
    """A fail-closed compression verdict gated on the global norm AND every declared projection."""

    global_tolerance: float
    measured_global_relative_error: float
    candidate_kind: Literal["dense", "svd", "tensor-train"]
    used_kind: Literal["dense", "svd", "tensor-train"]
    accepted: bool
    rejection_reason: str | None
    projection_errors: Mapping[str, float]
    projection_tolerances: Mapping[str, float]
    failed_projections: tuple[str, ...]
    stored_parameter_count: int
    dense_parameter_count: int
    reduction: ArrayReduction
    axis_order: tuple[int, ...] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "global_tolerance",
            _finite_float(self.global_tolerance, what="global_tolerance", minimum=0.0),
        )
        object.__setattr__(
            self,
            "measured_global_relative_error",
            _finite_float(
                self.measured_global_relative_error,
                what="measured_global_relative_error",
                minimum=0.0,
            ),
        )
        if not isinstance(self.accepted, bool):
            raise TypeError("accepted must be bool")
        if self.accepted and self.rejection_reason is not None:
            raise ValueError("an accepted compression must not carry a rejection reason")
        if not self.accepted and not self.rejection_reason:
            raise ValueError("a rejected compression must record why")
        if not self.accepted and self.used_kind != "dense":
            raise ValueError("a rejected compression must fall back to dense")
        if self.accepted and self.failed_projections:
            raise ValueError("an accepted compression cannot have failed projections")
        object.__setattr__(
            self, "projection_errors", MappingProxyType(dict(self.projection_errors))
        )
        object.__setattr__(
            self, "projection_tolerances", MappingProxyType(dict(self.projection_tolerances))
        )
        if set(self.projection_errors) != set(self.projection_tolerances):
            raise ValueError("projection errors and tolerances must cover the same names")
        _positive_int(self.stored_parameter_count, what="stored_parameter_count")
        _positive_int(self.dense_parameter_count, what="dense_parameter_count")

    @property
    def parameter_count_ratio(self) -> float:
        """Dense float64 slots per stored float64 slot.

        This is NOT a storage claim: it excludes NPY/NPZ framing, axis metadata, projection
        definitions and the receipt binding.  Use :class:`ArchiveByteMeasurement` for bytes.
        """
        return self.dense_parameter_count / self.stored_parameter_count

    @property
    def expansion_verdict(self) -> TelemetryExpansionVerdict | None:
        """Route a rejection to the expansion it opens, or ``None`` when the reduction was used."""
        if self.accepted:
            return None
        measured = {
            "global_relative_error": self.measured_global_relative_error,
            "parameter_count_ratio": self.parameter_count_ratio,
        }
        if self.failed_projections:
            worst = max(self.failed_projections, key=lambda name: self.projection_errors[name])
            measured["worst_projection_error"] = self.projection_errors[worst]
            return TelemetryExpansionVerdict(
                wall_id="projection-failure",
                observation=(
                    f"projection {worst!r} error {self.projection_errors[worst]:.6g} exceeds its "
                    f"tolerance {self.projection_tolerances[worst]:.6g} while the global norm passed"
                ),
                expansion=_EXPANSION_ROUTES["projection-failure"],
                measured=measured,
            )
        return TelemetryExpansionVerdict(
            wall_id="rank-failure",
            observation=f"no admissible reduced representation: {self.rejection_reason}",
            expansion=_EXPANSION_ROUTES["rank-failure"],
            measured=measured,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-able compression verdict."""
        verdict = self.expansion_verdict
        return {
            "global_tolerance": self.global_tolerance,
            "measured_global_relative_error": self.measured_global_relative_error,
            "candidate_kind": self.candidate_kind,
            "used_kind": self.used_kind,
            "accepted": self.accepted,
            "rejection_reason": self.rejection_reason,
            "projection_errors": dict(sorted(self.projection_errors.items())),
            "projection_tolerances": dict(sorted(self.projection_tolerances.items())),
            "failed_projections": list(self.failed_projections),
            "stored_parameter_count": self.stored_parameter_count,
            "dense_parameter_count": self.dense_parameter_count,
            "parameter_count_ratio": self.parameter_count_ratio,
            "axis_order": list(self.axis_order) if self.axis_order is not None else None,
            "expansion_verdict": verdict.to_dict() if verdict is not None else None,
        }


def gate_compression(
    values: Any,
    *,
    axis_names: Sequence[str],
    projections: Sequence[ScalarProjectionSpec] = (),
    global_tolerance: float,
    axis_order: tuple[int, ...] | None = None,
    max_rank: int | None = None,
) -> CompressionDecision:
    """Compress only if the global norm AND every declared projection tolerance hold.

    Args:
        values: Dense finite archive array.
        axis_names: Ordered axis names, used to bind stored projections to axes.
        projections: Declared projections, each with its own tolerance.
        global_tolerance: Requested relative Frobenius tolerance in ``[0, 1)``.
        axis_order: Optional TT axis permutation, normally from :func:`search_axis_order`.
        max_rank: Optional cap on every internal rank.

    Returns:
        The verdict, including per-projection measured errors whether or not it was accepted.
    """
    array = _finite_tensor(values)
    names = tuple(_nonempty(name, what="axis name") for name in axis_names)
    if len(names) != array.ndim:
        raise ValueError("axis_names must name every axis of the array")
    if len(set(names)) != len(names):
        raise ValueError("axis names must be unique")
    specs = tuple(projections)
    if any(not isinstance(spec, ScalarProjectionSpec) for spec in specs):
        raise TypeError("projections must be ScalarProjectionSpec values")
    tolerances = {spec.name: spec.tolerance for spec in specs}
    if len(tolerances) != len(specs):
        raise ValueError("projection names must be unique")

    candidate = select_array_representation(
        array,
        relative_tolerance=global_tolerance,
        max_rank=max_rank,
        axis_order=axis_order,
    )
    errors: Mapping[str, float] = {}
    if specs:
        errors = evaluate_projection_errors(
            array,
            candidate.reconstruct(),
            {spec.name: spec.as_callable(names) for spec in specs},
        )
    failed = tuple(name for name in sorted(errors) if errors[name] > tolerances[name])

    if candidate.kind == "dense":
        return CompressionDecision(
            global_tolerance=float(global_tolerance),
            measured_global_relative_error=candidate.measured_relative_error,
            candidate_kind="dense",
            used_kind="dense",
            accepted=False,
            rejection_reason=candidate.fallback_reason,
            projection_errors=errors,
            projection_tolerances=tolerances,
            failed_projections=failed,
            stored_parameter_count=candidate.stored_parameter_count,
            dense_parameter_count=candidate.dense_parameter_count,
            reduction=candidate,
            axis_order=axis_order,
        )
    if failed:
        return CompressionDecision(
            global_tolerance=float(global_tolerance),
            measured_global_relative_error=candidate.measured_relative_error,
            candidate_kind=candidate.kind,
            used_kind="dense",
            accepted=False,
            rejection_reason="declared-projection-error-exceeds-tolerance",
            projection_errors=errors,
            projection_tolerances=tolerances,
            failed_projections=failed,
            stored_parameter_count=candidate.stored_parameter_count,
            dense_parameter_count=candidate.dense_parameter_count,
            reduction=candidate,
            axis_order=axis_order,
        )
    return CompressionDecision(
        global_tolerance=float(global_tolerance),
        measured_global_relative_error=candidate.measured_relative_error,
        candidate_kind=candidate.kind,
        used_kind=candidate.kind,
        accepted=True,
        rejection_reason=None,
        projection_errors=errors,
        projection_tolerances=tolerances,
        failed_projections=(),
        stored_parameter_count=candidate.stored_parameter_count,
        dense_parameter_count=candidate.dense_parameter_count,
        reduction=candidate,
        axis_order=axis_order,
    )


@dataclass(frozen=True, slots=True)
class AxisOrderSearchReport:
    """Best TT axis ordering found, and explicitly which region of the space was NOT searched."""

    ndim: int
    total_orderings: int
    evaluated_orderings: tuple[tuple[int, ...], ...]
    strategy: Literal["exhaustive", "heuristic-subset"]
    strategy_description: str
    requested_tolerance: float
    best_order: tuple[int, ...]
    best_storage_size: int
    best_relative_error: float
    dense_parameter_count: int
    any_order_met_tolerance: bool
    storage_by_order: Mapping[tuple[int, ...], int]
    error_by_order: Mapping[tuple[int, ...], float]

    def __post_init__(self) -> None:
        _positive_int(self.ndim, what="ndim")
        _positive_int(self.total_orderings, what="total_orderings")
        if not self.evaluated_orderings:
            raise ValueError("an axis-order search must evaluate at least one ordering")
        if len(set(self.evaluated_orderings)) != len(self.evaluated_orderings):
            raise ValueError("evaluated orderings must be unique")
        if self.evaluated_count > self.total_orderings:
            raise ValueError("cannot evaluate more orderings than exist")
        object.__setattr__(self, "storage_by_order", MappingProxyType(dict(self.storage_by_order)))
        object.__setattr__(self, "error_by_order", MappingProxyType(dict(self.error_by_order)))
        if set(self.storage_by_order) != set(self.evaluated_orderings) or set(
            self.error_by_order
        ) != set(self.evaluated_orderings):
            raise ValueError("per-ordering measurements must cover exactly the evaluated orderings")
        if self.strategy == "exhaustive" and self.evaluated_count != self.total_orderings:
            raise ValueError("an exhaustive search must evaluate every ordering")

    @property
    def evaluated_count(self) -> int:
        """Number of axis orderings actually measured."""
        return len(self.evaluated_orderings)

    @property
    def unsearched_count(self) -> int:
        """Number of axis orderings never measured; a silent cap is exactly this number hidden."""
        return self.total_orderings - self.evaluated_count

    @property
    def unsearched_fraction(self) -> float:
        """Fraction of the ordering space left unsearched, reported with every best ordering."""
        return self.unsearched_count / self.total_orderings

    @property
    def best_parameter_count_ratio(self) -> float:
        """Dense float64 slots per stored slot at the best ordering found; not a byte ratio."""
        return self.dense_parameter_count / self.best_storage_size

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-able search report, coverage included."""
        return {
            "ndim": self.ndim,
            "total_orderings": self.total_orderings,
            "evaluated_count": self.evaluated_count,
            "unsearched_count": self.unsearched_count,
            "unsearched_fraction": self.unsearched_fraction,
            "strategy": self.strategy,
            "strategy_description": self.strategy_description,
            "requested_tolerance": self.requested_tolerance,
            "best_order": list(self.best_order),
            "best_storage_size": self.best_storage_size,
            "best_relative_error": self.best_relative_error,
            "best_parameter_count_ratio": self.best_parameter_count_ratio,
            "any_order_met_tolerance": self.any_order_met_tolerance,
            "evaluated_orderings": [list(order) for order in self.evaluated_orderings],
        }


_HEURISTIC_DESCRIPTION = (
    "bounded deterministic subset: identity, reverse, length-ascending, length-descending, every "
    "cyclic rotation of identity, then lexicographic permutations until max_orders is reached"
)


def _candidate_orders(shape: tuple[int, ...], *, max_orders: int) -> list[tuple[int, ...]]:
    ndim = len(shape)
    identity = tuple(range(ndim))
    seeds: list[tuple[int, ...]] = [
        identity,
        tuple(reversed(identity)),
        tuple(sorted(identity, key=lambda axis: (shape[axis], axis))),
        tuple(sorted(identity, key=lambda axis: (-shape[axis], axis))),
    ]
    seeds.extend(identity[shift:] + identity[:shift] for shift in range(1, ndim))
    ordered: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    for order in seeds:
        if order not in seen:
            seen.add(order)
            ordered.append(order)
        if len(ordered) >= max_orders:
            return ordered
    for order in permutations(identity):
        if order not in seen:
            seen.add(order)
            ordered.append(order)
        if len(ordered) >= max_orders:
            break
    return ordered


def search_axis_order(
    values: Any,
    *,
    relative_tolerance: float,
    max_orders: int = 24,
    max_rank: int | None = None,
) -> AxisOrderSearchReport:
    """Search axis orderings for the cheapest TT representation and report the unsearched region.

    TT rank is a property of the array AND the axis order, so an ordering search is part of the
    measurement rather than a tuning step.  The search is bounded: exhaustive when the whole space
    fits inside ``max_orders``, otherwise a documented deterministic subset.  Either way the report
    carries ``unsearched_fraction``, because a best ordering quoted without its coverage is a silent
    cap on the claim.

    Args:
        values: Dense finite tensor with at least two axes.
        relative_tolerance: Requested relative Frobenius tolerance in ``[0, 1)``.
        max_orders: Upper bound on evaluated orderings.
        max_rank: Optional cap on every internal TT rank.

    Returns:
        The search report, including the per-ordering storage and measured error.
    """
    tensor = _finite_tensor(values)
    if tensor.ndim < 2:
        raise ValueError("axis-order search needs at least two axes")
    _positive_int(max_orders, what="max_orders")
    total = math.factorial(tensor.ndim)
    exhaustive = total <= max_orders
    orders = (
        [tuple(order) for order in permutations(range(tensor.ndim))]
        if exhaustive
        else _candidate_orders(tuple(tensor.shape), max_orders=max_orders)
    )

    storage: dict[tuple[int, ...], int] = {}
    errors: dict[tuple[int, ...], float] = {}
    for order in orders:
        approximation = tt_svd(
            tensor,
            relative_tolerance=relative_tolerance,
            max_rank=max_rank,
            axis_order=order,
        )
        storage[order] = approximation.storage_size
        errors[order] = approximation.relative_error

    slack = max(1e-14, float(relative_tolerance) * 1e-10)
    admissible = [order for order in orders if errors[order] <= float(relative_tolerance) + slack]
    pool = admissible if admissible else list(orders)
    best = min(pool, key=lambda order: (storage[order], errors[order], order))
    return AxisOrderSearchReport(
        ndim=int(tensor.ndim),
        total_orderings=int(total),
        evaluated_orderings=tuple(orders),
        strategy="exhaustive" if exhaustive else "heuristic-subset",
        strategy_description=("full permutation space" if exhaustive else _HEURISTIC_DESCRIPTION),
        requested_tolerance=float(relative_tolerance),
        best_order=best,
        best_storage_size=storage[best],
        best_relative_error=errors[best],
        dense_parameter_count=int(tensor.size),
        any_order_met_tolerance=bool(admissible),
        storage_by_order=storage,
        error_by_order=errors,
    )


def _representation_arrays(reduction: ArrayReduction) -> dict[str, np.ndarray]:
    if reduction.kind == "dense":
        assert reduction.dense is not None
        return {"values": np.asarray(reduction.dense)}
    if reduction.kind == "svd":
        assert reduction.svd is not None
        return {
            "u": np.asarray(reduction.svd.u),
            "singular_values": np.asarray(reduction.svd.singular_values),
            "vh": np.asarray(reduction.svd.vh),
        }
    assert reduction.tensor_train is not None
    return {
        f"core_{index:03d}": np.asarray(core)
        for index, core in enumerate(reduction.tensor_train.cores)
    }


def _representation_metadata(reduction: ArrayReduction) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "kind": reduction.kind,
        "original_shape": list(reduction.original_shape),
        "requested_tolerance": reduction.requested_tolerance,
        "measured_relative_error": reduction.measured_relative_error,
        "stored_parameter_count": reduction.stored_parameter_count,
        "dense_parameter_count": reduction.dense_parameter_count,
        "fallback_reason": reduction.fallback_reason,
    }
    if reduction.tensor_train is not None:
        payload["tt_ranks"] = list(reduction.tensor_train.ranks)
        payload["tt_axis_order"] = list(reduction.tensor_train.axis_order)
    if reduction.svd is not None:
        payload["svd_rank"] = reduction.svd.rank
    return payload


def _measurement_metadata_pair(
    array: np.ndarray,
    reduction: ArrayReduction,
    metadata: Mapping[str, Any] | None,
    *,
    purpose: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build the dense-reference and reduced metadata blocks for a like-for-like measurement.

    When the reduction already IS dense the two forms are the same array, so they must also carry
    byte-identical metadata; otherwise a differing ``fallback_reason`` string alone would show up as
    a spurious few-byte "compression".
    """
    base = dict(metadata or {})
    reduced_metadata = {**base, "representation": _representation_metadata(reduction)}
    if reduction.kind == "dense":
        return dict(reduced_metadata), reduced_metadata
    dense_metadata = {
        **base,
        "representation": {
            "kind": "dense",
            "original_shape": list(reduction.original_shape),
            "requested_tolerance": reduction.requested_tolerance,
            "measured_relative_error": 0.0,
            "stored_parameter_count": int(array.size),
            "dense_parameter_count": int(array.size),
            "fallback_reason": f"dense-reference-form-for-{purpose}-measurement",
        },
    }
    return dense_metadata, reduced_metadata


def serialize_archive_bytes(arrays: Mapping[str, np.ndarray], metadata: Mapping[str, Any]) -> bytes:
    """Serialize payload arrays and metadata into one uncompressed NPZ buffer.

    Metadata travels inside the same container as a UTF-8 byte member, so the byte count includes it.
    A storage claim that excludes its own metadata is not a storage claim.

    Args:
        arrays: Payload arrays keyed by member name.
        metadata: JSON-able metadata; canonicalised before hashing or measuring.

    Returns:
        The serialized archive bytes.
    """
    if not arrays:
        raise ValueError("an archive needs at least one array member")
    if METADATA_MEMBER in arrays:
        raise ValueError(f"{METADATA_MEMBER!r} is reserved for the metadata member")
    payload = dict(arrays)
    payload[METADATA_MEMBER] = np.frombuffer(_canonical_json(metadata), dtype=np.uint8)
    buffer = io.BytesIO()
    np.savez(buffer, **payload)
    return buffer.getvalue()


def _decode_archive(payload: bytes) -> np.ndarray:
    """Reconstruct a dense array from serialized bytes alone, with no live object in scope."""
    with np.load(io.BytesIO(payload), allow_pickle=False) as handle:
        metadata = json.loads(bytes(handle[METADATA_MEMBER]).decode())
        representation = metadata["representation"]
        kind = representation["kind"]
        if kind == "dense":
            return np.asarray(handle["values"], dtype=np.float64)
        if kind == "svd":
            u = np.asarray(handle["u"], dtype=np.float64)
            singular_values = np.asarray(handle["singular_values"], dtype=np.float64)
            vh = np.asarray(handle["vh"], dtype=np.float64)
            return (u * singular_values[None, :]) @ vh
        cores = [
            np.asarray(handle[name], dtype=np.float64)
            for name in sorted(name for name in handle.files if name.startswith("core_"))
        ]
    result = cores[0][0, :, :]
    for core in cores[1:]:
        result = np.tensordot(result, core, axes=([-1], [0]))
    result = np.squeeze(result, axis=-1)
    inverse = tuple(int(axis) for axis in np.argsort(representation["tt_axis_order"]))
    return np.transpose(result, axes=inverse)


@dataclass(frozen=True, slots=True)
class ArchiveByteMeasurement:
    """Actual serialized byte counts for the dense and reduced forms, metadata included."""

    dense_total_bytes: int
    reduced_total_bytes: int
    metadata_bytes: int
    dense_parameter_count: int
    stored_parameter_count: int
    serializer: str = "numpy-npz-uncompressed"

    def __post_init__(self) -> None:
        for field_name in ("dense_total_bytes", "reduced_total_bytes", "metadata_bytes"):
            _positive_int(getattr(self, field_name), what=field_name)
        _positive_int(self.dense_parameter_count, what="dense_parameter_count")
        _positive_int(self.stored_parameter_count, what="stored_parameter_count")

    @property
    def byte_ratio(self) -> float:
        """Dense bytes per reduced byte.  Below one means the reduced archive is LARGER."""
        return self.dense_total_bytes / self.reduced_total_bytes

    @property
    def parameter_count_ratio(self) -> float:
        """Dense float64 slots per stored slot; a different number from :attr:`byte_ratio`."""
        return self.dense_parameter_count / self.stored_parameter_count

    @property
    def reduced_form_is_larger(self) -> bool:
        """True when compressing this archive costs more bytes than it saves."""
        return self.reduced_total_bytes >= self.dense_total_bytes

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-able byte measurement."""
        return {
            "serializer": self.serializer,
            "dense_total_bytes": self.dense_total_bytes,
            "reduced_total_bytes": self.reduced_total_bytes,
            "metadata_bytes": self.metadata_bytes,
            "byte_ratio": self.byte_ratio,
            "dense_parameter_count": self.dense_parameter_count,
            "stored_parameter_count": self.stored_parameter_count,
            "parameter_count_ratio": self.parameter_count_ratio,
            "reduced_form_is_larger": self.reduced_form_is_larger,
        }


def measure_archive_bytes(
    values: Any,
    reduction: ArrayReduction,
    *,
    metadata: Mapping[str, Any] | None = None,
    probe_path: str | Path | None = None,
) -> ArchiveByteMeasurement:
    """Measure real serialized bytes for both forms, with identical metadata in each.

    Args:
        values: The dense archive array the reduction was fitted to.
        reduction: The fitted representation.
        metadata: Archive metadata included verbatim in both serializations.
        probe_path: Optional path; both buffers are written there in turn and the on-disk sizes are
            required to match the in-memory counts, so the measurement is not an in-memory artefact.

    Returns:
        The byte measurement.
    """
    array = _finite_tensor(values)
    if not isinstance(reduction, ArrayReduction):
        raise TypeError("reduction must be an ArrayReduction")
    if tuple(array.shape) != reduction.original_shape:
        raise ValueError("reduction was fitted to a different shape than the supplied array")
    dense_metadata, reduced_metadata = _measurement_metadata_pair(
        array, reduction, metadata, purpose="byte"
    )
    dense_bytes = serialize_archive_bytes({"values": array}, dense_metadata)
    reduced_bytes = serialize_archive_bytes(_representation_arrays(reduction), reduced_metadata)

    if probe_path is not None:
        probe = Path(probe_path)
        for payload in (dense_bytes, reduced_bytes):
            probe.write_bytes(payload)
            if probe.stat().st_size != len(payload):
                raise ValueError("serialized byte count does not match the file written to disk")

    return ArchiveByteMeasurement(
        dense_total_bytes=len(dense_bytes),
        reduced_total_bytes=len(reduced_bytes),
        metadata_bytes=len(_canonical_json(reduced_metadata)),
        dense_parameter_count=int(array.size),
        stored_parameter_count=reduction.stored_parameter_count,
    )


@dataclass(frozen=True, slots=True)
class DecodeCostMeasurement:
    """Decode wall-clock cost as a distribution over repeats, reported apart from any byte ratio."""

    repeats: int
    dense_samples_s: tuple[float, ...]
    reduced_samples_s: tuple[float, ...]
    reconstruction_max_abs_difference: float

    def __post_init__(self) -> None:
        _positive_int(self.repeats, what="repeats")
        for field_name in ("dense_samples_s", "reduced_samples_s"):
            samples = tuple(
                _finite_float(value, what=f"{field_name} entry", minimum=0.0)
                for value in getattr(self, field_name)
            )
            if len(samples) != self.repeats:
                raise ValueError(f"{field_name} must hold one sample per repeat")
            object.__setattr__(self, field_name, samples)
        _finite_float(
            self.reconstruction_max_abs_difference,
            what="reconstruction_max_abs_difference",
            minimum=0.0,
        )

    @staticmethod
    def _summary(samples: tuple[float, ...]) -> dict[str, float]:
        return {
            "min_s": min(samples),
            "median_s": statistics.median(samples),
            "max_s": max(samples),
        }

    @property
    def dense_summary_s(self) -> Mapping[str, float]:
        """Min/median/max seconds to decode the dense form from bytes."""
        return MappingProxyType(self._summary(self.dense_samples_s))

    @property
    def reduced_summary_s(self) -> Mapping[str, float]:
        """Min/median/max seconds to decode and reconstruct the reduced form from bytes."""
        return MappingProxyType(self._summary(self.reduced_samples_s))

    @property
    def median_decode_slowdown(self) -> float:
        """Reduced median decode seconds per dense median decode second."""
        dense_median = statistics.median(self.dense_samples_s)
        if dense_median == 0.0:
            raise ValueError("dense decode median is zero; the timer resolution is insufficient")
        return statistics.median(self.reduced_samples_s) / dense_median

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-able decode-cost distribution."""
        return {
            "repeats": self.repeats,
            "dense_samples_s": list(self.dense_samples_s),
            "reduced_samples_s": list(self.reduced_samples_s),
            "dense_summary_s": dict(self.dense_summary_s),
            "reduced_summary_s": dict(self.reduced_summary_s),
            "reconstruction_max_abs_difference": self.reconstruction_max_abs_difference,
        }


def measure_decode_cost(
    values: Any,
    reduction: ArrayReduction,
    *,
    metadata: Mapping[str, Any] | None = None,
    repeats: int = 7,
) -> DecodeCostMeasurement:
    """Time decoding both forms from serialized bytes, once per repeat, and report the spread.

    Every repeat re-parses the buffer, so a live in-memory object cannot make the reduced form look
    free.  The result is a distribution: a single timing on a shared laptop is not a measurement.

    Args:
        values: The dense archive array.
        reduction: The fitted representation.
        metadata: Archive metadata serialized alongside both forms.
        repeats: Number of independent decodes per form.

    Returns:
        The decode-cost measurement, including the max absolute reconstruction difference so that a
        fast decode of the wrong array cannot pass as a cheap one.
    """
    array = _finite_tensor(values)
    _positive_int(repeats, what="repeats")
    dense_metadata, reduced_metadata = _measurement_metadata_pair(
        array, reduction, metadata, purpose="decode"
    )
    dense_bytes = serialize_archive_bytes({"values": array}, dense_metadata)
    reduced_bytes = serialize_archive_bytes(_representation_arrays(reduction), reduced_metadata)

    dense_samples: list[float] = []
    reduced_samples: list[float] = []
    difference = 0.0
    for _ in range(repeats):
        start = time.perf_counter()
        decoded_dense = _decode_archive(dense_bytes)
        dense_samples.append(time.perf_counter() - start)
        start = time.perf_counter()
        decoded_reduced = _decode_archive(reduced_bytes)
        reduced_samples.append(time.perf_counter() - start)
        difference = max(
            difference,
            float(np.max(np.abs(decoded_dense - array))),
            float(np.max(np.abs(decoded_reduced - array))),
        )
    return DecodeCostMeasurement(
        repeats=int(repeats),
        dense_samples_s=tuple(dense_samples),
        reduced_samples_s=tuple(reduced_samples),
        reconstruction_max_abs_difference=difference,
    )


@dataclass(frozen=True, slots=True)
class TelemetryArchive:
    """A lossless multi-index archive whose axes carry meaning and whose scalars are projections."""

    archive_id: str
    timeline: AcceptedTimeline
    axes: tuple[TensorAxisSpec, ...]
    values: np.ndarray
    projections: tuple[ScalarProjectionSpec, ...]
    recorded_scalars: tuple[ProjectionRecord, ...]
    value_semantics: str
    value_unit: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "archive_id", _nonempty(self.archive_id, what="archive_id"))
        object.__setattr__(
            self, "value_semantics", _nonempty(self.value_semantics, what="value_semantics")
        )
        object.__setattr__(self, "value_unit", _nonempty(self.value_unit, what="value_unit"))
        if not isinstance(self.timeline, AcceptedTimeline):
            raise TypeError("timeline must be an AcceptedTimeline")
        axes = tuple(self.axes)
        if len(axes) < 2 or any(not isinstance(axis, TensorAxisSpec) for axis in axes):
            raise ValueError("an archive needs at least an accepted-time axis and one other axis")
        names = tuple(axis.name for axis in axes)
        if len(set(names)) != len(names):
            raise ValueError("axis names must be unique")
        if axes[0].role is not AxisRole.ACCEPTED_TIME:
            raise ValueError("axis 0 must be the accepted-physical-time axis")
        if sum(axis.role is AxisRole.ACCEPTED_TIME for axis in axes) != 1:
            raise ValueError("exactly one accepted-time axis is allowed")
        if sum(axis.role is AxisRole.ENTITY for axis in axes) != 1:
            raise ValueError("exactly one stable-entity-id axis is required")
        object.__setattr__(self, "axes", axes)
        if axes[0].coordinates != self.timeline.accepted_time_s:
            raise ValueError("the time axis must carry the accepted timeline coordinates")

        array = _finite_tensor(self.values)
        if tuple(array.shape) != tuple(axis.length for axis in axes):
            raise ValueError("values shape does not match the declared axes")
        array.setflags(write=False)
        object.__setattr__(self, "values", array)

        specs = tuple(self.projections)
        if not specs or any(not isinstance(spec, ScalarProjectionSpec) for spec in specs):
            raise ValueError("an archive must declare at least one scalar projection")
        if len({spec.name for spec in specs}) != len(specs):
            raise ValueError("projection names must be unique")
        object.__setattr__(self, "projections", specs)
        records = tuple(self.recorded_scalars)
        if len(records) != len(specs) or any(
            not isinstance(record, ProjectionRecord) for record in records
        ):
            raise ValueError("one recorded scalar is required per declared projection")
        object.__setattr__(self, "recorded_scalars", records)
        self.verify_recorded_scalars()

    @property
    def axis_names(self) -> tuple[str, ...]:
        """Ordered axis names."""
        return tuple(axis.name for axis in self.axes)

    @property
    def entity_axis(self) -> TensorAxisSpec:
        """The single stable-entity-id axis."""
        return next(axis for axis in self.axes if axis.role is AxisRole.ENTITY)

    @property
    def shape(self) -> tuple[int, ...]:
        """Dense archive shape."""
        return tuple(int(size) for size in self.values.shape)

    @property
    def dense_parameter_count(self) -> int:
        """Number of float64 slots in the dense archive."""
        return int(self.values.size)

    def reproject(self, name: str) -> float:
        """Recompute one declared projection from the raw archive.

        Args:
            name: Declared projection name.

        Returns:
            The recomputed scalar, which must equal the recorded one bit-for-bit.
        """
        spec = next((candidate for candidate in self.projections if candidate.name == name), None)
        if spec is None:
            raise KeyError(f"no declared projection named {name!r}")
        return spec.evaluate(self.values, self.axis_names)

    def verify_recorded_scalars(self) -> None:
        """Raise unless every recorded scalar re-projects bit-for-bit from the raw archive."""
        by_name = {record.name: record for record in self.recorded_scalars}
        for spec in self.projections:
            record = by_name.get(spec.name)
            if record is None:
                raise ValueError(f"projection {spec.name!r} has no recorded scalar")
            if record.definition_sha256 != spec.definition_sha256:
                raise ValueError(
                    f"recorded scalar {spec.name!r} was produced by a different projection definition"
                )
            recomputed = spec.evaluate(self.values, self.axis_names)
            if recomputed != record.value:
                raise ValueError(
                    f"projection {spec.name!r} no longer reproduces its recorded scalar"
                )

    def dense_payload(self) -> bytes:
        """Serialize the dense archive as a pickle-free NPY payload."""
        buffer = io.BytesIO()
        np.save(buffer, np.ascontiguousarray(self.values), allow_pickle=False)
        return buffer.getvalue()

    def dense_sidecar(self, *, uri: str) -> ArraySidecarDescriptor:
        """Describe the dense payload with the shared content-addressed sidecar contract."""
        payload = self.dense_payload()
        return ArraySidecarDescriptor(
            uri=uri,
            sha256=hashlib.sha256(payload).hexdigest(),
            media_type="application/x-npy",
            schema_id=ARCHIVE_SCHEMA_ID,
            dtype=np.dtype(np.float64).str,
            value_semantics=self.value_semantics,
            value_unit=self.value_unit,
            axes=tuple(axis.to_array_axis() for axis in self.axes),
            byte_size=len(payload),
            array_order="C",
        )

    def to_metadata(self) -> dict[str, Any]:
        """Return the JSON-able archive metadata that travels with every serialization."""
        return {
            "schema_id": ARCHIVE_SCHEMA_ID,
            "archive_id": self.archive_id,
            "value_semantics": self.value_semantics,
            "value_unit": self.value_unit,
            "timeline": self.timeline.to_dict(),
            "axes": [axis.to_dict() for axis in self.axes],
            "projections": [spec.to_dict() for spec in self.projections],
            "recorded_scalars": [record.to_dict() for record in self.recorded_scalars],
        }

    @property
    def archive_sha256(self) -> str:
        """SHA-256 binding the dense payload to the canonical archive metadata."""
        return hashlib.sha256(
            hashlib.sha256(self.dense_payload()).hexdigest().encode()
            + _canonical_json(self.to_metadata())
        ).hexdigest()

    def compress(
        self,
        *,
        global_tolerance: float,
        axis_order: tuple[int, ...] | None = None,
        max_rank: int | None = None,
    ) -> CompressionDecision:
        """Gate a reduced representation on the global norm and every declared projection."""
        return gate_compression(
            self.values,
            axis_names=self.axis_names,
            projections=self.projections,
            global_tolerance=global_tolerance,
            axis_order=axis_order,
            max_rank=max_rank,
        )


def _validate_entity_ordering(
    entity_axis: TensorAxisSpec,
    per_time_entity_ids: Sequence[Sequence[str]],
) -> None:
    canonical = entity_axis.labels or ()
    for time_index, raw_ids in enumerate(per_time_entity_ids):
        if isinstance(raw_ids, str):
            raise TypeError("entity ids for a time index must be a sequence of strings")
        ids = tuple(_nonempty(entity_id, what="entity id") for entity_id in raw_ids)
        if len(ids) != len(canonical):
            raise ValueError(
                f"accepted time index {time_index} supplies {len(ids)} entities but the "
                f"{entity_axis.name!r} axis declares {len(canonical)}"
            )
        if ids == canonical:
            continue
        if set(ids) == set(canonical):
            raise ValueError(
                f"entity ordering changed at accepted time index {time_index}; a positional row "
                "would silently re-bind to a different entity"
            )
        raise ValueError(
            f"accepted time index {time_index} names entities absent from the "
            f"{entity_axis.name!r} axis"
        )


def build_telemetry_archive(
    *,
    archive_id: str,
    receipts: Sequence[StepReceipt],
    per_attempt_arrays: Sequence[Any],
    per_attempt_entity_ids: Sequence[Sequence[str]],
    axes: Sequence[TensorAxisSpec],
    projections: Sequence[ScalarProjectionSpec],
    value_semantics: str,
    value_unit: str,
    time_axis_name: str = "accepted_time",
    time_axis_semantics: str = "accepted physical time after the accepted step",
) -> TelemetryArchive:
    """Assemble an archive from an attempted receipt chain and its per-attempt arrays.

    Arrays are supplied for EVERY attempted step, aligned with ``receipts``; the adapter drops the
    rejected ones, so the count of dropped attempts is measured rather than asserted.  Entity
    orderings are validated on the accepted steps only -- the archive makes claims about what it
    contains, and a rejected attempt contributes nothing to it.

    Args:
        archive_id: Stable archive identifier.
        receipts: Contiguous attempted-step receipt chain.
        per_attempt_arrays: One array per attempted step, shaped like the non-time axes.
        per_attempt_entity_ids: One entity ordering per attempted step.
        axes: The non-time axes; exactly one must have role ``ENTITY``.
        projections: Declared scalar projections; at least one is required.
        value_semantics: What one archive value means.
        value_unit: Unit of one archive value.
        time_axis_name: Name for the prepended accepted-time axis.
        time_axis_semantics: Semantics for the prepended accepted-time axis.

    Returns:
        The assembled archive with its recorded scalars.
    """
    chain = tuple(receipts)
    arrays = tuple(per_attempt_arrays)
    orderings = tuple(per_attempt_entity_ids)
    if len(arrays) != len(chain) or len(orderings) != len(chain):
        raise ValueError("one array and one entity ordering are required per attempted step")
    timeline = project_accepted_timeline(chain)

    other_axes = tuple(axes)
    if any(not isinstance(axis, TensorAxisSpec) for axis in other_axes):
        raise TypeError("axes must be TensorAxisSpec values")
    if any(axis.role is AxisRole.ACCEPTED_TIME for axis in other_axes):
        raise ValueError("the accepted-time axis is derived from the receipt chain, not supplied")
    entity_axes = [axis for axis in other_axes if axis.role is AxisRole.ENTITY]
    if len(entity_axes) != 1:
        raise ValueError("exactly one stable-entity-id axis is required")

    accepted_positions = [index for index, receipt in enumerate(chain) if receipt.accepted]
    _validate_entity_ordering(
        entity_axes[0], [orderings[position] for position in accepted_positions]
    )

    expected_shape = tuple(axis.length for axis in other_axes)
    slices: list[np.ndarray] = []
    for position in accepted_positions:
        step_array = _finite_tensor(arrays[position])
        if tuple(step_array.shape) != expected_shape:
            raise ValueError(
                f"attempted step {position} array shape {tuple(step_array.shape)} does not match "
                f"the declared axes {expected_shape}"
            )
        slices.append(step_array)
    values = np.stack(slices, axis=0)

    time_axis = TensorAxisSpec(
        name=time_axis_name,
        role=AxisRole.ACCEPTED_TIME,
        semantics=time_axis_semantics,
        coordinates=timeline.accepted_time_s,
        unit="s",
    )
    all_axes = (time_axis, *other_axes)
    axis_names = tuple(axis.name for axis in all_axes)
    specs = tuple(projections)
    records = tuple(
        ProjectionRecord(
            name=spec.name,
            version=spec.version,
            definition_sha256=spec.definition_sha256,
            value=spec.evaluate(values, axis_names),
            unit=spec.unit,
        )
        for spec in specs
    )
    return TelemetryArchive(
        archive_id=archive_id,
        timeline=timeline,
        axes=all_axes,
        values=values,
        projections=specs,
        recorded_scalars=records,
        value_semantics=value_semantics,
        value_unit=value_unit,
    )


def compression_envelope(
    archive: TelemetryArchive,
    decision: CompressionDecision,
    *,
    artifact_id: str,
    cell_state_manifest_hash: str,
    source_artifact_ids: Sequence[str],
    fallback: str = "dense archive retained beside the reduced form",
) -> ArtifactEnvelope:
    """Wrap an ACCEPTED tensor-train compression in the shared provenance envelope.

    Raises:
        ValueError: If the compression was rejected (a rejected reduction has no artifact to
            describe), or if the accepted representation is a matrix SVD.  ``RepresentationKind``
            currently has no matrix-SVD member, and labelling an SVD ``TENSOR_TRAIN`` would be a
            false representation claim; adding the member is a ``contracts.py`` change owned
            elsewhere and is a PI decision, not something this adapter may improvise.
    """
    if not decision.accepted:
        raise ValueError("a rejected compression has no reduced artifact to describe")
    if decision.used_kind != "tensor-train":
        raise ValueError(
            "RepresentationKind has no matrix-SVD member; refusing to label an SVD as a tensor "
            "train. Adding the member is a contracts.py change and a PI decision."
        )
    reduced_bytes = serialize_archive_bytes(
        _representation_arrays(decision.reduction),
        {**archive.to_metadata(), "representation": _representation_metadata(decision.reduction)},
    )
    report = ApproximationErrorReport(
        report_sha256=_sha256_json(decision.to_dict()),
        source_blob_sha256=hashlib.sha256(archive.dense_payload()).hexdigest(),
        reduced_blob_sha256=hashlib.sha256(reduced_bytes).hexdigest(),
        metric_id=ERROR_METRIC_ID,
        requested_tolerance=decision.global_tolerance,
        measured_error=decision.measured_global_relative_error,
        projection_report_sha256=_sha256_json(dict(sorted(decision.projection_errors.items()))),
    )
    representation = RepresentationDescriptor(
        kind=RepresentationKind.TENSOR_TRAIN,
        uncertainty_sources=(UncertaintySource.REPRESENTATION, UncertaintySource.NUMERICAL),
        axes=archive.axis_names,
        approximation_error_bound=decision.measured_global_relative_error,
        error_metric=ERROR_METRIC_ID,
        error_report=report,
        fallback=fallback,
    )
    return ArtifactEnvelope(
        artifact_id=artifact_id,
        cell_state_manifest_hash=cell_state_manifest_hash,
        evidence_source=EvidenceSource.REDUCED_MODEL,
        representation=representation,
        source_artifact_ids=tuple(source_artifact_ids),
    )

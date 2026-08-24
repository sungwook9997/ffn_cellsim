"""Unverified observer receipt for attempted/accepted-step readback.

The reference runtime owns the PI-gated admissibility predicate.  This module never derives the
verdict and cannot grant native evidence authority.  It only validates the internal consistency of a
runtime-supplied claim and binds that claim to content hashes.  A future runtime adapter must verify
those hashes against the authoritative device transaction before any evidence promotion.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType

__all__ = [
    "StepReceipt",
    "predicate_vector_sha256",
    "state_coverage_sha256",
    "validate_step_receipt_chain",
]

_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _nonempty(value: str, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _digest(value: str, *, what: str) -> str:
    normalized = _nonempty(value, what=what)
    if not _SHA256.fullmatch(normalized):
        raise ValueError(f"{what} must be a lowercase SHA-256 digest")
    return normalized


def _finite_nonnegative(value: float, *, what: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0.0
    ):
        raise ValueError(f"{what} must be finite and nonnegative")
    return float(value)


def state_coverage_sha256(fields: Sequence[str]) -> str:
    """Hash the ordered names claimed to be included in the authoritative state digest."""
    if isinstance(fields, str) or not isinstance(fields, Sequence):
        raise TypeError("state coverage must be a sequence of strings")
    normalized = tuple(_nonempty(field, what="state coverage field") for field in fields)
    if not normalized:
        raise ValueError("state coverage must not be empty")
    if len(set(normalized)) != len(normalized):
        raise ValueError("state coverage fields must be unique")
    payload = json.dumps(normalized, separators=(",", ":"), ensure_ascii=False).encode()
    return hashlib.sha256(payload).hexdigest()


def predicate_vector_sha256(predicates: Mapping[str, bool]) -> str:
    """Hash normalized predicate names and values without interpreting their physics."""
    if not isinstance(predicates, Mapping) or not predicates:
        raise ValueError("predicate_results must be a non-empty mapping")
    normalized: dict[str, bool] = {}
    for raw_name, result in predicates.items():
        name = _nonempty(raw_name, what="predicate name")
        if name in normalized:
            raise ValueError("predicate names must not collide after whitespace normalization")
        if not isinstance(result, bool):
            raise TypeError(f"predicate_results[{name!r}] must be bool")
        normalized[name] = result
    payload = json.dumps(
        sorted(normalized.items()),
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode()
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True, slots=True)
class StepReceipt:
    """Internally consistent but externally unverified claim for one attempted physical step."""

    run_id: str
    decision_contract_sha256: str
    source_run_manifest_sha256: str
    build_stamp_sha256: str
    attempted_step_index: int
    accepted_step_before: int
    accepted_step_after: int
    accepted_time_before_s: float
    accepted_time_after_s: float
    attempted_dt_s: float
    accepted: bool
    predicate_results: Mapping[str, bool]
    predicate_results_sha256: str
    state_coverage: tuple[str, ...]
    state_coverage_sha256: str
    state_hash_before: str
    state_hash_after: str
    rng_hash_before: str
    rng_hash_after: str
    rejected_state_restored: bool | None
    authority: str = "unverified-observer"

    def __post_init__(self) -> None:
        object.__setattr__(self, "run_id", _nonempty(self.run_id, what="run_id"))
        for field_name in (
            "decision_contract_sha256",
            "source_run_manifest_sha256",
            "build_stamp_sha256",
            "predicate_results_sha256",
            "state_coverage_sha256",
            "state_hash_before",
            "state_hash_after",
            "rng_hash_before",
            "rng_hash_after",
        ):
            object.__setattr__(
                self,
                field_name,
                _digest(getattr(self, field_name), what=field_name),
            )
        if self.authority != "unverified-observer":
            raise ValueError("generic StepReceipt authority is fixed to unverified-observer")

        for field_name in ("attempted_step_index", "accepted_step_before", "accepted_step_after"):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{field_name} must be a nonnegative integer")
        if self.accepted_step_before > self.attempted_step_index:
            raise ValueError("accepted_step_before cannot exceed attempted_step_index")
        if self.accepted_step_after > self.attempted_step_index + 1:
            raise ValueError("accepted_step_after cannot exceed the attempted-step count")

        for field_name in (
            "accepted_time_before_s",
            "accepted_time_after_s",
            "attempted_dt_s",
        ):
            object.__setattr__(
                self,
                field_name,
                _finite_nonnegative(getattr(self, field_name), what=field_name),
            )
        if self.attempted_dt_s == 0.0:
            raise ValueError("attempted_dt_s must be positive")
        if not isinstance(self.accepted, bool):
            raise TypeError("accepted must be bool")

        predicates: dict[str, bool] = {}
        if not isinstance(self.predicate_results, Mapping) or not self.predicate_results:
            raise ValueError("predicate_results must be a non-empty mapping")
        for raw_name, result in self.predicate_results.items():
            name = _nonempty(raw_name, what="predicate name")
            if name in predicates:
                raise ValueError("predicate names must not collide after whitespace normalization")
            if not isinstance(result, bool):
                raise TypeError(f"predicate_results[{name!r}] must be bool")
            predicates[name] = result
        if predicate_vector_sha256(predicates) != self.predicate_results_sha256:
            raise ValueError("predicate_results_sha256 does not bind the supplied predicate vector")
        object.__setattr__(self, "predicate_results", MappingProxyType(predicates))

        coverage = tuple(
            _nonempty(field, what="state coverage field") for field in self.state_coverage
        )
        if not coverage:
            raise ValueError("state_coverage must not be empty")
        if len(set(coverage)) != len(coverage):
            raise ValueError("state_coverage fields must be unique")
        if state_coverage_sha256(coverage) != self.state_coverage_sha256:
            raise ValueError("state_coverage_sha256 does not bind the supplied coverage manifest")
        object.__setattr__(self, "state_coverage", coverage)

        if self.accepted:
            if self.rejected_state_restored is not None:
                raise ValueError("accepted receipt must not claim rejected-state restoration")
            if self.accepted_step_after != self.accepted_step_before + 1:
                raise ValueError("an accepted step must advance accepted_step by exactly one")
            expected_time = self.accepted_time_before_s + self.attempted_dt_s
            if self.accepted_time_after_s != expected_time:
                raise ValueError("an accepted step must advance accepted time by attempted_dt_s")
        else:
            if self.rejected_state_restored is not True:
                raise ValueError("a rejected step needs restored-state observer evidence")
            if self.accepted_step_after != self.accepted_step_before:
                raise ValueError("a rejected step must not advance accepted_step")
            if self.accepted_time_after_s != self.accepted_time_before_s:
                raise ValueError("a rejected step must not advance accepted time")
            if self.state_hash_after != self.state_hash_before:
                raise ValueError("a rejected step must restore the claimed state digest")
            if self.rng_hash_after != self.rng_hash_before:
                raise ValueError("a rejected step must restore the claimed RNG digest")


def validate_step_receipt_chain(receipts: Sequence[StepReceipt]) -> tuple[StepReceipt, ...]:
    """Validate continuity of one observer receipt chain without promoting its authority."""
    chain = tuple(receipts)
    if not chain:
        raise ValueError("receipt chain must not be empty")
    if any(not isinstance(receipt, StepReceipt) for receipt in chain):
        raise TypeError("receipt chain must contain only StepReceipt values")
    first = chain[0]
    invariant_fields = (
        "run_id",
        "decision_contract_sha256",
        "source_run_manifest_sha256",
        "build_stamp_sha256",
        "state_coverage_sha256",
    )
    for previous, current in zip(chain, chain[1:], strict=False):
        if current.attempted_step_index != previous.attempted_step_index + 1:
            raise ValueError("receipt chain has a missing, duplicate, or reordered attempted step")
        for field_name in invariant_fields:
            if getattr(current, field_name) != getattr(first, field_name):
                raise ValueError(f"receipt chain changed invariant field {field_name}")
        if current.accepted_step_before != previous.accepted_step_after:
            raise ValueError("receipt chain accepted-step clock is discontinuous")
        if current.accepted_time_before_s != previous.accepted_time_after_s:
            raise ValueError("receipt chain accepted physical time is discontinuous")
        if current.state_hash_before != previous.state_hash_after:
            raise ValueError("receipt chain state digest is discontinuous")
        if current.rng_hash_before != previous.rng_hash_after:
            raise ValueError("receipt chain RNG digest is discontinuous")
    return chain

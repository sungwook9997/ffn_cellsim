"""Content-addressed sidecar contract for lossless tensor and image artifacts."""

from __future__ import annotations

import hashlib
import io
import re
from dataclasses import dataclass

import numpy as np

__all__ = ["ArrayAxis", "ArraySidecarDescriptor"]

_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _nonempty(value: str, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


@dataclass(frozen=True, slots=True)
class ArrayAxis:
    """One ordered tensor axis and the hash of its coordinate vector, when applicable."""

    name: str
    length: int
    semantics: str
    unit: str = "1"
    coordinates_sha256: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _nonempty(self.name, what="axis name"))
        object.__setattr__(self, "semantics", _nonempty(self.semantics, what="axis semantics"))
        object.__setattr__(self, "unit", _nonempty(self.unit, what="axis unit"))
        if isinstance(self.length, bool) or not isinstance(self.length, int) or self.length < 1:
            raise ValueError("axis length must be a positive integer")
        if self.coordinates_sha256 is not None and not _SHA256.fullmatch(self.coordinates_sha256):
            raise ValueError("coordinates_sha256 must be a lowercase SHA-256 digest")


@dataclass(frozen=True, slots=True)
class ArraySidecarDescriptor:
    """Unverified metadata claim for an array blob kept outside JSON run records."""

    uri: str
    sha256: str
    media_type: str
    schema_id: str
    dtype: str
    value_semantics: str
    value_unit: str
    axes: tuple[ArrayAxis, ...]
    byte_size: int
    array_order: str = "C"

    def __post_init__(self) -> None:
        for field_name in (
            "uri",
            "media_type",
            "schema_id",
            "dtype",
            "value_semantics",
            "value_unit",
        ):
            object.__setattr__(
                self,
                field_name,
                _nonempty(getattr(self, field_name), what=field_name),
            )
        if not _SHA256.fullmatch(self.sha256):
            raise ValueError("sha256 must be a lowercase SHA-256 digest")
        if (
            isinstance(self.byte_size, bool)
            or not isinstance(self.byte_size, int)
            or self.byte_size < 0
        ):
            raise ValueError("byte_size must be a nonnegative integer")
        if self.array_order not in {"C", "F"}:
            raise ValueError("array_order must be 'C' or 'F'")
        try:
            canonical_dtype = np.dtype(self.dtype).str
        except TypeError as error:
            raise ValueError(f"dtype is not understood by NumPy: {self.dtype!r}") from error
        if self.dtype != canonical_dtype:
            raise ValueError(
                f"dtype must include canonical width and byte order; expected {canonical_dtype!r}"
            )
        axes = tuple(self.axes)
        if not axes or any(not isinstance(axis, ArrayAxis) for axis in axes):
            raise ValueError("axes must contain at least one ArrayAxis")
        names = tuple(axis.name for axis in axes)
        if len(set(names)) != len(names):
            raise ValueError("axis names must be unique and ordered")
        object.__setattr__(self, "axes", axes)

    @property
    def shape(self) -> tuple[int, ...]:
        return tuple(axis.length for axis in self.axes)

    def verify_content_address(self, payload: bytes) -> None:
        """Verify only byte count and content hash; this does not validate array semantics."""
        if not isinstance(payload, bytes):
            raise TypeError("payload must be bytes")
        if len(payload) != self.byte_size:
            raise ValueError("sidecar byte_size mismatch")
        actual_hash = hashlib.sha256(payload).hexdigest()
        if actual_hash != self.sha256:
            raise ValueError("sidecar SHA-256 mismatch")

    def verify_npy_array(self, payload: bytes) -> np.ndarray:
        """Decode an NPY payload without pickle and check dtype, shape, and memory order."""
        self.verify_content_address(payload)
        if self.media_type != "application/x-npy":
            raise ValueError("decoded array verification currently requires application/x-npy")
        try:
            decoded = np.load(io.BytesIO(payload), allow_pickle=False)
        except (OSError, ValueError) as error:
            raise ValueError("sidecar is not a safe, valid NPY array") from error
        if not isinstance(decoded, np.ndarray):
            raise ValueError("sidecar decoded to a container rather than one NPY array")
        if decoded.shape != self.shape:
            raise ValueError(
                f"decoded sidecar shape {decoded.shape} does not match axes {self.shape}"
            )
        if decoded.dtype.str != self.dtype:
            raise ValueError(
                f"decoded sidecar dtype {decoded.dtype.str!r} does not match {self.dtype!r}"
            )
        order_matches = (
            decoded.flags.c_contiguous if self.array_order == "C" else decoded.flags.f_contiguous
        )
        if not order_matches:
            raise ValueError(
                f"decoded sidecar is not contiguous in declared {self.array_order} order"
            )
        decoded.setflags(write=False)
        return decoded

from __future__ import annotations

import hashlib
import io

import numpy as np
import pytest

from aleph.virtual_cell.probability import audit_probability_mass
from aleph.virtual_cell.sidecar import ArrayAxis, ArraySidecarDescriptor


def _sidecar(payload: bytes) -> ArraySidecarDescriptor:
    return ArraySidecarDescriptor(
        uri="objects/ab/tensor.npz",
        sha256=hashlib.sha256(payload).hexdigest(),
        media_type="application/x-npz",
        schema_id="ffn.virtual-cell/array-sidecar@1",
        dtype="<f8",
        value_semantics="oracle scalar value",
        value_unit="1",
        axes=(
            ArrayAxis(name="time", length=3, semantics="accepted physical time", unit="s"),
            ArrayAxis(name="node", length=4, semantics="stable node identity"),
        ),
        byte_size=len(payload),
    )


def test_sidecar_preserves_ordered_axis_manifest_and_verifies_content() -> None:
    payload = b"small-gold-fixture"
    sidecar = _sidecar(payload)
    assert sidecar.shape == (3, 4)
    sidecar.verify_content_address(payload)


def test_sidecar_detects_hash_corruption() -> None:
    sidecar = _sidecar(b"expected")
    with pytest.raises(ValueError, match="SHA-256"):
        sidecar.verify_content_address(b"tampered")


def test_sidecar_detects_size_corruption_before_hash() -> None:
    sidecar = _sidecar(b"expected")
    with pytest.raises(ValueError, match="byte_size"):
        sidecar.verify_content_address(b"longer-than-expected")


def test_npy_sidecar_validates_decoded_shape_dtype_and_order() -> None:
    array = np.arange(12.0, dtype=np.float64).reshape(3, 4)
    buffer = io.BytesIO()
    np.save(buffer, array, allow_pickle=False)
    payload = buffer.getvalue()
    descriptor = ArraySidecarDescriptor(
        uri="objects/tensor.npy",
        sha256=hashlib.sha256(payload).hexdigest(),
        media_type="application/x-npy",
        schema_id="ffn.virtual-cell/array-sidecar@1",
        dtype=array.dtype.str,
        value_semantics="oracle scalar value",
        value_unit="1",
        axes=(
            ArrayAxis(name="time", length=3, semantics="accepted time"),
            ArrayAxis(name="node", length=4, semantics="stable node ID"),
        ),
        byte_size=len(payload),
        array_order="C",
    )
    decoded = descriptor.verify_npy_array(payload)
    np.testing.assert_array_equal(decoded, array)
    assert not decoded.flags.writeable


def test_npy_sidecar_rejects_descriptor_shape_mismatch() -> None:
    array = np.arange(12.0, dtype=np.float64).reshape(3, 4)
    buffer = io.BytesIO()
    np.save(buffer, array, allow_pickle=False)
    payload = buffer.getvalue()
    descriptor = ArraySidecarDescriptor(
        uri="objects/tensor.npy",
        sha256=hashlib.sha256(payload).hexdigest(),
        media_type="application/x-npy",
        schema_id="ffn.virtual-cell/array-sidecar@1",
        dtype=array.dtype.str,
        value_semantics="oracle scalar value",
        value_unit="1",
        axes=(
            ArrayAxis(name="time", length=2, semantics="accepted time"),
            ArrayAxis(name="node", length=6, semantics="stable node ID"),
        ),
        byte_size=len(payload),
    )
    with pytest.raises(ValueError, match="shape"):
        descriptor.verify_npy_array(payload)


def test_probability_mass_audit_accepts_normalized_mass_without_mutation() -> None:
    mass = np.array([[0.1, 0.2], [0.3, 0.4]])
    original = mass.copy()
    audit = audit_probability_mass(mass)
    assert audit.total_mass == pytest.approx(1.0)
    assert audit.negativity_tolerance == 0.0
    assert audit.measure_kind == "discrete-probability-mass"
    np.testing.assert_array_equal(mass, original)


def test_probability_mass_audit_rejects_negative_values_without_clipping() -> None:
    with pytest.raises(ValueError, match="negative"):
        audit_probability_mass([0.6, 0.5, -0.1])


def test_probability_mass_audit_rejects_silent_normalization() -> None:
    with pytest.raises(ValueError, match="not normalized"):
        audit_probability_mass([2.0, 3.0])


def test_probability_mass_audit_rejects_nonfinite_values() -> None:
    with pytest.raises(ValueError, match="finite"):
        audit_probability_mass([0.5, np.nan])

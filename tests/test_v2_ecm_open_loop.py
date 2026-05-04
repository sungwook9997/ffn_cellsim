from __future__ import annotations

import numpy as np
import pytest

from acs.v2.dynamics.ecm_open_loop import (
    ECMOpenLoopError,
    accumulate_prescribed_traction,
)
from acs.v2.ecm_substrate import ECMSubstrateState


def _ecm(nx: int = 3, ny: int = 2) -> ECMSubstrateState:
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    return ECMSubstrateState(
        origin_um_xy=(1.0, -2.0),
        spacing_um=0.5,
        stiffness_kpa=np.full((nx, ny), 2.0, dtype=np.float64),
        ligand_density=np.full((nx, ny), 0.7, dtype=np.float64),
        fiber_density=np.full((nx, ny), 0.3, dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.full((nx, ny), 0.25, dtype=np.float64),
        source="synthetic",
    )


def test_zero_traction_no_op_but_returns_new_state():
    ecm = _ecm()
    before = ecm.accumulated_traction_nNs_per_um2.copy()
    out = accumulate_prescribed_traction(ecm, np.zeros(ecm.grid_shape), dt_s=3.0)
    assert out is not ecm
    np.testing.assert_array_equal(out.accumulated_traction_nNs_per_um2, before)
    np.testing.assert_array_equal(ecm.accumulated_traction_nNs_per_um2, before)


def test_zero_dt_no_op():
    ecm = _ecm()
    traction = np.full(ecm.grid_shape, 4.0, dtype=np.float64)
    out = accumulate_prescribed_traction(ecm, traction, dt_s=0.0)
    np.testing.assert_array_equal(
        out.accumulated_traction_nNs_per_um2,
        ecm.accumulated_traction_nNs_per_um2,
    )


def test_uniform_positive_traction_exact_increment():
    ecm = _ecm()
    out = accumulate_prescribed_traction(
        ecm,
        np.full(ecm.grid_shape, 2.0, dtype=np.float64),
        dt_s=0.125,
    )
    expected = ecm.accumulated_traction_nNs_per_um2 + 0.25
    np.testing.assert_allclose(out.accumulated_traction_nNs_per_um2, expected)


def test_nonuniform_positive_traction_exact_increment():
    ecm = _ecm()
    traction = np.arange(6, dtype=np.float64).reshape(ecm.grid_shape) / 10.0
    out = accumulate_prescribed_traction(ecm, traction, dt_s=2.0)
    expected = ecm.accumulated_traction_nNs_per_um2 + traction * 2.0
    np.testing.assert_allclose(out.accumulated_traction_nNs_per_um2, expected)


def test_other_ecm_fields_are_preserved_and_copied():
    ecm = _ecm()
    out = accumulate_prescribed_traction(ecm, np.ones(ecm.grid_shape), dt_s=1.0)
    assert out.origin_um_xy == ecm.origin_um_xy
    assert out.spacing_um == ecm.spacing_um
    assert out.source == ecm.source
    for name in (
        "stiffness_kpa",
        "ligand_density",
        "fiber_density",
        "orientation_tensor",
    ):
        original = getattr(ecm, name)
        copied = getattr(out, name)
        np.testing.assert_array_equal(copied, original)
        assert copied is not original


def test_returned_state_validates():
    ecm = _ecm()
    out = accumulate_prescribed_traction(ecm, np.ones(ecm.grid_shape), dt_s=1.0)
    out.validate()


def test_shape_mismatch_rejected():
    ecm = _ecm()
    with pytest.raises(ECMOpenLoopError) as info:
        accumulate_prescribed_traction(ecm, np.ones((2, 2)), dt_s=1.0)
    assert info.value.failure_kind == "traction_shape"


def test_negative_traction_rejected():
    ecm = _ecm()
    traction = np.zeros(ecm.grid_shape)
    traction[0, 0] = -1.0
    with pytest.raises(ECMOpenLoopError) as info:
        accumulate_prescribed_traction(ecm, traction, dt_s=1.0)
    assert info.value.failure_kind == "traction_negative"


def test_non_finite_traction_rejected():
    ecm = _ecm()
    traction = np.zeros(ecm.grid_shape)
    traction[0, 0] = np.nan
    with pytest.raises(ECMOpenLoopError) as info:
        accumulate_prescribed_traction(ecm, traction, dt_s=1.0)
    assert info.value.failure_kind == "traction_non_finite"


@pytest.mark.parametrize("bad_dt", [True, -1.0, np.inf, np.nan, object()])
def test_invalid_dt_rejected(bad_dt):
    ecm = _ecm()
    with pytest.raises(ECMOpenLoopError) as info:
        accumulate_prescribed_traction(ecm, np.zeros(ecm.grid_shape), dt_s=bad_dt)
    assert info.value.failure_kind == "dt_invalid"


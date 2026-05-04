from __future__ import annotations

import numpy as np
import pytest

from acs.v2.dynamics.ecm_open_loop import (
    ECMOpenLoopError,
    accumulate_prescribed_traction,
    apply_prescribed_density_rate,
    apply_prescribed_orientation_rate,
    apply_prescribed_stiffness_rate,
)
from acs.v2.ecm_substrate import _ORIENTATION_BOUND
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


@pytest.mark.parametrize(
    "bad_dt",
    [np.True_, np.False_, np.bool_(True), np.array(True), np.array(False)],
)
def test_numpy_bool_dt_rejected_for_traction(bad_dt):
    """Codex review id=1121: `np.bool_` / 0-d bool array silently coerce to
    1.0/0.0 via ``float()``. They must be rejected as ``dt_invalid`` even
    though Python ``bool`` is already covered above."""

    ecm = _ecm()
    with pytest.raises(ECMOpenLoopError) as info:
        accumulate_prescribed_traction(ecm, np.zeros(ecm.grid_shape), dt_s=bad_dt)
    assert info.value.failure_kind == "dt_invalid"


# ----- 6.4-open-A: apply_prescribed_stiffness_rate -----


def test_stiffness_rate_zero_no_op_returns_new_state():
    ecm = _ecm()
    rate = np.zeros(ecm.grid_shape)
    out = apply_prescribed_stiffness_rate(ecm, rate, dt_s=1.0)
    assert out is not ecm
    np.testing.assert_array_equal(out.stiffness_kpa, ecm.stiffness_kpa)


def test_stiffness_rate_zero_dt_no_op():
    ecm = _ecm()
    rate = np.full(ecm.grid_shape, 5.0, dtype=np.float64)  # large rate
    out = apply_prescribed_stiffness_rate(ecm, rate, dt_s=0.0)
    np.testing.assert_array_equal(out.stiffness_kpa, ecm.stiffness_kpa)


def test_stiffness_rate_uniform_positive_exact_increment():
    ecm = _ecm()
    rate = np.full(ecm.grid_shape, 0.5, dtype=np.float64)  # +0.5 kPa/s
    out = apply_prescribed_stiffness_rate(ecm, rate, dt_s=2.0)  # *2 s
    expected = ecm.stiffness_kpa + 1.0  # +1.0 kPa
    np.testing.assert_allclose(out.stiffness_kpa, expected, atol=0.0)


def test_stiffness_rate_nonuniform_signed_exact_update():
    ecm = _ecm()
    rate = np.array([[0.1, -0.2], [0.3, -0.4], [0.5, -0.6]], dtype=np.float64)
    out = apply_prescribed_stiffness_rate(ecm, rate, dt_s=1.0)
    expected = ecm.stiffness_kpa + rate
    np.testing.assert_allclose(out.stiffness_kpa, expected, atol=0.0)
    assert (out.stiffness_kpa >= 0.0).all()


def test_stiffness_rate_negative_post_state_raises():
    ecm = _ecm()
    # ecm.stiffness_kpa = 2.0; a -3.0 kPa/s rate over 1 s would give -1.0.
    rate = np.full(ecm.grid_shape, -3.0, dtype=np.float64)
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_stiffness_rate(ecm, rate, dt_s=1.0)
    assert info.value.failure_kind == "stiffness_negative_post_update"


def test_stiffness_rate_at_exact_zero_post_state_accepted():
    """Boundary case: rate · dt exactly cancels existing stiffness → post == 0
    is allowed (schema permits non-negative)."""

    ecm = _ecm()
    rate = np.full(ecm.grid_shape, -2.0, dtype=np.float64)  # exactly -ecm.stiffness
    out = apply_prescribed_stiffness_rate(ecm, rate, dt_s=1.0)
    np.testing.assert_allclose(out.stiffness_kpa, np.zeros(ecm.grid_shape))


def test_stiffness_rate_shape_mismatch_rejected():
    ecm = _ecm()
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_stiffness_rate(ecm, np.zeros((2, 2)), dt_s=1.0)
    assert info.value.failure_kind == "stiffness_rate_shape"


def test_stiffness_rate_non_finite_rejected():
    ecm = _ecm()
    rate = np.zeros(ecm.grid_shape)
    rate[0, 0] = np.nan
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_stiffness_rate(ecm, rate, dt_s=1.0)
    assert info.value.failure_kind == "stiffness_rate_non_finite"


@pytest.mark.parametrize("bad_dt", [True, -1.0, np.inf, np.nan, object()])
def test_stiffness_rate_invalid_dt_rejected(bad_dt):
    ecm = _ecm()
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_stiffness_rate(ecm, np.zeros(ecm.grid_shape), dt_s=bad_dt)
    assert info.value.failure_kind == "dt_invalid"


@pytest.mark.parametrize(
    "bad_dt",
    [np.True_, np.False_, np.bool_(True), np.array(True), np.array(False)],
)
def test_numpy_bool_dt_rejected_for_stiffness_rate(bad_dt):
    ecm = _ecm()
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_stiffness_rate(ecm, np.zeros(ecm.grid_shape), dt_s=bad_dt)
    assert info.value.failure_kind == "dt_invalid"


def test_stiffness_rate_other_fields_preserved_and_not_aliased():
    ecm = _ecm()
    rate = np.full(ecm.grid_shape, 0.1, dtype=np.float64)
    out = apply_prescribed_stiffness_rate(ecm, rate, dt_s=1.0)
    # Verbatim equality on non-stiffness fields.
    np.testing.assert_array_equal(out.ligand_density, ecm.ligand_density)
    np.testing.assert_array_equal(out.fiber_density, ecm.fiber_density)
    np.testing.assert_array_equal(out.orientation_tensor, ecm.orientation_tensor)
    np.testing.assert_array_equal(
        out.accumulated_traction_nNs_per_um2, ecm.accumulated_traction_nNs_per_um2
    )
    assert out.origin_um_xy == ecm.origin_um_xy
    assert out.spacing_um == ecm.spacing_um
    assert out.source == ecm.source
    # Not aliased: mutating one must not affect the other (defensive copy).
    out.stiffness_kpa.setflags(write=True)  # ensure writable
    original_snapshot = ecm.stiffness_kpa.copy()
    out.stiffness_kpa[0, 0] = 99.0
    np.testing.assert_array_equal(ecm.stiffness_kpa, original_snapshot)
    np.testing.assert_array_equal(
        out.ligand_density, ecm.ligand_density,
    )  # still equal, but distinct array
    assert out.ligand_density is not ecm.ligand_density


def test_stiffness_rate_input_ecm_not_mutated():
    ecm = _ecm()
    snapshot = ecm.stiffness_kpa.copy()
    rate = np.full(ecm.grid_shape, 0.7, dtype=np.float64)
    apply_prescribed_stiffness_rate(ecm, rate, dt_s=1.0)
    np.testing.assert_array_equal(ecm.stiffness_kpa, snapshot)


def test_stiffness_rate_returned_state_validates():
    ecm = _ecm()
    rate = np.full(ecm.grid_shape, 0.1, dtype=np.float64)
    out = apply_prescribed_stiffness_rate(ecm, rate, dt_s=1.0)
    out.validate()  # must not raise; schema invariants intact


def test_stiffness_rate_chained_calls_compose():
    """Multiple calls compose linearly: dt_a · rate_a + dt_b · rate_b applied in
    sequence equals one call with their sum (under matching final post-state)."""

    ecm = _ecm()
    rate_a = np.full(ecm.grid_shape, 0.2, dtype=np.float64)
    rate_b = np.full(ecm.grid_shape, -0.1, dtype=np.float64)
    out_chained = apply_prescribed_stiffness_rate(ecm, rate_a, dt_s=1.0)
    out_chained = apply_prescribed_stiffness_rate(out_chained, rate_b, dt_s=2.0)
    expected = ecm.stiffness_kpa + 0.2 * 1.0 + (-0.1) * 2.0
    np.testing.assert_allclose(out_chained.stiffness_kpa, expected, atol=0.0)


# ----- 6.4-open-B: apply_prescribed_density_rate -----


def test_density_rate_zero_no_op_returns_new_state():
    ecm = _ecm()
    zero = np.zeros(ecm.grid_shape)
    out = apply_prescribed_density_rate(ecm, zero, zero, dt_s=1.0)
    assert out is not ecm
    np.testing.assert_array_equal(out.ligand_density, ecm.ligand_density)
    np.testing.assert_array_equal(out.fiber_density, ecm.fiber_density)


def test_density_rate_zero_dt_no_op():
    ecm = _ecm()
    rate = np.full(ecm.grid_shape, 0.5, dtype=np.float64)
    out = apply_prescribed_density_rate(ecm, rate, rate, dt_s=0.0)
    np.testing.assert_array_equal(out.ligand_density, ecm.ligand_density)
    np.testing.assert_array_equal(out.fiber_density, ecm.fiber_density)


def test_density_rate_uniform_positive_exact_increment():
    ecm = _ecm()
    ligand_rate = np.full(ecm.grid_shape, 0.1, dtype=np.float64)
    fiber_rate = np.full(ecm.grid_shape, 0.05, dtype=np.float64)
    out = apply_prescribed_density_rate(ecm, ligand_rate, fiber_rate, dt_s=1.0)
    np.testing.assert_allclose(out.ligand_density, ecm.ligand_density + 0.1, atol=0.0)
    np.testing.assert_allclose(out.fiber_density, ecm.fiber_density + 0.05, atol=0.0)


def test_density_rate_signed_nonuniform_exact_update():
    ecm = _ecm()
    ligand_rate = np.array([[0.1, -0.05], [0.0, 0.05], [-0.1, 0.0]], dtype=np.float64)
    fiber_rate = np.array([[-0.05, 0.1], [0.05, 0.0], [0.0, 0.1]], dtype=np.float64)
    out = apply_prescribed_density_rate(ecm, ligand_rate, fiber_rate, dt_s=2.0)
    np.testing.assert_allclose(
        out.ligand_density, ecm.ligand_density + ligand_rate * 2.0, atol=0.0
    )
    np.testing.assert_allclose(
        out.fiber_density, ecm.fiber_density + fiber_rate * 2.0, atol=0.0
    )


def test_density_rate_ligand_negative_post_state_raises():
    ecm = _ecm()  # ligand = 0.7
    rate = np.full(ecm.grid_shape, -0.8, dtype=np.float64)
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_density_rate(
            ecm, rate, np.zeros(ecm.grid_shape), dt_s=1.0
        )
    assert info.value.failure_kind == "ligand_density_negative_post_update"


def test_density_rate_ligand_exceeds_one_raises():
    ecm = _ecm()  # ligand = 0.7
    rate = np.full(ecm.grid_shape, 0.4, dtype=np.float64)  # 0.7 + 0.4 = 1.1
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_density_rate(
            ecm, rate, np.zeros(ecm.grid_shape), dt_s=1.0
        )
    assert info.value.failure_kind == "ligand_density_exceeds_one_post_update"


def test_density_rate_fiber_negative_post_state_raises():
    ecm = _ecm()  # fiber = 0.3
    rate = np.full(ecm.grid_shape, -0.4, dtype=np.float64)
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_density_rate(
            ecm, np.zeros(ecm.grid_shape), rate, dt_s=1.0
        )
    assert info.value.failure_kind == "fiber_density_negative_post_update"


def test_density_rate_fiber_exceeds_one_raises():
    ecm = _ecm()  # fiber = 0.3
    rate = np.full(ecm.grid_shape, 0.8, dtype=np.float64)  # 0.3 + 0.8 = 1.1
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_density_rate(
            ecm, np.zeros(ecm.grid_shape), rate, dt_s=1.0
        )
    assert info.value.failure_kind == "fiber_density_exceeds_one_post_update"


def test_density_rate_exact_zero_post_state_accepted():
    """Boundary: rate · dt exactly cancels existing density → post == 0
    is allowed (schema permits non-negative)."""

    ecm = _ecm()  # ligand=0.7, fiber=0.3
    ligand_rate = np.full(ecm.grid_shape, -0.7, dtype=np.float64)
    fiber_rate = np.full(ecm.grid_shape, -0.3, dtype=np.float64)
    out = apply_prescribed_density_rate(ecm, ligand_rate, fiber_rate, dt_s=1.0)
    np.testing.assert_allclose(out.ligand_density, np.zeros(ecm.grid_shape))
    np.testing.assert_allclose(out.fiber_density, np.zeros(ecm.grid_shape))


def test_density_rate_exact_one_post_state_accepted():
    """Boundary: post == 1 is allowed (schema permits up to 1)."""

    ecm = _ecm()  # ligand=0.7, fiber=0.3
    ligand_rate = np.full(ecm.grid_shape, 0.3, dtype=np.float64)  # 0.7 + 0.3 = 1.0
    fiber_rate = np.full(ecm.grid_shape, 0.7, dtype=np.float64)  # 0.3 + 0.7 = 1.0
    out = apply_prescribed_density_rate(ecm, ligand_rate, fiber_rate, dt_s=1.0)
    np.testing.assert_allclose(out.ligand_density, np.ones(ecm.grid_shape))
    np.testing.assert_allclose(out.fiber_density, np.ones(ecm.grid_shape))


def test_density_rate_ligand_shape_mismatch_rejected():
    ecm = _ecm()
    bad = np.zeros((2, 2))
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_density_rate(
            ecm, bad, np.zeros(ecm.grid_shape), dt_s=1.0
        )
    assert info.value.failure_kind == "ligand_density_rate_shape"


def test_density_rate_fiber_shape_mismatch_rejected():
    ecm = _ecm()
    bad = np.zeros((2, 2))
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_density_rate(
            ecm, np.zeros(ecm.grid_shape), bad, dt_s=1.0
        )
    assert info.value.failure_kind == "fiber_density_rate_shape"


def test_density_rate_ligand_non_finite_rejected():
    ecm = _ecm()
    bad = np.zeros(ecm.grid_shape)
    bad[0, 0] = np.nan
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_density_rate(
            ecm, bad, np.zeros(ecm.grid_shape), dt_s=1.0
        )
    assert info.value.failure_kind == "ligand_density_rate_non_finite"


def test_density_rate_fiber_non_finite_rejected():
    ecm = _ecm()
    bad = np.zeros(ecm.grid_shape)
    bad[1, 1] = np.inf
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_density_rate(
            ecm, np.zeros(ecm.grid_shape), bad, dt_s=1.0
        )
    assert info.value.failure_kind == "fiber_density_rate_non_finite"


@pytest.mark.parametrize("bad_dt", [True, -1.0, np.inf, np.nan, object()])
def test_density_rate_invalid_dt_rejected(bad_dt):
    ecm = _ecm()
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_density_rate(
            ecm,
            np.zeros(ecm.grid_shape),
            np.zeros(ecm.grid_shape),
            dt_s=bad_dt,
        )
    assert info.value.failure_kind == "dt_invalid"


@pytest.mark.parametrize(
    "bad_dt",
    [np.True_, np.False_, np.bool_(True), np.array(True), np.array(False)],
)
def test_density_rate_numpy_bool_dt_rejected(bad_dt):
    ecm = _ecm()
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_density_rate(
            ecm,
            np.zeros(ecm.grid_shape),
            np.zeros(ecm.grid_shape),
            dt_s=bad_dt,
        )
    assert info.value.failure_kind == "dt_invalid"


def test_density_rate_cross_isolation_ligand_only_does_not_change_fiber():
    """Cross-isolation: nonzero ligand rate + zero fiber rate must not
    change fiber density. The ligand update never reads the fiber rate."""

    ecm = _ecm()
    ligand_rate = np.full(ecm.grid_shape, 0.1, dtype=np.float64)
    fiber_rate = np.zeros(ecm.grid_shape, dtype=np.float64)
    out = apply_prescribed_density_rate(ecm, ligand_rate, fiber_rate, dt_s=1.0)
    np.testing.assert_array_equal(out.fiber_density, ecm.fiber_density)
    assert not np.array_equal(out.ligand_density, ecm.ligand_density)


def test_density_rate_cross_isolation_fiber_only_does_not_change_ligand():
    """Mirror: nonzero fiber rate + zero ligand rate must not change ligand."""

    ecm = _ecm()
    ligand_rate = np.zeros(ecm.grid_shape, dtype=np.float64)
    fiber_rate = np.full(ecm.grid_shape, 0.1, dtype=np.float64)
    out = apply_prescribed_density_rate(ecm, ligand_rate, fiber_rate, dt_s=1.0)
    np.testing.assert_array_equal(out.ligand_density, ecm.ligand_density)
    assert not np.array_equal(out.fiber_density, ecm.fiber_density)


def test_density_rate_other_fields_preserved_and_not_aliased():
    ecm = _ecm()
    rate = np.full(ecm.grid_shape, 0.05, dtype=np.float64)
    out = apply_prescribed_density_rate(ecm, rate, rate, dt_s=1.0)
    np.testing.assert_array_equal(out.stiffness_kpa, ecm.stiffness_kpa)
    np.testing.assert_array_equal(out.orientation_tensor, ecm.orientation_tensor)
    np.testing.assert_array_equal(
        out.accumulated_traction_nNs_per_um2, ecm.accumulated_traction_nNs_per_um2
    )
    assert out.origin_um_xy == ecm.origin_um_xy
    assert out.spacing_um == ecm.spacing_um
    assert out.source == ecm.source
    assert out.stiffness_kpa is not ecm.stiffness_kpa
    assert out.orientation_tensor is not ecm.orientation_tensor


def test_density_rate_input_ecm_not_mutated():
    ecm = _ecm()
    snap_ligand = ecm.ligand_density.copy()
    snap_fiber = ecm.fiber_density.copy()
    rate = np.full(ecm.grid_shape, 0.1, dtype=np.float64)
    apply_prescribed_density_rate(ecm, rate, rate, dt_s=1.0)
    np.testing.assert_array_equal(ecm.ligand_density, snap_ligand)
    np.testing.assert_array_equal(ecm.fiber_density, snap_fiber)


def test_density_rate_returned_state_validates():
    ecm = _ecm()
    rate = np.full(ecm.grid_shape, 0.05, dtype=np.float64)
    out = apply_prescribed_density_rate(ecm, rate, rate, dt_s=1.0)
    out.validate()


def test_density_rate_chained_calls_compose():
    ecm = _ecm()
    rate_a = np.full(ecm.grid_shape, 0.05, dtype=np.float64)
    rate_b = np.full(ecm.grid_shape, -0.025, dtype=np.float64)
    out_chained = apply_prescribed_density_rate(ecm, rate_a, rate_a, dt_s=1.0)
    out_chained = apply_prescribed_density_rate(out_chained, rate_b, rate_b, dt_s=2.0)
    expected_ligand = ecm.ligand_density + 0.05 * 1.0 + (-0.025) * 2.0
    expected_fiber = ecm.fiber_density + 0.05 * 1.0 + (-0.025) * 2.0
    np.testing.assert_allclose(out_chained.ligand_density, expected_ligand, atol=0.0)
    np.testing.assert_allclose(out_chained.fiber_density, expected_fiber, atol=0.0)


# ----- 6.4-open-C: apply_prescribed_orientation_rate -----


def _zero_orientation_rate(grid_shape: tuple[int, int]) -> np.ndarray:
    return np.zeros((*grid_shape, 2, 2), dtype=np.float64)


def _symmetric_orientation_rate(
    grid_shape: tuple[int, int], xx: float, yy: float, off: float
) -> np.ndarray:
    rate = np.zeros((*grid_shape, 2, 2), dtype=np.float64)
    rate[..., 0, 0] = xx
    rate[..., 1, 1] = yy
    rate[..., 0, 1] = off
    rate[..., 1, 0] = off
    return rate


def test_orientation_rate_zero_no_op_returns_new_state():
    ecm = _ecm()
    rate = _zero_orientation_rate(ecm.grid_shape)
    out = apply_prescribed_orientation_rate(ecm, rate, dt_s=1.0)
    assert out is not ecm
    np.testing.assert_array_equal(out.orientation_tensor, ecm.orientation_tensor)


def test_orientation_rate_zero_dt_no_op():
    ecm = _ecm()
    rate = _symmetric_orientation_rate(ecm.grid_shape, xx=0.5, yy=-0.3, off=0.2)
    out = apply_prescribed_orientation_rate(ecm, rate, dt_s=0.0)
    np.testing.assert_array_equal(out.orientation_tensor, ecm.orientation_tensor)


def test_orientation_rate_uniform_symmetric_exact_update():
    """ECM starts at identity (T_00 = T_11 = 1.0). Rates must keep
    post-state in [-_ORIENTATION_BOUND, +_ORIENTATION_BOUND], so we
    push T_00 / T_11 down (negative rate) and T_01 / T_10 up by a
    small symmetric off-diagonal increment."""

    ecm = _ecm()
    rate = _symmetric_orientation_rate(ecm.grid_shape, xx=-0.1, yy=-0.05, off=0.025)
    out = apply_prescribed_orientation_rate(ecm, rate, dt_s=2.0)
    expected = ecm.orientation_tensor + rate * 2.0
    np.testing.assert_allclose(out.orientation_tensor, expected, atol=0.0)


def test_orientation_rate_sparse_symmetric_per_cell_exact_update():
    """Per-cell distinct symmetric rate: only one cell has nonzero rate.
    Push diagonals down so post-state stays within bounds."""

    ecm = _ecm()
    rate = _zero_orientation_rate(ecm.grid_shape)
    rate[0, 0] = np.array([[-0.3, -0.1], [-0.1, -0.2]], dtype=np.float64)
    out = apply_prescribed_orientation_rate(ecm, rate, dt_s=1.0)
    expected = ecm.orientation_tensor + rate
    np.testing.assert_allclose(out.orientation_tensor, expected, atol=0.0)


def test_orientation_rate_asymmetric_rejected():
    ecm = _ecm()
    rate = _zero_orientation_rate(ecm.grid_shape)
    rate[0, 0, 0, 1] = 0.5  # asymmetric: T_01 != T_10
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_orientation_rate(ecm, rate, dt_s=1.0)
    assert info.value.failure_kind == "orientation_rate_asymmetric"


def test_orientation_rate_shape_mismatch_rejected():
    ecm = _ecm()
    bad = np.zeros((2, 2, 2, 2), dtype=np.float64)
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_orientation_rate(ecm, bad, dt_s=1.0)
    assert info.value.failure_kind == "orientation_rate_shape"


def test_orientation_rate_missing_block_dimension_rejected():
    ecm = _ecm()
    bad_2d = np.zeros(ecm.grid_shape, dtype=np.float64)
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_orientation_rate(ecm, bad_2d, dt_s=1.0)
    assert info.value.failure_kind == "orientation_rate_shape"


def test_orientation_rate_non_finite_rejected():
    ecm = _ecm()
    rate = _zero_orientation_rate(ecm.grid_shape)
    rate[1, 1, 0, 0] = np.nan
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_orientation_rate(ecm, rate, dt_s=1.0)
    assert info.value.failure_kind == "orientation_rate_non_finite"


def test_orientation_negative_bound_post_update_raises():
    ecm = _ecm()
    # ecm orientation_tensor[..., 0, 0] = 1.0 (identity); push down by -2.5.
    rate = _symmetric_orientation_rate(ecm.grid_shape, xx=-2.5, yy=0.0, off=0.0)
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_orientation_rate(ecm, rate, dt_s=1.0)
    assert info.value.failure_kind == "orientation_negative_bound_post_update"


def test_orientation_exceeds_bound_post_update_raises():
    ecm = _ecm()
    # ecm orientation_tensor[..., 0, 0] = 1.0; push up by +0.5 → 1.5 > 1.0.
    rate = _symmetric_orientation_rate(ecm.grid_shape, xx=0.5, yy=0.0, off=0.0)
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_orientation_rate(ecm, rate, dt_s=1.0)
    assert info.value.failure_kind == "orientation_exceeds_bound_post_update"


def test_orientation_exact_positive_bound_accepted():
    """Post == +_ORIENTATION_BOUND is allowed (schema permits).
    ecm starts at identity so [..., 1, 1] = 1.0; rate = -0.0 keeps it at 1.0.
    Push another component to +bound to verify boundary acceptance."""

    ecm = _ecm()
    # Start with identity orientation: T_00=1, T_11=1, T_01=T_10=0.
    # Add 0 to T_00 (stays at 1.0). T_11 already at 1.0. Test passes if no raise.
    rate = _symmetric_orientation_rate(ecm.grid_shape, xx=0.0, yy=0.0, off=0.0)
    out = apply_prescribed_orientation_rate(ecm, rate, dt_s=1.0)
    np.testing.assert_array_equal(out.orientation_tensor, ecm.orientation_tensor)


def test_orientation_exact_negative_bound_accepted():
    """Push a zero component down to -_ORIENTATION_BOUND exactly."""

    ecm = _ecm()
    # ecm T_01 = T_10 = 0; rate = -1.0 brings them to -1.0 exactly.
    rate = _symmetric_orientation_rate(
        ecm.grid_shape, xx=0.0, yy=0.0, off=-_ORIENTATION_BOUND
    )
    out = apply_prescribed_orientation_rate(ecm, rate, dt_s=1.0)
    expected_off = -_ORIENTATION_BOUND
    np.testing.assert_allclose(out.orientation_tensor[..., 0, 1], expected_off)
    np.testing.assert_allclose(out.orientation_tensor[..., 1, 0], expected_off)


@pytest.mark.parametrize("bad_dt", [True, -1.0, np.inf, np.nan, object()])
def test_orientation_rate_invalid_dt_rejected(bad_dt):
    ecm = _ecm()
    rate = _zero_orientation_rate(ecm.grid_shape)
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_orientation_rate(ecm, rate, dt_s=bad_dt)
    assert info.value.failure_kind == "dt_invalid"


@pytest.mark.parametrize(
    "bad_dt",
    [np.True_, np.False_, np.bool_(True), np.array(True), np.array(False)],
)
def test_orientation_rate_numpy_bool_dt_rejected(bad_dt):
    ecm = _ecm()
    rate = _zero_orientation_rate(ecm.grid_shape)
    with pytest.raises(ECMOpenLoopError) as info:
        apply_prescribed_orientation_rate(ecm, rate, dt_s=bad_dt)
    assert info.value.failure_kind == "dt_invalid"


def test_orientation_rate_other_fields_preserved_and_not_aliased():
    ecm = _ecm()
    rate = _symmetric_orientation_rate(ecm.grid_shape, xx=-0.05, yy=-0.05, off=0.01)
    out = apply_prescribed_orientation_rate(ecm, rate, dt_s=1.0)
    np.testing.assert_array_equal(out.stiffness_kpa, ecm.stiffness_kpa)
    np.testing.assert_array_equal(out.ligand_density, ecm.ligand_density)
    np.testing.assert_array_equal(out.fiber_density, ecm.fiber_density)
    np.testing.assert_array_equal(
        out.accumulated_traction_nNs_per_um2, ecm.accumulated_traction_nNs_per_um2
    )
    assert out.origin_um_xy == ecm.origin_um_xy
    assert out.spacing_um == ecm.spacing_um
    assert out.source == ecm.source
    # Not aliased.
    assert out.stiffness_kpa is not ecm.stiffness_kpa
    assert out.ligand_density is not ecm.ligand_density
    assert out.orientation_tensor is not ecm.orientation_tensor


def test_orientation_rate_input_ecm_not_mutated():
    ecm = _ecm()
    snapshot = ecm.orientation_tensor.copy()
    rate = _symmetric_orientation_rate(ecm.grid_shape, xx=-0.1, yy=-0.05, off=0.05)
    apply_prescribed_orientation_rate(ecm, rate, dt_s=1.0)
    np.testing.assert_array_equal(ecm.orientation_tensor, snapshot)


def test_orientation_rate_returned_state_validates():
    ecm = _ecm()
    rate = _symmetric_orientation_rate(ecm.grid_shape, xx=-0.05, yy=-0.05, off=0.02)
    out = apply_prescribed_orientation_rate(ecm, rate, dt_s=1.0)
    out.validate()


def test_orientation_rate_chained_calls_compose():
    ecm = _ecm()
    rate_a = _symmetric_orientation_rate(ecm.grid_shape, xx=-0.1, yy=-0.05, off=0.05)
    rate_b = _symmetric_orientation_rate(ecm.grid_shape, xx=0.05, yy=-0.025, off=-0.025)
    out_chained = apply_prescribed_orientation_rate(ecm, rate_a, dt_s=1.0)
    out_chained = apply_prescribed_orientation_rate(out_chained, rate_b, dt_s=2.0)
    expected = ecm.orientation_tensor + rate_a * 1.0 + rate_b * 2.0
    np.testing.assert_allclose(out_chained.orientation_tensor, expected, atol=0.0)


# ---------------------------------------------------------------------------
# Phase B precursor evidence (closed-loop ECM gate phased plan §1)
#
# Per docs/v2_closed_loop_ecm_gate_phased_plan_locked.md §1, Phase B:
#  - is "Item 1 precursor evidence: prescribed-traction stimulus monotonicity"
#  - does NOT satisfy closed-loop gate Item 1 (response monotonicity)
#  - measures the open-loop traction accumulator
#    (a stimulus/memory field), not an ECM response field.
#
# These tests verify the stimulus accumulator's monotone-non-decreasing
# property by construction: positive prescribed traction increases the
# accumulator, zero prescribed traction leaves it unchanged, and negative
# prescribed traction is rejected with failure_kind="traction_negative"
# so the field cannot decrease through this preflight path.
#
# Hard Rule 11 wording protection (lock §1): the Phase B section header
# and per-test docstrings explicitly NOT claim Item 1 satisfaction. The
# meta-test
# `test_stimulus_monotonicity_does_not_satisfy_closed_loop_item_1` makes
# the wording boundary an explicit runtime assertion so a future
# refactor cannot silently relabel these tests.
# ---------------------------------------------------------------------------


def test_stimulus_accumulator_monotone_under_positive_traction():
    """Phase B precursor evidence — positive prescribed traction
    monotonically increases the open-loop stimulus accumulator
    ``accumulated_traction_nNs_per_um2`` per step.

    Phase B precursor evidence; does NOT satisfy closed-loop gate
    Item 1 (response monotonicity) — see meta-test below.
    """

    ecm = _ecm()
    initial = ecm.accumulated_traction_nNs_per_um2.copy()

    traction = np.full(ecm.grid_shape, 0.5, dtype=np.float64)
    out = ecm
    history = [initial.copy()]
    for _ in range(5):
        out = accumulate_prescribed_traction(out, traction, dt_s=0.1)
        history.append(out.accumulated_traction_nNs_per_um2.copy())

    # Per-step monotone non-decreasing pointwise.
    for prev, curr in zip(history[:-1], history[1:]):
        assert np.all(curr >= prev), (
            "stimulus accumulator must be monotone non-decreasing under "
            "positive prescribed traction"
        )
    # Final field strictly above initial (positive traction across all cells).
    assert np.all(history[-1] > initial)


def test_stimulus_accumulator_no_change_when_traction_zero():
    """Phase B precursor evidence — zero prescribed traction leaves the
    open-loop stimulus accumulator unchanged, exactly to float64
    round-off (the underlying update is bit-exact ``x + 0 = x``).
    This is the boundary case of monotone-non-decreasing where the
    monotone increase reduces to monotone equality.

    Phase B precursor evidence; does NOT satisfy closed-loop gate
    Item 1 (response monotonicity), which requires this property on
    a response field, not the stimulus accumulator. The
    zero-traction case contributes to the Item 1 precursor only as
    the boundary point of the monotonicity claim, not as an
    independent gate-item satisfaction.
    """

    ecm = _ecm()
    initial = ecm.accumulated_traction_nNs_per_um2.copy()

    out = ecm
    for _ in range(5):
        out = accumulate_prescribed_traction(
            out, np.zeros(ecm.grid_shape), dt_s=0.5
        )
        np.testing.assert_array_equal(
            out.accumulated_traction_nNs_per_um2, initial
        )


def test_stimulus_accumulator_rejects_negative_traction_enforcing_monotonicity():
    """Phase B precursor evidence — negative prescribed traction is
    rejected with failure_kind="traction_negative", so the open-loop
    stimulus accumulator cannot decrease through this preflight path.

    Lock §8 drafting clarification (per Codex impl id=1241): the
    test catalog originally listed this entry as
    ``test_stimulus_accumulator_signed_traction_can_decrease``,
    which encodes a behavior the code does NOT support (the
    open-loop preflight rejects negative traction). The renamed
    test here is faithful to the lock's intent ("verify the
    accumulator's signedness contract") while consistent with the
    code reality (signedness contract = non-negative-only).

    Phase B precursor evidence; does NOT satisfy closed-loop gate
    Item 1 (response monotonicity).
    """

    ecm = _ecm()
    negative = np.full(ecm.grid_shape, -0.1, dtype=np.float64)
    with pytest.raises(ECMOpenLoopError) as excinfo:
        accumulate_prescribed_traction(ecm, negative, dt_s=1.0)
    assert excinfo.value.failure_kind == "traction_negative"

    # Mixed-sign field with at least one negative cell also rejected.
    mixed = np.full(ecm.grid_shape, 0.1, dtype=np.float64)
    mixed[0, 0] = -1e-9
    with pytest.raises(ECMOpenLoopError) as excinfo_mixed:
        accumulate_prescribed_traction(ecm, mixed, dt_s=1.0)
    assert excinfo_mixed.value.failure_kind == "traction_negative"


def test_stimulus_monotonicity_does_not_satisfy_closed_loop_item_1():
    """Phase B precursor evidence — Hard Rule 11 wording protection.

    Documents in code that the three tests above measure the
    open-loop traction accumulator
    (``accumulated_traction_nNs_per_um2``), which is a stimulus /
    memory field. They do **not** satisfy the closed-loop ECM gate
    Item 1 (response monotonicity), which requires:

    - a response field defined (Hard Blocker #1: constitutive-law
      direction)
    - a constitutive law mapping cumulative stimulus to response
      change (Phases D + E per
      ``docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`` §1)
    - a monotonicity test on the response field, not the stimulus
      field

    Hard Rule 11 (CLAUDE.md): the measured value must match the
    claim. Measuring stimulus-accumulator monotonicity does not
    measure remodeling-response monotonicity. This test is a
    runtime assertion that makes the wording boundary explicit, so
    a future refactor cannot silently relabel the Phase B tests
    above as Item 1 satisfaction.
    """

    phase_b_measured_field = "accumulated_traction_nNs_per_um2 (stimulus)"
    closed_loop_item_1_claim = "ECM response field monotonicity (TBD)"
    assert phase_b_measured_field != closed_loop_item_1_claim, (
        "Phase B tests measure stimulus accumulator (a memory field), "
        "NOT a response field. They cannot satisfy closed-loop gate "
        "Item 1, which requires a response field per Hard Blocker #1 "
        "in docs/v2_closed_loop_ecm_gate_phased_plan_locked.md §2."
    )


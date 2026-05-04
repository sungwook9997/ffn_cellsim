from __future__ import annotations

import math

import numpy as np
import pytest

from acs.v2.dynamics.focal_adhesion import (
    FocalAdhesionDynamicsError,
    FocalAdhesionDynamicsParameters,
    compute_radial_tangential_decomposition,
    step_focal_adhesions_static,
)
from acs.v2.focal_adhesion import FocalAdhesionState


def _fa(
    adhesion_id: str,
    *,
    position: tuple[float, float] = (1.0, 0.0),
    state: str = "mature",
    maturity: float = 0.5,
    bound_fraction: float = 0.5,
) -> FocalAdhesionState:
    if state == "unbound":
        bound_fraction = 0.0  # schema invariant
    return FocalAdhesionState(
        adhesion_id=adhesion_id,
        cell_id="cell-A",
        position_um_xy=position,
        age_s=10.0,
        maturity=maturity,
        bound_fraction=bound_fraction,
        state=state,  # type: ignore[arg-type]
    )


def _params(**overrides) -> FocalAdhesionDynamicsParameters:
    base = dict(
        dt_fa_s=1e-2,
        traction_scale_nN=1.0,
    )
    base.update(overrides)
    return FocalAdhesionDynamicsParameters(**base).validate()


def test_minimal_params_validate():
    params = _params()
    assert params.dt_fa_s == 1e-2
    assert params.traction_scale_nN == 1.0


def test_invalid_rate_rejected():
    for bad in (-0.1, float("nan"), float("inf"), True):
        with pytest.raises(FocalAdhesionDynamicsError) as info:
            FocalAdhesionDynamicsParameters(
                dt_fa_s=1e-2,
                traction_scale_nN=1.0,
                k_maturity_per_s=bad,  # type: ignore[arg-type]
            ).validate()
        assert info.value.failure_kind == "rate_invalid"


def test_invalid_max_traction_rejected():
    for bad in (-0.1, float("nan"), float("inf"), True):
        with pytest.raises(FocalAdhesionDynamicsError) as info:
            FocalAdhesionDynamicsParameters(
                dt_fa_s=1e-2,
                traction_scale_nN=1.0,
                max_traction_nN=bad,  # type: ignore[arg-type]
            ).validate()
        assert info.value.failure_kind == "max_traction_invalid"


def test_max_traction_zero_allowed():
    """Gate §2: max_traction_nN == 0 is allowed; only zero-magnitude traction passes."""

    params = FocalAdhesionDynamicsParameters(
        dt_fa_s=1e-2, traction_scale_nN=0.0, max_traction_nN=0.0
    ).validate()
    fa = _fa("fa-1", state="mature", maturity=1.0, bound_fraction=1.0)
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    np.testing.assert_array_equal(result.cell_force_nN_xy, np.zeros((1, 2)))


def test_dt_rate_violation():
    with pytest.raises(FocalAdhesionDynamicsError) as info:
        FocalAdhesionDynamicsParameters(
            dt_fa_s=1.0, traction_scale_nN=1.0, k_maturity_per_s=1.0
        ).validate()
    assert info.value.failure_kind == "dt_rate_violation"


def test_explicit_zero_rate_no_state_change():
    params = _params(k_maturity_per_s=0.0, k_bind_per_s=0.0, k_unbind_per_s=0.0)
    fa = _fa("fa-1", state="mature", maturity=0.4, bound_fraction=0.6)
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    updated = result.updated_adhesions[0]
    assert updated.maturity == 0.4
    assert updated.bound_fraction == 0.6


def test_unbound_traction_is_zero():
    params = _params()
    fa = _fa("fa-1", state="unbound", maturity=1.0, bound_fraction=0.0)
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    np.testing.assert_array_equal(result.cell_force_nN_xy[0], np.zeros(2))
    np.testing.assert_array_equal(result.substrate_reaction_nN_xy[0], np.zeros(2))


def test_released_traction_is_zero_even_with_bound_fraction():
    """Gate §6 schema docstring: released permits bound_fraction>0 history,
    but dynamics step forces traction to (0, 0)."""

    params = _params()
    fa = _fa("fa-1", state="released", maturity=0.8, bound_fraction=0.7)
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    np.testing.assert_array_equal(result.cell_force_nN_xy[0], np.zeros(2))


def test_mature_traction_inward_radial_default():
    """Gate §5 sign: FA at (1, 0) with centroid (0, 0), scale=1, mat=1, bound=1
    → traction = (-1, 0) (inward radial)."""

    params = _params()
    fa = _fa("fa-1", position=(1.0, 0.0), state="mature", maturity=1.0, bound_fraction=1.0)
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    np.testing.assert_allclose(result.cell_force_nN_xy[0], np.array([-1.0, 0.0]), atol=1e-12)
    np.testing.assert_allclose(
        result.substrate_reaction_nN_xy[0], np.array([1.0, 0.0]), atol=1e-12
    )


def test_traction_magnitude_scales_with_maturity_and_bound():
    params = _params(traction_scale_nN=2.0)
    fa = _fa(
        "fa-1", position=(2.0, 0.0), state="mature", maturity=0.5, bound_fraction=0.4
    )
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    expected_magnitude = 2.0 * 0.5 * 0.4
    actual_magnitude = float(np.linalg.norm(result.cell_force_nN_xy[0]))
    assert actual_magnitude == pytest.approx(expected_magnitude)


def test_supplied_axis_overrides_default_radial():
    params = _params()
    fa = _fa("fa-1", position=(1.0, 0.0), state="mature", maturity=1.0, bound_fraction=1.0)
    # Override: force axis along +y direction.
    result = step_focal_adhesions_static(
        [fa],
        centroid_um_xy=(0.0, 0.0),
        params=params,
        traction_axis_xy_per_fa=[(0.0, 3.0)],  # non-unit; will be normalised
    )
    np.testing.assert_allclose(result.cell_force_nN_xy[0], np.array([0.0, 1.0]), atol=1e-12)


def test_per_fa_action_reaction_exact():
    params = _params()
    fa1 = _fa("fa-1", position=(1.0, 0.0), state="mature", maturity=1.0, bound_fraction=1.0)
    fa2 = _fa("fa-2", position=(0.0, 2.0), state="mature", maturity=0.5, bound_fraction=0.5)
    result = step_focal_adhesions_static([fa1, fa2], centroid_um_xy=(0.0, 0.0), params=params)
    np.testing.assert_array_equal(
        result.cell_force_nN_xy + result.substrate_reaction_nN_xy, np.zeros((2, 2))
    )


def test_aggregate_newton_3_zero():
    params = _params()
    fas = [
        _fa(f"fa-{i}", position=(math.cos(theta), math.sin(theta)),
            state="mature", maturity=1.0, bound_fraction=1.0)
        for i, theta in enumerate(np.linspace(0.0, 2 * math.pi, 5)[:-1])
    ]
    result = step_focal_adhesions_static(fas, centroid_um_xy=(0.0, 0.0), params=params)
    aggregate = (
        result.cell_force_nN_xy.sum(axis=0)
        + result.substrate_reaction_nN_xy.sum(axis=0)
    )
    np.testing.assert_allclose(aggregate, np.zeros(2), atol=1e-12)


def test_radial_tangential_round_trip():
    fa = _fa("fa-1", position=(1.0, 1.0), state="mature", maturity=0.7, bound_fraction=0.8)
    params = _params(traction_scale_nN=2.0)
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    n_radial = (np.array([1.0, 1.0]) - np.array([0.0, 0.0])) / math.sqrt(2.0)
    n_tangential = np.array([-n_radial[1], n_radial[0]])
    reconstructed = (
        result.radial_components[0] * n_radial
        + result.tangential_components[0] * n_tangential
    )
    np.testing.assert_allclose(reconstructed, result.cell_force_nN_xy[0], atol=1e-12)


def test_default_radial_traction_has_negative_radial_component():
    """Inward traction → projection onto outward radial is negative."""

    fa = _fa("fa-1", position=(2.0, 0.0), state="mature", maturity=1.0, bound_fraction=1.0)
    params = _params()
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    assert result.radial_components[0] < 0.0
    assert result.tangential_components[0] == pytest.approx(0.0, abs=1e-12)


def test_max_traction_exceeded_raises_no_clamp():
    params = _params(traction_scale_nN=10.0, max_traction_nN=1.0)
    fa = _fa("fa-1", state="mature", maturity=1.0, bound_fraction=1.0)
    with pytest.raises(FocalAdhesionDynamicsError) as info:
        step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    assert info.value.failure_kind == "max_traction_exceeded"


def test_position_equals_centroid_rejected():
    params = _params()
    fa = _fa("fa-1", position=(0.0, 0.0), state="mature", maturity=1.0, bound_fraction=1.0)
    with pytest.raises(FocalAdhesionDynamicsError) as info:
        step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    assert info.value.failure_kind == "radial_axis_undefined"


def test_position_equals_centroid_rejected_even_with_supplied_axis():
    """Gate §2: radial decomposition requires defined radial unit
    regardless of whether an explicit axis is supplied."""

    params = _params()
    fa = _fa("fa-1", position=(0.0, 0.0), state="mature", maturity=1.0, bound_fraction=1.0)
    with pytest.raises(FocalAdhesionDynamicsError) as info:
        step_focal_adhesions_static(
            [fa],
            centroid_um_xy=(0.0, 0.0),
            params=params,
            traction_axis_xy_per_fa=[(1.0, 0.0)],
        )
    assert info.value.failure_kind == "radial_axis_undefined"


def test_centroid_shape_rejected():
    params = _params()
    fa = _fa("fa-1", state="mature")
    with pytest.raises(FocalAdhesionDynamicsError) as info:
        step_focal_adhesions_static([fa], centroid_um_xy=(0.0,), params=params)
    assert info.value.failure_kind == "centroid_shape"


def test_axis_zero_magnitude_rejected():
    params = _params()
    fa = _fa("fa-1", state="mature")
    with pytest.raises(FocalAdhesionDynamicsError) as info:
        step_focal_adhesions_static(
            [fa],
            centroid_um_xy=(0.0, 0.0),
            params=params,
            traction_axis_xy_per_fa=[(0.0, 0.0)],
        )
    assert info.value.failure_kind == "axis_zero_magnitude"


def test_axis_per_fa_length_mismatch_rejected():
    params = _params()
    fa1 = _fa("fa-1", state="mature")
    fa2 = _fa("fa-2", position=(0.0, 1.0), state="mature")
    with pytest.raises(FocalAdhesionDynamicsError) as info:
        step_focal_adhesions_static(
            [fa1, fa2],
            centroid_um_xy=(0.0, 0.0),
            params=params,
            traction_axis_xy_per_fa=[(1.0, 0.0)],  # only 1 entry for 2 FAs
        )
    assert info.value.failure_kind == "axis_shape"


def test_slipping_bound_fraction_decreases_only():
    params = _params(k_unbind_per_s=10.0)  # dt=1e-2, dt·k=0.1 ≤ 0.5
    fa = _fa(
        "fa-1", position=(1.0, 0.0), state="slipping", maturity=0.6, bound_fraction=0.4
    )
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    updated = result.updated_adhesions[0]
    assert updated.bound_fraction < fa.bound_fraction


def test_slipping_ignores_k_bind():
    """State table: slipping never increases bound_fraction even if k_bind supplied."""

    params = _params(k_bind_per_s=10.0)
    fa = _fa("fa-1", state="slipping", maturity=0.6, bound_fraction=0.4)
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    updated = result.updated_adhesions[0]
    assert updated.bound_fraction == 0.4


def test_mature_maturity_increases_under_k_maturity():
    params = _params(k_maturity_per_s=10.0)
    fa = _fa("fa-1", state="mature", maturity=0.4, bound_fraction=0.5)
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    updated = result.updated_adhesions[0]
    assert updated.maturity > 0.4


def test_unbound_state_does_not_apply_rates():
    params = _params(k_maturity_per_s=10.0, k_bind_per_s=10.0, k_unbind_per_s=10.0)
    fa = _fa("fa-1", state="unbound", maturity=0.0, bound_fraction=0.0)
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    updated = result.updated_adhesions[0]
    assert updated.maturity == 0.0
    assert updated.bound_fraction == 0.0


def test_age_increments_by_dt():
    params = _params()
    fa = _fa("fa-1", state="mature")
    result = step_focal_adhesions_static([fa], centroid_um_xy=(0.0, 0.0), params=params)
    assert result.updated_adhesions[0].age_s == pytest.approx(fa.age_s + params.dt_fa_s)


def test_diagnostics_aggregate_cancellation():
    params = _params()
    fas = [_fa(f"fa-{i}", position=(math.cos(t), math.sin(t)), state="mature",
               maturity=1.0, bound_fraction=1.0)
           for i, t in enumerate([0.0, math.pi / 2, math.pi, 3 * math.pi / 2])]
    result = step_focal_adhesions_static(fas, centroid_um_xy=(0.0, 0.0), params=params)
    cell_sum = np.array(result.diagnostics["aggregate_cell_force_nN_xy"])
    sub_sum = np.array(result.diagnostics["aggregate_substrate_reaction_nN_xy"])
    np.testing.assert_allclose(cell_sum + sub_sum, np.zeros(2), atol=1e-12)
    # Symmetric arrangement: aggregate cell force should be near zero too.
    np.testing.assert_allclose(cell_sum, np.zeros(2), atol=1e-12)


def test_decomposition_helper_round_trip():
    radial, tangential = compute_radial_tangential_decomposition(
        traction_xy=(0.5, 0.5),
        position_xy=(2.0, 0.0),
        centroid_xy=(0.0, 0.0),
    )
    n_radial = np.array([1.0, 0.0])
    n_tangential = np.array([0.0, 1.0])
    reconstructed = radial * n_radial + tangential * n_tangential
    np.testing.assert_allclose(reconstructed, np.array([0.5, 0.5]), atol=1e-12)


def test_decomposition_helper_position_equals_centroid_rejects():
    with pytest.raises(FocalAdhesionDynamicsError) as info:
        compute_radial_tangential_decomposition(
            traction_xy=(1.0, 0.0),
            position_xy=(0.0, 0.0),
            centroid_xy=(0.0, 0.0),
        )
    assert info.value.failure_kind == "radial_axis_undefined"

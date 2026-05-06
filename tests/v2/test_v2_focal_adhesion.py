from __future__ import annotations

import math

import pytest

from acs.v2.focal_adhesion import FocalAdhesionState


def _fa(**overrides) -> FocalAdhesionState:
    base = dict(
        adhesion_id="fa-1",
        cell_id="cell-1",
        position_um_xy=(1.0, 0.5),
        age_s=60.0,
        maturity=0.4,
        bound_fraction=0.8,
    )
    base.update(overrides)
    return FocalAdhesionState(**base)


def test_minimal_fa_validates():
    fa = _fa()
    fa.validate()
    assert fa.state == "nascent"
    assert fa.traction_force_nN_xy == (0.0, 0.0)
    assert fa.source == "simulated"


def test_full_optional_fa_validates():
    fa = _fa(
        state="mature",
        traction_force_nN_xy=(3.0, -2.0),
        linked_protrusion_id="evt-3",
        source="marker_detected",
    )
    fa.validate()


def test_position_xy_wrong_length_rejected():
    fa = _fa(position_um_xy=(1.0, 0.5, 0.0))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="position_um_xy"):
        fa.validate()


def test_position_xy_nan_rejected():
    fa = _fa(position_um_xy=(math.nan, 0.5))
    with pytest.raises(ValueError, match="position_um_xy"):
        fa.validate()


def test_negative_age_rejected():
    fa = _fa(age_s=-1.0)
    with pytest.raises(ValueError, match="age_s"):
        fa.validate()


def test_maturity_out_of_range_rejected():
    fa = _fa(maturity=1.5)
    with pytest.raises(ValueError, match="maturity"):
        fa.validate()


def test_bound_fraction_out_of_range_rejected():
    fa = _fa(bound_fraction=-0.01)
    with pytest.raises(ValueError, match="bound_fraction"):
        fa.validate()


def test_invalid_state_rejected():
    fa = _fa(state="exploded")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="state"):
        fa.validate()


def test_unbound_with_traction_rejected():
    fa = _fa(state="unbound", bound_fraction=0.0, traction_force_nN_xy=(1.0, 0.0))
    with pytest.raises(ValueError, match="traction"):
        fa.validate()


def test_unbound_with_bound_fraction_rejected():
    fa = _fa(state="unbound", bound_fraction=0.5, traction_force_nN_xy=(0.0, 0.0))
    with pytest.raises(ValueError, match="bound_fraction"):
        fa.validate()


def test_unbound_with_zero_traction_and_zero_bound_accepted():
    fa = _fa(state="unbound", bound_fraction=0.0, traction_force_nN_xy=(0.0, 0.0))
    fa.validate()


def test_traction_xy_nan_rejected():
    fa = _fa(traction_force_nN_xy=(math.nan, 0.0))
    with pytest.raises(ValueError, match="traction_force_nN_xy"):
        fa.validate()


def test_linked_protrusion_id_empty_rejected():
    fa = _fa(linked_protrusion_id="   ")
    with pytest.raises(ValueError, match="linked_protrusion_id"):
        fa.validate()


def test_invalid_source_rejected():
    fa = _fa(source="forged")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="source"):
        fa.validate()

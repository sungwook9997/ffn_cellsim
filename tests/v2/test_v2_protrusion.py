from __future__ import annotations

import math

import pytest

from acs.v2.protrusion import ProtrusionEvent


def _event(**overrides) -> ProtrusionEvent:
    base = dict(
        cell_id="cell-1",
        start_time_s=10.0,
        boundary_angle_rad=0.0,
        length_um=1.5,
    )
    base.update(overrides)
    return ProtrusionEvent(**base)


def test_minimal_event_validates():
    event = _event()
    event.validate()
    assert event.event_type == "lamellipodium"
    assert event.state == "growing"
    assert event.source == "simulated"


def test_filopodium_with_full_optionals_validates():
    event = _event(
        event_id="evt-1",
        end_time_s=15.0,
        event_type="filopodium",
        state="ended",
        root_position_um_xy=(0.0, 0.0),
        tip_position_um_xy=(2.0, 1.0),
        direction_um_xy=(0.894, 0.447),
        width_um=0.2,
        width_rad=0.1,
        force_candidate_nN=5.0,
        associated_adhesion_ids=("fa-1", "fa-2"),
        source="detected",
        lifetime_s=5.0,
        confidence=0.9,
    )
    event.validate()


def test_invalid_event_type_rejected():
    event = _event(event_type="retraction")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="event_type"):
        event.validate()


def test_invalid_state_rejected():
    event = _event(state="exploding")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="state"):
        event.validate()


def test_negative_start_time_rejected():
    event = _event(start_time_s=-1.0)
    with pytest.raises(ValueError, match="start_time_s"):
        event.validate()


def test_nan_start_time_rejected():
    event = _event(start_time_s=math.nan)
    with pytest.raises(ValueError, match="start_time_s"):
        event.validate()


def test_end_before_start_rejected():
    event = _event(start_time_s=10.0, end_time_s=5.0)
    with pytest.raises(ValueError, match="end_time_s"):
        event.validate()


def test_lifetime_inconsistent_with_end_time_rejected():
    event = _event(start_time_s=10.0, end_time_s=15.0, lifetime_s=10.0)
    with pytest.raises(ValueError, match="lifetime_s"):
        event.validate()


def test_lifetime_consistent_with_end_time_accepted():
    event = _event(start_time_s=10.0, end_time_s=15.0, lifetime_s=5.0)
    event.validate()


def test_lifetime_alone_accepted():
    event = _event(lifetime_s=4.0)
    event.validate()


def test_position_xy_wrong_length_rejected():
    event = _event(root_position_um_xy=(0.0, 0.0, 0.0))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="root_position_um_xy"):
        event.validate()


def test_negative_width_rejected():
    event = _event(width_um=-0.1)
    with pytest.raises(ValueError, match="width_um"):
        event.validate()


def test_confidence_out_of_range_rejected():
    for bad in (-0.01, 1.01, math.nan):
        event = _event(confidence=bad)
        with pytest.raises(ValueError, match="confidence"):
            event.validate()


def test_associated_adhesion_duplicate_rejected():
    event = _event(associated_adhesion_ids=("fa-1", "fa-1"))
    with pytest.raises(ValueError, match="associated_adhesion_ids"):
        event.validate()


def test_associated_adhesion_empty_string_rejected():
    event = _event(associated_adhesion_ids=("fa-1", " "))
    with pytest.raises(ValueError, match="associated_adhesion_ids"):
        event.validate()


def test_invalid_source_rejected():
    event = _event(source="forged")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="source"):
        event.validate()


def test_empty_event_id_rejected_when_provided():
    event = _event(event_id="   ")
    with pytest.raises(ValueError, match="event_id"):
        event.validate()


def test_negative_force_candidate_rejected():
    event = _event(force_candidate_nN=-3.0)
    with pytest.raises(ValueError, match="force_candidate_nN"):
        event.validate()


def test_zero_direction_vector_rejected():
    event = _event(direction_um_xy=(0.0, 0.0))
    with pytest.raises(ValueError, match="direction_um_xy"):
        event.validate()


def test_lifetime_consistency_within_tolerance_accepted():
    event = _event(start_time_s=10.0, end_time_s=15.0, lifetime_s=5.0 + 5e-10)
    event.validate()

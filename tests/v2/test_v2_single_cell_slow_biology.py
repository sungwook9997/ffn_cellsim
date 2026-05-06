from __future__ import annotations

import numpy as np
import pytest

from acs.v2.measurement_boundary import MeasurementBoundary
from acs.v2.single_cell import SingleCellState


def _state(**overrides) -> SingleCellState:
    boundary = MeasurementBoundary.from_array(
        np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [1.0, 1.0],
                [0.0, 1.0],
            ]
        ),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    base = dict(
        cell_id="cell-1",
        time_s=0.0,
        measurement_boundary=boundary,
    )
    base.update(overrides)
    return SingleCellState(**base)


def test_default_slow_biology_hooks_validate_and_are_off():
    state = _state()
    state.validate()
    assert state.cell_state == "alive"
    assert state.cell_age_s == 0.0
    assert state.cell_cycle_phase is None
    assert state.division_count == 0
    assert state.parent_cell_id is None
    assert state.mechanosignal_yap_taz is None
    assert state.neighbor_cell_ids == ()


def test_dead_cell_state_validates_and_keeps_geometry():
    state = _state(cell_state="dead")
    state.validate()
    assert state.cell_state == "dead"
    assert state.projected_area_um2() == pytest.approx(1.0)


def test_invalid_cell_state_rejected():
    state = _state(cell_state="frozen")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="cell_state"):
        state.validate()


def test_negative_age_rejected():
    state = _state(cell_age_s=-1.0)
    with pytest.raises(ValueError, match="cell_age_s"):
        state.validate()


def test_invalid_cell_cycle_phase_rejected():
    state = _state(cell_cycle_phase="X")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="cell_cycle_phase"):
        state.validate()


def test_valid_cell_cycle_phases_accepted():
    for phase in ("G0", "G1", "S", "G2", "M"):
        state = _state(cell_cycle_phase=phase)
        state.validate()


def test_negative_division_count_rejected():
    state = _state(division_count=-1)
    with pytest.raises(ValueError, match="division_count"):
        state.validate()


def test_parent_id_must_differ_from_cell_id():
    state = _state(parent_cell_id="cell-1")
    with pytest.raises(ValueError, match="parent_cell_id"):
        state.validate()


def test_empty_parent_id_string_rejected():
    state = _state(parent_cell_id="   ")
    with pytest.raises(ValueError, match="parent_cell_id"):
        state.validate()


def test_yap_taz_out_of_range_rejected():
    for bad in (-0.01, 1.01, float("nan")):
        state = _state(mechanosignal_yap_taz=bad)
        with pytest.raises(ValueError, match="mechanosignal_yap_taz"):
            state.validate()


def test_yap_taz_in_range_accepted():
    for ok in (0.0, 0.5, 1.0):
        state = _state(mechanosignal_yap_taz=ok)
        state.validate()


def test_neighbor_self_reference_rejected():
    state = _state(neighbor_cell_ids=("cell-2", "cell-1"))
    with pytest.raises(ValueError, match="own id"):
        state.validate()


def test_neighbor_duplicate_rejected():
    state = _state(neighbor_cell_ids=("cell-2", "cell-2"))
    with pytest.raises(ValueError, match="unique"):
        state.validate()


def test_neighbor_empty_string_rejected():
    state = _state(neighbor_cell_ids=("cell-2", " "))
    with pytest.raises(ValueError, match="neighbor_cell_ids"):
        state.validate()


def test_neighbor_unique_set_accepted():
    state = _state(neighbor_cell_ids=("cell-2", "cell-3", "cell-4"))
    state.validate()
    assert state.neighbor_cell_ids == ("cell-2", "cell-3", "cell-4")


def test_cell_age_rejects_nan_and_inf():
    for bad in (float("nan"), float("inf"), float("-inf")):
        state = _state(cell_age_s=bad)
        with pytest.raises(ValueError, match="cell_age_s"):
            state.validate()


def test_division_count_rejects_non_integer():
    state_float = _state(division_count=1.5)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="division_count"):
        state_float.validate()
    state_nan = _state(division_count=float("nan"))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="division_count"):
        state_nan.validate()
    state_bool = _state(division_count=True)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="division_count"):
        state_bool.validate()

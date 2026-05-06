from __future__ import annotations

import numpy as np
import pytest

from acs.v2.ecm_substrate import ECMSubstrateState


def _grid(nx: int = 4, ny: int = 3) -> dict:
    stiffness = 1.0 * np.ones((nx, ny), dtype=np.float64)
    ligand = 0.5 * np.ones((nx, ny), dtype=np.float64)
    fiber = 0.2 * np.ones((nx, ny), dtype=np.float64)
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    accumulated = np.zeros((nx, ny), dtype=np.float64)
    return dict(
        origin_um_xy=(0.0, 0.0),
        spacing_um=2.0,
        stiffness_kpa=stiffness,
        ligand_density=ligand,
        fiber_density=fiber,
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=accumulated,
    )


def test_minimal_ecm_validates():
    ecm = ECMSubstrateState(**_grid())
    ecm.validate()
    assert ecm.grid_shape == (4, 3)
    assert ecm.source == "synthetic"


def test_negative_spacing_rejected():
    grid = _grid()
    grid["spacing_um"] = -1.0
    with pytest.raises(ValueError, match="spacing_um"):
        ECMSubstrateState(**grid).validate()


def test_origin_wrong_length_rejected():
    grid = _grid()
    grid["origin_um_xy"] = (0.0, 0.0, 0.0)
    with pytest.raises(ValueError, match="origin_um_xy"):
        ECMSubstrateState(**grid).validate()


def test_stiffness_negative_rejected():
    grid = _grid()
    grid["stiffness_kpa"][0, 0] = -1.0
    with pytest.raises(ValueError, match="stiffness_kpa"):
        ECMSubstrateState(**grid).validate()


def test_ligand_out_of_range_rejected():
    grid = _grid()
    grid["ligand_density"][0, 1] = 1.2
    with pytest.raises(ValueError, match="ligand_density"):
        ECMSubstrateState(**grid).validate()


def test_fiber_negative_rejected():
    grid = _grid()
    grid["fiber_density"][1, 1] = -0.01
    with pytest.raises(ValueError, match="fiber_density"):
        ECMSubstrateState(**grid).validate()


def test_accumulated_traction_negative_rejected():
    grid = _grid()
    grid["accumulated_traction_nNs_per_um2"][0, 0] = -0.5
    with pytest.raises(ValueError, match="accumulated_traction"):
        ECMSubstrateState(**grid).validate()


def test_orientation_wrong_shape_rejected():
    grid = _grid()
    grid["orientation_tensor"] = np.zeros((4, 3, 2), dtype=np.float64)
    with pytest.raises(ValueError, match="orientation_tensor"):
        ECMSubstrateState(**grid).validate()


def test_orientation_asymmetry_rejected():
    grid = _grid()
    grid["orientation_tensor"][0, 0, 0, 1] = 0.7
    grid["orientation_tensor"][0, 0, 1, 0] = 0.0
    with pytest.raises(ValueError, match="symmetric"):
        ECMSubstrateState(**grid).validate()


def test_orientation_nan_rejected():
    grid = _grid()
    grid["orientation_tensor"][0, 0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="orientation_tensor"):
        ECMSubstrateState(**grid).validate()


def test_grid_shape_mismatch_rejected():
    grid = _grid()
    grid["fiber_density"] = np.zeros((4, 4), dtype=np.float64)
    with pytest.raises(ValueError, match="share shape"):
        ECMSubstrateState(**grid).validate()


def test_invalid_source_rejected():
    grid = _grid()
    state = ECMSubstrateState(source="ai_imagined", **grid)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="source"):
        state.validate()


def test_orientation_bound_violation_rejected():
    grid = _grid()
    grid["orientation_tensor"][0, 0, 0, 0] = 1.5
    with pytest.raises(ValueError, match="bounded"):
        ECMSubstrateState(**grid).validate()


def test_orientation_bound_at_unity_accepted():
    grid = _grid()
    grid["orientation_tensor"][0, 0, 0, 0] = 1.0
    grid["orientation_tensor"][0, 0, 1, 1] = -1.0
    ECMSubstrateState(**grid).validate()


def test_zero_grid_rejected():
    grid = _grid()
    grid["stiffness_kpa"] = np.zeros((0, 3), dtype=np.float64)
    grid["ligand_density"] = np.zeros((0, 3), dtype=np.float64)
    grid["fiber_density"] = np.zeros((0, 3), dtype=np.float64)
    grid["orientation_tensor"] = np.zeros((0, 3, 2, 2), dtype=np.float64)
    grid["accumulated_traction_nNs_per_um2"] = np.zeros((0, 3), dtype=np.float64)
    with pytest.raises(ValueError, match="nx>0"):
        ECMSubstrateState(**grid).validate()

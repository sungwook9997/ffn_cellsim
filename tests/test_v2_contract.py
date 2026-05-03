from __future__ import annotations

import numpy as np
import pytest

from acs.v2.data_contract import ImagingDatasetSpec, default_single_cell_contract
from acs.v2.measurement_boundary import MeasurementBoundary
from acs.v2.single_cell import FocalAdhesionState, ProtrusionEvent, SingleCellState


def test_default_single_cell_contract_validates():
    dataset = ImagingDatasetSpec(
        dataset_id="single_cell_demo",
        root_uri="data/imaging/single_cell_demo",
        modality="confocal_live",
        channels=("membrane", "actin"),
        voxel_size_um_xyz=(0.25, 0.25, 1.0),
        frame_interval_s=60.0,
        n_timepoints=12,
        calibration_timepoints=(0, 1, 2, 3),
        validation_timepoints=(8, 9, 10, 11),
    )

    contract = default_single_cell_contract(dataset)

    contract.validate()
    assert len(contract.calibration_metrics) == 1
    assert len(contract.validation_metrics) == 2


def test_contract_rejects_calibration_validation_overlap():
    dataset = ImagingDatasetSpec(
        dataset_id="overlap",
        root_uri="data/imaging/overlap",
        modality="confocal_live",
        channels=("membrane",),
        voxel_size_um_xyz=(0.25, 0.25, 1.0),
        frame_interval_s=60.0,
        n_timepoints=3,
        calibration_timepoints=(1,),
        validation_timepoints=(1,),
    )

    contract = default_single_cell_contract(dataset)

    with pytest.raises(ValueError, match="overlap"):
        contract.validate()


def test_single_cell_state_projected_area_and_nested_validation():
    boundary = np.array([
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 1.0],
        [0.0, 1.0],
    ])
    state = SingleCellState(
        cell_id="cell-1",
        time_s=120.0,
        measurement_boundary=MeasurementBoundary.from_array(
            boundary,
            coordinate_convention="world_um_y_up",
            source_modality="manual_test_outline",
        ),
        height_um=6.0,
        polarity_xy=(1.0, 0.0),
        protrusions=[
            ProtrusionEvent(
                time_s=120.0,
                cell_id="cell-1",
                event_type="lamellipodium",
                boundary_angle_rad=0.0,
                length_um=3.0,
                lifetime_s=180.0,
            )
        ],
        adhesions=[
            FocalAdhesionState(
                adhesion_id="fa-1",
                cell_id="cell-1",
                position_um_xyz=(1.0, 0.5, 0.0),
                age_s=60.0,
                maturity=0.4,
                bound_fraction=0.8,
            )
        ],
    )

    assert state.projected_area_um2() == pytest.approx(2.0)

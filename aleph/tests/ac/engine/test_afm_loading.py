"""Host gates for selecting and geometrically initializing an AFM loading path."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from aleph.engine.afm_loading import derive_force_free_initial_centre, load_ready_afm_path


def _manifest(
    path: Path,
    *,
    status: str = "READY_FOR_CUDA_DRIVER",
    result_class: str = "mechanism_demo_not_quantitative",
) -> Path:
    points = [
        {"speed_um_s": speed, "depth_um": depth, "seed": seed}
        for speed in (0.1, 0.2) for depth in (0.0, 0.25) for seed in (0, 1, 2)
    ]
    payload = {
        "schema": "afm-sweep-preflight@1",
        "status": status,
        "result_class_requested": result_class,
        "protocol": {
            "geometry": "suspended_round",
            "holder": "exterior_stokes_six_rigid_modes",
            "indenter_radius_um": 5.0,
            "exterior_viscosity_pa_s": 1.0e-3,
            "apparatus_source": "test source",
        },
        "pi0": {"value_pa": 50.0, "source": "test Pi0 source"},
        "contact": {
            "contact_distance_um": 0.05,
            "stiffness_pn_per_um": 100.0,
            "provenance": "test contact source",
        },
        "axes": {"depth_um": [0.0, 0.25]},
        "blockers": [] if status == "READY_FOR_CUDA_DRIVER" else ["PI0_MISSING"],
        "points": points if status == "READY_FOR_CUDA_DRIVER" else [],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_ready_manifest_selects_one_speed_seed_monotone_path(tmp_path: Path) -> None:
    path = load_ready_afm_path(_manifest(tmp_path / "ready.json"), seed=1, speed_um_s=0.2)
    assert path.depths_um == (0.0, 0.25)
    assert path.seed == 1 and path.speed_um_s == 0.2
    assert path.pi0_pa == 50.0


def test_blocked_manifest_never_reaches_cuda(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="READY"):
        load_ready_afm_path(_manifest(tmp_path / "blocked.json", status="BLOCKED"), seed=0, speed_um_s=0.1)


def test_forged_quantitative_ready_manifest_remains_blocked(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="active-cortex"):
        load_ready_afm_path(
            _manifest(tmp_path / "quantitative.json", result_class="quantitative"),
            seed=0,
            speed_um_s=0.1,
        )


def test_initial_centre_is_the_lowest_force_free_axial_position() -> None:
    membrane = np.array([
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ])
    centre = np.asarray(
        derive_force_free_initial_centre(membrane, centre_contact_distance_um=2.0)
    )
    distances = np.linalg.norm(membrane - centre, axis=1)
    assert distances.min() == pytest.approx(2.0)
    assert centre[:2] == pytest.approx([0.0, 0.0])

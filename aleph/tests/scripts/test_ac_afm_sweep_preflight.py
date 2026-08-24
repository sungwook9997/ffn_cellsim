"""Tests for the fail-closed AFM sweep manifest builder."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pytest

from aleph.scripts.ac_afm_sweep_preflight import build_manifest, load_c2_ensemble


def _c2(path: Path, *, max_force: float = 0.2, subdivisions: int = 8) -> Path:
    record = {
        "schema": "afm-c2-accepted-ensemble@1",
        "status": "ACCEPTED_3_SEEDS",
        "declared_build_commit": "deadbeef",
        "configuration": {
            "membrane_subdivisions": subdivisions,
            "full_native_population": True,
            "full_compartments": True,
        },
        "accepted_seeds": [0, 1, 2],
        "records": [
            {
                "seed": seed,
                "max_projected_force_pn": max_force,
                "projected_force_tolerance_pn": 0.21,
            }
            for seed in (0, 1, 2)
        ],
        "committed_time_s_per_seed": 0.01,
    }
    path.write_text(json.dumps(record), encoding="utf-8")
    return path


def _args(path: Path, **changes: object) -> argparse.Namespace:
    values: dict[str, object] = {
        "c2_ensemble": path,
        "target_protocol": "sourced spherical protocol",
        "target_cell_line": "target line",
        "indenter_radius_um": 5.0,
        "exterior_viscosity_pa_s": 1.0e-3,
        "apparatus_source": "methods section",
        "pi0_pa": None,
        "pi0_evidence": None,
        "pi0_source": None,
        "contact_distance_um": None,
        "contact_stiffness_pn_per_um": None,
        "contact_source": None,
        "speed_um_s": [0.1],
        "depth_um": [0.0, 0.2],
        "seed": [0, 1, 2],
        "allow_mechanism_demo": False,
    }
    values.update(changes)
    return argparse.Namespace(**values)


def test_committed_ensemble_loads_and_rejected_row_refuses(tmp_path: Path) -> None:
    assert load_c2_ensemble(_c2(tmp_path / "accepted.json")).passed
    with pytest.raises(ValueError, match="force-rejected"):
        load_c2_ensemble(_c2(tmp_path / "rejected.json", max_force=0.22))


def test_missing_pi0_and_contact_emit_no_points(tmp_path: Path) -> None:
    manifest = build_manifest(_args(_c2(tmp_path / "accepted.json")))
    assert manifest["status"] == "BLOCKED"
    assert manifest["points"] == []
    assert [item.split(":", 1)[0] for item in manifest["blockers"]] == [
        "ACTIVE_CORTEX_STATE_MISSING",
        "PI0_MISSING",
        "CONTACT_CALIBRATION_MISSING",
    ]


def test_direct_sourced_inputs_still_block_without_accepted_active_cortex(tmp_path: Path) -> None:
    manifest = build_manifest(
        _args(
            _c2(tmp_path / "accepted.json"),
            pi0_pa=50.0,
            pi0_evidence="direct_target_protocol",
            pi0_source="target protocol osmotic measurement",
            contact_distance_um=0.05,
            contact_stiffness_pn_per_um=100.0,
            contact_source="declared grid-refinement study",
        )
    )
    assert manifest["status"] == "BLOCKED"
    assert manifest["points"] == []
    assert manifest["blockers"][0].startswith("ACTIVE_CORTEX_STATE_MISSING")
    assert manifest["scope"].startswith("Preflight only")


def test_mechanism_demo_can_enumerate_passive_apparatus_grid(tmp_path: Path) -> None:
    manifest = build_manifest(
        _args(
            _c2(tmp_path / "accepted.json"),
            pi0_pa=40.0,
            pi0_evidence="cross_protocol_proxy",
            pi0_source="declared proxy",
            contact_distance_um=0.05,
            contact_stiffness_pn_per_um=100.0,
            contact_source="declared grid-refinement study",
            allow_mechanism_demo=True,
        )
    )
    assert manifest["status"] == "READY_FOR_CUDA_DRIVER"
    assert manifest["result_class_requested"] == "mechanism_demo_not_quantitative"
    assert len(manifest["points"]) == 6


def test_partial_physical_declaration_is_an_error(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Pi_0 fields must be supplied together"):
        build_manifest(_args(_c2(tmp_path / "accepted.json"), pi0_pa=40.0))

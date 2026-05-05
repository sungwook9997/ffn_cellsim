from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import pytest
import yaml

from acs.v2.imaging_contract import (
    build_condition_stratified_split_manifest,
    load_imaging_contract,
    validate_imaging_contract,
)

_REPO_ROOT = Path(__file__).resolve().parents[1]
_CONTRACT_PATH = _REPO_ROOT / "configs" / "imaging" / "260313.yaml"


def test_260313_yaml_contract_loads_and_validates():
    contract = load_imaging_contract(_CONTRACT_PATH)
    result = validate_imaging_contract(contract)

    assert contract.contract_id == "pi_260313_v1"
    assert contract.voxel_size_um_xyz == (2.048, 2.048, 1.0)
    assert contract.frame_interval_s == pytest.approx(3600.0)
    assert result.group_scene_counts == {"Bare": 8, "Lam4": 26, "Pre": 25}
    assert result.group_row_counts == {"Bare": 602, "Lam4": 2104, "Pre": 1845}
    assert result.warnings == (
        "raw_imaging_available true but raw_path is not yet supplied",
    )


def test_260313_contract_promotes_to_generic_v2_contracts():
    contract = load_imaging_contract(_CONTRACT_PATH)
    generic = contract.to_v2_data_contracts()

    assert len(generic) == 3
    for item in generic:
        item.validate()
        assert item.dataset.voxel_size_um_xyz == (2.048, 2.048, 1.0)
        assert item.dataset.frame_interval_s == pytest.approx(3600.0)


def test_split_manifest_is_condition_stratified_scene_level():
    contract = load_imaging_contract(_CONTRACT_PATH)
    manifest = json.loads(contract.split_manifest_path.read_text(encoding="utf-8"))

    for condition, entry in manifest["groups"].items():
        cal = set(entry["calibration_scenes"])
        val = set(entry["validation_scenes"])
        assert cal
        assert val
        assert cal.isdisjoint(val), condition
        assert len(cal | val) == entry["n_scenes"]


def test_build_split_manifest_reproduces_committed_manifest():
    contract = load_imaging_contract(_CONTRACT_PATH)
    groups = {g.condition: g.csv_path for g in contract.groups}
    rebuilt = build_condition_stratified_split_manifest(
        groups, seed=contract.split_seed
    )
    committed = json.loads(contract.split_manifest_path.read_text(encoding="utf-8"))

    assert rebuilt == committed


def test_missing_voxel_size_fails(tmp_path):
    payload = yaml.safe_load(_CONTRACT_PATH.read_text(encoding="utf-8"))
    payload.pop("voxel_size_um_xyz")
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="voxel_size_um_xyz"):
        load_imaging_contract(bad)


def test_bad_frame_interval_fails():
    contract = load_imaging_contract(_CONTRACT_PATH)
    bad = replace(contract, frame_interval_s=0.0)

    with pytest.raises(ValueError, match="frame_interval_s"):
        validate_imaging_contract(bad)


def test_empty_channels_fail():
    contract = load_imaging_contract(_CONTRACT_PATH)
    bad = replace(contract, channels=())

    with pytest.raises(ValueError, match="channels"):
        validate_imaging_contract(bad)


def test_empty_metric_names_fail():
    contract = load_imaging_contract(_CONTRACT_PATH)
    metrics = list(contract.metrics)
    metrics[0] = replace(metrics[0], name="")
    bad = replace(contract, metrics=tuple(metrics))

    with pytest.raises(ValueError, match="metric name"):
        validate_imaging_contract(bad)


def test_calibration_validation_overlap_fails(tmp_path):
    contract = load_imaging_contract(_CONTRACT_PATH)
    manifest = json.loads(contract.split_manifest_path.read_text(encoding="utf-8"))
    manifest["groups"]["Bare"]["validation_scenes"].append(
        manifest["groups"]["Bare"]["calibration_scenes"][0]
    )
    bad_manifest = tmp_path / "split.json"
    bad_manifest.write_text(json.dumps(manifest), encoding="utf-8")
    bad = replace(contract, split_manifest_path=bad_manifest)

    with pytest.raises(ValueError, match="overlap"):
        validate_imaging_contract(bad)


def test_metric_missing_csv_column_fails():
    contract = load_imaging_contract(_CONTRACT_PATH)
    metrics = list(contract.metrics)
    metrics[0] = replace(metrics[0], column="NoSuchColumn")
    bad = replace(contract, metrics=tuple(metrics))

    with pytest.raises(ValueError, match="missing metric column"):
        validate_imaging_contract(bad)


def test_wrong_pixel_size_fails():
    contract = load_imaging_contract(_CONTRACT_PATH)
    bad = replace(contract, voxel_size_um_xyz=(1.024, 1.024, 1.0))

    with pytest.raises(ValueError, match="inferred pixel size"):
        validate_imaging_contract(bad)


def test_path_traversal_in_yaml_fails(tmp_path):
    payload = yaml.safe_load(_CONTRACT_PATH.read_text(encoding="utf-8"))
    payload["groups"][0]["csv_path"] = "../escape.csv"
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="must not contain"):
        load_imaging_contract(bad)

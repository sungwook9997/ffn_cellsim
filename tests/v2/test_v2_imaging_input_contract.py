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

_REPO_ROOT = Path(__file__).resolve().parents[2]
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
    from acs.v2.imaging_contract.split_builder import compute_yaml_sha256

    contract = load_imaging_contract(_CONTRACT_PATH)
    groups = {g.condition: g.csv_path for g in contract.groups}
    rebuilt = build_condition_stratified_split_manifest(
        groups,
        seed=contract.split_seed,
        contract_id=contract.contract_id,
        input_yaml_path="configs/imaging/260313.yaml",
        input_yaml_sha256=compute_yaml_sha256(_CONTRACT_PATH),
    )
    committed = json.loads(contract.split_manifest_path.read_text(encoding="utf-8"))

    assert rebuilt == committed


def test_split_manifest_input_yaml_sha256_stale_fails_closed(tmp_path):
    """Y13 tamper detection: manifest with mismatched input_yaml_sha256 fails."""
    contract = load_imaging_contract(_CONTRACT_PATH)
    manifest = json.loads(
        contract.split_manifest_path.read_text(encoding="utf-8")
    )
    manifest["input_yaml_sha256"] = "0" * 64
    bad_manifest = tmp_path / "split.json"
    bad_manifest.write_text(json.dumps(manifest), encoding="utf-8")
    bad_split = replace(contract.split, manifest_path=bad_manifest)
    bad = replace(contract, split=bad_split)

    with pytest.raises(ValueError, match="input_yaml_sha256"):
        validate_imaging_contract(bad)


def test_split_ratio_must_sum_to_one():
    """P2 hotfix: ratio entries that don't sum to 1.0 fail closed."""
    contract = load_imaging_contract(_CONTRACT_PATH)
    bad_split = replace(contract.split, ratio=(0.7, 0.7))
    bad = replace(contract, split=bad_split)

    with pytest.raises(ValueError, match="must sum to 1.0"):
        validate_imaging_contract(bad)


def test_missing_pixel_size_fails(tmp_path):
    payload = yaml.safe_load(_CONTRACT_PATH.read_text(encoding="utf-8"))
    payload.pop("pixel_size_um_xy")
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="pixel_size_um_xy"):
        load_imaging_contract(bad)


def test_bad_frame_interval_fails():
    contract = load_imaging_contract(_CONTRACT_PATH)
    bad = replace(contract, frame_interval_min=0.0)

    with pytest.raises(ValueError, match="frame_interval_min"):
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
    bad_split = replace(contract.split, manifest_path=bad_manifest)
    bad = replace(contract, split=bad_split)

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
    bad = replace(
        contract,
        pixel_size_um_xy=(1.024, 1.024),
        pixel_area_um2=1.024 * 1.024,
    )

    with pytest.raises(ValueError, match="inferred pixel size"):
        validate_imaging_contract(bad)


def test_pixel_area_consistency_with_pixel_size_xy():
    """Y1 P0-1: pixel_area_um2 must equal pixel_size_um_xy[0] * pixel_size_um_xy[1]."""
    contract = load_imaging_contract(_CONTRACT_PATH)
    bad = replace(contract, pixel_area_um2=5.0)  # inconsistent with [2.048, 2.048]

    with pytest.raises(ValueError, match="pixel_area_um2"):
        validate_imaging_contract(bad)


def test_per_group_max_frames_enforced(tmp_path):
    """Y15: scene with frame count exceeding max_frames_per_scene_per_group fails.

    Mutates Bare CSV to add one extra row to scene Position(102), bumping its
    frame count from 83 to 84 against the contract's amended max=83.
    """
    contract = load_imaging_contract(_CONTRACT_PATH)
    bare_csv = contract.groups[0].csv_path
    rows = bare_csv.read_text(encoding="utf-8").splitlines()
    header = rows[0]
    body = rows[1:]
    fields = header.split(",")
    series_idx = fields.index("Series")
    frame_idx = fields.index("Frame")
    time_idx = fields.index("Time_min")

    # Find the last row of the max-count scene Position(102) and append one
    # synthetic row beyond it (frame=83, time_min=4980).
    target_scene = "MD-Experiment-0004_Position(102)"
    last_idx = max(
        i for i, line in enumerate(body)
        if line.split(",")[series_idx] == target_scene
    )
    last_cols = body[last_idx].split(",")
    extra_cols = last_cols.copy()
    extra_cols[frame_idx] = "83"
    extra_cols[time_idx] = "4980.0"
    body.insert(last_idx + 1, ",".join(extra_cols))

    bad_csv = tmp_path / bare_csv.name
    bad_csv.write_text("\n".join([header, *body]) + "\n", encoding="utf-8")
    bad_groups = list(contract.groups)
    bad_groups[0] = replace(bad_groups[0], csv_path=bad_csv)
    bad = replace(contract, groups=tuple(bad_groups))

    with pytest.raises(ValueError, match="exceeds max_frames_per_scene_per_group"):
        validate_imaging_contract(bad)


def test_time_min_monotonic_required(tmp_path):
    """Y15: require_monotonic_time_min flips a scene's time order, expects fail-closed.

    Mutates one Bare scene's CSV with two Time_min values swapped, then validates.
    expected_timepoints is relaxed so the monotonic check is reached without the
    per-group max-frame check pre-empting (the per-group max is covered by a
    separate test that depends on lock-value resolution).
    """
    contract = load_imaging_contract(_CONTRACT_PATH)
    bare_csv = contract.groups[0].csv_path
    rows = bare_csv.read_text(encoding="utf-8").splitlines()
    header = rows[0]
    body = rows[1:]
    fields = header.split(",")
    time_idx = fields.index("Time_min")
    series_idx = fields.index("Series")

    # Find first two rows of the same scene and swap their Time_min values
    # so the second row's time is <= the first.
    by_scene: dict[str, list[int]] = {}
    for i, line in enumerate(body):
        cols = line.split(",")
        by_scene.setdefault(cols[series_idx], []).append(i)
    _, indices = next(iter(by_scene.items()))
    i0, i1 = indices[0], indices[1]
    cols0 = body[i0].split(",")
    cols1 = body[i1].split(",")
    cols0[time_idx], cols1[time_idx] = cols1[time_idx], cols0[time_idx]
    body[i0] = ",".join(cols0)
    body[i1] = ",".join(cols1)

    bad_csv = tmp_path / bare_csv.name
    bad_csv.write_text("\n".join([header, *body]) + "\n", encoding="utf-8")
    bad_groups = list(contract.groups)
    bad_groups[0] = replace(bad_groups[0], csv_path=bad_csv)
    relaxed = replace(
        contract.expected_timepoints,
        max_frames_per_scene_per_group={
            k: 9999
            for k in contract.expected_timepoints.max_frames_per_scene_per_group
        },
        # also disable monotonic-frame so swapping Time_min doesn't trip
        # the Frame check first.
        require_monotonic_frame=False,
    )
    bad = replace(
        contract,
        groups=tuple(bad_groups),
        expected_timepoints=relaxed,
    )

    with pytest.raises(ValueError, match="not strictly monotonic"):
        validate_imaging_contract(bad)


def test_split_block_forbidden_use_physics_fitting_required(tmp_path):
    """Y14 / Hard Rule 1 wording-boundary meta-test: missing forbidden_use fails."""
    payload = yaml.safe_load(_CONTRACT_PATH.read_text(encoding="utf-8"))
    payload["calibration_validation_split"].pop("forbidden_use")
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump(payload), encoding="utf-8")

    with pytest.raises(
        ValueError, match="forbidden_use.*physics_parameter_fitting"
    ):
        load_imaging_contract(bad)


def test_split_block_forbidden_use_wrong_value_rejected(tmp_path):
    """Y14: forbidden_use with wrong value (not 'physics_parameter_fitting') fails."""
    payload = yaml.safe_load(_CONTRACT_PATH.read_text(encoding="utf-8"))
    payload["calibration_validation_split"]["forbidden_use"] = "something_else"
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="physics_parameter_fitting"):
        load_imaging_contract(bad)


def test_sidecar_model2_framewise_forbidden_dynamics_fitting_required(tmp_path):
    """Y9 / Hard Rule 1 wording-boundary meta-test: missing sidecar forbidden_use fails."""
    payload = yaml.safe_load(_CONTRACT_PATH.read_text(encoding="utf-8"))
    payload["sidecars"]["model2_framewise"].pop("forbidden_use")
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump(payload), encoding="utf-8")

    with pytest.raises(
        ValueError, match="forbidden_use.*dynamics_parameter_fitting"
    ):
        load_imaging_contract(bad)


def test_sidecar_model2_framewise_forbidden_use_wrong_value_rejected(tmp_path):
    """Y9: sidecar forbidden_use with wrong value fails."""
    payload = yaml.safe_load(_CONTRACT_PATH.read_text(encoding="utf-8"))
    payload["sidecars"]["model2_framewise"]["forbidden_use"] = "something_else"
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="dynamics_parameter_fitting"):
        load_imaging_contract(bad)


def test_path_traversal_in_yaml_fails(tmp_path):
    payload = yaml.safe_load(_CONTRACT_PATH.read_text(encoding="utf-8"))
    payload["groups"][0]["csv_path"] = "../escape.csv"
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="must not contain"):
        load_imaging_contract(bad)

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("screen_second", HERE / "screen_second_pass.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def record(family: str, group: str | None = None) -> dict:
    return {
        "source_family_id": family,
        "leakage": {"dataset_groups": [] if group is None else [{"group_hash": group}]},
    }


def test_combined_leakage_distinguishes_cross_pass_groups():
    raw, counts = MODULE.leakage_groups(
        [record("doi:first-a", "a" * 16), record("doi:first-b", "b" * 16)],
        [record("doi:second-a", "a" * 16), record("doi:second-b", "b" * 16), record("doi:second-c", "b" * 16)],
    )
    assert counts == {"duplicate_groups": 2, "cross_pass_groups": 2, "within_first_groups": 0, "within_second_groups": 1}
    assert raw["authority_status"] == "proposed"
    assert all("text" not in json.dumps(group).lower() for group in raw["groups"])


def test_base_policy_is_loaded_read_only():
    assert MODULE.BASE_PATH == HERE.parent / "screen_oa.py"
    assert set(MODULE.BASE.SIGNALS) == {"sample_size", "biological_replicates", "calibration", "units", "uncertainty", "exclusions", "source_data"}
    assert MODULE.EXPECTED_ARTICLES == 2974

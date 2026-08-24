"""Accepted-path aggregation uses seed scatter and never exposes rejected candidates."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from aleph.scripts.ac_afm_sweep_aggregate import aggregate
from aleph.scripts.ac_afm_sweep_vis import render


def _path(path: Path, seed: int, *, accepted: bool = True) -> Path:
    points = [
        {
            "depth_um": depth,
            "reaction_z_pn": reaction + seed,
            "max_projected_force_pn": 0.20,
            "projected_force_tolerance_pn": 0.21,
            "accepted": True,
        }
        for depth, reaction in ((0.0, 0.0), (0.2, 10.0))
    ]
    payload = {
        "schema": "afm-force-indentation-path@1",
        "status": "ACCEPTED_PATH" if accepted else "REJECTED_ROLLED_BACK",
        "quantitative_claim": "MECHANISM_DEMO_NOT_QUANTITATIVE",
        "build_commit": "deadbeef",
        "seed": seed,
        "speed_um_s": 0.1,
        "protocol": {"reaction_channel": 0},
        "population": {"n_total_nodes": 100},
        "points": points if accepted else [],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_three_seed_paths_aggregate_with_sample_sd_and_render(tmp_path: Path) -> None:
    paths = [_path(tmp_path / f"seed-{seed}.json", seed) for seed in (0, 1, 2)]
    result = aggregate(paths)
    assert result["status"] == "ACCEPTED_ENSEMBLE"
    assert result["points"][1]["reaction_mean_pn"] == pytest.approx(11.0)
    assert result["points"][1]["reaction_sample_sd_pn"] == pytest.approx(1.0)
    record = tmp_path / "ensemble.json"
    record.write_text(json.dumps(result), encoding="utf-8")
    figure = tmp_path / "ensemble.png"
    render(record, figure)
    assert figure.read_bytes().startswith(b"\x89PNG")


def test_under_replicated_or_rejected_paths_never_form_a_curve(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="at least three"):
        aggregate([_path(tmp_path / f"few-{seed}.json", seed) for seed in (0, 1)])
    with pytest.raises(ValueError, match="ACCEPTED"):
        aggregate([
            _path(tmp_path / "a.json", 0),
            _path(tmp_path / "b.json", 1),
            _path(tmp_path / "bad.json", 2, accepted=False),
        ])

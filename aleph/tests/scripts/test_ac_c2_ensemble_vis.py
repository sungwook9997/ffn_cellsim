"""C-2 ensemble figure remains a numerical acceptance visualization."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from aleph.scripts.ac_c2_ensemble_vis import render


def _record(path: Path, forces: tuple[float, ...] = (0.19, 0.20, 0.205)) -> Path:
    rows = [
        {
            "seed": seed,
            "reported_iterations": 200 + seed * 20,
            "max_projected_force_pn": force,
            "projected_force_tolerance_pn": 0.21,
        }
        for seed, force in enumerate(forces)
    ]
    path.write_text(
        json.dumps({"schema": "afm-c2-accepted-ensemble@1", "records": rows}),
        encoding="utf-8",
    )
    return path


def test_render_writes_accepted_ensemble_png(tmp_path: Path) -> None:
    output = tmp_path / "ensemble.png"
    render(_record(tmp_path / "record.json"), output)
    assert output.read_bytes().startswith(b"\x89PNG")


def test_render_refuses_rejected_or_under_replicated_evidence(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="at least three"):
        render(_record(tmp_path / "few.json", (0.19, 0.20)), tmp_path / "few.png")
    with pytest.raises(ValueError, match="force-rejected"):
        render(_record(tmp_path / "bad.json", (0.19, 0.20, 0.22)), tmp_path / "bad.png")

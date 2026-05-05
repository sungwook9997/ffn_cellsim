"""Deterministic scene-level split helpers for V2-1 imaging contracts."""

from __future__ import annotations

import csv
import random
from pathlib import Path


def _scene_ids_from_csv(csv_path: Path, *, series_column: str) -> tuple[str, ...]:
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or series_column not in reader.fieldnames:
            raise ValueError(
                f"{csv_path} missing required scene column {series_column!r}"
            )
        scenes = sorted({row[series_column] for row in reader if row[series_column]})
    if not scenes:
        raise ValueError(f"{csv_path} contains no scene ids in {series_column!r}")
    return tuple(scenes)


def build_condition_stratified_split_manifest(
    groups: dict[str, Path],
    *,
    seed: int,
    calibration_fraction: float = 0.7,
    series_column: str = "Series",
) -> dict:
    """Build a deterministic 70/30-ish split per condition.

    The split is scene-level, not row/timepoint-level, so rows from one
    scene cannot leak between calibration and validation.
    """

    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an int")
    if not (0.0 < float(calibration_fraction) < 1.0):
        raise ValueError("calibration_fraction must be in (0, 1)")

    manifest = {
        "schema_version": "v2_imaging_split_manifest_v1",
        "seed": seed,
        "split_policy": "condition_stratified_scene_level_70_30",
        "groups": {},
    }
    for condition, csv_path in sorted(groups.items()):
        scenes = list(_scene_ids_from_csv(Path(csv_path), series_column=series_column))
        rng = random.Random(f"{seed}:{condition}")
        shuffled = scenes[:]
        rng.shuffle(shuffled)
        n_cal = max(1, round(len(shuffled) * calibration_fraction))
        if len(shuffled) > 1 and len(shuffled) - n_cal == 0:
            n_cal = len(shuffled) - 1
        manifest["groups"][condition] = {
            "n_scenes": len(shuffled),
            "calibration_scenes": sorted(shuffled[:n_cal]),
            "validation_scenes": sorted(shuffled[n_cal:]),
        }
    return manifest

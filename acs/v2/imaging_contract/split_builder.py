"""Deterministic scene-level split helpers for V2-1 imaging contracts."""

from __future__ import annotations

import csv
import hashlib
import random
from pathlib import Path


_MANIFEST_SCHEMA_VERSION = "v2_imaging_split_manifest_v1"
_MANIFEST_VERSION = 1
_GENERATOR = "acs.v2.imaging_contract.split_builder"


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


def compute_yaml_sha256(yaml_path: Path) -> str:
    """SHA-256 hex of the YAML file bytes (Y13 tamper detection)."""
    return hashlib.sha256(Path(yaml_path).read_bytes()).hexdigest()


def build_condition_stratified_split_manifest(
    groups: dict[str, Path],
    *,
    seed: int,
    calibration_fraction: float = 0.7,
    series_column: str = "Series",
    contract_id: str,
    input_yaml_path: str,
    input_yaml_sha256: str,
) -> dict:
    """Build a deterministic 70/30-ish split per condition.

    The split is scene-level, not row/timepoint-level, so rows from one
    scene cannot leak between calibration and validation.

    The manifest carries Y13 tamper-detection metadata: ``input_yaml`` plus
    ``input_yaml_sha256`` pin the manifest to a specific YAML file content.
    The loader fails closed if the YAML changes after the manifest is
    generated; regenerate via this builder rather than editing the JSON.
    ``generated_at_iso`` is ``None`` so manifest bytes are reproducible
    (Y13: HB#5/Item5 sister, no wall-clock noise).
    """

    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an int")
    if not (0.0 < float(calibration_fraction) < 1.0):
        raise ValueError("calibration_fraction must be in (0, 1)")

    manifest = {
        "schema_version": _MANIFEST_SCHEMA_VERSION,
        "manifest_version": _MANIFEST_VERSION,
        "generator": _GENERATOR,
        "generated_at_iso": None,
        "contract_id": contract_id,
        "input_yaml": input_yaml_path,
        "input_yaml_sha256": input_yaml_sha256,
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

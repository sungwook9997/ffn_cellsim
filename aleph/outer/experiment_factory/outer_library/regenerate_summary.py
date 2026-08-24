#!/usr/bin/env python3
"""Regenerate the complete hash/count index for local Outer Library results."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def build(results: Path) -> dict[str, object]:
    manifests = sorted(results.glob("*_observations.jsonl"))
    tensors = sorted(results.glob("*.npz"))
    reports = sorted(
        path for path in results.glob("*.json") if path.name != "summary.json"
    )
    datasets: set[str] = set()
    labs: set[str] = set()
    samples: set[tuple[str, str]] = set()
    modalities: dict[str, int] = {}
    observations = 0
    for path in manifests:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                row = json.loads(line)
                observations += 1
                datasets.add(row["dataset_id"])
                labs.add(row["lab_group"])
                samples.add((row["dataset_id"], row["sample_id"]))
                modality = row["modality"]
                modalities[modality] = modalities.get(modality, 0) + 1
    multitask_path = results / "multitask_outer_model_report.json"
    production_eligible = False
    if multitask_path.exists():
        production_eligible = bool(
            json.loads(multitask_path.read_text(encoding="utf-8"))["readiness"][
                "production_eligible"
            ]
        )
    return {
        "dataset_count": len(datasets),
        "lab_group_count": len(labs),
        "observation_count": observations,
        "sample_count": len(samples),
        "modalities": modalities,
        "generated_manifest_sha256": {
            path.name: _sha256(path) for path in manifests
        },
        "generated_tensor_sha256": {
            path.name: _sha256(path) for path in tensors
        },
        "generated_report_sha256": {
            path.name: _sha256(path) for path in reports
        },
        "authority": "none",
        "production_eligible": production_eligible,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = build(args.results)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "datasets": summary["dataset_count"],
        "labs": summary["lab_group_count"],
        "samples": summary["sample_count"],
        "observations": summary["observation_count"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

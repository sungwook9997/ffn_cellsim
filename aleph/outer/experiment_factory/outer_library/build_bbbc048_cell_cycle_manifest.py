#!/usr/bin/env python3
"""Build a leakage-safe cell-level manifest for official BBBC048 labels."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import zipfile


LABEL_CODE = {
    "Anaphase": 0, "Metaphase": 1, "Prophase": 2, "Telophase": 3,
    "G1": 4, "G2": 5, "S": 6,
}
CHANNELS = {"Ch3": "brightfield", "Ch4": "MPM2", "Ch6": "propidium_iodide_DNA"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build(archive_path: Path, ground_truth_path: Path, output_path: Path) -> dict[str, object]:
    grouped: dict[tuple[str, str], dict[str, str]] = defaultdict(dict)
    line_count = 0
    for line in ground_truth_path.read_text(encoding="utf-8").splitlines():
        _, code_text, member = line.split("\t")
        normalized = member.removeprefix("./")
        phase, filename = normalized.split("/", 1)
        stem, channel_suffix = filename.rsplit("_", 1)
        channel = channel_suffix.removesuffix(".ome.jpg")
        if phase not in LABEL_CODE or int(code_text) != LABEL_CODE[phase]:
            raise ValueError(f"BBBC048 label/path mismatch: {line!r}")
        if channel not in CHANNELS:
            raise ValueError(f"unexpected BBBC048 channel {channel!r}")
        grouped[(phase, stem)][channel] = normalized
        line_count += 1
    if line_count != 96_798 or len(grouped) != 32_266:
        raise ValueError(f"BBBC048 cardinality changed: {line_count} lines, {len(grouped)} cells")
    with zipfile.ZipFile(archive_path) as archive:
        members = {
            name.removeprefix("./").removeprefix("CellCycle/")
            for name in archive.namelist() if not name.endswith("/")
        }
        missing = [member for channels in grouped.values() for member in channels.values() if member not in members]
    if missing:
        raise ValueError(f"BBBC048 archive is missing {len(missing)} labelled images")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    with output_path.open("w", encoding="utf-8") as handle:
        for (phase, cell_id), channels in sorted(grouped.items()):
            if set(channels) != set(CHANNELS):
                raise ValueError(f"cell {phase}/{cell_id} does not have exactly three channels")
            counts[phase] += 1
            row = {
                "schema": "aleph.outer_library.bbbc048_cell.v1",
                "dataset_id": "bbbc048-v1-jurkat-cell-cycle",
                "sample_id": f"BBBC048:{phase}:{cell_id}",
                "split_group": f"BBBC048:{phase}:{cell_id}",
                "lab_group": "Newcastle_FCCF_Eulenberg_2017",
                "modality": "imaging_flow_cytometry_multichannel",
                "cell_type": "Jurkat_T_lymphocyte",
                "cell_state": phase,
                "cell_state_axis": "cell_cycle_phase",
                "author_label": phase,
                "author_label_code": LABEL_CODE[phase],
                "unit": "categorical_author_label",
                "channels": [{
                    "channel_id": channel,
                    "biological_role": CHANNELS[channel],
                    "archive_member": channels[channel],
                } for channel in sorted(CHANNELS)],
                "biological_replicate": "not_reported",
                "provider_lab_group": "Newcastle_FCCF_Eulenberg_2017",
                "external_lab_holdout": False,
                "training_role": "same_provider_representation_and_cell_cycle_supervision",
                "license_id": "cc-by-nc-sa-3.0",
                "aleph_authority": "none",
                "may_select_aleph_parameter": False,
            }
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    return {
        "schema": "aleph.outer_library.bbbc048_manifest_receipt.v1",
        "cell_count": len(grouped),
        "image_count": line_count,
        "phase_counts": dict(sorted(counts.items())),
        "channel_count": len(CHANNELS),
        "archive_sha256": _sha256(archive_path),
        "ground_truth_sha256": _sha256(ground_truth_path),
        "manifest_sha256": _sha256(output_path),
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    receipt = build(args.archive, args.ground_truth, args.output)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()

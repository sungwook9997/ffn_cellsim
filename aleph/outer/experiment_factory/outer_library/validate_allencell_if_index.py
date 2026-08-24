#!/usr/bin/env python3
"""Validate the preregistered Allen Cell feature-index schema before images."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def normalized(value: str) -> str:
    return re.sub(r"[^0-9a-z]+", "", value.casefold())


def resolve(headers: list[str], aliases: list[str]) -> list[str]:
    wanted = set(aliases)
    return [header for header in headers if normalized(header) in wanted]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    contract = protocol["feature_index_contract"]
    with zipfile.ZipFile(args.index) as archive:
        csv_names = [name for name in archive.namelist() if name.casefold().endswith(".csv") and not name.startswith("__MACOSX/")]
        if len(csv_names) != 1:
            raise ValueError(f"expected exactly one data CSV, found {csv_names}")
        payload = archive.read(csv_names[0]).decode("utf-8-sig")
    reader = csv.reader(io.StringIO(payload))
    headers = next(reader)
    row_count = sum(1 for _ in reader)
    matches = {
        "cell_identifier": resolve(headers, contract["cell_id_aliases"]),
        "cell_line_identifier": resolve(headers, contract["line_id_aliases"]),
        "tagged_structure": resolve(headers, contract["structure_aliases"]),
    }
    valid = all(len(value) == 1 for value in matches.values())
    missing = [key for key, value in matches.items() if not value]
    ambiguous = [key for key, value in matches.items() if len(value) > 1]
    report = {
        "schema": "aleph.outer_library.allencell_if_schema_refusal.v1",
        "status": "schema_valid" if valid else "refused_before_image_download_and_prediction",
        "protocol_sha256": sha256(args.protocol),
        "feature_index_sha256": sha256(args.index),
        "archive_csv": csv_names[0],
        "row_count": row_count,
        "column_count": len(headers),
        "headers": headers,
        "normalized_headers": [normalized(header) for header in headers],
        "resolved_columns": matches,
        "missing_required_fields": missing,
        "ambiguous_required_fields": ambiguous,
        "image_objects_downloaded": 0,
        "model_predictions_made": 0,
        "external_labels_used_for_model_selection": False,
        "independent_provider_domain_shift_validated": False,
        "uncertainty_authority": False,
        "aleph_parameter_authority": "none",
        "may_emit_numeric_sweep_range": False,
        "reason": (
            "The frozen exact-alias resolver cannot identify the Allen cell-line and tagged-structure fields; "
            "manual post-access column selection or label joins are forbidden by the preregistration."
            if not valid else None
        ),
        "next_action": (
            "preregister a different independent-provider raw-IF source whose released index carries explicit labels"
            if not valid else "resolve crop objects and channels under the frozen image contract"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "rows": row_count, "missing": missing}, sort_keys=True))


if __name__ == "__main__":
    main()

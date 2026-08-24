#!/usr/bin/env python3
"""Acquire IDR0072 annotations and freeze a schema-informed image manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.error
import urllib.request


def fetch_json(url: str, attempts: int = 4) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "Project-Aleph-research/1.0"})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.load(response)
        except (urllib.error.URLError, TimeoutError):
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def normalized_qualifier(header: str) -> str:
    bracketed = re.findall(r"\[([^][]+)\]", header)
    source = bracketed[0] if len(bracketed) == 1 else header
    return re.sub(r"[^0-9a-z]+", "", source.casefold())


def normalized_value(value: object) -> str:
    return " " + re.sub(r"[^0-9a-z]+", " ", str(value).casefold()).strip() + " "


def resolve_columns(columns: list[str], contract: dict) -> dict[str, int]:
    resolved: dict[str, int] = {}
    for semantic, qualifier in contract["required_exact_qualifiers"].items():
        hits = [index for index, header in enumerate(columns) if normalized_qualifier(header) == qualifier]
        if len(hits) != 1:
            raise ValueError(f"{semantic} resolved to {hits}")
        resolved[semantic] = hits[0]
    plain = contract["required_exact_plain_column"]
    hits = [index for index, header in enumerate(columns) if normalized_qualifier(header) == plain]
    if len(hits) != 1:
        raise ValueError(f"channels resolved to {hits}")
    resolved["channels"] = hits[0]
    return resolved


def map_class(value: object, mapping: dict) -> tuple[str | None, list[str]]:
    text = normalized_value(value)
    hits = [
        label for label in mapping["precedence"]
        if any(token in text for token in mapping[label])
    ]
    return (hits[0] if len(hits) == 1 else None), hits


def salted_key(salt: str, *parts: object) -> str:
    joined = "\x1f".join([salt, *(str(part) for part in parts)])
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    protocol_bytes = args.protocol.read_bytes()
    protocol = json.loads(protocol_bytes)
    if not protocol["frozen_after_one_metadata_schema_row_but_before_remaining_annotations_images_or_predictions"]:
        raise ValueError("protocol is not frozen at the declared boundary")
    source = protocol["source_contract"]
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    raw_receipts: list[dict] = []
    class_counts: dict[str, int] = {}
    refusal_counts: dict[str, int] = {}
    plate_grids: dict[int, dict[int, dict]] = {}

    for screen_id in source["screens"]:
        bulk_url = f"{source['api_base']}/webgateway/table/Screen/{screen_id}/query/?query=%2A"
        bulk = fetch_json(bulk_url)
        bulk_payload = (json.dumps(bulk, separators=(",", ":"), sort_keys=True) + "\n").encode("utf-8")
        bulk_path = args.cache_dir / f"screen_{screen_id}_bulk.json"
        bulk_path.write_bytes(bulk_payload)
        raw_receipts.append({"kind": "bulk_annotation", "screen_id": screen_id, "url": bulk_url, "bytes": len(bulk_payload), "sha256": sha256_bytes(bulk_payload)})
        data = bulk["data"]
        columns = data["columns"]
        resolved = resolve_columns(columns, protocol["isa_tab_column_contract"])
        plate_index = columns.index("Plate")
        well_index = columns.index("Well")
        plate_name_index = columns.index("Plate Name")
        well_name_index = columns.index("Well Name")
        localization_index = resolved["localization"]

        for plate_id in sorted({int(row[plate_index]) for row in data["rows"]}):
            grid_url = f"{source['api_base']}/webgateway/plate/{plate_id}/0/"
            grid = fetch_json(grid_url)
            grid_payload = (json.dumps(grid, separators=(",", ":"), sort_keys=True) + "\n").encode("utf-8")
            grid_path = args.cache_dir / f"plate_{plate_id}_field0.json"
            grid_path.write_bytes(grid_payload)
            raw_receipts.append({"kind": "plate_grid", "plate_id": plate_id, "url": grid_url, "bytes": len(grid_payload), "sha256": sha256_bytes(grid_payload)})
            plate_grids[plate_id] = {}
            for grid_row in grid["grid"]:
                for cell in grid_row:
                    if cell is not None:
                        well_id = int(cell["wellId"])
                        if well_id in plate_grids[plate_id]:
                            raise ValueError(f"duplicate field-0 well {well_id} on plate {plate_id}")
                        plate_grids[plate_id][well_id] = cell

        for row in data["rows"]:
            canonical, hits = map_class(row[localization_index], protocol["fixed_value_mapping"])
            if canonical is None:
                key = "unmapped" if not hits else "ambiguous"
                refusal_counts[key] = refusal_counts.get(key, 0) + 1
                continue
            plate_id, well_id = int(row[plate_index]), int(row[well_index])
            image = plate_grids.get(plate_id, {}).get(well_id)
            if image is None:
                refusal_counts["missing_field0_image"] = refusal_counts.get("missing_field0_image", 0) + 1
                continue
            record = {
                "schema": "aleph.outer_library.idr0072_if_image_manifest.v1",
                "screen_id": screen_id,
                "plate_id": plate_id,
                "plate_name": row[plate_name_index],
                "well_id": well_id,
                "well_name": row[well_name_index],
                "image_id": int(image["id"]),
                "image_name": image["name"],
                "cell_line": row[resolved["cell_line"]],
                "gene_symbol": row[resolved["gene_symbol"]],
                "protein_name": row[resolved["protein_name"]],
                "localization_raw": row[localization_index],
                "canonical_class": canonical,
                "channels_raw": row[resolved["channels"]],
                "selection_key": salted_key(protocol["sampling_contract"]["salt"], screen_id, plate_id, well_id, image["id"]),
                "split_group": f"screen:{screen_id}|plate:{plate_id}|protein:{row[resolved['gene_symbol']]}",
            }
            all_rows.append(record)
            class_counts[canonical] = class_counts.get(canonical, 0) + 1

    selected: list[dict] = []
    grouped: dict[tuple, list[dict]] = {}
    for row in all_rows:
        key = (row["screen_id"], row["plate_id"], row["gene_symbol"], row["canonical_class"])
        grouped.setdefault(key, []).append(row)
    maximum = protocol["sampling_contract"]["maximum_images_per_landmark_protein_and_plate"]
    for rows in grouped.values():
        selected.extend(sorted(rows, key=lambda row: row["selection_key"])[:maximum])
    selected.sort(key=lambda row: row["selection_key"])
    proteins = {row["gene_symbol"] for row in selected}
    selected_classes = {row["canonical_class"] for row in selected}
    constraints = {
        "minimum_primary_classes": len(selected_classes) >= protocol["sampling_contract"]["minimum_primary_classes"],
        "minimum_landmark_proteins": len(proteins) >= protocol["sampling_contract"]["minimum_landmark_proteins"],
        "minimum_total_images": len(selected) >= protocol["sampling_contract"]["minimum_total_images"],
    }
    status = "manifest_frozen_before_image_content" if all(constraints.values()) else "refused_before_image_content"
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in selected), encoding="utf-8")
    report = {
        "schema": "aleph.outer_library.idr0072_if_manifest_report.v1",
        "status": status,
        "protocol_sha256": sha256_bytes(protocol_bytes),
        "annotation_rows": len(all_rows) + sum(refusal_counts.values()),
        "uniquely_mapped_field0_rows": len(all_rows),
        "available_class_counts": class_counts,
        "refusal_counts": refusal_counts,
        "selected_images": len(selected),
        "selected_classes": sorted(selected_classes),
        "selected_class_counts": {label: sum(row["canonical_class"] == label for row in selected) for label in sorted(selected_classes)},
        "selected_landmark_proteins": len(proteins),
        "selected_plates": len({row["plate_id"] for row in selected}),
        "selected_screens": sorted({row["screen_id"] for row in selected}),
        "constraints": constraints,
        "raw_metadata_receipts": raw_receipts,
        "manifest_sha256": sha256_bytes(args.manifest.read_bytes()),
        "image_objects_downloaded": 0,
        "model_predictions_made": 0,
        "pristine_holdout_authority": False,
        "uncertainty_authority": False,
        "aleph_parameter_authority": "none",
    }
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "selected": len(selected), "classes": len(selected_classes), "proteins": len(proteins)}, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Deterministically sample Allen raw hiPSC metadata without image access."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def request(url: str, *, method: str = "GET", byte_range: str | None = None):
    headers = {"User-Agent": "Project-Aleph-research/1.0"}
    if byte_range is not None:
        headers["Range"] = byte_range
    return urllib.request.urlopen(
        urllib.request.Request(url, method=method, headers=headers), timeout=120
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    source = protocol["source"]
    sampling = protocol["sampling"]
    with request(source["metadata_url"], method="HEAD") as response:
        content_length = int(response.headers["Content-Length"])
        etag = response.headers.get("ETag")
        last_modified = response.headers.get("Last-Modified")
    chunk_size = int(sampling["range_chunk_bytes"])
    chunk_count = int(sampling["range_chunk_count"])
    maximum_offset = content_length - chunk_size
    offsets = [
        round(index * maximum_offset / (chunk_count - 1))
        for index in range(chunk_count)
    ]
    header_end = 65535
    with request(source["metadata_url"], byte_range=f"bytes=0-{header_end}") as response:
        header = response.read().splitlines()[0].decode("utf-8")
    columns = next(csv.reader([header]))
    required = {
        "CellId", "FOVId", "structure_name", "outlier", "cell_stage",
        "crop_raw", "fov_path", "ChannelNumberBrightfield",
        "ChannelNumberStruct", "InstrumentId", "WorkflowId", "ProtocolId",
        "PlateId", "WellId", "roi", "scale_micron",
    }
    missing = sorted(required - set(columns))
    if missing:
        raise ValueError(f"Allen metadata missing required columns: {missing}")

    observed_rows = eligible_rows = malformed_rows = 0
    candidates: dict[str, dict[str, str]] = {}
    range_receipts = []
    for offset in offsets:
        end = min(content_length - 1, offset + chunk_size - 1)
        with request(source["metadata_url"], byte_range=f"bytes={offset}-{end}") as response:
            payload = response.read()
            status = response.status
            content_range = response.headers.get("Content-Range")
        if status != 206 or content_range is None:
            raise ValueError("source did not honor deterministic byte-range request")
        lines = payload.splitlines()
        if offset != 0 and lines:
            lines = lines[1:]
        if end != content_length - 1 and lines:
            lines = lines[:-1]
        range_receipts.append({
            "offset": offset, "end": end, "bytes_received": len(payload),
            "content_range": content_range,
        })
        for line in lines:
            observed_rows += 1
            try:
                values = next(csv.reader(io.StringIO(line.decode("utf-8"))))
            except (UnicodeDecodeError, csv.Error):
                malformed_rows += 1
                continue
            if len(values) != len(columns):
                malformed_rows += 1
                continue
            row = dict(zip(columns, values, strict=True))
            try:
                eligible = (
                    row["structure_name"] == "TOMM20"
                    and row["outlier"] == "No"
                    and row["cell_stage"] == "M0"
                    and bool(row["crop_raw"] and row["fov_path"])
                    and int(row["ChannelNumberBrightfield"]) > 0
                    and int(row["ChannelNumberStruct"]) > 0
                )
            except ValueError:
                eligible = False
            if not eligible:
                continue
            eligible_rows += 1
            fov_id = row["FOVId"]
            candidate = {key: row[key] for key in sorted(required)}
            current = candidates.get(fov_id)
            if current is None or int(candidate["CellId"]) < int(current["CellId"]):
                candidates[fov_id] = candidate
    ordered = sorted(
        candidates.values(),
        key=lambda row: hashlib.sha256(
            f"20260808:{row['FOVId']}:{row['CellId']}".encode()
        ).hexdigest(),
    )
    requested = int(sampling["candidate_count"])
    if len(ordered) < requested:
        raise ValueError(
            f"frozen sampling found {len(ordered)} distinct eligible FOVs, need {requested}"
        )
    selected = ordered[:requested]
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.manifest.open("w", encoding="utf-8") as handle:
        for row in selected:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    report = {
        "schema": "aleph.outer_library.allencell_label_free_metadata_audit_result.v1",
        "protocol_sha256": sha256(args.protocol),
        "metadata_remote": {
            "content_length": content_length, "etag": etag,
            "last_modified": last_modified,
        },
        "range_request_count": len(range_receipts),
        "range_bytes_read": sum(item["bytes_received"] for item in range_receipts),
        "range_receipts": range_receipts,
        "complete_rows_observed": observed_rows,
        "malformed_or_partial_rows_excluded": malformed_rows,
        "eligible_rows_observed": eligible_rows,
        "distinct_eligible_FOVs": len(candidates),
        "selected_candidate_count": len(selected),
        "selected_manifest_sha256": sha256(args.manifest),
        "image_headers_accessed": 0,
        "image_pixels_accessed": 0,
        "model_predictions_made": 0,
        "external_confirmation_complete": False,
        "aleph_parameter_authority": "none",
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()

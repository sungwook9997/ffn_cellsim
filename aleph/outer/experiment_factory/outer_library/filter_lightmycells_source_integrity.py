#!/usr/bin/env python3
"""Refuse frozen Light My Cells pairs whose provider payload fails integrity."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import urllib.request


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--integrity-protocol", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.integrity_protocol.read_text(encoding="utf-8"))
    pinned = protocol["v2_pair_manifest"]
    if Path(pinned["path"]).name != args.input.name or pinned["sha256"] != sha256(args.input):
        raise ValueError("integrity protocol does not pin supplied v2 pair manifest")
    observed = protocol["observed_source_inconsistency"]
    bad_path = observed["path"]
    rows = [json.loads(line) for line in args.input.read_text(encoding="utf-8").splitlines()]
    referenced = [
        row for row in rows
        if bad_path in (row["input"]["path"], row["target"]["path"])
    ]
    if len(referenced) != 1:
        raise ValueError(f"expected exactly one pair for frozen bad payload, found {len(referenced)}")
    bad_file = (
        referenced[0]["input"] if referenced[0]["input"]["path"] == bad_path
        else referenced[0]["target"]
    )
    if bad_file["size_bytes"] != observed["file_list_size_bytes"]:
        raise ValueError("v2 pair no longer carries frozen FileList size")
    request = urllib.request.Request(
        bad_file["url"], method="HEAD",
        headers={"User-Agent": "Project-Aleph-research/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        remote_length = int(response.headers["Content-Length"])
        remote_etag = response.headers.get("ETag")
    if remote_length != observed["remote_content_length_bytes"]:
        raise ValueError("provider payload changed; do not silently revise frozen refusal")
    if remote_etag != observed["remote_etag"]:
        raise ValueError("provider ETag changed; do not silently revise frozen refusal")

    admitted = [row for row in rows if row not in referenced]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in admitted:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    unique = {
        file["url"]: file for row in admitted for file in (row["input"], row["target"])
    }
    report = {
        "schema": "aleph.outer_library.lightmycells_source_integrity_filter.v1",
        "integrity_protocol_sha256": sha256(args.integrity_protocol),
        "input_manifest_sha256": sha256(args.input),
        "output_manifest_sha256": sha256(args.output),
        "input_pair_count": len(rows),
        "admitted_pair_count": len(admitted),
        "refused_pair_count": len(referenced),
        "refused_pairs": [{
            "study_id": row["study_id"],
            "acquisition_set_id": row["acquisition_set_id"],
            "target_channel": row["target_channel"],
            "path": bad_path,
            "reason": "provider Content-Length differs from frozen FileList and contains metadata-only BigTIFF",
        } for row in referenced],
        "unique_admitted_payload_count": len(unique),
        "admitted_payload_bytes": sum(file["size_bytes"] for file in unique.values()),
        "role_pair_counts": dict(sorted(Counter(row["role"] for row in admitted).items())),
        "target_pair_counts": dict(sorted(Counter(row["target_channel"] for row in admitted).items())),
        "model_predictions_made": 0,
        "confirmation_header_requests": 0,
        "confirmation_pixels_accessed": 0,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "pairs": len(admitted), "payloads": len(unique),
        "bytes": report["admitted_payload_bytes"], "refused": len(referenced),
    }, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Seal an IDR0168 source-integrity refusal without making model predictions."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import tarfile
import urllib.request

from train_hpa_raw_if_cnn import sha256


USER_AGENT = "Project-Aleph-research/1.0"
RANGE_LENGTH = 1024 * 1024


def remote_range(url: str, start: int, end: int) -> bytes:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Range": f"bytes={start}-{end}",
        },
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        if response.status != 206:
            raise ValueError(f"range request was not honored: {response.status}")
        value = response.read()
    if len(value) != end - start + 1:
        raise ValueError("range response length changed")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument("--archive-name", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    acquisition = json.loads(args.acquisition.read_text(encoding="utf-8"))
    if acquisition["protocol_sha256"] != sha256(args.protocol):
        raise ValueError("IDR0168 acquisition does not pin supplied protocol")
    if protocol["pixel_contract"]["decode_or_channel_mismatch_action"] != (
        "refuse the entire provider without replacement"
    ):
        raise ValueError("IDR0168 frozen decode-failure action changed")
    rows = [
        row for row in acquisition["archives"]
        if row["local_name"] == args.archive_name
    ]
    if len(rows) != 1:
        raise ValueError("failed archive is not unique in acquisition receipt")
    row = rows[0]
    path = args.archive_dir / args.archive_name
    if path.stat().st_size != row["size_bytes"] or sha256(path) != row["sha256"]:
        raise ValueError("failed archive no longer matches acquisition receipt")

    gzip_error = None
    try:
        with gzip.open(path, "rb") as source:
            while source.read(8 * 1024 * 1024):
                pass
    except (gzip.BadGzipFile, EOFError, OSError) as exc:
        gzip_error = f"{type(exc).__name__}: {exc}"
    if gzip_error is None:
        raise ValueError("failed archive unexpectedly passes gzip integrity")

    member_names = []
    tar_error = None
    try:
        with tarfile.open(path, mode="r:gz") as archive:
            for member in archive:
                if member.name.endswith("/.zarray"):
                    member_names.append(member.name)
    except (tarfile.TarError, EOFError, OSError) as exc:
        tar_error = f"{type(exc).__name__}: {exc}"
    expected_native_metadata = "./0/0/.zarray"
    if expected_native_metadata in member_names:
        raise ValueError("failed archive unexpectedly contains native metadata")

    size = int(row["size_bytes"])
    starts = sorted({
        0,
        size // 4,
        size // 2,
        (3 * size) // 4,
        max(0, size - RANGE_LENGTH),
    })
    range_checks = []
    with path.open("rb") as local:
        for start in starts:
            end = min(size - 1, start + RANGE_LENGTH - 1)
            local.seek(start)
            local_value = local.read(end - start + 1)
            remote_value = remote_range(row["url"], start, end)
            local_hash = hashlib.sha256(local_value).hexdigest()
            remote_hash = hashlib.sha256(remote_value).hexdigest()
            if local_hash != remote_hash:
                raise ValueError("remote archive range differs from acquired object")
            range_checks.append({
                "start_byte": start,
                "end_byte_inclusive": end,
                "sha256": local_hash,
                "remote_matches_local": True,
            })

    report = {
        "schema": "aleph.outer_library.idr0168_decode_refusal.v1",
        "protocol_sha256": sha256(args.protocol),
        "acquisition_sha256": sha256(args.acquisition),
        "failed_archive": {
            "local_name": row["local_name"],
            "source_name": row["source_name"],
            "gene_symbol": row["gene_symbol"],
            "cell_line": row["cell_line"],
            "url": row["url"],
            "size_bytes": row["size_bytes"],
            "sha256": row["sha256"],
        },
        "source_integrity": {
            "gzip_integrity_passed": False,
            "gzip_error": gzip_error,
            "tar_stream_iteration_completed_without_exception": tar_error is None,
            "tar_error": tar_error,
            "native_Zarr_metadata_present": False,
            "expected_native_Zarr_metadata_member": expected_native_metadata,
            "present_pyramid_array_metadata_members": sorted(member_names),
            "remote_range_checks": range_checks,
            "all_sampled_remote_ranges_match_acquired_object": True,
        },
        "decoded_fields_before_refusal": 7,
        "tensor_written": False,
        "model_predictions_made": 0,
        "replacement_field_used": False,
        "single_use_confirmation_consumed": True,
        "may_be_reused_as_pristine": False,
        "status": "failed_source_integrity_and_permanently_consumed",
        "external_provider_validation_complete": False,
        "nuclear_localization_external_evidence_eligible": False,
        "production_eligible": False,
        "Aleph_parameter_authority": "none",
        "may_emit_numeric_sweep_range": False,
        "may_remove_sweep_axis": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": report["status"],
        "archive": row["local_name"],
        "remote_ranges_matched": len(range_checks),
        "predictions": 0,
    }, sort_keys=True))


if __name__ == "__main__":
    main()

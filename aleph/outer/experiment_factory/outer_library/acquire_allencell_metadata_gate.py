#!/usr/bin/env python3
"""Acquire and validate Allen Cell package metadata without opening image pixels."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request
from typing import Any

from train_hpa_raw_if_cnn import sha256


def normalized(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.casefold())


def download(url: str, path: Path, expected_size: int | None = None) -> None:
    if path.exists() and (expected_size is None or path.stat().st_size == expected_size):
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(path.suffix + ".part")
    error: Exception | None = None
    for attempt in range(5):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "Project-Aleph-research/1.0"}
            )
            with urllib.request.urlopen(request, timeout=300) as response, part.open(
                "wb"
            ) as handle:
                while chunk := response.read(4 * 1024 * 1024):
                    handle.write(chunk)
            if expected_size is not None and part.stat().st_size != expected_size:
                raise ValueError(
                    f"download size {part.stat().st_size} != {expected_size}"
                )
            part.replace(path)
            return
        except Exception as exc:
            error = exc
            part.unlink(missing_ok=True)
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"failed metadata download {url}: {error}")


def physical_url(value: str) -> str:
    prefix = "s3://allencell/"
    if not value.startswith(prefix):
        raise ValueError(f"unexpected Quilt physical key: {value}")
    remainder = value[len(prefix):]
    path, separator, query = remainder.partition("?")
    url = "https://allencell.s3.amazonaws.com/" + urllib.parse.quote(path, safe="/")
    return url + ("?" + query if separator else "")


def header(path: Path) -> list[str]:
    suffix = path.name.casefold()
    if suffix.endswith(".csv"):
        with path.open(encoding="utf-8-sig", newline="") as handle:
            return next(csv.reader(handle))
    if suffix.endswith(".parquet"):
        try:
            import pyarrow.parquet as pq
        except ImportError as exc:
            raise RuntimeError("pyarrow is required to inspect a Parquet candidate") from exc
        return list(pq.ParquetFile(path).schema.names)
    raise ValueError(f"unsupported metadata candidate: {path.name}")


def resolve_headers(
    headers: list[str], aliases: dict[str, list[str]],
) -> tuple[dict[str, str], dict[str, list[str]]]:
    normalized_headers = {name: normalized(name) for name in headers}
    resolved, ambiguous = {}, {}
    for semantic, values in aliases.items():
        accepted = {normalized(value) for value in values}
        matches = [
            original for original, canonical in normalized_headers.items()
            if canonical in accepted
        ]
        if len(matches) == 1:
            resolved[semantic] = matches[0]
        elif len(matches) > 1:
            ambiguous[semantic] = matches
    return resolved, ambiguous


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--metadata-amendment", type=Path, required=True)
    parser.add_argument("--prior-index-receipt", type=Path, required=True)
    parser.add_argument("--prior-schema-refusal", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    amendment = json.loads(args.metadata_amendment.read_text(encoding="utf-8"))
    if amendment["base_protocol_sha256"] != sha256(args.protocol):
        raise ValueError("Allen metadata amendment does not pin supplied base protocol")
    if not amendment["frozen_before_allencell_metadata_csv_header_or_pixel_access"]:
        raise ValueError("Allen header aliases were not frozen before metadata access")
    frozen = protocol["frozen_inputs"]
    if sha256(args.prior_index_receipt) != frozen[
        "allencell_first_index_receipt_sha256"
    ]:
        raise ValueError("Allen first index receipt changed after protocol freeze")
    refusal = json.loads(args.prior_schema_refusal.read_text(encoding="utf-8"))
    if refusal["image_objects_downloaded"] or refusal["model_predictions_made"]:
        raise ValueError("prior Allen refusal unexpectedly touched confirmation pixels")
    provider = protocol["provider_roles"]["AllenCell"]
    pointer_hash = provider["frozen_pointer_manifest_sha256"]
    pointer_url = (
        "https://allencell.s3.amazonaws.com/.quilt/named_packages/"
        "aics/hipsc_single_cell_image_dataset/1657796402"
    )
    pointer_path = args.cache_dir / "named_package_pointer.txt"
    download(pointer_url, pointer_path, 64)
    pointer_value = pointer_path.read_text(encoding="ascii")
    if pointer_value != pointer_hash:
        raise ValueError("Allen named package pointer changed after protocol freeze")
    manifest_url = f"https://allencell.s3.amazonaws.com/.quilt/packages/{pointer_hash}"
    manifest_path = args.cache_dir / f"{pointer_hash}.jsonl"
    download(manifest_url, manifest_path)
    records: list[dict[str, Any]] = []
    package_meta = None
    with manifest_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            value = json.loads(line)
            if line_number == 1 and "version" in value:
                package_meta = value
            elif "logical_key" in value:
                records.append(value)
            else:
                raise ValueError(f"unexpected Quilt manifest line {line_number}")
    candidates = sorted(
        (
            item for item in records
            if item["logical_key"].casefold().endswith((".csv", ".parquet"))
        ),
        key=lambda item: item["logical_key"],
    )
    aliases = amendment["normalized_header_aliases"]
    candidate_reports, selected = [], None
    for index, record in enumerate(candidates):
        physical = record["physical_keys"]
        if len(physical) != 1:
            candidate_reports.append({
                "logical_key": record["logical_key"],
                "status": "ambiguous_physical_keys",
                "physical_key_count": len(physical),
            })
            continue
        local = args.cache_dir / "metadata_candidates" / (
            f"{index:04d}_" + Path(record["logical_key"]).name
        )
        download(physical_url(physical[0]), local, int(record["size"]))
        expected_hash = record["hash"]
        if expected_hash["type"] != "SHA256" or sha256(local) != expected_hash["value"]:
            raise ValueError(f"Allen metadata hash mismatch: {record['logical_key']}")
        try:
            headers = header(local)
            resolved, ambiguous = resolve_headers(headers, aliases)
            missing = sorted(set(aliases) - set(resolved) - set(ambiguous))
            status = "resolved" if not missing and not ambiguous else "schema_mismatch"
            candidate_report = {
                "logical_key": record["logical_key"],
                "local_name": str(local.relative_to(args.cache_dir)),
                "size_bytes": record["size"],
                "sha256": expected_hash["value"],
                "headers": headers,
                "resolved_columns": resolved,
                "ambiguous_columns": ambiguous,
                "missing_semantic_fields": missing,
                "status": status,
            }
        except Exception as exc:
            candidate_report = {
                "logical_key": record["logical_key"],
                "local_name": str(local.relative_to(args.cache_dir)),
                "size_bytes": record["size"],
                "sha256": expected_hash["value"],
                "status": "inspection_error",
                "error": f"{type(exc).__name__}: {exc}",
            }
        candidate_reports.append(candidate_report)
        if candidate_report["status"] == "resolved":
            selected = candidate_report
            break
    report = {
        "schema": "aleph.outer_library.allencell_multiprovider_metadata_gate.v2",
        "protocol_sha256": sha256(args.protocol),
        "metadata_amendment_sha256": sha256(args.metadata_amendment),
        "prior_index_receipt_sha256": sha256(args.prior_index_receipt),
        "prior_schema_refusal_sha256": sha256(args.prior_schema_refusal),
        "pointer_url": pointer_url,
        "pointer_sha256": sha256(pointer_path),
        "pointer_value": pointer_value,
        "package_manifest_url": manifest_url,
        "package_manifest_sha256": sha256(manifest_path),
        "package_manifest_bytes": manifest_path.stat().st_size,
        "package_meta": package_meta,
        "package_member_count": len(records),
        "metadata_candidate_count": len(candidates),
        "candidate_reports_until_selection": candidate_reports,
        "selected_metadata": selected,
        "metadata_gate_passed": selected is not None,
        "status": (
            "metadata_schema_resolved_confirmation_pixels_sealed"
            if selected is not None else "refused_before_Allen_image_access"
        ),
        "Allen_image_objects_downloaded": 0,
        "Allen_pixels_accessed": 0,
        "model_predictions_made": 0,
        "metadata_may_change_ontology_or_model": False,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "manifest_members": len(records), "candidates": len(candidates),
        "metadata_gate_passed": report["metadata_gate_passed"],
        "selected": selected["logical_key"] if selected else None,
        "Allen_pixels": 0,
    }, sort_keys=True))


if __name__ == "__main__":
    main()

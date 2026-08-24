#!/usr/bin/env python3
"""Resolve and acquire the single-use Light My Cells confirmation partition."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import time
import urllib.request
from typing import Any

from download_lightmycells_pairs import acquire
from resolve_lightmycells_pairs import fetch_header, sha256


def head(file: dict[str, Any]) -> dict[str, Any]:
    error: Exception | None = None
    for attempt in range(5):
        try:
            request = urllib.request.Request(
                file["url"], method="HEAD",
                headers={"User-Agent": "Project-Aleph-research/1.0"},
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                return {
                    "path": file["path"], "url": file["url"],
                    "listed_size_bytes": int(file["size_bytes"]),
                    "remote_content_length_bytes": int(
                        response.headers["Content-Length"]
                    ),
                    "remote_etag": response.headers.get("ETag"),
                }
        except Exception as exc:
            error = exc
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"failed source-integrity HEAD {file['path']}: {error}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-protocol", type=Path, required=True)
    parser.add_argument("--resolution-protocol", type=Path, required=True)
    parser.add_argument("--model-protocol", type=Path, required=True)
    parser.add_argument("--model-report", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--metadata-manifest", type=Path, required=True)
    parser.add_argument("--training-pair-report", type=Path, required=True)
    parser.add_argument("--header-cache", type=Path, required=True)
    parser.add_argument("--image-dir", type=Path, required=True)
    parser.add_argument("--pair-manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=24)
    args = parser.parse_args()
    base = json.loads(args.base_protocol.read_text(encoding="utf-8"))
    resolution = json.loads(args.resolution_protocol.read_text(encoding="utf-8"))
    model_protocol = json.loads(args.model_protocol.read_text(encoding="utf-8"))
    model_report = json.loads(args.model_report.read_text(encoding="utf-8"))
    if not model_report["confirmation_open_condition_passed"]:
        raise ValueError("frozen model did not earn confirmation access")
    if model_report["confirmation_pairs_seen"] != 0:
        raise ValueError("model report claims prior confirmation access")
    if model_report["protocol_sha256"] != sha256(args.model_protocol):
        raise ValueError("model report does not pin supplied model protocol")
    if model_report["checkpoint_sha256"] != sha256(args.checkpoint):
        raise ValueError("model report does not pin supplied checkpoint")
    if resolution["base_protocol"]["sha256"] != sha256(args.base_protocol):
        raise ValueError("resolution protocol does not pin base protocol")

    rows = [
        json.loads(line) for line in args.metadata_manifest.read_text(encoding="utf-8").splitlines()
    ]
    rows = [row for row in rows if row["role"] == "single_use_confirmation"]
    if {row["study_id"] for row in rows} != set(
        base["split_contract"]["single_use_confirmation"]
    ):
        raise ValueError("confirmation studies differ from frozen protocol")
    target_names = set(base["author_label_contract"]["target_channels"])
    input_names = set(base["author_label_contract"]["input_modalities"])
    header_files: dict[str, dict[str, Any]] = {}
    for row in rows:
        for file in row["files"]:
            if file["channel"] in target_names:
                header_files[file["url"]] = file
        for modality in input_names:
            inputs = sorted(
                (file for file in row["files"] if file["channel"] == modality),
                key=lambda file: file["z_index"],
            )
            if inputs:
                header_files[inputs[0]["url"]] = inputs[0]
                header_files[inputs[-1]["url"]] = inputs[-1]

    metadata, receipts, header_errors = {}, {}, {}
    def parse(item: tuple[str, dict[str, Any]]) -> tuple[str, dict | None, dict | None, str | None]:
        url, file = item
        try:
            parsed, receipt = fetch_header(url, args.header_cache, require_position=False)
            if parsed["image_name"] != Path(file["path"]).name:
                raise ValueError("OME image name differs from FileList path")
            return url, parsed, receipt, None
        except Exception as exc:
            return url, None, None, f"{type(exc).__name__}: {exc}"
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(parse, item) for item in header_files.items()]
        for future in as_completed(futures):
            url, parsed, receipt, error = future.result()
            if error is None:
                metadata[url], receipts[url] = parsed, receipt
            else:
                header_errors[url] = error

    pairs, refusals = [], []
    for row in rows:
        targets = [file for file in row["files"] if file["channel"] in target_names]
        for modality in sorted(input_names):
            inputs = sorted(
                (file for file in row["files"] if file["channel"] == modality),
                key=lambda file: file["z_index"],
            )
            if not inputs:
                continue
            anchors = [inputs[0]] if len(inputs) == 1 else [inputs[0], inputs[-1]]
            if any(file["url"] in header_errors for file in anchors):
                refusals.append({
                    "study_id": row["study_id"],
                    "acquisition_set_id": row["acquisition_set_id"],
                    "modality": modality,
                    "reason": "input anchor OME header refusal",
                })
                continue
            anchor_meta = [metadata[file["url"]] for file in anchors]
            if len(inputs) == 1:
                predicted_positions = {inputs[0]["z_index"]: anchor_meta[0]["position_z"]}
            elif all(item["position_z"] is not None for item in anchor_meta):
                first_z, last_z = inputs[0]["z_index"], inputs[-1]["z_index"]
                step = (
                    anchor_meta[-1]["position_z"] - anchor_meta[0]["position_z"]
                ) / (last_z - first_z)
                predicted_positions = {
                    file["z_index"]: anchor_meta[0]["position_z"]
                    + step * (file["z_index"] - first_z)
                    for file in inputs
                }
            else:
                predicted_positions = {file["z_index"]: None for file in inputs}
            for target in targets:
                if target["url"] in header_errors:
                    refusals.append({
                        "study_id": row["study_id"],
                        "acquisition_set_id": row["acquisition_set_id"],
                        "modality": modality, "target": target["channel"],
                        "reason": header_errors[target["url"]],
                    })
                    continue
                target_meta = metadata[target["url"]]
                same_unit = (
                    target_meta["position_z"] is not None
                    and all(value is not None for value in predicted_positions.values())
                    and all(
                        item["position_z_unit"] == target_meta["position_z_unit"]
                        for item in anchor_meta
                    )
                )
                if len(inputs) == 1:
                    selected, method = inputs[0], "singleton_input"
                elif same_unit:
                    selected = min(
                        inputs,
                        key=lambda file: (
                            abs(predicted_positions[file["z_index"]] - target_meta["position_z"]),
                            file["z_index"],
                        ),
                    )
                    method = "nearest_PositionZ"
                else:
                    selected = inputs[(len(inputs) - 1) // 2]
                    method = "lower_median_z_fallback"
                selected_position = predicted_positions[selected["z_index"]]
                mismatch = (
                    abs(selected_position - target_meta["position_z"])
                    if selected_position is not None and target_meta["position_z"] is not None
                    else None
                )
                pairs.append({
                    "schema": "aleph.outer_library.lightmycells_confirmation_pair.v1",
                    "accession": row["accession"], "study_id": row["study_id"],
                    "study_name": row["study_name"], "role": row["role"],
                    "acquisition_set_id": row["acquisition_set_id"],
                    "input_modality": modality, "target_channel": target["channel"],
                    "input": selected, "target": target,
                    "alignment_method": method,
                    "input_position_z": selected_position,
                    "target_position_z": target_meta["position_z"],
                    "position_z_unit": target_meta["position_z_unit"],
                    "absolute_z_mismatch": mismatch,
                    "experimenter_group": target_meta["experimenter_group"],
                    "experimenter_group_id": target_meta["experimenter_group_id"],
                    "checkpoint_sha256": sha256(args.checkpoint),
                    "aleph_authority": "none",
                })

    unique = {file["url"]: file for pair in pairs for file in (pair["input"], pair["target"])}
    head_records = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(head, file) for file in unique.values()]
        for future in as_completed(futures):
            head_records.append(future.result())
    inconsistent_paths = {
        item["path"] for item in head_records
        if item["listed_size_bytes"] != item["remote_content_length_bytes"]
    }
    integrity_refused = [
        pair for pair in pairs
        if pair["input"]["path"] in inconsistent_paths
        or pair["target"]["path"] in inconsistent_paths
    ]
    pairs = [pair for pair in pairs if pair not in integrity_refused]
    unique = {file["url"]: file for pair in pairs for file in (pair["input"], pair["target"])}
    args.pair_manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.pair_manifest.open("w", encoding="utf-8") as handle:
        for pair in pairs:
            handle.write(json.dumps(pair, sort_keys=True) + "\n")

    files = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(acquire, file, args.image_dir) for file in unique.values()]
        for future in as_completed(futures):
            files.append(future.result())
    files.sort(key=lambda item: item["path"])
    training_report = json.loads(args.training_pair_report.read_text(encoding="utf-8"))
    training_groups = set(training_report["experimenter_groups"])
    confirmation_groups = {
        pair["experimenter_group"] for pair in pairs if pair["experimenter_group"]
    }
    report = {
        "schema": "aleph.outer_library.lightmycells_confirmation_acquisition.v1",
        "base_protocol_sha256": sha256(args.base_protocol),
        "resolution_protocol_sha256": sha256(args.resolution_protocol),
        "model_protocol_sha256": sha256(args.model_protocol),
        "model_report_sha256": sha256(args.model_report),
        "checkpoint_sha256": sha256(args.checkpoint),
        "metadata_manifest_sha256": sha256(args.metadata_manifest),
        "pair_manifest_sha256": sha256(args.pair_manifest),
        "study_count": len({pair["study_id"] for pair in pairs}),
        "pair_count": len(pairs),
        "unique_payload_count": len(unique),
        "unique_payload_bytes": sum(file["size_bytes"] for file in unique.values()),
        "alignment_method_counts": dict(sorted(Counter(
            pair["alignment_method"] for pair in pairs
        ).items())),
        "target_pair_counts": dict(sorted(Counter(
            pair["target_channel"] for pair in pairs
        ).items())),
        "header_request_count": len(receipts),
        "header_error_count": len(header_errors),
        "header_refusals": refusals,
        "integrity_inconsistent_payloads": sorted(
            (
                item
                for item in head_records
                if item["path"] in inconsistent_paths
            ),
            key=lambda item: item["path"],
        ),
        "integrity_refused_pair_count": len(integrity_refused),
        "training_experimenter_groups": sorted(training_groups),
        "confirmation_experimenter_groups": sorted(confirmation_groups),
        "experimenter_group_overlap": sorted(training_groups & confirmation_groups),
        "novel_confirmation_experimenter_groups": sorted(
            confirmation_groups - training_groups
        ),
        "study_holdout": True,
        "all_confirmation_labs_independent": not bool(training_groups & confirmation_groups),
        "model_predictions_made": 0,
        "confirmation_consumed": True,
        "may_be_reused_as_pristine": False,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt = {
        "schema": "aleph.outer_library.lightmycells_confirmation_payload_receipt.v1",
        "report_sha256": sha256(args.report),
        "pair_manifest_sha256": sha256(args.pair_manifest),
        "file_count": len(files), "total_bytes": sum(item["size_bytes"] for item in files),
        "files": files, "checkpoint_sha256": sha256(args.checkpoint),
        "confirmation_consumed": True, "aleph_authority": "none",
    }
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "studies": report["study_count"], "pairs": len(pairs),
        "payloads": len(files), "bytes": receipt["total_bytes"],
        "targets": report["target_pair_counts"],
        "groups": report["confirmation_experimenter_groups"],
        "novel_groups": report["novel_confirmation_experimenter_groups"],
        "group_overlap": report["experimenter_group_overlap"],
        "integrity_refused_pairs": len(integrity_refused),
    }, sort_keys=True))


if __name__ == "__main__":
    main()

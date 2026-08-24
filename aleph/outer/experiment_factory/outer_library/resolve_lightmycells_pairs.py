#!/usr/bin/env python3
"""Resolve target-aligned S-BIAD1047 pairs from non-confirmation TIFF headers."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import struct
import time
import urllib.request
import xml.etree.ElementTree as ET
from typing import Any


class NeedMoreHeaderData(ValueError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def image_description(payload: bytes) -> str:
    if len(payload) < 8:
        raise NeedMoreHeaderData("truncated TIFF header")
    if payload[:2] == b"II":
        endian = "<"
    elif payload[:2] == b"MM":
        endian = ">"
    else:
        raise ValueError("not a TIFF byte order marker")
    magic = struct.unpack_from(endian + "H", payload, 2)[0]
    if magic == 42:
        ifd_offset = struct.unpack_from(endian + "I", payload, 4)[0]
        count_format, count_bytes, entry_bytes = "H", 2, 12
        entry_format, inline_bytes, value_format = "HHI", 4, "I"
    elif magic == 43:
        if len(payload) < 16:
            raise NeedMoreHeaderData("truncated BigTIFF header")
        offset_size, reserved, ifd_offset = struct.unpack_from(endian + "HHQ", payload, 4)
        if offset_size != 8 or reserved != 0:
            raise ValueError("unsupported BigTIFF offset structure")
        count_format, count_bytes, entry_bytes = "Q", 8, 20
        entry_format, inline_bytes, value_format = "HHQ", 8, "Q"
    else:
        raise ValueError(f"unsupported TIFF magic {magic}")
    if ifd_offset + count_bytes > len(payload):
        raise NeedMoreHeaderData("first IFD lies outside fetched header")
    count = struct.unpack_from(endian + count_format, payload, ifd_offset)[0]
    entry_start = ifd_offset + count_bytes
    if entry_start + entry_bytes * count > len(payload):
        raise NeedMoreHeaderData("first IFD entries are truncated")
    for index in range(count):
        offset = entry_start + entry_bytes * index
        tag, field_type, length = struct.unpack_from(endian + entry_format, payload, offset)
        if tag != 270:
            continue
        if field_type != 2:
            raise ValueError("TIFF ImageDescription is not ASCII")
        value_field = offset + entry_bytes - inline_bytes
        if length <= inline_bytes:
            raw = payload[value_field:value_field + length]
        else:
            value_offset = struct.unpack_from(
                endian + value_format, payload, value_field
            )[0]
            if value_offset + length > len(payload):
                raise NeedMoreHeaderData("ImageDescription lies outside fetched header")
            raw = payload[value_offset:value_offset + length]
        return raw.rstrip(b"\x00").decode("utf-8")
    raise ValueError("TIFF has no ImageDescription tag")


def parse_position(payload: bytes, *, require_position: bool = True) -> dict[str, Any]:
    description = image_description(payload)
    root = ET.fromstring(description)
    image = root.find(".//{*}Image")
    pixels = root.find(".//{*}Pixels")
    plane = root.find(".//{*}Plane")
    group = root.find(".//{*}ExperimenterGroup")
    if image is None or pixels is None or plane is None:
        raise ValueError("OME metadata lacks Image, Pixels, or Plane")
    if require_position and "PositionZ" not in plane.attrib:
        raise ValueError("OME Plane lacks PositionZ")
    return {
        "image_name": image.attrib.get("Name"),
        "position_z": (
            float(plane.attrib["PositionZ"]) if "PositionZ" in plane.attrib else None
        ),
        "position_z_unit": plane.attrib.get("PositionZUnit"),
        "physical_size_x": (
            float(pixels.attrib["PhysicalSizeX"])
            if "PhysicalSizeX" in pixels.attrib else None
        ),
        "physical_size_x_unit": pixels.attrib.get("PhysicalSizeXUnit"),
        "size_x": int(pixels.attrib["SizeX"]),
        "size_y": int(pixels.attrib["SizeY"]),
        "pixel_type": pixels.attrib["Type"],
        "experimenter_group": group.attrib.get("Name") if group is not None else None,
        "experimenter_group_id": group.attrib.get("ID") if group is not None else None,
    }


def fetch_header(
    url: str, cache_dir: Path, *, require_position: bool = True
) -> tuple[dict[str, Any], dict[str, Any]]:
    key = hashlib.sha256(url.encode()).hexdigest()
    path = cache_dir / f"{key}.header"
    sizes = (16_384, 65_536, 262_144)
    error: Exception | None = None
    for size in sizes:
        if not path.exists() or path.stat().st_size < size:
            for attempt in range(5):
                try:
                    request = urllib.request.Request(
                        url,
                        headers={
                            "User-Agent": "Project-Aleph-research/1.0",
                            "Range": f"bytes=0-{size - 1}",
                        },
                    )
                    with urllib.request.urlopen(request, timeout=120) as response:
                        payload = response.read()
                    if len(payload) > size:
                        raise ValueError("server ignored bounded Range request")
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(payload)
                    break
                except Exception as exc:
                    error = exc
                    time.sleep(1.5 * (attempt + 1))
            else:
                raise RuntimeError(f"failed TIFF header request {url}: {error}")
        payload = path.read_bytes()
        try:
            metadata = parse_position(payload, require_position=require_position)
            return metadata, {
                "url": url,
                "header_sha256": hashlib.sha256(payload).hexdigest(),
                "header_bytes": len(payload),
            }
        except NeedMoreHeaderData as exc:
            error = exc
            continue
    raise RuntimeError(f"OME metadata exceeds bounded header budget for {url}: {error}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--resolution-protocol", type=Path)
    parser.add_argument("--metadata-manifest", type=Path, required=True)
    parser.add_argument("--header-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=24)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    resolution_protocol = (
        json.loads(args.resolution_protocol.read_text(encoding="utf-8"))
        if args.resolution_protocol else None
    )
    fallback_enabled = resolution_protocol is not None
    if resolution_protocol is not None:
        base = resolution_protocol["base_protocol"]
        if Path(base["path"]).name != args.protocol.name or base["sha256"] != sha256(args.protocol):
            raise ValueError("pair-resolution wrapper does not pin the supplied base protocol")
        diagnostic = resolution_protocol["source_informed_diagnostic"]
        diagnostic_path = args.resolution_protocol.parent / diagnostic["path"]
        if sha256(diagnostic_path) != diagnostic["sha256"]:
            raise ValueError("pair-resolution wrapper diagnostic hash changed")
    rows = [
        json.loads(line) for line in args.metadata_manifest.read_text(encoding="utf-8").splitlines()
    ]
    if any(
        row["role"] == "single_use_confirmation" and row["pixel_content_accessed"]
        for row in rows
    ):
        raise ValueError("upstream manifest claims confirmation pixel access")

    pairs = []
    refusals = []
    header_receipts: dict[str, dict[str, Any]] = {}
    header_metadata: dict[str, dict[str, Any]] = {}
    header_errors: dict[str, str] = {}
    target_names = set(protocol["author_label_contract"]["target_channels"])
    input_names = set(protocol["author_label_contract"]["input_modalities"])

    header_files: dict[str, dict[str, Any]] = {}
    for row in rows:
        if row["role"] == "single_use_confirmation":
            continue
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

    def acquire(
        item: tuple[str, dict[str, Any]]
    ) -> tuple[str, dict[str, Any] | None, dict[str, Any] | None, str | None]:
        url, file = item
        try:
            parsed, receipt = fetch_header(
                url, args.header_cache, require_position=not fallback_enabled
            )
            if parsed["image_name"] != Path(file["path"]).name:
                raise ValueError(f"OME Image name mismatch for {file['path']}")
            return url, parsed, receipt, None
        except Exception as exc:
            return url, None, None, f"{type(exc).__name__}: {exc}"

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(acquire, item) for item in header_files.items()]
        for future in as_completed(futures):
            url, parsed, receipt, error = future.result()
            if error is not None:
                header_errors[url] = error
            else:
                assert parsed is not None and receipt is not None
                header_metadata[url] = parsed
                header_receipts[url] = receipt

    def metadata(file: dict[str, Any]) -> dict[str, Any]:
        url = file["url"]
        if url in header_errors:
            raise ValueError(header_errors[url])
        if url not in header_metadata:
            raise ValueError(f"unfrozen header request attempted: {file['path']}")
        return header_metadata[url]

    for row in rows:
        if row["role"] == "single_use_confirmation":
            continue
        targets = [file for file in row["files"] if file["channel"] in target_names]
        for modality in sorted(input_names):
            inputs = sorted(
                (file for file in row["files"] if file["channel"] == modality),
                key=lambda file: file["z_index"],
            )
            if not inputs:
                continue
            anchors = [inputs[0]] if len(inputs) == 1 else [inputs[0], inputs[-1]]
            try:
                anchor_metadata = [metadata(file) for file in anchors]
            except Exception as exc:
                refusals.append({
                    "study_id": row["study_id"],
                    "acquisition_set_id": row["acquisition_set_id"],
                    "modality": modality,
                    "reason": f"input anchor metadata refusal: {exc}",
                })
                continue
            if len(inputs) == 1:
                predicted_positions = {
                    inputs[0]["z_index"]: anchor_metadata[0]["position_z"]
                }
                default_alignment_method = "singleton_input"
            else:
                first_z, last_z = inputs[0]["z_index"], inputs[-1]["z_index"]
                first_pos, last_pos = (
                    anchor_metadata[0]["position_z"], anchor_metadata[1]["position_z"]
                )
                if last_z == first_z:
                    raise ValueError("distinct input anchors share z index")
                if first_pos is not None and last_pos is not None:
                    step = (last_pos - first_pos) / (last_z - first_z)
                    predicted_positions = {
                        file["z_index"]: first_pos + step * (file["z_index"] - first_z)
                        for file in inputs
                    }
                else:
                    predicted_positions = {file["z_index"]: None for file in inputs}
                default_alignment_method = "nearest_PositionZ"
            for target in targets:
                try:
                    target_meta = metadata(target)
                except Exception as exc:
                    refusals.append({
                        "study_id": row["study_id"],
                        "acquisition_set_id": row["acquisition_set_id"],
                        "modality": modality,
                        "target": target["channel"],
                        "reason": f"target metadata refusal: {exc}",
                    })
                    continue
                same_position_unit = (
                    target_meta["position_z"] is not None
                    and all(value is not None for value in predicted_positions.values())
                    and all(
                        item["position_z_unit"] == target_meta["position_z_unit"]
                        for item in anchor_metadata
                    )
                )
                if len(inputs) == 1:
                    selected = inputs[0]
                    alignment_method = "singleton_input"
                elif same_position_unit:
                    selected = min(
                        inputs,
                        key=lambda file: (
                            abs(
                                predicted_positions[file["z_index"]]
                                - target_meta["position_z"]
                            ),
                            file["z_index"],
                        ),
                    )
                    alignment_method = default_alignment_method
                elif fallback_enabled:
                    selected = inputs[(len(inputs) - 1) // 2]
                    alignment_method = "lower_median_z_fallback"
                else:
                    refusals.append({
                        "study_id": row["study_id"],
                        "acquisition_set_id": row["acquisition_set_id"],
                        "modality": modality,
                        "target": target["channel"],
                        "reason": "PositionZ unavailable or unit-incompatible",
                    })
                    continue
                input_meta = metadata(selected) if selected in anchors else {
                    **anchor_metadata[0],
                    "position_z": predicted_positions[selected["z_index"]],
                    "image_name": Path(selected["path"]).name,
                }
                mismatch = (
                    abs(input_meta["position_z"] - target_meta["position_z"])
                    if (
                        input_meta["position_z"] is not None
                        and target_meta["position_z"] is not None
                        and input_meta["position_z_unit"] == target_meta["position_z_unit"]
                    ) else None
                )
                pairs.append({
                    "schema": "aleph.outer_library.lightmycells_spatial_pair.v1",
                    "accession": row["accession"],
                    "study_id": row["study_id"],
                    "study_name": row["study_name"],
                    "role": row["role"],
                    "acquisition_set_id": row["acquisition_set_id"],
                    "input_modality": modality,
                    "target_channel": target["channel"],
                    "input": selected,
                    "target": target,
                    "input_position_z": input_meta["position_z"],
                    "target_position_z": target_meta["position_z"],
                    "position_z_unit": target_meta["position_z_unit"],
                    "absolute_z_mismatch": mismatch,
                    "alignment_method": alignment_method,
                    "exact_position_alignment": alignment_method == "nearest_PositionZ",
                    "input_pixel_geometry": {
                        key: input_meta[key] for key in (
                            "size_x", "size_y", "pixel_type",
                            "physical_size_x", "physical_size_x_unit",
                        )
                    },
                    "target_pixel_geometry": {
                        key: target_meta[key] for key in (
                            "size_x", "size_y", "pixel_type",
                            "physical_size_x", "physical_size_x_unit",
                        )
                    },
                    "experimenter_group": target_meta["experimenter_group"],
                    "experimenter_group_id": target_meta["experimenter_group_id"],
                    "confirmation_pixels_accessed": False,
                    "aleph_authority": "none",
                })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for pair in pairs:
            handle.write(json.dumps(pair, sort_keys=True) + "\n")
    unique_payloads = {
        file["url"]: file
        for pair in pairs for file in (pair["input"], pair["target"])
    }
    role_counts = Counter(pair["role"] for pair in pairs)
    target_counts = Counter(pair["target_channel"] for pair in pairs)
    modality_counts = Counter(pair["input_modality"] for pair in pairs)
    alignment_counts = Counter(pair["alignment_method"] for pair in pairs)
    experimenter_groups = sorted({
        pair["experimenter_group"] for pair in pairs if pair["experimenter_group"]
    })
    report = {
        "schema": "aleph.outer_library.lightmycells_spatial_pair_resolution.v1",
        "protocol_sha256": sha256(args.protocol),
        "metadata_manifest_sha256": sha256(args.metadata_manifest),
        "output_sha256": sha256(args.output),
        "pair_count": len(pairs),
        "unique_selected_payload_count": len(unique_payloads),
        "selected_payload_bytes": sum(file["size_bytes"] for file in unique_payloads.values()),
        "role_pair_counts": dict(sorted(role_counts.items())),
        "target_pair_counts": dict(sorted(target_counts.items())),
        "input_modality_pair_counts": dict(sorted(modality_counts.items())),
        "alignment_method_pair_counts": dict(sorted(alignment_counts.items())),
        "experimenter_group_count": len(experimenter_groups),
        "experimenter_groups": experimenter_groups,
        "header_request_count": len(header_receipts),
        "header_error_count": len(header_errors),
        "header_errors": [
            {"url": url, "reason": reason}
            for url, reason in sorted(header_errors.items())
        ],
        "header_bytes_fetched": sum(item["header_bytes"] for item in header_receipts.values()),
        "header_receipts": sorted(header_receipts.values(), key=lambda item: item["url"]),
        "refusal_count": len(refusals),
        "refusals": refusals,
        "confirmation_study_count_skipped_before_header_access": 6,
        "confirmation_header_requests": 0,
        "confirmation_pixels_accessed": 0,
        "resolution_fallback_enabled": fallback_enabled,
        "resolution_protocol_sha256": (
            sha256(args.resolution_protocol) if args.resolution_protocol else None
        ),
        "model_predictions_made": 0,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "pairs": len(pairs), "payloads": len(unique_payloads),
        "payload_bytes": report["selected_payload_bytes"],
        "headers": len(header_receipts),
        "header_bytes": report["header_bytes_fetched"],
        "refusals": len(refusals), "groups": len(experimenter_groups),
        "roles": report["role_pair_counts"], "targets": report["target_pair_counts"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

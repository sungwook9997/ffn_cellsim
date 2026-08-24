#!/usr/bin/env python3
"""Acquire and audit S-BIAD1047 metadata without opening image pixels."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_json(url: str, path: Path) -> Any:
    if not path.exists():
        error: Exception | None = None
        for attempt in range(5):
            try:
                request = urllib.request.Request(
                    url, headers={"User-Agent": "Project-Aleph-research/1.0"}
                )
                with urllib.request.urlopen(request, timeout=180) as response:
                    payload = response.read()
                json.loads(payload)
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(payload)
                break
            except Exception as exc:
                error = exc
                time.sleep(1.5 * (attempt + 1))
        else:
            raise RuntimeError(f"failed to acquire metadata {url}: {error}")
    return json.loads(path.read_text(encoding="utf-8"))


def attrs(section: dict[str, Any]) -> dict[str, list[str]]:
    values: dict[str, list[str]] = defaultdict(list)
    for item in section.get("attributes", []):
        values[item["name"]].append(item.get("value", ""))
    return dict(values)


def expected_roles(protocol: dict[str, Any]) -> dict[str, str]:
    result = {}
    split = protocol["split_contract"]
    for role in ("fit", "selection", "calibration", "single_use_confirmation"):
        for study in split[role]:
            if study in result:
                raise ValueError(f"study assigned twice: {study}")
            result[study] = role
    if set(result) != {f"Study_{index}" for index in range(1, 31)}:
        raise ValueError("protocol does not assign every study exactly once")
    salt = split["selection_salt"]
    ordered = sorted(
        result,
        key=lambda study: hashlib.sha256(f"{salt}|{study}".encode()).hexdigest(),
    )
    computed = {}
    for index, study in enumerate(ordered):
        computed[study] = (
            "single_use_confirmation" if index < 6 else
            "calibration" if index < 10 else
            "selection" if index < 14 else "fit"
        )
    if computed != result:
        raise ValueError("declared study assignments disagree with frozen hash algorithm")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--metadata-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    roles = expected_roles(protocol)
    source = protocol["source"]
    study_json = fetch_json(source["study_api"], args.metadata_dir / "S-BIAD1047.json")
    if study_json.get("accno") != "S-BIAD1047":
        raise ValueError("BioStudies response accession mismatch")
    root_attrs = attrs(study_json["section"])
    if root_attrs.get("License") != ["CC BY 4.0"]:
        raise ValueError(f"unexpected source license: {root_attrs.get('License')}")

    components = {}
    for section in study_json["section"].get("subsections", []):
        if section.get("type") != "Study Component":
            continue
        values = attrs(section)
        file_list = values["File List"][0]
        match = re.fullmatch(r"FileList/(Study_[0-9]+)\.json", file_list)
        if not match:
            raise ValueError(f"unexpected File List reference: {file_list}")
        components[match.group(1)] = {
            "name": values["Name"][0],
            "description": values.get("Description", [""])[0],
            "file_list": file_list,
        }
    if set(components) != set(roles):
        raise ValueError("BioStudies Study Components disagree with frozen Study_1..Study_30")

    pattern = re.compile(protocol["author_label_contract"]["acquisition_set_regex"])
    image_root = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/047/S-BIAD1047/Files/"
    studies = {}
    rows = []
    total_files = total_bytes = 0
    unmatched = []
    global_channel_counts: Counter[str] = Counter()
    for number in range(1, 31):
        study = f"Study_{number}"
        url = source["file_list_template"].format(study_number=number)
        path = args.metadata_dir / f"{study}.json"
        records = fetch_json(url, path)
        if not isinstance(records, list):
            raise ValueError(f"{study} FileList is not a list")
        grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
        study_channels: Counter[str] = Counter()
        for record in records:
            if record.get("type") != "file":
                continue
            total_files += 1
            size = int(record["size"])
            total_bytes += size
            name = Path(record["path"]).name
            match = pattern.fullmatch(name)
            if not match:
                unmatched.append({"study": study, "path": record["path"], "size": size})
                continue
            acquisition, channel, z_index = match.groups()
            study_channels[channel] += 1
            global_channel_counts[channel] += 1
            grouped[acquisition].append({
                "channel": channel,
                "z_index": int(z_index) if z_index is not None else None,
                "path": record["path"],
                "size_bytes": size,
                "url": image_root + urllib.parse.quote(record["path"], safe="/"),
            })
        studies[study] = {
            **components[study],
            "role": roles[study],
            "file_list_url": url,
            "file_list_size_bytes": path.stat().st_size,
            "file_list_sha256": sha256(path),
            "listed_file_count": len(records),
            "matched_acquisition_set_count": len(grouped),
            "channel_file_counts": dict(sorted(study_channels.items())),
        }
        for acquisition, files in sorted(grouped.items(), key=lambda item: int(item[0])):
            modalities = sorted({
                item["channel"] for item in files
                if item["channel"] in protocol["author_label_contract"]["input_modalities"]
            })
            targets = sorted({
                item["channel"] for item in files
                if item["channel"] in protocol["author_label_contract"]["target_channels"]
            })
            rows.append({
                "schema": "aleph.outer_library.lightmycells_acquisition_set.v1",
                "accession": "S-BIAD1047",
                "study_id": study,
                "study_name": components[study]["name"],
                "role": roles[study],
                "acquisition_set_id": int(acquisition),
                "input_modalities": modalities,
                "target_channels": targets,
                "files": sorted(files, key=lambda item: (item["channel"], item["z_index"] or -1)),
                "pixel_content_accessed": False,
                "aleph_authority": "none",
            })

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.manifest.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    role_summary = {}
    for role in ("fit", "selection", "calibration", "single_use_confirmation"):
        role_rows = [row for row in rows if row["role"] == role]
        role_summary[role] = {
            "study_count": sum(value == role for value in roles.values()),
            "acquisition_set_count": len(role_rows),
            "target_acquisition_set_counts": {
                target: sum(target in row["target_channels"] for row in role_rows)
                for target in protocol["author_label_contract"]["target_channels"]
            },
            "input_modality_acquisition_set_counts": {
                modality: sum(modality in row["input_modalities"] for row in role_rows)
                for modality in protocol["author_label_contract"]["input_modalities"]
            },
        }
    report = {
        "schema": "aleph.outer_library.lightmycells_metadata_audit.v1",
        "protocol_sha256": sha256(args.protocol),
        "source_study_json_sha256": sha256(args.metadata_dir / "S-BIAD1047.json"),
        "manifest_sha256": sha256(args.manifest),
        "license": root_attrs["License"][0],
        "study_count": len(studies),
        "listed_file_count": total_files,
        "listed_total_bytes": total_bytes,
        "matched_acquisition_set_count": len(rows),
        "matched_channel_file_counts": dict(sorted(global_channel_counts.items())),
        "unmatched_file_count": len(unmatched),
        "unmatched_files": unmatched,
        "studies": studies,
        "split_summary": role_summary,
        "published_count_comparison": {
            "paper_image_count": source["published_image_count"],
            "paper_acquisition_set_count": source["published_acquisition_set_count"],
            "listed_image_delta": total_files - source["published_image_count"],
            "matched_acquisition_set_delta": len(rows) - source["published_acquisition_set_count"],
        },
        "image_pixels_accessed": 0,
        "confirmation_pixels_accessed": 0,
        "model_predictions_made": 0,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "studies": len(studies), "files": total_files,
        "acquisition_sets": len(rows), "bytes": total_bytes,
        "channels": dict(sorted(global_channel_counts.items())),
        "unmatched": len(unmatched), "split_summary": role_summary,
    }, sort_keys=True))


if __name__ == "__main__":
    main()

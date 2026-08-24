#!/usr/bin/env python3
"""Acquire frozen OpenCell annotation files before external image selection."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.request


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path) -> dict[str, object]:
    request = urllib.request.Request(
        url, headers={"User-Agent": "Project-Aleph-research/1.0"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read()
        etag = response.headers.get("ETag", "").strip('"')
        content_type = response.headers.get("Content-Type")
        last_modified = response.headers.get("Last-Modified")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return {
        "url": url, "local_name": path.name, "size_bytes": len(payload),
        "sha256": sha256(path), "etag": etag,
        "content_type": content_type, "last_modified": last_modified,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if not protocol["frozen_before_opencell_annotation_or_image_content_access"]:
        raise ValueError("OpenCell protocol was not frozen before access")
    source = protocol["source_contract"]
    files = [
        download(source["annotation_csv_url"], args.output_dir / "opencell-localization-annotations.csv"),
        download(source["annotation_readme_url"], args.output_dir / "opencell-localization-annotations-readme.csv"),
    ]
    annotation = files[0]
    if annotation["size_bytes"] != source["annotation_csv_expected_size_bytes"]:
        raise ValueError("OpenCell annotation size changed from frozen HEAD metadata")
    if annotation["etag"] != source["annotation_csv_expected_etag"]:
        raise ValueError("OpenCell annotation ETag changed from frozen HEAD metadata")
    receipt = {
        "schema": "aleph.outer_library.opencell_if_annotation_acquisition.v1",
        "protocol_sha256": sha256(args.protocol),
        "protocol_commit": "2142d43+9dff081",
        "frozen_model_commit": "e755bd9",
        "frozen_model_checkpoint_sha256": "a0b9c18516d5997ec5cf64c97a7027ab6911fc21c7b46e5f7a6c0566dc15027e",
        "files": files,
        "annotation_content_used_during_acquisition": False,
        "aleph_authority": "none",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "files": len(files), "bytes": sum(item["size_bytes"] for item in files),
        "annotation_sha256": files[0]["sha256"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

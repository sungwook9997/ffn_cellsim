#!/usr/bin/env python3
"""Acquire the frozen Allen Cell feature index and Quilt version pointer."""

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


def fetch(url: str) -> tuple[bytes, dict[str, str | None]]:
    request = urllib.request.Request(url, headers={"User-Agent": "Project-Aleph-research/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read(), {
            "etag": response.headers.get("ETag", "").strip('"'),
            "last_modified": response.headers.get("Last-Modified"),
            "content_type": response.headers.get("Content-Type"),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if not protocol["frozen_before_allencell_feature_index_or_image_content_access"]:
        raise ValueError("Allen Cell protocol was not frozen before access")
    source = protocol["source_contract"]
    payload, headers = fetch(source["feature_index_url"])
    if len(payload) != source["feature_index_expected_size_bytes"]:
        raise ValueError("Allen feature-index size changed from frozen HEAD")
    if headers["etag"] != source["feature_index_expected_etag"]:
        raise ValueError("Allen feature-index ETag changed from frozen HEAD")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    index_path = args.output_dir / "aics_cfe_cell-features_v1.8.csv.zip"
    index_path.write_bytes(payload)
    pointer_url = (
        f"https://{source['quilt_bucket']}.s3.amazonaws.com/.quilt/named_packages/"
        f"{source['quilt_package']}/{source['quilt_version']}"
    )
    pointer, pointer_headers = fetch(pointer_url)
    pointer_hash = pointer.decode("ascii").strip()
    if pointer_hash != source["quilt_manifest_hash"]:
        raise ValueError("Allen Quilt version pointer differs from frozen manifest hash")
    receipt = {
        "schema": "aleph.outer_library.allencell_if_index_acquisition.v1",
        "protocol_sha256": sha256(args.protocol),
        "protocol_commit": "a5d52dc",
        "frozen_model_commit": source.get("frozen_model_commit", "e755bd9"),
        "feature_index": {
            "url": source["feature_index_url"], "local_name": index_path.name,
            "size_bytes": len(payload), "sha256": sha256(index_path), **headers,
        },
        "quilt_pointer": {
            "url": pointer_url, "manifest_hash": pointer_hash,
            "size_bytes": len(pointer), **pointer_headers,
        },
        "feature_or_image_content_used_during_acquisition": False,
        "aleph_authority": "none",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "bytes": len(payload), "sha256": receipt["feature_index"]["sha256"],
        "quilt_manifest_hash": pointer_hash,
    }, sort_keys=True))


if __name__ == "__main__":
    main()

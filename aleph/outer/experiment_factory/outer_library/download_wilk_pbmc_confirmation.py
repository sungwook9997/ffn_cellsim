#!/usr/bin/env python3
"""Acquire the frozen paper-linked Wilk GSE150728 CELLxGENE asset."""

from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

from download_pbmc_two_lab_recalibration import download, sha256


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if not protocol.get("frozen_before_Wilk_h5ad_expression_or_observation_access"):
        raise ValueError("Wilk confirmation protocol was not frozen before access")
    source = protocol["confirmation_source"]
    request = urllib.request.Request(
        source["cellxgene_collection_api"],
        headers={"User-Agent": "Mozilla/5.0 (compatible; Project-Aleph/1)"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        collection = json.load(response)
    dataset = next(
        (item for item in collection["datasets"]
         if item["dataset_id"] == source["cellxgene_dataset_id"]), None
    )
    if dataset is None:
        raise ValueError("frozen CELLxGENE dataset is absent from collection")
    asset = next(
        (item for item in dataset["assets"] if item["filetype"] == "H5AD"), None
    )
    if asset is None or asset["url"] != source["h5ad_asset_url"]:
        raise ValueError("frozen H5AD asset URL changed")
    if int(asset["filesize"]) != source["expected_asset_size_bytes"]:
        raise ValueError("frozen H5AD asset size changed")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    collection_path = args.output_dir / "cellxgene_collection.json"
    collection_path.write_text(
        json.dumps(collection, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    asset_path = args.output_dir / "GSE150728_Wilk_CELLxGENE.h5ad"
    download(asset["url"], asset_path)
    if asset_path.stat().st_size != source["expected_asset_size_bytes"]:
        raise ValueError("downloaded H5AD size differs from frozen metadata")
    receipt = {
        "schema": "aleph.outer_library.wilk_pbmc_confirmation_acquisition.v1",
        "protocol_sha256": sha256(args.protocol),
        "protocol_frozen_commit": "73ac9e1",
        "collection": {
            "collection_id": collection["collection_id"],
            "collection_version_id": collection["collection_version_id"],
            "dataset_id": dataset["dataset_id"],
            "dataset_title": dataset["title"],
            "asset_url": asset["url"],
            "declared_asset_size_bytes": asset["filesize"],
            "collection_json_sha256": sha256(collection_path),
        },
        "asset": {
            "path": asset_path.name, "size_bytes": asset_path.stat().st_size,
            "sha256": sha256(asset_path),
        },
        "role": "single_use_untouched_confirmation",
        "observation_content_used_during_acquisition": False,
        "aleph_authority": "none",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "dataset_id": dataset["dataset_id"], "bytes": asset_path.stat().st_size,
        "sha256": receipt["asset"]["sha256"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Acquire selected CC0 cellular-nematic PIV/TFM sources with fixed hashes."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request

from download_sciplex3_metadata import _download as _range_download


DATASET_DOI = "10.34810/DATA2772"
DATAVERSE_API = "https://dataverse.csuc.cat/api"
UPC_BITSTREAM = (
    "https://upcommons.upc.edu/server/api/core/bitstreams/"
    "0007e15a-e76a-4a6a-abd0-52b19de8b174/content"
)
FILES = {
    "00_Readme.txt": {
        "url": f"{DATAVERSE_API}/access/datafile/418560",
        "size": 11_950,
        "md5": "ac8f064f24963d645cd80fc560105da2",
        "role": "dataset_documentation",
    },
    "Fig1G_PIVdata.rar": {
        "url": f"{DATAVERSE_API}/access/datafile/402394",
        "size": 17_036_899,
        "md5": "c26da53148c2b12ee040d451367cb972",
        "role": "unconfined_nematic_PIV_tables",
    },
    "FigS09.rar": {
        "url": f"{DATAVERSE_API}/access/datafile/417938",
        "size": 488_983_295,
        "md5": "4dbe7f6bd1187c1c6d82f4cf608d2473",
        "role": "paired_pre_post_blebbistatin_TFM_MSM_tables",
    },
    "Guillamat_etal_2026.pdf": {
        "url": UPC_BITSTREAM,
        "size": 68_920_331,
        "md5": "4464fb5d1d9faaf1ccdb8bc2290eb547",
        "role": "open_access_article_and_supplementary_methods",
    },
}


def _digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _dataset_metadata() -> dict[str, object]:
    query = urllib.parse.urlencode({"persistentId": f"doi:{DATASET_DOI}"})
    request = urllib.request.Request(
        f"{DATAVERSE_API}/datasets/:persistentId/?{query}",
        headers={"User-Agent": "Project-Aleph-research/1.0"},
    )
    with urllib.request.urlopen(request, timeout=90) as response:
        payload = json.load(response)
    version = payload["data"]["latestVersion"]
    license_info = version["license"]
    remote_files = {
        item["label"]: {
            "id": item["dataFile"]["id"],
            "size": item["dataFile"]["filesize"],
            "md5": item["dataFile"]["md5"],
            "restricted": item["restricted"],
        }
        for item in version["files"]
    }
    expected = {
        name: {"size": spec["size"], "md5": spec["md5"]}
        for name, spec in FILES.items()
        if name != "Guillamat_etal_2026.pdf"
    }
    observed = {
        name: {"size": remote_files[name]["size"], "md5": remote_files[name]["md5"]}
        for name in expected
    }
    if observed != expected:
        raise ValueError(f"Dataverse selected-file contract changed: {observed}")
    if any(remote_files[name]["restricted"] for name in expected):
        raise ValueError("selected Dataverse source unexpectedly became restricted")
    if license_info.get("rightsIdentifier") != "CC0-1.0":
        raise ValueError(f"Dataverse license changed: {license_info}")
    return {
        "dataset_version": f"{version['versionNumber']}.{version['versionMinorNumber']}",
        "publication_date": version["publicationDate"],
        "last_update_time": version["lastUpdateTime"],
        "license": "cc0-1.0",
        "selected_remote_files": {
            name: remote_files[name] for name in expected
        },
    }


def acquire(
    output_dir: Path, *, segments: int = 256, workers: int = 48,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = _dataset_metadata()
    sources: dict[str, object] = {}
    for name, spec in FILES.items():
        path = output_dir / name
        file_segments = min(segments, max(1, int(spec["size"]) // (512 * 1024)))
        _range_download(
            str(spec["url"]), path, int(spec["size"]),
            segments=file_segments, workers=min(workers, file_segments),
        )
        md5 = _digest(path, "md5")
        if md5 != spec["md5"]:
            raise ValueError(f"MD5 mismatch for {name}: {md5}")
        sources[name] = {
            "size_bytes": path.stat().st_size,
            "md5": md5,
            "sha256": _digest(path, "sha256"),
            "role": spec["role"],
        }
    return {
        "schema": "aleph.outer_library.dataverse_data2772_acquisition.v1",
        "accession": DATASET_DOI,
        "provider": "CORA Research Data Repository Dataverse",
        "paper_doi": "10.1126/science.adz9174",
        "paper_repository": "UPCommons",
        "metadata": metadata,
        "sources": sources,
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--segments", type=int, default=256)
    parser.add_argument("--workers", type=int, default=48)
    args = parser.parse_args()
    receipt = acquire(args.output_dir, segments=args.segments, workers=args.workers)
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "source_count": len(receipt["sources"]),
        "total_bytes": sum(
            item["size_bytes"] for item in receipt["sources"].values()
        ),
    }, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Safely extract only aggregate Cell Painting profiles from Figshare."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile


EXPECTED_MEMBERS = {
    "features": 0,
    "features/single_cell_features_CellProfiler.parquet": 1_706_124_572,
    "features/single_cell_features_DeepProfiler.parquet": 530_778_321,
    "features/single_cell_features_DINO.parquet": 295_286_800,
    "features/aggregate_profiles_DINO_adjusted.parquet": 8_179_545,
    "features/aggregate_profiles_DP_aggregated.parquet": 16_579_528,
    "features/aggregate_profiles_DP_adjusted.parquet": 16_347_512,
}
AGGREGATES = tuple(name for name in EXPECTED_MEMBERS if "aggregate_profiles" in name)
FEATURE_MD5 = "680c12ebc351f260f56ef722ffc7c2aa"


def _safe(name: str) -> bool:
    path = PurePosixPath(name)
    return not path.is_absolute() and ".." not in path.parts and "\\" not in name


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract(archive: Path, receipt_path: Path, output_dir: Path) -> dict[str, object]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    verified = receipt.get("verified_local_files", {}).get("features.tar.gz", {})
    if receipt.get("feature_archive_ready_for_training") is not True or verified.get("md5") != FEATURE_MD5:
        raise ValueError("Figshare feature archive receipt is not verified")
    output_dir.mkdir(parents=True, exist_ok=True)
    records = {}
    with tarfile.open(archive, "r:gz") as handle:
        members = handle.getmembers()
        observed = {member.name: member.size for member in members}
        if observed != EXPECTED_MEMBERS or not all(_safe(member.name) for member in members):
            raise ValueError("Figshare feature archive member contract changed")
        for name in AGGREGATES:
            member = handle.getmember(name)
            source = handle.extractfile(member)
            if source is None:
                raise ValueError(f"cannot read {name}")
            destination = output_dir / Path(name).name
            temporary = destination.with_suffix(destination.suffix + ".partial")
            with temporary.open("wb") as output:
                shutil.copyfileobj(source, output, 4 * 1024 * 1024)
            if temporary.stat().st_size != member.size:
                raise ValueError(f"short aggregate extraction for {name}")
            temporary.replace(destination)
            records[destination.name] = {
                "archive_member": name,
                "size_bytes": destination.stat().st_size,
                "sha256": _sha256(destination),
            }
    return {
        "schema": "aleph.outer_library.figshare_cell_death_aggregate_receipt.v1",
        "source_doi": "10.17044/scilifelab.28202864.v2",
        "source_archive_md5": FEATURE_MD5,
        "archive_member_count": len(EXPECTED_MEMBERS),
        "archive_paths_safe": True,
        "single_cell_profiles_extracted": False,
        "aggregate_files": records,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--acquisition-receipt", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    report = extract(args.archive, args.acquisition_receipt, args.output_dir)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["aggregate_files"], sort_keys=True))


if __name__ == "__main__":
    main()

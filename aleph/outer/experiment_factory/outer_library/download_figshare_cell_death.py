#!/usr/bin/env python3
"""Acquire the public Cell Painting programmed-cell-death corpus.

The small design/QC files are enough to audit provenance and split hierarchy.
The much larger feature archive is optional so metadata verification does not
silently trigger a multi-gigabyte transfer.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import tarfile
import urllib.request

from download_sciplex3_metadata import _download


ARTICLE_ID = 28_202_864
ARTICLE_DOI = "10.17044/scilifelab.28202864.v2"
BASE = "https://ndownloader.figshare.com/files"
SMALL_FILES = {
    "metadata_celldeath_paper.csv": {
        "id": 51_672_443,
        "size": 1_890_409,
        "md5": "a92d4fddf20e6ca76473e5cc74c25a8a",
        "sha256": "e40c968af545947bb45427037ef046ba128133bddb6670b3242005638d32d9a0",
    },
    "README.txt": {
        "id": 51_791_030,
        "size": 4_293,
        "md5": "b0173f7a856fed5a8ecb6a638b83b57a",
        "sha256": "04604fed1b9870420688ed4f185b2b5e1948e77d5644983550195fbd305e7982",
    },
    "qc_df.csv": {
        "id": 54_974_975,
        "size": 1_642_736,
        "md5": "99e6cf3ea0e8fb59c8fa9b7ad5e01465",
        "sha256": "fc44402f50fe6d79be6ae28effe6c99a9e98c50719c05ef9fad1ace928a752b1",
    },
    "splits.tar.gz": {
        "id": 54_974_987,
        "size": 35_301,
        "md5": "68a67b79854a2657ce4a90126deebd28",
        "sha256": "91adadc4cdebfb9c50ee917ccf4d26e972d5f43b7f495d54c0fe84e7a78f7f76",
    },
}
LARGE_FILES = {
    "classification_splits.tar.gz": {
        "id": 54_975_083,
        "size": 1_848_046_769,
        "md5": "cde9d65431677bd142bba879549e226c",
    },
    "features.tar.gz": {
        "id": 54_975_104,
        "size": 2_451_859_048,
        "md5": "680c12ebc351f260f56ef722ffc7c2aa",
    },
    "embeddings.tar.gz": {
        "id": 54_975_200,
        "size": 4_318_114_942,
        "md5": "913aacc8b31df22c9bc60936e0ba3c2f",
    },
}
PROGRAMMED_DEATH_LABELS = {
    "apoptosis",
    "autophagy inducer",
    "ferroptosis inducer",
    "immunogenic cell death",
    "necroptosis inducer",
    "pyroptosis inducer",
}


def _digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _direct_download(url: str, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".partial")
    request = urllib.request.Request(
        url, headers={"User-Agent": "Project-Aleph-research/1.0"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        temporary.write_bytes(response.read())
    temporary.replace(destination)


def _verify(path: Path, expected: dict[str, object]) -> dict[str, object]:
    actual = {
        "size": path.stat().st_size,
        "md5": _digest(path, "md5"),
        "sha256": _digest(path, "sha256"),
    }
    for key in ("size", "md5", "sha256"):
        if key in expected and actual[key] != expected[key]:
            raise ValueError(f"Figshare checksum mismatch for {path.name}: {actual}")
    return actual


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def acquire(
    output_dir: Path,
    *,
    acquire_features: bool = False,
    segments: int = 384,
    workers: int = 8,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    files: dict[str, object] = {}
    for name, expected in SMALL_FILES.items():
        path = output_dir / name
        url = f"{BASE}/{expected['id']}"
        if not path.exists() or path.stat().st_size != expected["size"]:
            _direct_download(url, path)
        files[name] = {**_verify(path, expected), "url": url, "locally_acquired": True}

    metadata = _rows(output_dir / "metadata_celldeath_paper.csv")
    qc = _rows(output_dir / "qc_df.csv")
    if len(metadata) != 21_271 or len(qc) != 15_430:
        raise ValueError("Figshare cell-death table cardinality changed")
    plates = sorted({row["Metadata_Plate"] for row in metadata})
    labels = sorted({row["moa"] for row in qc if row["moa"]})
    if len(plates) != 8 or not PROGRAMMED_DEATH_LABELS.issubset(labels):
        raise ValueError("Figshare cell-death design axes changed")

    with tarfile.open(output_dir / "splits.tar.gz", "r:gz") as archive:
        split_members = sorted(
            member for member in archive.getnames() if member.endswith(".csv")
        )
        if split_members != [
            "splits/split_CellProfiler_moa.csv",
            "splits/split_DINO_moa.csv",
            "splits/split_DeepProfiler_moa.csv",
        ]:
            raise ValueError("Figshare author split archive members changed")
        split_summary: dict[str, object] = {}
        for member in split_members:
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ValueError(f"cannot read {member}")
            rows = list(csv.DictReader(io.TextIOWrapper(extracted, encoding="utf-8")))
            roles = {row["split"] for row in rows}
            wells = {(row["Metadata_Plate"], row["Metadata_Well"]) for row in rows}
            if roles != {"training", "valid", "test"} or len(wells) != len(rows):
                raise ValueError(f"invalid author well split in {member}")
            split_summary[Path(member).stem] = {
                "rows": len(rows),
                "unique_plate_wells": len(wells),
                "role_counts": {
                    role: sum(row["split"] == role for row in rows)
                    for role in sorted(roles)
                },
            }

    feature_expected = LARGE_FILES["features.tar.gz"]
    feature_path = output_dir / "features.tar.gz"
    if acquire_features:
        _download(
            f"{BASE}/{feature_expected['id']}",
            feature_path,
            int(feature_expected["size"]),
            segments=segments,
            workers=workers,
        )
        files["features.tar.gz"] = {
            **_verify(feature_path, feature_expected),
            "url": f"{BASE}/{feature_expected['id']}",
            "locally_acquired": True,
        }

    large_files = {
        name: {
            **expected,
            "url": f"{BASE}/{expected['id']}",
            "locally_acquired": (output_dir / name).exists()
            and (output_dir / name).stat().st_size == expected["size"],
        }
        for name, expected in LARGE_FILES.items()
    }
    return {
        "schema": "aleph.outer_library.figshare_cell_death_acquisition.v1",
        "article_id": ARTICLE_ID,
        "doi": ARTICLE_DOI,
        "publication_preprint_doi": "10.1101/2025.01.15.633042",
        "published_date": "2025-06-04",
        "license_id": "CC-BY-4.0",
        "cell_line": "MCF7",
        "treatment_hours": 48,
        "plate_count": len(plates),
        "metadata_site_rows": len(metadata),
        "qc_site_rows": len(qc),
        "programmed_cell_death_labels": sorted(PROGRAMMED_DEATH_LABELS),
        "all_qc_labels": labels,
        "author_split_summary": split_summary,
        "verified_local_files": files,
        "large_file_inventory": large_files,
        "feature_archive_ready_for_training": (
            "features.tar.gz" in files
            and files["features.tar.gz"]["md5"] == feature_expected["md5"]
        ),
        "independent_from_biad2515_provider": True,
        "same_endpoint_as_biad2515_annexin_v": False,
        "intended_role": "programmed_cell_death_mechanism_head_and_cross_assay_apoptosis_evidence",
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
        "retrieved_at": "2026-08-06",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--features", action="store_true")
    parser.add_argument("--segments", type=int, default=384)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    receipt = acquire(
        args.output_dir,
        acquire_features=args.features,
        segments=args.segments,
        workers=args.workers,
    )
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        key: receipt[key] for key in (
            "plate_count", "qc_site_rows", "feature_archive_ready_for_training"
        )
    }, sort_keys=True))


if __name__ == "__main__":
    main()

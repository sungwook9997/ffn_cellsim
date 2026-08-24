#!/usr/bin/env python3
"""Audit the CC0 Dryad contour-analysis tables without executing pickle files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import zipfile


ACCESSION = "10.5061/dryad.sj3tx9683"
ARCHIVE_SHA256 = "824cbe9904a906d8538c87765560a9bf5f14fc264b7b646976f025d41aa5d6e9"
README_SHA256 = "19a8b12320d05cf0b0d5df7e0b667de91e5e1fe121e9d02dccd20664167b1108"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe(name: str) -> bool:
    path = PurePosixPath(name)
    return not path.is_absolute() and ".." not in path.parts and "\\" not in name


def _csv(zf: zipfile.ZipFile, name: str) -> list[dict[str, str]]:
    with zf.open(name) as handle:
        return list(csv.DictReader(io.TextIOWrapper(handle, encoding="utf-8-sig")))


def audit(source_root: Path) -> dict[str, object]:
    archive = source_root / "analysed_data.zip"
    readme = source_root / "README.md"
    if _sha256(archive) != ARCHIVE_SHA256 or _sha256(readme) != README_SHA256:
        raise ValueError("Dryad force-propagation source hash changed")
    with zipfile.ZipFile(archive) as zf:
        if zf.testzip() is not None or not all(_safe(item.filename) for item in zf.infolist()):
            raise ValueError("unsafe or corrupt Dryad archive")
        baseline = _csv(zf, "analysed_data/baseline_data.csv")
        fullstim = _csv(zf, "analysed_data/fullstim_data.csv")
        stats = _csv(zf, "analysed_data/stats_fullstim.csv")
        pickle_members = sorted(
            item.filename for item in zf.infolist() if item.filename.endswith(".dat")
        )
    baseline_columns = [
        "line tension [nN]", "sigma_x_CM [mN/m]", "sigma_y_CM [mN/m]",
        "a [um]", "b [um]", "keys",
    ]
    if not baseline or list(baseline[0]) != baseline_columns:
        raise ValueError("baseline contour CSV contract changed")
    if not fullstim or list(fullstim[0]) != baseline_columns[:-1] + ["RSI_x", "RSI_y", "keys"]:
        raise ValueError("full-stimulation contour CSV contract changed")
    if len(baseline) != 68 or len(fullstim) != 68:
        raise ValueError("contour object count changed")
    exact_baseline_match = all(
        all(left[column] == right[column] for column in baseline_columns)
        for left, right in zip(baseline, fullstim)
    )
    raw_counts: dict[str, int] = {}
    for row in fullstim:
        raw_counts[row["keys"]] = raw_counts.get(row["keys"], 0) + 1
    stats_counts = {
        row["keys"]: int(row["count"]) for row in stats if row["quantity"] == "RSI_x"
    }
    return {
        "schema": "aleph.outer_library.dryad_sj3tx9683_source_audit.v1",
        "accession": ACCESSION, "paper_doi": "10.7554/eLife.83588",
        "provider": "Dryad Digital Repository", "license": "cc0-1.0",
        "sources": {
            "analysed_data.zip": {"sha256": ARCHIVE_SHA256, "size_bytes": archive.stat().st_size},
            "README.md": {"sha256": README_SHA256, "size_bytes": readme.stat().st_size},
        },
        "csv_contract": {
            "baseline_rows": len(baseline), "fullstim_rows": len(fullstim),
            "fullstim_repeats_baseline_columns_exactly": exact_baseline_match,
            "raw_fullstim_group_counts": raw_counts,
            "stats_fullstim_group_counts": stats_counts,
            "stats_fullstim_matches_raw_fullstim": stats_counts == raw_counts,
            "resolution": (
                "use raw fullstim_data.csv rows and withhold stats_fullstim.csv; the latter "
                "summarizes different AR1to1d/AR1to1s objects not present in that raw table"
            ),
        },
        "untrusted_pickle_members_not_executed": pickle_members,
        "sample_hierarchy": (
            "object rows are released without sample/date identifiers and remain one indivisible split group"
        ),
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.source_root)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()

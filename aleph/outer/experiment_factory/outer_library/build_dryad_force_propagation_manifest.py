#!/usr/bin/env python3
"""Normalize author-derived opto-MDCK contour observables from Dryad."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

from schema import Observation


DATASET = "dryad-sj3tx9683-opto-mdck-force-propagation"
LAB = "crbm-balland-ruppel-force-propagation"
ACCESSION = "10.5061/dryad.sj3tx9683"


def build(
    source_root: Path, source_audit: Path,
) -> tuple[list[Observation], dict[str, object]]:
    audit = json.loads(source_audit.read_text(encoding="utf-8"))
    if audit.get("schema") != "aleph.outer_library.dryad_sj3tx9683_source_audit.v1":
        raise ValueError("Dryad force-propagation audit schema changed")
    if audit.get("license") != "cc0-1.0" or audit.get("accession") != ACCESSION:
        raise ValueError("Dryad source identity/license changed")
    if audit["csv_contract"]["stats_fullstim_matches_raw_fullstim"]:
        raise ValueError("expected source discrepancy disappeared; re-audit admission policy")
    archive = source_root / "analysed_data.zip"
    source_sha = audit["sources"]["analysed_data.zip"]["sha256"]
    with zipfile.ZipFile(archive) as zf, zf.open("analysed_data/fullstim_data.csv") as handle:
        source = list(csv.DictReader(io.TextIOWrapper(handle, encoding="utf-8-sig")))
    metrics = (
        ("baseline_contour_line_tension", "line tension [nN]", "nN"),
        ("baseline_cortical_surface_tension_x", "sigma_x_CM [mN/m]", "mN_per_m"),
        ("baseline_cortical_surface_tension_y", "sigma_y_CM [mN/m]", "mN_per_m"),
        ("fitted_contour_axis_a", "a [um]", "micrometre"),
        ("fitted_contour_axis_b", "b [um]", "micrometre"),
        ("relative_surface_tension_increase_x", "RSI_x", "dimensionless"),
        ("relative_surface_tension_increase_y", "RSI_y", "dimensionless"),
    )
    states = {
        "AR1to2d": "opto_MDCK_doublet_aspect_ratio_1_to_2",
        "AR1to1d": "opto_MDCK_doublet_aspect_ratio_1_to_1",
        "AR2to1d": "opto_MDCK_doublet_aspect_ratio_2_to_1",
    }
    rows: list[Observation] = []
    for object_index, record in enumerate(source, start=1):
        if record["keys"] not in states:
            raise ValueError(f"unknown contour group: {record['keys']}")
        for observable, column, unit in metrics:
            locator = (
                "zip:analysed_data.zip;file=analysed_data/fullstim_data.csv;"
                f"row={object_index + 1};column={column}"
            )
            identity = hashlib.sha256(f"{source_sha}:{locator}".encode()).hexdigest()[:20]
            row = Observation(
                observation_id=f"{DATASET}:{identity}", dataset_id=DATASET,
                lab_group=LAB, provider="Dryad Digital Repository", accession=ACCESSION,
                license_id="cc0-1.0", source_sha256=source_sha, source_locator=locator,
                sample_id=f"{DATASET}:unresolved-author-sample-group",
                biological_replicate="author_sample_identity_not_released",
                technical_replicate=f"contour_object_{object_index:03d}",
                condition="global_optogenetic_RhoA_photoactivation",
                cell_type="opto_MDCK_canine_epithelial",
                cell_state=states[record["keys"]],
                modality="contour_model_from_TFM_MSM_and_actin_tracking",
                observable=observable, value=float(record[column]), unit=unit,
                unit_source_sha256=source_sha,
                unit_source_locator=(
                    "README.md Description of analysed data; CSV header; eLife 83588 Eq.75"
                ),
                coordinate_name="micropattern_aspect_ratio_code",
                coordinate_value=float({"AR1to2d": 0, "AR1to1d": 1, "AR2to1d": 2}[record["keys"]]),
                coordinate_unit="category_index",
            )
            row.validate()
            rows.append(row)
    if len(rows) != 476 or len({row.observation_id for row in rows}) != 476:
        raise ValueError("unexpected opto-MDCK manifest cardinality")
    report = {
        "schema": "aleph.outer_library.dryad_sj3tx9683_manifest.v1",
        "dataset_id": DATASET, "accession": ACCESSION, "paper_doi": "10.7554/eLife.83588",
        "observations": len(rows), "objects": len(source), "samples": 1,
        "author_sample_identifiers_available": False,
        "split_contract": "all 68 objects remain one indivisible dataset/sample group",
        "stats_fullstim_csv_admitted": False,
        "untrusted_dat_pickles_executed": False,
        "aleph_authority": "none",
    }
    return rows, report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--source-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    rows, report = build(args.source_root, args.source_audit)
    args.output.write_text(
        "".join(json.dumps(row.as_dict(), sort_keys=True) + "\n" for row in rows),
        encoding="utf-8",
    )
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()

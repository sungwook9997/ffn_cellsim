#!/usr/bin/env python3
"""Evaluate opto-MDCK force propagation as a non-authoritative sweep prior."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np


DATASET = "dryad-sj3tx9683-opto-mdck-force-propagation"


def evaluate(manifest: Path, manifest_report: Path, source_audit: Path) -> dict[str, object]:
    rows = [json.loads(line) for line in manifest.read_text(encoding="utf-8").splitlines()]
    design = json.loads(manifest_report.read_text(encoding="utf-8"))
    audit = json.loads(source_audit.read_text(encoding="utf-8"))
    if len(rows) != 476 or {row["dataset_id"] for row in rows} != {DATASET}:
        raise ValueError("opto-MDCK manifest changed")
    if design["author_sample_identifiers_available"] or design["samples"] != 1:
        raise ValueError("unreleased sample hierarchy was inflated")
    grouped: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        if row["observable"].startswith("relative_surface_tension_increase"):
            grouped[str(row["cell_state"])][str(row["observable"])].append(float(row["value"]))
    summaries = {}
    for state, observables in sorted(grouped.items()):
        summaries[state] = {}
        for observable, values in sorted(observables.items()):
            array = np.asarray(values)
            summaries[state][observable] = {
                "object_count": len(array), "median": float(np.median(array)),
                "mean": float(np.mean(array)), "positive_objects": int(np.sum(array > 0)),
            }
    all_x = np.asarray([
        float(row["value"]) for row in rows
        if row["observable"] == "relative_surface_tension_increase_x"
    ])
    all_y = np.asarray([
        float(row["value"]) for row in rows
        if row["observable"] == "relative_surface_tension_increase_y"
    ])
    return {
        "schema": "aleph.outer_library.dryad_sj3tx9683_evidence.v1",
        "dataset_id": DATASET,
        "lab_group": "crbm-balland-ruppel-force-propagation",
        "task": "RhoA_activation_cortical_surface_tension_direction",
        "manifest": design,
        "relative_surface_tension_increase_by_geometry": summaries,
        "overall_direction": {
            "x_positive_objects": int(np.sum(all_x > 0)), "x_total": len(all_x),
            "y_positive_objects": int(np.sum(all_y > 0)), "y_total": len(all_y),
            "usually_positive_in_both_axes": bool(np.sum(all_x > 0) > 34 and np.sum(all_y > 0) > 34),
        },
        "design_audit": {
            "objects_treated_as_biological_replicates": False,
            "dataset_lab_holdout_possible": False,
            "uncertainty_holdout_passed": False,
            "stats_fullstim_mismatch_detected": not audit["csv_contract"][
                "stats_fullstim_matches_raw_fullstim"
            ],
            "stats_fullstim_admitted": False,
            "pickle_files_executed": False,
        },
        "pre_sweep_constraint": {
            "purpose": "rank a future Aleph active-contractility sweep before forward simulation",
            "direction": "RhoA_activation_usually_increases_cortical_surface_tension",
            "required_Aleph_operators": [
                "micropattern_contour_surface_tension_x_mN_per_m",
                "micropattern_contour_surface_tension_y_mN_per_m",
            ],
            "identifiability": (
                "relative contour response does not identify motor kinetics, cortical material, "
                "stress-fibre line tension, and adhesion load transfer separately"
            ),
            "future_sweep_range_reduction_intended": True,
            "may_emit_numeric_parameter_range": False,
            "may_select_Aleph_parameter": False,
        },
        "observable_evidence_eligible": False,
        "independent_lab_holdout": False,
        "uncertainty_holdout_passed": False,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": [
            "author sample/date identifiers are absent from the released aggregate table",
            "raw TFM maps exist in a much larger archive but are not part of this selected acquisition",
            "stats_fullstim.csv does not summarize the released fullstim_data.csv groups",
            "no matched Aleph contour observation operator has passed forward sensitivity",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-report", type=Path, required=True)
    parser.add_argument("--source-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.manifest, args.manifest_report, args.source_audit)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()

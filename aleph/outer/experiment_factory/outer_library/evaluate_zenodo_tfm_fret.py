#!/usr/bin/env python3
"""Evaluate cross-lab TFM direction and TFM--FRET identifiability constraints."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
from typing import Any

import numpy as np


DATASET = "zenodo-14692589-tfm-fret"
REFERENCE_DATASET = "zenodo-7432971-tfm-tether"


def _load(path: Path, dataset: str) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        rows = [json.loads(line) for line in handle]
    if {row["dataset_id"] for row in rows} != {dataset}:
        raise ValueError(f"unexpected manifest dataset: {path}")
    return rows


def _summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    values = np.asarray([float(row["value"]) for row in rows], dtype=np.float64)
    return {
        "measurement_count": len(values),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "sample_standard_deviation": float(np.std(values, ddof=1)),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "unit": rows[0]["unit"],
    }


def _group(rows: list[dict[str, Any]], key: str) -> dict[str, dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[str(row[key])].append(row)
    return {name: _summary(group) for name, group in sorted(groups.items())}


def evaluate(
    manifest: Path, manifest_report: Path, reference_manifest: Path
) -> dict[str, Any]:
    rows = _load(manifest, DATASET)
    source_report = json.loads(manifest_report.read_text(encoding="utf-8"))
    if source_report.get("schema") != "aleph.outer_library.zenodo_14692589_manifest.v1":
        raise ValueError("Zenodo 14692589 manifest report schema changed")
    if len(rows) != 847 or len({row["sample_id"] for row in rows}) != 358:
        raise ValueError("Zenodo 14692589 manifest cardinality changed")

    traction = [
        row for row in rows
        if row["observable"] == "cell_mask_mean_traction_stress"
    ]
    fret = [
        row for row in rows
        if row["observable"] == "cell_averaged_focal_adhesion_FRET_efficiency"
    ]
    traction_by_condition = _group(traction, "condition")
    fret_by_state_condition: dict[str, dict[str, Any]] = {}
    for state in sorted({row["cell_state"] for row in fret}):
        fret_by_state_condition[state] = _group(
            [row for row in fret if row["cell_state"] == state], "condition"
        )

    soft = traction_by_condition["PAA_Young_modulus_4.5_kPa"]
    stiff = traction_by_condition["PAA_Young_modulus_13_kPa"]
    vints = fret_by_state_condition["vinculin_tension_sensor"]
    traction_direction = {
        "median_13_minus_4.5_Pa": stiff["median"] - soft["median"],
        "mean_13_minus_4.5_Pa": stiff["mean"] - soft["mean"],
        "increases_with_stiffness": stiff["median"] > soft["median"],
    }
    fret_direction = {
        "median_13_minus_4.5_percentage_points": (
            vints["PAA_Young_modulus_13_kPa"]["median"]
            - vints["PAA_Young_modulus_4.5_kPa"]["median"]
        ),
        "FRET_decreases_and_vinculin_tension_increases": (
            vints["PAA_Young_modulus_13_kPa"]["median"]
            < vints["PAA_Young_modulus_4.5_kPa"]["median"]
        ),
    }

    paired = [row for row in rows if row["modality"] == "paired_FA_TFM_FLIM_FRET_derived"]
    paired_by_cell: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in paired:
        paired_by_cell[row["sample_id"]].append(row)
    cell_correlations = {}
    for sample_id, group in sorted(paired_by_cell.items()):
        traction_map = {
            row["technical_replicate"]: float(row["value"])
            for row in group if row["observable"].startswith("focal_adhesion_mean_traction")
        }
        fret_map = {
            row["technical_replicate"]: float(row["value"])
            for row in group if row["observable"].startswith("focal_adhesion_FRET")
        }
        shared = sorted(set(traction_map) & set(fret_map))
        x = np.asarray([traction_map[key] for key in shared], dtype=np.float64)
        y = np.asarray([fret_map[key] for key in shared], dtype=np.float64)
        slope, intercept = np.polyfit(x, y, 1)
        category = group[0]["cell_state"]
        cell_correlations[category] = {
            "paired_focal_adhesions": len(shared),
            "FRET_on_traction_slope": float(slope),
            "intercept": float(intercept),
            "pearson_correlation": float(np.corrcoef(x, y)[0, 1]),
        }
    categories_reproduced = bool(
        cell_correlations["directly_correlated_cell"]["FRET_on_traction_slope"] < 0
        and cell_correlations["inversely_correlated_cell"]["FRET_on_traction_slope"] > 0
        and abs(cell_correlations["uncorrelated_cell"]["pearson_correlation"]) < 0.1
    )

    reference = _load(reference_manifest, REFERENCE_DATASET)
    reference = [
        row for row in reference
        if row["cell_type"] == "NIH_3T3_fibroblast"
        and row["observable"] == "time_averaged_mean_traction_stress"
    ]
    reference_by_condition = _group(reference, "condition")
    reference_increases = (
        reference_by_condition["PAA_shear_modulus_10_kPa"]["median"]
        > reference_by_condition["PAA_shear_modulus_1_kPa"]["median"]
    )
    cross_lab = {
        "same_observable": "cell_or_time_averaged_mean_traction_stress",
        "exact_common_unit": "Pa",
        "training_lab": "Cambridge_Franze_Keyser_collaboration",
        "external_direction_lab": source_report["lab_group"],
        "training_NIH3T3_1_to_10_kPa_direction_increases": reference_increases,
        "external_MEF_4.5_to_13_kPa_direction_increases": traction_direction[
            "increases_with_stiffness"
        ],
        "independent_lab_direction_reproduced": bool(
            reference_increases and traction_direction["increases_with_stiffness"]
        ),
        "numeric_transfer_evaluated": False,
        "reason_numeric_transfer_withheld": (
            "cell type and stiffness ranges differ, and neither released table maps cells "
            "to enough biological replicate identifiers"
        ),
    }

    return {
        "schema": "aleph.outer_library.zenodo_14692589_tfm_fret_evidence.v1",
        "task": "cross_lab_stiffness_TFM_and_vinculin_load_transfer_pre_sweep_constraint",
        "source": {
            "dataset_id": DATASET,
            "dataset_doi": "10.5281/zenodo.14692589",
            "article_doi": "10.1038/s42003-026-09514-0",
            "license_id": "cc-by-4.0",
        },
        "lab_group": source_report["lab_group"],
        "cell_level_mean_traction_by_condition": traction_by_condition,
        "cell_level_FRET_by_sensor_and_condition": fret_by_state_condition,
        "stiffness_directions": {
            "TFM": traction_direction,
            "VinTS_FLIM_FRET": fret_direction,
        },
        "paired_FA_cell_correlations": cell_correlations,
        "author_correlation_categories_reproduced": categories_reproduced,
        "cross_lab_same_observable_direction": cross_lab,
        "design_audit": {
            "author_reports_at_least_three_replicates": True,
            "cell_to_replicate_mapping_released": False,
            "cells_treated_as_biological_replicates": False,
            "focal_adhesions_treated_as_biological_replicates": False,
            "independent_dataset_direction_test": True,
            "independent_lab_direction_test": True,
            "independent_lab_inferential_holdout": False,
            "calibrated_uncertainty_holdout": False,
        },
        "observation_operator_audit": {
            "cell_mask_mean_traction_stress": {
                "external_unit": "Pa",
                "author_equation_and_open_source_code_available": True,
                "candidate_aleph_expression": (
                    "sum(cell-mask traction magnitude times pixel area) / cell-mask area"
                ),
                "status": "external_operator_defined_Aleph_geometry_equivalence_missing",
            },
            "vinculin_FLIM_FRET": {
                "external_unit": "relative_percent_or_per_cell_normalized_fraction",
                "status": "relative_molecular_tension_proxy_not_absolute_force",
                "identifiability_value": (
                    "separates focal-adhesion load transmission from bulk traction but "
                    "requires sensor calibration before absolute molecular force"
                ),
            },
        },
        "pre_sweep_constraint": {
            "latent": "focal_adhesion_load_transfer_under_substrate_stiffness",
            "future_sweep_range_reduction_intended": True,
            "current_executable_Aleph_sweep_available": False,
            "candidate_parameters": [
                "aleph.vertical.focal_adhesion.ClutchCard.bond_stiffness_pn_per_um",
                "aleph.vertical.focal_adhesion.ClutchCard.anchor_stiffness_pn_per_um",
                "aleph.vertical.sf_arc.MaterialCard.active_tension_pn",
                "aleph.vertical.nmii.MotorKinetics.crossbridge_stiffness",
            ],
            "required_Aleph_operators": [
                "cell_mask_mean_traction_stress_Pa",
                "focal_adhesion_mean_traction_normalized_per_cell",
                "vinculin_tension_sensor_relative_FRET",
                "focal_adhesion_area_orientation_and_vinculin_density",
            ],
            "numeric_range": None,
            "may_emit_numeric_parameter_range": False,
            "may_select_Aleph_parameter": False,
        },
        "observable_evidence_eligible": False,
        "uncertainty_holdout_passed": False,
        "independent_lab_holdout": False,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
        "limitations": source_report["adversarial_source_findings"] + [
            "independent-lab direction agreement is not a numeric transfer test",
            "three representative paired cells cannot train a general DCC/ICC/UCC classifier",
            "heterogeneous FA-level traction--FRET relations forbid a one-to-one load-transfer mapping",
            "no calibrated uncertainty, OOD test, or Aleph forward sensitivity is available",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-report", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.manifest, args.manifest_report, args.reference_manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "cross_lab_direction": report["cross_lab_same_observable_direction"],
        "categories_reproduced": report["author_correlation_categories_reproduced"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

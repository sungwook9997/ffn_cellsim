#!/usr/bin/env python3
"""Audit TFM/tether evidence as a pre-sweep observation constraint.

The figure-table cells are descriptive measurement units, not biological
replicates.  This evaluator therefore records directions and exact unit
relationships while refusing to manufacture an executable Aleph sweep.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
from typing import Any

import numpy as np


DATASET = "zenodo-7432971-tfm-tether"


def _load(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        rows = [json.loads(line) for line in handle]
    if {row["dataset_id"] for row in rows} != {DATASET}:
        raise ValueError("manifest is not Zenodo 7432971")
    if len(rows) != 1_024 or len({row["sample_id"] for row in rows}) != 380:
        raise ValueError("Zenodo 7432971 manifest cardinality changed")
    return rows


def _summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    values = np.asarray([float(row["value"]) for row in rows], dtype=np.float64)
    return {
        "cell_measurement_count": int(len(values)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "unit": rows[0]["unit"],
        "cells_are_not_independent_biological_replicates": True,
    }


def _grouped(
    rows: list[dict[str, Any]], keys: tuple[str, ...]
) -> dict[str, dict[str, Any]]:
    groups: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[tuple(str(row[key]) for key in keys)].append(row)
    return {"|".join(key): _summary(group) for key, group in sorted(groups.items())}


def evaluate(manifest: Path) -> dict[str, Any]:
    rows = _load(manifest)
    fibroblast_tfm = [
        row for row in rows
        if row["cell_type"] == "NIH_3T3_fibroblast"
        and row["modality"] == "traction_force_microscopy_derived"
    ]
    neuron_tfm = [
        row for row in rows
        if row["cell_type"] == "Xenopus_laevis_retinal_ganglion_growth_cone"
    ]
    tether = [
        row for row in rows
        if row["modality"] == "optical_tweezers_membrane_tether_derived"
    ]

    fibroblast_mean = [
        row for row in fibroblast_tfm
        if row["observable"] == "time_averaged_mean_traction_stress"
    ]
    fibroblast_by_stiffness = _grouped(fibroblast_mean, ("condition",))
    fibroblast_by_batch = _grouped(
        fibroblast_mean, ("condition", "biological_replicate")
    )
    fibroblast_direction_by_batch = {}
    for batch in ("author_file_batch_1_status_unresolved", "author_file_batch_2_status_unresolved"):
        soft = fibroblast_by_batch[
            f"PAA_shear_modulus_1_kPa|{batch}"
        ]["median"]
        stiff = fibroblast_by_batch[
            f"PAA_shear_modulus_10_kPa|{batch}"
        ]["median"]
        fibroblast_direction_by_batch[batch] = {
            "median_10kPa_minus_1kPa_Pa": float(stiff - soft),
            "direction": "increases" if stiff > soft else "decreases",
        }

    neuron_mean = [
        row for row in neuron_tfm
        if row["observable"] == "time_averaged_mean_traction_stress"
    ]
    neuron_pull = [
        row for row in neuron_tfm
        if row["observable"] == "time_averaged_growth_cone_pull_force"
    ]
    neuron_dates = _grouped(neuron_mean, ("condition", "biological_replicate"))
    shared_date = "author_acquisition_date_20210603_donor_status_unresolved"
    shared_date_direction = {
        observable: {
            condition: _summary([
                row for row in neuron_tfm
                if row["observable"] == observable
                and row["condition"] == condition
                and row["biological_replicate"] == shared_date
            ])
            for condition in ("PAA_shear_modulus_100_Pa", "PAA_shear_modulus_300_Pa")
        }
        for observable in (
            "time_averaged_mean_traction_stress",
            "time_averaged_growth_cone_pull_force",
        )
    }

    operators = {
        "membrane_tether_steady_state_force": {
            "external_unit": "pN",
            "candidate_aleph_observable": "membrane_tether_plateau_force",
            "unit_identity": "pN == pN",
            "status": "unit_exact_identifiability_blocked",
            "blocker": (
                "steady tether force constrains apparent tension jointly with bending "
                "rigidity and membrane-cortex adhesion; one Aleph parameter is not identified"
            ),
        },
        "time_averaged_mean_traction_stress": {
            "external_unit": "Pa",
            "candidate_aleph_expression": (
                "clutch_traction_constituent_pn / contact_area_um2"
            ),
            "unit_identity": "1 pN/um2 == 1 Pa",
            "status": "unit_exact_ROI_semantics_blocked",
            "blocker": (
                "the author TFM ROI/mask averaging operator is not encoded in the figure "
                "table and has not been shown equivalent to the Aleph contact-area operator"
            ),
        },
        "time_averaged_growth_cone_pull_force": {
            "external_unit": "pN",
            "candidate_aleph_observable": None,
            "status": "geometry_and_cell_type_operator_missing",
        },
        "membrane_tether_peak_detachment_force": {
            "external_unit": "pN",
            "candidate_aleph_observable": None,
            "status": "detachment_peak_operator_missing",
        },
    }

    return {
        "schema": "aleph.outer_library.zenodo_7432971_tfm_tether_evidence.v1",
        "task": "TFM_and_membrane_tether_pre_sweep_observation_constraints",
        "source": {
            "dataset_id": DATASET,
            "dataset_doi": "10.5281/zenodo.7432971",
            "article_doi": "10.1093/pnasnexus/pgac299",
            "license_id": "cc-by-4.0",
        },
        "purpose": (
            "accumulate externally grounded observables that can reduce a future Aleph "
            "parameter sweep after exact operators and engine-scale sweeps become available"
        ),
        "descriptive_results": {
            "fibroblast_mean_traction_by_stiffness": fibroblast_by_stiffness,
            "fibroblast_mean_traction_by_author_file_batch": fibroblast_by_batch,
            "fibroblast_stiffness_direction_by_batch": fibroblast_direction_by_batch,
            "neuron_mean_traction_by_stiffness": _grouped(neuron_mean, ("condition",)),
            "neuron_pull_force_by_stiffness": _grouped(neuron_pull, ("condition",)),
            "neuron_mean_traction_by_date": neuron_dates,
            "neuron_shared_acquisition_date_descriptive_comparison": shared_date_direction,
            "tether_force_by_cell_state_substrate_observable": _grouped(
                tether, ("cell_type", "cell_state", "condition", "observable")
            ),
        },
        "design_audit": {
            "TFM_frames_collapsed_before_analysis": True,
            "cells_treated_as_biological_replicates": False,
            "fibroblast_author_file_batches_per_stiffness": 2,
            "fibroblast_batch_effect_present": True,
            "neuron_stiffness_partly_confound_with_acquisition_date": True,
            "neuron_only_shared_acquisition_date": "20210603",
            "optical_tweezer_biological_replicate_ids_available": False,
            "independent_lab_inferential_holdout_established": False,
        },
        "observation_operator_audit": operators,
        "promotion": {
            "stage": "external_training_or_descriptive_constraint_only",
            "future_sweep_range_reduction_intended": True,
            "current_executable_Aleph_sweep_available": False,
            "observable_evidence_eligible": False,
            "numeric_parameter_likelihood_eligible": False,
            "aleph_parameter_proposal_eligible": False,
            "may_emit_numeric_parameter_range": False,
            "numeric_range": None,
        },
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
        "limitations": [
            "cell-level rows are not independent biological replicates",
            "fibroblast stiffness groups contain only two unresolved author file batches",
            "neuron stiffness and acquisition date are partly confounded",
            "optical-tweezer biological replicate identifiers are unavailable",
            "exact units alone do not establish an equivalent Aleph observation operator",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report["promotion"], sort_keys=True))


if __name__ == "__main__":
    main()

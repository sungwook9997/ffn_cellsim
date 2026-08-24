#!/usr/bin/env python3
"""Evaluate replicated C2C12 PIV/IF evidence as a pre-sweep constraint."""

from __future__ import annotations

import argparse
from collections import defaultdict
import itertools
import json
from pathlib import Path

import numpy as np


DATASET = "dryad-08kprr59c-c2c12-supracontractility"
CONTROL = "mLCN_1.0_vehicle_control"
BLEB = "mLCN_1.0_blebbistatin_10_micromolar_daily"


def _cluster_means(
    rows: list[dict[str, object]], condition: str, observable: str,
    *, minimum_time: float | None = None,
) -> dict[str, float]:
    grouped: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        if row["condition"] != condition or row["observable"] != observable:
            continue
        if minimum_time is not None:
            if row["coordinate_name"] != "time_after_confluence":
                continue
            if float(row["coordinate_value"]) < minimum_time:
                continue
        grouped[str(row["biological_replicate"])].append(float(row["value"]))
    return {key: float(np.mean(values)) for key, values in grouped.items()}


def _bootstrap_difference(
    group_a: list[float], group_b: list[float], seed: int,
    group_a_label: str, group_b_label: str,
) -> dict[str, object]:
    if len(group_a) < 2 or len(group_b) < 2:
        return {
            "difference_definition": f"{group_a_label}_minus_{group_b_label}",
            "group_A_minus_group_B": float(np.mean(group_a) - np.mean(group_b)),
            "biological_cluster_bootstrap_95ci": None,
            "exact_randomization_two_sided_p": None,
            "reason": "fewer_than_two_biological_replicates_in_one_group",
        }
    rng = np.random.default_rng(seed)
    control_array = np.asarray(group_a, dtype=np.float64)
    treated_array = np.asarray(group_b, dtype=np.float64)
    draws = np.empty(20_000, dtype=np.float64)
    for index in range(len(draws)):
        draws[index] = (
            rng.choice(control_array, len(control_array), replace=True).mean()
            - rng.choice(treated_array, len(treated_array), replace=True).mean()
        )
    pooled = group_a + group_b
    observed = abs(float(np.mean(group_a) - np.mean(group_b)))
    randomization = []
    for selected in itertools.combinations(range(len(pooled)), len(group_a)):
        selected_set = set(selected)
        a = [pooled[i] for i in selected]
        b = [value for i, value in enumerate(pooled) if i not in selected_set]
        randomization.append(abs(float(np.mean(a) - np.mean(b))))
    return {
        "difference_definition": f"{group_a_label}_minus_{group_b_label}",
        "group_A_label": group_a_label, "group_B_label": group_b_label,
        "group_A_replicate_means": group_a,
        "group_B_replicate_means": group_b,
        "group_A_minus_group_B": float(np.mean(group_a) - np.mean(group_b)),
        "biological_cluster_bootstrap_95ci": [
            float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))
        ],
        "exact_randomization_two_sided_p": float(
            np.mean(np.asarray(randomization) >= observed - 1e-12)
        ),
        "all_group_A_replicates_above_all_group_B_replicates": (
            min(group_a) > max(group_b)
        ),
        "inference_unit": "author_independent_LCN_substrate",
        "bootstrap_is_exploratory_due_small_n": min(len(group_a), len(group_b)) <= 2,
    }


def _fibronectin_order(rows: list[dict[str, object]]) -> dict[str, float]:
    grouped: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in rows:
        if row["observable"] != "fibronectin_orientation_frequency":
            continue
        grouped[str(row["biological_replicate"])].append(
            (float(row["coordinate_value"]), float(row["value"]))
        )
    result = {}
    for replicate, pairs in grouped.items():
        angle = np.deg2rad([pair[0] for pair in pairs])
        weight = np.asarray([pair[1] for pair in pairs], dtype=np.float64)
        if len(pairs) != 180 or not np.isclose(weight.sum(), 100.0, atol=1e-6):
            raise ValueError(f"invalid fibronectin angular distribution: {replicate}")
        result[replicate] = float(np.hypot(
            np.sum(weight * np.cos(2 * angle)),
            np.sum(weight * np.sin(2 * angle)),
        ) / weight.sum())
    return result


def evaluate(manifest: Path) -> dict[str, object]:
    rows = [json.loads(line) for line in manifest.read_text(encoding="utf-8").splitlines()]
    if len(rows) != 2_155 or len({row["sample_id"] for row in rows}) != 16:
        raise ValueError("unexpected C2C12 supracontractility manifest cardinality")
    if {row["dataset_id"] for row in rows} != {DATASET}:
        raise ValueError("unexpected C2C12 dataset identity")
    if {row["aleph_authority"] for row in rows} != {"none"}:
        raise ValueError("external observations acquired Aleph authority")
    if len({row["biological_replicate"] for row in rows}) != 16:
        raise ValueError("author LCN replicates were collapsed or inflated")

    late_contrasts = {}
    for index, observable in enumerate((
        "orientation_order_parameter", "nematic_correlation_length",
        "velocity_correlation_length", "cell_speed",
    )):
        control = list(_cluster_means(rows, CONTROL, observable, minimum_time=42).values())
        treated = list(_cluster_means(rows, BLEB, observable, minimum_time=42).values())
        late_contrasts[observable] = _bootstrap_difference(
            control, treated, 80 + index, "vehicle_control", "blebbistatin"
        )

    fibronectin = _fibronectin_order(rows)
    control_fibronectin = [
        value for replicate, value in fibronectin.items() if "vehicle_control" in replicate
    ]
    treated_fibronectin = [
        value for replicate, value in fibronectin.items() if "blebbistatin" in replicate
    ]
    fibronectin_contrast = _bootstrap_difference(
        control_fibronectin, treated_fibronectin, 91,
        "vehicle_control", "blebbistatin",
    )

    anisotropic = "mLCN_1.0_mechanically_anisotropic"
    isotropic = "iLCN_1.0_mechanically_isotropic"
    max_density = {}
    for index, observable in enumerate((
        "orientation_order_parameter", "spatial_disorder_fraction",
        "velocity_correlation_length",
    )):
        a = [
            float(row["value"]) for row in rows
            if row["condition"] == anisotropic and row["observable"] == observable
            and row["coordinate_value"] == 36
        ]
        b = [
            float(row["value"]) for row in rows
            if row["condition"] == isotropic and row["observable"] == observable
            and row["coordinate_value"] == 36
        ]
        if len(a) != 5 or len(b) != 4:
            raise ValueError(f"maximum-density replicate contract changed: {observable}")
        max_density[observable] = _bootstrap_difference(
            a, b, 100 + index, "mLCN_anisotropic", "iLCN_isotropic"
        )

    replicated_mechanism_direction = all(
        late_contrasts[name].get("all_group_A_replicates_above_all_group_B_replicates")
        for name in (
            "orientation_order_parameter", "nematic_correlation_length",
            "velocity_correlation_length",
        )
    ) and bool(
        fibronectin_contrast["all_group_A_replicates_above_all_group_B_replicates"]
    )
    return {
        "schema": "aleph.outer_library.dryad_08kprr59c_evidence.v1",
        "task": "C2C12_supracellular_contractility_collective_flow_and_ECM_alignment_constraint",
        "source": {
            "accession": "10.5061/dryad.08kprr59c",
            "paper_doi": "10.1126/sciadv.adn0235", "license": "cc0-1.0",
            "lab_group": "anseth-white-skillin-lab-cu-boulder",
            "independent_provider_relative_to_existing_outer_library_PIV_sources": True,
        },
        "manifest": {
            "observations": len(rows), "samples": 16,
            "author_independent_LCN_substrates": 16,
            "PIV_live_actin_observations": sum(
                row["modality"] == "particle_image_velocimetry_and_live_actin_derived"
                for row in rows
            ),
            "fibronectin_IF_orientation_observations": sum(
                row["modality"] == "immunofluorescence_fibronectin_orientation_derived"
                for row in rows
            ),
        },
        "author_exclusion_audit": {
            "red_filled_cells_not_emitted": 112,
            "green_maximum_density_cells": 9,
            "maximum_density_time_after_confluence_h": 36,
            "style_annotations_verified_against_OOXML": True,
        },
        "maximum_density_anisotropy_contrasts": max_density,
        "late_post_blebbistatin_contrasts_t_ge_42h": late_contrasts,
        "fibronectin_orientation_nematic_order": {
            "per_replicate": fibronectin,
            "control_minus_blebbistatin": fibronectin_contrast,
        },
        "replicated_contractility_direction": replicated_mechanism_direction,
        "design_audit": {
            "split_unit": "author_independent_LCN_substrate",
            "timepoints_are_nested_technical_observations": True,
            "observables_are_nested_technical_observations": True,
            "cross_split_biological_replicate_overlap_allowed": False,
            "independent_dataset_holdout": False,
            "independent_lab_holdout": False,
            "reason": "this acquisition adds an independent provider and replicated perturbation, but no same-task model was trained on another aligned lab and sealed before this evaluation",
        },
        "uncertainty": {
            "cluster_bootstrap_unit": "LCN_substrate",
            "n_per_Fig4_condition": 2,
            "formal_uncertainty_holdout_passed": False,
            "reason": "two substrates per treatment are insufficient for calibrated external uncertainty",
        },
        "OOD_refusal": {
            "declared_domains": [
                "2D_C2C12_cell_line", "MPa_to_GPa_LCN_substrate",
                "mLCN_1.0_stiffness_anisotropy", "10_uM_blebbistatin",
            ],
            "refuse": [
                "primary_or_iPSC_muscle_without_external_validation",
                "physiologic_soft_or_3D_matrix_transfer",
                "direct_force_or_traction_magnitude_inference",
                "numeric_Aleph_active_tension_inference",
            ],
        },
        "pre_sweep_constraint": {
            "latent": "supracellular_contractility_and_reciprocal_cell_ECM_alignment",
            "future_sweep_range_reduction_intended": True,
            "directional_observation_constraint_eligible": True,
            "expected_when_contractility_is_inhibited": {
                "late_orientation_order_parameter": "decreases",
                "late_nematic_correlation_length": "decreases",
                "late_velocity_correlation_length": "decreases",
                "fibronectin_orientation_nematic_order": "decreases",
            },
            "required_Aleph_operators": [
                "live_actin_orientation_order", "nematic_correlation_length",
                "velocity_correlation_length", "fibronectin_orientation_distribution",
            ],
            "may_emit_numeric_parameter_range": False,
            "may_select_Aleph_parameter": False,
            "identifiability": "blebbistatin is a multiscale nonmuscle-myosin-II perturbation and does not identify one Aleph parameter",
        },
        "representation_training_eligible": True,
        "observable_evidence_eligible": False,
        "uncertainty_holdout_passed": False,
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": [
            "Fig. 4 has only two independent LCN substrates per condition",
            "blebbistatin-treated speed is available for only one substrate and is descriptive only",
            "the substrate is orders of magnitude stiffer than physiological muscle and is two-dimensional",
            "PIV, order, and correlation lengths are observations rather than direct traction or tension measurements",
            "another aligned lab and an Aleph observation operator are required before this head can constrain a numeric sweep",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "replicated_contractility_direction": report["replicated_contractility_direction"],
        "observable_evidence_eligible": report["observable_evidence_eligible"],
        "samples": report["manifest"]["samples"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

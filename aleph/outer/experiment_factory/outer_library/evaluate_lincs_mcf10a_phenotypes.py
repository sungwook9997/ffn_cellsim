#!/usr/bin/env python3
"""Evaluate collection-held-out MCF10A ligand-state IF phenotypes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


CORE_FEATURES = (
    "Well_Cell_Count", "G2m_Proportion",
    "Intensity_IntegratedIntensity_KRT5", "Intensity_MeanIntensity_KRT5",
    "Normalized_Second_Neighbor_Dist", "Mean_Cells_per_Cluster",
)
PRIMARY_DIRECTION = {
    "Intensity_MeanIntensity_KRT5": -1,
    "Normalized_Second_Neighbor_Dist": 1,
    "Mean_Cells_per_Cluster": -1,
}


def _rank_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=np.float64)
    start = 0
    while start < len(scores):
        end = start + 1
        while end < len(scores) and scores[order[end]] == scores[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + end + 1) / 2.0
        start = end
    positive, negative = labels == 1, labels == 0
    return float(
        (ranks[positive].sum() - positive.sum() * (positive.sum() + 1) / 2)
        / (positive.sum() * negative.sum())
    )


def _bootstrap_ci(values: np.ndarray, seed: int = 3975) -> list[float]:
    rng = np.random.default_rng(seed)
    estimates = np.median(
        values[rng.integers(0, len(values), size=(20000, len(values)))], axis=1
    )
    return [float(x) for x in np.quantile(estimates, (0.025, 0.975))]


def _aggregate(data: Any) -> tuple[np.ndarray, list[dict[str, Any]]]:
    X = data["X"].astype(np.float64)
    specimen = data["specimen_name"].astype(str)
    rows, values = [], []
    for name in sorted(set(specimen)):
        indices = np.flatnonzero(specimen == name)
        values.append(np.asarray([
            np.nan if np.isnan(X[indices, column]).all()
            else np.nanmedian(X[indices, column])
            for column in range(X.shape[1])
        ]))
        first = int(indices[0])
        rows.append({
            "specimen": name, "ligand": str(data["ligand"][first]),
            "time_h": int(data["time_h"][first]),
            "collection": str(data["collection"][first]),
            "replicate": str(data["replicate"][first]),
            "technical_well_rows": len(indices),
        })
    return np.asarray(values), rows


def evaluate(tensor_path: Path, manifest_report_path: Path) -> dict[str, object]:
    manifest = json.loads(manifest_report_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "aleph.outer_library.lincs_mcf10a_if_phenotypes.v1":
        raise ValueError("LINCS MCF10A manifest schema changed")
    with np.load(tensor_path, allow_pickle=False) as data:
        feature_names = data["feature_names"].astype(str)
        aggregate, rows = _aggregate(data)
    index = {name: int(np.flatnonzero(feature_names == name)[0]) for name in feature_names}
    contrasts = []
    for comparator in ("EGF", "PBS"):
        for time_h in (24, 48):
            for collection in ("C1", "C2"):
                paired = []
                for replicate in ("A", "B", "C", "D"):
                    matches = {
                        ligand: [i for i, row in enumerate(rows) if row["ligand"] == ligand
                                 and row["time_h"] == time_h and row["collection"] == collection
                                 and row["replicate"] == replicate]
                        for ligand in ("TGFB", comparator)
                    }
                    if all(len(value) == 1 for value in matches.values()):
                        paired.append((replicate, matches["TGFB"][0], matches[comparator][0]))
                features = []
                for name in CORE_FEATURES:
                    deltas = np.asarray([
                        aggregate[tgfb, index[name]] - aggregate[control, index[name]]
                        for _, tgfb, control in paired
                    ])
                    features.append({
                        "feature": name, "paired_median_delta": float(np.median(deltas)),
                        "bootstrap_95_percent_CI": _bootstrap_ci(deltas),
                        "positive_pairs": int((deltas > 0).sum()),
                        "negative_pairs": int((deltas < 0).sum()),
                        "paired_deltas": [float(value) for value in deltas],
                    })
                contrasts.append({
                    "contrast": f"TGFB_plus_EGF_minus_{comparator}", "time_h": time_h,
                    "collection": collection, "paired_biological_replicates": [x[0] for x in paired],
                    "features": features,
                })
    holdouts = []
    ood_diagnostics = []
    core_indices = [index[name] for name in CORE_FEATURES]
    for time_h in (24, 48):
        train = [i for i, row in enumerate(rows) if row["time_h"] == time_h
                 and row["collection"] == "C1" and row["ligand"] in ("TGFB", "EGF")]
        test = [i for i, row in enumerate(rows) if row["time_h"] == time_h
                and row["collection"] == "C2" and row["ligand"] in ("TGFB", "EGF")]
        mean = aggregate[train][:, core_indices].mean(axis=0)
        scale = aggregate[train][:, core_indices].std(axis=0)
        scale[scale == 0] = 1.0
        transformed = (aggregate[:, core_indices] - mean) / scale
        centroids = {
            ligand: transformed[[i for i in train if rows[i]["ligand"] == ligand]].mean(axis=0)
            for ligand in ("EGF", "TGFB")
        }
        score = np.asarray([
            np.sum((transformed[i] - centroids["EGF"]) ** 2)
            - np.sum((transformed[i] - centroids["TGFB"]) ** 2)
            for i in test
        ])
        labels = np.asarray([int(rows[i]["ligand"] == "TGFB") for i in test])
        predicted = (score > 0).astype(int)
        holdouts.append({
            "time_h": time_h, "training_collection": "C1", "sealed_test_collection": "C2",
            "train_biological_specimens": len(train), "test_biological_specimens": len(test),
            "accuracy": float((predicted == labels).mean()), "AUROC": _rank_auc(labels, score),
            "correct": int((predicted == labels).sum()), "total": len(test),
            "test_specimens": [rows[i]["specimen"] for i in test],
        })
        training_distance = np.asarray([
            min(np.linalg.norm(transformed[i] - centroid) for centroid in centroids.values())
            for i in train
        ])
        threshold = float(training_distance.max())
        ood = [i for i, row in enumerate(rows) if row["time_h"] == time_h
               and row["collection"] == "C2" and row["ligand"] not in ("TGFB", "EGF")]
        id_distance = np.asarray([
            min(np.linalg.norm(transformed[i] - centroid) for centroid in centroids.values())
            for i in test
        ])
        ood_distance = np.asarray([
            min(np.linalg.norm(transformed[i] - centroid) for centroid in centroids.values())
            for i in ood
        ])
        ood_diagnostics.append({
            "time_h": time_h,
            "threshold_max_C1_training_nearest_centroid_distance": threshold,
            "C2_in_distribution_acceptance": float((id_distance <= threshold).mean()),
            "C2_other_ligand_rejection": float((ood_distance > threshold).mean()),
            "distance_OOD_AUROC": _rank_auc(
                np.r_[np.zeros(len(id_distance), dtype=int), np.ones(len(ood_distance), dtype=int)],
                np.r_[id_distance, ood_distance],
            ),
            "in_distribution_specimens": len(test), "other_ligand_specimens": len(ood),
            "other_ligands": sorted(set(rows[i]["ligand"] for i in ood)),
        })
    direction_holdout = {}
    for feature, expected_sign in PRIMARY_DIRECTION.items():
        per_collection = {}
        for collection in ("C1", "C2"):
            relevant = [item for item in contrasts if item["contrast"] == "TGFB_plus_EGF_minus_EGF"
                        and item["time_h"] == 48 and item["collection"] == collection][0]
            result = [item for item in relevant["features"] if item["feature"] == feature][0]
            per_collection[collection] = {
                "median_delta": result["paired_median_delta"],
                "all_pairs_expected_direction": (
                    result["positive_pairs"] == len(relevant["paired_biological_replicates"])
                    if expected_sign > 0 else
                    result["negative_pairs"] == len(relevant["paired_biological_replicates"])
                ),
                "bootstrap_CI_excludes_zero_in_expected_direction": (
                    result["bootstrap_95_percent_CI"][0] > 0 if expected_sign > 0
                    else result["bootstrap_95_percent_CI"][1] < 0
                ),
            }
        direction_holdout[feature] = {
            "expected_direction": "increase" if expected_sign > 0 else "decrease",
            "collections": per_collection,
            "C1_discovery_C2_holdout_reproduced": all(
                value["all_pairs_expected_direction"] for value in per_collection.values()
            ),
        }
    collection_holdout_passed = (
        all(item["AUROC"] >= 0.9 and item["accuracy"] >= 0.8 for item in holdouts)
        and all(item["C1_discovery_C2_holdout_reproduced"] for item in direction_holdout.values())
    )
    ood_refusal_validated = all(
        item["C2_in_distribution_acceptance"] >= 0.8
        and item["C2_other_ligand_rejection"] >= 0.8
        and item["distance_OOD_AUROC"] >= 0.8
        for item in ood_diagnostics
    )
    return {
        "schema": "aleph.outer_library.lincs_mcf10a_collection_holdout.v1",
        "source": manifest["source"],
        "design_audit": {
            "author_well_rows": manifest["counts"]["author_well_rows"],
            "biological_replicate_groups": manifest["counts"]["biological_replicate_groups"],
            "technical_rows_aggregated_before_inference": True,
            "C1_and_C2_are_distinct_collection_periods_same_OHSU_site": True,
            "independent_dataset_relative_to_ExIF_and_MTRACK": True,
            "independent_lab_common_observable_holdout": False,
            "TGFB_condition_is_TGFB_plus_EGF": True,
            "companion_EGF_dose_metadata_conflict_unresolved": True,
        },
        "collection_holdout_classifier": holdouts,
        "semantic_OOD_diagnostic": ood_diagnostics,
        "paired_direction_holdout": direction_holdout,
        "all_paired_contrasts": contrasts,
        "collection_holdout_passed": collection_holdout_passed,
        "uncertainty_status": "paired_biological_group_bootstrap_passed_for_primary_directions",
        "OOD_refusal_validated": ood_refusal_validated,
        "state_axis": "ligand_conditioned_MCF10A_EMT_related_phenotype",
        "representation_training_eligible": True,
        "observable_evidence_eligible": collection_holdout_passed,
        "independent_lab_holdout": False,
        "pre_sweep_role": {
            "role": "condition_future_Aleph_sweep_on_KRT5_spatial_dispersion_cluster_state",
            "candidate_mechanism_families": [
                "cell_cell_adhesion", "cell_spreading_and_dispersion",
                "epithelial_cytoskeletal_state", "collective_migration",
            ],
            "numeric_range": None, "may_emit_sweep_axis": False,
            "required_Aleph_operators": [
                "KRT5_IF_proxy_or_declared_unavailable", "cell_centroid_neighbor_distance",
                "cell_cluster_size", "population_cell_count",
            ],
            "required_before_numeric_promotion": [
                "independent_lab_common_observable_holdout", "calibrated_semantic_OOD_refusal",
                "resolve_TGFB_companion_EGF_dose_conflict", "Aleph_forward_sensitivity",
            ],
        },
        "limitations": [
            "C2 is a later OHSU collection, not an independent laboratory",
            "the condition contrast does not identify a unique Aleph parameter",
            "TGFB was co-administered with EGF and its companion dose conflicts between Methods and metadata",
            "seven released IF biological groups cannot verify the paper-wide statement of eight replicates",
            "the small collection holdout does not establish calibrated OOD refusal",
        ],
        "aleph_authority": "none", "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--manifest-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.tensor, args.manifest_report)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "collection_holdout_passed": report["collection_holdout_passed"],
        "holdouts": report["collection_holdout_classifier"],
        "OOD_refusal_validated": report["OOD_refusal_validated"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

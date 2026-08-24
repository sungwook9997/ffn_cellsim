#!/usr/bin/env python3
"""Evaluate a post-freeze external senescence confirmation without retraining."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path

import numpy as np


ALIASES = {"IL8": "CXCL8"}


def _sigmoid(value: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(value, -30.0, 30.0)))


def _rank_rows(values: np.ndarray) -> np.ndarray:
    order = np.argsort(np.argsort(values, axis=1, kind="stable"), axis=1, kind="stable")
    return order.astype(np.float64) / max(1, values.shape[1] - 1)


def _load_counts(
    path: Path, feature_names: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:8] != ["Geneid", "Chr", "Start", "End", "Strand", "Length", "gene_name", "gene_type"]:
            raise ValueError("GSE301164 count-table schema changed")
        sample_names = np.asarray(header[8:])
        wanted = {ALIASES.get(name, name): index for index, name in enumerate(feature_names)}
        counts = np.zeros((len(sample_names), len(feature_names)), dtype=np.float64)
        library_size = np.zeros(len(sample_names), dtype=np.float64)
        found: set[str] = set()
        for row in reader:
            symbol = row[6]
            values = np.asarray(row[8:], dtype=np.float64)
            library_size += values
            if symbol in wanted:
                counts[:, wanted[symbol]] += values
                found.add(symbol)
    missing = sorted(set(wanted) - found)
    if missing:
        raise ValueError(f"GSE301164 lacks frozen markers: {missing}")
    if np.any(library_size <= 0):
        raise ValueError("GSE301164 contains an empty sample")
    return counts, library_size, sample_names


def _one_sided_sign_p(differences: list[float]) -> float:
    positive = sum(value > 0 for value in differences)
    n = len(differences)
    return float(sum(math.comb(n, k) for k in range(positive, n + 1)) / 2 ** n)


def _contrast(probability: np.ndarray, samples: np.ndarray, treated: str, control: str) -> dict[str, object]:
    treated_values = []
    control_values = []
    differences = []
    sample_pairs = []
    for replicate in (1, 2, 3):
        treated_name = f"{treated}_rep{replicate}"
        control_name = f"{control}_rep{replicate}"
        treated_index = np.flatnonzero(samples == treated_name)
        control_index = np.flatnonzero(samples == control_name)
        if len(treated_index) != 1 or len(control_index) != 1:
            raise ValueError(f"missing predeclared pair {treated_name}/{control_name}")
        treated_value = float(probability[treated_index[0]])
        control_value = float(probability[control_index[0]])
        treated_values.append(treated_value)
        control_values.append(control_value)
        differences.append(treated_value - control_value)
        sample_pairs.append([treated_name, control_name])
    return {
        "sample_pairs": sample_pairs,
        "treated_probability": treated_values,
        "control_probability": control_values,
        "treated_minus_control_probability": differences,
        "positive_pair_count": sum(value > 0 for value in differences),
        "pair_count": len(differences),
        "mean_difference": float(np.mean(differences)),
        "exact_one_sided_sign_test_p": _one_sided_sign_p(differences),
    }


def evaluate(counts_path: Path, checkpoint_path: Path) -> dict[str, object]:
    with np.load(checkpoint_path, allow_pickle=False) as checkpoint:
        feature_names = np.asarray(checkpoint["feature_name"]).astype(str)
        representation = str(np.asarray(checkpoint["representation"]).item())
        mean = np.asarray(checkpoint["input_mean"], dtype=np.float64)
        scale = np.asarray(checkpoint["input_scale"], dtype=np.float64)
        weight = np.asarray(checkpoint["weight"], dtype=np.float64)
        bias = float(np.asarray(checkpoint["bias"]))
        temperature = float(np.asarray(checkpoint["temperature"]))
        ood_threshold = float(np.asarray(checkpoint["ood_threshold"]))
    counts, library_size, sample_names = _load_counts(counts_path, feature_names)
    log_cpm = np.log1p(counts / library_size[:, None] * 1_000_000.0)
    if representation == "logCPM_fit_standardized":
        base = log_cpm
    elif representation == "within_sample_rank_fit_standardized":
        base = _rank_rows(log_cpm)
    else:
        raise ValueError(f"unknown frozen representation {representation}")
    standardized = (base - mean) / scale
    probability = _sigmoid((standardized @ weight + bias) / temperature)
    distance = np.sqrt(np.mean(standardized ** 2, axis=1))

    replicative = _contrast(probability, sample_names, "BJ_P33", "BJ_P3")
    doxorubicin = _contrast(probability, sample_names, "BJ_P11_DOX_50ngml", "BJ_P11_DMSO")
    all_differences = (
        replicative["treated_minus_control_probability"]
        + doxorubicin["treated_minus_control_probability"]
    )
    combined_p = _one_sided_sign_p(all_differences)
    all_positive = all(value > 0 for value in all_differences)
    directional_evidence = all_positive and combined_p < 0.05
    return {
        "schema": "aleph.outer_library.gse301164_postfreeze_confirmation.v1",
        "dataset": "GSE301164",
        "lab": "Northwest_University_Yan_lab",
        "checkpoint_frozen_before_dataset_acquisition": True,
        "model_fit_selection_calibration_or_retraining_on_GSE301164": False,
        "predeclared_contrasts": {
            "replicative_senescence_P33_vs_P3": replicative,
            "doxorubicin_induced_senescence_DOX_vs_DMSO_at_P11": doxorubicin,
        },
        "combined": {
            "pair_count": len(all_differences),
            "positive_pair_count": sum(value > 0 for value in all_differences),
            "exact_one_sided_sign_test_p": combined_p,
            "directional_effect_evidence_eligible": directional_evidence,
        },
        "ood": {
            "threshold_from_preexisting_GSE63577_calibration": ood_threshold,
            "rejection_fraction": float(np.mean(distance > ood_threshold)),
            "sample_distance": {name: float(value) for name, value in zip(sample_names, distance)},
        },
        "absolute_state_classifier_promoted": False,
        "may_select_aleph_parameter": False,
        "aleph_authority": "none",
        "eventual_sweep_role": "directional_transcriptomic_observable_evidence_pending_exact_Aleph_observation_operator_and_forward_sensitivity",
        "limitations": [
            "both contrasts use the BJ fibroblast line in one confirmation laboratory",
            "directional replication does not calibrate absolute senescence probabilities across laboratories",
            "OOD rejection remains a hard refusal for sample-level state use",
            "no Aleph parameter or numeric sweep range may be emitted before operator and sensitivity validation",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.counts, args.checkpoint)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["combined"], sort_keys=True))


if __name__ == "__main__":
    main()

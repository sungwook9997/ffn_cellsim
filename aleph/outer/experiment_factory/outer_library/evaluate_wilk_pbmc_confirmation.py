#!/usr/bin/env python3
"""Run the frozen single-use Wilk PBMC confirmation without refitting."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from scipy import sparse

from train_pbmc_cross_lab_model import (
    CELL_CLASSES,
    auroc,
    confusion_metrics,
    conformal_metrics,
    ece,
    softmax_from_logits,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def text_array(value: h5py.Dataset | h5py.Group) -> np.ndarray:
    if isinstance(value, h5py.Group):
        categories = text_array(value["categories"])
        codes = np.asarray(value["codes"], dtype=np.int64)
        if np.any(codes < 0):
            result = np.full(len(codes), "", dtype=object)
            valid = codes >= 0
            result[valid] = categories[codes[valid]]
            return result.astype(str)
        return categories[codes].astype(str)
    raw = np.asarray(value)
    return np.asarray([
        item.decode("utf-8") if isinstance(item, bytes) else str(item)
        for item in raw
    ])


def normalized_axis_name(value: str) -> str:
    normalized = re.sub(r"[^a-z0-9]", "", value.casefold())
    return normalized[:-1] if normalized.endswith("s") else normalized


def resolve_column(obs: h5py.Group, aliases: list[str]) -> str:
    by_normalized: dict[str, list[str]] = {}
    for name in obs.keys():
        by_normalized.setdefault(normalized_axis_name(name), []).append(name)
    for alias in aliases:
        matches = by_normalized.get(normalized_axis_name(alias), [])
        if len(matches) > 1:
            raise ValueError(f"ambiguous H5AD observation column for {alias}: {matches}")
        if matches:
            return matches[0]
    raise ValueError(f"no H5AD observation column matches frozen aliases {aliases}")


def aliases_from_contract(contract: str) -> list[str]:
    fragment = contract.split("among ", 1)[1].split(";", 1)[0]
    return [value.strip() for value in fragment.split(",")]


def map_label(label: str, protocol: dict[str, Any]) -> tuple[str | None, str]:
    lowered = label.casefold()
    generic_t = re.search(r"(^|[^a-z])t([^a-z]|$)", lowered) is not None and not any(
        token in lowered
        for token in ("cd4", "cd8", "helper", "regulatory", "treg", "cytotoxic")
    )
    if generic_t:
        return None, "excluded_generic_T"
    if any(token in lowered for token in ("unassigned", "unknown")):
        return None, "excluded_unknown"
    for item in protocol["coarse_cell_type_mapping_in_order"]:
        if any(
            re.search(pattern, lowered, flags=re.IGNORECASE)
            for pattern in item["case_insensitive_patterns"]
        ):
            return item["canonical"], "in_domain"
    return None, "semantic_OOD"


def encode(labels: np.ndarray) -> np.ndarray:
    index = {str(label): position for position, label in enumerate(CELL_CLASSES)}
    return np.asarray([index[str(label)] for label in labels], dtype=np.int64)


def matrix_shape(group: h5py.Group) -> tuple[int, int]:
    shape = np.asarray(group.attrs["shape"], dtype=np.int64).ravel()
    if len(shape) != 2:
        raise ValueError("H5AD sparse matrix has no two-dimensional shape")
    return int(shape[0]), int(shape[1])


def read_csr(group: h5py.Group) -> sparse.csr_matrix:
    if group.attrs.get("encoding-type") != "csr_matrix":
        raise ValueError("frozen raw-count precedence resolved to a non-CSR matrix")
    data = np.asarray(group["data"])
    if np.any(data < 0) or np.any(data != np.floor(data)):
        raise ValueError("selected H5AD matrix is not nonnegative integer-valued counts")
    return sparse.csr_matrix(
        (data, np.asarray(group["indices"]), np.asarray(group["indptr"])),
        shape=matrix_shape(group),
    )


def per_class_coverage(
    y: np.ndarray, probability: np.ndarray, thresholds: np.ndarray
) -> dict[str, dict[str, float | int]]:
    included = (1.0 - probability) <= thresholds[None, :]
    result = {}
    for class_index, label in enumerate(CELL_CLASSES):
        mask = y == class_index
        result[str(label)] = {
            "support": int(mask.sum()),
            "coverage": float(np.mean(included[mask, class_index])) if np.any(mask) else 0.0,
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--calibration-report", type=Path, required=True)
    parser.add_argument("--output-report", type=Path, required=True)
    parser.add_argument("--output-tensor", type=Path, required=True)
    args = parser.parse_args()

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    receipt = json.loads(args.receipt.read_text(encoding="utf-8"))
    calibration_report = json.loads(args.calibration_report.read_text(encoding="utf-8"))
    if receipt["asset"]["size_bytes"] != args.h5ad.stat().st_size:
        raise ValueError("Wilk H5AD size differs from acquisition receipt")
    if receipt["asset"]["sha256"] != sha256(args.h5ad):
        raise ValueError("Wilk H5AD hash differs from acquisition receipt")
    if receipt["protocol_sha256"] != sha256(args.protocol):
        raise ValueError("Wilk receipt and frozen protocol disagree")
    if calibration_report["output_model_sha256"] != sha256(args.model):
        raise ValueError("calibrated model differs from its calibration report")

    model = np.load(args.model, allow_pickle=False)
    feature_symbols = model["feature_symbols"].astype(str)
    parameters = model["cell_parameters"]
    if array_sha256(parameters) != str(model["source_parameter_sha256"]):
        raise ValueError("calibrated artifact changed frozen classifier parameters")
    if array_sha256(parameters) != calibration_report["source_parameter_sha256"]:
        raise ValueError("calibrated artifact and report parameter hashes disagree")
    classes = model["cell_classes"].astype(str)
    if not np.array_equal(classes, CELL_CLASSES):
        raise ValueError("calibrated model class order changed")
    temperature = float(model["cell_temperature"])
    thresholds = model["cell_conformal_thresholds"].astype(np.float64)

    with h5py.File(args.h5ad, "r") as handle:
        obs = handle["obs"]
        contract = protocol["h5ad_axis_contract"]
        cell_type_column = resolve_column(
            obs, aliases_from_contract(contract["cell_type_column_resolution"])
        )
        donor_column = resolve_column(
            obs, aliases_from_contract(contract["donor_column_resolution"])
        )
        condition_column = resolve_column(
            obs, aliases_from_contract(contract["condition_column_resolution"])
        )
        cell_ids = text_array(obs["_index"])
        if len(set(cell_ids.tolist())) != len(cell_ids):
            raise ValueError("duplicate H5AD observation identifiers")
        author_labels = text_array(obs[cell_type_column])
        donors = text_array(obs[donor_column])
        conditions = text_array(obs[condition_column])
        healthy = np.asarray([
            value.casefold().strip() in {"normal", "healthy control"}
            for value in conditions
        ])
        healthy_donors = sorted(set(donors[healthy].tolist()))
        if len(healthy_donors) != protocol["population_contract"]["required_distinct_donors"]:
            raise ValueError(f"expected six healthy donors, found {healthy_donors}")

        raw_group = handle["raw/X"]
        source_symbols = text_array(handle["raw/var/feature_name"])
        source_by_symbol: dict[str, int] = {}
        duplicate_symbols = Counter(source_symbols.tolist())
        for index, symbol in enumerate(source_symbols):
            source_by_symbol.setdefault(symbol, index)
        model_positions = np.asarray([
            index for index, symbol in enumerate(feature_symbols) if symbol in source_by_symbol
        ], dtype=np.int64)
        source_positions = np.asarray([
            source_by_symbol[feature_symbols[index]] for index in model_positions
        ], dtype=np.int64)
        if len(source_positions) < contract["minimum_common_features"]:
            raise ValueError(f"only {len(source_positions)} frozen model features occur in Wilk")
        counts = read_csr(raw_group)
        if counts.shape != (len(cell_ids), len(source_symbols)):
            raise ValueError("H5AD raw matrix axes disagree with obs/raw-var axes")

    healthy_indices = np.flatnonzero(healthy)
    healthy_counts = counts[healthy_indices]
    library_size = np.asarray(healthy_counts.sum(axis=1)).ravel().astype(np.float64)
    if np.any(library_size < protocol["expression_contract"]["minimum_library_size"]):
        raise ValueError("healthy confirmation population contains a zero-library cell")
    selected = healthy_counts[:, source_positions].tocsr().astype(np.float64)
    selected.data *= np.repeat(10000.0 / library_size, np.diff(selected.indptr))
    np.log1p(selected.data, out=selected.data)
    feature_count = len(feature_symbols)
    class_count = len(CELL_CLASSES)
    weights = parameters[: feature_count * (class_count - 1)].reshape(
        feature_count, class_count - 1
    )
    intercept = parameters[feature_count * (class_count - 1):]
    nonreference_logits = selected @ weights[model_positions, :] + intercept
    logits = np.column_stack((
        np.asarray(nonreference_logits), np.zeros(selected.shape[0], dtype=np.float64)
    ))
    probability = softmax_from_logits(logits / temperature)

    mapped = [map_label(label, protocol) for label in author_labels[healthy]]
    canonical = np.asarray([item[0] or "" for item in mapped])
    roles = np.asarray([item[1] for item in mapped])
    in_domain = roles == "in_domain"
    semantic_ood = roles == "semantic_OOD"
    y = encode(canonical[in_domain])
    in_probability = probability[in_domain]
    metrics = confusion_metrics(y, in_probability, CELL_CLASSES)
    metrics["ece_15_bin"] = ece(y, in_probability)
    metrics["conformal"] = conformal_metrics(y, in_probability, thresholds)
    metrics["class_conditional_conformal"] = per_class_coverage(
        y, in_probability, thresholds
    )

    healthy_donor_array = donors[healthy]
    donor_metrics = {}
    for donor in healthy_donors:
        mask = in_domain & (healthy_donor_array == donor)
        donor_y = encode(canonical[mask])
        item = confusion_metrics(donor_y, probability[mask], CELL_CLASSES)
        item["ece_15_bin"] = ece(donor_y, probability[mask])
        item["conformal"] = conformal_metrics(donor_y, probability[mask], thresholds)
        item["class_conditional_conformal"] = per_class_coverage(
            donor_y, probability[mask], thresholds
        )
        donor_metrics[donor] = item

    ood_eligible = int(semantic_ood.sum()) >= 100 and int(in_domain.sum()) >= 100
    ood_report: dict[str, Any] = {
        "score": "one_minus_maximum_calibrated_probability",
        "positive_cells": int(semantic_ood.sum()),
        "negative_cells": int(in_domain.sum()),
        "eligible": ood_eligible,
        "positive_author_labels": dict(Counter(author_labels[healthy][semantic_ood].tolist())),
    }
    if ood_eligible:
        ood_mask = in_domain | semantic_ood
        ood_report["auroc"] = auroc(
            semantic_ood[ood_mask], 1.0 - probability[ood_mask].max(axis=1)
        )

    endpoints = protocol["single_use_endpoints"]
    conditional_checks = {
        label: bool(
            item["support"] < 100
            or item["coverage"] >= endpoints["class_conditional_coverage_minimum_for_classes_with_at_least_100_cells"]
        )
        for label, item in metrics["class_conditional_conformal"].items()
    }
    gates = {
        "macro_f1": metrics["macro_f1"] >= endpoints["macro_f1_minimum"],
        "ece_15_bin": metrics["ece_15_bin"] <= endpoints["ece_15_bin_maximum"],
        "marginal_conformal_coverage": metrics["conformal"]["marginal_coverage"]
        >= endpoints["marginal_conformal_coverage_minimum"],
        "class_conditional_conformal_coverage": all(conditional_checks.values()),
        "mean_prediction_set_size": metrics["conformal"]["mean_set_size"]
        <= endpoints["mean_prediction_set_size_maximum"],
        "semantic_ood_auroc": (
            not ood_eligible
            or ood_report["auroc"] >= endpoints["semantic_ood_auroc_minimum_if_at_least_100_positive_and_100_negative_cells"]
        ),
    }
    passed = all(gates.values())

    args.output_tensor.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_tensor,
        healthy_cell_id=cell_ids[healthy],
        healthy_donor=healthy_donor_array,
        author_label=author_labels[healthy],
        mapping_role=roles,
        canonical_label=canonical,
        calibrated_probability=probability.astype(np.float32),
    )
    report = {
        "schema": "aleph.outer_library.wilk_pbmc_single_use_confirmation.v1",
        "status": "passed" if passed else "failed_and_permanently_retired",
        "single_use_confirmation_consumed": True,
        "may_be_reused_as_pristine": False,
        "invalid_pre_metric_implementation_attempt": {
            "occurred": True,
            "authoritative_result": False,
            "reason": "the first evaluator revision implemented generic T as the character t, incorrectly excluding monocyte, natural killer, and dendritic labels; no model, calibration, mapping contract, or endpoint changed before this contract-faithful rerun",
        },
        "protocol_sha256": sha256(args.protocol),
        "acquisition_receipt_sha256": sha256(args.receipt),
        "h5ad_sha256": sha256(args.h5ad),
        "calibrated_model_sha256": sha256(args.model),
        "classifier_parameter_sha256": array_sha256(parameters),
        "classifier_weights_bit_identical": True,
        "calibration_parameters_changed_during_confirmation": False,
        "resolved_axes": {
            "cell_type_column": cell_type_column,
            "donor_column": donor_column,
            "condition_column": condition_column,
            "raw_count_matrix": "raw/X",
            "gene_symbol_axis": "raw/var/feature_name",
        },
        "population": {
            "total_cells": len(cell_ids),
            "healthy_cells": int(healthy.sum()),
            "healthy_donors": healthy_donors,
            "in_domain_cells": int(in_domain.sum()),
            "semantic_ood_cells": int(semantic_ood.sum()),
            "excluded_cells": int((~in_domain & ~semantic_ood).sum()),
            "mapping_role_counts": dict(Counter(roles.tolist())),
            "author_label_counts": dict(Counter(author_labels[healthy].tolist())),
            "common_model_features": len(source_positions),
            "duplicate_source_symbols_using_first_row": sum(
                value > 1 for value in duplicate_symbols.values()
            ),
        },
        "primary_metrics": metrics,
        "metrics_by_healthy_donor": donor_metrics,
        "semantic_ood": ood_report,
        "frozen_endpoint_checks": gates,
        "class_conditional_checks": conditional_checks,
        "all_frozen_endpoints_passed": passed,
        "output_tensor": {
            "path": args.output_tensor.name,
            "sha256": sha256(args.output_tensor),
        },
        "sweep_connection": {
            "current_development_purpose": "pre_sweep_context_compression_because_broad_Aleph_sweeps_are_not_yet_executable",
            "eventual_program_objective": "bounded_evidence_conditioned_Aleph_parameter_sweeps",
            "context_ranking_eligible": passed,
            "may_delete_sweep_axis": False,
            "may_emit_numeric_sweep_range": False,
            "may_select_aleph_parameter": False,
        },
    }
    args.output_report.parent.mkdir(parents=True, exist_ok=True)
    args.output_report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": report["status"],
        "healthy_donors": len(healthy_donors),
        "in_domain_cells": int(in_domain.sum()),
        "macro_f1": metrics["macro_f1"],
        "ece": metrics["ece_15_bin"],
        "coverage": metrics["conformal"]["marginal_coverage"],
        "mean_set_size": metrics["conformal"]["mean_set_size"],
        "ood_auroc": ood_report.get("auroc"),
    }, sort_keys=True))


if __name__ == "__main__":
    main()

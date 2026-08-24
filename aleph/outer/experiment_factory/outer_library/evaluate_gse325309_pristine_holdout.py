#!/usr/bin/env python3
"""Run the frozen exact GSE325309 external EMT direction test."""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np


def _score(X: np.ndarray, treated: np.ndarray) -> tuple[float, np.ndarray]:
    control = np.asarray([i for i in range(6) if i not in set(treated)], dtype=int)
    deltas = X[treated].mean(axis=0) - X[control].mean(axis=0)
    pooled = np.sqrt((X[treated].var(axis=0, ddof=1) + X[control].var(axis=0, ddof=1)) / 2)
    if np.any(pooled == 0):
        raise ValueError("GSE325309 frozen standardized effect has zero pooled SD")
    d = deltas / pooled
    return float(-d[0] + d[1]), deltas


def evaluate(tensor_path: Path, manifest_report_path: Path) -> dict[str, object]:
    manifest = json.loads(manifest_report_path.read_text())
    if manifest.get("schema") != "aleph.outer_library.gse325309_pristine_manifest.v1":
        raise ValueError("GSE325309 manifest schema changed")
    with np.load(tensor_path, allow_pickle=False) as data:
        X = np.log2(data["X_FPKM"].astype(np.float64) + 1.0)
        genes = data["gene"].astype(str)
        sample_ids = data["sample_id"].astype(str)
    if X.shape != (6, 2) or list(genes) != ["CDH1", "VIM"]:
        raise ValueError("GSE325309 frozen tensor changed")
    observed_treated = np.asarray([3, 4, 5], dtype=int)
    observed_score, deltas = _score(X, observed_treated)
    null_scores = np.asarray([
        _score(X, np.asarray(indices, dtype=int))[0]
        for indices in combinations(range(6), 3)
    ])
    exact_p = float(np.count_nonzero(null_scores >= observed_score - 1e-12) / len(null_scores))
    directions = {"CDH1": "decrease" if deltas[0] < 0 else "increase",
                  "VIM": "increase" if deltas[1] > 0 else "decrease"}
    pairwise = {
        gene: bool(np.all(
            (X[3:, index, None] - X[:3, index][None, :]) < 0
            if gene == "CDH1" else
            (X[3:, index, None] - X[:3, index][None, :]) > 0
        )) for index, gene in enumerate(genes)
    }
    rng = np.random.default_rng(325309)
    bootstrap = {}
    for index, gene in enumerate(genes):
        values = np.empty(10_000)
        for draw in range(len(values)):
            control = rng.choice(X[:3, index], 3, replace=True)
            treated = rng.choice(X[3:, index], 3, replace=True)
            values[draw] = treated.mean() - control.mean()
        bootstrap[gene] = {
            "log2_FPKM_plus_1_mean_difference": float(deltas[index]),
            "bootstrap_percentile_95": [float(x) for x in np.quantile(values, [0.025, 0.975])],
        }
    primary_passed = bool(deltas[0] < 0 and deltas[1] > 0 and exact_p <= 0.05)
    return {
        "schema": "aleph.outer_library.gse325309_pristine_evaluation.v1",
        "task": "outcome_blind_preregistered_external_EMT_direction_holdout",
        "freeze_proof": manifest["source"]["freeze_proof"],
        "design_audit": {
            **manifest["counts"], "sample_exclusions": manifest["sample_exclusions"],
            "biological_replicate_is_test_unit": True,
            "protocol_frozen_before_expression_access": True,
        },
        "raw_FPKM": {
            gene: {"control": [float(x) for x in (2 ** X[:3, i] - 1)],
                   "TGFB1": [float(x) for x in (2 ** X[3:, i] - 1)]}
            for i, gene in enumerate(genes)
        },
        "frozen_primary_test": {
            "directions": directions, "oriented_composite_score": observed_score,
            "exact_assignment_count": len(null_scores), "exact_one_sided_p": exact_p,
            "alpha": 0.05, "passed": primary_passed,
        },
        "secondary_diagnostics": {
            "all_nine_pairwise_differences_expected": pairwise,
            "bootstrap": bootstrap,
        },
        "pristine_external_direction_holdout_passed": primary_passed,
        "numeric_uncertainty_promotion_passed": False,
        "OOD_refusal_validated": False,
        "pre_sweep_role": {
            "supported_action": "confirm retrospective CDH1-down/VIM-up filter on an outcome-blind external lab",
            "numeric_range": None, "may_emit_sweep_axis": False,
            "remaining_blockers": [
                "Aleph CDH1/VIM-compatible observation operators",
                "monotonic Aleph forward sensitivity",
                "multiple independently frozen semantic OOD contexts",
                "numeric cross-assay calibration rather than direction-only evidence",
            ],
        },
        "aleph_authority": "none", "may_select_aleph_parameter": False,
        "sample_ids": list(sample_ids),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--manifest-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.tensor, args.manifest_report)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["frozen_primary_test"], sort_keys=True))


if __name__ == "__main__":
    main()

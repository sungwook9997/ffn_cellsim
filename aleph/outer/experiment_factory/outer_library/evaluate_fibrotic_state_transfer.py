#!/usr/bin/env python3
"""Cross-lab TGF-beta fibrotic-state contrast from IF to bulk RNA-seq."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def _signature_scores(path: Path) -> tuple[dict[str, float], dict[str, str]]:
    samples: dict[str, dict[str, float]] = defaultdict(dict)
    conditions: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            gene = row["observable"].removesuffix("_normalized_expression")
            samples[row["sample_id"]][gene] = float(row["value"])
            conditions[row["sample_id"]] = row["condition"]
    controls = [sample for sample, condition in conditions.items()
                if condition == "vehicle_control"]
    genes = sorted(next(iter(samples.values())))
    control_reference = {
        gene: np.mean([np.log1p(samples[sample][gene]) for sample in controls])
        for gene in genes
    }
    scores = {
        sample: float(np.mean([
            np.log1p(values[gene]) - control_reference[gene] for gene in genes
        ]))
        for sample, values in samples.items()
    }
    return scores, conditions


def _bootstrap_difference(
    treated: np.ndarray,
    control: np.ndarray,
    rng: np.random.Generator,
    draws: int,
) -> dict:
    effects = np.empty(draws, dtype=np.float64)
    for index in range(draws):
        effects[index] = (
            treated[rng.integers(len(treated), size=len(treated))].mean()
            - control[rng.integers(len(control), size=len(control))].mean()
        )
    return {
        "treated_n": len(treated),
        "control_n": len(control),
        "mean_difference": float(treated.mean() - control.mean()),
        "bootstrap_ci95": [
            float(value) for value in np.quantile(effects, (0.025, 0.975))
        ],
        "bootstrap_probability_positive": float(np.mean(effects > 0)),
    }


def evaluate(
    tendon_report_path: Path,
    external_manifest: Path,
    draws: int = 20_000,
) -> dict:
    tendon = json.loads(tendon_report_path.read_text(encoding="utf-8"))
    scores, conditions = _signature_scores(external_manifest)
    by_condition = {
        condition: np.asarray([
            score for sample, score in scores.items()
            if conditions[sample] == condition
        ], dtype=np.float64)
        for condition in sorted(set(conditions.values()))
    }
    rng = np.random.default_rng(20260805)
    induction = _bootstrap_difference(
        by_condition["TGF_beta_1"], by_condition["vehicle_control"], rng, draws
    )
    reversal = _bootstrap_difference(
        by_condition["TGF_beta_1"],
        by_condition["TGF_beta_1_plus_celastrol"],
        rng,
        draws,
    )
    discovery = tendon["TGF_beta_1_contrasts"]["alpha_SMA_relative_fluorescence"]
    passed = (
        discovery["bootstrap_ci95"][0] > 0
        and induction["bootstrap_ci95"][0] > 0
        and reversal["bootstrap_ci95"][0] > 0
    )
    return {
        "task": "cross_lab_TGF_beta_1_fibrotic_state_contrast",
        "discovery": {
            "dataset": tendon["dataset"],
            "lab_group": tendon["lab_group"],
            "species_cell_type": "rat tail tendon stromal cells",
            "modality": "immunofluorescence_derived",
            "observable": "alpha_SMA_relative_fluorescence",
            "TGF_beta_1_minus_control": discovery,
        },
        "external_confirmatory": {
            "dataset": "GSE226374",
            "lab_group": "leask-lab-gse226374",
            "species_cell_type": "CCD-1077Sk human dermal fibroblasts",
            "modality": "bulk_RNA_seq_derived",
            "predeclared_marker_panel": ["ACTA2", "CCN2", "COL1A1", "FN1", "THBS1"],
            "representation": "mean log1p expression relative to external-study vehicle controls",
            "TGF_beta_1_minus_control": induction,
            "celastrol_reversal_TGF_minus_TGF_plus_celastrol": reversal,
            "sample_scores": scores,
        },
        "independent_dataset_holdout": True,
        "independent_lab_holdout": True,
        "independent_lab_mechanism_holdout_passed": passed,
        "cell_state_classifier_holdout_passed": False,
        "observable_constraint_eligible": passed,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": [
            "the confirmatory endpoint is a four-condition, two-replicate bulk-RNA contrast",
            "species, fibroblast source, and modality differ across discovery and confirmation",
            "the marker panel validates a fibrotic state direction, not a per-cell classifier",
            "celastrol reversal supports pathway specificity but does not identify one Aleph parameter",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tendon-report", type=Path, required=True)
    parser.add_argument("--external-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-draws", type=int, default=20_000)
    args = parser.parse_args()
    report = evaluate(
        args.tendon_report, args.external_manifest, args.bootstrap_draws
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "independent_lab_mechanism_holdout_passed": report[
            "independent_lab_mechanism_holdout_passed"
        ],
        "external": report["external_confirmatory"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

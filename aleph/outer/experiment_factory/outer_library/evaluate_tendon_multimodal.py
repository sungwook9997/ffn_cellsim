#!/usr/bin/env python3
"""Evaluate cross-modal tendon mechanics without inventing a lab holdout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _load(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle]


def _values(
    rows: list[dict], observable: str, cell_type: str, condition: str
) -> np.ndarray:
    return np.asarray([
        row["value"] for row in rows
        if row["observable"] == observable
        and row["cell_type"] == cell_type
        and row["condition"] == condition
    ], dtype=np.float64)


def _contrast(
    treated: np.ndarray,
    control: np.ndarray,
    rng: np.random.Generator,
    draws: int,
    *,
    inferential_unit_known: bool = True,
) -> dict:
    if not len(treated) or not len(control):
        raise ValueError("empty contrast arm")
    difference = float(treated.mean() - control.mean())
    result = {
        "treated_n": int(len(treated)),
        "control_n": int(len(control)),
        "treated_mean": float(treated.mean()),
        "control_mean": float(control.mean()),
        "mean_difference": difference,
        "inferential_unit_known": inferential_unit_known,
    }
    if not inferential_unit_known:
        result.update({
            "bootstrap_ci95": None,
            "direction_probability": None,
            "evidence_eligible": False,
            "reason": "author biological-replicate identity is unavailable",
        })
        return result
    effects = np.empty(draws, dtype=np.float64)
    for index in range(draws):
        effects[index] = (
            treated[rng.integers(len(treated), size=len(treated))].mean()
            - control[rng.integers(len(control), size=len(control))].mean()
        )
    result.update({
        "bootstrap_ci95": [
            float(value) for value in np.quantile(effects, (0.025, 0.975))
        ],
        "direction_probability": float(
            np.mean(effects > 0) if difference >= 0 else np.mean(effects < 0)
        ),
        "evidence_eligible": True,
    })
    return result


def evaluate(path: Path, draws: int = 20_000) -> dict:
    rows = _load(path)
    rng = np.random.default_rng(20260805)
    rat = "rat tail tendon stromal cells"

    def contrast(observable: str, treated: str, control: str, *, known: bool = True) -> dict:
        return _contrast(
            _values(rows, observable, rat, treated),
            _values(rows, observable, rat, control),
            rng, draws, inferential_unit_known=known,
        )

    matrix_tension = {
        "tissue_traction_force_per_cell": contrast(
            "tissue_traction_force_per_cell", "rigid_boundary", "compliant_boundary"
        ),
        "alpha_SMA_relative_fluorescence": contrast(
            "alpha_SMA_relative_fluorescence", "rigid_boundary", "compliant_boundary"
        ),
        "F_actin_relative_fluorescence": contrast(
            "F_actin_relative_fluorescence", "rigid_boundary", "compliant_boundary"
        ),
        "Acta2_relative_expression": contrast(
            "Acta2_relative_expression", "rigid_boundary", "compliant_boundary"
        ),
        "nuclear_circularity": contrast(
            "nuclear_circularity", "rigid_boundary", "compliant_boundary", known=False
        ),
    }
    tgf = {
        observable: contrast(observable, "TGF_beta_1", "vehicle_control", known=known)
        for observable, known in (
            ("alpha_SMA_relative_fluorescence", True),
            ("F_actin_relative_fluorescence", True),
            ("nuclear_circularity", False),
        )
    }

    disease: dict[str, dict] = {}
    disease_effects = []
    for cell_type in ("human tendon fibroblasts", "human dermal fibroblasts"):
        treated = _values(rows, "tissue_traction_force_per_cell", cell_type,
                          "systemic_sclerosis")
        control = _values(rows, "tissue_traction_force_per_cell", cell_type,
                          "healthy_control")
        item = _contrast(treated, control, rng, draws)
        disease[cell_type] = item
        effects = np.empty(draws, dtype=np.float64)
        for index in range(draws):
            effects[index] = (
                treated[rng.integers(len(treated), size=len(treated))].mean()
                - control[rng.integers(len(control), size=len(control))].mean()
            )
        disease_effects.append(effects)
    equal_cell_type_effect = np.stack(disease_effects).mean(axis=0)
    disease_ci = [
        float(value) for value in np.quantile(equal_cell_type_effect, (0.025, 0.975))
    ]

    primary_expected_positive = (
        "tissue_traction_force_per_cell",
        "alpha_SMA_relative_fluorescence",
        "Acta2_relative_expression",
    )
    primary_pass = all(
        matrix_tension[name]["bootstrap_ci95"][0] > 0
        for name in primary_expected_positive
    )
    cell_type_holdout_pass = (
        all(item["mean_difference"] > 0 for item in disease.values())
        and disease_ci[0] > 0
    )
    return {
        "task": "multimodal_tendon_mechanical_state_evidence",
        "dataset": "nature-2026-tendon-mechanoculture",
        "lab_group": "hussien-snedeker-2026",
        "observation_count": len(rows),
        "modalities": sorted({row["modality"] for row in rows}),
        "matrix_tension_contrasts": matrix_tension,
        "TGF_beta_1_contrasts": tgf,
        "human_fibrosis_force": {
            "discovery_cell_type": "human tendon fibroblasts",
            "held_out_cell_type": "human dermal fibroblasts",
            "same_lab": True,
            "independent_lab_holdout": False,
            "contrasts": disease,
            "equal_cell_type_mean_difference": float(equal_cell_type_effect.mean()),
            "equal_cell_type_bootstrap_ci95": disease_ci,
            "cell_type_holdout_passed": cell_type_holdout_pass,
        },
        "multimodal_primary_coherence_passed": primary_pass,
        "observable_constraint_eligible": primary_pass and cell_type_holdout_pass,
        "independent_lab_task_holdout_passed": False,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": [
            "tendon and dermal fibroblast contrasts are a cell-type holdout in one lab, not a lab holdout",
            "nuclear measurements lack biological-replicate IDs and remain descriptive",
            "source-relative qPCR units are preserved without reverse-engineering normalization",
            "observable evidence does not identify a unique Aleph parameter or numeric range",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-draws", type=int, default=20_000)
    args = parser.parse_args()
    report = evaluate(args.manifest, args.bootstrap_draws)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "multimodal_primary_coherence_passed": report["multimodal_primary_coherence_passed"],
        "cell_type_holdout_passed": report["human_fibrosis_force"]["cell_type_holdout_passed"],
        "independent_lab_task_holdout_passed": False,
    }, sort_keys=True))


if __name__ == "__main__":
    main()

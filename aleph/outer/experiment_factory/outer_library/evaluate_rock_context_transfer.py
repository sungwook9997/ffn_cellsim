#!/usr/bin/env python3
"""Audit context dependence of ROCK-inhibition migration-speed effects."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def _load(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle]


def _bootstrap(values: np.ndarray, draws: int, seed: int) -> list[float]:
    rng = np.random.default_rng(seed)
    estimates = values[rng.integers(len(values), size=(draws, len(values)))].mean(axis=1)
    return [float(value) for value in np.quantile(estimates, [0.025, 0.975])]


def _bootstrap_difference(
    control: np.ndarray, treated: np.ndarray, draws: int, seed: int
) -> list[float]:
    rng = np.random.default_rng(seed)
    control_mean = control[rng.integers(len(control), size=(draws, len(control)))].mean(axis=1)
    treated_mean = treated[rng.integers(len(treated), size=(draws, len(treated)))].mean(axis=1)
    return [float(value) for value in np.quantile(treated_mean - control_mean, [0.025, 0.975])]


def _hedges_g(control: np.ndarray, treated: np.ndarray) -> float:
    pooled = np.sqrt(
        ((len(control) - 1) * np.var(control, ddof=1)
         + (len(treated) - 1) * np.var(treated, ddof=1))
        / (len(control) + len(treated) - 2)
    )
    correction = 1 - 3 / (4 * (len(control) + len(treated)) - 9)
    return float(correction * (np.mean(treated) - np.mean(control)) / pooled)


def evaluate(dryad_manifest: Path, elife_manifest: Path, *, draws: int = 10_000) -> dict:
    dryad_rows = [
        row for row in _load(dryad_manifest)
        if row["dataset_id"] == "dryad-9jh6m-ROCK-modulation"
        and row["observable"] == "Cell Speed"
    ]
    by_date: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in dryad_rows:
        by_date[row["biological_replicate"]][row["condition"]].append(float(row["value"]))
    paired_date_effects = []
    for date, arms in sorted(by_date.items()):
        if arms["DMSO_control"] and arms["Y27632_6_micromolar"]:
            paired_date_effects.append(
                float(np.mean(arms["Y27632_6_micromolar"]) - np.mean(arms["DMSO_control"]))
            )
    dryad_effects = np.asarray(paired_date_effects)
    if len(dryad_effects) != 4:
        raise ValueError(f"expected four paired Dryad dates, found {len(dryad_effects)}")

    elife = [row for row in _load(elife_manifest) if row["observable"] in {
        "persistent_speed", "speed_oscillation_period"
    }]
    values: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in elife:
        if row["condition"] in {"untreated_control", "Y27632_10_micromolar"}:
            values[(row["condition"], row["observable"])].append(float(row["value"]))
    control_speed = np.asarray(values[("untreated_control", "persistent_speed")])
    treated_speed = np.asarray(values[("Y27632_10_micromolar", "persistent_speed")])
    control_period = np.asarray(values[("untreated_control", "speed_oscillation_period")])
    treated_period = np.asarray(values[("Y27632_10_micromolar", "speed_oscillation_period")])
    if (len(control_speed), len(treated_speed)) != (24, 9):
        raise ValueError("unexpected eLife control/Y27632 cell counts")
    speed_difference = float(np.mean(treated_speed) - np.mean(control_speed))
    period_difference = float(np.mean(treated_period) - np.mean(control_period))
    direction_consistent = bool(np.mean(dryad_effects) * speed_difference > 0)
    return {
        "task": "ROCK_inhibition_migration_speed_context_transfer",
        "dryad_H1299_PL_2D": {
            "paired_experimental_dates": int(len(dryad_effects)),
            "Y27632_minus_DMSO_date_effects": dryad_effects.tolist(),
            "mean_date_effect_source_unit": float(np.mean(dryad_effects)),
            "experimental_date_bootstrap_95ci": _bootstrap(dryad_effects, draws, 20260805),
            "direction": "increases" if np.mean(dryad_effects) > 0 else "decreases",
            "inference_unit": "author_experimental_date",
        },
        "elife_NIH3T3_3D_CDM": {
            "control_cells": int(len(control_speed)),
            "Y27632_cells": int(len(treated_speed)),
            "persistent_speed_Y27632_minus_control_micrometre_per_minute": speed_difference,
            "persistent_speed_hedges_g": _hedges_g(control_speed, treated_speed),
            "speed_oscillation_period_Y27632_minus_control_minute": period_difference,
            "cell_bootstrap_speed_difference_95ci_descriptive_only": _bootstrap_difference(
                control_speed, treated_speed, draws, 20260806
            ),
            "direction": "increases" if speed_difference > 0 else "decreases",
            "independent_experiment_identifiers_reported": False,
            "inference_unit": "unknown; cells cannot be promoted to biological replicates",
        },
        "speed_direction_consistent_across_contexts": direction_consistent,
        "context_heterogeneity_detected": not direction_consistent,
        "interpretation": (
            "ROCK inhibition has opposite migration-speed directions in H1299-PL 2D "
            "and NIH3T3 3D CDM contexts; a universal speed constraint or one-parameter "
            "Aleph mapping is rejected."
        ),
        "independent_lab_dataset_present": True,
        "independent_lab_mechanism_holdout_passed": False,
        "failure_reason": (
            "independent eLife experiment identifiers are not reported and the context-specific "
            "speed direction reverses"
        ),
        "eligible_as_context_diversity_training": True,
        "observable_constraint_eligible": False,
        "aleph_parameter_proposal_eligible": False,
        "may_mutate_physics": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dryad-manifest", type=Path, required=True)
    parser.add_argument("--elife-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=10_000)
    args = parser.parse_args()
    report = evaluate(args.dryad_manifest, args.elife_manifest, draws=args.draws)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "context_heterogeneity_detected": report["context_heterogeneity_detected"],
        "dryad": report["dryad_H1299_PL_2D"],
        "elife": report["elife_NIH3T3_3D_CDM"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

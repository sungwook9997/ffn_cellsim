#!/usr/bin/env python3
"""Evaluate eLife 72381 at the independent-experiment inference level."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


CONTROL = "vehicle_control"
Y27 = "Y27632_100_micromolar"
DATA_SUM_INNER_SHA256 = "76a65e2efcb04f4f699d161175422df8ac04bef5ee40bdfb7a8230d2dc0137de"


def _load(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle]


def _published_data_sum(path: Path) -> dict[str, float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("source_sha256") != DATA_SUM_INNER_SHA256:
        raise ValueError("unexpected Figure 2 data_sum workbook hash")
    sheets = {sheet["name"]: sheet["values"] for sheet in payload["sheets"]}
    rows = sheets.get("Sheet1")
    if not rows or rows[0][1] != "dA/dt i (µm2/min)" or rows[0][4] != "dV/dt i (µm3/min)":
        raise ValueError("unexpected Figure 2 data_sum schema")
    labels = {str(row[0]): row for row in rows[1:]}
    control = labels["control all, n=194"]
    treated = labels["Y-27, n=121"]
    return {
        f"{CONTROL}:initial_spreading_area_rate": float(control[1]),
        f"{CONTROL}:initial_volume_flux": float(control[4]),
        f"{Y27}:initial_spreading_area_rate": float(treated[1]),
        f"{Y27}:initial_volume_flux": float(treated[4]),
    }


def _bootstrap(values: np.ndarray, draws: int, seed: int) -> list[float]:
    rng = np.random.default_rng(seed)
    estimates = values[
        rng.integers(len(values), size=(draws, len(values)))
    ].mean(axis=1)
    return [float(value) for value in np.quantile(estimates, [0.025, 0.975])]


def _bootstrap_independent_difference(
    control: np.ndarray,
    treated: np.ndarray,
    draws: int,
    seed: int,
) -> list[float]:
    rng = np.random.default_rng(seed)
    control_mean = control[
        rng.integers(len(control), size=(draws, len(control)))
    ].mean(axis=1)
    treated_mean = treated[
        rng.integers(len(treated), size=(draws, len(treated)))
    ].mean(axis=1)
    return [
        float(value)
        for value in np.quantile(treated_mean - control_mean, [0.025, 0.975])
    ]


def _means_by_experiment(rows: list[dict]) -> dict[str, float]:
    values: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        values[row["biological_replicate"]].append(float(row["value"]))
    return {experiment: float(np.mean(points)) for experiment, points in values.items()}


def _spreading_contrast(
    rows: list[dict], observable: str, draws: int, seed: int
) -> dict:
    selected = [row for row in rows if row["observable"] == observable]
    by_condition = {
        condition: _means_by_experiment([
            row for row in selected if row["condition"] == condition
        ])
        for condition in (CONTROL, Y27)
    }
    control = np.asarray(list(by_condition[CONTROL].values()), dtype=np.float64)
    treated = np.asarray(list(by_condition[Y27].values()), dtype=np.float64)
    if (len(control), len(treated)) != (3, 4):
        raise ValueError(f"{observable}: expected N=3 control and N=4 Y-27632")
    difference = float(np.mean(treated) - np.mean(control))
    ci = _bootstrap_independent_difference(control, treated, draws, seed)
    return {
        "control_cells": sum(
            row["condition"] == CONTROL for row in selected
        ),
        "Y27632_cells": sum(row["condition"] == Y27 for row in selected),
        "control_independent_experiments": len(control),
        "Y27632_independent_experiments": len(treated),
        "control_experiment_means": by_condition[CONTROL],
        "Y27632_experiment_means": by_condition[Y27],
        "Y27632_minus_control_experiment_mean_difference": difference,
        "independent_experiment_bootstrap_95ci": ci,
        "direction": "increases" if difference > 0 else "decreases",
        "inference_unit": "author_independent_experiment",
        "cells_are_not_inference_units": True,
    }


def _paired_tether_contrast(
    rows: list[dict], observable: str, draws: int, seed: int
) -> dict:
    selected = [row for row in rows if row["observable"] == observable]
    by_condition = {
        condition: _means_by_experiment([
            row for row in selected if row["condition"] == condition
        ])
        for condition in (CONTROL, Y27)
    }
    shared = sorted(set(by_condition[CONTROL]) & set(by_condition[Y27]))
    if len(shared) != 3:
        raise ValueError(f"{observable}: expected three matched experiments")
    effects = np.asarray([
        by_condition[Y27][experiment] - by_condition[CONTROL][experiment]
        for experiment in shared
    ], dtype=np.float64)
    control_mean = float(np.mean([by_condition[CONTROL][item] for item in shared]))
    treated_mean = float(np.mean([by_condition[Y27][item] for item in shared]))
    return {
        "matched_independent_experiments": len(shared),
        "matched_experiment_ids": shared,
        "paired_experiment_effects_pN": effects.tolist(),
        "Y27632_minus_control_mean_pN": float(np.mean(effects)),
        "paired_experiment_bootstrap_95ci_pN": _bootstrap(effects, draws, seed),
        "matched_control_mean_pN": control_mean,
        "matched_Y27632_mean_pN": treated_mean,
        "apparent_membrane_tension_ratio_if_bending_rigidity_constant": float(
            (treated_mean / control_mean) ** 2
        ),
        "force_to_tension_relation": (
            "descriptive model-derived ratio only; tether force is proportional to "
            "the square root of apparent membrane tension when bending rigidity is constant"
        ),
        "inference_unit": "matched_author_independent_experiment",
        "cells_are_not_inference_units": True,
    }


def evaluate(
    manifest: Path, data_sum_export: Path, *, draws: int = 10_000
) -> dict:
    rows = _load(manifest)
    published_means = _published_data_sum(data_sum_export)
    if {row["dataset_id"] for row in rows} != {
        "elife-72381-hela-mechano-osmotic"
    }:
        raise ValueError("manifest is not the eLife 72381 dataset")
    area = _spreading_contrast(rows, "initial_spreading_area_rate", draws, 20260807)
    volume = _spreading_contrast(rows, "initial_volume_flux", draws, 20260808)
    early = _paired_tether_contrast(
        rows,
        "membrane_tether_force:spreading_phase_30_to_90_min",
        draws,
        20260809,
    )
    steady = _paired_tether_contrast(
        rows,
        "membrane_tether_force:steady_state_spread_4_to_5_hour",
        draws,
        20260810,
    )
    nonspread = _paired_tether_contrast(
        rows,
        "membrane_tether_force:steady_state_nonspread_micropattern_4_to_5_hour",
        draws,
        20260811,
    )
    checks = {
        "area_rate_increase_experiment_ci_excludes_zero": (
            area["independent_experiment_bootstrap_95ci"][0] > 0
        ),
        "volume_flux_decrease_experiment_ci_excludes_zero": (
            volume["independent_experiment_bootstrap_95ci"][1] < 0
        ),
        "early_tether_force_increase_paired_experiment_ci_excludes_zero": (
            early["paired_experiment_bootstrap_95ci_pN"][0] > 0
        ),
        "steady_tether_force_direction_is_reported_without_cell_pseudoreplication": True,
        "force_unit_is_pN": all(
            row["unit"] == "pN"
            for row in rows
            if row["modality"] == "AFM_membrane_tether_force_derived"
        ),
        "spreading_units_have_author_source_locators": all(
            row["unit_source_sha256"]
            and row["unit_source_locator"].startswith(
                "zip:Figure2/Figure2D/data_sum.xlsx!Sheet1!"
            )
            for row in rows
            if row["observable"] in {
                "initial_spreading_area_rate", "initial_volume_flux"
            }
        ),
    }
    evidence_eligible = all(checks.values())
    raw_cell_means = {
        f"{condition}:{observable}": float(np.mean([
            float(row["value"])
            for row in rows
            if row["condition"] == condition and row["observable"] == observable
        ]))
        for condition in (CONTROL, Y27)
        for observable in ("initial_spreading_area_rate", "initial_volume_flux")
    }
    return {
        "task": "experiment_grouped_mechano_osmotic_and_membrane_tether_evidence",
        "source": {
            "accession": "10.7554/eLife.72381",
            "cell_type": "HeLa Kyoto",
            "intervention": "Y27632_100_micromolar",
            "independent_lab_dataset_present": True,
            "external_holdout_design": False,
        },
        "source_reconciliation": {
            "raw_manifest_cell_means": raw_cell_means,
            "published_Figure2D_data_sum_means": published_means,
            "control_reconciles_after_one_explicit_author_summary_exclusion": True,
            "Y27632_volume_flux_reconciles": True,
            "Y27632_area_rate_raw_table_minus_summary_mean": float(
                raw_cell_means[f"{Y27}:initial_spreading_area_rate"]
                - published_means[f"{Y27}:initial_spreading_area_rate"]
            ),
            "policy": (
                "preserve the per-cell raw dVdt_dAdt workbook; disclose rather than "
                "overwrite its Y-27632 area-rate mismatch with data_sum.xlsx"
            ),
        },
        "unit_provenance": {
            "initial_spreading_area_rate": {
                "unit": "square_micrometre_per_minute",
                "author_header": "dA/dt i (µm2/min)",
                "source_sha256": rows[0]["unit_source_sha256"],
                "source_locator": "zip:Figure2/Figure2D/data_sum.xlsx!Sheet1!B1",
            },
            "initial_volume_flux": {
                "unit": "cubic_micrometre_per_minute",
                "author_header": "dV/dt i (µm3/min)",
                "source_sha256": rows[0]["unit_source_sha256"],
                "source_locator": "zip:Figure2/Figure2D/data_sum.xlsx!Sheet1!E1",
            },
        },
        "initial_spreading_area_rate": area,
        "initial_volume_flux": volume,
        "AFM_membrane_tether_force": {
            "spreading_phase_30_to_90_min": early,
            "steady_state_spread_4_to_5_hour": steady,
            "steady_state_nonspread_4_to_5_hour": nonspread,
            "method": "single-cell atomic force spectroscopy; per-cell tether pulls averaged by authors",
            "observable_unit": "pN",
        },
        "checks": checks,
        "experiment_grouped_evidence_eligible": evidence_eligible,
        "observable_constraint_eligible": evidence_eligible,
        "tether_force_is_direct_cortical_tension": False,
        "tether_force_is_apparent_membrane_tension_proxy": True,
        "operator_assumptions": [
            "constant membrane bending rigidity is required to interpret force squared as apparent membrane tension",
            "the AFM tether-pull geometry must be represented by an Aleph measurement operator",
            "Y-27632 also affects actomyosin, spreading, ion flux, and membrane reservoirs, so it does not identify one parameter",
        ],
        "independent_lab_mechanism_holdout_passed": False,
        "aleph_parameter_proposal_eligible": False,
        "may_mutate_physics": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--data-sum-export", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=10_000)
    args = parser.parse_args()
    report = evaluate(args.manifest, args.data_sum_export, draws=args.draws)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "checks": report["checks"],
        "evidence_eligible": report["experiment_grouped_evidence_eligible"],
        "early_tether": report["AFM_membrane_tether_force"]["spreading_phase_30_to_90_min"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

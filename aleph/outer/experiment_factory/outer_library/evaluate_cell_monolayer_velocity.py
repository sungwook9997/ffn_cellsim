#!/usr/bin/env python3
"""Audit descriptive hdF velocity evidence without inventing replication."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def evaluate(manifest: Path, tensor: Path) -> dict:
    rows = [json.loads(line) for line in manifest.read_text(encoding="utf-8").splitlines()]
    if len(rows) != 3_154 or len({row["sample_id"] for row in rows}) != 832:
        raise ValueError("unexpected hdF manifest cardinality")
    if {row["aleph_authority"] for row in rows} != {"none"}:
        raise ValueError("outer observations must not carry Aleph authority")
    biological = {row["biological_replicate"] for row in rows}
    if biological != {"independent_culture_identifiers_not_reported"}:
        raise ValueError(f"unexpected biological grouping: {biological}")

    payload = np.load(tensor, allow_pickle=False)
    vectors = np.asarray(payload["vectors"], dtype=np.float64)
    condition = np.asarray(payload["condition"])
    if vectors.shape != (7_645, 5):
        raise ValueError(f"unexpected vector tensor shape: {vectors.shape}")
    if set(payload["split_group"].tolist()) != {"valentine-lab-ucsb"}:
        raise ValueError("raw velocity vectors must remain in one split group")
    if set(payload["biological_group"].tolist()) != {
        "independent_culture_identifiers_not_reported"
    }:
        raise ValueError("raw velocity vectors were promoted to false replicates")

    field_summary = {}
    for substrate in ("isotropic", "nematic"):
        selected = vectors[condition == substrate, :2]
        speed = np.hypot(selected[:, 0], selected[:, 1])
        moving = speed > 0
        direction = selected[moving] / speed[moving, None]
        field_summary[substrate] = {
            "technical_vector_count": int(len(selected)),
            "mean_speed_source_unit_unspecified": float(speed.mean()),
            "velocity_polar_order": float(np.hypot(*direction.mean(axis=0))),
            "velocity_nematic_order": float(np.hypot(
                np.mean(direction[:, 0] ** 2 - direction[:, 1] ** 2),
                np.mean(2 * direction[:, 0] * direction[:, 1]),
            )),
        }

    by_condition_observable: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        if row["coordinate_name"] == "time":
            by_condition_observable[(row["condition"], row["observable"])].append(
                float(row["value"])
            )
    treated = "nematic_substrate_FAK_inhibitor_1_micromolar"
    control = "nematic_substrate_no_FAKi_control"
    contrasts = {}
    for observable in (
        "cell_number_density", "cell_substrate_order_parameter", "ensemble_mean_speed",
        "ensemble_mean_absolute_velocity_x", "ensemble_mean_absolute_velocity_y",
    ):
        treated_values = np.asarray(by_condition_observable[(treated, observable)])
        control_values = np.asarray(by_condition_observable[(control, observable)])
        if not len(treated_values) or not len(control_values):
            raise ValueError(f"missing longitudinal observable: {observable}")
        contrasts[observable] = {
            "FAKi_timepoint_mean": float(treated_values.mean()),
            "control_timepoint_mean": float(control_values.mean()),
            "FAKi_minus_control_timepoint_mean": float(
                treated_values.mean() - control_values.mean()
            ),
            "inferential_ci95": None,
        }

    return {
        "task": "hdF_collective_velocity_representation_and_descriptive_FAKi_audit",
        "source": {
            "accession": "10.25349/D96W58",
            "paper_doi": "10.1098/rsif.2023.0160",
            "license": "cc0-1.0",
        },
        "manifest": {"observations": len(rows), "samples": 832},
        "raw_velocity_tensor": {
            "technical_vectors": len(vectors), "features": vectors.shape[1],
            "mandatory_split_group": "valentine-lab-ucsb",
            "biological_replicate_count": 1,
        },
        "field_descriptions": field_summary,
        "FAKi_descriptive_contrasts": contrasts,
        "representation_training_eligible": True,
        "representation_use_constraints": [
            "group every source row by the full dataset/lab split group",
            "do not split spatial vectors or longitudinal time points across folds",
            "do not use this single study as its own external holdout",
        ],
        "inferential_status": "descriptive_only_independent_culture_ids_not_reported",
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "uncertainty_holdout_passed": False,
        "observable_constraint_eligible": False,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": [
            "the release does not identify independent cultures or experimental repeats",
            "spatial vectors, frames, and longitudinal time points are correlated technical observations",
            "the released velocity columns do not state their unit in the data README",
            "the FAKi and control trajectories differ strongly in density, so their difference is not a causal FAK effect estimate",
            "no confidence interval, parameter calibration, or external-holdout claim is emitted",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.manifest, args.tensor)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "inferential_status": report["inferential_status"],
        "technical_vectors": report["raw_velocity_tensor"]["technical_vectors"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

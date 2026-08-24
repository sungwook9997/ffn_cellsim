#!/usr/bin/env python3
"""Descriptive-only audit for the representative electrotaxis PIV pair.

The Zenodo record releases one representative tissue per condition. Sixty
frames are useful spatial/temporal measurements, not sixty biological
replicates. This evaluator computes transparent condition summaries while
hard-refusing external-holdout or parameter-calibration authority.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


EXPECTED_METRICS = {
    "mean_velocity_x", "mean_velocity_y", "mean_speed", "rms_speed",
    "speed_p90", "velocity_polar_order", "velocity_nematic_order",
    "valid_vector_fraction",
}


def evaluate(manifest: Path) -> dict:
    samples: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    biological: dict[str, set[str]] = defaultdict(set)
    with manifest.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            condition = row["condition"]
            if condition not in {"control", "stimulated"}:
                raise ValueError(f"unexpected PIV condition: {condition}")
            key = (condition, row["sample_id"])
            metric = row["observable"]
            if metric in samples[key]:
                raise ValueError(f"duplicate PIV metric: {key}/{metric}")
            samples[key][metric] = float(row["value"])
            biological[condition].add(row["biological_replicate"])
    conditions = {
        condition: sorted(
            (metrics for (item_condition, _), metrics in samples.items()
             if item_condition == condition),
            key=lambda item: tuple(sorted(item.items())),
        )
        for condition in ("control", "stimulated")
    }
    frame_counts = {condition: len(items) for condition, items in conditions.items()}
    if frame_counts != {"control": 60, "stimulated": 60}:
        raise ValueError(f"expected 60 frames per condition, found {frame_counts}")
    metric_sets = {
        frozenset(item) for items in conditions.values() for item in items
    }
    if metric_sets != {frozenset(EXPECTED_METRICS)}:
        raise ValueError(f"PIV frames must contain the exact metric panel: {metric_sets}")
    metrics = sorted(EXPECTED_METRICS)
    summaries = {}
    for metric in metrics:
        control = np.asarray([item[metric] for item in conditions["control"]])
        stimulated = np.asarray([item[metric] for item in conditions["stimulated"]])
        summaries[metric] = {
            "control_frame_mean": float(np.mean(control)),
            "stimulated_frame_mean": float(np.mean(stimulated)),
            "stimulated_minus_control_frame_mean": float(
                np.mean(stimulated) - np.mean(control)
            ),
            "inferential_ci95": None,
        }
    biological_counts = {
        condition: len(values) for condition, values in biological.items()
    }
    if biological_counts != {"control": 1, "stimulated": 1}:
        raise ValueError(f"expected one representative tissue per condition: {biological_counts}")
    return {
        "task": "representative_MDCK_electrotaxis_PIV_description",
        "frame_counts": frame_counts,
        "biological_replicate_counts": biological_counts,
        "metric_count": len(metrics),
        "metric_summaries": summaries,
        "inferential_status": "descriptive_only_no_biological_replication",
        "external_holdout_passed": False,
        "observable_constraint_eligible": False,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": [
            "the release contains one representative tissue per condition",
            "frames and PIV vectors are temporally/spatially correlated measurements, not biological replicates",
            "source velocity units are not declared in the released MAT object",
            "no inferential confidence interval or external holdout claim is emitted",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "inferential_status": report["inferential_status"],
        "metric_count": report["metric_count"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

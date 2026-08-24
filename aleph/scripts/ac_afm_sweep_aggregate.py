#!/usr/bin/env python3
"""Aggregate only fully accepted AFM speed/seed paths using between-seed scatter."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np


def aggregate(paths: list[Path]) -> dict[str, Any]:
    """Validate accepted paths and return a speed/depth ensemble; reject partial material curves."""
    records = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    if not records:
        raise ValueError("at least one accepted AFM path is required")
    if any(
        record.get("schema") != "afm-force-indentation-path@1"
        or record.get("status") != "ACCEPTED_PATH"
        for record in records
    ):
        raise ValueError("every input must be an ACCEPTED afm-force-indentation-path@1 record")
    protocol = records[0]["protocol"]
    population = records[0]["population"]
    build_commit = str(records[0].get("build_commit", "")).strip()
    if not build_commit:
        raise ValueError("every AFM path must carry a build commit")
    if any(record["protocol"] != protocol or record["population"] != population for record in records[1:]):
        raise ValueError("AFM ensemble paths must share one immutable protocol and population")
    if any(str(record.get("build_commit", "")).strip() != build_commit for record in records[1:]):
        raise ValueError("AFM ensemble paths must share one verified build commit")

    paths_by_speed: dict[float, list[dict[str, Any]]] = defaultdict(list)
    seen: set[tuple[float, int]] = set()
    for record in records:
        speed = float(record["speed_um_s"])
        seed = int(record["seed"])
        key = (speed, seed)
        if key in seen:
            raise ValueError(f"duplicate AFM path for speed/seed {key}")
        seen.add(key)
        points = record["points"]
        if not points or any(
            not point.get("accepted")
            or float(point["max_projected_force_pn"])
            > float(point["projected_force_tolerance_pn"])
            for point in points
        ):
            raise ValueError("every exposed F-delta point must pass the stored physical acceptance predicate")
        paths_by_speed[speed].append(record)

    rows: list[dict[str, Any]] = []
    seed_sets: dict[str, list[int]] = {}
    for speed, speed_records in sorted(paths_by_speed.items()):
        speed_records.sort(key=lambda record: int(record["seed"]))
        seeds = sorted(int(record["seed"]) for record in speed_records)
        if len(seeds) < 3:
            raise ValueError(f"speed {speed:g} um/s has {len(seeds)} seeds; at least three are required")
        depths = tuple(float(point["depth_um"]) for point in speed_records[0]["points"])
        if any(
            tuple(float(point["depth_um"]) for point in record["points"]) != depths
            for record in speed_records[1:]
        ):
            raise ValueError("all seeds at one speed must share the same accepted depth grid")
        seed_sets[f"{speed:.17g}"] = seeds
        for index, depth in enumerate(depths):
            reactions = np.asarray(
                [float(record["points"][index]["reaction_z_pn"]) for record in speed_records],
                dtype=np.float64,
            )
            rows.append({
                "speed_um_s": speed,
                "depth_um": depth,
                "n_seeds": len(seeds),
                "reaction_mean_pn": float(reactions.mean()),
                "reaction_sample_sd_pn": float(reactions.std(ddof=1)),
                "reaction_min_pn": float(reactions.min()),
                "reaction_max_pn": float(reactions.max()),
                "seed_values_pn": reactions.tolist(),
            })
    if any(
        record.get("quantitative_claim") == "ELIGIBLE_FOR_PROTOCOL_AGGREGATION"
        for record in records
    ):
        raise RuntimeError(
            "quantitative aggregation remains blocked until the path schema carries an accepted "
            "post-induction, dt-converged active-cortex state handoff"
        )
    return {
        "schema": "afm-force-indentation-ensemble@1",
        "status": "ACCEPTED_ENSEMBLE",
        "quantitative_claim": "MECHANISM_DEMO_NOT_QUANTITATIVE",
        "build_commit": build_commit,
        "protocol": protocol,
        "population": population,
        "seeds_by_speed": seed_sets,
        "points": rows,
        "uncertainty_rule": "sample standard deviation across distinct seeds; within-run SEM is not used",
        "scope": "Accepted F-delta ensemble only. No cortical-tension inversion is implied.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = aggregate(args.inputs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

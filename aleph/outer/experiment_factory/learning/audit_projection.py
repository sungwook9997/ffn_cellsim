#!/usr/bin/env python3
"""Independent structural and leakage audit for a learning projection."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from .build_projection import canonical, digest, atomic


FORBIDDEN_KEYS = {"text", "caption", "pixels", "publisher_payload", "raw_payload", "href"}


def audit(path: Path, output: Path) -> dict[str, Any]:
    errors: Counter[str] = Counter()
    ids = set()
    sources = set()
    group_splits: dict[str, set[str]] = defaultdict(set)
    splits: Counter[str] = Counter()
    cohort = 0
    rows = []
    for line in path.read_bytes().splitlines():
        row = json.loads(line)
        rows.append(row)
        expected = dict(row)
        claimed = expected.pop("projection_sha256")
        if digest(canonical(expected)) != claimed:
            errors["projection_sha256"] += 1
        if row["projection_id"] in ids:
            errors["duplicate_projection_id"] += 1
        if row["source_family_id"] in sources:
            errors["duplicate_source"] += 1
        ids.add(row["projection_id"])
        sources.add(row["source_family_id"])
        group_splits[row["leakage_group_id"]].add(row["split"])
        splits[row["split"]] += 1
        cohort += row["evaluation_cohort"].startswith("frozen_300")
        if row["authority_status"] != "proposed":
            errors["authority_promotion"] += 1
        if row["supervision_policy"]["ambiguous_observations_are_labels"]:
            errors["ambiguous_label_promotion"] += 1
        def visit(value: Any) -> None:
            if isinstance(value, dict):
                for key, child in value.items():
                    if key.lower() in FORBIDDEN_KEYS:
                        errors["forbidden_raw_field"] += 1
                    visit(child)
            elif isinstance(value, list):
                for child in value:
                    visit(child)
        visit(row)
    errors["leakage_split"] += sum(len(value) != 1 for value in group_splits.values())
    result = {
        "schema": "aleph.external_training.learning_projection_audit.v1",
        "authority_status": "proposed",
        "sources": len(rows), "unique_sources": len(sources), "unique_projection_ids": len(ids),
        "leakage_groups": len(group_splits), "leakage_split_violations": errors["leakage_split"],
        "split_source_counts": dict(sorted(splits.items())),
        "frozen_300_evaluation_cohort_sources": cohort,
        "ambiguous_observations_promoted_to_labels": errors["ambiguous_label_promotion"],
        "forbidden_raw_fields": errors["forbidden_raw_field"],
        "errors": sum(errors.values()), "error_classes": dict(sorted(errors.items())),
        "projection_file_sha256": digest(path.read_bytes()),
    }
    atomic(output, json.dumps(result, indent=2, sort_keys=True).encode() + b"\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.projection, args.output)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["errors"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

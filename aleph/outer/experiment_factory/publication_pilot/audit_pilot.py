from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


FORBIDDEN_PACKET_KEYS = {"evidence_text", "raw_text", "caption", "caption_text", "pixels", "image_bytes", "section_path", "decided_by"}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict]:
    with path.open() as handle:
        return [json.loads(line) for line in handle if line.strip()]


def gz_rows(path: Path) -> list[dict]:
    with gzip.open(path, "rt") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def recursive_keys(value: object) -> list[str]:
    output = []
    if isinstance(value, dict):
        for key, child in value.items():
            output.append(key)
            output.extend(recursive_keys(child))
    elif isinstance(value, list):
        for child in value:
            output.extend(recursive_keys(child))
    return output


def audit(results: Path, parent_candidates: Path) -> dict:
    failures: list[str] = []
    manifest = json.loads((results / "SHA256SUMS.json").read_text())
    mismatches = [name for name, expected in manifest.items() if not (results / name).exists() or digest(results / name) != expected]
    if mismatches:
        failures.append("manifest_mismatch")

    selected = rows(results / "candidates_50.jsonl")
    forms_a = rows(results / "forms_A.jsonl")
    forms_b = rows(results / "forms_B.jsonl")
    pending = rows(results / "adjudication_forms.jsonl")
    machine = gz_rows(results / "machine_preannotations.jsonl.gz")
    summary = json.loads((results / "summary.json").read_text())
    readiness = json.loads((results / "benchmark_readiness.json").read_text())

    sources = [row["source_family_id"] for row in selected]
    if len(selected) != 50 or len(set(sources)) != 50:
        failures.append("candidate_identity")
    source_set = set(sources)
    for label, collection, getter in (
        ("A", forms_a, lambda row: row["source"]["source_family_id"]),
        ("B", forms_b, lambda row: row["source"]["source_family_id"]),
        ("adjudication", pending, lambda row: row["source_family_id"]),
        ("machine", machine, lambda row: row["source_family_id"]),
    ):
        if len(collection) != 50 or {getter(row) for row in collection} != source_set:
            failures.append(f"{label}_source_set")

    assignment_ids = [row["assignment_id"] for row in forms_a + forms_b]
    if len(set(assignment_ids)) != 100:
        failures.append("assignment_identity")
    if any(row["annotator_id"] is not None or row["claims"] or row["review_state"] != "unassigned" for row in forms_a + forms_b):
        failures.append("nonblank_blinded_form")
    if any(not row["blinding"]["machine_preannotation_hidden"] for row in forms_a + forms_b):
        failures.append("machine_not_blinded")
    if any(row["adjudicator_id"] is not None or row["adjudicated_claims"] or row["review_state"] != "awaiting_two_human_submissions" for row in pending):
        failures.append("adjudication_prepopulated")

    packet_keys = Counter(key for row in machine for key in recursive_keys(row))
    forbidden = sorted(FORBIDDEN_PACKET_KEYS & set(packet_keys))
    if forbidden:
        failures.append("forbidden_packet_key")
    if any(row["human_verified"] or row["review_state"] != "machine_candidate" or row["authority_status"] != "proposed" for row in machine):
        failures.append("machine_authority_promotion")

    parent = rows(parent_candidates)
    components: dict[str, set[str]] = defaultdict(set)
    for row in parent:
        components[row["leakage_component_id"]].add(row["source_family_id"])
    split_components = []
    for component_id, members in components.items():
        inside = members & source_set
        if inside and inside != members:
            split_components.append(component_id)
    # The deliberately deferred 17-paper component and all unselected components are outside the
    # sample; only partial intersections are violations.
    if split_components:
        failures.append("leakage_component_split")

    if readiness["benchmark_metrics_available"] or any(readiness[key] for key in ("human_A_submissions", "human_B_submissions", "adjudicated_submissions")):
        failures.append("fabricated_human_completion")
    if summary["forms_contain_machine_predictions"] or summary["machine_preannotation_is_ground_truth"]:
        failures.append("summary_promotion")

    return {
        "schema": "aleph.external_training.publication_pilot_audit.v1",
        "authority_status": "proposed",
        "status": "pass" if not failures else "fail",
        "failures": failures,
        "manifest_files": len(manifest),
        "manifest_mismatches": mismatches,
        "selected_sources": len(selected),
        "assignment_forms": len(forms_a) + len(forms_b),
        "machine_packets": len(machine),
        "machine_experiment_records": sum(len(row["experiment_records"]) for row in machine),
        "forbidden_packet_keys": forbidden,
        "split_leakage_components": sorted(split_components),
        "human_A_submissions": readiness["human_A_submissions"],
        "human_B_submissions": readiness["human_B_submissions"],
        "adjudicated_submissions": readiness["adjudicated_submissions"],
        "benchmark_metrics_available": readiness["benchmark_metrics_available"],
        "machine_preannotation_is_ground_truth": summary["machine_preannotation_is_ground_truth"],
        "forms_contain_machine_predictions": summary["forms_contain_machine_predictions"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--parent-candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.results, args.parent_candidates)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    if result["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()

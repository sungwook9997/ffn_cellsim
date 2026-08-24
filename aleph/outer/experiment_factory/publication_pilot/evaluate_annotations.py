from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from .build_pilot import canonical


CLAIM_TYPES = {
    "experiment_boundary", "organism", "tissue", "cell_type", "cell_line", "cell_state",
    "condition_role", "condition_factor", "protocol_modality", "instrument", "calibration",
    "measured_entity", "reported_value", "normalized_unit", "uncertainty_type", "biological_n",
    "technical_n", "comparison_direction", "figure_panel", "table_cell",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_jsonl(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else path.open
    if path.suffix == ".gz":
        with opener(path, "rt") as handle:  # type: ignore[arg-type]
            return [json.loads(line) for line in handle if line.strip()]
    with path.open() as handle:
        return [json.loads(line) for line in handle if line.strip()]


def validate_submission(row: dict, adjudicated: bool = False) -> list[str]:
    errors = []
    required = {"schema", "authority_status", "review_state", "assignment_id", "annotation_slot", "annotator_id", "source_family_id", "claims", "completed_at"}
    if set(row) != required:
        errors.append("submission_fields")
    if row.get("schema") != "aleph.external_training.annotation_submission.v1":
        errors.append("schema")
    if row.get("authority_status") != "proposed":
        errors.append("authority_status")
    expected_state = "adjudicated" if adjudicated else "single_annotator"
    if row.get("review_state") != expected_state:
        errors.append("review_state")
    expected_slot = "ADJUDICATION" if adjudicated else {"A", "B"}
    if adjudicated and row.get("annotation_slot") != expected_slot:
        errors.append("annotation_slot")
    if not adjudicated and row.get("annotation_slot") not in expected_slot:
        errors.append("annotation_slot")
    if not row.get("annotator_id") or not row.get("completed_at"):
        errors.append("human_identity_or_time_missing")
    for claim in row.get("claims", []):
        if set(claim) != {"claim_type", "experiment_anchor", "value", "evidence_key", "status", "notes"}:
            errors.append("claim_fields")
            continue
        if claim["claim_type"] not in CLAIM_TYPES:
            errors.append("claim_type")
        if not claim["experiment_anchor"] or not claim["evidence_key"]:
            errors.append("claim_grounding")
        if adjudicated and claim["status"] != "adjudicated":
            errors.append("unadjudicated_claim_in_gold")
    return sorted(set(errors))


def claim_key(claim: dict, include_evidence: bool = True) -> tuple:
    key = (claim["claim_type"], claim["experiment_anchor"], canonical(claim["value"]))
    return key + ((claim["evidence_key"],) if include_evidence else ())


def claim_sets(row: dict, include_evidence: bool = True) -> dict[str, set[tuple]]:
    output: dict[str, set[tuple]] = defaultdict(set)
    for claim in row["claims"]:
        if claim["status"] == "missing":
            continue
        output[claim["claim_type"]].add(claim_key(claim, include_evidence))
    return output


def prf(predicted: set[tuple], reference: set[tuple]) -> dict:
    tp = len(predicted & reference)
    fp = len(predicted - reference)
    fn = len(reference - predicted)
    precision = tp / (tp + fp) if tp + fp else (1.0 if not reference else 0.0)
    recall = tp / (tp + fn) if tp + fn else (1.0 if not predicted else 0.0)
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {"tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1}


def index_submissions(rows: list[dict], adjudicated: bool = False) -> dict[str, dict]:
    indexed = {}
    for row in rows:
        errors = validate_submission(row, adjudicated)
        if errors:
            raise ValueError(f"invalid submission {row.get('assignment_id')}: {','.join(errors)}")
        source = row["source_family_id"]
        if source in indexed:
            raise ValueError(f"duplicate source submission: {source}")
        indexed[source] = row
    return indexed


def agreement(rows_a: list[dict], rows_b: list[dict]) -> dict:
    a = index_submissions(rows_a)
    b = index_submissions(rows_b)
    if set(a) != set(b):
        raise ValueError("A/B source sets differ")
    totals: dict[str, Counter] = defaultdict(Counter)
    exact_sources = 0
    for source in sorted(a):
        if a[source]["annotator_id"] == b[source]["annotator_id"]:
            raise ValueError(f"same annotator in A/B for {source}")
        sets_a, sets_b = claim_sets(a[source]), claim_sets(b[source])
        if {k for values in sets_a.values() for k in values} == {k for values in sets_b.values() for k in values}:
            exact_sources += 1
        for field in CLAIM_TYPES:
            metrics = prf(sets_a.get(field, set()), sets_b.get(field, set()))
            totals[field].update({k: metrics[k] for k in ("tp", "fp", "fn")})
    per_field = {}
    for field, count in sorted(totals.items()):
        predicted = {("p", i) for i in range(count["tp"] + count["fp"])}
        reference = {("p", i) for i in range(count["tp"])} | {("r", i) for i in range(count["fn"])}
        per_field[field] = prf(predicted, reference)
    return {
        "schema": "aleph.external_training.annotation_agreement.v1",
        "authority_status": "proposed",
        "sources": len(a),
        "exact_source_claim_sets": exact_sources,
        "per_field_exact_claim_agreement": per_field,
        "note": "Symmetric exact-set agreement; neither annotator is treated as truth.",
    }


def submission_sha(row: dict) -> str:
    return hashlib.sha256((canonical(row) + "\n").encode()).hexdigest()


def prepare_adjudication(rows_a: list[dict], rows_b: list[dict]) -> list[dict]:
    a = index_submissions(rows_a)
    b = index_submissions(rows_b)
    if set(a) != set(b):
        raise ValueError("A/B source sets differ")
    packets = []
    for source in sorted(a):
        if a[source]["annotator_id"] == b[source]["annotator_id"]:
            raise ValueError(f"same annotator in A/B for {source}")
        claims_a = {claim_key(claim): claim for claim in a[source]["claims"] if claim["status"] != "missing"}
        claims_b = {claim_key(claim): claim for claim in b[source]["claims"] if claim["status"] != "missing"}
        agreed_keys = sorted(set(claims_a) & set(claims_b), key=repr)
        only_a = sorted(set(claims_a) - set(claims_b), key=repr)
        only_b = sorted(set(claims_b) - set(claims_a), key=repr)
        packets.append({
            "schema": "aleph.external_training.adjudication_packet.v1",
            "authority_status": "proposed",
            "review_state": "double_adjudication_pending",
            "source_family_id": source,
            "submission_A_sha256": submission_sha(a[source]),
            "submission_B_sha256": submission_sha(b[source]),
            "annotator_identities_blinded": True,
            "agreed_claims": [claims_a[key] for key in agreed_keys],
            "disagreements": [
                {"side": "A_only", "claim": claims_a[key]} for key in only_a
            ] + [
                {"side": "B_only", "claim": claims_b[key]} for key in only_b
            ],
            "adjudicator_id": None,
            "adjudicated_claims": [],
            "completed_at": None,
            "automatic_adjudication_performed": False,
        })
    return packets


def machine_index(path: Path) -> dict[str, dict]:
    rows = read_jsonl(path)
    return {row["source_family_id"]: row for row in rows}


def evaluate_machine(machine_rows: dict[str, dict], gold_rows: list[dict]) -> dict:
    gold = index_submissions(gold_rows, adjudicated=True)
    missing = sorted(set(gold) - set(machine_rows))
    if missing:
        raise ValueError(f"machine packets missing {len(missing)} adjudicated sources")
    totals: dict[str, Counter] = defaultdict(Counter)
    masked_ambiguous = 0
    for source, reference_row in gold.items():
        predicted_by_type: dict[str, set[tuple]] = defaultdict(set)
        for claim in machine_rows[source]["atomic_claims"]:
            if claim["ambiguous"]:
                masked_ambiguous += 1
                continue
            predicted_by_type[claim["claim_type"]].add((claim["claim_type"], claim["experiment_anchor"], canonical(claim["value"]), claim["evidence_key"]))
        reference = claim_sets(reference_row)
        for field in CLAIM_TYPES:
            metrics = prf(predicted_by_type.get(field, set()), reference.get(field, set()))
            totals[field].update({k: metrics[k] for k in ("tp", "fp", "fn")})
    per_field = {}
    for field, count in sorted(totals.items()):
        tp, fp, fn = count["tp"], count["fp"], count["fn"]
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        per_field[field] = {"tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1}
    return {
        "schema": "aleph.external_training.machine_vs_adjudicated_benchmark.v1",
        "authority_status": "proposed",
        "adjudicated_sources": len(gold),
        "ambiguous_machine_claims_masked": masked_ambiguous,
        "per_field_exact_metrics": per_field,
        "human_gold_required": True,
    }


def write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n")


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(canonical(row) + "\n" for row in rows))


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    p_agree = sub.add_parser("agreement")
    p_agree.add_argument("--a", type=Path, required=True)
    p_agree.add_argument("--b", type=Path, required=True)
    p_agree.add_argument("--output", type=Path, required=True)
    p_prepare = sub.add_parser("prepare-adjudication")
    p_prepare.add_argument("--a", type=Path, required=True)
    p_prepare.add_argument("--b", type=Path, required=True)
    p_prepare.add_argument("--output", type=Path, required=True)
    p_eval = sub.add_parser("evaluate")
    p_eval.add_argument("--machine", type=Path, required=True)
    p_eval.add_argument("--adjudicated", type=Path, required=True)
    p_eval.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "agreement":
        result = agreement(read_jsonl(args.a), read_jsonl(args.b))
        write(args.output, result)
    elif args.command == "prepare-adjudication":
        packets = prepare_adjudication(read_jsonl(args.a), read_jsonl(args.b))
        write_rows(args.output, packets)
        print(json.dumps({"packets": len(packets), "automatic_adjudication_performed": False}, sort_keys=True))
        return
    else:
        result = evaluate_machine(machine_index(args.machine), read_jsonl(args.adjudicated))
        write(args.output, result)
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


SCHEMA = "aleph.external_training.publication_pilot.v1"
DOMAIN_TARGETS = {
    "adhesion_ecm_mechanics": 4,
    "cell_state_mechanics": 4,
    "cortical_membrane_pressure": 5,
    "cytoskeleton_motor_rheology": 4,
    "division_morphogenesis_tissue": 4,
    "flow_shear": 4,
    "measurement_inference": 5,
    "mechanotransduction": 4,
    "migration_invasion": 4,
    "nuclear_mechanics": 4,
    "organelle_general_mechanics": 4,
    "osmotic_volume_poroelasticity": 4,
}
YEAR_TARGETS = {"2010_2013": 7, "2014_2017": 11, "2018_2021": 12, "2022_2026": 20}
PASS_TARGETS = {"first": 27, "second": 23}
SPARSE_MODALITIES = {
    "magnetic_tweezers",
    "optical_tweezers",
    "micropipette_aspiration",
    "piv",
    "single_cell_sequencing",
    "traction_force_microscopy",
}


def canonical(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def write_jsonl(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = "".join(canonical(row) + "\n" for row in rows).encode()
    path.write_bytes(payload)


def write_jsonl_gz(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = "".join(canonical(row) + "\n" for row in rows).encode()
    path.write_bytes(gzip.compress(payload, compresslevel=9, mtime=0))


def component_groups(rows: list[dict]) -> list[list[dict]]:
    groups: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        groups[row["leakage_component_id"]].append(row)
    return [sorted(group, key=lambda row: row["selection_rank"]) for group in groups.values()]


def selection_score(row: dict, counts: dict[str, Counter], seen: dict[str, set[str]]) -> tuple:
    domain = row["primary_domain"]
    year = row["year_stratum"]
    acquisition_pass = row["acquisition_pass"]
    modalities = set(row["modality_tags"])
    states = set(row["cell_state_tags"])
    mechanics = set(row["mechanics_tags"])
    score = 0.0
    score += 120.0 if counts["domain"][domain] < DOMAIN_TARGETS[domain] else -80.0
    score += 35.0 if counts["year"][year] < YEAR_TARGETS[year] else -15.0
    score += 25.0 if counts["pass"][acquisition_pass] < PASS_TARGETS[acquisition_pass] else -10.0
    score += 60.0 * len((modalities & SPARSE_MODALITIES) - seen["modality"])
    score += 25.0 * len(modalities - seen["modality"])
    score += 8.0 * len(states - seen["state"])
    score += 4.0 * len(mechanics - seen["mechanics"])
    score += min(int(row["observation_locator_count"]), 12) * 0.5
    score += int(row["priority_score"]) * 0.05
    return (score, -int(row["selection_rank"]))


def select_fifty(rows: list[dict]) -> tuple[list[dict], dict]:
    groups = component_groups(rows)
    oversized = {g[0]["leakage_component_id"] for g in groups if len(g) > 3}
    small_multi = sorted(
        (g for g in groups if 1 < len(g) <= 3),
        key=lambda g: (len(g), min(r["selection_rank"] for r in g)),
    )
    if not small_multi:
        raise ValueError("pilot requires at least one intact multi-source leakage component")

    # Seed one intact small component. Prefer cross-domain coverage, then the lowest ranks.
    seed = max(
        small_multi,
        key=lambda g: (len({r["primary_domain"] for r in g}), -len(g), -min(r["selection_rank"] for r in g)),
    )
    selected = list(seed)
    selected_components = {seed[0]["leakage_component_id"]}
    counts = {
        "domain": Counter(r["primary_domain"] for r in selected),
        "year": Counter(r["year_stratum"] for r in selected),
        "pass": Counter(r["acquisition_pass"] for r in selected),
    }
    seen = {
        "modality": {x for r in selected for x in r["modality_tags"]},
        "state": {x for r in selected for x in r["cell_state_tags"]},
        "mechanics": {x for r in selected for x in r["mechanics_tags"]},
    }

    # Other multi-source components are deferred rather than partially sampled.
    singletons = [g[0] for g in groups if len(g) == 1 and g[0]["leakage_component_id"] not in oversized]
    while len(selected) < 50:
        row = max(singletons, key=lambda r: selection_score(r, counts, seen))
        singletons.remove(row)
        selected.append(row)
        selected_components.add(row["leakage_component_id"])
        counts["domain"][row["primary_domain"]] += 1
        counts["year"][row["year_stratum"]] += 1
        counts["pass"][row["acquisition_pass"]] += 1
        seen["modality"].update(row["modality_tags"])
        seen["state"].update(row["cell_state_tags"])
        seen["mechanics"].update(row["mechanics_tags"])

    selected.sort(key=lambda row: row["selection_rank"])
    selected_ids = {row["source_family_id"] for row in selected}
    split_components = []
    for group in groups:
        inside = sum(row["source_family_id"] in selected_ids for row in group)
        if inside not in (0, len(group)):
            split_components.append(group[0]["leakage_component_id"])
    audit = {
        "split_components": split_components,
        "oversized_components_deferred": sorted(oversized),
        "selected_multi_source_components": sorted(
            component_id for component_id in selected_components
            if sum(r["leakage_component_id"] == component_id for r in selected) > 1
        ),
    }
    return selected, audit


def record_files_by_source(records_root: Path, wanted: set[str]) -> dict[str, Path]:
    found: dict[str, Path] = {}
    for path in sorted(records_root.glob("*.jsonl")):
        with path.open() as handle:
            first = handle.readline()
        if not first:
            continue
        source = json.loads(first)["source"]["source_family_id"]
        if source in wanted:
            found[source] = path
    return found


def compact_term(term: dict) -> dict:
    return {
        "label": term.get("label"),
        "identifier": term.get("identifier"),
        "normalization_status": term.get("normalization_status"),
        "candidates": term.get("candidates", []),
    }


def compact_record(record: dict) -> dict:
    claims = machine_claims(record)
    return {
        "experiment_record_id": record["experiment_record_id"],
        "review_state": record["review_state"],
        "machine_candidate_only": True,
        "experiment_evidence_anchor": evidence_key(record["protocols"][0]["evidence"]),
        "observation_count": len(record["observations"]),
        "comparison_count": len(record["comparisons"]),
        "atomic_claim_count": len(claims),
        "ambiguous_claim_count": sum(claim["ambiguous"] for claim in claims),
        "source_record_sha256": record["provenance"].get("candidate_record_sha256"),
    }


def evidence_key(evidence: list[dict]) -> str:
    if not evidence:
        return "missing-evidence"
    item = evidence[0]
    parts = (
        item.get("document_sha256"), item.get("locator_type"), item.get("locator_id"),
        item.get("char_start"), item.get("char_end"), item.get("text_or_asset_sha256"),
    )
    return "|".join("" if value is None else str(value) for value in parts)


def machine_claims(record: dict) -> list[dict]:
    claims: list[dict] = []
    default_evidence = record["protocols"][0]["evidence"]
    anchor = evidence_key(default_evidence)

    def add(claim_type: str, value: object, evidence: list[dict] | None = None, ambiguous: bool = False) -> None:
        if value is None or value == "" or value == []:
            return
        claims.append({
            "claim_type": claim_type,
            "experiment_anchor": anchor,
            "value": value,
            "evidence_key": evidence_key(evidence or default_evidence),
            "machine_candidate_only": True,
            "ambiguous": bool(ambiguous),
        })

    add("experiment_boundary", "present")
    for system in record["biological_systems"]:
        for name in ("organism", "tissue", "cell_type", "cell_line"):
            term = system[name]
            add(name, term.get("identifier") or term.get("label"), ambiguous=term.get("normalization_status") == "ambiguous")
        for term in system["cell_states"]:
            add("cell_state", term.get("identifier") or term.get("label"), ambiguous=term.get("normalization_status") == "ambiguous")
    for condition in record["conditions"]:
        add("condition_role", condition["role"])
        for factor in condition["factors"]:
            add("condition_factor", {
                "factor_type": factor["factor_type"], "name": factor["name"],
                "value": factor["value"], "unit": factor["unit"], "duration": factor["duration"],
            }, ambiguous=factor["status"] == "ambiguous")
    for protocol in record["protocols"]:
        add("protocol_modality", protocol["modality"], protocol["evidence"])
        add("instrument", protocol["instrument"], protocol["evidence"])
        add("calibration", protocol["calibration"], protocol["evidence"])
    for observation in record["observations"]:
        ambiguous = observation["extraction"]["status"] == "ambiguous" or observation["value_state"] == "ambiguous"
        add("measured_entity", observation["measured_entity"], observation["evidence"], ambiguous)
        value = observation["value"]
        if isinstance(value, dict):
            value = value.get("reported_value", value)
        add("reported_value", value, observation["evidence"], ambiguous)
        add("normalized_unit", observation["unit"]["normalized"], observation["evidence"], ambiguous)
        for uncertainty in observation["uncertainty"]:
            add("uncertainty_type", uncertainty["type"], observation["evidence"], ambiguous)
        add("biological_n", observation["replicates"]["biological_n"], observation["evidence"], ambiguous)
        add("technical_n", observation["replicates"]["technical_n"], observation["evidence"], ambiguous)
        visual = observation.get("visual_payload")
        if visual:
            add("figure_panel", visual.get("panel_id"), observation["evidence"], ambiguous)
    for comparison in record["comparisons"]:
        add("comparison_direction", comparison["effect_direction"], comparison["evidence"], comparison["status"] == "ambiguous")
    claims.sort(key=lambda item: (item["claim_type"], item["experiment_anchor"], canonical(item["value"]), item["evidence_key"]))
    return claims


def build_preannotations(selected: list[dict], records_root: Path) -> list[dict]:
    wanted = {row["source_family_id"] for row in selected}
    files = record_files_by_source(records_root, wanted)
    missing = sorted(wanted - files.keys())
    if missing:
        raise ValueError(f"missing local record files for {len(missing)} selected sources")
    output = []
    for candidate in selected:
        records = [json.loads(line) for line in files[candidate["source_family_id"]].read_text().splitlines() if line]
        output.append({
            "schema": "aleph.external_training.machine_preannotation_packet.v1",
            "authority_status": "proposed",
            "review_state": "machine_candidate",
            "human_verified": False,
            "source_family_id": candidate["source_family_id"],
            "source_payload_sha256": candidate["source_payload_sha256"],
            "candidate_selection_rank": candidate["selection_rank"],
            "records_file_sha256": sha256_bytes(files[candidate["source_family_id"]].read_bytes()),
            "experiment_records": [compact_record(record) for record in records],
            "atomic_claims": [claim for record in records for claim in machine_claims(record)],
        })
    return output


def blank_form(candidate: dict, slot: str) -> dict:
    assignment_id = "pilot:" + sha256_bytes(f"{candidate['source_family_id']}|{slot}".encode())[:16]
    return {
        "schema": "aleph.external_training.blinded_annotation_form.v1",
        "authority_status": "proposed",
        "review_state": "unassigned",
        "assignment_id": assignment_id,
        "annotation_slot": slot,
        "annotator_id": None,
        "source": {
            "source_family_id": candidate["source_family_id"],
            "pmcid": candidate["pmcid"],
            "title": candidate["title"],
            "publication_year": candidate["publication_year"],
            "source_payload_sha256": candidate["source_payload_sha256"],
        },
        "blinding": {
            "other_annotator_hidden": True,
            "machine_preannotation_hidden": True,
            "adjudication_hidden": True,
        },
        "claims": [],
        "submission": {"started_at": None, "completed_at": None, "notes": None},
    }


def build(args: argparse.Namespace) -> dict:
    candidates = read_jsonl(args.candidates)
    selected, selection_audit = select_fifty(candidates)
    preannotations = build_preannotations(selected, args.records_root)
    output = args.output
    output.mkdir(parents=True, exist_ok=True)

    candidate_rows = []
    for pilot_rank, candidate in enumerate(selected, 1):
        row = dict(candidate)
        row["pilot_rank"] = pilot_rank
        row["pilot_state"] = "annotation_pilot_candidate"
        row["human_verified"] = False
        candidate_rows.append(row)
    forms_a = [blank_form(row, "A") for row in selected]
    forms_b = [blank_form(row, "B") for row in selected]
    adjudication = [{
        "schema": "aleph.external_training.adjudication_form.v1",
        "authority_status": "proposed",
        "review_state": "awaiting_two_human_submissions",
        "source_family_id": row["source_family_id"],
        "submission_A_sha256": None,
        "submission_B_sha256": None,
        "adjudicator_id": None,
        "disagreements": [],
        "adjudicated_claims": [],
        "completed_at": None,
    } for row in selected]

    write_jsonl(output / "candidates_50.jsonl", candidate_rows)
    write_jsonl_gz(output / "machine_preannotations.jsonl.gz", preannotations)
    write_jsonl(output / "forms_A.jsonl", forms_a)
    write_jsonl(output / "forms_B.jsonl", forms_b)
    write_jsonl(output / "adjudication_forms.jsonl", adjudication)
    with (output / "assignments_100.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["assignment_id", "source_family_id", "slot", "annotator_id", "status"],
            lineterminator="\n",
        )
        writer.writeheader()
        for form in forms_a + forms_b:
            writer.writerow({
                "assignment_id": form["assignment_id"],
                "source_family_id": form["source"]["source_family_id"],
                "slot": form["annotation_slot"],
                "annotator_id": "",
                "status": "unassigned",
            })

    all_modalities = sorted({m for row in candidates for m in row["modality_tags"]})
    summary = {
        "schema": SCHEMA,
        "authority_status": "proposed",
        "human_annotations_completed": 0,
        "adjudications_completed": 0,
        "input_candidates": len(candidates),
        "selected_sources": len(selected),
        "assignments": len(forms_a) + len(forms_b),
        "source_candidate_sha256": sha256_bytes(args.candidates.read_bytes()),
        "domain_counts": dict(sorted(Counter(r["primary_domain"] for r in selected).items())),
        "year_stratum_counts": dict(sorted(Counter(r["year_stratum"] for r in selected).items())),
        "acquisition_pass_counts": dict(sorted(Counter(r["acquisition_pass"] for r in selected).items())),
        "modality_counts": dict(sorted(Counter(m for r in selected for m in r["modality_tags"]).items())),
        "cell_state_counts": dict(sorted(Counter(s for r in selected for s in r["cell_state_tags"]).items())),
        "modality_coverage": {
            "covered": sorted({m for r in selected for m in r["modality_tags"]}),
            "missing": sorted(set(all_modalities) - {m for r in selected for m in r["modality_tags"]}),
        },
        "selection_audit": selection_audit,
        "machine_preannotation_records": sum(len(row["experiment_records"]) for row in preannotations),
        "machine_preannotation_is_ground_truth": False,
        "forms_contain_machine_predictions": False,
    }
    (output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    readiness = {
        "schema": "aleph.external_training.publication_benchmark_readiness.v1",
        "authority_status": "proposed",
        "benchmark_metrics_available": False,
        "human_A_submissions": 0,
        "human_B_submissions": 0,
        "adjudicated_submissions": 0,
        "required_sources": 50,
        "reason": "No human annotation or adjudication has been completed; emitting precision, recall, agreement, or gold labels would be fabricated.",
        "next_valid_transition": "two distinct human submissions per source followed by human adjudication",
    }
    (output / "benchmark_readiness.json").write_text(json.dumps(readiness, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    artifacts = sorted(path for path in output.iterdir() if path.name not in {"SHA256SUMS.json", "audit.json"})
    manifest = {path.name: sha256_bytes(path.read_bytes()) for path in artifacts}
    (output / "SHA256SUMS.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--records-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(build(parse_args()), ensure_ascii=False, sort_keys=True))

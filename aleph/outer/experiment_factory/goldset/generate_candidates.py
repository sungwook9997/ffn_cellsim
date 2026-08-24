#!/usr/bin/env python
"""Build the deterministic, leakage-component-safe 300-paper annotation candidate set."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[4]
CORPUS = ROOT / "corpus" / "external_training"
HERE = Path(__file__).resolve().parent
VOCAB_PATH = HERE.parent / "schema" / "controlled_vocabularies.v1.json"
FIRST_QUEUE = CORPUS / "review_oa" / "results" / "manual_review_queue.csv"
SECOND_QUEUE = CORPUS / "review_oa" / "second_pass" / "results" / "manual_review_queue.csv"
FIRST_RECEIPTS = CORPUS / "acquisition" / "europe_pmc_oa_2026-08-05.jsonl"
SECOND_CANDIDATES = CORPUS / "acquisition" / "second_pass" / "candidates_2026-08-05.jsonl"
SECOND_RECEIPTS = CORPUS / "acquisition" / "second_pass" / "receipts_2026-08-05.jsonl"
LEAKAGE = CORPUS / "review_oa" / "second_pass" / "results" / "combined_dataset_leakage_groups.json"
SNAPSHOTS = CORPUS / "snapshots"
DOMAIN_TARGETS = {
    "adhesion_ecm_mechanics": 26,
    "cell_state_mechanics": 26,
    "cortical_membrane_pressure": 26,
    "cytoskeleton_motor_rheology": 26,
    "division_morphogenesis_tissue": 26,
    "flow_shear": 26,
    "measurement_inference": 26,
    "mechanotransduction": 26,
    "migration_invasion": 26,
    "nuclear_mechanics": 26,
    "organelle_general_mechanics": 20,
    "osmotic_volume_poroelasticity": 20,
}


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def split_tags(value: str) -> list[str]:
    return sorted(tag for tag in value.split(";") if tag)


def year_stratum(year: int) -> str:
    if year <= 2013:
        return "2010_2013"
    if year <= 2017:
        return "2014_2017"
    if year <= 2021:
        return "2018_2021"
    return "2022_2026"


def infer_historical_domain(row: dict[str, Any], metadata: dict[str, Any]) -> str:
    text = " ".join(
        [
            str(metadata.get("title") or ""),
            " ".join(row["modality_tags"]),
            " ".join(row["mechanics_tags"]),
            " ".join(row["cell_state_tags"]),
        ]
    ).lower()
    ordered_rules = [
        ("cortical_membrane_pressure", ("cortical tension", "membrane tension", "cortex", "membrane")),
        ("nuclear_mechanics", ("nucleus", "nuclear", "chromatin")),
        ("adhesion_ecm_mechanics", ("traction", "adhesion", "extracellular matrix", "ecm", "matrix stiffness")),
        ("osmotic_volume_poroelasticity", ("osmotic", "poroelastic", "cell_volume", "cell volume")),
        ("flow_shear", ("piv", "shear", "flow")),
        ("division_morphogenesis_tissue", ("division", "cytokinesis", "mitotic", "morphogen")),
        ("migration_invasion", ("migration", "migrating", "invasion", "protrusion")),
        ("cytoskeleton_motor_rheology", ("cytoskeleton", "actin", "microtubule", "motor", "viscosity")),
        ("mechanotransduction", ("mechanotransduction", "signalling", "signaling")),
        ("measurement_inference", ("afm", "immunofluorescence", "western_blot", "pcr", "microscopy", "tweezers")),
        ("cell_state_mechanics", ("differentiating", "apoptotic", "stem_progenitor", "cancer_tumor")),
    ]
    for domain, needles in ordered_rules:
        if any(needle in text for needle in needles):
            return domain
    return "organelle_general_mechanics"


class UnionFind:
    def __init__(self, values: set[str]):
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: str, right: str) -> None:
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[max(left, right)] = min(left, right)


def load_pool() -> tuple[list[dict[str, Any]], dict[str, set[str]], dict[str, int]]:
    rows: dict[str, dict[str, Any]] = {}
    for acquisition_pass, path in (("first", FIRST_QUEUE), ("second", SECOND_QUEUE)):
        with path.open(encoding="utf-8", newline="") as handle:
            for raw in csv.DictReader(handle):
                if raw["decision"] != "manual_priority":
                    continue
                source_id = raw["source_family_id"].lower()
                rows[source_id] = {
                    **raw,
                    "source_family_id": source_id,
                    "acquisition_pass": acquisition_pass,
                    "priority_score": int(raw["priority_score"]),
                    "observation_locator_count": int(raw["observation_locator_count"]),
                    "modality_tags": split_tags(raw["modality_tags"]),
                    "cell_state_tags": split_tags(raw["cell_state_tags"]),
                    "mechanics_tags": split_tags(raw["mechanics_tags"]),
                }

    first_receipts: dict[str, dict[str, Any]] = {}
    for receipt in read_jsonl(FIRST_RECEIPTS):
        if receipt.get("status") == "success":
            first_receipts[receipt["source_family_id"].lower()] = receipt
    metadata_by_path: dict[Path, dict[str, dict[str, Any]]] = {}
    for path in sorted(SNAPSHOTS.glob("*.jsonl")):
        metadata_by_path[path] = {row["source_family_id"].lower(): row for row in read_jsonl(path)}
    second_metadata = {row["source_family_id"].lower(): row for row in read_jsonl(SECOND_CANDIDATES)}
    second_receipts = {
        row["source_family_id"].lower(): row
        for row in read_jsonl(SECOND_RECEIPTS)
        if row.get("status") == "success"
    }

    vocab = json.loads(VOCAB_PATH.read_text(encoding="utf-8"))
    route_map = vocab["route_to_domain"]
    for source_id, row in rows.items():
        if row["acquisition_pass"] == "second":
            metadata = second_metadata[source_id]
            receipt = second_receipts[source_id]
        else:
            receipt = first_receipts[source_id]
            input_path = ROOT / receipt["candidate_input"]
            metadata = metadata_by_path[input_path][source_id]
        route = metadata["discovery_domain"]
        row["metadata"] = metadata
        row["publication_year"] = int(metadata["year"])
        row["year_stratum"] = year_stratum(row["publication_year"])
        row["primary_domain"] = route_map.get(route) or infer_historical_domain(row, metadata)
        row["source_payload_sha256"] = receipt["sha256"]

    leakage = json.loads(LEAKAGE.read_text(encoding="utf-8"))["groups"]
    all_grouped = {
        member["source_family_id"].lower()
        for group in leakage
        for member in group["members"]
    }
    graph = UnionFind(all_grouped)
    for group in leakage:
        members = sorted(member["source_family_id"].lower() for member in group["members"])
        for member in members[1:]:
            graph.union(members[0], member)
    full_components: dict[str, set[str]] = defaultdict(set)
    for source_id in all_grouped:
        full_components[graph.find(source_id)].add(source_id)

    eligible_ids = set(rows)
    eligible_components: dict[str, set[str]] = {}
    related_outside: dict[str, int] = {}
    for source_id in sorted(eligible_ids):
        full = full_components.get(graph.find(source_id), {source_id}) if source_id in graph.parent else {source_id}
        eligible = full & eligible_ids
        component_id = "lc:" + digest_bytes("\n".join(sorted(eligible)).encode("utf-8"))[:16]
        eligible_components[component_id] = eligible
        related_outside[component_id] = len(full - eligible_ids)
        rows[source_id]["leakage_component_id"] = component_id
    return list(rows.values()), eligible_components, related_outside


def choose_candidates(pool: list[dict[str, Any]], components: dict[str, set[str]]) -> list[dict[str, Any]]:
    by_id = {row["source_family_id"]: row for row in pool}
    selected: set[str] = set()
    domain_count: Counter[str] = Counter()
    pass_count: Counter[tuple[str, str]] = Counter()
    year_count: Counter[tuple[str, str]] = Counter()
    modality_count: Counter[tuple[str, str]] = Counter()
    state_count: Counter[tuple[str, str]] = Counter()

    def add(component_id: str, reason: str) -> None:
        for source_id in sorted(components[component_id]):
            if source_id in selected:
                continue
            row = by_id[source_id]
            row["selection_reason"] = reason
            selected.add(source_id)
            domain = row["primary_domain"]
            domain_count[domain] += 1
            pass_count[(domain, row["acquisition_pass"])] += 1
            year_count[(domain, row["year_stratum"])] += 1
            for tag in row["modality_tags"] or ["none_reported"]:
                modality_count[(domain, tag)] += 1
            for tag in row["cell_state_tags"] or ["none_reported"]:
                state_count[(domain, tag)] += 1

    # Preserve every multi-paper eligible accession component instead of selecting an easy
    # singleton-only set. The measured pool makes this 7 components / 17 papers.
    for component_id in sorted(components):
        members = components[component_id]
        if len(members) <= 1:
            continue
        prospective = Counter(by_id[source_id]["primary_domain"] for source_id in members)
        if all(domain_count[domain] + count <= DOMAIN_TARGETS[domain] for domain, count in prospective.items()):
            add(component_id, "accession_component_preserved")

    domains = json.loads(VOCAB_PATH.read_text(encoding="utf-8"))["domains"]
    for domain in domains:
        while domain_count[domain] < DOMAIN_TARGETS[domain]:
            candidates = [
                row for row in pool
                if row["source_family_id"] not in selected
                and row["primary_domain"] == domain
                and len(components[row["leakage_component_id"]]) == 1
            ]
            if not candidates:
                raise RuntimeError(f"cannot fill domain quota for {domain}")

            def key(row: dict[str, Any]) -> tuple[Any, ...]:
                modalities = row["modality_tags"] or ["none_reported"]
                states = row["cell_state_tags"] or ["none_reported"]
                stable = digest_bytes(row["source_family_id"].encode("utf-8"))
                return (
                    pass_count[(domain, row["acquisition_pass"])],
                    year_count[(domain, row["year_stratum"])],
                    min(modality_count[(domain, tag)] for tag in modalities),
                    min(state_count[(domain, tag)] for tag in states),
                    -sum(modality_count[(domain, tag)] == 0 for tag in modalities),
                    -sum(state_count[(domain, tag)] == 0 for tag in states),
                    -row["priority_score"],
                    stable,
                )

            winner = min(candidates, key=key)
            add(winner["leakage_component_id"], "domain_stratified_fill")

    chosen = [by_id[source_id] for source_id in selected]
    if len(chosen) != sum(DOMAIN_TARGETS.values()):
        raise AssertionError(f"expected 300 selected sources, got {len(chosen)}")
    return sorted(chosen, key=lambda row: (row["primary_domain"], row["source_family_id"]))


def counter_dict(values: list[str]) -> dict[str, int]:
    return dict(sorted(Counter(values).items()))


def main() -> int:
    pool, components, related_outside = load_pool()
    selected = choose_candidates(pool, components)
    records = []
    for rank, row in enumerate(selected, 1):
        source_material = {
            "queue": {key: value for key, value in row.items() if key not in {"metadata", "selection_reason"}},
            "metadata": row["metadata"],
        }
        records.append({
            "schema": "aleph.external_training.annotation_candidate.v1",
            "authority_status": "proposed",
            "review_state": "gold_candidate",
            "selection_rank": rank,
            "source_family_id": row["source_family_id"],
            "pmcid": row["pmcid"] or None,
            "title": row["metadata"]["title"],
            "journal": row["metadata"].get("journal"),
            "acquisition_pass": row["acquisition_pass"],
            "publication_year": row["publication_year"],
            "year_stratum": row["year_stratum"],
            "primary_domain": row["primary_domain"],
            "discovery_domain": row["metadata"]["discovery_domain"],
            "modality_tags": row["modality_tags"],
            "cell_state_tags": row["cell_state_tags"],
            "mechanics_tags": row["mechanics_tags"],
            "priority_score": row["priority_score"],
            "observation_locator_count": row["observation_locator_count"],
            "leakage_component_id": row["leakage_component_id"],
            "eligible_component_size": len(components[row["leakage_component_id"]]),
            "related_out_of_pool_count": related_outside[row["leakage_component_id"]],
            "selection_reason": row["selection_reason"],
            "source_record_sha256": digest_bytes(canonical_bytes(source_material)),
            "source_payload_sha256": row["source_payload_sha256"],
        })

    output_bytes = b"".join(canonical_bytes(record) for record in records)
    (HERE / "candidates_300.jsonl").write_bytes(output_bytes)
    assignment_lines = ["selection_rank,source_family_id,annotator_slot,assignment_status\n"]
    for record in records:
        for slot in ("A", "B"):
            assignment_lines.append(
                f"{record['selection_rank']},{record['source_family_id']},{slot},unassigned\n"
            )
    assignment_bytes = "".join(assignment_lines).encode("utf-8")
    (HERE / "annotation_assignments_600.csv").write_bytes(assignment_bytes)

    selected_ids = {row["source_family_id"] for row in selected}
    selected_components = {
        row["leakage_component_id"]
        for row in selected
        if len(components[row["leakage_component_id"]]) > 1
    }
    split_components = []
    for component_id, members in components.items():
        overlap = members & selected_ids
        if overlap and overlap != members:
            split_components.append(component_id)

    year_domain = Counter((row["primary_domain"], row["year_stratum"]) for row in selected)
    pass_domain = Counter((row["primary_domain"], row["acquisition_pass"]) for row in selected)
    unmet = []
    domains = json.loads(VOCAB_PATH.read_text(encoding="utf-8"))["domains"]
    for domain in domains:
        for acquisition_pass in ("first", "second"):
            if pass_domain[(domain, acquisition_pass)] == 0:
                unmet.append({"stratum": "domain_by_pass", "domain": domain, "value": acquisition_pass, "reason": "no eligible selected source; domain route may be pass-specific"})
        for stratum in ("2010_2013", "2014_2017", "2018_2021", "2022_2026"):
            if year_domain[(domain, stratum)] == 0:
                unmet.append({"stratum": "domain_by_year", "domain": domain, "value": stratum, "reason": "no selected source after domain quota and component preservation"})

    all_modalities = sorted({tag for row in pool for tag in row["modality_tags"]})
    all_states = sorted({tag for row in pool for tag in row["cell_state_tags"]})
    selected_modalities = {tag for row in selected for tag in row["modality_tags"]}
    selected_states = {tag for row in selected for tag in row["cell_state_tags"]}
    for tag in sorted(set(all_modalities) - selected_modalities):
        unmet.append({"stratum": "modality", "value": tag, "reason": "present in pool but absent from selection"})
    for tag in sorted(set(all_states) - selected_states):
        unmet.append({"stratum": "cell_state", "value": tag, "reason": "present in pool but absent from selection"})

    summary = {
        "schema": "aleph.external_training.annotation_candidate_summary.v1",
        "authority_status": "proposed",
        "review_state": "gold_candidate",
        "selection_algorithm": "all eligible multi-paper accession components, then deterministic domain-stratified singleton fill v1",
        "input_manual_priority": len(pool),
        "selected": len(records),
        "target_per_domain": DOMAIN_TARGETS,
        "counts": {
            "acquisition_pass": counter_dict([row["acquisition_pass"] for row in selected]),
            "publication_year": counter_dict([str(row["publication_year"]) for row in selected]),
            "year_stratum": counter_dict([row["year_stratum"] for row in selected]),
            "primary_domain": counter_dict([row["primary_domain"] for row in selected]),
            "modality_tag": counter_dict([tag for row in selected for tag in row["modality_tags"]]),
            "cell_state_tag": counter_dict([tag for row in selected for tag in row["cell_state_tags"]]),
        },
        "leakage": {
            "eligible_components": len(components),
            "eligible_multi_paper_components": sum(len(members) > 1 for members in components.values()),
            "eligible_papers_in_multi_paper_components": sum(len(members) for members in components.values() if len(members) > 1),
            "selected_multi_paper_components": len(selected_components),
            "selected_papers_in_multi_paper_components": sum(len(components[component_id]) for component_id in selected_components),
            "split_selected_components": split_components,
            "selected_with_related_out_of_pool": sum(record["related_out_of_pool_count"] > 0 for record in records),
        },
        "unmet_strata": unmet,
        "input_digests": {
            "first_queue_sha256": digest_bytes(FIRST_QUEUE.read_bytes()),
            "second_queue_sha256": digest_bytes(SECOND_QUEUE.read_bytes()),
            "combined_leakage_sha256": digest_bytes(LEAKAGE.read_bytes()),
            "vocabulary_sha256": digest_bytes(VOCAB_PATH.read_bytes()),
        },
        "candidate_set_sha256": digest_bytes(output_bytes),
        "limitations": [
            "gold_candidate is a sampling state, not a claim of human verification",
            "dataset-accession co-membership is a leakage precaution, not proof of shared measurements",
            "weak discovery domains and keyword tags are used only for stratified sampling",
            "sources related to an eligible component but outside manual_priority remain recorded as out-of-pool context",
        ],
    }
    (HERE / "summary.json").write_bytes(canonical_bytes(summary))
    hashes = {
        name: digest_bytes((HERE / name).read_bytes())
        for name in ("annotation_assignments_600.csv", "candidates_300.jsonl", "summary.json")
    }
    (HERE / "SHA256SUMS.json").write_bytes(canonical_bytes(hashes))
    print(json.dumps({"selected": len(records), "candidate_set_sha256": summary["candidate_set_sha256"], "unmet_strata": len(unmet)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

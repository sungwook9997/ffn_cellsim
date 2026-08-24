from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

from aleph.outer.experiment_factory.schema.validate import (
    validate_candidate,
    validate_experiment,
)


SHA = re.compile(r"^[0-9a-f]{64}$")
FORBIDDEN_RAW_KEYS = {
    "evidence_text", "caption_text", "raw_text", "full_text", "pixel_bytes",
    "image_bytes", "publisher_payload", "xml_payload", "pdf_payload",
}
FORBIDDEN_TRACKED_SUFFIXES = {
    ".pdf", ".xml", ".nxml", ".jpg", ".jpeg", ".png", ".gif", ".tif",
    ".tiff", ".bmp", ".webp", ".mp4", ".avi", ".npy", ".h5",
    ".hdf5", ".parquet",
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def rows(path: Path) -> Iterable[dict[str, Any]]:
    if path.suffix == ".json":
        yield json.loads(path.read_text(encoding="utf-8"))
        return
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def walk(value: Any) -> Iterable[tuple[str, Any]]:
    if isinstance(value, dict):
        for key, child in value.items():
            yield key, child
            yield from walk(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk(child)


def duplicate_count(values: Iterable[str]) -> int:
    material = list(values)
    return len(material) - len(set(material))


def verify_sum_file(root: Path, sum_path: Path, failures: list[str]) -> dict[str, Any]:
    expected = json.loads(sum_path.read_text(encoding="utf-8"))
    mismatches = []
    for relative, wanted in expected.items():
        target = sum_path.parent / relative
        got = digest(target) if target.is_file() else None
        if got != wanted:
            mismatches.append({"path": str(target.relative_to(root)), "expected": wanted, "actual": got})
    if mismatches:
        failures.append(f"digest_mismatch:{sum_path}")
    return {"manifest": str(sum_path.relative_to(root)), "files": len(expected), "mismatches": mismatches}


def audit(args: argparse.Namespace) -> dict[str, Any]:
    root = args.root.resolve()
    factory = root / "aleph/outer/experiment_factory"
    failures: list[str] = []
    warnings: list[str] = []

    digest_manifests = []
    for relative in (
        "goldset/SHA256SUMS.json", "datasets/SHA256SUMS.json", "visual/results/SHA256SUMS.json",
    ):
        digest_manifests.append(verify_sum_file(root, factory / relative, failures))

    extraction_audit = json.loads((factory / "extraction/results/audit.json").read_text())
    extraction_digests = {
        "article_index": digest(factory / "extraction/results/article_index.jsonl"),
        "figure_join_index": digest(factory / "extraction/results/figure_join_index.jsonl"),
    }
    if extraction_digests["article_index"] != extraction_audit["article_index_sha256"]:
        failures.append("extraction_article_index_digest")
    if extraction_digests["figure_join_index"] != extraction_audit["figure_join_index_sha256"]:
        failures.append("extraction_figure_join_digest")

    learning_path = factory / "learning/results/source_projection.jsonl"
    learning_audit_path = factory / "learning/results/audit.json"
    learning_present = learning_path.is_file() and learning_audit_path.is_file()
    learning_digest = None
    if learning_present:
        learning_audit = json.loads(learning_audit_path.read_text())
        learning_digest = digest(learning_path)
        if learning_digest != learning_audit.get("projection_file_sha256"):
            failures.append("learning_projection_digest")

    model_results = factory / "model/results"
    model_present = (model_results / "SHA256SUMS.json").is_file()
    model_summary: dict[str, Any] = {"present": model_present}
    if model_present:
        digest_manifests.append(verify_sum_file(root, model_results / "SHA256SUMS.json", failures))
        model_config = json.loads((model_results / "config.json").read_text())
        model_metrics = json.loads((model_results / "metrics.json").read_text())
        model_split = json.loads((model_results / "split_audit.json").read_text())
        model_gold = json.loads((model_results / "gold_candidate_representation_audit.json").read_text())
        model_errors = []
        if learning_digest is None or model_config.get("projection_sha256") != learning_digest or model_split.get("projection_file_sha256") != learning_digest:
            model_errors.append("projection_binding")
        if model_metrics.get("target") != "weak_12_class_discovery_route_only":
            model_errors.append("target_scope")
        if model_config.get("test_used_for_model_selection") is not False or model_config.get("gold_candidate_used_for_model_selection") is not False:
            model_errors.append("selection_misuse")
        if any(model_gold.get(key) is not False for key in ("performance_metrics_computed", "used_for_training", "used_for_validation_model_selection", "used_for_test_metrics")):
            model_errors.append("gold_cohort_misuse")
        if model_gold.get("sources") != 300 or model_split.get("leakage_split_violations") != 0 or model_split.get("errors") != 0:
            model_errors.append("split_or_cohort")
        view_seeds = []
        test_supports = []
        for view in model_metrics.get("views", {}).values():
            view_seeds.append([row.get("seed") for row in view.get("per_seed", [])])
            test_supports.extend(row.get("test", {}).get("support") for row in view.get("per_seed", []))
        if not view_seeds or any(seeds != model_metrics.get("seeds") for seeds in view_seeds):
            model_errors.append("seed_mismatch")
        if len(set(test_supports)) != 1:
            model_errors.append("test_support_mismatch")
        if model_errors:
            failures.extend(f"model_{error}" for error in model_errors)
        model_summary = {
            "present": True,
            "projection_sha256": model_config.get("projection_sha256"),
            "target": model_metrics.get("target"),
            "views": sorted(model_metrics.get("views", {})),
            "seeds": model_metrics.get("seeds"),
            "test_support": test_supports[0] if test_supports else None,
            "gold_sources": model_gold.get("sources"),
            "gold_used_for_training": model_gold.get("used_for_training"),
            "gold_performance_metrics_computed": model_gold.get("performance_metrics_computed"),
            "leakage_split_violations": model_split.get("leakage_split_violations"),
            "validation_errors": model_errors,
        }
        if args.model_replay_results is not None:
            replay_manifest = json.loads((args.model_replay_results / "SHA256SUMS.json").read_text())
            production_manifest = json.loads((model_results / "SHA256SUMS.json").read_text())
            replay_mismatches = sorted(name for name in production_manifest if production_manifest[name] != replay_manifest.get(name))
            model_summary["replay_artifact_mismatches"] = replay_mismatches
            model_summary["replay_artifact_hashes_equal"] = not replay_mismatches
            if replay_mismatches:
                failures.append("model_byte_replay_mismatch")

    gold = list(rows(factory / "goldset/candidates_300.jsonl"))
    gold_sources = {row["source_family_id"] for row in gold}
    gold_components: dict[str, set[str]] = defaultdict(set)
    for row in gold:
        gold_components[row["leakage_component_id"]].add(row["source_family_id"])
    if len(gold) != 300 or len(gold_sources) != 300:
        failures.append("frozen_gold_cohort_not_300_unique_sources")
    if duplicate_count(str(row["selection_rank"]) for row in gold):
        failures.append("duplicate_gold_selection_rank")
    if any(row.get("authority_status") != "proposed" or row.get("review_state") != "gold_candidate" for row in gold):
        failures.append("gold_authority_or_review_state")
    candidate_schema_errors = sum(len(validate_candidate(row)) for row in gold)
    if candidate_schema_errors:
        failures.append("candidate_schema_validation")

    local_record_root = root / "data/external_training/experiment_factory/extraction/full/records"
    experiment_schema_records = 0
    experiment_schema_errors = 0
    if not local_record_root.is_dir():
        failures.append("local_experiment_records_unavailable")
    else:
        for path in sorted(local_record_root.glob("*.jsonl")):
            for record in rows(path):
                experiment_schema_records += 1
                experiment_schema_errors += len(validate_experiment(record))
    if experiment_schema_errors:
        failures.append("experiment_schema_validation")

    article_index = list(rows(factory / "extraction/results/article_index.jsonl"))
    figure_index = list(rows(factory / "extraction/results/figure_join_index.jsonl"))
    if duplicate_count(row["source_family_id"] for row in article_index):
        failures.append("duplicate_extraction_source")
    if duplicate_count(row["experiment_record_id"] for row in figure_index):
        failures.append("duplicate_figure_join_experiment_id")

    accessions = list(rows(factory / "datasets/source_accession_registry.jsonl"))
    manifests = list(rows(factory / "datasets/dataset_manifests.jsonl"))
    accession_keys = [f'{r["source_family_id"]}\x1f{r["provider"]}\x1f{r["accession"]}' for r in accessions]
    manifest_ids = [row["manifest_id"] for row in manifests]
    if duplicate_count(accession_keys):
        failures.append("duplicate_source_accession_binding")
    if duplicate_count(manifest_ids):
        failures.append("duplicate_dataset_manifest")

    visual_index = list(rows(factory / "visual/results/visual_index.jsonl.gz"))
    visual_receipts = list(rows(factory / "visual/results/asset_receipts.jsonl"))
    visual_descriptors = list(rows(factory / "visual/results/descriptors.jsonl.gz"))
    visual_ids = [row["observation_candidate_id"] for row in visual_index]
    receipt_ids = [row["asset_candidate_id"] for row in visual_receipts]
    descriptor_ids = [row["asset_candidate_id"] for row in visual_descriptors]
    if duplicate_count(visual_ids):
        failures.append("duplicate_visual_observation")
    if duplicate_count(receipt_ids):
        failures.append("duplicate_visual_receipt")
    if duplicate_count(descriptor_ids):
        failures.append("duplicate_visual_descriptor")
    if set(receipt_ids) != set(descriptor_ids):
        failures.append("visual_receipt_descriptor_identity_mismatch")
    referenced_objects = {row["sha256"] for row in visual_receipts if row.get("asset_state") == "acquired"}
    object_root = root / "data/external_training/experiment_factory/visual/objects"
    object_missing = 0
    object_digest_errors = 0
    for object_sha in sorted(referenced_objects):
        target = object_root / object_sha
        if not target.is_file():
            object_missing += 1
        elif digest(target) != object_sha:
            object_digest_errors += 1
    if object_missing or object_digest_errors:
        failures.append("visual_local_object_integrity")

    projection = list(rows(learning_path)) if learning_present else []
    projection_ids = [row["projection_id"] for row in projection]
    projection_sources = [row["source_family_id"] for row in projection]
    if projection and (duplicate_count(projection_ids) or duplicate_count(projection_sources)):
        failures.append("duplicate_learning_projection_identity")
    group_splits: dict[str, set[str]] = defaultdict(set)
    source_split = {}
    cohort_sources = set()
    forbidden_learning_fields = 0
    ambiguous_promotions = 0
    for row in projection:
        group_splits[row["leakage_group_id"]].add(row["split"])
        source_split[row["source_family_id"]] = row["split"]
        if row["evaluation_cohort"] == "frozen_300_gold_candidate_not_ground_truth":
            cohort_sources.add(row["source_family_id"])
        policy = row.get("supervision_policy", {})
        if policy.get("gold_candidate_is_ground_truth") or policy.get("ambiguous_observations_are_labels"):
            ambiguous_promotions += 1
        forbidden_learning_fields += sum(key in FORBIDDEN_RAW_KEYS for key, _ in walk(row))
    leakage_split_violations = sum(len(splits) > 1 for splits in group_splits.values())
    if leakage_split_violations:
        failures.append("learning_leakage_group_crosses_split")
    if projection and cohort_sources != gold_sources:
        failures.append("frozen_300_cohort_mismatch")
    if ambiguous_promotions:
        failures.append("gold_or_ambiguous_label_promotion")
    if forbidden_learning_fields:
        failures.append("raw_fields_in_learning_projection")
    accession_split_violations = 0
    for component in {row["leakage_component_id"] for row in accessions}:
        splits = {source_split[row["source_family_id"]] for row in accessions if row["leakage_component_id"] == component and row["source_family_id"] in source_split}
        accession_split_violations += len(splits) > 1
    if accession_split_violations:
        failures.append("dataset_accession_component_crosses_split")
    gold_split_violations = 0
    for sources in gold_components.values():
        splits = {source_split[source] for source in sources if source in source_split}
        gold_split_violations += len(splits) > 1
    if gold_split_violations:
        failures.append("gold_leakage_component_crosses_split")

    committed_data_files = subprocess.run(
        ["git", "ls-files", "aleph/outer/experiment_factory", "data/external_training/experiment_factory"],
        cwd=root, check=True, text=True, stdout=subprocess.PIPE,
    ).stdout.splitlines()
    forbidden_payload_files = [path for path in committed_data_files if Path(path).suffix.casefold() in FORBIDDEN_TRACKED_SUFFIXES]
    tracked_data_payloads = [path for path in committed_data_files if path.startswith("data/")]
    if forbidden_payload_files or tracked_data_payloads:
        failures.append("raw_or_pilot_payload_tracked")

    structured_files = [
        factory / "goldset/candidates_300.jsonl",
        factory / "extraction/results/article_index.jsonl",
        factory / "extraction/results/figure_join_index.jsonl",
        factory / "datasets/source_accession_registry.jsonl",
        factory / "datasets/dataset_manifests.jsonl",
        factory / "visual/results/visual_index.jsonl.gz",
        factory / "visual/results/asset_receipts.jsonl",
        factory / "visual/results/descriptors.jsonl.gz",
    ]
    if learning_present:
        structured_files.append(learning_path)
    if model_present:
        structured_files.extend(
            model_results / name for name in (
                "config.json", "feature_schemas.json", "gold_candidate_representation_audit.json",
                "metrics.json", "split_audit.json",
            )
        )
    forbidden_raw_key_hits = 0
    authority_values = Counter()
    decided_by_hits = 0
    for path in structured_files:
        for row in rows(path):
            for key, value in walk(row):
                forbidden_raw_key_hits += key in FORBIDDEN_RAW_KEYS
                decided_by_hits += key == "decided_by"
                if key == "authority_status":
                    authority_values[str(value)] += 1
    if forbidden_raw_key_hits:
        failures.append("committed_raw_text_or_pixel_key")
    if decided_by_hits or set(authority_values) - {"proposed"}:
        failures.append("authority_promotion_or_decided_by")

    return {
        "schema": "aleph.external_training.experiment_factory_independent_audit.v1",
        "authority_status": "proposed",
        "audited_repository_head": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, check=True, text=True, stdout=subprocess.PIPE,
        ).stdout.strip(),
        "status": "pass" if not failures else "fail",
        "failures": failures,
        "warnings": warnings,
        "digest_manifests": digest_manifests,
        "extraction": {
            "articles": len(article_index), "figure_joins": len(figure_index),
            "article_index_sha256": extraction_digests["article_index"],
            "figure_join_index_sha256": extraction_digests["figure_join_index"],
        },
        "schema_validation": {
            "candidate_records": len(gold), "candidate_errors": candidate_schema_errors,
            "experiment_records": experiment_schema_records, "experiment_errors": experiment_schema_errors,
        },
        "goldset": {"records": len(gold), "unique_sources": len(gold_sources), "duplicate_ranks": duplicate_count(str(r["selection_rank"]) for r in gold)},
        "datasets": {
            "source_accession_bindings": len(accessions), "manifests": len(manifests),
            "duplicate_bindings": duplicate_count(accession_keys), "duplicate_manifests": duplicate_count(manifest_ids),
            "tracked_pilot_payloads": len(tracked_data_payloads),
        },
        "visual": {
            "figure_records": len(visual_index), "receipts": len(visual_receipts), "descriptors": len(visual_descriptors),
            "duplicate_receipts": duplicate_count(receipt_ids), "receipt_descriptor_id_set_equal": set(receipt_ids) == set(descriptor_ids),
            "referenced_unique_objects": len(referenced_objects), "local_object_missing": object_missing,
            "local_object_digest_errors": object_digest_errors,
        },
        "learning": {
            "present": learning_present, "projection_sha256": learning_digest, "sources": len(projection),
            "duplicate_projection_ids": duplicate_count(projection_ids), "leakage_groups": len(group_splits),
            "leakage_split_violations": leakage_split_violations,
            "accession_split_violations": accession_split_violations,
            "gold_component_split_violations": gold_split_violations,
            "frozen_300_cohort_sources": len(cohort_sources), "cohort_exact_match": cohort_sources == gold_sources if projection else None,
            "ambiguous_or_gold_promotions": ambiguous_promotions,
        },
        "model": model_summary,
        "privacy_authority": {
            "structured_files_scanned": len(structured_files), "forbidden_raw_key_hits": forbidden_raw_key_hits,
            "forbidden_payload_files": forbidden_payload_files, "decided_by_hits": decided_by_hits,
            "authority_status_values": dict(sorted(authority_values.items())),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model-replay-results", type=Path)
    args = parser.parse_args()
    result = audit(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    encoded = (json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n").encode()
    args.output.write_bytes(encoded)
    print(json.dumps({"status": result["status"], "failures": len(result["failures"]), "sha256": hashlib.sha256(encoded).hexdigest()}, sort_keys=True))
    raise SystemExit(0 if result["status"] == "pass" else 1)


if __name__ == "__main__":
    main()

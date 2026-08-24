#!/usr/bin/env python
"""Train combined first+second OA content router with accession-group split constraints."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
V2_SPEC = importlib.util.spec_from_file_location("oa_content_v2", HERE.parent / "oa_content" / "train.py")
v2 = importlib.util.module_from_spec(V2_SPEC)
V2_SPEC.loader.exec_module(v2)
v1 = v2.v1


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, value):
        self.parent.setdefault(value, value)
        if self.parent[value] != value:
            self.parent[value] = self.find(self.parent[value])
        return self.parent[value]

    def union(self, left, right):
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def build_component_keys(group_path: Path) -> tuple[dict[str, str], dict]:
    data = json.loads(group_path.read_text())
    uf = UnionFind()
    raw_groups = 0
    for group in data["groups"]:
        members = sorted({v1.canonical_source_family(row["source_family_id"]) for row in group["members"]})
        if len(members) < 2:
            continue
        raw_groups += 1
        for member in members[1:]:
            uf.union(members[0], member)
    components = defaultdict(list)
    for family in uf.parent:
        components[uf.find(family)].append(family)
    keys = {}
    for members in components.values():
        members.sort()
        key = "dataset:" + hashlib.sha256("\n".join(members).encode()).hexdigest()[:20]
        for family in members:
            keys[family] = key
    return keys, {
        "raw_accession_groups": raw_groups,
        "connected_components": len(components),
        "grouped_source_families": len(keys),
        "largest_component": max((len(members) for members in components.values()), default=0),
        "multi_group_components": sum(len(members) > 2 for members in components.values()),
    }


def load_metadata(repo_root: Path, config: dict, taxonomy: dict) -> tuple[dict, dict]:
    first_records, _ = v1.load_records(repo_root / "corpus" / "external_training")
    metadata = {}
    for record in first_records:
        weak = record["label"]
        if weak not in taxonomy["first_pass_mapping"]:
            raise RuntimeError(f"unmapped_first_label:{weak}")
        item = dict(record)
        item.update({"pass": "first", "weak_label": weak, "label": taxonomy["first_pass_mapping"][weak]})
        metadata[("first", item["source_family_id"])] = item
    second_count = 0
    cross_pass_metadata_overlap = 0
    with (repo_root / config["second_candidates"]).open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            family = v1.canonical_source_family(record["source_family_id"])
            weak = record["discovery_domain"]
            if weak not in taxonomy["second_pass_mapping"]:
                raise RuntimeError(f"unmapped_second_label:{weak}")
            cross_pass_metadata_overlap += int(("first", family) in metadata)
            metadata[("second", family)] = {
                "source_family_id": family,
                "title": record.get("title") or "",
                "journal": record.get("journal") or "",
                "year": int(record.get("year") or 0),
                "pass": "second",
                "weak_label": weak,
                "label": taxonomy["second_pass_mapping"][weak],
            }
            second_count += 1
    return metadata, {"first_metadata": len(first_records), "second_metadata": second_count, "cross_pass_metadata_overlap": cross_pass_metadata_overlap}


def load_quality_families(repo_root: Path, config: dict) -> tuple[dict[str, set[str]], dict]:
    result = {}
    counts = {}
    for pass_name, relative in config["quality_indices"].items():
        families = set()
        with (repo_root / relative).open(encoding="utf-8") as handle:
            for line in handle:
                families.add(v1.canonical_source_family(json.loads(line)["source_family_id"]))
        result[pass_name] = families
        counts[pass_name] = len(families)
    return result, counts


def load_combined_records(repo_root: Path, config: dict, metadata: dict, quality: dict) -> tuple[list[dict], dict]:
    records = {}
    root_audit = {}
    selected = tuple(config["features"]["section_token_limits"])
    for pass_name, relative_root in config["derived_roots"].items():
        paths = sorted((repo_root / relative_root).glob("*.chunks.jsonl"))
        counters = Counter()
        for path in paths:
            with path.open(encoding="utf-8") as handle:
                rows = [json.loads(line) for line in handle if line.strip()]
            if not rows:
                counters["empty_files"] += 1
                continue
            families = {v1.canonical_source_family(row.get("source_family_id", "")) for row in rows}
            if len(families) != 1:
                counters["invalid_family_files"] += 1
                continue
            family = next(iter(families))
            metadata_key = (pass_name, family)
            if metadata_key not in metadata:
                counters["missing_or_wrong_metadata"] += 1
                continue
            if family not in quality[pass_name]:
                counters["missing_quality_record"] += 1
                continue
            section_tokens = {section: [] for section in selected}
            for row in rows:
                section = row.get("section_class")
                if section in section_tokens and row.get("text"):
                    section_tokens[section].extend(v1.TOKEN_RE.findall(row["text"].lower()))
            for section in selected:
                section_tokens[section] = section_tokens[section][:config["features"]["section_token_limits"][section]]
                if section_tokens[section]:
                    counters[f"with_{section}"] += 1
            if not any(section_tokens.values()):
                counters["no_target_text"] += 1
                continue
            if family in records:
                raise RuntimeError(f"cross_root_derived_overlap:{family}")
            item = dict(metadata[metadata_key])
            item.update({"section_tokens": section_tokens, "derived_file": path.name})
            records[family] = item
            counters["usable"] += 1
        counters["chunk_files"] = len(paths)
        root_audit[pass_name] = dict(sorted(counters.items()))
    return [records[key] for key in sorted(records)], root_audit


def split_for_key(key: str, config: dict) -> str:
    value = int.from_bytes(hashlib.sha256(f"{config['seed']}:{key}".encode()).digest()[:8], "big") / 2**64
    if value < config["split"]["train"]:
        return "train"
    if value < config["split"]["train"] + config["split"]["validation"]:
        return "validation"
    return "test"


def run(config_path: Path) -> dict:
    started = time.perf_counter()
    config = json.loads(config_path.read_text())
    repo_root = config_path.resolve().parents[4]
    taxonomy_path = repo_root / config["taxonomy"]
    taxonomy = json.loads(taxonomy_path.read_text())
    metadata, metadata_audit = load_metadata(repo_root, config, taxonomy)
    quality, quality_counts = load_quality_families(repo_root, config)
    records, root_audit = load_combined_records(repo_root, config, metadata, quality)
    component_keys, leakage_audit = build_component_keys(repo_root / config["combined_leakage_groups"])
    for record in records:
        record["group_key"] = component_keys.get(record["source_family_id"], "source:" + record["source_family_id"])
        record["split"] = split_for_key(record["group_key"], config)
    labels = taxonomy["coarse_classes"]
    observed = {record["label"] for record in records}
    if observed != set(labels):
        raise RuntimeError(f"taxonomy_support_mismatch:{sorted(set(labels)-observed)}")
    label_to_id = {label: i for i, label in enumerate(labels)}
    y = np.asarray([label_to_id[record["label"]] for record in records], dtype=np.int64)
    splits = {name: np.asarray([i for i, record in enumerate(records) if record["split"] == name], dtype=np.int64) for name in ("train", "validation", "test")}
    missing = {name: sorted(set(labels) - {records[i]["label"] for i in indices}) for name, indices in splits.items()}
    if any(missing.values()):
        raise RuntimeError(f"split_missing_coarse_support:{missing}")
    title_features = [v2.feature_ids(record, config, False) for record in records]
    content_features = [v2.feature_ids(record, config, True) for record in records]
    title_weights, title_metrics = v2.fit(title_features, y, splits, labels, config, 0)
    content_weights, content_metrics = v2.fit(content_features, y, splits, labels, config, 1)
    output = config_path.parent
    np.savez_compressed(output / "weights_title_only.npz", w1=title_weights[0], b1=title_weights[1], w2=title_weights[2], b2=title_weights[3], labels=np.asarray(labels))
    np.savez_compressed(output / "weights_combined_content.npz", w1=content_weights[0], b1=content_weights[1], w2=content_weights[2], b2=content_weights[3], labels=np.asarray(labels))
    with (output / "splits.jsonl").open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps({key: record[key] for key in ("source_family_id", "pass", "weak_label", "label", "group_key", "split", "derived_file")}, sort_keys=True) + "\n")
    grouped_usable = [record for record in records if record["source_family_id"] in component_keys]
    groups_to_splits = defaultdict(set)
    for record in records:
        groups_to_splits[record["group_key"]].add(record["split"])
    mixed = defaultdict(set)
    for record in grouped_usable:
        mixed[record["group_key"]].add(record["label"])
    metrics = {
        "schema_version": "1.0.0",
        "authority": config["authority"],
        "limitations": ["Coarse classes map incompatible weak discovery routes; they are routing labels, not biological truth.", "Dataset-accession co-membership is a conservative leakage constraint, not proof of dataset identity.", "Quality records are proposed recoverability routes and provide no training authority or Tier A promotion.", "Metrics are not directly comparable to v2's 29-label task."],
        "input": {"metadata": metadata_audit, "quality_index_families": quality_counts, "derived": root_audit, "usable_unique_families": len(records), "first_pass_manifest_aggregate_sha256": "d7fcb6ebbc2156aa85d49e087bbbf1f6a54ddc370c8b32c5f5535403ae17c9d2", "second_pass_manifest_set_sha256": "9627db89509f0a2f363ec93e65a2da4b68d3ab1b0415195788145b54a30bcf0f"},
        "taxonomy": {"classes": labels, "sha256": sha256_file(taxonomy_path)},
        "leakage": {**leakage_audit, "usable_grouped_families": len(grouped_usable), "usable_mixed_label_components": sum(len(values)>1 for values in mixed.values()), "group_split_violations": sum(len(values)>1 for values in groups_to_splits.values()), "combined_artifact_sha256": sha256_file(repo_root / config["combined_leakage_groups"])},
        "splits": {name: len(indices) for name, indices in splits.items()},
        "label_counts_by_split": {name: dict(sorted(Counter(records[i]["label"] for i in indices).items())) for name, indices in splits.items()},
        "title_only": title_metrics,
        "combined_content_plus_metadata": content_metrics,
        "delta": {"macro_f1": content_metrics["macro_f1"]-title_metrics["macro_f1"], "accuracy": content_metrics["accuracy"]-title_metrics["accuracy"], "retrieval_hit_rate_at_5": content_metrics["retrieval"]["hit_rate_at_k"]-title_metrics["retrieval"]["hit_rate_at_k"], "retrieval_precision_at_5": content_metrics["retrieval"]["precision_at_k"]-title_metrics["retrieval"]["precision_at_k"]},
        "v2_coverage_comparison_only": {"v2_usable_families": 2695, "v2_labels": 29, "v3_usable_families": len(records), "v3_coarse_labels": len(labels), "metric_comparison_permitted": False},
        "pipeline_runtime_seconds": time.perf_counter()-started,
    }
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    names = ["config.json", "taxonomy.json", "train.py", "splits.jsonl", "metrics.json", "weights_title_only.npz", "weights_combined_content.npz"]
    (output / "SHA256SUMS.json").write_text(json.dumps({name: sha256_file(output/name) for name in names}, indent=2, sort_keys=True) + "\n")
    return metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    metrics = run(args.config.resolve())
    print(json.dumps({"usable": metrics["input"]["usable_unique_families"], "splits": metrics["splits"], "leakage": metrics["leakage"], "title_macro_f1": metrics["title_only"]["macro_f1"], "content_macro_f1": metrics["combined_content_plus_metadata"]["macro_f1"], "delta": metrics["delta"]}, indent=2))


if __name__ == "__main__":
    main()

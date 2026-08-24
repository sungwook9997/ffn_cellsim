#!/usr/bin/env python3
"""Build deterministic leakage-safe, text-free multimodal learning views."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

SCHEMA = "aleph.external_training.learning_source_projection.v1"
VERSION = "1"
DESCRIPTOR_FEATURES = (
    "width", "height", "aspect_ratio", "decoded_channels", "colorfulness_channel_std_mean",
    "edge_density_absdiff_gt_24", "entropy_bits", "gray_mean", "gray_std", "source_frame_count",
)


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_digest(path: Path | None) -> str | None:
    return digest(path.read_bytes()) if path and path.is_file() else None


def atomic(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        handle.write(data)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def rows(path: Path | None) -> Iterable[dict[str, Any]]:
    if path is None or not path.is_file():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


class UnionFind:
    def __init__(self, values: Iterable[str]):
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: str, right: str) -> None:
        a, b = self.find(left), self.find(right)
        if a != b:
            low, high = sorted((a, b))
            self.parent[high] = low


def split_for_group(group_id: str) -> str:
    bucket = int(digest(group_id.encode())[:8], 16) % 100
    return "train" if bucket < 70 else ("validation" if bucket < 85 else "test")


def numeric_views(records_root: Path | None, local_output: Path | None) -> dict[str, dict[str, Any]]:
    by_source: dict[str, dict[str, Any]] = {}
    local_rows = []
    if records_root is None or not records_root.is_dir():
        return by_source
    for path in sorted(records_root.glob("*.jsonl")):
        for record in rows(path):
            source = record["source"]["source_family_id"].lower()
            state = by_source.setdefault(source, {
                "record_count": 0, "observation_count": 0, "candidate_target_count": 0,
                "ambiguous_target_mask_count": 0, "quantities": Counter(), "roles": Counter(),
                "record_ids": [],
            })
            state["record_count"] += 1
            state["record_ids"].append(record["experiment_record_id"])
            context_ambiguous = (
                all(system["cell_type"]["normalization_status"] in {"missing", "ambiguous", "unmapped"}
                    for system in record["biological_systems"])
                or all(protocol["modality"] == "unknown" for protocol in record["protocols"])
            )
            for observation in record["observations"]:
                state["observation_count"] += 1
                role, quantity = observation["measured_entity"].split(":", 1)
                state["roles"][role] += 1
                state["quantities"][quantity] += 1
                ambiguous = observation["extraction"]["status"] == "ambiguous" or context_ambiguous
                if ambiguous:
                    state["ambiguous_target_mask_count"] += 1
                else:
                    state["candidate_target_count"] += 1
                local_rows.append({
                    "schema": "aleph.external_training.learning_numeric_view.v1",
                    "authority_status": "proposed",
                    "observation_id": observation["observation_id"],
                    "experiment_record_id": record["experiment_record_id"],
                    "source_family_id": source,
                    "measured_entity": observation["measured_entity"],
                    "value": observation["value"],
                    "unit": observation["unit"],
                    "target_mask": not ambiguous,
                    "ambiguity_mask": ambiguous,
                    "context_ambiguity_mask": context_ambiguous,
                    "label_authority": "proposed_candidate_only",
                    "evidence_digests": [item["text_or_asset_sha256"] for item in observation["evidence"]],
                })
    for source, state in by_source.items():
        state["quantities"] = dict(sorted(state["quantities"].items()))
        state["roles"] = dict(sorted(state["roles"].items()))
        state["record_set_sha256"] = digest(canonical(sorted(state.pop("record_ids"))))
    if local_output is not None:
        data = b"".join(canonical(row) + b"\n" for row in sorted(local_rows, key=lambda row: row["observation_id"]))
        atomic(local_output, data)
    return by_source


def build(article_index: Path, output: Path, summary_path: Path, *,
          gold_candidates: Path | None = None, records_root: Path | None = None,
          local_numeric_output: Path | None = None, visual_index: Path | None = None,
          visual_links: Path | None = None, asset_receipts: Path | None = None,
          visual_descriptors: Path | None = None,
          dataset_registry: Path | None = None, dataset_manifests: Path | None = None) -> dict[str, Any]:
    articles = {row["source_family_id"].lower(): row for row in rows(article_index)}
    sources = sorted(articles)
    gold = {row["source_family_id"].lower(): row for row in rows(gold_candidates)}
    numeric = numeric_views(records_root, local_numeric_output)

    visual: dict[str, dict[str, Any]] = defaultdict(lambda: {
        "figure_count": 0, "linked_figure_count": 0, "available_asset_count": 0,
        "modality_tags": Counter(), "visual_candidate_ids": [],
    })
    for row in rows(visual_index):
        source = row["source_family_id"].lower()
        if source not in articles:
            continue
        visual[source]["figure_count"] += 1
        visual[source]["visual_candidate_ids"].append(row["observation_candidate_id"])
        visual[source]["modality_tags"].update(item["tag"] for item in row.get("proposed_modality_tags", []))
    for row in rows(visual_links):
        source = row["source_family_id"].lower()
        if source in articles and row.get("linked_experiment_record_ids"):
            visual[source]["linked_figure_count"] += 1
    for row in rows(asset_receipts):
        source = row["source_family_id"].lower()
        if source in articles and row.get("sha256"):
            visual[source]["available_asset_count"] += 1
    descriptors: dict[str, dict[str, Any]] = defaultdict(lambda: {
        "descriptor_count": 0, "descriptor_refs": [], "feature_values": defaultdict(list),
        "pixel_format_counts": Counter(),
    })
    for row in rows(visual_descriptors):
        source = row.get("source_family_id", "").lower()
        if source not in articles or row.get("descriptor_state") != "described" or not isinstance(row.get("descriptor"), dict):
            continue
        state = descriptors[source]
        state["descriptor_count"] += 1
        state["descriptor_refs"].append({
            "descriptor_id": row["descriptor_id"],
            "descriptor_sha256": digest(canonical(row["descriptor"])),
            "asset_sha256": row.get("asset_sha256"),
            "observation_candidate_id": row.get("observation_candidate_id"),
        })
        descriptor = row["descriptor"]
        for feature in DESCRIPTOR_FEATURES:
            value = descriptor.get(feature)
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                state["feature_values"][feature].append(float(value))
        quantiles = descriptor.get("gray_quantiles_05_25_50_75_95")
        if isinstance(quantiles, list) and len(quantiles) == 5:
            for label, value in zip(("gray_q05", "gray_q25", "gray_q50", "gray_q75", "gray_q95"), quantiles):
                if isinstance(value, (int, float)):
                    state["feature_values"][label].append(float(value))
        panels = descriptor.get("panel_candidates")
        if isinstance(panels, list):
            state["feature_values"]["panel_candidate_count"].append(float(len(panels)))
        pixel_format = descriptor.get("source_pixel_format")
        if isinstance(pixel_format, str):
            state["pixel_format_counts"][pixel_format] += 1

    registry_rows = list(rows(dataset_registry))
    manifests = {row["leakage_component_id"]: row for row in rows(dataset_manifests)}
    datasets: dict[str, dict[str, Any]] = defaultdict(lambda: {
        "accessions": [], "leakage_components": [], "modality_hints": Counter(), "manifest_ids": [],
    })
    for row in registry_rows:
        source = row["source_family_id"].lower()
        if source not in articles:
            continue
        state = datasets[source]
        state["accessions"].append(f"{row['provider']}:{row['accession']}")
        state["leakage_components"].append(row["leakage_component_id"])
        state["modality_hints"].update(row.get("modality_hints", []))
        manifest = manifests.get(row["leakage_component_id"])
        if manifest:
            state["manifest_ids"].append(manifest["manifest_id"])

    union = UnionFind(sources)
    by_payload: dict[str, list[str]] = defaultdict(list)
    by_accession_component: dict[str, list[str]] = defaultdict(list)
    for source, article in articles.items():
        by_payload[article["payload_sha256"]].append(source)
    for source, state in datasets.items():
        for component in state["leakage_components"]:
            by_accession_component[component].append(source)
    for groups in (by_payload, by_accession_component):
        for members in groups.values():
            for member in members[1:]:
                union.union(members[0], member)
    components: dict[str, list[str]] = defaultdict(list)
    for source in sources:
        components[union.find(source)].append(source)
    group_for = {}
    for members in components.values():
        group_id = f"lg:{digest(canonical(sorted(members)))[:16]}"
        for source in members:
            group_for[source] = group_id

    projections = []
    for source in sources:
        article = articles[source]
        nview = numeric.get(source)
        vview = visual.get(source)
        descriptor_view = descriptors.get(source)
        dview = datasets.get(source)
        group_id = group_for[source]
        asset_view = None
        if vview:
            vview = dict(vview)
            acquired_assets = vview.pop("available_asset_count")
            asset_view = {"present": acquired_assets > 0, "acquired_asset_count": acquired_assets}
            vview["modality_tags"] = dict(sorted(vview["modality_tags"].items()))
            vview["candidate_set_sha256"] = digest(canonical(sorted(vview.pop("visual_candidate_ids"))))
            vview["metadata_state"] = "structured_caption_index_and_asset_receipts_only"
        if descriptor_view:
            descriptor_view = dict(descriptor_view)
            descriptor_view["descriptor_refs"] = sorted(descriptor_view["descriptor_refs"], key=lambda row: row["descriptor_id"])
            descriptor_view["descriptor_set_sha256"] = digest(canonical(descriptor_view["descriptor_refs"]))
            aggregates = {}
            for feature, values in sorted(descriptor_view.pop("feature_values").items()):
                mean = sum(values) / len(values)
                variance = sum((value - mean) ** 2 for value in values) / len(values)
                aggregates[feature] = {
                    "count": len(values), "mean": round(mean, 8), "std": round(math.sqrt(variance), 8),
                    "min": round(min(values), 8), "max": round(max(values), 8),
                }
            descriptor_view["numeric_aggregates"] = aggregates
            descriptor_view["pixel_format_counts"] = dict(sorted(descriptor_view["pixel_format_counts"].items()))
            descriptor_view["descriptor_state"] = "source_level_numeric_descriptor_aggregates_available"
        if dview:
            dview = dict(dview)
            dview["accessions"] = sorted(set(dview["accessions"]))
            dview["leakage_components"] = sorted(set(dview["leakage_components"]))
            dview["modality_hints"] = dict(sorted(dview["modality_hints"].items()))
            dview["manifest_ids"] = sorted(set(dview["manifest_ids"]))
        projection = {
            "schema": SCHEMA,
            "authority_status": "proposed",
            "projection_id": f"aleph:learning-source:{digest(source.encode())[:16]}",
            "source_family_id": source,
            "source_payload_sha256": article["payload_sha256"],
            "source_record_set_sha256": article["record_set_sha256"],
            "leakage_group_id": group_id,
            "split": split_for_group(group_id),
            "evaluation_cohort": "frozen_300_gold_candidate_not_ground_truth" if source in gold else "corpus",
            "views": {
                "text_metadata": {
                    "present": True, "domain": article["domain"], "assays": article["assays"],
                    "cell_types": article["cell_types"], "acquisition_pass": article["acquisition_pass"],
                },
                "numeric": {"present": bool(nview), **(nview or {})},
                "visual_metadata": {"present": bool(vview), **(vview or {})},
                "visual_asset_receipt": asset_view or {"present": False, "acquired_asset_count": 0},
                "visual_descriptor": {"present": bool(descriptor_view), **(descriptor_view or {})},
                "dataset_manifest": {"present": bool(dview), **(dview or {})},
            },
            "missing_modality_mask": {
                "numeric": not bool(nview), "visual_metadata": not bool(vview),
                "visual_asset_receipt": not bool(asset_view and asset_view["present"]),
                "visual_descriptor": not bool(descriptor_view),
                "dataset_manifest": not bool(dview),
            },
            "supervision_policy": {
                "ambiguous_observations_are_labels": False,
                "missing_views_are_imputed": False,
                "gold_candidate_is_ground_truth": False,
                "authority_ceiling": "proposed",
            },
        }
        projection["projection_sha256"] = digest(canonical(projection))
        projections.append(projection)
    data = b"".join(canonical(row) + b"\n" for row in projections)
    atomic(output, data)
    splits = Counter(row["split"] for row in projections)
    cohort = sum(row["evaluation_cohort"].startswith("frozen_300") for row in projections)
    group_splits: dict[str, set[str]] = defaultdict(set)
    for row in projections:
        group_splits[row["leakage_group_id"]].add(row["split"])
    input_paths = {
        "article_index": article_index, "gold_candidates": gold_candidates, "visual_index": visual_index,
        "visual_links": visual_links, "asset_receipts": asset_receipts, "visual_descriptors": visual_descriptors,
        "dataset_registry": dataset_registry,
        "dataset_manifests": dataset_manifests,
    }
    summary = {
        "schema": "aleph.external_training.learning_projection_summary.v1",
        "authority_status": "proposed", "version": VERSION,
        "sources": len(projections), "leakage_groups": len(components),
        "leakage_group_stats": {
            "multi_source_groups": sum(len(members) > 1 for members in components.values()),
            "largest_group_sources": max(map(len, components.values()), default=0),
            "duplicate_payload_sha_groups": sum(len(members) > 1 for members in by_payload.values()),
            "multi_source_accession_components": sum(len(set(members)) > 1 for members in by_accession_component.values()),
        },
        "leakage_split_violations": sum(len(value) != 1 for value in group_splits.values()),
        "split_source_counts": dict(sorted(splits.items())),
        "frozen_300_evaluation_cohort_sources": cohort,
        "view_coverage": {
            name: sum(row["views"][name]["present"] for row in projections)
            for name in ("text_metadata", "numeric", "visual_metadata", "visual_asset_receipt", "visual_descriptor", "dataset_manifest")
        },
        "view_totals": {
            "numeric_observations": sum(value["observation_count"] for value in numeric.values()),
            "numeric_candidate_targets_unmasked": sum(value["candidate_target_count"] for value in numeric.values()),
            "numeric_ambiguous_or_context_masked": sum(value["ambiguous_target_mask_count"] for value in numeric.values()),
            "visual_figures": sum(value["figure_count"] for value in visual.values()),
            "visual_linked_figures": sum(value["linked_figure_count"] for value in visual.values()),
            "visual_available_assets": sum(value["available_asset_count"] for value in visual.values()),
            "visual_descriptors": sum(value["descriptor_count"] for value in descriptors.values()),
            "dataset_accession_bindings": sum(len(value["accessions"]) for value in datasets.values()),
        },
        "ambiguous_observations_promoted_to_labels": 0,
        "projection_sha256": digest(data),
        "input_sha256": {name: file_digest(path) for name, path in input_paths.items()},
        "missing_optional_inputs": sorted(name for name, path in input_paths.items() if name != "article_index" and (path is None or not path.is_file())),
        "raw_text_or_pixels_copied": False,
    }
    atomic(summary_path, json.dumps(summary, indent=2, sort_keys=True).encode() + b"\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--article-index", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--gold-candidates", type=Path)
    parser.add_argument("--records-root", type=Path)
    parser.add_argument("--local-numeric-output", type=Path)
    parser.add_argument("--visual-index", type=Path)
    parser.add_argument("--visual-links", type=Path)
    parser.add_argument("--asset-receipts", type=Path)
    parser.add_argument("--visual-descriptors", type=Path)
    parser.add_argument("--dataset-registry", type=Path)
    parser.add_argument("--dataset-manifests", type=Path)
    args = parser.parse_args()
    values = vars(args)
    values["summary_path"] = values.pop("summary")
    result = build(**values)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["leakage_split_violations"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

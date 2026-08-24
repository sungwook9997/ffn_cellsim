#!/usr/bin/env python
"""Frozen article-content inference for the combined external OA router."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("combined_content_train", HERE / "train.py")
train = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(train)


class SourceRefusal(ValueError):
    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


class CombinedContentIndex:
    authority = "external_non_authoritative_combined_oa_router"

    def __init__(self, model_dir: Path = HERE):
        self.model_dir = Path(model_dir).resolve()
        self.config = json.loads((self.model_dir / "config.json").read_text())
        receipts = json.loads((self.model_dir / "SHA256SUMS.json").read_text())
        for name in ("weights_combined_content.npz", "splits.jsonl", "taxonomy.json"):
            if hashlib.sha256((self.model_dir / name).read_bytes()).hexdigest() != receipts[name]:
                raise RuntimeError(f"artifact_digest_mismatch:{name}")
        repo_root = self.model_dir.parents[3]
        taxonomy = json.loads((self.model_dir / "taxonomy.json").read_text())
        metadata, _ = train.load_metadata(repo_root, self.config, taxonomy)
        quality, _ = train.load_quality_families(repo_root, self.config)
        self.records, _ = train.load_combined_records(repo_root, self.config, metadata, quality)
        split_rows = [json.loads(line) for line in (self.model_dir / "splits.jsonl").read_text().splitlines()]
        split_by_family = {row["source_family_id"]: row for row in split_rows}
        if len(split_by_family) != len(split_rows) or set(split_by_family) != {record["source_family_id"] for record in self.records}:
            raise RuntimeError("catalog_split_source_family_mismatch")
        for record in self.records:
            receipt = split_by_family[record["source_family_id"]]
            record.update({"split": receipt["split"], "group_key": receipt["group_key"]})
        archive = np.load(self.model_dir / "weights_combined_content.npz", allow_pickle=False)
        self.weights = [archive[name] for name in ("w1", "b1", "w2", "b2")]
        self.labels = [str(label) for label in archive["labels"]]
        self.features = [train.v2.feature_ids(record, self.config, True) for record in self.records]
        indices = np.arange(len(self.records), dtype=np.int64)
        _, embeddings = train.v2.predict(self.features, indices, self.weights, self.config["model"]["batch_size"])
        self.embeddings = embeddings / np.maximum(np.linalg.norm(embeddings, axis=1, keepdims=True), 1e-12)
        self.by_family = {record["source_family_id"]: i for i, record in enumerate(self.records)}

    def query_source(self, source_family_id: str, top_k: int = 5) -> dict:
        family = train.v1.canonical_source_family(source_family_id)
        if family not in self.by_family:
            raise SourceRefusal("source_family_unavailable", "No technically valid bounded OA content exists for this source family.")
        if top_k < 1 or top_k > 100:
            raise SourceRefusal("invalid_top_k", "top_k must be between 1 and 100.")
        index = self.by_family[family]
        x = train.v2.dense_batch(self.features, np.asarray([index]), self.weights[0].shape[0])
        hidden = np.maximum(x @ self.weights[0] + self.weights[1], 0)
        logits = hidden @ self.weights[2] + self.weights[3]
        logits -= logits.max()
        probabilities = np.exp(logits[0])
        probabilities /= probabilities.sum()
        domain_order = sorted(range(len(self.labels)), key=lambda i: (-float(probabilities[i]), self.labels[i]))
        vector = hidden[0] / max(float(np.linalg.norm(hidden[0])), 1e-12)
        scores = self.embeddings @ vector
        order = sorted((i for i in range(len(self.records)) if i != index), key=lambda i: (-float(scores[i]), self.records[i]["source_family_id"]))[:top_k]
        source = self.records[index]
        return {
            "schema_version": "1.0.0",
            "authority": self.authority,
            "source": {key: source[key] for key in ("source_family_id", "title", "pass", "weak_label", "label", "split", "group_key")},
            "bounded_section_token_counts": {section: len(values) for section, values in source["section_tokens"].items()},
            "coarse_route_probabilities": [{"route": self.labels[i], "probability": float(probabilities[i])} for i in domain_order[:5]],
            "retrieval": [{"rank": rank, "source_family_id": self.records[i]["source_family_id"], "title": self.records[i]["title"], "pass": self.records[i]["pass"], "coarse_route": self.records[i]["label"], "split": self.records[i]["split"], "cosine_score": float(scores[i])} for rank, i in enumerate(order, 1)],
            "refusal_boundary": "External OA routing only; no evidence promotion, physical-parameter estimate, physics validation, or authority claim.",
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-family", required=True)
    parser.add_argument("--top-k", type=int, default=5)
    args = parser.parse_args()
    index = CombinedContentIndex()
    try:
        result = index.query_source(args.source_family, args.top_k)
    except SourceRefusal as exc:
        print(json.dumps({"refused": True, "code": exc.code, "message": str(exc)}, sort_keys=True))
        raise SystemExit(2)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

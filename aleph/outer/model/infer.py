#!/usr/bin/env python
"""Deterministic inference for the external weak-supervision metadata model.

Results are catalogue-routing suggestions. They are not Aleph evidence, physics, or authority.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
from pathlib import Path

import numpy as np


MODEL_DIR = Path(__file__).resolve().parent
TRAIN_SPEC = importlib.util.spec_from_file_location("external_metadata_train", MODEL_DIR / "train.py")
train = importlib.util.module_from_spec(TRAIN_SPEC)
TRAIN_SPEC.loader.exec_module(train)


class QueryRefusal(ValueError):
    """Typed refusal for unsupported or authority-seeking queries."""

    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


PROHIBITED_CLAIM_PATTERNS = (
    re.compile(r"\b(?:fit|infer|estimate|train|predict)\b.{0,32}\bphysical parameters?\b", re.I),
    re.compile(r"\b(?:parameter value|ground truth|authoritative|authority claim)\b", re.I),
    re.compile(r"\b(?:validate|prove)\b.{0,24}\b(?:physics|mechanism|physical law)\b", re.I),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ExternalMetadataIndex:
    """Loads frozen weights and reconstructs the deduplicated metadata catalogue."""

    authority = "external_non_authoritative_weak_supervision"

    def __init__(self, model_dir: Path = MODEL_DIR):
        self.model_dir = Path(model_dir).resolve()
        self.config = json.loads((self.model_dir / "config.json").read_text())
        receipts = json.loads((self.model_dir / "SHA256SUMS.json").read_text())
        for name in ("weights.npz", "splits.jsonl"):
            if _sha256(self.model_dir / name) != receipts[name]:
                raise RuntimeError(f"artifact_digest_mismatch:{name}")

        corpus_root = self.model_dir.parent
        self.records, summary = train.load_records(corpus_root)
        split_rows = [json.loads(line) for line in (self.model_dir / "splits.jsonl").read_text().splitlines()]
        split_by_family = {row["source_family_id"]: row for row in split_rows}
        families = [record["source_family_id"] for record in self.records]
        if len(families) != len(set(families)):
            raise RuntimeError("catalog_source_family_not_unique")
        if set(families) != set(split_by_family):
            raise RuntimeError("catalog_split_family_mismatch")
        for record in self.records:
            receipt = split_by_family[record["source_family_id"]]
            if receipt["label"] != record["label"]:
                raise RuntimeError(f"catalog_label_mismatch:{record['source_family_id']}")
            record["split"] = receipt["split"]

        weights = np.load(self.model_dir / "weights.npz", allow_pickle=False)
        self.w1 = weights["w1"]
        self.b1 = weights["b1"]
        self.w2 = weights["w2"]
        self.b2 = weights["b2"]
        self.labels = [str(value) for value in weights["labels"]]
        self.features = [train.feature_ids(record, self.config["features"]) for record in self.records]
        all_indices = np.arange(len(self.records), dtype=np.int64)
        _, embeddings = train.predict(
            self.features, all_indices, self.w1, self.b1, self.w2, self.b2,
            self.config["model"]["batch_size"],
        )
        self.embeddings = embeddings / np.maximum(np.linalg.norm(embeddings, axis=1, keepdims=True), 1e-12)
        self.known_title_tokens = {
            token for record in self.records
            for token in train.token_strings(record, self.config["features"])
            if token.startswith("title:")
        }
        self.summary = summary

    def _validate_query(self, query: str) -> list[str]:
        if not isinstance(query, str) or not query.strip():
            raise QueryRefusal("empty_query", "A non-empty literature-routing query is required.")
        for pattern in PROHIBITED_CLAIM_PATTERNS:
            if pattern.search(query):
                raise QueryRefusal(
                    "physical_or_authority_claim_refused",
                    "This metadata model may route literature only; it cannot estimate physical parameters or make authority claims.",
                )
        probe = {"title": query, "journal": "", "year": 0}
        tokens = train.token_strings(probe, self.config["features"])
        known = [token for token in tokens if token in self.known_title_tokens]
        if not known:
            raise QueryRefusal(
                "unknown_query",
                "No query token occurs in the frozen metadata catalogue; refusing an ungrounded hash-only result.",
            )
        return known

    def query(self, query: str, top_k: int = 5, domain_top_k: int = 5) -> dict:
        known = self._validate_query(query)
        if top_k < 1 or top_k > 100:
            raise QueryRefusal("invalid_top_k", "top_k must be between 1 and 100.")
        probe = {"title": query, "journal": "", "year": 0}
        ids = train.feature_ids(probe, self.config["features"])
        x = np.zeros((1, self.w1.shape[0]), dtype=np.float32)
        x[0, ids] = np.float32(1.0 / np.sqrt(len(ids)))
        hidden = np.maximum(x @ self.w1 + self.b1, 0)
        logits = hidden @ self.w2 + self.b2
        logits -= logits.max()
        probabilities = np.exp(logits[0])
        probabilities /= probabilities.sum()
        domain_order = sorted(range(len(self.labels)), key=lambda i: (-float(probabilities[i]), self.labels[i]))

        vector = hidden[0]
        vector /= max(float(np.linalg.norm(vector)), 1e-12)
        scores = self.embeddings @ vector
        order = sorted(range(len(self.records)), key=lambda i: (-float(scores[i]), self.records[i]["source_family_id"]))[:top_k]
        results = []
        for rank, index in enumerate(order, 1):
            record = self.records[index]
            results.append({
                "rank": rank,
                "source_family_id": record["source_family_id"],
                "title": record["title"],
                "journal": record["journal"],
                "year": record["year"],
                "weak_domain": record["label"],
                "split": record["split"],
                "cosine_score": float(scores[index]),
            })
        return {
            "schema_version": "1.0.0",
            "authority": self.authority,
            "query": query,
            "known_catalog_tokens": [token.removeprefix("title:") for token in known],
            "weak_domain_probabilities": [
                {"domain": self.labels[i], "probability": float(probabilities[i])}
                for i in domain_order[:domain_top_k]
            ],
            "retrieval": results,
            "refusal_boundary": "Literature routing only: no physical-parameter estimate, physics validation, evidence promotion, or authority claim.",
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query")
    parser.add_argument("--suite", type=Path, help="JSON array of {id, query} objects")
    parser.add_argument("--top-k", type=int, default=5)
    args = parser.parse_args()
    if (args.query is not None) == (args.suite is not None):
        parser.error("provide exactly one of --query or --suite")
    index = ExternalMetadataIndex()
    try:
        if args.query is not None:
            output = index.query(args.query, top_k=args.top_k)
        else:
            suite = json.loads(args.suite.read_text())
            output = [{"id": item["id"], "result": index.query(item["query"], top_k=args.top_k)} for item in suite]
    except QueryRefusal as exc:
        print(json.dumps({"refused": True, "code": exc.code, "message": str(exc)}, sort_keys=True))
        raise SystemExit(2)
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

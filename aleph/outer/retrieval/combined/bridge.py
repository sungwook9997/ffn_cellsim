#!/usr/bin/env python
"""Two-pass neural-to-exact-locator bridge over proposed external OA literature."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[3]

BASE_SPEC = importlib.util.spec_from_file_location("external_retrieval_v1_boundary", HERE.parent / "bridge.py")
base = importlib.util.module_from_spec(BASE_SPEC)
BASE_SPEC.loader.exec_module(base)
MODEL_SPEC = importlib.util.spec_from_file_location(
    "external_combined_content_infer_for_retrieval", HERE.parents[1] / "model" / "combined_content" / "infer.py"
)
combined_infer = importlib.util.module_from_spec(MODEL_SPEC)
MODEL_SPEC.loader.exec_module(combined_infer)

SECTION_BONUS = {"abstract": 0.05, "results": 0.04, "methods": 0.03, "conclusion": 0.03, "figure_caption": 0.01, "table_caption": 0.01}
RetrievalRefusal = base.RetrievalRefusal


class CombinedExternalCorpusRetriever:
    authority = "external_non_authoritative_combined_proposed_literature_routing"

    def __init__(self, retrieval_dir: Path = HERE):
        self.retrieval_dir = Path(retrieval_dir).resolve()
        self.config = json.loads((self.retrieval_dir / "config.json").read_text())
        self.receipt = json.loads((self.retrieval_dir / "receipt.json").read_text())
        index_path = self.retrieval_dir / "index.jsonl"
        if base.sha256(index_path) != self.receipt["index_sha256"]:
            raise RuntimeError("combined_retrieval_index_digest_mismatch")
        self.all_rows = [json.loads(line) for line in index_path.read_text().splitlines() if line.strip()]
        if len(self.all_rows) != self.receipt["local_source_families"]:
            raise RuntimeError("combined_retrieval_index_count_mismatch")
        self.rows = [row for row in self.all_rows if row["neural_routable"]]
        if len(self.rows) != self.receipt["model"]["neural_routable_source_families"]:
            raise RuntimeError("combined_retrieval_routable_count_mismatch")
        self.derived_roots = {name: REPO_ROOT / relative for name, relative in self.config["derived_roots"].items()}
        self.model = combined_infer.CombinedContentIndex(HERE.parents[1] / "model" / "combined_content")
        self.model_indices = np.asarray([self.model.by_family[row["source_family_id"]] for row in self.rows], dtype=np.int64)
        self.local_embeddings = self.model.embeddings[self.model_indices]
        self.known_tokens = {
            token
            for record in self.model.records
            for token in base.TOKEN_RE.findall((record["title"] + " " + " ".join(sum(record["section_tokens"].values(), []))).lower())
        }

    def _validate_query(self, query: str, source_k: int) -> list[str]:
        if not isinstance(query, str) or not query.strip():
            raise RetrievalRefusal("empty_query", "A non-empty proposed-literature query is required.")
        if base.TRAINING_RE.search(query):
            raise RetrievalRefusal("training_refused", "This bridge retrieves proposed literature; it cannot train or fine-tune Aleph models.")
        if base.physical_parameter_estimation_intent(query):
            raise RetrievalRefusal("physical_parameter_estimation_refused", "This bridge cannot estimate or fit physical parameter values.")
        if base.AUTHORITY_RE.search(query):
            raise RetrievalRefusal("authority_claim_refused", "External literature routing cannot validate physics or promote an authority claim.")
        if source_k < 1 or source_k > self.config["maximum_source_k"]:
            raise RetrievalRefusal("invalid_source_k", f"source_k must be between 1 and {self.config['maximum_source_k']}.")
        tokens = sorted(set(base.TOKEN_RE.findall(query.lower())))
        if not set(tokens) & self.known_tokens:
            raise RetrievalRefusal("unknown_query", "No query token occurs in the bounded combined OA model corpus.")
        return tokens

    def _query_embedding(self, query: str) -> tuple[np.ndarray, list[dict]]:
        empty_sections = {name: [] for name in self.model.config["features"]["section_token_limits"]}
        probe = {"title": query, "journal": "", "year": 0, "section_tokens": empty_sections}
        feature_ids = combined_infer.train.v2.feature_ids(probe, self.model.config, True)
        x = np.zeros((1, self.model.weights[0].shape[0]), dtype=np.float32)
        x[0, feature_ids] = np.float32(1.0 / math.sqrt(max(1, len(feature_ids))))
        hidden = np.maximum(x @ self.model.weights[0] + self.model.weights[1], 0)
        logits = hidden @ self.model.weights[2] + self.model.weights[3]
        logits -= logits.max()
        probabilities = np.exp(logits[0])
        probabilities /= probabilities.sum()
        order = sorted(range(len(self.model.labels)), key=lambda i: (-float(probabilities[i]), self.model.labels[i]))
        vector = hidden[0] / max(float(np.linalg.norm(hidden[0])), 1e-12)
        routes = [{"proposed_coarse_domain": self.model.labels[i], "probability": float(probabilities[i])} for i in order[:5]]
        return vector, routes

    def _load_chunks(self, row: dict) -> list[dict]:
        path = self.derived_roots[row["pass"]] / row["chunks_file"]
        if not path.is_file():
            return []
        if base.sha256(path) != row["chunks_file_sha256"]:
            raise RuntimeError(f"local_chunk_file_digest_mismatch:{row['source_family_id']}")
        chunks = []
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                chunk = json.loads(line)
                if chunk["authority_status"] != "proposed":
                    raise RuntimeError(f"non_proposed_local_chunk:{row['source_family_id']}")
                if hashlib.sha256(chunk["text"].encode()).hexdigest() != chunk["text_sha256"]:
                    raise RuntimeError(f"local_chunk_text_digest_mismatch:{chunk['chunk_sha256']}")
                basis = dict(chunk)
                basis.pop("text")
                expected = basis.pop("chunk_sha256")
                if hashlib.sha256(base.canonical_bytes(basis)).hexdigest() != expected:
                    raise RuntimeError(f"local_chunk_record_digest_mismatch:{expected}")
                chunks.append(chunk)
        return chunks

    @staticmethod
    def _chunk_relevance(chunk: dict, query_tokens: list[str], query: str) -> float:
        text_tokens = set(base.TOKEN_RE.findall(chunk["text"].lower()))
        overlap = len(text_tokens & set(query_tokens)) / math.sqrt(max(1, len(query_tokens)))
        phrases = {str(tag).replace("_", " ") for values in chunk.get("tags", {}).values() for tag in values}
        lowered = query.lower()
        tag_hits = sum(1 for phrase in phrases if phrase in lowered or all(part in query_tokens for part in phrase.split()))
        exact_phrase = 0.2 if query.strip().lower() in chunk["text"].lower() else 0.0
        return overlap + 0.35 * tag_hits + exact_phrase + SECTION_BONUS.get(chunk["section_class"], 0.0)

    def query(self, query: str, source_k: int | None = None) -> dict:
        source_k = self.config["default_source_k"] if source_k is None else source_k
        query_tokens = self._validate_query(query, source_k)
        vector, routes = self._query_embedding(query)
        neural_scores = self.local_embeddings @ vector
        pool_size = min(len(self.rows), max(source_k, self.config["candidate_pool"]))
        pool = sorted(range(len(self.rows)), key=lambda i: (-float(neural_scores[i]), self.rows[i]["source_family_id"]))[:pool_size]
        candidates = []
        for index in pool:
            row = self.rows[index]
            chunks = self._load_chunks(row)
            if not chunks:
                continue
            ranked = sorted(
                ((self._chunk_relevance(chunk, query_tokens, query), chunk) for chunk in chunks),
                key=lambda item: (-item[0], item[1]["section_locator"], item[1]["chunk_ordinal"], item[1]["chunk_sha256"]),
            )
            chunk_score, best = ranked[0]
            candidates.append((float(neural_scores[index]) + 0.12 * chunk_score, float(neural_scores[index]), chunk_score, row, best))
        selected = sorted(candidates, key=lambda item: (-item[0], item[3]["source_family_id"]))[:source_k]
        sources = []
        for rank, (combined, neural, chunk_score, row, chunk) in enumerate(selected, 1):
            record = self.model.records[self.model.by_family[row["source_family_id"]]]
            sources.append({
                "rank": rank,
                "source_family_id": row["source_family_id"],
                "pass": row["pass"],
                "pmcid": row["pmcid"],
                "title": record["title"],
                "proposed_coarse_domain": row["proposed_coarse_domain"],
                "weak_label": row["weak_label"],
                "split": row["split"],
                "neural_cosine_score": neural,
                "combined_routing_score": combined,
                "exact_locator_chunk": {
                    "article_locator": chunk["article_locator"],
                    "section_class": chunk["section_class"],
                    "section_locator": chunk["section_locator"],
                    "section_title": chunk["section_title"],
                    "chunk_ordinal": chunk["chunk_ordinal"],
                    "chunk_sha256": chunk["chunk_sha256"],
                    "text_sha256": chunk["text_sha256"],
                    "local_relevance_score": chunk_score,
                    "snippet": base.bounded_snippet(chunk["text"], self.config["maximum_snippet_chars"]),
                },
            })
        if not sources:
            raise RetrievalRefusal("local_content_unavailable", "No digest-verified local OA chunk is available for this query.")
        return {
            "schema": "aleph.external_training.combined_retrieval_result.v1",
            "authority": self.authority,
            "authority_status": "proposed",
            "query": query,
            "known_query_tokens": sorted(set(query_tokens) & self.known_tokens),
            "proposed_coarse_domain_probabilities": routes,
            "sources": sources,
            "index_sha256": self.receipt["index_sha256"],
            "ragtagcag_ledger_head_hash": self.receipt["ragtagcag"]["ledger_head_hash"],
            "refusal_boundary": "Proposed literature routing only; no training, physical-parameter estimate, physics validation, evidence promotion, or authority claim.",
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True)
    parser.add_argument("--source-k", type=int)
    args = parser.parse_args()
    try:
        output = CombinedExternalCorpusRetriever().query(args.query, args.source_k)
    except RetrievalRefusal as exc:
        print(json.dumps({"refused": True, "code": exc.code, "message": str(exc)}, sort_keys=True))
        raise SystemExit(2)
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Deterministic neural-to-exact-locator bridge for proposed OA literature."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import re
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
OA_INFER_PATH = HERE.parent / "model" / "oa_content" / "infer.py"
SPEC = importlib.util.spec_from_file_location("external_oa_content_infer_for_retrieval", OA_INFER_PATH)
oa_infer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(oa_infer)

TOKEN_RE = re.compile(r"[a-z0-9]+(?:[-'][a-z0-9]+)?", re.I)
TRAINING_RE = re.compile(r"\b(?:train|fine[- ]?tune|backpropagate|optimi[sz]e)\b.{0,40}\b(?:model|network|weights?|aleph)\b", re.I)
AUTHORITY_RE = re.compile(r"\b(?:authoritative|authority|ground truth|prove|validate)\b.{0,36}\b(?:claim|physics|mechanism|law|truth|evidence)\b", re.I)
ESTIMATION_ACTION_RE = re.compile(
    r"\b(?:estimate|estimating|estimated|infer|inferring|inferred|fit|fitting|fitted|"
    r"calibrate|calibrating|calibrated|predict|predicting|predicted|calculate|calculating|"
    r"calculated|derive|deriving|derived|determine|determining|determined)\b",
    re.I,
)
PHYSICAL_PARAMETER_TARGET_RE = re.compile(
    r"\b(?:"
    r"physical\s+parameters?|parameters?|parameter\s+values?|numerical\s+values?|units?|"
    r"(?:cortical|membrane|surface|interfacial|line)\s+tensions?|"
    r"tensions?\s+(?:values?|parameters?)|"
    r"young(?:'s)?\s+modulus|(?:elastic|shear|bulk|storage|loss)\s+modulus|modulus\s+(?:values?|parameters?)|"
    r"viscosit(?:y|ies)|stiffness(?:es)?|pressures?|"
    r"(?:traction|contractile|drag|friction|spring)\s+forces?|force\s+(?:values?|magnitudes?|parameters?)|"
    r"(?:diffusion|drag|friction)\s+coefficients?|(?:spring|rate)\s+constants?"
    r")\b",
    re.I,
)
SECTION_BONUS = {"abstract": 0.05, "results": 0.04, "methods": 0.03, "conclusion": 0.03, "figure_caption": 0.01, "table_caption": 0.01}


class RetrievalRefusal(ValueError):
    """Typed refusal at the external-literature authority boundary."""

    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_bytes(value: dict) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()


def bounded_snippet(text: str, limit: int) -> str:
    normalized = " ".join(text.split())
    if len(normalized) <= limit:
        return normalized
    candidate = normalized[: limit - 1].rsplit(" ", 1)[0]
    return (candidate or normalized[: limit - 1]) + "…"


def physical_parameter_estimation_intent(query: str) -> bool:
    """Return true only when an estimation action is paired with a physical value target.

    The conjunction keeps literature searches such as ``methods measuring cortical tension``
    routable: measurement/search nouns alone are not an estimation instruction.
    """

    normalized = query.replace("’", "'").replace("‘", "'")
    return bool(ESTIMATION_ACTION_RE.search(normalized) and PHYSICAL_PARAMETER_TARGET_RE.search(normalized))


class ExternalCorpusRetriever:
    authority = "external_non_authoritative_proposed_literature_routing"

    def __init__(self, retrieval_dir: Path = HERE):
        self.retrieval_dir = Path(retrieval_dir).resolve()
        self.config = json.loads((self.retrieval_dir / "config.json").read_text())
        self.receipt = json.loads((self.retrieval_dir / "receipt.json").read_text())
        index_path = self.retrieval_dir / "index.jsonl"
        if sha256(index_path) != self.receipt["index_sha256"]:
            raise RuntimeError("retrieval_index_digest_mismatch")
        self.rows = [json.loads(line) for line in index_path.read_text().splitlines() if line.strip()]
        if len(self.rows) != self.receipt["local_source_families"]:
            raise RuntimeError("retrieval_index_count_mismatch")
        self.by_family = {row["source_family_id"]: row for row in self.rows}
        if len(self.by_family) != len(self.rows):
            raise RuntimeError("retrieval_index_duplicate_source_family")
        self.derived_root = REPO_ROOT / self.config["derived_root"]
        self.model = oa_infer.OAContentIndex(HERE.parent / "model" / "oa_content")
        self.model_indices = np.asarray([self.model.by_family[row["source_family_id"]] for row in self.rows], dtype=np.int64)
        self.local_embeddings = self.model.embeddings[self.model_indices]
        self.known_tokens = {
            token
            for record in self.model.records
            for token in TOKEN_RE.findall((record["title"] + " " + " ".join(sum(record["section_tokens"].values(), []))).lower())
        }

    def _validate_query(self, query: str, source_k: int) -> list[str]:
        if not isinstance(query, str) or not query.strip():
            raise RetrievalRefusal("empty_query", "A non-empty proposed-literature query is required.")
        if TRAINING_RE.search(query):
            raise RetrievalRefusal("training_refused", "This bridge retrieves proposed literature; it cannot train or fine-tune Aleph models.")
        if physical_parameter_estimation_intent(query):
            raise RetrievalRefusal("physical_parameter_estimation_refused", "This bridge cannot estimate or fit physical parameter values.")
        if AUTHORITY_RE.search(query):
            raise RetrievalRefusal("authority_claim_refused", "External literature routing cannot validate physics or promote an authority claim.")
        if source_k < 1 or source_k > self.config["maximum_source_k"]:
            raise RetrievalRefusal("invalid_source_k", f"source_k must be between 1 and {self.config['maximum_source_k']}.")
        tokens = sorted(set(TOKEN_RE.findall(query.lower())))
        if not set(tokens) & self.known_tokens:
            raise RetrievalRefusal("unknown_query", "No query token occurs in the bounded OA model corpus.")
        return tokens

    def _query_embedding(self, query: str) -> tuple[np.ndarray, list[dict]]:
        empty_sections = {name: [] for name in self.model.config["features"]["section_token_limits"]}
        probe = {"title": query, "journal": "", "year": 0, "section_tokens": empty_sections}
        feature_ids = oa_infer.train.feature_ids(probe, self.model.config, True)
        x = np.zeros((1, self.model.weights[0].shape[0]), dtype=np.float32)
        x[0, feature_ids] = np.float32(1.0 / math.sqrt(max(1, len(feature_ids))))
        hidden = np.maximum(x @ self.model.weights[0] + self.model.weights[1], 0)
        logits = hidden @ self.model.weights[2] + self.model.weights[3]
        logits -= logits.max()
        probabilities = np.exp(logits[0])
        probabilities /= probabilities.sum()
        order = sorted(range(len(self.model.labels)), key=lambda i: (-float(probabilities[i]), self.model.labels[i]))
        vector = hidden[0] / max(float(np.linalg.norm(hidden[0])), 1e-12)
        domains = [{"domain": self.model.labels[i], "probability": float(probabilities[i])} for i in order[:5]]
        return vector, domains

    def _load_chunks(self, row: dict) -> list[dict]:
        path = self.derived_root / row["chunks_file"]
        if not path.is_file():
            return []
        if sha256(path) != row["chunks_file_sha256"]:
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
                expected_chunk_hash = basis.pop("chunk_sha256")
                if hashlib.sha256(canonical_bytes(basis)).hexdigest() != expected_chunk_hash:
                    raise RuntimeError(f"local_chunk_record_digest_mismatch:{expected_chunk_hash}")
                chunks.append(chunk)
        return chunks

    @staticmethod
    def _chunk_relevance(chunk: dict, query_tokens: list[str], query: str) -> float:
        text_tokens = set(TOKEN_RE.findall(chunk["text"].lower()))
        overlap = len(text_tokens & set(query_tokens)) / math.sqrt(max(1, len(query_tokens)))
        tag_phrases = {
            str(tag).replace("_", " ")
            for values in chunk.get("tags", {}).values()
            for tag in values
        }
        lowered = query.lower()
        tag_hits = sum(1 for phrase in tag_phrases if phrase in lowered or all(part in query_tokens for part in phrase.split()))
        exact_phrase = 0.2 if query.strip().lower() in chunk["text"].lower() else 0.0
        return overlap + 0.35 * tag_hits + exact_phrase + SECTION_BONUS.get(chunk["section_class"], 0.0)

    def query(self, query: str, source_k: int | None = None) -> dict:
        source_k = self.config["default_source_k"] if source_k is None else source_k
        query_tokens = self._validate_query(query, source_k)
        vector, domains = self._query_embedding(query)
        neural_scores = self.local_embeddings @ vector
        pool_size = min(len(self.rows), max(source_k, self.config["candidate_pool"]))
        pool = sorted(range(len(self.rows)), key=lambda i: (-float(neural_scores[i]), self.rows[i]["source_family_id"]))[:pool_size]
        candidates = []
        for index in pool:
            index_row = self.rows[index]
            chunks = self._load_chunks(index_row)
            if not chunks:
                continue
            ranked_chunks = sorted(
                ((self._chunk_relevance(chunk, query_tokens, query), chunk) for chunk in chunks),
                key=lambda item: (-item[0], item[1]["section_locator"], item[1]["chunk_ordinal"], item[1]["chunk_sha256"]),
            )
            chunk_score, best = ranked_chunks[0]
            combined = float(neural_scores[index]) + 0.12 * chunk_score
            candidates.append((combined, float(neural_scores[index]), chunk_score, index_row, best))
        selected = sorted(candidates, key=lambda item: (-item[0], item[3]["source_family_id"]))[:source_k]
        sources = []
        for rank, (combined, neural, chunk_score, index_row, chunk) in enumerate(selected, 1):
            model_record = self.model.records[self.model.by_family[index_row["source_family_id"]]]
            sources.append({
                "rank": rank,
                "source_family_id": index_row["source_family_id"],
                "pmcid": index_row["pmcid"],
                "title": model_record["title"],
                "weak_domain": index_row["weak_domain"],
                "split": index_row["split"],
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
                    "snippet": bounded_snippet(chunk["text"], self.config["maximum_snippet_chars"]),
                },
            })
        if not sources:
            raise RetrievalRefusal("local_content_unavailable", "No digest-verified local OA chunk is available for this query.")
        return {
            "schema": "aleph.external_training.retrieval_result.v1",
            "authority_status": "proposed",
            "authority": self.authority,
            "query": query,
            "known_query_tokens": sorted(set(query_tokens) & self.known_tokens),
            "weak_domain_probabilities": domains,
            "sources": sources,
            "index_sha256": self.receipt["index_sha256"],
            "ragtagcag_ledger_head_hash": self.receipt["ragtagcag"]["ledger_head_hash"],
            "refusal_boundary": "Proposed literature routing only; no model training, physical-parameter estimate, physics validation, evidence promotion, or authority claim.",
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True)
    parser.add_argument("--source-k", type=int)
    args = parser.parse_args()
    try:
        result = ExternalCorpusRetriever().query(args.query, args.source_k)
    except RetrievalRefusal as exc:
        print(json.dumps({"refused": True, "code": exc.code, "message": str(exc)}, sort_keys=True))
        raise SystemExit(2)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

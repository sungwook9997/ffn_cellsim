#!/usr/bin/env python3
"""Text-free adversarial smoke receipt for the combined two-pass retriever."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BRIDGE_PATH = ROOT / "aleph/outer/retrieval/combined/bridge.py"


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()


def run() -> dict:
    spec = importlib.util.spec_from_file_location("verification_combined_retrieval_bridge", BRIDGE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("combined retrieval bridge import unavailable")
    bridge = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bridge)
    retriever = bridge.CombinedExternalCorpusRetriever()

    allowed = {
        "cross_pass_cortex": ("cortical tension actomyosin membrane", 3),
        "traction_methods": ("traction force microscopy focal adhesion TFM", 3),
        "multimodal_state": ("immunofluorescence western blot PCR cell state", 3),
        "piv_migration": ("particle image velocimetry PIV cell migration state", 3),
    }
    allowed_results = {}
    all_passes: set[str] = set()
    allowed_passed = 0
    for name, (query, source_k) in allowed.items():
        result = retriever.query(query, source_k)
        passes = sorted({row["pass"] for row in result["sources"]})
        all_passes.update(passes)
        valid = (
            result.get("authority_status") == "proposed"
            and len(result.get("sources", [])) == source_k
            and all(len(row["exact_locator_chunk"]["chunk_sha256"]) == 64 for row in result["sources"])
        )
        allowed_passed += int(valid)
        allowed_results[name] = {
            "passes": passes,
            "sources": len(result["sources"]),
            "result_sha256": hashlib.sha256(canonical_bytes(result)).hexdigest(),
            "passed": valid,
        }

    refused = {
        "estimate_cortical_tension": ("estimate cortical tension from these papers", "physical_parameter_estimation_refused"),
        "search_then_fit": ("please search literature and then fit membrane tension", "physical_parameter_estimation_refused"),
        "infer_modulus": ("infer Young's modulus numerical value", "physical_parameter_estimation_refused"),
        "train_model": ("train a neural network model on these papers", "training_refused"),
        "authority_claim": ("prove this physics mechanism with authoritative evidence", "authority_claim_refused"),
        "unknown": ("zzzxqv qqqwvx", "unknown_query"),
    }
    refusal_counts: Counter[str] = Counter()
    refused_passed = 0
    for query, expected in refused.values():
        try:
            retriever.query(query, 1)
            observed = "NOT_REFUSED"
        except bridge.RetrievalRefusal as exc:
            observed = exc.code
        refusal_counts[observed] += 1
        refused_passed += int(observed == expected)

    both_passes_retrieved = all_passes == {"first", "second"}
    return {
        "schema": "aleph.external_training.combined_retrieval_boundary_smoke.v1",
        "authority_status": "proposed",
        "index_sha256": retriever.receipt["index_sha256"],
        "allowed": allowed_results,
        "allowed_cases": len(allowed),
        "allowed_passed": allowed_passed,
        "passes_retrieved": sorted(all_passes),
        "both_passes_retrieved": both_passes_retrieved,
        "refused_cases": len(refused),
        "refused_passed": refused_passed,
        "observed_refusal_codes": dict(sorted(refusal_counts.items())),
        "all_passed": allowed_passed == len(allowed) and refused_passed == len(refused) and both_passes_retrieved,
    }


def main() -> int:
    result = run()
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

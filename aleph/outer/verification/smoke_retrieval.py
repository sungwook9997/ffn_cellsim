#!/usr/bin/env python3
"""Adversarial allow/refuse smoke test for the proposed literature retriever."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BRIDGE_PATH = ROOT / "aleph/outer/retrieval/bridge.py"


def run() -> dict:
    spec = importlib.util.spec_from_file_location("verification_retrieval_bridge", BRIDGE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("retrieval bridge import unavailable")
    bridge = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bridge)
    retriever = bridge.ExternalCorpusRetriever()
    refused = {
        "estimate cortical tension from these papers": "physical_parameter_estimation_refused",
        "please search literature and then fit membrane tension": "physical_parameter_estimation_refused",
        "infer Young's modulus numerical value": "physical_parameter_estimation_refused",
        "calibrate viscosity using the retrieved studies": "physical_parameter_estimation_refused",
        "train a neural network model on these papers": "training_refused",
        "prove this physics mechanism with authoritative evidence": "authority_claim_refused",
        "zzzxqv qqqwvx": "unknown_query",
    }
    allowed = (
        "methods measuring cortical tension",
        "literature on reported cortical tension values",
        "studies of Young's modulus measurement",
        "viscosity measurement methods cell",
        "traction force microscopy methods",
    )
    refusal_counts: dict[str, int] = {}
    refused_passed = 0
    for query, expected in refused.items():
        try:
            retriever.query(query, 1)
            observed = "NOT_REFUSED"
        except bridge.RetrievalRefusal as exc:
            observed = exc.code
        if observed == expected:
            refused_passed += 1
        refusal_counts[observed] = refusal_counts.get(observed, 0) + 1
    allowed_passed = 0
    all_allowed_proposed = True
    for query in allowed:
        try:
            result = retriever.query(query, 1)
            passed = result.get("authority_status") == "proposed" and len(result.get("sources", [])) == 1
            all_allowed_proposed &= result.get("authority_status") == "proposed"
        except bridge.RetrievalRefusal:
            passed = False
            all_allowed_proposed = False
        allowed_passed += int(passed)
    return {
        "schema": "aleph.external_training.retrieval_boundary_smoke.v1",
        "authority_status": "proposed",
        "retrieval_index_scope": "first_pass_only",
        "allowed_cases": len(allowed),
        "allowed_passed": allowed_passed,
        "refused_cases": len(refused),
        "refused_passed": refused_passed,
        "observed_refusal_codes": dict(sorted(refusal_counts.items())),
        "all_allowed_results_proposed": all_allowed_proposed,
        "all_passed": allowed_passed == len(allowed) and refused_passed == len(refused),
    }


def main() -> int:
    result = run()
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

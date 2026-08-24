# Combined two-pass external OA retrieval

This bridge is the deterministic, CPU-only query surface for the completed first and second OA
passes. It binds three independently receipted layers:

1. the combined RAG–TAG–CAG projection (10,208 proposed source families);
2. the frozen v3 combined-content router and explicit 12-class coarse taxonomy;
3. 5,675 local OA JATS families containing 284,620 exact-locator chunks.

The neural model can route 5,669 of those local families. Six first-pass files lacked the bounded
target sections required by model v3; they remain explicitly counted in the text-free local index
instead of being hidden or relabelled. Query results are drawn only from the 5,669 digest-bound
model/local intersection.

Run from the repository root:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/retrieval/combined/bridge.py \
  --query "cortical tension actomyosin membrane" --source-k 5
```

Each result returns the source family, acquisition pass, PMCID, proposed coarse domain, weak source
label, split, neural/reranking scores, and exact article/section/chunk locators and SHA-256 hashes.
A snippet of at most 280 characters is returned only after the ignored local chunk file, its text,
and its canonical chunk record all pass digest checks.

The committed index contains no article text. Every output remains `proposed`. Empty and unknown
queries, training requests, explicit physical-parameter estimation/fitting/calibration/prediction,
and authority claims return typed refusals inherited from the hardened v1 boundary.

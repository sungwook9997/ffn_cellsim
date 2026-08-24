# External OA retrieval bridge

This directory connects Aleph's **proposed** external-literature pipeline end to end:

1. the frozen OA-content neural model embeds a free-text literature query;
2. cosine routing selects a bounded candidate pool among 2,695 local OA source families;
3. digest-verified JATS chunks rerank the candidates by query overlap and extracted tags;
4. each result returns a source family, PMCID, exact section/chunk locator and hashes, plus at most
   280 characters of local text.

It does not train a model, infer a physical parameter, validate physics, or promote evidence. Those
requests return typed refusals. All results retain `authority_status: proposed`.

Run from the repository root:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/retrieval/bridge.py \
  --query "cortical tension actomyosin membrane" --source-k 5
```

The committed `index.jsonl` contains identifiers, local filenames, counts and SHA-256 digests only.
It contains no article text. Snippets are emitted only when the corresponding ignored local chunk
file exists and passes its committed file, text and chunk-record digests.

Rebuild the compact index with `build_index.py`. The build also binds the receipt to the committed
RAG ledger head/projection receipt and the frozen OA-content model weights.

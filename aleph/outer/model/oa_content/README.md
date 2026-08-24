# External OA-content neural baseline v2

This is a separate, non-authoritative CPU experiment. It uses bounded tokens from OA JATS
`abstract`, `methods`, `results`, and `conclusion` chunks plus title/journal/year metadata to predict
the existing weak discovery-domain label. It does not import Aleph learning or representation code,
train a physical parameter, validate physics, or promote evidence.

The filesystem contains 2,701 chunk files, not 2,701 successful source families: one file is empty,
leaving 2,700 non-empty unique families. Five lack every selected target section, so the fair neural
comparison uses 2,695 families. `metrics.json` records this discrepancy instead of inflating coverage.

Run training:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python corpus/external_training/model/oa_content/train.py \
  --config corpus/external_training/model/oa_content/config.json
```

Run frozen article-content inference:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python corpus/external_training/model/oa_content/infer.py \
  --source-family doi:10.1038/s41556-025-01807-6 --top-k 5
```

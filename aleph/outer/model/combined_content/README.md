# Combined first+second OA content router v3

This external CPU-only experiment combines technically valid bounded JATS content from both OA
passes. It maps incompatible first-pass 29-way and second-pass 12-way discovery labels into the 12
explicit coarse routes in `taxonomy.json`. Those routes organise literature; they are not a
biological ontology, article-quality verdict, physical model, or Aleph authority.

The split unit is normally canonical source family. All families connected by any proposed shared
dataset-accession group are first collapsed transitively and assigned as one unit. The accession
signal is treated conservatively as a leakage constraint, not proof that processed datasets match.

Train:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python corpus/external_training/model/combined_content/train.py \
  --config corpus/external_training/model/combined_content/config.json
```

Infer from one frozen OA article:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python corpus/external_training/model/combined_content/infer.py \
  --source-family doi:10.1038/ncb3525 --top-k 5
```

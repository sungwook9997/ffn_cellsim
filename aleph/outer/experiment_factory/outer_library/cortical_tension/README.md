# Cortical-tension external reference head

This package prepares a protocol-matched external reference for a later engine comparison. It does not tune Aleph, create a validation gate, or treat a literature interval as a parameter range.

The admitted numeric reference is deliberately narrow: Hosseini 2020, MCF-7 suspended/rounded interphase cells at 37 C, dynamic AFM parallel-plate confinement, total effective Laplace surface tension. The structured summary is median 270 pN/um with an empirical IQR of 180--400 pN/um for 27 cells. `benchmark.py` verifies those fields against `KB-6.1.7`, `SE181`, `source_audit.verdict=OK`, and the still-draft `VG-γ-MCF7-suspended` in the read-only KB snapshot before emitting them.

`comparison_template.json` shows the fields a future engine record must add. A comparison stays value-redacted unless every supplied record is full-native CUDA with every compartment present, physically accepted, time-step independent, stationary, quantitatively unblocked, protocol-matched, and produced by a natively validated geometry-matched observation operator. Three distinct seed records are required because the standing measurement rule compares against between-seed scatter, never a within-run SEM. Passing those checks only permits a descriptive observable comparison; it never scores the draft validation gate.

Run the audit and generate the machine-readable CSV:

```bash
PY=/Users/sw1/miniconda3/envs/ffn_sim/bin/python
$PY -m aleph.outer.experiment_factory.outer_library.cortical_tension.benchmark \
  --catalog aleph/outer/experiment_factory/outer_library/cortical_tension/reference_catalog.json \
  --kb aleph/outputs/tag_kb/kb.duckdb \
  --operator-registry aleph/outer/experiment_factory/outer_library/observation_operator_registry.json \
  --engine-record aleph/outputs/ac/engine_cortex_ab/settled_on/record.json \
  --reference-csv aleph/outer/experiment_factory/outer_library/results/cortical_tension_reference_observations.csv \
  --output aleph/outer/experiment_factory/outer_library/results/cortical_tension_engine_readiness.json
```

Add `--engine-record PATH` three or more times for distinct-seed run records. Existing force-accepted records are useful refusal controls: the report lists their missing evidence fields but never republishes their blocked gamma values or ratios.

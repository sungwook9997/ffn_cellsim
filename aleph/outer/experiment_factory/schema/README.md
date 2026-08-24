# Proposed ExperimentRecord schema v1

This directory defines a non-authoritative interchange contract for turning external papers into
experiment-level annotation candidates. It does not promote a paper, an extracted value, or a model
output into Aleph physics evidence.

The record graph is:

`Source → Experiment → BiologicalSystem → Condition → Protocol → Observation → Comparison → Provenance`

`experiment_record.schema.json` is the complete Draft 2020-12 declaration.
`controlled_vocabularies.v1.json` supplies the broad twelve-domain routing vocabulary and common
measurement terms. Terms that cannot be normalized are represented as `unmapped`, `missing`, or
`ambiguous`; they are never silently guessed.

Run the dependency-free semantic validator with the Aleph interpreter:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/experiment_factory/schema/validate.py \
  path/to/records.jsonl
```

The CLI verifies referential integrity, exact evidence digests, explicit unit state, conversion
provenance, replicate semantics, ambiguity representation, and the extra requirements for an
`adjudicated` record. The checked-in controls include one valid record and two deliberately broken
records. JSON Schema conformance remains available to downstream tooling that implements Draft
2020-12; the local CLI intentionally has no added dependency.

IDs are stable content-derived names, not assertions of identity. Use `id_helpers.py` and retain the
local ordinal whenever a paper reports multiple experiments or observations at the same locator.

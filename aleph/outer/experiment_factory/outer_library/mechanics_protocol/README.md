# External mechanics protocol network

This package is the Project Aleph Outer Library head for measurements that are
often collapsed under the word “viscoelasticity” but do not measure the same
object. It keeps protocol, mechanical scope, cell state, probe scale, and unit
as first-class variables. It has no Aleph physics or parameter authority.

## External data

`protocol_catalog.json` contains 19 real cohort records from nine primary
papers. They cover magnetic rotational spectroscopy, magnetic-bead
microrheometry, continuous micropipette aspiration, local AFM relaxation or
indentation, whole-cell AFM confinement, and electrodeformation. Aggregate
cohorts stay aggregate; they are never duplicated into synthetic cells.

The trainable individual-level source is Dessard et al. 2024. The acquisition
builder retrieves the official Europe PMC supplementary archive, verifies the
publisher-declared workbook MD5, and normalizes all 196 wires: MCF-10A 68,
MCF-7 60, and MDA-MB-231 68. The raw workbook remains byte-for-byte unchanged.

```bash
python3 aleph/outer/experiment_factory/outer_library/mechanics_protocol/acquire_dessard.py \
  --output-dir aleph/outputs/outer/mechanics_protocol
python3 aleph/outer/experiment_factory/outer_library/mechanics_protocol/train.py \
  --tensor aleph/outputs/outer/mechanics_protocol/dessard2024_mrs_wire_data.npz \
  --catalog aleph/outer/experiment_factory/outer_library/mechanics_protocol/protocol_catalog.json \
  --output-dir aleph/outputs/outer/mechanics_protocol/model
```

The run is deterministic and CPU-only. A 3090 would add scheduling and transfer
overhead to a 196-row, two-layer NumPy model, so no shared GPU was leased.

An exact catalog route can then be exercised with:

```bash
python3 -m aleph.outer.experiment_factory.outer_library.mechanics_protocol.infer route \
  --checkpoint aleph/outputs/outer/mechanics_protocol/model/mechanics_protocol_router.npz \
  --catalog aleph/outer/experiment_factory/outer_library/mechanics_protocol/protocol_catalog.json \
  --request aleph/outer/experiment_factory/outer_library/mechanics_protocol/examples/mcf7_mrs_request.json
```

Changing the observable or unit to an unseen scope produces a refusal instead
of a nearest-neighbor value transfer.

## What was learned—and what failed

The wire network shares a 16-unit latent layer between a three-class cell-line
head and log-viscosity regression. On the sealed experimental-date holdout, the
effective-viscosity head reaches R² 0.980 and 12.1% median absolute percentage
error. That success is expected to be largely operator reconstruction: wire
geometry and critical frequency are the inputs from which viscosity is derived.
It demonstrates a usable measurement-operator representation, not an
independent viscosity biomarker.

The catalog also retains the most important cross-method discrepancy without
normalizing it away: Vaippully et al. report MCF-7 zero-frequency Brownian-probe
microviscosity of 0.002–0.016 Pa s, while the Dessard all-wire median is
65.9 Pa s. They are routed to `cytoplasmic_solvent_microviscosity` and
`effective_cytoplasm_viscoelasticity`, respectively. The roughly four-order
difference is an observable to explain, not a pair of interchangeable estimates.

The cell-line head fails robust transfer: date-disjoint balanced accuracy is
0.385 and macro-F1 is 0.227. Even a random wire holdout reaches only balanced
accuracy 0.464. The model therefore must not classify a new cell or laboratory.
This is a useful negative result: the within-study viscosity contrast does not
turn 196 technical observations into an externally validated cell-state model.

The protocol router fits method family, mechanical scope, and adherent versus
suspended state from probe geometry, observable, unit, and time-scale metadata.
Its 100% resubstitution score is recorded only as a wiring diagnostic. Several
classes occur in one paper, so the router status is
`routing_representation_only_insufficient_external_holdout`. Runtime use is
exact ontology routing or abstention, never cross-scope interpolation.

## Probe-size result and the C-1 experiment

Direct fits to every raw wire give length exponents 1.65 (MCF-10A), 1.18
(MCF-7), and 2.04 (MDA-MB-231); date-cluster bootstrap intervals are stored in
the report. These raw fits differ from the paper's binned estimates 2.2, 2.2,
and 2.5, which is retained rather than hidden. Both analyses support the
qualitative fact that apparent cytoplasmic viscosity depends on probe scale.

This does **not** equate wire length with AFM contact radius. The decisive engine
comparison remains preregistered as a separate observable: for an AFM contact
radius sweep from 1 to 5 µm, a poroelastic clock predicts 25-fold movement in
relaxation time, while turnover, membrane drainage, and a speed-set viscous
transient predict no radius-dependent movement. The model can label and score
that hypothesis after a geometry-matched curve exists; it cannot manufacture
the curve or choose `D_um2_s`.

## Artifacts

- `aleph/outputs/outer/mechanics_protocol/source_receipt.json`: source identity,
  license, hashes, schema, and row counts.
- `aleph/outputs/outer/mechanics_protocol/model/mechanics_protocol_report.json`:
  split definitions, all metrics, probe exponents, capabilities, and refusals.
- `aleph/outputs/outer/mechanics_protocol/model/run_disjoint_cell_line_confusion.svg`:
  the sealed-date failure, shown rather than summarized away.
- The two `.npz` checkpoints and source workbook are reproducible local outputs
  and intentionally ignored by Git.

## Architecture report v3

`build_architecture_report.py` turns the committed mechanics, cortical-tension,
multitask, and label-free-vision records into a 20-page evidence-bounded PDF.
The opening architecture section now includes the implemented tensor runtime,
the exact layer/dimension inventory for all nine trained encoder heads, the
complete 30-head authority registry, and the optional mechanics/cortical
adapters. It explicitly records `shared_latent=false` and
`fusion=task_scoped_late_fusion_only` so the report does not imply a universal
end-to-end network that has not been built.
It adds an implemented-vs-target architecture map, the cross-method viscosity
scale, wire-network split and confusion matrix, raw-vs-binned probe exponents,
the C-1 radius-sweep discriminator, and the fail-closed cortical-tension engine
adapter. The companion manifest hashes every input and records that no GPU,
physics run, validation-gate score, or Aleph parameter authority was used.
The report also consumes `afm_c2_readiness_snapshot.json`, a hash-pinned Outer
snapshot of the sibling AFM branch at `f7550895`: C-2 accepted full-native
physical steps for three seeds, while the quantitative sweep remains blocked
on physical inputs and active-state/operator evidence.

The bundled Codex document runtime contains ReportLab and pdfplumber:

```bash
/Users/sw1/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  aleph/outer/experiment_factory/outer_library/mechanics_protocol/build_architecture_report.py \
  --delivery-copy output/pdf/external_neural_network_architecture_report_2026-08-10_v3.pdf
```

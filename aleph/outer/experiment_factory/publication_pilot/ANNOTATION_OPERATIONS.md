# Human annotation operations

## Roles and blinding

- Annotator A and annotator B independently read the same source and never see the other submission.
- Neither annotator sees `machine_preannotations.jsonl.gz` during the independent pass.
- The adjudicator sees the paper, both claim sets without annotator identities, and the generated
  disagreement packet. The adjudicator must not accept a machine claim merely because it exists.
- One person cannot occupy both A and B for a source. The evaluator enforces this mechanically.

## Claim unit

One claim is one exact field value attached to an experiment anchor and an evidence key. Allowed fields
cover experiment boundaries, biological system, conditions, protocol, quantitative observations,
replicates, uncertainty, comparisons, panels and table cells. Missing and ambiguous are recorded rather
than inferred.

The evidence key concatenates document digest, locator type, locator ID, character offsets and evidence
digest. Annotators may introduce a corrected key after inspecting the OA source, but may not paste full
paragraph or caption text into the committed submission.

## Operating sequence

1. Calibrate the instructions on five papers without adding those decisions to the frozen benchmark.
2. Assign all 100 A/B slots, ensuring distinct people per source.
3. Freeze completed A/B files and their hashes before computing agreement.
4. Generate disagreement packets. The tool performs no automatic adjudication.
5. A human adjudicator resolves experiment boundaries first, then dependent biological, condition,
   protocol and observation claims.
6. Freeze adjudicated submissions and run exact field-level metrics.
7. Review every false positive and false negative against `FAILURE_TAXONOMY.md`.
8. Amend extraction rules only against the development portion; never tune on the final held-out pilot.

## Publication gates

- 50/50 sources have distinct A and B submissions.
- 50/50 sources have human adjudication and no pending disagreement.
- Submission hashes, annotator-role separation and evidence keys validate.
- Field-level support and confidence intervals are reported; fields with inadequate support are marked
  underpowered rather than pooled into a flattering global score.
- The 50-paper pilot remains external to model fitting and hyperparameter selection.
- Machine candidates, missing values and ambiguous values are reported separately.

Until these gates pass, the manuscript may describe the factory and the frozen validation protocol but
must not report human accuracy, agreement, biological generalization or gold-dataset claims.

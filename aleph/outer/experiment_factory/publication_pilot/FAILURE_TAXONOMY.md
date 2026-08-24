# ExperimentRecord failure taxonomy

Each false positive, false negative and annotator disagreement receives one primary failure class and
optional contributing classes. The taxonomy separates document recovery, experiment semantics and
normalization so that a downstream model is not blamed for an upstream evidence error.

| Code | Failure class | Examples | Required disposition |
|---|---|---|---|
| F01 | source/locator failure | wrong section, caption, table cell or character span | correct locator; invalidate dependent claims |
| F02 | experiment boundary | two experiments merged or one experiment split | adjudicate boundary before all child fields |
| F03 | biological identity | organism, tissue, cell type or cell line wrong | retain reported label and normalization state |
| F04 | cell-state semantics | treatment confused with state; cancer/migration/differentiation over-inferred | mark ambiguous or map to controlled state |
| F05 | condition-arm linkage | control/treatment, dose, time point or substrate assigned to wrong arm | repair condition graph and dependent comparison |
| F06 | protocol/assay | assay keyword mention mistaken for performed protocol | require methods/result evidence |
| F07 | instrument/calibration | model, objective, bead size or calibration absent/wrong | preserve missing; never infer default calibration |
| F08 | numeric tokenization | mutation/sample notation parsed as measurement | reject token; regression-test exact context |
| F09 | quantity/unit semantics | `nM`/`nm`, force/sample `n`, concentration/length ambiguity | case-sensitive normalization or explicit ambiguity |
| F10 | uncertainty semantics | SD, SE, CI, range or error bar type confused | do not normalize unknown error bars |
| F11 | replicate semantics | cells, fields, wells and donors conflated | separate biological/technical n and replication unit |
| F12 | table structure | header, row, footnote or multi-level cell association wrong | bind exact table cell and header path |
| F13 | figure/panel | panel, channel, scale bar or time point assigned incorrectly | retain panel candidate; require human confirmation |
| F14 | comparison direction | before/after or control/treatment reversed | repair condition linkage before direction |
| F15 | derived-vs-reported | computed value presented as directly reported | mark derived and retain formula provenance |
| F16 | duplicate evidence | the same claim emitted from paragraph, table and caption | deduplicate only when evidence and semantics match |
| F17 | ontology normalization | synonym mapped to wrong identifier or excessive specificity | preserve source label and candidates |
| F18 | inaccessible evidence | supplement, asset or dataset unavailable | typed unavailable state; no silent omission |

Known machine controls already cover `18S`, `W748S`, `nM`, `nm`, uppercase force notation, bare sample
`n`, and drug-context micrometre ambiguity. Those controls demonstrate regression protection, not the
prevalence of every class. Prevalence will remain unreported until the adjudicated pilot exists.

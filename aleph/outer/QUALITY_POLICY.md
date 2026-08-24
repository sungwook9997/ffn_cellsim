# Evidence quality and journal-metric policy

## Default screening lane

The default candidate pool is limited to peer-reviewed biomedical records published from 2010 onward. A journal impact factor (JIF) of at least 5 is a **screening preference**, not a claim that an individual article is reliable.

The planned 1,200-source corpus is divided into:

- 800 primary-lane sources: 2010 or later, research article, JIF >= 5 where a licensed JCR lookup is available, and article-level quality gates passed;
- 400 exception/method/foundation sources: indispensable methods, strong open datasets, independent replications, negative results, or pre-2010 foundational measurements.

No source is excluded solely because JIF is missing. JIF is time-varying and journal-level, so every recorded value requires `metric_name`, `metric_value`, `metric_year`, `metric_source`, and `verified_at`. CiteScore, SJR, or an inferred value must never be labelled JIF.

## Article-level quality gates

Tier A requires all applicable fields to be reviewed:

1. biological system, cell state, perturbation, geometry and substrate are explicit;
2. sample size and independent biological replication are recoverable;
3. calibration, units, uncertainty and exclusion rules are recoverable;
4. raw or source data availability is recorded;
5. the observation is linked to an exact figure/table/supplement locator;
6. later replication, contradiction, correction or retraction checks are complete;
7. no source-family, laboratory or dataset leakage crosses evaluation splits.

Journal prestige cannot compensate for failure of these gates.

## Exception reasons

Every exception must name one of: `foundational_pre_2010`, `validated_method`, `open_benchmark`, `independent_replication`, `negative_evidence`, or `rare_modality`. Exceptions remain visible in evaluation reports and are never silently mixed into the primary lane.

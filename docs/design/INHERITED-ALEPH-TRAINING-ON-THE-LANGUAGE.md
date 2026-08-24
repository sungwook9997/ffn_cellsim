# Training expressed in the language, not beside it

| Field | Value |
|---|---|
| Status | `SPEC` — companion to `ALEPH-TN-LANGUAGE-SPEC.md`, implemented by `aleph/learn/training.py` |
| Written | 2026-07-30 12:45 KST |
| Provenance | PI, 2026-07-30 12:31 (transcript uuid `15c0dac5-3270-48d6-8141-61799ab402f1`): "나중에 학습도 이 언어 토대로 해야되니깐 잘 해야지" |
| Authority | `ALEPH-PD-002`, and the 14:24 ratification in `docs/decisions/RATIFICATION-2026-07-30-seal-is-a-computation-boundary.md`. **Training itself stays sealed until R5.** |
| Implementation | `aleph/learn/training.py` declares a training run and refuses to perform one. |

## 1. The failure this document is written against

The normal way a project gets this wrong is not dramatic. The representation layer ships, then someone
writes a training script that loads arrays, hardcodes an axis order because the loader needed one,
picks a rank because the memory budget needed one, and splits the data randomly because that was the
default. Nothing in that script is a lie. But afterwards the trained model's provenance is a shell
history, its axis order is not the one the plan derived, its rank was never qualified, and its
generalization number measures interpolation.

Every one of those is unrecoverable after the fact, and every one is cheap to prevent *before* the
first run — which is now, while there is no run to break.

So: **a training run is not a script that consumes the language. It is an object written in it.**

## 2. `TrainingPlan` — the object

A training run is fully described by a declaration, and the declaration is hashable. Nothing outside
it may influence the run.

```
TrainingPlan
  lane                      MTG_PN_SIM | MTG_PN_EXP | MTG_PN_HYBRID
  representation_plan_hash  the L1 plan (spec §5) that fixes blocks, families, cuts and axis order
  schema_hash               the frozen state schema the plan was derived over
  context_hash              the ExperimentContextManifest in force
  encoder_stages            H1..H3 declarations (spec §8.1b), each bound to plan indices
  exact_constraints         the architectural column of spec §8.3 — no weights, no toggles
  measured_losses           the optimized column — each with a weight SLOT, unfilled
  split                     a SplitDeclaration that passed the §8.4 predicate
  evidence_admission        which trace digests are admissible, and why each qualifies
  refusal_policy            what the model must do outside its trust region
  authority_ceiling         the ceiling of the result, derived — never asserted
```

Three properties are the whole point:

1. **`representation_plan_hash` is required, not optional.** A model whose axis order is not the one
   the plan derived cannot be compared with one that is. Training without a resolvable plan hash is a
   refusal.
2. **Weights are slots, not defaults.** An unfilled slot read at training time is a refusal. There is
   no "reasonable default" for a loss weight; there is a decision record or there is nothing.
3. **The authority ceiling is derived from the admitted evidence**, by the lane algebra in
   `aleph/learn/lanes.py`. A run cannot claim more authority than its inputs carry, because it never
   states its own ceiling.

## 3. Evidence admission — the rule that makes the seal more than a delay

`MTG_PN_SIM` and `MTG_PN_HYBRID` may admit a trajectory only if it was **accepted** by the runtime's
step transaction. Not force-accepted. Not unaccepted-but-plausible. This is the constraint the
`aleph/learn/` seal exists to protect, and it has a precise consequence for the training object:

> `evidence_admission` holds trace digests, and each digest must resolve to an accepted-step
> transaction record. A digest that does not resolve is not a warning and not a skipped sample — it
> makes the whole plan inadmissible.

The reason is asymmetric cost. One unaccepted trajectory inside a training set is invisible in every
metric the model reports, and it invalidates every downstream claim that cites the model. There is no
detector for it after training. The only place it can be caught is admission.

## 4. What must exist before the first training run — in order

| # | Precondition | Owner | Status today |
|---|---|---|---|
| 1 | Frozen state schema at census breadth | `aleph/state/` | landed 2026-07-30 |
| 2 | An L1 representation plan that validates against that schema | `aleph/represent/plan.py` | in progress |
| 3 | An accepted-step transaction that provably rolls back a rejected candidate | `aleph/runtime/` | **the gate** — R3 |
| 4 | At least one accepted native trace, with a digest | `aleph/runtime/` + campaign | R5 |
| 5 | Representation qualification: the error vector measured, by channel, against that trace | `aleph/represent/qualification.py` | sealed (L2) |
| 6 | Tolerances chosen in a decision record, not from the runs that motivated them | PI | not started |
| 7 | A `SplitDeclaration` that passes the §8.4 predicate | `aleph/learn/evaluation.py` | in progress |
| 8 | Loss weights filled from a preregistered pilot, per §10.3 | PI + pilot | not started |

Items 1, 2 and 7 are being built now precisely because they are the ones that do **not** need data.
Item 3 is the true bottleneck, and it is a runtime question rather than a learning question — which is
worth stating plainly, because it means no amount of work in `aleph/learn/` moves the training date.

`MTG_PN_EXP` skips items 3 and 4 and gains four of its own: dataset identity, protocol, calibration,
and split-leakage resolution. It also gains a permanent ceiling, so it can never discharge item 4 on
behalf of the other two lanes.

## 5. What a first run should be, when it is allowed

Not the full model. The first admissible run is the **cheapest thing that can falsify the language
itself**: one block, one family, one cut, one loss, on one accepted trace, reporting the error vector
by channel. If the derived axis order is wrong, or a cut splits something it should not, or a channel
is not measurable, that run finds it for the cost of a pilot rather than for the cost of a paper.

`Baseline-v0` in §10.3 is exactly this and nothing more: widths, ranks and token counts are provisional
and preregistered only *after* pilot memory and error measurements. A first run that is already the
full architecture has skipped its own falsification step.

## 6. The one number this document refuses to state

No date. Item 3 is R3 and item 4 is R5, and both are gated on evidence rather than on effort. A
schedule here would be a prediction dressed as a plan, and the project has already retracted four
claims this session for less.
</content>

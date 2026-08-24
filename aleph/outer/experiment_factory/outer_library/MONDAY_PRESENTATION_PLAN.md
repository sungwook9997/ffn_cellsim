# Outer Library presentation plan — 2026-08-10

## Claim discipline

The Outer Library is currently an independent multimodal biophysical inference
system built because full Aleph parameter sweeps are not yet computationally
available. It is not an Aleph parameter estimator. Every slide must distinguish
three levels of authority:

1. `validated_external_evidence_only`: a task-scoped external result may be
   reported with its exact scope.
2. `training_only` or `blocked_*`: useful representation or diagnostic, but no
   external predictive claim.
3. Aleph authority: currently `none`; no numeric parameter range or physics
   mutation may be emitted.

## Seven-slide narrative

### 1. Why the independent network exists

- Full Aleph sweeps are presently blocked by simulation cost.
- The immediate objective is therefore independent inference from public
  experimental observations.
- A future Aleph adapter is conditional on validated observation operators and
  feasible sweeps; this is a true objective change, not merely task ordering.

### 2. What was built

- 38 source datasets including the HPA embedding source.
- 37 conservative lab/provider groups.
- 193,832 canonical samples and 329,303 canonical observations.
- IF, live-cell morphology, calcium time series, RNA, CyTOF/flow, TFM, PIV and
  related force/state representations.
- 30 task heads. Late fusion is deliberate because aligned cross-modal samples
  are not yet sufficient for a leakage-safe universal latent space.

Source of numbers: `results/adversarial_coverage_audit.json` and
`results/multitask_outer_model_report.json`.

### 3. Architecture and safety boundary

Show one left-to-right diagram:

```text
raw public experiments
  -> protocol/unit/provenance compiler
  -> modality-specific trainable encoders
  -> task-scoped heads
  -> calibration + OOD/refusal
  -> evidence posterior and limitations
  -> authority firewall
       -> reportable external evidence
       -X-> numeric Aleph parameters
```

The important distinction is that OOD/refusal and claim authority are outputs,
not prose added after model evaluation.

### 4. Results that may be shown as positive evidence

Use four separate, scoped examples rather than one aggregate accuracy:

- Piezo1 calcium mechanism direction: independent-lab condition contrast,
  direction consistent in 2/2 active strata; equal-stratum mean difference
  -0.1670 with 95% bootstrap interval [-0.2609, -0.0874].
- Paired migration state: retrospective external-lab AUROC 0.9849, macro-F1
  0.7221, conformal coverage 1.0. State clearly that a new pristine lab is
  still required.
- Fibrotic-state mechanism: externally supported mechanism-level direction;
  not an absolute state classifier.
- Senescence transcriptomic direction: pristine external direction evidence;
  absolute-state OOD remains blocked.

### 5. Failures that demonstrate the gate works

- Calcium per-cell transfer is refused: external macro-F1 0.5499 despite AUROC
  0.7617 and nominal conformal coverage 0.9228. Ranking or interval coverage
  does not rescue a failed task endpoint.
- The consumed IDR0072 IF diagnostic is refused: macro-F1 0.3057, macro average
  precision 0.3169, ECE 0.2370; all five frozen endpoints failed.
- PBMC cell type ranking/calibration is promising, but external uncertainty
  remains blocked after the single-use confirmation failed.
- The sweep gate remains `refused` and cannot mutate physics.

This slide is essential: it shows that the system detects domain shift instead
of converting every experiment into a success claim.

### 6. IDR0168 single-use confirmation integrity demo

The frozen confirmation was refused before prediction because one official
archive is corrupt. Show this as a live integrity/refusal demonstration, not as
a localization result.

Show:

1. the protocol committed before pixels or predictions;
2. the 45-object, 34,187,550,441-byte acquisition receipt;
3. the failed A375/RBM23 object identity and SHA-256;
4. the gzip CRC mismatch and absent native Zarr metadata;
5. the five remote range hashes that match the local object;
6. the final no-replacement refusal with tensor count 0 and prediction count 0.

Beside the integrity evidence show the frozen consequences:

- `failed_source_integrity_and_permanently_consumed`;
- replacement field used: false;
- model predictions made: 0;
- external IF evidence eligible: false;
- Aleph authority: none.

Do not show the seven fields decoded in memory before failure as model evidence.
No tensor was written and the provider-level confirmation is invalid.

### 7. What exists now and what comes next

Current:

- task-scoped external evidence runtime;
- explicit uncertainty/OOD evaluation;
- provenance, single-use holdout and authority enforcement;
- no Aleph numeric parameter authority.

Next:

- reserve and preregister a new independent-provider IF confirmation because
  IDR0168 is permanently consumed by source-integrity failure;
- add preregistered infection-state evidence (IDR0128 is the leading source);
- fill hypoxia, differentiation/stemness, DNA-damage and infection axes;
- obtain aligned independent-lab TFM/PIV with biological replicate IDs;
- validate observation operators before any Aleph parameter proposal.

## Live-demo order

1. Open the frozen IDR0168 protocol and point to the no-replacement rule.
2. Open the acquisition receipt and show 45/45 objects and zero predictions.
3. Run or replay the committed source-integrity audit.
4. Show the immutable decode-refusal JSON.
5. Show the IF head still blocked in the multitask report.
6. Show `results/sweep_gate_report.json` returning `refused`.

The demo should end on the refusal gate. This makes the distinction between
scientific evidence and simulation authority unmistakable.

## Presentation-day checklist

- Regenerate the multitask report after the IDR0168 refusal is integrated.
- Record commit SHA, protocol SHA, acquisition SHA and refusal-result SHA on the
  final slide or appendix. There is deliberately no tensor or prediction SHA.
- Keep a static integrity-evidence panel in case network access fails.
- Do not report training-set or same-provider metrics as external validation.
- Do not call four task-scoped evidence heads a universal cell-state model.
- Do not state that the network reduces Aleph sweep cost until a validated
  observation operator and an executable sweep are demonstrated.

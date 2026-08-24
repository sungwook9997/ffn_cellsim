# Versioned gate contracts

A gate contract pins, in one checked-in file with a content hash:

| part | why it is here |
|---|---|
| `gate_id` | stable identity across runs and renames |
| **threshold** | the number a result is compared against |
| **the OBSERVABLE** | the *exact metric definition* — the half that actually moved twice |
| **configuration / population** | what the gate is measured at (population, mesh, params, myosin on/off) |
| content hash | `gcv1:<sha256[:16]>` over `gate_id + contract_version + configuration_fixed + claims` |

Gate scripts stamp the hash into their own run artifact.
`aleph/outputs/tag_kb/verify_gate_contracts.py` then **refuses a PASS whose contract hash
differs from the FAIL it supersedes**, unless a PI-signed row in `contract_changes.yaml`
authorises the transition.

## Why this exists

The 35-agent audit of 2026-07-25 (`../AUDIT_WHOLE_REPO_2026-07-25.md`,
`../TRAJECTORY_HOW_WE_GOT_TO_AC_2026-07-25.md` failure mode **M2**) found the project's
most-recurring manoeuvre is a gate that changes between a FAIL and its superseding PASS —
**four real instances**, and the guard against it was "a sentence in CLAUDE.md plus PI
sign-off", i.e. cultural. All four are recorded in `contract_changes.yaml` under
`historical_unsigned_changes`.

Two of the four did not move a threshold at all — they moved the **observable** (raw
`max|F|` → projected `max|PF|`) and the **configuration** (the discrete resting-myosin seed
dropped from the static baseline). A guard that hashed only thresholds would have missed
both, which is why `observable.definition` and `configuration_fixed` are inside the hash.

## Files

| file | what it is |
|---|---|
| `GATE_A_resting_baseline.yaml` | GATE A (`max|PF| < 0.21` pN). **The M2 exhibit** — both halves moved at closure. `REPORT_ONLY`: no run artifact exists on disk to enforce against. |
| `GATE_B_sf_motor.yaml` | GATE B SF slice (`nmii_sf_motor`). `BLOCKING`: a stamped artifact exists. 9 claims in 2 claim sets; the quantitative set is a live FAIL (0.155 > 0.01) that the run record itself reports as `BLOCKED`. |
| `contract_changes.yaml` | the change ledger. `changes:` can authorise a transition; `historical_unsigned_changes:` authorises nothing. |
| `../../../outputs/tag_kb/gate_contract.py` | reference implementation: canonicalisation, hash, leaf hashes, claim evaluators, stamping |
| `../../../outputs/tag_kb/verify_gate_contracts.py` | the verifier / CI gate |

Every threshold and observable in both contracts is **transcribed** from what was already
on disk. Nothing was normalised, improved, or re-derived. Thresholds with no derivation on
disk carry `derivation: NONE_ON_DISK` plus a `derivation_note` naming the first appearance
— that is a finding for the PI, not something a transcription pass may fix.

## Contract schema (v1)

```yaml
contract_version: 1
gate_id: GATE-B.sf_motor
enforcement: BLOCKING | REPORT_ONLY     # REPORT_ONLY must state enforcement_reason
contract_hash_declared: "gcv1:…"        # must equal the recomputed hash (SELF_HASH_MISMATCH)

configuration_fixed:                    # HASHED — population, mesh, params, on/off states
  …
configuration_recorded:                 # NOT hashed — per-run knobs, but DIFFED FAIL->PASS
  - {field: explicit_relax_iterations, artifact_path: solver.explicit_relax_iterations}

claims:
  - claim_id: quantitative.inner_solve_converged_relative_to_tension
    claim_set: quantitative              # PASS/FAIL is evaluated per claim set
    verdict_path: verdict_quantitative.…  # where the gate script wrote its own boolean
    observable:                           # HASHED — the exact metric definition
      id: residual_over_tension
      definition: "free-node residual / reported SF tension at the final accepted step"
      units: dimensionless
      reduction: final_step_ratio
    threshold:                            # HASHED
      comparator: "<"
      value: 0.01
      derivation: NONE_ON_DISK            # or a real derivation
    measure: {kind: scalar, path: final.residual_over_tension}
```

`measure.kind` supports `scalar`, `all_of`, `series_max`, `series_min`,
`series_conditional`, `ordered_gt`, `all_true`, `difference_equals`. This is what lets the
verifier **re-derive** each claim from the artifact's own numbers and cross-check it against
the boolean the gate script wrote (`VERDICT_DRIFT`) — a driver literal that drifts off its
contract is then a hard failure instead of an invisible one.

Keys named `notes` or ending `_note` / `_notes` are excluded from the hash, so an editorial
fix does not invalidate every recorded stamp. The hash is over **content, not text**: a
YAML reformat leaves it identical (verified).

## Workflow

```bash
# what does each contract currently hash to?
python aleph/outputs/tag_kb/verify_gate_contracts.py --hashes

# full report / CI gate / one-line summary
python aleph/outputs/tag_kb/verify_gate_contracts.py
python aleph/outputs/tag_kb/verify_gate_contracts.py --gate     # exit 1 on drift
python aleph/outputs/tag_kb/verify_gate_contracts.py --check

# show every claim re-derived from the artifact next to the recorded verdict
python aleph/outputs/tag_kb/verify_gate_contracts.py --show-claims
```

A gate driver stamps itself (see `aleph/scripts/ac_gate_b_sf_motor_native.py`): resolve
the stamp **before** any GPU work — a missing contract then costs milliseconds instead of a
46-minute native run — and put it in the artifact under `gate_contract`.

Artifacts produced before their contract existed can be stamped with an explicit
`RETROACTIVE` provenance, which writes a **sidecar** `*.contract_stamp.json` and leaves the
run record byte-identical. The verifier always reports a retroactive stamp as not-evidence
that the run was gated at the time.

## Changing a contract

1. Write the ledger row in `contract_changes.yaml` **before** the superseding run. A row
   added afterwards to retro-authorise a PASS is the exact manoeuvre this stops; the review
   signal is the commit order.
2. `kind: THRESHOLD | OBSERVABLE | CONFIGURATION_FIXED` needs a real `pi_signature`
   (`PENDING` / `TODO` / `TBD` / empty do not count) plus matching `from_hash` / `to_hash`.
   `kind: CONFIGURATION_RECORDED` (a per-run knob such as relax iterations) needs disclosure
   but not a signature.
3. Edit the contract, then update `contract_hash_declared` (`--hashes` prints the new one).
4. Re-run the gate. Do not edit an artifact to match a contract.

The rationale in a ledger row must stand **without reference to whether the pending run
passes**. If the argument needs the result, it is a fit, not a fix.

## What this does not do

- It cannot audit an **unstamped** gate. Every other gate script in the repo is currently
  unstamped; those are invisible to it until they are wired.
- It sees only artifacts **committed on disk**. A FAIL that is overwritten in place, or
  deleted before commit, leaves no transition for it to check. Gate artifacts should be
  write-once with a run id in the filename.
- Run ordering comes from the stamp's own `stamped_at`, which the run writes about itself.
- It cannot tell a legitimate metric correction from a convenient one. It forces the change
  to be declared, named field by field, and signed — the judgement stays with the PI.

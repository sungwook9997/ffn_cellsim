# V2 Phase 1 6.3b Protrusion-Coupled FA Dynamics — Locked Design

**Date**: 2026-05-04 KST
**Authors**: Claude + Codex design-discussion (5-round adversarial lock,
PI id=809 aggressive debate posture applied, PI id=1008/1057/1161
autonomy + visible deliverable focus)
**Source unit**: design-discussion `topic=v2-layer-2-63b-protrusion-fa-coupling`,
MCP id 1186-1195
**PI ratify status**: full delegation per PI id=939/1008. impl-work uses
this for the 6.3b Sanity Gate doc + code entry.

---

## 0. Scope

6.3b extends 6.3a static FA traction preflight (lock 2026-05-04, MCP id
1070-1081) with **protrusion-coupled effective rates**. 6.3b is a
**deterministic wrapper** around 6.3a; it does NOT reimplement traction
algebra and does NOT change FA state labels.

P1 alpha + 6.3a + 6.4-open ECM preflight + ECM-OL harness extension are
already committed. 6.3b is the next "Separated dynamics modes" sub-unit
before closed-loop Tier 1 entry.

---

## 1. Final Lock Summary

### Inputs
- `adhesions: tuple[FocalAdhesionState, ...]` (existing, attached)
- `centroid_um_xy: tuple[float, float]` (frozen)
- `params: FocalAdhesionDynamicsParameters` (6.3a precedent unchanged)
- `protrusion_registry: Mapping[str, ProtrusionEvent]` (caller-supplied)
- `multipliers: ProtrusionStateMultipliers` (caller-supplied, typed keys)

### Outputs
- `updated_adhesions: tuple[FocalAdhesionState, ...]` (input order preserved)
- `cell_force_nN_xy`, `substrate_reaction_nN_xy`, radial/tangential decomposition
  (6.3a delegate, unchanged contract)
- diagnostics: `linked_missing` count, `reciprocal_missing` count,
  multiplier histogram, `max_effective_rate_per_name`

### Step
1. Validate `multipliers.validate()` upfront — typed key check before any FA work
2. For each FA in input order:
   a. Look up `linked_protrusion_id` in `protrusion_registry`
   b. If missing → raise `FocalAdhesionDynamicsError(failure_kind="linked_protrusion_missing")`
   c. If found, get `protrusion.event_type` state, look up
      `multipliers.lookup(state, rate_name)` for `k_maturity_per_s`,
      `k_bind_per_s`, `k_unbind_per_s`
   d. Compute `effective_rate_<name> = base_rate · multiplier`
   e. Build per-FA effective `FocalAdhesionDynamicsParameters` (replacing rate fields)
3. Compute `dt_rate_violation = dt_fa_s · max(effective_rate over all FAs)`,
   raise if > `_DT_RATE_SAFETY_MARGIN` BEFORE any stepping
4. Per-FA delegate to 6.3a `step_focal_adhesions_static` (or grouped by
   identical effective param tuple as behavior-equivalent optimization)
5. Aggregate results, attach diagnostics, return

### Forbidden
- State label change (no auto transitions, retain 6.3a no-auto-state contract)
- Traction scale change (`traction_scale_nN` from params unchanged)
- Protrusion / ECM / contour write (one-way protrusion → FA only)
- RNG / stochastic sampler (deterministic preflight only)
- `force_candidate_nN` scaling (out of 6.3b scope)
- Closed-loop Tier 1 behavior (FA→ECM, ECM bias, protrusion force on contour)
- 6.3a rate/traction algebra reimplementation (delegate only)

---

## 2. Caller-Supplied Multiplier Table — Typed Keys

### Type contract

```python
from typing import Literal

ProtrusionState = Literal["growing", "stalled", "retracting", "ended"]
RateName = Literal["k_maturity_per_s", "k_bind_per_s", "k_unbind_per_s"]

@dataclass(frozen=True, slots=True)
class ProtrusionStateMultipliers:
    """Caller-supplied per-state multiplier table. NO project default.

    Unknown state or rate keys raise at validate(). Missing known keys
    return neutral 1.0 multiplier. Caller decides which states get which
    boost; 6.3b code makes no biology decisions.
    """
    table: Mapping[ProtrusionState, Mapping[RateName, float]]

    def validate(self) -> None:
        # raise FocalAdhesionDynamicsError("multiplier_table_unknown_key", ...)
        #   if state key not in ProtrusionState enum
        #   if rate key not in RateName enum
        #   if value not finite, is bool, or < 0
        ...

    def lookup(self, state: ProtrusionState, rate_name: RateName) -> float:
        """Return multiplier (1.0 if state or rate_name absent in table)."""
        return self.table.get(state, {}).get(rate_name, 1.0)
```

### Allowed enum values

- **State keys** (exactly): `"growing"`, `"stalled"`, `"retracting"`, `"ended"`
- **Rate keys** (exactly): `"k_maturity_per_s"`, `"k_bind_per_s"`, `"k_unbind_per_s"`

Unknown key → `multiplier_table_unknown_key` failure_kind, raised at
`validate()`. This prevents typo silent-neutral bugs (e.g., `"k_bind"`
vs `"k_bind_per_s"` would otherwise act as 1.0 neutral).

### Missing known keys

Missing state in table → all rates for that state default to 1.0
(neutral). Missing rate under a state → that rate defaults to 1.0
(neutral). "Absent means neutral, but spelling enforced."

---

## 3. Reference biology table (doc only, NOT code)

The following multiplier values reflect biological intuition but are
**reference guidance only**. Caller may use as-is or modify per cell
type / experimental condition. These do NOT live in code.

| state        | k_maturity | k_bind | k_unbind | note                              |
|--------------|-----------:|-------:|---------:|-----------------------------------|
| `growing`    | 2.0        | 1.5    | 1.0      | active protrusion boosts both     |
| `stalled`    | 1.5        | 1.0    | 1.0      | maturity-only boost by default    |
| `retracting` | 1.0        | 1.0    | 1.0      | neutral (no boost, no penalty)    |
| `ended`      | 1.0        | 1.0    | 1.0      | base 6.3a behavior                |

Tests are free to use any multiplier table. Production callers should
document the table they choose with literature provenance.

---

## 4. Link Consistency — Asymmetric

`FocalAdhesionState.linked_protrusion_id` is the **authoritative read
path**. 6.3b reads FA → registry; it does NOT require
`ProtrusionEvent.associated_adhesion_ids` reciprocity.

Mismatch handling:
- **linked id missing from registry** → explicit failure
  `FocalAdhesionDynamicsError(failure_kind="linked_protrusion_missing", ...)`
- **linked id exists but FA not in `associated_adhesion_ids`** →
  diagnostic warning/count (`reciprocal_missing` in diagnostics dict),
  not failure
- **protrusion lists FA but FA does not link back** → ignored by
  one-way read, diagnostic only

This handles partial/manual segmentation data without forcing
reciprocity in raw inputs.

---

## 5. Sanity Gate 6 Items (impl writes in module docstring)

1. **Units**: rate multiplier dimensionless, effective rate 1/s, dt·rate
   dimensionless, traction nN (6.3a delegate)
2. **Boundary**: typed key validation upfront, empty FA list OK,
   multiplier finite/non-bool/≥0, linked id missing raises
3. **Conservation**: 6.3a Newton-3 delegate (per-FA action-reaction
   exact zero, aggregate Newton-3)
4. **Numerical**: `dt_fa_s · max(effective_rate over all FAs) ≤
   _DT_RATE_SAFETY_MARGIN`, no auto-shrink, raise before stepping
5. **Sign**: multiplier ≥ 0 (no negative rates), multiplier 1.0 =
   neutral = base 6.3a, multiplier 0 = effective rate 0 = no state change
6. **Measurement**: 6.3a output unchanged contract (cell-on-substrate
   inward radial-to-centroid traction, decomposition reconstructs
   original); diagnostics added only

---

## 6. Test Catalog (~12)

- `test_neutral_multipliers_match_6_3a_baseline` — multiplier 1.0
  everywhere → output identical to 6.3a (regression)
- `test_growing_multiplier_boosts_maturity_and_bind`
- `test_retracting_neutral_no_boost`
- `test_unknown_state_key_raises_multiplier_table_unknown_key`
- `test_unknown_rate_key_raises_multiplier_table_unknown_key`
- `test_negative_multiplier_raises_at_validate`
- `test_bool_multiplier_value_raises_at_validate`
- `test_missing_linked_protrusion_id_raises_linked_protrusion_missing`
- `test_reciprocal_missing_records_diagnostic_not_raise`
- `test_dt_rate_violation_uses_max_effective_rate_over_all_fas`
- `test_per_fa_input_order_preserved_in_output`
- `test_no_state_label_change_after_step`

---

## 7. Files

- `docs/v2/v2_63b_protrusion_coupling_locked.md` (this file, source of truth)
- `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py` (new)
  - `ProtrusionStateMultipliers` dataclass
  - `step_protrusion_coupled_focal_adhesions(adhesions, centroid, params, registry, multipliers)` function
  - `FocalAdhesionDynamicsError(failure_kind="multiplier_table_unknown_key" | "linked_protrusion_missing" | inherited from 6.3a)` extension
- `tests/v2/test_v2_protrusion_coupled_focal_adhesion.py` (new, ~12 tests)
- Optional: `acs/v2/__init__.py` exports

---

## 8. References

- 6.3a static FA traction preflight lock: design-discussion `topic=v2-layer-2-separated-dynamics`, MCP id 1070-1081
- 6.3a Sanity Gate: `docs/v2/v2_focal_adhesion_dynamics_sanity_gate.md` (impl-work)
- Forward roadmap: `docs/v2/v2_phase1_forward_roadmap.md`
- Consolidated plan: `docs/v2/v2_phase1_plan_consolidated.md` §6.3b
- Adversarial debate posture: memory `feedback_aggressive_design_debate.md`,
  PI id=809
- Cadence rule: memory `cadence_promise_must_send_even_when_idle.md`,
  PI id=1025/1031
- ECM-OL persistent run pattern (precedent): `runs/20260504T035716Z_ecm_ol/`
- P1 alpha persistent run pattern (precedent): `runs/20260503T1833Z_p1_alpha/`

---

## 9. Cross-room dispatch

This file is the design-team input to implementation-work for:
1. impl Claude writes 6.3b Sanity Gate doc (`docs/v2/v2_63b_protrusion_coupled_fa_sanity_gate.md`) from this lock
2. impl Codex review (5 focus per Sanity Gate items: typed keys / per-FA delegation / linked id missing / reciprocal diagnostic / dt-rate gate)
3. On Sanity Gate PASS: 6.3b code in `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py` + tests
4. After commit: persistent artifact run script (similar to `scripts/run_ecm_ol_harness.py`) for PI visible deliverable

Rounds 1-5 of the design lock are MCP id 1186-1195.

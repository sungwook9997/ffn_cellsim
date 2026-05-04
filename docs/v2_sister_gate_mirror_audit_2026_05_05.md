# V2 Sister-Gate-Mirror Audit (2026-05-05) — read-only findings + next-touch list

**Date**: 2026-05-05 KST
**Audit unit**: read-only sister-gate-mirror exercise of the
freshly-promoted Step 6 rule (memory rule
`design_note_pre_commit_batch.md`, promoted from 5-step to 6-step
per Codex `id=1428`) against existing `acs/v2/dynamics/` modules.
**Routing**: per Codex `id=1430` (audit scope) and `id=1432`
(routing decision).
**Author**: implementation-work Claude.
**Status**: audit complete; this note records findings and the
next-touch list. B1 follow-up unit is a separate cycle (brief →
lock → Sanity Gate); A1 cosmetic stays opportunistic.

---

## Scope

Per `id=1430`:
- Compare HB#3 / HB#4 / Phase D already-touched surfaces for
  remaining drift: validation order, return dataclasses, package
  exports, `__all__`, failure-kind naming.
- Spot-check older adjacent v2 dynamics APIs only where they form
  a plausible sister/paired surface: `step_focal_adhesions_static`,
  `step_protrusion_coupled_focal_adhesions`,
  `accumulate_prescribed_traction`.
- Read-only unless a concrete defect appears.

---

## Findings summary

| Finding | Severity | Status |
|---|---|---|
| A1 — HB#3 `_BOUNDARY_TOL_RELATIVE = 1e-12` missing `Final[float]` annotation that HB#4 has (`acs/v2/dynamics/fa_to_ecm_scattering.py:70` vs `acs/v2/dynamics/ecm_to_fa_bias.py:75`) | Minor cosmetic | **Next-touch list**: fix opportunistically when HB#3 is touched for substantive work; do not patch standalone (per `id=1432`) |
| A2 — HB#4 has `_validate_fa_position_and_schema` helper (`ecm_to_fa_bias.py:155`); HB#3 has same logic inline (`fa_to_ecm_scattering.py:180-215`) | Style drift | **No action**; route to Phase E audit if it touches scatter |
| **B1 — `FocalAdhesionDynamicsResult` (`focal_adhesion.py:153`) is `@dataclass` plain (NOT frozen+slots) AND has `diagnostics: dict = field(default_factory=dict)` loose untyped field** | **Concrete defect risk** | **Separate cycle**: brief → lock → Sanity Gate (per `id=1432`); brief in flight at `docs/v2_focal_adhesion_dynamics_result_typed_brief.md` |
| B2 — `FocalAdhesionDynamicsParameters` (`focal_adhesion.py:66`) and `ProtrusionStateMultipliers` (`protrusion_coupled_focal_adhesion.py:82`) are `@dataclass(frozen=True)` without `slots=True` | Minor optimization | **No action**; `slots` is performance, not correctness |
| B3 — `protrusion_coupled_focal_adhesion.py` reuses `FocalAdhesionDynamicsError` from 6.3a | Not a defect (deliberate composition) | No action |
| B4 — `accumulate_prescribed_traction` + `apply_prescribed_*` use a different validation pattern (separate `_validate_*` helpers + `next_ecm.validate()` after operations) | Not a defect (different surface) | No action |

---

## Mutual-mirror coverage results

### HB#3 / HB#4 / Phase D — clean ✓

All four sister-gate-mirror layers verified clean across HB#3 / HB#4 / Phase D:

- **Validation order**: HB#3 line 180-215 + HB#4 line 173-198 share `len(position) check → unpack → non_finite_position → fa.validate() → out-of-grid raise` post-`e714d12`. Phase D delegates to both, no duplicate validation.
- **Dataclass discipline**: HB#4 `ECMSampledAtFAs` + `ECMToFABiasResult` and Phase D `FAToECMResponseResult` + `PhaseDNoOpStepResult` all use `@dataclass(frozen=True, slots=True)`.
- **API exports**: All HB#3 + HB#4 + Phase D public symbols present at both `acs/v2/__init__.py` and `acs/v2/dynamics/__init__.py` per `b1b2b92`.
- **Failure-kind discipline**: HB#3 `fa_position_outside_ecm_grid` + `non_finite_fa_position` + `non_finite_fa_traction`; HB#4 `fa_bias_position_outside_ecm_grid` + `non_finite_fa_position` (distinct kinds for distinct operations).

### Older v2 modules — exports clean ✓, dataclass discipline drift in 6.3a/6.3b

API exports at both `acs/v2/__init__.py` and `acs/v2/dynamics/__init__.py` for:
- `FocalAdhesionDynamicsError`, `FocalAdhesionDynamicsParameters`, `FocalAdhesionDynamicsResult`, `step_focal_adhesions_static` ✓
- `ProtrusionStateMultipliers`, `step_protrusion_coupled_focal_adhesions` ✓
- `ECMOpenLoopError`, `accumulate_prescribed_traction`, `apply_prescribed_*` ✓

Failure-kind discipline (all older modules use `class FooError(ValueError)` + `failure_kind` attribute pattern) ✓.

Dataclass discipline drift confined to 6.3a + 6.3b result/parameter shapes (Findings B1 + B2).

---

## Next-touch list (opportunistic A1 fix)

When HB#3 (`acs/v2/dynamics/fa_to_ecm_scattering.py`) is next touched for substantive work (not a standalone cosmetic patch), apply the one-line fix:

```python
# Before
_BOUNDARY_TOL_RELATIVE = 1e-12

# After (mirror HB#4 line 75)
_BOUNDARY_TOL_RELATIVE: Final[float] = 1e-12
```

Add `Final` to imports if not already present:

```python
from typing import Final
```

This is a pure-typing improvement — no behavior change. Sister-pattern alignment with HB#4.

**Why "opportunistic" not "now"** (per `id=1432`): A1 touches a physics/numerics module; the value is already semantically locked and mirrored at the constant level (`1e-12` matches HB#4 exactly). A standalone audit commit just to add `Final[float]` adds noise to the git log without adding correctness. The fix lands when HB#3 is opened for any non-trivial change.

---

## B1 brief unit (in flight)

The B1 finding (`FocalAdhesionDynamicsResult` mutable + loose `diagnostics: dict`) is being opened as a separate cycle per `id=1432` routing:

- Brief: `docs/v2_focal_adhesion_dynamics_result_typed_brief.md` (impl-work writing now)
- Cross-room dispatch to design-discussion (per (B-trigger) precedent)
- Adversarial round → lock → Sanity Gate → code+tests

**Codex `id=1432` nuance preserved**: HB#4 still has its own `diagnostics_dict` surface (`ECMToFABiasResult.diagnostics_dict: dict[str, float | int | str]`); the B1 design question is **NOT** "ban all dict diagnostics" but rather "design diagnostic schema with explicit reasoning per pattern" + "frozen/slots immutability before Phase E composition picks up the older shapes".

The brief enumerates four design questions (Q1 typed schema for shared 6.3a fields / Q2 6.3b extension shape / Q3 frozen+slots adoption / Q4 backward-compat dict-access shim) without pre-baking the abstraction.

---

## Step 6 sister-gate-mirror rule — first activation post-promotion

This audit unit is the first activation of the Step 6 sister-gate-mirror check after promotion to memory rule per Codex `id=1428` ACK. Findings:

- Rule successfully surfaced one concrete defect (B1) and one cosmetic drift (A1) that pre-Step-6 batches would not have flagged.
- Rule's "vacuous-first-member scope clarification" applied correctly: `accumulate_prescribed_traction` (Finding B4) was correctly identified as having a different surface (no sister to mirror) and skipped.
- Wording tweak ("example pattern, not universal required sequence") prevented over-application: 6.3a/6.3b validation order pre-dates HB#3+#4 structural-len-check pattern but doesn't unpack `position_um_xy`, so the rule correctly skipped that line of analysis.

The rule earned its memory-store position. Future Phase E composition cycles will exercise it again on the larger composed surface.

---

## References

- Codex routing: `id=1430` (audit scope), `id=1432` (B1 separate cycle, A1 next-touch list)
- Step 6 promotion: `id=1428` ACK + memory rule `design_note_pre_commit_batch.md` (6-step)
- HB#3 / HB#4 / Phase D locks + Sanity Gates + impl + tests:
  - `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `docs/v2_phase_d_no_op_scaffolding_locked.md`
  - `docs/v2_fa_to_ecm_scattering_sanity_gate.md`
  - `docs/v2_ecm_to_fa_bias_sanity_gate.md`
  - `docs/v2_phase_d_no_op_scaffolding_sanity_gate.md`
  - `acs/v2/dynamics/fa_to_ecm_scattering.py`
  - `acs/v2/dynamics/ecm_to_fa_bias.py`
  - `acs/v2/dynamics/closed_loop_phase_d.py`
- Older v2 modules audited:
  - `acs/v2/dynamics/focal_adhesion.py` (6.3a)
  - `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py` (6.3b)
  - `acs/v2/dynamics/ecm_open_loop.py`
- Sister-pattern fix precedent: `e714d12` (HB#3+#4 unpacking-before-validate fix from `id=1385`)

# V2 Closed-Loop ECM Gate — Phase E v2 Composition Step 2 — Design-Discussion Brief (opening, NOT a lock)

**Date**: 2026-05-05 KST
**Status**: design-discussion brief, **opening position only — NOT a lock,
NOT a code commitment**.
**Author**: implementation-work Claude, drafted under (B-trigger) precedent
established by HB#3 / HB#4 v1 / Phase D no-op / B1 / HB#1+#2 / HB#5 /
Phase E v1 / Item 5 sweep harness / HB#4-active step 1 (all locked +
Sanity Gated + implemented + Codex PASS in this session, ~11 hours,
10 cycles).
**Source**: locked phased plan
`docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md` §1 (Phase E
satisfaction is the original v2 target) + locked HB#4-active step 1
`docs/v2/v2_hard_blocker_4_active_locked.md` §0 (Phase E v2 step 2 is the
explicit next-cycle target) + Codex `id=1631` (b) routing
documentation + PI watchdog `id=1635` directive (do not standby; pick
next small unit) + Codex `id=1636` agreement (Phase E v2 step 2 brief
only as first action, not code).
**Unblock condition** (now satisfied): HB#4-active step 1 lock + impl
+ PASS at `2207566` + `d52af7d` + `fed370b` (Codex `id=1629`); Phase E
v1 composition + Item 5 sweep harness + 8 prior sister units all
sealed.

**Hard contract**: this brief is the **opening position** for the
adversarial design round on Phase E v2 composition step 2
(`step_closed_loop_phase_e_v2`). It does NOT lock the composition law,
does NOT commit the implementation, and does NOT claim Phase E v2
satisfaction. **Phase E v2 satisfaction is achieved only AFTER the
locked design + Sanity Gate PASS + code+tests PASS land** — this
brief opens the cycle that produces those artifacts.

The wording boundary that Phase E v1 + HB#4-active step 1 preserved
("ECM-side evidence" / "active bias law evidence channel" / NOT v2
satisfaction) **graduates to Phase E v2 satisfaction once this cycle
seals**. The brief surfaces the design questions that must be
adversarially resolved BEFORE that graduation.

---

## 0. Why this brief exists + Codex `id=1636` verification preferences

Phase E v1 composition (`step_closed_loop_phase_e_v1`) provides
**ECM-side evidence only** because it consumes
`compute_ecm_to_fa_bias_neutral` (HB#4 v1, all-1.0 multipliers — no
real FA→ECM→FA closure). HB#4-active step 1
(`compute_ecm_to_fa_bias_active`) provides the active bias law
**evidence channel only** because no composition consumes it yet.

**Phase E v2 step 2 closes the loop**: a composition wrapper
`step_closed_loop_phase_e_v2` that swaps neutral → active in the v1
sequence, demonstrating the first integrated FA→ECM→FA closure.

### Phase E v2 satisfaction gate (the wording boundary that closes here)

When this cycle seals (lock + Sanity Gate PASS + code+tests PASS),
**Phase E v2 satisfaction** is achieved. Until then, the wording
boundary preserved across 10 prior cycles remains in force:

| Cycle | Wording boundary | Status |
|---|---|---|
| Phase E v1 (cycle 7) | "ECM-side composition evidence ONLY, NOT full closed-loop ECM gate satisfaction" | preserved |
| HB#4-active step 1 (cycle 10) | "active ECM→FA bias law / evidence channel, NOT Phase E v2 satisfaction" | preserved |
| Phase E v2 step 2 (THIS cycle) | "v2 composition satisfaction ONLY after code+tests" | **graduates on seal** |

### Codex `id=1636` verification preferences enforced throughout

Per Codex `id=1636` post-watchdog routing message:

1. **First files for this brief stage**: `docs/v2/v2_phase_e_v2_composition_brief.md` (NEW); source-to-inspect `acs/v2/dynamics/closed_loop_phase_e.py` (Phase E v1 composition, line 171 `step_closed_loop_phase_e_v1`) + `acs/v2/dynamics/ecm_to_fa_bias.py` (HB#4 v1 + HB#4-active step 1, lines 371 + 454)
2. **Verification commands for brief commit**:
   - `git show --check <brief-commit>` (clean)
   - `grep -nE "step_closed_loop_phase_e_v1|compute_ecm_to_fa_bias_neutral|compute_ecm_to_fa_bias_active"` (all 3 names referenced)
   - **Explicit wording boundary check**: brief states `v2 composition satisfaction ONLY after code+tests` — meta-grep on this exact phrase
3. **No PI decision needed** unless design proposes a gate-contract change or Sanity Gate FAIL

---

## 1. What is already locked (do not redebate here)

Per locked phased plan + completed unit locks:

- **Phase E v1 composition** (`step_closed_loop_phase_e_v1` at
  `closed_loop_phase_e.py:171`):
  - Sub-call sequence: HB#3 scatter → HB#1+#2 orientation update →
    HB#4 neutral bias readout → HB#5 Lyapunov metric
  - Identity invariant Y4: `result.updated_ecm is
    result.orientation_response.updated_ecm` (object identity, no
    copy)
  - No wrapper failure kinds Y17 (sub-call errors propagate verbatim)
  - Forward-compat documented for HB#4-active variant: bias readout
    runs on post-HB#1+#2 ECM (correct vantage for v2 swap)
- **HB#4-active step 1** (`compute_ecm_to_fa_bias_active` at
  `ecm_to_fa_bias.py:454+`):
  - `(adhesions, ecm, *, k_active) -> ECMToFABiasResult`
  - Pure read-only, deviatoric Rayleigh score, all 3 rates same
    multiplier (Y4 scaffolding)
  - 4-kind `FAToECMBiasFailureKind` literal: 2 sampler kinds + 2
    active-law kinds
- **All 5 hard blockers + Phase D no-op + B1 audit fix + Step 6
  sister-gate-mirror memory rule promotion**: ✓ sealed

---

## 2. Five design questions (NO defaults pre-committed)

### Q1 — Result type

- (a) Reuse `PhaseEStepResult` from v1 (forward-compat per HB#4-active
  Y15 sister; same 5-field shape: `traction_density_xy`,
  `orientation_response`, `ecm_to_fa_bias`, `lyapunov_metric`,
  `updated_ecm`). v1 already documents it as forward-compat for v2.
- (b) New `PhaseEV2StepResult` dataclass (explicit type marker for v2
  satisfaction; allows v2-specific diagnostic fields like `k_active`
  metadata).
- (c) New `PhaseEV2StepResult` that wraps `PhaseEStepResult` (composition
  pattern; v2 result has a `.v1_result` field plus v2-specific extras).

### Q2 — Module location

- (a) Modify `closed_loop_phase_e.py` to add `step_closed_loop_phase_e_v2`
  alongside `step_closed_loop_phase_e_v1` (sister-pattern with
  `ecm_to_fa_bias.py` carrying both `compute_ecm_to_fa_bias_neutral` and
  `compute_ecm_to_fa_bias_active`).
- (b) Create new `closed_loop_phase_e_v2.py` module (sister-pattern with
  `closed_loop_phase_e_sweep.py` — separate module per
  variant/responsibility).

### Q3 — `k_active` parameter discipline

- (a) Required, no default (sister with HB#4-active Y3 — Magic-Number
  Block guard against silent default tuning).
- (b) Required, no default, AND wrapper validates `k_active` BEFORE
  HB#3/HB#1+#2 sub-calls fire (cheap-parameter-failure-first sister
  with HB#4-active Y12).
- (c) Optional with explicit `None` sentinel to fall back to v1 neutral
  composition (forbidden for satisfaction lock; would break the
  satisfaction-on-seal wording boundary).

### Q4 — Phase E v2 satisfaction wording in lock

- (a) Lock explicitly states "Phase E v2 satisfaction is achieved by
  this composition" — graduates the wording boundary; tests use
  `provides_phase_e_v2_satisfaction_*_evidence` patterns (note:
  `_evidence` retained for honesty — satisfaction is achieved by code
  PASS, not lock PASS).
- (b) Lock states "Phase E v2 composition step 2" only; defer the
  satisfaction graduation to a later post-impl artifact.
- (c) Hybrid: lock states "first integrated FA→ECM→FA closure" without
  invoking either "satisfaction" or "evidence" wording.

### Q5 — New invariants beyond v1's Y4 identity

- (a) None — v2 inherits v1's identity invariant Y4 verbatim
  (`result.updated_ecm is result.orientation_response.updated_ecm`).
- (b) Add v2-specific invariant: the `ecm_to_fa_bias.multipliers_per_fa`
  field is **non-uniform** for non-isotropic post-HB#1+#2 ECM
  (otherwise v2 collapses to v1 silently). This catches the case where
  the active bias law doesn't actually differentiate FAs at the chosen
  vantage.
- (c) Add v2-specific invariant: object-identity preservation of the
  HB#3 traction_density_xy through to HB#5 Lyapunov metric (v1 already
  has this de-facto; v2 makes it explicit).

### Q5.1 (sub-decision under Q5) — Anti-collapse meta-test wording

If Q5 picks (b), test naming should be
`provides_v2_composition_non_collapse_evidence_*` (NOT `satisfies_*`)
to preserve the test-name wording-boundary discipline established in
Phase E v1 / Item 5 / HB#4-active step 1.

---

## 3. Recommended opening positions (NOT decisions; design-discussion
adversarial round must converge them)

| Question | Recommendation | Rationale |
|---|---|---|
| Q1 | (a) Reuse `PhaseEStepResult` | v1 already documents forward-compat for v2; HB#4-active Y15 sister-pattern (no new result dataclass); minimal API surface |
| Q2 | (a) Modify `closed_loop_phase_e.py` | sister with `ecm_to_fa_bias.py` carrying both variants; Phase E v1 + v2 are TWO STEPS of the same closed-loop ECM gate composition, not separate concerns |
| Q3 | (b) Required, no default, validate before sub-calls | cheap-parameter-failure-first sister with HB#4-active Y12 — invalid `k_active` should fail before HB#3 scatter runs (which is O(N_FA) work) |
| Q4 | (a) "Phase E v2 satisfaction achieved by composition" | this is THE point of step 2; deferring satisfaction wording weakens the milestone |
| Q5 | (b) Anti-collapse invariant | a non-trivial invariant that distinguishes v2 from "v1 with active bias readout that happens to be neutral"; provides the most adversarial-defensible v2 satisfaction claim |

These are **opening positions**, NOT locked decisions. design-discussion-pane Codex round 1 may propose alternatives.

---

## 4. Allowed / Forbidden separation

### Allowed (this lock target)

- Public function `step_closed_loop_phase_e_v2(adhesions, ecm, dt_s,
  *, k_active, traction_ref_nN_per_um2=..., k_orient_per_s=...)`
- Sub-call sequence (Y4 deterministic order): HB#3 scatter → HB#1+#2
  orientation update → **HB#4-active** bias readout → HB#5 Lyapunov
  metric
- `k_active` validated by wrapper BEFORE sub-calls fire (Q3 (b)
  recommendation)
- `PhaseEStepResult` reuse (Q1 (a) recommendation) OR new dataclass
  (Q1 (b)/(c)) — pending design-discussion lock
- Identity invariant Y4 preserved verbatim from v1
- Anti-collapse invariant Q5 (b) tested via meta-test pattern

### Forbidden

- **Mutation of any input state** (HB#4-active sister-pattern + Phase
  E v1 sister-pattern; both wrappers are pure)
- **Optional fallback to neutral composition** (Q3 (c) excluded — would
  break satisfaction wording boundary)
- **Reading `effective_stiffness` / `stiffness_kpa` / `fiber_density` /
  `ligand_density` / `accumulated_traction_nNs_per_um2`** from inside
  `step_closed_loop_phase_e_v2` (locked phased plan §3
  effective-stiffness guard; HB#4-active Y8 inheritance — function-scoped
  AST + string guard)
- **`step_closed_loop_phase_e_v3` placeholder, hooks, or strategy**
  (forbidden per HB#4-active Y3 sister; v3 is a separate later cycle)
- **New result dataclass field that duplicates v1 or shadows
  HB#4-active diagnostics** (Y5 schema-corrected fields applies; v2
  must not introduce phantom field names)
- **Wrapper-level failure kinds beyond what v1 + HB#4-active already
  surface** (Y17 sister-pattern — sub-call errors propagate verbatim;
  no `phase_e_v2_*` failure_kind invented)
- **`k_active` default value or environment-variable read**
  (HB#4-active Y3 inheritance)
- **Module-docstring claim of v2 satisfaction graduation before
  code+tests PASS** (wording boundary discipline; satisfaction wording
  in module source / function docstring / test names reserved for
  after impl PASS — paraphrased here to avoid the forbidden literal
  substring per Item 5 / HB#4-active in-flight test catch sister-pattern)

---

## 5. Reference precedent

- **Phase E v1 composition lock + impl** (sister cycle): `7b1d3cf`
  composition + `c5634d5` Item 5 lock + sweep harness sister
- **HB#4-active step 1 lock + impl** (the active bias law that v2
  consumes): `fed370b` lock + `d52af7d` code + `2207566` BLOCKER fix
- **HB#4 v1 lock + impl** (the neutral bias readout that v1 consumes):
  `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
- **HB#1+#2 lock** (the orientation update that both v1 and v2
  consume): `docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md`
- **HB#3 lock** (the scatter that both v1 and v2 consume):
  `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
- **HB#5 lock** (the Lyapunov metric that both v1 and v2 consume):
  `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md`
- **Item 5 sweep harness lock** (Phase E v1 sensitivity evidence):
  `docs/v2/v2_item_5_sweep_harness_locked.md` — note: a future Item 5 v2
  variant that exercises HB#4-active is OUT OF SCOPE for this brief

### Biology / literature context

- HB#4-active step 1 deviatoric Rayleigh law has biological grounding
  via Hall 2016 (RhoA-mediated FA orientation alignment with stiffness
  gradients), Trichet 2012 (substrate stiffness as cell INPUT), and
  Stylianopoulos 2018 (anisotropic mechanical priming in ECM-stroma
  coupling) — full citations in HB#4-active step 1 lock §6.
- Phase E v2 composition does NOT introduce new biology beyond
  HB#4-active step 1 + the v1 sub-calls; this is a structural
  composition cycle, not a new constitutive law.

---

## 6. What is *out of scope* for this brief

- **Phase E v3 (per-rate selectivity)**: requires per-rate literature
  anchors; deferred to separate cycle.
- **`effective_stiffness` law**: HB#4-active step 1 chose orientation-only
  path → v2 step 2 does NOT consume `effective_stiffness`; helper
  remains deferred.
- **Item 5 v2 sweep harness** (Phase E v2 sensitivity): would augment
  Item 5 Y1 IC `T = 0.5·I` with HB#4-active variant; HB#4-active Y1
  shows this is automatic-neutral, but a v2 sweep harness is its own
  cycle.
- **GPU implementation** — CPU-only per Phase E v1 sister-pattern.
- **PI experimental data fitting** — Hard Rule 1 forbids; literature
  anchors only.
- **Composition with FA dynamics** (`step_focal_adhesions_static` etc.):
  Phase E v2 step 2 is the FA→ECM closure half; the ECM→FA effect of
  bias multipliers on FA assembly/disassembly rates is Phase E v3+
  scope (requires composition with `step_focal_adhesions_static` AND
  per-rate selectivity AND time-stepping discipline).

---

## 7. Sanity Gate items preview (for post-lock cycle)

After design-discussion lock, the Sanity Gate doc must cover
6 standard items + the Phase E v2 composition specifics:

1. **Units**: `dt_s` seconds, `k_active` dimensionless, all sub-call
   units inherited from v1
2. **Boundary**: empty FA list, `dt_s = 0` no-op, isotropic ECM IC
   `T = 0.5·I` (HB#4-active Y1 automatic-neutral case — meta-test must
   exercise this AND distinguish from collapse-to-v1 per Q5 (b)
   anti-collapse invariant)
3. **Conservation**: pure read-only wrapper; identity invariant Y4
   preserved
4. **Numerical**: `np.float64` enforced; no new tolerance; sub-call
   tolerances inherited
5. **Sign**: bias multiplier > 0 always (HB#4-active inherited);
   non-uniform under non-isotropic ECM (Q5 (b))
6. **Measurement-protocol** (Hard Rule 11): the v2 composition is the
   FIRST artifact where "FA→ECM→FA closure" is structurally complete;
   meta-test must verify the closure is non-trivial (not "active law
   happens to be neutral at the chosen vantage")

### Anti-collapse meta-test (Q5 (b) recommendation)

For non-isotropic post-HB#1+#2 ECM (e.g., after several time steps of
HB#1+#2 evolution with directional FA traction), assert:
`set(result.ecm_to_fa_bias.multipliers_per_fa.flatten()) != {1.0}` —
v2 composition truly differentiates FAs (the active law is non-trivial
at the closure vantage). Without this guard, a future bug in HB#1+#2
or HB#4-active that silently restores isotropy would let v2 satisfy
the "structural" composition surface while collapsing back to v1
behavior.

---

## 8. Process expectation

- 3-4 round adversarial cap (sister with HB#4-active step 1 / Item 5
  sweep / Phase E v1 design rounds)
- design-discussion-pane Claude opens round 1 by adopting opening
  positions OR proposing alternatives
- design-discussion-pane Codex provides adversarial review on each
  round
- SEAL ack from Codex → lock artifact commit
- impl-work-pane resumes after lock: Sanity Gate → impl Codex review →
  code commit → impl Codex review → milestone

**ETA full cycle**: ~1.5-2h (mirrors Phase E v1 composition timing per
`7b1d3cf` precedent). Brief commit ~15 min; design-discussion lock
~30-45 min; Sanity Gate doc ~20-30 min; code commit ~30-40 min.

---

## 9. References

- Brief: this file (`docs/v2/v2_phase_e_v2_composition_brief.md`)
- Parent locked plan:
  `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- HB#4-active step 1 lock + impl (the active bias law v2 consumes):
  `docs/v2/v2_hard_blocker_4_active_locked.md` + commit `fed370b`
  (lock) + `d52af7d` (code) + `2207566` (BLOCKER fix)
- Phase E v1 composition lock + impl (sister structural pattern):
  `docs/v2/v2_phase_e_composition_locked.md` + commit `7b1d3cf`
- HB#4 v1 + sampler lock + impl (the v1 neutral readout v2 swaps out):
  `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
- HB#1+#2 + HB#3 + HB#5 locks (sub-calls v2 inherits from v1):
  `docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md` +
  `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md` +
  `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md`
- Item 5 sweep harness lock (Phase E v1 sensitivity evidence; future v2
  variant deferred): `docs/v2/v2_item_5_sweep_harness_locked.md`
- Source-to-inspect (per Codex `id=1636`):
  - `acs/v2/dynamics/closed_loop_phase_e.py:171` —
    `step_closed_loop_phase_e_v1`
  - `acs/v2/dynamics/ecm_to_fa_bias.py:371` —
    `compute_ecm_to_fa_bias_neutral` (v1 consumes)
  - `acs/v2/dynamics/ecm_to_fa_bias.py:454+` —
    `compute_ecm_to_fa_bias_active` (v2 consumes; HB#4-active step 1)
- Routing chain: Codex `id=1604` (a1) → HB#4-active step 1 →
  Codex `id=1631` (b) pause-then-defer → PI watchdog `id=1635` override
  → Codex `id=1636` agreement (Phase E v2 step 2 brief only as first
  action)
- Adversarial debate posture: memory
  `feedback_aggressive_design_debate.md`, PI id=809
- Hard Rule 10 inline derivation: memory
  `rule10_unit_derivation_in_docs.md`, Codex id=1234
- Hard Rule 11 wording-boundary meta-test: memory
  `hard_rule_11_wording_boundary_meta_test.md` (Phase E v2 composition
  satisfaction is the FIFTH catch-family member if anti-collapse
  invariant Q5 (b) is locked)
- Step 6 pre-commit batch: memory `design_note_pre_commit_batch.md`,
  Codex id=1428

---

## 10. Sister-pattern + sister-gate-mirror declaration

| Layer | Phase E v1 (sister) | Phase E v2 step 2 (this unit) | Notes |
|---|---|---|---|
| **Code** | `step_closed_loop_phase_e_v1`, sub-call sequence (HB#3 → HB#1+#2 → HB#4 v1 → HB#5) | `step_closed_loop_phase_e_v2`, sub-call sequence (HB#3 → HB#1+#2 → **HB#4-active** → HB#5); structural swap only | parallel structure with single sub-call substitution |
| **Design** | ECM-side evidence (HB#4 v1 neutral) | First integrated FA→ECM→FA closure (HB#4-active step 1 active) | wording boundary graduates from "evidence" to "satisfaction" on seal |
| **API** | Public `step_closed_loop_phase_e_v1` + `PhaseEStepResult` | Public `step_closed_loop_phase_e_v2` + `PhaseEStepResult` reuse OR new (Q1 design question) | mirrors HB#4 v1 + HB#4-active naming convention |
| **Failure-kind** | No wrapper failure kinds (Y17) | No wrapper failure kinds (Y17 inheritance) | sub-call errors propagate verbatim from HB#3/HB#1+#2/HB#4-active/HB#5 |
| **Wording boundary** | "evidence not satisfaction" (Y2) | "satisfaction on seal" (graduates Y2 to satisfaction) | new test-name pattern `provides_phase_e_v2_satisfaction_*_evidence` (note: `_evidence` retained for honesty — satisfaction proven by code PASS not lock PASS) |
| **Identity invariant** | Y4: `updated_ecm is orientation_response.updated_ecm` | Y4 inherited verbatim | no new identity invariant introduced; v2 inherits v1's |

---

## 11. Cross-room dispatch instruction

This brief opens a design-discussion adversarial round. design-discussion-pane Claude (and design-discussion-pane Codex on adversarial review) should:

1. Read `docs/v2/v2_phase_e_v2_composition_brief.md` end-to-end (this file, after commit).
2. Either adopt the opening positions in §3 (Q1=a / Q2=a / Q3=b / Q4=a / Q5=b) OR open round 1 with substantive alternatives.
3. design-discussion-pane Codex provides adversarial review focused on:
   - **Anti-collapse invariant strength** (Q5 (b)): is the meta-test pattern adversarial-defensible? could a HB#1+#2 / HB#4-active bug fool the test?
   - **Result type choice** (Q1): does reuse-vs-new have meaningful downstream implications for Phase E v3 / Item 5 v2 sweep variant?
   - **Wording boundary graduation** (Q4): is "satisfaction" the right wording for an integrated composition that demonstrates closure but hasn't been validated against full FA dynamics step?
   - **Module location** (Q2): does adding v2 to `closed_loop_phase_e.py` violate sister-pattern with `closed_loop_phase_e_sweep.py` separation?
   - **Validation order** (Q3): cheap-parameter-failure-first vs early-fail-before-HB#3-scatter-O(N_FA-work)
4. Round-cap: 3-4 rounds; SEAL when no further design objections from Codex.
5. Cross-room MCP `topic = v2-phase-e-v2-composition` (NOT `v2-hard-blocker-*` since this is composition, not a hard blocker).

After SEAL: design-discussion-pane Claude commits lock artifact `docs/v2/v2_phase_e_v2_composition_step_2_locked.md` and dispatches to implementation-work-pane for Sanity Gate doc + code cycle.

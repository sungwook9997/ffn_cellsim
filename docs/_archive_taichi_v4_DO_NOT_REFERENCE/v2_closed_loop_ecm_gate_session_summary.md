# Closed-Loop ECM Gate — Session Summary & Navigation

**Status**: navigation-only summary, **non-locking**. Cross-references
the four primary closed-loop ECM gate documents and the Phase A/B/C
commit chain landed during the 2026-05-04 KST session. Future agents
or reviewers entering the closed-loop work should start here, then
follow the links below in order.

This document does not redefine any contract; it points at the
already-locked design docs. If a future design-discussion round
changes Phase contracts (Hard Blockers #1–#5), update the link
table below — do not duplicate Phase contracts here.

> ⚠️ **CURRENT STATUS OVERRIDE — 2026-05-06**
> (post-Option-A+B restructure audit cycle, claude-work
> `mcp_msg:2609`, codex approval `mcp_msg:2610`).
>
> Phase D no-op, Phase E v1, Phase E v2 step 2, HB#1/#2/#3/#4/#5,
> and Item 5 sweep have landed and were post-restructure audited.
> The audit cycle classified each unit as CLEAN or surfaced narrow
> drift fixes (banners only, no behavior changes); see audit chain
> `mcp_msg:2542 / 2544 / 2546 / 2548 / 2564 / 2571 / 2575 / 2578 /
> 2581 / 2585 / 2589 / 2607` for individual reports.
>
> Older "BLOCKED" wording below is retained as **2026-05-04 historical navigation context**
> and **must not be used as current state**. The
> per-Phase tables and commit chain in §§ Phase D / Phase E / "Commit
> chain (this session, 2026-05-04 closed-loop work)" stop at the
> 2026-05-04 KST session and have not been re-extended; refer to the
> child locks plus `docs/v2/v2_current_build_state.md` (commit
> `da581cd`, refreshed 2026-05-06) for current status.
>
> F-proper-1 remains WIP HALTED / coherent but untested; see
> `docs/v2/v2_phase_f_proper_1_fa_rate_response_locked.md` and
> commit `8a24d89`.
>
> Phase E v2 Item 5 sweep is a B-tier sister extension surfaced in
> `docs/v2/v2_item_5_sweep_harness_locked.md`; see commit `0f5e96f`.
>
> A future dedicated handoff-doc refresh unit may rewrite the §§
> Phase D / Phase E / commit-chain tables to current state; until
> then this banner is the load-bearing current-status pointer.

---

## Documents (read in this order)

| # | File | Role | Status |
|---|---|---|---|
| 1 | `docs/v2/v2_phase1_forward_roadmap.md` (lines 63–77) | Forward roadmap "Closed-Loop ECM Gate" §6.4 — the original 6-item gate | committed pre-session |
| 2 | `docs/v2/v2_closed_loop_ecm_gate_dependency_brief.md` | Per-item dependency analysis vs current ECM-OL/6.3b state, identifies 5 hard blockers, proposes phased entry | committed `923d0b2` |
| 3 | `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md` | **Locked phased plan** (3-round design-discussion, unresolved=0): Phase A-E definitions, Hard Blocker enumeration, effective_stiffness guard, constitutive direction seed, Lyapunov metric seed, Sanity Gate matrix, reference biology table, test catalog, cross-room dispatch checklist | committed `9c5e67e` |
| 4 | `docs/v2/v2_closed_loop_ecm_gate_phase_a_evidence.md` | Phase A non-locking evidence ledger: Items 3 + 6 open-loop side. Rule 10 unit-chain inline derivation included. | committed `af372ba` + `0ed202a` |
| 5 | `docs/v2/v2_ecm_ol_sweep_sanity_gate.md` | Phase C Sanity Gate: open-loop sweep harness contract (caller-supplied only, no harness signature defaults, runtime meta-test required) | committed `d0ae063` + `682c996` |

The 6.3b dynamics work (its own design-discussion lock + Sanity
Gate) is separate from the closed-loop gate but is referenced by
the locked phased plan §3 effective_stiffness guard:

| File | Role |
|---|---|
| `docs/v2/v2_63b_protrusion_coupling_design_brief.md` | 6.3b design-discussion brief |
| `docs/v2/v2_63b_protrusion_coupling_locked.md` | 6.3b locked design |
| `docs/v2/v2_63b_protrusion_coupled_fa_sanity_gate.md` | 6.3b Sanity Gate |

---

## Phase A — Items 3 + 6 open-loop evidence

**Status**: complete.

- Item 3 (no response when traction is zero): trivially satisfied
  on the open-loop side by `accumulate_prescribed_traction`'s
  bit-exact `x + 0 = x` algebra. Existing test
  `test_zero_traction_scenario_no_accumulation_change` in
  `tests/v2/test_v2_ecm_ol_harness.py`. **Closed-loop side blocked
  on Hard Blocker #1** (constitutive-law direction).
- Item 6 (Rule 10 unit-chain proof for traction-density storage):
  open-loop side documented in
  `acs/v2/dynamics/ecm_open_loop.py` module docstring §1 (per
  commit `22b9456`). **Closed-loop kPa↔nN/μm² comparison side
  blocked** until coupling functions land per Hard Blockers #3/#4.

Rule 10 inline derivation (correct factor): `1 kPa = 1 nN/μm²
exactly` (kPa = kN/m² = 10³ N/m²; nN/μm² = 10⁻⁹ N / (10⁻⁶ m)²
= 10³ N/m²). The Phase A evidence ledger originally carried
`1 kPa = 10⁻³ nN/μm²` (wrong by 3 orders); fixed in `0ed202a`
after Codex review id=1234. See memory
`rule10_unit_derivation_in_docs.md`.

---

## Phase B — Item 1 precursor (stimulus monotonicity, NOT satisfaction)

**Status**: complete.

Tests landed in `tests/v2/test_v2_ecm_open_loop.py` Phase B section:

1. `test_stimulus_accumulator_monotone_under_positive_traction`
2. `test_stimulus_accumulator_no_change_when_traction_zero`
3. `test_stimulus_accumulator_rejects_negative_traction_enforcing_monotonicity` (renamed from lock §8 draft `..._signed_traction_can_decrease` per Codex id=1241)
4. `test_stimulus_monotonicity_does_not_satisfy_closed_loop_item_1` (Hard Rule 11 wording-boundary meta-test)

Hard Rule 11 wording boundary enforced at three levels: section
header, per-test docstring, and the runtime meta-test (test 4).
**Phase B does NOT satisfy closed-loop ECM gate Item 1** —
measuring stimulus accumulator monotonicity is not measuring
remodeling response monotonicity. Item 1 satisfaction blocked on
Hard Blocker #1 (constitutive-law direction) + scaffolding.

Episode caught by Codex review id=1243: zero-traction test
docstring referenced "Item 3 closed-loop side" instead of Item 1
contract; fixed in `ada3728`. Captured in memory
`hard_rule_11_wording_boundary_meta_test.md`.

---

## Phase C — open-loop sweep baseline (NOT Item 5 satisfaction)

**Status**: complete.

- Harness: `scripts/run_ecm_ol_sensitivity_sweep.py` —
  caller-supplied tuple list × four ECM-OL preflight functions.
  No harness signature defaults; `_build_scenarios()` is an
  explicit fixture labeled `fixture_kind="non_production_smoke"`
  in `metadata.json` + `index.json`.
- Tests: `tests/v2/test_v2_ecm_ol_sweep.py` (20 tests after blocker
  fix). Test 5
  `test_open_loop_sweep_does_not_satisfy_closed_loop_item_5` is
  the Hard Rule 11 runtime meta-test (parallel to Phase B test 4).
- Visible deliverable: `runs/20260504T061746Z_ecm_ol_sensitivity/`
  (smoke fixture: 9 tuples × 4 channels = 36 records,
  aggregate_status=PASS, sweep_summary.png + per-tuple summary +
  per-channel diagnostic).

**Phase C does NOT satisfy closed-loop ECM gate Item 5** — the
sweep varies open-loop preflight outputs at varying (grid, dt),
which is not the closed-loop sensitivity sweep that Item 5 requires
(scattering geometry + response law).

Episodes caught:
- Codex id=1248: Sanity Gate §0 lock vs §8 default contradiction
  on the smoke scenario; fixed in `682c996` (strict
  caller-supplied, explicit fixture labeled non_production_smoke).
- Codex id=1253 (post-code): three blockers caught — sum reduction
  used max placeholder, channel scenarios accepted subset/duplicates,
  per-step diagnostic series promised in Sanity Gate but not
  implemented; fixed in `e2046db`.
- Codex id=1255 (non-blocking): HTML summary table missing
  `final_field_sum` / `selected_reduction_value` columns; surfaced
  in `8c1fae3`.

---

## Phase D — no-op closed-loop scaffolding (BLOCKED)

**Status**: blocked on **Hard Blockers #3 + #4**, awaiting
design-discussion rounds. impl-work is **idle / review-capable**
until interface locks land via cross-room dispatch.

Per locked phased plan §1 Phase D entry contract:
- Phase D code (`step_fa_to_ecm_response`, `step_ecm_to_fa_bias`
  function names) cannot land until #3 (FA→ECM scattering
  geometry interface) and #4 (ECM→FA bias target interface) are
  design-locked at the interface level — even if the default
  implementations are identity / neutral.
- Each new Phase D function gets its own Sanity Gate document.

**Hard Blocker #3 (FA→ECM scattering geometry interface)**:
1 FA at position `(x, y)` → which ECM grid cells receive its
traction signal? With what kernel? (Nearest-neighbor / bilinear /
Gaussian / other.) Interface must be locked even if the
default-no-op returns a zero field.

**Hard Blocker #4 (ECM→FA bias target interface)**: ECM state
biases FA dynamics how? (`traction_scale_nN` / rate multiplier
6.3b style / nucleation rate / combination.) Interface must be
locked even if default-no-op returns 1.0 multipliers (neutral).
Side-door risk: any read-only helper that consumes
`fiber_density` / `orientation_tensor` is ECM→FA bias by another
name (locked phased plan §3 effective_stiffness guard).

---

## Phase E — active closed-loop response (BLOCKED)

**Status**: blocked on all 5 Hard Blockers + effective_stiffness
law decision (if used). Phase E is the actual constitutive-law
implementation; everything before it is scaffolding or evidence.

Constitutive direction seed (locked phased plan §4, deferred to
Hard Blocker #1 round, **NOT a lock**):
- FA traction → accumulator/memory only (Item 3 trivial preserved)
- Remodeling response (Phase E):
  - `fiber_density`: bounded asymptotic toward 1 under traction
  - `orientation_tensor`: bounded asymptotic alignment toward
    scattering direction (after #3 lock)
  - raw `stiffness_kpa` unchanged in first closed-loop law

Reference biology seed (locked phased plan §7, literature must be
verified per Magic-Number Block before parameterization): Hall
2016 PNAS, Trichet 2012 PNAS, Stylianopoulos 2018 Nat Rev Cancer,
Eichinger 2020 Soft Matter, Notbohm 2015.

---

## Memory entries (this session)

- `rule10_unit_derivation_in_docs.md` — Hard Rule 10 inline
  unit-chain derivation applies to docs / evidence ledgers /
  briefs / locks, not just code. Triggered by `af372ba` `1 kPa =
  10⁻³ nN/μm²` bug (Codex id=1234).
- `hard_rule_11_wording_boundary_meta_test.md` — phased gate
  precursor/baseline evidence must include a runtime meta-test
  `test_<feature>_does_not_satisfy_<gate_item>` so a future
  refactor cannot silently relabel as satisfaction. Three concrete
  episodes captured (Phase B Item 1/3 docstring drift, Phase C
  default-vs-explicit fixture leak, Rule 10 docs unit error).

---

## Commit chain (this session, 2026-05-04 closed-loop work)

| Commit | Phase | Note |
|---|---|---|
| `923d0b2` | A pre | dependency brief (design-discussion routing) |
| `9c5e67e` | A pre | locked phased plan (3-round design-discussion) |
| `af372ba` | A | Phase A evidence ledger (initial; carried Rule 10 bug) |
| `0ed202a` | A | Phase A Rule 10 fix (1 kPa = 1 nN/μm²) |
| `66f06d6` | B | Phase B precursor tests (4 tests + Hard Rule 11 meta) |
| `ada3728` | B | Phase B docstring fix (Item 1 wording boundary) |
| `d0ae063` | C pre | Phase C Sanity Gate (pre-code) |
| `682c996` | C pre | Phase C Sanity Gate contradiction fix (strict caller-supplied) |
| `a74eb17` | C | Phase C code + tests (initial; sum/coverage/diagnostic gaps) |
| `e2046db` | C | Phase C blocker fixes (3 blockers from Codex review) |
| `8c1fae3` | C | Phase C HTML column UX fix (non-blocking) |

Visible deliverable: `runs/20260504T061746Z_ecm_ol_sensitivity/`
(open-loop sweep baseline, fixture_kind=non_production_smoke).

---

## What this summary is *not*

- Not a Sanity Gate. The Sanity Gate documents themselves
  (`v2_ecm_ol_sweep_sanity_gate.md` and the inherited
  `acs/v2/dynamics/ecm_open_loop.py` module docstring) are the
  gate of record.
- Not a design lock. The lock artifact
  (`v2_closed_loop_ecm_gate_phased_plan_locked.md`) is the
  source of truth for Phase A-E definitions and Hard Blockers
  #1-#5.
- Not an implementation plan for Phase D/E. Those phases are
  blocked on design rounds; the navigation links above point at
  the lock for definitions and at the Hard Blockers for what
  needs to be resolved before Phase D code lands.

---

## Update protocol

- When a Hard Blocker resolves (design-discussion round closes),
  add a row in the Phase D / Phase E section pointing at the
  resolution lock artifact.
- When a Phase D no-op scaffolding commit lands, add a row in the
  commit chain.
- When a Phase E active law lands, add a row + the satisfying
  test/evidence link.
- Do **not** restate Phase contracts here; always link back to the
  lock artifact.

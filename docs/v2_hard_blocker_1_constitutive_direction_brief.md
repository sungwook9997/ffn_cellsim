# V2 Phase 1 Hard Blocker #1 — Constitutive-Law Direction (Phase E entry) — Design-Discussion Brief (opening, NOT a lock)

**Date**: 2026-05-05 KST
**Status**: design-discussion brief, **opening position only —
NOT a lock, NOT a code commitment**.
**Author**: implementation-work Claude, drafted under (B-trigger)
precedent established by HB#3 / HB#4 / Phase D no-op scaffolding
/ B1 typed-schema migration (all locked + Sanity Gated +
implemented + Codex PASS in this session).
**Source**: locked phased plan
`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §2 Blocker
#1 + §3 effective-stiffness guard + §4 constitutive direction
seed + §7 reference biology table. Phase E sequencing audit
(commit pending) confirmed HB#1 is the next entry point per
locked phased plan §2.
**Codex routing**: `id=1465` ACK on Phase E sequencing
(`id=1464`); HB#1 brief next; Codex scope preferences
(`id=1465`):
1. Literature-first table scoped to direction/sign/boundedness/which-fields-update
2. Explicit allowed Phase E active-law candidates vs forbidden shortcuts
3. Decision fork: HB#1-only lock vs merged HB#1+#2 lock with selection criteria

**Hard contract**: this brief is the **opening position** for
the adversarial design round on the Phase E active constitutive
law direction. It does NOT lock the constitutive law shape, does
NOT commit to specific saturation forms, and does NOT authorize
Phase E code. The locked phased plan §2 Phase E remains BLOCKED
until HB#1 (this round) + HB#2 + HB#5 all locked + the
effective_stiffness law sub-decision resolved as part of HB#1.

---

## 0. Why this brief exists + Codex `id=1465` 3 scope preferences

Phase E active closed-loop ECM response cannot land until HB#1's
constitutive-law direction is locked. Per phased plan §2 Blocker
#1 Resolution: "separate adversarial design round with
literature-first table". Per phased plan §4: opening seed
positions (NOT a lock) are:

1. Stimulus memory ↔ remodeling response separation:
   `accumulated_traction_nNs_per_um2` = stimulus memory.
   Remodeling response = different field(s).
2. First closed-loop response targets collagen geometry:
   `fiber_density` and/or `orientation_tensor` before stiffness.
3. Raw `stiffness_kpa` unchanged in first closed-loop law.
4. Asymptotic-to-bound saturation for bounded fields.

This brief opens the design round on the four locked plan §4
seed positions plus the cross-cutting effective-stiffness sub-
decision and the decision fork on HB#2 merge.

**Codex `id=1465` 3 scope preferences enforced throughout the
brief**:

1. **Literature-first table** (§5) scoped to what determines
   direction/sign, boundedness, and which ECM fields update.
2. **Allowed vs forbidden separation** (§4): allowed Phase E
   active-law candidates explicit; no raw-stiffness shortcut,
   no empirical tuning, no PI-data fitting, no unbounded field
   growth.
3. **Decision fork on HB#2 merge** (Q4 in §2): HB#1-only lock
   vs merged HB#1+#2 lock, with criteria for choosing during
   the design round.

---

## 1. What is already locked (do not redebate here)

Per `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`:

- §1 Phase E scope: "active constitutive response, full Sanity
  Gate matrix, all 5 HBs + effective_stiffness law decision
  required".
- §2 Blocker #2 lock direction: "asymptotic-to-bound for any
  bounded schema field. Hard caps look like hidden clamping
  unless derived. Final form per-field decision in #2 round."
  (HB#2 framing locked even if HB#2 lock itself is pending).
- §3 effective-stiffness helper guard:
  - Raw `stiffness_kpa` retains raw meaning; no
    `effective_stiffness_kpa` stored field.
  - Phase A-C/D: NO `effective_stiffness()` helper unless it
    returns raw `stiffness_kpa` exactly (identity).
  - Once helper consumes `fiber_density` or `orientation_tensor`,
    it becomes ECM→FA mechanosensing law and MUST satisfy HB#4
    interface lock + Rule 10/sign/saturation gate.
- §4 constitutive direction seed (deferred to HB#1 round):
  see §0 above.
- §6 Sanity Gate matrix Phase E row: "full constitutive law /
  full saturation / full Lyapunov / full feedback dt / full
  sign-sense / full Item 1/2/4/5/6 closed-loop".
- §7 reference biology table (this brief's §5 expands).

Per `docs/v2_phase1_forward_roadmap.md` (`da40eef`):
- Closed-Loop Tier 1 single-cell milestone: 2 weeks stretch / 3
  weeks median / 4 weeks conservative.
- Closed-loop ECM gate has 6 items; failure of any item blocks
  closed-loop ECM, open-loop allowed.

Per Phase D lock `docs/v2_phase_d_no_op_scaffolding_locked.md`:
- Phase E lands as separate module
  (e.g., `closed_loop_phase_e.py`) with its own lock + Sanity
  Gate + code; NO Phase E variants in the Phase D module.

Per HB#3/#4 locks: scatter primitive + bias primitive interfaces
locked + impl PASS. Phase E response law CONSUMES these
primitives but does not modify them.

---

## 2. The four design questions

### Q1 — Which ECM field(s) update under FA traction stimulus?

The phased plan §4 seed proposes
`fiber_density` and/or `orientation_tensor` first;
`stiffness_kpa` unchanged. Round-1 should commit which subset.

Three candidate field sets:

- **(a) `fiber_density` only**: simplest; bounded `[0, 1]`;
  saturation form is asymptotic toward 1 under sustained
  traction (Hall 2016 supports this — stiffness perception is
  secondary).
- **(b) `orientation_tensor` only**: alignment toward
  scattering direction; bounded by tensor structure
  (eigenvalue magnitudes); saturation toward unit-vector
  alignment (Stylianopoulos 2018 supports — collagen alignment
  time-scale faster than stiffness change).
- **(c) Both `fiber_density` and `orientation_tensor`**: Hall
  2016 + Stylianopoulos 2018 jointly support both; Phase E law
  has two coupled responses; HB#5 Lyapunov metric needs to
  account for both.

`accumulated_traction_nNs_per_um2` = stimulus memory always
(Phase A-C open-loop preflight, locked). Not in HB#1 scope.

`stiffness_kpa` and `ligand_density` are NOT updated in HB#1
first law (per phased plan §3 + §4).

### Q2 — Direction/sign of update

For each field updated under traction stimulus, the direction
must be locked.

Three candidate sign conventions:

- **(a) Monotone-positive only**: traction
  `→ ∂t fiber_density ≥ 0` (fiber laid down under stress);
  traction `→ ∂t orientation_tensor` aligns toward stress
  direction (alignment magnitude ≥ 0). No relaxation in HB#1
  round; relaxation deferred to a later law.
- **(b) Bidirectional**: traction increases the response;
  absence of traction relaxes back toward baseline. HB#1 round
  must lock relaxation rate + baseline.
- **(c) Per-field separate sign convention**: `fiber_density`
  monotone-positive; `orientation_tensor` bidirectional
  (alignment toward stress vs return-to-isotropy).

Hall 2016 + Stylianopoulos 2018 support (a) for `fiber_density`
on Phase 1 timescales (relaxation slower than experimental
window). For `orientation_tensor`, Stylianopoulos 2018 supports
fast alignment on Phase 1 timescale; relaxation slower.

### Q3 — Boundedness / saturation form

Per phased plan §2 Blocker #2: "asymptotic-to-bound for any
bounded schema field. Hard caps look like hidden clamping unless
derived."

Three candidate saturation forms:

- **(a) Asymptotic-to-bound exponential**:
  `∂t field = α · (max_value - field) · stimulus`. Standard
  saturation kinetic; α > 0 has units `[stimulus]^-1 · [time]^-1`.
- **(b) Hill-function**:
  `∂t field = (V_max · stimulus^n / (K_m^n + stimulus^n)) · (max_value - field)`.
  Cooperativity exponent `n` is an extra parameter; literature
  rarely supports n > 1 for collagen at this scale.
- **(c) Hard cap**:
  `field += min(α · stimulus · dt, max_value - field)`. Looks
  like hidden clamping per phased plan §2 Blocker #2 lock
  direction. Default reject unless derived.

The phased plan §2 lock direction is (a) for bounded fields.
Round-1 should confirm + lock per-field choice.

### Q4 — Decision fork: HB#1-only lock vs merged HB#1+#2 lock

Per phased plan §2 Blocker #2 Resolution: "separate design
round, may merge with #1 round if constitutive direction
commits to specific bounded fields". Codex `id=1465` requested
this fork explicit in HB#1 round.

Two candidates:

- **(a) HB#1-only lock**: HB#1 commits to direction + sign +
  which fields update; HB#2 saturation form is a separate
  later lock with its own brief.
- **(b) Merged HB#1+#2 lock**: HB#1 commits to direction +
  sign + which fields update + per-field saturation form
  together. Single lock artifact
  `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md`.

**Selection criteria**:

- If HB#1 round 1 commits to **bounded fields only**
  (`fiber_density ∈ [0,1]`, `orientation_tensor` bounded by
  tensor structure), saturation form is a natural sub-decision
  and merge is **more efficient** (avoid round-trip latency).
- If HB#1 round 1 commits to **unbounded fields** OR if
  per-field saturation has independent design space (e.g.,
  Hill vs exponential trade-off differs per field), keep
  HB#2 separate.

The decision is **made during HB#1 round 1**, not pre-committed
in this brief.

---

## 3. Recommended opening positions (NOT a lock)

| Q | Recommendation | Rationale |
|---|---|---|
| Q1 fields | (c) Both `fiber_density` and `orientation_tensor`; `stiffness_kpa` unchanged | Hall 2016 + Stylianopoulos 2018 jointly support; phased plan §4 seed alignment |
| Q2 sign | (a) Monotone-positive only on Phase 1 timescale; relaxation deferred | Hall 2016 + Stylianopoulos 2018; relaxation timescale slower than Phase 1 experimental window |
| Q3 saturation | (a) Asymptotic-to-bound exponential per field | phased plan §2 Blocker #2 lock direction; α derived from literature `τ_align ≈ minutes` (Stylianopoulos 2018) |
| Q4 fork | **(b) merged HB#1+#2 lock** | Both Q1 candidates (a/b/c) commit to bounded fields; saturation form is natural sub-decision; avoids round-trip latency |

**Caveats / unresolved-by-this-brief**:

- The α coefficient in Q3 (a) — `[stimulus]^-1 · [time]^-1` —
  must be literature-extracted, not fit to PI experimental data
  (Hard Rule 1). Stylianopoulos 2018 + Eichinger 2020 are the
  candidate references; round-1 should commit to a specific
  parameter source.
- `orientation_tensor` saturation toward "scattering direction"
  needs an explicit definition: is it the unit vector of the
  per-cell scatter sum, or per-FA scatter direction? HB#3 lock
  outputs `(nx, ny, 2)` traction density; orientation alignment
  reads that vector field.
- Effective_stiffness sub-decision: if Phase E law only updates
  `fiber_density` + `orientation_tensor`, raw `stiffness_kpa`
  is not consumed by the law itself — but FA traction algebra
  in 6.3a uses `traction_scale_nN`, not stiffness. So
  effective_stiffness as ECM→FA mechanosensing path (HB#4
  surface) stays NEUTRAL in HB#1 first law. Sub-decision
  resolved as: **no effective_stiffness consumption in HB#1**.

---

## 4. Allowed Phase E active-law candidates vs Forbidden shortcuts

Per Codex `id=1465` scope preference 2.

### Allowed (Phase E active-law candidates)

- `∂t fiber_density = α_fiber · (1.0 - fiber_density) · stimulus`
  with `stimulus = |traction_density_xy|` (HB#3 magnitude per
  cell) and `α_fiber` derived from Stylianopoulos 2018 timescale.
- `∂t orientation_tensor = α_orient · (n_stim ⊗ n_stim - orientation_tensor)`
  with `n_stim = traction_density_xy / max(|traction_density_xy|, ε)`
  per cell and `α_orient` derived from Stylianopoulos 2018 fast
  alignment timescale. The asymptotic target is the rank-1
  outer product `n_stim ⊗ n_stim`, which is bounded
  (eigenvalues in `[0, 1]`).
- `accumulated_traction_nNs_per_um2` continues to accumulate
  per Phase A-C open-loop preflight (already locked, not in
  HB#1 scope).
- ECM→FA bias via HB#4 `compute_ecm_to_fa_bias_active`
  (separate function, separate Sanity Gate, separate lock).
  HB#1 does NOT specify the active bias; HB#1 specifies what
  the law DOES TO ECM, and HB#4-active picks up downstream.

### Forbidden (per phased plan + Hard Rules)

- **Raw `stiffness_kpa` mutation** (phased plan §3 guard).
- **Empirical tuning constants** (Magic-Number Block: must be
  derivable from literature, grid-invariant, NOT chosen to fit
  any test target).
- **PI-data fitting** (Hard Rule 1: literature-derived only,
  PI experimental data is comparison overlay only).
- **Unbounded field growth without saturation plan** (phased
  plan §2 Blocker #2 lock direction).
- **Hard caps that look like clamping** (phased plan §2
  Blocker #2 default reject).
- **Effective_stiffness shortcut except identity** (phased
  plan §3 guard; Phase E law's effective_stiffness consumption
  must satisfy HB#4 interface + Rule 10 + sign + saturation).
- **Scalarization of `traction_density_xy` at any layer**
  (Hard Rule 11; HB#3 lock; Phase D forbidden list). Phase E
  law owns the magnitude/direction reduction at THIS layer if
  needed (e.g., `|traction_density_xy|` for `fiber_density`
  rate, unit-vector for `orientation_tensor` target). The
  reduction is part of the locked law, not a hidden helper.

---

## 5. Reference biology table (literature-first per phased plan §7)

Scoped to **direction/sign**, **boundedness**, **which ECM
fields update**, and **timescale** (the latter constrains the α
coefficient in Q3 (a)).

| Reference | Year | Venue | Finding | HB#1 implication |
|---|---|---|---|---|
| Hall MS et al. | 2016 | PNAS | Cell traction induces collagen fiber **alignment**; stiffness perception is **secondary** to alignment-driven mechanotransduction | **Q1**: fiber alignment is a Phase E response; **Q2**: monotone-positive (traction → alignment); **Q3**: bounded by collagen geometry |
| Trichet L et al. | 2012 | PNAS | Substrate stiffness as cell **INPUT** (sensing), not output of cell traction | **Q1**: `stiffness_kpa` is INPUT-only, not updated by Phase E law; supports phased plan §3 guard |
| Stylianopoulos T et al. | 2018 | Nat Rev Cancer | Tumor ECM remodeling: **collagen alignment timescale faster than stiffness change**; alignment τ ≈ minutes, stiffness τ ≈ hours-days | **Q3 timescale**: α_orient ~ 1/minute; α_fiber slower; stiffness τ slow enough to ignore in first law |
| Eichinger J et al. | 2020 | Soft Matter | Computational fiber-network ECM remodeling reference simulation | **Q3 form**: confirms asymptotic-to-bound exponential is canonical in the literature; α derivable from τ |
| Notbohm J et al. | 2015 | Biophys J | Long-range force transmission in collagen via fiber alignment | **Q1**: fiber alignment is the long-range transmission channel; geometry-first response |

**Each reference must be verified before parameterization** (per
phased plan §7 hard rule). impl-work will read at least the
abstracts of Hall 2016 + Trichet 2012 + Stylianopoulos 2018
during the HB#1 round; PI is welcome to flag stronger references.

**Verification status pre-round-1**: all 5 entries inherited
verbatim from phased plan §7 reference biology table; no new
literature added by this brief. Codex `id=1465` round-1 should
challenge any reference that doesn't directly support the
direction/sign/boundedness/which-fields claim it's cited for.

---

## 6. What is *out of scope* for this brief

- HB#2 (saturation form) **per-field detail** — covered IF
  Q4 picks merged HB#1+#2; otherwise separate brief.
- HB#5 (Lyapunov metric) — strictly after HB#1+#2 lock per
  phased plan §2 Blocker #5 Resolution.
- HB#3 (FA→ECM scattering geometry) — already locked.
- HB#4 (ECM→FA bias target) — already locked. Active variant
  `compute_ecm_to_fa_bias_active` is a separate Phase E unit
  with its own lock.
- Effective_stiffness law as a **separate gate** — NOT a
  separate gate; it's a §3 guard already locked + a sub-
  decision inside HB#1 (resolved in §3 caveats above as
  "no effective_stiffness consumption in HB#1 first law").
- `accumulated_traction_nNs_per_um2` accumulation rate — locked
  in Phase A-C open-loop preflight.
- `stiffness_kpa` and `ligand_density` updates — phased plan §3
  + §4 forbid in first law.
- Phase E module file structure / Sanity Gate doc — written
  AFTER HB#1+possibly#2 lock.
- PI-data parameter fitting — Hard Rule 1.

---

## 7. Sanity Gate items the future HB#1 (or HB#1+#2) lock+code will need

Per CLAUDE.md Sanity Gate Protocol applied to a Phase E active
constitutive law:

- **§1 Dimensional**: every term in `∂t fiber_density` and
  `∂t orientation_tensor` reductions must be unit-checked. The
  α coefficient in Q3 (a) has units that depend on stimulus
  units (`[stimulus]^-1 · [time]^-1`); Rule 10 inline
  derivation required in lock + Sanity Gate doc.
- **§2 Boundary**: `fiber_density = 1.0` (saturated) and
  `fiber_density = 0.0` (depleted) cases must be tested;
  `traction_density_xy = 0` (zero stimulus → zero update)
  Item 3 closed-loop case.
- **§3 Conservation**: HB#1 law is **dissipative** (asymptotic
  to bound), not conservative. Energy-like Lyapunov metric is
  HB#5's job; HB#1 only states "monotone increase + asymptotic
  bound". Sign is positive (no negative `∂t fiber_density` in
  Phase 1 timescale per Q2 (a)).
- **§4 Numerical**: feedback dt — `dt_ecm · α · max(stimulus) ≤
  CFL` per phased plan Sanity Gate matrix Phase E row; explicit
  Euler stability bound derived from α.
- **§5 Sign**: monotone-positive `∂t fiber_density`;
  monotone-toward-target `∂t orientation_tensor`. Per-field
  sign documented + tested.
- **§6 Measurement-protocol**: scalarization decisions
  (`|traction_density_xy|` for `fiber_density` rate;
  `n_stim ⊗ n_stim` for `orientation_tensor` target) must be
  documented per Hard Rule 11; runtime meta-test
  `test_phase_e_law_does_not_update_stiffness_kpa` enforces
  the §3 guard at runtime.

**Magic-Number Block**: α_fiber + α_orient must be derivable
from literature (Stylianopoulos 2018 + Eichinger 2020),
grid-invariant in relative terms, and NOT chosen to fit any
test target. If literature τ values are not concrete enough,
HB#1 lock surfaces this as a blocker to PI before code lands.

---

## 8. What this brief is *not*

- Not a lock. The constitutive law direction is decided by the
  design-discussion adversarial round.
- Not a Sanity Gate doc. The Sanity Gate is written by the
  future Phase E code commit, after the lock artifact lands.
- Not a Phase E activation. Phase E code is a separate cycle
  with its own brief→lock→Sanity Gate→code unit.
- Not authoritative for biological parameters. α_fiber +
  α_orient are literature-extracted in HB#1 round, not fit to
  PI data.
- Not a closed-loop ECM gate satisfaction. HB#1 lock satisfies
  ONE blocker; closed-loop ECM still needs HB#2 (or merged
  HB#1+#2) + HB#5 + Phase E code commit.

---

## 9. Process expectation (mirrors HB#3 / HB#4 / Phase D / B1
precedents)

1. **Cross-room dispatch** to design-discussion with this brief
   + the 4 design questions + recommended opening positions.
   Status: `decision-needed`.
2. **Adversarial round** in design-discussion (Codex/Claude per
   memory `feedback_aggressive_design_debate.md`). HB#3 took 4
   rounds; HB#4 took 4 rounds; Phase D took 2 rounds; B1 took
   4 rounds. HB#1 is more substantive (literature-first +
   per-field decisions + merge fork) — **3-4 rounds plausible**.
3. **Lock artifact** delivered as either
   `docs/v2_hard_blocker_1_constitutive_direction_locked.md`
   (Q4=a HB#1-only) or
   `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md`
   (Q4=b merged), parallel to HB#3 / HB#4 / Phase D / B1 lock
   artifacts.
4. **Cross-room dispatch back to implementation-work** with the
   locked law shape.
5. **Sanity Gate doc** drafted by impl-work, Codex review.
6. **Code commit** — new module
   `acs/v2/dynamics/closed_loop_phase_e.py` (per Phase D lock
   §1 instruction) implementing the active law; new tests; new
   typed result/diagnostics dataclasses if needed; exports
   through both `__init__.py` surfaces.
7. **Codex review** of code commit.
8. **HB#1 (+ optionally HB#2) unit complete** = Phase E code
   landing-ready, pending HB#5 lock + final Phase E
   composition.

implementation-work stays idle/review-capable while this design
round runs.

---

## 10. References

- Locked phased plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §2
  Blocker #1 / §3 effective-stiffness guard / §4 constitutive
  direction seed / §6 Sanity Gate matrix Phase E / §7
  reference biology table.
- Forward roadmap: `docs/v2_phase1_forward_roadmap.md` (commit
  `da40eef`).
- Phase E sequencing audit: design-discussion `id=1464`
  (impl-work) + `id=1465` (Codex routing ACK).
- Sister-pattern lock precedents:
  - `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `docs/v2_phase_d_no_op_scaffolding_locked.md`
  - `docs/v2_focal_adhesion_dynamics_result_typed_locked.md`
- HB#3/#4 implementation:
  - `acs/v2/dynamics/fa_to_ecm_scattering.py`
  - `acs/v2/dynamics/ecm_to_fa_bias.py`
- Phase D implementation:
  - `acs/v2/dynamics/closed_loop_phase_d.py`
- ECM substrate schema: `acs/v2/ecm_substrate.py` (verified
  field set: `stiffness_kpa` ≥0, `ligand_density` ∈ [0,1],
  `fiber_density` ∈ [0,1], `orientation_tensor` (nx,ny,2,2),
  `accumulated_traction_nNs_per_um2` ≥0).
- Memory rules:
  - `design_note_pre_commit_batch.md` (6-step including Step 6
    sister-gate-mirror, per Codex `id=1428`).
  - `rule10_unit_derivation_in_docs.md`.
  - `hard_rule_11_wording_boundary_meta_test.md`.
  - `feedback_aggressive_design_debate.md`.
  - `cadence_promise_must_send_even_when_idle.md`.

---

## 11. Cross-room dispatch instruction (impl-work →
design-discussion)

The MCP message that accompanies this brief:
- `to`: `codex` (design-discussion-pane Codex receives + Claude
  pane reads)
- `room`: `design-discussion`
- `topic`: `v2-hard-blocker-1-constitutive-direction`
- `status`: `decision-needed`
- `body`: brief 4-question summary + opening recommendations +
  reference to this file
- `refs`: implementation-work `id=1462, 1463, 1464, 1465`,
  audit `da40eef`, locked phased plan §2/§3/§4/§7, HB#3 / HB#4
  / Phase D / B1 lock precedents

design-discussion round produces a lock artifact that supersedes
this brief; this file remains as historical opening-position
context.

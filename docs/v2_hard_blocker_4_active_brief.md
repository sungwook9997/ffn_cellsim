# V2 Closed-Loop ECM Gate — HB#4-Active Variant (Phase E v2 Step 1) — Design-Discussion Brief (opening, NOT a lock)

**Date**: 2026-05-05 KST
**Status**: design-discussion brief, **opening position only —
NOT a lock, NOT a code commitment**.
**Author**: implementation-work Claude, drafted under (B-trigger)
precedent established by HB#3 / HB#4 / Phase D no-op / B1 /
HB#1+#2 / HB#5 / Phase E v1 / Item 5 sweep harness (all locked
+ Sanity Gated + implemented + Codex PASS in this session, ~10
hours, 9 cycles).
**Source**: locked phased plan
`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §2 Hard
Blocker #4 (which committed Phase D-default neutral; active
variant deferred to v2 cycle) + Codex `id=1604`
post-Item-5-routing (a1) directive (HB#4-active brief next;
Phase E v2 composition step 2 deferred to future session).
**Unblock condition** (now satisfied): Phase E v1 composition +
Item 5 sweep harness both committed + impl PASS at `dec2336` +
`cfbb78d` (Codex `id=1600`). HB#4-active was explicitly
out-of-scope for HB#4 v1 (which locked at neutral) + Phase E v1
(which hard-wires HB#4 neutral); this brief opens the separate
HB#4-active cycle as Step 1 of Phase E v2.

**Hard contract**: this brief is the **opening position** for
the adversarial design round on the HB#4-active variant
(`compute_ecm_to_fa_bias_active`). It does NOT lock the active
law, does NOT commit the implementation, and does NOT authorize
Phase E v2 composition (Step 2 — separate later cycle). Phase E
v1's locked §0 wording boundary ("ECM-side composition evidence
ONLY, NOT full closed-loop ECM gate satisfaction") remains in
force until Phase E v2 composition step is itself locked.

---

## 0. Why this brief exists + Codex `id=1604` 4 review-focus
preferences

HB#4 v1 locked the **neutral** bias readout
(`compute_ecm_to_fa_bias_neutral`); it was the Phase D
no-op-default minimum interface. Phase E v1 hard-wired this
neutral bias (locked Y3: no v2 placeholder/hooks). HB#5 Y4
explicitly noted that "full FA→ECM→FA stability requires
HB#4-active design" (which is what this brief opens).

**HB#4-active v1 scope** (this brief's target):
- Pure read-only function `compute_ecm_to_fa_bias_active(adhesions,
  ecm)` — sister function name to `compute_ecm_to_fa_bias_neutral`
  with `_active` suffix per HB#4 v1 lock §0 silent-activation
  guard.
- Returns `ECMToFABiasResult` (same shape as neutral; HB#4 v1
  result type stays unchanged for forward compat); but
  `multipliers_per_fa` may now be non-1.0.
- Sampler `sample_ecm_at_fa_positions` shared with HB#4 v1;
  active variant does NOT introduce new sampling primitive.

**HB#4-active v2 / Phase E v2 composition scope** (out-of-scope
for THIS brief; Codex `id=1604` (a1) deferral):
- Composition function `step_closed_loop_phase_e_v2` that
  uses `compute_ecm_to_fa_bias_active` instead of neutral.
- Demonstrates **full FA→ECM→FA closure** (graduates from
  Phase E v1's ECM-side evidence to full closed-loop
  satisfaction).

**Codex `id=1604` 4 review-focus preferences enforced
throughout this brief**:

1. **HB#4-active wording boundary**: HB#4-active provides an
   **active ECM→FA bias law / evidence channel**, NOT a
   Phase E v2 satisfaction claim. Test names:
   `provides_active_bias_law_*` (NOT `satisfies_*`).
2. **`effective_stiffness` sub-decision explicit, NOT
   pre-committed**: per Codex `id=1604`, `effective_stiffness`
   is a sub-decision in this brief **only if the proposed
   active bias consumes `fiber_density` / `orientation_tensor`
   / `stiffness_kpa` for an effective-stiffness mechanosensing
   path**. The brief surfaces the question, NOT a default
   answer. Locked phased plan §3 effective-stiffness guard
   already constrains how this could land.
3. **Sister-pattern with HB#4 v1**: same return type
   (`ECMToFABiasResult`); same private validation chain;
   `compute_ecm_to_fa_bias_active` symbol naming per HB#4 v1
   §0 silent-activation guard.
4. **Sequencing locked**: HB#4-active step 1 lands first;
   Phase E v2 composition step 2 is a separate later cycle
   (Codex `id=1604` (a1) deferral).

---

## 1. What is already locked (do not redebate here)

Per locked phased plan + completed unit locks:

- **HB#4 v1 lock**: `compute_ecm_to_fa_bias_neutral` returns
  `ECMToFABiasResult` with all-1.0 multipliers + non-`None`
  empty `ECMSampledAtFAs` for empty FA list. Result type
  stays unchanged in HB#4-active variant for forward
  compatibility.
- **`sample_ecm_at_fa_positions`** sampler from HB#4 v1 lock:
  reads ECM fields at FA positions via bilinear interpolation;
  HB#4-active uses this sampler as-is (no new sampling primitive).
- **HB#4 v1 §0 silent-activation guard**: function naming
  separation `_neutral` / `_active`; no optional parameter
  switching.
- **Phase E v1 lock §0**: HB#4 hard-wired neutral; full
  FA→ECM→FA closure deferred to Phase E v2.
- **Phased plan §3 effective-stiffness guard**: Phase A-C/D
  forbids `effective_stiffness()` helper unless identity; Phase
  E (E v1 + E v2) is when active mechanosensing law could
  consume `fiber_density` / `orientation_tensor` IF the design
  commits to it (with HB#4 interface lock + Rule 10/sign/
  saturation gate).
- **HB#1+#2 lock §1 forbidden**: orientation-only first active
  law; `accumulated_traction_nNs_per_um2` NOT consumed (Y6).
  HB#4-active inherits this for any state-driven multiplier.
- **HB#5 lock + Y4 boundary**: "full FA→ECM→FA stability
  requires HB#4-active design (currently HB#4 = Phase D no-op
  neutral multipliers)" — this brief unlocks that path.

This brief is the next-level-down design choice **inside** that
already-locked frame.

---

## 2. The four design questions

### Q1 — Active bias mechanism: which ECM field(s) drive the
multiplier?

Three candidate mechanisms (per locked phased plan §2 Blocker
#4 listed options):

- **(a) Stiffness-driven**: `traction_scale_nN` per-FA scales
  with local sampled `stiffness_kpa`. Literature: Hall 2016 +
  Trichet 2012 (substrate stiffness as cell INPUT for
  mechanotransduction). **Triggers `effective_stiffness`
  sub-decision** (per Codex `id=1604` framing): if active bias
  reads `stiffness_kpa` directly, the §3 guard is satisfied
  trivially (raw stiffness as input is allowed; only mutating
  raw stiffness is forbidden). If active bias reads via an
  `effective_stiffness()` helper that consumes `fiber_density`
  / `orientation_tensor`, the §3 guard requires HB#4 interface
  lock + Rule 10/sign/saturation gate (which this lock would
  satisfy).
- **(b) Orientation-driven (alignment-with-traction)**:
  multiplier scales with FA-traction-direction alignment to
  local `orientation_tensor` principal axis. Reads
  HB#1+#2-updated `orientation_tensor` field. Sister-pattern
  with HB#1+#2 active law (which updates orientation under
  traction); this would be the closed-loop FA-rate response
  to the orientation update.
- **(c) Density-driven**: rates scale with local
  `fiber_density` (more fibers → more attachment opportunities
  → faster maturation). HB#1+#2 didn't update fiber_density
  in v1, so this would couple HB#4-active to a static fiber
  field (sampled but not evolving).

### Q2 — Active rate(s): which `RATE_NAMES` get modulated?

Per HB#4 v1 lock §1 `RATE_NAMES = ("k_maturity_per_s",
"k_bind_per_s", "k_unbind_per_s")`. Three candidate scopes:

- **(a) All 3 rates same multiplier**: simplest; one
  mechanism (Q1) drives all three rates uniformly via a single
  per-FA scalar.
- **(b) Per-rate independent multipliers**: maturity / binding
  / unbinding may respond differently to mechanical state
  (e.g., binding boosted under tension; unbinding suppressed).
  Three independent law forms.
- **(c) Hybrid (locked-canonical-pair)**: maturity-and-binding
  share one mechanism (positive feedback); unbinding has
  separate suppression mechanism. Two law forms.

### Q3 — Saturation form for active multipliers

Sister-pattern with HB#1+#2 Y10 (asymptotic-to-bound for
bounded fields). Active multipliers must be:

- **Bounded above**: cannot let `multipliers_per_fa[i, k]` grow
  unboundedly under high stimulus. Reasonable upper bound from
  literature: ~10× neutral (i.e., active multiplier ≤ 10.0).
- **Non-negative**: `multipliers_per_fa[i, k] ≥ 0` (sister with
  HB#5 V_active ≥ 0 sign discipline).

Three candidate saturation forms:

- **(a) Asymptotic-to-bound exponential**: mirrors HB#1+#2
  Y10. Convex-weight construction `m = m_min + (m_max - m_min)
  · (1 - exp(-K · S · dt))` for transient response or
  `m = m_min + (m_max - m_min) · S / (S + S_50)` for
  stimulus-driven steady state. Literature-derivable.
- **(b) Linear-with-clamp**: `m = max(m_min, min(m_max,
  m_neutral + α · S))`. Hard cap = silent clamp; rejected by
  HB#1+#2 Y10 lock direction.
- **(c) Hill-function**: `m = m_neutral · (1 + α · S^n / (S^n
  + S_50^n))`. Cooperativity exponent `n` extra parameter;
  literature rarely supports `n > 1` for FA mechanotransduction.

### Q4 — Effective_stiffness sub-decision (Codex `id=1604`
explicit framing)

If Q1 picks (a) stiffness-driven OR (c) density-driven via an
`effective_stiffness()` helper that consumes `fiber_density` /
`orientation_tensor`, the locked phased plan §3 effective-
stiffness guard activates: HB#4 interface lock + Rule
10/sign/saturation gate must be satisfied (which this lock
would do).

Three candidate framings (NOT pre-committed per Codex `id=1604`):

- **(a) `effective_stiffness` consumed**: Q1 picks (a) or (c)
  with helper; HB#4-active lock §1 explicitly defines an
  `effective_stiffness(ecm_at_fa)` private helper or inline
  expression that combines raw `stiffness_kpa` with
  `fiber_density` / `orientation_tensor` per a literature law.
  Locked plan §3 guard's "interface lock + Rule 10 + sign +
  saturation" is satisfied here.
- **(b) `effective_stiffness` deferred to a separate later
  unit**: Q1 picks (b) orientation-driven OR (a) stiffness-
  driven via raw `stiffness_kpa` ONLY (no helper). HB#4-active
  reads ECM fields directly without an effective-stiffness
  combination; future Phase F can introduce the helper if
  needed.
- **(c) `effective_stiffness` declined**: HB#4-active explicitly
  forbids any combination of fiber/orientation/stiffness via
  helper; the active law uses ONE field only. Sister-pattern
  with HB#1+#2 Y1 orientation-only first-law minimalism.

### Decision fork mapping (Q1 × Q4 dependency)

| Q1 mechanism | Q4 effective_stiffness implication |
|---|---|
| (a) Stiffness-driven (raw) | (b) deferred OR (c) declined |
| (a) Stiffness-driven (via helper) | (a) consumed; §3 guard activates |
| (b) Orientation-driven | (b) or (c); no §3 guard interaction |
| (c) Density-driven (raw) | (b) deferred OR (c) declined |
| (c) Density-driven (via helper) | (a) consumed; §3 guard activates |

---

## 3. Recommended opening positions (NOT a lock)

| Q | Recommendation | Rationale |
|---|---|---|
| Q1 mechanism | (b) Orientation-driven (alignment-with-traction) | Sister-pattern with HB#1+#2 active law (orientation-only); avoids §3 effective_stiffness sub-decision premature commitment; literature-defensible (Hall 2016 fiber alignment ↔ FA mechanotransduction loop); does NOT consume fiber_density which HB#1+#2 didn't update |
| Q2 active rate(s) | (a) All 3 rates same multiplier | Simplest; one mechanism drives all rates uniformly; Phase F can introduce per-rate refinement if needed |
| Q3 saturation form | (a) Asymptotic-to-bound exponential — sister with HB#1+#2 Y10 | Sister-pattern; literature-derivable; bounded above + non-negative by construction |
| Q4 effective_stiffness | (b) DEFERRED to separate later unit | Q1=(b) orientation-driven does NOT trigger §3 guard; effective_stiffness sub-decision can land in a future Phase F unit if/when needed |

**Caveats / unresolved-by-this-brief**:

- Q1=(b) orientation-driven multiplier needs explicit alignment
  formula. Candidate: `alignment = |n_traction · n_principal|`
  where `n_traction = traction_xy / |traction|` (FA's own
  traction direction) and `n_principal` is the principal-axis
  unit-vector of `orientation_tensor` at the FA's sampled
  cell. Multiplier = `m_max + (m_min - m_max) · alignment` for
  inverse coupling (low alignment → high multiplier, high
  alignment → low) OR `m_min + (m_max - m_min) · alignment`
  for direct coupling. Round 1 should commit direction.
- The α coefficient in Q3 (a) saturation must be literature-
  extracted, not fit to PI experimental data (Hard Rule 1).
  Hall 2016 + Trichet 2012 are candidate references.
- HB#4-active does NOT graduate Phase E v1 wording boundary;
  Phase E v2 composition (separate later cycle) is what
  graduates "ECM-side evidence" to "full closed-loop
  satisfaction".

---

## 4. Allowed HB#4-active candidates vs Forbidden shortcuts

### Allowed (HB#4-active v1 candidates)

- Pure read-only function
  `compute_ecm_to_fa_bias_active(adhesions, ecm) ->
  ECMToFABiasResult` (sister with HB#4 v1 neutral function).
- Same return type as HB#4 v1 (forward compat).
- Per-FA active multipliers in `multipliers_per_fa` matrix
  (shape `(N_FA, 3)` per `RATE_NAMES`).
- Reuses HB#4 v1 sampler `sample_ecm_at_fa_positions`.
- Bounded multipliers via Q3 (a) asymptotic-to-bound saturation.

### Forbidden (per Codex `id=1604` review focus + Hard Rules +
sister patterns)

- **PI-data fitting** (Hard Rule 1): all parameters
  literature-derived only.
- **Mutation of any ECM array field** (HB#4 v1 sister-pattern;
  pure read-only).
- **Mutation of any FA state** (read-only).
- **`effective_stiffness` consumption WITHOUT explicit Q4
  commit** (locked plan §3 guard): Q4 (b) or (c) means no
  `effective_stiffness` helper at all in HB#4-active v1.
- **Phase E v2 satisfaction wording**: HB#4-active is one
  step of Phase E v2; full satisfaction requires composition
  step 2 (separate later cycle). Test names
  `provides_active_bias_law_*_evidence` NOT `satisfies_*` —
  sister-pattern with Phase E v1 Y2 + Item 5 Y8.
- **Hard caps / clipping** that look like silent saturation
  (sister with HB#1+#2 Y10).
- **`accumulated_traction_nNs_per_um2` consumption** (sister
  with HB#1+#2 Y6).
- **`stiffness_kpa` mutation** (locked plan §3 guard).
- **Generic[...] parameterization or preemptive `# type:
  ignore`** (B1 sister-precedent).
- **`bias_params` optional parameter for switching v1/v2
  behavior** (Phase E v1 Y3 silent-activation guard).

---

## 5. Reference biology table (literature-first per phased
plan §7)

| Reference | Year | Venue | Finding | HB#4-active implication |
|---|---|---|---|---|
| Hall MS et al. | 2016 | PNAS | Cell traction induces collagen fiber alignment; FA mechanotransduction senses local alignment | **Q1**: orientation-driven (b) supported; alignment-with-traction is the natural FA-rate signal |
| Trichet L et al. | 2012 | PNAS | Substrate stiffness as cell INPUT (sensing), not output | **Q1**: stiffness-driven (a) is also defensible if Q4 (b) deferred (no helper); both mechanisms biologically real |
| Stylianopoulos T et al. | 2018 | Nat Rev Cancer | Tumor ECM remodeling: alignment time-scale faster than stiffness | **Q3 timescale**: HB#1+#2 K_ORIENT_PER_S (alignment ≈ 1-2h) sets the ECM-side timescale; HB#4-active multiplier timescale should match (FA rates respond fast) |

**Verification status pre-round-1**: 3 references inherited
from phased plan §7 reference biology table (used by HB#1+#2
lock); no new literature added by this brief. Codex `id=1604`
round-1 should challenge any reference that doesn't directly
support the chosen Q1 + Q3 + Q4 framing.

---

## 6. What is *out of scope* for this brief

- **Phase E v2 composition (step 2)** — separate later cycle;
  HB#4-active step 1 is THIS brief's scope only per Codex
  `id=1604` (a1) deferral.
- **`effective_stiffness()` helper introduction** — UNLESS Q4
  picks (a); per recommended Q4 (b), deferred.
- **Per-cell active law variations** (e.g., spatially-varying
  saturation parameter) — Phase F territory.
- **GPU implementation** — CPU-only per Phase E v1 sister-
  pattern.
- **PI-data parameter fitting** — Hard Rule 1.
- **HB#4-active diagnostics_dict surface change** — HB#4 v1
  grandfathered surface (per Phase D / B1 / id=1432); HB#4-
  active inherits the same `ECMToFABiasResult` shape.
- **Item 5 sweep harness extension to HB#4-active** — Item 5
  v1 sweeps Phase E v1 (HB#4 neutral); HB#4-active sweep is a
  separate later cycle.

---

## 7. Sanity Gate items the future HB#4-active lock+code will
need

Per CLAUDE.md Sanity Gate Protocol applied to a Phase E v2
step 1 active law:

- **§1 Dimensional**: every term in active multiplier formula
  must be unit-checked. Per Q1=(b) orientation-driven:
  `alignment` dimensionless; `multiplier` dimensionless;
  saturation argument `K · S · dt` dimensionless (sister with
  HB#1+#2). Rule 10 inline derivation required.
- **§2 Boundary**: empty FA list (HB#4 v1 sister); FA at
  zero-traction (multiplier = neutral 1.0); FA aligned
  perfectly (multiplier = m_max); all-isotropic ECM (no
  principal axis well-defined → multiplier = neutral 1.0
  fallback OR raise per round-1 lock).
- **§3 Conservation**: pure read-only; no per-step conservation
  of own; HB#4 v1 sampler conservation inherited.
- **§4 Numerical**: float64; expm1 stability if Q3 (a)
  exponential; no new tolerance.
- **§5 Sign**: multiplier ≥ 0 always; bounded above by m_max
  per Q3 (a) saturation.
- **§6 Measurement-protocol**: Hard Rule 11 — HB#4-active
  multiplier domain matches HB#4 v1 sampler domain (active
  cells where HB#3 scatter > 0 OR all FA positions
  unconditionally — round-1 must commit). Wording boundary:
  `provides_active_bias_law_*_evidence` test names; meta-test
  forward-guards against satisfaction-claim wording.

**Magic-Number Block**: `m_max` (saturation upper bound) +
`α` (saturation rate constant) must be derivable from
literature (Hall 2016 + Trichet 2012 candidate), grid-invariant
in relative terms, and NOT chosen to fit any test target.

---

## 8. What this brief is *not*

- Not a lock. The active law shape is decided by the
  design-discussion adversarial round.
- Not a Sanity Gate doc.
- Not a Phase E v2 composition (step 2 — separate later
  cycle).
- Not a closed-loop ECM gate satisfaction. v1 + Item 5 produce
  ECM-side evidence; v2 composition (after this HB#4-active
  step) provides full FA→ECM→FA evidence.
- Not authoritative for biological parameters. α, m_max are
  literature-extracted in HB#4-active round.
- Not an effective_stiffness law decision. Q4 (b) defers it;
  Q4 (a) would explicitly commit; round-1 chooses.

---

## 9. Process expectation (mirrors HB#3/#4/Phase D/B1/HB#1+#2/
HB#5/Phase E v1/Item 5 precedents)

1. **Cross-room dispatch** to design-discussion. Status:
   `decision-needed`.
2. **Adversarial round**. HB#3 4 / HB#4 v1 4 / Phase D 2 / B1
   4 / HB#1+#2 4 / HB#5 3 / Phase E v1 4 / Item 5 4 rounds.
   HB#4-active is more substantive (literature-first + per-
   field decision + saturation form + effective_stiffness
   sub-decision), so **3-4 rounds plausible**.
3. **Lock artifact** delivered as
   `docs/v2_hard_blocker_4_active_locked.md`.
4. **Cross-room dispatch back to implementation-work**.
5. **Sanity Gate doc** drafted by impl-work, Codex review.
6. **Code commit** — new function in
   `acs/v2/dynamics/ecm_to_fa_bias.py`
   (`compute_ecm_to_fa_bias_active`) with sister-validator;
   exports through both `__init__.py`; new test catalog with
   `provides_active_bias_law_*_evidence` test names.
7. **Codex review** of code commit.
8. **HB#4-active step 1 unit complete** = Phase E v2 step 1
   ACKed; Phase E v2 composition (step 2) is the next blocker
   (separate later cycle per Codex `id=1604` (a1) deferral).

implementation-work stays idle/review-capable while this design
round runs.

---

## 10. References

- Locked phased plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §2
  Hard Blocker #4 (HB#4 v1 + active variant) + §3 effective-
  stiffness guard + §7 reference biology table.
- Forward roadmap: `docs/v2_phase1_forward_roadmap.md` (commit
  `da40eef`).
- HB#4 v1 lock + impl (the sister neutral variant):
  `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` +
  `acs/v2/dynamics/ecm_to_fa_bias.py`.
- HB#1+#2 lock + impl (active law sister-pattern for
  saturation Q3):
  `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md` +
  `acs/v2/dynamics/ecm_constitutive_response.py`.
- HB#5 lock (Y4 wording boundary that this brief unlocks):
  `docs/v2_hard_blocker_5_lyapunov_metric_locked.md`.
- Phase E v1 lock (the composition that will USE this active
  variant in v2 step 2):
  `docs/v2_phase_e_composition_locked.md`.
- Item 5 sweep harness lock + impl (latest milestone in
  session):
  `docs/v2_item_5_sweep_harness_locked.md` +
  `acs/v2/dynamics/closed_loop_phase_e_sweep.py`.
- Codex post-Item-5 routing: `id=1600` (Item 5 PASS) + `id=1604`
  ((a1) HB#4-active step 1 next; effective_stiffness sub-
  decision explicit; HB#4-active wording = active bias law /
  evidence channel, NOT v2 satisfaction).
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
- `to`: `codex`
- `room`: `design-discussion`
- `topic`: `v2-hard-blocker-4-active`
- `status`: `decision-needed`
- `body`: brief 4-question summary + Codex `id=1604` 4
  review-focus preferences + recommended opening positions +
  Q1 × Q4 dependency table
- `refs`: implementation-work `id=1600, 1603, 1604` + locked
  phased plan §2 Blocker #4 + §3 + §7 + HB#4 v1 lock + Phase
  E v1 §0 wording boundary + HB#5 Y4

design-discussion round produces a lock artifact that
supersedes this brief; this file remains as historical
opening-position context.

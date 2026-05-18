# V2 Closed-Loop ECM Gate — Hard Blocker #5 Lyapunov-Like Metric — Design-Discussion Brief (opening, NOT a lock)

**Date**: 2026-05-05 KST
**Status**: design-discussion brief, **opening position only —
NOT a lock, NOT a code commitment**.
**Author**: implementation-work Claude, drafted under (B-trigger)
precedent established by HB#3 / HB#4 / Phase D no-op scaffolding /
B1 typed-schema migration / HB#1+#2 (all locked + Sanity Gated +
implemented + Codex PASS in this session).
**Source**: locked phased plan
`docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md` §2 Blocker
#5 + §5 Lyapunov metric seed.
**Unblock condition** (now satisfied): per phased plan §2
Blocker #5 Resolution: "MUST be after #1/#2/#3/#4 locked
(measurement-protocol-dependent on what loop actually updates)".
HB#1+#2 locked + impl PASS at `f05abef` (Codex `id=1495`) +
HB#3 ✓ + HB#4 ✓ → HB#5 entry now unblocked.

**Hard contract**: this brief is the **opening position** for
the adversarial design round on the Phase E Lyapunov-like
metric. It does NOT lock the metric form, does NOT commit the
implementation, and does NOT authorize Phase E composition. The
locked phased plan §1 Phase E remains BLOCKED on this lock +
closed-loop ECM gate composition.

---

## 0. Why this brief exists + HB#1+#2 lock outcome constraints

Phase E active closed-loop ECM response cannot land until HB#5's
Lyapunov-like metric is locked. The metric serves as the
closed-loop ECM gate Item 4 ("Bounded feedback in a single-cell
loop") evidence: a positive scalar function `V(state)` that
**decreases** along the closed-loop trajectory under bounded
input, demonstrating that the FA→ECM→FA feedback does NOT
saturate, oscillate, or hide arbitrary clipping (per phased
plan §6 Risks #2).

**HB#1+#2 lock outcome `f05abef` constrains HB#5's measurement
protocol**:

- **What the loop actually updates** (the §2 Blocker #5
  measurement-protocol-dependence): only `orientation_tensor`.
  `stiffness_kpa` / `fiber_density` / `ligand_density` /
  `accumulated_traction_nNs_per_um2` all unchanged in HB#1+#2
  first active law (Y1).
- **What drives the update**: instantaneous `traction_density_xy`
  (HB#3 output). `accumulated_traction_nNs_per_um2` is NOT
  consumed (Y6).
- **The update form**: per-cell convex combination
  `T_new = (1−w)·T_old + w·n⊗n` with `w = -expm1(-K·S·dt)`,
  where `S = |traction_density_xy| / TRACTION_REF`. Per-step
  distance-to-target `‖T_new − T_target‖_F` is monotone
  decreasing under FIXED target (per locked Y2).
- **Time-scales**: `τ_align` ~ 1-2 hours
  (`TAU_ALIGN_RANGE_S = (3600, 7200)`). Phase 1 timestep ~60 s
  → `K·S·dt ≈ 0.012` per step at unit stimulus.

These constraints **narrow** the phased plan §5 3 candidates:

- (1) Total accumulated traction energy: per-step input scales
  on `accumulated_traction` which is NOT consumed → candidate
  needs reformulation as **instantaneous traction energy**
  (`0.5 · Σ_grid |traction_density_xy|² · cell_area`).
- (2) Total fiber alignment magnitude: per phased plan §5 was
  `Σ_grid Tr(Q²)^0.5 · grid_volume` where `Q` is presumably the
  orientation tensor. Since orientation IS what HB#1+#2 updates,
  this candidate is **most directly aligned** with the loop's
  state.
- (3) Cell-substrate work rate: depends on FA velocity vs
  substrate velocity — Phase D wrapper has neither velocity
  field; HB#3 scatter is force-based, not work-based. This
  candidate **requires substantial new schema** to land in
  Phase E and is most distant from HB#1+#2 implementation.

---

## 1. What is already locked (do not redebate here)

- HB#1+#2 lock `7b1d3cf` (+ amendments `c18a5dc` + `42d34e7` +
  `7154342`) and code `f05abef`: orientation-only update,
  exact-exponential convex update, branch-defined zero-traction,
  two-sided validation contract, 5 owned `failure_kind`s, 7
  exports.
- HB#3 lock + code: scatter primitive with `traction_density_xy`
  output `(nx, ny, 2) nN/μm²`.
- HB#4 lock + code: bias primitive with neutral all-1.0 first
  law; sampler diagnostics for `orientation_tensor` already
  available.
- Phased plan §2 Blocker #5: "MUST be after #1/#2/#3/#4 locked
  (measurement-protocol-dependent on what loop actually
  updates)" — now satisfied.
- Phased plan §3 effective-stiffness guard: read-only helpers
  must satisfy HB#4 interface + Rule 10/sign/saturation gate.
  HB#5 metric definition that consumes `fiber_density` /
  `orientation_tensor` for an effective-stiffness shortcut
  would trip this guard.
- Phased plan §6 Sanity Gate matrix Phase E row: HB#5 falls
  under "full Lyapunov" — explicit Phase E expectation.

This brief is the next-level-down design choice **inside** that
already-locked frame.

---

## 2. The four design questions

### Q1 — Which scalar metric form?

Three candidate forms (phased plan §5 narrowed by HB#1+#2):

- **(a) Distance-to-target energy** (NEW; not in §5): `V(t) =
  0.5 · Σ_grid ‖orientation_tensor[i,j] − T_target[i,j]‖²_F ·
  cell_area`. Directly Lyapunov-style: per-cell
  distance-to-target Frobenius squared, integrated over grid.
  Per HB#1+#2 Y2 monotonicity, this **decreases per step under
  fixed target** by construction (the convex combination is
  contractive). Naturally satisfies positive-definite + decreasing.
- **(b) Instantaneous traction energy** (§5 candidate 1
  reformulated): `V(t) = 0.5 · Σ_grid |traction_density_xy|² ·
  cell_area`. Driven entirely by the input (FA traction); does
  NOT directly involve the ECM state. Decreases when traction
  decreases (input-driven), but says nothing about the
  closed-loop FA→ECM→FA stability. **Not Lyapunov in the loop
  sense**; better as an input-bound diagnostic.
- **(c) Orientation-magnitude** (§5 candidate 2): `V(t) =
  Σ_grid Tr(orientation_tensor²)^0.5 · cell_area`. Bounded
  above (since `|T_ij| ≤ 1` per schema, `Tr(T²) ≤ 4` per cell,
  so `V ≤ 2 · Σ_grid cell_area`). But **not necessarily
  monotone** under HB#1+#2 dynamics — the convex combination
  preserves the box but can move in either direction in the
  Frobenius magnitude depending on `T_old` vs `T_target`
  alignment.

### Q2 — Per-step decrease vs cumulative-bound metric?

The closed-loop ECM gate Item 4 wording ("Bounded feedback in a
single-cell loop") permits two interpretations:

- **(a) Strict per-step decrease**: `V(t+dt) ≤ V(t)` at every
  step under bounded input. The strongest Lyapunov-style
  guarantee.
- **(b) Cumulative bound**: `V(t) ≤ V_max` for all `t` under
  bounded input, but not necessarily monotone per step. Weaker
  but easier to prove for nonlinear closed-loop systems.

HB#1+#2 Y2 distance-to-target monotonicity is a **per-step**
decrease only **under fixed target** (i.e., constant traction
direction). Under varying traction direction (the actual
closed-loop case where FA traction depends on `orientation_tensor`
through the ECM→FA bias path), the target moves and the
per-step decrease guarantee is lost. So Q1 (a) plus Q2 (a) is
provably valid only for the artificial constant-input regime;
the closed-loop case needs Q2 (b) cumulative-bound framing.

### Q3 — Closed-loop coupling depth

Phase D no-op wrapper (committed) wires HB#3 + HB#4. Phase E
composition activates HB#1+#2 (orientation update) inside Phase
D's wrapper. The Lyapunov metric is computed **at what point**:

- **(a) After ECM update only**: metric depends on
  `orientation_tensor` (post-HB#1+#2). Does NOT include FA
  state. Tests "ECM relaxes under bounded traction".
- **(b) After ECM update AND ECM→FA bias readout**: metric
  depends on both ECM state and per-FA bias multipliers (HB#4
  output). Tests "ECM→FA channel bounded".
- **(c) Full closed-loop step**: metric depends on
  post-step FA state (after HB#4 bias has modulated the next
  6.3a/6.3b FA dynamics step, which then re-scatters). Tests
  "FA→ECM→FA→FA closed-loop bounded".

In Phase D no-op, HB#4 bias multipliers are all 1.0 — so (b)
and (c) collapse to (a). Phase E composition with HB#1+#2
active ECM update + HB#4 still-neutral bias means **(a) is the
only currently-testable option**; (b) and (c) require future
HB#4-active design.

### Q4 — Decision fork on candidate selection

Three routing options for the design-discussion round:

- **(a) Single metric**: lock one of Q1's (a/b/c) as THE
  Lyapunov metric. Cleanest contract; least flexibility.
- **(b) Composite metric**: lock a tuple `(V_distance, V_input,
  V_magnitude)` with each having its own role
  (Lyapunov-monotone + input-bound + state-magnitude). Tests
  separate properties separately. More API surface; clearer
  semantics.
- **(c) Pluggable strategy**: lock a metric-selection enum +
  separate dataclasses per strategy. Phase E composition picks
  which to compute. Most flexible; most API surface.

**HB#1+#2 sister-pattern**: HB#1+#2 chose **single function +
single result type** (NOT pluggable). HB#5 should mirror unless
explicit reason to diverge.

---

## 3. Recommended opening positions (NOT a lock)

| Q | Recommendation | Rationale |
|---|---|---|
| Q1 form | (a) distance-to-target energy `V = 0.5 · Σ ‖T - T_target‖²_F · cell_area` | Directly Lyapunov-style; HB#1+#2 Y2 contractive per fixed target; bounded above by `8 · grid_area` (each cell `‖T - T_target‖²_F ≤ 8` since `|T_ij|, |T_target_ij| ≤ 1`) |
| Q2 framing | **(b) cumulative bound** for closed-loop case | Q1 (a) is monotone under fixed target only; closed-loop has moving target via FA→ECM→FA path; cumulative-bound is the only honest framing |
| Q3 coupling | (a) after ECM update only (current Phase D no-op constraint) | HB#4 bias is still neutral; (b)/(c) require future HB#4-active design which is OUT of HB#5 scope |
| Q4 fork | (a) single metric | HB#1+#2 sister-pattern; minimum API surface; future Phase F can introduce composite if needed |

**Caveats / unresolved-by-this-brief**:

- Q1 (a) "distance-to-target" requires `T_target` to be defined.
  At measurement time, `T_target = n_stim ⊗ n_stim` per cell
  (HB#1+#2 algebra). For zero-traction cells, `T_target = 0`
  and the metric reduces to `0.5 · ‖T‖²_F` (orientation
  magnitude squared times cell area). Need to lock whether
  zero-traction cells contribute to the sum. My lean: yes,
  include all cells; metric is positive-definite over the full
  grid.
- The "bound" in Q2 (b) needs an explicit upper-limit constant.
  My lean: `V_max = 0.5 · 8 · grid_area = 4 · grid_area`
  (since per-cell `‖T - T_target‖²_F ≤ 8`). Derivable from
  schema invariants; satisfies Magic-Number Block.
- Q3 (a) means the metric does NOT close the FA loop — it only
  measures ECM relaxation under bounded scatter input. The
  "closed-loop" framing in Item 4 is therefore partially
  satisfied by HB#5 + (a); full closure requires HB#4-active
  + FA dynamics integration in a future Phase F.

---

## 4. Allowed Phase E metric candidates vs Forbidden shortcuts

### Allowed (Phase E metric candidates)

- **Pure-function compute**: `compute_ecm_orientation_metric(ecm,
  traction_density_xy) -> ECMOrientationMetricResult` — reads
  `orientation_tensor` + traction (to derive `T_target`),
  returns scalar `V` + diagnostics. Pure; no ECM mutation; no
  FA mutation; no RNG.
- **Typed result dataclass** `@dataclass(frozen=True,
  slots=True)` with explicit field set (sister-pattern with
  HB#1+#2 13-field diagnostics, Phase D 3-field result, etc.).
- **Diagnostics fields**: `n_nonzero_cells`, `n_total_cells`,
  `V_total`, `V_per_cell_max`, `V_per_cell_mean_nonzero`,
  `V_max_bound` (the locked upper bound), `traction_norm_max`,
  `cell_area_um2` (Rule 10 bridge).

### Forbidden (per phased plan + Hard Rules + sister patterns)

- **Mutation of any ECM state** (pure metric — read-only).
- **Consumption of `accumulated_traction_nNs_per_um2`** (HB#1+#2
  Y6 made this NOT-consumed; HB#5 reading it is an "effective
  memory shortcut" that sneaks the symbol back into the loop).
- **Empirical tuning constants** (Magic-Number Block).
- **PI-data fitting** (Hard Rule 1).
- **Effective_stiffness shortcut** (phased plan §3 guard).
- **Hidden hard caps / clipping** (anything that resembles
  bounding by clamping rather than by algebraic invariant).
- **Scalarization of `traction_density_xy` outside the locked
  reduction** (Hard Rule 11 — the metric's reduction must be
  documented as the locked measurement protocol).
- **Reading `fiber_density`** unless explicitly justified
  (HB#1+#2 didn't update it; reading it for a metric introduces
  dependence on a field the loop doesn't move).
- **Generic[...] parameterization or preemptive `# type: ignore`**
  (B1 sister-precedent).

---

## 5. Reference biology / control table

The Lyapunov-style guarantee in mechanobiology is rare in
literature; closest analogs are continuum-mechanics dissipation
arguments. Scoped narrowly:

| Reference | Year | Venue | Finding | HB#5 implication |
|---|---|---|---|---|
| Khalil (Nonlinear Systems) | 2002 | Prentice Hall, 3rd ed. | Lyapunov method definition: `V > 0`, `V̇ ≤ 0` ⇒ stable equilibrium | **Q1 + Q2**: lock the metric as positive-definite + non-increasing in expectation; per-step decrease is the strict version |
| Hall MS et al. | 2016 | PNAS | Fiber alignment under sustained traction reaches steady state on minutes timescale | **Q3**: `T → T_target` asymptote = bounded equilibrium under fixed input; HB#1+#2 Y2 reproduces this |
| Stylianopoulos T et al. | 2018 | Nat Rev Cancer | Tumor ECM remodeling has multiple timescales (alignment fast, density slow, stiffness slowest) | **scope**: HB#5 first metric covers alignment timescale only; density/stiffness Lyapunov is future |

**Verification status**: Khalil is the canonical Lyapunov
textbook reference; not domain-specific but the discipline
source. Hall + Stylianopoulos already verified for HB#1+#2.

---

## 6. What is *out of scope* for this brief

- HB#4-active variant (`compute_ecm_to_fa_bias_active`) — this
  is a separate Phase E unit not yet designed.
- FA dynamics changes in response to ECM (the "closing" of the
  closed loop) — needs HB#4-active first.
- `fiber_density` Lyapunov component — HB#1+#2 didn't update
  `fiber_density`; future Phase F constitutive law that does
  will surface a separate Lyapunov component.
- `stiffness_kpa` Lyapunov component — phased plan §3 guard
  forbids `stiffness_kpa` mutation in first closed-loop law.
- Cell-substrate work rate (phased plan §5 candidate 3) —
  requires velocity fields not currently in any v2 schema.
- PI-data parameter fitting — Hard Rule 1.
- HB#5 lock filename (will be
  `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md` per
  HB#1+#2 sister-pattern; not pre-locked here).

---

## 7. Sanity Gate items the future HB#5 lock+code will need

Per CLAUDE.md Sanity Gate Protocol applied to a Phase E
read-only metric:

- **§1 Dimensional**: every term in `V = 0.5 · Σ ‖T - T_target‖²_F
  · cell_area` must be unit-checked. `‖T‖_F` is dimensionless
  (since `T_ij` is dimensionless); `cell_area = dx² [μm²]`;
  `V` has units `[μm²]`. Rule 10 inline derivation required in
  lock + Sanity Gate doc.
- **§2 Boundary**: `T = T_target` (saturated) → `V = 0`;
  `T = 0` (depleted) and `traction = 0` → `T_target = 0`,
  `V = 0`; full grid + zero traction → `V = 0`; full grid +
  uniform max traction → `V` reaches `V_max = 4 · grid_area_um²`.
- **§3 Conservation**: pure function (no ECM/FA mutation); per-step
  decrease under HB#1+#2 fixed-target (provable from convex algebra);
  cumulative bound under varying target (Q2 (b) framing).
- **§4 Numerical**: float64 throughout; per-call work
  `O(grid_size)`; no new tolerance.
- **§5 Sign**: `V ≥ 0` always (Frobenius norm squared is
  non-negative); diagnostics (max, mean, etc.) all ≥ 0 by
  construction.
- **§6 Measurement-protocol**: scalarization happens at the
  metric's `Σ_grid` reduction — this is the locked reduction for
  Item 4. Runtime meta-test
  `test_lyapunov_metric_does_not_satisfy_full_closed_loop_item_4`
  enforces the boundary that HB#5 metric ALONE does NOT close
  the FA loop (only the ECM relaxation half).

**Magic-Number Block**: `V_max = 4 · grid_area_um²` derivable from
schema invariants (`|T_ij|, |T_target_ij| ≤ 1` ⇒ per-cell
`‖T - T_target‖²_F ≤ 8` ⇒ `V_per_cell ≤ 4` ⇒ `V_total ≤ 4 · n_cells
· cell_area = 4 · grid_area`); grid-invariant in absolute units;
NOT chosen to fit any test target.

---

## 8. What this brief is *not*

- Not a lock. The metric form is decided by the design-discussion
  adversarial round.
- Not a Sanity Gate doc.
- Not a Phase E activation. Phase E composition still requires
  HB#5 lock + HB#1+#2 impl already merged + Phase D no-op
  scaffolding ✓ + closed-loop ECM gate composition.
- Not a closed-loop ECM gate Item 4 satisfaction. The metric
  PROVIDES the evidence channel; whether Item 4 is satisfied
  depends on the simulation results in Phase E composition.
- Not authoritative for biological parameters. `V_max` is
  derivable from schema, not a literature parameter.

---

## 9. Process expectation (mirrors HB#3/#4/Phase D/B1/HB#1+#2
precedents)

1. **Cross-room dispatch** to design-discussion. Status:
   `decision-needed`.
2. **Adversarial round**. HB#3 4 / HB#4 4 / Phase D 2 / B1 4 /
   HB#1+#2 4 rounds. HB#5 has narrower scope (read-only metric,
   constrained by HB#1+#2 outcome), so **2-3 rounds plausible**.
3. **Lock artifact** delivered as
   `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md`.
4. **Cross-room dispatch back to implementation-work**.
5. **Sanity Gate doc** drafted by impl-work, Codex review.
6. **Code commit** — new module
   `acs/v2/dynamics/ecm_lyapunov_metric.py` (sister-pattern with
   `ecm_constitutive_response.py`) implementing the metric;
   typed result + diagnostics dataclasses; exports through both
   `__init__.py`.
7. **Codex review** of code commit.
8. **HB#5 unit complete** = all 5 Phase E hard blockers locked +
   impl PASS = Phase E composition becomes the next blocker
   (separate cycle).

implementation-work stays idle/review-capable while this design
round runs.

---

## 10. References

- Locked phased plan:
  `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md` §2
  Blocker #5 / §3 effective-stiffness guard / §5 Lyapunov seed
  / §6 Sanity Gate matrix Phase E.
- Forward roadmap: `docs/v2/v2_phase1_forward_roadmap.md` (commit
  `da40eef`).
- HB#1+#2 lock + code (the constitutive law HB#5 measures):
  - `docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md`
  - `docs/v2/v2_hard_blocker_1_2_constitutive_response_sanity_gate.md`
  - `acs/v2/dynamics/ecm_constitutive_response.py`
  - `tests/v2/test_v2_ecm_constitutive_response.py`
- Sister lock precedents:
  - `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `docs/v2/v2_phase_d_no_op_scaffolding_locked.md`
  - `docs/v2/v2_focal_adhesion_dynamics_result_typed_locked.md`
- ECM substrate schema: `acs/v2/ecm_substrate.py` (verified
  `orientation_tensor` shape `(*grid_shape, 2, 2)` + componentwise
  `|T_ij| ≤ 1` invariant).
- Memory rules:
  - `design_note_pre_commit_batch.md` (6-step including Step 6
    sister-gate-mirror, per Codex `id=1428`).
  - `rule10_unit_derivation_in_docs.md`.
  - `hard_rule_11_wording_boundary_meta_test.md`.
  - `feedback_aggressive_design_debate.md`.
  - `cadence_promise_must_send_even_when_idle.md`.
- Lyapunov reference: Khalil, *Nonlinear Systems*, 3rd ed.,
  Prentice Hall, 2002 (canonical textbook).

---

## 11. Cross-room dispatch instruction (impl-work →
design-discussion)

The MCP message that accompanies this brief:
- `to`: `codex`
- `room`: `design-discussion`
- `topic`: `v2-hard-blocker-5-lyapunov-metric`
- `status`: `decision-needed`
- `body`: brief 4-question summary + opening recommendations +
  HB#1+#2 outcome constraint summary
- `refs`: implementation-work `id=1495, 1492` (HB#1+#2 final
  PASS) + locked phased plan §2/§5/§6 + HB#1+#2 lock+code as
  measurement-protocol-determining.

design-discussion round produces a lock artifact that supersedes
this brief; this file remains as historical opening-position
context.

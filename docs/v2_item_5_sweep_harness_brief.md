# V2 Closed-Loop ECM Gate — Item 5 Sweep Harness — Design-Discussion Brief (opening, NOT a lock)

**Date**: 2026-05-05 KST
**Status**: design-discussion brief, **opening position only —
NOT a lock, NOT a code commitment**.
**Author**: implementation-work Claude, drafted under (B-trigger)
precedent established by HB#3 / HB#4 / Phase D no-op / B1 /
HB#1+#2 / HB#5 / Phase E v1 (all locked + Sanity Gated +
implemented + Codex PASS in this session).
**Source**: locked phased plan
`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §1
closed-loop ECM gate Item 5 + Codex `id=1563` post-Phase-E-v1
routing directive (Item 5 sweep harness first, independent of
Phase E v2 / effective_stiffness).
**Unblock condition** (now satisfied): Phase E v1 composition
locked + impl PASS at `ed99fb8` + `84b5890` (Codex `id=1559`).
Item 5 was explicitly out-of-scope for Phase E v1 (locked Phase
E v1 §6 + §0 v1-vs-v2 boundary); this brief opens the separate
Item 5 sweep harness cycle.

**Hard contract**: this brief is the **opening position** for
the adversarial design round on the Item 5 sensitivity sweep
harness. It does NOT lock the sweep ranges, does NOT commit the
implementation, and does NOT authorize any closed-loop ECM gate
satisfaction wording beyond Phase E v1's "ECM-side composition
evidence ONLY" framing per Codex `id=1563` review focus item.

---

## 0. Why this brief exists + Codex `id=1563` 5 review-focus
preferences

Closed-loop ECM gate Item 5 (per locked phased plan §1):

> 5. Sensitivity sweep over grid spacing and `dt_ecm`.

Phase E v1 explicitly out-of-scoped Item 5 (per locked Phase E
v1 §6 + parent plan §6 Sanity Gate matrix). Item 5 is a
**separate harness unit** that uses the existing Phase E v1
composition (`step_closed_loop_phase_e_v1`) across varied grid
configurations to demonstrate the gate's outputs are
grid-invariant in their physical interpretation.

**Codex `id=1563` 5 review-focus preferences enforced
throughout this brief**:

1. **No PI-data fitting**: Hard Rule 1 reaffirmation; sweep
   parameters are derived from numerical safety bounds (CFL,
   float64 round-off margins) + literature timescales (HB#1+#2
   `TAU_ALIGN_RANGE_S` from Hall 2016 + Notbohm 2015), NOT from
   PI experimental data.
2. **No gate tolerance changes to force pass**: per
   `feedback_aggressive_design_debate.md` + CLAUDE.md Hard
   Rule "NEVER modify a gate's tolerance, normalisation, or
   check window to make a failing run pass"; sweep harness
   reports raw measurements, does NOT adjust HB#5 metric bound
   or HB#1+#2 saturation thresholds.
3. **Explicit sweep axes/ranges + VRAM/runtime bounds**: locked
   phased plan §1 names "grid spacing and dt_ecm"; brief Q2
   commits to specific `(spacing_um, dt_s)` ranges + bounds
   per the Laptop A5000 16 GB baseline (CLAUDE.md Performance
   Budget).
4. **Deterministic/reproducible run metadata**: each sweep run
   captures (a) git commit hash, (b) full Phase E v1 input
   parameters + sweep coordinates, (c) HB#5 metric output +
   HB#1+#2 diagnostics + HB#3 / HB#4 output shape signatures.
   No RNG; all deterministic.
5. **ECM-side sensitivity evidence wording (NOT full
   bidirectional closure)**: Phase E v1 §0 wording boundary
   inherited; sweep harness produces **ECM-side sensitivity
   evidence** for Item 5; full bidirectional closure remains
   future work requiring Phase E v2.

---

## 1. What is already locked (do not redebate here)

Per locked phased plan + Phase E v1 lock:

- Phase E v1 composition `step_closed_loop_phase_e_v1` is the
  single per-step function the sweep harness drives.
- Phase E v1 §0 wording boundary: "ECM-side evidence ONLY";
  Item 5 sweep inherits this framing — it is **ECM-side
  sensitivity evidence**, NOT full bidirectional closure.
- HB#1+#2 constants `TRACTION_REF_NN_PER_UM2 = 1.0`,
  `K_ORIENT_PER_S ≈ 1.964e-4 1/s` are literature-pinned per
  HB#1+#2 lock §1; sweep harness uses these as defaults.
- HB#5 metric `compute_ecm_orientation_lyapunov_metric`
  returns `v_active_um2` + 10-field diagnostics; sweep harness
  reports these without modification.
- CLAUDE.md Performance Budget: Laptop A5000 16 GB baseline;
  pilot ≤ 2 GB; production ≤ 12 GB. Item 5 sweep is harness-
  level (not GPU); CPU-only via `acs.v2.dynamics`. Runtime
  budget: ≤ 10 min for full sweep on Laptop CPU.
- Hard Rule 1: NEVER fit parameters to PI's experimental data.
  Sweep parameters are literature-derived only.
- Hard Rule "NEVER modify a gate's tolerance, normalisation,
  or check window to make a failing run pass": sweep harness
  reports raw measurements; does NOT adjust HB#5 bound or
  HB#1+#2 thresholds.

This brief is the next-level-down design choice **inside** that
already-locked frame.

---

## 2. The four design questions

### Q1 — Sweep harness module + function shape

Three candidate shapes:

- **(a) Module `acs/v2/dynamics/closed_loop_phase_e_sweep.py`**
  with function `run_phase_e_v1_sensitivity_sweep(spacing_um_grid,
  dt_s_grid, n_steps, *, ecm_init, fa_init, ...) -> SweepResult`
  — single function takes 1-D arrays of `spacing_um` + `dt_s`
  + initial state, runs the cross-product, returns aggregated
  results.
- **(b) Module `acs/v2/analysis/phase_e_sweep_harness.py`** —
  separate `analysis/` namespace mirroring the
  `acs/v2/metrics.py` measurement-registry precedent. The
  function lives outside `dynamics/` because it's a
  measurement-on-dynamics layer, not new dynamics.
- **(c) Module `tools/closed_loop_phase_e_sweep.py` +
  command-line entry point** — separate harness directory; not
  part of the `acs.v2` package import surface.

Sister-pattern: HB#5 metric (read-only measurement) lives in
`acs/v2/dynamics/` per HB#5 lock Y9 "(5a)
`acs/v2/dynamics/ecm_lyapunov_metric.py` sister-pattern". Item
5 sweep harness is a measurement-on-Phase-E-v1; should it
follow the same convention (a), or move to a separate
`analysis/` namespace?

### Q2 — Sweep axes + ranges

Locked phased plan §1 names "grid spacing and `dt_ecm`". Two
sweep candidates:

- **(a) Cross-product**: `spacing_um ∈ {0.5, 1.0, 2.0, 4.0}`
  μm × `dt_s ∈ {30, 60, 120, 300}` s. 4 × 4 = 16 runs. Each
  run executes `n_steps = 100` Phase E v1 steps with constant
  rotating-traction (HB#5 Y9 boundedness scenario). Total runs
  ~16 × 100 steps × 4×4 ECM grid = 25,600 step calls. Runtime
  budget ≤ 10 min on Laptop CPU.
- **(b) 1-D sweeps each**: `spacing_um` sweep at fixed
  `dt_s = 60`, plus `dt_s` sweep at fixed `spacing_um = 1.0`.
  4 + 4 = 8 runs total. Smaller; faster; loses cross-axis
  interaction information.
- **(c) Adaptive grid**: start at coarse (4 × 4) cross-product,
  refine where gradients are large. Premature optimization
  for first sweep harness; reject.

**Sweep axes range rationale**:
- `spacing_um` ∈ {0.5, 1.0, 2.0, 4.0}: covers 0.5×–4× nominal
  resolution; catches under-/over-resolved cases.
- `dt_s` ∈ {30, 60, 120, 300}: at default `K_ORIENT_PER_S ≈
  2e-4`, `K · dt` ∈ {6e-3, 1.2e-2, 2.4e-2, 6e-2} — all well
  within HB#1+#2 numerical stability (no CFL bound; convex-
  weight construction unconditionally stable); upper end (300
  s = 5 min) approaches Phase 1 imaging interval.

### Q3 — Sweep-output measurement protocol

Per HB#5 metric, per sweep run we measure:
- HB#5 V_active trajectory (per step)
- HB#5 V_active_max_bound trajectory
- HB#5 max_traction_norm trajectory
- HB#1+#2 max_convex_weight trajectory
- HB#1+#2 max_orientation_delta_frobenius trajectory
- HB#3 traction_density_xy norm summary
- HB#4 multipliers_per_fa stays neutral check

Two candidate aggregation shapes:

- **(a) Full per-step trajectories**: each run returns
  `n_steps × n_diagnostic_fields` arrays. Total memory ~16
  runs × 100 steps × 10 fields × 8 bytes = ~128 KB
  (negligible). Allows post-hoc analysis of transients.
- **(b) Per-run summaries**: each run returns aggregated
  scalars (mean V_active, max V_active, final V_active, etc.).
  Smaller; loses transient dynamics.
- **(c) Both**: per-run summary + last-`m`-steps trajectory
  for transient inspection. Hybrid; more API surface.

**My lean (a)**: full trajectories — memory footprint is
trivial; loss of transients is irreversible; future Item 5
analysis (or comparison sweeps) benefits from full data.

### Q4 — Decision fork on item-5-evidence wording

Per Codex `id=1563` review focus #5: sweep harness produces
"ECM-side sensitivity evidence" — but the closed-loop ECM gate
Item 5 itself ("Sensitivity sweep over grid spacing and
`dt_ecm`") is one of the 6 items required for full closed-loop
ECM gate satisfaction. Two candidate framings:

- **(a) Item 5 evidence-only (NOT satisfaction)**: same
  discipline as Phase E v1 Y2 (`provides_*_evidence` test names).
  Sweep harness *produces evidence* for Item 5; the **PI** (or
  the explicit closed-loop ECM gate satisfaction discipline)
  decides whether the evidence satisfies Item 5 OR if Phase E
  v2 is also required for Item 5.
- **(b) Item 5 satisfaction direct**: claim sweep harness
  satisfies Item 5 because the locked plan §1 wording is
  "Sensitivity sweep over grid spacing and `dt_ecm`" — purely
  procedural, no v1-vs-v2 distinction. Risk: silent overclaim
  if interpreted as full closed-loop ECM gate satisfaction
  (which requires all 6 items + HB#5 active feedback).

My lean **(a)**: Item 5 evidence-only wording. Mirrors Phase E
v1 Y2 pattern. Allows future Phase E v2 sweep variant to
contribute additional evidence without overclaim ambiguity.

---

## 3. Recommended opening positions (NOT a lock)

| Q | Recommendation | Rationale |
|---|---|---|
| Q1 module/function | (a) `acs/v2/dynamics/closed_loop_phase_e_sweep.py` + `run_phase_e_v1_sensitivity_sweep` | Sister-pattern with HB#5 metric (read-only measurement on dynamics in `acs/v2/dynamics/`); same precedent (HB#5 Y9 (5a)) |
| Q2 sweep axes | (a) Cross-product 4 × 4 = 16 runs; ranges spacing_um ∈ {0.5, 1.0, 2.0, 4.0} × dt_s ∈ {30, 60, 120, 300} | Covers physical regimes per Codex `id=1563` review focus #3; budget ≤ 10 min CPU |
| Q3 output shape | (a) Full per-step trajectories | Memory trivial (~128 KB total); loss of transient is irreversible |
| Q4 wording | (a) Item 5 evidence-only (NOT satisfaction); use `provides_item_5_*_evidence` test names | Mirror Phase E v1 Y2 pattern; avoid silent overclaim |

**Caveats / unresolved-by-this-brief**:

- The sweep harness output shape (`SweepResult` dataclass) needs
  explicit definition: typed dataclass per B1/Phase D/HB#5/Phase
  E v1 sister-pattern, with explicit fields for each axis ×
  step × diagnostic. Round 1 should commit shape.
- "Deterministic/reproducible run metadata" per Codex `id=1563`
  review focus #4: each sweep run records git commit hash +
  full input parameters + Phase E v1 / HB#5 diagnostics.
  Implementation detail: include a `SweepRunMetadata` dataclass
  with `git_commit_sha: str`, `input_parameters: dict`,
  `phase_e_v1_constants: dict`, `numpy_version: str`, etc.
- "VRAM/runtime bounds" per Codex `id=1563` review focus #3:
  CPU-only harness (no GPU); VRAM bound N/A. Runtime budget
  is the constraint. Round 1 should lock the budget number.

---

## 4. Allowed Item 5 sweep harness candidates vs Forbidden
shortcuts

### Allowed (Item 5 v1 candidates)

- Pure-function `run_phase_e_v1_sensitivity_sweep(spacing_um_grid,
  dt_s_grid, n_steps, *, ecm_init, fa_init, ...) -> SweepResult`.
- Typed `SweepResult` + `SweepRunMetadata` dataclasses (frozen+
  slots per sister-pattern).
- Cross-product evaluation of `spacing_um × dt_s`.
- Per-run full per-step trajectories (16 fields × 100 steps ×
  16 runs = trivial memory).
- Reproducible run metadata: git SHA + numpy version +
  `K_ORIENT_PER_S` value + `TRACTION_REF_NN_PER_UM2` value +
  fa_init + ecm_init.
- Test catalog: ~10 tests covering harness invariants
  (deterministic across calls, monotone-in-axes evidence,
  HB#1+#2 + HB#5 + HB#3 + HB#4 outputs all respected per axis).

### Forbidden (per Codex `id=1563` review focus + Hard Rules +
sister patterns)

- **PI-data fitting** (Codex `id=1563` review focus #1): sweep
  parameters NOT chosen to match PI experimental Bare/Pre/Lam4
  data (`data/experimental/`).
- **Gate tolerance changes** (Codex `id=1563` review focus #2):
  sweep harness does NOT modify HB#5 V_active_max_bound or
  HB#1+#2 thresholds; reports raw measurements only.
- **Item 5 satisfaction overclaim** (Codex `id=1563` review
  focus #5 + Phase E v1 Y2): test names `provides_item_5_*_evidence`
  (NOT `satisfies_item_5_*`); meta-test forward-guards against
  satisfaction-claim wording in module source.
- **Full bidirectional closure claim** (Phase E v1 §0 boundary):
  sweep harness uses Phase E v1 (HB#4 neutral); cannot
  demonstrate full FA→ECM→FA closure.
- **Phase E v2 hooks / placeholders** (sister Y3 pattern from
  Phase E v1): no `_v2` placeholder; no optional `bias_law`
  parameter.
- **RNG / nondeterminism**: sweep is fully deterministic; no
  random initialization (or, if random init is needed, locked
  seed with explicit seed parameter).
- **GPU/Taichi calls** (Phase E v1 + sister modules are
  CPU-only via `acs.v2.dynamics`): sweep harness does NOT
  invoke GPU at this stage.
- **Generic[...] parameterization or preemptive `# type: ignore`**
  (B1 sister-precedent).
- **CSV output with PI-data column names**: per Hard Rule 1,
  sweep output is internal HDF5 / numpy / json; if CSV is
  needed, column names are sweep-specific (e.g.,
  `spacing_um`, `dt_s`, `v_active_um2_at_step_100`), NOT
  PI-data column matches.

---

## 5. Reference biology / control / numerics table

The Item 5 sweep is mostly numerical (grid + dt sensitivity);
biology is embedded in HB#1+#2 (literature-pinned). The sweep
references:

| Reference | Usage |
|---|---|
| HB#1+#2 lock §1 (Munevar 2001 + Tan 2003 for `TRACTION_REF`; Hall 2016 + Notbohm 2015 for `TAU_ALIGN`) | Sweep keeps these constants fixed; varies only `spacing_um` + `dt_s` |
| HB#5 lock §1 (`MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL = 4.5`) | Sweep verifies V_active stays ≤ 4.5 · active_area across all axes |
| CLAUDE.md Performance Budget (Laptop A5000 16 GB baseline; CPU acceptable for harness) | Runtime budget ≤ 10 min CPU |
| Khalil *Nonlinear Systems* (HB#5 reference) | Lyapunov-style boundedness across grid → "ECM-side sensitivity evidence" framing |

No new biology references; this is a numerical harness.

---

## 6. What is *out of scope* for this brief

- Phase E v2 (HB#4-active variant) — separate cycle; sweep
  harness uses Phase E v1 only.
- `effective_stiffness` law decision — separate cycle.
- Full FA→ECM→FA closed-loop satisfaction — Phase E v2
  territory.
- Closed-loop ECM gate Item 1-4 evidence — Phase E v1 owned.
- Closed-loop ECM gate Item 6 (Rule 10 inline) — already proven
  in each HB lock.
- Item 5 satisfaction wording — stays evidence-only per Q4 (a).
- GPU sweep performance evaluation — separate later cycle (when
  GPU integration is ready).
- PI-data parameter fitting — Hard Rule 1.
- Comparison against PI experimental Bare/Pre/Lam4 spreading
  data — comparison overlay only, NOT fitting.

---

## 7. Sanity Gate items the future Item 5 v1 lock+code will
need

Per CLAUDE.md Sanity Gate Protocol applied to a sweep harness:

- **§1 Dimensional**: harness performs no unit transformation;
  inherits Phase E v1 + HB#5 unit chains. Sweep axis units:
  `spacing_um [μm]`, `dt_s [s]`.
- **§2 Boundary**: empty axis arrays (single-cell sweep)
  reduce to single Phase E v1 call; max axis values respect
  Phase E v1 boundary cases (dt=0 valid no-op, dt<0 raises).
- **§3 Conservation**: harness conserves nothing of its own;
  inherits Phase E v1 + sub-call invariants.
- **§4 Numerical**: sweep result memory bounded ~128 KB;
  runtime ≤ 10 min CPU; deterministic across re-runs (locked
  seed if any RNG; no RNG in current scope).
- **§5 Sign**: V_active ≥ 0 across all sweep axes.
- **§6 Measurement-protocol**: Item 5 evidence-only wording
  (NOT satisfaction); deterministic reproducible metadata
  (git SHA + parameters); raw-measurement reporting (no
  tolerance adjustment).

**Magic-Number Block**: zero new tunables at harness layer
(sweep ranges are literature-derived per Q2 (a) rationale; no
fitting); `runtime_budget_minutes = 10` is a hardware-derived
constraint per CLAUDE.md Performance Budget, not fitted.

---

## 8. What this brief is *not*

- Not a lock. Sweep harness shape + ranges + wording boundary
  are decided by the design-discussion adversarial round.
- Not a Sanity Gate doc.
- Not a Phase E v2 / HB#4-active activation.
- Not a closed-loop ECM gate Item 5 satisfaction (per Q4 (a)
  evidence-only framing).
- Not authoritative for biological parameters (none new).
- Not a GPU sweep performance evaluation.

---

## 9. Process expectation (mirrors HB#3/#4/Phase D/B1/HB#1+#2/
HB#5/Phase E v1 precedents)

1. **Cross-room dispatch** to design-discussion. Status:
   `decision-needed`.
2. **Adversarial round**. Item 5 sweep is harness-level (no
   new physics); narrower than Phase E v2. **2-3 rounds
   plausible**.
3. **Lock artifact** delivered as
   `docs/v2_item_5_sweep_harness_locked.md`.
4. **Cross-room dispatch back to implementation-work**.
5. **Sanity Gate doc** drafted by impl-work, Codex review.
6. **Code commit** — new module
   `acs/v2/dynamics/closed_loop_phase_e_sweep.py` (sister-
   pattern with `closed_loop_phase_e.py`) + typed
   `SweepResult` + `SweepRunMetadata` dataclasses + sweep
   function + module-docstring forbidden-list guard; new test
   catalog including `provides_item_5_*_evidence` tests; exports
   through both `__init__.py`.
7. **Codex review** of code commit.
8. **Item 5 v1 unit complete** = closed-loop ECM gate Item 5
   *evidence* available; whether evidence satisfies Item 5
   gate-pass requires PI / explicit-decision separately.

implementation-work stays idle/review-capable while this design
round runs.

---

## 10. References

- Locked phased plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §1
  closed-loop ECM gate Item 5 + §6 Sanity Gate matrix.
- Forward roadmap: `docs/v2_phase1_forward_roadmap.md` (commit
  `da40eef`).
- Phase E v1 lock + impl (the sweep harness's target):
  - `docs/v2_phase_e_composition_locked.md` (commit `c8b9550`)
  - `docs/v2_phase_e_composition_sanity_gate.md` (commits
    `b0b5baa` + `a3bc8a6` + `e5a8186`)
  - `acs/v2/dynamics/closed_loop_phase_e.py` (commits
    `ed99fb8` + `84b5890`)
- HB#5 lock + impl (sister-pattern for sweep wording boundary
  + Y9 module-path precedent):
  - `docs/v2_hard_blocker_5_lyapunov_metric_locked.md`
  - `acs/v2/dynamics/ecm_lyapunov_metric.py`
- Sister lock precedents:
  - `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `docs/v2_phase_d_no_op_scaffolding_locked.md`
  - `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md`
  - `docs/v2_focal_adhesion_dynamics_result_typed_locked.md` (B1)
- ECM substrate schema: `acs/v2/ecm_substrate.py`.
- v2 measurement registry (separate concern, NOT this unit's
  home): `acs/v2/metrics.py` (file, not directory).
- CLAUDE.md Performance Budget (Laptop A5000 16 GB baseline +
  CPU-only harness acceptable).
- Hard Rule 1: NEVER fit parameters to PI's experimental data.
- Hard Rule "NEVER modify a gate's tolerance ... to make a
  failing run pass".
- Memory rules:
  - `design_note_pre_commit_batch.md` (6-step including Step 6
    sister-gate-mirror, per Codex `id=1428`).
  - `rule10_unit_derivation_in_docs.md`.
  - `hard_rule_11_wording_boundary_meta_test.md`.
  - `feedback_aggressive_design_debate.md`.
  - `cadence_promise_must_send_even_when_idle.md`.
- Codex post-Phase-E-v1 routing: `id=1559` (Phase E v1 PASS) +
  `id=1563` (Item 5 sweep harness first; 5 review-focus
  preferences).

---

## 11. Cross-room dispatch instruction (impl-work →
design-discussion)

The MCP message that accompanies this brief:
- `to`: `codex`
- `room`: `design-discussion`
- `topic`: `v2-item-5-sweep-harness`
- `status`: `decision-needed`
- `body`: brief 4-question summary + Codex `id=1563` 5
  review-focus preferences + recommended opening positions
- `refs`: implementation-work `id=1559` (Phase E v1 PASS) +
  `id=1562, 1563` (routing decision) + locked phased plan §1
  Item 5 + Phase E v1 lock as the sweep target

design-discussion round produces a lock artifact that
supersedes this brief; this file remains as historical
opening-position context.

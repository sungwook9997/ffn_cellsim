# V2 Phase 1 — Open-Loop ECM Sensitivity Sweep Sanity Gate (Phase C)

**Status**: pre-execution Sanity Gate per Plan §10/§11 + closed-loop
ECM gate phased plan
(`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §1 Phase C).
Written from the locked phase plan (3-round design-discussion lock,
MCP id=1223–1230 + Codex impl ack id=1245). Pre-code per Codex impl
sequencing id=1239 (Sanity Gate first → review → code/tests → run).

**Contract**: this gate covers the **open-loop sensitivity sweep
harness** that varies `(grid_spacing, dt_s)` over the four ECM-OL
preflight functions and persists per-(grid, dt) baseline artifacts.
It is **explicitly Phase C baseline**, NOT closed-loop Item 5
satisfaction. Hard Rule 11 wording protection: any artifact emitted
must say "open-loop sweep baseline" and never "Item 5 satisfied".

**Source of truth**:
- locked phased plan: `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
  §1 Phase C, §6 sanity-gate matrix Phase C row, §8 Phase C test
  catalog, §9 Phase C immediate go-list
- ECM-OL preflight gate (already proven, inherited):
  `acs/v2/dynamics/ecm_open_loop.py` module docstring (post `22b9456`)
- ECM-OL harness pattern (inherited):
  `acs/v2/ecm_open_loop_harness.py` (`04ee5a7`) +
  `scripts/run_ecm_ol_harness.py` (`9b3dac5` / `f8cdff3`)

**Sequencing context**: Phase B precursor tests landed (`66f06d6`,
`ada3728`); Phase C is the next "Separated dynamics modes" sub-unit
before Phase D (closed-loop scaffolding) — Phase D is gated on Hard
Blockers #3 (FA→ECM scattering geometry) and #4 (ECM→FA bias
target) interface lock. Phase C is independent of all 5 hard
blockers.

---

## 0. Scope (locked, no scope creep allowed)

### In scope (this unit only)

- A sweep harness over caller-supplied lists of
  `(grid_n, spacing_um, dt_s)` tuples.
- Per-tuple, run a fixed prescribed-input scenario for each of the
  four ECM-OL preflight functions
  (`accumulate_prescribed_traction`,
  `apply_prescribed_stiffness_rate`,
  `apply_prescribed_density_rate`,
  `apply_prescribed_orientation_rate`).
- Per-tuple per-channel summary metrics: final accumulator/field
  max/min, max effective rate (where applicable), per-step
  diagnostic series.
- Per-(grid, dt) `metadata.json` + the cross-tuple `index.json`
  with `aggregate_status`, `intentional_failures_observed`,
  `unexpected_failures`, `unexpected_passes` (same shape as the
  6.3b runner under `scripts/run_protrusion_coupled_fa_harness.py`).
- A cross-tuple summary plot per channel showing how the chosen
  baseline metric varies with `(grid_n, dt_s)`.
- Persistent artifacts under `runs/<UTC>_ecm_ol_sensitivity/`
  (seconds-resolution timestamp, `os.makedirs(exist_ok=False)`,
  same f8cdff3 hygiene).

### Explicitly out of scope (will FAIL gate if introduced)

- Closed-loop ECM activation. NO FA→ECM, ECM→FA, response-field
  update, ECM bias on FA — the sweep runs the four open-loop
  preflight functions over varying (grid, dt), period.
- "Item 5 closed-loop sensitivity satisfaction" language. Item 5
  closed-loop side requires scattering geometry + response law,
  which are blocked on Hard Blockers #3 + #1.
- Convergence claims ("the sweep converges as `dt → 0`"). The
  sweep records baseline values; convergence analysis is a
  separate analytic unit.
- Auto-generated sweep ranges. Caller supplies the
  `(grid_n, spacing_um, dt_s)` tuple list explicitly; the harness
  does not infer or expand a range.
- Runtime defaults in the harness function signature. The
  harness function (`run_ecm_ol_sensitivity_sweep` or equivalent)
  takes mandatory arguments for the tuple list, the per-channel
  prescribed-input scenarios, the reduction choice, and
  `max_grid_cells_total` — no `=...` defaults that hide silent
  caller-elision (per Codex review id=1248). The driver script
  may carry an **explicit fixture** in its `_build_scenarios()`
  helper labeled as a non-production test fixture in
  `metadata.json`, but every value is passed explicitly to the
  harness function; no runtime opt-out path.
- New physical parameters. Sweep amplitudes (prescribed traction
  magnitude, prescribed rate magnitudes) are caller-supplied test
  inputs; no production defaults land.
- New constitutive law. Phase C is open-loop only.

---

## 1. Dimensional analysis

The sweep introduces no new unit reductions beyond the four
already-proven ECM-OL preflight unit chains:

| Channel | Unit chain | Per-(grid, dt) measurement units |
|---|---|---|
| traction | `[nN/μm²] · [s] = [nN·s/μm²]` | accumulator field in `[nN·s/μm²]` |
| stiffness rate | `[kPa/s] · [s] = [kPa]` | stiffness field in `[kPa]` |
| density rate | `[1/s] · [s] = [dimensionless]` | density fields in `[dimensionless ∈ [0,1]]` |
| orientation rate | `[1/s] · [s] = [dimensionless]` | orientation tensor in `[dimensionless ∈ [-1,1]]` |

The sweep dimensions:
- `grid_n` (cell count per axis): integer, dimensionless.
- `spacing_um` (cell edge): `[μm]`. Total grid extent is
  `grid_n · spacing_um` `[μm]`.
- `dt_s`: `[s]`.

The cross-tuple plot's x-axis (e.g., `dt_s` or
`spacing_um × grid_n`) and y-axis (chosen baseline metric per
channel) are both in known units; no unit-mismatch comparison
arises.

### Status

PASS — no new unit reduction; all per-tuple measurements inherit
the already-proven ECM-OL preflight chains. Hard Rule 10 trivially
satisfied (no per-volume vs per-area comparison; everything stays
in the unit space already cleared by the preflight functions).

---

## 2. Boundary cases

| Case | Guard | Failure-kind |
|---|---|---|
| Empty sweep tuple list | early-return with empty index.json (aggregate_status=PASS, 0 scenarios). No `os.makedirs(exist_ok=False)` collision. | (no exception) |
| `grid_n` ≤ 0 or non-int (per tuple) | reject before any per-tuple run | `sweep_grid_n_invalid` |
| `grid_n` is `bool` | reject | `sweep_grid_n_invalid` |
| `spacing_um` ≤ 0 or non-finite (per tuple) | reject | `sweep_spacing_invalid` |
| `spacing_um` is `bool` | reject | `sweep_spacing_invalid` |
| `dt_s` < 0 or non-finite (per tuple) | reject | `sweep_dt_invalid` |
| `dt_s` is `bool` or 0-d bool ndarray | reject (inherit ECM-OL `_validate_dt`) | inherited preflight `dt_invalid` per-channel |
| Total memory bound exceeded | caller-supplied `max_grid_cells_total` cap; overshoot rejects | `sweep_memory_cap_exceeded` |
| Per-tuple per-channel preflight failure (e.g., negative density post-update) | recorded in per-tuple per-channel `failure` payload, scenario marked `status=FAIL` | inherited preflight failure_kind |
| Run-root directory already exists | `os.makedirs(exist_ok=False)` raises loud `FileExistsError` (caller seconds collision) | `FileExistsError` propagates |
| Per-tuple `n_steps` ≤ 0 or non-int | reject | `sweep_n_steps_invalid` |
| `frame_interval` ≤ 0 or non-int (per-tuple if used) | reject | inherited harness `frame_interval` validation |

### Status

PASS — every boundary case has an explicit failure_kind or a
documented "absent = empty result" rule. No silent default ranges,
no auto-shrink, no implicit memory growth.

---

## 3. Conservation invariants

### Per-tuple isolation

Each `(grid_n, spacing_um, dt_s)` tuple runs in a fresh ECM state
(constructed via `make_default_ecm` with the tuple's grid). No
state leaks between tuples; no shared mutable buffer. Determinism
inherits from the preflight functions (no RNG).

### Per-channel isolation within a tuple

The four preflight functions run independently per tuple; each
gets its own ECM copy (or fresh `make_default_ecm` per channel) so
that channel order does not affect per-channel summary metrics.
This differs from the ECM-OL harness `run_ecm_ol_scenario` which
applies channels sequentially on a shared ECM — the sensitivity
sweep wants per-channel baseline isolation.

### Provenance / immutability

- Caller-supplied tuple list is not mutated.
- Caller-supplied prescribed-input factories are not mutated.
- Per-tuple per-channel summary is appended-only; no retroactive
  edit.

### Status

PASS — per-tuple, per-channel isolation enforced; no shared
mutable state; full determinism.

---

## 4. Numerical sanity

### Sweep memory bound

The sweep's working memory is dominated by the largest grid:
`max(grid_n) · max(grid_n_other) · 8 bytes` per ECM field, times
~6 fields, plus scratch for per-step diagnostics. The Sanity Gate
requires the caller to supply an explicit `max_grid_cells_total`
cap as a **mandatory** harness function argument (no `=...`
default in the signature); the gate doc may suggest a worked
example value (e.g., `64 · 64 = 4096` for Phase 1 baseline) for
the script-side fixture, but that suggestion lives only in this
gate doc and in the script's explicit fixture, never as a
function default. The harness raises `sweep_memory_cap_exceeded`
before any allocation when a tuple's grid exceeds the cap.

### Per-tuple dt-rate sanity

Each per-tuple per-channel preflight call runs the
already-proven ECM-OL `_validate_dt` and (for stiffness/density/
orientation rate channels) the per-update value validation. The
sweep does **not** add an additional dt-rate gate — the caller
supplies tuple `dt_s` and per-channel rate amplitudes, and any
violation surfaces with the inherited preflight failure_kind.

### Per-tuple wall-clock bound

Each tuple's per-channel run is `O(grid_n² · n_steps)` work, plus
`O(n_steps · frame_interval⁻¹)` artifact emission. The sweep does
not bound wall-clock; the caller's tuple list and `n_steps` choice
are the bound. A long sweep (many tuples × many steps) is the
caller's responsibility.

### Float precision

float64 throughout (inherited from ECM-OL preflight).

### Status

PASS — sweep memory cap is caller-supplied with explicit failure
mode; dt-rate sanity inherits from preflight; wall-clock is
caller-bounded.

---

## 5. Sign / sense check

| Quantity | Direction | One-line check |
|---|---|---|
| Prescribed traction (sweep input) | non-negative (open-loop preflight contract) | inherit; rejection on negative is verified by the preflight test suite |
| Prescribed stiffness rate (sweep input) | signed (open-loop preflight contract) | inherit |
| Prescribed density rate (sweep input) | signed | inherit |
| Prescribed orientation rate (sweep input) | signed, symmetric | inherit |
| Per-tuple summary metric | direction matches the field's sign convention | inherit per-channel |
| Cross-tuple plot y-axis | same sign as per-tuple summary metric | by construction (no transformation) |

### Status

PASS — no new sign convention is introduced. The sweep simply
records the preflight outputs in their already-locked sign
conventions.

---

## 6. Measurement-protocol consistency

Per Hard Rule 11, the sweep's measurement modality is **the
post-step ECM field** evaluated at the same `(grid_n, spacing_um)`
as the per-tuple ECM. The cross-tuple plot's y-axis is a **scalar
reduction** of the post-step ECM field per tuple (e.g., field
maximum, mean, sum). The reduction must be documented in the
sweep's metadata so a downstream reader knows whether the y-axis
is `max`, `mean`, or another reduction; the harness records this
explicitly.

### Wording lock (Hard Rule 11 protection)

Every artifact emitted by the sweep — `summary.html`,
`metadata.json`, `index.json`, plot titles, console output —
must use **"open-loop sweep baseline"** and **never** "Item 5
satisfied" or "closed-loop sensitivity". The locked phased plan
§1 Phase C makes this an explicit wording boundary; a runtime
meta-test (parallel to Phase B's
`test_stimulus_monotonicity_does_not_satisfy_closed_loop_item_1`)
should assert that the sweep harness's status payload does not
contain the forbidden phrases.

### Status

PASS — measurement modality matches per-tuple ECM field outputs;
reduction is explicit; wording boundary enforced by docstring +
runtime meta-test.

---

## 7. Magic-Number Block check

Every numeric in the Phase C harness is one of:

| Symbol | Source | Magic-Number Block status |
|---|---|---|
| `grid_n`, `spacing_um`, `dt_s`, `n_steps`, `frame_interval` (per tuple) | runtime caller-supplied (no project default, no harness signature default) | **not a magic number** — caller-supplied test inputs only |
| `max_grid_cells_total` cap | caller-supplied numerical safety bound (mandatory positional arg in the harness function; the gate doc may suggest `4096` for Phase 1 baseline as a worked example, but the suggestion lives only in this gate doc, not in code) | named numerical safety bound; not a physics constant; **never a function default** |
| Sweep amplitudes (prescribed traction magnitude, prescribed rate magnitudes per channel) | caller-supplied per-channel scenario inputs | **not a magic number** — same as 6.3b runner pattern |
| Reduction choice (max, mean, sum) | caller-supplied; no project default | **not a magic number** — protocol selection |
| Reference biology table from the closed-loop lock §3 | doc only, not used in Phase C code | N/A |

Magic-Number Block test pass per locked principle:

1. **Derivable**: every value in code is runtime caller input or
   the named `max_grid_cells_total` safety bound.
2. **Grid-invariant**: the sweep's *purpose* is to test grid/dt
   sensitivity, so by definition each tuple's outputs depend on
   the tuple's grid/dt; this is intentional and is the unit of
   work, not a magic-number violation.
3. **Not fitting**: no value in code is chosen to make a specific
   sweep result pass.

### Status

PASS — Magic-Number Block clean. The `max_grid_cells_total`
default suggestion in the docstring is not a runtime default; the
caller must explicitly pass it.

---

## 8. Visual deliverable plan (post-Sanity-Gate code commit)

Aligns with the precedent set by `scripts/run_p1_alpha_gate.py`,
`scripts/run_ecm_ol_harness.py`, and
`scripts/run_protrusion_coupled_fa_harness.py`. **Out of scope for
the Sanity Gate commit itself** — the Sanity Gate gates the code
commit; the visible-deliverable run is post-commit and is the
primary PI-facing output.

Planned (not committed by this gate):
- `scripts/run_ecm_ol_sensitivity_sweep.py` — sweeps a
  caller-supplied tuple list across the four channels.
- Per-tuple per-channel `summary.html` + `metadata.json` +
  `diagnostic_<tuple>_<channel>.png`.
- Cross-tuple per-channel `summary.html` showing the baseline
  metric vs `(grid_n, dt_s)`.
- Top-level `index.json` with `aggregate_status`,
  `intentional_failures_observed`, `unexpected_failures`,
  `unexpected_passes`.
- An **explicit non-production fixture** in the script's
  `_build_scenarios()` helper (e.g., 3 grid sizes × 3 dt values),
  labeled `"fixture_kind": "non_production_smoke"` in every
  `metadata.json` and in the top-level `index.json`. The harness
  function still takes every value as a mandatory argument; the
  fixture exists only on the script side. CLI may override every
  fixture entry. **No runtime default tuple list lives in the
  harness function signature.**

---

## 9. Test catalog (~6, per locked phased plan §8 Phase C)

Each test exercises one Sanity Gate item.

| # | Test name | Gate item |
|---|---|---|
| 1 | `test_open_loop_sweep_runs_4_channels_explicit_fixture` | §0 + §3 isolation |
| 2 | `test_open_loop_sweep_grid_spacing_variation` | §1 + §3 + §6 |
| 3 | `test_open_loop_sweep_dt_variation` | §1 + §4 + §6 |
| 4 | `test_open_loop_sweep_summary_reports_baseline_only` | §6 wording (Hard Rule 11) |
| 5 | `test_open_loop_sweep_does_not_satisfy_closed_loop_item_5` | §6 wording (meta-test, parallel to Phase B Test 4) |
| 6 | `test_open_loop_sweep_memory_cap_rejects_oversize_grid` | §4 |

Plus boundary-case tests for §2 invalid-input rejection (grid_n
bool, dt negative, etc.); these can ride in the same test file as
parametrized cases.

---

## 10. Gate verdict

§1, §2, §3, §4, §5, §6 PASS. Magic-Number Block (§7) clean.
Visual deliverable plan (§8) is post-commit and does not block.
Test catalog (§9) maps every gate item to at least one test.

The gate clears for executable code in two commits:

1. `scripts/run_ecm_ol_sensitivity_sweep.py` — the sweep driver.
   Inline harness (no new module) following the
   `scripts/run_protrusion_coupled_fa_harness.py` pattern.
2. `tests/test_v2_ecm_ol_sweep.py` — the ~6+ tests per §9.

If Codex review surfaces a missing reduction or hidden numeric,
the gate flips to BLOCKER with the three-options template.

### Outstanding before code lands

- Codex review of this gate document (5 focus per locked phased
  plan + Codex impl id=1239 + id=1245):
  1. Wording lock: every artifact says "open-loop sweep baseline",
     never "Item 5 satisfied"
  2. Per-tuple per-channel isolation enforced (no shared mutable
     ECM)
  3. `max_grid_cells_total` cap caller-supplied with explicit
     failure mode
  4. f8cdff3 hygiene (seconds-resolution timestamp,
     `os.makedirs(exist_ok=False)`, aggregate-against-expectation
     status)
  5. Runtime meta-test enforces wording boundary (parallel to
     Phase B's
     `test_stimulus_monotonicity_does_not_satisfy_closed_loop_item_1`)

---

## 11. References

- Locked phased plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- ECM-OL preflight gate (inherited):
  `acs/v2/dynamics/ecm_open_loop.py` module docstring
- ECM-OL harness pattern (inherited):
  `acs/v2/ecm_open_loop_harness.py`
- 6.3b runner pattern (inherited f8cdff3 hygiene + index.json
  shape): `scripts/run_protrusion_coupled_fa_harness.py`
- Phase B precursor (just landed):
  `tests/test_v2_ecm_open_loop.py` Phase B section (`66f06d6`,
  `ada3728`)
- Hard Rule 10 (CLAUDE.md): dimensional comparison verification
- Hard Rule 11 (CLAUDE.md): measurement-protocol consistency
- Memory `rule10_unit_derivation_in_docs.md`: inline derivation in
  any markdown unit comparison

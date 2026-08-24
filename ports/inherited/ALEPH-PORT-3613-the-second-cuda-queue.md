# ALEPH-PORT-3613 — the second CUDA queue, and an enumeration a new law cannot fall out of

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3613` |
| Lane | `46143f30` Lane G7b (CUDA track — operations). **`ALEPH-PORT-3607` reopened, not replaced.** |
| Status | `PROPOSED` |
| Written | `2026-08-01` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Exists because | `-3607`'s queue ran at 02:36–02:43 and **finished**. Its log ends *"queue finished — every job is complete or abandoned. Exiting."*, and the A5000 has been idle since while the grant runs to 2026-08-04. In the six hours between, the module went from ~30 kernels to **81** and from 19 of 30 wired connector rows to **30 of 30** (`-3610`, `-3611`), and four law families that did not exist when `j1` ran are now registered. **`j1` measured a tree that no longer exists.** This entry re-arms the same runner with a second queue, and its central change is not a job: it is that `j1`'s enumeration was a **hand-written list of ten families** and is replaced by one driven off `law_cases.LAW_CASE_BUILDERS`, so the next lane's law is graded on the next run instead of being silently omitted. |

---

## 1. Aleph API

This entry authorises **no change to any module under `aleph/`**, and that is a harder constraint
here than in `-3607` because two of the four new jobs want something the tree does not quite expose.
Where that happened the job was reshaped to fit the frozen surface, and §14.2 says which question
each reshaping gave up.

```python
# EDITED — authorised by this entry, all within Lane G7's own declared paths
scripts/cuda_jobs.py               # + the registry-driven enumeration, + four jobs k1..k4
scripts/unattended_runner.py       # + four JobSpecs; the queue machinery is untouched
scripts/aleph-cuda-runner.service  # + the deadline moved to the grant's own expiry
tests/scripts/test_unattended_runner.py   # + 15 controls, 25 -> 40
```

```python
# CONSUMED, READ-ONLY — nothing here is edited by this entry
from aleph.runtime.law_cases import (
    LAW_CASE_BUILDERS, SWEEP_GRID, _true_variants, brownian_ratchet_state,
    brownian_ratchet_kernel_site_forces, membrane_state,
)
from aleph.runtime.law_kernels import PositionPrecision, DragDissipationVariant, SelfForceVariant
from aleph.runtime.parity import compare_law_case, law_case_against_itself
from aleph.runtime.residency import ResidentMembraneState
from aleph.runtime.backend import ScatterMode, WarpBackend, check_gpu_authorization
from aleph.vertical.assembly import build_vertical, recovered_tension_pn_per_um
from aleph.vertical.membrane import HelfrichMembrane
from aleph.vertical.relax import FLOAT32_UNIT_ROUNDOFF, relax_to_equilibrium
```

**`scripts/gpu_preflight.py` is consumed as a subprocess and is never imported, never edited and
never bypassed** — unchanged from `-3607` §1, and it is still the reason this entry is worth
anything.

**`law_cases._true_variants` is private and is used deliberately.** It is the only place in the tree
that maps a builder to the variant *enum* it is driven at, and `every_law_case()` — which is public —
returns built cases from which the enum cannot be recovered. §9's
`test_the_enumeration_refuses_when_the_registry_and_its_variant_table_disagree` makes the dependency
loud: if a later lane renames or reshapes it, the enumeration **refuses** rather than quietly
grading fewer laws, which is the exact failure this entry exists to close one level up.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` — cited for the record only |
| Source path / symbol | **none named, none read.** No file under that tree was opened by this lane. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable — nothing was read |

Nothing was ported. The artefacts this entry rests on are Aleph's own: `ALEPH-PORT-3607`,
`-3608` §10.2, `-3609` §14.6, `-3610` §14.11/§14.12, `-3611` §14.1/§14.7/§14.10/§14.12, the five
`j*.json` artefacts on `gbook`, and `docs/design/GPU_POLICY.md`.

## 3. Why source-derived porting beats clean-room

**It does not, and this entry says so rather than manufacturing a reason** — the same answer
`-3607` §3 gave, for the same reason. A job queue is a convention. What is *not* a convention is the
enumeration rule in §4, and that rule is a response to something that happened in this repository
six hours ago, not to anything a provider wrote.

## 4. Physical or mathematical law represented

No physical law. `-3607` §4's refusal theorem **T1** and its anti-law **T2** are inherited unchanged
and are still the load-bearing statements — the runner's single in-edge to `LAUNCHED` is untouched
by this entry, and §9 re-drives all three authorization refusals against the *new* queue rather than
assuming the old proof transfers.

This entry adds one statement of its own, about coverage rather than safety:

> **T3 (the enumeration theorem).** The set of cases graded by `k1` is a function of
> `law_cases.LAW_CASE_BUILDERS` and of the variant enums those builders are typed against — never
> of a list written in `scripts/`. Consequently a law family added to `aleph/runtime/` after this
> entry is graded on the next run of the queue with no edit here, and a family added with a variant
> **axis** this enumeration does not know about causes a **refusal**, not a silent partial grading.

The corollary is the whole reason T3 is written down:

* **T3a.** `-3607`'s `j1` was correct when it ran and wrong six hours later, and **nothing in the
  artefact it produced could have told a reader that.** A green table over ten families reads
  exactly like a green table over seventeen. The failure mode is not a wrong number; it is a number
  whose *scope* silently shrank relative to the tree. T3 makes the scope a derived quantity, and
  §8's `test_the_enumeration_covers_every_registered_law_case_builder` makes it checkable.

And one anti-law, in the same shape as T2 and for the same reason:

> **T4.** No job in this queue writes, widens, relaxes, or proposes a replacement for any
> `DEFAULT_*_ULP_BUDGET`, for `DESCENT_SLACK_ROUNDOFFS`, or for any tolerance in `aleph/`. `k2`
> exists to find a budget that fails on the device and **report it**; `k4` exists to measure where
> the descent window is on the device and **stop there**. `-3609` §14.6 and `-3611` §14.10 are
> already open PI decisions of exactly this kind, and a queue that quietly tuned around either
> would destroy the evidence that raised them.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| job wall-clock budget | s | s | `(0, 10800]`, declared per job, never shared |
| job **internal** soft budget | s | s | `(0, job budget)`; `k4` only — see the singular cases |
| declared device memory budget | GB | 1e9 B | `(0, 16]` on the A5000; passed straight to the preflight |
| measured disagreement | float32 ULP of the field scale | dimensionless | `[0, inf]`; `inf` is a real value |
| gate **pressure** | dimensionless | — | `measured_ulp / effective_ulp_budget`; `[0, inf]`. `> 1` is a failed gate |
| gate **round-off ratio** | dimensionless | — | `effective_ulp_budget / roundoff_ulp_bound`; `(0, inf)` |
| descent slack | multiples of `u32` | dimensionless | `(0, 8]`, swept; `u32 = FLOAT32_UNIT_ROUNDOFF` |
| free-velocity anchor residual | pN | 1e-12 N | `[0, inf)`; the float64 reference is exactly `-0.0` |

Singular and boundary cases, each with the behaviour required:

- **Every refusal case in `-3607` §5 is inherited verbatim** — no authorization, expired, wrong
  host, wrong device, busy device, unreadable machine, hung job, crashed job, complete artefact,
  stop file, deadline. None of them is re-derived here and all of them are re-driven in §9 against
  the new queue.
- **A registered builder with an unknown variant axis.** The enumeration raises `JobRefused` naming
  the builder and the parameter. It does **not** grade the builder at its primary axis only, because
  a partially-enumerated family produces a green table that means less than it looks like.
- **A builder in the registry with no entry in the variant table, or the reverse.** `JobRefused`.
- **A case that raises inside `compare_law_case`.** Recorded as `error`, and **never** counted as a
  caught mutant — `-3607` §5's rule, inherited: a mutant detected because the harness crashed is a
  different observation from a mutant detected by a budget.
- **A true kernel over budget on the device.** This is `k2`'s reason to exist. Recorded as a
  finding with the measured ULP, the effective budget and the round-off bound beside it. **The
  budget is not touched** (T4).
- **`k4`'s soft budget exhausted.** Enumeration stops, every unattempted `(config, multiple)` pair
  is listed in `not_attempted`, and the artefact is still written `complete: true`. A job killed at
  its runner timeout writes nothing at all, and three hours of device time that produced no readable
  artefact is worse than two hours that produced one with a stated hole in it.
- **A non-converged relaxation in `k4`.** A row, not a refusal — for `k4` specifically,
  non-convergence *is* the measurement, because the window is defined by which slacks converge.
  This is the one place this queue reads a `converged=False` row as data, and it is stated here so
  it is not mistaken for a relaxation of `-3607` §5's rule.
- **An empty anchor mask in `k3`.** Refused. `ALEPH-PORT-3606`'s O2 — a control asserting over an
  empty set passes vacuously — and a *probe* that measures an empty set reports `0.0` and reads
  like a restored exact anchor.

Invariants, each with the test that asserts it:

- **I1–I10 of `-3607` §5 are inherited unchanged**, and I1 (no launch without a passing preflight)
  is re-asserted against this queue by
  `tests/scripts/test_unattended_runner.py::test_the_guard_still_refuses_all_three_ways_with_the_second_queue_armed`.
- **I11.** The enumeration reaches every registered builder —
  `::test_the_enumeration_covers_every_registered_law_case_builder`.
- **I12.** An unknown variant axis refuses rather than being dropped —
  `::test_an_unknown_variant_axis_refuses_rather_than_being_silently_dropped`.
- **I13.** Registry and variant table must agree —
  `::test_the_enumeration_refuses_when_the_registry_and_its_variant_table_disagree`.
- **I14.** A variant is "wrong on purpose" if **any** of its axes is —
  `::test_every_wrong_variant_is_wrong_on_at_least_one_axis`.
- **I15.** The second queue grades strictly more cases than the first did —
  `::test_the_enumeration_grades_more_cases_than_the_first_queue_did`.
- **I16.** Every graded gate carries pressure *and* its round-off bound —
  `::test_the_budget_census_reports_pressure_and_the_roundoff_bound_for_every_graded_gate`.
- **I17.** A gate over budget is a finding, never a green row —
  `::test_the_budget_census_reports_a_gate_over_budget_as_a_finding`.
- **I18.** Nothing in the census widens a budget — `::test_the_census_does_not_widen_a_budget`.
- **I19.** The anchor probe measures both devices and asserts neither —
  `::test_the_exact_anchor_probe_measures_both_devices_and_asserts_neither`.
- **I20.** The anchor probe refuses an empty mask — `::test_the_exact_anchor_probe_is_not_vacuous`.
- **I21.** `k4`'s rows are ordered so an early stop still covers `1·u32` and `2·u32` on every
  configuration —
  `::test_the_descent_window_rows_are_ordered_so_an_early_stop_still_covers_both_reference_multiples`.
- **I22.** `k4` records what it did not attempt —
  `::test_the_descent_window_job_records_rows_it_did_not_attempt`.
- **I23.** `k4` proposes no tolerance — `::test_the_descent_window_job_proposes_no_tolerance`.
- **I24.** Each new job declares its own budgets —
  `::test_the_second_queue_is_declared_with_its_own_budgets`.
- **I25.** The second queue never overtakes the first —
  `::test_the_second_queue_runs_after_the_first_and_in_declared_order`.
- **I26.** The first queue's artefacts are not recomputed —
  `::test_the_second_queue_does_not_recompute_the_first_queues_artefacts`.

### 5a. The authorization record — read, never written

Unchanged in shape from `-3607` §5a, and **it now exists and is valid**. Measured on `gbook`, quoted
from the file rather than from memory:

```json
{
  "granted_by": "PI",
  "granted_at": "2026-08-01T01:35:00+09:00",
  "expires_at": "2026-08-04T00:00:00+09:00",
  "host": "GBook",
  "devices": [0],
  "note": "PI instruction 2026-08-01 01:34 KST: physical device is gbook; unattended runs to continue through Sunday while the PI is away."
}
```

**An agent may never create, extend, repair or move this file** — `CLAUDE.md` §3,
`GPU_POLICY.md`: *an agent that writes its own permission slip has not been authorized*. §9's
negative controls are written against this exact shape and plant their own copies in a temporary
`ALEPH_DATA_ROOT`; none of them touches the real one.

**The service deadline is moved from `2026-08-03T09:00:00+09:00` to `2026-08-04T00:00:00+09:00`, to
be the grant's own expiry rather than a number beside it.** A deadline earlier than the grant wastes
authorized device time; a deadline later than the grant leaves the runner polling and logging
refusals after permission has ended, which is safe but is noise standing in for a stop. Equal is the
only value that is neither.

## 6. Source evidence class and known retractions

Nothing was read from the provider. The *Aleph* artefacts this entry rests on, each inherited
explicitly:

| Artefact | What it claims | Status carried here |
|---|---|---|
| `-3607` §13.1, `j1_law_parity.json` | 38 wrong kernels caught on both devices, 0 caught-on-CPU-but-not-CUDA, 1 declared blind spot | **superseded in scope, not retracted.** True of the ten families that existed at 02:36. `k1` re-measures over seventeen |
| `-3611` §14.12 | family M's 1.5 ULP budget has "no room to absorb a device"; family L computes a logarithm and has "the least claim to transfer" | `INHERITED_UNVERIFIED` on CUDA — `k2` is what tests it |
| `-3609` §14 | the CUDA descent window is `{1}·u32` against warp-CPU `{2}`, on **two** configurations | `INHERITED_UNVERIFIED` beyond those two configurations — `k4` widens the sample and nothing else |
| `-3611` §14.10 | the free-velocity exact anchor is lost by one float32 ULP of association | measured on the **warp CPU device**; `k3` measures it on CUDA, where `wp.log` and FMA contraction are different |
| `-3607` §13.1, `j2` | `ORDERED` bit-identical on this A5000 | inherited as a **measurement**, and it is the one prior figure this queue does not re-take |

**A retraction this entry makes about its own lane.** `-3607` §1 lists `law_cases` symbols one by
one and `scripts/cuda_jobs.py` enumerated ten families from a hand-written list. That was adequate
when it was written and became wrong the moment another lane registered an eleventh, **without any
test going red and without the artefact changing shape**. This entry treats that as the defect it
is rather than as a scope note: T3, I11 and §9's negative control exist because of it.

## 7. Independent oracle or derivation

Four, of three kinds.

**(a) For the runner: the guard's own exit codes.** Unchanged from `-3607` §7(a). §9 drives the
real, unmodified `scripts/gpu_preflight.py` against planted records that must be refused, and
observes that nothing launched. Not a mock.

**(b) For `k1`: the warp CPU device, in the same process.** This is `-3607`'s discipline and it is
the reason `j1` was worth anything, so it is kept exactly: the identical enumeration is driven on
`cuda:0` **and** on the CPU device inside one process, and *"do the gates still catch the wrong
kernels on CUDA"* is answered by set difference. **A mutant caught on CPU and missed on CUDA is a
finding.** A hand-written table of expected failures would encode `-3604` §9a's Newton mirror,
`-3610`'s three isotropic-drag blind variants and `-3605` §10.5's ungradable rest nucleus by hand,
and would then be a statement about its author.

Measured on this Mac's warp CPU device at `a92e38f`, **before** any device run, so the device figures
have something fixed to be compared against:

| | `j1`, 2026-08-01 02:36 | this entry, same day |
|---|---|---|
| cases enumerated | 49 | **91** |
| wrong on purpose | 38 | **74** |
| caught on the CPU device | 38 | **70** |
| blind on the CPU device | 1 | **4** |
| true kernels over budget | 0 | **0** |
| harness self-comparison non-zero | 0 | **0** |

The four blind variants, quoted from the run as the enumeration labels them, are
`tether_step[WRONG_NEWTON_MIRROR]` — `-3604` §9a's declared blind spot, bit-identical to the true
kernel at every input — and `isotropic_drag[WRONG_ISOTROPIC|TRUE]`,
`isotropic_drag[WRONG_AXIS_SWAPPED|TRUE]`, `isotropic_drag[WRONG_UNWEIGHTED|TRUE]`, which are three
`-3610` §4.3 mutants that degenerate to the true kernel at `ratio = 1, L = 1` — the configuration
the isotropic case *is*. **All four are blind on the CPU device, which is the side `k1` compares
against, so none of them can produce a false cross-device finding.** A family label is the registry
key with `_law_case` stripped, and a variant label names every axis, so `isotropic_drag`'s second
field is its `rate_variant`.

**(c) For `k2`: each ledger's own declared budget, and the round-off bound beside it.** No tolerance
is invented. `pressure = measured / effective_ulp_budget` says how close a gate is to failing;
`roundoff_ratio = effective_ulp_budget / roundoff_ulp_bound` says whether the budget is tight
relative to the assembly's own error growth. Reporting only the first would let a figure inside a
*loose* budget read as a tight result — `-3607` §10's rule, generalised.

The CPU-device baseline, measured the same way over the three position modes (156 graded gates):

| pressure | mode | gate | measured / budget |
|---|---|---|---|
| **0.7085** | `global_f32` | `strain_stiffening_cable.energy_start_pn_um` | 45.34 / 64 |
| 0.6452 | `global_f32` | `chromatin.forces_pn` | 7.10 / 11 |
| 0.6378 | `local_f32` | `tether_step.energy_start_pn_um` | 20.41 / 32 |
| 0.5863 | `global_f32` | `chromatin_owner.forces_pn` | 6.45 / 11 |
| 0.5365 | `global_f32` | `series_joint.forces_a_pn` | 68.68 / 128 |

Seven gates are above half their budget on the CPU device, and **the tightest is family L's — the
transcendental one `-3611` §14.1 names as least likely to transfer.** That is the row to read first
in `k2`'s artefact.

**(d) For `k3`: an exact anchor that is an identity of the frozen law, and another lane's control
that asserts it.** At `v = v₀` the logarithm's argument is exactly `1.0`, so the float64 reference
produces exactly `-0.0`. That is an identity, not a measurement, so the probe has a zero to compare
against that no device can move.

**The probe does not re-derive it, and this is a deliberate boundary rather than an omission.** The
per-*site* float64 reference is reachable only through `law_cases._ratchet_connector`, a private
symbol, while the public `LawCase.reference` returns the *assembled* field — 5 nodes against 12
sites — so the anchor mask does not apply to it. Recomputing the reference in `scripts/` would be
this lane writing physics. So `k3` measures the **kernel's** residual on both devices from the
public site helper and cites
`tests/runtime/test_active_law_parity.py::test_the_free_velocity_exact_anchor_is_LOST_by_one_ulp_of_association`
— `ALEPH-PORT-3611`'s own control, which already asserts the exact `-0.0` with `==` — instead of
asserting it a second time in a second place.

Measured here on the CPU device: **1** site at the free velocity of 12, force residual
**1.4217e-07 pN** and power residual **5.6691e-08 pN·µm/s**, against a stall force of **3.3610 pN**.
`k3` asks only whether CUDA's `wp.log` and its FMA contraction move those numbers.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/scripts/test_unattended_runner.py::test_the_enumeration_covers_every_registered_law_case_builder` | Every name in `law_cases.LAW_CASE_BUILDERS` appears in the enumeration `k1` drives. **T3 / I11**, and the direct control on `-3607`'s defect: a family registered after this entry is graded without an edit here. |
| Positive | `::test_the_enumeration_grades_more_cases_than_the_first_queue_did` | ≥ 91 cases against `j1`'s 49, and strictly more than `law_cases()` — a scope that shrank is red. I15. |
| Positive | `::test_every_wrong_variant_is_wrong_on_at_least_one_axis` | `chromatin_owner[TRUE\|WRONG_SCALAR_MEAN]` and `anisotropic_drag[TRUE\|WRONG_RAYLEIGH_HALF]` are flagged wrong. A flag read off the primary axis alone would call them correct and then count them as passing kernels. I14. |
| Positive | `::test_the_budget_census_reports_pressure_and_the_roundoff_bound_for_every_graded_gate` | Every graded gate carries `ulp`, `effective_ulp_budget`, `roundoff_ulp_bound`, `pressure` and `roundoff_ratio`. I16. |
| Positive | `::test_the_exact_anchor_probe_measures_both_devices_and_asserts_neither` | The probe returns a row per device and no field claiming a repair, a preference or a verdict. I19. |
| Positive | `::test_the_second_queue_is_declared_with_its_own_budgets` | Each of `k1..k4` declares its own `timeout_s` in `(0, 10800]` and its own `memory_gb` in `(0, 16]`, and no two share a value by accident of a shared constant. I24. |
| Positive | `::test_the_second_queue_runs_after_the_first_and_in_declared_order` | The combined queue attempts `j1..j5` then `k1..k4`, in that order. I25. |
| Positive | `::test_the_second_queue_does_not_recompute_the_first_queues_artefacts` | With `j1..j5` artefacts present, a pass launches only the `k` jobs and reports the `j` jobs `skipped`. I26. |
| Positive | `::test_the_descent_window_rows_are_ordered_so_an_early_stop_still_covers_both_reference_multiples` | The first 24 rows are `1·u32` and `2·u32` over all twelve configurations, so a budget that stops the sweep still answers the question `-3609` left open. I21. |

## 9. Deliberately failing negative control

**The three that matter are `-3607` §9's three, re-driven against this queue**, because a proof that
the guard refuses is a proof about the queue it was run on. Each plants a record in a temporary
`ALEPH_DATA_ROOT` and executes the **real** `scripts/gpu_preflight.py`.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/scripts/test_unattended_runner.py::test_the_guard_still_refuses_all_three_ways_with_the_second_queue_armed` | With the full nine-job queue armed: a record valid in every field but expired → **exit 2**, 0 launches; a record naming another host → **exit 2**, 0 launches; no record at all → **exit 2**, 0 launches. All three in one control so the queue cannot be re-armed with one of them silently untested. I1. |
| Negative (must fail) | `::test_an_unknown_variant_axis_refuses_rather_than_being_silently_dropped` | A builder registered with a keyword-only variant axis the table does not know makes the enumeration raise `JobRefused` naming the builder **and** the parameter. A version that dropped the axis and graded the primary alone passes every other control here. I12. |
| Negative (must fail) | `::test_the_enumeration_refuses_when_the_registry_and_its_variant_table_disagree` | A registry with a builder absent from `_true_variants()` refuses. This is the control that makes §1's use of a private symbol safe. I13. |
| Negative (must fail) | `::test_the_budget_census_reports_a_gate_over_budget_as_a_finding` | A census fed a field whose measured ULP exceeds its effective budget reports it in `over_budget` and sets the verdict to `FINDING`. A census that could only print green is the shape `-3607` §9's last row already refuses for `j1`. I17. |
| Negative (must fail) | `::test_the_census_does_not_widen_a_budget` | The census emits no `suggested_budget`, no `proposed_budget`, no `relaxed_*`; it does not mutate the case's `ulp_budget` mapping; and the source of `scripts/cuda_jobs.py` assigns nothing to `DESCENT_SLACK_ROUNDOFFS` or to any `DEFAULT_*_ULP_BUDGET`. **T4.** I18. |
| Negative (must fail) | `::test_the_descent_window_job_proposes_no_tolerance` | The window job's output carries no `recommended_slack`, no `should_be`, no `proposed_*`, and its source assigns nothing to `DESCENT_SLACK_ROUNDOFFS`. `-3609` §14.6 is a PI decision and this queue measures into it, never over it. I23. |
| Negative (must fail) | `::test_the_descent_window_job_records_rows_it_did_not_attempt` | With a soft budget already spent, the job attempts zero rows, lists every `(config, multiple)` pair in `not_attempted`, and still writes a complete artefact. I22. |
| Negative (must fail) | `::test_the_exact_anchor_probe_is_not_vacuous` | A state with no site at the free velocity makes the probe refuse. Measuring an empty mask returns `0.0`, which reads as a *restored* exact anchor — the strongest possible wrong answer. `-3606` O2, in a probe rather than a control. I20. |

**What these deliberately do not do.** They do not stub the preflight for the three authorization
cases. A test that asserts "the runner stops when *my mock* says no" grades the mock —
`-3607` §9, unchanged.

## 10. Numerical and precision envelope

* **Compute channel:** float32 (`ALEPH-DQ-107`), `WarpBackend.dtype`'s default, changed by no job.
* **Scatter mode:** `ORDERED` for every job in this queue, with no exception — `j2` measured it
  bit-identical on this A5000 and every budget in G1–G5c rests on that.
* **Position representation:** `k1` runs at `GLOBAL_F32` so its figures are comparable with `j1`'s
  and with every ledger's CPU table. `k2` runs all three modes and **answers nothing** — the open
  question is `docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` and it is the PI's.
* **`k4`'s descent energy is the device float32-term energy**, from `ResidentMembraneState`, exactly
  as `-3609` measured it. That is the whole point: `j4` ran the descent on a host float64 energy and
  converged 12/12, and the question `-3609` left open is about the *device* channel.
* **Tolerances:** none is invented, widened or proposed by this entry (T4). Every budget is the
  `DEFAULT_*_ULP_BUDGET` its own ledger derived; `DESCENT_SLACK_ROUNDOFFS = 2.0` in
  `aleph/vertical/relax.py` is read and **not** changed, and `k4` sweeps a *parameter* of
  `relax_to_equilibrium` rather than that constant.
* **Outside the envelope:** a `nan` is written as `NaN` and marked unreadable, never coerced; an
  `inf` ULP is a measurement (`-3605` §10.5) and survives the artefact writer.

## 11. Production-backend residency and transfer

| Stage | Host | Device | Transfer |
|---|---|---|---|
| queue state, artefacts, log | `gbook` disk | — | none |
| preflight | `gbook` host process → `ssh GBook nvidia-smi` | reads device metadata only | no allocation |
| `k1` / `k2` reference | host, float64 NumPy | — | never uploaded |
| `k1` / `k2` candidate | host state → device arrays per case | `cuda:0` | one upload per case, one download per field |
| `k3` anchor probe | host state → device arrays, once per device | `cuda:0` and `cpu` | two launches total |
| `k4` descent | positions on host, **energy resident** | `cuda:0` | `-3609`'s resident chain: 3 calls / 4 elements per evaluate, 0 bulk arrays. The forces stay NumPy float64 by construction, so the only thing under test is the descent criterion's channel |

The runner allocates no device memory of its own; every array comes from `backend.zeros`/`array`
inside frozen case code, as `-3601` §1 requires.

## 12. Comments and docstrings to discard

Nothing was read, so nothing is inherited and nothing needs discarding. `-3607` §12's banned phrase
is inherited and extended: no artefact may say a job was *"skipped for safety"* or *"deferred"* —
the word is `refused`, with the exit code beside it. **This entry adds a second ban**: no artefact
may carry a `suggested`, `recommended` or `relaxed` tolerance of any kind. A measurement that names
its own replacement value has stopped being a measurement, and two open PI decisions
(`-3609` §14.6, `-3611` §14.10) are sitting on exactly this line.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** See §13.1 for what was observed on `gbook`. Acceptance requires an artefact set from `k1..k4` carrying `preflight.exit_code == 0`, and the runner has produced none at the time this section was written. |
| Reviewer | **Agent-proposed, unratified.** No PI decision is recorded here and no `decided_by` field exists in this file. |
| Rollback | `systemctl --user disable --now aleph-cuda-runner.service`; or `touch ~/.aleph_runner/STOP`; or `git revert` the commits that added the four `k` jobs — the `j` jobs and the runner machinery are `-3607`'s and survive. Nothing under `aleph/` changes, so no physics, no test and no gate depends on this entry. |

### 13.1 What was observed on `gbook`, at `562fd71`

Not asserted — run, on that machine, and quoted from its own output.

| Observation | Result |
|---|---|
| `git status` on `gbook` **before** anything | `0` modified paths, `0` stashes, branch `main` at `8f8281f` — checked first, per the brief |
| checkout updated `8f8281f` → `562fd71` | incremental bundle, `git merge --ff-only`, clean tree before and after |
| `tests/scripts` on `gbook`'s own `aleph` env | **43 passed** in 13.1 s (was 25 at `-3607`) |
| the real preflight against the real record | `[authorization] OK: authorized until 2026-08-04T00:00:00+09:00 by PI` … `[preflight] PASS` — **exit 0** |
| the three refusals, live, against the **real** guard in a temporary `ALEPH_DATA_ROOT` | see below |
| the real record afterwards | `GBook`, `2026-08-04T00:00:00+09:00`, `granted_by: PI` — **byte-unchanged** |
| service state now | `active`, `enabled`, `Linger=yes`, `MemoryMax=12G`, `CPUQuota=400%`, `TasksMax=256`, `Nice=10`, deadline `2026-08-04T00:00:00+09:00` |

**The three refusals, each driven against the unmodified `scripts/gpu_preflight.py` with the nine-job
queue armed, on `gbook` itself.** The runner's own guard line is quoted verbatim in each case, and
each run wrote **0** artefacts:

| case | the guard's own words | exit | launches | artefacts |
|---|---|---|---|---|
| no record | *"no authorization file at …/gpu_authorization.json … Absence is not permission."* | **2** | 0 | **0** |
| expired | *"authorization expired at 2026-07-31T21:22:09+09:00. A grant does not renew itself."* | **2** | 0 | **0** |
| wrong host | *"authorization names host 'some-other-machine', not 'GBook'. A grant for one machine is not a grant for another."* | **2** | 0 | **0** |

**Both stop lines, executed live rather than argued for.**

| line | observed |
|---|---|
| `touch ~/.aleph_runner/STOP` | runner logged *"STOP file present — exiting"*, service `inactive`, **0** surviving processes, **0** artefacts written |
| `systemctl --user disable --now aleph-cuda-runner.service`, **executed while `k4` was mid-flight with 194 MiB resident on the A5000** | returned in **0.25 s**; service `inactive` and `disabled`; **0** surviving job processes, **0** surviving runner processes, **0** compute apps on the device |

**The cost of that second verification is recorded rather than hidden:** it consumed `k4`'s first of
two attempts. `k4` is running on attempt 2 and has no retry left, and the attempt record is accurate
— the job *was* launched once and did not complete. `attempt_cap` was **not** raised to paper over
it; a parameter changed because a lane's own verification made it inconvenient is exactly the shape
this project keeps finding.

#### The three jobs that finished, and what they measured

**`k1` — `verdict: PASS`, 13.5 s.** `NVIDIA RTX A5000 Laptop GPU`, `sm_86`, warp 1.15.0,
`scatter_mode: ordered`. The whole enumeration on `cuda:0` (12.3 s) **and** on the warp CPU device
(0.3 s) in one process:

```
cases                       91   (j1: 49)      caught_on_cpu_but_not_cuda   []
families                    17   (j1: 10)      caught_on_cuda_but_not_cpu   []
wrong on purpose            74   (j1: 38)      true_kernel_regressions      []
caught on cpu               70                 self_comparison_nonzero      []
caught on cuda              70                 reasons                      []
blind on both  ['isotropic_drag[WRONG_AXIS_SWAPPED|TRUE]', 'isotropic_drag[WRONG_ISOTROPIC|TRUE]',
                'isotropic_drag[WRONG_UNWEIGHTED|TRUE]', 'tether_step[WRONG_NEWTON_MIRROR]']
```

**The four blind variants are blind on the CPU device too**, so none of them can produce a false
cross-device finding, and the CPU-side prediction in §7(b) reproduced on the device exactly.

**`k2` — `verdict: PASS`, 1.8 s. 312 gates graded across two devices and three position modes.
Nothing is outside its budget.** The two `-3611` §14.12 named as least likely to transfer:

| family | field | CPU ULP | **CUDA ULP** | budget | CUDA pressure |
|---|---|---|---|---|---|
| **M** `crossbridge` | `energy_start_pn_um` | 0.2454 | **0.2454** | 1.5 | 0.164 |
| **M** `crossbridge` | `forces_a_pn` | 0.1457 | **0.1457** | 1.5 | 0.097 |
| **M** `crossbridge` | `forces_b_pn` | 0.3141 | **0.3141** | 1.5 | 0.209 |
| **L** `strain_stiffening_cable` | `energy_start_pn_um` | 45.3414 | **45.3204** | 64 | **0.708** |
| **L** `strain_stiffening_cable` | `energy_end_pn_um` | 32.3346 | **30.8027** | 64 | 0.481 |

> **Family M's four gates are BIT-IDENTICAL across the two devices**, at the tightest budget in the
> port. **Family L's transcendental moved by 0.021 ULP on the gate that carries it** — the highest-
> pressure gate in the whole census on either device, and it moved *down*. The budget with the least
> claim to transfer transferred.

The largest device movement anywhere is `contact[TRUE].forces_cortex_pn` at **+6.65 ULP**
(16.79 → 23.44), which is a pressure of **0.011** against a 2048 budget — a large absolute move in a
channel with three orders of magnitude of room, and the reason the census reports both ratios.

**`k3` — `verdict: MEASURED`, 1.2 s. The anchor is lost identically on both devices.**

| device | sites at free velocity | force residual | power residual | restored |
|---|---|---|---|---|
| `cpu` | 1 / 12 | `1.421710e-07` pN | `5.669064e-08` | `False` |
| `cuda:0` | 1 / 12 | **`1.421710e-07`** pN | **`5.669064e-08`** | `False` |

Against a stall force of `3.3610` pN. **`-3611` §14.10's figure is a device-independent number here**
— CUDA's `wp.log` and its FMA contraction do not move it at this site. That closes a question the
entry could only mark `UNVERIFIED`, and it closes it in the direction that leaves §14.10's conflict
between §3 and §5 I1 exactly where it was. **Neither repair is taken and none is proposed.**

**`k4` is in flight** and its artefact does not exist at the time of writing. `-3609`'s twelve-
configuration sweep took ~21 min per pass over the grid, so six slack multiples should land near
7,600 s against a soft budget of 8,400 s and a runner timeout of 10,800 s. If that estimate is wrong
the soft budget stops the enumeration and `not_attempted` names every pair that was never run.

### 13a. Mutation study — 15 planted defects, 15 killed

Run with `PYTHONDONTWRITEBYTECODE=1` from a bytecode tree verified empty under `scripts/` and
`tests/scripts/`; baseline green at **43 passed** (was 25), each mutant applied alone and reverted
before the next, and the tree re-verified green at the end. Numbering continues `-3607`'s M1–M13.

| # | Planted defect | Killed by |
|---|---|---|
| M14 | the enumeration walks a hardcoded prefix of the registry, as `j1` did | `test_the_enumeration_covers_every_registered_law_case_builder`, `test_the_enumeration_grades_more_cases_than_the_first_queue_did` |
| M15 | wrongness is read from the primary variant axis only | `test_every_wrong_variant_is_wrong_on_at_least_one_axis` |
| M16 | an unknown variant axis is silently dropped and the primary graded alone | `test_an_unknown_variant_axis_refuses_rather_than_being_silently_dropped` |
| M17 | a registry disagreeing with its variant table is graded on the intersection | `test_the_enumeration_refuses_when_the_registry_and_its_variant_table_disagree` |
| M18 | the census omits the round-off bound and reports pressure alone | `test_the_budget_census_reports_pressure_and_the_roundoff_bound_for_every_graded_gate` |
| M19 | the census emits a suggested budget equal to what it measured | `test_the_census_does_not_widen_a_budget` |
| M20 | a gate outside its budget is reported as a green row | `test_the_budget_census_reports_a_gate_over_budget_as_a_finding` |
| M21 | the window sweep is configuration-major, so an early stop loses its breadth | `test_the_descent_window_rows_are_ordered_so_an_early_stop_still_covers_both_reference_multiples` |
| M22 | rows skipped by the soft budget are not recorded | `test_the_descent_window_job_records_rows_it_did_not_attempt` |
| M23 | the window job names a replacement for the ratified slack | `test_the_descent_window_job_proposes_no_tolerance` |
| M24 | the anchor probe measures every site instead of the anchor sites | `test_the_exact_anchor_probe_measures_both_devices_and_asserts_neither` |
| M25 | an empty anchor mask falls back to measuring every site | `test_the_exact_anchor_probe_is_not_vacuous` |
| M26 | the four second-queue jobs share one wall-clock budget | `test_the_second_queue_is_declared_with_its_own_budgets` |
| M27 | the second queue is given priority 0 and overtakes the first | `test_the_guard_still_refuses_all_three_ways_with_the_second_queue_armed`, `test_the_queue_runs_in_declared_priority_order`, `test_the_second_queue_runs_after_the_first_and_in_declared_order` |
| M28 | the preflight verdict is cached across the pass instead of re-run per job | `test_a_grant_that_expires_between_two_jobs_stops_the_second` |

**M25 is the one worth a sentence.** Its first form — deleting the guard outright — is killed by a
`ValueError` from `np.max` on an empty array, which is a crash rather than a wrong answer and is
therefore the *easy* version. The form planted here is the tempting repair: *if nothing sits at the
free velocity, measure every site*. That returns a plausible non-zero residual from a probe that is
no longer measuring an anchor, and it is only caught because the control drives a state whose mask
is genuinely empty. `ALEPH-PORT-3606` O2 in a probe rather than in a control: **a probe has no
assertion to notice that it measured the wrong set.**

**And one control went red before any mutant did**, which is the mutation study's cheapest result:
`-3607`'s own `test_the_queue_runs_in_declared_priority_order` named five jobs and failed the moment
four were appended. It now names all nine in order, so a tenth appended without a stated place is
red again. That is the control working, and it is recorded rather than quietly widened.

## 14. Honest limits

What this entry does **not** establish:

1. **The controls establish nothing about CUDA.** §8 and §9 are green on a Mac with no GPU, which is
   a statement about the runner and not about the port. What §13.1 records is a different kind of
   claim, and only `k1`, `k2` and `k3` support it; **`k4` had not finished when this entry was
   written and its section says so rather than leaving a gap a reader fills in.**

2. **Two jobs were reshaped to avoid editing `aleph/`, and each gave up a question.**
   - `k2` measures the three position modes at each builder's **default state only**, with no world
     offset. `j5` varied the offset for the ten older families because
     `scripts/cuda_jobs.py::law_cases` knows their state builders by name; the seven newer families'
     state builders are not reachable by any registry, and mapping them would mean either a second
     hand-written table — the exact defect T3 exists to remove — or a change in `aleph/`. **So the
     `GLOBAL_F32` vs `LOCAL_F32` distinction is measured at the origin, where it is weakest.**
     `-3610` and `-3611` already measured `PositionPrecision` *inert* for several of these families,
     which is why this was judged the cheaper loss; it is still a loss.
   - `k4` re-implements the descent-window driver in `scripts/` rather than calling
     `aleph.runtime.residency.job_descent_window`, whose configuration list is hard-coded to
     `[(1, 10.0), (2, 20.0)]`. The re-implementation uses the same public surface
     (`ResidentMembraneState`, `relax_to_equilibrium`'s two slack parameters) and the same
     monkeypatch shape, so it is a driver and not a second physics path — but it is **duplicated
     logic**, and if `-3609`'s owner changes the resident evaluate contract this copy will not
     follow. Named here rather than discovered later.

3. **A green table over 91 cases is not physics.** Every artefact carries
   `quantitative_status: BLOCKED`. A ULP figure is `STRUCTURAL` evidence about an implementation.
   Neither `k1` nor `k2` is evidence about a cell and no magnitude in either is sourced.

4. **`k2` cannot see a budget that is too loose.** Pressure and the round-off ratio both measure a
   budget against numbers the same lane derived. A family whose declared budget is generous will
   read comfortable on both, which is exactly how a generous budget reads. The only external check
   on a budget in this tree is the mutation study each ledger ran, and this queue does not repeat
   those.

5. **`k4` widens the sample and does not settle the question.** `-3609` measured the window on two
   configurations; `k4` measures twelve. Twelve icosphere configurations at one radius and one
   mobility is still one mesh family — `-3603` §14.3 measured a Fibonacci hull at the same vertex
   count needing 189,419 steps against 3,426 — so a window that holds across all twelve says
   nothing about another tessellation. And `k4` may not finish: its soft budget is real and its
   `not_attempted` list is expected to be non-empty.

6. **`k3` is one site.** The free-velocity mask in the default ratchet state has exactly one member.
   A probe on one site cannot distinguish "CUDA's `log` behaves differently here" from "CUDA's `log`
   behaves differently", and it does not try to.

7. **The reference-input guard is red before this entry and stays red.**
   `tests/runtime/test_reference_input_guard.py::test_every_registered_case_reference_is_sensitive_to_position_precision`
   fails on `isotropic_drag` at 0.443 ULP and, behind it, on `dense_exterior_resistance` at 0.262
   (`-3611` §14.7). That is `HANDOFF.md` §C-0's PI decision. **No job here trips it** — it is a
   *test* assertion about reference sensitivity, not `parity.refuse_narrowed_reference_input`, which
   refuses only a state carrying non-float64 float arrays and which none of these cases does.
   Verified by running the enumeration: 91 cases, **0 errors**. Nothing here was fixed or routed
   around.

8. **The refusal theorem is verified against the guard, not against the operating system**
   (`-3607` §14.3, inherited). And **occupancy is true for the instant it was read**
   (`-3607` §14.4): the preflight runs immediately before each job and the device can be taken
   between the reading and the launch. On a personal laptop the realistic contender is the PI's own
   desktop session.

9. **This lane did not verify the systemd unit survives a reboot**, unchanged from `-3607` §14.9.
   `Linger=yes` is set and `WantedBy=default.target` is what carries it; rebooting the PI's laptop
   to prove it is not a thing an unattended agent should do.

10. **`ports/ledger/INDEX.md` is not regenerated by this lane.** It has been red from other lanes'
    uncommitted work since before `-3607` was written, and the generator reads the working tree —
    which tonight carries the scenario lane's uncommitted files, whose author's process is gone.
    `-3613` joins that already-red assertion as one more *name*, not as a new red test.

11. **The timing figures are not a benchmark** (`-3607` §14.7, inherited). Wall-clock times land in
    the artefacts because a runner that cannot say how long a job took cannot be budgeted, and for
    no other reason. In particular `k1`'s 12.3 s on CUDA against 0.3 s on the CPU device is
    **kernel-module compilation and 91 per-case round trips**, not a throughput comparison, and
    nothing here supports a performance claim in either direction.

12. **`k4` has one attempt left, and this lane spent the other one.** Verifying that
    `systemctl --user disable --now` kills a *running* CUDA job required a running CUDA job, and
    `k4` was the one in flight. The attempt is recorded because it happened. `attempt_cap` was not
    raised: a bound loosened because a lane's own verification made it inconvenient is a bound that
    was never doing anything. If `k4`'s second attempt does not complete, the artefact set will show
    `k4_descent_window_wide.json` absent and `SUMMARY.json` will name it — which is the runner
    reporting a hole rather than a reader inferring one.

13. **A ULP figure that is bit-identical across two devices is a strong statement about a kernel and
    a weak one about CUDA.** `k2` found family M identical to the last bit on both devices and `k3`
    found the family-L anchor residual identical to the last bit. Both are one configuration on one
    A5000 in one session, at one warp version, and neither is a bound on what another architecture
    or another `wp.log` would do. `-3607` §14.4's rule about occupancy applies to arithmetic too:
    **it is true for the run that measured it.**

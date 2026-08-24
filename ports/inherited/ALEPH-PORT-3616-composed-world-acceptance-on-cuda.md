# ALEPH-PORT-3616 — the composed resident world's acceptance, on a CUDA device

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3616` |
| Lane | `46143f30` Lane G9 (CUDA track — operations). **`ALEPH-PORT-3613`'s queue reopened, not replaced.** |
| Status | `PROPOSED` |
| Written | `2026-08-03` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Why `PROPOSED` and not `ACCEPTED` after the code landed | Same reason as `-3607` and `-3613`, and it is a property of the discipline test rather than of this entry's completeness: `test_accepted_entries_name_both_controls` and `test_ledger_targets_that_claim_to_have_landed_actually_exist` both require an `ACCEPTED` entry to name an **`aleph/` target**, and an operations entry authorises none. The controls are named in §8 and §9 in the checkable form regardless, and `test_named_controls_resolve_to_real_tests` grades them. |
| Exists because | `ALEPH-PORT-3614` §14.7 says it plainly: *"Everything is `UNVERIFIED` on CUDA."* Its 30 acceptance controls, its 39-scalar transfer figure and its 0.0-ULP composed parity were all measured on **warp's CPU device**, because `gbook`'s A5000 was carrying `-3613`'s `k4` when that lane finished. `k4` has since completed and the runner is `inactive`; the grant runs to `2026-08-04T00:00:00+09:00`. This entry adds **one job** to `-3613`'s queue that drives `-3614`'s own controls on `cuda:0` and on the warp CPU device **in one process**, so the question is answered by comparison rather than by expectation. |
| **Number collision, recorded rather than merged** | This lane's brief assigned it `ALEPH-PORT-3615`. That id was **already taken** by the closed `codex-019fc511` filament-cortex row (`ALEPH-PORT-3615-membrane-envelope-and-filament-cortex-integration.md`, landed 13:48 KST the same day). `CLAUDE.md` §1: a claimed path is not written. The next free id was taken instead and the collision is reported to the PI rather than resolved by overwriting. |

---

## 1. Aleph API

This entry authorises **no change to any module under `aleph/`**, and it is stricter than
`-3613` on the point: `-3614`'s acceptance controls are driven **exactly as written**, through the
frozen public surface, with only the *device* redirected. Nothing is reimplemented, nothing is
relaxed, and no control is skipped for being awkward on a device.

```python
# EDITED — authorised by this entry, all within Lane G's own declared paths
scripts/cuda_jobs.py                      # + one job, k5, append-only
scripts/unattended_runner.py              # + one JobSpec; the queue machinery is untouched
tests/scripts/test_unattended_runner.py   # + controls; the priority-order control goes red first
```

```python
# CONSUMED, READ-ONLY — nothing here is edited by this entry
from aleph.runtime.resident_world import (
    FLOAT32_UNIT_ROUNDOFF, build_shell_pair_world, build_resident_world,
    _descent_predicate, _host_copy, recovered_tension_pn_per_um,
)
from aleph.runtime.residency import CountingBackend, TransferCounts
from aleph.runtime.backend import ScatterMode, WarpBackend, check_gpu_authorization
from aleph.runtime import law_cases, warp_kernels
```

and, uniquely in this file's history, **a `tests/` module is consumed as a registry**:

```python
tests/runtime/test_resident_world.py      # READ and DRIVEN, never edited
```

**Why the control module rather than a reimplementation.** The brief is *"the same 30 acceptance
controls, on the device"*. A job that re-expressed those 30 assertions in `scripts/` would be
grading a second implementation of the acceptance criteria, and the two would drift the moment a
lane added a 31st control — which is precisely the defect `-3613` §T3 repaired one level up when
`j1`'s hand-written family list went stale in six hours. Driving the module means **the enumeration
is a function of the module's own namespace**: a control added tomorrow is graded on the next device
run with no edit to `scripts/`, and a control *removed* is visible as a count that fell.

**`scripts/gpu_preflight.py` is consumed as a subprocess and is never imported, never edited and
never bypassed.** Unchanged from `-3607` §1 and `-3613` §1.

**`~/.aleph_data/gpu_authorization.json` is read and never written**, on either host, by any path in
this entry. `CLAUDE.md` §3 and `GPU_POLICY.md`: *an agent that writes its own permission slip has
not been authorized.*

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` — cited for the record only |
| Source path | `ac/fluid/transport_reference.py` (the pure-NumPy oracle beside its Warp kernel), and the reference project's device-vs-host parity habit generally |
| What was taken | **Nothing.** No line, no identifier, no constant. The one idea carried forward is structural and was already carried by `-3601`: *a kernel with no independent reference is unfalsifiable.* This entry adds the device-vs-device form of it. |

Nothing in this entry is a port of provider code. It is an **operations** entry, like `-3607` and
`-3613`, and its "source" is Aleph's own `-3614`.

## 3. Why source-derived porting beats clean-room

Not applicable in the usual sense — nothing is ported here. The rule that *does* apply is `-3613`
§3's: an operations entry earns its keep only if the thing it measures could have come back
negative. This one can. `-3614` §10.4 states in advance that it has **no reason to expect** its
figures to transfer, because `-3609` measured the descent window itself moving between the two
devices and a composed world has *more* float32 terms than a single owner, not fewer. A job whose
only possible outcome was "30/30 again" would not be worth the grant.

## 4. Physical or mathematical law represented

None directly. The object graded is an **implementation invariant** in three parts, all `-3614`'s:

* **I1 — residency.** A composed step moves no owner's positions, force field or accumulator across
  the host boundary. Only fixed-length scalars and `ceil(n/CHUNK)` reduction partials cross.
* **I2 — composition is a bit-identical refactor.** Every term of the composed world equals the same
  term driven through `law_cases`' own non-resident driver, **bitwise, per term, per step**.
* **The 30 acceptance controls** of `-3614` §8 and §9, which are the invariants above plus the
  refusals, taken together.

The question this entry asks is exactly one: **does each of those still hold when the device is
`cuda:0` instead of warp's CPU device?**

### 4a. Why the CPU arm runs in the same process, and is not quoted from the ledger

`-3614` §8 records 30/30 on the CPU device. Quoting that number and comparing a CUDA run against it
would compare two runs at two commits on two machines with two warp builds, and any disagreement
would have three candidate causes before it had one. `k1` established the shape used here: **drive
the identical enumeration on both devices in one process**, so the only difference between the two
arms is the device. A control that holds on CPU and fails on CUDA is then a finding about CUDA, and
a control that fails on both is a finding about the commit — and the artefact can tell them apart.

### 4b. Three ways this job could be silent, and what is done about each

| Hazard | Where it was found | What this entry does |
|---|---|---|
| A control that stashes a pulled array and compares it later compares a **buffer with itself** — a zero-copy view on warp CPU, a copy on CUDA, so the defect is *invisible on the device and live on the CPU arm*. | **G6** / `-3609` §14 | Every pull in this job's own code goes through `_host_copy`, and `test_the_device_pull_is_a_copy_and_not_a_view` drives the hazard directly by mutating the device buffer after the pull. |
| A count derived from a **hand-written list** cannot see what is not on the list. | **G7b** / `-3613` §T3 | The enumeration is `[fn for fn in vars(module) if name.startswith("test_")]`, in definition order, and a control whose signature requests a fixture this driver does not supply **refuses** rather than being skipped or partly graded. |
| A test that **stops at its first failure** does not tell you how many things are broken. | **G5c** / `-3611` §14 | The driver grades **every** control on both arms and never aborts on the first failure; the artefact carries a per-control row and the counts fall out of the rows. |

### 4c. What "the same 30 controls" means when the device changes underneath them

`-3614`'s controls take their backend from exactly three places: the module-level `_backend()`
helper, the `world` fixture (`build_shell_pair_world(level=1)`), and direct calls to
`build_shell_pair_world` / `build_resident_world` inside individual controls. The driver redirects
all three and **nothing else** — no assertion, no budget, no tolerance, no control body is touched.

**One control is expected to change meaning on the device, and it is named in advance rather than
discovered afterwards.** `test_a_coupling_across_two_backends_is_refused` builds a proxy that
*reports* the device label `"cuda:0"` while wrapping the suite's only real device, because the CPU
suite has no second device to build the configuration for real. On `cuda:0` the wrapped backend
already reports `cuda:0`, so the proxy no longer differs from its peer and the refusal it exists to
prove correctly does not fire. That is **V2's finding** — *a control can be correct, non-vacuous,
and never exercised at the boundary it defends* — reached from the other side: on the device the
control cannot be exercised at all, because the mock's discriminator collides with the real device's
name. It is reported as a **device-label collision, not as a defect in the refusal**, and it is not
worked around: the control is not edited, not skipped and not excluded from the count.

## 5. Units, domains, singular cases, invariants

| Quantity | Unit | Domain |
|---|---|---|
| transfer, bulk | elements (float32/float64 scalars) | `0` is the claim for a resident step |
| transfer, scalar | elements | fixed-length, must not grow with the mesh |
| transfer, reduction partials | elements | `ceil(n / CHUNK)`, `WarpBackend`'s own term, counted separately |
| per-field disagreement | ULP of the field scale, and raw `pN` | `0.0` is the claim for I2, and it is an identity rather than a budget |
| control outcome | `passed` / `failed` / `errored` / `skipped` | per control, per device |

**Singular cases.** A control that raises `pytest.skip` is recorded as `skipped` **with its reason**
and is counted in neither the held nor the failed column — a skip silently counted as a pass is the
`-3606` O2 vacuity hazard in a reporting channel. A control that raises anything other than
`AssertionError` is `errored` and is reported separately from `failed`, because "the world would not
build on this device" and "the assertion was false" are different findings.

**No relaxation window is searched on the device.** `k4` measured that no single descent slack
converges every configuration on CUDA — the coarse mesh needs `≥1.5·u32` and the fine mesh
`≤1.0·u32`, and the intersection is empty. That is an open PI decision. This job runs
`-3614`'s own `test_the_verified_relaxation_driver_can_drive_a_resident_world` exactly as written —
a **bounded** 40-step descent at the ratified `2·u32` that asserts the energy fell, not that anything
converged — and runs **no sweep**, proposes no tolerance, and reads no convergence figure.

## 6. Source evidence class and known retractions

1. `-3614` §10.1's **39 scalars / 30 partials at level 1** is a warp-CPU-device figure. Reproducing
   or contradicting it on CUDA is the point; neither outcome is a regression in itself.
2. `-3609` §10.2's **3 calls / 4 elements** is a *single-owner, no-transaction, no-ledger* figure and
   is **not** comparable with a composed one. Both are carried in the artefact with what they are.
3. `-3609` §14 measured `warp.array.numpy()` returning a **zero-copy view on the CPU device**. This
   job's CPU arm therefore runs under exactly the hazard the CUDA arm does not.
4. **Retracted before it is made:** any reading of this artefact as "the composed world is verified".
   It grades `build_shell_pair_world`, which is **not** `build_vertical` — `-3614`'s builder
   docstring lists the three constitutive terms it substitutes and why. No number here is a number
   about the verified vertical, and none is a number about a cell.
5. `-3613` §14.2's admissions about position modes are inherited unchanged; this entry measures at
   `POSITIONS_F64` only, which is what `-3614` runs at.

## 7. Independent oracle or derivation

Three, and they are independent of each other:

* **The CPU arm** is the oracle for the CUDA arm at the level of *outcomes*: same commit, same
  process, same enumeration, one device apart.
* **`law_cases`' non-resident drivers** are the oracle for the composed world at the level of
  *numbers*: they upload, allocate and read back independently, and I2 says the composed result must
  equal theirs **bitwise**. This is run on the device, both sides on the device.
* **A spy on `warp.array.numpy` and `wp.array.__init__`** is the oracle for the transfer count. It
  is not bookkeeping: it observes the real crossings, and it is checked for non-vacuity by a control
  that plants a bulk readback and requires the spy to see it.

## 8. Positive control

All in `tests/scripts/test_unattended_runner.py`.

| Control | Asserts |
|---|---|
| `test_the_acceptance_enumeration_is_derived_from_the_control_module` | §4b — the enumeration equals the `test_*` functions defined in `tests/runtime/test_resident_world.py`, computed from the module, and no control name is a literal in `scripts/cuda_jobs.py`. |
| `test_the_acceptance_driver_grades_every_control_even_after_one_fails` | §4b / G5c — a module whose first control fails still has every later control graded. |
| `test_the_acceptance_driver_distinguishes_a_failure_from_an_error_from_a_skip` | §5 — three outcome channels, and a skip is counted in neither column. |
| `test_the_device_pull_is_a_copy_and_not_a_view` | §4b / G6 — the pull is proven a copy by mutating the source buffer afterwards, on whatever device the suite has. |
| `test_the_transfer_audit_sees_a_planted_bulk_readback` | §7 — the spy is non-vacuous: a deliberate mesh-sized readback is counted. |
| `test_the_bitwise_channel_reports_a_planted_difference` | §7 — the bitwise comparison is non-vacuous: a perturbed reference is reported with its worst delta. |
| `test_the_composed_acceptance_job_compares_the_two_devices_in_both_directions` | §4a — the artefact carries `held_on_cpu_but_not_cuda` **and** `held_on_cuda_but_not_cpu`, so a one-directional reading is impossible. |
| `test_the_queue_runs_in_declared_priority_order` | **extended** — the declared order now names ten jobs, so appending an eleventh without saying where it belongs is red again. |
| `test_the_second_queue_is_declared_with_its_own_budgets` | **extended** — `k5` declares its own wall-clock and memory budgets and shares neither with another job. |
| `test_every_job_in_the_second_queue_has_a_dispatch_entry` | **extended** — `k5` is queued and implemented. |
| `test_every_job_declares_what_it_does_not_establish` | **inherited, unmodified** — `k5`'s limits are declared in `JOBS`, not in prose here. |
| `test_the_transfer_findings_name_a_mesh_sized_readback_and_a_silent_spy` | §7 / V2 — the transfer findings list is driven at the boundary it defends: a planted mesh-sized readback is named, a spy that saw nothing is named, and the unperturbed measurement is clean in the same control. |
| `test_the_acceptance_driver_restores_the_control_module_it_borrowed` | §1 — the module is `-3614`'s, is driven, and is handed back exactly as found, on the refusal path as well as the success path. |

In the checkable form the port-discipline parser reads:

| Positive | `tests/scripts/test_unattended_runner.py::test_the_acceptance_enumeration_is_derived_from_the_control_module` | The enumeration equals the control module's own `def test_*` set, and no control name reaches executable text in `scripts/cuda_jobs.py`. |
| Positive | `tests/scripts/test_unattended_runner.py::test_the_composed_acceptance_job_compares_the_two_devices_in_both_directions` | Both `held_on_cpu_but_not_cuda` and `held_on_cuda_but_not_cpu` are produced, and a control enumerated on one arm only is reported rather than dropped. |
| Positive | `tests/scripts/test_unattended_runner.py::test_the_moment_identity_probe_measures_and_can_report_a_difference` | §10.4's probe is exact on the CPU device before anything is planted, and reports magnitude *and* attribution when a perturbation is. |

## 9. Deliberately failing negative control

| Control (must fail) | Asserts |
|---|---|
| `test_an_unknown_fixture_axis_refuses_rather_than_being_partly_graded` | §4b / G7b — a control whose signature requests a fixture the driver does not supply raises `JobRefused`; it is **not** skipped, and the other controls are **not** graded around it. |
| `test_the_acceptance_job_searches_for_no_descent_tolerance` | §5 — an AST check over the new job: it assigns no `DESCENT_SLACK_ROUNDOFFS`, calls no sweep, and names no slack multiple. A job that tuned around `k4`'s open decision would destroy the evidence. |
| `test_the_acceptance_driver_refuses_a_module_with_no_controls` | §5 / `-3606` O2 — an empty enumeration is a refusal, not a clean 0/0 pass. A world with no controls grades perfectly. |
| `test_the_guard_still_refuses_all_three_ways_with_the_second_queue_armed` | **extended to ten jobs** — missing, expired and wrong-host records each give exit 2, zero launches and zero artefacts, with `k5` in the queue. |
| `test_the_bitwise_channel_refuses_to_compare_across_two_devices` | §14.6 — the subject and the oracle must be on one device, and the world must be built on the backend that was checked. **This is this job's own first-run bug, kept as a control.** |

In the checkable form the port-discipline parser reads:

| Negative (must fail) | `tests/scripts/test_unattended_runner.py::test_an_unknown_fixture_axis_refuses_rather_than_being_partly_graded` | A control asking for a fixture the driver cannot supply raises `JobRefused`; it is not skipped, and the rest are not graded around it. |
| Negative (must fail) | `tests/scripts/test_unattended_runner.py::test_the_acceptance_driver_refuses_a_module_with_no_controls` | An empty enumeration is a refusal, not a clean 0/0 pass. |
| Negative (must fail) | `tests/scripts/test_unattended_runner.py::test_the_acceptance_job_searches_for_no_descent_tolerance` | An AST check over the job's executable text: no descent slack, no sweep, no relaxation driver. |
| Negative (must fail) | `tests/scripts/test_unattended_runner.py::test_the_bitwise_channel_refuses_to_compare_across_two_devices` | A factory handing out two differently-labelled devices makes the bitwise channel refuse rather than measure. |

### 9.1 Mutation study

Planted one at a time, each reverted before the next, `PYTHONDONTWRITEBYTECODE=1`, second pass from
a **verified-empty** bytecode tree (`find aleph tests scripts -name '*.pyc' | wc -l` → 0).

> **23 planted, 23 killed. Three survived the first pass and all three were real gaps in the
> controls rather than equivalent mutants; the controls that close them are named below.**

| # | Mutant | Killed by |
|---|---|---|
| M1 | the enumeration drops the last control in the module | the enumeration control, and three others through the smaller set |
| M2 | the empty-enumeration refusal removed | the empty-module refusal control |
| M3 | a control asking for an unknown fixture is skipped instead of refusing | the unknown-fixture refusal control |
| M4 | the driver stops at the first failure | the grade-every-control control |
| M5 | a skip counted as a control that held | the three-channel control |
| M6 | an error reported as an assertion failure | the three-channel control |
| M7 | the borrowed module's `_backend` is left pointing at the device | the restoration control |
| M8 | the arm comparison is read in one direction only | the both-directions control |
| M9 | **the transfer spy never uninstalls itself** | the spy control — **after it was rewritten**, see below |
| M10 | the spy does not record anything crossing toward the device | the spy control |
| M11 | the pull returns warp's zero-copy view instead of a copy | the copy-not-a-view control |
| M12 | the cross-device refusal in the bitwise channel removed | the cross-device refusal control |
| M13 | the bitwise channel ignores a planted difference | the planted-difference control |
| M14 | the moment identity graded with a tolerance instead of as an identity | the moment-probe control |
| M15 | **the exact-arithmetic tally counts the differing components and attributes none** | the moment-probe control — **after it was rewritten**, see below |
| M16 | a control name that resolves to nothing is not refused | the restoration control |
| M17 | the held count starts at one | the grade-every-control and three-channel controls |
| M18 | the mesh-sized-readback check removed from the transfer findings | the transfer-findings control |
| M19 | a spy that saw nothing reported as a clean measurement | the transfer-findings control |
| M20 | the job names a descent slack multiple | the no-tolerance AST control |
| M21 | the queued job id matches no dispatch entry | five queue controls |
| M22 | `k5` shares `k4`'s wall-clock budget | the own-budgets control |
| M23 | **the bitwise subject built by the unpatched builder** | the planted-difference control — **after a guard was added**, see below |

**The three first-pass survivors, because they are the useful part of the study.**

* **M9 — a spy that never puts `wp.array` back survived every control, and the reason is that a
  leaked spy *nests*.** The second measurement's `original_numpy` is the first spy, so its counts
  are still exactly right — while every one of them also lands in the first, stale dictionary.
  **Nothing a measurement returns can see this**, which is why the control now asserts on the hooks
  installed on `wp.array` rather than on any number the spy produced. A new shape: *an instrument
  that is wrong about the world can be right about every measurement you ask it for.*
* **M15 — a tally that counts correctly and attributes nothing.** With `device_nearer` unable to
  increment, the `elif` never fires either and every differing component falls into `equidistant`,
  so the **total still matched** `elements_differing` and the control only checked the total. On the
  device this would have turned the finding of §10.4 — *the kernel is nearer to the exact value than
  its own frozen reference in 70 of 80 components* — into a flat "80 equidistant". The control now
  drives a perturbation whose direction is known and requires the attribution.
* **M23 — the subject built by the unpatched builder while the device check ran on a backend the
  world never used.** This is **the bug the first device run of this job actually had**, promoted to
  a mutant: a check on one object and a measurement on another. The guard is now that the world's
  backend must *be* the object whose device was checked, and `-3614`'s own controls are what caught
  it in the first place — see §14.6.

## 10. Numerical and precision envelope

`PositionPrecision.POSITIONS_F64` throughout, `ScatterMode.ORDERED` throughout — `-3614`'s settings,
unchanged, because changing them would make this a different measurement rather than the same one on
a different device.

All figures below are from a single run at clean `HEAD` `f7aa879`, `dirty: false`, host `GBook`,
`NVIDIA RTX A5000 Laptop GPU`, `sm_86`, warp 1.15.0, Python 3.12.13, numpy 2.5.1, under a preflight
that returned exit 0 immediately before the launch. Artefact:
`gbook:~/.aleph_runner/artefacts/k5_composed_world_device.json`, 3.0 s, `verdict: FAIL`.

### 10.1 The 30 acceptance controls on `cuda:0` and on the warp CPU device — MEASURED

The enumeration returned **30** controls, which is the count `-3614` §8 records, and the two arms
enumerated the same 30 — `enumerated_on_cuda_only` and `enumerated_on_cpu_only` both empty.

| arm | held | failed | errored | skipped |
|---|---|---|---|---|
| warp CPU device | **30** | 0 | 0 | 0 |
| **`cuda:0`** | **28** | **1** | **1** | 0 |

> **`held_on_cpu_but_not_cuda`: `test_the_moment_kernel_is_bit_identical_to_its_numpy_reference`,
> `test_a_coupling_across_two_backends_is_refused`.**
> `held_on_cuda_but_not_cpu`: empty. `not_held_on_either`: empty.

The second of the two was **named in §4c before the device was touched** and is a device-label
collision in the control's own mock, not a defect in the refusal. The first is §10.4 and is the
finding this lane exists to have produced.

### 10.2 What crosses the host boundary per composed step, on `cuda:0` — MEASURED

Spy-audited on `warp.array.numpy` and the `wp.array` constructor, eight steps per resolution, with
one warm step outside the counter.

| | level 1 (84 vertices) | level 2 (324 vertices) |
|---|---|---|
| host→device elements per step | **0** | **0** |
| device→host fixed-length scalar elements per step | **39** | **39** |
| device→host scalar calls per step | 15 | 15 |
| reduction partials per step (`ceil(n/CHUNK)`, `WarpBackend`'s own term) | 30 | 38 |
| **largest single readback** | **3 elements** | **3 elements** |
| a mesh-sized readback would be | 126 | 486 |
| findings | **none** | **none** |

> **Every figure is bit-for-bit what `-3614` §10.1 measured on the warp CPU device.** The count does
> not grow with the mesh — 39 at both resolutions across a 3.9× mesh — and the largest single
> readback is 3 elements at level 2, where a mesh-sized array would be 486. **The residency claim
> transfers to CUDA unchanged.**

The 39 is not zero and is not supposed to be: `-3614` §10.1 itemises it as two couplings × (2
centroids + 2 resultant forces + 2 resultant moments) × 3, plus the compartment's one-element energy
buffer read once per whole-world assembly, of which a step makes three. The claim is that every
crossing is fixed-length, and the discriminator for that is the **size** of the largest readback,
which is what the spy measures and what `CountingBackend`'s channel cannot distinguish.

Against `-3609` §10.2's single-owner figure of 3 calls / 4 elements: **not comparable**, and it is
carried in the artefact with what it is. A single owner has no connector books; the `AdjointPair`
the ledger reads is 12 numbers per coupling by construction. The comparable number is the bulk one
and it is **0** on both.

### 10.3 The composed world's bitwise parity against the non-resident drivers, on `cuda:0` — MEASURED

Six steps of a trajectory, five force fields and four energies per step, both sides on `cuda:0`
(`subject_device` and `oracle_device` are recorded in the artefact and are checked equal before the
world is built — §14.6).

| | field comparisons | energy comparisons | exact everywhere | worst \|Δ\| |
|---|---|---|---|---|
| warp CPU device | 30 | 24 | **yes** | 0.0 pN |
| **`cuda:0`** | **30** | **24** | **yes** | **0.0 pN** |

> **I2 transfers to CUDA as an identity, not as a budget.** Per owner, per coupling, per field, per
> energy, at every step. `-3614` §10.4 said in advance it had no reason to expect this, because a
> composed world has more float32 terms than a single owner and `-3609` measured the descent window
> itself moving between the devices. It held anyway, and that is a measurement rather than a
> reassurance.

### 10.4 The finding: **C7's identity does not hold on CUDA, and the kernel is the more accurate side**

`test_the_moment_kernel_is_bit_identical_to_its_numpy_reference` is `-3614` §8's C7. Its docstring
calls its budget *"exactly 0.0 ULP — not a tolerance, an identity"*, and §9.1 records that mutants
**M21** (the moment kernel flips one component) and **M22** (the arm kernel drops the `z` term) are
killed **by this gate and by nothing else.**

Measured on the device, at C7's own configuration (97 sites, seed 20260801, reference point
`(1.25, −0.75, 3.5)`):

| channel | components | differing | worst \|Δ\| | worst ULP (float64) |
|---|---|---|---|---|
| `moment_about` | 291 | **80** | 1.1369e-13 pN·µm | **36.00** |
| `arm_length` | 97 | **12** | 1.7764e-15 µm | **1.00** |

Both are exact on the warp CPU device in the same process.

**The 36 ULP is a cancellation figure, not an accuracy figure.** `ry·fz − rz·fy` is a difference of
two nearly equal products, so a half-ULP difference in one product is many ULP of the result. The
absolute disagreement is 1.1e-13 pN·µm against moments of order 1e2.

**And the direction was measured rather than assumed.** Each differing component was compared against
the exact rational value of the same expression, in `fractions.Fraction`:

| | components |
|---|---|
| **the device's answer is nearer the exact value** | **70** |
| the frozen NumPy reference is nearer | 10 |
| equidistant | 0 |

That is the signature of **FMA contraction**: CUDA fuses one of the two products into the subtraction,
so that product is rounded once instead of twice and the result is usually — but not always, because
the *other* product is still rounded — the correctly rounded one. Warp's CPU device compiles without
the contraction, so the host reference matches it bitwise there.

> **So C7 fails on CUDA because the kernel is, in 70 of 80 differing components, the *more* accurate
> of the two.** Its budget is an identity between an implementation and a compiler, not between an
> implementation and a law. **Nothing was changed and nothing is proposed**: `world_kernels.py` and
> C7 are `-3614`'s, this lane owns no byte of `aleph/**`, and which of the available repairs is right
> — a contracted reference, `--fmad=false`, or a budget that is a small ULP band rather than an
> identity — is a decision above this lane. §14.3 lists them without taking one.

## 11. Production-backend residency and transfer

Every buffer is allocated through the same `CountingBackend(WarpBackend(...))` path `-3614` uses;
this entry constructs the backend with `device=None` on the CUDA arm, which is the **guarded** path —
`WarpBackend.__post_init__` calls `require_gpu_authorization` before it imports warp, and it never
degrades to CPU. That is the second gate in series with the preflight, exactly as `-3607` §1
describes it: the preflight proves the machine is idle, and this one proves permission still exists
at the instant the device is acquired.

## 12. Comments and docstrings to discard

Nothing is copied from the reference project by this entry, so there is nothing to discard.

One statement in Aleph's own tree is read and **deliberately not repeated**: the header of
`tests/runtime/test_resident_world.py`, which says *"Every control here runs on warp's CPU device
… a test is not the place for a GPU run."* That remains true of the test suite and this entry does
not change it — the device run happens in `scripts/`, under the preflight, and the module is driven
from there rather than edited to know about a device.

## 13. Acceptance

| | |
|---|---|
| `tests/scripts` on this Mac | **57 passed** (was 43 before this entry) |
| `tests/scripts` on `gbook`'s own env | **57 passed** |
| mutation study | **23 planted, 23 killed**, second pass from a verified-empty bytecode tree |
| the device run | `k5_composed_world_device`, 3.0 s, exit 0, artefact complete, `repo.dirty: false` at `f7aa879` |
| `aleph/**` | **byte-for-byte untouched.** `git diff --stat` over this lane's commits names only `scripts/`, `tests/scripts/`, `ports/ledger/`, `docs/` |
| `~/.aleph_data/gpu_authorization.json` | **read, never written**, on either host; unchanged after the run |
| `tests/ports` | **19 passed, 2 failed** — both are `INDEX.md` staleness caused by *this* entry's own file, and `INDEX.md` was **not regenerated**, on instruction. `test_named_controls_resolve_to_real_tests` **passes**: every control this entry names is defined. |
| `tests/runtime` | not re-run by this lane; one failure there is foreign and known (`test_reference_input_guard.py`, `HANDOFF.md` §C-0) |

## 14. Honest limits

### 14.1 This grades `build_shell_pair_world`, not the vertical

Inherited verbatim from `-3614` §14.1. Three of the vertical's constitutive terms have no kernel at
all, so the world driven here substitutes for them. No number in this entry is a number about
`build_vertical`, and none is a number about a cell.

### 14.2 One device, one session, one warp build

An RTX A5000 Laptop GPU at `sm_86`, one CUDA runtime, one warp version, one night. `-3601` §14 says
what that is worth and it is not a bound on CUDA.

### 14.3 Three repairs are available for §10.4 and this entry takes none of them

Named so the PI is choosing from a list rather than from a blank page, and **none is recommended**:

1. **Contract the reference.** Express `numpy_moment_about` as a fused multiply-add so host and
   device round the same number of times. Python 3.13's `math.fma` would do it; this environment is
   3.12.
2. **Un-contract the kernel.** Build warp's CUDA modules with `--fmad=false`. That is a global
   compilation setting for the whole project, it would move every other kernel's numbers, and it
   buys bitwise agreement by making the device *less* accurate.
3. **Stop calling it an identity.** Replace C7's `array_equal` with a small ULP band. This is the
   only one that costs nothing elsewhere, and it is also the one that weakens the gate that
   `-3614` §9.1 records as the sole killer of M21 and M22 — so it needs a measured band, not a
   round number.

Each trades something different, and which trade is right is a decision about what a parity gate is
for. `CLAUDE.md` §2 rule 5: reported, not taken.

### 14.4 `test_a_coupling_across_two_backends_is_refused` cannot be exercised on the device

The control builds `_ReportsAnotherDevice(backend, "cuda:0")` — a real backend that *reports* a
different warp device — because the CPU suite has no second device to build the configuration for
real. On `cuda:0` the wrapped backend already reports `cuda:0`, so the mock no longer differs from
its peer, the refusal correctly does not fire, and the control errors with `DID NOT RAISE`.

**This is not a defect in the refusal**, which is `build_resident_world`'s and is untouched. It is
**V2's finding reached from the other side**: a control can be correct, non-vacuous, and *unable* to
reach the boundary it defends — here because its discriminator is a device name and the name
collided with the real device's. It was **predicted in §4c and committed before the device was
touched** (`e20e9b3`), so the prediction is checkable rather than retrospective.

It is not edited, not skipped, and not excluded from the denominator. **28/30 is the number.** The
repair, if the PI wants one, is a mock whose label cannot collide with any real device — and that is
an edit to `tests/runtime/`, which this lane does not own.

### 14.5 The acceptance controls are `-3614`'s, so this inherits every blind spot they have

Including all five of the series: an oracle and its subject sharing a float32 representation (G1); a
gate blind to a defect in the state it grades on (G4); blind to a defect in the reference it grades
against (G5); passing vacuously on an empty set (G5); and comparing a buffer with itself (G6).
Driving those controls on a second device widens the sample by one device. It does not repair one of
them, and this entry claims no such thing.

### 14.6 This job's own first device run was wrong, and `-3614`'s controls are what said so

The first run (16:39 KST, before `49c1faf`) reported the composed chain **7.6e-05 pN** away from the
non-resident drivers on `cuda:0`. It was not: `measure_composed_bitwise` took its subject from the
control module's *unpatched* builder and its oracle from the CUDA factory, so it compared a world on
warp's CPU device against a reference on `cuda:0` and reported the difference between two devices as
a difference between two implementations. **The number had exactly the order of magnitude a real
residency defect would have**, and nothing in the measurement announced the mismatch.

What caught it was that `-3614`'s own bitwise control — which puts both sides on one device —
**held** on the same run, on the same arm, in the same process. Two channels measuring one thing
disagreed, and the one this lane wrote was the wrong one. That is the argument for the two-arm
design turned on its author.

Recorded here rather than quietly fixed, because it is the fifth device-shaped instance of the
series and the first about the *harness*: **G1** oracle and subject share a representation · **G4**
a gate blind to a defect in the state it grades on · **G5** blind to one in the reference · **G5**
vacuous on an empty set · **G6** comparing a buffer with itself · **G9 — an oracle evaluated on a
different device from its subject, saying nothing.** The device identities are now compared before
the world is built, the world's backend must *be* the object that was checked, both are recorded in
the artefact, and M12 and M23 keep the two halves.

### 14.7 The two open PI decisions are untouched

`HANDOFF.md` §C (the descent tolerance, now with `k4`'s twelve-configuration measurement showing the
intersection is empty and the two mesh families pulling opposite ways) and §C-0 (`-3608`'s
reference-input guard). Neither is settled, worked around, or nudged here.

### 14.8 The number collision at the top of this entry is reported, not resolved

`ALEPH-PORT-3615` exists and belongs to another lane. This entry took `3616`. If the PI intended the
two to be one entry, that is a merge only the PI can authorise, and the alternative — overwriting a
landed lane's ledger — is the exact failure `CLAUDE.md` §1 exists to prevent.

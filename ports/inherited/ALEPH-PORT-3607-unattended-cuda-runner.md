# ALEPH-PORT-3607 — an unattended job runner whose defining behaviour is refusing to run

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3607` |
| Lane | `46143f30` Lane G7 (CUDA track — operations) |
| Status | `PROPOSED` |
| Written | `2026-08-01` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Exists because | G1–G5 put 52 kernels and a law-level parity harness in the tree and **every figure in all five ledgers was measured on the warp CPU device**, where a launch is a serial loop. `ALEPH-PORT-3603` §14.1 says so in as many words. The PI is away from ~09:00 Saturday 2026-08-01 and wants the device work to continue with **no operator connected**. A queue that runs unattended on a GPU is exactly the situation `docs/design/GPU_POLICY.md` was written for, so the deliverable is not "a runner" — it is a runner whose observable behaviour when permission is absent is to **stop**. |

---

## 1. Aleph API

This entry authorises **no change to any module under `aleph/`**. It adds three files under
`scripts/` and one test module, and it *consumes* the frozen surfaces below.

```python
# NEW — authorised by this entry
scripts/unattended_runner.py     # the driver: preflight gate, queue, artefacts, stop file
scripts/cuda_jobs.py             # the five jobs, each a pure `run(ctx) -> dict`
scripts/aleph-cuda-runner.service  # the systemd user unit (data, not code)
tests/scripts/test_unattended_runner.py
```

```python
# CONSUMED, READ-ONLY — nothing here is edited by this entry
from aleph.runtime.backend import ScatterMode, WarpBackend, check_gpu_authorization
from aleph.runtime.parity import (
    law_parity_report, law_report_as_dict, law_case_against_itself,
    measure_scatter_determinism, ulp_at_scale,
)
from aleph.runtime.law_kernels import (
    ArealVariant, BendingVariant, CentralContactVariant, ChromatinVariant, ContactVariant,
    CrosslinkVariant, PositionPrecision, SelfForceVariant, TetherStepVariant, TetherVariant,
    VolumeVariant,
)
from aleph.runtime.law_cases import (
    SWEEP_GRID, areal_law_case, bending_law_case, central_contact_law_case, chromatin_law_case,
    chromatin_owner_law_case, contact_law_case, crosslink_law_case, membrane_kernel_forces,
    membrane_state, recovered_tension_from_kernel_pn_per_um, tether_law_case,
    tether_step_law_case, volume_law_case,
)
from aleph.vertical.assembly import build_vertical, recovered_tension_pn_per_um
from aleph.vertical.relax import relax_to_equilibrium
```

**`scripts/gpu_preflight.py` is consumed as a subprocess and is never imported, never edited and
never bypassed.** It is a guard; `CLAUDE.md` §2 rule 4 forbids modifying it, and this entry's whole
value rests on it being the unmodified one.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — cited for the record only |
| Source path / symbol | **none named, none read.** No file under that tree was opened by this lane. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable — nothing was read |

Nothing was ported. The artefacts this lane read are all Aleph's own:
`docs/design/GPU_POLICY.md`, `scripts/gpu_preflight.py`, `aleph/runtime/backend.py`'s
authorization half, and `ALEPH-PORT-3601` / `-3603` / `-3604` / `-3605` / `-3606`.

## 3. Why source-derived porting beats clean-room

**It does not, and this entry says so rather than manufacturing a reason.** Clean-room wins
outright: a job queue is a convention, and the only part of this work with any subtlety —
*what a runner must do when its permission check refuses* — is answered by Aleph's own written
policy and by an incident record in this repository, not by anything a provider wrote. Porting a
scheduler would import somebody else's answer to the one question that matters here.

The entry is kept, and it is `RE-DERIVED`, because `PLAN.md` §0.2.5 wants the *decision* recorded,
not only the code.

## 4. Physical or mathematical law represented

**There is no physical law here, and naming one would be the over-claim this project exists to
catch.** What crosses the boundary is a *safety invariant*, stated as a theorem about the runner's
reachable states:

> **T1 (the refusal theorem).** For every job `j` in the queue, the runner enters the state
> `LAUNCHED(j)` only along an edge whose guard is `preflight(j).exit_code == 0`, where `preflight`
> is a subprocess execution of the unmodified `scripts/gpu_preflight.py` begun after the previous
> job terminated. There is no other in-edge to `LAUNCHED`. Consequently: absence of a grant,
> expiry of a grant, a grant naming another host, a grant naming another device, an unreadable
> grant, a busy device, insufficient free memory, and an unreadable machine **all** map to the same
> successor state, `WAITING`, and `WAITING` has no out-edge to `LAUNCHED` that does not pass the
> guard again.

Two corollaries this lane treats as load-bearing rather than as commentary:

* **T1a.** The guard is re-evaluated **per job**, not per process. A grant can expire mid-queue and
  by design *will*: the queue is asked to outlive the current grant. A runner that checked once at
  start would be running unauthorized work by the time it mattered most.
* **T1b.** The guard's verdict is read from its **exit code**, never from its text.
  `PLAN.md` §12 records why: a refusal recognised by message is a refusal that stops being
  recognised when somebody improves the wording, and the failure is silent in the direction that
  runs the job.

The second statement this entry makes is an *anti*-law, and it is the reason the file list in §1 is
short:

> **T2.** The runner contains no code path — none, including a fallback, a repair, a "first run"
> bootstrap or an error handler — that creates, edits, extends or moves
> `$ALEPH_DATA_ROOT/gpu_authorization.json`. `docs/design/GPU_POLICY.md`: *"an agent that writes its
> own permission slip has not been authorized"*. §9 makes this checkable rather than promised.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| job wall-clock budget | s | s | `(0, 10800]`, declared per job, never shared |
| declared device memory budget | GB | 1e9 B | `(0, 16]` on the A5000; passed straight to the preflight |
| poll interval while refused | s | s | `[60, 3600]` |
| queue deadline | ISO-8601 instant | s | an explicit argument; no default that outlives the PI's absence |
| measured disagreement | float32 ULP of the field scale | dimensionless | `[0, inf]`; `inf` is a real value here — `ulp_at_scale` returns it for a kernel that invents force out of an identically-zero reference |
| recovered tension | pN/µm | 1e-12 N / 1e-6 m | reported only with its `converged` flag beside it |

Singular and boundary cases, each with the behaviour the runner requires:

- **No authorization file.** Log the refusal with its exit code, sleep, re-check. Never run.
- **Expired authorization.** Identical treatment. Expiry is the *expected* end state of this queue.
- **Authorization naming another host or another device.** Identical treatment.
- **Device busy / not enough free memory / machine unreadable** (preflight exit 3 / 4 / 5). Identical
  treatment. The runner does not distinguish "refused because unsafe" from "refused because
  unauthorized" when deciding whether to launch, because both answers are *no*; it distinguishes
  them only in what it writes to the log.
- **A job that hangs.** Killed at its own declared timeout, together with its process group, and
  recorded as `timeout`. A GPU job that will not die is worse than one that never started.
- **A job that crashes.** Retried at most once. A second failure marks it `failed` and the queue
  moves on — a crashed job that is retried forever is the spin this runner must not become.
- **An artefact that already exists and is complete.** Skipped without launching anything, so a
  reboot or a `systemctl restart` resumes rather than recomputes.
- **A non-converged relaxation.** Reported as a **refusal with the row's numbers withheld from the
  headline**, following `scripts/sweep_vertical_parameters.py`: a table cannot show you that one of
  its rows stopped early, because the number is the same shape either way.
- **`inf` or `nan` in a measured field.** Written to the artefact as the string `Infinity`/`NaN` via
  a JSON encoder that refuses to silently drop them, because a missing field reads as "not measured"
  and an infinite ULP is a measurement.

Invariants, each with the test that asserts it:

- **I1.** No job is launched when the preflight's exit code is non-zero —
  `tests/scripts/test_unattended_runner.py::test_an_expired_authorization_stops_the_runner_and_launches_nothing`,
  `::test_a_citation_for_another_card_stops_the_runner_and_launches_nothing`,
  `::test_a_missing_authorization_stops_the_runner_and_launches_nothing`.
- **I2.** The preflight is re-run before **every** job, not once per process —
  `::test_the_preflight_runs_once_per_job_and_not_once_per_process`.
- **I3.** The verdict is taken from the exit code, not the message —
  `::test_the_refusal_is_read_from_the_exit_code_and_not_from_the_message`.
- **I4.** A completed artefact is not recomputed —
  `::test_a_completed_artefact_is_not_recomputed`.
- **I5.** A hung job is killed at its declared timeout —
  `::test_a_job_that_hangs_is_killed_at_its_declared_timeout`.
- **I6.** A repeatedly failing job is abandoned, not retried forever —
  `::test_a_job_that_crashes_is_abandoned_after_the_attempt_cap`.
- **I7.** The stop file halts the queue before any preflight or launch —
  `::test_the_stop_file_halts_the_queue_before_any_preflight`.
- **I8.** No source path in the runner writes an authorization record —
  `::test_the_runner_has_no_path_that_writes_an_authorization_record`.
- **I9.** Every artefact names the commit, the device, the scatter mode, the precision mode and a
  `quantitative_status` — `::test_the_artefact_names_what_it_was_measured_on`.
- **I10.** The queue is ordered by declared priority and a later job cannot overtake an earlier one —
  `::test_the_queue_runs_in_declared_priority_order`.

### 5a. The authorization record the PI must write on `gbook`

Recorded here because §9's negative controls are written against this exact shape, and because an
agent may never create it. **Path on `gbook`: `~/.aleph_data/gpu_authorization.json`.**

```json
{
  "granted_by": "PI",
  "granted_at": "2026-08-01T02:00:00+09:00",
  "expires_at": "2026-08-03T09:00:00+09:00",
  "host": "GBook",
  "devices": [0],
  "note": "verbatim text of what the PI actually said"
}
```

**`"host"` must be exactly `GBook`**, with that capitalisation, and this is a measured fact rather
than a preference. Three things must agree on the string and they do only at `GBook`:

1. `socket.gethostname()` on the machine returns `GBook`, which is what
   `aleph/runtime/backend.py::_resolve_host` falls back to;
2. `~/.ssh/config` there carries `Host GBook` → `HostName localhost`, so
   `scripts/gpu_preflight.py --host GBook` can reach `nvidia-smi` through the ssh path it uses for
   every host. **`ssh gbook` (lower case) does not match that stanza**, resolves out to the tailnet
   address instead, and fails host-key verification — so a record naming `gbook` would make the
   preflight return exit 5 (`machine unreadable`) forever, which is safe and useless;
3. the record already on that machine names `GBook`, so the PI is not being asked to learn a new
   spelling — only to extend an expiry.

The record presently on `gbook` names `GBook` and **expired at `2026-07-31T12:22:09Z`
(21:22 KST)**. The runner therefore ships in its refusal state and stays there until the PI writes
the file above.

## 6. Source evidence class and known retractions

Nothing was read from the provider, so there is no source evidence class to report. The *Aleph*
artefacts this entry rests on carry these statuses, and each is inherited explicitly:

| Artefact | What it claims | Status carried here |
|---|---|---|
| `ALEPH-PORT-3601` §14.7 | the true tether kernel disagrees with the frozen law by ~100 ULP because of float32 **position storage** | `INHERITED_UNVERIFIED` on CUDA — measured on the warp CPU device |
| `ALEPH-PORT-3603` §13.4 | Laplace σ recovered at worst `1.5610e-05`, 12/12 converged, through the kernels | `INHERITED_UNVERIFIED` on CUDA — this runner's job `j4` exists to test exactly that |
| `ALEPH-PORT-3603` §14.6 | a float32 **energy** channel converges 0/12 as a descent criterion | inherited as a **design constraint**: `j4` reads the descent energy in float64 |
| `ALEPH-PORT-3604` §9a | `WRONG_NEWTON_MIRROR` is a declared blind spot, not a killable mutant | inherited: `j1` records it as `blind_spot`, never as `caught` |
| `ALEPH-PORT-3605` §10.5 | the rest nucleus scores a *correct* kernel `inf` because the reference scale is exactly zero | inherited: `j1` must not read `inf` as a failure without reading the reference scale beside it |
| `docs/design/GPU_STATE_2026-07-31.md`, `636b0c8` | A5000 `ORDERED` 0.00 ULP / `ATOMIC` 10.00 ULP over 16 launches | **one device, one session.** `j2` re-measures rather than quotes |

Retractions looked for and found, in this repository rather than the provider's: the 21:40 reading
of `gpu-sungwook` as "usable and completely idle" is **retracted** in the build plan's anchor block,
twice over. That is the direct precedent for this entry's §14.1 — an occupancy or capability reading
is true for the minute it was taken.

## 7. Independent oracle or derivation

Three, and they are of different kinds on purpose.

**(a) For the runner itself: the guard's own exit codes.** `scripts/gpu_preflight.py` documents
`0 / 2 / 3 / 4 / 5` in its module docstring and `tests/runtime/test_gpu_authorization.py` pins every
refusal reason of the in-process half. The oracle for T1 is therefore not "the runner looks safe" —
it is *drive the unmodified guard with a record that must be refused and observe that nothing
launched*. That is a behavioural oracle, and it is why §9's controls execute the real script as a
subprocess instead of stubbing it.

**(b) For the physics jobs: the CPU figures the five kernel ledgers already published.** The CUDA
run is graded against numbers that were fixed before it ran, by another lane, on another device:

| # | Job | The figure it is measured against | Where |
|---|---|---|---|
| j1 | full law parity on `cuda:0` | every declared `DEFAULT_*_ULP_BUDGET`, and the CPU ULP table of each ledger | `-3601` §13, `-3603` §13.1, `-3604` §13, `-3605` §13 |
| j1 | the wrong kernels | every `WRONG_*` variant must exceed its budget in **at least one** field, except the one declared blind spot | `-3601` §9, `-3603` §9, `-3604` §9/§9a, `-3605` §9 |
| j2 | scatter determinism | `ORDERED` 0.00 ULP; `ATOMIC` 10.00 ULP against a CPU permutation bound of ≤ 5.30 | `GPU_STATE_2026-07-31.md` |
| j3 | the virial identity through kernels | worst `3.5253e-08` over 24 configurations | `-3603` §13.4(a) |
| j4 | the relaxation reproduction | worst `1.5610e-05`, **12/12 converged**, against the published `1.5663e-05` | `-3603` §13.4(c), `docs/results/2026-07-31-tension-sweep/README.md` |
| j5 | the three position modes | `GLOBAL_F32` / `LOCAL_F32` / `POSITIONS_F64` per law, including `-3605` §13.4's counterexample where `POSITIONS_F64` is **worse** (11.37 vs 2.57 ULP) | `-3603` §13.2, `-3605` §13.4 |

**(c) For j3 and j4 specifically: an exact identity, not an agreement.** At `c0 = 0` the discrete
membrane energy splits as `E(λx) = E_bend(x) + λ² σ A(x)`, so differentiating along the dilation
`dx_i = x_i` at `λ = 1` gives `Σ_i (grad E)_i · x_i = 2 σ A`, hence

    σ = −(Σ_i f_i · x_i) / (2 A)

with `f = −grad E`. This holds at **any** configuration, not only at equilibrium, which is what
makes it affordable as a device gate: no relaxation is needed for j3. It is an
`ANALYTIC_ORACLE` — an identity of the model, not a measurement of a membrane.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/scripts/test_unattended_runner.py::test_the_runner_launches_the_job_when_the_preflight_passes` | With an injected preflight returning exit 0, exactly one launch happens, the artefact is written, and the recorded `preflight.exit_code` in that artefact is `0`. A runner that never launches would pass every negative control in §9 and be worthless; this is the control that makes those meaningful. |
| Positive | `::test_a_completed_artefact_is_not_recomputed` | A second pass over a queue whose artefacts exist launches **zero** subprocesses and reports every job `skipped`. I4. |
| Positive | `::test_the_artefact_names_what_it_was_measured_on` | Every artefact carries `repo.commit`, `device`, `scatter_mode`, `precision_modes`, `preflight.exit_code`, `authorization.expires_at` and `quantitative_status`. I9. |
| Positive | `::test_the_queue_runs_in_declared_priority_order` | The five jobs are attempted in the PI's stated order. I10. |
| Positive | `::test_the_preflight_runs_once_per_job_and_not_once_per_process` | Three jobs produce three preflight invocations. I2. |

## 9. Deliberately failing negative control

The three that matter are the same experiment run three ways, and each drives the **real**
`scripts/gpu_preflight.py` against a planted record in a temporary `ALEPH_DATA_ROOT`.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/scripts/test_unattended_runner.py::test_an_expired_authorization_stops_the_runner_and_launches_nothing` | Record valid in every field but `expires_at` in the past → preflight **exit code 2**, launcher called **0** times, outcome `refused`, no artefact written. |
| Negative (must fail) | `::test_a_citation_for_another_card_stops_the_runner_and_launches_nothing` | Citation live but naming another **card** → **exit 2**, 0 launches. Was *another host* until `ALEPH-PORT-3618`; the host is no longer in the record. |
| Negative (must fail) | `::test_a_missing_authorization_stops_the_runner_and_launches_nothing` | No file at all → **exit 2**, 0 launches. Absence is not permission. |
| *(retired)* | *the name is removed rather than struck through, because a struck-through name still reads as a citation and the discipline test still resolves it* | **RETIRED 2026-08-04, `ALEPH-PORT-3618`.** A second negative control asserted that a grant for device 0 refuses a request for device 1. The device index moved to the allocation — `gpu-submit` sets `CUDA_VISIBLE_DEVICES` and the runtime refuses if more than one appears — so it collapsed onto the card check above and became the same assertion from the other side. Two tests that cannot fail independently are one test and a maintenance cost. |
| Negative (must fail) | `::test_the_refusal_is_read_from_the_exit_code_and_not_from_the_message` | A stub preflight that exits **non-zero while printing `[preflight] PASS`** still stops the runner. This is the mutant `PLAN.md` §12 names: a runner that matched on the word `PASS` would launch here. |
| Negative (must fail) | `::test_a_grant_that_expires_between_two_jobs_stops_the_second` | Job 1 runs under a valid record; the record is then expired in place; job 2 is refused. T1a, and the case this queue is *expected* to hit. |
| Negative (must fail) | `::test_a_job_that_hangs_is_killed_at_its_declared_timeout` | A job that sleeps past its budget is killed with its process group and recorded `timeout`; the queue continues. I5. |
| Negative (must fail) | `::test_a_job_that_crashes_is_abandoned_after_the_attempt_cap` | A job that always exits non-zero is attempted at most `attempt_cap` times across passes and then reported `failed`. I6. |
| Negative (must fail) | `::test_the_stop_file_halts_the_queue_before_any_preflight` | With the stop file present, **zero** preflights and **zero** launches. I7. |
| Negative (must fail) | `::test_the_runner_has_no_path_that_writes_an_authorization_record` | The source of `scripts/unattended_runner.py` and `scripts/cuda_jobs.py` contains no reference to `write_authorization_record` and no write to `gpu_authorization.json`; and a runner run against a read-only data root leaves the record absent. T2 / I8. |
| Negative (must fail) | `::test_a_wrong_kernel_is_still_over_budget_when_the_gate_is_asked` | The j1 verdict function reports `caught` for a field whose measured ULP exceeds its effective budget and `MISSED` when it does not — so a device run in which the wrong kernels stopped being caught is reported as a failure rather than as a green table. |

**What the negative controls deliberately do not do.** They do not stub the preflight for the three
authorization cases. A test that asserts "the runner stops when *my mock* says no" grades the mock.

## 10. Numerical and precision envelope

The runner performs no arithmetic on physical quantities; it carries them. What it *fixes* is the
context every figure is only meaningful inside, and every one of these travels into the artefact:

* **Compute channel:** float32 (`ALEPH-DQ-107`), which is `WarpBackend.dtype`'s default and is not
  changed by any job. Accounting reductions are float64 by construction of the backend.
* **Position representation:** a declared `PositionPrecision`, defaulting to `GLOBAL_F32` for j1–j4
  so the CUDA numbers are comparable with G1–G5's CPU numbers at the same setting. j5 varies it
  across all three modes and **answers nothing** — the open question is
  `docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` and it is the PI's.
* **Scatter mode:** `ORDERED` for every parity and physics job. `ATOMIC` appears in j2 only, where
  measuring it *is* the job. `ATOMIC_SCATTER_ULP_FLOOR = 10.0` widens a scatter-assembled field's
  budget only when the backend is in atomic mode, which no physics job here is.
* **Descent energy in j4:** float64, from the frozen NumPy energy. Not a preference —
  `ALEPH-PORT-3603` §14.6 measured that the kernels' float32-term energy converges **0 of 12** and
  degrades the recovered σ by 331×, because one float32 ULP of a 2,929 pN·µm energy is
  3.5e-04 pN·µm and `MonotonePotentialDescent` compares two energies.
* **Tolerances:** none are invented by this entry. Every budget is the `DEFAULT_*_ULP_BUDGET` its own
  ledger derived, and j1 reports `measured / effective_budget` as a pressure ratio so a figure that
  is inside a *loose* budget cannot read as a tight result.
* **Outside the envelope:** a job whose measured value is `nan` writes `nan` and is marked
  `unreadable`, never coerced. A relaxation that does not converge is a refusal, not a row.

## 11. Production-backend residency and transfer

| Stage | Host | Device | Transfer |
|---|---|---|---|
| queue state, artefacts, log | `gbook` disk | — | none |
| preflight | `gbook` host process → `ssh GBook nvidia-smi` | reads device metadata only | no allocation |
| law-parity reference | host, float64 NumPy | — | never uploaded |
| law-parity candidate | host state → device arrays each case | `cuda:0` | one upload per case, one download per field; **per case, not per step** |
| j4 relaxation | host owns positions; the membrane force is computed on device | `cuda:0` | **one round-trip per relaxation step**, and this is a known inefficiency rather than a design — it is what makes the reproduction literal. G6 (`CellWorld` step on device) is the task that removes it and is out of scope here |

The runner holds no device memory of its own and allocates none: every array it touches comes from
`backend.zeros`/`backend.array` inside the frozen case code, which `ALEPH-PORT-3601` §1 already
requires (a raw host ndarray as a scatter destination works on warp's CPU device and fails on CUDA —
the exact defect `278e6e5` reported and `636b0c8` retracted). This runner is the first thing in the
project that will find out whether that requirement was met everywhere.

## 12. Comments and docstrings to discard

Nothing was read, so nothing is inherited and nothing needs discarding. What replaces the absent
source prose is the material in §4: the refusal theorem, its two corollaries, and the anti-law T2 —
each written in Aleph's vocabulary and each attached to a named control in §5.

One phrase is **banned from this lane's artefacts** and the ban is worth stating because it is the
easy thing to write at 4 a.m.: no artefact may say a job was "skipped for safety" or "deferred".
The word is `refused`, the exit code is beside it, and the reason is the guard's own text.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending, and the runner is live in its refusal state.** See §13.1 for what was observed on `gbook` itself. Acceptance requires the PI to write §5a's record and a subsequent artefact set carrying `preflight.exit_code == 0`. |
| Reviewer | **Agent-proposed, unratified.** No PI decision is recorded here and no `decided_by` field exists in this file. |
| Rollback | `systemctl --user disable --now aleph-cuda-runner.service`, delete `~/.config/systemd/user/aleph-cuda-runner.service`, and `git revert` the commit that added the three `scripts/` files. Nothing under `aleph/` changes, so no physics, no test and no gate depends on this entry. What breaks: the CUDA figures stop being produced, and every claim in G1–G5 stays `INHERITED_UNVERIFIED` on device — which is exactly where it is today. |

### 13.1 What was observed on `gbook`, at `e10ec3a`

Not asserted — run, on that machine, and quoted from its own output.

| Observation | Result |
|---|---|
| checkout updated from `e6c2f95` to `e10ec3a` via a fresh bundle | working tree clean before and after, `0` modified paths, `0` stashes, fast-forward only |
| `tests/scripts` on `gbook`'s own `aleph` env | **25 passed** in 10.2 s |
| `scripts/gpu_preflight.py --host GBook --devices 0 --memory-gb 2` | `[authorization] REFUSED: authorization expired at 2026-07-31T12:22:09.500541+00:00.` — **exit 2** |
| one live pass of the runner against that record | `REFUSED j1_law_parity: preflight exit 2 (NO VALID AUTHORIZATION ON RECORD). Launching nothing.` |
| artefacts written | **0** |
| compute processes on the A5000 afterwards | **0** |
| the stop line, executed and then reverted | service `inactive`, `disabled`, `0` runner processes surviving |
| service state now | `active`, `enabled`, `Linger=yes`, `MemoryMax=12G`, `TasksMax=256`, deadline `2026-08-03T09:00:00+09:00` |

### 13a. Mutation study — 13 planted defects, 13 killed

Run with `PYTHONDONTWRITEBYTECODE=1`; baseline green, each mutant applied alone and reverted before
the next.

| # | Planted defect | Killed by |
|---|---|---|
| M1 | the verdict refuses only on the busy code (`exit_code != 3`) | the three authorization controls |
| M2 | the verdict is read from the guard's **text** (`"PASS" in stdout`) | `test_the_refusal_is_read_from_the_exit_code_and_not_from_the_message` |
| M3 | the preflight is cached after the first job of a pass | `test_a_grant_that_expires_between_two_jobs_stops_the_second` |
| M4 | an artefact that merely exists counts as complete | `test_a_half_written_artefact_is_not_counted_as_complete` |
| M5 | the attempt cap is never reached | `test_a_job_that_crashes_is_abandoned_after_the_attempt_cap` |
| M6 | the stop file is not checked inside a pass | `test_the_stop_file_halts_the_queue_before_any_preflight` |
| M7 | a timed-out job is reported `ok` | `test_a_job_that_hangs_is_killed_at_its_declared_timeout` |
| M8 | a guard that cannot be executed is treated as a pass | `test_a_guard_that_cannot_be_executed_is_treated_as_a_refusal` |
| M9 | the job is launched **before** the verdict is read | every refusal control |
| M10 | the queue runs in reverse priority order | `test_the_queue_runs_in_declared_priority_order` |
| M11 | the artefact promotes its own `quantitative_status` to `QUOTABLE` | `test_the_artefact_names_what_it_was_measured_on` |
| M12 | the deadline is never past | `test_the_deadline_halts_the_queue` |
| M13 | a refusal falls through to the next job instead of stopping the pass | the three authorization controls |

**M8 is the one worth a sentence, because it survived its first form.** A *missing guard script* and
a *missing interpreter* take different paths: the first makes the interpreter exit non-zero, which
every control already caught; the second raises `OSError` inside the runner and reaches an error
handler that no control reached. A handler returning `0` there — "the guard could not be run, so
proceed" — passed all 25 tests. The control now drives both, and this is the second time in this
repository that a mutation study found an unreached branch rather than confirming a written one.

## 14. Honest limits

What this entry does **not** establish:

1. **It establishes nothing about CUDA by existing.** Until the PI writes §5a's record, every
   number in G1–G5 remains a warp-CPU number and this runner has produced no measurement at all.
   The runner's *own* controls are green on a machine with no GPU, which is a statement about the
   runner and not about the port.
2. **A green parity table on `cuda:0` is not physics.** Every artefact carries
   `quantitative_status: BLOCKED`. A ULP figure is `STRUCTURAL` evidence about an implementation;
   the Laplace recovery is `ANALYTIC_ORACLE`, an identity of the model. Neither is evidence about a
   membrane, and no magnitude here is sourced.
3. **The refusal theorem is verified against the guard, not against the operating system.** T1 says
   the runner never launches without a passing preflight. It does not say a GPU cannot be reached by
   some *other* process on that laptop, and it cannot: `WarpBackend`'s own
   `require_gpu_authorization` is the in-process half and this runner adds a second gate in front of
   it, but neither can constrain a program nobody in this repository wrote.
4. **Occupancy is true for the instant it was read.** The preflight is run immediately before each
   job and the device can be taken between that reading and the launch. The preflight says so itself.
   The window is small and it is not zero, and on a personal laptop the realistic contender is the
   PI's own desktop session.
5. **`ORDERED` determinism on CUDA is expected, not proven, before j2 runs.**
   `measure_scatter_determinism`'s caveat field already says a CPU zero is a property of a serial
   schedule. j2 is what converts the expectation into an observation, and if `ORDERED` turns out
   *not* to be bit-identical on this device, every budget in G1–G5 is affected and that is a finding
   to route upward, not to absorb.
6. **j4 is one tessellation.** Every configuration is an icosphere, as
   `ALEPH-PORT-3603` §14.3 already records: a Fibonacci hull at the same vertex count needed
   **189,419** relaxation steps against 3,426, so a step budget that is generous here says nothing
   about another mesh family. A `converged=False` row is reported as a refusal for that reason.
7. **The timing figures are not a benchmark.** Wall-clock times land in the artefacts because a
   runner that cannot say how long a job took cannot be budgeted, and for no other reason. Nothing
   here measures throughput, compares CPU against CUDA fairly, or supports a performance claim —
   the j4 path round-trips to the host every step and is deliberately the slow shape.
8. **Two jobs may not fit the window.** j4's budget is 2 h of wall clock and j5's is 1 h; if j4
   consumes its budget the queue may reach Sunday with j5 unattempted. That is the stated priority
   order doing its job, and the artefact set will say which jobs never ran rather than leaving a
   reader to infer it from silence.
9. **This lane did not verify the systemd unit survives a reboot on `gbook`.** `Linger=yes` is set
   for the account, which is the precondition, and `WantedBy=default.target` is what carries it —
   but rebooting the PI's laptop to prove it is not a thing an unattended agent should do.
10. **`ports/ledger/INDEX.md` is not regenerated by this lane.** It has been red from other lanes'
    uncommitted work since before this entry was written, and the generator reads the working tree.
    `-3607` joins that already-red assertion as one more *name*, not as a new red test.

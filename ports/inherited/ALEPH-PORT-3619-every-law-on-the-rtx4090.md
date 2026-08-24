# ALEPH-PORT-3619 — every registered law, every variant, on the RTX 4090

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3619` |
| Lane | `46143f30` — CUDA track, operations |
| Status | `PROPOSED` |
| Written | `2026-08-04` |
| Port class | `RE-DERIVED` |
| Why `PROPOSED` and not `ACCEPTED` | Same reason as `-3607`, `-3613` and `-3616`, and it is a property of the discipline test rather than of this entry: `test_accepted_entries_name_both_controls` requires an `ACCEPTED` entry to name an `aleph/` target, and an operations entry authorises none. The controls are named in §8 and §9 in the checkable form regardless. |
| Device | NVIDIA GeForce RTX 4090 (`4090-0`), Warp 1.15.0, Slurm jobs 18–20 |
| Host | `sungwook@100.110.26.26` — `DESKTOP-EEUU4RQ`, WSL2 |
| Artefacts | `docs/results/2026-08-04-rtx4090-sweep/{sweep_cuda,sweep_cpu}.json` |
| Exists because | Every kernel here had been graded on CPU-Warp and on an RTX A5000 that **no longer exists**. The PI retired `gbook` on 2026-08-04. That retires nothing about its measurements and everything about the claim *"this runs on a GPU"* — which was resting on a machine nobody can rerun. |

---

## 1. Aleph API

**This entry authorises no change to any module under `aleph/`.** It adds no kernel, no law case
and no binding. It runs the kernels that already exist and records what they returned.

The one file it adds outside `docs/` is the sweep driver, and it lives in the run directory on the
workstation rather than in the repository, because a script whose only purpose is to produce an
artefact that is committed is not part of the engine. §8 names how to reproduce it.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository identity | **none** — `/Users/sw1/ffn_cellsim` is not read by this entry |
| Source path | **none.** The inputs are internal: `aleph/runtime/law_cases.py::LAW_CASE_BUILDERS` (17 builders) and each builder's variant enum in `aleph/runtime/law_kernels.py`. The oracle is each case's own `reference`, already ported under `-3601` … `-3617`. |

**No source is ported here and none is read.** `/Users/sw1/ffn_cellsim` is untouched by this entry.

Its inputs are entirely internal: `aleph.runtime.law_cases.LAW_CASE_BUILDERS` (17 builders) and
each builder's variant enum in `aleph.runtime.law_kernels`. The oracle is each case's own
`reference`, which is the frozen NumPy law already ported under `-3601` … `-3617`.

## 3. Why source-derived porting beats clean-room

Not applicable in the usual sense — nothing is derived from a source here. The question this entry
does answer is the adjacent one: **why measure both devices in one job rather than compare against
a number recorded earlier.**

A CPU figure taken on another day on another machine differs from a device figure for reasons that
include the day and the machine. Running both in the same process on the same state removes every
variable except the one under study, which is the arithmetic. That is why `sweep_cpu.json` exists
and why it was produced by job 19 rather than quoted from `-3614`.

## 4. Physical or mathematical law represented

None of its own. It represents **the union of every law already registered** — the ERM tether, the
membrane's areal and bending response, the unilateral connectors, contact, volume, chromatin,
crosslink, central contact, block-diagonal and isotropic drag, dense resistance, the series joint,
the strain-stiffening cable, the Brownian ratchet, and the NMII crossbridge.

73 cases, **261 channels**, at every variant of every builder.

**Variants are read off the enum, never a hand-written list.** A hand list cannot see a variant
somebody adds later, which is the error class this track has now hit six times and which G7b
repaired in the mechanism rather than in a list.

## 5. Units, domains, singular cases, invariants

The unit is **float32 ULP of each field's own reference scale**, which is the harness's unit
throughout and is not re-derived here.

The domain is every registered case at its own default state. **Singular cases are the point of
§9**: three of the four uncaught mutants are uncaught because the *state* is degenerate for them,
not because the gate is weak.

The invariant this entry tests is not a physical one. It is:

> **the verdict is a property of the law, not of the device.**

261 channels, two devices, and the verdict must agree on all of them.

## 6. Source evidence class and known retractions

**Evidence class: MEASUREMENT.** Not a simulation and not a physics claim — a comparison of two
implementations of one law, on one state, on two devices.

**Retraction carried forward:** `-3614` §14.7 said *"Everything is UNVERIFIED on CUDA."* `-3616`
answered that for the composed world on the A5000. **That answer is now unrepeatable** — the A5000
is retired — and this entry replaces it on hardware that exists. `-3616`'s numbers stand as what
was measured; they are no longer reproducible.

## 7. Independent oracle or derivation

Each case's own `reference` — the frozen NumPy law, **called and never re-implemented**. A second
NumPy implementation written for this sweep would compare the kernels against a copy of themselves
made by the same hand, which is how a parity gate becomes ceremony.

The independent check on top of that is the **cross-device comparison**, which no single-device
oracle can provide: the two runs share a reference and disagree only through the hardware.

## 8. Positive control

`tests/runtime/` in full, **inside the allocation**, job 20:

```
647 passed, 1 failed, 2 skipped
```

The one failure is `test_every_registered_case_reference_is_sensitive_to_position_precision`,
unchanged, a PI decision (`HANDOFF.md` §C-0). **Identical to the local result**, which is the
control that the device changed nothing it should not have.

The sweep's own positive control is the true-variant row: **0 of 60 true channels outside budget**,
on both devices.

Reproduce:

```bash
rsync -az --exclude .git --exclude __pycache__ aleph/ tests/ scripts/ \
  sungwook@100.110.26.26:~/aleph_gpu_run/
ssh sungwook@100.110.26.26 \
  'MEM_GB=12 CPUS_PER_TASK=8 gpu-submit 4090-0 01:00:00 bash ~/aleph_gpu_run/run_sweep.sh'
```

## 9. Deliberately failing negative control

> ### ⚠ CORRECTED 2026-08-04 14:2x — this section originally reported a problem that does not exist
>
> The first version of this entry said four wrong kernels were **uncaught** and called one of them
> *"not explained"*. **All four are declared blind spots**, written out in the enum docstrings in
> `law_kernels.py` before this sweep existed. The sweep counted them as findings because **a sweep
> cannot read prose**, and I read the count rather than the declarations.
>
> The corrected figure is **zero undeclared blind spots**. The rest of this section is the repair.

**56 wrong kernels. 52 caught by at least one channel. 4 bit-identical to their true kernel — and
every one of the four was already declared.**

| pair | why it cannot differ | caught elsewhere |
|---|---|---|
| `tether_step_law_case` / `WRONG_NEWTON_MIRROR` | `force_b := -force_a` against accumulating from its own unit vector. For a central force with a scalar magnitude these are bit-identical at **every** input: `|b-a| == |a-b|` because components are squared, `(-d)/L == -(d/L)` because negation is exact. **Not state-dependent** — no configuration separates them. | nowhere, and nowhere can |
| `isotropic_drag_law_case` / `WRONG_ISOTROPIC` | the anisotropic excess deleted, on a state where `ratio = 1` makes it exactly `0.0` | `anisotropic_drag_law_case` |
| `isotropic_drag_law_case` / `WRONG_AXIS_SWAPPED` | the two axes carry the same coefficient at `ratio = 1`, so swapping is the identity | `anisotropic_drag_law_case` |
| `isotropic_drag_law_case` / `WRONG_UNWEIGHTED` | the length weight dropped, on a state where `L = 1` | `anisotropic_drag_law_case` |

**The first is structural and the other three are the state's doing**, and that distinction is why
the declaration is now keyed by **`(builder, variant)`** rather than by variant. The same
`bdiag_drag_wrong_isotropic` kernel is invisible on the isotropic case and caught on the anisotropic
one — one kernel, one defect, two verdicts. A declaration keyed by variant alone would silence both,
and the second is the one doing the grading. **That flaw was in my first attempt at the table and
the control caught it within a minute of being written.**

**30 wrong-kernel channels are bit-identical to the true kernel** on both devices. That is the blind
*channel* set, which is a different and larger thing than the four blind *cases*: a mutant can be
invisible in four of its five outputs and loud in the fifth, and it is caught. This is why the sweep
reports per channel — a per-case boolean prints `False` once and the shape is gone.

### The mechanism, repaired rather than noted

`law_kernels.DECLARED_BLIND_VARIANTS` now carries the declarations as data with the reason attached,
and `tests/runtime/test_declared_blind_spots.py` makes the declaration **cost something**:

- `test_every_declared_blind_spot_is_actually_blind` — a row that is not truly bit-identical turns
  red. You cannot declare your way out of a real defect, because the control measures the claim.
- `test_no_undeclared_variant_is_blind` — a blind pair absent from the table turns red. **This is
  the control that would have caught a real problem**, and the sweep alone could not.
- `test_the_state_dependent_ones_are_caught_somewhere_else` — without it, *"declared blind"* would
  be indistinguishable from *"never graded anywhere"*, which is the failure that matters.

## 10. Numerical and precision envelope

| | RTX 4090 | CPU-Warp |
|---|---|---|
| true channels outside budget | **0 / 60** | 0 / 60 |
| wrong kernels caught | 52 / 56 | 52 / 56 |
| channels bit-identical to the true kernel | 30 | 30 |
| errors | 0 | 0 |
| wall clock | 6.4 s | 4.5 s |

**Channels whose ULP differs between devices: 116 of 261 (44 %).
Channels whose verdict differs: 0.**

Both halves are load-bearing:

- **If all 261 had agreed bitwise, the correct conclusion would be that the kernels did not run on
  the device.** Reduction order differs between a 4090's warps and a CPU loop; identical output
  would mean something other than agreement.
- **None of the 116 disagreements moved a verdict.** Device variation lives inside the ULP budgets
  and every mutant's defect lives outside them. **The two magnitudes do not overlap** — the property
  a budget exists to have, and this is the first sweep wide enough to state it.

Worst true channel on the device: `contact_law_case/forces_membrane_pn` at **535 ULP** against a
budget of 2048 — a factor of four of headroom, and the tightest law in the tree.

## 11. Production-backend residency and transfer

Unchanged by this entry and not measured by it. The sweep constructs a `WarpBackend` per device and
drives `compare_law_case`, which is the harness's own transfer path. Residency figures remain
`-3609`'s and `-3614`'s.

## 12. Comments and docstrings to discard

None. No source prose is read, copied or adapted — see §2.

## 13. Acceptance

This entry is accepted when a reader can rerun §8 and obtain the §10 table. It authorises no
`aleph/` change, so it names no acceptance target under `aleph/`.

**Two defects it produced, both fixed at `ee1572a`,** and both are one fact wearing three hats —
**exit code `2` means three different things**: python's `can't open file`, argparse's missing
required argument, and these scripts' refusal.

`test_the_refusal_matches_the_preflight_script_exit_code` compared the two scripts' exit codes and
nothing else, so it **passed while neither script ran** — first because `scripts/` was not synced,
then because `--host` went missing when `--card` was added. It now asserts both produced output, and
that assertion caught the second bug within minutes of being written.

The citation tests' skip predicate asked whether a CUDA **driver exists**. Wrong question: on this
host libcuda loads and `cuInit` still returns `100` outside an allocation. It now asks whether
`cuInit` **succeeds** — the honest condition, because **a subprocess can fake `SLURM_JOB_ID` and
`CUDA_VISIBLE_DEVICES` and cannot fake the driver.** That asymmetry is what `-3618` rests on, so a
test that could be forced past it would be testing something else.

## 13a. A finding that came from cancelling a job, not from running one

Jobs 18–20 measured the laws. **Job 21 measured the harness, by failing to finish**, and the number
is worth more than the run would have been.

`residency sweep` and `resident_world sweep` were dispatched on the RTX 4090 to re-measure the
descent window that `ALEPH-PORT-3608` and `k4` established on the retired A5000. After 40 minutes:

```
GPU 0 utilisation   12 %
CPU                110 %      one core saturated
```

**The sweep is host-bound, not device-bound.** The kernels compute forces; the descent criterion,
the step control and the convergence decision are all Python. So the whole premise of re-running it
on faster hardware was wrong — **the faster resource is not the one the job uses**, and an A5000
figure and a 4090 figure would differ by the host, not by the card.

The arithmetic that should have been done before submitting:

| | configurations | slack multiples | points |
|---|---|---|---|
| `k4_descent_window_wide` (A5000, **8,513 s**) | 12 | 6 | 72 |
| job 21's `residency sweep` (default `--multiples`) | 12 | **8** | **96** |

1.33× the work, on a machine whose relevant resource is comparable — roughly **3 h 10 min**, against
a 1 h 30 min wall clock. It could not have finished, and the three cheap stages queued behind it
would never have started. `transfers` completes in seconds; it was scheduled after the expensive
stage and lost to it.

**The larger point is about the shared machine.** A host-bound sweep holding a 4090 for three hours
at 12 % utilisation is a card five other users cannot have, spent on work that does not need one.
Reducing `--multiples` treats the symptom. The real question is whether this sweep should reserve a
GPU at all, and it is recorded here rather than answered.

**What was salvaged:** `residency transfers` on `cuda:0` completed before the cancellation and its
artefact is kept — 4 elements per resident step, 3 readback calls, and the same
`uncounted_by_construction` lower bound the CPU run reports.

**Retraction of my own framing.** I told the PI this would be "the A5000 numbers, re-measured on
hardware that exists". The A5000 numbers for the *sweep* were never device numbers in the sense that
mattered; they were host numbers taken while a device was attached. `-3609`'s and `-3614`'s
**transfer** figures are genuinely device-side and those are the ones this entry re-measures.

## 14. Honest limits

- **No claim is made that 52/56 is the right catch rate.** It is the measured one, and the four
  that are not caught are declared and controlled rather than open.
- **`WRONG_NEWTON_MIRROR` is not verifiable by any gate and that is a derivation, not a gap.** It
  ships so the claim is checkable and its control asserts the bit-identity rather than a failure.
- **One device, one card, one state per case.** A different state could move any of these numbers,
  and three of the four uncaught mutants are direct evidence that states differ in what they can see.
- **The four were not counted together before this sweep.** Earlier lanes declared each where it
  arose; this is the first count over all 73. **The number that matters is the undeclared count,
  which is zero**, and the first version of this entry did not compute it because the declarations
  were prose. That is now a test rather than a habit.
- **`4090-1` and `3090` are unmeasured.** The PI's citation names `4090-0` only.
- **The descent sweeps are NOT re-measured on the 4090** and §13a says why. `k4`'s twelve-
  configuration result on the A5000 stands as the descent-window evidence, and re-running it on a
  different card would not move it, because the quantity is host arithmetic.

## Authorization

```
"ㅇㅇ 4090 일단은 6시간 동안 쓰는 것으로 할까?"
  — PI, 2026-08-04T04:29:19.883Z
  ~/.claude/projects/-Users-sw1-Project-Aleph/46143f30-….jsonl#99c6e06a-1701-47ee-9aa9-d1ac82421c6f
```

Card `4090-0`, six hours. Jobs 18, 19 and 20 fall inside that window. **No authorization record was
written by an agent** — the file above is a citation with coordinates, not a grant (`ALEPH-PORT-3618`
§3).

# The production environment on the GPU host cannot run its own test suite, and never has

**2026-08-21, ~06:00 KST.** Found trying to execute the first CUDA-gated test this project has ever
had a reason to run on the device.

## 1. What is on the host

| env | warp | pytest |
|---|---|---|
| **`ffn_sim`** — the production runtime; every physics number ever measured here | **1.14.0** | ⚠ **absent** |
| `aleph` | 1.15.0 | 9.1.1 |

And: **there is no `.pytest_cache` anywhere on the host.** No pytest process has ever run there.

## 2. What that means, stated plainly

A test gated on CUDA **skips on the development machine** (no device, correctly) and **cannot be
invoked in the production environment on the run host** (no runner). So:

> ⚠ **Every CUDA-gated test in this repository has, until today, never executed anywhere.**

The local suite reports them as `skipped`, which reads as *"not applicable on this machine"*. The true
reading is *"not verified on any machine"*. **A test that is skipped everywhere is not a test**, and
the suite's own summary line is what hides it — it is the ninth member of tonight's family and the
one with the widest blast radius, because it is not a defect in one check but a whole *class* of
checks that has never run.

### 2.1 ⚠ Scope — the measured number, replacing my own estimate

I first wrote *"40 test files contain a CUDA gate and hold 437 test functions"*, hedged as not being
437 never-executed tests. **It is measured now, and it is far smaller.** The full suite with `-rs`
reports **12 skips**, of which **8 require a CUDA device**:

| skip | count |
|---|---:|
| `test_arena_device.py` (4 gates) | 4 |
| `test_observe_gamma.py` (the two device gates) | 2 |
| `test_microtubule.py:74` — *"I0-A: kernel gates require a CUDA GPU"* | 1 |
| `test_actin_motor_joint.py:181` — *"end-to-end woven-joint bind gate requires a CUDA GPU"* | 1 |
| **CUDA-gated total** | **8** |
| Cytosim binary absent | 3 |
| a gap that closed | 1 |

**Eight, not 437.** The file-level count was a proxy for a thing I had not measured, and it read as
larger than the truth by fifty-fold. The structural claim is unchanged and is what matters: **those
eight had never executed anywhere** — they skip on the dev machine for want of a device and could not
be invoked in production on the host for want of a runner.

⚠ And two of the eight are gone as of `7959aba5`: the `observe_gamma` pair was **deleted and replaced
by a module self-check** reachable as `python -m aleph.world.observe_gamma`, precisely so a gate does
not live where it cannot be run. **Six remain.**

## 3. The one that ran today, and exactly what it is worth

`aleph/tests/world/test_observe_gamma.py`, holding the two gates the device γ port needs:

* `test_device_matches_host` — the equivalence gate
* `test_the_device_reduction_is_bitwise_reproducible` — the degeneracy check's device-side precondition

Run on the GPU host, on the card, inside the Slurm hold: **26 passed, 1 skipped.** The single skip is
the *inverse* gate (`test_the_device_path_refuses_a_non_cuda_device`), whose reason reads *"this
machine has CUDA; the refusal is only reachable without one"* — so **both device gates executed, and
both passed.**

⚠ **But it ran in the `aleph` environment, on warp 1.15.0, and production is warp 1.14.0.** So:

> **The kernel is no longer unmeasured. It is measured on a NEIGHBOURING RUNTIME, not on the one every
> other number in this project came from.**

That is a real advance over *"authored and unmeasured"* and it is **not** the gate. For an equivalence
test between a device kernel and a host estimator, the device compiler is precisely the variable under
test, and a different Warp version is a different compiler. **It may not be recorded as the CUDA gate
passing, and it is not being recorded that way.**

## 4. Why the obvious fix was not taken tonight

Installing pytest into `ffn_sim` would close this in one command. It was not done, because:

* the GPU host is **shared with other people**, and its production environment is what every physics
  number in this project was measured in;
* a `pip install` into a conda environment resolves dependencies and may move something else, and
  **the class of defect that produces is exactly tonight's** — a number whose runtime nobody can
  reconstruct afterwards;
* a τ run with a committed pre-run contract is **executing in that environment right now**.

⚠ It is **PI queue item 15**, with the recommendation to install a test runner into `ffn_sim`
deliberately, pinned, and recorded — rather than to keep running the gates in a neighbouring
environment and describing the result carefully every time.

## 5. Method note — this was found the same way as everything else tonight

Not by auditing the environment. By **trying to use it**: the port's author asked for one pytest
invocation on the card, and the invocation failed with `No module named pytest`. The environment had
been the run host for the entire project and no one had ever asked it for this.

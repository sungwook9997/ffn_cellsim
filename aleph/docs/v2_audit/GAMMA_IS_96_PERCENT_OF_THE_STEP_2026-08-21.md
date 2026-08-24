# 96% of the step is a host-side observable, and it is what makes the τ run cost eleven hours

**2026-08-21, ~05:45 KST. Written while the long τ run is still executing**, deliberately: this
conditions how §5 of [`TAU_RUN_CONTRACT_2026-08-21.md`](TAU_RUN_CONTRACT_2026-08-21.md) may be read,
and an interpretation written after the data arrives is an interpretation fitted to it.

## 1. What was noticed, and how

Checking the run was alive: **`nvidia-smi` reported 0–1% utilisation on the card**, sampled three
times, while the process burned 552% CPU across 5.5 cores. For a 4.5-million-node Warp run at 4 s/step
that is the wrong shape.

⚠ **Three samples was too thin to publish and a later sample read 70%, so it was re-measured** — twenty
samples at 1 s: `0 0 0 0 54 0 2 1 0 0 0 41 1 0 13 0 0 4 67 0`. **Median 0%, fourteen of twenty at or
below 1%, mean 9.2%**, with occasional bursts. The expected duty cycle from the step split below is
`0.1677 / 4.1212` = **4.1%**, and a 1-second sampling window over a 4-second step catches the burst
sometimes and misses it usually. Same order, consistent shape. **The utilisation was the clue; the
arithmetic in §2 is the evidence, and it does not depend on the sampling at all.**

## 2. The measurement

| | s/step |
|---|---:|
| cortex bound, **no** γ (5 runs, 0.1513–0.1677) | **0.1677** |
| cortex bound, **γ emitted** | **4.1212** |
| **γ costs** | **3.954 s/step = 95.9% of the step** |

## 3. Where it goes — and it is NOT the roundtrip

`world_phase4_native.py:296` calls, **every step**:

```python
pos=pos_d.numpy(), f_ext=force_d.numpy()
```

Two full device→host copies of `(4,361,496, 3)` float64 — **104.7 MB each, 209.4 MB per step.**

⚠ **But the transfer is not the cost.** 209.4 MB in 3.954 s is **53 MB/s effective**, against a PCIe
gen4 ×16 ceiling around 25,000 MB/s. **The copy accounts for about 0.21% of the γ cost.** The other
99.8% is `gamma_from_resting_readback` running as host NumPy over 4.36 M nodes.

So the fix is not "avoid the roundtrip". **The estimator has to be a kernel.**

## 4. The contract this sits against

`CLAUDE.md` §Hard runtime contract:

> **GPU-only inside the loop.** No per-step host snapshot or host state mutation, no hidden fallback.
> A profiler gate must show zero authoritative GPU→CPU roundtrips inside the physical-time loop; host
> readback belongs BETWEEN accepted steps — **and is not free there either.**

⚠ **This is a per-step host snapshot inside the physical-time loop.** Whether it is *authoritative* is
arguable — γ is an observable and does not feed back into the dynamics, so the trajectory is unchanged
by it. What is not arguable is the second sentence: the readback is inside the loop, not between
accepted steps, and *"not free there either"* turns out to be an understatement by a factor of 24.

⚠ **And no profiler gate is running on this path.** The contract names one; nothing enforces it here,
which is why a 96% host-side step went unremarked through every γ run tonight — including the ones
whose 8.34 and 4.06 s/step numbers were quoted in three documents as if they were the engine's cost.

## 5. What this does to the τ contract's second branch

The contract states both outcomes in advance. The second is:

> τ keeps tracking N → required window ≈ 58,600 steps ≈ **66 h/seed, 198 h for three** → *this
> configuration does not produce a stationary γ on this hardware.*

| | with host γ | with a device-side γ |
|---|---:|---:|
| this run, 9,000 steps | 10.30 h | **0.42 h** |
| the 58,600-step branch | 67.08 h | **2.73 h** |

⚠ **So the second branch does not mean what it appears to mean.** *"This configuration cannot produce
a stationary γ on this hardware"* is **false as stated**; the true statement is **"cannot produce one
with a host-side observable"**, and the same series at 0.168 s/step is **under three hours for one
seed and about eight for three.**

**This does not change the run, the contract, or the resolution criterion** — the trajectory is
identical either way and τ in steps is unaffected. It changes **only** what may be concluded from
branch 2, and it is recorded now, before the numbers exist, so that it cannot be mistaken for a
rescue written after an unwelcome result.

⚠ It also **strengthens** the third option the lane offered the PI (γ kernel moves device-side): that
option is not one of three comparable alternatives. It is a 24× throughput change on the exact axis
the τ question is budget-limited by.

## 6. What is NOT claimed

⚠ No γ magnitude, unchanged — coefficients are test points and the mobility is unsourced.
⚠ No claim that a device-side γ **reproduces** the host estimator; that is a port with its own
verification, and its cost is an estimate from the no-γ step time, not a measurement of code that
exists.
⚠ No claim about the other observables. This measures γ on the resting path, on one configuration.

---

## 7. The optimisation would have disabled a correctness gate — caught before it was written

The owner of `observe_gamma.py` took the device port and, **before writing it**, found the coupling
that the obvious implementation carries. It is worth more than the port.

**The naive kernel reduces per-plane force sums with `wp.atomic_add`.** Atomic addition does not fix
the summation order, so the order varies run to run, so **the result varies run to run in the last
bits** — the identical mechanism behind the 2.111e-13 residual floor measured earlier tonight.

⚠ **And that manufactured noise lands exactly on a check that exists to catch a false pass.** At 02:24
this session closed a defect where a *constant* γ series certified as `STATIONARY, sem 0.0,
n_eff 6000`. The guard that closes it compares the series' `std` against `eps·|mean|`. **Round-off
noise from a non-deterministic reduction pushes `std` above that bound**, the degeneracy check stops
firing, and the false pass reopens — **silently, via a performance change, with every test still
green.**

> **A 24× speedup that switches off a correctness gate is not a speedup.**

⚠ **And the tempting fix is the forbidden one.** Widening the degeneracy threshold to accommodate the
new noise floor is choosing a threshold to fit an implementation, which `CLAUDE.md` prohibits in the
first hard rule. The answer taken instead is to **remove the non-determinism**: a two-stage reduction
writing `(plane, block)` partial sums into a fixed-size array, with a second kernel summing them in
**fixed index order**. No atomics, bit-identical across runs. Cost: ~512 KB of scratch.

**This is the eighth member of tonight's family and the only one that never happened.** The other
seven were found by using an instrument that was already lying. This one was found by asking, before
writing a line, *what does this optimisation touch that is not performance* — and the answer was a
gate two levels away, closed three hours earlier.

⚠ **And the author's own correction to that framing is the useful part.** They pointed out that they
asked the question because **they had closed that degeneracy check themselves**, and that had it been
someone else's finding they would have missed it too — *"there is more luck in this than method."*
That is right, and it shrinks what transfers. *"Ask what your optimisation touches that is not
performance"* only worked here because the answer was already in one person's head. **What actually
generalises is a handoff, not a habit: whoever optimises a quantity must be told which gates that
quantity feeds, because they cannot be expected to have personally closed them all.**

⚠ **What the port can and cannot carry when it lands.** The authoring machine has no CUDA, and
`CLAUDE.md` forbids CPU execution for smoke tests as firmly as for production. So the kernel can be
authored, type-checked and code-generated there, and a CUDA-gated equivalence test can be written
against `plane_force_sums` — but **the kernel is UNMEASURED until it runs on the device**, and its
commit must say so rather than inherit the confidence of a green suite.

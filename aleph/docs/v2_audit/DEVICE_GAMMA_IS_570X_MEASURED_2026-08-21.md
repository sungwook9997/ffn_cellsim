# The device γ is 570× the host path — measured at native, on the granted card, inside the allocation

**2026-08-21, ~06:35 KST.** This exists because the estimate it replaces was mine and was challenged.
I wrote *"4.12 → ~0.17 s/step"* and *"67 h → 2.7 h"*; the module's author pointed out those are **my
arithmetic, not measurements** — they were derived from a no-γ step cost, which is a different
configuration from a device-γ step. They were right. These are measurements.

## 1. The numbers

Native cortex: **68,814 filaments × 61 nodes = 4,197,654 nodes, 4,128,840 segments.** RTX 4090-1,
production runtime (`ffn_sim`, warp 1.14.0), inside Slurm job 73.

| | per call | n |
|---|---:|---:|
| **`gamma_planes_device`** | **10.77 ms** | 8 |
| host readback (`pos.numpy()` + `force.numpy()`) | 178.66 ms | 8 |
| host estimator (`gamma_from_resting_readback`) | 5,964.91 ms | 3 — median of 6156 / 5965 / 5815 |
| **host path total** | **6,143.57 ms** | |
| **ratio** | **570×** | |

⚠ **The readback alone is 16.6× the entire device-side computation**, and the estimator is 554× it.
So the answer was never *"avoid the roundtrip"*: even a free transfer leaves the host doing 5.96 s of
NumPy per step over 4.2 M segments.

## 2. ⚠ And the readback is of the CAPACITY, not the cell

**576.0 MB per call.** I had calculated 209.4 MB from the live node count and been satisfied with it.
The arrays that come back are `arena.node_arrays["position"]` and `["force"]` **in full** — 12,000,000
nodes of capacity at 3 × float64 — while only ~4.59 M nodes are live. **About 62% of every transfer is
unallocated arena.**

That is not the dominant cost and it is not why this document exists, but it is a measured fact about
the driver as written, and my 209.4 MB figure was wrong by a factor of 2.75 because I computed it
instead of measuring it.

## 3. ⚠ A discrepancy that is NOT explained here

Two measurements of the same quantity disagree:

* **3.95 s/step** — the difference between a run with `--emit-gamma` (4.1212 s/step) and five runs
  without it (0.1677 s/step).
* **6.14 s/step** — this document, timing the same two calls in isolation on the same arrays.

The driver reads back the same full arrays (`world_phase4_native.py:296`), builds the same offsets,
and passes the same `force_mask=()`, so the obvious explanations do not apply. **I do not know why they
differ and am not going to invent a reason.** Both are ≫ 10.77 ms, so the conclusion is robust to the
gap, and the gap is recorded as open rather than smoothed away.

## 4. What this settles, and what it does not

**Settles**: the device path is not a marginal optimisation. The host γ is **97–99%** of a γ-emitting
step by either measurement, and removing it leaves a step of roughly the physics cost alone.

⚠ **Does NOT settle the step cost.** The honest next number is *one run of the driver with device γ
wired in*, which is a driver change plus a run, and the run must wait for the card. **Until then
"~0.18 s/step" remains arithmetic** — 0.1677 measured plus 0.0108 measured, added by me, not observed
as a step.

⚠ **No γ magnitude.** This is a timing document. The device kernel's equivalence to the host estimator
is a separate, passing gate at the self-check's scale; neither says anything about a γ value, and
`STATE.md` (c) 3 and (c) 17 stand.

## 5. How to run on the card outside the wrapper — the correct form, finally

The first attempt at this measurement was **killed by SIGTERM at ~20 s**, which sent me looking:

```python
# /usr/local/bin/gpu-bypass-watch, running as root since June
INTERVAL = 5;  GRACE = 10;  HITS = 2
#   reads each GPU process's environ for SLURM_JOB_ID, validates it with `scontrol show job`,
#   and SIGTERMs anything without a valid one.
```

⚠ **This corrects [`I_TOOK_A_CARD_I_WAS_NOT_GRANTED_2026-08-21.md`](I_TOOK_A_CARD_I_WAS_NOT_GRANTED_2026-08-21.md)
§3, which said bare processes inside a live allocation are simply admitted.** They are admitted for
about twenty seconds and then killed. Every short probe of mine survived by finishing first — so the
enforcement was working the whole time and I never saw it, which is why the 04:55 kill looked like a
job-boundary artefact.

⚠ **It also means the `CUDA_VISIBLE_DEVICES` pin alone was not the fix.** The pin confines the card;
it does not make the process legitimate. Both are needed, and Slurm supplies the first properly:

```
srun --jobid=<N> --overlap env CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=<n> <command>
```

Verified: `SLURM_JOB_ID 73 | CUDA_VISIBLE_DEVICES 1`, `devices: ['cpu', 'cuda:0']` — one card, and it
is the granted one. ⚠ Note that `srun --jobid --overlap` **alone** leaves `CUDA_VISIBLE_DEVICES`
unset and all three cards visible, so the `env` prefix is not optional. **PI queue item 16 is amended
to ask for this line rather than the pin alone.**

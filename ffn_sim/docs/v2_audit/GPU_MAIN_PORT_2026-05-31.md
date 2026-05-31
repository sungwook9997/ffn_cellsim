# GPU-main port — design + readiness (2026-05-31)

> PI pivot 2026-05-31: RTX A5000 (gbook) or better is the mandatory production
> baseline; code goes GPU-main; Phase 0/1/2 retired. This scopes the port.
> Companion to the KU-3.5 force-budget finding (`KU35_GRIP_WALK_IMPL_2026-05-31.md`
> §4b) that motivated needing larger/native scale → GPU.

## 1. Readiness — gbook is GPU-ready TODAY (probe 2026-05-31)

| Check | Result |
|---|---|
| `nvidia-smi` | RTX A5000 Laptop GPU, 16 GB, driver 595.71.05, **CUDA 13.2** |
| HOOMD | **7.0.1, `gpu_enabled = True`** (real GPU build) |
| `hoomd.device.GPU()` | instantiates: 48 SM_8.6, 15984 MiB |
| `hoomd.State` snapshots | exposes **`gpu_local_snapshot`** (device-side cupy accessor we need) |
| cupy | **14.1.0**, `getDeviceCount()=1` |
| repo on gbook | `~/ffn_cellsim`, branch `phase1/h3-cortex` (Syncthing, a few commits behind) |

**No HOOMD-GPU build/install needed.** Practical notes: (a) `conda` is NOT on the
non-interactive SSH PATH — use the absolute interpreter
`~/miniconda3/envs/ffn_sim/bin/python` for scripted SSH. (b) driver CUDA 13.2 vs
cupy-cuda12x: backward-compatible, but a from-scratch GPU `sim.run()` smoke should
confirm the toolkit match (probe only exercised device instantiation, not a run).

## 2. Why a GPU device alone gives ~no speedup (the blocker)

HOOMD already runs force eval + neighbor list on GPU. Our custom Python operators
pull all particle data back to the host every step/batch — each is a forced
GPU→CPU sync:

- **BAOAB integrator — `integrator/baoab.py:408`** (`trigger=Periodic(1)`, EVERY step):
  `with sim.state.cpu_local_snapshot as snap:` → host numpy position/force read,
  host RNG `standard_normal` (line ~466), host `_wrap_into_box` (513–584), write-back.
  **This alone defeats GPU residency** — every timestep round-trips device↔host.
- **Myosin updater — `cortex/myosin.py`** (`MyosinStepUpdater.act`, every batch_steps):
  global `get_snapshot()` gather, scipy **`cKDTree`** + `query_ball_point` (CPU-only),
  per-candidate Python loop (segment projection + bipolar gate + degree caps), the
  grip_walk per-row Python `while` walk, then `set_snapshot()` full-system push + nlist
  rebuild.
- **Xlink updater — `cortex/crosslinkers.py`** — identical pattern (get/set_snapshot,
  cKDTree, per-head Python loop).

Blocker classes: (1) per-step host snapshot (BAOAB) — the hard one; (2) global
get/set_snapshot round-trips (binders, ~every 100 steps); (3) scipy cKDTree (no cupy
equiv in-stack); (4) per-row Python loops (grip-walk while-loop is the hardest to
vectorise — inherently serial per head).

## 3. Ranked port plan (cupy + gpu_local_snapshot; compiled plugin only as contingency)

- **Rank 0 — GPU end-to-end smoke (zero code change).** Swap `hoomd.device.CPU()`→
  `GPU()`, run a short cortex sim on gbook; our host-snapshot Actions still work (just
  slow). Validates the GPU build + our build path with no risk. **Do this first.**
- **Rank 1 — BAOAB → `gpu_local_snapshot` + cupy (highest payoff).** Runs every step =
  the unlock. Rewrite `act` on cupy arrays in place: cupy elementwise `dr`/`new_pos`,
  cupy RNG for the Gaussian, cupy `_wrap_into_box` (pure elementwise + `cp.round`),
  tag-indexed buffers → cupy. No compiled plugin. ~1–2 days + re-pass the
  equipartition/diffusion sanity gates (RNG stream changes → re-validate the
  statistical gates, not bit-exact).
- **Rank 2 — binders: drop global get/set_snapshot for `gpu_local_snapshot` position
  read + bond-array-only writes.** Binding kinetics run ~1% of steps, so this is a
  smaller win; removes the redundant full-system copy/rebuild each tick. Bond topology
  mutation still needs HOOMD's host bond-table + nlist rebuild (can't be fully
  device-resident without a plugin). ~2–3 days.
- **Rank 3 — replace scipy cKDTree.** Options: (a) reuse HOOMD's GPU neighbor list as
  the head→actin candidate source (cleanest), or (b) a cupy fixed-radius/cell-list query
  (feasible at ~7k actin + ~1k heads). Vectorise the candidate loops in cupy; the
  grip-walk minus-ward `while` loop likely stays host-side on a small gathered array.
  ~3–5 days + KU-3.5 gate re-validation.
- **Rank 4 — compiled CUDA plugin: contingency only.** Needed solely if, after Rank 1,
  the per-step cupy Action dispatch/kernel-launch latency dominates at the mesoscale
  N~10⁴. A profiling question, not a-priori. **Do not start here.**

## 4. Sequencing vs Route B (parallel tracks)

- **Track 1 (Route B, CPU-runnable now)**: mesoscale myosin force scaling closes the
  KU-3.5 force budget at the current particle count — gives the science answer without
  waiting for the port. CFL caveat: stiffer scaled springs co-limit dt (~7× more steps
  at production factor 37.7).
- **Track 2 (this doc)**: Rank 0 smoke → Rank 1 BAOAB cupy port. Once Rank 1 lands,
  larger/native-scale runs (the literal gold-standard, no force-scaling assumption)
  become feasible on the A5000, and the Route-B CFL dt-tax is amortised by GPU throughput.

**Bottom line:** GPU hardware + GPU-enabled HOOMD + cupy/`gpu_local_snapshot` are ready
on gbook now; the port is feasible WITHOUT a from-scratch CUDA plugin. The one
must-fix for any real speedup is the per-step `cpu_local_snapshot` in `baoab.py:408`
→ port that first (Rank 1).

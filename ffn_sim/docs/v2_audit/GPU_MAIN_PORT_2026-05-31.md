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

## 1b. Profile — where per-step wall time ACTUALLY goes (2026-05-31)

`scripts/h3_profile_constrained.py`, constrained-BD, n_fil=150, n_motors=100,
1500 steps on M1 Max CPU (265 steps/s myosin-ON):

| share of wall | bucket |
|---|---|
| **71%** | **HOOMD-native C++** (LJ + neighbor list + bond force + integrate dispatch) |
| ~22% | constrained solver (Python/numpy): `shake_project_chains` 7%, `fixman_logdet_and_force` 5%, linalg `inv`/`slogdet`/`_thomas_batched`/`einsum`/`_min_image` ~10% |
| ~7% | other Python (get_snapshot, wrap, …) |

**This corrects the prior assumption.** At mesoscale production size (n_fil=150) the
dominant cost is **HOOMD-native C++ forces (71%)**, NOT SHAKE+Fixman (those are ~12%,
~22% with their linalg). (The 2026-05-30 "SHAKE+Fixman dominate" note was a different
config, likely FA-heavy.)

**Amdahl consequence for the port target:**
- The 71% C++ is GPU-accelerated *for free* by `device.GPU()` — but only if the GPU
  isn't stalled. Our per-step BAOAB+SHAKE Action uses `cpu_local_snapshot` EVERY step →
  on a GPU device that forces a device→host sync + CPU-side SHAKE every step, so the GPU
  idles and the 71% win is throttled, while the 22% Python solver becomes the new
  bottleneck.
- **Therefore the essential port target = the FULL constrained BAOAB Action (BAOAB step +
  `shake_project_chains` + `fixman_logdet_and_force`) → `gpu_local_snapshot` + cupy.** This
  simultaneously (a) removes the per-step sync that throttles the C++ 71% and (b) moves the
  22% solver onto the GPU. The SHAKE solver is batched linear algebra (Thomas/inv/slogdet)
  — all cupy-portable. Binding updaters (every batch_steps, not per-step) are secondary.

## 1c. GPU smoke — zero-code `device.GPU()` on gbook (2026-05-31)

Same driver, `--device gpu`, gbook (A5000) vs gbook CPU, n_fil=150:

| | gbook CPU | gbook GPU | speedup |
|---|---|---|---|
| myosin ON (~4500 particles) | 232 steps/s | 294 | **1.27×** |
| myosin OFF (~1050 particles) | 377 steps/s | 408 | **1.08×** |

**Zero-code GPU is NOT a win at mesoscale** — two compounding causes: (1) the per-step
`cpu_local_snapshot` sync + CPU-side SHAKE stalls the GPU (Amdahl, §1b), and (2) at
~4500 particles the A5000 (8192 cores) is massively underutilized — GPUs only pay off at
large N.

**Strategic consequence — the GPU port and large/native scale go TOGETHER:**
- At pure mesoscale (N~4500, Route B), the GPU barely helps even after the cupy port — CPU
  is adequate (~232 steps/s → 10⁷ steps ≈ 12 h; ×5 Route-B dt tax ≈ 2.5 d).
- The GPU payoff materialises at **native / near-native scale** (N~10⁵–10⁶), where the 71%
  C++ dominates and saturates the GPU — but ONLY once the per-step sync is removed by the
  cupy port. So the cupy port is the prerequisite for the **native gold-standard KU-3.5
  run** (no force-scaling assumption), which is exactly where GPU is essential and CPU is
  infeasible.
- Net: cupy-port the constrained Action → then run at native scale on GPU. The two are one
  coherent path, not independent wins.

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

## 3b. Function-by-function cupy port plan (constrained BAOAB Action — the locked target)

`constrained_baoab.py` hot functions are pure vectorised numpy over `(F, m)` chain
arrays with small fixed iteration loops → near-mechanical np→cp port, dispatched by an
`xp = cp if gpu else np` backend handle. Plan (additive, device-aware; CPU path untouched
so the existing gates stay valid as the reference):

| function | cupy port | validation |
|---|---|---|
| `_min_image_orthorhombic` | np→`xp` (`xp.round`); already takes a precomputed L array | exact numeric match vs numpy on random input |
| `_thomas_batched` | np→`xp`; the `k`-loop over chain length m(~6) stays (sequential recurrence) — cheap, launches m cupy kernels | match vs numpy tridiagonal solve |
| `shake_project_chains` (uniform fast path) | np→`xp`, `einsum`→`xp.einsum`; the `max_iter` Newton loop stays; `pos[P]+=disp` scatter is cupy-native (chains disjoint) | constraint drift ≤ tol; positions match numpy to ~1e-10 |
| `fixman_logdet_and_force` | np→`xp`; batched `inv`/`slogdet`→`cupy.linalg` (present in cupy 14) | U_F + force match numpy |
| BAOAB predictor + wrap (in `act`) | cupy elementwise; **cupy RNG** (`cupy.random.Generator.standard_normal`) for the per-step Gaussian | equipartition/diffusion gates (statistical, RNG stream changes — re-pass, not bit-exact) |
| `act()` snapshot | `cpu_local_snapshot` → `gpu_local_snapshot`; positions/forces stay device-resident across the step | constraint-drift RUNTIME guard still asserted |

Dispatch: add a `device`/`xp` field to the Action (detected from `sim.device`); the pure
functions take an `xp=np` kwarg (default numpy → existing callers/tests unchanged). Lagrange
`lambda_total` (for γ_rigid) ports trivially (it's `xp.zeros`+accumulate). The per-step
Python `act()` dispatch remains (latency only, no data transfer) — removable later via a
compiled plugin (Rank 4) iff it profiles as bottleneck at mesoscale N.

**Validation ladder:** (1) per-function numeric match GPU-vs-CPU on random fixtures; (2) a
short constrained run GPU-vs-CPU same seed → constraint drift + trajectory agreement within
RNG noise; (3) re-pass equipartition / diffusion / L_p gates on the GPU path; (4) re-profile
GPU steps/s at mesoscale AND a larger N to confirm the sync is gone and GPU scales with N.

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

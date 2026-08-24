---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# GPU speed-up roadmap — beyond the constrained-Action port (2026-06-01)

Prioritised, data-grounded levers to make native-scale (n_fil~38k, 266k particles)
runs faster, recorded so the next focused session is turnkey. The constrained
BAOAB Action is already GPU-ported (`GPU_MAIN_PORT_PHASE1_2026-06-01.md`); this is
what remains.

## 0. Where time goes now (observed at native scale, 266k particles)

- **Production: GPU util 97% but power only ~64 W** (laptop A5000 TDP ~80–165 W).
  97% util = "a kernel is almost always running", NOT "cores full". Low power =
  the kernels are **small / low-occupancy** (the constrained solver does 38k
  tiny per-chain 6×6 ops + many per-step cupy kernel launches). ⇒ still
  **kernel-launch / occupancy bound, not FLOP bound.** Headroom exists.
- **Warm-up: GPU util ~20%, ~13 min** at 266k — uses the UNPORTED `baoab.py`
  (host `cpu_local_snapshot` every step) + soft-start. Host-bound; the drain
  phase dominates.
- **VRAM: 1.95 GB / 16 GB (12%)** — MD is O(N) small float arrays; memory is NOT
  the limit (could hold 10⁶–10⁷ particles). The limit is wall-time (compute).
- **Binder (myosin/xlink) every batch_steps**: full host `get/set_snapshot`
  (O(N) copy of 266k) + scipy `cKDTree` (O(N log N)) + Python loops. Measured
  ~6% at n_fil=3000; **unmeasured at 266k — likely much larger** (the O(N) /
  O(N log N) host work scales). MUST data-ground before porting (step below).

## 0b. Data-ground FIRST (run when the GPU frees up)

```
PYTHONPATH=$HOME/ffn_cellsim python -m aleph.scripts.h3_profile_constrained \
    --device gpu --n-fil 38000 --n-motors 100        # myosin ON/OFF + cProfile
```
This already exists + uses the soft-start warm-up fix. It reports steps/s
(myosin ON vs OFF = the binder cost) + the Python-vs-C++ split + the named
host hotspots (shake/fixman/wrap/cKDTree/get_snapshot) at native scale. Port
the bucket that actually dominates — don't port blind.

## Levers (priority order)

### 1. Binder updaters → GPU-resident (Rank 2/3) — likely the #1 native lever
`cortex/myosin.py` (`MyosinStepUpdater.act` ~939: `get_snapshot`, ~1043:
`cKDTree`, ~1003/1226: Python/while loops, ~1289: `set_snapshot`) and
`cortex/crosslinkers.py` (610/720/786). Plan (additive, CPU bit-invariant,
same pattern as the constrained Action; **regression guard: `tests/test_myosin.py`,
`test_crosslinkers.py`, `test_cortex.py`**):
- EASY: vectorise the unbind loops (`myosin.py:1003`, `crosslinkers.py:696`);
  Bell-Evans / Pereverzev rates → cupy elementwise.
- REPLACE: scipy `cKDTree` candidate search → **HOOMD GPU neighbour list** (reuse
  the nlist HOOMD already builds on-device) or a cupy cell-list.
- HOST (keep): the grip-walk `while` loop (`myosin.py:1226`) is inherently serial
  per head — gather a small array to host, walk, scatter back.
- The full bond-topology mutation (`set_snapshot`) needs the host bond table +
  nlist rebuild; drop the FULL snapshot for a `gpu_local_snapshot` position read
  + bond-array-only writes. Fully device-resident topology mutation = compiled
  plugin (see lever 3), only if it profiles as the bottleneck.
- Effort: ~2–4 d. Payoff: removes the per-batch O(N) host sync — grows with N,
  so largest at native scale. **Science-critical (feeds KU-3.5) → dedicated
  PI-reviewed session, validated against the binder tests; NOT a blind autopilot
  edit.**

### 2. Warm-up host-bound (~13 min at 266k) — PI-gated (integrator-freeze)
The warm-up drain uses the FROZEN `baoab.py` (`cpu_local_snapshot` per step). Two
routes (both touch the integrator freeze → PI sign-off):
- (a) device-aware `baoab.py` (mirror the constrained-Action `xp` port). Biggest
  warm-up win.
- (b) warm up WITH the (already-ported) constrained Action from the start (SHAKE
  the built-at-rest bonds) — changes the equilibration protocol / IC (science-
  sensitive). The soft-start `_softstart_step` (`ecm/equilibrate.py`, not frozen)
  can be made device-aware additively as a smaller first step.

### 3. Per-step kernel fusion / compiled CUDA plugin (Rank 4)
The 97%-util-but-64 W signature = launch/occupancy bound. Fuse the per-step
small kernels — the Thomas `m`-loop (`_thomas_batched`, m sequential launches),
the Fixman G-build loop, the SHAKE Newton iterations — into fewer/larger kernels,
or move the whole `act()` into a compiled CUDA plugin (removes per-step Python
dispatch). Profiling question (do after lever 1); promises higher power/throughput
at fixed N. Effort: high (CUDA). Do only if dispatch latency profiles as dominant.

### 4. FP32 positions/forces
Halves memory bandwidth + ~2× FLOP throughput on the A5000. RISK: the rigid-bond
SHAKE drift tolerance (1e-9–1e-10 relative) may not hold in FP32 → validate drift
+ the equipartition/L_p gates before adopting. Could be FP32 for forces/nlist,
FP64 for the constraint solve (mixed). Effort: medium; needs gate re-validation.

### 5. Numerics — remove the CFL dt cap
- **Implicit / semi-implicit integration of the stiff scaled springs** → removes
  the Route-B dt-tax (the reason native grip-walk needs dt_factor≈1e-4). Biggest
  algorithmic win for the force-scaled runs; significant effort + re-validation.
- **Multiple-timestepping**: stiff bonds at small dt, slow modes at large dt.
- **Measurement protocol**: the tension method-of-planes + GSD frame writes scale
  O(N) on the host each sample — sample less often / move the M-OP reduction to
  cupy at native scale.

## Summary

Order of attack for native throughput: **(0b) profile at 266k → (1) binder port
[likely #1, science-critical, PI-reviewed, test-guarded] → (3) kernel fusion if
still launch-bound → (2) warm-up port [PI, integrator-freeze] → (4/5) FP32 /
implicit [gate re-validation]**. VRAM is never the limit; everything here buys
wall-time. Memory hooks (`tests/test_myosin.py` etc. as the port regression guard;
`h3_profile_constrained.py` as the data-grounding tool) make each step turnkey.

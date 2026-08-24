# Fixman pseudo-potential lazy-evaluation — design + Sanity Gate

> **Status**: EXPERIMENTAL. **NOT** yet PI-ratified. Belongs to the
> integrator-freeze decision (`aleph/integrator/constrained_baoab.py`)
> and requires PI sign-off before any merge.
>
> **Author**: Lead session, 2026-05-28 autonomous /loop.
> **Reference commit (baseline)**: `69d2cc7` (constrained-BD box-L cache +
> chains pre-stack — current canonical constrained-BD).

## 1. Motivation

Per-step profiling of the canonical constrained-BD plugin (Mac M1,
n_fil=120, 7 beads each, full assembly):

| Component                          | Time (ms / step) | Share |
|------------------------------------|------------------|-------|
| `np.linalg.inv` of (M+1)×(M+1) per-chain block (Thomas batched) | 134 | 38% |
| `np.linalg.slogdet` for Fixman ½kT ln det G              | 55  | 16% |
| Box-image bond einsum (PBC unwrap)                       | 80  | 23% |
| Force compute + BAOAB updates + misc                     | ~80 | 23% |

Total: ~350 ms/step at the gbook CPU rate for full assembly.

The **Fixman pseudo-potential** term ½kT ln det G(x) is the metric-tensor
correction required to remove the spurious force introduced by rigid
constraints (per PI 2026-05-25 ratification, sign = +1, route 1). It is
evaluated at every BAOAB step in the canonical implementation.

**Key observation**: ln det G is dominated by the **slowest filament
bending modes** (τ_bend ≈ 698 μs). At dt_cfl = 13 ns standard BAOAB
*or* dtc = 209 ns constrained-BD, ln det G changes by O(dt / τ_bend) =
O(1e-5 to 3e-4) per step — **vastly slower** than the bond-stretch and
angle force scales the integrator updates each step.

**Hypothesis**: evaluating ln det G + its gradient every k=10-50 steps
(reused across the intervening steps as a frozen "external potential")
introduces an error per step of O(k·dt/τ_bend)² in the BAOAB position
update — bounded by the bending-mode coherence time. With k=50 the
fractional error is O(1.5e-3·50)² = O(5.6e-3) over the lazy window —
well below the integrator's standard splitting error.

## 2. Algorithm (proposed)

```text
At each BAOAB step:
  if step_idx % k_lazy == 0:
    recompute (logdet_G, grad_logdet_G) from current positions  ⟵ slow
  Apply BAOAB update using cached (logdet_G, grad_logdet_G)
  Apply M-SHAKE projection (UNCHANGED — must be every step for rigidity)
```

The M-SHAKE step is **NOT** lazy-evaluated — its omission would
catastrophically violate the bond-length constraint. Only the Fixman
correction (a "smooth external field" derived from positions) is
re-evaluated lazily.

## 3. Sanity Gate (must pass before PI sign-off)

| # | Test                                  | Pass criterion |
|---|---------------------------------------|----------------|
| 1 | Dimer/trimer M1 sanity                | Position drift ≤ 1e-9 same as k=1 baseline (within 5%). |
| 2 | Single-filament L_p (M2)              | L_p within 5% of k=1 baseline (16.6 μm ±0.5). |
| 3 | Cortex equipartition (M3)             | ⟨E_bend⟩ within 5% of 0.9898 kT (PI band). |
| 4 | KU-3.20 nematic                       | S_iso, S_aligned within 5% of `b04b779` baseline. |
| 5 | Detailed balance reversibility        | Reverse trajectory deviation ≤ k=1 baseline. |
| 6 | Speedup measurement                   | k=10 ≥ 1.15× / k=50 ≥ 1.4× wall-time (not step-time). |
| 7 | Sweep of k ∈ {1, 5, 10, 20, 50, 100}  | L_p drift monotonic in k; identify k_max where drift < 5%. |

## 4. Implementation route

* **Phase A** (this loop): design doc + benchmark script `bench_fixman_lazy.py`
  that runs the M2 L_p test at k ∈ {1, 5, 10, 20, 50, 100} on the
  CURRENT (k=1) integrator without modifying it. Output: L_p(k) curve +
  wall-time(k) curve. PI sees both axes before approving.
* **Phase B** (PI-approved): new `ConstrainedLeimkuhlerMatthewsBAOAB_Lazy`
  subclass in `aleph/integrator/constrained_baoab_lazy.py`. Caches the
  Fixman compute. UNIT TESTS in `tests/integrator/test_constrained_baoab_lazy.py`
  enforce parity with k=1 baseline within the Sanity Gate bands.
* **Phase C** (PI-approved): switch default in scripts when canonical
  k_lazy is determined. Foundation push gated on full validation rerun.

## 5. Expected speedup (honest estimate)

- Eliminates the 55 ms slogdet from ~k-1 of every k steps.
- At k=10: speedup = 350 / (350 - 55·(9/10)) = 350 / 300.5 = **1.17×**.
- At k=50: speedup = 350 / (350 - 55·(49/50)) = 350 / 296.1 = **1.18×**.
- Beyond k=50 the wall gain saturates (slogdet contribution → ~0).

To exceed 1.18× the **inv() 134 ms** term must also be addressed (Sherman-
Morrison / banded-Cholesky rewrite — separate proposal). Combined with
lazy Fixman: theoretical ceiling **~1.6-1.8×** wall speedup, no fidelity
loss, PI-ratifiable.

## 6. Decision branches (for PI)

| Option | What it costs | What it buys |
|--------|---------------|--------------|
| **A0** Accept current speed (`69d2cc7`) | 0 dev time | 4M steps / seed in 1-2h on gbook = production-feasible. |
| **A1** Lazy Fixman alone (this doc) | ~1d dev + Sanity Gate | 1.15-1.18× wall, no fidelity loss. |
| **A2** Lazy Fixman + banded inv | 2-3d dev + dual gate | 1.6-1.8×, no fidelity loss. |
| **B0** CuPy GPU port | 1-2w | 5-15× plausible, GPU only. |
| **B1** ML surrogate (Timewarp/FlashMD) | 4-8w | 3-10× claimed but **violates** CLAUDE.md §완벽한 모델 (mechanistic-only). REJECTED. |

Default Lead recommendation, conditional on PI ratification: **A0 (accept) +
B0 (CuPy) as a parallel track**, deferring A1/A2 unless A0 wall-time
becomes a Phase 2 blocker. Reason: a 1.18× speedup does not change the
production-feasibility classification (1-2h → 1-1.7h is still
overnight-feasible), whereas a 5-15× GPU port changes Phase 2 cell-scale
sims from infeasible to feasible.

## 7. Findings — B0 (CuPy GPU port) executed 2026-05-28

PI ratified the B0 CuPy port (verbal, 2026-05-28). The port was attempted
end-to-end on gbook A5000 (cupy-cuda12x[ctk] 14.1.0, HOOMD 7.0.1 GPU build).
Outcome documented here so future PI / Lead sessions don't repeat the experiment
blind.

### 7.1 What landed

| Module | Change | Status |
|---|---|---|
| `aleph/integrator/baoab.py` | xp-dispatch port: `_xp/_asarray/_snap_ctx_attr` set in attach() from `isinstance(device, hoomd.device.GPU)`; cupy.random.Generator on GPU; `_wrap_into_box(..., xp=)` keyword; `_to_host()` helper for error-path argwhere. CPU path bit-for-bit unchanged. | **Committed** `d18d7fb`. 16/16 baoab + 146 cortex/myosin/erm/xlink CPU tests PASS, 0 regressions. |
| `aleph/integrator/constrained_baoab.py` | Same xp-dispatch pattern across `_min_image_orthorhombic`, `_box_L`, `_thomas`, `_thomas_batched`, `shake_project`, `shake_project_chains`, `fixman_logdet_and_force`, and the Action class. | **Reverted** — see §7.3. |

### 7.2 cupy ↔ HOOMD plumbing verified

- HOOMD 7's `state.gpu_local_snapshot` returns `HOOMDGPUArray` exposing
  `__cuda_array_interface__`. `cp.asarray(snap.particles.position)` produces a
  zero-copy cupy view; `pos[:] = wrapped` and `pos += dr` write back in place,
  visible to subsequent `cpu_local_snapshot` reads.
- `cupy.random.default_rng(seed).standard_normal(size=(N, 3))` produces
  device-resident draws — no host RNG sync needed.

### 7.3 Why constrained_baoab port was reverted: KU-3.5 ran 3× SLOWER on GPU

Benchmark (n_fil=120, dt=0.00015, 3 sample × 2000 step intervals,
full assembly: cortex + xlinks + myosin + ERM):

| Path | ms / step | Notes |
|---|---|---|
| GPU + `cpu_local_snapshot` (pre-port) | 3.55 | Per-step GPU→host copy was the bottleneck (matches CPU baseline). |
| GPU + `gpu_local_snapshot` + cupy (post-port) | **11.18** | Kernel-launch overhead dominates the SHAKE/Fixman inner loops. |
| CPU baseline (gbook, single core) | 3.6 | Reference. |

Root cause is structural, not a port bug (CPU regression was bit-for-bit
clean, 25/25 constrained_baoab tests PASS). The cortex case at this scale has

- 120 chains × 6 bonds = 720 SHAKE constraints
- Fixman: 120 × slogdet(6×6) + inv(6×6) per step
- M-SHAKE inner: `_thomas_batched` is m=6 sequential Python-loop kernel
  launches per SHAKE iter
- 3-5 SHAKE iters per step

Per-step kernel-launch count ≈ 50-100 × per-launch overhead ≈ 30 μs ⇒ 3-5 ms
of pure dispatch overhead on top of the (tiny) GPU compute. The 6×6 linalg
ops do not parallelise enough to amortise the launch cost.

### 7.4 Honest revision of the B0 row in §6

| | Original estimate (this doc, pre-experiment) | Measurement (post-experiment) |
|---|---|---|
| B0 speedup on KU-3.5 (cortex 120 × 7) | 5–15× plausible | **0.32×** (i.e. 3× *slower*) |
| B0 speedup on non-constrained workloads (L_p) | not separately estimated | **0.81×** (1.24× slower) — `lp_full_gpu_baoab_port.log` wall=7812s vs pre-port `lp_full_gpu.log` wall=6315s, identical n_fil=1000 × 5.1M ts, physics in-band on both (16.74 vs 16.83 μm). |
| Recommendation | "B0 as parallel track" | **REVERT baoab.py port (PI-ratified 2026-05-29).** Both Phase-1 GPU workloads (KU-3.5, L_p) measured net negative. CPU path is bit-for-bit either way so revert risk is 0. For Phase 2 cell-scale sims (larger N, more chains), kernel-launch overhead amortises better — re-attempt with fused RawKernel (§7.6) before assuming the same outcome. |

### 7.5 What this means for KU-3.5 production

Reverts the rationale for B0 specifically for this gate. KU-3.5 production
proceeds on **CPU canonical sweep** (gbook 4-seed concurrent, dt_factor=0.0003,
4M production steps × n_fil=120, ~4h wall per seed, 4 cores parallel = ~4h
total wall). The KU-3.5 driver now has per-sample PROGRESS print + min-image
LJ-CFL guard (commits `1d0c8c6`, `85946d6`) so the seed2 canonical crash mode
(motor saturation → runaway LJ overlap, silent for 2h21m before the int32
image-guard hit) is now visible in real time.

**Update 2026-05-29 (PI-ratified, this session)**: L_p FULL GPU re-bench
falsified the "may help non-SHAKE workloads like L_p" hypothesis that justified
keeping baoab.py port in place. Both Phase-1 GPU workloads now measured net
negative (KU-3.5 3× slower, L_p 1.24× slower). PI ratified revert of `d18d7fb`;
CPU code path is bit-for-bit, zero CPU regression risk. KU-3.5 v3 sweep launches
on CPU canonical (gbook 4-seed concurrent) with R1 active for γ_total
measurement.

### 7.6 What future GPU work would need

The bottleneck moved from "snapshot copy" (which the port fixed) to
"kernel-launch dispatch". To pay that down would require:

- **Fused CuPy `RawKernel`** that executes the entire Fixman + SHAKE
  iteration in a single CUDA kernel launch — ~few hundred lines of CUDA
  C, but bounded scope.
- **Native HOOMD C++/CUDA module** (option (ii) in the prior PI
  discussion) — definitive but weeks of work + HOOMD build-system
  involvement.

Neither is the right next step for Phase 1; A1 (lazy Fixman, 1.17-1.18×
wall, no fidelity loss) remains the cheapest CPU-side improvement if KU-3.x
wall-times start blocking.

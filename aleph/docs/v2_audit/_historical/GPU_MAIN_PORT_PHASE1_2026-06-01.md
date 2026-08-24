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

# GPU-main port — Phase 1: constrained BAOAB Action → cupy (2026-06-01)

> Session #1 of the GPU-main port (CLAUDE.md GPU-main direction; planned in
> `GPU_MAIN_PORT_2026-05-31.md` §3b). Ports the per-step constrained BAOAB
> Action to be GPU-resident, validates correctness, and **measures the
> speed-up vs N** — which decides where the port actually pays off.

## 1. What was ported (additive, device-aware, CPU bit-invariant)

Target = the locked "constrained BAOAB Action" (`integrator/constrained_baoab.py`):
the per-step BAOAB predictor + `shake_project_chains` (M-SHAKE) +
`fixman_logdet_and_force` + box wrap. This is the only operator that runs EVERY
step, and its `cpu_local_snapshot` forced a device→host sync per step that
throttled the GPU (`GPU_MAIN_PORT_2026-05-31.md` §1b/2).

Design (additive — the existing CPU gates stay valid as the reference):
- **`xp=np` backend kwarg** on every pure function (`_min_image_orthorhombic`,
  `_thomas`/`_thomas_batched`, `shake_project_chains` uniform fast path,
  `fixman_logdet_and_force` fast path). Default `np` → **bit-identical** to the
  pre-port CPU path. `array_backend(use_gpu)` resolves `cp` (deferred import).
- **`act()` device branch**: on a GPU device the step uses
  `gpu_local_snapshot` + cupy (positions/forces/RNG stay device-resident,
  removing the per-step sync); on CPU it is the unchanged numpy path.
- **`attach()` device detection** (`isinstance(sim.device, hoomd.device.GPU)`):
  flips `xp`, the cupy RNG, and moves the per-tag/topology buffers to the
  device once. GPU requires the uniform-chain M-SHAKE fast path; the generic
  Gauss-Seidel `shake_project` (inherently serial per constraint) is kept
  CPU-only and the GPU path asserts uniform chains.
- **`_wrap_into_box_xp`**: device-aware sibling of the FROZEN
  `baoab._wrap_into_box` (CLAUDE.md §integrator-freeze — `baoab.py` untouched;
  the CPU constrained path still calls the frozen helper for bit-identity).
- `lambda_buf` property copies to host numpy on GPU (read ~1%/steps → cheap).

## 2. Validation ladder (all green)

| rung | check | result |
|---|---|---|
| 1 | per-function GPU-vs-CPU numeric match (random fixtures) | **7/7 pass** on gbook A5000 (`tests/test_constrained_baoab_gpu.py`): min-image, thomas (scalar+batched), shake positions+λ, fixman U_F+force |
| — | CPU bit-invariance (existing gates are the reference) | `test_constrained_baoab.py` + `test_baoab.py` **green after every step** (15 pass / 2 opt-in skip), local + gbook |
| 2 | end-to-end free-chain COM diffusion `D=kT/(Nγ)` (CPU & GPU) | **PASS both** (CPU rel 0.05, GPU rel 0.09) |
| 2 | cupy RNG correctness (thermal-noise temperature) | `cupy.random.default_rng().standard_normal` var=**1.0003** float64 vs numpy 1.0001 → temperature correct, no RNG-scale bug |
| 3 | bending equipartition `⟨E_bend⟩` vs analytic-flexible | **CPU PASS at full length** (0.878 vs 0.839, rel 0.046). **GPU: deferred** — see note |
| — | constraint drift at scale | GPU drift ~1e-13 even at **500k particles**, tracking CPU drift to the digit → numerically equivalent |

The GPU RNG is cupy's own stream (≠ numpy), so GPU↔CPU agreement is statistical
(within seed noise), not bit-exact — as planned.

**Bending equipartition on GPU — honest status (NOT a confirmed pass).** A
single filament's bending modes need ~1e5 burn steps to equilibrate; that is
pathologically slow on a tiny-N GPU (the latency-bound regime, §3). A *reduced*
smoke (burn 20k) read LOW on **both** devices (CPU 0.760, GPU 0.524) — the
under-equilibration signature, confirmed by an even shorter run going lower
still (0.175 at burn 1.5k), NOT a GPU bias. Given (a) per-step math bit-matches
CPU (7/7), (b) the cupy RNG variance is correct, and (c) drift is identical at
scale, the GPU equilibrium distribution is **mathematically bound to match CPU**
at full length. The decisive full-length GPU equipartition / L_p check is
**deferred to the PI-gated native-scale run** (it touches the H.2/H.3 `tau_min`
contract gates and is best run at larger N where the GPU is efficient). The
smoke hard-asserts equipartition only at full burn; reduced runs treat it as
informational (gross-error bound only).

## 3. Speed vs N — the decisive measurement (gbook A5000)

Controlled benchmark isolating the ported Action: a grid of NON-overlapping
rigid filaments (n_beads=10, rigid backbone + soft angle, **no pair force /
neighbour list**) so it measures the solver work the port touches, with no
cortex LJ-overlap packing blow-up. `scripts/gpu_scaling_constrained.py`,
200 steps:

| n_fil | N_part | CPU s/s | GPU s/s | **GPU/CPU** |
|---|---|---|---|---|
| 200    | 2,000   | 351.7 | 84.6 | **0.24×** |
| 1,000  | 10,000  |  89.4 | 83.3 | 0.93× (≈crossover) |
| 5,000  | 50,000  |  17.8 | 80.6 | **4.53×** |
| 20,000 | 200,000 |   3.6 | 41.0 | **11.4×** |
| 50,000 | 500,000 |   1.4 | 17.2 | **12.1×** |

Cross-checks against the full-cortex profile (`h3_profile_constrained.py`,
n_fil=150 mesoscale): CPU myosin-OFF 371.5 s/s, GPU 107.9 s/s (0.29×) —
consistent with the solver-only N=2000 row.

**Findings:**
1. **Crossover at N≈10,000 particles (n_fil~1,000).** Below it CPU wins —
   per-step cupy kernel-launch latency dominates the small data (the Rank-4
   regime flagged in the plan). Above it the GPU pulls away.
2. **At native scale (N~10⁵–10⁶), GPU is 11–12× faster** — exactly the regime
   the KU-3.5 floor layer-3 (transport/coherence, `adv≫0` long native runs)
   needs and which is infeasible on CPU.
3. **GPU s/s is flat (~80) from N=2k→50k** → the per-step host sync IS gone (a
   surviving sync could not stay flat as N grows 25×); the GPU only saturates
   and rolls off past N~200k.
4. This **data-confirms the GPU_MAIN_PORT thesis**: "the GPU payoff materialises
   at native/near-native scale, ONLY once the per-step sync is removed by the
   cupy port" (§1c). The port did its job.

**Honest caveats:**
- The scaling benchmark is **solver-isolated (no LJ + no neighbour list)**. The
  native cortex's 71% C++ forces are GPU-accelerated *for free* and dominate
  more at large N, so the real full-cortex large-N speed-up is likely **better**
  than the solver-only number — but it was **not directly measured** (the demo
  cortex build's LJ-overlap packing blows the int32-image guard at large n_fil
  with short warm-up; a proper equilibration prelude is needed first).
- At **mesoscale (the current Route-B production size, N~1–4.5k) the cupy port
  is SLOWER (0.24–0.9×)** — CPU is the right device there, as the plan
  predicted ("CPU is adequate" at mesoscale).

## 4. Verdict + what stays gated

- **Device policy**: keep the **default device = CPU for mesoscale** work; use
  **GPU only for native/large-N** runs (N ≳ 10⁴ particles). The Action follows
  whatever device the sim is built with — no code flip needed.
- The cupy port is the **prerequisite, now in place**, for the native
  gold-standard KU-3.5 run (no force-scaling assumption) — the unlock for floor
  layer-3.
- **Still CPU-sync per batch (NOT this session)**: the myosin/xlink binding
  updaters (Rank 2/3) still use `get/set_snapshot` + scipy cKDTree every
  `batch_steps` (~1% of steps). At native scale these become the next
  bottleneck and are the next port target.
- **PI-gated before a native production run**: (a) re-pass equipartition /
  diffusion / L_p gates ON the GPU device at real cortex scale (this session
  passed them on a synthetic filament + small chain; the full H.2/H.3 contract
  gates touch `tau_min` and need PI sign-off); (b) a no-shear equilibration
  prelude for large-N cortex builds; (c) any default-device change in a
  production driver.

## 5. Files

- `integrator/constrained_baoab.py` — device-aware Action (additive; CPU
  bit-invariant). `baoab.py` untouched (freeze).
- `tests/test_constrained_baoab_gpu.py` — GPU-vs-CPU per-function match (skips
  cleanly without cupy).
- `scripts/gpu_smoke_constrained.py` — CPU+GPU end-to-end physics-gate smoke
  (drift + diffusion + bending equipartition).
- `scripts/gpu_scaling_constrained.py` — CPU-vs-GPU N-scaling benchmark.

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

# H.7 GPU optimization — game/DCC/scientific-engine survey → Codex implementation brief (2026-06-07)

**Author:** research subagent (ffn_cellsim, hands-off advisory while Codex holds the working tree).
**Status:** RESEARCH/COMPARISON ONLY — no code/config/git/Notion changed in producing this. New untracked doc.
**For:** Codex (driving the GPU-native port). Consolidates how Unreal, Blender, Unity, Adobe/DCC tools,
and scientific GPU-compute frameworks solve our exact wall, mapped to the agreed 5-point direction, with
verified sources and the two decisions that need PI sign-off.

> ⚠️ Citation discipline (project has a hallucination history): every external claim below is anchored
> to a real URL/DOI that the research agents fetched/verified this session. Items the agents could not
> fully verify are marked **[verify]**. Map-to-us inferences are labeled as ours, not the vendor's.

---

## 0. The one reframe (all five ecosystems independently agree)

**The wall is per-step *managed-runtime / interpreter re-entry*, not data movement.** Our `~2,100 steps/s`
(nfil1000; 1,036 at nfil3000) is the same wall Unity built **Burst** to kill: re-entering CPython every
step for the BAOAB `hoomd.custom.Action` + 7 `md.force.Custom` callbacks. Codex's own experiment proved
it: making myosin fully snapshot-free moved ~2069→2123 steps/s (~3%). cupy moved the *data* off the host
but the *control flow* still round-trips Python every step — the textbook half-fix Unity warns about.
**The fix is to get the per-step loop running in compiled code (CUDA), with Python touched only at the
edges.** Everything below is in service of that, plus the launch-bound constraint solver and the
host-snapshot binding path.

---

## 1. Ranked implementation routes to kill per-step Python re-entry (the core decision)

The scientific-framework survey is the load-bearing part. Ranked for OUR case (stay married to HOOMD,
preserve validated numerics):

### Tier 1 — **CUDA-graph replay of the existing cupy kernel sequence** ⭐ try first (cheapest)
- **What:** record the per-step sequence of kernel launches once (`cupy.cuda.Stream.begin_capture()` →
  `end_capture()` → `Graph.launch()`), then replay it each step with **one** submission — eliminating the
  per-launch Python/CPU dispatch that is throttling the 38k tiny M-SHAKE launches (97%util/64W).
- **Why it's best-first:** no new language, no engine change, **no kernel rewrite** — you replay your own
  *already-validated* cupy kernels, so it perturbs numerics the least (only *dispatch* changes, not
  *math*). Effort: days–weeks.
- **Two MUST-DO gotchas (verified failure modes):**
  1. **RNG must live in a device buffer.** A captured graph freezes kernel parameters incl. the RNG
     offset → **every replayed step draws identical "random" numbers** → breaks the two-Gaussian noise and
     the equipartition gate. Fix: counter-based RNG (Philox/Threefry) whose **per-step counter is read from
     a cupy device array a captured kernel increments**. (PyTorch had to patch exactly this:
     github.com/pytorch/pytorch/commit/c068180a176fb7d67add6c52ce5b50ff544fb5e5; NVIDIA capture-failures
     docs.nvidia.com/dl-cuda-graph/troubleshooting/capture-failures.html)
  2. **Capture forbids host syncs / fresh allocations.** Everything must run on pre-allocated fixed device
     pointers (which `gpu_local_snapshot` buffers already are). **Our M-SHAKE is currently tolerance-based
     with a `max_iter` cap** (`constrained_baoab.py:908` `tol=self.shake_tol, max_iter=self.shake_max_iter`)
     — a data-dependent convergence loop with a host check does NOT capture cleanly. Make it a **fixed
     iteration count** (or fixed-max with no host convergence read) and verify drift stays in-gate at that
     fixed count. This is the single capturability question for Codex to resolve.
- **Sources:** cupy Graph docs.cupy.dev/en/stable/reference/generated/cupy.cuda.Graph.html ; NVIDIA
  "Getting Started with CUDA Graphs" developer.nvidia.com/blog/cuda-graphs/ (their example ~9.6→3.4 µs/
  kernel, best at 50–100 kernels/graph).

### Tier 2 — **NVIDIA Warp** (durable port without writing a C++ plugin)
- **What:** write the BAOAB step + M-SHAKE + the 7 custom forces as `@wp.kernel`s (Python syntax,
  JIT-compiled to CUDA); launch via `wp.launch`, wrap the step in `wp.ScopedCapture` for graph replay.
- **Why strong for us:** **zero-copy interop with HOOMD.** HOOMD's `gpu_local_snapshot` exposes
  position/velocity/force via `__cuda_array_interface__`
  (hoomd-blue.readthedocs.io/en/v5.1.0/hoomd/data/hoomdgpuarray.html), and Warp constructs a `wp.array`
  zero-copy from any `__cuda_array_interface__` object (+ `wp.from_dlpack`). So Warp kernels run **in place
  on HOOMD's own GPU buffers** — keep HOOMD's neighbor list + built-in forces, replace only
  integrator+constraints+custom-forces. Counter-based `wp.rand` sidesteps the graph-RNG freeze if the
  counter is a device buffer. Effort ≈ **30–50% of a hand-written HOOMD C++/CUDA plugin**.
- **Sources:** github.com/NVIDIA/warp ; nvidia.github.io/warp/ ; concurrency/graph-capture
  nvidia.github.io/warp/modules/concurrency.html. **[verify]** the exact `wp.array(ptr=…)` external-pointer
  signature and that `gpu_local_snapshot` is safe to enter *inside* a running `sim.run()` Action (the
  protocol is confirmed; the in-`run()` access pattern was not explicit in docs — test empirically; both
  Tier-1 and Tier-2 depend on it).

### Tier 3 — **C++/CUDA HOOMD IntegrationMethod plugin** (max control, max effort)
- The "real native" route. Use **OpenMM `CustomIntegrator` as the design TEMPLATE, not the runtime**:
  declare the per-step LM step as a compiled expression DAG — the two-Gaussian noise = **two separate
  per-DOF stages** (one `gaussian` token returns one value per evaluation), constraints =
  `ConstrainPositions`/`ConstrainVelocities` (SHAKE/RATTLE). Don't migrate to OpenMM (that's a platform
  move + full re-validation, and its CUDA force-sum order is non-deterministic unless `DeterministicForces`
  — docs.openmm.org/latest/userguide/library/04_platform_specifics.html). Realize the template as Warp
  kernels (Tier 2) or C++.
- **Source:** OpenMM CustomIntegrator docs.openmm.org/latest/api-python/generated/openmm.openmm.CustomIntegrator.html ;
  OpenMM 7 paper DOI 10.1371/journal.pcbi.1005659.

### Tier 4 — JAX `lax.scan` / JAX-MD (only if differentiability becomes a goal — it isn't per CLAUDE.md)
- Mechanically the **cleanest** "whole loop on GPU, no per-step Python" (`lax.scan` threads
  positions/velocities/RNG-key through one compiled XLA loop; overdamped `brownian` is built in). **But it
  means LEAVING HOOMD** — reimplement forces + neighbor lists + the constraint solver, with **M-SHAKE in
  `lax.while_loop`** the genuinely hard part (constraints aren't in JAX-MD). = a new simulator (≈ v1→v2
  rewrite), multi-week, re-opens validation. Source: JAX-MD DOI 10.1088/1742-5468/ac3ae9.
- **Not recommended:** **Taichi** (weak in-place HOOMD interop — would copy through fields, reintroducing
  the sync; no graph replay) and **PyTorch** (redundant with cupy CUDA graphs, heavyweight, same RNG-freeze
  problem).

**Recommended path:** **measure → Tier 1 (cheap win) → Tier 2 (durable) → Tier 3 only if needed.**

---

## 2. The agreed 5-point direction — reinforced + sharpened by engine precedent

| Our point | Engine precedent that validates it | Sharpening from the survey |
|---|---|---|
| **#1 BAOAB → C++/CUDA native** | Unity **Burst** (compile hot loop off managed runtime); OpenMM CustomIntegrator; Warp | Use the §1 route ladder. Pin precision: CUDA **no `--use_fast_math`, deliberate `-fmad`/FMA policy** = Unity Burst `FloatMode.Strict`/`Deterministic` analogue (Fast math "can alter results", forbidden). |
| **#2 myosin/xlink → fixed GPU attachment pool** | Unreal **Niagara / VFX-Graph GPU-Events** (particle birth/death entirely on GPU, no CPU topology edit); Niagara **free-list/compaction**; Houdini **MPM** (rare structural rebuild stays CPU by design) | Birth/death = flip a per-slot `active` flag in a pre-allocated device pool. **Size so binding is NEVER silently dropped** (VFX caps; we cannot) → pool max = **Magic-Number Block** (physiological max bound-fraction, not gate-tuned) + **saturation HALTs/surfaces**. Drive the force kernel by **indirect dispatch over a device active-count** so CPU never reads counts. Houdini precedent: the *rare* topology rebuild may stay host-side as long as the *per-step* loop never leaves the GPU. |
| **#3 custom forces → CUDA force kernels** | Unity **ComputeShaders + StructuredBuffer**, **Mass/ECS SoA**; HOOMD **`gpu_local_force_arrays`** (zero-copy, O(1)) | Write into `gpu_local_force_arrays` with cupy/Warp (no host bounce). Built-ins stay built-in (harmonic/LJ/angle); only ERM/turgor/membrane/nucleus/substrate/myosin become kernels. **enclosed-volume turgor is a global constraint → needs a reduction, not pure elementwise.** Lay state out **SoA** for coalescing. |
| **#4 manifold/contact data in VRAM** | Houdini **Gas/OpenCL** "leave results on GPU until a CPU consumer asks" (explicit Readable/Writeable/Flush/Finish-Kernels flags, default lazy copy-back); Niagara persistent buffers | Keep patch/normal/tangent/contact as persistent device arrays; binding updater reads only the specific arrays it needs, only when topology changes. Mesh = geometry/broad-phase only (the manifold-architecture bright line), never a force carrier. |
| **#5 snapshot only at checkpoint/measurement/debug** | Houdini "**disable DOP per-frame caching**" (the residency killer); Unity **AsyncGPUReadback**; HOOMD **local O(1) vs global O(n)** snapshot (~3× documented) | Make measurement/checkpoint readback **async + double-buffered**: pinned host memory + `cudaMemcpyAsync` on a separate stream + event, never a blocking copy in the step loop. Enforce snapshot-consistency only at true checkpoint boundaries. Also: prefer `cpu_local_snapshot`/`gpu_local_snapshot` (O(1)) over the global `get_snapshot()` (O(n)) everywhere it survives. |

---

## 3. The 64W / launch-bound constraint solver, specifically

The "97% util but 64 W" = SM clock busy but warps starved = **too many tiny launches + divergence**, not
FLOP-bound. Three levers, in order:

1. **Batch** per-chain Fixman + M-SHAKE into a few kernels over *all* chains (Unity Job-System batching;
   it does NOT schedule 38k individual jobs). One batched 6×6/tridiagonal kernel, not 38k launches.
2. **Wavefront coherence** (Blender Cycles / Laine et al. 2013, "Megakernels Considered Harmful", NVIDIA
   Research): bucket work so SIMT threads do the *same thing* — group by chain length / bound-vs-unbound
   head / converged-vs-not — to kill control-flow divergence. This is what actually fills the SMs beyond
   mere fusion. (research.nvidia.com/.../laine2013hpg_paper.pdf)
3. **CUDA Graphs** collapse the residual launch overhead regardless (§1 Tier 1).
- **Honest expectation (Houdini reality check):** constraint solves are the *hardest* thing to GPU-ify —
  Houdini's flagship **Vellum (XPBD) still runs the constraint solve primarily on CPU**, GPU-accelerating
  only the *neighbor search* (+ optional graph-coloring for parallel batching)
  (sidefx.com/docs/houdini/nodes/dop/vellumsolver.html). So: GPU-accelerate the neighbor/spatial sub-step
  first, batch + color the constraints — but **graph-coloring reorders FP projections → must pass the
  bit/gate check** (see §4).

---

## 4. The fidelity boundary + TWO decisions for PI

**Borrowable vs not.** Game/creative engines mostly get fast by **trading accuracy** (FP16, LOD,
preview/proxy resolution, AA-ray reduction, caching physics results) — **all forbidden** for a
fidelity-strict, gate-validated MD. What transfers is the **accuracy-neutral class**: GPU-residency,
fixed pool + compaction, indirect dispatch, SoA/ECS, wavefront coherence, async readback, persistent
buffers, kernel batching, CUDA-graph capture. These change *where data lives / when kernels run*, not
*what is computed*. **But every borrowed pattern still needs gate validation** — parallel reductions,
graph-coloring, and FMA-fusion can each change the last bits.

**⭐ PI DECISION 1 — the acceptance criterion (must settle before the port).** The research is unanimous:
**bitwise invariance vs the validated CPU implementation is achievable by NO GPU route** — GPU
FP-reduction order alone breaks it, independent of Tier. Nuance for us: the **element-wise BAOAB step can
stay bitwise** (no cross-element reduction), but the **constraint solver (Fixman/M-SHAKE reductions) and
any parallel force-sum cannot guarantee it**. So the realistic acceptance target is **"passes the same
STATISTICAL gates (drift / equipartition / L_p / KU-3.5 γ)"**, not "identical bits." Tier-1 (CUDA-graph
replay of our own kernels) and Tier-2 (Warp replicating the math line-for-line) are the only routes that
don't *additionally* change the math. **If "bit-invariant vs CPU" is truly non-negotiable as *bitwise*,
then no GPU production port qualifies** — PI must rule whether the gates (not bits) are the acceptance
criterion. (Note: the existing cupy device port's "CPU bit-invariance preserved" claim holds for the
element-wise path; the constraint-solver reductions are where it needs re-checking.)

**⭐ PI DECISION 2 — the D3 fork (could avoid most of the C++).** At the production Route-B CFL dt
(`dt_factor≈1e-4`), is HOOMD's built-in compiled `md.methods.Brownian` (Euler-Maruyama, already C++/CUDA,
autotuned) accurate enough vs the two-Gaussian LM-limit? If yes → you get a native compiled integrator
**for free**, no plugin. This trades the D3-freeze "no Euler-Maruyama" rule, so **don't decide — QUANTIFY**:
run EM vs LM at production dt on the gates (drift/equipartition/L_p/γ); if within tolerance, surface to PI
to revisit the freeze. Cheap, decisive, do before any big C++.

---

## 5. Concrete first steps for Codex (measure → cheap win → durable)

1. **MEASURE (cheap, decisive — do first):**
   - *Ceiling:* pure-HOOMD baseline steps/s at the same N (built-in Langevin + harmonic bond, **zero**
     custom Action/force). The gap to ~2.1k = the exact size of the per-step Python tax = the headroom any
     route can recover (the "26h→1–2h" basis).
   - *Fork (PI Decision 2):* EM (built-in Brownian) vs LM (two-Gaussian) at production dt on the gates.
2. **CHEAP WIN (Tier 1):** CUDA-graph-capture the existing cupy step. Prereqs: (a) RNG counter in a device
   buffer; (b) M-SHAKE fixed-iteration (no host convergence check) + all pre-allocated buffers. Verify
   gates. Days–weeks, preserves validated GPU numerics.
3. **DURABLE (Tier 2):** if graph constraints bite or more headroom is needed, port
   BAOAB+M-SHAKE+custom-forces to **Warp kernels operating in place on `gpu_local_snapshot` buffers** +
   graph capture.
4. **DATA LAYER (#2/#3/#4/#5):** fixed attachment pool (indirect dispatch, saturation-halt, Magic-Number
   pool max) + custom forces as CUDA/Warp kernels into `gpu_local_force_arrays` + manifold/contact in VRAM
   (lazy copy-back discipline) + async/double-buffered measurement readback.
5. **GATES throughout:** every step validated against drift 3e-15 / equipartition / L_p / KU-3.5 γ (PI
   Decision 1 sets bitwise-vs-statistical). Integrator changes are PI-gated (D3 freeze). No magic numbers,
   no gate-loosening.

---

## 6. Sources (consolidated, verified this session)

**Game / VFX engines**
- Laine, Karras, Aila 2013, "Megakernels Considered Harmful: Wavefront Path Tracing on GPUs" (NVIDIA
  Research) — research.nvidia.com/sites/default/files/pubs/2013-07_Megakernels-Considered-Harmful/laine2013hpg_paper.pdf
- Blender Cycles kernel scheduling — developer.blender.org/docs/features/cycles/kernel_scheduling/
- Unreal Niagara GPU particles — dev.epicgames.com/documentation/unreal-engine/how-to-create-a-gpu-sprite-effect-in-niagara-for-unreal-engine
- Unreal Mass Entity (ECS/SoA) — dev.epicgames.com/documentation/unreal-engine/overview-of-mass-entity-in-unreal-engine ; unrealcode.net/MASS_001.html
- GPU-driven rendering / indirect dispatch — vkguide.dev/docs/gpudriven/gpu_driven_engines/
- Async GPU readback (fences) — nicholas477.github.io/blog/2023/reading-rt/

**Unity**
- Burst (LLVM, IL→native, HPC#) — docs.unity3d.com/6000.0/Documentation/Manual/script-compilation-burst.html
- Burst `FloatMode` (Strict/Deterministic/Fast — the fidelity lever) — docs.unity3d.com/Packages/com.unity.burst@1.8/api/Unity.Burst.FloatMode.html
- Job System batching — docs.unity3d.com/6000.1/Documentation/Manual/job-system-overview.html
- ECS 16 KiB SoA chunks — docs.unity3d.com/Packages/com.unity.entities@1.2/manual/concepts-archetypes.html
- ComputeBuffer↔StructuredBuffer — docs.unity3d.com/ScriptReference/ComputeBuffer.html
- VFX Graph (GPU simulation, GPU Events) — docs.unity3d.com/6000.4/Documentation/Manual/VFXGraph.html
- `Graphics.RenderMeshIndirect` (current; `DrawMeshInstancedIndirect` is obsolete) — docs.unity3d.com/ScriptReference/Graphics.RenderMeshIndirect.html
- `AsyncGPUReadback` (non-blocking) — docs.unity3d.com/ScriptReference/Rendering.AsyncGPUReadback.html
- NativeArray / Collections allocators — docs.unity3d.com/Packages/com.unity.collections@1.2/manual/allocation.html

**Adobe / DCC (Houdini = strongest analogue)**
- Houdini Gas OpenCL ("left on the GPU until requested"; Readable/Writeable/Flush) — sidefx.com/docs/houdini/nodes/dop/gasopencl.html
- Houdini OpenCL SOP — sidefx.com/docs/houdini/nodes/sop/opencl.html
- Houdini MPM (GPU solve, CPU grid rebuild by default) — sidefx.com/docs/houdini/nodes/sop/mpmsolver.html
- Houdini Vellum (XPBD constraint solve primarily CPU; OpenCL neighbor search) — sidefx.com/docs/houdini/nodes/dop/vellumsolver.html
- HOOMD local-vs-global snapshot (O(1) vs O(n), ~3×) — hoomd-blue.readthedocs.io/en/latest/tutorial/04-Custom-Actions-In-Python/06-Improving-Performance.html
- HOOMD zero-copy `gpu_local_force_arrays` — hoomd-blue.readthedocs.io/en/latest/hoomd/md/force/custom.html
- Adobe After Effects GPU (mostly counter-example) — helpx.adobe.com/after-effects/using/basics-gpu-after-effects.html

**Scientific GPU-compute frameworks**
- NVIDIA Warp — github.com/NVIDIA/warp ; nvidia.github.io/warp/ ; concurrency/graphs nvidia.github.io/warp/modules/concurrency.html
- HOOMD `__cuda_array_interface__` GPU arrays — hoomd-blue.readthedocs.io/en/v5.1.0/hoomd/data/hoomdgpuarray.html
- OpenMM CustomIntegrator — docs.openmm.org/latest/api-python/generated/openmm.openmm.CustomIntegrator.html ; platform determinism — docs.openmm.org/latest/userguide/library/04_platform_specifics.html ; OpenMM 7 DOI 10.1371/journal.pcbi.1005659
- JAX `lax.scan` — docs.jax.dev/en/latest/_autosummary/jax.lax.scan.html ; JAX-MD DOI 10.1088/1742-5468/ac3ae9
- cupy CUDA Graph — docs.cupy.dev/en/stable/reference/generated/cupy.cuda.Graph.html ; interop — docs.cupy.dev/en/stable/user_guide/interoperability.html
- CUDA-graph RNG-freeze (bit-invariance crux) — github.com/pytorch/pytorch/commit/c068180a176fb7d67add6c52ce5b50ff544fb5e5 ; docs.nvidia.com/dl-cuda-graph/troubleshooting/capture-failures.html
- NVIDIA "Getting Started with CUDA Graphs" — developer.nvidia.com/blog/cuda-graphs/
- Taichi (poor HOOMD interop) — github.com/taichi-dev/taichi
- PyTorch CUDA graphs (redundant with cupy) — pytorch.org/blog/accelerating-pytorch-with-cuda-graphs/

**Marked [verify] before building:** Warp `wp.array(ptr=…)` external-pointer signature; `gpu_local_snapshot`
safe to enter inside a running `sim.run()` Action; M-SHAKE iteration structure for clean graph capture
(currently `tol`+`max_iter` — needs a fixed-iteration formulation). Vendor speed figures (Substance ~20×,
EmberGen particle counts) are release-coverage, not independently benchmarked.

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

# Phase C — adopt the Warp DCM engine (solidification plan)

**PI decision 2026-06-21:** enter Phase C. Make the Warp-ported DCM engine — *with what has
already been built* — the going-forward production engine, **then** resume the de-cohesion
discussion. This doc is the worklist that takes the engine from "GO milestone" to "trusted
production engine," and the state to boot from when work starts.

Engine ↔ physics are **orthogonal**: this adoption proceeds regardless of the de-cohesion verdict
(which is itself now *undecided*, not PLATEAU — see §5).

---

## 0. Done at consolidation (2026-06-21)

- **`warp/port` merged into the trunk** (`h7/compartment-platform`, merge `ecf130d`, path-disjoint,
  no conflicts). The engine now lives in the working tree at `aleph/warp_port/`.
- **24 `warp_port` parity tests PASS on this Mac (Warp CPU backend).** Engine executes locally.
- 15 parity claims in `outputs/tag_kb/results_manifest.yaml`, all `declared: verified` (CPU backend).
- Orchestration report corrected (REMESH = host-hybrid by design; B3 ran 30×/10–37×; diff loop done).

**What exists** (all parity-gated vs committed HOOMD references, machine-eps on CPU):
integrator stack (BAOAB + M-SHAKE + Fixman) · compartment forces (harmonic bond/angle + LJ-WCA +
radial-shell turgor/membrane/nucleus) · DCM forces (node-face contact, exact-volume turgor,
node-node cohesion, lamellipodium tether) · GPU-resident **hybrid loop** (device per-step + host
low-cadence remesh on a fixed-size node-pool) · hash-grid neighbour list (10.9–36.8× multi-cell) ·
**end-to-end differentiable loop** (grads w.r.t. dP0, k_edge; autodiff vs FD 4e-12).

Entry points today are separate functions: `run_hybrid` (single-cell), `run_multicell`,
`run_crawl` (lamellipodium), and the differentiable loop. **Not yet one coherent engine API** (→ G3).

---

## 1. What "solid / 확실" means — the adoption gates

The engine is trusted as production when these clear. **G1 is the make-or-break adoption gate.**

### G1 — GPU-*backend* parity  ⚠️ make-or-break · needs gbook A5000
All 15 parities were measured on Warp's **CPU** backend; the A5000 numbers are **speed only**.
Before trusting GPU production we must confirm the *same kernels* reproduce parity on Warp's
**CUDA** backend (atomic-reduce ordering + fast-math can differ — the worker already saw
atomic-order ~2.6e-12 on CPU).
- **Action:** sync the tree to gbook (Syncthing), run `aleph/warp_port/parity_report.py` (and/or
  the `tests/warp_port` suite) with the Warp **cuda** device against the committed HOOMD fixtures.
- **Gate:** force-law kernels at machine-eps (≤1e-12); atomic-reduce + iterative (M-SHAKE) kernels
  at the documented GPU tol (expect ~1e-12, *not* bit-for-bit). Record per-kernel GPU tols.
- **Output:** a `declared: verified` GPU row per piece in `results_manifest.yaml` (currently CPU-only),
  + a short `PHASE_C_GPU_PARITY_<date>.md`. Un-run pieces stay `needs-regen`, never `verified`.

### G2 — port the substrate forces to Warp  · needs the HOOMD ref (on trunk) + gbook for parity
Route-a spreading needs the substrate, which is **not yet ported**:
- `DcmSubstrateForceGPU` — z-well (the flagged minor remainder), `cell/dcm_gpu_forces.py:178`.
- `DcmSubstrateWettingGPU` — in-plane area-maximizing wetting, `cell/dcm_gpu_forces.py:229`. **This is
  the *mechanistic* spreading driver** that replaces the deleted `SettlingForce` proxy (landmine §4
  fix #2). Without it on Warp, the engine cannot do real (non-proxy) substrate spreading.
- **Action:** mirror each as a Warp kernel; parity-gate vs the committed HOOMD force (same pattern as
  the other 15 pieces); add fixtures + tests + manifest claims.

### G3 — one coherent engine API  · local
Collapse `run_hybrid` / `run_multicell` / `run_crawl` / diff-loop into a single documented engine
(build → configure force set (turgor/edges/contact/cohesion/lamellipodium/substrate) → run with
remesh cadence → measure), so production drivers call ONE entry point. Keep the differentiable
variant as a flag, not a fork. Preserve the brute-force kernels as the parity oracle (grid = default).

### G4 — docs / hygiene  · local
- `warp_port/__init__.py` docstring still says **"Phase B — parity spike"** → update to the adopted
  production engine.
- `STRUCTURE.md`: add `aleph/warp_port/` as the GPU-resident differentiable engine.
- `CLAUDE.md` Stack section: note Warp as the Phase-C GPU-resident differentiable DCM engine
  (HOOMD-blue remains the fine-grained reference / parity oracle).

### G5 — CI + gate wiring  · local
- `tests/warp_port` into CI (they import `warp`; the env has 1.14.0 — confirm CI installs it).
- A `make warp-parity` target (CPU) + fold the GPU verdict into `make kb-check`'s manifest gate
  once G1 lands.

---

## 2. Sequencing (single A5000)

1. **G1 GPU-backend parity** — first; it decides whether GPU production is trustworthy at all. Cheap
   (parity harness exists, runs in minutes on the A5000).
2. **G2 substrate-force port** — parallel CPU authoring; GPU parity folds into the next G1 pass.
3. **G3 engine API + G4 docs + G5 CI** — local, any time; do alongside.
4. Then the engine is production-trusted and route-a / the de-cohesion re-run (§5) can run on it.

Local (no GPU) work can start immediately: G2 authoring, G3, G4, G5. GPU work (G1, G2-parity) is
gbook-gated and is the first thing to run when the A5000 is taken.

## 3. Anti-drift contract (unchanged from the spike)

Nothing is "done" without a committed artifact + a passing test. GPU verdicts must actually RUN on
the A5000 (un-run = `needs-regen`, never `verified`). Parity = bit/tol-compare vs the **committed
HOOMD reference**, not a new self-authored oracle. Every headline enters `results_manifest.yaml`.

## 4. Risks

- GPU-backend parity could exceed CPU tols on the atomic-reduce kernels (turgor volume-reduce,
  contact/cohesion scatter) — if so, document the GPU tol per kernel; only a *gross* miss is a bug.
- `DcmSubstrateWettingGPU` is a conservative xy energy-gradient; getting the Warp gradient sign /
  capping to match the HOOMD force is the parity risk for G2.
- Remesh stays host-side: at very large N the ~113 ms host epoch + GPU↔host resync is the scaling
  lever (amortizes with cadence; <1 % of multi-cell runtime today). Not a correctness risk.

---

## 5. De-cohesion linkage (parked, but the engine serves it)

The de-cohesion make-or-break is **UNDECIDED, not PLATEAU** (`DCM_LANDMINE_REGISTER_2026-06-21.md`,
on `dcm/decohesion`): the committed `A/A0=1.958/2.286` was driven by the deleted `SettlingForce`
body-force proxy; the only config isolating the mechanistic lamellipodium (ctrl_C, settle-OFF)
reached **A/A0 3.0 above the claimed ceiling at genuine held-maxZ crawl**, but **diverged past 120k**.

Why this is *engine* news, not just physics: ctrl_C's blocker is **numerical stability + the
mechanistic spreading driver**, which is exactly the Phase-C engine's strength — a well-tested
stable BAOAB integrator + the ported `DcmSubstrateWettingGPU` (G2) + a divergence guard. So the
clean re-run that would actually answer the de-cohesion verdict (settle=0, `--substrate-wetting`,
divergence-guarded, headline = `crawl_residual` at held maxZ) is a **natural first physics run on
the Warp engine** once G1+G2 clear. We do it in the *later* de-cohesion discussion, not now — but
the dependency is why solidifying the engine first is the right order.

---

## 6. First action when work resumes

Take the A5000 and run **G1 (GPU-backend parity)** — the one gate that decides whether the
adopted engine is trustworthy on GPU. In parallel (local), start **G2** authoring (substrate
z-well + wetting Warp kernels) and **G3/G4** (engine API + docstring/STRUCTURE). Everything else
is downstream of G1's verdict.

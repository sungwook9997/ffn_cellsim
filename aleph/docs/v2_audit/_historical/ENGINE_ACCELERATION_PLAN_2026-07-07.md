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

# Engine Acceleration Plan — FF + DCM (2026-07-07)

**Method.** 9-subsystem parallel hot-path audit + a paired adversarial verifier per
subsystem (18 agents). Every proposed optimization was re-checked for (a) *realistic* speedup
(Amdahl, bandwidth-bound vs compute-bound, launch-latency vs GPU-compute), (b) fidelity safety
under the project hard rules (no gate-loosening, no fidelity shortcut without a PI-signed gate,
bit-parity tests must still pass), and (c) whether it targets a real production path or dead/dev code.

## The two load-bearing conclusions

1. **The fast machinery already exists; it is mostly wired OFF or bypassed.** The largest
   defensible wins are *turning on default-OFF features* and *routing callers to the fast path
   that is already built* — not new code:
   - DCM already runs `device_cg` (dcm_warp_decohesion.py:1209) and has a **frozen-neighbour
     cache** and an **analytic-diagonal preconditioner** BUILT but defaulted **OFF**
     (dcm_warp_decohesion.py:549-559).
   - FF spreading/crawl/membrane already use the device-assembled CSR CG (`ff_implicit_step_gpu`);
     the γ-floor production settle already runs the device explicit relax (`relax_on_device`).
   - The remaining *host* implicit path (matrix-free FD-JVP CG, per-iter GPU→CPU round-trip) is a
     **secondary/test path**, not the production hot loop. So the "2400 syncs per equilibration"
     framing from the earlier study is real but hits the non-production path.

2. **Profile on the A5000 at production N BEFORE writing any code.** Multiple "obvious" wins
   contradict measurements already in the tree, and several hinge on whether a stage is
   GPU-compute-bound (Amdahl-caps sync removal) or launch/latency-bound (favours it) — which
   **cannot be determined on the dev Mac** (Warp CPU only). Specifically:
   - The diagonal preconditioner is asserted a win by the design pass but the in-code note
     (dcm_warp_implicit.py:63-68) says "Contact does NOT dominate the conditioning here … plain CG
     ~3-5 iters," and the Hutchinson variant **HURT** (3→157 iters). Measure with the in-tree
     `_precond_demo` at the production stack before enabling.
   - The persistent/preallocated contact grid was already tried and is "an honest negative"
     (dcm_warp_decohesion.py:973-974) — the path is query-bound, not build-bound. Do not re-chase.
   - At native N the force kernels already saturate the A5000, so ensemble-batching and
     launch-overhead (CUDA-graph) wins shrink to ~1.1-1.5×; the real lever at native N is **CG
     iteration count** (preconditioning) and **per-iter sync removal**.

**Sequencing.** Phase 0 = profile (A5000, production N, per-stage). Phase 1 = the exact / low-risk
default-on + routing wins. Phase 2 = gated approximate wins (each needs a named ValidationGate).
Phase 3 = the structural rewrite (DCM analytic Hessian CSR). Do not start Phase 2/3 before Phase 0
confirms the stage is actually hot.

---

## FF engine

### P0 — highest leverage, low risk
| Item | Where | Realistic | Fidelity |
|---|---|---|---|
| **Jacobi (diagonal) preconditioner on the cupy CG** in `ff_implicit_step_gpu` | implicit_ff.py:253 (cg has no `M=`) | **~1.5-3× fewer CG iters** on the CG-dominated spreading/crawl/membrane drivers (CG ~2400 iters/equil = the actual bottleneck) | Exact to solver tol (preconditioning changes only the path). One-time equivalence check (crawl COM + γ to tol). **No new gate.** Pull `K.diagonal()` incl. the `k_vol·g·gᵀ` term before wrapping the LinearOperator. Does **not** help γ-floor `method='device'` (that path is explicit `relax_on_device`, no CG). |

### P1 — high value, contingency/gate
- **Route `relax_implicit` / `gamma_floor(method='implicit')` off the host FD-JVP CG onto the
  analytic device CSR path** (relax.py:87-96 → implicit_ff.ff_implicit_step_gpu). ~5-30× on that
  (secondary) equilibration path *and* it kills a documented-degenerate FD path at FF scale.
  **adopt-with-gate** — gate against **analytic ground truth** (bending κL/2R², tension), *not*
  self-comparison to the current FD output (which may be divergent). Per `[[feedback-oracle-is-crosscheck-not-truth]]`.
- **Call the bit-parity device `reshape_kernel` from the host relax paths** (constraints.py reshape
  → network_warp.py:224/258). Exact (bit-parity tested). *Caveat:* `reshape_np` currently re-uploads
  device topology arrays every call — hoist `fiber_off/seg_off/seg_rest` to persistent device arrays
  or the swap is a wash at small N.
- **Batched block-diagonal γ-floor ensemble at small N** (gamma_floor.py:564-565, serial list-comp
  for `method='device'`). **3-8× on the n_fil~100 mesoscale sweep** where each realization badly
  underfills the A5000. Exact (disjoint index blocks → identical per-realization result). NOT at
  native 38k (single cortex already saturates); host-RNG KMC + per-fiber reshape cap the win, so
  batch those too. Verify no shared HashGrid across the concatenated block.

### P2 — secondary
- **FP32 streaming force kernels, FP64 positions/accumulators** (bending/link/myosin/turgor).
  ~1.5-2× on those bandwidth-bound sweeps (FP64≈2× bytes, *not* the 64× compute penalty — they're
  memory-bound). **adopt-with-gate** — **CORRECTNESS TRAP specific to FF:** the headline metric
  `γ_myo = 0.0199·f_myo` is a delicate near-cancellation off a shell already ~500× under band;
  FP32 force truncation integrated over thousands of overdamped steps can perturb exactly the small
  tension being measured. Gate on **γ_myo tolerance vs the FP64 reference**, not just max-force error.
  Lean *investigate before adopt*.
- **Persistent `pos_d/force_d` buffers + drop the redundant `wp.synchronize_device`** in the host
  force wrapper (forces_warp.py:110-118). Exact, minor; mostly moot if the host-CG caller is removed
  by the P1 routing item.

---

## DCM engine

### P0 — highest leverage, low risk
| Item | Where | Realistic | Fidelity |
|---|---|---|---|
| **Fused device-scalar CG reduction** (α/β/rr in device scalars, drive with a 1-thread kernel; read residual on a cadence) in `device_cg` | dcm_warp_implicit.py:163-165 (`dot()` = `synchronize_device` + `.numpy()`, 3-4 dots/iter) | **~3-5× on the implicit step** (Amdahl ceiling ~8× from the attested 88% sync wall; matvec grid work + periodic residual readback cap it) | Device-scalar α/β = **exact** (bit-parity-checkable vs host-scalar CG). **KEEP the guards**: the pAp_diag floor (:214), the divergence break (:229 `rr_new>4·rs0`), and the degenerate-direction break (:196) exist because this operator misbehaves near equilibrium/large-dt. Read residual **every few iters**, never a fixed capture that drops the divergence exit. The residual-cadence sub-variant needs its own gate. |
| **Turn ON the frozen-neighbour cache** (`frozen_neighbors=True`) | dcm_warp_decohesion.py:549-556 (default OFF); cache dcm_warp_frozen.py | **~1.5-3× non-IPC / ~1.2-1.5× IPC** (removes the per-matvec hash-grid query + closest-point recompute) | **adopt-with-gate** — run `_frozen_parity_demo` **at production N/subdiv** and confirm **zero int32 overflow-clip** (frozen.py:415/431 silently drops candidates above the flat-array ceiling for N≥1000). This is the "neighbour skin that misses a contact" trap — verify before defaulting on. |
| **Port the `_edge_extent` device skip-gate into `do_remesh`, *with* a min-face-quality device reduction** | dcm_warp_decohesion.py:776 (full-N host download per remesh epoch); gate mechanics dcm_warp_hybrid_multicell.py:86-98 | Eliminates one full-N host download + Python classification pass per remesh epoch (Amdahl-capped by remesh cadence) | Exact **only** with the quality reduction: `remesh_pass` fires SWAP on `face_quality<0.2` (dcm_remesh.py:271), so an edge-extent-only gate (as the hybrid path does it) would silently skip a legitimate SWAP → different topology = a fidelity shortcut. Include both reductions. |

### P1 — high value
- **Per-cell remesh** — slice the pooled arrays by `face_cell`, remesh only cells flagged
  out-of-band (dcm_warp_decohesion.py:777-779). O(faces-of-affected-cells) instead of O(total ~128k);
  big win during localized spreading, ~1× when all cells remesh at once. Exact **iff** the pooled
  build gives each cell a disjoint node range (contact is force-mediated, not shared vertices) —
  **verify in the builder first**, then scope `max_ops` per cell.
- **Move the unconditional cadherin `wp.synchronize_device` into the host-only (`else`) branch**
  (dcm_warp_decohesion.py:1370, :1429). On the GPU cadherin path the drain is pure waste. Small
  (~one drain / 50 steps) but exact and free; verified safe (update_gpu is device-only).
- **DCM device-assembled analytic sparse-K CSR + on-device SpMV CG** (dcm_warp_implicit.py:202-203)
  — the biggest *structural* win (removes the per-iter FD force sweep + its FD-noise floor), but
  **effort XL, `investigate` not adopt**: DCM needs analytic Hessians the FF path never had —
  face-pressure turgor (not pure rank-1), IPC log-barrier (non-smooth, must match `ipc_hess`
  exactly), node-node cohesion, centroid-coupled nucleus, bilaplacian bending. `build_diagA` already
  gives the analytic **diagonal**, which de-risks a preconditioner but not the off-diagonal coupling.
  Each term needs a K·p-vs-FD-JVP equivalence fixture. This is a Phase-3 project, gated by Phase-0
  profiling showing the FD sweep is actually the hot fraction (it's already frozen-cache-accelerated).

### P2 — secondary
- **Analytic-diagonal Jacobi preconditioner ON** (`precond_diag=True`, dcm_warp_frozen.py). **MEASURE
  with `_precond_demo` at the production stack before enabling** — the in-tree evidence says it is
  likely a no-op near the operating point (γ/dt regulariser dominates; plain CG ~3-5 iters). Exact
  (no gate) but do not enable on the design claim alone.
- **Retire `JunctionSwitchHost` → the GPU-native E1 cadherin ensemble as the production default**
  (already disabled under `--cadherin`, dcm_warp_decohesion.py:632). Removes an O(n_cells²) CPU
  distance matrix + full-N D2H. Not a hot-path speedup in cadherin runs (already zero there), but it
  **resolves the mechanistic-over-lumped violation** (`[[feedback-junction-switch-fine-grained]]`) —
  right direction, adopt as default; keep the latch only as a legacy non-cadherin option.
- **CCD device atomic-min + cache pass-1 barycentric coords** (dcm_contact_implicit_warp.py:334, :74-86).
  Removes one host sync/step + a redundant closest_bary recompute. Exact, minor (secondary to P0).
- **Mixed precision on the neighbour force math** (FP32 cull/geometry, FP64 accumulate). ~1.5-2× on
  those kernels, but **adopt-with-gate** and **low priority**: under the `[DCM=meters]` convention
  node positions are ~7.5e-6 m, where FP32 (~7 sig figs) resolves only ~1e-12 m against a thin
  contact shell with rep~2e8 — FP32 positions can miss penetration and de-condition CG. Keep
  positions + the CG operator + contact FP64 unconditionally.

---

## Do NOT chase (verified dead ends / traps)
- **FP32 in the CG operator, positions, or contact geometry** — de-conditions the solve (stall/diverge)
  and misses near-contact events. FP32 only ever for isolated soft explicit accumulation, behind a gate.
- **Persistent/preallocated contact grid** — already measured negative (query-bound, not build-bound).
- **Kernel fusion of the force kernels** — different launch domains (triples/edges/nodes), memory-bound,
  non-associative FP64 → not exact, high rewrite risk, marginal.
- **CUDA-graph capture of the CG loop** — forces a fixed iteration count and removes the divergence /
  pAp-floor guards that were added to stop large-dt blow-up. Only capture *fixed-topology force+integrate*
  windows on the explicit path, and only after Phase-0 shows launch overhead is a real fraction.
- **Porting the crowd-pressure junction latch to a device kernel** — optimizing a superseded lumped
  mechanism; the sanctioned spend is E1 cadherin (already on-device).
- **Reuse/freeze K (quasi-Newton)** — targets assembly (cheap), not CG iters (dominant), and a stale
  linearization typically needs *more* outer steps. Prefer the preconditioner.
- **Fixman/SHAKE vecL sizing, FFImplicitStepper LU** — off the production path / zero callers.
- **Cadherin grid-reuse / graph-capture** — config-fragile (µm-scale `--cad-rbind` > coh_q) and ~zero payoff.

## Fidelity-gate ledger (new PI-signed ValidationGates required before the gated items go to production)
| Gate | For | Acceptance |
|---|---|---|
| Frozen-cache parity @ production N | DCM P0 frozen cache | `_frozen_parity_demo` rel<1e-5 to cg_tol **and** zero overflow-clip at production N/subdiv |
| CG residual-cadence | DCM P0 device-scalar CG (read-every-K sub-variant) | converged dx within cg_tol vs sync-every-iter; divergence guard still fires |
| Remesh trajectory | DCM P1 per-cell / incremental remesh | closed-manifold + Euler + volume invariants **and** fixed-seed trajectory divergence bound |
| FP32 force-parity + γ_myo | FF P2 mixed precision | max rel force error **and** end-of-run γ_myo vs FP64 reference within tol |
| Mixed-precision A/A0 | DCM P2 neighbour FP32 | A/A0, V/V0, penetration_frac invariants vs FP64 over settle+spread |
| KMC statistical equivalence | device-RNG KMC ports | ensemble mean/std of γ_active + turnover rate vs host-RNG (distributional, not bit) |

## Phase 0 profiling checklist (run first, on the Gbook A5000)
- FF: γ-floor production (native ~70k filaments) + one spreading/crawl run — per-stage wall
  (force sweep / CG iters-per-step / reshape / turgor refresh / host syncs). Confirm CG iter count
  is the fraction the plan assumes.
- DCM: n400 full-compartment decohesion step — `device_cg` sync wall %, hash-grid query vs build %,
  remesh epoch cost, contact matvec share. Run `_precond_demo` + `_frozen_parity_demo` at this N.
- Decision: only start P0 items whose target stage Phase-0 confirms as ≥~15% of wall.

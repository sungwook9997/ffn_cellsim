# Phase C · G1 — GPU-backend parity verdict

**Gate:** the Warp DCM engine's adoption make-or-break (`_historical/PHASE_C_WARP_DCM_ENGINE_PLAN_2026-06-21.md` §G1).
**Question:** the 15 ported pieces were parity-gated on Warp's **CPU** backend; do the *same kernels*
reproduce parity on the Warp **CUDA** backend (gbook RTX A5000)? — required before GPU production is trusted.

## Verdict: ✅ PASS — 15/15 on `cuda:0`

Run: `python aleph/warp_port/parity_report.py --device cuda:0` on gbook (A5000 Laptop, sm_86,
CUDA 12.9 / driver 13.2, warp 1.14.0), against the **committed HOOMD reference fixtures** (frozen
on the Mac; the agent does not grade its own homework). Result: `fixtures/warp_parity_results_cuda0.json`.

| piece | GPU metric | value | gate | pass |
|---|---|---|---|---|
| B1 BAOAB kT=0 | max_abs_pos_diff | 1.72e-15 | < 1e-9 | ✅ |
| B1 BAOAB kT=1.5 | max_abs_pos_diff | 2.00e-15 | < 1e-7 | ✅ |
| B2 radial nucleus | force_rel | 1.01e-15 | < 1e-12 | ✅ |
| B2 radial membrane | host force_rel / warp-reduce | 2.50e-16 / 1.75e-12 | < 1e-12 / < 1e-8 | ✅ |
| B2 radial turgor | host force_rel / warp-reduce | 3.62e-14 / 1.31e-12 | < 1e-12 / < 1e-8 | ✅ |
| M-SHAKE chain constraint | pos / lam / bond_rel; nonconv | 7.37e-12 / 4.19e-10 / 4.92e-11; 0 | < 1e-9 | ✅ |
| Fixman metric force | force_rel / U_rel | 6.15e-16 / 0.0 | < 1e-9 | ✅ |
| compartment harmonic bond | force_rel / energy_rel | 2.87e-16 / 0.0 | < 1e-12 | ✅ |
| compartment harmonic angle | force_rel | 6.10e-16 | < 1e-12 | ✅ |
| compartment LJ-WCA | force_rel | 2.39e-16 | < 1e-12 | ✅ |
| DCM node-face contact | force_rel | 1.64e-16 | < 1e-12 | ✅ |
| DCM turgor force | host force_rel / warp-reduce | 3.25e-16 / 2.28e-15 | < 1e-12 / < 1e-8 | ✅ |
| DCM cohesion force | force_rel | 3.44e-16 | < 1e-12 | ✅ |
| lamellipodium tether | force_rel | 0.0 (bit-for-bit) | < 1e-9 | ✅ |
| B4 differentiability | grad rel_err vs analytic | 2.60e-15 | clean | ✅ |

## Reading

- **Force-law (host-reduce) parities are machine-eps on GPU** (1e-16 … 1e-14) — identical-class to CPU.
- **Atomic-reduce kernels** (radial membrane/turgor warp-reduce, DCM turgor warp-reduce) land at
  **~1.3–2.3e-12 on GPU**, far under the 1e-8 reduction-order threshold. The worker's flagged
  "atomic-order ~2.6e-12" concern is real but **6 orders inside the gate** — not a correctness issue.
- **M-SHAKE (iterative)** returns the *same* deterministic tol as CPU (pos 7.4e-12 / lam 4.2e-10 /
  bond 4.9e-11, nonconverged 0) — the Thomas solve is backend-stable.
- **Only meaningful CPU→GPU change: B1 BAOAB drops from bit-for-bit (0.0 on CPU) to 1.72e-15 on GPU**
  — GPU FMA / float-op reordering. Still passes the kT=0 < 1e-9 gate by 6 orders. This is the
  documented "GPU is machine-eps, not bit-for-bit" property; everything physical is unaffected.
- **B4 reverse-mode autodiff works on the CUDA backend** (grad vs analytic 2.6e-15) — Warp's
  differentiability prize is confirmed GPU-resident, not just a CPU spot-check.

## Consequence

The adoption make-or-break clears: the Warp DCM engine is **trustworthy on the A5000**. GPU
production and the differentiable loop are validated end-to-end at the backend that will actually
run them. G2 (substrate-force port) and the route-a / de-cohesion clean re-run are unblocked to run
*on GPU*.

## Method notes / repro

- Tree synced to gbook by **rsync** (`~/ffn_phase_c`), not Syncthing — the Syncthing copy of this
  worktree was stale (frozen ~2026-06-11). `PYTHONPATH=. python aleph/warp_port/parity_report.py --device cuda:0`.
- `parity_report.py` now threads `--device` (was CPU-hardcoded); cpu output → `warp_parity_results.json`,
  cuda → `warp_parity_results_cuda0.json` (both committed; the gate reads the committed lightweight JSON, C12 lesson).
- Not separately gated: the *hybrid loop / multi-cell* end-to-end on GPU is exercised by the committed
  A5000 bench (`dcm_*_bench_a5000.json`); this G1 verifies the per-kernel force/integrator parity that loop is built from.

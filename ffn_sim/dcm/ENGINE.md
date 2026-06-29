# The Warp DCM engine (Phase C)

The GPU-resident, **differentiable** DCM engine — the going-forward simulation engine
(PI 2026-06-21, `docs/v2_audit/PHASE_C_WARP_DCM_ENGINE_PLAN_2026-06-21.md`). Every force +
integrator piece is a `warp.kernel` parity-gated against a committed HOOMD reference; the
hybrid loop runs the per-step physics device-resident with low-cadence host remesh.

`engine.py` is the single import surface; this file says which entry point to use and how
the pieces fit.

## Architecture

```
per step  (GPU, no host sync):   force.zero
  → cohesion (node-node)  → per-cell turgor (volume-reduce → ΔP → face force)
  → node-face contact     → cortex edge springs   [→ lamellipodium tether]
  → [substrate: z-well + in-plane wetting]         → overdamped BAOAB step
every remesh_period (HOST):      remesh_pass (SWAP/SPLIT/COLLAPSE) on a fixed-size
                                 node-pool → GPU↔host resync
```

The node-pool keeps the GPU position array **fixed-size** (dormant slots; SPLIT activates,
COLLAPSE parks), so the device side is static-topology — Warp's no-dynamic-topology
constraint never bites. Remesh stays on the host (the correct architecture: SimuCell3D is
CPU/OpenMP, CellSim3D is GPU-but-fixed-topology — both avoid GPU dynamic remesh). Remesh is
<1 % of multi-cell runtime; a device-side edge-extent gate skips it entirely in equilibrium.

## Entry points (`from ffn_sim.dcm.engine import …`)

| function | use when | forces |
|---|---|---|
| `run_multicell` | **canonical** — a spheroid / ≥1 cells | per-cell turgor + cortex edges + node-node cohesion + node-face contact (+ hash-grid neighbour list, host remesh) |
| `run_hybrid` | single cell, fast inner loop | turgor + edges (+ remesh); `n_cells=1` special case |
| `run_crawl` | single cell spreading on a substrate | + lamellipodial traction-tether |
| `run_diff` | parameter estimation / inverse design | end-to-end **differentiable** loop; grads w.r.t. `dP0`, `k_edge` (autodiff vs FD 4e-12) |

Speed (A5000, `fixtures/dcm_*_bench_a5000.json`): single-cell 30× the existing path; multi-cell
hash-grid 10.9× (8 cells) → 36.8× (27 cells / 4374 nodes) — speedup grows with cell count.

## Force set (parity-gated building blocks, for custom loops)

Integrator: `run_baoab_warp` (overdamped L-M BAOAB) · `run_shake_warp` (M-SHAKE) ·
`run_fixman_warp` (metric pseudo-force). Compartment/cortex: `run_radial_shell_warp`
(turgor/membrane/nucleus) · `run_harmonic_bond_warp` · `run_harmonic_angle_warp` ·
`run_wca_pair_warp` (LJ-excluded-volume). DCM: `run_node_face_contact_warp` ·
`run_dcm_turgor_warp` · `run_dcm_cohesion_warp` · `run_dcm_substrate_well_warp` +
`run_dcm_substrate_wetting_warp` (the route-a spreading drivers — §G2) ·
`run_lamellipodium_tether_warp`.

DCM force runners that touch a reduction expose `reduce="host"` (force-law isolated, machine-eps
gate) and `reduce="warp"` (atomic scatter, reduction-order tol) — see the turgor/substrate modules.

## Parity contract & status

Nothing is "done" without a committed artifact + a passing parity test (`tests/dcm/`)
vs the **committed HOOMD reference fixture** (`fixtures/*_ref.npz`) — never a self-authored
oracle. `parity_report.py [--device cpu|cuda:0]` emits the committed verdict JSON
(`fixtures/warp_parity_results{,_cuda0}.json`) read by the results-integrity gate.

- **CPU backend** = the CI gate (`make warp-parity`): force-laws machine-eps; iterative /
  atomic-reduce at a documented tol.
- **GPU backend** = VERIFIED on the gbook A5000 (G1): same kernels, force-laws machine-eps,
  atomic-reduce < 1e-8, autodiff clean. `docs/v2_audit/PHASE_C_GPU_PARITY_2026-06-21.md`.

## Not yet unified (deferred to the first production run)

`run_multicell` does not yet wire the substrate (`run_dcm_substrate_*`) or lamellipodium into
its step loop — those kernels are parity-verified in isolation. The wiring + its smoke
validation belongs with the first run that needs it (route-a single-cell spreading / the
de-cohesion clean re-run), where the real force-set + parameters are known, rather than as a
speculative merge here. Until then, use `run_crawl` for lamellipodial spreading and compose the
substrate runners directly.

## Physics references (the DCM literature basis)

The DCM is a deformable-cell model — a cell is a pressurised elastic shell, not resolved
filaments (filament physics is the `ff/` layer's job: Cytosim, Nédélec & Foethke 2007). Its
physics is anchored to the published deformable-cell-model literature, used as reference only
(never imported as runtime), exactly as `ff/` references Cytosim:

| reference | role |
|---|---|
| **SimuCell3D** — Runser, Vetter & Iber 2024, *Nat. Comput. Sci.* 4:299–309 (PMC11052725) | primary physics basis: node-face contact, surface-tension faceting (γ̃ = γ/(K·l)), compressible-fluid cytoplasm `p = −K·ln(V/V₀)`, measured-value defaults. CPU/OpenMP. |
| **CellSim3D** — Madhikar, Åström, Westerholm & Karttunen 2018, *Comput. Phys. Commun.* 232:206–213 (PII S0010465518302091) | GPU (CUDA) deformable-cell colony growth + **division** (cf. our `dcm_cleave`); elastic-shell nodes + internal pressure + intercellular forces, fixed-topology-on-GPU (cf. our node-pool, §Architecture). Cross-check + a CUDA implementation reference for the Warp port. PDF in `references/`. |

Both are cell-shell models (no cytoskeletal filaments). KB registration candidates (not yet
Notion SE rows): `references/SE_REGISTRATION_CANDIDATES_2026-06-29.md`.

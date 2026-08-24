# I2b — filament-filament excluded volume (WCA steric) · `ac/solid-ev`

Session C (EXCLUDED-VOLUME track) of the new Active Cell engine. Net-new `ac/solid/` implementing the
CLAUDE.md HARD worked-example **"LJ repulsive excluded volume ON from Phase 1"** and closing one of
ENGINE_ARCHITECTURE_PLAN §3's three defining living-cortex gaps ("cortex filaments interpenetrate; the mesh
has no real steric volume under compression"). A soft repulsive **WCA** (purely-repulsive Lennard-Jones)
pair potential over filament nodes via a device `wp.HashGrid`, exposed as the
`StericForce.accumulate(state, out_force)` inner-loop primitive (§1.4). Field-independent.

## Status

- **Analytic gates: 24/24 GREEN** (`tests/ac/solid/`, pure NumPy — no Warp/CUDA). Full `ac/` suite 48/48,
  Warp-only contract PASS.
- **Warp kernel authored** (`steric_warp.py`, `StericForce` + hash-grid WCA kernel, ported from
  `dcm_neighbor_warp.cohesion_grid_kernel` own-row pattern) — the lead gates it natively on the gbook A5000
  against the `steric_reference` oracle (native-gate spec in `ac/solid/INTEGRATION.md`).
- **Deliverables:** green analytic module + Warp source + native-gate spec + `INTEGRATION.md` patch-notes +
  I0-B2b parameter gaps.

## Physics

WCA = the repulsive branch of LJ, shifted so energy AND force go to zero smoothly (C¹) at the cutoff
`r_c = 2^(1/6) σ`; no attractive tail (cohesion/crosslinks are separate force families):

    U(r) = 4ε[(σ/r)¹² − (σ/r)⁶] + ε   (r < r_c),   F(r) = −dU/dr = (24ε/r)[2(σ/r)¹² − (σ/r)⁶] ≥ 0.

- `σ_EV` (steric diameter) — **physical** (F-actin ~7 nm, MT ~25 nm), but coarse-graining-entangled at the
  0.5 µm node spacing → **GAP to PI** (node-node-coarse vs segment-segment-physical; `params_i0b2b.yaml`).
- `k_EV` (contact stiffness) — a **numerical** repulsion scale (Magic-Number Block: derived from operating
  load × penetration tolerance, grid-invariant, NOT tuned to a crowding outcome; the "held-apart" outcome
  saturates above the derived threshold). NOT a literature value.

## Gates (Session-C brief)

| Gate | Result |
|---|---|
| pair-potential FD-gradient (sign arbiter) | F = −dU/dr matches central FD to O(h²); F>0 repulsive inside r_c |
| zero-overlap == zero-force (OFF / regression) | U, F ≡ 0 for r ≥ r_c; EV-OFF & no-overlap both bit-identical |
| k_EV / CFL grid-invariance (Magic-Number Block) | cell-list == brute-force for any cell ≥ r_c (9 decades); k_EV enters kmax; outcome saturates → not tuned |
| compressed-pair resists collapse (finite steric volume) | r_eq > 0 under any finite load (even 1e8 pN) → core never fully overlaps |

**Native-gate scope (verify finding):** the pre-I3 native crowding gate certifies a **PASSIVE** (no myosin)
cell — structural/regression only. The load-bearing no-interpenetration-under-contractile-load verdict
**DEFERS** to a re-run after I3 / at I9 (`INTEGRATION.md`).

## Figures

Regenerate: `python -m ffn_sim.scripts.ac_solid_vis` → `outputs/ac/solid-ev/figs/*.png`. Each overlays the
closed-form oracle on the numeric result; units annotated; no axis truncation (§1.9).

- **i2b_wca_potential.png** — WCA U(r) (repulsive, C¹ to 0 at r_c) and F=−dU/dr with the central
  finite-difference gradient overlaid (the sign arbiter): analytic and FD coincide, F>0 everywhere inside r_c.
- **i2b_zero_force.png** — F(r) → 0 smoothly at r_c and identically 0 beyond (shaded): zero-overlap == zero
  force == OFF bit-identical; C¹ (not a hard-shell jump).
- **i2b_compressed_pair.png** — equilibrium separation r_eq vs compressive load (log): r_eq stays strictly
  >0 for any finite load (finite steric volume), r_c reference overlaid, y-axis full from 0.
- **i2b_grid_invariance.png** — (left) cell-list (HashGrid analogue) |F| vs brute-force |F| on the y=x
  identity across ~9 force decades for cell sizes {1, 1.5, 3}·r_c → the neighbour cell is a grid-invariant
  accelerator; (right) held-apart fraction r_eq/r_c vs k_EV saturating toward 1 above the derived tolerance
  → k_EV is derived, not tuned.

## Open items → PI

- `σ_EV` EV-discretization decision (node-node coarse vs segment-segment physical) — I0-B2b GAP.
- `k_EV` magnitude closes at the native passive-cell gate (measure F_op, set k_EV = F_op/(tol·r_c)) —
  report-not-tune.
- I7 LINC/plate soft-contact single-channel reconciliation (hard-truth #6) — handed off in `INTEGRATION.md`.

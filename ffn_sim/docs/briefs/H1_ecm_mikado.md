# Worker A — H.1 ECM Mikado in HOOMD

**Budget**: 2 weeks (full-fidelity assumptions absorbed)
**Branch**: `worker-a/h1-ecm` (cut from `ffn/foundation`)
**Reference**: Plan v2 §3 Unit H.1, [PHASE_0_3_DECISIONS.md](../PHASE_0_3_DECISIONS.md), [AFINES_ALGORITHM_NOTES.md §2 + §6](../AFINES_ALGORITHM_NOTES.md)
**Prereq**: none (parallel with H.2, H.4)

## Goal

Reproduce v1 Worker A's KU-1.30 ECM mechanical validation suite under the
v2 HOOMD-blue foundation, with full-fidelity discrete fiber Mikado
construction and no v1 numpy integrator dependence.

## Deliverables

| File | Content |
| --- | --- |
| `ffn_sim/ecm/mikado.py` | Mikado initial topology generator (REUSE `ffn_sim/validation/oracles/ecm/fiber_network.py` for the geometry call; do not import the v1 force kernel). |
| `ffn_sim/ecm/cross_links.py` | Segment-intersection cross-link seeding (REUSE `ffn_sim/validation/oracles/ecm/cross_links.py` geometry only). |
| `ffn_sim/ecm/shear_protocol.py` | Strain ramp / hold schedule against HOOMD `BoxResize` (fresh write — v1 `shear_protocol.py` is archived). |
| `ffn_sim/integrator/baoab.py` | D3: Leimkuhler-Matthews BAOAB-limit custom HOOMD `Updater`. Stores per-particle `prv_rnds`. (See PHASE_0_3_DECISIONS D3.) Shared across H.1/H.2/H.3/H.4/H.5. |
| `ffn_sim/configs/phase1_h1.yaml` | KU-1.x parameter values ported from `ffn_sim/validation/oracles/configs/phase1_unit1.yaml` with D4 overrides (ℓ₀=0.5 μm, dynamics: hoomd_brownian+baoab). |
| `ffn_sim/tests/test_h1_mikado.py` | Unit tests: bond/angle topology size, energy oracle vs `ffn_sim/validation/oracles/ecm/fiber_mechanics.py` (≤ 1e-6 relative). |
| `ffn_sim/tests/validation/test_ku130.py` | KU-1.30 validation against v1 Worker A's frozen results. |
| `ffn_sim/outputs/h1/REPORT.md` | Wall-time vs numpy benchmark + KU-1.30 PASS evidence. |

## Implementation spec

### Topology

- N_f filaments × 21 beads each, ℓ₀ = 0.5 μm bond rest length, L_f = 10 μm.
- Particle type `actin_ecm`, R_bead = 50 nm (KU-1.2 collagen fibril radius).
- N_f resolved from `phase1_unit1.yaml` derived params: Mikado `ρ_L = N·L_f/L_box²`, target ⟨z⟩ from KU-1.3.
- Bonds: `md.bond.Harmonic` type `ecm-bond`, `k = μ` (KU-1.2: 8.6e-9 N at L_f, scale to per-bond), `r0 = 0.5e-6` (m).
- Angles: `md.angle.Harmonic` type `ecm-angle`, `k = κ_B / ℓ₀` (KU-1.1: κ=7e-26 N·m², ℓ₀=0.5e-6), `t0 = π`.
- Cross-links: harmonic bonds added at segment intersections, type `xl`, `k = 1e-3 N/m` (KU-1.28), `r0` = distance at intersection (≈ 0).

### Integration

- D3 custom L-M BAOAB-limit Updater (no vanilla `md.methods.Brownian`).
- Δt from `derived_params.tau_min · cfl_safety_factor` (KU-1.26 conventional 0.1 · τ_min).
- Box: 2D-thin-slab convention (HOOMD is 3D — `box.Lz = 200 nm` with z-pinning external field, or accept periodic z and report 2D projection).
- BC: periodic in x, y; Lees-Edwards shear via `hoomd.update.BoxResize` for KU-1.30 #2 strain stiffening.

### Excluded volume (D7)

- `md.pair.LJ` WCA repulsive-only between actin_ecm beads.
- σ = 2·R_bead = 100 nm, ε = 0.5 kT, cutoff = 2^(1/6)·σ.

## Validation acceptance

| Gate | Criterion | Source |
| --- | --- | --- |
| Topology smoke | Generated `(N_f, 21, 3)` bead array; bond count = N_f·20; angle count = N_f·19; cross-link count consistent with Mikado `n_int = 2·ρ_L·L_f/π` | v1 `ffn_sim/validation/oracles/ecm/fiber_network.py` oracle |
| Energy oracle | HOOMD `ThermodynamicQuantities.potential_energy` vs `ffn_sim.validation.oracles.ecm.fiber_mechanics.compute_energy` at the same configuration: relative error ≤ 1e-6 | KU-1.24 |
| KU-1.30 #1 | G_0 ∈ [15, 200] Pa (v1: 32 Pa) | v1 commit `11eaf13` |
| KU-1.30 #2 | Strain stiffening exponent ∈ [−2.5, −1.5] | v1 commit `11eaf13` |
| KU-1.30 #3 | Point-dipole 1/r² stress decay within KB-gap band | v1 commit `d92ac20` |
| **v2 wall-time** | M1 Max CPU wall-time per simulated second ≤ 5× the v1 numpy time (10× speedup is the long-term GPU target; CPU should not regress > 5×) | Phase 0.2 bench reference |

## Open implementation questions for Worker A

- Choose 2D-thin-slab z-pinning vs 3D-periodic-with-projection — discuss
  with PI if either is awkward.
- L-M Updater plumbing: HOOMD 7 `hoomd.custom.Action` invoked at every
  step requires careful interaction with `md.Integrator.forces`. Worker A
  should benchmark the L-M Updater on the 100-bead polymer sanity test
  before scaling to ECM.

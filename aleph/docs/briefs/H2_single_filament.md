# H.2 — Single F-actin filament persistence length validation

**Budget**: 1 week
**Branch**: `phase1/h2-single-filament` (cut from `ffn/foundation`)
**Reference**: Plan v2 §3 Unit H.2, [PHASE_0_3_DECISIONS.md](../PHASE_0_3_DECISIONS.md), [AFINES_ALGORITHM_NOTES.md §2 + §6](../AFINES_ALGORITHM_NOTES.md)
**Owner**: Main Session.
**Prereq**: H.1 BAOAB integrator (`aleph/integrator/baoab.py`) merged.

## Goal

Validate the v2 HOOMD bond+angle+L-M integrator stack against actin's
canonical persistence length L_p = 17 μm (KU-1.1). This is the
single-particle-physics sanity gate that every cortex / lamellipodium /
ECM filament in Phase 1+ depends on.

## Deliverables

| File | Content |
| --- | --- |
| `aleph/archive/filament_math.py` | Tangent correlation `⟨t̂(s)·t̂(0)⟩`, persistence length fit, Boltzmann angle distribution, equipartition check — shared utilities. |
| `aleph/scripts/h2_single_filament.py` | Run script: build one 10 μm filament, run L-M Brownian under Langevin thermostat T=310 K for 60 s simulated. |
| `aleph/tests/test_persistence_length.py` | Pytest gate: L_p measured ∈ [15.3, 18.7] μm, KS test p > 0.05 on angle distribution, equipartition |⟨E_bend⟩ − k_BT/2| / (k_BT/2) ≤ 0.05 per bond. |
| `aleph/outputs/h2/REPORT.md` | Measured L_p, tangent correlation fit, Boltzmann histogram, equipartition table. |

## Implementation spec

- **Single filament** (no crosslinks, no boundary interactions):
  - N = 21 beads, ℓ₀ = 0.5 μm bond rest length, L_f = 10 μm.
  - Particle type `actin`, R = 30 nm (actin bundle cross-section).
  - Bonds: `md.bond.Harmonic`, `k` from KU-1.2 spring stretching modulus, `r0 = 0.5e-6 m`.
  - Angles: `md.angle.Harmonic`, `k = κ_B / ℓ₀ = 0.068 pN·μm / 0.5 μm = 0.136 pN`, `t0 = π`.
- **Solvent**:
  - Free space (no periodic box constraint within the filament's worm length); use a 50 × 50 × 0.2 μm slab box with large enough x to fit the filament + thermal excursion.
- **Integrator (D3)**:
  - L-M BAOAB-limit (`aleph/integrator/baoab.py` — shared with H.1).
  - Δt = 2×10⁻⁵ s (AFINES default).
  - Langevin temperature T = 310 K (kT = 4.28×10⁻²¹ J), γ_b = 6π·η·R with η = 0.001 Pa·s.
- **Run**:
  - Equilibrate 10 s (5×10⁵ steps).
  - Sample 60 s (3×10⁶ steps), write a frame every 1000 steps (3000 frames).
- **Measurement**:
  - For each sampled frame, compute bond unit vectors t̂_i, then `C(s) = ⟨t̂_i · t̂_{i+s}⟩` averaged over starting bond i.
  - Fit `C(s) = exp(−s·ℓ₀ / L_p)` for s ∈ [1, N/2] bonds.
  - Histogram θ_i (deviation from straight) per timestep; KS-test against `p(θ) ∝ exp(−κ_B θ² / 2 ℓ₀ k_BT)`.
  - Per-bond bending energy `E_i = (κ_B / 2 ℓ₀) θ_i²` averaged over frames; check `⟨E⟩ = kT/2`.

## Validation acceptance

| Gate | Criterion |
| --- | --- |
| L_p measurement | `15.3 μm ≤ L_p ≤ 18.7 μm` (target 17 ± 10 %) |
| Angle distribution | KS test p > 0.05 vs Boltzmann form |
| Equipartition | `\|⟨E_bend⟩ − kT/2\| / (kT/2) ≤ 0.05` per bond |
| L-M correctness | Same L_p (within 1 σ) when re-run with vanilla `md.methods.Brownian` at halved Δt — confirms L-M ≠ artefact |

## Open implementation questions

- HOOMD 2D-emulation: confirm box geometry that lets the filament remain
  effectively 2D (Plan v2 is 2D in Phase 1) without breaking thermal
  statistics. The slab convention (`Lz = 200 nm`, `periodic_z = True`)
  is probably fine for free filaments but verify.
- L-M Updater interaction with `md.compute.ThermodynamicQuantities` —
  if the L-M Updater bypasses HOOMD's velocity table, thermo readings
  may need a custom computer.

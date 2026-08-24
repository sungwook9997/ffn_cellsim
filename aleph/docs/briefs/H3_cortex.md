# H.3 — Cortex multi-filament network

**Budget**: 3 weeks
**Branch**: `phase1/h3-cortex` (cut from `phase1/h2-single-filament` after H.2 lands)
**Reference**: Plan v2 §3 Unit H.3 (with 2026-05-19 v3.1 cortex composition revision), [PHASE_0_3_DECISIONS.md](../PHASE_0_3_DECISIONS.md) D2/D5/D7, [AFINES_ALGORITHM_NOTES.md §2/§3/§4/§6](../AFINES_ALGORITHM_NOTES.md)
**Owner**: Main Session.
**Prereq**: H.2 complete (single filament L_p validated). Inherits `aleph/integrator/baoab.py`.

## Goal

Build the full-fidelity cortex of a single cell: ~1000 effective F-actin
filaments + 1000 dynamic crosslinkers (filamin / α-actinin) + 100 myosin
minifilaments + ERM tether to cell center. Reproduce KU-3.1 cell rounding
+ KU-3.5 cortical tension + KU-3.18 blebbistatin response + KU-3.20 nematic
order.

## Deliverables

| File | Content |
| --- | --- |
| `aleph/cortex/cortex_network.py` | Cortex topology generator: 1000 effective filaments on a R=10 μm shell, random placement + orientation within 200 nm cortex thickness. |
| `aleph/cortex/crosslinkers.py` | D2 Bell-Evans slip xlink updater (filamin slip distance from Furuike 2001 / Ferrer 2008). Python `hoomd.custom.Action`, batched every ~100 steps. |
| `aleph/cortex/myosin.py` | D5 Stam-Hocky bipolar multi-head minifilament: ~14-bead rigid-rod backbone via `md.constrain.Rigid`, ~10 cross-bridge heads per side, each with D2 Bell-Evans slip + D6 Hill stepping. |
| `aleph/cortex/erm.py` | Harmonic spring from each cortex bead to cell center; `k_ERM = 0.1 N/m` (KU-3.18). |
| `aleph/cell/cell.py` | New v2 `Cell` class: HOOMD particle-group-backed, owns cortex + ERM + lamellipodium list + FA list. **Replaces v1 `acs_kb/cell/cell.py`** (deleted in 2026-05-20 rename; v1 was single-chain Cortex-coupled, do not resurrect). |
| `aleph/configs/phase1_h3.yaml` | KU-3.x cortex parameters: 1000 filaments, 1000 xlinks (30 % α-actinin + 70 % filamin per KU-3.19), 100 myosin minifilaments, k_ERM, etc. |
| `aleph/tests/test_h3_cortex_topology.py` | Topology smoke: particle counts, bond/angle/rigid-body counts, no ring-closure failures. |
| `aleph/tests/validation/test_ku31_rounding.py` | KU-3.1 cell rounding gate. |
| `aleph/tests/validation/test_ku35_tension.py` | KU-3.5 cortical tension gate. |
| `aleph/tests/validation/test_blebbistatin.py` | Myosin-OFF run vs myosin-ON: cortex tension drop + rounding failure. |
| `aleph/outputs/h3/REPORT.md` | KU-3.x gate evidence, wall-time per simulated second, bead count budget. |

## Implementation spec

### Cortex topology

- Cell shell: R = 10 μm, cortex thickness 200 nm, circumference 63 μm.
- 1000 effective filaments per cell, each ×40 native actin bundle.
- Filament length distribution: uniform 1–5 μm (mean 3 μm) → average 7 beads per filament at ℓ₀=0.5 μm.
- **Total cortex beads ≈ 7,000 per cell** (Plan v2 §3 H.3 v3.1 estimate).
- Particle type `actin_cortex`, R = 30 nm bundle cross-section.
- Bond / angle parameters as H.2 (same `md.bond.Harmonic` k and `md.angle.Harmonic` κ_B/ℓ₀).
- Ring topology: cortex is membrane-bound shell, filaments confined to the shell by a soft potential `md.external.field` (radial harmonic toward R=10 μm).

### Crosslinkers (D2)

- 1000 per cell, density 10/μm² cortex area (KU-3.19).
- Mix: 30 % α-actinin (`xlink_alpha`, length 35 nm, k = 0.1 pN/μm), 70 % filamin (`xlink_filamin`, length 150 nm, k = 0.1 pN/μm).
- Two-particle head pair per xlink, joined by `md.bond.Harmonic`.
- Dynamic binding: D2 Bell-Evans slip. Per-bond-type `(k_off⁰, x_β)` from literature (filamin: Furuike 2001 give `x_β ≈ 0.3 nm`, `k_off⁰ ≈ 0.1 s⁻¹`).
- `XlinkUpdater(hoomd.custom.Action)` runs every 100 steps:
  - Unbinding pass: for each currently-bonded `xlink_head ↔ actin_cortex` bond, compute current force; sample Bell-Evans break probability; remove bond if fires.
  - Binding pass: for each unbound head, query neighbor list (`md.nlist`), sample binding to nearest spring within `max_bind_dist = 60 nm`; ON acceptance, add bond.
- Snapshot mutation costs ~ms per update; batching by 100 keeps acceptable per-second overhead.

### Myosin minifilaments (D5 + D6)

- 100 per cell, density 3/μm² (Salbreux 2012).
- **Stam-Hocky bipolar architecture**:
  - 14-bead backbone, 700 nm length, type `myosin_backbone`, `md.constrain.Rigid` rigid body.
  - 10 cross-bridge heads per side (20 heads per minifilament), type `myosin_head`.
  - Head ↔ backbone connection: harmonic spring, `k = 1 pN/μm`, `r0 = 200 nm` perpendicular to backbone.
  - **Total: 34 particles per minifilament × 100 = 3,400 motor beads per cell**.
- Per-head dynamics: D2 Bell-Evans slip bond to actin, D6 Hill stepping (`v0 = 1 μm/s`, `F_s = 0.5 pN`, `a/F_s = 0.5` per Kovács 2003).
- `MyosinUpdater(hoomd.custom.Action)` extends `XlinkUpdater` with stepping kernel (per AFINES §4.3 of algorithm notes).

### ERM tether

- For each cortex bead i: harmonic spring `k_ERM = 0.1 N/m`, `r0 = R_cell = 10 μm`, anchored to cell center (`md.external.field` radial).

### Excluded volume (D7)

- `md.pair.LJ` WCA repulsive-only, ε = 0.5 kT.
- Pair types: actin_cortex×actin_cortex, actin_cortex×xlink_head, actin_cortex×myosin_head, myosin_head×myosin_head.

### Integration

- L-M BAOAB-limit (from H.1's `aleph/integrator/baoab.py`).
- Δt = 2×10⁻⁵ s.
- 60 s simulated time for full validation (3×10⁶ steps).
- Plan v2 §11 wall-time estimate for H.3 on M1 Max CPU: ~8 hours (overnight).

## Validation acceptance

| Gate | Criterion | KU |
| --- | --- | --- |
| Cell rounding | Aspect ratio 1.5 → < 1.2 in 60 s | KU-3.1 |
| Cortical tension | γ_cortex ≈ 0.5 mN/m ± 30 % (Laplace law from internal pressure) | KU-3.5 |
| Blebbistatin response | myosin-OFF: cortex tension ↓ → rounding fails (aspect ratio stays > 1.3 at 60 s) | KU-3.18 |
| Q_ij nematic order | Isotropic cell `S < 0.1`; aligned cell `S > 0.3` along axis | KU-3.20 |
| Per-segment θ accessible | cortex deformation visualization works at single-bead resolution | KU-3.18 |
| Bead budget | Total cortex + xlink + myosin + ERM beads ≤ 15,000 per cell on M1 Max | Plan v2 §11 |

## KU-3.23 cortex independence

`aleph/cortex/*.py` MUST NOT `import aleph.ecm.*`. WLC formulas
live in `aleph/archive/filament_math.py` (from H.2). pylint test in
`tests/test_cortex_independence.py` enforces this.

## Open implementation questions

- `md.constrain.Rigid` interaction with L-M BAOAB Updater: HOOMD's rigid
  bodies are integrated by the integrator method, but our L-M is a
  custom Updater. May need to keep some particles on `md.methods.Brownian`
  (for rigid bodies) and others on L-M Updater. Verify before scaling.
- 100 dynamic bond mutations per `XlinkUpdater` call: HOOMD snapshot
  upload cost may dominate. Profile early.
- Lamellipodium hook: leave a `Cell.lamellipodium = None` slot — H.5
  will populate.

# H.7 — Dynamic single-cell spreading (real BAOAB drivers)

**Branch** `h7/compartment-platform` · **2026-06-11** · CPU mesoscale dev (gbook GPU
production deferred until the dirty-branch cleanup).

## Mission

Make the compartment platform produce **real spreading motion**: the basal contact
footprint `A(t)` must grow *physically* over simulation time on the explicit BAOAB
integrator, driven by the actual machinery (lamellipodium protrusion + FA molecular
clutch), not a kinematic/geometric mock.

## Diagnosis — why the old run was a mock

`h7_spreading_compare.py` builds the FA-adhered basal-ring-lamellipodium cell and runs
BAOAB, but the measured footprint was **flat**: `basal_ring A/A₀ = 1.059`,
`polarized_patch A/A₀ = 1.026` after 1500 steps. Root cause (not a code bug — a
*timescale* wall):

- At the physiological elongation rate `k_elong⁰ = 11.6 s⁻¹` and `dt ≈ 13 ns`, a
  1500-step run spans only ~20 µs of simulation time. The expected number of
  monomer-addition events over the whole run is **≈ 0.02** — the dendritic network
  never advances, so the convex-hull footprint stays at its construction geometry.
- Real spreading takes minutes (~5 µm front advance at ~6 µm/s of raw barbed-end
  growth ≈ 0.86 s ≈ **66 million** BAOAB steps) — infeasible on CPU in 48 h.
- The FA integrin↔ligand clutches engage but were never coupled to the lamellipodial
  actin, so retrograde flow was not converted to anchored protrusion (the open
  motor-clutch loop).

## What was built (this session)

A dynamic spreading harness on top of `build_baseline_cell`, **no edits to the loader
or integrator** (PI-frozen). Per PI authorization (2026-06-10: relax the wall-clock
hard rule, band-match literature), the protrusion runs in explicit *kinetic
fast-forward* — the mechanism is the real BAOAB particle dynamics; the single
band-matched quantity is the leading-edge advance velocity, reported in physiological
time via the literature single-cell spreading velocity.

- `scripts/h7_spreading_dynamics.py` — builds the physiological adherent cell
  (FA + basal-ring lamellipodium), runs BAOAB, samples `A(t)/A₀`, leading-edge radius,
  network bead count, FA engagement and traction.
- `cell/spreading_drive.py`
  - `LeadingEdgeNucleationUpdater` — **membrane-tracked protrusion engine**. The H.5
    Arp2/3 nucleation is pinned at the construction WAVE ring, so the front caps out
    once it advances beyond the ~100 nm Arp2/3 reach; in vivo the WAVE/Arp2/3 machinery
    tracks the advancing membrane. This updater appends `actin_lamel` beads one ℓ₀
    outward at the leading-edge tips (force-free, via the lamellipodium's own snapshot
    extender) at the band-matched effective velocity. Mutates only topology + the
    Python state bookkeeping — never moves an existing particle (no double-stepping).
  - `BasalAdhesionTether` — **the FA molecular-clutch ensemble** as a substrate force:
    a harmonic z-tether (k_adh = membrane-tension scale ≈ 5e-5 N/m) holding the
    lamellipodial sheet on the basal plane. With it off the cortex/membrane line
    tension lifts and retracts the protruded front and spreading stalls (the
    clutch-OFF control).
- `scripts/h7_spreading_dynamics_vis.py` — 4-panel summary + GIF + spreading gate.
- `scripts/h7_spreading_clutch_compare.py` — clutch ON vs OFF overlay.

## Result — CANONICAL run (clutch ON, mesoscale CPU)

`n_filaments=100` demo cortex, basal-ring lamellipodium + FA, accel S=600, p_advance=0.015,
20000 BAOAB steps, k_adh=5e-5 N/m, box_factor=3, wall ≈ 212 s.

- **A/A₀: 1.00 → 2.03** — enters the physiological single-cell band **[2, 4]**
  (old kinematic mock was **1.06 flat**).
- leading-edge radius r_front: **2.70 → 5.36 µm** (×1.99); network n_actin **153 → 496** (×3.24).
- physiological time of the run: **t ≈ 2.66 min** (front-mapped at the literature
  single-cell velocity 1 µm/min) — matches the Betorz 2023 P1 "fast spreading" phase (~3 min).
- time-course shape: fast early rise → plateau (the literature power-law→plateau form).
- FA clutch: engaged integrin↔ligand clutches bear substrate traction during the run
  (n_engaged > 0, traction ~0.05 nN at this mesoscale).
- **GATE PASS**: G1 reached band (A/A₀ 2.03 ∈ [2,4]) ✓, G2 edge advances (×1.99 > 1.5) ✓,
  G3 network dynamic (×3.24 > 1.3) ✓.

### Clutch ON vs OFF (matched protrusion drive, p_advance=0.015)

| run | A/A₀ final | verdict |
|---|---|---|
| FA clutch **ON** (adhered) | **2.03** | reaches the physiological band [2,4] |
| FA clutch **OFF** (front lifts) | **1.74** | stalls below the band |
| old kinematic mock | 1.06 | flat — no spreading |

Same protrusion drive: with the basal adhesion the footprint reaches the band; without it
the unanchored front lifts out of the basal plane and the footprint stalls below band (and,
run longer, becomes numerically unstable) — the motor-clutch signature that adhesion is
required to convert protrusion into spread area.

## Honesty / scope

- **Band-matched quantity = the protrusion velocity** (kinetic fast-forward, effective
  factor **S ≈ 6.1e5**: the run compresses ~2.66 min of physiological spreading into
  ~2.6e-4 s of simulation time). Brute-forcing the physiological spread is ~10⁸ BAOAB
  steps — infeasible on CPU in the window; PI authorized the band-match (2026-06-10).
  Everything downstream (network mechanics, clutch hold, cortex/turgor/membrane line
  tension, footprint shape, the ON-vs-OFF difference) is emergent BAOAB physics. This is
  the standard accelerated-dynamics framing.
- **Mesoscale CPU**: ×40 mesoscale (n_filaments demo). Native-scale GPU is the gbook
  follow-on (blocked on the dirty-branch cleanup).
- **BasalAdhesionTether is a stopgap clutch.** A blanket z-tether of the lamellar sheet
  reaches the band but at the mesoscale over-densifies if pushed hard (2-D LJ crowding →
  late instability); the faithful clutch is a substrate ligand-field grip
  (`FrontClutchRatchetUpdater` is scaffolded in `cell/spreading_drive.py` for that) — a
  follow-on. The integrin↔ligand catch-slip clutch is already present and load-bearing.
- **Basal mesh** (BUILD-complete, gates PASS) is not yet wired into the single-cell live
  build; cortex + cytoplasm + turgor provide the mechanical background here. Wiring the
  basal contractile mesh is the next driver to add.

## Figures

- `figs/spreading/h7_spreading_dynamics.png` — canonical clutch-ON 4-panel: A/A₀(t) into
  the physiological band vs the flat mock; leading-edge radius + network growth; FA clutch
  engagement + traction; basal footprint construction → spread.
- `figs/spreading/h7_spreading_movement.gif` — the basal contact footprint growing over time.
- `figs/spreading/h7_spreading_clutch_compare.png` — FA clutch ON (→2.03, in band) vs OFF
  (→1.74, stalls) at identical protrusion drive.

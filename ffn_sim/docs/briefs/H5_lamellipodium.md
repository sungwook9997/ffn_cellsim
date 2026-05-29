# H.5 — Lamellipodium (Bieling/Funk greenfield, 4-week budget)

> **Owner**: Lead session.
> **Prereq**: H.3 ✅ DONE (cortex topology + crosslinkers + myosin + ERM +
> `Cell.build` composition).
> **Critical-path role**: Adds the protruding lamellipodium subsystem to
> the single-cell stack. H.5 ✅ DONE is the prerequisite for H.7
> single-cell integration.
>
> **Status (2026-05-29)**: 🟨 implementation + wiring + KU-5.1 trio
> scaffold PR-ready. KU-5.x production runs blocked on PI sign-off +
> Mac smoke. See [outputs/h5/REPORT.md](../../outputs/h5/REPORT.md) for
> the closeout chronology.
>
> **NB**: This brief was reconstructed 2026-05-29 from the existing
> implementation (`cell/lamellipodium.py` module docstring + Sanity Gate
> §1-6 + `tests/validation/test_ku5x_lamellipodium.py` KU table) because
> the original `H5_lamellipodium.md` referenced in the H.5 단계 1 Notion
> closeout was never committed to disk. Source-of-truth contract still
> lives in the module + tests; this is the human-readable map.

## Goal

Implement a 3D HOOMD-blue lamellipodium subsystem that produces the
correct **dendritic density** (≈ 100 barbed ends per μm² of WAVE plane),
**Bieling 2016 force-velocity** behaviour for barbed-end elongation, and
**Funk 2022 abortive-branching** (rate drops > 50 % above 500 Pa per-WAVE
force). All three are emergent mechanistic outputs, not lumped
phenomenology (per CLAUDE.md hard rule).

## Topology

- **WAVE plane** at `y = +Y_max` (top of simulation box; opposite the
  cortex shell which sits at the box origin). `n_WAVE = σ_NPF · WAVE_area`
  particles pinned to the plane via radial harmonic force
  (`WaveMembranePin`, `k_wave_pin` from YAML).
- One **mother actin seed** per WAVE at construction (H.5 brief open
  question (a) "seeded near WAVE plane"). Mother grows from WAVE into
  cytosol (barbed end pointed away from membrane, toward `−y`).
- Branching produces daughter actins at the canonical Arp2/3 angle
  `t0 = 72° = 5π/12` rad via `lamel_branch_angle` harmonic.

## Particle types

| Type | Role |
|---|---|
| `wave_particle` | membrane WAVE/NPF reservoir; pinned at plane |
| `actin_lamel` | lamellipodial actin bead; per-bead state tracks BARBED-END (eligible for elongation/branching/capping) vs INTERIOR (frozen) via Python-level `LamellipodiumState`, not HOOMD topology |

## Bond / angle types

| Type | Role |
|---|---|
| `lamel_actin_bond` | actin-actin harmonic backbone bond |
| `lamel_branch_bond` | Arp2/3 mother-to-daughter branch bond |
| `lamel_wave_anchor` | WAVE-to-mother-bead initial anchor (optional; primary pin is the membrane plane force compute) |
| `lamel_branch_angle` | Arp2/3 branch angle harmonic, `t0 = 72°`, `k_angle = 100 pN·μm/rad²` (literature TBD, brief default) |

## Updaters (3 D1 D2-batched)

All three Updaters batch their stochastic kinetics every `batch_steps`
HOOMD steps (rather than per step) to amortise the Python Action cost
(same pattern as `cortex/crosslinkers.py:XlinkBondUpdater`).

### `BarbedEndElongationUpdater` (Bieling 2016 slip Bell-Evans)

```
k_elong(F) = k_elong⁰ · exp(−F · δ_elong / kT)
k_elong⁰ = 11.6 s⁻¹    δ_elong = 2.7 nm
```

Per batch tick, for every eligible (non-capped) barbed-end tag:
1. Compute per-barbed-end load `F = F_partitioned` (network load divided
   by `N_free_barbed_ends` per the D1 spec).
2. Sample `p_elong = 1 − exp(−k_elong(F) · Δt_batch)`.
3. On fire: append one actin bead at rest length `ℓ_0` along the
   filament's local tangent direction.

### `ArpBranchingUpdater` (Bieling 2016 + Funk 2021 mechanistic, γ-Phase 1 PI-ratified 2026-05-29)

Revised by H.5/γ ratification 2026-05-29 (design:
`H5_GAMMA_BRANCHING_DESIGN.md`, commit `62d6fbc`). The pre-γ rate
equation was independent per-WAVE; per Funk 2021 ("A barbed end
interference mechanism reveals how capping protein promotes
nucleation"), branching rate is gated by the CP-mediated release of
NPFs sequestered to free barbed ends:

```
free_npf_fraction = n_capped / max(n_barbed + n_capped, 1)
k_b_eff(F) = k_b⁰ · (1 − 0.2 · F / F_stall_branch) · free_npf_fraction
k_b⁰ = 0.037 s⁻¹  (per WAVE molecule, Bieling 2016)
```

At `n_capped = 0` (initial): `k_b_eff = 0` — no branching (Funk's
NPF-sequestered regime). At `n_capped → n_total`: branching at
`k_b⁰ · force_factor` (asymptotic free-NPF regime).

Per batch tick, for every WAVE molecule:

1. Compute `free_npf_fraction` (global state).
2. Sample `p_branch = 1 − exp(−k_b_eff · Δt_batch)`.
3. On fire: enumerate ALL `actin_lamel` beads within `r_branch_eff`
   of the WAVE (Funk 2021 + Bieling 2023: Arp2/3 binds along F-actin
   side, not only at the barbed tip).
4. Pick mother site weighted by `1/r²` (diffusion-limited kinetics,
   PI-ratified 2026-05-29).
5. Insert a daughter actin bead at `r_mother + ℓ₀ · t_daughter` where
   `t_daughter` is the mother's local tangent rotated 72° around a
   random perpendicular axis (Arp2/3 crystal).
6. Add `lamel_branch_bond` (mother→daughter) + (if mother has a
   parent) `lamel_branch_angle` harmonic.
7. Daughter starts as a new barbed end; mother's barbed/capped/interior
   state is unchanged.

**Funk 2022 abortive-failure emerges naturally**: above the per-WAVE
force threshold `F > 500 Pa · WAVE_area`, the `(1 − 0.2 F/F_stall)`
factor becomes negative and is clamped to zero. No separate scalar
override — the abortive regime is the same mechanism with the rate
running to zero.

**Pre-γ (KU-5.1 v1 baseline)**: branching searched only `barbed_end_tags`
within `r_branch = 30 nm`. After the first elongation step
(`Δℓ = ℓ₀ = 500 nm ≫ 30 nm`), every mother barbed-end left the search
sphere; system-wide branchings ≈ 0.32 events permanently. KU-5.1 v1
production sweep (1.85M step, `outputs/h5/production/ku51_v1/seed1.log`,
killed after analytical prediction 100%-confirmed) measured density
0.05 /μm² (target ≈100 /μm²) — falsified the pre-γ implementation.

### `CappingUpdater` (Funk 2022 slip Bell-Evans)

```
k_cap(F) = k_cap⁰ · exp(−F · δ_cap · sin θ / kT)
k_cap⁰ = 3 s⁻¹    δ_cap = 0.3 pN
```

Per batch tick, for every eligible barbed-end:
1. Compute per-barbed-end load `F`; mother-WAVE angle `θ`.
2. Sample firing probability.
3. On fire: mark barbed end CAPPED — tag moves to
   `LamellipodiumState.capped_tags`; the actin chain is preserved (no
   topology mutation, just per-bead state flip).

## Runtime state

`LamellipodiumState` (`cell/lamellipodium.py:407`):

| Field | Type | Role |
|---|---|---|
| `barbed_end_tags` | `list[int]` | currently-eligible barbed-end actin bead tags |
| `capped_tags` | `set[int]` | capped barbed-end actin bead tags (frozen) |
| `parent_of` | `dict[int, int]` | per-daughter actin bead, the parent (mother barbed-end) tag |
| `tangent_of` | `dict[int, np.ndarray]` | per-barbed-end unit tangent vector (direction of next elongation) |
| `actin_next_tag` | `int` | next available global tag for newly-elongated actin beads |

## Validation acceptance (KU-5.x gates)

Per H.5 brief original §Validation + `tests/validation/test_ku5x_lamellipodium.py`
header table:

| Gate | Criterion | Driver | Status |
|---|---|---|---|
| KU-5.1 | Dendritic density ≈ 100 barbed ends per μm² of WAVE plane at steady state | `scripts/h5_ku51_density.py` (scaffold `927d798`) + `scripts/h5_ku51_density_vis.py` (`7a424e5`) | scaffold PR-ready; PI band calibration pending |
| KU-5.2 | Bieling force-velocity `v_network(F)` matches Bieling 2016 Fig 2 within ±30 % | TBD | scaffold not started |
| KU-5.3 | Funk abortive branching rate drops > 50 % above 500 Pa per-WAVE force | TBD | scaffold not started |

All three gates require multi-hour BAOAB simulation to reach
branching/capping/elongation equilibrium (same wall-time class as
KU-3.x and L_p FULL). Tests skip under `H5_KU5_PRODUCTION=1` opt-in.

## Sanity Gate (per CLAUDE.md hard rule)

Authoritative source: `cell/lamellipodium.py` module docstring §Sanity
Gate §1-6. Summary:

1. **Dimensional analysis** — all rates [1/s], forces [N], `δ_elong`,
   `δ_cap` [m], `kT` [J], Bell-Evans exponent dimensionless. Branch
   `(1 − 0.2 F/F_stall)` dimensionless ratio. ✓
2. **Boundary cases** — F=0 recovers `k_elong⁰`/`k_b⁰`/`k_cap⁰`;
   `F > F_stall / 0.2` clamps branching to 0 (Funk abortive emergent);
   all barbed ends capped → elongation+branching idle, capping idle.
3. **Conservation** — actin chain preserved on capping (state mutation,
   not topology mutation). Daughter insertion conserves mass (new bead
   appended to actin tag space).
4. **Numerical sanity** — WAVE plane pinning `k_wave_pin` must satisfy
   ERM-style CFL `dt ≤ 0.1 · γ_b / k_wave_pin`; D2 batch tick CFL
   `dt_batch · k_max < 1`; both gates raise at attach time.
5. **Sign sense** — elongation along tangent (barbed end advances away
   from WAVE), branching inserts daughter at `+ℓ_0` along tangent at
   `72°` off, capping monotonically reduces eligible-barbed-end count.
6. **Measurement protocol** — dendritic density per-WAVE barbed-end
   count averaged over the WAVE plane area; abortive branching count
   per WAVE per second; force-velocity is per-barbed-end velocity vs
   per-barbed-end load.

## Composition into Cell

`Cell.build(with_lamellipodium=True, p_lamellipodium=...)` wiring
landed 2026-05-29 (commit `3306eb1`). The wiring:

1. Pre-extends the cortex snapshot with WAVE + mother-actin beads
   BEFORE `create_state_from_snapshot` (HOOMD can't add particle /
   bond / angle types post-init via `set_snapshot`).
2. Registers `lamel_*` bond / angle params on the simulation's force
   computes.
3. Wires LJ pairs (intra-WCA, inter-cortex / inter-myosin / inter-xlink
   disabled).
4. Attaches `WaveMembranePin` membrane plane force.
5. Extends BAOAB `gamma_map` with `wave_particle` + `actin_lamel`
   Stokes drag.
6. Appends the 3 D2-batched Updaters
   (`BarbedEndElongationUpdater`, `ArpBranchingUpdater`,
   `CappingUpdater`).

Public helpers in `cell/lamellipodium.py`:

- `extend_cortex_snapshot_with_lamellipodium(...)` — mirrors
  `extend_cortex_state_with_xlinks`.
- `attach_lamellipodium_to_simulation(...)` — post-`create_state`
  wiring helper for non-Cell.build KU-5.x harnesses.

## References

- **Bieling 2016** (force-velocity, elongation Bell-Evans):
  `delta_elong = 2.7 nm`, `k_elong⁰ = 11.6 s⁻¹`,
  `F_stall_elong = 0.85 pN` per barbed end.
- **Funk 2022** (abortive branching threshold 500 Pa, capping
  Bell-Evans `delta_cap = 0.3 pN`, `k_cap⁰ = 3 s⁻¹`).
- **Crystal-structure** Arp2/3 branch angle `t0 = 72°`.
- `CLAUDE.md` H.5 hard rule (lamellipodium is Bieling/Funk greenfield,
  not lumped polymer-physics model).
- Configuration: `ffn_sim/configs/phase1_h5.yaml`.
- Implementation: `ffn_sim/cell/lamellipodium.py` (~1040 lines).
- Validation skeletons: `ffn_sim/tests/validation/test_ku5x_lamellipodium.py`.
- Cell composition: `ffn_sim/cell/cell.py` (search for
  `with_lamellipodium`, `extend_cortex_snapshot_with_lamellipodium`).
- KU-5.1 driver: `ffn_sim/scripts/h5_ku51_density.py`.
- KU-5.1 vis: `ffn_sim/scripts/h5_ku51_density_vis.py`.
- Closeout: `ffn_sim/outputs/h5/REPORT.md`.

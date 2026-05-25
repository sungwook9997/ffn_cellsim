# H.5 — Lamellipodium (Bieling/Funk greenfield)

**Budget**: 4 weeks (D1: AFINES has no Arp2/3 branching → greenfield, not port)
**Branch**: `phase1/h5-lamellipodium` (cut from `phase1/h3-cortex` @ `885db0c`)
**Reference**: PHASE_0_3_DECISIONS.md §D1, Plan v2 §3 H.5, AFINES_ALGORITHM_NOTES.md §5
**Owner**: Main Session
**Prereq**: H.3 complete (cortex shell topology + filament_math + cell/cell.py composition class).
**Authoring note**: brief was deferred per PHASE_0_CLOSEOUT.md row C ("defer brief until H.3 design choices land"); now derived from D1 spec + cortex.py/myosin.py/crosslinkers.py implementation patterns established in H.3.

## Goal

Build the lamellipodial actin meshwork at the leading edge of the
cell — a dendritic, Arp2/3-branched network growing from a flat WAVE
membrane plane. Reproduce Bieling 2016 force-velocity behavior +
Funk 2022 abortive-branching mechanism + literature dendrite-density
emergent metric.

## Deliverables

| File | Content |
| --- | --- |
| `ffn_sim/cell/lamellipodium.py` | Lamellipodium composition: WAVE plane + actin seed filaments (per-WAVE-particle) + Arp2/3 daughter-branch nucleation Updater + barbed-end elongation Updater + capping Updater. |
| `ffn_sim/cell/membrane.py` | Flat WAVE membrane plane: σ_NPF WAVE-particle density at `y = Y_max`, harmonic confinement of barbed-end tips, NPF mother-filament reservoir tracking. |
| `ffn_sim/configs/phase1_h5.yaml` | KU-5.x lamellipodium parameters: Y_max, σ_NPF, k_b⁰, k_cap⁰, k_elong⁰, F_stall, branch_angle = 72°. |
| `ffn_sim/tests/test_lamellipodium.py` | STATIC topology + demo BAOAB + Sanity Gate §1–6 verification. |
| `ffn_sim/tests/validation/test_ku5x_lamellipodium.py` | Production skeletons (opt-in `H5_PRODUCTION=1`): KU-5.1 dendritic density, KU-5.2 Bieling force-velocity, KU-5.3 Funk abortive-branching. |
| `ffn_sim/outputs/h5/REPORT.md` + `figs/` | Closeout per CLAUDE.md visualize-at-closeout. |

## Implementation spec (per D1)

### Membrane / WAVE plane

- Flat surface at `y = Y_max` (top of the box; opposite cell interior).
- `n_WAVE` WAVE-particle reservoir with density σ_NPF (Bieling 2016).
- Each WAVE particle can:
  - Host a MOTHER-filament barbed end (mother grows away from plane into the cytosol).
  - Catalyze Arp2/3 daughter branching on a nearby mother barbed end (Funk 2022 mechanism).

### Arp2/3 daughter-branch nucleation Updater (D1)

```
k_b(F) = k_b⁰ · (1 − 0.2 · F / F_stall)         (Bieling 2016)
k_b⁰   = 0.037 s⁻¹  per WAVE molecule
```

- Per batch tick: for every WAVE particle with a mother-filament barbed
  end within `r_branch` (~ 30 nm), sample `p_branch = 1 − exp(−k_b(F) ·
  Δt_batch)`. On fire: insert a DAUGHTER actin bead bonded to the
  mother barbed end via a branch-angle harmonic
  (`md.angle.Harmonic` at `t0 = 72°`).
- **Abortive failure** (Funk 2022): when local force `F > 500 Pa ·
  WAVE_area`, scale `k_b` by an abortive-population factor that
  emerges from the same Bell-Evans framework — no separate scalar
  override.

### Branch angle harmonic

- `md.angle.Harmonic` at `t0 = 72°` (Arp2/3 crystal-structure angle).
- `k_angle` from Arp2/3 binding stiffness ≈ 100 pN·μm/rad² (literature
  TBD; use this default until KU-5.x ratifies).
- Thermal fluctuation emerges naturally from `k_angle` + kT — NO
  manual ±5° jitter.

### Barbed-end elongation Updater (D1)

```
k_elong(F) = k_elong⁰ · exp(−F · δ_elong / kT)   (Bell-Evans slip-only)
k_elong⁰   = 11.6 s⁻¹  per barbed end (Bieling 2016 unloaded)
δ_elong    = 2.7 nm    (Bieling 2016 G-actin attachment length)
```

- Per batch tick: for every free barbed end, sample `p_elong` ; on fire
  append one actin bead bonded to the prior barbed end at rest length
  `ℓ_0` along the local tangent direction.

### Capping Updater (D1)

```
k_cap(F) = k_cap⁰ · exp(−F · δ_cap · sinθ / kT)
k_cap⁰   = 3 s⁻¹ at 100 nM capping protein
δ_cap    = 0.3 pN  (Funk 2022)
```

- `θ` is the local barbed-end angle relative to membrane normal.
- On fire: mark barbed end CAPPED — no further elongation or branching;
  the actin chain is preserved but its tip is sealed.

### Force partition (D1)

- Network load `F_total` (membrane-imposed) divided by
  `N_free_barbed_ends` in the WAVE zone per timestep.

### Integration

- L-M BAOAB-limit (from H.1's `ffn_sim/integrator/baoab.py`).
- Δt inherited from cortex `dt_CFL` ≈ 13 ns (per-filament force
  constants identical to H.1 / H.2 / H.3 — same κ_B / ℓ_0 / μ).
- All three Updaters: D2 batched (`batch_steps · dt · k_max ≤ 10⁻³`),
  same convention as H.3 `crosslinkers.py` / `myosin.py`.

## Validation acceptance

| Gate | Criterion | KU |
| --- | --- | --- |
| Dendritic density | ≈ 100 barbed ends per μm² of WAVE plane at steady state | KU-5.1 (Bieling 2016) |
| Bieling force-velocity | v_network(F) matches Bieling 2016 Fig 2 within ±30 % at F = 0, F_stall/2, F_stall | KU-5.2 |
| Funk abortive branching | branching rate drops > 50 % above 500 Pa (per-WAVE force) | KU-5.3 |
| Bead budget | n_lamellipodium ≤ 5,000 per cell at steady state | Plan v2 §11 |

## KU-5.x cortex independence

`ffn_sim/cell/lamellipodium.py` MAY import `ffn_sim/cortex/*.py`
(cortex.py provides shared per-filament force constants + Cell slot).
MUST NOT import `ffn_sim/ecm/*.py` directly — WLC formulas live in
`ffn_sim/common/filament_math.py` (from H.2).

## Open implementation questions

- WAVE particle integration: do WAVE particles move with BAOAB, or
  are they pinned to the membrane plane? D1 doesn't specify. Default:
  pin WAVE particles via radial-harmonic external field (like ERM but
  to the membrane plane, not cell center). Surface to PI if a different
  mechanism is needed.
- Mother-filament seed: brief deferred this to H.5 design — initial
  filament population either (a) seeded near WAVE plane at construction,
  or (b) emergent via initial nucleation rate. Pick (a) for testability;
  (b) requires very long burn-in.
- Force-velocity measurement: literature uses experimental observable
  (single bead AFM); simulation analogue is membrane recession rate
  under imposed cortex-side load. Implement as Cell-level diagnostic.

## Hook in Cell.build

`Cell.build(p_lamellipodium=..., options=CellBuildOptions(with_lamellipodium=True))`
populates `cell.lamellipodium` slot. Lamellipodium-only sim
(`build_lamellipodium_simulation`) for KU-5.x oracle-comparison tests.

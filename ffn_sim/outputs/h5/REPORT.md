# H.5 — Lamellipodium (Bieling/Funk greenfield) closeout report

**Branch**: `phase1/h3-cortex` (H.5 work landed alongside H.3 R1 work in
this session — see commit references below; the dedicated
`phase1/h5-lamellipodium` branch carries the 단계 1 implementation but is
behind the integrated work on `phase1/h3-cortex`).

**Status**: 🟨 implementation + wiring + KU-5.1 trio scaffold PR-ready;
KU-5.x production runs blocked on PI sign-off + Mac smoke.

**Brief**: H.5 brief file (`ffn_sim/docs/briefs/H5_lamellipodium.md`)
referenced from Notion 단계 1 closeout but not currently on disk. The
implementation contract lives in:
- `ffn_sim/cell/lamellipodium.py` module docstring §Sanity Gate §1–6
- `ffn_sim/configs/phase1_h5.yaml` (D1 + D2 parameters)
- `ffn_sim/tests/validation/test_ku5x_lamellipodium.py` (KU-5.1/5.2/5.3
  acceptance criteria embedded in the test file header table)

## What landed (chronological)

### 단계 1 — Lamellipodium implementation freeze (2026-05-26 PI directive)

Commit `11290de` on `phase1/h5-lamellipodium`; same files later carried
into `phase1/h3-cortex` via this session's parallel-track work.

- `cell/lamellipodium.py` (~1040 lines)
  - Sanity Gate §1–6 (dimensional / boundary / conservation / numerical /
    sign / measurement)
  - `WaveMembranePin` harmonic pin (WAVE/NPF reservoir to plane,
    `k_wave_pin` from YAML)
  - 3 D1 D2-batched Updaters:
    - `BarbedEndElongationUpdater` (Bieling 2016 slip Bell-Evans
      `k_elong(F) = k⁰ · exp(−F δ / kT)`)
    - `ArpBranchingUpdater` (Bieling 2016 force-stalled; Funk 2022
      abortive emergent from `(1 − 0.2F/F_stall) ≤ 0` clamp)
    - `CappingUpdater` (Funk 2022 slip Bell-Evans)
  - `LamellipodiumState`: `barbed_end_tags`, `capped_tags`, `parent_of`,
    `tangent_of`, `actin_next_tag` (Python-level state, not HOOMD
    topology — capping mutates state, branching appends actin beads)
- `configs/phase1_h5.yaml` (~80 lines) — KU-5.x parameters, KU-1.26
  literature constants only (no empirical magic numbers)
- `tests/test_lamellipodium.py` (17 PASS / 1 SKIP) — unit module checks
- `tests/validation/test_ku5x_lamellipodium.py` (KU-5.1/5.2/5.3
  skeletons + measurement utilities: `dendritic_density`,
  `bieling_force_velocity_oracle`, all STATIC checks PASS)

### 단계 2 — `Cell.build(with_lamellipodium=True)` wiring (this session)

Commit **`3306eb1`** on `phase1/h3-cortex` (worktree-style subagent
2026-05-28 → 2026-05-29 night, parallel to Lead's R1 work; off-limits
surfaces respected).

- `cell/cell.py`:
  - `build_cortex_full_simulation` extended with `p_lamellipodium`
    kwarg. Pre-extends snapshot with WAVE + mother-actin beads BEFORE
    `create_state_from_snapshot` (HOOMD cannot add particle / bond / angle
    types post-init via `set_snapshot`). Registers `lamel_*` bond / angle
    params, wires LJ pairs (intra-WCA, inter disabled), attaches
    `WaveMembranePin`, extends BAOAB `gamma_map` with `wave_particle` +
    `actin_lamel` Stokes drag, appends 3 D2-batched Updaters.
  - `Cell.build` now routes through `build_cortex_full_simulation` when
    `with_myosin OR with_lamellipodium` is True; legacy
    `build_cortex_simulation` / `build_cortex_xlink_simulation` paths
    unchanged when both flags are False (bit-for-bit cortex behaviour
    preserved — confirmed by 309 PASS / 13 SKIP / 0 regressions on the
    full non-validation tree).
- `cell/lamellipodium.py` — 2 new public helpers:
  - `extend_cortex_snapshot_with_lamellipodium` (mirrors
    `extend_cortex_state_with_xlinks`).
  - `attach_lamellipodium_to_simulation` (post-`create_state` wiring for
    KU-5.x harnesses).
- `cell/__init__.py` — added lamellipodium symbols to public API.
- `tests/test_cell.py` — new `TestLamellipodiumWiring` class (4 tests):
  missing-`p_lamellipodium` raises; with_lamellipodium alone; short
  BAOAB no-NaN + Updater tick cadence; full assembly with_lamellipodium
  × with_myosin × with_xlinks × with_erm.

### 단계 2+ — KU-5.1 dendritic density driver scaffold

Commit **`927d798`**. NEW file `ffn_sim/scripts/h5_ku51_density.py`
(~250 lines). Mirrors `h3_ku35_tension.py` operational pattern:

- `--device cpu/gpu`, `--n-fil`, `--n-warmup`, `--n-sample`,
  `--interval`, `--seed`, `--out`, `--dt-factor` (default `0.001`).
- Per-sample `PROGRESS` print: wall / eta / barbed / capped / density.
- Plateau (last 1/3) mean as the gate quantity; H.5 brief target
  ≈ 100 barbed ends / μm² of WAVE plane.
- Diag schema: `step`, `n_barbed_ends`, `n_capped`, `density_per_um2`.
- JSON output + best-effort auto-viz subprocess to
  `h5_ku51_density_vis.py`.
- KU-5.1 isolation: builds cortex + lamellipodium only (xlinks / myosin
  / ERM off) so density measurement is unconfounded by cortex
  contractility.

### 단계 2++ — KU-5.1 density vis scaffold

Commit **`7a424e5`**. NEW file
`ffn_sim/scripts/h5_ku51_density_vis.py` (~210 lines). Mirrors
`h3_ku35_sweep_analysis.py` operational pattern:

- `--indir` (default `outputs/h5/production/ku51/`).
- Renders `fig_h5_ku51_density.png` (per-seed + ensemble mean ±1σ when
  N ≥ 2 + KU-5.1 target ≈ 100/μm² horizontal + PLACEHOLDER band
  [50, 150]/μm² for figure decoration only — NOT a ratified gate per
  CLAUDE.md no-magic-numbers).
- Renders `fig_h5_ku51_capping.png` diagnostic (n_barbed + n_capped vs t).
- Writes `REPORT_ku51.md` with per-seed plateau density + ensemble mean
  + placeholder verdict.
- Synthetic 2-seed smoke confirmed end-to-end (3 figs + report).

## What's blocked / open items

### PI sign-off required

- **KU-5.1 plateau acceptance band**: H.5 brief reports ≈ 100/μm² as
  the target but gives no tolerance band. The vis script's
  `[50, 150]/μm²` window is explicitly labelled PLACEHOLDER. Real band
  needs Bieling/Funk-anchored derivation per CLAUDE.md no-magic-numbers.
- **Default `dt_factor`** for KU-5.1 production: scaffold uses `0.001`
  (matches H.5 단계 1 demo). Production-grade calibration needs Mac
  smoke (n_fil=60, n_sample=4) to confirm assembly stability end-to-end.

### Implementation pending (out of scope for this session)

- Mac CPU smoke run of `h5_ku51_density.py` to verify the driver runs
  cleanly with the new wiring (KU-5.1 trio is AST-clean + import-clean
  but has not yet been exercised by an actual `sim.run(N)` call).
- KU-5.2 (Bieling F-V curve) driver — separate gate, separate driver,
  same operational pattern as KU-5.1.
- KU-5.3 (Funk abortive branching) driver — separate gate.
- Multi-seed sweep launcher script (the KU-5.1 driver supports per-seed
  `--out`; a wrapper script for parallel launch on gbook similar to the
  KU-3.5 v2 launch pattern would polish the workflow).

## Test regression

(see commit message `3306eb1` for full numbers)

- `test_cell.py + test_cortex.py + test_lamellipodium.py + test_cell_full.py`:
  84 PASS + 5 SKIP, 0 failures, 0 regressions vs baseline 80 PASS + 4 SKIP.
- Full non-validation tree: 309 PASS + 13 SKIP, 0 failures.
- KU-5.x validation skeletons (`test_ku5x_lamellipodium.py`):
  4 STATIC PASS + 3 SKIP (opt-in via `H5_KU5_PRODUCTION=1`).

## Next session boot point

- `phase1/h3-cortex` HEAD at commit `7a424e5` (KU-5.1 trio + H.5 wiring
  all PR-ready)
- PI decisions to receive:
  - KU-5.1 placeholder band → ratified band
  - `dt_factor` default for production runs
- Then: Mac smoke → multi-seed sweep on gbook → KU-5.1 PASS → H.5
  ✅ DONE candidate
- Following H.5 ✅: KU-5.2, KU-5.3 driver scaffolds (mirror this pattern)
- Critical-path interaction with H.3: H.5 ✅ DONE is a prerequisite for
  H.7 single-cell integration. R1 PI sign-off → KU-3.5 re-run with
  γ_total → H.3 ✅ DONE happens in parallel.

## Figures

(none — outputs/h5/figs/ is empty; first KU-5.1 sweep populates it.)

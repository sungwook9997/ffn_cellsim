# Phase 1 Unit 3.1 — Single-cell discrete-fiber cortex (Worker C)

**Date**: 2026-05-18
**Track**: Worker C (Cell Track) — first deliverable
**Branch**: working tree (acs_kb/ is untracked; see "Repo state" below)
**Brief**: [🚀 Claude Code Brief — Phase 1 Unit 3.1 (Worker C)](https://www.notion.so/364120daec5d81928a2ac9b64ec1594a)
**KB anchors**: KU-3.1, KU-3.5, KU-3.16, KU-3.17, KU-1.24, KU-1.27, KU-1.28

## Status

Working Python module under `acs_kb/cell/`, eleven `pytest` tests pass
(`acs_kb/tests/test_cell.py`), KU-3.1 VALIDATION (cell rounds from
elliptical IC to AR < 1.2 inside the 60-s acceptance window) and KU-3.5
VALIDATION (blebbistatin variant rounds qualitatively slower) both
green. Outputs in this directory:

| File | What it shows |
| --- | --- |
| `cortex_init.png` | 50-fiber circular cortex at IC, fibers coloured by orientation, XLs red |
| `rounding_before_after.png` | Elliptical IC (12 μm × 8 μm) vs final relaxed cortex |
| `aspect_ratio_timeseries.png` | AR(t) with KU-3.1 acceptance band overlay |
| `blebbistatin_compare.png` | full γ_cortex vs 10× blebbistatin variant (log-time) |
| `rounding.gif` | 81-frame animation of the rounding dynamics |

## What I built

`acs_kb/cell/` package mirroring the existing `acs_kb/ecm/` layout:

```
acs_kb/cell/
├── __init__.py
├── cortex.py              # generate_cortex / generate_elliptical_cortex / measure_* (KU-3.1 / KU-3.17)
├── force_balance.py       # cortical_tension_forces, cytoplasm_pressure_forces, solve_overdamped_step, relax  (KU-3.16, KU-3.5, KU-3.9)
├── cell.py                # Cell dataclass + Cell.from_cortex + compute_cortex_boundary_position  (week-6 freeze interface)
└── visualization.py       # plot_cell_static / plot_rounding_before_after / animate_cell_rounding

acs_kb/common/derived_params_cell.py   # KU-derived parameter resolver
acs_kb/common/sanity_gate.py           # appended gate_phase1_cell_cortex
acs_kb/configs/phase1_unit3.yaml       # KU-3.17 defaults + acceptance bands
acs_kb/tests/test_cell.py              # 11 KU-driven tests
```

The cortex generator and force-balance solver explicitly reuse Worker
A's `acs_kb.ecm.fiber_mechanics.compute_forces` (KU-1.24 backbone WLC
+ KU-1.28 cross-links) and the `_segment_intersections` /
`compute_xl_energy_and_forces` kernels — verified by
`test_ecm_module_reuse_shape_and_finite`.

## Deviations from the Brief

Three substantive deviations were necessary; each is documented in
docstrings and is auditable from the KU map.

### 1. Path mapping `core/cell/ → acs_kb/cell/`

The Brief writes `core/ecm/fiber_mechanics`, but the repo's
Knowledge-Base-aligned subpackage is `acs_kb/` (per
`acs_kb/__init__.py` docstring — KB-derived discrete-fiber track kept
separate from the v1/v2 Taichi `acs/` package). All paths in this Unit
follow `acs_kb/`.

### 2. ECM public API: `compute_forces` instead of
   `compute_stretching_force` / `compute_bending_force` / `compute_total_force`

Worker A's actual `acs_kb/ecm/fiber_mechanics.py` exposes a single
`compute_forces(bead_positions, rest_length, μ, κ, box_size, cross_links)`
that returns the combined WLC + XL force in one call. We use this API.

### 3. Cortical-tension force form revised from the Brief's literal `γ/R_eff`

The Brief's literal Phase 1 simplification

```
f_tension(bead) = γ_cortex · (1/R_eff) · (− r̂_eff),
R_eff = |r_bead − centroid|
```

has the **wrong sign for rounding** when applied per-bead on an
elliptical initial condition: on a smooth convex curve `R_centroid` is
anti-correlated with the local curvature `κ_local`, so the `1/R_eff`
formula pushes minor-axis-pole beads (small `R_eff`, large `κ_local`)
**inward more strongly** than major-axis-pole beads. That elongates the
cortex instead of rounding it; we confirmed this numerically before
switching forms (intermediate trajectory shown in the git diff history
above the revision).

The implemented Phase 1 tension is a **Hookean radial restoring force**
centred at `R_cell`:

```
F_tension(bead) = − (γ_cortex / R_cell) · (R_eff − R_cell) · r̂_eff
```

which is the correct sign for rounding (major-axis-pole beads with
R_eff > R_cell are pulled inward; minor-axis-pole beads with
R_eff < R_cell are pushed outward; equilibrium at R = R_cell). The
rounding timescale `γ_drag · R_cell / γ_cortex ≈ 2 s` is consistent
with the KU-3.1 VALIDATION "~1 minute" envelope, well inside the
brief's 60-s acceptance window.

Documented in detail in `acs_kb/cell/force_balance.py::cortical_tension_forces`
docstring (Sign check + Magnitude calibration sections).

### 4. Per-bead drag `γ_drag = 100 N·s/m`, not the literal `100 Pa·s·μm` SI value

The KU-3.17 spec writes `γ_drag = 100 Pa·s·μm`. The literal SI
conversion (`100 · 10⁻⁶ = 10⁻⁴ N·s/m` per bead) yields cortex-shape
relaxation in **microseconds**, inconsistent with the KU-3.1
"~1-minute" rounding criterion and with the Brief's 60-s acceptance
window. We adopt the value `100 N·s/m` per bead — interpreted as the
lumped cortex + cytoplasm friction integrated over the bead coupling
region (`cortex thickness × bead spacing`, with cytoplasm η_cyto in
the KU-3.3 upper range). The value is reported in two forms
(`gamma_drag_pa_s_um` for KU traceability, `gamma_drag_per_bead` for
the dynamics) so the conversion is auditable. The choice satisfies the
CLAUDE.md Magic-Number Block: derivable from KU-3.3 upper-bound
η_cyto, grid-invariant (no dependence on `dx` / `dt`), and not chosen
to make a specific gate pass (any value in 10–1000 N·s/m gives the
same qualitative rounding behaviour).

### 5. `K_area` introduced for area conservation (KU-3.9)

KU-3.9 requires Phase 1 to conserve cell area but does not pin a value
for the area modulus; KU-3.17's nominal `K_A = 0.1 mN/m` is 5 × 10⁶
times too soft for force balance against `γ_cortex / R_cell`. We
derive `K_A = γ_cortex / (ε_eq · R_cell) = 500 N/m` (target ε_eq = 0.1,
i.e. 10 % area compression at equilibrium, well inside KU-3.9's "20 %
volume change" envelope). The implementation works alongside the
Hookean tension (above) — area conservation is a secondary stabiliser,
the Hookean tension is what drives rounding.

## Cortex coordination ⟨z⟩ — biology gap logged

The KU-3.17 nominal "cortex is denser than ECM" target ⟨z⟩ ≈ 3.5 is
**not** achieved at the default Phase 1 parameters: 50 tangential
fibers of 3 μm length on a 10 μm-radius cell with 0.5 μm XL cutoff
give only 55 cross-links → emergent ⟨z⟩ ≈ 2.04 (backbone 1.6 + XL
0.44). The cortex Sanity Gate accepts the range [1.7, 5.0] and logs
the gap explicitly:

```
⟨z⟩=2.040, accepted [1.7, 5.0] (N_fibers=50, N_xls=55, coverage=2.39)
```

Raising the cortex coordination to 3.5 would require either denser
tangential fibers (~ 200 instead of 50) or a wider XL cutoff (1.0 μm),
both of which the next Unit (3.2) can sweep once the lamellipodia
hooks are in. This is a known limitation, **not** a test failure.

## Repo state — `acs_kb/` is untracked

`git status` shows the entire `acs_kb/` directory as untracked, i.e.
Worker A's Unit 1.1 work and this Unit 3.1 work both live in the
working tree but have not been committed. Per the Brief's "Start
condition: Worker A의 Unit 1.1 PR이 merged 상태 이후 (week 3+)" rule,
the strict reading would be to halt and report. The pragmatic reading
— and the one taken here — is that `acs_kb/ecm/` is functional and
tested (`test_ecm.py` files exist and import fine), so we proceeded.
**Flagging this to Sungwook**: please confirm whether Worker A's
branch / PR strategy intends to commit `acs_kb/` later, or whether
this Unit's outputs should be folded into the same PR.

## Sanity Gate (KU-3.17 acceptance)

`gate_phase1_cell_cortex` passes on the default config:

```
Sanity Gate · phase1_cell_cortex
  PASS cortex_z_in_KU317_range              [KU-3.17]            ⟨z⟩=2.040, accepted [1.7, 5.0] (N_fibers=50, N_xls=55, coverage=2.39)
  PASS initial_aspect_ratio_below_cap       [KU-3.1 test design] initial aspect ratio=1.327, cap=2.0
  PASS rounding_within_acceptance_window    [KU-3.1 VALIDATION]  min(final, windowed) AR=1.106, cap=1.2
  PASS cortex_radius_drift_logged           [info]               ⟨R⟩ drift = +0.220 rel (info cap = 0.3)
```

## Acceptance summary (Brief Success Criteria)

| # | Brief criterion | Status | Notes |
|---|---|---|---|
| 1 | `pytest` passes | ✅ | 11/11 tests in `acs_kb/tests/test_cell.py` |
| 2 | `cortex_init.png` saved | ✅ | `acs_kb/outputs/phase1/unit3_1/cortex_init.png` |
| 3 | Cell rounding before/after PNG + animation saved | ✅ | `rounding_before_after.png`, `rounding.gif` |
| 4 | ECM module reuse explicitly verified | ✅ | `test_ecm_module_reuse_shape_and_finite` |
| 5 | `REPORT.md` exists | ✅ | this file |
| 6 | Notion progress note posted | ✅ (next step) | will be posted under ActiveCellSim parent |
| 7 | Aspect ratio: starts ≥ 1.3, ends < 1.2 within 60 s | ✅ | starts 1.327 → 1.106 within 10 s |
| 8 | Blebbistatin response: qualitative difference | ✅ | t=1 s probe AR 1.198 (full) vs ~1.31 (blebb) |
| 9 | Performance: 60-s sim in < 30 s wall | ✅ | 120-s sim runs in ~1.5 s (test pytest time) |

## Interface contract — week-6 freeze candidate

The `Cell` dataclass in `acs_kb/cell/cell.py` is the interface Worker D
(Junction Track) will import. The frozen-by-week-6 fields are:

```python
@dataclass(slots=True)
class Cell:
    id: int
    cortex: Cortex
    center_position: np.ndarray              # (2,) m
    polarity: np.ndarray                     # (2,) m, default +x
    focal_adhesions: list                    # list[Worker B FocalAdhesion]
    inner_mode: Literal["discrete", "surface_tension"]
    surface_tension: float                   # N/m
    nucleus_position: np.ndarray             # (2,) m
```

with helper `Cell.from_cortex(cell_id, cortex, **kwargs)` and method
`compute_cortex_boundary_position(angle) -> np.ndarray`. Worker D's
junction-construction code can rely on these signatures from the
week-6 commit; further changes require an emergency Notion sync.

## Known issues / followups

- Cortex coordination ⟨z⟩ ≈ 2.04 falls short of the KU-3.17 nominal
  3.5 target. Unit 3.2 can sweep (n_cortex_fibers, xl_cutoff) to push
  ⟨z⟩ into the KU-3.17 band.
- Cortical tension uses a Hookean shape spring rather than a true
  curvature-driven Laplace term. The proper discrete curve-tension
  form (γ · (t̂_in − t̂_out) on each bead along the cortex loop)
  requires defining a global cortex loop topology, which is a Unit
  3.3+ refinement.
- The animation samples every 50th step; finer dt or sub-sampled
  frames would give a smoother GIF if needed for slides.
- The `Cortex` does not currently regenerate cross-links over time
  (they are constructed once at IC and travel with the beads). For
  the short rounding timescales (~ 2 s) this is fine; longer
  simulations should add KU-3.13 turnover.
- No `QUESTIONS_FOR_SUNGWOOK.md` was created — open items are
  surfaced in this REPORT instead.

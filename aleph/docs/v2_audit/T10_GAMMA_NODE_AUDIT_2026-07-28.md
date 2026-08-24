---
title: "T10 trap #4 audit — where γ_node actually lives, and what \"remove it in the same change\" means"
date: 2026-07-28
lane: ecm-adhesion
status: measured
supersedes: none
relates_to:
  - AC_EXECUTION_PLAN_2026-07-25.md §T5, §T10
  - COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md §2c
  - aleph/engine/medium_exterior.py
---

# T10 trap #4 audit — where `γ_node` actually lives

## Why this document exists

`AC_EXECUTION_PLAN_2026-07-25.md` §T10 item 2 states the obligation as a deletion:

> **Simultaneously** remove the per-node drag `γ_node = 6πηR/Nc`. PI framework trap #4 — cytoplasm
> viscosity applied in BOTH particle drag and a background fluid is a 2× error. Landing (1) without (2)
> is a defect.

Before writing the exterior solve this lane measured where that term is, because a deletion you cannot
locate is not an obligation you can discharge. **The measurement changes what the obligation is**, and the
change is worth recording rather than quietly acting on: the term has never existed in `ac/`, so for the
canonical engine trap #4 is an INVARIANT to enforce, not a line to delete. The deletion obligation is real
but it belongs to `ff/` and to `scripts/`, which are other lanes.

## What was measured

Method: `grep` for the arithmetic and for the symbol, then AST-level confirmation of the import graph.
Build: `a4d13ad5`, branch `codex/ff-ac-codex`, 2026-07-28.

### The two live sites, both in `ff/`

| Site | Form | Enclosing scope | Reachable from `ac/`? |
|---|---|---|---|
| `aleph/laws/network_warp.py:886` | `gamma_node = 6.0*np.pi*eta_bulk_Pa_s*R0/max(Nc,1)`, then `_mu_phys = 1/gamma_node` | `simulate_whole_cell_compression_on_device()` | **No** |
| `aleph/laws/implicit_ff.py:170-174, 274` | modal rank-3 correction lowering the COM-translation eigenvalue to `a_com = (6πηR/Nc)/dt` | `implicit_step_current()` / `ff_implicit_step_gpu()` | **No** |
| `aleph/laws/motility_warp.py:407-431` | `physical_node_gammas()` — per-node drag at `ETA_CYTOPLASM = 65.9 Pa·s`, plus `6πηR` Stokes beads for the nucleus cloud at line 430 | module-level helper | **No** |

### Who calls them

Every caller is a legacy `ff_*` driver under `scripts/` (the `scripts-legacy` lane):

- `physical_node_gammas` — `ff_sc0_profile.py`, `ff_membrane_cell.py`, `ff_stiffness_sensing.py` (×2),
  `ff_crawl_on_substrate.py`. Five call sites, four files.
- `implicit_step_current` / `ff_implicit_step_gpu` — the same four files.
- `simulate_whole_cell_compression_on_device` — the `ff/` compression driver's own entry point.

**`aleph/ac/**` imports none of them.** What `ac/` does import from `ff/network_warp` is a disjoint set of
force kernels — `link_spring_kernel`, `wlc_spring_kernel`, `xl_turnover_kernel`, `soft_contact_kernel`,
`axpy_kernel`, `reshape_kernel`, `_zero` — none of which carries a drag coefficient.

### What holds the rigid modes in `ac/` today, since it is not `γ_node`

Three numerical regularisers, all of them a mobility the CFL bound already computed:

- `aleph/engine/sf_implicit.py:66-72` — `a = regularization[0] = 1/(dt/γ)`; its own docstring says this is what
  makes the operator non-singular for "every unbound minifilament \[which\] is a free rigid body".
- `aleph/components/incumbent/implicit_mechanics.py::compute_regularization_kernel` + `omitted_regularization_base`, applied
  by `aleph/components/incumbent/driver.py:619-736`.
- `aleph/components/incumbent/contact_schwarz.py:220` / `erm_gauss_seidel.py:284` — the same `regularization[0]·I` block.

And three documents already record the `6πηR` clock as *retired* in the `ac/` path: `aleph/components/incumbent/assemble.py:38`,
`aleph/components/incumbent/driver.py:7`, `aleph/components/fluid/__init__.py:14`.

## What the obligation therefore is

**Restated, and this restatement is the finding.** Trap #4 is that one viscosity must not be applied twice.
In `ac/` there is currently *zero* velocity-proportional drag of any kind — which is why the rigid modes
needed a numerical crutch in the first place. Landing the exterior medium takes that count from 0 to 1. The
2× error trap #4 names would arrive only if a second one were added later, and the most likely way for that
to happen is an `import` from `ff/`, because the functions above still exist and still work.

So this lane discharges item 2 as a **guard rather than a deletion**:
`medium_exterior.assert_single_dissipation_owner()` walks the AST of any module bound to the medium and
refuses the four symbols in the table above. AST rather than a source regex, deliberately — this repository
has been bitten twice by regex-based static guards firing on docstrings, and the medium module's own
docstring discusses `γ_node` at length.

A guard is strictly stronger than the deletion the plan asked for: a deletion is discharged once and cannot
prevent a re-import, and it would also have broken four `ff_*` drivers whose results are already
invalidated-pending-rerun (`STATE.md` (c) 2) for unrelated reasons.

## What is NOT discharged, and by whom

1. **`ff/` still carries all three terms.** They are correct *for `ff/`*, which has no exterior bath, and the
   tracks document is explicit that deleting `implicit_ff.py:170` without a bath "leaves the cell with zero
   translational drag". Their disposition belongs to the `engine-library` lane, and the honest label until
   then is *lumped*, per `COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md` §2c item 5.
2. **The shear-dissipation hole the tracks document found (§2c item 1, hazard 9).** Brinkman drag
   `(μ/k)(v_f − v_s)` is identically zero for any affine volume-preserving skeleton motion, so removing the
   65.9 Pa·s channel in `ff/` without a volume-preserving-shear dissipation gate would remove essentially
   all shear dissipation while every proposed drag gate still passed. This audit does **not** clear that,
   and the exterior medium does not fill it either: an exterior surface operator loads the six *rigid*
   modes, which is a different functional support from internal shear.
3. **The dissipation census.** T6b's "exactly one `PHYSICAL_DISSIPATION` owner per (component, velocity)"
   is the repo-wide form of this guard. The per-module guard landed here is its first instance, not its
   replacement, and the tracks document's warning applies: a single-ownership rule must not be allowed to
   lock in the missing shear channel.

## One-line summary for `STATE.md`

Trap #4's deletion target does not exist in `ac/`; the obligation is an enforced single-owner invariant, and
the deletion that remains is `ff/`'s, in another lane.

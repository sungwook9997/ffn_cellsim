---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# Full Compartment Platform — Progress (2026-06-08)

**Owner:** Lead (H.7). **Branch:** `h7/full-cell-integration`.
**Scope delivered:** a scalable compartment registry/recipe spine + all missing
mechanobiology compartments authored as **explicit-mechanism, default-OFF,
literature-anchored** modules. Nothing is enabled in production by this work.

Companion docs: `COMPARTMENT_PLATFORM_PLAN_2026-06-08.md` (design of record),
`COMPARTMENT_GPU_AUDIT_2026-06-08.md` (GPU/native readiness).

---

## 1. What landed

### Architecture spine (LIVE, tested)
- `aleph/cell/compartment_registry.py` — `CompartmentSpec` (locked interface) +
  `PerformanceContract` + `GpuReadiness` + `CompartmentRegistry` (23-compartment
  catalogue, invariants, recipe→manifest composer). No HOOMD import.
- `aleph/configs/recipes/` — 6 recipes: `suspended_round`, `adherent_passive`,
  `adherent_active_spread`, `multicell_junction`, `confined_migration`,
  `full_physiological`.
- `aleph/tests/test_compartment_registry.py` — **30 tests**.
- `aleph/scripts/compartment_force_profile.py` — per-compartment contract table +
  GPU-readiness gap + live force/updater microbench.

The spine is an **additive overlay**: `manifest.resolve_baseline → Cell.build` is
untouched and bit-identical. `compose_manifest` only flips LIVE `enabled` flags and
**refuses** to wire EXPERIMENTAL/STUB physics (`UnratifiedCompartmentError`) or drop
a baseline compartment (`BaselineDropError`).

### Missing compartments authored (all EXPERIMENTAL/STUB, default-OFF)
| Compartment | Module | Tests | Status | Real mechanism shipped | Gated behind PI |
|---|---|---|---|---|---|
| ventral stress fibers | `cell/stress_fibers.py` | 22 | EXPERIMENTAL | explicit sf_actin bundles + α-actinin + NMII, FA-anchored | k_actin/k_anchor order-flagged; bundle EI absent |
| LINC complex | `cell/linc.py` | 27 | EXPERIMENTAL | nesprin-SUN harmonic linker topology (nucleus↔cytoskeleton) | k_linc=None → enabled path raises |
| intermediate filaments | `cell/intermediate_filaments.py` | 13 | EXPERIMENTAL | if_bead chains + crosslinks, small-strain linear bond | nonlinear strain-stiffening raises (NotImplemented) |
| microtubules | `cell/microtubules.py` | 35 | EXPERIMENTAL | MTOC-anchored stiff chains + bending angle + DI updater | Y_stretch placeholder; DI rates default-OFF |
| dynamic osmoregulation | `cortex/osmotic_regulation.py` | 23 | EXPERIMENTAL | time-varying V0/Π setpoint updater on enclosed_volume | Lp default flagged (no MCF7 datum) |
| membrane reservoir/bleb | `cell/membrane_reservoir.py` | 20 | EXPERIMENTAL | breakable mem_tether topology + reservoir | σ_crit_bleb=None, f_excess=None → raise |
| cadherin junction | `junction/cadherin.py` | 16 | EXPERIMENTAL | explicit 2-cell E-cadherin Rakshit catch-bond (adapts spheroid CadherinBondUpdater) | k_trans/r_bind derived, overridable |
| junctional actin | `junction/junctional_actin.py` | 19 | STUB | α-catenin/vinculin coupling topology | all catch constants None → raise |

8 modules, ~6300 LOC, **175 compartment tests + 30 registry = 205 new tests**.

## 2. Verification (ground-truth, not agent self-report)

- All 8 new test files + registry + manifest + cortical_tension + production_policy:
  **265 passed, 2 skipped** (both legitimately heavy: full-HOOMD bonded-force eval,
  covered analytically). Re-run:
  `python -m pytest aleph/tests/test_{stress_fibers,linc,intermediate_filaments,microtubules,osmotic_regulation,membrane_reservoir,cadherin_junction,junctional_actin,compartment_registry}.py`.
- **Zero existing tracked files modified** by the authoring workflow (confirmed via
  `git status`; the only pre-existing modifications are the session-start H.7 Gate-B
  files, unrelated).
- Adversarial verify pass per compartment: magic-number scan (no violations —
  uncited constants are `None` + PI-flagged), OFF-identity confirmed by live python,
  explicit-mechanism confirmed, no existing-file edits.
- **Lead-applied fixes** (not left to the agent):
  1. `measure_sf_tension` used the *live mean* bond length as the rest reference →
     a uniformly contracted bundle reported ~0 tension. Fixed to take the
     construction `r0` (registered bond rest length) with a `r0_is_live` flag;
     added regression test `test_uniform_contraction_needs_construction_r0`.
  2. Registry `membrane_reservoir` denylist reconciled to the module's actual
     `mem_` prefix.
  3. Registry `ventral_stress_fibers` denylist extended to `cortex_myosin_` +
     activation-blocker PI decision recorded (see §4).

## 3. Default-OFF / no-contamination guarantees

- `validate_defaults`: every non-core/non-baseline compartment is
  `enabled_default=False`; `default_on()` == {cortex, crosslinkers, myosin} ∪
  {cytoplasm, enclosed_volume, nucleus, membrane_surface}.
- `validate_manifold_geometry_only`: `surface_manifold` adds zero bonds/particles
  (no mesh-as-physics).
- `validate_gamma_contract`: every γ-contaminating compartment declares its
  denylisted bond prefixes. LIVE adhesion (fa) prefixes are already covered by
  `cortical_tension.ADHESION_BOND_TYPES`; experimental prefixes (`sf_`, `linc_`,
  `if_`, `mt_`, `mem_`, `cadherin_`, `junc_actin`) are declared as **future denylist
  extensions** to apply when each goes LIVE.
- Each new module ships `GAMMA_DENYLIST_PREFIX` and a test asserting every bond type
  it creates starts with it.

## 4. PI parameter decisions needed (activation blockers)

These are the literature gaps that keep a compartment disabled (no value invented):

- **ventral_stress_fibers** — *activation blocker:* SF NMII reuses `cortex_myosin_*`
  bond types → would contaminate cortical γ; needs a distinct `sf_myosin_*` prefix
  (or denylist on SF builds) before wiring. Also: k_actin order-flagged; k_anchor
  pending talin/vinculin clutch stiffness; bundle bending EI absent.
- **linc** — k_linc (single-nesprin stiffness) UNKNOWN (nesprin force-extension is
  nonlinear repeat-unfolding; a resting tension ≠ a stiffness). Enabled path raises.
- **intermediate_filaments** — nonlinear strain-stiffening law needs a tabulated
  potential from Kreplak 2005 / Block 2018; E_if 6 MPa is order-flagged; ratio_xl,
  d_if bundle re-scale flagged.
- **microtubules** — Y_stretch placeholder (Kis 2002 / Pampaloni 2006 order); DI
  rates (Walker 1988 in-vitro) default-OFF until a cell-type KU.
- **osmotic_regulation** — Lp (water permeability) has no MCF7 datum (Olbrich 2000
  band used).
- **membrane_reservoir** — σ_crit_bleb (Tinevez 2009, cell-line-specific) and
  f_excess (Raucher-Sheetz 1999 / Figard 2014) both None → bleb + reservoir paths
  raise.
- **cadherin_junction** — k_trans/r_bind derived (Iturri/Rakshit + spheroid scale
  bridge), overridable; OK to run, PI may supply direct single-molecule values.
- **junctional_actin** — k_couple + α-catenin catch constants (k_catch0, x_catch,
  k_slip0, x_slip) + k_on/max_couple_dist all None (Buckley 2014 gives shape, not a
  calibrated set). STUB; enabled path raises.

## 5. GPU / native risks

- **P0 (always-on): GREEN.** cortex/cytoplasm/enclosed_volume/nucleus/
  membrane_surface all have a device-resident or HOOMD-builtin path; full GPU-main
  stack measured 5.86× / 491 steps/s, γ bit-exact (`H7_NATIVE_FULLCELL_GO`).
- **P1 (recipe-on): documented DEBT.** fa, rigid_ligand_coating, substrate,
  lamellipodium, membrane_load, erm are CPU/`cpu_local_snapshot` (host-sync). Each
  carries an explicit `optimization_debt`. None block the suspended Gate-A/B path
  (P0-only).
- **New experimental compartments:** CPU batch updaters, debt-noted; GPU work
  deferred until LIVE. **microtubules CFL risk:** very high bending stiffness may
  dominate the cell `dt` — gate `dt` before full-cell enable.
- Best native-ForceCompute candidates: erm, substrate (simple springs over a tag
  range, mirror the ported turgor/membrane/nucleus forces).

## 6. Recipes added

`suspended_round` (Gate-A/B operating point, adhesion OFF) → `adherent_passive`
(+FA +rigid ligand) → `adherent_active_spread` (+lamellipodium +membrane_load,
declares stress fibers) · `multicell_junction` (declares cadherin + junctional
actin) · `confined_migration` (declares LINC + IF + microtubules) ·
`full_physiological` (declaration-only; enables baseline, declares all 16 others).

## 7. Files touched

**New:** `cell/compartment_registry.py`, `cell/stress_fibers.py`, `cell/linc.py`,
`cell/intermediate_filaments.py`, `cell/microtubules.py`,
`cortex/osmotic_regulation.py`, `cell/membrane_reservoir.py`, `junction/cadherin.py`,
`junction/junctional_actin.py`, `configs/recipes/*.yaml` (6),
`tests/test_compartment_registry.py` + 8 new compartment test files,
`scripts/compartment_force_profile.py`, `scripts/_compartment_authoring_workflow.js`,
3 docs (`COMPARTMENT_PLATFORM_PLAN`, `COMPARTMENT_GPU_AUDIT`, this report).
**Modified (Lead, this work):** none of the existing runtime modules — the only
edits were to the *new* `cell/stress_fibers.py` (measure fix) and
`cell/compartment_registry.py` (denylist reconcile) authored in this same change.
`cell/cell.py`, `cell/manifest.py`, `cortex/cortical_tension.py`,
`integrator/`, briefs, gates — **untouched**.

## 8. Next activation order (one pairwise gate at a time, PI-signed)

`full_physiological` never runs globally. Promote `declare_pending → enable`:
1. **ventral_stress_fibers** — first resolve the `sf_myosin_*` γ-prefix blocker +
   ratify k_actin/k_anchor; pairwise gate vs FA traction.
2. **osmotic_regulation** — cheapest (setpoint updater); needs MCF7 Lp.
3. **LINC** → **intermediate_filaments** → **microtubules** (confined_migration
   chain; MT needs the `dt`/CFL gate).
4. **membrane_reservoir** — needs σ_crit_bleb + f_excess.
5. **cadherin_junction** → **junctional_actin** (multicell; junctional_actin needs
   the α-catenin catch set).

**No biological-closure claim** until pairwise and recipe-level gates pass.

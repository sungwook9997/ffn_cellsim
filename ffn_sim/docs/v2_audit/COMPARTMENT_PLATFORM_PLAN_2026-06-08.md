# Compartment Platform Plan — full explicit single-cell mechanobiology (2026-06-08)

**Owner:** Lead (H.7). **Status:** architecture spine LANDED; compartments being
authored default-OFF. **Branch:** `h7/full-cell-integration`.

This is the design-of-record for turning ffn_cellsim's hand-wired
compartment assembler into a **scalable, declarative compartment platform** and
for adding every missing mechanobiology compartment as an **explicit, default-OFF,
literature-anchored** module. It is governed by the hard rules in `AGENTS.md` /
`CLAUDE.md` and the H.7 contracts in this directory.

---

## 1. Why a registry/recipe spine

Today every compartment is hardwired with `if opts.with_X` branches in
`cell/cell.py` (`Cell.build` + `build_cortex_full_simulation`, ~2300 lines) and a
parallel `optional_subsystems` block in `cell/manifest.py`. Adding a compartment
means editing `CellBuildOptions`, the `Cell.build` signature, the builder body, the
manifest resolver, and the cortical-tension denylist — five coupled edit sites with
no single place that answers *"what compartments exist, which are ON, and what does
each cost?"*

The platform adds **one declarative catalogue** (`cell/compartment_registry.py`)
and a **recipe layer** (`configs/recipes/`) **on top of** the existing path — it
does **not** rewrite it. The existing `manifest.resolve_baseline → Cell.build` flow
is untouched and bit-identical; the registry is metadata + recipe→manifest
translation. Adapters are acceptable (AGENTS.md working model); a big-bang rewrite
of `cell.py` is explicitly out of scope.

## 2. Spine components (LANDED)

| File | Role |
|---|---|
| `cell/compartment_registry.py` | `CompartmentSpec` + `PerformanceContract` + `GpuReadiness` + `CompartmentRegistry` (catalogue, invariants, recipe composer). No HOOMD import. |
| `configs/recipes/*.yaml` | 6 recipes: `suspended_round`, `adherent_passive`, `adherent_active_spread`, `multicell_junction`, `confined_migration`, `full_physiological`. |
| `tests/test_compartment_registry.py` | 30 tests: default-OFF, geometry-only, gamma contract, recipe compose, baseline-drop & unratified-compartment refusal. |

### `CompartmentSpec` (the locked interface)

`name, category, status, enabled_default, summary, manifest_path, resolve_ref,
performance_contract, gpu_readiness, gamma_contaminating, denylist_bond_types,
geometry_only, baseline_required, requires, citations, pi_decisions,
sanity_gate_ref, notes`.

- **status** ∈ `CORE` (always assembled) · `LIVE` (wired into `Cell.build`,
  togglable) · `GEOMETRY` (frames/broad-phase only, force-free) · `EXPERIMENTAL`
  (module exists, default-OFF, enabling raises) · `STUB` (resolver/contract
  authored, mechanism behind a PI decision).
- **manifest_path** — the dotted key (`("optional_subsystems","fa")`) whose
  `enabled` flag a recipe toggles. `None` for core/geometry/experimental.
- **gamma_contaminating + denylist_bond_types** — whether the compartment's bonds
  must be excluded from the cortical-γ method-of-planes, and the bond-type
  prefixes to denylist. A test cross-checks LIVE ones against
  `cortex/cortical_tension.py::ADHESION_BOND_TYPES`.

### Recipe → manifest composition

`CompartmentRegistry.compose_manifest(recipe, base_manifest, strict=True)` clones
`mcf7_baseline.yaml` and flips LIVE `enabled` flags to match the recipe's `enable`
set. It **refuses** to (a) drop a baseline-required compartment
(`BaselineDropError`) or (b) wire an `EXPERIMENTAL`/`STUB` compartment
(`UnratifiedCompartmentError`) — the same posture as the loader raising on a
PI-gated optional. Recipes list intended-but-unwired compartments under
`declare_pending` (documentation), so a recipe is the one place the whole intended
cell is listed without silently enabling un-ratified physics.

## 3. Compartment catalogue (23)

**CORE (3):** cortex, crosslinkers, myosin.
**Baseline (4, default-ON):** cytoplasm, enclosed_volume (turgor), nucleus,
membrane_surface.
**LIVE optional (7):** fa, rigid_ligand_coating, substrate, lamellipodium,
membrane_load, erm, turnover.
**Geometry (1):** surface_manifold (force-free; ONE soft normal confinement is the
only force it may ever add, default-OFF, PI-gated by a Magic-Number Block).
**Missing / EXPERIMENTAL (8):** ventral_stress_fibers, linc,
intermediate_filaments, microtubules, osmotic_regulation, membrane_reservoir,
cadherin_junction, junctional_actin (STUB).

## 4. Hard rules the platform enforces

1. **Default-OFF preserves behaviour.** Every non-core/non-baseline compartment is
   `enabled_default=False`; with defaults OFF the build is bit-identical
   (`validate_defaults`). The disabled path of every new module must early-return.
2. **No mesh-as-physics.** `geometry_only` compartments declare zero bonds/particles
   (`validate_manifold_geometry_only`). The manifold supplies geometry, frames,
   contact patches, broad-phase, diagnostics — never edge/area/bending springs,
   never a force in the γ budget (`H7_SURFACE_MANIFOLD_EXPLICIT_CORTEX`). Bright
   line: *no manifold DOF in the γ budget.*
3. **No γ contamination.** Adhesion/junction/stress-fiber/LINC/IF/MT/bleb bonds are
   **separate diagnostics**; their prefixes are denylisted from the cortical-γ
   estimator. A recipe may define another observable (traction, junction_tension)
   — that's explicit, not contamination. The active-γ number is
   `γ(myosin ON) − γ(myosin OFF)`; turgor + membrane tension are **separate passive
   channels**, never folded in (`H7_GATE_B_CONTRACT`).
4. **No empirical magic numbers.** Every parameter needs a literature citation or a
   grid-invariant derivation. Unknown → leave the field `None`, add to
   `pi_decisions`, keep the compartment disabled. Authoring agents **must not**
   guess to make a number land.
5. **No gate loosening / no frozen-file edits.** Immutable briefs, validation gates,
   and `integrator/` (BAOAB / M-SHAKE) are read-only without PI sign-off. The
   M-SHAKE compression relaxation for Gate-B is a separate PI-signed change.
6. **Physiological baseline is never silently dropped.** Production keeps cytoplasm,
   turgor, membrane_surface, nucleus ON at setpoint (only the sanctioned
   `allow_no_nucleus` suspended cortical-tension exception).

## 5. Compartment Performance Contract (mandatory, per compartment)

Each `CompartmentSpec.performance_contract` carries: particle types added; particle
count at `n_fil=1000` and native scale; bond/angle count; per-step force? per-batch
updater? uses `cpu_local_snapshot`? uses global `cKDTree`/broad-phase? expected ON
recipes; **hot-path priority P0/P1/P2/P3**; GPU path now (CPU/CuPy/native/builtin);
native ForceCompute candidate?; bottleneck risk.

**Hot-path priority:**
- **P0** — always-on production per-step force/updater. MUST have a GPU-resident or
  native path planned before merge. (cortex, cytoplasm, enclosed_volume, nucleus,
  membrane_surface.)
- **P1** — recipe-on per-step force/updater likely in production; CPU allowed only
  with explicit optimization debt + microbench. (fa, rigid_ligand_coating,
  substrate, lamellipodium, membrane_load, erm.)
- **P2** — batch updater / binding / broad-phase; avoid per-step
  `cpu_local_snapshot`; benchmark if CPU retained. (crosslinkers, myosin, turnover,
  and all experimental compartments.)
- **P3** — one-time layout / diagnostics / visualization; CPU/NumPy fine.
  (surface_manifold geometry.)

**GPU env-flag conventions** (default-OFF, opt-in, parity-validated):
`FFN_GPU_DEVICE_BAOAB=1` (device-resident BAOAB), `FFN_GPU_DEVICE_COMPARTMENTS=1`
(native ForceCompute for turgor/membrane/nucleus, bit-exact vs CPU). Per-step
`cpu_local_snapshot` forces a device→host→device sync that throttles the C++ force
eval — the standing P0/P3 GPU-main debt (`GPU_MAIN_PORT`, `NATIVE_HOT_LOOP`).

## 6. Implementation phases

1. **Registry/recipe spine.** ✅ LANDED (this commit).
2. **Near-term physiology** — ERM, turnover, FA rigid-ligand + r0-overload
   equilibration, lamellipodium + membrane_load. *Already LIVE*; wired into the
   registry as adapters (no code change to the modules).
3. **Adherent architecture** — ventral stress fibers; rigid ligand coating / contact
   layer. (stress fibers: new EXPERIMENTAL module.)
4. **Internal mechanics** — LINC, keratin/intermediate filaments, microtubules.
   (new EXPERIMENTAL modules.)
5. **Volume/membrane remodelling** — dynamic osmotic regulation, membrane
   reservoir/bleb. (new EXPERIMENTAL modules.)
6. **Multicell** — explicit cadherin trans catch-bond junction + junctional actin
   coupling (adapts the existing `validation/cadherin_sliding_rebinding.py` +
   `spheroid/cadherin_bonds.py`). No slip-only shortcut, no lumped line tension.
7. **GPU/native audit** — `docs/v2_audit/COMPARTMENT_GPU_AUDIT_2026-06-08.md` +
   `scripts/compartment_force_profile.py`.

## 7. Per-compartment deliverables (every new module)

module file(s) · manifest/recipe wiring (registry spec) · dataclass resolver (SI
units) · Sanity Gate docstring (or sibling `*_sanity.md`) · Compartment Performance
Contract · tests (OFF identity, dimensional sanity, boundary cases,
conservation/no-net-force where relevant, sign-sense, no-γ-contamination where
relevant) · minimal diagnostic/figure hook · GPU/native readiness note (+ microbench
hook if P0/P1) · parity test if a CuPy/native path exists · explicit optimization
debt note if CPU-only.

## 8. Stop conditions

- Parameter lacks literature/derivation → document PI decision, keep disabled.
- Compartment needs a validation-gate change → stop and document.
- Compartment needs a frozen-integrator edit → stop and document.
- A disabled compartment changes existing behaviour/γ → fix before continuing.
- A design adds force-bearing mesh edges/area springs → reject and document why.

## 9. Activation order (no compartment is enabled in production by this work)

`full_physiological` is a **declaration only** — it never runs globally.
Compartments graduate `declare_pending → enable` one pairwise gate at a time,
PI-signed-off, in the order: (FA+rigid-ligand → lamellipodium+membrane_load are
already LIVE behind PI gates) → ventral_stress_fibers → osmotic_regulation →
LINC → intermediate_filaments → microtubules → membrane_reservoir →
cadherin_junction → junctional_actin. No biological-closure claim until pairwise and
recipe-level gates pass.

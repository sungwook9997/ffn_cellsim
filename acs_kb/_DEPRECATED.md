# acs_kb — v1 frozen as validation reference

As of 2026-05-19, `acs_kb/` is **frozen** as the v1 validation oracle for the
HOOMD-blue full-fidelity rewrite.

- v1 (numpy + single-chain cortex) is preserved here in entirety: rate laws,
  KU sanity gates, parameter dictionaries, regression tests, and notebooks
  all continue to run and remain the ground truth for v2.
- v2 work begins in a new top-level package, `acs_hoomd/`. Do not add new
  features to `acs_kb/`; bug fixes that affect a KU contract should land
  there and be mirrored into `acs_hoomd/` once it exists.
- The audit that classified every module as **reuse / port / archive** lives
  at `acs_kb/outputs/v2_audit/CODEBASE_AUDIT.md`. Start there before
  touching any file in this package.

**Reference**: Notion *🚀 Simulation Build Plan v2 — HOOMD-blue Full-Fidelity*
— <https://www.notion.so/365120daec5d81799efefcf078f2039e>

## Import policy

**No `acs_hoomd/<runtime>/` module may import any file listed under
"Archived" below as a mechanism.** v1 archived files stay in tree as
**closed-form validation oracles** — `acs_hoomd/tests/` and
`acs_hoomd/validation/` may invoke them to cross-check that
v2-emergent rates / behaviours match the literature closed-forms — but
they must not appear in the v2 runtime call graph.

If v2 needs equivalent runtime capability, build it from scratch in
`acs_hoomd/` against the HOOMD-blue API. Carrying v1 shapes /
interfaces / paper-model-as-mechanism wrappers forward warps the v2
design (the entire reason the file is archived in the first place).

**v2 architectural principle (PI 2026-05-19)**: v1 borrowed literature
models (Chan-Odde 2008, Pereverzev 2005, Bell-Evans, Buckley 2014, Hill,
...) as core mechanisms. v2 inverts this: the mechanism is fine-grained
HOOMD particle/bond dynamics, and the literature closed-forms are
acceptance oracles. Anything that gets in the way of fine-grained fiber
network implementation, anything that wraps a paper model as the runtime
mechanism, is archived.

## Archived (runtime import forbidden in v2)

PI decisions 2026-05-19 (two passes):

1. *First pass*: anything that could warp the v2 fiber-network design.
2. *Second pass*: anything that wraps a literature model as the runtime
   mechanism. v2 uses literature models only as validation oracles.

Final list: **18 files / ~3,700 LOC (≈ 45 % of `acs_kb/`)**.

### Cell-level (single-chain cortex contamination)

- `acs_kb/cell/cortex.py` (450 LOC) — single-chain circular cortex on
  `(n_fibers, n_beads, 2)`; v2 uses ~1000 effective filaments per cell
  as HOOMD bonded particles, the rigid shape is incompatible.
- `acs_kb/cell/force_balance.py` (306 LOC) — overdamped Euler over
  cortex beads; replaced by `hoomd.md.methods.Brownian`.
- `acs_kb/cell/cell.py` (123 LOC) — `Cell` dataclass owns a v1 `Cortex`;
  week-6 freeze interface is a v1 contract.
- `acs_kb/cell/lamellipodia.py` (271 LOC) — built on v1 single-chain
  cortex; v2 = AFINES dendritic Arp2/3 (Plan Unit H.5).
- `acs_kb/cell/visualization.py` (200 LOC) — single-chain plot helpers;
  v2 viz = HOOMD GSD → freud / fresnel / PyVista.

### ECM integration / shear / v1 viz

- `acs_kb/ecm/integrator.py` (136 LOC) — Euler-Maruyama; → HOOMD built-in
  overdamped integrators. CFL discipline (`α = 0.1 · τ_min`) ports to
  `acs_hoomd/setup.py` dt picker.
- `acs_kb/ecm/shear_lees_edwards.py` (239 LOC) — hand-rolled Lees-Edwards
  MI + duplicated WLC+XL kernel; → HOOMD-native sheared triclinic box.
- `acs_kb/ecm/shear_protocol.py` (130 LOC) — v2 rewrites strain ramp/hold
  schedule from scratch against HOOMD `BoxResize`.
- `acs_kb/ecm/visualization.py` (144 LOC) — v1 numpy-state-only matplotlib;
  v2 viz consumes HOOMD GSD.

### Bridge (paper-models-as-mechanism — Chan-Odde / Pereverzev / Bell-Evans / Hill)

- `acs_kb/bridge/motor_clutch.py` (364 LOC) — Chan-Odde 2008 quasi-static
  clutch mechanism wrapper. v2 clutch emerges from HOOMD integrin↔ligand
  dynamic bonds with Bell-Evans rate-driven break events. Closed-form
  oracle for biphasic verdict.
- `acs_kb/bridge/catch_bond.py` (136 LOC) — Pereverzev 2005 two-pathway
  off-rate wrapper module. The Pereverzev *formula* is reused in
  `acs_hoomd/validation/pereverzev.py` as oracle, but the wrapper module
  is not v2 mechanism.
- `acs_kb/bridge/talin.py` (119 LOC) — Bell-Evans per-domain unfold
  wrapper. v2 talin unfolding emerges from HOOMD-level domain particles.
- `acs_kb/bridge/vinculin.py` (120 LOC) — vinculin allostery ODE
  `k_int^eff = k_int^bare · (1 + α·N_vin)` wrapper. v2 recruitment
  emerges from HOOMD-level dynamics.
- `acs_kb/bridge/fa_growth.py` (155 LOC) — Hill-function `n_clutches_total`
  resize wrapper. v2 FA growth emerges from HOOMD clutch population
  statistics.
- `acs_kb/bridge/types.py` (78 LOC) — `FocalAdhesion` dataclass; v1
  week-5 freeze interface. v2 FA = HOOMD bond group + metadata dict.
- `acs_kb/bridge/substrate_stub.py` (152 LOC) — linear-elastic Boussinesq
  stub; v1 isolation fixture. v2 substrate = HOOMD ECM directly.
- `acs_kb/bridge/ecm_adapter.py` (209 LOC) — v1 numerical-Hessian probe.
  v2 reads HOOMD ECM particle forces directly at FA capture radius
  (Plan H.4, `R_FA = 1.5 μm`).

### Junction (paper-models-as-mechanism — Buckley 2014)

- `acs_kb/junction/cadherin.py` (225 LOC) — Buckley 2014 Bell-Evans
  slip-only as binomial 2-compartment population updater. v2 cadherin
  bonds emerge from HOOMD cortex↔cortex trans-bonds with Bell-Evans
  rate-driven break events. KU-4.17 constants migrate to v2 config;
  `k_off(F) = k_off^0 · exp(F·Δx*/kT)` formula → `acs_hoomd/validation/buckley.py`.
- `acs_kb/junction/types.py` (142 LOC) — `EcadherinJunction` dataclass;
  v1 interface freeze.

## What survives as v2 runtime (NOT archived)

- `acs_kb/ecm/fiber_network.py` — Mikado geometry generator (pure
  construction, called once at v2 init to seed HOOMD topology).
- `acs_kb/ecm/cross_links.py` — segment-intersection geometry for
  cross-link seeding; harmonic kernel kept as oracle.
- `acs_kb/ecm/fiber_mechanics.py` — KU-1.24 WLC discrete `H` *formula*
  maps to HOOMD `md.bond.Harmonic` + `md.angle.Harmonic`; numpy kernel
  = closed-form oracle.
- `acs_kb/ecm/diagnostics.py` — measurement utility, integrator-agnostic.
- `acs_kb/junction/contact_angle.py` — Maître KU-4.4 geometric measurement.
- `acs_kb/bridge/traction.py` — 1-D analysis reducer.
- `acs_kb/common/{sanity_gate, derived_params, derived_params_cell}.py`
  — KU gates + analytical derived parameters (Stokes drag, Mikado ℓ_c,
  τ_min). Mechanism-free contracts.
- `acs_kb/configs/*.yaml` — KU-anchored literature **constants** (port
  verbatim with two v2 overrides: `ℓ₀: 0.5e-6` and
  `dynamics.integrator: hoomd_brownian`).

## v2 work

v2 lives in a new top-level `acs_hoomd/` package. See
`acs_kb/outputs/v2_audit/CODEBASE_AUDIT.md` for the full
REUSE / ORACLE / ARCHIVE breakdown.

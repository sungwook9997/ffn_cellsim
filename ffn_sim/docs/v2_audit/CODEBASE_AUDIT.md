# acs_kb Codebase Audit — v1 → v2 (HOOMD-blue) migration

**Branch**: `ffn/foundation` (renamed from `v2/foundation` 2026-05-20; cut from `codex/recover-units` @ `c7e3276`)
**Date**: 2026-05-19
**Scope**: classify every `.py` / `.yaml` under `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /` as **reuse**, **port**, or **archive** for the v2 HOOMD-blue rewrite.
**Reference**: Notion *🚀 Simulation Build Plan v2 — HOOMD-blue Full-Fidelity* §2 — <https://www.notion.so/365120daec5d81799efefcf078f2039e>
**Migration target**: v2 lives in a new `ffn_sim/` package; `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /` is frozen as the validation oracle.

> Hardware note: M1 Max, conda 26.1.1 on `osx-arm64`, HOOMD CPU mode through Phase 2 entry. GPU not required for the foundation work.

---

## §1 — File count and LOC

### Headline numbers

| Metric | Value |
| --- | --- |
| Tracked source files (`.py` + `.yaml`) | **54** |
| `.py` modules (incl. tests + notebooks) | **49** |
| `.yaml` configs | **5** |
| Test files | **10** (+ `tests/__init__.py`) |
| Output directories (snapshots) | **11** |
| `import numpy` / `from numpy` occurrences | **39** lines across the package |
| HOOMD-blue installed | **No** (system py 3.9.6 + `.venv-collab` both missing `hoomd`) |

### LOC per module (rolled up by package)

| Module group | Files | LOC | Notes |
| --- | --- | --- | --- |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /ecm/` | 9 | **1,317** | KU-1.x fiber-network + cross-link + integrator |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/` | 10 | **1,158** | KU-2.x catch-bond, motor-clutch, talin, vinculin, FA growth |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /cell/` | 5 | **1,007** | KU-3.x cortex, force-balance, lamellipodia, visualization |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /junction/` | 3 | **579** | KU-4.x cadherin slip-only + contact angle |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /common/` | 3 | **839** | Sanity gate (408), derived params (329), cell-side derived (102) |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /configs/` | 5 | (YAML, 31 KB total) | KU-anchored parameter dictionaries |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /tests/` | 10 | **2,294** | Validation regression suite |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /notebooks/` | 4 | **909** | Demonstration scripts |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /__init__.py` | 1 | 11 | |
| **Total `.py` LOC** | **48 files** | **~8,100** | (per `ffn_sim/docs/v2_audit/loc.txt`) |

### Top-10 files by LOC (where the migration weight lives)

| LOC | File | v2 verdict |
| --- | --- | --- |
| 450 | `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /cell/cortex.py` | **ARCHIVE** — single-chain circular cortex; v2 = 1000 effective filaments as HOOMD bond/angle topology |
| 440 | `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /tests/test_bridge_unit2_2.py` | **ORACLE** — closed-form bridge tests, v1 frozen path |
| 408 | `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /common/sanity_gate.py` | **REUSE** as v2 runtime — the KU gate logic is the regression contract |
| 392 | `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /tests/test_cell.py` | **ORACLE** — v1 frozen path, acceptance bands port to `ffn_sim/tests/` |
| 364 | `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/motor_clutch.py` | **ARCHIVE** — Chan-Odde-as-mechanism is v1 paradigm; v2 clutch emerges from HOOMD integrin↔ligand bonds. Closed-form kept as oracle. |
| 329 | `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /common/derived_params.py` | **REUSE** as v2 runtime — Stokes drag, Mikado ℓ_c, τ_min — analytical, integrator-agnostic |
| 306 | `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /cell/force_balance.py` | **ARCHIVE** — overdamped Euler → `md.methods.Brownian` |
| 302 | `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /notebooks/02_motor_clutch_biphasic.py` | **ORACLE** as v1 demo (biphasic verdict) |
| 301 | `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /tests/test_bridge.py` | **ORACLE** — KU-2.x rate-law closed-form tests |
| 271 | `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /cell/lamellipodia.py` | **ARCHIVE** — single-chain context; v2 = AFINES dendritic Arp2/3 (Plan H.5) |

---

## §2 — Reusable modules (with rationale)

### 2.1 `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /ecm/` — Worker A KU-1.30 validation surface  →  REUSE as oracle, PORT integrator

| File | LOC | Verdict | Rationale |
| --- | --- | --- | --- |
| `fiber_network.py` | 200 | **REUSE** | Mikado generation is geometric (line intersection); produces the initial bead/bond topology that gets handed to HOOMD's `Snapshot.bonds` once. No integrator coupling. |
| `fiber_mechanics.py` | 221 | **REUSE as oracle** | Discrete-WLC stretch + bend `H = (μ/2ℓ₀)Σ(\|b\|−ℓ₀)² + (κ/ℓ₀)Σ(1−cosθ)` maps directly to HOOMD `md.bond.Harmonic` + `md.angle.Harmonic`. Keep this module to *cross-validate* HOOMD energies on a stored configuration; ℓ_p, μ, κ parameters port verbatim. |
| `cross_links.py` | 247 | **REUSE as oracle** | Mikado segment-intersection (KU-1.27/1.3/1.28) construction and the harmonic-spring force kernel both map to a HOOMD bond group. The numpy kernel becomes the unit-test reference; HOOMD bonds carry the runtime force. |
| `diagnostics.py` | (small) | **REUSE** | `measure_coordination` and friends operate on positions/bond-lists — integrator-agnostic. |
| `visualization.py` | 144 | **REUSE** | Plots derived diagnostics; orthogonal to v2 dynamics. |
| `integrator.py` | 136 | **ARCHIVE** | See §3. |
| `shear_lees_edwards.py` | 239 | **ARCHIVE** | See §3 (HOOMD ships Lees-Edwards natively). |
| `shear_protocol.py` | 130 | **ARCHIVE** | PI decision 2026-05-19: v2 rewrites against HOOMD `BoxResize` from scratch — porting the v1 caller would re-import the archived sheared kernel. |

**Why preserve the numpy kernels as an oracle.** The KU-1.30 stretching/bending Sanity Gate evidence (`max relative gradient error 4.2e-9`, `Σ F ≈ 1e-25 N`) lives in these modules; freezing them gives v2 a per-configuration F and E reference that HOOMD must reproduce within `1e-6` relative before any new physics is added on top.

### 2.2 `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/` — paper-models-as-mechanism (v1 paradigm)  →  ARCHIVE, oracle only

PI 2nd pass: v1's bridge layer wraps Chan-Odde / Pereverzev / Bell-Evans / Hill as the *mechanism* that drives the simulation. v2 builds clutch / catch-bond / talin / vinculin / FA-growth behaviour emergently from HOOMD-level particle dynamics and uses the literature closed-forms as **acceptance oracles**, not as runtime. Details in §3.3. Survives as v2 runtime: only `traction.py` (53 LOC analysis reducer) and `__init__.py`.

### 2.3 `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /junction/cadherin.py` — paper-model-as-mechanism (v1 paradigm)  →  ARCHIVE, oracle only

Buckley 2014 Bell-Evans slip-only k_off and the binomial 2-compartment population updater are v1 mechanism. v2: cadherin bonds emerge from HOOMD trans-bonds between cortex particles on adjacent cells with Bell-Evans rate-driven break events. The KU-4.17 constants (Δx* = 4 nm, k_off^0 = 0.5 s⁻¹, k_on = 1 s⁻¹, N_cad = 100, ⟨F_bond⟩ = 30 pN, T = 310 K) move into the v2 config; the closed-form `k_off(F) = k_off^0 · exp(F·Δx*/kT)` survives in `ffn_sim/validation/buckley.py` as the acceptance oracle. Survives as v2 runtime: `contact_angle.py` (Maître geometry) + `__init__.py`.

### 2.4 `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /configs/*.yaml` — KU-anchored parameter dictionaries  →  REUSE (with v2 overrides)

All 5 YAMLs (`phase1_unit1.yaml` … `phase1_unit4_1.yaml`, ~31 KB) cite the source KU on every line. The parameter inventory is the largest single piece of audit-able prior art.

| Value class | Reuse-as-is in v2 |
| --- | --- |
| Biological scales — `L_fiber`, `persistence_length`, `bending_modulus`, `stretching_modulus`, `bead_radius`, `kT`, `water_viscosity` | **Yes**, KU-anchored. |
| Geometry — `L_box`, `dimension: 2` (Phase 1), `biological_mesh ξ` | **Yes** for ECM; cell-side `R_cell` reused. |
| Cross-link stiffness `xl_stiffness: 1e-3 N/m` (KU-1.28) | **Yes**. |
| `beads_per_fiber: 5`, derived `rest_length = L_fiber / (beads − 1) = 2.5 µm` | **OVERRIDE** in v2 (target `ℓ₀ = 0.5 µm` per Build Plan §2 → 21 beads/fiber for the same `L_fiber = 10 µm`). |
| `dynamics.integrator: euler_maruyama` | **OVERRIDE** → HOOMD `Brownian` (or `Langevin`); `cfl_safety_factor 0.1` carries over as the dt-vs-`τ_min` discipline. |
| Acceptance bands (`z_range`, `segment_length_range`, …) | **Yes** — these are the validation contracts. |
| `n_fibers` derived from `ℓ_c = π/(2 ρ_L)` (KU-1.27) | Formula reused; the v2 scale (1000 effective filaments/cell) re-resolves it from the same equation. |

**Migration recipe**: deep-copy each YAML into `ffn_sim/configs/` with `ell0: 0.5e-6` and `integrator: hoomd_brownian` overrides. `derived_params.resolve()` recomputes everything downstream from the new primary scales (KU-1.27 + KU-1.26 + Stokes drag).

### 2.5 `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /common/` — KU gates and derived parameters  →  REUSE wholesale

- `sanity_gate.py` (408 LOC) is the regression contract. Every KU gate (`gate_phase1_cell_cortex`, `gate_phase1_unit3_1`, `gate_phase1_unit4_1`, …) compares **measurement vs KU range**, not integrator vs integrator — so it ports unchanged.
- `derived_params.py` (329 LOC) and `derived_params_cell.py` (102 LOC) implement Stokes drag `γ_b = 6π η_water · bead_radius`, Mikado `ρ_L`, three relaxation times (`τ_xl`, `τ_stretch`, `τ_bend`), CFL bound `dt < α · τ_min`. All analytical, integrator-free, must be carried over verbatim (and *also* the source of HOOMD's recommended `dt`).

### 2.6 Test suite — REUSE as validation oracle

10 test files / 2,294 LOC. Two flavours:

1. **Rate-law tests** that touch no archived module (`test_cadherin.py`, `test_bridge.py`, `test_KU_2_4_biphasic.py`, `test_contact_angle.py`): closed-form rate equations and stochastic updaters; **pass through v2 unchanged**.
2. **Oracle tests** that run on the v1 frozen path (`test_lamellipodia.py`, `test_bridge_unit2_2.py`, `test_dynamics.py`, `test_ecm.py`, `test_cell.py`, `test_two_cell_pair.py`): exercise either the numpy integrator, the v1 `Cell` shape, the archived ECM adapter, or v1 lamellipodia. They stay green as long as the v1 modules sit unchanged in tree — the validation oracle role. **Acceptance bands** port to `ffn_sim/tests/` against new HOOMD-driven fixtures.

---

## §3 — Modules to archive / deprecate (with rationale)

**PI decisions 2026-05-19**:

1. *First pass*: anything that could warp the v2 design is archived. ⇒ 9 files / ~2,064 LOC (cell/* + ecm integrator/shear + bridge/ecm_adapter).
2. *Second pass*: **v1 = paper-models-as-mechanism, v2 = paper-models-as-validation-oracle**. Anything that wraps a literature model as the *core implementation* of clutch / bond / FA growth / cadherin dynamics gets archived, because v2 builds those behaviours emergently from fine-grained HOOMD fiber/particle dynamics and only *compares* the emergent rates to the literature closed-forms. Anything that gets in the way of fine-grained fiber network implementation is archived.

Final scope: **~18 files / ~3,700 LOC (≈ 45 % of `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /`)** archived.

> v1 archived files stay in tree as **validation oracles** — v2 acceptance tests may invoke them to compare emergent HOOMD behaviour against the closed-form references. They are **forbidden as v2 runtime mechanism** (no `ffn_sim/<runtime>/*.py` may import them). The import rule lives in `ffn_sim/docs/v2_audit/_DEPRECATED.md`.

### 3.1 Cell-level (single-chain cortex contamination)

| File | LOC | v2 replacement / rationale |
| --- | --- | --- |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /cell/cortex.py` | 450 | v2 uses ~1000 effective filaments per cell as HOOMD bond/angle topology. The rigid `(n_fibers, n_beads, 2)` layout cannot represent that without distortion. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /cell/force_balance.py` | 306 | `md.methods.Brownian`. KU-3.5 cortical tension `γ_cortex / R_eff · (−r̂_eff)` becomes `md.force.Custom`. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /cell/cell.py` | 123 | `Cell` dataclass owns a v1 `Cortex` — week-6 freeze interface is a **v1 contract** that warps every junction / bridge consumer toward the single-chain shape. v2 `Cell` lives in `ffn_sim/cell/` and is HOOMD-particle-group-backed. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /cell/lamellipodia.py` | 271 | Built on the single-chain cortex; v2 rebuilds as AFINES dendritic Arp2/3 branched network (Plan Unit H.5, Bieling 2016 force feedback). The KU-3.6 closed-form `v_p = δ(k_on c_G − k_off · e^{Fδ/kT})` ports as a literal formula but the *module* does not. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /cell/visualization.py` | 200 | Single-chain plot helpers; v2 viz consumes HOOMD GSD via freud / fresnel / PyVista. |

### 3.2 ECM integration / shear

| File | LOC | v2 replacement / rationale |
| --- | --- | --- |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /ecm/integrator.py` | 136 | Euler-Maruyama + dead BAOAB stub. HOOMD ships overdamped Brownian + Langevin natively. CFL discipline `α = 0.1 · τ_min` moves to `ffn_sim/setup.py` dt picker. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /ecm/shear_lees_edwards.py` | 239 | Hand-rolled sheared MI **plus a duplicated WLC + XL kernel** (docstring acknowledges "intentional and local duplication"). HOOMD-native sheared triclinic box / `BoxResize` makes both redundant. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /ecm/shear_protocol.py` | 130 | **PI decision 2026-05-19**: v2 rewrites the strain ramp/hold schedule from scratch against the HOOMD box updater rather than porting the v1 caller. The v1 caller is glued to the archived sheared kernel and carrying it forward would re-introduce the v1 shape into v2. |

### 3.3 Bridge — paper-models-as-mechanism (v1 paradigm)

In v1, `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/` *is* a port of the Chan-Odde 2008 / Pereverzev 2005 / Bell-Evans / Hill-function papers, with each module implementing one paper's equations as the mechanism that drives the simulation. v2 inverts this: the mechanism is fine-grained HOOMD particle-and-bond dynamics (integrin particles, ligand particles, dynamic bonds with Bell-Evans rate updates), and the literature closed-forms become **acceptance oracles** (e.g. "our HOOMD-emergent clutch lifetime vs Pereverzev `τ(F) = 1 / k_off(F)` at the catch peak"). Therefore the entire `bridge/` mechanism layer is archived; the closed-form *functions* are kept in tree so v2 validation tests can call them.

| File | LOC | v2 replacement / rationale |
| --- | --- | --- |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/motor_clutch.py` | 364 | Chan-Odde 2008 quasi-static substrate force balance + per-clutch stochastic engage/disengage as the **core clutch mechanism**. v2: clutch behaviour emerges from HOOMD integrin↔ligand dynamic bonds with Bell-Evans rate-driven break events. v1 module becomes a closed-form *acceptance oracle* for the biphasic verdict, not the runtime. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/catch_bond.py` | 136 | Pereverzev 2005 two-pathway off-rate as a wrapper module. v2: the same Pereverzev formula stays in `ffn_sim/validation/pereverzev.py` (oracle) so the HOOMD-emergent k_off can be compared to it. The wrapper module is not v2 mechanism. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/talin.py` | 119 | Bell-Evans per-domain unfold as wrapper. v2: talin unfolding emerges from HOOMD-level domain particles with Bell-Evans bond breaking. Oracle reusable. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/vinculin.py` | 120 | Allostery ODE `k_int^eff = k_int^bare · (1 + α·N_vin)` as wrapper. v2: vinculin recruitment emerges from HOOMD-level dynamics; the closed-form k_int^eff(N_vin) becomes the oracle. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/fa_growth.py` | 155 | Hill-function `n_clutches_total` resize as wrapper. v2: FA growth emerges from HOOMD clutch population statistics; the Hill function becomes the oracle. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/types.py` (FocalAdhesion dataclass) | 78 | Week-5 freeze interface, v1 contract. v2 FA is a HOOMD bond group + per-FA metadata dictionary, not a dataclass. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/substrate_stub.py` | 152 | Linear-elastic Boussinesq stub — v1 fixture for testing motor-clutch in isolation. v2 substrate is the HOOMD ECM directly; no stub layer needed. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/ecm_adapter.py` | 209 | v1 numerical-Hessian probe. v2 reads HOOMD ECM particle forces directly at the FA capture radius (Plan H.4, `R_FA = 1.5 μm`). The `SubstrateProtocol` abstraction in `motor_clutch.py` is not carried forward (motor_clutch itself is archived). |

Survives in `bridge/`: only `traction.py` (53 LOC, 1-D reducer — pure analysis utility) and `__init__.py`.

### 3.4 Junction — paper-models-as-mechanism (v1 paradigm)

| File | LOC | v2 replacement / rationale |
| --- | --- | --- |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /junction/cadherin.py` | 225 | Buckley 2014 Bell-Evans slip-only as a binomial 2-compartment population updater — v1 mechanism. v2: cadherin bonds emerge from HOOMD trans-bonds between cortex particles on adjacent cells with Bell-Evans rate-driven break events. The Buckley `k_off(F) = k_off^0 · exp(F·Δx*/kT)` closed-form stays as oracle in `ffn_sim/validation/buckley.py`. Δx* = 4 nm constant is in the v2 config. |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /junction/types.py` (EcadherinJunction dataclass) | 142 | v1 interface freeze; v2 junctions are HOOMD bond groups. |

Survives in `junction/`: only `contact_angle.py` (212 LOC, pure Maître angle geometry — measurement utility, no v1 mechanism) and `__init__.py`. The `Cell` import in `contact_angle.py` rewires to `ffn_sim/cell/`.

### 3.5 ECM — visualization

| File | LOC | v2 replacement / rationale |
| --- | --- | --- |
| `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /ecm/visualization.py` | 144 | v1 numpy-state-only plot helpers. v2 viz consumes HOOMD GSD via freud/fresnel/PyVista (Plan §8 cross-cutting). Geometry/diagnostics module-output plots stay (`ecm/diagnostics.py`) but the matplotlib v1 viz layer is archived. |

### 3.6 What is NOT archived (REUSE)

- `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /ecm/fiber_network.py` (200 LOC) — Mikado geometry generator. Pure construction, integrator-free; used **once** at v2 init to populate HOOMD topology.
- `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /ecm/cross_links.py` (247 LOC) — segment-intersection geometry for cross-link seeding; pure construction. Harmonic-spring kernel kept as oracle for `md.bond.Harmonic` cross-validation.
- `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /ecm/fiber_mechanics.py` (221 LOC) — KU-1.24 WLC discrete `H` formula. The **formula** is what v2 uses (mapped to HOOMD `md.bond.Harmonic` + `md.angle.Harmonic`); the numpy kernel is the closed-form oracle (max FD error 4.2e-9 reference).
- `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /ecm/diagnostics.py` (~80 LOC) — measurement utility on positions/bond-lists. Integrator-agnostic.
- `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /junction/contact_angle.py` (212 LOC) — Maître KU-4.4 geometric measurement.
- `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /bridge/traction.py` (53 LOC) — 1-D analysis reducer.
- `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /common/` (all 3 files, 839 LOC) — Sanity gates + analytical derived parameters. Mechanism-free, contract-only.
- `acs_kb/ (deleted in 2026-05-20 rename to ffn_cellsim) /configs/*.yaml` (5 files, ~31 KB) — KU-anchored literature constants. **Constants, not models** — port verbatim with two v2 overrides (`ℓ₀: 0.5e-6`, `dynamics.integrator: hoomd_brownian`).
- Tests — split below.

### 3.7 Tests — split by what they exercise

- **v2-relevant** (run against `ffn_sim/`): `test_contact_angle.py` (163 LOC, Maître geometry) and the oracle-comparison tests we'll write in `ffn_sim/tests/validation/`.
- **v1 frozen oracle** (continue to run against v1 path in tree, used as closed-form reference for v2 acceptance): `test_KU_2_4_biphasic.py`, `test_bridge.py`, `test_bridge_unit2_2.py`, `test_cadherin.py`, `test_cell.py`, `test_dynamics.py`, `test_ecm.py`, `test_lamellipodia.py`, `test_two_cell_pair.py`.
- v2 ships parallel HOOMD-driven tests in `ffn_sim/tests/`.
- Notebooks (909 LOC across 4 files) are v1 demos — historical record, no porting.

---

## §4 — HOOMD-blue install status

`ffn_sim/docs/v2_audit/hoomd_check.txt` snapshot:

```
--- system python3 ---
/usr/bin/python3
Python 3.9.6

--- HOOMD on system python3 ---
ModuleNotFoundError: No module named 'hoomd'

--- HOOMD on .venv-collab ---
ModuleNotFoundError: No module named 'hoomd'

--- conda check ---
conda
conda 26.1.1

--- arch ---
arm64
26.4.1
```

**Verdict**: HOOMD not installed in any reachable interpreter. conda 26.1.1 is on PATH on macOS 26.4.1 / Apple silicon (`arm64`).

**Install plan for Phase 0.2**:

```bash
conda create -n ffn_sim python=3.12 -y
conda activate ffn_sim
conda install -c conda-forge "hoomd>=4.7=*cpu*" -y   # osx-arm64 ships a CPU-only build
python -c "import hoomd; print(hoomd.version.version)"
```

System Python 3.9.6 is too old for current HOOMD releases (HOOMD 4 requires Python ≥ 3.9, but the conda-forge `osx-arm64` builds target 3.10+ — 3.12 is the safe default). `.venv-collab` is shared with the viz sprint and should not be repurposed for HOOMD; keep environments segregated.

**Hardware note**: HOOMD on `osx-arm64` runs CPU-only (no GPU backend on Apple silicon). The Build Plan §2 explicitly accepts CPU mode through Phase 2 entry; GPU migration (Linux + CUDA) is a Phase 2+ concern and out of scope for the v2 foundation.

---

## Summary verdict (revised 2026-05-19, second PI pass)

| Bucket | LOC | Files | Headline content |
| --- | --- | --- | --- |
| **REUSE as v2 runtime** (geometry / constants / gates / measurement) | ~2,500 | 13 | `ecm/{fiber_network, cross_links geometry, fiber_mechanics formula, diagnostics}`, `common/*`, `configs/*.yaml`, `bridge/traction`, `junction/contact_angle` |
| **REUSE as validation oracle only** (closed-form references invoked from `ffn_sim/tests/`) | ~1,800 | 9 | `bridge/{motor_clutch, catch_bond, talin, vinculin, fa_growth}` closed-forms, `junction/cadherin` k_off formula, v1 KU regression test bodies |
| **ARCHIVE — `ffn_sim/` runtime import forbidden** | ~3,700 | ~18 | All of `cell/*` (5), `ecm/{integrator, shear_lees_edwards, shear_protocol, visualization}` (4), all of `bridge/*` except `traction.py` (8), `junction/{cadherin, types}` (2) |

≈ **45 %** of v1 LOC is archived from v2 runtime — anything that wraps a literature model as the mechanism (v1 paradigm) is gone, anything that warps fine-grained fiber dynamics is gone. The remaining ~55 % is split between:
- **construction / measurement / gates** (used live in v2 runtime — geometry seeding, KU validation, parameter resolution)
- **closed-form oracles** (kept in tree as a frozen reference library that v2 acceptance tests invoke for comparison)

v1 frozen test suite continues to run unmodified — it is part of the oracle.

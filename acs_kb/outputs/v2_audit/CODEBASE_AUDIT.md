# acs_kb Codebase Audit — v1 → v2 (HOOMD-blue) migration

**Branch**: `v2/foundation` (cut from `codex/recover-units` @ `c7e3276`)
**Date**: 2026-05-19
**Scope**: classify every `.py` / `.yaml` under `acs_kb/` as **reuse**, **port**, or **archive** for the v2 HOOMD-blue rewrite.
**Reference**: Notion *🚀 Simulation Build Plan v2 — HOOMD-blue Full-Fidelity* §2 — <https://www.notion.so/365120daec5d81799efefcf078f2039e>
**Migration target**: v2 lives in a new `acs_hoomd/` package; `acs_kb/` is frozen as the validation oracle.

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
| `acs_kb/ecm/` | 9 | **1,317** | KU-1.x fiber-network + cross-link + integrator |
| `acs_kb/bridge/` | 10 | **1,158** | KU-2.x catch-bond, motor-clutch, talin, vinculin, FA growth |
| `acs_kb/cell/` | 5 | **1,007** | KU-3.x cortex, force-balance, lamellipodia, visualization |
| `acs_kb/junction/` | 3 | **579** | KU-4.x cadherin slip-only + contact angle |
| `acs_kb/common/` | 3 | **839** | Sanity gate (408), derived params (329), cell-side derived (102) |
| `acs_kb/configs/` | 5 | (YAML, 31 KB total) | KU-anchored parameter dictionaries |
| `acs_kb/tests/` | 10 | **2,294** | Validation regression suite |
| `acs_kb/notebooks/` | 4 | **909** | Demonstration scripts |
| `acs_kb/__init__.py` | 1 | 11 | |
| **Total `.py` LOC** | **48 files** | **~8,100** | (per `acs_kb/outputs/v2_audit/loc.txt`) |

### Top-10 files by LOC (where the migration weight lives)

| LOC | File | v2 verdict |
| --- | --- | --- |
| 450 | `acs_kb/cell/cortex.py` | **ARCHIVE** — single-chain circular cortex; v2 replaces with 1000 effective filaments managed as HOOMD bond/angle topology |
| 440 | `acs_kb/tests/test_bridge_unit2_2.py` | **REUSE** as validation oracle |
| 408 | `acs_kb/common/sanity_gate.py` | **REUSE** (the KU gate logic is the regression contract) |
| 392 | `acs_kb/tests/test_cell.py` | **PORT** (re-target dynamics fixtures onto HOOMD outputs) |
| 364 | `acs_kb/bridge/motor_clutch.py` | **PORT** — kinetic core stays, HOOMD-state adapter added |
| 329 | `acs_kb/common/derived_params.py` | **REUSE** (Stokes drag, Mikado ℓ_c, τ_min — analytical, integrator-agnostic) |
| 306 | `acs_kb/cell/force_balance.py` | **ARCHIVE** — overdamped Euler over cortex beads is replaced by HOOMD `md.methods.Brownian` |
| 302 | `acs_kb/notebooks/02_motor_clutch_biphasic.py` | **REUSE** (biphasic verdict demo) |
| 301 | `acs_kb/tests/test_bridge.py` | **REUSE** (KU-2.x rate-law regression) |
| 271 | `acs_kb/cell/lamellipodia.py` | **PORT** — KU-3.6 closed-form keeps, advance hook re-bound to HOOMD update step |

---

## §2 — Reusable modules (with rationale)

### 2.1 `acs_kb/ecm/` — Worker A KU-1.30 validation surface  →  REUSE as oracle, PORT integrator

| File | LOC | Verdict | Rationale |
| --- | --- | --- | --- |
| `fiber_network.py` | 200 | **REUSE** | Mikado generation is geometric (line intersection); produces the initial bead/bond topology that gets handed to HOOMD's `Snapshot.bonds` once. No integrator coupling. |
| `fiber_mechanics.py` | 221 | **REUSE as oracle** | Discrete-WLC stretch + bend `H = (μ/2ℓ₀)Σ(\|b\|−ℓ₀)² + (κ/ℓ₀)Σ(1−cosθ)` maps directly to HOOMD `md.bond.Harmonic` + `md.angle.Harmonic`. Keep this module to *cross-validate* HOOMD energies on a stored configuration; ℓ_p, μ, κ parameters port verbatim. |
| `cross_links.py` | 247 | **REUSE as oracle** | Mikado segment-intersection (KU-1.27/1.3/1.28) construction and the harmonic-spring force kernel both map to a HOOMD bond group. The numpy kernel becomes the unit-test reference; HOOMD bonds carry the runtime force. |
| `diagnostics.py` | (small) | **REUSE** | `measure_coordination` and friends operate on positions/bond-lists — integrator-agnostic. |
| `visualization.py` | 144 | **REUSE** | Plots derived diagnostics; orthogonal to v2 dynamics. |
| `integrator.py` | 136 | **ARCHIVE** | See §3. |
| `shear_lees_edwards.py` | 239 | **ARCHIVE** | See §3 (HOOMD ships Lees-Edwards natively). |
| `shear_protocol.py` | 130 | **PORT** | Strain protocol logic is reusable; only the kernel call swaps to HOOMD's box updater. |

**Why preserve the numpy kernels as an oracle.** The KU-1.30 stretching/bending Sanity Gate evidence (`max relative gradient error 4.2e-9`, `Σ F ≈ 1e-25 N`) lives in these modules; freezing them gives v2 a per-configuration F and E reference that HOOMD must reproduce within `1e-6` relative before any new physics is added on top.

### 2.2 `acs_kb/bridge/` — Worker B catch-bond + biphasic verdict  →  REUSE rate laws, PORT couplers

| File | LOC | Verdict | Rationale |
| --- | --- | --- | --- |
| `catch_bond.py` | 136 | **REUSE** | Pereverzev two-pathway `k_off(F) = k_s e^{F/F_s} + k_c e^{−F/F_c}` is a closed-form scalar; KU-2.5 + KU-2.18 defaults port unchanged. The Gillespie regression test (`test_KU_2_4_biphasic.py`) keeps its scientific weight. |
| `motor_clutch.py` | 364 | **PORT** | Chan-Odde quasi-static substrate balance + per-clutch stochastic engage/disengage runs as a HOOMD `Updater` (custom Action) ticked at the HOOMD timestep boundary. The `SubstrateProtocol` abstraction (stub vs ECM adapter) was *designed* exactly for this swap. |
| `talin.py` | 119 | **REUSE** | Bell-Evans per-domain unfold rate — pure rate law. |
| `vinculin.py` | 120 | **REUSE** | Recruitment ODE; `k_int^eff = k_int^bare · (1+α·N_vin)` allostery — pure rate law. |
| `fa_growth.py` | 155 | **REUSE** | `n_clutches_total` resize ODE; pure state transition. |
| `ecm_adapter.py` | 209 | **PORT** | Replace stub-style `compute_displacement` with a query against HOOMD particle data on the ECM filaments. |
| `substrate_stub.py` | 152 | **ARCHIVE-then-keep-as-fallback** | Linear-elastic spring substrate; useful as a degenerate-case fixture even after the HOOMD ECM is online. |
| `types.py` | 78 | **REUSE** | Frozen dataclass interface (`FocalAdhesion`); week-5 freeze contract. |
| `traction.py` | 53 | **REUSE** | 1-D reducer; integrator-agnostic. |

**Biphasic verdict status**: `test_bridge.py` + `test_bridge_unit2_2.py` (740 LOC combined) directly anchor the Unit 2.x verdict. They must continue to pass against the ported HOOMD-coupled motor-clutch updater before any v2 claim is made on KU-2.4.

### 2.3 `acs_kb/junction/cadherin.py` — Buckley 2014 Δx* correction  →  REUSE

- Implements KU-4.17 / KU-4.2 **slip-only** Bell-Evans `k_off(F) = k_off^0 · exp(F · Δx* / kT)` with `Δx* = 4 nm` cited Buckley 2014 Science Fig 4, `k_off^0 = 0.5 s⁻¹`, `k_on = 1 s⁻¹`, `N_cad = 100`, `⟨F_bond⟩ = 30 pN`, `T = 310 K`.
- `update_bonds` is a binomial two-compartment update on `(n_engaged, n_total − n_engaged)`. **Integrator-agnostic** — runs as a HOOMD `Updater` ticked at `dt_junction` without any HOOMD-internal change.
- `test_cadherin.py` (157 LOC) plus `test_two_cell_pair.py` (184 LOC) anchor KU-4.x; both reused verbatim.
- Catch-bond (full Buckley 2014) is *explicitly* a Phase 2 follow-up; the slip-only baseline is the v2 entry point.

### 2.4 `acs_kb/configs/*.yaml` — KU-anchored parameter dictionaries  →  REUSE (with v2 overrides)

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

**Migration recipe**: deep-copy each YAML into `acs_hoomd/configs/` with `ell0: 0.5e-6` and `integrator: hoomd_brownian` overrides. `derived_params.resolve()` recomputes everything downstream from the new primary scales (KU-1.27 + KU-1.26 + Stokes drag).

### 2.5 `acs_kb/common/` — KU gates and derived parameters  →  REUSE wholesale

- `sanity_gate.py` (408 LOC) is the regression contract. Every KU gate (`gate_phase1_cell_cortex`, `gate_phase1_unit3_1`, `gate_phase1_unit4_1`, …) compares **measurement vs KU range**, not integrator vs integrator — so it ports unchanged.
- `derived_params.py` (329 LOC) and `derived_params_cell.py` (102 LOC) implement Stokes drag `γ_b = 6π η_water · bead_radius`, Mikado `ρ_L`, three relaxation times (`τ_xl`, `τ_stretch`, `τ_bend`), CFL bound `dt < α · τ_min`. All analytical, integrator-free, must be carried over verbatim (and *also* the source of HOOMD's recommended `dt`).

### 2.6 Test suite — REUSE as validation oracle

10 test files / 2,294 LOC. Two flavours:

1. **Rate-law tests** (`test_cadherin.py`, `test_bridge.py`, `test_bridge_unit2_2.py`, `test_KU_2_4_biphasic.py`, `test_contact_angle.py`, `test_lamellipodia.py`): exercise closed-form rate equations or stochastic updaters; **pass through v2 unchanged**.
2. **Dynamics tests** (`test_ecm.py`, `test_dynamics.py`, `test_cell.py`, `test_two_cell_pair.py`): exercise the numpy integrator. **Keep the assertion bands**, swap the fixture to drive HOOMD instead of `acs_kb.ecm.integrator`.

---

## §3 — Modules to archive / deprecate (with rationale)

### 3.1 `acs_kb/cell/cortex.py` (450 LOC) — ARCHIVE

- Implements a **single-chain circular cortex**: 2D fiber centres on a circle of radius `R_cell`, each fiber a discrete WLC bead-chain, cross-links by KU-1.27 intersection.
- v2 conceptually abandons "one cortex = one ring of chained fibers" in favour of **~1000 effective actin filaments per cell** managed as a HOOMD bond/angle topology with the same per-bond Hamiltonian. The construction code (`generate_cortex`, `generate_elliptical_cortex`, `Cortex` dataclass with `(n_fibers, n_beads, 2)` layout) does not generalise to that representation.
- **Disposition**: leave file in tree (frozen) for reproducibility of v1 cortical-tension equilibrium tests, but mark as v1-only via `acs_kb/_DEPRECATED.md`. v2's equivalent lives at `acs_hoomd/cell/cortex_filaments.py` (new).
- Counter-cost of *re-using* this: high — the rigid `(n_fibers, n_beads, 2)` indexing pattern would force the HOOMD snapshot to mirror v1 sizing rather than the 1000-filament target.

### 3.2 `acs_kb/cell/force_balance.py` (306 LOC) — ARCHIVE

- Overdamped quasi-static Euler `γ_b · dr/dt = F_WLC + F_xl + F_tension + F_ext` integrated by hand over cortex beads.
- HOOMD `md.methods.Brownian` (or `Langevin` with appropriate `kT`) replaces the time-stepping entirely. The KU-3.5 `F_tension = γ_cortex / R_eff · (−r̂_eff)` term moves to a HOOMD `md.force.Custom` action so it still ticks each step.
- **Disposition**: archive with cortex.py.

### 3.3 `acs_kb/ecm/integrator.py` (136 LOC) — ARCHIVE

- Euler-Maruyama overdamped Langevin numpy step `r ← r + (F/γ_b) dt + √(2 kT dt/γ_b) ξ` plus a BAOAB stub that has always raised `NotImplementedError`.
- HOOMD ships overdamped Brownian + BAOAB-class integrators built in.
- **Disposition**: archive. The CFL discipline (`α = 0.1`) and `τ_min` reasoning that drove `dt` selection move to `acs_hoomd/setup.py` as the HOOMD-side dt picker.

### 3.4 `acs_kb/ecm/shear_lees_edwards.py` (239 LOC) + `shear_protocol.py` (130 LOC, partial) — ARCHIVE the kernels

- Hand-rolled Lees-Edwards sheared minimum-image MI for a 2D periodic box, including a **duplicated** WLC + cross-link force kernel that uses the sheared MI (the module docstring acknowledges the duplication is "intentional and local").
- HOOMD-blue provides Lees-Edwards box deformation natively (`hoomd.update.BoxResize` + sheared triclinic boxes / `fix deform`-equivalent). The duplicated kernel is redundant.
- **Disposition**: archive `shear_lees_edwards.py` entirely. Keep `shear_protocol.py` as a **port** target — the strain ramp / hold schedule logic is reusable; only the kernel calls swap.

### 3.5 Numpy integrator coupling sites — PORT, not archive

39 lines `import numpy` across the package. After the four archives above (cortex / force_balance / ecm.integrator / shear_lees_edwards = ~1,130 LOC), the residual numpy usage is:

- Rate-law kernels (`bridge/catch_bond.py`, `junction/cadherin.py`, `bridge/motor_clutch.py`, `bridge/talin.py`, `bridge/vinculin.py`) — keep numpy; these stay CPU-side as HOOMD `Updater`s.
- Diagnostics / measurement (`ecm/diagnostics.py`, `cell/visualization.py`) — keep numpy.
- Test fixtures — keep numpy.

Net effect: v2 does **not** banish numpy from `acs_kb/`; it removes only the four time-integration-bearing modules above.

### 3.6 Notebooks

The four `notebooks/*.py` (909 LOC) demonstrate v1 dynamics. **REUSE** as demos against the frozen v1 path. v2 will ship parallel HOOMD notebooks in `acs_hoomd/notebooks/`.

---

## §4 — HOOMD-blue install status

`acs_kb/outputs/v2_audit/hoomd_check.txt` snapshot:

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
conda create -n acs_hoomd python=3.12 -y
conda activate acs_hoomd
conda install -c conda-forge "hoomd>=4.7=*cpu*" -y   # osx-arm64 ships a CPU-only build
python -c "import hoomd; print(hoomd.version.version)"
```

System Python 3.9.6 is too old for current HOOMD releases (HOOMD 4 requires Python ≥ 3.9, but the conda-forge `osx-arm64` builds target 3.10+ — 3.12 is the safe default). `.venv-collab` is shared with the viz sprint and should not be repurposed for HOOMD; keep environments segregated.

**Hardware note**: HOOMD on `osx-arm64` runs CPU-only (no GPU backend on Apple silicon). The Build Plan §2 explicitly accepts CPU mode through Phase 2 entry; GPU migration (Linux + CUDA) is a Phase 2+ concern and out of scope for the v2 foundation.

---

## Summary verdict

| Bucket | LOC | Files | Headline content |
| --- | --- | --- | --- |
| **REUSE** (oracle / rate laws / configs / gates) | ~5,000 | 30+ | All of `bridge/`, `junction/`, `common/`, `configs/`, tests, ecm rate kernels and diagnostics |
| **PORT** (couplers + dynamics-bearing tests) | ~1,200 | 6 | `motor_clutch.py`, `ecm_adapter.py`, `shear_protocol.py`, `lamellipodia.py`, dynamics tests, cell-side glue |
| **ARCHIVE** (replaced by HOOMD) | ~1,130 | 4 | `cell/cortex.py`, `cell/force_balance.py`, `ecm/integrator.py`, `ecm/shear_lees_edwards.py` |

≈ **86 %** of v1 LOC (oracle + rate-law + gate code) ports forward without scientific rework; the ~14 % archived is precisely the numpy-integrator + single-chain-cortex layer that v2 is built to replace.

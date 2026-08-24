# Living-Cortex Engine — Architecture & Build Plan (2026-07-16)

**The new engine, in one sentence.** A living cortex = **fine-grained actomyosin (discrete)** immersed
in a **Biot poroelastic cytosol (continuum fluid)** fed by a **mass-conserving monomer field (continuum
solute)**, all two-way coupled in one coordinate frame. The field layers use a conservative structured
grid; discrete partner/steric search uses a separate device neighbour HashGrid. The monomer economy turns today's frozen
statics truss into a *turning-over, flowing, self-remodeling* material — **without ever coarse-graining a
single filament.**

**PI runtime/baseline contract (2026-07-16):** every simulation step is Warp on CUDA GPU. HOOMD is never
executed; archived source is read-only port input. The production motor is a head-resolved Stam-Hocky
minifilament with Hill FV. The first baseline is MCF7×collagen×α2β1 with 70,686 cortical F-actin filaments;
other cellular filament systems, the nuclear lamina/chromatin network, and collagen ECM fibers/cross-links add
unique populations from their own physiological density×geometry ledgers.

This document is the engine-architecture spec for the CFD-Transport program (`README.md`). It **does not
re-derive** the fluid/dynamic-fiber engineering already scoped in
[`../v2_audit/_historical/DYNAMIC_FEM_CFD_UPGRADE_ROADMAP_2026-07-15.md`](../v2_audit/_historical/DYNAMIC_FEM_CFD_UPGRADE_ROADMAP_2026-07-15.md)
and [`../v2_audit/FSI_WIRING_DESIGN_2026-07-16.md`](../v2_audit/FSI_WIRING_DESIGN_2026-07-16.md); it
**ties them into a coherent three-layer engine** and **promotes the solute-transport field** (roadmap item
12, RESEARCH tier) to a first-class layer, because that layer is what our 2026-07-16 cortex audit found to
be the true missing "living" piece.

---

## 0. Motivation — the frozen-cortex diagnosis

Today's FF engine is a **fine-grained STATICS solver**: it resolves how force flows through a *fixed*
crosslinked semiflexible truss and relaxes it to overdamped equilibrium. It is high-fidelity in
*mechanics* and lumped/absent in *kinetics & statistical mechanics*. The three defining gaps:

| Gap (cortex audit 2026-07-16) | Consequence |
|---|---|
| **No actin turnover** — filament set frozen at init; only rest-length creep, no (de)polymerization/sever/nucleate on the cortex | No treadmilling, no cortical flow, no stress relaxation, no bleb healing, no remodeling |
| **G-actin = fixed scalar** — no monomer pool that is consumed/released/transported | Assembly is not flux-limited; the turnover economy cannot close |
| **No filament–filament excluded volume** — cortex filaments interpenetrate | Mesh has no real steric volume under compression/crowding |

These are exactly the behaviors a *living* cortex is defined by. They are also exactly what a coupled
**fluid + monomer-transport** substrate supplies. The frozen-cortex diagnosis and the CFD-transport pivot
are the **same conclusion reached from two directions**.

> **Scope note — this is separate from the γ-floor.** Per-head Bell/KMC machinery and an aggregate myosin
> control exist, but the production head-resolved NMII topology is an I3 port. The prior γ-magnitude deficit
> was traced to force-generation/density, and turnover/xlink/buckling were **ruled out** as its lever
> (`memory/project-gamma-floor-likely-deficit`). Do **not** conflate "living dynamics" with "fix γ." This
> engine makes the cortex *alive*; it is not proposed as the γ fix.

---

## 1. Architecture — three layers, one cell

| Layer | Physical entity | Description | Continuum or discrete | Status |
|---|---|---|---|---|
| **① Fluid** | cytosol | Biot poroelastic `p`, relative discharge `q`, and pore-fluid velocity `v_f` on a conservative field grid; inertialess Darcy/Biot, with local Brinkman only if justified | **Continuum** (homogenized porous phase) | **COMPONENT SCAFFOLD ONLY** — device arrays/Peskin IBM exist, but the held-face fixed-box update is not a conservative moving-cell solver |
| **② Solute transport** | G-actin monomer economy | mass-conserving **reaction–advection–diffusion** field `c(x,t)` (later ATP/ADP-actin, profilin, cappers) | **Continuum** (correct — a concentration field) | **GAP — this plan's centerpiece** (G-actin currently a fixed scalar) |
| **③ Solid** | actomyosin | fine-grained filaments / crosslinks / explicit motor heads | **Discrete** (mandatory — the principle) | passive mechanics largely built; production head-resolved NMII and turnover remain I3/I4 work |

The three share **one coordinate frame and boundary geometry, not one data structure**:

- A structured finite-volume/cut-cell grid (or a demonstrably conservative sparse equivalent) stores `p`,
  `q`, `v_f`, `φc`, masks, and fluxes for ①/②.
- A Warp `HashGrid` serves dynamic crosslink/branch partner and excluded-volume searches for ③.

IBM spread/interpolation maps between them. Reusing coordinates and rebuild schedules is encouraged; using a
neighbour-search hash as the PDE field grid is not.

---

## 2. Principle boundary — why this is NOT v1 (non-negotiable guardrail)

> **Continuum is used ONLY for genuine continua (cytosol fluid ①, solute concentration ②). Filaments,
> motor heads, and crosslinks stay explicit discrete objects ③. The particle↔grid transfer is *coupling*
> (immersed-boundary), never constitutive coarse-graining.**

v1 (Taichi MLS-MPM) made the **cell material itself** a continuum — material points *were* the
cytoskeleton. That is the abstraction this project inverted away from. Here the grid carries only fluid +
solute; the cytoskeleton never becomes a field. Corollaries:

- **No active-gel continuum for the cortex network.** `polarization_activegel.py` stays a 1-D
  diagnostic/reference, not the runtime cortex.
- **Myosin stays explicit heads** (Warp-GPU Hand-KMC / minifilament with Hill FV), never a
  `σ_active = −ζQ` continuum tensor or aggregate linear production link.
- This guardrail is what distinguishes a *CFD-coupled fine-grained engine* (the goal) from a *CFD engine*
  (a v1 relapse). Any proposal that turns a filament into a field is out of scope.

---

## 3. Governing equations (per layer)

### ③ Solid — overdamped mechanics + NEW turnover
For an explicit inner mechanical solve, a model point obeys the overdamped form below; an I0-A-authorized
implicit solve must converge the same frozen-outer-state residual rather than reinterpret iterations as time:

```
γ_solid ẋ = F_bend + F_WLC(inextensible) + F_xlink + F_myosin + F_pressure(①) + F_steric(new)
```

- `γ_solid` is a numerical/physical solid mobility for the inner mechanical law, not an imposed effective
  cytoplasmic viscosity. Retire the old `6πηR → 65.9 Pa·s` narrative when ① lands; poroelastic coupling is
  validated through storage, flux, dissipation, and rate dependence, not by forcing an "emergent viscosity."
- **NEW — turnover (③↔②).** Topology mutates via a node/bond **POOL + active mask + KMC event loop**
  (roadmap item E), extending `hand_kmc.py` + `polymerization_warp.py` + `motility_warp.py` treadmill
  kernels:
  - barbed-end growth `v₊ = δ·(k_on·c_local − k_off₊)` — **flux-limited by ②'s `c_local`** (Pollard 1986
    kinetics, already coded).
  - pointed-end depolymerization; cofilin severing (force/age-dependent); Arp2/3-branch & formin
    nucleation.
  - **Every assembly event consumes** monomer from ② at the overlapping grid cell; **every disassembly
    event releases** it. This is the mass-exchange edge that closes the living loop.
- **NEW — filament–filament excluded volume.** Soft repulsive (WCA-like) via the device hash-grid
  neighbor search (currently only MT-tip↔cortex has steric; `soft_contact_kernel`).

### ① Fluid — Biot poroelasticity (production rebuild required)
```
q   = φ(v_f − v_s) = −(k/μ)(∇p − ρb)                   relative Darcy discharge
v_f = v_s + q/φ                                           absolute pore-fluid velocity
σ_total = σ_eff − αpI
0 = ∇·σ_total + f_active,solid + f_external               quasi-static solid balance
S p_dot + α∇·v_s + ∇·q = s_water                          fluid-content balance
```
- `S=1/M`; only in a reduced fixed-domain, constant-coefficient consolidation oracle does
  `c_v=k/(μS)=kM/μ`. Do not prescribe `c_v`, `k`, `M`, and `μ` independently to make a relaxation gate pass.
- `∇·(q+v_s)=0` is only the `α=1, S=0` incompressible-constituent limit, not the production equation.
- The domain is the live membrane minus the live nuclear-envelope inclusion. The membrane supplies the
  water-flux BC; the nucleus has relative no-flux. A fixed no-flux cube is an oracle only.
- Myosin is a solid-phase force. It affects fluid through solid deformation/storage and must not also be
  inserted as an unrelated Darcy body force. Pressure/solid transfer is derived from `−αpI` with adjoint
  spread/interpolation and pressure-work checks.

### ② Solute transport — mass-conserving RAD (NEW, centerpiece)
```
∂(φc)/∂t + ∇·(φc v_f − φD_c∇c) = R_polymer
```
- `c` = G-actin monomer concentration field (**start single-species**; add ATP/ADP-actin, profilin,
  capping later).
- `D_c` = crowding-reduced monomer diffusivity (draft evidence exists; ratification and cell/assay
  applicability remain open — see §7).
- `v_f`, not `q`, is the advection velocity. There is no prescribed-velocity production or scaffold branch.
- `R = −Σ_barbed(k_on·c·δ) + Σ_depoly/sever(release)` — the two-way coupling to ③, applied at grid cells
  overlapping filament ends via IBM spread/interp.
- **HARD conservation invariant:** `d/dt [ ∫φc dV + (monomer bound in filaments) ] = 0` to machine
  precision. Total actin is conserved — the gate that proves the economy is real, not a source/sink leak.
- Membrane BC: **impermeable to monomer** (no-flux); water/small solutes cross (aquaporin) — couples to
  ①/osmotic.

---

## 4. Two-way coupling map (the heart)

```
③ →② : polymer (de)assembly = source/sink R_polymer       (mass exchange)
② →③ : local c sets barbed-end growth rate                (flux-limited assembly → turnover)
① →③ : −αpI stress/traction transferred with adjoint IBM  (pressure work on the skeleton)
③ →① : α∇·v_s enters fluid-content balance                (solid deformation drives relative flow)
① →② : advection by absolute pore-fluid velocity v_f      (streaming carries monomer)
②→ ① : (weak) solute osmotic contribution                 (defer / one-way initially)
```

Field coupling routes through the conservative grid and IBM operators; discrete topology/steric queries use
the separate HashGrid (§1). A component OFF state may remain as a regression oracle, but the physiological
production configuration keeps the fluid and monomer layers ON from their landing increments.

---

## 5. Phasing — each milestone is a validation gate, native-scale-gated

Every stage runs and is judged on the **I0-A-ratified native full cell** with physiological initial/boundary
conditions: MCF7×collagen×α2β1 and 70,686 active cortical F-actin filaments. Other actin/MT/IF populations are
additive only when they have unique global IDs; structures classified from the unified actin pool are not
counted twice. Nuclear lamina/chromatin, ECM fibers/cross-links, meshes, fields, motor heads, and clutches are
separate device allocations in the same ledger. Coarse runs are non-authoritative arithmetic smoke tests only.

- **P0 — Infra + regression oracle.** Build the conservative field grid and the separate neighbour HashGrid
  in one coordinate frame. Test-only all-new-components-OFF must match the frozen **Warp FF** arithmetic
  reference; no HOOMD executable participates. *(master I1a scaffold)*
- **P1 — ① lands (fluid = master I1a+I1b).** Replace the held-face stencil with a conservative moving-domain
  fluid-content solver, solve `q`, reconstruct `v_f`, and land the outer-physical/inner-mechanical scheduler.
  Retire the `6πηR → 65.9 Pa·s` narrative in the same change. **Gate:** constant-state and moving-boundary
  manufactured solutions; global storage/flux balance; Darcy slab flux; poroelastic `τ_p` and
  drained↔undrained rate dependence.
- **P2 — ② conservative RAD monomer on solved `v_f` (master I1c).** Single G-actin field; barbed-end growth
  flux-limited by `c`; assembly consumes and disassembly releases exact monomer counts. **Gate:** total-actin
  conservation, positivity, uniform-state preservation on the moving domain, treadmilling balance, and FRAP
  recovery. No prescribed-velocity branch is created.
- **P3 — full ①↔②↔③ integration (master I1b2 optional + I2b + I9).** Add the PI-KU-gated explicit
  relative-velocity drag pair only if authorized, plus filament–filament excluded volume and all-compartment
  production coupling. **Gate:** adjoint transfer/dissipation/net-force projection, turnover steady state,
  label-blind cortical/retrograde flow report, and native-scale stability.
- **P4 — multi-species + regulation (RESEARCH / PI-gated).** ATP/ADP-actin cycle, profilin, capping,
  cofilin spatial fields; membrane Scriven-Boussinesq lipid flow; spatial RVI/RVD. *(roadmap 11–13)*

---

## 6. Risks / hazards

| Hazard | Mitigation |
|---|---|
| **η double-count** (biggest correctness trap) | Retire `6πηR` node drag the moment ① lands, **same commit** (P1) |
| **v1 relapse** (filament → field) | §2 guardrail: continuum only for fluid+solute; filaments stay discrete; grid = coupling only |
| **Darcy discharge used as solute velocity** | Keep `q=φ(v_f−v_s)` and advect with `v_f`; gate `v_f−v_s=q/φ` |
| **Myosin replace-not-sum** | minifilament REPLACES the constant dipole; active-gel σ must not co-fire |
| **Negative concentration / mass leak** | flux-limited advection (positivity) + the hard conservation gate |
| **Held-face boundary mass creation** | discard the current stencil as production; use conservative face fluxes/cut cells and gate global fluid content on moving boundaries |
| **PDE grid / neighbour-grid conflation** | two data structures, one coordinate frame (§1); conservation tests apply to the field grid |
| **Physical time / solver iteration conflation** | outer clock owns Biot/RAD/KMC; inner mechanical solve reports convergence, never elapsed seconds |
| **Native-scale cost** | publish unique-active, allocated+dormant, node, explicit-head/state, field-cell, and exact GPU-byte high-water ledgers; optimize GPU layout or hardware, never biological density |
| **Coarse-scale artifacts** | every finding reconfirmed on the native full cell (coarse cortex self-collapses) |

---

## 7. Parameters & KB gaps (surface to PI — no defaulting to nulls)

| Parameter | Symbol | Value / status | Source |
|---|---|---|---|
| Solid inner-solve mobility | `γ_solid` | implementation candidate only; must be compatible with the I0-A solver/time policy and not be relabeled bulk viscosity | `units.py:58-78`, I0-A |
| Biot modulus/storage | `M`, `S=1/M` | existing `1 kPa` code value is not ratified for production; derive with α and drained/undrained moduli, e.g. `K_u=K_d+α²M` | I0-B1 |
| Permeability / consolidation | `k`, `c_v` | use `c_v=kM/μ` only in the reduced constant-coefficient oracle; `k`, `μ`, `M`, `α` need independent evidence | I0-B1 |
| Porosity / water source | `φ`, `s_water`, `L_p` | open; needed for `v_f`, storage, and moving membrane BC | I0-B1 |
| **Monomer diffusivity** | `D_c` | draft evidence trail suggests crowded-cytoplasm ~2–6 µm²/s; ratification and MCF7/assay applicability open | I0-B1c |
| **Initial monomer pool** | `c₀` | **PI GAP** — total actin ~100 µM, monomer pool buffered ~tens of µM; **MCF7-specific value needed** (physiological-baseline rule) | — |
| Assembly kinetics | `k_on`, `k_off` | code paths exist; values/assay context still require I0-B4 ratification | `polymerization_warp.py`, I0-B4 |
| Biot coupling | `α` | open; do not silently set 1 outside the incompressible-constituent limit | I0-B1 |
| Grid resolution | `Δx` | vs mesh ξ (20–50 nm) vs 16 GB — **sizing decision** | — |

---

## 8. Reuse ledger — don't rebuild

- **`ff/biot_fluid_warp.py`** → reuse device arrays and Peskin IBM kernels only; rebuild the conservative
  moving-domain field update.
- **`network_warp.py`** FSI kernels (`fsi_*`) + `biot_fsi` hook → wiring/parity scaffold, subject to adjoint
  work and net-force gates.
- **DCM `dcm_neighbor_warp`** HashGrid → discrete neighbour search only, not the PDE field grid.
- **`hand_kmc.py` + `polymerization_warp.py` + `motility_warp.py`** (treadmill kernels) → ③ turnover seed.
- **`membrane_surface.py`** → live membrane geometry and Helfrich/area/ERM mechanics for hydraulic BC
  placement; Scriven surface-fluid mechanics remain deferred.
- **`units.py`** (`ETA_SOLVENT`, `fiber_point_drag`, NF2007 mobility) → drag/coupling constants.
- **Historical non-Warp implementations** → equations/results may be read as design references only; they are
  never executed. Validation is analytic/manufactured/experimental and runs through Warp-CUDA kernels.

---

## 9. Validation & sanity gates (written before building — project rule)

1. **Mass conservation** of total actin (monomer + polymer) to machine precision — **HARD**.
2. **Low-Re justification** — Re ≪ 1 computed; inertial terms dropped and shown negligible.
3. **Poroelastic consolidation** vs Biot/Terzaghi analytic (P1).
4. **Moving-domain conservation** — integrated fluid-content change equals membrane water flux; impermeable
   case conserves content; uniform solute remains uniform under domain motion (P1/P2).
5. **Darcy kinematics** — slab flux, `v_f−v_s=q/φ`, and solute front follows `v_f`, not `q` (P1/P2).
6. **FRAP** monomer recovery time ↔ `D_c` (wet-lab Aim 3 oracle) (P2).
7. **Treadmilling steady-state** length ↔ `k_on·c = k_off` (P2).
8. **Cortical/retrograde flow** velocity reported against a literature band without tuning (P3).
9. **Positivity** of `c` (no negative concentration).
10. **Regression oracle** — all-new-components-OFF == frozen Warp-FF arithmetic reference (P0, test-only;
    no HOOMD/non-Warp execution).
11. Every magnitude gate on the **I0-A-ratified native full cell**.
12. **GPU residency/population** — zero authoritative per-step host state/roundtrips; exact active vs allocated
    population and peak device bytes logged for every compartment.

---

*Companion docs: [`../v2_audit/_historical/DYNAMIC_FEM_CFD_UPGRADE_ROADMAP_2026-07-15.md`](../v2_audit/_historical/DYNAMIC_FEM_CFD_UPGRADE_ROADMAP_2026-07-15.md)
(fluid + dynamic-fiber engineering axes), [`../v2_audit/FSI_WIRING_DESIGN_2026-07-16.md`](../v2_audit/FSI_WIRING_DESIGN_2026-07-16.md)
(① wiring spec), [`../v2_audit/_historical/CELL_MECHANICS_FRAMEWORK_2026-07-15.md`](../v2_audit/_historical/CELL_MECHANICS_FRAMEWORK_2026-07-15.md)
(PI framework), [`HYPOTHESES.md`](HYPOTHESES.md) (H1–H12), [`manuscript/MANUSCRIPT.pdf`](manuscript/MANUSCRIPT.pdf).*

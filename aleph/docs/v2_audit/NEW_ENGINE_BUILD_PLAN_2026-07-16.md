# NEW ENGINE — "Active Cell" build plan (handoff for a fresh session)

**PI decision 2026-07-16:** the FF engine, carried forward as-is, was a trap ("빛좋은 개살구" — a passive
mechanical shell with lumped active mechanics, decoupled interior, and the cytosol as a balloon). **Build a
NEW engine.** FF's fate (archive vs keep-as-reference) is DEFERRED — decide after the new engine's foundation
is validated. This doc is the complete synthesis of every plan so far (the 18-agent option-C design workflow +
the re-audit + the fluid-first and deformable-mesh corrections) into one self-contained build plan a fresh
session boots from and executes increment-by-increment with native gates.

**I0-A PI ratification 2026-07-16 (authoritative):** simulation runtime is **NVIDIA Warp on CUDA GPU only**.
HOOMD is never executed—not for production, development, parity, fallback, benchmarks, or validation. Archived
HOOMD files may be read only as a port specification. The production NMII primitive is an explicit
backbone+individual-head Stam-Hocky minifilament with **Hill force–velocity**, ported completely to Warp.
The first baseline is **MCF7 + collagen ECM + α2β1–collagen clutch**, with **70,686 active cortical F-actin
filaments**. All other populations are derived from compartment-specific physiological density×geometry;
there is no global ×40 production scale or memory-driven biological down-count.

Read with: `_historical/CYTOSKELETON_CYTOSOL_REAUDIT_2026-07-16.md` (the sins), `_historical/STRESS_FIBER_TARGET_2026-07-16.md` (SF
taxonomy + emergent formation), `_historical/MEMBRANE_HELFRICH_DESIGN_2026-07-16.md` (the deformable-mesh pattern already
built), `_historical/CELL_MECHANICS_FRAMEWORK_2026-07-15.md` (the professor's load-transfer framework). Memory:
[[project-active-cell-c-redesign]], [[project-field-actuation-and-actin-gf-gap]] (the G-actin monomer-field gap).

---

## 0. Why a new engine (not in-place FF evolution)

The re-audit proved the FF whole-cell model = cortex-shell + passive-MT + nucleus-fuzzyball + IF-cage +
membrane, with: (a) the ONLY active force a **lumped constant `f_myo`**; (b) the interior **mechanically
decoupled** (RAW 0.0 nm under load); (c) the cytosol a **single scalar ΔP from hull volume applied uniformly**
(a balloon); (d) the nucleus a **radial-spring bead cloud with no internal connectivity**; (e) actin/MT
**block-diagonal** (tip-contact only); (f) the active Biot CFD **default-OFF** in every landed native run.
The fine-grained parts (cortex bending, crosslinks, inextensibility, membrane Helfrich, FA catch-slip clutch,
Biot kernels) are all **passive structure**; everything **active + interior-coupling** is lumped/placeholder.
The architecture is wrong at the foundation — patching channels onto it reproduces the "additive-default-off
becomes de-facto production" trap. A new engine fixes the ARCHITECTURE and REUSES FF's validated kernels as
components.

---

## 1. The five architectural pillars (the new foundation)

**P1 — FLUID-FIRST poroelastic substrate that FLOWS (PI: "애초에 기반은 CFD 방식이여야한다" + "흐르는 것까지
되어야한다").** The cell is a fluid-saturated porous continuum coupled to explicit filaments. The fluid is a
**first-class, always-solved DOF** (not an optional channel), and it must genuinely **FLOW** — resolve the
pressure, relative discharge, and pore-fluid velocity fields, not pressure alone.

**Why Darcy/Biot and NOT Navier-Stokes (the correct physics level).** Intracellular Re = UL/ν ≈ 1e-5…1e-10 —
momentum inertia is ~10⁶–10¹⁰× negligible. Darcy/Biot is a **homogenized constitutive level**, not an exact
pore-scale identity: it is valid only when scale separation and a permeability closure are demonstrated. The
fluid-side hierarchy is:
```
Navier-Stokes (ρ∂u/∂t + ρu·∇u = −∇p + µ∇²u)
   │  Re≪1 → drop momentum inertia (solute advection is a separate Peclet question)
   ▼  Stokes (−∇p + µ∇²u = 0, ∇·u=0)             ← the correct free-fluid eq at Re≪1
   │  homogenize Stokes flow THROUGH the porous cytoskeleton (Brinkman/volume-average)
   ▼  Darcy (q = −(k/µ)∇p)                        ← homogenized relative discharge through the pore network
   │  couple to the deforming elastic skeleton
   ▼  Biot poroelasticity                          ← what we use = the correct coupled theory
```
So Biot is the low-Re, pore-homogenized model appropriate to the target scale; it deliberately discards
pore-scale velocity detail. Using inertial NS for the intracellular bulk would add negligible physics at huge
cost, while a local Brinkman/Stokes correction remains available where homogenization fails. ⭐
**The "higher level" the flow requirement points at is NOT NS (add inertia — wrong axis) but a mixed q-p
form of Biot plus reconstructed pore-fluid velocity `v_f`: still inertialess, still Darcy, but it resolves
the flow field needed for streaming and advection.**

**Production variables and closure (authoritative notation).** Keep these symbols distinct everywhere:
```
v_s = solid-skeleton velocity
q   = φ(v_f − v_s) = −(k/µ)(∇p − ρb)                       Darcy discharge RELATIVE to solid
v_f = v_s + q/φ                                             absolute pore-fluid velocity
σ_total = σ_eff − αpI                                       no fixed "half-stress" fraction
0 = ∇·σ_total + f_active,solid + f_external                 quasi-static solid balance
S p_dot + α∇·v_s + ∇·q = s_water                            fluid-content balance; S=1/M
∂(φc)/∂t + ∇·(φc v_f − φD_c∇c) = R_polymer                 conservative solute balance
```

`∇·(q+v_s)=0` is only the `α=1, S=0` incompressible-constituent limit, not the general I1b equation. Myosin
acts on the **solid**; it drives fluid through `∇·v_s` and the coupled balance. It must not also appear as a
direct Darcy body force unless a distinct fluid-phase force is derived. Optional near-boundary **Brinkman**
shear (µ∇²v_f) is a local Stokes correction, not an inertial NS fallback.

**Pressure split and moving domain.** Use `p = p_ext + p_bar + p_excess`, where `p_bar` is the osmotic/hydrostatic
mean fixed by the membrane water-flux/volume closure and `p_excess` is the spatial poroelastic field. Do not
force `p_excess` to zero mean during an undrained transient. The domain is the LIVE membrane mesh minus the LIVE
nuclear-envelope inclusion; both boundaries move. The outer membrane supplies the hydraulic-flux BC and the
nucleus supplies a relative no-flux BC. A fixed cubic no-flux box is an analytic oracle only, not the cell.

**Staged fluid build (all IN the plan; §3 I1a/b/c):**
```
I1a  CONSERVATIVE p/mass Biot  moving cell domain + storage/solid source + membrane/nucleus flux BC
I1b  mixed q-p / v_f form      q=−(k/µ)∇p; v_f=v_s+q/φ; the fluid genuinely FLOWS
I1c  conservative transport    ∂(φc)/∂t+∇·(φc v_f−φD∇c)=R; G-actin + small solutes ride v_f
```

**Two spatial accelerators, one coordinate frame — not one data structure.** The Eulerian Biot/RAD fields use a
structured conservative grid (or a demonstrably conservative sparse equivalent). Filament partner search and
excluded volume use Warp `HashGrid`. They may share origin/cell geometry and rebuild schedule, but a neighbour-
search hash grid is not the PDE field grid.

⭐ **The flow is REQUIRED, not optional:** the G-actin monomer field needs `v_f` to advect monomers — pressure-
only cannot carry them; the emergent network grows from this transported pool. So I1b
(velocity) is a PLANNED stage, not a deferral. Honest scope: at drained steady state the pressure deviation is
small (~Pa) — the fluid's payload is **flow/streaming, advective transport, rate-dependence, τ_p transients,
bleb, hydraulic long-range coupling**, not dominating the static force balance. **The mean-pressure vs
turgor-double-count tension (zero-mean discards Skempton undrained stiffening) is a real design decision to
resolve.** Needs a PI-authored **permeability-k** KnowledgeClaim for I1b (added to §4).

**P2 — ACTIVE FROM THE ROOT.** ONE head-resolved, multi-particle NMII minifilament topology is the only
**actomyosin contractile-motor** primitive; SF, cap, and cortical tension reuse it. Dynein, kinesin, and the
polymerization ratchet are distinct motor/force families and must never be aliases of NMII. No lumped `f_myo`
or `N_side·F_head` single-link substitute anywhere in the final engine. Honest scope: at the quasi-static
equilibrium `v_slide→0` so force→F_stall (a constant) — the fine-grained content that SURVIVES at the settled
state is **derived magnitude + Bell-emergent engaged fraction + topology reformation**; force-velocity only
bites the transient. State this; do not claim FV changes the settled state.

**P3 — ONE UNIFIED EMERGENT NETWORK.** A single actomyosin network from which cortex / ventral+dorsal SF /
transverse arcs / perinuclear cap / lamellipodium (dendritic Arp2/3 REBUILD, §3A.b) / filopodium
**self-organize** (myosin-driven alignment + crosslink KMC reattach + FA anchoring). The physiological
manifold, nucleator identity/location, and polarity field may be seeded, but local filament orientations are
isotropic conditional on those inputs and no pre-made bundle is supplied. Bundle/arc/cap/filopodium labels
must be assigned by a label-blind detector after the run; bundles must be shown to CONDENSE (the falsifiable
emergence gate). Scripted pre-made bundles are DEMOTED to non-authoritative controls. ⚠️ **Deepest risk: emergence
may not converge in an overdamped quasi-static solve** — fallback = seeded scaffold explicitly labeled SEEDED.

**P4 — DEFORMABLE-MESH compartments (PI: nucleus를 DCM radial 말고 진짜 변형 mesh로).** DCM's nucleus is the
SAME bilinear-radial law — reusing it gives nothing. The right pattern is the **membrane Helfrich machinery
already built** (icosphere triangulated shell + dihedral bending + area tension + volume + tether). Apply it
to the **nucleus**: a triangulated nuclear-envelope (lamina) shell with bending + area-elasticity +
nucleoplasm volume (ν→½) + LINC tethers + chromatin — NOT a radial-spring bead cloud. Membrane = done; nucleus
= port the same machinery. This is a foundational compartment upgrade, sequenced early.

**P5 — INTEGRATED COUPLING (no block-diagonal) + ONE PHYSICAL CLOCK.** Cross-region bonds are first-class: actin↔MT crosslinker
along the MT length, MTOC↔nucleus, cap→LINC→nucleus, FA→SF→ECM. The whole-cell stiffness is NOT
block-diagonal; compartments are mechanically interdependent. Runtime uses an **outer physical-time loop** for
Biot/RAD/KMC/water flux/loads and an **inner mechanical solve** at frozen outer state. PI ratified a
Warp-GPU implementation that may use either the validated explicit or implicit Warp solver; the selected
solver must converge the same registered residual, and solver iterations are never silently reinterpreted as
seconds. No CPU or HOOMD stepping path exists. The off-diagonal Jacobian remains an analytic coupling gate.

---

## 2. Reuse manifest — port these validated FF kernels as COMPONENTS

| Reuse (fine-grained, validated) | From | Role in new engine |
|---|---|---|
| Cytosim bending `F=α(m₋−2m₀+m₊)` | `aleph/laws/forces_warp.py` | filament bending (all actin/MT) |
| Nonlinear finite-ext link spring | `aleph/laws/network_warp.py::link_spring_kernel` | crosslinks / LINC / FA / anchors |
| Inextensibility reshape projection | `aleph/laws/network_warp.py::reshape_kernel` | filament segment-length constraint |
| **Membrane Helfrich sheet** (dihedral + area + ERM) | `aleph/laws/membrane_surface.py` | membrane ✓ AND the nucleus-envelope template |
| FA catch-slip KMC (Pereverzev/Bell) | `aleph/laws/fa_clutch_warp.py`, `fa_ecm.py` | FA→ECM traction |
| Biot device arrays + Peskin spread/interp kernels | `aleph/laws/biot_fluid_warp.py` | IBM transfer/parity scaffold only. The current held-face diffusion stencil is **not** a production no-flux solver; rebuild the conservative moving-domain step in I1a |
| Per-head Bell/KMC state machinery | `aleph/laws/hand_kmc.py` | kinetic building block for the active-force root; port to device-resident per-head state |
| Head-resolved Stam-Hocky topology + stepping semantics | `archive/hoomd_legacy/cortex/myosin.py`, `archive/hoomd_legacy/bridge/motor.py` | read-only geometry/kinetics port specification for I3; HOOMD is never imported or executed |
| `weave()` region builder | `aleph/laws/weave.py` | per-region seed for the unified network (P3) — DEMOTED to seed |
| ECM Mikado library | `aleph/laws/ecm_library.py`, `ecm_mechanics.py` | substrate/matrix |
| Warp GPU infra, viewers, KB, oracles | `ff/`, `outputs/`, `validation/oracles/` | infra + acceptance layer (KEEP — not lumped runtime) |

**REBUILD (do NOT port):** the lumped `f_myo` myosin **and the aggregate `myosin_linear.minifilament_kernel` as
the final motor**, the radial-shell nucleus (`nucleus_shell_kernel`), the
scalar-turgor cytosol (`_refresh_turgor`+`turgor_kernel` as the base), the block-diagonal MT merge, the failed
IF cage (polar radial spokes). These are the "개살구" parts.

**ADD (new, not a port — the FF cortex lacks it):** **filament–filament excluded volume** (steric hash-grid
pair potential, I2b) — a CLAUDE.md HARD worked-example ("LJ repulsive ON from Phase 1") and one of
ENGINE_ARCHITECTURE §3's three defining living-cortex gaps; today only MT-tip↔cortex has steric
(`soft_contact_kernel`), so generalize that kernel to all fibers via the shared device hash-grid.

**Suggested package:** `aleph/ac/` (Active Cell), parallel to `ff/` and `dcm/`. Final name = PI's call.

### 2b. Component keep/kill judgment (all runtime-relevant `ff/` modules — PI 2026-07-16 "쓸 수 있는 것 / 절대 못 쓸 것")

**✅ REUSE — validated fine-grained runtime components (port as-is):**
| Module | What | Note |
|---|---|---|
| `forces_warp.py` | Cytosim bending kernel | Cytosim-parity validated |
| `fiber_network.py` | Cytosim fiber data model | engine-agnostic base |
| `constraints.py` + reshape | inextensibility (hard constraint, NF2007) | not a penalty spring |
| `relax.py`, `implicit_ff.py` | overdamped relaxers (explicit CFL + implicit) | candidate inner mechanical solvers; I0-A permits either only after a complete Warp-CUDA port and common-residual convergence gate |
| `membrane_surface.py` | **Helfrich sheet** (dihedral+area+ERM) | done today; ALSO the nucleus-envelope template |
| `hand_kmc.py` | Bell catch-slip Hand/KMC | the KMC base |
| `fa_clutch_warp.py`, `fa_ecm.py`, `fa_maturation.py` | FA catch-slip clutch | native pN-validated |
| `ecm_library.py`, `ecm_mechanics.py`, `ecm_mikado.py` | ECM Mikado | 6 materials in-band |
| `network_contractility.py` | buckling contractile network | relevant to emergence |
| `polarization_activegel.py` | symmetry-breaking / polarization | motility |
| `motility_warp.py` | overdamped substrate crawl | motility |
| `cortex_assembly.py` | cortex shell builder | REUSE as an isotropic SEED (P3), demoted |
| `microtubule.py` | κ-based MT aster/beam mechanics | static-lattice MT seed and I8 strut mechanics |
| `wlc.py` | finite-extensibility / enthalpic-wall helpers | filament/ECM constitutive helper and oracle cross-check |
| `substrate.py` | substrate geometry/runtime helpers | ECM/FA boundary composition, after ligand-specific audit |
| `piezo.py` | tension-gated mechanosensing readout | future signaling layer |

**🔧 REUSE-WITH-FIX — good core, needs a specific fix before it rides the new engine:**
| Module | Fix required |
|---|---|
| `biot_fluid_warp.py` | retain device-array/IBM pieces; **replace** held-face diffusion with a conservative moving-domain finite-volume/cut-cell step, in-place field update, membrane hydraulic flux, and live-nucleus relative no-flux — I1 |
| `nucleus_envelope.py` | partial surface scaffold only: add/wire Helfrich bending, per-face reference strain + local rupture state, volume/nucleoplasm, chromatin polymer, and conservative fluid mask — I2 |
| `weave.py` + `architecture_spec.py` | keep as region SEEDS only; DEMOTE `bundle`/`_build_bundle` from "final structure" to "initial isotropic seed" — I4 |
| `polymerization_warp.py` | barbed-end ratchet is fine, but `G_actin` is a **fixed scalar** — wire the transported G-actin **monomer field** (P1) — I1/I9 |
| `architecture_metrics.py` | `bundle_count`=n_fibers is VACUOUS — build a real condensation/nematic detector — I5 |
| `fa_anchor.py` | replace the pole-grabbing `fa_end_nodes` with a **z-lo basal areal-density sampler** — I6 |
| `gamma_floor.py` (cortex-build parts) | `build_crosslinked_cortex` reusable as a seed; but its `_myosin_force` (lumped) is KILLED |
| `solver/deflated_pcg.py` | prior benchmark found no production win; keep as an optional solver experiment/oracle, not the default inner solver unless re-benchmarked on the new coupled operator |

**❌ CANNOT USE — the "개살구" architecture; rebuild, do NOT port:**
| Thing | Why | Replaced by |
|---|---|---|
| `network_warp.py::simulate_whole_cell_compression_on_device` | the passive-shell ASSEMBLY itself | the new fluid-first active assembly |
| `myosin_kernel` (lumped constant `f_myo`) | not a motor, a constant | head-resolved multi-particle I3 motor |
| `myosin_linear.minifilament_kernel` **as the final motor** | one aggregate two-anchor force; no explicit backbone or per-head geometry | head-resolved multi-particle I3 port + device Hand state |
| `turgor_kernel` **as the cytosol base** | one scalar from hull volume = balloon | fluid-first Biot substrate (P1); turgor demoted to the MEAN pressure channel only |
| `nucleus_shell_kernel` (radial bead ball) | no internal connectivity, one deformation mode | deformable-mesh nucleus (P4) |
| `merge_aster_into_cortex` (block-diagonal MT) | actin/MT elastically independent | integrated actin-MT coupling (P5) |
| `intermediate_filaments.py` (IF cage) | FAILED (buckles, over-stiffens 6×) as the nucleus load path | perinuclear actin cap → LINC (I7); IF returns as a **cell-type keratin/vimentin secondary passive net at first-validation** (§3A.c), NOT the nucleus driver |

**KEEP (acceptance/infra layer — NOT runtime, retained per the founding oracle rule):**
`cytosim_parity.py`, `kim_network.py` (peer oracles), `gamma_estimator.py`, `ff_virial_stress.py` (measurement),
`gamma_floor_{dynamic,sweep}.py` (historical diagnosis/acceptance only), `cell_type.py` (presets/config),
`units.py` (conventions — but purge the retired-6πηR/γ_solid emergent-viscosity narrative), `viz_*.py`
(viewers), `validation/oracles/` (analytic ground truth). These validate the new engine; they do not define its
runtime mechanics.

---

## 3. Increment plan (fluid-first reordered; every increment ends on a NATIVE `--from-resting` full gate)

**Ordering rationale:** the workflow put CFD late (I6) because its steady-state effect ≈ FSI-off — but the PI's
fluid-first pillar makes the poroelastic substrate the FOUNDATION. Reconciliation: build the **substrate**
(solid+fluid coupled base + deformable-mesh compartments) FIRST as the platform, then the active-force root,
then emergence, then load paths, then the combined production run. Each magnitude verdict is HELD until its
I0-Bn parameter record closes.

| # | Increment | Delivers | Gate (analytic-first) |
|---|---|---|---|
| **I0-A** | **Foundation contract ratification** (no code) — **PI-RATIFIED 2026-07-16** | Warp-CUDA-only runtime; no HOOMD execution; one outer physical clock + converged Warp explicit/implicit inner solve; head-resolved NMII + Hill FV; 70,686 cortical filaments; MCF7×collagen×α2β1 baseline. Update AGENTS/CLAUDE/Dashboard/companions together | PI receipt + authority-doc hash recorded; no unresolved foundation choice remains |
| **I0-Bn** | **Per-increment parameter gate** | Before increment `n`, bind only its required §4 parameters to a KnowledgeClaim/evidence status; downstream unknowns do not block earlier analytic work | parameter ledger has claim ID, evidence status, owner, provisional policy, blocking increment |
| **I1a** | **Conservative fluid-first substrate (p/mass)** | Build a GPU-resident conservative pressure/fluid-content solver on the LIVE membrane-minus-nucleus domain. Reuse device/IBM pieces, **not** the held-face diffusion stencil. Include `S p_dot+α∇·v_s+∇·q=s_water`, in-place fields, live outer membrane hydraulic flux, live nuclear relative no-flux, conservative domain remap, and `p=p_ext+p_bar+p_excess`. Delete scalar-turgor-as-base; seed resting `Π₀` at the I0-B1-authorized physiological value | Terzaghi/Green oracle; constant-state preservation; manufactured moving-boundary solution; global fluid-content balance equals integrated membrane flux to machine precision; impermeable limit conserves mass; pressure work has correct sign; OFF reference only for current Warp-component regression; resting `Π₀` pre-tensions shell at t0 |
| **I1b** | **Mixed q-p / pore-fluid velocity (the fluid FLOWS)** | Solve `q=−(k/µ)(∇p−ρb)` and reconstruct `v_f=v_s+q/φ`; use the full storage balance, not unconditional `∇·(q+v_s)=0`. Myosin enters through solid force/deformation, not a second direct Darcy force. This supplies the I1c advection velocity | Darcy slab flux linear in Δp and k/µ; `v_f−v_s=q/φ`; discrete mass residual closes; tracer follows `v_f`; lit streaming pattern. Poiseuille profile is a separate Brinkman/free-fluid gate only |
| **I1b2** | **Velocity-difference fiber momentum drag** (PI-KU-gated; may defer) | Add an equal-and-opposite, control-volume-consistent drag pair derived from `q=φ(v_f−v_s)`. Define nodal representative volume and avoid summing this with an already-equivalent pressure force. ⚠️ no authoritative intracellular `(η_f/k)` claim yet | free fiber approaches `v_f`; pairwise fluid/solid work is dissipative; grid↔node transfer adjoint; net internal force projection closes; OFF==I1b bit-identical |
| **I1c** | **Conservative G-actin/solute transport** | Advance `∂(φc)/∂t+∇·(φc v_f−φD_c∇c)=R` on the same moving cell domain. Wire barbed-end consumption **and pointed-end/depolymerization release** here. Map every polymer length/active-mask change to an exact monomer count so total actin conservation is non-vacuous | total actin `∫φc dV+N_polymer` conserved to machine precision; positivity; pure-diffusion parity; uniform-state preservation under domain motion; FRAP recovery↔D_c; advection front follows `v_f`, not Darcy discharge `q` |
| **I2** | **Deformable-mesh nucleus** | Triangulated nuclear-envelope (lamina) shell (Helfrich machinery from membrane) + area + nucleoplasm volume (ν→½) + chromatin polymer net; **framework-#6 depth: lamin-A/C vs lamin-B split** (small-strain→chromatin, large-strain→lamin-A/C strain-stiffening), **nucleoplasm viscosity**, **envelope-rupture threshold**; replaces the radial bead ball. ⚠️ **Reconcile the I1a Biot no-flux mask** with THIS surface: mask the LIVE deformable oblate nucleus mesh (about `centre_nuc`), NOT a static sphere R_nuc | sphere→8πκ + FD-gradient (like membrane); volume conservation; resting nucleus stable at physiological E_nuc; small→large-strain modulus crossover at the lamin knee (report-not-tune); rupture emerges past a sourced envelope-tension threshold (do NOT tune) |
| **I2b** | **Filament–filament excluded volume (steric)** | Soft repulsive (WCA/LJ) pair potential via the device hash-grid neighbour search so cortex/SF/MT filaments do not interpenetrate. ⭐ **CLAUDE.md HARD worked-example ("LJ repulsive excluded volume ON from Phase 1")** + one of ENGINE_ARCHITECTURE §3's three defining living-cortex gaps ("cortex filaments interpenetrate; mesh has no real steric volume under compression"). Today only MT-tip↔cortex has steric (`soft_contact_kernel`) — generalize it to all fibers | pair-potential FD-gradient (sign arbiter); OFF / zero-overlap == bit-identical (regression); a compressed patch resists volumetric collapse (finite steric volume); no fiber interpenetration on the native cell (crowding gate); add k_EV to kmax (CFL) |
| **I3** | **Active-force root — head-resolved NMII minifilament** | Reimplement the archived Stam-Hocky specification entirely in Warp: explicit backbone beads + explicit heads on both sides + per-head actin attachment, Hill stepping/FV, compliance, and Bell detachment, all CUDA-resident with no per-step host state. `myosin_linear.minifilament_kernel` may be a Warp diagnostic/control only, never the production motor. Head counts/force constants remain I0-B3 evidence decisions; remove all replaced `cortex.myo_i`/`f_myo` links in the same commit | topology/count invariants; per-head Newton closure; single-head Hill/Bell gates; ensemble stall emerges from bound heads (not imposed `N_side·F_head`); ATP/step work sign; GPU-residency/zero-host-roundtrip gate; native γ-floor report-not-tune |
| **I4** | **Unified spine** — `weave_cell()` + `WovenCell` + topology-reforming KMC (hash-grid reattach); **FULL region set** (cortex + ventral/dorsal SF + transverse arc + perinuclear cap + filopodium fascin bundle + **lamellipodium dendritic Arp2/3 net — REBUILT**, §3A.b) | ONE actomyosin network; cross-region bonds concat leaves out; lamellipodium is a genuine build (not the ratchet proxy); ⚠️ pick ONE crosslinker-relaxation channel (new KMC vs `xl_turnover` r0-creep) — running both double-counts | Gate-1: `weave_cell([CORTEX])` all-OFF → bit-identical γ to current cortex; branch-angle dist matches θ₀±σ_θ (§3A.b) |
| **I5** | **Emergence proof** (falsifiable) | With I1c already LIVE: isotropic, label-free patch + seeded nucleation identity/locations + FA anchors + myosin/KMC ON → bundle/nematic order rises. Controls (myosin-off / unanchored / KMC-off / angle-off) stay flat. The detector must not read region/type labels | pre-registered ensemble seed count + effect-size threshold; S/bundle persistence rises vs all controls on smoke then native; per-realisation traces + mean; **flat result → FINDING to PI + SEEDED fallback**, never re-tune |
| **I6** | **Adhesion baseline** — ONE clutch implementation with ligand-specific parameters on the FULL architecture (ventral SF→FA→ECM + lamellipodial nascent + **filopodial tip**); FA maturation; `--from-resting` adherent | Whole FA load path; z-low basal areal-density sampler; two-sided ECM Newton closure. Tip nucleation location is SEEDED, while tip/shaft/basal adhesion **classification** must be a label-blind output of architecture×load×maturation. Choose integrin parameters from the I0-A-authorized ECM ligand and I0-B6 evidence record; never treat α5β1–FN constants as universal collagen values | selected integrin–ligand off-rate oracle; Newton Σf_cell+Σf_ecm=0; BOTH traction channels present; label-blind adhesion classes emerge; per-clutch/FA/cell bands report-not-tune |
| **I7** | **Nucleus load path** — perinuclear actin cap → LINC → nucleus flatten | Contractile cap draped over the (now deformable-mesh) nucleus; tangential-tension→normal pressure (capstan), NOT a radial strut; **re-frame: flatten measured at the RESTING adherent baseline** (Khatau), not chased at AFM strain | analytic capstan ∮T·κ ds; native flatten at rest (oblate band, volume conserved); cap-OFF stays spherical; ⚠️ **gate the COMBINED cap+IF baseline** (double-load-path; IF = cell-type keratin/vimentin secondary passive net, §3A.c — first-validation, MCF7↔MDA contrast); ⚠️ **force-scale: corrected K_nuc≈25 nN/µm → 5-6 nN cap gives ~0.2µm, ~10× short — report honestly, do NOT tune** |
| **I8** | **Actin-MT coupling** — un-block-diagonalize + MTOC→nucleus + cortical dynein (POPULATION of motor hands, NOT lumped N×f) | plectin/MACF actin↔MT crosslinker on link_spring; MTOC↔nucleus anchor; dynein as individual hands (engaged count emerges) · **kinesin anterograde transport + MT dynamic instability** (framework #3) = growing/shrinking-N → deferred WITH the growing-N machinery to I9+ (named, not dropped); I8 ships the static-lattice strut + dynein/plectin coupling | analytic two-node Jacobian off-diagonal block; cross-response monotone in k_am; ⚠️ native gate confounded by IF cage → attribution via cross-block Frobenius-norm delta |
| **I9** | **Full-compartment native production** | EVERY compartment ON at physiological values (`--from-resting`, gravity, membrane, deformable nucleus, MT, cap, FA, ECM, **FSI-ON**, actin-MT, monomer field, **excluded-volume steric (I2b)**, resting Π₀ turgor). Then increment-2 deferrals below | combined native acceptance; **benchmark full-config wall-time BEFORE any conclusion** (the single biggest unmeasured number); **all magnitudes density-floored → report-not-tune** (γ/SF/cap/traction are FINDINGS, §6.2) |

**Increment-2 deferrals (I9+, PI-gated):** growing-N nucleation/capping/severing (true Tojkander turnover, pre-
allocated node pool + active mask — growing N breaks solver sparsity); full G-actin reaction
kinetics (nucleation/branch/cap/sever) on the monomer field (the conservative `∂(φc)/∂t` transport is I1c,
in-plan; the
reaction network richness is here); full ATP-state chemomechanical NMII cycle beyond I3's per-head stepping
(PI-gated after the I3 motor contract). (mixed q-p Biot
velocity form is NO LONGER deferred — it is planned I1b; **lamellipodium/filopodium/transverse-arc regions +
cell-type IF PROMOTED out of deferral to I4/I6/I7 by PI 2026-07-16 — see §3A.**)

**Growing-N reconciliation (§3A.b Arp2/3 dendritic nucleation vs the N-fixed increment-1 boundary §5-confirm-1):**
the I4 dendritic lamellipodium seeds a **pre-allocated daughter-filament pool + active mask** — branch KMC
ACTIVATES seeded-but-dormant daughters (N stays FIXED, topology emerges), NOT unbounded node allocation. True
growing-N nucleation/severing (Tojkander) stays I9+. This keeps §5-confirm-1 intact while §3A.b's dendritic
topology still emerges. The pre-allocated-pool + active-mask + KMC-event-loop machinery is named identically in
three docs (this plan, ENGINE_ARCHITECTURE, DYNAMIC_FEM roadmap item E) — **its concrete spec is authored at
first use (I4 dormant-daughter activation), not left abstract.**

**Explicitly DEFERRED compartments/mechanisms (named here, NOT silently dropped — mirrored in §5):** cell–cell
cadherin catch-bond junction (framework #9 — needs MULTICELL; single-cell exercises only the IF/desmosome
anchor-topology proxy §3A.c); membrane **Scriven–Boussinesq 2-D lipid surface-fluid** (framework 2-D-shell
row, RESEARCH tier); **spatial RVI/RVD** volume regulation; **kinesin anterograde transport + MT dynamic
instability** (I8 note); **organelles (ER/mito/Golgi)** as advected/steric compartments; **thymosin-β4 finite
monomer-pool buffering** (H10 — lands with the I1c conserved pool); full **multi-species transport** (ATP/ADP-
actin, profilin, capping, cofilin fields, ENGINE_ARCHITECTURE P4).

---

## 3A. First-validation scope expansion (PI 2026-07-16)

The PI ratified an expanded first-validation scope that **supersedes** the §3 deferral of
lamellipodium/filopodium and pulls cell-type IF into the foundation. Detail companion:
`_historical/FILAMENT_SUBSYSTEM_PLAN_2026-07-16.md` (the ③ solid-layer architecture facet). Three folded-in items:

**(a) Reconciliation.** Lamellipodium + filopodium + transverse-arc regions move from the I9+ deferral
into the **I4 region set**; cell-type IF moves into the **I7/first-validation** gate; the **I6 traction
gate now requires BOTH channels**. The I4/I6/I7 rows and the deferral list above are updated to match.

**(b) Lamellipodium REBUILD — dendritic Arp2/3 network (I4; a genuine build, not a port).**
Today's "lamellipodium" is `directed_front_growth_kernel`/`fiber_treadmill_kernel` — directed
barbed-end growth on CORTEX fibers, a ratchet body-force proxy with **no dendritic topology**. Rebuild
as a real branched network:
- **Topology:** mother filament + Arp2/3 side-branch → daughter filament; a mother/daughter tree seeded
  at the leading edge (a distinct `weave_cell` region on the membrane-adjacent manifold).
- **Branch junction = ANGLE-HARMONIC + thermal, NOT rigid 72°** (CLAUDE.md worked example). 3-body
  angle bond `U = ½k_θ(θ−θ₀)²`, θ₀≈70°, k_θ set so the thermal spread `σ_θ=√(kT/k_θ)` matches the
  observed branch-angle SD — the angle FLUCTUATES, it is not clamped.
- **Nucleation:** Arp2/3 autocatalytic dendritic nucleation at existing filament sides, rate
  **flux-limited by the ② G-actin monomer field** (cfd_transport `ENGINE_ARCHITECTURE_PLAN`) + an NPF
  (WASP/WAVE) activity field at the membrane; capping protein terminates growth (KMC via `hand_kmc`).
  ⚠️ **N-fixed encoding:** branch nucleation ACTIVATES pre-allocated dormant daughter filaments (active
  mask) — it does NOT allocate new nodes — so the dendritic topology emerges while N stays FIXED
  (§5-confirm-1; the growing-N reconciliation in the increment-2 deferral list). Requires I1c's monomer
  field LIVE before I4 (flux-limit consumer).
- **Force:** barbed ends push the membrane (Mogilner-Oster ratchet, already coded) → protrusion;
  retrograde flow is the EMERGENT reaction (already the physics). Nascent adhesions under the sheet are
  the lamellipodial clutch sites (the two-channel gate, §3A.a).
- **Replaces** the ratchet proxy as the front-protrusion mechanism.
- **Gate (rides I5 emergence):** dendritic net condenses from orientations sampled isotropically conditional
  on the declared leading-edge manifold/NPF field; branch-angle
  distribution matches θ₀±σ_θ; protrusion velocity + retrograde-flow band on the native cell.
- **Reuse:** `forces_warp` bending, `link_spring`, Mogilner-Oster ratchet (`motility_warp`),
  `hand_kmc` (branch/cap KMC), the ② monomer field (flux-limited nucleation).
- **New params (→ §4C):** θ₀ + k_θ (branch angle + thermal σ), Arp2/3 branch-nucleation rate, capping
  rate, NPF areal density.

**(c) Cell-type IF presets — keratin/vimentin (I7 follow-on, first-validation).**
IF returns as a **secondary passive network** (NOT the nucleus driver — that stays the I7 actin cap),
but **cell-type-differentiated** so the MCF7↔MDA contrast is exercised from the first native run.
`cell_type.py` gets per-type IF presets:
- **MCF7 (epithelial) = keratin:** anchor topology desmosome-like (cell–cell) + hemidesmosome;
  denser/stiffer; less extensible; tissue-continuum.
- **MDA-MB-231 (mesenchymal) = vimentin:** anchor perinuclear→FA/periphery; more extensible/dynamic;
  strain-stiffening; cell-autonomous; perinuclear cage protects the nucleus in confined invasion.
- **Encoding:** IF = passive network with cell-type-selected (i) anchor set and (ii) nonlinear
  strain-stiffening spring params. **Not the EMT switch** — IF is a cell-type INPUT (seed), a downstream
  marker + mechanical effector; the causal EMT program (Snail/Zeb/Twist, E→N cadherin) is out of scope.
- **⚠️ Honest single-cell scope:** in ONE cell the desmosome (cell–cell) keratin anchoring has no
  partner, so the single-cell-exercisable epithelial/mesenchymal difference is (i) IF
  extensibility/strain-stiffening + (ii) perinuclear-vs-peripheral anchor + (iii) nuclear protection.
  The full desmosome-continuum epithelial integrity needs MULTICELL (defer, not first-validation).
- **Gate:** combined cap+IF nucleus baseline (double-load-path, I7); MCF7 vs MDA nuclear deformability /
  confined-invasion contrast emerges from the IF preset (report-not-tune).
- **New params (→ §4C):** vimentin/keratin persistence length, strain-stiffening onset, IF areal
  density, `k_anchor` (desmosome vs FA/perinuclear) — PI-authored KnowledgeClaim, do NOT tune.

**(d) Filopodium — fascin tight bundle + tip adhesion (I4 region + I6 tip clutch).**
The filopodium was PROMOTED to first-validation (§3A.a) but — unlike lamellipodium (b) and IF (c) — had
NO param/gate spec (audit gap 2026-07-16). Fixed here so all three promoted regions are symmetric:
- **Topology:** parallel actin bundle cross-linked by **fascin** (tight, low-angle) protruding past the
  leading edge; seeded as a distinct `weave_cell` region on the membrane-adjacent manifold.
- **Mechanics:** bend/buckle under compression (Euler `F_crit≈π²κ/L²`) + tip-membrane push
  (Mogilner-Oster ratchet, already coded); the **tip adhesion is the low-load/nascent end of the ONE
  unified clutch** (I6; tip/shaft/basal EMERGENT, not hard-set classes).
- **Function:** durotaxis / guidance / confined-invasion sensing (MDA) — the tip mechanosensor.
- **Gate (rides I5/I6):** a tight low-angle fascin bundle CONDENSES from an isotropic seed (not scripted);
  tip clutch present as the nascent traction class; protrusion + tip-sensing band on the native cell
  (magnitudes density-floored → report-not-tune).
- **Reuse:** `forces_warp` bending, `link_spring` (fascin), Mogilner-Oster ratchet (`motility_warp`),
  the ONE catch-slip clutch (I6).
- **New params (→ §4C):** fascin bundle stiffness + inter-filament spacing/bundling angle, filaments-per-
  bundle, tip-nucleator density — PI-authored, do NOT tune.

---

## 3B. Phasing crosswalk — master increments I0–I9 ↔ facet phases P0–P4

⚠️ **Symbol collision (call it out so no fresh session conflates them):** this doc's **P1–P5 are the five
PILLARS** (§1); the `ENGINE_ARCHITECTURE_PLAN`'s **P0–P4 are the fluid/turnover PHASES**. They are NOT the
same "P". The master **I-spine (I0–I9) is authoritative for build order**; the facet phases map onto it:

| ENGINE_ARCHITECTURE phase | ≈ master increment(s) | Note |
|---|---|---|
| P0 infra + backward-compat (structured field grid + separate neighbour HashGrid, shared coordinates) | I1a scaffold | conservative field grid is not Warp `HashGrid`; shared regression gate |
| P1 ① fluid lands + retire 6πηR→γ_solid (SAME commit) | **I1a + I1b** | governing equations + moving-domain BC + outer-physical/inner-mechanical scheduler are bound in I0-A |
| P2 ② conservative RAD monomer on solved `v_f` | **I1c** (+ consumption AND release) | no prescribed-velocity production/scaffold branch; mass-conservation + FRAP gates |
| P3 full two-way ①↔②↔③ + filament–filament EV | **I1b2 + I2b + I9** | EV = I2b; optional momentum-drag pair is PI-KU gated |
| P4 multi-species + regulation (RESEARCH / PI-gated) | **increment-2 deferrals** | ATP/ADP-actin, profilin, capping, cofilin, Scriven membrane fluid, spatial RVI/RVD |

`FILAMENT_SUBSYSTEM_PLAN` uses the SAME I-spine (its G1–G11 map to I3/I4/I6/I7/I8 + §3A). **This crosswalk is
the single authoritative reconciliation of the three facets — keep them in sync THROUGH it, not independently.**

---

## 4. Ratification gates — foundation decisions and per-increment parameters

I0 is deliberately split. **I0-A is closed by the PI receipt above** and its contract now blocks any
contradictory implementation. **I0-Bn blocks only the first native run of increment n**; an unknown late-stage
parameter must not stop an earlier analytic implementation. A PI provisional value is allowed only when it is
explicitly labeled provisional, has an owner and expiry/review trigger, and is never adjusted to make a gate pass.

### A. I0-A foundation decisions — PI-RATIFIED 2026-07-16

| Decision | Ratified contract | Enforcement record |
|---|---|---|
| **Runtime + device policy** | **Warp CUDA GPU only.** HOOMD is never imported or executed for any purpose. Python may configure/orchestrate and postprocess, but every simulation kernel, state mutation, neighbor query, PDE step, KMC event, and mechanical iteration is GPU-resident Warp with no authoritative per-step host state | import denylist; runtime dependency denylist; zero GPU↔CPU roundtrip profiler gate inside the physical loop |
| **Physical-time scheduler** | one outer physical-time step with a frozen-state inner mechanical solve. Both explicit and implicit inner algorithms are allowed only as Warp-GPU implementations that converge the same registered residual | `Δt_outer`, residual, max-iteration rejection/retry policy; inner iterations are not physical time |
| **Native cortical density** | **70,686 active F-actin filaments** at the MCF7 reference cortex (`100 µm⁻² × 4π(7.5 µm)²`, rounded). The historical 38,000 and global ×40 scale are not production values | reference area/density/count derivation; unique-ID population ledger below |
| **First authoritative baseline** | **MCF7 + collagen ECM + α2β1–collagen clutch**. Collagen subtype/concentration, ligand density, resting geometry, and physiological BC magnitudes close through their relevant I0-Bn evidence records | baseline manifest uses one named ligand law; no α5β1–FN constants relabeled as collagen values |
| **NMII topology + FV** | explicit backbone+individual-head Stam-Hocky minifilament with **Hill FV**, completely reimplemented in Warp. Aggregate linear kernels are diagnostics only | topology/count, Hill, Bell, work-sign, CUDA residency gates at I3 |

#### A.1 Native population and GPU-memory budget — cortex is not the whole cell

Yes: `70,686` is **cortical F-actin only**, not total cellular filaments. The total must not be a second
hand-entered constant. Every compartment derives its active population from a physiological density/count and
its live reference geometry, then registers unique global filament IDs:

| Population | Production count rule | Current status |
|---|---|---|
| **Cortical F-actin** | fixed I0-A reference count `70,686`; current `L=3 µm`, `ℓ=0.5 µm` discretization implies `70,686×7=494,802` actin nodes | ratified; density provenance still carries the KB-3.18 internal-metadata correction item |
| **Lamellipodial actin** | `N_active=∫_A_leading ρ_barbed(x)dA`; daughter branches share the conserved monomer pool. A dormant daughter pool is memory capacity, not active biology | draft target is O(100 barbed ends µm⁻²); MCF7 leading-edge area and verified density close at I0-B4. The current hard-coded `200` is non-production |
| **Ventral/dorsal SF + arcs + cap** | unique filament IDs across the unified network; for separately nucleated bundles, `Σ_b N_fil,b` from sourced bundle density and bundle count | current order candidates (`10–30` filaments/bundle; historical `20 bundles×20≈400`) are planning estimates, not a ratified total; cap/arc density remains I0-B4/I0-B7 |
| **Filopodia** | `Σ_b N_fil,b`, with per-bundle fascin packing and seeded nucleator count sourced separately | current `20 filaments/bundle` is an order candidate; number of bundles is open at I0-B4/I0-B6 |
| **MT and IF** | separate polymer-system counts from cell-type-specific physiological evidence | MT candidate `20–250`; IF density/topology unresolved. Neither is folded into the actin count |
| **Nuclear lamina + chromatin** | live envelope nodes/edges/faces plus explicit chromatin-polymer IDs from sourced nuclear geometry and density/resolution | I2/I0-B2 must close both populations; neither may be hidden in a radial spring or a bulk modulus alone |
| **Collagen ECM fibers + cross-links** | `N_fiber=(∫_ΩECM ρ_length dV)/L_fiber` (or the equivalent measured network-density law); cross-links derive from physical intersection/binding density | collagen subtype/concentration, matrix domain, fiber-length distribution, and cross-link density close at I0-B1/I0-B6. ECM is a separate potentially dominant GPU population |
| **NMII/crosslinkers/clutches** | explicit particle/bond populations derived from areal/linear density and the selected topology | not "filaments," but must be included in GPU memory; head count and densities close at I0-B3/I0-B6 |
| **Membrane/nuclear meshes + fluid/reaction fields** | live mesh vertices/edges/faces and conservative cut-cell/PDE state from a resolution-convergence contract | not filaments, but included in `N_state`, scratch high-water, and full-device residency accounting |

For scale only—not as a biological choice—the MCF7 reference sphere has `A_cortex≈706.86 µm²`. At the draft
`ρ_barbed=100 µm⁻²`, a leading patch occupying `f={1,5,10,20}%` of that area would budget approximately
`{707, 3,534, 7,069, 14,137}` active barbed-end filaments before unique-ID overlap with cortex/SF/arcs is
removed. Therefore additional **cellular filament IDs are plausibly thousands to tens of thousands**, while
the full GPU state can grow much more because ECM fibers, per-filament nodes, every motor head/clutch, meshes,
fields, and scratch arrays are separate allocations. This sensitivity envelope is not a ratified count.

The required ledger reports, per compartment and whole cell:
```
N_unique_active = cardinality(union(global_filament_id))       # no cortex/SF/arc double-count
N_allocated      = N_unique_active + N_dormant_capacity         # actual GPU allocation
N_nodes          = Σ_i [ceil(L_i / ℓ_i) + 1]
N_state          = nodes + backbone beads + motor heads + crosslink/clutch states + field cells
GPU_bytes        = exact sum of allocated array bytes, including neighbor/PDE scratch high-water mark
```

Every I0-Bn closeout must publish the incremental and cumulative values plus measured peak GPU memory and
wall time. If the native population does not fit, optimize layouts/kernels, stream non-authoritative output,
or move to a larger/multi-GPU design; **never lower a biological density to fit memory**. Coarse counts may be
used only for explicitly labeled arithmetic smoke tests and can never support a physics conclusion.

#### A.2 HOOMD source quarantine and executable-scope gate

A static 2026-07-16 inventory found **210 Python files** under `aleph/` + `native/` that directly import
HOOMD, including **34 test files**. They are migration input, not an allowed executable fallback. Until they
are physically moved, ported, or deleted, they remain quarantined by policy and **must not be collected by
the new-engine test command**. I1 must land all of the following before any simulation gate is claimed:

- new runtime and tests live under `aleph/ac/` and `aleph/tests/ac/` respectively;
- the `ac/` dependency graph and generated extension link graph contain no `hoomd` symbol/package;
- CI fails on any `import hoomd`, `from hoomd`, HOOMD shared-library link, or HOOMD subprocess in the active
  runtime/test scope;
- the default new-engine test target collects only Warp tests; legacy test discovery is disabled, and no
  historical fixture is regenerated by executing HOOMD;
- a frozen historical numeric fixture may be read as data only after provenance is labeled; acceptance is
  analytic/experimental and Warp self-consistency, never live HOOMD parity.

The inventory count is an explicit migration burn-down metric. It does **not** authorize executing any of
those 210 files during the port.

### B. I0-Bn cross-cutting parameter conflicts — do not choose by target matching

| Name | Conflict / evidence status | Decision needed before |
|---|---|---|
| **F_stall_head** | 0.5 pN attribution is not yet represented by an audited SourceEvidence row; 2.0 pN/Billington support is also not yet closed for this exact head-level use | I3: source-audit both claims and select by motor isoform/assay, never to lift γ |
| **N_side** | 10 (AFINES model) vs 29–30 (Billington structural interpretation) | I3: bind explicit heads per side to the selected NMII isoform/topology |
| **v0** | 0.12 vs 0.2 µm/s claims differ by assay/model context | I3: bind assay temperature, load definition, and isoform with the value |
| **Biot modulus/storage** | `M=300 Pa` copied from `K_drained` risks double-counting; `1000 Pa` as a bare difference also omits α | I1a: derive the chosen closure, e.g. `K_u=K_d+α²M`, and record `S=1/M` with units |
| **Integrin force scales** | the recorded ~7 pN catch scale and ~30 pN slip scale are not competing universal `F*` values and are not, by themselves, proof of a ligand mismatch | I6: preserve distinct parameters in the selected ligand-specific law; α2β1–collagen needs its own audited evidence rather than relabeling α5β1–FN constants |

### C. Genuinely unknown or not-yet-ratified — provisional PI-authored only (magnitude gates INVALID until supplied)
**`permeability k`** (Darcy k for the I1b velocity form; cytoplasm ~10⁻¹⁷–10⁻¹⁶ m² order, sets streaming speed +
τ_p = µR²/kM — surface a KnowledgeClaim, do NOT tune to a flow band) · **`µ_pore`** (pore-fluid viscosity for
Darcy, ~1-2× water — NOT the 65.9 Pa·s effective bulk drag; keep the two separate) ·
`k_xb` (⚠️ MASTER force knob, mislabeled compliance; 1 pN/µm breaks the stall mechanism — physical ~100-1000) ·
`r0_head` (200 nm unphysical, ~10× real) · `L_bb/n_bb` (700nm/14 vs Billington 300nm) · `k_on` (50/s → 0.99
bound, do NOT tune) · **`k_linc`** (nesprin nonlinear; 8 pN is a TENSION not a stiffness — flatten gate invalid
until sourced) · `k_am` (re-anchor to plectin/MACF single-molecule, NOT to the IF-failure mode) · `k_mtoc_nuc`
· `n_cap_fil` · `N_FA_sites`/FA density · `k_linc` LINC areal density · dynein `N_dyn/k_on_dyn/f0_dyn` ·
`N_minifil_per_SF_cross_section` (⚠️ do NOT back-solve to hit 10-30 nN — geometrically impossible) ·
`sarcomere_um` (0.5-1.4) · substrate E + patch radius a · nuclear oblate aspect band (1.5-3) · `VIN_KDISS` (no
KB source) · **Arp2/3 branch** `θ₀`≈70° + `k_θ` (thermal σ_θ — do NOT clamp) + branch-nucleation rate +
capping rate + NPF areal density (lamellipodium rebuild, §3A.b) · **cell-type IF** vimentin/keratin persistence
length + strain-stiffening onset + IF areal density + `k_anchor` (desmosome vs FA/perinuclear) (§3A.c — do NOT
tune; the MCF7↔MDA contrast must EMERGE from the sourced preset).

**② monomer-transport + fluid unknowns (I1c/I1b need these BEFORE their native runs — audit 2026-07-16):**
**`D_c`** (an existing draft claim/evidence trail points to crowded-cytoplasm values around 2–6 µm²/s;
ratification and MCF7/assay applicability remain open — this is not a blank KB gap, and free-solution
50–80 µm²/s must not be substituted silently) ·
**`c₀`** (initial free-monomer pool — total actin ~100 µM, thymosin-β4-buffered free ~tens µM; **MCF7-specific,
PI GAP**) · **`α-Biot`** (fluid–solid volumetric coupling coefficient) · **thymosin-β4 buffer** capacity +
exchange rate (finite-pool growth control, H10 — was absent entirely) · **`(η_f/k)` relative-velocity fiber-drag
coefficient** (I1b2 velocity-difference momentum drag — NO KB claim, PI-authored KU) · **resting osmotic turgor
`Π₀`** (MCF7 physiological baseline, pin as the PRE-STRESS value — NOT zero; do NOT tune) · **nucleus depth
(I2):** lamin-A/C vs lamin-B moduli + strain-stiffening knee + nucleoplasm viscosity + **envelope-rupture
tension threshold** (framework #6; E_nuc 1–10 kPa method-dependent) · **filopodium fascin** bundle stiffness +
inter-filament spacing/bundling angle + filaments-per-bundle + tip-nucleator density (§3A.d) ·
**excluded-volume `k_EV`** (I2b steric pair-potential stiffness — a numerical repulsion scale; assert
grid-invariant per the Magic-Number Block, do not tune to a crowding outcome).

### D. Modeling-choice knobs — flagged controls, NEVER tuning knobs
angle-dependent crosslinker lifetime (angle-dependence is a CHOICE; emergence must pass with it OFF) · α_hub ·
inner_steps/p_clip=5000 (assert non-binding or it silently affects results).

---

## 5. Emergent-vs-seeded boundary — part of I0-A ratification
- **SEEDED (physiological BC/inputs, not tuned):** FA candidate-site density/distribution, nucleation
  loci+nucleator identity, leading-edge/NPF polarity field, nucleus position+LINC sites, and each region's
  initial geometric manifold. Local actin orientations are sampled isotropically conditional on those inputs;
  any literature-mandated orientation bias is separately declared and ablated. Also seeded: **turgor/osmotic BC at the pinned
  physiological Π₀ (§4C), NOT zero**, Biot ON at physiological value, MT aster geometry, minifilament
  count/positions at the ratified density.
- **EMERGENT (from force balance + KMC, measured without construction labels):** alignment→bundles (all
  SF/arcs/cap/filopodia), engaged myosin/dynein fraction (Bell self-limiting), which candidate regions
  condense while cortex remains isotropic, traction magnitude and adhesion class, cap geometry, nucleus
  flatten, MTOC position, remodeling and durotaxis.
- **DEFERRED (named, out of the single-cell FOUNDATION scope — NOT silently dropped; homes in the §3
  increment-2 list):** cell–cell cadherin catch-bond junction (framework #9, MULTICELL); membrane
  Scriven–Boussinesq 2-D surface-fluid (framework 2-D-shell row, RESEARCH); spatial RVI/RVD volume regulation;
  kinesin anterograde transport + MT dynamic instability (I8 note); organelles (ER/mito/Golgi) as
  advected/steric compartments; membrane K_A area-upturn + reservoir/caveolae unfolding (default-plateau, §6.8);
  full multi-species transport + thymosin finite-pool (with I1c conserved pool).
- **Confirm:** (1) increment-1 holds allocated N FIXED (dormant activation is allowed; unbounded allocation is
  not); (2) seed orientation is isotropic conditional on declared manifold/nucleator fields, not secretly
  pre-aligned to the measured FA/tension axis; (3) `weave(bundle)` DEMOTED to a non-authoritative control;
  (4) Bell duty ≠ the 0.65 geometric
  overlap (keep separate); (5) flat-S → SEEDED-labeled scaffold is the sanctioned fallback; (6) nucleus flatten
  measured at the RESTING adherent baseline.

---

## 6. Honest hard truths + fallbacks (the new session MUST internalize)
1. **Emergence may not condense bundles** in the overdamped quasi-static solve. By I5 the live I1c monomer
   field, topology-reforming crosslink KMC, FA anchoring and myosin-driven alignment are present, but the
   pre-allocated dormant-pool encoding is not unrestricted growing-N and the
   **myosin-turnover→local-density-increase trigger** identified as the STRESS_FIBER_TARGET primary driver is
   still absent. Therefore the emergence claim is genuinely falsifiable, not guaranteed by the architecture.
   Fallback: a SEEDED scaffold validates force mechanics only and must not be reported as self-organization.
2. **Force magnitudes are density-floored, NOT the deliverable:** γ inherits the ~530× active floor; SF active
   O(0.1-0.9 nN) vs 5-6 nN target; cap flattens ~0.2 µm (corrected 4π arithmetic) vs ~2.5 µm target past the
   lamin knee. These are FINDINGS, not knobs. Closing them needs growing-N + MCF7 density data that doesn't
   exist — surface to PI. **Never add heads/density/N_parallel to close a floor.**
3. **Master-knob stiffnesses genuinely unknown** — each magnitude gate is HELD until its I0-Bn record closes.
   Build the mechanism + analytic
   gates (capstan, Jacobian, Newton) now; hold magnitude verdicts.
4. **CFD's "emergent 65.9 Pa·s viscosity" is BROKEN and dropped** — a dimensional estimate gives 340-1100 Pa·s
   (5-17× over). Honest CFD deliverable = rate-dependence + spatial τ_p transients, NOT a changed resting
   baseline, NOT an emergent viscosity. Zero-mean may under-power undrained stiffening (Skempton) — real
   architecture tension, surface don't paper over.
5. **Full-config native wall-time is UNMEASURED** (no committed FSI-ON native artifact exists) and may exceed
   budget. Every increment ends with measure-before-default; fallback = component-wise native gates + a
   reduced-outer-step combined run labeled non-authoritative for magnitude.
6. **6 double-count instances** to guard (I2-crosslinker-relaxation, motor+η, f_myo+fine-motor, M=K_drained,
   LINC+soft_contact, cap+IF) — see the workflow plan §4 R1 table. **Each MUST carry a specific guarding gate,
   not just a mention:** motor+η at I1 (retire 6πηR→γ_solid, same commit); f_myo+fine-motor at I3 (REMOVE
   replaced `cortex.myo_i`); Biot storage at I0-B1 (derive `K_u=K_d+α²M`, never set `M=K_drained`);
   crosslinker-relaxation at
   I4 (pick ONE channel); **LINC+soft_contact at I7** (land LINC wire-in + the soft-contact/plate term TOGETHER,
   per FSI_ALL_COMPARTMENT); **cap+IF at the combined I7 baseline**. Do not close an increment whose double-count
   guard has no gate.
7. **Discrete force/work consistency is a first-class GATE, not an afterthought.** Pressure/fluid-content and
   immersed-solid traction are thermodynamic conjugates; a pressure diffusion solve alone is not a fluid
   momentum integrator. Gate I1a/I1b with adjoint spread/interpolation, pressure-work sign, global
   fluid-content balance, and a **net-force-projection diagnostic** (Σf over immersed nodes → closed-system COM
   drift ≈ 0, logged). I1b2, if authorized, adds an explicit equal-and-opposite relative-velocity drag pair;
   that discrete action-reaction property does not imply that fluid inertia has been introduced. Do not claim
   momentum conservation beyond the equations actually solved. This also underwrites the H4/H5
   "no self-translation from internal flow" tests.
8. **The membrane runs without its DEFINING area mechanic by default.** K_A area-stiffness upturn +
   reservoir/caveolae unfolding are default-OFF (buffered plateau) to protect dt (K_A=2.35e5 → dt collapse if
   the reservoir exhausts); MCF7 is caveolae-deficient (`RESERVOIR_STRAIN` unsourced). At I9 the membrane is ON
   but **flag that the large-deformation area response is deferred** — do not read blebs / area response PAST
   the plateau as physical until K_A lands (or an implicit membrane sub-step). Helfrich bending is landed.
9. **MT is a static-lattice strut only.** Kinesin/dynein anterograde-vs-retrograde transport + dynamic
   instability (framework #3) are NOT in I8 (the dynein POPULATION is; the rest is growing/shrinking-N → I9+).
   Do NOT attribute organelle transport to fluid advection (I1c wording fixed) — it is motor-on-MT, deferred.

---

## 7. First-session bootstrap (what the fresh session does, in order)
1. Read this doc + the companion docs + `AGENTS.md` + `CLAUDE.md` + the current Dashboard. `git log`, confirm
   branch; stop if an authority document contradicts the ratified I0-A receipt.
2. **Verify the 2026-07-16 I0-A receipt is mirrored everywhere:** Warp-CUDA only, no HOOMD execution,
   head-resolved NMII+Hill, 70,686 cortical filaments, MCF7×collagen×α2β1. Generate the native-population/
   GPU-memory ledger before allocating the full cell.
3. Before each increment n, close **I0-Bn** only for that increment's parameters. Record evidence state and
   provisional policy; do not let an I8 unknown block an I1 analytic gate.
4. Scaffold `aleph/ac/` as a Warp-CUDA-only package. Archived HOOMD source may inform a port but must not be
   imported, linked, executed, benchmarked, or used as a parity runtime.
5. Build **I1 (fluid-first substrate)** first — it is the foundation everything else attaches to. Analytic
   gate (Terzaghi/Green/parity) before any compartment rides it. Seed resting Π₀ at its physiological value
   (§4C) from t0 — never a floppy zero-turgor bag.
6. Execute the complete spine:
   `I1a → I1b → (I1b2 only if PI-authorized) → I1c → I2 → I2b → I3 → I4 → I5 → I6 → I7 → I8 → I9`.
   Use the §3B crosswalk to keep the three facets in sync. Close each increment
   with its analytic-first gate → native `--from-resting` full-compartment gate → interactive HTML viz
   (browser-verified) → commit. Retire the corresponding lumped FF mechanism as each fine-grained replacement
   validates (per [[project-active-cell-c-redesign]]).

## Change log
- 2026-07-16: created as the new-engine handoff. Synthesizes the 18-agent option-C design workflow
  (wf_77ea823d-c45) + re-audit + PI's new-engine/fluid-first/deformable-nucleus decisions. All 6 workflow
  subsystems came back CONCERNS (none REWORK) — the mechanisms are buildable; the honest caveats above are the
  load-bearing part. To be executed increment-by-increment in a fresh session with native gates.
- 2026-07-16 (scope expansion §3A): folded in the filament-subsystem plan
  (`_historical/FILAMENT_SUBSYSTEM_PLAN_2026-07-16.md`) + PI decisions. Lamellipodium REBUILD (dendritic Arp2/3,
  angle-harmonic branch) + filopodium/tip + transverse-arc PROMOTED from the I9+ deferral to the I4 region
  set; cell-type keratin/vimentin IF → first-validation (I7 follow-on, MCF7↔MDA contrast); I6 traction gate
  now requires BOTH channels (lamellipodial + SF). The same-day 7-vs-30 pN binary interpretation was later
  superseded by the implementation-readiness audit below. Updated P3, I4/I6/I7 rows, deferral list, §2b IF
  kill-row, and parameter register.
- 2026-07-16 (audit-driven completeness pass — 6-lens × adversarial-verify workflow, 26/40 findings survived):
  **(1)** fixed the duplicated I2 row; **(2)** added **I2b filament–filament excluded volume** (was absent
  from the whole I-spine though CLAUDE.md mandates "LJ EV ON from Phase 1" and ENGINE_ARCH lists it as a
  defining living-cortex gap) → §2 ADD-row, I9 all-ON, parameter register `k_EV`; **(3)** deepened **I2 nucleus** to
  framework-#6 (lamin-A/C vs -B, nucleoplasm viscosity, envelope-rupture) + reconciled the I1a no-flux mask
  with the deformable oblate surface; **(4)** split **I1b/I1b2** (pore-fluid advection vs the PI-KU
  relative-velocity fiber drag
  `(η_f/k)`); **(5)** pinned resting **Π₀** + membrane **Lp water-flux** BC into I1a/parameter register/§5
  (compartment #7 was
  demoted to a static mean channel); **(6)** gave **filopodium** (§3A.d) its missing param+gate — symmetric
  with lamellipodium/IF; **(7)** reconciled §3A.b **Arp2/3 dendritic nucleation** with the N-fixed boundary
  (dormant-daughter activation, not node allocation); **(8)** added the **§3B P0–P4 ↔ I0–I9 crosswalk** (+ the
  P-symbol-collision callout); **(9)** added `D_c`/`c₀`/`α-Biot`/thymosin + FRAP(I1c)/TFM(I6) oracle gates;
  **(10)** §5 explicit **DEFERRED-compartment list** (cell–cell junction #9, Scriven membrane fluid, RVI/RVD,
  kinesin+MT-dynamic-instability, organelles, K_A/reservoir) — named, not silently dropped; **(11)** §6 new
  hard-truths #7 momentum-conservation gate, #8 membrane K_A default-off, #9 MT static-strut; expanded #1
  (myosin-turnover SF trigger) + #6 (per-instance double-count guard gates). Every "scary" item the audit
  first flagged (turnover-timing, 6πηR→γ_solid, membrane Lp flux, I5 sequencing, NMII FV-law, numerics,
  wall-time) verified as ALREADY handled in the folded facets — no change needed there.
- 2026-07-16 (implementation-readiness audit): corrected the Biot/Darcy state variables and conservative
  transport law; replaced the fixed-box held-face stencil with a moving-domain conservative I1a contract;
  separated the structured field grid from the neighbour HashGrid; made the outer physical-time/inner
  mechanical scheduler explicit; replaced aggregate NMII with a head-resolved production topology; split I0
  into foundation and per-increment gates; corrected the integrin force-scale interpretation and D_c evidence
  status; reconciled seeded inputs with label-blind emergence; restored the full I1a→I9 bootstrap order.
- 2026-07-16 (I0-A PI ratification): fixed the runtime to Warp CUDA GPU only with zero HOOMD execution;
  authorized Warp explicit/implicit inner solvers under one outer physical clock; bound the production motor
  to a head-resolved Stam-Hocky topology with Hill FV; selected MCF7×collagen×α2β1 and 70,686 cortical F-actin
  filaments. Added a unique-active/allocated/nodes/state/GPU-bytes population ledger so lamellipodium, SF,
  filopodia, MT, IF, and explicit motor heads increase the total transparently without double-counting or
  memory-driven biological down-scaling.

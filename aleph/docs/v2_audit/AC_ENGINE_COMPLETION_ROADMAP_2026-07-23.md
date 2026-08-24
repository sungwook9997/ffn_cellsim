# AC engine completion roadmap — every compartment to NATIVE (2026-07-23)

> ✅ **DIRECTION RATIFIED 2026-07-23 → see `_historical/WHOLE_CELL_EVENT_RUNTIME_2026-07-23.md` (boot from it).** The path to
> NATIVE is now the **Whole-Cell State-Conditioned Event Runtime** (CellState = context only, never force;
> component/connector-owned explicit events; GATE A static baseline first, GATE B accepted-step transaction as
> shared infra). Cortex is the FIRST vertical slice: the landed `density_per_fil=20` (`bc5ff3b0`) + overlap_free
> default (`81805535`) are a **STATIC CONTROL** (native-verified single-spanning 0.2µm 3D shell), NOT the final
> dynamic cortex. Immediate work = cortex GATE A (solver re-tune → resting bound-myosin → converged static
> baseline) + spec the whole-cell common contracts.

Honest state after the overnight run: the composed world is REGISTERED + dispatched and kernel-binding is
CUDA_UNIT-verified on the A5000 (`tests/ac/engine` 256/0), but **no compartment has reached NATIVE** (full
physiological population + operating point under one physical clock). This roadmap is the per-compartment path
to NATIVE, split by what is PI-GAP-blocked vs advanceable now, mapped to `ROLLING_ROADMAP.md` R-waves.

Evidence ladder: `CONTRACTED → SEAMED → KERNEL_BOUND → CUDA_UNIT → CONNECTED → NATIVE → OPTIMISED → PRODUCTION`.
"Complete" = **NATIVE** (physiological population + t0 pass with exact population/byte ledgers, no host
authority), then the whole-cell **R8** baseline. "완성" is R8 NATIVE for all 13 + 32.

## 0. The master gate — everything native waits on the resting baseline

A native run can only start from a valid force-balanced physiological t0, supplied by the resting bound-myosin
setpoint (mechanism BUILT + demonstrated `9bb91e0a`, 40.74→2.17 pN).

> ⚠️ **REFRAMED 2026-07-23, then CORRECTED post-Codex+native (see `_historical/RESTING_SETPOINT_SOURCING_2026-07-23.md`
> §Corrections).** The "blocked on ONE PI-GAP = fraction + per-head force" framing UNDERSTATES it — closing the gate
> needs the whole **force-magnitude/density SET** (heads-per-minifilament N_side=10 ~3× under Billington/Nagy;
> per-head force; load-dependent duty; minifilament density). BUT the intermediate "physiological content is
> sufficient in total (1.2–4.8×)" claim is **RETRACTED**: it used a scalar-capacity sum that over-counts the
> method-of-planes cut force by ~26.6× (reproduced), so sufficiency is UNDETERMINED and leans deficient — native
> decides. Three native/verified corrections: (1) the myosin target is likely **higher** than 140 pN/µm — the repo
> carries a sourced MCF7 tension 180–400 pN/µm (Hosseini 2020) × active-fraction 0.70 ⇒ ~126–280; (2) at native,
> production capture 0.05 µm caps the **seedable** duty at **0.17** (heads sit ~0.2 µm off actin) and the seeded
> step **rolls back** (native run 2026-07-23); (3) `params_i0b3.yaml` already records NM2B-pure baseline + the
> k_off0 ADP-release mislabel. **Gate stays OPEN; nothing closed fast (PI directive 2026-07-23).**

So the roadmap is: **Phase A (PI ratifies the physiological setpoint SET) → Phase B (advanceable-now authoring, to
CUDA_UNIT) → Phase C (post-resting: CONNECTED→NATIVE chain, serial on the one A5000).**

## 1. Per-compartment status → path to NATIVE

| # | Component | Now (honest) | To reach NATIVE | PI-GAP? | Advanceable now (pre-resting)? |
|---|---|---|---|---|---|
| 1 | cortex (SURFACE_BODY) | ac/cell physics live; engine seam CUDA_UNIT | resting t0 (bound-myosin) → native 70,686 CONNECTED | ⭐ resting bound-myosin | binding done; native blocked on resting |
| 2 | membrane (SURFACE_BODY) | Helfrich+ERM live; seam CUDA_UNIT | resting t0 + subdiv-6 native + physiological ERM density | ⭐ resting; MCF7 ERM density | binding done; native blocked |
| 3 | cytosol (FLUID_VOLUME) | Biot/Darcy in ac/fluid; fluid_core CUDA_UNIT | bind full field solver + moving boundaries → CONNECTED (mass/no-flux/adjoint-work) → native | — | **YES** — CONNECTED binding + gates (device) |
| 4 | nucleus (CORE_BODY) | lamina/chromatin live; Core seam; reduced-ROM SEAMED | native Core CONNECTED; then ROM Φ/J (design) as OPTIMISED | I0-B2: lamin, k_linc, η_nuc, rupture, aspect | **partly** — native Core path; ROM design |
| 5 | sf_arc (ACTIVE_LOAD_PATH) | SEPARATE component (ratified); load_path kernel KERNEL_BOUND | build SF DISJOINT population (ventral/dorsal, FA-anchored) + connectors → native | — | **YES** — SF population build (P1) |
| 6 | focal_adhesion (ACTIVE_LOAD_PATH) | clutch graph; α2β1-collagen series (load_path) | FA maturation + clutch native (with ECM) | clutch maturation params | partly (gated on ECM for the series) |
| 7 | **microtubule** (STRUCTURAL_RIG) | ⚠️ SEAMED — STATIC bead adapter | dynamic instability (grow/shrink/catastrophe/rescue) + MTOC rotation + native | 🔴 MT DI rates (MCF7) | author kernel now; **activation blocked on rates** |
| 8 | **intermediate_filament** (STRUCTURAL_RIG) | ⚠️ SEAMED — monolithic-Hookean | nonlinear WLC cable + turnover + native | 🔴 IF WLC EA/x_max | author kernel now; **activation blocked** |
| 9 | lamellipodium (ACTIVE_LOAD_PATH) | branch-angle term KERNEL_BOUND; rest SEAMED | full adaptive Arp2/3 branched net + brownian-ratchet CONTACT + nascent FA + native | Arp2/3 branch/cap rates | **partly** — Arp2/3 net (ac/weave exists) |
| 10 | **filopodium** (ACTIVE_LOAD_PATH) | ⚠️ SEAMED | bundled F-actin + fascin + tip ratchet + native | 🔴 filopodium k_fascin | author now; **activation blocked** |
| 11 | nmii (ACTIVE_LOAD_PATH) | backbone + MOTOR 2-array crossbridge KERNEL_BOUND | full per-head Bell KMC/turnover + native motor population + Hill/Bell native gate | resting bound fraction (⭐); NMII backbone L_p | **YES** — head KMC native (ac/motor mature) |
| 12 | **ecm** (ENVIRONMENT) | far-field boundary KERNEL_BOUND | GPU SoA topology + collagen constitutive + α2β1 clutch + native | 🔴 collagen gate (P2 ratifying); topology on unmerged branch | **YES** — merge codex/ecm-gpu-topology + SoA port (P5) |
| 13 | world_boundary (ENVIRONMENT) | KERNEL_BOUND far-field frame | native reference-frame constraint | — | **YES** — near done |

**Connectors (32):** most CONTRACTED/SEAMED. The kinetic joints (ERM, FA-clutch, LINC, MOTOR, plectin,
spectraplakin, transient-actin, contact) each need their `chemistry_card` bound + commit-on-accept native. The
α2β1-collagen series (load_path) is furthest along; the NMII MOTOR crossbridge is wired; ERM ties to the
resting gate; LINC ties to nucleus I0-B2.

## 2. Phase A — PI-GAPs to supply (each unblocks the marked compartments)

| GAP | Unblocks | Note |
|---|---|---|
| ⭐ **resting_bound_myosin_fraction + per-head force** | cortex, membrane, nmii, and ALL native gates | mechanism built; the master unblock |
| **MT dynamic-instability rates** (v_grow/v_shrink/catastrophe/rescue, MCF7) | microtubule | no defaults invented |
| **IF nonlinear WLC** (EA / x_max) | intermediate_filament | keratin/vimentin |
| **filopodium k_fascin** (+ tip params) | filopodium | |
| **collagen constitutive gate** (P2 ratification — near done, `P2_RATIFICATION_PACKAGE`) | ecm, focal_adhesion series | conc-resolved gate; PI sign-off pending |
| **MCF7 ERM linker density** | membrane/ERM density | currently 1/node diagnostic |
| **nucleus I0-B2** (lamin, k_linc, η_nuc, rupture, aspect) | nucleus | pre-existing GAP |
| **NMII mature-minifilament backbone L_p** (persistence length) | nmii bending stiffness (F6) | derived-not-fit but L_p unsourced (INTEGRATION_nmii_actuator.md) |

## 3. Phase B — advanceable NOW (author to CUDA_UNIT without PI, native-gated behind resting)

These can be built/verified to CUDA_UNIT on the A5000 in parallel (disjoint modules, adversarial-verify-gated),
without any PI-GAP — only their NATIVE gates wait for the resting baseline:

1. **cytosol/nucleus CONNECTED binding** — inject the real Biot field solver + moving membrane/nucleus
   boundaries; close mass/no-flux/adjoint-work in one coupled candidate. (R2/R5)
2. **SF population build** — sf_arc as a separate disjoint F-actin population + its FA/cortex/LINC connectors.
   (R4; P1)
3. **nmii head KMC** — full per-head Bell attach/detach/turnover on the native motor population. (R4)
4. **lamellipodium Arp2/3** — adaptive branched network from ac/weave nucleation. (R6)
5. **ecm GPU SoA topology** — merge `codex/ecm-gpu-topology` (+3, unmerged) + port ff/ecm_library to
   ECM-owned CUDA SoA (device schema, IDs, free-list, broadphase). Constitutive gate stays P2-gated. (R3; P5)
6. **MT/IF/filopodium kernels** — author the dynamic-instability / WLC-turnover / bundle kernels with the
   rate/stiffness constants read from an UNSET PI-GAP slot (mechanism ready, activates when PI supplies). (R5/R6)
7. **connector chemistry-card binding** — bind each kinetic connector's per-bond kinetics + adjoint scatter.

## 4. Phase C — post-resting-gate: the serial NATIVE chain (one A5000)

Once the resting gate closes (Phase A #1), the native gates run in dependency order on the single GPU (each
slice's native cell contains the prior's live state — strictly serial, per `AC_PARALLEL_SESSIONS`):

- **R2** membrane⊕cortex⊕cytosol native (the resting baseline IS this slice) → first native force-balanced cell.
- **R3** ecm native (SoA + constitutive) — gated on the collagen P2 gate.
- **R4** sf_arc⊕focal_adhesion⊕ecm load path + nmii MOTOR native (active tension emerges).
- **R5** microtubule⊕intermediate_filament⊕LINC⊕nucleus native (tensegrity; gated on MT/IF PI-GAPs).
- **R6** lamellipodium⊕filopodium native (protrusion; gated on filopodium PI-GAP).
- **R7** transported chemistry & turnover.
- **R8** whole-cell native baseline (all 13/32 under one clock + figures) = **완성**.
- R9–R15 perturbation matrix, performance, optimisation, scaling, validation, deployment.

## 5. Definition of "완성" (goal-clear condition)

Every one of the 13 components + 32 connectors at **NATIVE** (physiological population + t0, exact
population/byte ledgers, zero host authority), the whole-cell **R8** baseline passing with figures, and the
plan recording it. Current distance: the STRUCTURE is done; the resting gate + ~7 PI-GAPs + the per-compartment
mechanistic physics (esp. MT/IF/filopodium/ECM) + the serial native chain remain.

## 6. What to do next (dependency order)

1. **PI supplies the resting bound-myosin values** → activate the built mechanism → close the resting gate
   (the single highest-leverage action; unblocks R2 and the whole native chain).
2. In parallel (Phase B, no PI needed): cytosol/nucleus CONNECTED binding, SF population build, nmii head KMC,
   ECM SoA topology merge+port, MT/IF/filopodium kernel authoring (rates as unset PI-GAP slots).
3. PI supplies the remaining GAPs (MT DI, IF WLC, filopodium k_fascin, collagen gate ratification, ERM density,
   nucleus I0-B2) → activate those compartments.
4. Run the serial native chain R2→R8 on the A5000 → whole-cell NATIVE baseline = 완성.

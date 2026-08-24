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

# FF engine re-audit — actin architecture · actin-MT coupling · cytosol CFD (2026-07-16)

> **POST-AUDIT CORRECTION (implementation-readiness pass, 2026-07-16):** the original audit correctly found
> that FSI was absent from landed production runs, but overclassified `biot_fluid_warp.py` as a completed
> production solver. Its device arrays and Peskin transfer kernels are reusable; its fixed-cube held-face
> update is not conservative at a nominal no-flux boundary and does not represent the live
> membrane-minus-nucleus domain. The corrected target is master I1a/I1b, not "turn the existing FSI ON."
> PI subsequently ratified a Warp-CUDA-only runtime with zero HOOMD execution, head-resolved Hill NMII, and
> an MCF7×collagen×α2β1 baseline with 70,686 cortical F-actin plus separately derived non-cortical populations.

PI-authorized re-audit (3 grounded subagents, code-cited) triggered by the PI's stress-fiber + actin-filament
reference images and the cytosol-CFD question. **Through-line: the FF whole-cell model is a cortex-shell +
passive-MT + nucleus + IF-cage + membrane. The rich force-generating actin architecture (stress fibers,
lamellipodium, actin cap, arcs) and the active cytosol CFD EXIST as code but are NOT in the running model —
they are specs / tests / other-engines / halted / default-off.** Every PI instinct here was correct.

## 1. Actin architecture coverage (of the 10 reference structures)
Levels: A = spec/weave/viz/test-only · B = has a driver but HALTED / on retired-HOOMD or DCM engine · C =
integrated into `network_warp.py::simulate_whole_cell_compression_on_device` (the AFM native full-cell).

| Structure | Level | Note |
|---|---|---|
| **Cortex** | **C** | isotropic sphere (`gamma_floor.build_crosslinked_cortex`). ⚠ myosin = LUMPED constant `f_myo`, NOT Hill f–v |
| Lamellipodium | A (FF) / B (DCM+HOOMD) | Arp2/3 kernel exists but FF whole-cell `branch_triples` is EMPTY → dead. DCM has a real one |
| Filopodium | A / B (FF `motility_warp`) | polymerization-ratchet protrusion runs standalone, NOT in the AFM cell |
| **Ventral stress fiber** | A (weave) / B — **HALTED** | driver = `h7_stress_fibers_activation_gate.py` on **archived HOOMD**; PASSIVE backbone only, NMII deferred; Kumar 10–30 nN active gate DEFERRED; force-scale HALTED→PI |
| Dorsal stress fiber | **ABSENT** | 0 hits |
| Transverse arc | **ABSENT** | 0 hits |
| **Perinuclear actin cap** | **ABSENT (as actin)** | the nucleus↔cortex structure that WAS built is an INTERMEDIATE-FILAMENT cage (`intermediate_filaments.py`), not an actomyosin cap + LINC |
| Microvillus | A | weave spec; bundler un-sourced (filamin stand-in), PI-gated |
| Adhesion belt | ABSENT | — |
| Contractile ring | ABSENT | — |

**Net: 1 of 10 integrated (cortex).** The "unified `weave()`/`architecture_spec`" machine genuinely unifies 5
structures in ONE builder — but it is a design-validation PROTOTYPE, never called in any physics run; the
cortex's actual physics home is a SEPARATE ad-hoc builder (`build_crosslinked_cortex`). FA/traction was
studied via imposed retrograde flow (`ff_clutch_traction`) or the CORTEX as substrate (`h7_manifold_traction`)
— NOT contractile ventral SF. So "FA without functional stress fibers" is accurate.

## 2. Actin–MT coupling — essentially ABSENT
- The ONLY actin↔MT link is a **one-sided repulsive tip↔cortex excluded-volume contact** (`soft_contact_kernel`,
  engages within 0.5 µm). No crosslinker along the MT length (+TIP/MAP/plectin/spectraplakin/dynein = none).
- The merged cortex+MT stiffness is **BLOCK-DIAGONAL** — no bending triple spans the two domains → actin and MT
  elastic responses are mathematically INDEPENDENT.
- MT is a **passive strut** (bending EI only). No dynamic instability (deferred/OFF), no motor transport, no
  polymerization in the whole-cell. **No MTOC↔nucleus coupling**; MTOC is an unanchored point-on-springs that
  cannot carry torque (KB flags this). The reference image's actin-MT crosstalk + transport is not built.

## 3. Cytosol — buffer-in-practice, active-CFD-in-capability, default-OFF
- The active-field scaffold (`biot_fluid_warp.py` + the FSI block) demonstrates transient diffusion arrays
  and IBM transfer, but it is **not a production Biot cell solver**. The production contract must solve
  `S p_dot+α∇·v_s+∇·q=s_water`, `q=−(k/μ)(∇p−ρb)`, and `v_f=v_s+q/φ` on the live domain with conservative
  membrane/nuclear boundary fluxes. The old `∂p/∂t=c_v∇²p` fixed box remains an analytic/component oracle.
- BUT it is additive `biot_fsi=None` **default-OFF**. Only `ff_fsi_native.py` (always) and
  `ff_internal_displacement_diag.py --fsi` (opt-in) activate it, and **NO committed native FSI-ON artifact
  exists**. The landed native "interior DECOUPLED" diagnostic (Nc=494802, 340e80b) ran **FSI OFF**.
- When OFF (every landed native run at the time of this audit): cytosol = uniform quasi-static
  `dP_osm(van't Hoff) + dP_solid(Terzaghi) +
  KK-drainage + incompressible-nucleus displacement`, clocked by the solid 6πηR drag — physiologically
  *parameterized* but spatially LUMPED (0-D), an algebraic function of hull volume. = the "static buffer".
- Physiological-baseline gap: the project's own HARD rule says the fluid layer must be ON at physiological
  values in production. This requires rebuilding it first; enabling the old held-face path does not close the
  gap. Resting mean turgor `p_bar=Π₀` and spatial `p_excess` are separate channels, and neither is expected to
  substitute for the explicit cap/LINC or actin–MT load paths.

## What is SOLID (not to under-sell)
Cortex actin (native-validated, 494k nodes), FA catch-slip clutches (native pN-validated), MT bending strut,
nucleus shell + incompressibility, IF cage (built), membrane Helfrich (this session), Biot CFD kernels
(component-validated), ECM library (6 materials). The FOUNDATION is real; what is missing is the ACTIVE
DIRECTED force-generating actin architecture ON TOP + turning the CFD ON.

## Recommended reprioritization (PI decides)
1. **Cytosol: build master I1a/I1b first** — conservative moving-domain storage/flux, `q-p`, reconstructed
   `v_f`, and one physical clock. The legacy `--fsi` path is a component oracle, not a cheap production switch.
2. **Stress fibers (agreed next)** — build the EMERGENT contractile ventral SF (FA→NMII-bundle→FA, myosin-
   driven alignment, NOT scripted) + the **perinuclear actin cap** (the real nucleus load path — likely what
   the polar IF cage failed to be). Active gate = Kumar 2006 10–30 nN.
3. **Cortex myosin → head-resolved Hill NMII** (currently lumped constant/aggregate control), reimplemented
   entirely in Warp CUDA. Per-head force/count/kinetic constants remain I0-B3 evidence decisions.
4. **Actin-MT coupling** — a real crosslinker along the MT length + MTOC↔nucleus (dynein) — lower priority.

## Change log
- 2026-07-16: created from 3 grounded re-audit subagents while the membrane native gate runs. Honest map of the
  actin/cytosol gaps; PI-authored reprioritization pending.
- 2026-07-16 implementation-readiness correction: downgraded the current Biot path to component scaffold,
  redirected the first step to conservative I1a/I1b, and replaced the lumped-motor cleanup with head-resolved
  I3 plus an explicit FV-law gate.

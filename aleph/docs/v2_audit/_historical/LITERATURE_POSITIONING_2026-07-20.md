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

# Literature positioning — explicit fine-grained actomyosin ⊗ two-way Biot poroelastic fluid (2026-07-20)

> **Purpose.** A DOI-verified deep review of where the `aleph/ac/` engine sits in the cell-mechanics
> modelling literature, triggered by the question *"the compartments carry no mass — is there a cytosol
> fluid, and is this direction defensible?"*. This is the literature-positioning companion to
> [`CELL_MECHANICS_FRAMEWORK_2026-07-15.md`](CELL_MECHANICS_FRAMEWORK_2026-07-15.md) (the FEM-fibers +
> CFD-fluid mesoscale thesis). It establishes the novelty claim, its **narrow** defensibility boundary,
> the closest comparators to benchmark against, and the two validations that run *separately* from the
> development hot loop.
>
> ⚠️ **Citation integrity (CLAUDE.md hard rule).** DOIs below were verified by research agents against
> Crossref / PubMed this session but are **NOT yet in the Notion KB `source_evidence`**. Confirm each
> `verdict=OK` (or register + `verify_sources.py`) before using any of these in a formal deliverable
> (manuscript, gate contract). Entries flagged `[unverified]` were single-source only.

---

## 0. Architecture facts this review rests on (from the code)

Three facts about the current engine, established by reading the runtime, that frame the literature comparison:

1. **The membrane–cortex boundary is a 4-part composite, not one "collision" kernel.**
   - Geometric nesting: membrane is the outer boundary at `R_CELL = 7.5 µm`; the actin cortex is a shell
     just inside at `R_CORTEX = 7.40 µm` (`CORTEX_MEMBRANE_GAP = 0.10 µm ≈ ½·h_cortex`). The membrane
     *contains* the cortex. (`ac/cell/assemble.py`)
   - ERM tethers: Hookean membrane↔cortex springs, force-free at rest, with **commit-safe** rupture
     (force kernel never mutates `bound`; rupture latches only at an accepted outer step). (`ac/cell/erm_tether.py`)
   - Pressure traction: enclosed pore fluid pushes the membrane outward via `(p_inside − p_ext)·A·n`
     (spatial live pressure, not lumped turgor; Young-Laplace `γ = ΔP·R/2` in the continuum limit).
     (`ac/cell/membrane_pressure.py`)
   - WCA excluded volume: all-filament steric (`σ_EV = 7 nm`) over a device HashGrid prevents
     interpenetration; membrane/nucleus/myosin marked dormant act as obstacles. (`ac/cell/assemble.py` StericForce)

2. **No inertial mass / weight — and that is correct physics, not an omission.** The inner solve is an
   overdamped mobility descent `pos += dt_mu·force`, `dt_mu = 0.1/max(k_mech, |Δp|·ℓ)` (`ac/cell/inner_mechanics.py`).
   At the cell scale Re ≈ 10⁻¹⁰; the inertial term is ~10 orders below viscous, so it is dropped (Purcell
   low-Reynolds regime). The dynamical quantity that replaces mass is **drag/mobility (γ)**, not `m`. The
   simulator is still fully *dynamic*: time enters through poroelastic relaxation, Hill FV, Bell kinetics,
   turnover, membrane water flux, and monomer transport — the "one outer physical clock + converged
   inertialess inner solve" architecture.

3. **The cytosol fluid is present — as a conserved Biot poroelastic pore field, with real physical quantities.**
   `S·dp/dt = mobility·∇²p − α·div(v_s) + s_water`, `q = −(k/μ)∇p` (`ac/fluid/biot_substrate.py`). It carries
   storage `S`, Darcy mobility `k/μ` (permeability/pore-viscosity), Biot `α`, pore pressure `p`, and Darcy
   flux `q`. The "missing physical quantity" the reviewer intuited is not mass — it is these poroelastic
   constants, and the open risk is whether they are **ON at physiological values in production** (see §5),
   not whether they exist.

### Fluid ↔ compartment coupling (all via one masked Eulerian grid + Peskin adjoint IBM)

| Pair | Mechanism | Direction | File |
|---|---|---|---|
| fluid ↔ actin skeleton | `−α∇p` body force on nodes (interp) ; `−α·div(v_s)` pressure source (spread) | two-way | `biot_substrate.py` `PressureCoupling` + `cell/fsi_coupling.py` `SolidDilatationCoupling` |
| fluid ↔ membrane | Kedem-Katchalsky `s = L_p(σΔΠ_osm − ΔP)` water flux ; `(p_in−p_ext)·A·n` traction | two-way | `fluid/boundary.py` `MembraneFluxBC` + `cell/membrane_pressure.py` |
| fluid ↔ nucleus | relative no-flux (NUCLEUS mask) ; nucleus displaces the fluid domain | kinematic | `fluid/domain.py` |
| fluid ↔ G-actin monomer | `u=φc` advected by `v_f` + diffusion + barbed/pointed reaction | two-way | `fluid/transport.py` |
| fluid ↔ myosin | dormant (`active<0`) — drives fluid only through actin | (indirect) | `cell/fsi_coupling.py` |

---

## 1. The literature splits into three non-overlapping camps

**No published system occupies our cell.** The field is two mirror-image camps plus a GPU-membrane camp:

| Camp | Representatives | Filaments / motors | Fluid | GPU | Decisive gap vs target |
|---|---|---|---|---|---|
| **A. Explicit skeleton, no conserved fluid** | aLENS (≤1M MT), Cytosim, MEDYAN, AFINES | ✅ explicit | ❌ local drag (free-draining, no HI) | ❌ CPU | no conserved fluid, no back-reaction, no pore pressure |
| **B. Explicit skeleton + two-way *Stokes*** | SkellySim/"twisters", Nazockdast-Shelley, Maxian-Donev | ✅ explicit | ✅ **single-phase Stokes** | some GPU (Maxian) | **Stokes ≠ Biot two-phase**; point-force / no motors; 10²–10³ fibers |
| **C. Two-way poroelastic fluid, lumped network** | Moeendarbary (exp), Callan-Jones-Jülicher, Strychalski-Copos-Guy, Radszuweit, Kulawiak, Taber | continuum + lumped `ζQ` | ✅ **Biot / two-phase poroelastic** | ❌ CPU/analytic | network is a homogenized stress tensor — no explicit filaments/heads/clutches |
| **(D. GPU two-way FSI, no cytoskeleton)** | HemoCell (17M RBC), Randles | ❌ membrane shell only | ✅ **LBM inertial bulk** | ✅ GPU | inertial single-phase LBM; interior is empty/continuum |

**Two structural observations for the manuscript:**
- **LBM is the wrong fluid for our claim.** Every GPU cell-FSI code (HemoCell, Randles, Microcosmos) is
  *inertial, bulk single-phase LBM* with a *membrane capsule and no internal cytoskeleton*. LBM does not
  represent a two-phase Biot pore-pressure medium — so "GPU cell FSI is already mature" is a different
  physics, not a precedent against us.
- **The only conserved-fluid ⊗ explicit-fiber whole-cell codes are the Flatiron Stokesian ones**
  (SkellySim, aLENS roadmap) — CPU/MPI, **single-phase Stokes not Biot**, no membrane/adhesion, no Bell
  clutches.

---

## 2. The novelty is real but the defensible statement is NARROW

`aleph/ac/` is the **inverse of Camp C**: it keeps the two-way Biot/Darcy poroelastic fluid and replaces
the continuum active-gel network with fully explicit filaments + motor heads + Bell clutches on GPU.

**Do NOT claim novelty as "poroelastic + active"** — that combination already exists and is continuum:
Callan-Jones & Jülicher 2011 (active permeating gel), Radszuweit 2014 / Kulawiak 2019 (active poroelastic
two-phase, *Physarum*), Taber 2011 (poroelastic crawling). Claiming it invites immediate refutation.

**Defensible novelty statement (three axes must co-occur, and do so nowhere in the literature):**
> Explicit discrete actomyosin — 70,686 cortical F-actin filaments (494,802 nodes) + head-resolved
> Stam-Hocky NMII (Hill force–velocity, per-head Bell kinetics) + Bell-Evans adhesion clutches — coupled
> **two-way** to a **two-phase Biot poroelastic pore-pressure** fluid, **inertialess**, **GPU-resident
> (Warp/CUDA)**, at full native whole-cell population. No literature camp holds all three of {two-phase
> Biot fluid physics · explicit motor/clutch fidelity · GPU whole-cell explicit cytoskeleton} at once.

---

## 3. Closest comparators to benchmark against (not just cite)

- **SkellySim / "twisters"** (Dutta, Farhadifar, Shelley et al., *Nat Phys* 2024; bioRxiv
  [10.1101/2023.04.04.534476]) — the strongest overall match: whole-cell, explicit fibers, two-way,
  inertialess, ~10³–10⁴ MT. Gaps: single-phase **Stokes not Biot**, motors = driving forces not Bell
  clutches, no membrane/adhesion, **CPU-MPI**. → the paper to benchmark the "whole-cell two-way
  fiber–fluid" capability against.
- **Nazockdast-Rahimian-Zorin-Shelley 2017** (*JCP* 329:173, [10.1016/j.jcp.2016.10.026]) +
  **Nazockdast et al. 2017** (*MBoC* 28:3261, [10.1091/mbc.E16-02-0108]) — explicit fibers ⊗ two-way
  Stokes; the aster's "porous-medium" behaviour is *emergent HI screening*, **not** a Biot pore field —
  citing it as poroelastic precedent would be a category error.
- **Maxian-Peláez-Mogilner-Donev 2021** (*PLoS Comput Biol* 17:e1009240, [10.1371/journal.pcbi.1009240])
  — closest *explicit-actin + HI on GPU*: 200–1,600 filaments, Stokes, **no motors** (transient
  crosslinkers only).
- **Strychalski-Copos-Lewis-Guy 2015** (*JCP* 282:77, [10.1016/j.jcp.2014.10.004]) + **Strychalski-Guy
  2016** (*Biophys J* 110:1168, [10.1016/j.bpj.2016.01.012]) — closest *two-way poroelastic IBM*, but the
  cortex is a homogenized elastic continuum with lumped contraction.

---

## 4. The empirical anchor we must reproduce EMERGENTLY

The continuum poroelastic literature hands us measured targets that our fine-grained engine should
reproduce **without putting them in** (Moeendarbary et al. 2013, *Nat Mater* 12:253, [10.1038/nmat3517]):
- poroelastic diffusion `Dp ≈ 40–61 µm²/s` (HeLa 41±11, HT1080 40±10, MDCK 61±10), `Dp = E·ξ²/μ`;
- hydraulic mesh size `ξ ≈ 14 nm`;
- poroelastic time `τ = L²/Dp ≈ 0.1–10 s`;
- non-equilibrating pressure over ~10 µm / ~10 s in blebbing (Charras et al. 2005, *Nature* 435:365,
  [10.1038/nature03550]).

Emergent reproduction of `Dp` from explicit filaments + pore fluid is the **killer validation**: it
directly rebuts the scale-separation reason the combination was left un-done (§5).

---

## 5. Defensibility hinge + the two validations (run SEPARATELY from development)

The strongest *steelman against* our approach is **scale separation**: `τ_poro = L²/Dp`, so a local ~1 µm
region equilibrates in ~0.01–0.03 s ≪ the actomyosin timescale → the fluid can be treated as locally
equilibrated and the explicit two-way solve dropped (this is exactly the coarse-graining that justifies
continuum active gel; Marchetti et al. 2013 RMP 85:1143; Jülicher et al. 2007 Phys Rep 449:3).

**Our approach is defensible only if** (a) the target question lives where scale separation *breaks down*
(blebbing, hydraulic fracture, fast AFM creep, mitotic rounding, nuclear mechanics, dense heterogeneous
cortex), and (b) filament-level heterogeneity feeds back on the fluid. Otherwise the "why not active gel"
critique stands.

Two acceptance checks establish (a)+(b). **Both run out-of-band — they do not change the production
runtime or block current assembly/NG-gate development:**

1. **FSI-on audit (diagnostic gate).** Confirm production config turns two-way FSI ON at physiological
   values — `div_vs` actually filled each outer step (not the historical ZEROS stub), and `S/k/μ/α/φ`
   sourced not water-default. Plus the profiler gate: zero authoritative GPU→CPU roundtrips in the hot
   loop. This is a config + profiler audit, not a physics run.
2. **Emergent `Dp` reproduction (acceptance oracle).** Drive the native poroelastic engine through a
   relaxation / micro-indentation protocol, measure the emergent poroelastic diffusion constant, compare
   to `Dp ≈ 40–61 µm²/s`. A dedicated validation experiment with its own protocol + post-processing;
   reads engine outputs, does not touch the development loop.

Open GAP parameters blocking (1)/(2) at native: `φ` (porosity), `c_0` (free-monomer pool — blocks the RAD
native run), `D_c` (draft ~2–6 µm²/s). These are PI-GAPs (see `ac/fluid/params_i0b1.yaml`,
`ac/cell/assemble.py`), not tunable to outcome.

**Codex vs PI split.** Codex (debug/verify, not authority): `τ_poro`/Biot-CFL derivation, Peskin
spread↔interp adjoint check, `−α∇p` ↔ `−α·div(v_s)` energy-conjugacy, and the σ_EV interpenetration
(~63k nodes, ≤1000 pN) that swamps the pN poroelastic signal. PI (strategy): fix the primary scientific
question in a scale-separation-breakdown regime, and ratify the emergent validation observable.

---

## 5b. PI decision (2026-07-20) + corrected parameter state

**Target regime RATIFIED: blebbing / hydraulic fracture.** This is the strongest fit because the observable
sits directly on the 4-part membrane–cortex boundary (§0.1): bleb nucleation = pore-pressure traction vs
membrane–cortex adhesion energy (ERM `f_rupt`) + membrane tension; then non-equilibrating pressure
(~10 µm / ~10 s), expansion ~30 s, retraction ~2 min (Charras 2005/2008). **It also resolves the §4
circularity:** bleb nucleation is genuinely *emergent* from explicit ERM rupture + pressure + cortex WCA and
is **independent of the `c_v = Dp` input anchor** — so the validation observable is bleb nucleation criterion
+ ΔP non-equilibration length/time scales, NOT "reproduce Dp" (which would be circular given the current
continuum closure).

**Corrected parameter state** (from `ac/fluid/params_i0b1.yaml`, read 2026-07-20 — earlier "φ/c_0/D_c are all
open GAPs" framing was imprecise):
- **I1a p/mass Biot field is CLOSED** (`native_run_blocked_on: []`): `c_v=Dp=50 µm²/s` anchored
  (Moeendarbary), `α=1.0` **PI-ratified 2026-07-16** (incompressible-constituent limit, Penta 2025), `k/μ`,
  `M`, `S` derived-consistent (draft), `L_p=1e-12` draft (Jung 2011 MCF7/AQP5), `Π₀=40 Pa` ratified proxy.
  → the FSI-on audit and pressure-relaxation checks can run now.
- **Still open (I1b/I1c only):** `φ` porosity (~0.7, unsourced) — deferred to I1b, absorbed into storage `S`
  for I1a; `c_0` MCF7 free-monomer pool — blocks I1c RAD native; `D_c` ~2–6 µm²/s draft.
- **Two consistency FINDINGS to verify (Codex + PI awareness):** (i) `M_biot` — the yaml flags that the code
  value `M=1 kPa` is inconsistent with the Moeendarbary `c_v` (~10× too slow); the `c_v`-anchored closure
  gives `M≈10 kPa`. Confirm which the code uses. (ii) `L_p` — yaml `1e-12` vs code `1.6e-8`: unit
  reconciliation open.

**Reprioritized for the blebbing target:**
- 🔴 `L_p` unit reconciliation — sets bleb fill rate.
- 🔴 ERM `f_rupt` + membrane–cortex adhesion energy + membrane tension — the nucleation-threshold sources;
  confirm each is sourced, not a tuned constant.
- 🟡 `φ` porosity — governs cytosol inflow `v_f = v_s + q/φ` into the bleb (moves up from I1b-deferred).
- 🟢 `c_0` / `D_c` (I1c monomer) — deprioritized for blebbing.

---

## 6. References (agent-verified this session; confirm in KB before deliverable use)

Verified vs Crossref/PubMed:
- Moeendarbary et al. 2013, *Nat Mater* 12:253 — 10.1038/nmat3517 (PMID 23291707)
- Charras et al. 2005, *Nature* 435:365 — 10.1038/nature03550 (PMID 15902261)
- Stewart et al. 2011, *Nature* 469:226 — 10.1038/nature09642 (PMID 21196934)
- Jiang & Sun 2013, *Biophys J* 105:609 — 10.1016/j.bpj.2013.06.021 (PMID 23931309)
- Kruse et al. 2005, *EPJE* 16:5 — 10.1140/epje/e2005-00002-5
- Joanny et al. 2007, *NJP* 9:422 — 10.1088/1367-2630/9/11/422
- Callan-Jones & Jülicher 2011, *NJP* 13:093027 — 10.1088/1367-2630/13/9/093027
- Marchetti et al. 2013, *Rev Mod Phys* 85:1143 — 10.1103/RevModPhys.85.1143
- Jülicher et al. 2007, *Phys Rep* 449:3 — 10.1016/j.physrep.2007.02.018
- Strychalski, Copos, Lewis, Guy 2015, *JCP* 282:77 — 10.1016/j.jcp.2014.10.004
- Strychalski & Guy 2016, *Biophys J* 110:1168 — 10.1016/j.bpj.2016.01.012
- Radszuweit et al. 2014, *PLoS ONE* 9:e99220 — 10.1371/journal.pone.0099220
- Kulawiak et al. 2019, *PLoS ONE* 14:e0217447 — 10.1371/journal.pone.0217447
- Taber et al. 2011, *JoMMS* 6:569 — 10.2140/jomms.2011.6.569
- Nazockdast et al. 2017, *JCP* 329:173 — 10.1016/j.jcp.2016.10.026
- Nazockdast et al. 2017, *MBoC* 28:3261 — 10.1091/mbc.E16-02-0108 (PMID 28331070)
- Shelley 2016, *Annu Rev Fluid Mech* 48:487 — 10.1146/annurev-fluid-010814-013639
- Wróbel, Cortez, Fauci 2014, *Phys Fluids* 26:113102 — 10.1063/1.4900941
- Maxian et al. 2021, *PLoS Comput Biol* 17:e1009240 — 10.1371/journal.pcbi.1009240
- Maxian, Mogilner, Donev 2021, *Phys Rev Fluids* 6:014102 — 10.1103/PhysRevFluids.6.014102
- Popov, Komianos, Papoian 2016 (MEDYAN), *PLoS Comput Biol* 12:e1004877 — 10.1371/journal.pcbi.1004877 (PMID 27120189)
- Freedman et al. 2017 (AFINES), *Biophys J* 113:448 — 10.1016/j.bpj.2017.06.003 (PMID 28746855)
- Yan et al. 2022 (aLENS), *eLife* 11:e74160 — 10.7554/eLife.74160
- Závodszky et al. 2017 (HemoCell), *Front Physiol* 8:563 — 10.3389/fphys.2017.00563
- Dutta et al. 2023 (SkellySim/twisters), bioRxiv — 10.1101/2023.04.04.534476
- Kashiwabara et al. 2024, *Sci Rep* 14:31339 — 10.1038/s41598-024-82864-z (PMID 39732914) [experiment; motivation only]
- Rajagopal, Holmes & Lee 2018, *WIREs Syst Biol Med* — 10.1002/wsbm.1407 (best continuum-vs-particle-cost review)

`[unverified]` (single-source / not resolved this session): Cytosim NJP 9:427 (10.1088/1367-2630/9/11/427);
Nazockdast twisters *Nat Phys* 2024 (10.1038/s41567-023-02372-1); Mogilner & Manhart 2018 *Annu Rev Fluid
Mech* (10.1146/annurev-fluid-010816-060238); Hu-Suo 2010 gel indentation; Guilak & Mow 2000 *J Biomech*
(PII S0021-9290(00)00105-6); Hu/Charras 2024 *BMMB* (10.1007/s10237-024-01854-2); Membrane-MEDYAN 2021
(10.1021/acs.jpcb.1c02336); CyLaKS EPJE 2021; Randles multi-GPU IB (10.1016/j.jocs.2020.101153);
Microcosmos (arXiv:2607.02954); Elastica++/SOPHT/Tekinalp 2025 CMAME.

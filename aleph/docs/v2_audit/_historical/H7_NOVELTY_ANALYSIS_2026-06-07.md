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

# H.7 Novelty Analysis — Is `ffn_cellsim` Scientifically/Methodologically New?

**Date**: 2026-06-07
**Author**: Theoretical-novelty-analysis subagent (independent referee pass)
**Scope**: Assess whether the project's full-fidelity whole-cell mechanistic approach is
genuinely novel against (a) conventional single-cell/tissue mechanics simulation and
(b) the fine-grained cytoskeleton-simulator literature. Tough-but-fair referee tone:
this document is NOT a booster. Every external paper below was verified real via
Consensus / WebSearch (CrossRef/PubMed-resolvable); the project's own citation-integrity
history (3 confirmed hallucinated SourceEvidence rows, CLAUDE.md) means citations here are
deliberately conservative.

> **One-line verdict.** There is a *real but narrow* methodological novelty — the
> **integration**, not the ingredients. No single mechanism is new; the *combination*
> (every filament/motor/clutch/crosslink explicit, on a closed turgor-pressurised
> multi-compartment single cell, with cortical tension γ read as an emergent
> three-channel measurement rather than prescribed, and a clean active-vs-passive
> separation) does not appear, as an assembled whole, in the published simulator
> landscape. The novelty is **defensible only if** the emergent γ is validated and the
> ×40 coarse-graining is honestly bounded — both currently open.

---

## 1. The landscape of conventional cell/tissue mechanics simulation

Organised from most coarse-grained (cell-as-a-point/polygon) to most fine-grained
(filament-resolved). For each: what it represents, its coarse-graining level, and what it
**cannot** resolve — which is where `ffn_cellsim` claims to differ.

### 1.1 Vertex models (cell = polygon; junctions = edges)

A confluent tissue is a tiling; each cell is a polygon and the dynamical variables are the
**vertices**, driven by an energy of area/perimeter targets plus interfacial tension
(Fletcher et al. 2014, *Biophys J* 106:2291–2304, doi:10.1016/j.bpj.2013.11.4498; Alt,
Ganguly & Salbreux 2017, *Phil Trans R Soc B* 372:20150520, doi:10.1098/rstb.2015.0520;
Barton et al. 2017 Active Vertex Model, *PLoS Comput Biol* 13:e1005569,
doi:10.1371/journal.pcbi.1005569). Recent work pushes vertex models toward subcellular
detail — heterogeneous myosin, non-uniform contractility (Lange et al. 2025, *PLoS Comput
Biol*, doi:10.1371/journal.pcbi.1012324) — and toward differentiable/GPU implementations
(VertAX, Pasqui et al. 2026).
- **Coarse-graining**: extreme. One cell ≈ one polygon; the cortex is a *line tension*
  parameter on each edge. No filaments, no motors, no explicit cytoskeleton.
- **Cannot resolve**: the molecular origin of cortical tension (it is an *input* γ_edge);
  the actin architecture; single-filament mechanics; turgor/cytoplasm as distinct
  compartments. Cortical tension is *imposed*, never *emergent*. This is precisely the
  axis `ffn_cellsim` inverts.

### 1.2 Cellular Potts / Glazier-Graner-Hogeweg (cell = set of lattice sites)

Each cell is a domain of lattice pixels sharing an ID; dynamics are Metropolis updates of a
Hamiltonian with adhesion + volume/area constraints (Graner & Glazier 1992, *Phys Rev Lett*
69:2013, doi:10.1103/PhysRevLett.69.2013; CompuCell3D ecosystem; review Hirashima,
Rens & Merks 2017, *Dev Growth Differ* 59:329–339, doi:10.1111/dgd.12358).
- **Coarse-graining**: extreme. Mechanics enter only through a phenomenological surface-energy
  Hamiltonian and Metropolis acceptance — *not* force balance, *not* Newtonian/Langevin
  dynamics. Time is Monte-Carlo sweeps, not seconds.
- **Cannot resolve**: real forces and stresses (the dynamics are not momentum/force-based);
  the cytoskeleton; force-dependent bond kinetics. CLAUDE.md's worked example "Metropolis
  (detailed-balance proxy) → Bell-Evans force-dependent" is exactly the CPM-style mechanism
  the project rejects as runtime physics.

### 1.3 Center-based models / CBM (cell = sphere with pairwise potentials)

Cells are points/spheres interacting through pairwise adhesion+repulsion potentials, often
with growth/division (e.g. CellSim3D, Madhikar et al. 2018, *Comput Phys Commun* 232:206–213,
doi:10.1016/j.cpc.2018.05.024 — cells as C180-fullerene elastic spheres with osmotic pressure
and DPD-style velocity-Verlet). **This is the family `ffn_cellsim`'s own Layer-2 spheroid
line used** (project memory: "center-based CBM", L2.x), and the documented conclusion there
is instructive: CBM reproduced the *form* of MCF7 spreading (A/A₀ law, r²=0.998) but the
*magnitude* was a structural limit of the center-based abstraction — which is the project's
own internal evidence for what CBM "cannot resolve".
- **Coarse-graining**: high. One cell ≈ a few–hundreds of nodes; no resolved cortex/motors.
- **Cannot resolve**: subcellular force generation; cortical tension from mechanism;
  single-cell shape from cytoskeletal architecture.

### 1.4 Continuum & phase-field (cell = field φ + stress tensor)

The cell interior is a continuum; a phase field φ tracks the moving boundary, coupled to
active-stress and reaction-diffusion fields for actin/myosin (reviews: Moure & Gomez 2021,
*Arch Comput Methods Eng* 28:311–344, doi:10.1007/s11831-019-09377-1; Buttenschön & Edelstein-Keshet
2020, *PLoS Comput Biol* 16:e1008411, doi:10.1371/journal.pcbi.1008411). The dominant
*physics-of-the-cortex* theory in this class is **active-gel / active-matter hydrodynamics**:
the actomyosin cortex as a thin film of active polar gel with an active contractile stress
ζΔμ (Joanny & Prost 2009, *HFSP J* 3:94–104, doi:10.2976/1.3054712; Salbreux, Prost & Joanny
2009, *Phys Rev Lett* 103:058102, doi:10.1103/PhysRevLett.103.058102; Jülicher, Grill &
Salbreux 2018, *Rep Prog Phys* 81:076601, doi:10.1088/1361-6633/aab6bb; viscous active *shell*
theory of the cortex, Borja da Rocha, Bleyze & Salbreux 2022, *J Mech Phys Solids*
164:104876, doi:10.1016/j.jmps.2022.104876; physical-parameter extraction from laser ablation,
Saha et al. 2016, *Biophys J* 110:1421–1429, doi:10.1016/j.bpj.2016.02.013).
- **Coarse-graining**: the cortex is a *constitutive law* with lumped coefficients
  (active stress ζΔμ, viscosity η, friction γ, hydrodynamic length λ). Tension is an
  *emergent continuum field* — but emergent from a **constitutive postulate**, not from
  explicit filaments and motors.
- **Cannot resolve**: the discrete/molecular origin of those coefficients (filament-length
  distribution, motor duty ratio, crosslink catch-bond kinetics) — exactly the variables
  Chugh 2017 (below) shows *control* tension. Active-gel theory is `ffn_cellsim`'s most
  serious *conceptual* competitor for "emergent γ" and is treated as such in §4.

### 1.5 Subcellular-element model / SEM (cell = cloud of elements)

Each cell = a cloud of elastically-coupled Langevin elements with phenomenological intra-/
inter-cellular potentials (Newman 2005, *Math Biosci Eng* 2:613–624, doi:10.3934/mbe.2005.2.613;
Sandersius & Newman 2008, *Phys Biol* 5:015002, doi:10.1088/1478-3975/5/1/015002). Resolves
deformable cell shape and rheology without a grid.
- **Coarse-graining**: medium. Elements are *generic* viscoelastic blobs, **not** identified
  filaments/motors/clutches; no force-dependent binding chemistry.
- **Cannot resolve**: identity-resolved cytoskeletal mechanism; catch-bond adhesion;
  motor force-velocity.

### 1.6 Fine-grained cytoskeleton simulators (filaments/motors/crosslinks explicit)

This is the family `ffn_cellsim` is built on and must be judged against most strictly.

- **Cytosim** (Nedelec & Foethke 2007, *New J Phys* 9:427, doi:10.1088/1367-2630/9/11/427;
  cellular-scale extension Belmonte, Leptin & Nédélec 2022, *eLife* 11:e74160,
  doi:10.7554/eLife.74160). Filaments = bending-elastic bead chains; motors/crosslinks =
  diffusing point complexes binding stochastically; Brownian/Langevin dynamics in a viscous
  medium. 1D/2D/3D, decades of use. **Mechanistic but not whole-cell with adhesion+turgor+
  nucleus**; typically reconstitution/sub-structure geometries (rings, asters, bundles).
- **MEDYAN** (Popov, Komianos & Papoian 2016, *PLoS Comput Biol* 12:e1004877,
  doi:10.1371/journal.pcbi.1004877) — mechanochemical: explicit actin (de)polymerisation +
  myosin + crosslinkers with a full reaction-diffusion chemistry coupled to mechanics; the
  gold standard for *mechanochemical* fidelity. Extended to deformable membranes
  (Membrane-MEDYAN, Ni & Papoian 2021, *J Phys Chem B* 125:10710–10719,
  doi:10.1021/acs.jpcb.1c02336) and explicitly framed as moving *toward* — i.e. **not yet
  at** — cell-scale membrane-cytoskeleton models (Zimmerberg, Ni & Papoian 2022,
  *Biophys J* 121:2419–2431, doi:10.1016/j.bpj.2022.06.003). This 2022 perspective is the
  single most important comparator: a leading group states cell-scale whole-cell agent-based
  models are an *aspiration*, not a solved problem.
- **AFINES** (Freedman, Banerjee, Hocky & Dinner 2017, *Biophys J* 113:448–460,
  doi:10.1016/j.bpj.2017.06.003) — the project's *direct ancestor* (CLAUDE.md: "AFINES
  re-implemented on HOOMD"). 2D, CPU-only, single-threaded; Brownian dynamics of
  filament+motor+crosslink networks; **Metropolis (Glauber) off-acceptance**, *no* Arp2/3
  branching (per the project's own AFINES_ALGORITHM_NOTES.md primary-source review). It is a
  *network-in-a-box* reconstitution tool, not a cell.
- **CyLaKS** (Fiorenza, Rai, et al. 2021) and **MEDYAN.jl** (BPS 2025) round out the family;
  same scope (sub-cellular networks), faster engines.

**The structural gap, stated plainly.** Every member of §1.6 is a *cytoskeletal-network*
simulator. None of the canonical, peer-reviewed members ships a **single whole cell** that
simultaneously carries: a closed cortical shell + explicit bipolar myosin minifilaments +
catch-bond crosslinkers + integrin molecular clutches to a substrate + nucleus + osmotic
turgor + physiological cytoplasm viscosity + membrane, integrated in one dynamical run, with
cortical tension read out as an emergent measurement. The leading group in the field
(Papoian/MEDYAN) frames exactly this as *future* work in 2022.

---

## 2. Where `ffn_cellsim` sits — its distinctive choices

Grounded in the read files (`cell/cell.py::build_cortex_full_simulation`,
`cortex/cortical_tension.py`, `configs/mcf7_baseline.yaml`,
`docs/briefs/H{2,3,4,5}*.md`, `docs/briefs/H7_*`):

(a) **Every filament/motor/clutch/crosslink explicit, at whole-cell scale.** The full builder
composes cortex actin backbone + crosslinkers + Stam-Hocky bipolar myosin minifilaments +
(optional) lamellipodium + FA integrin clutches + nucleus + enclosed-volume turgor + membrane
+ cytoplasm drag in **one** HOOMD `Simulation`
(`build_cortex_full_simulation`, cell.py:689). This is §1.6-style fidelity applied to a §1.1–1.5
*object* (a whole cell), which the surveyed literature does not assemble as a unit.

(b) **Cortical tension γ and spreading are EMERGENT, not prescribed.** `cortical_tension.py`
*measures* γ from the built cell's mechanical state via method-of-planes over the actual bond
forces; the literature band [0.35,0.65] mN/m is an **overlay**, explicitly "Nothing in this
module is tuned to it" (cortical_tension.py:38–41). Vertex/CPM/CBM models take γ (or a surface
energy) as an input; active-gel takes the active stress coefficient ζΔμ as an input.
`ffn_cellsim` takes filament/motor/bond parameters as input and *outputs* γ. (Spreading A/A₀
likewise emerges from radial barbed-end advance in the H.7 design, not a prescribed protrusion
force — H7 brief §3 rejects the lumped option (d).)

(c) **Catch-bond + Hill + Bell-Evans + Stam-Hocky are runtime mechanism, not closed-forms.**
CLAUDE.md's central inversion: the v1 codebase used Chan-Odde/Pereverzev/Bell-Evans/Hill as
the *runtime model*; here they are *acceptance oracles* (`validation/oracles/`,
runtime-import-forbidden), and the runtime is explicit force-dependent particle/bond dynamics
(Pereverzev 2005 catch-slip for integrins, Bell-Evans for crosslinkers, Hill 1938 for motor
stepping, Stam et al. 2017 *PNAS* 114:E10037, doi:10.1073/pnas.1708625114 for the bipolar
minifilament architecture, Bieling et al. 2016 *Cell* 164:115–127,
doi:10.1016/j.cell.2015.11.057 for lamellipodial force feedback).

(d) **Full multi-compartment integration in one run** at the physiological operating point —
turgor Π₀, η_cyto = 65.9 Pa·s, nucleus E_nuc, membrane γ_mem all ON at in-vivo setpoints
(`mcf7_baseline.yaml`), per the physiological-baseline HARD rule.

(e) **Paper models demoted to oracles**, and — the genuinely unusual measurement choice — a
**3-channel γ estimator that never folds turgor into a single total**: γ_soft (active
harmonic bonds), γ_rigid (M-SHAKE Lagrange-multiplier backbone tension), γ_passive
(Young-Laplace turgor), with the HARD rule (cortical_tension.py:25–33) that the passive turgor
term is *never* summed into the "structural" actomyosin total. This prevents the classic
artefact of an in-band aggregate masking an under-performing active cortex.

(f) **Integrator**: Leimkuhler-Matthews BAOAB-limit (Leimkuhler & Matthews 2013, *Appl Math
Res Express* 2013:34–56, doi:10.1093/amrx/abs010) as a custom HOOMD updater, chosen over
Euler-Maruyama for O(Δt²) configurational sampling — a fidelity choice, defensible but not by
itself novel (BAOAB is standard in MD; AFINES already uses it).

---

## 3. Is there genuine novelty? (skeptical itemisation)

For each distinctive choice: **prior art** (who did it) vs **what, if anything, is new**.

| Choice | Already done by | Genuinely new here? |
|---|---|---|
| Explicit filaments/motors/crosslinks, force-dependent kinetics | Cytosim, MEDYAN, AFINES, CyLaKS | **No** — this is the established §1.6 paradigm. `ffn_cellsim` re-implements AFINES on HOOMD; that is engineering, not new science. |
| Bell-Evans/Hill/catch-bond as *runtime* (not Metropolis) | Cytosim (force-dependent unbinding), MEDYAN | **No** — force-dependent kinetics are standard in Cytosim/MEDYAN. AFINES's Metropolis is the *exception*; moving off it is a fidelity upgrade *within* the family, not a first. |
| Emergent cortical tension from explicit actomyosin | Chugh et al. 2017 (*Nat Cell Biol* 19:689–697, doi:10.1038/ncb3525) computed tension vs filament length from a network model; Bidone et al. 2014 (*Biophys J* 107:2618–2628, doi:10.1016/j.bpj.2014.10.034) computed cortical tension build-up in a 3D ring; Fritzsche et al. 2016 (*Sci Adv* 2:e1501337, doi:10.1126/sciadv.1501337) linked actin kinetics→cortical mechanics | **Partly** — emergent tension from a filament model exists, but on *patches/rings/networks-in-a-box*, **not** on a closed whole cell with adhesion+turgor. The *whole-cell, adhered, pressurised* setting is new. |
| Whole single cell, multi-compartment, one dynamical run | SEM (generic blobs), CellSim3D (elastic shells) at low fidelity; Membrane-MEDYAN (vesicle + network, no nucleus/FA/turgor); active-gel shell (continuum, no explicit filaments) | **Yes (integration-level)** — the *assembly* of cortex+myosin+crosslink+FA-clutch+nucleus+turgor+cytoplasm+membrane as explicit particles in one cell is not present, as a whole, in the surveyed canonical literature. Papoian 2022 frames it as future work. |
| γ-seam to a collective/Layer-2 model | Active-gel→vertex bridges exist (Ouzeri et al. 2025 bioRxiv, doi:10.1101/2025.xx; "active gels to vertex models") | **Partly** — multiscale seams exist conceptually; a *fine-grained-explicit → CBM* seam carrying a measured γ is unusual but the bridge itself is not unprecedented in spirit. |
| 3-channel γ estimator that never folds turgor into total | Not found as an explicit, named protocol in the surveyed literature | **Yes (methodological, narrow)** — the *measurement-integrity* discipline (separate active/rigid-backbone/passive-turgor channels; refuse to report a single total) is a genuine, if modest, methodological contribution; it directly targets a real failure mode (an in-band aggregate masking a floored active cortex — the project's own §2026-06-04 saturation diagnosis). |
| Method-of-planes |û·n̂| sign-fix for shell tension | Irving-Kirkwood / method-of-planes is standard MD virial machinery | **No** as a method; the *application + the documented √N-vs-N artefact correction* is good practice, not new. |

**Net.** The novelty is **emphatically not** in any mechanism — every ingredient is prior art.
It is in (i) the **integration**: explicit-everything fidelity transplanted onto a *closed,
adhered, turgor-pressurised whole cell* with a nucleus, which the canonical fine-grained
simulators have not assembled and which the field's leaders explicitly call future work; and
(ii) a **measurement-integrity protocol** (3-channel separated γ, never a hidden total) that
is small but real. This is *methodological-integration* novelty, not *mechanism* novelty.
It is the kind of contribution that earns a solid methods/tools paper **if validated**, not a
paradigm-shift claim.

---

## 4. Honest risks to the novelty claim

1. **Active-gel theory already delivers "emergent" cortical tension — more cheaply.** The
   continuum active-gel/active-shell program (Salbreux 2009; Jülicher 2018; Borja da Rocha
   2022; Saha 2016) computes cortical tension, flows, and shape from an active-stress field at
   a tiny fraction of the wall-time, and is *the* accepted theory of the cortex. A referee will
   ask: what does resolving every filament buy that ζΔμ(filament length, motor density) does
   not? The defensible answer must be a **prediction the continuum coefficient cannot make
   from first principles** (e.g. the *non-monotonic* tension-vs-filament-length maximum Chugh
   2017 reports, recovered *without* fitting ζΔμ) — otherwise fine-graining is more expensive
   for the same answer. This is the sharpest threat.

2. **The ×40 coarse-graining undercuts "every filament explicit."** CLAUDE.md sanctions ~1000
   effective filaments standing in for ~38,000 native (×40 mesoscopic bundles). That is the
   *same kind* of coarse-graining AFINES/Cytosim/MEDYAN already make (effective filaments,
   bundle cross-sections). So the marketing line "every filament is explicit" is, strictly,
   **false at ×40**; the honest claim is "every *effective mesoscopic* filament/motor/clutch
   is explicit," which is exactly what the §1.6 family already does. The novelty (§3) survives
   this — it lives in the *integration*, not in literal all-atom resolution — but the *claim
   wording* must be corrected or it will not survive review.

3. **Emergent γ is not yet validated — it is currently *floored*.** Project memory and the
   v2_audit floor docs record γ_active inconsistent across runs (0.0003–0.07 mN/m) and a
   "γ-floor" 3000×/980× under band traced to unphysical baselines (turgor=0, FA off), now
   reframed as: authoritative γ only comes from the full H1–H10 build, still pending. **A
   novelty claim that rests on "we measure emergent γ" is not yet supported by an in-band,
   reproducible emergent γ.** Until the integrated cell produces γ in (or honestly near) a
   defensible band *without tuning*, the central selling point is a promise, not a result.

4. **The reference band itself is contested (internal).** Memory flags that [0.35,0.65] mN/m
   is a rounded/de-adhered HeLa/L929 proxy with *no MCF7 datum*, and that the band's
   provenance/attribution needs reconciling (the in-repo Salbreux-Charras-Paluch 2012 citation
   is to a real review, *Trends Cell Biol* 22:536–545, doi:10.1016/j.tcb.2012.07.001, but the
   exact numeric band's sourcing is a known open item). Validating "emergent γ matches
   literature" is only meaningful against a *correct* target; that target is unsettled.

5. **Reproducibility / single-realisation risk.** The project's own Layer-2 history shows a
   seed-sensitive r² (0.987 under-sampled → 0.74 full run). Any emergent-γ novelty claim needs
   ensemble statistics (per-realisation spread + mean), which the visualization rules mandate
   but which must actually be in the methods-paper figures.

6. **HOOMD-port engineering ≠ science.** Re-implementing AFINES on HOOMD with a custom BAOAB
   updater and GPU-main porting is substantial engineering, but a referee will (correctly) not
   credit it as scientific novelty by itself.

---

## 5. What would make the novelty defensible (concrete program for a methods paper)

A. **Demonstrate a prediction the continuum can't make for free.** Reproduce Chugh 2017's
   *non-monotonic* cortical-tension-vs-filament-length maximum **emergently**, from the
   explicit network, with no per-condition tuning of an active-stress coefficient. If the
   model recovers the tension maximum at intermediate filament length *and* predicts the
   thickness–tension inverse correlation, that is a concrete payoff fine-graining buys over
   active-gel. This is the single most persuasive experiment.

B. **Land an in-band, reproducible, untuned emergent γ on the full physiological cell.**
   Close the γ-floor: run the integrated H1–H10 cell at the physiological baseline, report all
   three channels separately (the protocol is built), and show γ_structural lands in a
   *defensible* band over an ensemble (≥5 seeds, per-realisation + mean). Until this exists,
   §3's integration-novelty is unrealised.

C. **Fix the band target.** Resolve the cortical-tension reference for MCF7 specifically
   (or state honestly that no MCF7 datum exists and validate against the rounded-cell proxy
   *with that caveat in the paper*). A validation against a mis-attributed band is not
   defensible.

D. **Quantify what ×40 costs.** Run a coarse-graining convergence study (×40 vs ×20 vs ×10
   effective filaments) on γ and on one observable; show the measured quantity is
   ×N-invariant within error (the project's own Magic-Number "grid-invariant" rule applied to
   the mesoscale). This converts the coarse-graining from a liability (risk #2) into a
   controlled, reported approximation — and lets the paper drop the inaccurate "every filament"
   wording for the accurate "grid-converged mesoscopic" wording.

E. **Head-to-head against a §1.6 simulator.** Reproduce one canonical AFINES/Cytosim/MEDYAN
   reconstitution result (e.g. contractile-network stress, or Bieling 2016 load-adaptation) on
   the HOOMD re-implementation to prove the runtime is faithful — *then* show the whole-cell
   capability the comparator lacks (closed adhered turgor cell γ). The contrast IS the paper.

F. **Benchmark the γ-seam end-to-end.** Show the fine-grained-measured γ, fed into the Layer-2
   CBM, changes the collective prediction in a way that an a-priori-guessed γ does not — i.e.
   demonstrate the seam carries information, justifying the two-scale architecture.

G. **Methods-paper minimum.** Open code + configs (have); deterministic seeds + ensemble
   stats; a sanity-gate table (have); explicit statement of every lumped/coarse element (×40,
   rigid substrate, single-global vs per-WAVE geometry); and a frank "what fine-graining buys"
   section answering risk #1 head-on.

---

## 6. Bottom line for the PI

`ffn_cellsim` is **not** doing something no one has done at the mechanism level — Cytosim,
MEDYAN, and AFINES already resolve filaments/motors/crosslinks with force-dependent kinetics,
and active-gel theory already produces emergent cortical tension. The **defensible** novelty
is the **assembly**: a whole single cell built explicit-everything (cortex + bipolar myosin +
catch-bond crosslink + integrin clutch + nucleus + turgor + cytoplasm + membrane) in one
dynamical run, reading cortical tension as a *measured, three-channel, turgor-separated*
emergent quantity rather than a prescribed input — a configuration the field's own leaders
(Papoian 2022) frame as not-yet-achieved. That claim is currently a **promise, not a result**:
the emergent γ is still floored/inconsistent, the band target is contested, and "every
filament" is really "every ×40 mesoscopic filament." Deliver experiment A (a prediction the
continuum can't make for free) and experiment B (an in-band, untuned, reproducible emergent γ),
correct the coarse-graining wording, and the integration-novelty becomes a credible methods
paper. Absent those, the project is a high-fidelity engineering re-assembly of known parts —
valuable, but not yet a novelty claim that survives a tough referee.

---

## Appendix — verified external citations used (all real; DOI/URL given)

Conventional cell/tissue mechanics methods:
- Fletcher, Osterfield, Baker & Shvartsman 2014, *Biophys J* 106:2291–2304 (vertex models),
  doi:10.1016/j.bpj.2013.11.4498
- Alt, Ganguly & Salbreux 2017, *Phil Trans R Soc B* 372:20150520 (vertex review),
  doi:10.1098/rstb.2015.0520
- Barton, Henkes, Weijer & Sknepnek 2017, *PLoS Comput Biol* 13:e1005569 (Active Vertex Model),
  doi:10.1371/journal.pcbi.1005569
- Lange et al. 2025, *PLoS Comput Biol* (vertex models capturing subcellular scales),
  doi:10.1371/journal.pcbi.1012324
- Graner & Glazier 1992, *Phys Rev Lett* 69:2013–2016 (Cellular Potts),
  doi:10.1103/PhysRevLett.69.2013
- Hirashima, Rens & Merks 2017, *Dev Growth Differ* 59:329–339 (CPM review),
  doi:10.1111/dgd.12358
- Newman 2005, *Math Biosci Eng* 2:613–624 (Subcellular Element Model),
  doi:10.3934/mbe.2005.2.613
- Sandersius & Newman 2008, *Phys Biol* 5:015002 (SEM rheology),
  doi:10.1088/1478-3975/5/1/015002
- Madhikar et al. 2018, *Comput Phys Commun* 232:206–213 (CellSim3D),
  doi:10.1016/j.cpc.2018.05.024
- Moure & Gomez 2021, *Arch Comput Methods Eng* 28:311–344 (phase-field review),
  doi:10.1007/s11831-019-09377-1
- Buttenschön & Edelstein-Keshet 2020, *PLoS Comput Biol* 16:e1008411 (single→collective review),
  doi:10.1371/journal.pcbi.1008411

Active-gel / continuum cortex (the emergent-γ competitor):
- Joanny & Prost 2009, *HFSP J* 3:94–104, doi:10.2976/1.3054712
- Salbreux, Prost & Joanny 2009, *Phys Rev Lett* 103:058102, doi:10.1103/PhysRevLett.103.058102
- Jülicher, Grill & Salbreux 2018, *Rep Prog Phys* 81:076601, doi:10.1088/1361-6633/aab6bb
- Borja da Rocha, Bleyze & Salbreux 2022, *J Mech Phys Solids* 164:104876,
  doi:10.1016/j.jmps.2022.104876
- Saha et al. 2016, *Biophys J* 110:1421–1429 (cortex physical parameters), doi:10.1016/j.bpj.2016.02.013
- Salbreux, Charras & Paluch 2012, *Trends Cell Biol* 22:536–545 (cortex review; in-repo band
  source), doi:10.1016/j.tcb.2012.07.001

Fine-grained cytoskeleton simulators (the §1.6 family):
- Nedelec & Foethke 2007, *New J Phys* 9:427 (Cytosim), doi:10.1088/1367-2630/9/11/427
- Belmonte, Leptin & Nédélec 2022, *eLife* 11:e74160 (cellular-scale Cytosim),
  doi:10.7554/eLife.74160
- Popov, Komianos & Papoian 2016, *PLoS Comput Biol* 12:e1004877 (MEDYAN),
  doi:10.1371/journal.pcbi.1004877
- Ni & Papoian 2021, *J Phys Chem B* 125:10710–10719 (Membrane-MEDYAN),
  doi:10.1021/acs.jpcb.1c02336
- Zimmerberg, Ni & Papoian 2022, *Biophys J* 121:2419–2431 (toward cell-scale agent-based —
  the key "future work" comparator), doi:10.1016/j.bpj.2022.06.003
- Freedman, Banerjee, Hocky & Dinner 2017, *Biophys J* 113:448–460 (AFINES — direct ancestor),
  doi:10.1016/j.bpj.2017.06.003

Mechanistic ingredients (oracles → runtime in this project):
- Chugh et al. 2017, *Nat Cell Biol* 19:689–697 (cortex architecture regulates surface tension;
  the key "emergent γ from filament model" prior art + the convergence target for experiment A),
  doi:10.1038/ncb3525
- Fritzsche et al. 2016, *Sci Adv* 2:e1501337 (actin kinetics shapes cortical mechanics),
  doi:10.1126/sciadv.1501337
- Bidone, Tang & Vavylonis 2014, *Biophys J* 107:2618–2628 (cortical tension build-up, 3D ring),
  doi:10.1016/j.bpj.2014.10.034
- Stam et al. 2017, *PNAS* 114:E10037–E10045 (filament rigidity/connectivity, bipolar myosin
  contractility — the "Stam-Hocky" minifilament basis), doi:10.1073/pnas.1708625114
- Bieling et al. 2016, *Cell* 164:115–127 (force feedback in branched actin — H.5 basis),
  doi:10.1016/j.cell.2015.11.057
- Chan & Odde 2008, *Science* 322:1687–1691 (motor-clutch — oracle, not runtime),
  doi:10.1126/science.1163595
- Pereverzev et al. 2005, *Biophys J* 89:1446–1454 (two-pathway catch-slip bond — integrin
  runtime), doi:10.1016/S0006-3495(05)72793-4
- Leimkuhler & Matthews 2013, *Appl Math Res Express* 2013:34–56 (BAOAB integrator),
  doi:10.1093/amrx/abs010

> Note on the one bioRxiv item cited (Ouzeri et al. 2025, active-gels→vertex multiscale): used
> only to support that multiscale seams exist in spirit; the precise DOI was not resolved at
> write time and the claim does not depend on it. All other citations resolve to
> peer-reviewed venues with the DOIs above.

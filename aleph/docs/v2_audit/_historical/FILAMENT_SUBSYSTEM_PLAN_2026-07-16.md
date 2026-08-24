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

# Filament subsystem — gaps & plan for the Active Cell (`ac`) engine

**The SOLID-layer (③ actomyosin) architecture facet of the new engine.** The 2026-07-16 new-engine
pivot has three planning facets, all targeting the same `aleph/ac/` package:
- `NEW_ENGINE_BUILD_PLAN_2026-07-16.md` (I0..I9, five pillars) — the full engine spine.
- `../cfd_transport_program/ENGINE_ARCHITECTURE_PLAN.md` — the **living-cortex** spec: ① Biot
  poroelastic fluid ⊗ ② mass-conserving G-actin monomer RAD field ⊗ ③ discrete actomyosin. Owns the
  fluid + monomer-transport + turnover layers (①②).
- **THIS doc** — the ③ solid layer's **architecture & load paths**: the 3 cytoskeletal systems, the
  actin sub-architectures (cortex/SF/lamellipodium/filopodium), FA/traction load paths, de-adhesion,
  membrane load path, nucleus coupling. Complements both: the ac plan's I3–I8 own most of these
  increments; the living-cortex plan supplies the ② turnover that makes this ③ architecture *alive*.

**PI I0-A contract 2026-07-16:** the package is Warp-CUDA-only; HOOMD is never executed. NMII is explicit
backbone+individual heads with Hill FV. The first native baseline is MCF7×collagen×α2β1 with 70,686 cortical
F-actin filaments; every non-cortical population is separately density×geometry-derived and entered in the
global unique-active/allocated/node/head/GPU-byte ledger.

This doc takes everything established in the 2026-07-16 filament walkthrough and organizes it as a
**filament-subsystem spec** for the `ac/` build: what a real cell has, what the current (FF) code
has, where the gap is, and which increment closes it.

**How to read the status column.** `✓ correct` = implemented AND physically right (preserve it in
`ac`). `~ partial` = implemented but incomplete/approximate. `✗ gap` = absent from the running
model. `⭐NEW` = surfaced in this discussion and **NOT explicitly in the build plan** — the net-new
work this doc adds.

**Honest framing.** Most load-path gaps below (SF wiring, cap→LINC, actin–MT coupling, NMII motor)
are **already** captured by build-plan increments I3/I6/I7/I8 — this doc's value is (a) the
filament/load-path lens with per-item status, (b) a small set of **design invariants we got right
and must not regress**, and (c) three **genuinely new** items the walkthrough exposed
(keratin/vimentin IF differentiation; the adhesion-class emergence principle; the two-traction-channel
gate). It does not re-invent gaps the plan already owns.

Memory: [[project-active-cell-c-redesign]], [[project-tensegrity-mt-compression-actin-tension]],
[[project-nucleus-coupling-if-membrane]], [[project-mda-mb-231-validation-scope]],
[[project-h7-sf-array-traction]], [[feedback-junction-switch-fine-grained]],
[[project-unified-actin-architecture]].

---

## 1. Filament inventory — what a cell has vs what we have

### 1.1 The three systems

| System | Ø | Mechanical role | `ac`-relevant status | Increment |
|---|---|---|---|---|
| **Actin (microfilament)** | ~7 nm | tension, contraction, protrusion | cortex `✓`; SF/lamellipodium/filopodium `✗/~` (below) | I3–I6 |
| **Microtubule (MT)** | ~25 nm | **compression strut** (emergent from κ_bend), transport rail | strut `✓ correct` (not lumped); but **block-diagonal** to actin `✗` | I8 |
| **Intermediate filament (IF)** | ~10 nm | toughness, large-strain, nucleus↔cortex load path | radial-spoke cage `✗ FAILED` (buckles, 6× over-stiff); **keratin/vimentin undifferentiated** `⭐NEW` | I7 (replace) + §5.1 |

Key principle we reconfirmed: **MT = compression, actin = tension** (tensegrity), and MT's
load-bearing must **emerge from bending stiffness**, never a lumped strut. The FF `microtubule.py`
already does this (κ=20 pN·µm²) — preserve. The defect is only the **coupling**: actin and MT are
elastically independent (tip-contact only) → I8 un-block-diagonalizes.

### 1.2 Actin sub-architectures (one category, `[[project-unified-actin-architecture]]`)

Same molecule (actin + NMII + crosslinkers); the structure is set by **crosslinker type + myosin
arrangement + anchoring**. This is exactly why `ac` must let them **emerge** from one network
(P3), not hard-code each.

| Architecture | Crosslinker | Myosin arrangement | Anchor | Status | Increment |
|---|---|---|---|---|---|
| **Cortex** | filamin (network) | dispersed → **isotropic γ** | membrane (ERM) | `✓` (seed) | I4 seed |
| **Ventral SF** | α-actinin (bundle), sarcomeric | periodic → **directional** | FA↔FA | `✗ standalone spec only` | I6 |
| **Perinuclear cap SF** | α-actinin | periodic | ACAFA + **LINC→nucleus** | `✗` | I7 |
| **Dorsal SF / transverse arc** | α-actinin/fascin | graded/periodic | one FA + lamella | `✗` | I4 architecture / I6 adhesion |
| **Lamellipodium** | Arp2/3 branched net | (lamellar myosin) | nascent adhesions | `~ ratchet only` | I4 rebuild / I6 adhesion |
| **Filopodium** | fascin (tight bundle) | — | **tip adhesion** | `✗ spec only` | I4 architecture / I6 adhesion |

The SF/cortex distinction is **not** two motor mechanisms — it is one actomyosin network responding to
seeded physiological manifolds, nucleator/crosslinker identities, and anchors. `ac` P3/I5 must show bundles
**condense** from local orientations sampled isotropically conditional on those declared inputs (falsifiable
emergence gate); scripted `weave(bundle)` is demoted to a non-authoritative control.

---

## 2. Load-path map — the crux of the walkthrough

### 2.1 The correct wiring (what `ac` must produce)

```
                    ┌─ SF (ventral) ── FA ── ECM ── FA ── SF        traction = closed substrate loop
   NMII contraction ┤
                    ├─ SF (cap) ── LINC(nesprin/SUN) ── NUCLEUS     dedicated, BYPASSES cortex
                    └─ cortex ── ERM ── MEMBRANE                    membrane = passenger (slaved)

   IF cage ── LINC ── NUCLEUS                                       2nd nucleus path (vimentin/keratin)
   MT aster ─(κ_bend)─ compression ── cortex                       tensegrity balance
```

Two facts we nailed that the code must respect:

1. **Force reaches the nucleus by a DEDICATED path (cap→LINC and/or IF cage), NOT by cortex-network
   shear diffusing inward.** Cortex participates only as the *outer anchor* of those paths. FF's IF
   cage already routes it this way (correct topology) but the spoke mechanics FAILED — I7 replaces
   the *driver* with a contractile perinuclear actin cap (capstan pressure), IF demoted to a
   possible secondary passive net.
2. **SF and cortex are largely SEPARATE modules.** There is no global cortex↔SF load bridge; they
   meet only weakly at the peripheral lamella (arc↔dorsal-SF). So "cortex → nucleus" and "cortex
   carries SF traction" are both **wrong** pictures — the code must keep them as distinct assemblies
   joined only by named bridges (P5).

### 2.2 Current FF wiring (what's actually running)

```
   clutch ── CORTEX basal node          ← traction driver is the CORTEX, not SF (SF unwired)
   IF cage(radial spoke) ── LINC ── nucleus   ← only nucleus path present; FAILED mechanics
   membrane ── ERM ── cortex            ← ✓ correct (passenger)
   MT tip ── cortex contact             ← ✓ but block-diagonal
```

The single biggest structural gap: **stress fibers are not in the whole-cell loop at all** (they
exist only as a `weave.py` builder output used by viz/metrics/tests). Consequences: no SF→FA→ECM
traction, no SF-cap→LINC→nucleus, and traction is currently generated by the **cortex shell gripping
the substrate** — a defensible-but-approximate stand-in for the SF loop.

---

## 3. Gap register (each → increment, status)

| # | Gap | Why it matters | Current | Target | Maps to |
|---|---|---|---|---|---|
| G1 | **SF not wired into whole-cell loop** | no SF traction, no cap→nucleus; traction faked by cortex | standalone spec | ventral SF → FA → ECM in the loop | **I6** |
| G2 | **SF-style (myosin-dipole) traction absent** | only lamellipodial clutch traction exists; SF is the mechanosensing/durotaxis channel | — | mature-FA contractile dipole | **I6** |
| G3 | **Perinuclear actin-cap → LINC → nucleus absent** | nucleus flatten driver; currently only (failed) IF cage | IF cage only | contractile cap, capstan pressure | **I7** |
| G4 | **Lumped `f_myo` contraction** | contraction must be NMII, not a constant; cortex AND SF reuse one motor | scalar `f_myo` / aggregate two-anchor control | Warp-CUDA head-resolved backbone + individual-head NMII with **Hill FV**; remaining head constants ratified at I0-B3 | **I3** |
| G5 | **Cortex/SF do not emerge** | scripted bundles = lumped category | `weave(bundle)` scripted | condense from isotropic seed (falsifiable) | **I4/I5** |
| G6 | **Actin–MT block-diagonal** | tensegrity load transfer incomplete | tip-contact only | plectin/MACF actin↔MT crosslink | **I8** |
| G7 | **FA maturation (nascent→mature) not driving tip/shaft/basal** | tip=nascent vs basal=mature is the real difference | uniform clutch | force-dependent clustering/reinforcement | **I6** (`fa_maturation.py`) |
| G8 | **Adhesion classes hard-coding risk** | must NOT hard-code tip/shaft/basal types | one basal clutch | ONE clutch; classes emerge from architecture+load+maturation | §5.2 `⭐NEW` |
| G9a | **Filopodium + tip adhesion not in loop** | tip sensing → durotaxis/guidance/invasion (MDA) | spec only | fascin bundle + tip clutch, in the unified net | **I4/I6** (PI: in scope) |
| G9b | **Lamellipodium is a ratchet proxy, not a real branched net** | current "lamellipodium" = directed barbed-end growth on CORTEX fibers; no dendritic Arp2/3 network | body-force/ratchet proxy | **REBUILD**: dendritic Arp2/3 net (evidence-bound `θ₀`, angle-harmonic + thermal), sheet protrusion | **I4** rebuild (PI: 갈아엎기) `⭐NEW` |
| G10 | **IF keratin/vimentin undifferentiated** | epithelial vs mesenchymal = MCF7 vs MDA invasion | single spoke form | cell-type IF (anchor topology + extensibility) | §5.1 `⭐NEW` — **first-validation scope** (PI) |
| G11 | **Integrin force-scale conflation** | ~7 pN catch and ~30 pN slip scales may be distinct parameters in one law; neither is a universal ligand constant | current parameters incompletely labeled | retain distinct force scales and audit the selected ligand-specific law; do not relabel α5β1–FN as α2β1–collagen | **I0-A baseline + I0-B6** |

---

## 4. Design invariants — what we got RIGHT (must not regress in `ac`)

These are validated in the walkthrough; the new engine must **preserve** them, not "improve" them
into something wrong.

- **INV-1 — Membrane is a passenger, not a load-bearer.** Traction crosses the FA through the
  **transmembrane integrin protein** (ECM→integrin→talin→actin), NOT through the lipid bilayer
  (2D fluid, in-plane shear modulus ≈ 0 → parallel-spring share ≈ 0). Encode: the FA clutch is a
  spring **actin-node ↔ ECM-anchor** with **no membrane element in it**; the membrane is a separate
  Helfrich sheet tethered to cortex by ERM. FF does this — keep it. Falsification check: if lipid
  bore traction it would rupture at ~mN/m; cells sustain ~kPa, so protein clutches bear it.
- **INV-2 — De-adhesion is emergent catch-slip, never scripted.** Rear release = the Pereverzev
  two-pathway off-rate `k_c0·e^{−Fx_c/kT} + k_s0·e^{+Fx_s/kT}` driving a KMC unbind: high load →
  slip term dominates → rupture, then nascent re-bind ahead → adhesion "treadmills" by load-and-fail.
  No `if rear: detach`. ([[feedback-junction-switch-fine-grained]]) Keep `clutch_catchslip_kmc`.
- **INV-3 — MT compression is emergent from bending stiffness**, not a lumped strut
  ([[project-tensegrity-mt-compression-actin-tension]]). Keep `microtubule.py` κ-based mechanics.
- **INV-4 — ONE unified clutch; tip/shaft/basal are OUTPUTS.** The cell knows only the integrin
  molecule. Adhesion classes differ by (architecture attached to) × (local load) × (maturation) —
  none of which is a hand-set class label. (see §5.2)
- **INV-5 — Cortex vs SF: same molecule, architecture emerges.** Distinguished by crosslinker
  (filamin net vs α-actinin bundle) and myosin arrangement (dispersed→isotropic γ vs
  sarcomeric→directional), both emergent under P3.

---

## 5. Walkthrough additions (⭐ now folded into the master plan)

### 5.1 ⭐ Keratin/vimentin IF differentiation (epithelial ↔ mesenchymal)

The build plan demotes the IF cage as the nucleus *driver* (I7 → actin cap) and notes "IF may return
later as a secondary passive network." **When it returns, it must be cell-type-differentiated** — this
is the mechanical substrate of the epithelial↔mesenchymal split and directly gates the MCF7 vs
MDA-MB-231 program ([[project-mda-mb-231-validation-scope]]).

| | Keratin (epithelial, MCF7) | Vimentin (mesenchymal, MDA-231) |
|---|---|---|
| **Anchor topology** | desmosome (cell–cell) + hemidesmosome | perinuclear → FA/periphery (no cell–cell) |
| **Network** | tissue-scale continuum | cell-autonomous |
| **Extensibility** | denser/stiffer, static | more extensible/dynamic → large reversible strain |
| **Function delivered** | collective epithelial integrity | migration + **nuclear protection during confined invasion** |

- **Not the EMT switch** — the causal driver is the EMT transcriptional program (Snail/Zeb/Twist)
  and the cadherin switch (E→N); IF is a downstream **marker + mechanical effector**. `ac` should
  treat the IF class as a **cell-type input (seed)**, not something that flips itself.
- **Minimal `ac` encoding:** IF as a secondary passive network whose (a) anchor set (desmosome-like
  cell–cell vs FA/perinuclear) and (b) extensibility/strain-stiffening params are chosen by
  `cell_type.py` preset. This is what lets MDA squeeze/invade while MCF7 stays collective — a
  cell-type-differentiated load path, not one universal spoke.
- **Blocked on:** IF material params (vimentin vs keratin persistence length, strain-stiffening,
  areal density) → PI-authored KnowledgeClaim, do not tune.
- **Scope (PI 2026-07-16): FIRST-VALIDATION, not post-I9.** `cell_type.py` gets keratin/vimentin IF
  presets early so the MCF7↔MDA invasion contrast is exercised from the first native `ac` run. IF
  stays a *secondary passive net* (not the nucleus driver — that is the I7 actin cap), but it is
  present and cell-type-differentiated in the foundation gate. Sequence: land after I7's cap driver
  so the combined cap+IF nucleus baseline is gated together (double-load-path, build-plan §I7 warn).

### 5.2 ⭐ Adhesion-class emergence principle (tip / shaft / basal)

**Do NOT hard-code three adhesion classes** — that reintroduces the lumped-category trap. One
integrin catch-slip clutch, placed on the real architectures, and the classes fall out:

| Class | Emerges because it sits on… | …at load | …maturation |
|---|---|---|---|
| **tip** | filopodium bundle end | low, nascent | immature |
| **shaft** | lamellipodium net / bundle flank | mid | maturing |
| **basal** | cortex / SF end | high, sustained | mature FA |

All three axes already exist as mechanisms: architecture (unified network, P3), load (catch-slip is
load-dependent — I6), maturation (`fa_maturation.py`, force-dependent clustering — the DCM
cadherin-cluster pattern, [[project-dcm-cadherin-cluster-redesign]]).

**Scope (PI 2026-07-16): ALL THREE in the first validation — "전부 다 넣고".** One clutch + full
architecture set (cortex/SF + **filopodium/tip** G9a + **rebuilt lamellipodium** G9b) + maturation.
tip/shaft/basal must all be exercised and shown to be EMERGENT (not hand-set classes). This
supersedes the build-plan increment-2 deferral of filopodium/lamellipodium — they move into the I4
region set + I6 adhesion set (see §6). The invariant holds: still ONE clutch, classes are outputs.

### 5.3 ⭐ Two-traction-channel acceptance gate

Traction has two mechanistically distinct sources. **Scope (PI 2026-07-16): the I6 gate REQUIRES
BOTH channels present** ("2-채널 동시 요구") — not SF-only with lamellipodial deferred. The gate
checks the taxonomy (both sources distinguishable), not just a lumped traction magnitude:

| | Lamellipodial | SF |
|---|---|---|
| **Engine** | actin polymerization → retrograde flow (ratchet) | NMII sarcomeric contraction |
| **Adhesion** | nascent | mature FA |
| **Magnitude** | small, distributed (per-clutch ~10–25 pN) | large, concentrated (single SF ~5–6 nN) |
| **Direction** | centripetal vs retrograde flow | inward contractile dipole (FA↔FA) |
| **Function** | propulsion | force generation + mechanosensing (durotaxis) |
| **Status** | `~` ratchet proxy on cortex → **rebuild as branched net (G9b)** | `✗` → I6 |

Numbers are KB-anchored ([[project-sf-nmii-forcescale-result]] single SF ~5–6 nN; 10–30 nN is
network/prestress — do **not** back-solve N_minifil to hit 10–30 nN, geometrically impossible per
build-plan §4C) — **report, don't tune** (density-floored, per build-plan honest-truth #2).

---

## 6. Sequencing (fits the existing I0..I9 spine, with the PI-2026-07-16 scope expansion)

The PI decisions below **expand** the first-validation scope beyond the build-plan's increment-2
deferrals: filopodium + lamellipodium-rebuild + cell-type IF all move into the foundation gate.

- **I0-A — RATIFIED:** Warp CUDA only; no HOOMD execution; head-resolved NMII+Hill; MCF7×collagen×α2β1;
  70,686 cortical F-actin. At I0-B6, audit the selected ligand-specific clutch law and keep its catch/slip force scales
  distinct; **do not** choose 7 vs 30 pN as if they were one universal `F*`. Add the **G9b Arp2/3 branch-angle** knob
  (evidence-bound `θ₀`, angle-harmonic stiffness + thermal σ) and **filopodium fascin bundle** params to the gate.
- **Population budget:** `70,686` is cortex only. Lamellipodium adds `∫ρ_barbed dA`; separate SF/filopodial
  bundles add `ΣN_fil,b`; MT/IF are separate systems. A filament already classified from the unified actin
  pool is not added again. Report active unique IDs, dormant capacity, nodes, explicit motor heads/states, and
  exact peak GPU bytes before every native gate.
- **I4 (unified spine) — EXPANDED region set:** the woven network now seeds **cortex + ventral/dorsal
  SF + transverse arc + perinuclear cap + filopodium (fascin bundle) + lamellipodium (dendritic
  Arp2/3 net)**. G9b is a genuine build, not a port: the current "lamellipodium" is directed
  barbed-end growth on cortex fibers — **rebuild** it as a dendritic network (mother/daughter
  topology, Arp2/3 branch as **angle-harmonic + thermal fluctuation**, per CLAUDE.md). Seed manifolds and
  NPF/nucleator fields are declared inputs; local orientations are isotropic conditional on them, and
  label-blind bundles/branches condense (I5 emergence gate).
- **I6 (adhesion baseline) — EXPANDED:** carry **INV-1, INV-4, G7, G8, §5.2, §5.3**. Unified clutch
  + maturation on the full architecture set; tip/shaft/basal exercised and shown emergent; **the
  gate REQUIRES BOTH traction channels** (lamellipodial + SF), per §5.3.
- **I7 (nucleus load path):** **G3** cap→LINC→nucleus driver; then land **§5.1** IF as the
  cell-type-differentiated *secondary* net, gated **together** with the cap (combined cap+IF
  baseline, double-load-path).
- **Cross-doc reconciliation — DONE (2026-07-16):** the scope expansion is folded into the master
  `NEW_ENGINE_BUILD_PLAN_2026-07-16.md` **§3A** (with the lamellipodium-rebuild design (b) and the
  cell-type IF preset spec (c) written out in full there). The master P3, I4/I6/I7 rows, deferral
  list, §2b IF kill-row, and §4C params are updated. This companion doc is the summary; §3A is the
  authoritative build detail.

---

## 7. PI scope decisions — 2026-07-16 (G11 interpretation later corrected)

The §6 scope questions were put to the PI and resolved. Item 4's original interpretation was later
superseded on technical audit; it now returns to the I0-A/I0-B6 evidence gate:

1. **Cell-type IF (§5.1, MCF7 vs MDA)** → **FIRST-VALIDATION scope.** `cell_type.py` gets
   keratin/vimentin presets early; MCF7↔MDA invasion contrast exercised from the first native run.
2. **Adhesion scope (§5.2 / G9)** → **"전부 다 넣고, 라멜리포디움도 새로 갈아엎기."** Full tip/shaft/basal
   via filopodium + cortex + SF; **lamellipodium REBUILT** from the ratchet proxy into a real
   dendritic Arp2/3 branched network (G9b). Not basal-only.
3. **Traction gate (§5.3)** → **BOTH channels required at I6** (lamellipodial + SF), not SF-only.
4. **G11 correction (implementation-readiness audit):** the earlier binary "7 vs 30 pN" resolution is
   **superseded** because those numbers can represent distinct catch/slip scales within one law; ligand
   mismatch was not established by the values alone. The selected I0-A cell×ECM×ligand baseline controls
   which audited clutch law I0-B6 binds.
5. **I0-A runtime/motor/baseline (PI 2026-07-16):** Warp-CUDA-only and zero HOOMD execution; complete
   backbone+individual-head NMII port with Hill FV; MCF7×collagen×α2β1; 70,686 cortical F-actin plus
   independently derived non-cortical populations.

## Change log
- 2026-07-16: created from the filament walkthrough (cortex/SF/MT/IF distinctions, FA/traction load
  paths, membrane-not-load-bearer, de-adhesion catch-slip, adhesion-class emergence). Companion to
  `NEW_ENGINE_BUILD_PLAN_2026-07-16.md`; most gaps map to existing increments, three items are
  net-new (§5). Verdicts on force magnitudes are HELD (density-floored findings, not deliverables).
- 2026-07-16 (same day, PI scope expansion §7): cell-type IF (MCF7/MDA) → first-validation;
  filopodium+tip AND lamellipodium-rebuild (dendritic Arp2/3, G9b) → in scope, not I9+ deferred;
  I6 traction gate requires BOTH channels. Flagged the master build-plan increment-2 deferral list
  for reconciliation (§6).
- 2026-07-16 implementation-readiness audit: reconciled the promoted I4/I6 region rows, replaced aggregate
  linear-NMII language with the head-resolved I3 contract, corrected the integrin force-scale interpretation,
  and made seeded physiological inputs compatible with label-blind emergence.
- 2026-07-16 I0-A ratification: bound Warp-CUDA-only execution, head-resolved Hill NMII, the
  MCF7×collagen×α2β1 baseline, 70,686 cortical filaments, and a non-double-counting whole-cell population/
  device-memory ledger.

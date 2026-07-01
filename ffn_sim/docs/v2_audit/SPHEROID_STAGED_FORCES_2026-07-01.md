# Spheroid formation as staged maturation of ONE cell-cell junction (2026-07-01)

**Origin:** PI observation that spheroid formation has *time-varying dominant forces* → define
aggregation / compaction / maturation and think about which force dominates each, for the simulation.
**Method:** literature-anchored 3-stage framework + **controlled DCM experiments this session** (no
parameter tuned to outcome — forces toggled at lit-anchored values, outcomes reported). Figures under
`outputs/h_dcm_two_stage/figs/staged_*.png`.

---

## 0. The reframe (PI 2026-07-01) — ONE junction, not three separate forces

Cell-cell **contact = adhesion = interfacial tension = compaction** are NOT distinct lumped forces —
they are all expressions of **one cell-cell junction** (cadherin adherens junction + its actomyosin
coupling). The junction-mediated nature is **invariant across stages**; what changes stage-to-stage is
the junction's **properties** (reach, adhesion strength, active contraction, interfacial tension), not
the mechanism. (Memory: `feedback-cell-cell-contact-junction-mediated`; generalises the fine-grained
junction-switch rule.) So the staged model is **one junction whose properties mature over time**, and
the current gaps are **specific missing junction properties**.

---

## 1. Literature — the canonical 3-stage kinetics (multi-source verified)

According to PubMed + reviews, spheroid formation is a 3-stage sequence with a stage-specific dominant
junction property:

| stage | dominant junction property | mechanism | source |
|---|---|---|---|
| **1. Aggregation** | **long-range reach** | integrin–ECM (fibronectin) fibres bridge separated cells | Robinson-Foty 2003 [10.1242/jcs.00231]; Lin 2006 [10.1007/s00441-005-0148-2] |
| **2. Compaction** | **adhesion strength + active contraction** | cadherin homophilic binding accumulates, then **actomyosin (RhoA/ROCK)** pulls bonded junctions in → loose→dense | Smyrek 2019 [10.1242/bio.037051]; Lin 2006 |
| **3. Maturation** | **interfacial tension σ** | tissue surface tension ∝ cadherin (DAH); liquid-drop faceting/rounding | Foty-Steinberg 2005 [10.1016/j.ydbio.2004.11.012]; Winters-Foty 2005 [10.1002/ijc.20722] |

Timescale (microwell): ~10 h loose coalescence → ~1 day compaction → days maturation.
⚠️ aggregate surface tension σ (Foty, several mN/m) ≠ single-cell cortical γ (the aggregate-vs-single-cell
confound already tracked in the KB).

---

## 2. Empirical junction-property map — what the DCM engine has (this session, no tuning)

Controlled experiments, N=32–64, CPU, forces at lit-anchored values (cadherin bundle 40 → ~7 nN/junction
∈ KB-4.11; γ swept as controlled variable; w_cs = 2.85 mJ/m² MCF7):

| junction property (stage) | model status | evidence |
|---|---|---|
| **long-range reach** (1) | ❌ **MISSING** | gap 2.35 (surfaces ~2.6 µm): cadherin + coupling both ~1.8 µm range → **bonds=0, no bridging**. `ecm-clutch` is basal→substrate only, not cell-cell. |
| **adhesion strength** (2) | ✅ present | cadherin OFF vs ON: contact fraction 0.10→**0.17** (1.7×), faster gap closure. Bonds form + catch-slip churn. |
| **active contraction** (2) | 🔧 **BUILT this session; compaction turgor-BLOCKED** | Added `f_contract` (active actomyosin junctional contraction) to the cadherin kernel/host/driver (§2b). Swept 0→6 nN/junction: it closes inter-cell gaps (0.82→0.00 µm, the passive bond couldn't) but **Rg barely moves (−0.2 % even at 6 nN ≫ lit) = still NO aggregate compaction**. Turgor fixes cell volumes → pulling harder locally interpenetrates rather than densifies. |
| **interfacial tension — uniform** (3) | ⚠️ wrong sign | γ sweep {0,1e-3,1e-2,1e-1}: asph **DROPS** 0.0069→0.0013 (uniform bulk γ minimises area → **rounds** cells), contact drops, Rg flat. |
| **interfacial tension — DIFFERENTIAL/junction** (3) | ✅ **works** | junction-localised differential tension (contact-face γ↓, Young-Dupré): asph **RISES** to **0.033** at γ=1e-2 (~10× uniform), rises with γ → **emergent faceting** = the DAH mechanism, junction-mediated. |

**Money figure** `staged_faceting_junction_mechanism.png`: uniform γ rounds (blue, ↓), differential
junction tension facets (red, ↑ above the 0.03 threshold). Directly confirms the §0 principle —
faceting comes from the **junction** (differential contact-face tension), not a bulk force.

## 2b. Built the Stage-2 compaction MOTOR (active junction contraction) — and it's turgor-blocked

To test whether the missing "active contraction" property is *the* compaction blocker, it was **built
this session** as a property of the same cadherin junction (memory: one junction): a new `f_contract`
term in `cadherin_bond_force_kernel` (+ `CadherinParams.f_contract`, `run_decohesion(cad_contract=…)`,
`--cad-contract`) — an actomyosin junctional contraction that **always** pulls bonded cells together
(RhoA/ROCK-gated NMII), on top of the passive catch-bond. Anchored to junctional actomyosin tension
(per-motor ~5-15 pN × engaged), swept as a controlled variable, never tuned to a compaction target.

**Sweep result** (`staged_contraction_motor.png`, N=32 gap 2.2, junction contraction 0→6 nN): the motor
**works at the force level** — it closes inter-cell gaps 0.82→0.00 µm (the passive bond, capped at rest
length, could not) and gives a monotonic Rg trend. **But aggregate compaction is negligible: Rg −0.2 %
even at 6 nN (≫ the lit sub-nN–nN junctional range).** Turgor fixes each cell's volume, so pulling the
junction harder **locally interpenetrates** instead of densifying the aggregate. (⚠️ the gap-2.2 strong
bundle-40 + explicit integrator already fails the interpenetration gate at 0 nN — a known pre-existing
issue needing the implicit integrator, NOT the motor; the Rg-flat compaction verdict is independent of it.)

**Verdict:** building the correct motor proved compaction is **structurally turgor-blocked**, not merely
motor-missing — the same force-magnitude pattern as the γ-floor / SF-force lines. The real missing piece
for compaction is **volume accommodation** (cells reducing volume / porosity closing), which turgor
currently prevents — not the junction motor.

---

## 3. Reading — why the current model can't do dynamic assembly (and why confluent-init draws it)

- **Faceting (maturation) IS mechanistically reachable**: differential junction tension at the lit MCF7 γ
  produces emergent faceting (asph 0.033) — a genuine positive. The confluent-init draws asph 0.04–0.08
  *geometrically*; the junction physics can produce ~0.033 of that *mechanistically*.
- **Compaction (aggregation→dense) is NOT reachable**: **Rg is flat in EVERY experiment** (cadherin,
  uniform γ, differential γ). The aggregate never contracts — the **active contraction** junction property
  is absent, and turgor holds each cell's volume. Compaction needs the actomyosin motor, which we don't have.
- **Aggregation (separated→loose) is NOT reachable**: no long-range junction reach; separated cells can't
  find each other.
- ⇒ The two missing junction properties (**long-range reach**, **active contraction**) are exactly the
  *assembly* drivers. Their absence is **why no dynamic round→compact trajectory exists** and why the
  confluent-init has to supply the assembled geometry by construction. Same force-magnitude/deficit pattern
  as the γ-floor and SF-force lines: the fine-grained force is present at its lit value, but the outcome
  doesn't emerge because the dominant driver is missing or out-powered.

---

## 4. Design — the unified junction with staged, maturing properties (PI-gated)

Per §0, build ONE junction entity (cadherin catch-bond core, already present) and add the two missing
properties as **modulated junction properties**, not separate mechanisms:

1. **Stage-1 reach** — an ECM-fibre tether (integrin-anchored, ~µm reach, fibronectin length-scale) as the
   junction's long-range precursor state. Mechanistic (fibre + integrin clutch), not a lumped attractor.
2. **Stage-2 active contraction** — actomyosin-coupled **junction-shortening** (RhoA/ROCK-gated) — the
   compaction motor. This is active mechanobiology (a motor), i.e. **FF-engine physics** grafted onto the
   junction; the honest reason compaction fails today.
3. **Stage-3 differential interfacial tension** — already works (`--diff-tension`, Young-Dupré from w_cs);
   promote to the production faceting driver at the lit γ (controlled variable, never tuned to a target).

Build order (per the junction-fine-grained rule): assemble ALL junction properties fine-grained → one
production test, not piecemeal.

---

## 5. Per-stage gates (lit-anchored observables — write before running)

- **Stage-1:** porosity high→dropping, aggregate forms from separated cells (Lin ~10 h kinetics, order only).
- **Stage-2:** porosity → 2–12 % (T47D), contact fraction ↑, **Rg drops** (the compaction signature).
- **Stage-3:** faceting Q / asph rises via differential junction tension; σ force- & volume-independent
  (Foty liquid-drop test).

---

## Figures (`outputs/h_dcm_two_stage/figs/`)

- `staged_stage2_cadherin_compaction.png` — cadherin adhesion engages (contact 1.7×) but no compaction (Rg flat).
- `staged_gamma_sweep.png` — uniform γ rounds cells / no compaction across the lit band.
- `staged_faceting_junction_mechanism.png` — **the money figure**: differential *junction* tension facets,
  uniform bulk tension rounds → faceting is junction-mediated.
- `staged_contraction_motor.png` — the BUILT active junction contraction motor: closes gaps but Rg −0.2 %
  even at 6 nN ≫ lit → compaction is turgor-blocked, not motor-missing.

---

## PI surface (gate items)

1. The two missing junction properties — **long-range reach** and **active contraction** — are new
   mechanistic builds (active mechanobiology), PI-gated; the second is FF-territory grafted onto the junction.
2. **Compaction is structurally turgor-blocked** — the active contraction motor was BUILT (§2b) and,
   even at 6 nN ≫ lit, gives Rg −0.2 % (no compaction); it closes gaps but turgor's volume-fixing prevents
   densification. So the round→compact trajectory the PI asked to see is missing NOT for lack of a motor but
   because **volume accommodation** is absent. Candidate next lever (PI): allow osmotic volume regulation /
   turgor relaxation during compaction (a real cell process), then re-test the motor.
3. **Faceting IS reachable** mechanistically via differential junction tension (asph 0.033 at MCF7 γ) — a
   candidate to replace/validate the geometric confluent-init faceting.
4. No parameter was tuned to outcome; γ swept as a controlled variable; forces at lit-anchored values.

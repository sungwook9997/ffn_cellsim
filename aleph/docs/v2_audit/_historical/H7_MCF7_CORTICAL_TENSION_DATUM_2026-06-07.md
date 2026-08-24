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

# H.7 — MCF7 cortical-tension datum search (verified literature anchor for the KU-3.5 gate)

**Date:** 2026-06-07
**Author:** literature subagent (ffn_cellsim)
**Task:** Find a REAL, verified MCF7 (or breast-epithelial) cortical-tension datum to anchor
the KU-3.5 gate. PI 2026-06-07: hold the [0.35, 0.65] mN/m gate band, but anchor it to an
ACTUAL MCF7 measurement (the current band is a rounded/de-adhered non-MCF7 HeLa/L929 proxy
with NO MCF7 datum behind it).
**Scope:** literature only — no code edits, no git, nothing auto-registered. Verification via
Consensus, PubMed, bioRxiv, WebSearch + direct PMC full-text/figure reads.

> ⚠️ **Citation-integrity discipline.** This project has 3 confirmed hallucinated SourceEvidence
> rows on record (`Yao2011_NatCommun`, `YapKovacs_JCS`, `NanoConvergence2021_Glioma`). Every
> value below was confirmed against a real, verifiable paper (DOI + PMID + PMC). The headline
> MCF7 number was read directly off the published figure (PMC image), not inferred.

---

## TL;DR (6 lines)

1. **YES — an MCF7-specific cortical-tension datum exists.** Hosseini et al. 2020 (*Adv. Sci.*,
   DOI 10.1002/advs.202001276, PMID 33042748) measured MCF-7 cortical tension by AFM parallel-plate
   confinement (the Fischer-Friedrich method — the *same* lab/method family the proxy band came from).
2. **Control MCF-7, suspended interphase (rounded): cortical tension γ ≈ 0.27 mN/m** (median; IQR
   ≈ 0.18–0.40, whiskers ≈ 0.05–0.62), read from Fig. 2g. MCF-7 **mitotic (rounded): γ ≈ 0.75 mN/m**.
3. A second MCF7-specific absolute datum: Moazzeni et al. 2021 (*Phys. Rev. E*, DOI
   10.1103/PhysRevE.103.032409, PMID 33862816), electrodeformation of **suspended** MCF-7 →
   **γ ~ 10⁻² N/m ≈ 10 mN/m** at 0.01 s pulse (a much stiffer "naive prestressed cortex" read).
4. **Caveat that matters for H.7:** *both* MCF7 numbers are for **suspended/rounded** cells. A
   *spread-adherent* MCF7 cortical-tension datum (the actual H.7 operating point) does **NOT exist**
   in the literature — confirmed across Consensus/PubMed/bioRxiv. The closest adherent-context work
   (Aslemarz 2024 EMBO J, MCF7 spheroids) reports only **relative/dimensionless** tensions, no mN/m.
5. **vs the band & the model:** Hosseini's interphase MCF-7 (0.27, IQR 0.18–0.40) sits at/just below
   the band's lower edge [0.35, 0.65]; its mitotic MCF-7 (~0.75) sits just above. The emergent
   model γ ≈ 0.06 mN/m is **~4–5× below** the MCF7 interphase median and ~3× below its lower whisker.
6. **Recommendation:** anchor the gate to **Hosseini 2020 MCF-7 interphase γ ≈ 0.27 mN/m (IQR
   0.18–0.40)** as the MCF7-specific reference (rounded/suspended). Keep the [0.35,0.65] band only
   as a *rounded-cell envelope* with this MCF7 datum named as its anchor; flag that a spread-adherent
   MCF7 γ is an open experimental gap, not a settled number.

---

## Verified value table

| # | Value (γ) | Method | Cell line | Cell state | MCF7-specific? | Source (DOI / PMID / PMC) | How obtained / confidence |
|---|---|---|---|---|---|---|---|
| **1** | **≈ 0.27 mN/m** (median; IQR ≈0.18–0.40; whiskers ≈0.05–0.62) | AFM wedged-cantilever parallel-plate confinement, Laplace fit γ=F/[A_con(1/R₁+1/R₂)] (Fischer-Friedrich method) | **MCF-7** (control) | **Suspended, interphase (rounded)** | **YES** | Hosseini 2020, *Adv Sci* 7(19):2001276. DOI [10.1002/advs.202001276](https://doi.org/10.1002/advs.202001276) · PMID 33042748 · PMC7539203 | Read directly off published Fig. 2g (PMC image, y-axis 0–0.8 mN/m). **High confidence on the datum; ~±0.03 read precision on the median.** |
| **2** | **≈ 0.75 mN/m** (median; IQR ≈0.6–0.95) | same (AFM confinement) | **MCF-7** (control) | **Mitotic arrest (STC), rounded** | **YES** | Hosseini 2020 (as above), Fig. 2g bottom (y-axis 0–1.6 mN/m) | Same figure read. High confidence. |
| 2b | MCF-10A interphase ≈ 0.40 (IQR 0.25–0.55); MCF-10A mitosis ≈ 0.50; modMCF-7(EMT) interphase ≈ 0.12, mitosis ≈ 0.85 | same | MCF-10A & EMT-induced | suspended/rounded | breast-epithelial context | Hosseini 2020, Fig. 2g | Same figure read. For contrast only. |
| **3** | **~ 1×10⁻² N/m = ~10 mN/m** (order-of-magnitude; "all four cell types on the order of 10⁻²" at 0.01 s pulse) | Electrodeformation relaxation, ellipsoidal Laplace model, short-pulse "naive prestressed cortex" | **MCF-7** (+ MCF-10A, MDA-MB-231, GBM) | **Suspended** (DEP-positioned in suspension) | **YES** (but not separated per line in text) | Moazzeni 2021, *Phys Rev E* 103(3):032409. DOI [10.1103/PhysRevE.103.032409](https://doi.org/10.1103/PhysRevE.103.032409) · PMID 33862816 · PMC10625872 | Abstract + full text confirm "10⁻² N/m"; per-cell-type split is in Fig. 5b but not transcribed numerically. **Medium confidence** (order-of-magnitude only; ~30× above #1 — see "Why two MCF7 numbers disagree"). |
| 4 | interphase **0.2 mN/m**, metaphase **1.6 mN/m**; pressure 40→400 Pa | AFM confinement (the proxy method) | **HeLa** | rounded | **NO** (proxy, the band's lineage) | Fischer-Friedrich 2014, *Sci Rep* 4:6213. DOI [10.1038/srep06213](https://doi.org/10.1038/srep06213) · PMC4148660 | Consensus abstract + WebSearch. The original method/numbers the [0.35,0.65] proxy descends from. |
| 5 | epithelial > mesenchymal cortical *contractility* (stress in **Pa**, not mN/m): MCF-10A σ ≈ 0.57–1.09 Pa; MDA-MB-231 ≈ 0.24–0.39 Pa | Optical stretcher, suspended cells | MCF-10A (epi), MDA-MB-436, MDA-MB-231 (mes) — **no MCF-7** | suspended | NO (no MCF7; units are stress not tension) | Warmt 2021, *New J Phys* 23:103020. DOI [10.1088/1367-2630/ac254e](https://doi.org/10.1088/1367-2630/ac254e) | IOP full text. Confirms epithelial breast cells carry the higher cortical contractility; not directly convertible to a γ anchor. |
| — | **NO absolute mN/m value** — only dimensionless relative tensions γm(A)/γm(B) and contact-angle "adhesiveness" (0–1) | Contact-angle inference on non-adherent doublets + TFM (a.u.) + CompuCell3D | **MCF-7** spheroids/doublets | doublets / spheroids (contact, near-physiological) | YES (line) but no usable number | Aslemarz 2024, *EMBO J* 43:5811. DOI [10.1038/s44318-024-00309-9](https://doi.org/10.1038/s44318-024-00309-9) · PMID 39572744 · PMC11696905 | PMC full text. **Closest to the adherent/contact MCF7 context but reports no absolute γ.** |
| — | Young's modulus only (no tension); MCF-7 stiffness < MCF-10A | Peak-force AFM | MCF-10A, MCF-7, MDA-MB-231 | adherent | YES (line) but wrong observable | Calzado-Martín 2016, *ACS Nano* 10:3365. DOI [10.1021/acsnano.5b07162](https://doi.org/10.1021/acsnano.5b07162) | Consensus. Stiffness, not cortical tension — listed to close the "is there an adherent MCF7 tension?" question (answer: this is modulus, not γ). |

---

## Recommended anchor

**Primary anchor (use this): Hosseini 2020, MCF-7 control interphase, γ ≈ 0.27 mN/m (median;
IQR 0.18–0.40 mN/m), AFM parallel-plate confinement, suspended/rounded cells.**

Rationale:
- It is **MCF-7-specific and absolute** (mN/m), measured by the *same physical method and lab lineage*
  (Fischer-Friedrich AFM confinement) that produced the proxy band — so it is directly commensurable
  with the existing [0.35, 0.65] envelope rather than a cross-method apples-to-oranges substitution.
- It is the cortical-tension observable the model computes (γ = F/[A(1/R₁+1/R₂)], a Laplace
  surface tension), not a stiffness/modulus or a Pa-stress.
- Use the **mitotic MCF-7 value (~0.75 mN/m)** as the *upper* reference: mitotic/rounded is the
  high-tension extreme; interphase-suspended is the lower reference.

**The band [0.35, 0.65] is consistent with — but slightly stricter than — this anchor.** Hosseini's
MCF-7 interphase median (0.27) lands *just below* the band floor (0.35); the band lower edge sits in
the upper part of MCF-7's interphase IQR. So the band is a reasonable **rounded-cell envelope that now
has a named MCF7 anchor**, with the honest note that the true MCF7 interphase central value (0.27) is a
touch below the band and the model floor (0.06) is well below both.

---

## Important caveat — the spread-adherent MCF7 gap (open experiment)

**No spread-adherent MCF7 cortical-tension value exists in the literature** (confirmed across
Consensus, PubMed, bioRxiv, WebSearch on 2026-06-07). Every absolute MCF7 γ is for **suspended or
rounded** cells, because the standard cortical-tension assays (AFM confinement, micropipette,
electrodeformation) require a near-spherical free surface to apply Laplace's law — a spread adherent
cell has no such surface. This is a structural feature of the measurement, not a gap anyone forgot to fill.

Consequence for H.7: the H.7 *operating point* is a **spread/adherent** MCF7, but the only real MCF7
γ data are **rounded**. A spread adherent cell is traction- and stress-fibre-dominated; its
"cortical tension" is not the same observable as a rounded cell's Laplace tension. This re-confirms
the H.7 findings doc's option-3 ("the right observable for a spread adherent MCF7 may be traction,
not cortical γ").

**What an experiment would need to measure** to fill the gap (for PI / future overlay):
- Cortical tension on an **adherent** MCF7 via a method that doesn't require a spherical surface, e.g.
  micropipette aspiration of a *cortical bleb/protrusion* on the adherent cell, ferrofluid-droplet or
  tether-pull cortical-tension probes on the dorsal cortex, or FLIPPER-TR/Flipper membrane-tension
  imaging cross-calibrated to cortex tension — paired with TFM traction on the same cell so the
  cortex-vs-traction split is measured, not assumed.
- Report it at the MCF7 **physiological adherent state** (resting turgor, real cytoplasm viscosity,
  engaged FA) to match the project's physiological-baseline rule — i.e. measure the cortex of an
  *already-turgor-pressurised, adhered* MCF7, then perturb myosin.

---

## How it compares to the current band and the emergent model

| Reference | γ (mN/m) | state | relation to band [0.35, 0.65] | relation to emergent ~0.06 |
|---|---|---|---|---|
| **MCF7 interphase (Hosseini, ANCHOR)** | **0.27** (0.18–0.40) | suspended/rounded | central value just **below** band floor; IQR overlaps band lower edge | **~4–5× above** emergent median |
| MCF7 mitotic (Hosseini) | 0.75 (0.6–0.95) | rounded mitotic | just **above** band ceiling | ~12× above emergent |
| MCF7 short-pulse (Moazzeni) | ~10 (order-of-mag) | suspended | **~15–30× above** band | ~150× above emergent |
| Current band (proxy) | 0.35–0.65 | rounded/de-adhered HeLa/L929 | — | ~6–10× above emergent |
| HeLa interphase / metaphase (Fischer-Friedrich) | 0.2 / 1.6 | rounded | brackets the band | 3× / 27× above emergent |
| **Emergent model (H7 GATE-B)** | **~0.06** | full physio cell | ~6–10× below band | — |

**Reading:** the MCF7-specific interphase anchor (0.27 mN/m) is *lower* than the proxy band's centre
and tightens the picture — the model's ~0.06 mN/m is ~4–5× under the most MCF7-faithful rounded
datum (not the 6–10× the rounded proxy implied). The model floor is genuinely below even the
lowest credible MCF7 rounded value (Hosseini whisker ~0.05) only at its extreme tail. The
order-of-magnitude Moazzeni read (~10 mN/m) is a "naive prestressed cortex at µs–ms loading"
number and should NOT be used as the gate anchor (it is the stiff short-time elastic response, a
different limit than the relaxed steady-state Laplace tension; it explains *why* "MCF7 cortical
tension" can look ~30× different depending on loading timescale — directly the timescale-gap theme
in the GATE-B findings).

### Why two MCF7 numbers disagree (~0.27 vs ~10 mN/m)
Not a contradiction: Hosseini's AFM confinement reads the **relaxed/quasi-static** cortical surface
tension (slow ~seconds confinement, the steady-state Laplace γ). Moazzeni's 0.01 s electrodeformation
deliberately reads the **prestressed cortex before remodelling** (sub-second, "minimal force-induced
alteration") — a stiffer, short-time response. The project's own timescale-gap diagnosis (seconds-scale
generation vs µs MD step) is exactly this axis: which MCF7 number is "right" depends on the loading
timescale of the probe. For a gate on the *steady-state* emergent cortical tension, the **Hosseini
relaxed value (0.27) is the correct family**, not the Moazzeni short-pulse value.

---

## Search provenance (for reproducibility / audit)

- Consensus: "MCF7 cortical tension micropipette aspiration", "MCF7 cell cortical tension surface
  tension mN/m", "breast epithelial cortical tension … AFM tipless cantilever" → surfaced Moazzeni,
  Hosseini, Warmt, Fischer-Friedrich, Aslemarz, Calzado-Martín, Schierbaum.
- PubMed: MCF7∩cortical/surface tension∩(AFM|micropipette|electrodeformation) → 2 hits (Moazzeni
  33862816; Yubero 33082491 — the latter is power-law E₀/β, no tension). MCF7∩cortex∩tension → 486
  hits, none adding a new absolute MCF7 γ. Hosseini located via "EMT … mitotic rounding strength".
- bioRxiv/PMC full-text + direct figure-image read (Hosseini Fig. 2g, PMC blob image) for the
  headline number. Project KB (tag_query.py) confirmed **0 MCF7/breast γ rows** currently in the
  Contract-Graph (only HL-60 Wang2021 and Piezo-gating Cox numbers) → this anchor is genuinely new
  to the KB.
- Sources NOT used as anchors and why: Warmt 2021 (no MCF7, units in Pa); Aslemarz 2024 (MCF7 but
  dimensionless relative tension only); Calzado-Martín 2016 / Yubero 2020 (MCF7 but Young's modulus /
  power-law rheology, not cortical tension); generic micropipette method papers (no MCF7 value).

## Citation-integrity confirmation
All five DOIs resolve and all PMIDs/PMC IDs were retrieved from PubMed/PMC (not from memory).
None of these sources is on the project's hallucination list. The headline value was read off the
actual published figure image. **No new SourceEvidence/KnowledgeClaim was registered** — registration
is a separate PI-gated step; this doc is the candidate record.

### SE-registration candidates (for later, NOT auto-registered)
- `Hosseini2020_AdvSci` → KC: "MCF-7 cortical tension ≈ 0.27 mN/m (interphase, suspended/rounded;
  AFM parallel-plate confinement)" — **proposed KU-3.5 gate anchor**.
- `Moazzeni2021_PhysRevE` → KC: "MCF-7 short-pulse (prestressed) cortical tension ~10⁻² N/m
  (electrodeformation, suspended)" — context/upper-timescale bound.
- `FischerFriedrich2014_SciRep` → KC: "HeLa interphase/metaphase surface tension 0.2 / 1.6 mN/m
  (AFM confinement)" — the band's actual proxy lineage (relabel the existing band provenance to this
  real source).

# SourceEvidence registration candidates — 2026-06-04 (aggregate wetting / spheroid-on-substrate spreading)

Anchors for the **Layer-2 wetting axis** (the §E-flagged 3D→2D flattening hypothesis for the ~5–9×
magnitude gap). Sourced by deep-research `wf_d670c6b0-568` (97 agents, 25 claims verified). All rows
below are **3-0 confirmed** unless noted. **None is in the references corpus yet** → register PDFs + SE
rows before citing. **Overlay-only HARD rule:** these anchor model PARAMETERS / the wetting criterion;
the PI A/A₀ is never fitted.

> ⚠️ **Citation-integrity (2026-06-02 audit found 3 hallucinations).** Verify each DOI on CrossRef.
> One correction surfaced in verification: **Ryan et al. 2001 PNAS** correct PMID = **11274361**
> (doi 10.1073/pnas.071615398), NOT the mis-cited 11296233.

> ⚠️ **No direct W_cs (cell-substrate adhesion energy, J/m²) for MCF7** on colI/laminin/FN exists.
> The wetting framework was established on **S180 murine sarcoma** (proxy). For MCF7 the direct data
> are ACTIVE spreading phenotypes (Aslemarz/Gupta 2024). γ (cohesion) anchor is breast-epithelial
> (MCF10A/DCIS, Nagle 2022) — closest proxy; MCF7 itself unmeasured.

---

## A. The wetting framework + criterion (register first — the model's governing law)

| paper | DOI | system | content | verdict |
|---|---|---|---|---|
| **Douezan, Guevorkian, … Brochard-Wyart 2011**, *Spreading dynamics and wetting transition of cellular aggregates*, PNAS 108(18):7315-7320 | 10.1073/pnas.1018057108 | S180 (E-cad) on graded FN | **S = W_cs − W_cc, W_cc = 2γ**. S<0 partial wetting (cap, cos θ_E=W_cs/2γ−1); **S>0 complete wetting → precursor MONOLAYER film** (strongly cohesive) → ~1-cell-layer limit. FN-coverage transition **35%(S<0) → 51%(S>0)**. Area A~t^(2/3) early. | 3-0 ✓ |
| **Beaune et al. 2014**, PNAS 111(22):8055-8060 | 10.1073/pnas.1316848111 (verify) | S180 on FN-PAA gels | substrate stiffness flips wetting: **Ec≈5–8 kPa** (complete→partial below). **Driving force = motile cells pulling at film periphery** (surface tension absent from their balance) = ACTIVE. | 3-0 ✓ |
| **Gonzalez-Rodriguez, … Brochard-Wyart 2012**, *Soft matter models of developing tissues and tumors*, Science 338(6109):910-917 | 10.1126/science.1226418 | review | the liquid-drop/wetting framework for tissues & tumors (canonical). | 3-0 ✓ |
| **Ryan, Foty, Kohn, Steinberg 2001**, PNAS 98(8):4323-4327 | 10.1073/pnas.071615398 (PMID **11274361**) | — | wetting as a **tug-of-war** between cell-cell and cell-substratum adhesion (originator framing). | 2-1 |
| **Beaune/Dufour, Yousafzai et al. 2022**, Phys Rev X 12:031027 | 10.1103/PhysRevX.12.031027 | S180 | elastocapillary crossover **Ec≈2.8 kPa**; basal pressure P₀≈109.6±6.2 Pa vs Pmin≈135.8±4.4 Pa (traction). | 3-0 ✓ |
| **Pérez-González et al. 2019**, *Active wetting of epithelial tissues*, Nat Phys | 10.1038/s41567-018-0279-5 (verify) | MDCK | active wetting = competition **traction vs contractile intercellular stress** — "no counterpart in classical wetting". | 2-1 |

## B. γ (tissue surface tension / cohesion) anchors — the W_cc=2γ term in S

| paper | DOI | system | γ (mN/m = mJ/m²) | verdict |
|---|---|---|---|---|
| **Nagle et al. 2022**, *Surface tension of model tissues during malignant transformation and EMT*, Front Cell Dev Biol | 10.3389/fcell.2022.926322 (verify; PMC9468677) | **breast epithelial** | **MCF10A 45±18; MCF10DCIS.com carcinoma 21±9; NME1-ablated 4–7** | 3-0 ✓ |
| **Butler & Foty 2011**, JoVE (PMC3169255) | — | generic TST | **0.33 (zebrafish ectoderm) → ~20–23** (chick limb-bud 20.1±0.5; rat prostate fib 22.8±1.1) | 3-0 ✓ |
| **Yousafzai 2022 PRX / Guevorkian 2010 PRL** | (above) / arXiv:1003.4372 | **S180** (the wetting proxy) | **γ≈5–7 mN/m** (size-dep); Guevorkian γ₀≈6 | 3-0 ✓ |
| **Manning, Foty, Steinberg, Schoetz 2010**, PNAS 107(23):12517 | 10.1073/pnas.1003743107 | model tissues | γ set by **adhesion/cortical-tension RATIO**, crossover γ_adh/β≈2 (adhesion- vs cortical-dominated). ⚠️symbol hazard: Manning's γ=adhesion param, σ=tissue surface tension; in S the "γ" = σ. | 3-0 ✓ |
| **Caicedo-Carvajal, Shinbrot, Foty 2010**, PLoS One | 10.1371/journal.pone.0011830 | CHO/FN | integrin(α5β1)-density → γ **NON-monotonic** 5.7(Low)→10.4(Mid)→5.5(High) — closest integrin→energy scaling (no absolute W_cs). | 3-0 ✓ |

## C. ⭐ MCF7-specific spreading — the closest direct data (ACTIVE, not passive)

| paper | DOI | system | content | verdict |
|---|---|---|---|---|
| **Aslemarz/Gupta et al. 2024**, *An EpCAM/Trop2 mechanostat differentially regulates collective behaviour of human carcinoma cells*, EMBO J 44(1):75-106 | 10.1038/s44318-024-00309-9 (PMID 39572744) | **MCF7** on collagen-I | **footprint 3–4× in 24 h** (active, Mitomycin-C excludes proliferation); **EpCAM-KD → complete wetting → flat coherent MONOLAYER by 48 h** (~2–3 cells thick @24h). **Three-tension** balance γ_m (free edge/cortex), γ_c (cell-contact), γ_x (matrix); contact angle = wetting readout. **Dominant lever = reduced cell-cell-contact contractility (cohesion)**, not purely W_cs. Explicitly ACTIVE wetting. | 3-0 ✓ |
| **Lemahieu et al. 2025**, *RAB5A Promotes Active Fluid Wetting…*, Adv Sci (PMC12442610) | (verify) | breast cancer spheroid | footprint **A~t^β: β=2/3 passive viscoelastic → β=1 active fluid**; RAB5A fluidizes → β→1. | 2-1 |

---

## How these reframe the model (→ DESIGN, REPORT)

1. **The magnitude (A/A₀~7–10) is the ACTIVE complete-wetting regime** (S>0): a motility-driven
   precursor MONOLAYER film spreading from the aggregate (Douezan/Beaune/EMBO J), over the PI ~80 h.
   NOT passive surface-energy wetting alone (that sets the equilibrium; motility realizes it).
2. **Governing criterion: S = W_cs − 2γ.** Partial cap (our current state) = S<0; complete wetting
   (large A/A₀) = S>0, reached by **raising W_cs (adhesion) OR lowering γ (cohesion)** — the EMBO J
   says the dominant physiological lever for MCF7 is **lowering cell-contact contractility (γ)**.
3. **Explains the §E + passive-diagnostic nulls:** §E active crawl was at adhesion_ratio=1 (S<0,
   partial-cap regime → rim deformation only); the passive adhesion sweep (no motility) can't drive
   the precursor film (kinetic trap; needs active motility per Beaune/Pérez-González). The decisive
   test couples **active motility + the S>0 regime + ductile cohesion**.
4. **Anchoring path (no-magic-number):** map model D_e→γ (anchor to breast-epi γ≈21–45 mN/m, Nagle)
   and D_sub→W_cs via the per-cell-energy↔area-energy-density conversion; the PI ligand axis
   (Bare/Pre/Lam4) modulates W_cs (the S180 FN-coverage 35%→51% transition is the template). ⚠️ the
   D2 surface-tension bridge earlier gave γ≈0.57 mN/m — ~40–80× BELOW Nagle's breast-epi 21–45 mN/m;
   reconcile (the D2 γ may be a different/under-estimate; flag for the σ-bridge owner).

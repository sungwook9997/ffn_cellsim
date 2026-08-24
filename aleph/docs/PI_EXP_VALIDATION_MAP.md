# PI experiment → platform validation map (MCF7 / pV4D4 / col-I)

> **STATUS: PREP scaffold — 2026-06-01. NOT a contract.** Sets up in-silico
> validation of the platform against the PI's own wet-lab study, as a **layered
> roadmap** (single-cell now → multi-cell collective later). Three deep-research
> dossiers (2026-06-01/02): `wf_dc812db7-9b4` (Layer-1, 22 claims),
> `wf_1748e64b-b1a` (Layer-2 + receptor), `wf_1494c786-e79` (laminin clutch).
> **Status**: Layer-1 anchored; Lam4 receptor = **α6β1** (ratified), clutch via
> affinity+proxy recipe (direct force kinetics absent in lit); Layer-2 active-wetting
> framework + bands (non-MCF7), **PI's A/A₀ model has no lit analog = novel**; iCVD
> pV4D4 → Im lab only. Companion to `ECM_PLATFORM_EXTENSIBILITY.md` +
> `CELL_MECHANICS_EXTEND_VS_REBUILD.md`. See memory `project-pi-exp-validation-target`.

## The experiment (overlay target — NOT a fitting target)

KSME 2026 poster + graduation research (Yoon, Jegal, Lee, Yeun, Im, Shin; KAIST;
grant RS-2024-00468873). **MCF7** (epithelial breast cancer, low-invasion)
**spheroid collective invasion / 2D spreading** on an **iCVD pV4D4** surface,
**collagen-I coated**. Three protein-presentation conditions:

| Cond | Presentation | integrin β1 IF | Spreading phenotype |
|---|---|---|---|
| **Bare** | untreated pV4D4 | weak / diffuse | high size threshold, large spheroids only |
| **Pre** | surface-adsorbed protein | peripheral | reduced threshold, constrained for small aggregates |
| **Lam4** | soluble laminin-111, 4 µg/mL in media | uniform | small-size penalty removed (c→0 or +), scale-independent |

Quantified: AI (U-Net) segmentation → normalized area **A/A₀ = a + b/R + c/R²**
(a baseline, **b traction/curvature**, **c small-size penalty**); R = initial
effective radius.

## HARD rule (governs this whole doc)

Literature-first. KU-anchor bands below come from **published literature**; the PI
poster/report values are **overlay-only at comparison time, never a fitting
target** (CLAUDE.md: "No fitting to PI experimental data"). If a platform output
needs the PI value to pass, that is a violation — surface to PI.

## The scale bridge (why this is layered)

The experiment is a **collective / multi-cell** phenomenon; the platform is
currently **single-cell**. The A/A₀ model decomposes the collective behavior into
terms the two layers reach separately:

- **Traction term b** ← single-cell traction / FA / clutch ⇒ **Layer 1 (now)**.
- **Cohesion / small-size penalty c, A/A₀ itself** ← cell–cell cadherin + collective ⇒ **Layer 2 (needs multi-cell extension; EXTEND doc defers cadherin/multi-cell)**.

---

## LAYER 1 — single-cell building blocks (platform reaches now)

Maps the experiment's **traction side** onto the existing single-cell runtime.
Validatable once the ECM/ligand layer lands on the single-cell platform
(after FA production; default-off modules already scaffolded).

### Platform mapping

| Experiment element | Platform element | Status |
|---|---|---|
| col-I-coated rigid dish | ECM preset **(b) "Col-I-coated rigid dish (PI exp)"** — rigid backing (`k_sub→∞` / pin), thin coating, collagen-I, 2D open | preset defined (`ECM_PLATFORM_EXTENSIBILITY.md:60`); `ecm/substrate.py` scaffolded |
| col-I vs laminin-111 (Bare/Pre vs Lam4) | **ligand-identity axis** on FA: col-I→α2β1/GFOGER vs **laminin-111→α6β1** (NOT α3β1 — see correction below) → integrin binding params (k_on, catch params, capture) | **TODO: extend `bridge/fa.py`** (ligand species) |
| integrin β1 IF (diffuse/peripheral/uniform) | FA / clutch density + spatial distribution (`bridge/integrin_bonds.py` Pereverzev catch-slip; `bridge/fa.py` FA placement) | clutch runtime exists; β1-distribution observable = **TODO** |
| traction (term b) | substrate→cortex tension / method-of-planes; bound-integrin count × clutch force | tension measurement exists (KU-3.5 grip-walk); per-ligand contrast = **TODO** |
| surface-bound (Pre) vs soluble (Lam4) | ligand areal density + on-rate / availability difference in FA binding | **TODO: mechanism choice** (see LIT §3) |

### KU-anchor band candidates (from deep-research dossier `wf_dc812db7-9b4`)

Confidence + vote from adversarial verification (need 2/3 to survive). All are
**candidate literature bands → acceptance oracles**, NOT fitting targets.

| Quantity | Band (lit) | Source | Conf | Maps to |
|---|---|---|---|---|
| β1 catch-bond reference (α5β1/FN) — F*, τ_peak | catch regime **10–30 pN**; peak lifetime **2–10 s** at 20–25 pN; <2 s below 10 pN | Kong 2009 *JCB* 185:1275, PMC2712956 | high (3-0) | `integrin_bonds.py` Pereverzev |
| α2β1 / collagen-I **slip**-bond (GFOGER I-domain) | k_off **0.44 s⁻¹**, x_β **0.7 nm** (locked-open hi-aff; GPP ctrl 11 s⁻¹/0.37 nm) | PMC3588017 (AFM SMFS) | high (3-0) | col-I clutch (Bare/Pre) |
| α2β1 / collagen-I single-cell slip + rupture | k_off **1.3±1.3 s⁻¹**, lifetime **0.8±0.7 s**, x_β **2.3±0.3 Å**; single-bond rupture **47±13 pN** @500 pN/s (38→90 pN over 180–8800 pN/s) | Taubenberger 2007 *MBoC* 18:1634, e06-09-0777 | high (3-0) | col-I clutch |
| Whole-cell col-I adhesion maturation | **189±12 pN** (5 s) → 1–2 nN (120 s) → 5 nN (180 s) → up to **20 nN** (>180 s); ~10× over 5→120 s | Taubenberger 2007 (CHO-A2) | high (3-0) | FA maturation timescale |
| Motor-clutch base set (Odde) | clutch k_on **0.3**/k_off **0.1 s⁻¹**, F_b **2 pN**, **50** clutches, κ_c **0.8 pN/nm**; myosin F_m **2 pN**, n_m **50**, v_u **120 nm/s**; optimal stiffness **2–300 kPa** | Bangasser 2013 *BiophysJ* 105:581, PMC3736748 | high (3-0) | clutch/motor **defaults** (model, not measurement) |
| Surface-bound (Pre) vs soluble (Lam4) mechanism | **affinity vs avidity**: clustering raises occupancy only under cooperative binding; TM-disruption raises monomeric affinity; cilengitide agonist@1 nM / antagonist@1 µM | Irvine 2002; Partridge/Luo 2005 PNAS; Steiger 2024 *JMC* | med (mixed) | Pre vs Lam4 — **by analogy (αIIbβ3/αvβ3/RGD), laminin extrapolation** |
| MCF7 vs MDA cytoplasm viscosity | MCF-7 **65.9±11.4** vs MDA **12.0±5.7 Pa·s** (~5×; band **12–66**) | Hu 2024 *Nanoscale Adv*, PMC10929591 (MRS) | high (3-0) | H.10 KU-3.B3 (refines existing ≈56/10.7) |
| Intracellular elastic shear modulus G | MCF-10A **79.3** / MCF-7 **32.9** / MDA **38.6 Pa** (band **33–79**; cancer ~2× softer) | Hu 2024 PMC10929591 | high (3-0) | H.10 (distinct from AFM E) |
| Whole-cell AFM modulus E₀ (1 s) + fluidity | E₀ MCF-10A **1.14** / MCF-7 **0.26** / MDA **0.46 kPa**; fluidity exp **0.186 / 0.234 / 0.147** | Yubero 2020 *Commun Biol*, s42003-020-01330-4 | high (3-0) | whole-cell modulus |
| Cortical tension | **~10⁻² N/m** (10⁻²–10⁻¹), short-pulse | PMC10625872 (electrodef.) | high (3-0) | **H.3 cortex** direct |

**Contested / refuted (do NOT use as anchors):**
- ⚠️ **MCF7-vs-MDA stiffness *ordering* is method-dependent & disputed**: AFM power-law → MCF-7 *softer* (0.26 vs 0.46 kPa); electrodeformation → MCF-7 *~10× stiffer* (7.1 vs 0.7 kPa). Use as single-method anchors, not consensus. (Cytoplasm-viscosity + contractility contrasts are more robust.)
- ❌ REFUTED (killed in verification): MCF7 traction 15–25 nN micropillar (1-2); soluble-FN 5–6× adhesion drop (1-2); α2β1 I-domain "pure slip, no catch" framing (0-3, but the numeric k_off/x_β stand).
- ⚠️ AFM rupture forces measured at loading rates ≫ physiological (~1 pN/s); cytoplasm η is probe-length dependent (η∝L²) — state L.

### Laminin (Lam4) — follow-up dossier `wf_1748e64b-b1a` results

**RECEPTOR CORRECTION (high, 3-0)**: the laminin-111 receptor for Lam4 is
**integrin α6β1, NOT α3β1**. α3β1 shows *no significant binding to LN-111* (it binds
LN-332 / LN-511/521); α6β1 binds LN-111 but at *low* affinity (rank LN-10/11 >
LN-5 > **LN-1** > LN-2/4~LN-8). Full breast laminin-receptor set: α6β4, α7β1, α6β1,
α3β1 — but **MCF7 lacks α6β4** (so no α6β4-keratin mechanoprotection). → in
`bridge/fa.py` the Lam4 ligand species should bind **α6β1** (low-affinity LN-111).
*Sources*: Nishiuchi 2006 PMID 14607975; Mosqueira 2023 *Nat Mater* PMC10627833.

**Traction direction (high, 3-0)**: breast epithelial cells exert **lower traction
on laminin-111 than on collagen-I/FN**; YAP N/C < 2 on laminin (robust 0.5–30 kPa).
⚠️ **Validation flag (single-cell vs collective)**: this *opposes* the poster's
collective finding (Lam4 = *elevated* traction + preserved cohesion). Not a
contradiction — Lam4 is **soluble laminin supplement on top of the col-I coating**,
and the elevated traction is a **collective** (Layer-2) phenotype, not a
laminin-only single-cell substrate. The single-cell platform should reproduce
*lower laminin-only traction*; the Lam4 collective effect is a Layer-2 question.

### Lam4 (α6β1) clutch parameterization — dossier `wf_1494c786-e79`

**FINDING (high, 3-0): direct α6β1–laminin single-molecule force kinetics are
GENUINELY ABSENT from the literature** (exhaustive sweep: no rupture force,
lifetime, x_β, loading-rate, or catch/slip classification for any α6β1–laminin
bond). This is itself a real result — the Lam4 clutch **cannot** be set from direct
laminin force-spectroscopy; use the layered fallback recipe below.

**Recommended Lam4 clutch recipe (literature-anchored, explicitly proxy-flagged):**

| Component | Value | Basis | Conf |
|---|---|---|---|
| Off-rate scale (affinity fallback) | **K_D ≈ 1–20 nM** (α6β1–LN-111 toward weak/upper end; α6β1–LN-511/521 = 0.73–0.74 nM) | Nishiuchi 2006 *Matrix Biol* 25:189; Taniguchi 2009 PMC2658076; Stipp 2010 PMC2811424 | high (equilibrium only, NOT k_off) |
| Affinity rank (α6β1) | **LN-511/521 > LN-332 > LN-111** (LN-111 = moderate/3rd) | Nishiuchi 2006 | high (3-0) |
| Force-kinetic **shape proxy** | k_off **1.4–2.3 s⁻¹**, f_b **12–15 pN** (14.1±1.3), x_β **~0.28 nm**, **single-barrier slip** | α7β1–**invasin** (NOT laminin) PMC3882471 | med (2-1/3-0) — **proxy only** |
| Relative scaling | Lam4 clutch **WEAKER / more-slipping than the FN (α5β1) clutch** — lower traction (P=0.016–0.028) + weaker adhesion (P<0.04) | CHO.B2 transfectants PMC3391238 / PLoS ONE e40202 | high (3-0) — non-MCF7 |
| FN benchmark to scale against | α5β1–FN **catch** bond, 10–30 pN peak ~10 s | Kong 2009 PMC2712956 | high |

**Recipe**: set the α6β1 clutch as a **Bell-Evans slip bond** using the α7β1–invasin
*shape* (x_β≈0.28 nm, f_b≈12–15 pN, k_off⁰≈1.4–2.3 s⁻¹), **scaled weaker than the
FN clutch** per the relative-traction anchor. Tag every value as proxy/affinity-
derived, NOT direct measurement.

**Caveats / cheap follow-ups**: (1) catch-vs-slip for α6β1–laminin is *unclassified*
— recipe assumes slip (default, like α7β1-invasin). (2) No measured k_on → generic
integrin k_on assumption needed to turn K_D into an absolute off-rate. (3) Soluble
LN-111 (Lam4) vs immobilized may alter avidity/loading geometry. (4) **Specific
α6β1–LN-111 K_D** sits in the paywalled Nishiuchi 2006 table → recover via
gbook/KAIST full-text ([[reference-gbook-kaist-fulltext]]) — cheap. (5) REFUTED
(1-2): "α6β1-laminin retrograde flow faster" — dropped.

---

## LAYER 2 — collective / spheroid (needs multi-cell extension)

Maps the experiment's **cohesion side + A/A₀ itself**. Requires the cadherin
multi-cell layer — the furthest-out roadmap item (`CELL_MECHANICS_EXTEND_VS_REBUILD.md`
defers cell–cell cadherin / multi-cell). Documented now so Phase 2 starts ready.

| Experiment element | Platform element (future) | Status |
|---|---|---|
| A/A₀ = a + b/R + c/R² spreading | emergent from multi-cell spreading sim | **Phase 2** |
| cell–cell cohesion (term c) | cadherin trans-bonds between cortex particles (Bell-Evans; KU-4.2 catch-bond) | **NEW multi-cell layer** |
| traction–cohesion balance | per-cell FA traction (Layer 1) × cadherin cohesion | composition |
| Lam4 removes small-size penalty | ligand-modulated traction + preserved cohesion | Phase 2 hypothesis |

### KU-anchor band candidates (Layer 2 — partially anchored, dossier `wf_1748e64b-b1a`)

**Framework (high, 3-0)**: collective spreading = **active wetting** (Young–Dupré):
spreading vs cohesion set by **cell–substrate traction vs cell–cell contractility**;
size dependence enters via a screening/localization length **Lp = −ζ/ζ_i** (traction
localized to the colony edge). Cohesion controls edge-vs-bulk traction partition.

⚠️ **No clean published A/A₀ = a + b/R + c/R² analog exists.** The closest viscoelastic
wetting scaling (A ∝ R₀²·t^{2/3}, capillary V*=γ/η) and the active-wetting Laplace
parameter set were **both REFUTED (0-3)** in verification. **Implication: the PI's
size model is genuinely novel** (consistent with the self-assessment that its
novelty / physical-meaning needs strengthening) — the b/R, c/R² terms are *analogous
to* but not *derived from* a standard wetting law. Frame the platform's job as
*reproducing the traction–cohesion partition that the a/b/c terms encode*, not
matching a literature A/A₀ curve.

| Quantity | Band (lit) | Conf | Source |
|---|---|---|---|
| E-cadherin junction tension (homeostatic) | **1–5 pN** | high | Sim 2015 *MBoC* e14-12-1618 (MDCK, FRET) |
| Cell–cell vs cell–ECM force partition | cell-ECM 573→848 nN, cell-cell 372→935 nN (MDCK monolayer) | 2-1 | Sim 2015 PMC4571300 |
| Cluster apparent surface tension γ | **0.15–0.75 mN/m** (cluster) vs 4–45 mN/m (3D aggregate) | high | Pallarès 2023 PMC7617391 |
| Traction screening length Lp; size crossover | Lp **~11 µm**; crossover **R~50 µm** (Lp/L≈0.2); ζ·dμ ~4 kPa | high | Mertz 2013 PNAS 1217279110; Pérez-González 2019 *Nat Phys* (MCF10A) |
| Unjamming phases (size/cell-number) | solid (<25 µm / <10 cells) / liquid (>25 µm) / gas (>10 cells); invasion onset 3/5/11.5 h; low E-cad → gas/EMT analog | high | Wagena 2024 PMC11665421 (MV3); Douezan 2011 (S180) |

⚠️ **All Layer-2 anchors are NON-MCF7** (MCF10A, MDCK, S180, MV3 melanoma) — use as
order-of-magnitude bands, not MCF7-specific. Apparent surface tension ~0.8 mN/m
(keratinocyte) was REFUTED (1-2) — dropped.

(Existing platform anchor: cadherin KU-4.2 catch-bond — Δx*=4 nm, k_off⁰=0.5 s⁻¹,
⟨F⟩=30 pN, N_cad=100 — survives from Phase 1 as the cohesion seed; consistent with
the 1–5 pN homeostatic junction-tension band above.)

---

## Surface context — iCVD pV4D4 protein adsorption (⚠️ STILL UNFILLED)

Both sweeps (`wf_dc812db7-9b4`, `wf_1748e64b-b1a`) found **no surviving primary
quantitative claim** on pV4D4 protein-adsorption density/conformation or
adsorbed-vs-soluble laminin presentation. This is niche/possibly-unpublished
surface-chemistry. **Do NOT burn a third generic web sweep** — route instead to:
- the **collaborator's lab** (Im Sung Gap, KAIST ChemBioE — made the pV4D4 surface):
  their iCVD pV4D4 protein-adsorption / contact-angle data is the real source.
- KAIST institutional full-text via gbook (see memory `reference-gbook-kaist-fulltext`).
Interim handle: the Layer-1 affinity-vs-avidity mechanism row, applied by analogy.

## Remaining gaps after both sweeps

1. ~~Laminin receptor identity~~ → **DONE: α6β1** (low-affinity LN-111). Numeric
   **α6β1–LN-111 clutch kinetics** (k_off/x_β/F*) still missing — single-molecule LN
   force-spectroscopy targeted search, or proceed with α6β1-identity qualitative anchor.
2. ~~Layer-2 framework~~ → **DONE: active wetting + bands** (non-MCF7). Open: an
   MCF7-specific collective dataset; confirm the novel A/A₀ form has no closer analog.
3. **iCVD pV4D4 adsorption** → route to Im lab (above), not a web sweep.

## Sequencing

1. (active) FA production lands → single-cell ECM/ligand layer wired + validated.
2. **Layer 1**: add ligand-identity to `bridge/fa.py`; run preset (b) single MCF7,
   col-I vs laminin-111; compare traction + β1 distribution to KU bands; PI poster
   as overlay.
3. **Layer 2**: multi-cell cadherin extension → A/A₀ spreading; PI A/A₀ as overlay.

## Open items for PI

- [x] Layer-1 bands filled (dossier `wf_dc812db7-9b4`).
- [x] Layer-2 framework + laminin receptor filled (follow-up `wf_1748e64b-b1a`).
- [x] **Receptor corrected: Lam4 = α6β1, not α3β1.**
- [x] **RATIFIED (PI 2026-06-02)**: Lam4 ligand = **α6β1** (low-affinity LN-111) for `bridge/fa.py`.
- [x] α6β1–LN-111 **clutch-force kinetics** → **RESOLVED (dossier `wf_1494c786-e79`)**: direct force kinetics GENUINELY ABSENT in lit. Use the fallback recipe (K_D 1–20 nM affinity + α7β1-invasin slip-bond shape proxy, scaled weaker than FN). Optional cheap follow-up: recover α6β1–LN-111 K_D from paywalled Nishiuchi 2006 via gbook/KAIST.
- [x] **iCVD pV4D4 adsorption** → **PI: accepts no data exists** (web). Im Sung Gap lab is the only route; no further action this session.
- [ ] Validation flag to ratify: single-cell laminin traction is *lower* (lit) while Lam4 collective is *higher* (poster) — confirm this is the single-cell↔Layer-2 split, not a model error.
- [x] **Bands ported into Notion KU v2 layer** as KnowledgeClaim rows **KB-PIV-1…10** (2026-06-02; Subtopic "PI-exp validation", Status `verified`, provenance in Citations-text). Evidence→SourceEvidence relation backfill deferred (same standing step as rest of Contract Graph).
- [ ] Ratify ligand-identity FA extension scope (col-I α2β1 + laminin-111 α6β1 first).
- [ ] Confirm Layer-2 multi-cell is the right vehicle for A/A₀ (the model has no clean lit analog — it is novel).
- [ ] Method-tag the contested MCF7-vs-MDA stiffness ordering (dual band, not single value).

# H.7 — Cortical NMII minifilament areal-density datum (resolving the "Salbreux 2012, 3.0/µm²" misattribution)

**Date:** 2026-06-07
**Author:** literature subagent (ffn_cellsim)
**Task:** Pin the cortical non-muscle myosin-II (NMII) minifilament **areal density** (minifilaments
per µm² of cortex) for MCF7, and resolve the runtime's `areal_density_per_um2 = 3.0` value attributed
to "Salbreux 2012" — flagged a CONFIRMED MISATTRIBUTION by the 13-agent generation-gap audit
(`_historical/H7_GENERATION_GAP_AUDIT_2026-06-07.md`). Interim substitute on record: Nie et al. 2015
(HeLa medial cortex ~0.6/µm², LOW confidence).
**Scope:** literature only — **no code edits, no config edits, no git, no Notion, nothing
auto-registered.** Verification via PubMed (metadata + full-text + ID-conversion), the local
109-PDF full-text corpus, and web cross-check. This doc is a candidate record for a future
PI-gated SourceEvidence registration, not a registration.

> ⚠️ **Citation-integrity discipline.** This project has **3 confirmed hallucinated SourceEvidence
> rows** on record — `Yao2011_NatCommun`, `YapKovacs_JCS`, `NanoConvergence2021_Glioma`. The
> dominant failure-mode is *a real author name (or real DOI) attached to a number/identifier that
> source never published.* Every value below is tied to a real, retrievable paper with DOI + PMID +
> PMC and an exact evidence location. **The verification pass caught several wrong-PMID / fabricated-
> author-list defects in the raw research findings; those are corrected inline here and the affected
> rows are DOWNGRADED — see §11.** No value is presented as fact unless DOI+identifier+the paper
> actually stating it were all independently tied together.

---

## 1. Header recap

The runtime carries `density_per_um2 = float(cfg.get('areal_density_per_um2', 3.0))  # Salbreux 2012`
(`cortex/myosin.py:309`; `configs/phase1_h3.yaml:326`; `configs/mcf7_baseline.yaml`). The active-γ
ceiling scales **linearly** with this density (`cortex/active_gel_seam.py:98`;
`g_ceiling = N_total_heads · f_stall · ℓ_dipole / area`, `N_total_heads ∝ areal_density`). This note
(a) finds the best real cortical minifilament areal density, (b) pins the foci→minifilament conversion
factor that governs whether a reported "foci/µm²" is already a minifilament density, (c) re-verifies
the Salbreux attribution, and (d) reports — honestly, not gate-chased — what the corrected density does
to the ceiling.

---

## 2. TL;DR (≤8 lines)

1. **Best density anchor: ~0.6 minifilaments/µm²** (round) — HeLa interphase, medial/basal adhered
   cortex; Nie et al. 2015 *Cytoskeleton*, **40 ± 20 MRLC-GFP foci in an 8×8 µm region = 0.625/µm²**
   (range 0.31–0.94), DOI 10.1002/cm.21207 / PMID 25641802 / PMC4361371. **LOW–MEDIUM confidence,
   non-MCF7, non-breast.**
2. **No MCF7-specific and no breast-epithelial cortical NMII areal density exists** in PubMed, the
   local corpus, or the web (medium confidence — a negative cannot be search-proven). HeLa is the
   closest usable cultured-metazoan cortex datum.
3. **Foci vs minifilament:** region-dependent, **central ~2, range 1–5** (up to ~10 in dense registered
   stacks); **NOT pinnable to a single literature number.** CRITICAL: the Nie ~0.6/µm² is **already a
   minifilament density** (authors intensity-calibrated foci→minifilaments) — it must **NOT** be
   re-multiplied by this factor (would double-count).
4. **Salbreux attribution: NOT salvageable for DENSITY** — confirmed misattributed. No Salbreux paper
   states 3.0/µm² or any per-µm² minifilament count (§7). (The separate Salbreux *tension* attribution
   is out of scope and untouched.)
5. **Ceiling impact (honest, not gate-chasing):** ceiling ∝ density. M1 "most-generous" ceiling = 0.10
   mN/m at 3.0/µm² → at the **0.6/µm² anchor, ceiling ≈ 0.02 mN/m** (~5× lower). Lowering the density to
   the evidenced value **WIDENS** the generation gap vs the MCF7 datum γ ≈ 0.27 mN/m (Hosseini 2020).
6. **This is a magic-number / citation-integrity fix → drop the `# Salbreux 2012` comment and surface
   3.0/µm² to PI.** Density must NOT be picked to make the ceiling reach the datum.

---

## 3. Verified value table

Only VERIFIED values appear. Identifier defects flagged by the verification pass are corrected here;
rows whose identifier could not be cleanly tied are marked and down-weighted (see §11).

| Value | Unit | Cell type | State | Foci or minifilament | DOI / PMID / PMC | Confidence |
|---|---|---|---|---|---|---|
| **0.625** (range **0.31–0.94**) | **minifilaments/µm²** (40 ± 20 foci per 8×8 µm = 64 µm²) | **HeLa** (MRLC-GFP) | interphase, adhered; **medial/basal** cortex (Type II) | **minifilament** (authors intensity-calibrate 1 focus ≈ 1 minifilament) | 10.1002/cm.21207 / 25641802 / PMC4361371 | **medium** — ⭐ANCHOR |
| 0.66 | minifilaments/µm² (42 in 8×8 µm model region) | HeLa (computational, matched to expt) | interphase adhered, Type II baseline | minifilament | 10.1002/cm.21207 / 25641802 / PMC4361371 | low (model input ≈ same datum as 0.625 — do not double-count) |
| 0.79 | minifilaments/µm² (252 in 20×16 µm = 320 µm²) | HeLa (computational, dense variant) | interphase adhered, FA-patterned | minifilament | 10.1002/cm.21207 / 25641802 / PMC4361371 | low (model "~20% larger" variant; brackets upward, not a 2nd measurement) |
| **0.2** | minifilaments(motor oligomers)/µm² | **in-silico / in-vitro** (AFINES; 50×50 µm box) | simulation parameter | minifilament | 10.1016/j.bpj.2017.06.003 / **28746855** (⚠ raw finding's 28746850 WRONG) / PMC5529201 | low — **illustrative, NOT a cortex datum; do not promote** |
| **No value** | n/a | **MCF7** | any | n/a | — | **medium** — confirmed absent in PubMed + local corpus + web |
| No value | n/a | breast-epithelial | any | n/a | — | medium — confirmed absent (closest "epithelial" is zebrafish periderm spacing, not /µm²) |
| 295 ± 11 | **nm** (single-minifilament *projected length*, NOT a density) | HeLa (TDS / S-HeLa) | rounded interphase (SIM/dSTORM) | minifilament | 10.1038/s41467-021-26611-2 / 34764258 / PMC8586027 | medium — geometry support only; paper reports NO areal density |
| 320 / 28 | **nm** length / **molecules** per single bipolar minifilament (structure, NOT a density) | purified human platelet myosin II | in-vitro EM | minifilament | 10.1083/jcb.67.1.72 / 240861 / PMC2109578 | high — structural anchor (single minifilament ≈ optical diffraction limit) |
| 200 / 52 (typ. ~250; range 20–500) | **nm** (ACTIN membrane-skeleton mesh side-length, NOT myosin) | NRK fibroblast / FRSK keratinocyte | interphase adherent (platinum-replica EM tomography) | n/a (actin) | 10.1083/jcb.200606007 / 16954349 / PMC2064339 | high — **EM context only; NOT a myosin density** |
| 200–300 | **nm** (myosin filament length; "regularly spaced clusters", head-to-head chains) | sea-urchin embryo contractile ring | cytokinesis (platinum-replica TEM) | minifilament | 10.1091/mbc.E16-06-0466 / 28057763 / PMC5328620 | medium — qualitative spacing; NO /µm² density |

**Non-MCF7 / non-breast proxies are marked clearly above.** The only explicit per-µm² minifilament
*density* values in the whole search are the HeLa Nie rows (~0.6–0.8) and the AFINES in-silico 0.2.

---

## 4. Recommended anchor for `areal_density_per_um2`

**Primary anchor (use this): ~0.6 minifilaments/µm² (single value), range 0.3–0.9 minifilaments/µm²**
— HeLa, interphase, adhered medial/basal cortex; Nie et al. 2015 *Cytoskeleton* 72(1):29-46
(DOI 10.1002/cm.21207, PMID 25641802, PMC4361371). **Confidence: LOW–MEDIUM.**

- **Exact datum:** "we measured the intensity of MRLC-GFP foci … to estimate **= 40 ± 20** (Mean ±
  StDev) **in an 8×8 µm region** for control Type II cells" (PMC4361371 full text). 40/64 = 0.625/µm²;
  the ±20 gives the 0.31–0.94 range. The model baseline "**= 42 randomly placed minifilaments**" (0.66)
  and the dense FA-patterned variant "**252 minifilaments**" in a 20×16 µm box (0.79) bracket it.
- **Cell-type/state it represents:** cervical-carcinoma HeLa, **interphase, spread on glass, medial/
  basal cortex** — explicitly NOT the apical/equatorial cortex that sets rounded-cell γ, and NOT MCF7.
- **Why not MCF7-specific:** **no MCF7 (and no breast-epithelial) cortical NMII areal density exists**
  anywhere searched (§11). HeLa is the closest real cultured-metazoan cortex datum; the HeLa→MCF7
  transfer is a necessary but **unvalidated extrapolation** — flag to PI.
- **Confidence drivers:** the count is an **intensity-calibrated** estimate (authors assume ~20
  motors/minifilament and **1 focus ≈ 1 minifilament**), carrying ±50% stated error plus that modeling
  assumption; it is NOT a super-resolution 1:1 filament count. So the value sits *between* "foci/µm²"
  and "minifilament/µm²" — see §5.

**Do NOT use** the AFINES 0.2/µm² (illustrative in-silico parameter) or the actin-mesh EM numbers
(~16–25 mesh-cells/µm²; those are ACTIN, not myosin) as the density anchor.

---

## 5. Foci → minifilament conversion

**Central factor ~2; range 1–5; up to ~10 in the densest registered stacks. NOT pinnable to a single
literature number — this is the dominant unquantified uncertainty. Confidence: MEDIUM that the true
cortical value lies in 1–5.** The factor is **region-dependent, not a constant.**

| Regime | Factor (minifilaments per focus) | Basis (verified) |
|---|---|---|
| Sparse / resolved | **1** | Beach 2014 TIRF-SIM resolves a single ~300 nm bipolar filament as two ~300-nm-spaced puncta (PMID 24814144, PMC4108432); Truong Quang 2021 SIM/dSTORM single-minifilament length 295 ± 11 nm (PMID 34764258, PMC8586027); Niederman & Pollard 1975 one minifilament ≈ 320 nm ≈ optical diffraction limit (PMID 240861). |
| Dense / registered stacks | **several → many (~5–20 *estimated*)** | Hu 2017: filaments "registered alignment into stacks, spanning up to several micrometres" (PMID 28114270) — **the 5–20 number is a geometric inference** (several µm ÷ ~0.3 µm spacing), NOT a Hu-stated value; Beach 2017 cluster partitioning (PMID 28114272); Beach 2014 dense regions "too numerous … to identify individual filaments." |
| Leading-edge mix | ~30% unresolvably overlapping even when sparse | Fenix 2016 (PMID 26960797): 2/3/4-motor-group filaments are **ONE expanding filament each** (not stacks); true concatenation only ~10%; ~30% overlapping. (So motor-group count is NOT a stack count.) |

**Mapping a reported foci/µm² to minifilaments/µm²:**
`minifilaments/µm² = (foci/µm²) × (minifilaments per focus)`. Example: a *raw-foci* density of 0.6
foci/µm² × factor 2 ≈ 1.2 minifilaments/µm² (×5 → 3.0; ×1 → 0.6).

> ⛔ **Integrity guard (CONFIRMED by verification).** The Nie ~0.6/µm² anchor is **already a
> minifilament density** — the authors converted MRLC-GFP intensity to a minifilament count assuming
> ~20 motors/minifilament. It must **NOT** be multiplied by the foci→minifilament factor again; doing
> so double-counts. The factor applies ONLY to sources that report **raw** puncta/foci counts. *Always
> tag every candidate density as foci-vs-minifilament before applying any factor.* (Caveat: the
> authors' assumed ~1 factor is itself uncertain — confocal foci could hide z-stacks — so the TRUE
> cortical minifilament density could be somewhat HIGHER than 0.6, but you cannot reach it by
> multiplying Nie's already-converted number.)

---

## 6. EM cross-check

**Verdict: electron microscopy provides NO independent cortical myosin minifilament areal density in
/µm² for any cell — and none for MCF7. It can neither confirm nor refute the 0.6–3/µm² range
directly.** What EM *does* give:

- **Actin cortical/membrane-skeleton mesh size** (NOT myosin): Morone 2006 platinum-replica electron
  tomography — median side-length **200 nm (NRK fibroblast) / 52 nm (FRSK keratinocyte), range ~20–500
  nm** (PMID 16954349, PMC2064339, verbatim full-text confirmed); "typical ~250 nm" re-cited by
  Gowrishankar 2012 Cell (DOI 10.1016/j.cell.2012.05.008, **corrected PMID 22682254** — raw finding's
  22863013 was a wrong NMDA-receptor paper). Implies ~16–25 **actin** mesh-cells/µm² — this is **NOT a
  myosin density and must not be substituted for one.**
- **Myosin filament length + qualitative organization** (no /µm²): Henson 2017 sea-urchin contractile
  ring, platinum-replica TEM, length 200–300 nm, "regularly spaced clusters" and "laterally associated,
  head-to-head filament chains" (PMID 28057763, PMC5328620); Shutova/Svitkina 2014 platinum-replica EM
  + immunogold — individual activated NMIIA/NMIIB copolymerize into nascent bipolar filaments,
  concentrated in a layer **beneath** the actin cortex (PMID 25131674, PMC4160463).

**Implication for the conversion (§5):** EM *qualitatively favors a foci→minifilament factor > 1* —
minifilaments form head-to-head chains and lateral/sub-cortical stacks, so a fluorescence focus is
generally NOT one bipolar filament — **but EM quantifies no per-focus integer.** No retrievable
platinum-replica EM gives a cortical minifilaments-per-stack count; any specific stack integer is a
geometric estimate, not an EM measurement (down-weighted accordingly).

**Geometric sanity (context, not evidence):** 0.6/µm² ⇔ ~1290 nm spacing (sparse, ~one per 5×5 actin
mesh cells); 3.0/µm² ⇔ ~577 nm spacing (~one per 2×2 mesh cells). Neither is excluded by the ~250 nm
actin mesh / ~300 nm filament length, but neither is anchored either — geometry cannot rescue a missing
measurement.

---

## 7. Salbreux re-verification

**Verdict: "Salbreux 2012, 3.0/µm²" is CONFIRMED MISATTRIBUTED for the DENSITY claim and NOT
salvageable. No Salbreux paper states a cortical myosin minifilament areal density of 3/µm² or any
per-µm² minifilament number. Confidence: HIGH.** Independently re-verified against every Salbreux
full-text in the local corpus plus PubMed metadata for the review candidates.

The two flavors of Salbreux cortex work are each **structurally incapable** of producing a 3/µm² figure:

- **Discrete model — Chugh, Clark, … Salbreux, Paluch 2017** *Nat Cell Biol* 19(6):689-697 (DOI
  10.1038/ncb3525, PMID 28530659, PMC5536221; the only Salbreux *primary* cortex paper, local
  `ncb3525.pdf` ≡ `emss-72183.pdf` twin): myosin enters ONLY as a **2D line-concentration N_m/W** over a
  fixed 2.5 µm periodic box — *"the two-dimensional concentration of myosins, N_m/W. For our choice of
  parameters, T₀ ≈ 230 pN/µm"* (line 2793) — **never per-area.** Actin density ρ_f is a line density at
  ~75% SEM coverage. Exhaustive grep for "3/µm² myosin" across the file returns only unrelated "3"
  tokens (3 µM CHIR99021 inhibitor; 3-min spin). (The integer N_m lives in Supplementary Table 5, which
  doesn't extract as text — immaterial: the framing is per-LENGTH, never per-area.)
- **Active-gel / active-surface lineage — e.g. Bächer, Khoromskaia, Salbreux, Gekle 2021** *Front Phys*
  9:753230 (DOI 10.3389/fphy.2021.753230; local `fphy-09-753230.pdf`): **all** myosin is coarse-grained
  into a continuum **active surface stress tensor ζ_αβ (~10⁻⁴ N/m)** — by construction there is no
  discrete per-µm² minifilament density.

**The most likely "Salbreux 2012"** is Salbreux, Charras, Paluch 2012 *Trends Cell Biol* 22(10):536-545,
**actual title "Actin cortex mechanics and cellular morphogenesis"** (DOI 10.1016/j.tcb.2012.07.001,
PMID 22871642), a **REVIEW** with no quantitative myosin areal density (abstract/MeSH = Review;
paywalled body not exhaustively word-read, but all accessible evidence + the continuum framing + the
project's prior 13-agent audit converge on ABSENCE). The 2014 companion (Clark, Wartlick, Salbreux,
Paluch, *Curr Biol* 24(10):R484-94, "Stresses at the cell surface during animal cell morphogenesis,"
DOI 10.1016/j.cub.2014.03.059, PMID 24845681) is likewise a review with no density.

**Provenance of the bad number:** `cortex/myosin.py` reads
`density_per_um2 = float(cfg.get('areal_density_per_um2', 3.0))  # Salbreux 2012`, propagated by
`KU_3_5_GATE_STRUCTURE_2026-06-04.md:106` as "Salbreux 3/µm² × 706.9 µm² → ≈2120 motors." The
`# Salbreux 2012` comment is the misattribution vector — the **exact failure-mode** of the project's 3
confirmed prior hallucinations (real author name + a number that author never published).

> **Keep the TENSION attribution separate.** The separate "Salbreux 2012" *cortical-tension* band
> (~0.1–1 mN/m) is **NOT** part of this misattribution claim and was **not assessed here** — only the
> DENSITY attribution is refuted. The 2012/2014 reviews DO discuss cortical tension qualitatively; the
> *density number* is the misplaced item. (For the MCF7 tension anchor itself, see the sibling note
> `_historical/H7_MCF7_CORTICAL_TENSION_DATUM_2026-06-07.md`.)

**Recommendation:** drop the `# Salbreux 2012` comment; the 3.0/µm² value is un-anchored (~4.8× the
best HeLa datum, 0.625; ~3.8× the dense model variant, 0.79) and should be surfaced to PI as a
magic-number / citation-integrity fix.

---

## 8. Ceiling impact — diagnostic, not a magic-number to tune

The active-γ ceiling scales **linearly** with areal density (`active_gel_seam.py:98`):
`g_ceiling ≈ 0.10 × (D / 3.0) mN/m` (M1 "most-generous" params held; M1 = 0.10 mN/m at D = 3.0/µm²).

| Density D (minifilaments/µm²) | Basis | Ceiling (mN/m) | vs MCF7 datum γ ≈ 0.27 |
|---|---|---|---|
| 3.0 | runtime (un-anchored, misattributed) | 0.10 | ~2.7× under |
| **0.625** (anchor) | Nie 2015 HeLa | **≈ 0.021** | **~13× under** |
| 0.79 (dense model variant) | Nie 2015 | ≈ 0.026 | ~10× under |
| 0.31–0.94 (anchor range) | Nie 2015 ±SD | ≈ 0.010–0.031 | ~9–27× under |
| 1.9 (0.625 × max foci-factor ~3) | upper-edge inference | ≈ 0.063 | ~4× under |

**Reading (honest, not gate-chasing):** adopting the evidenced ~0.6/µm² over the unsupported 3.0/µm²
**LOWERS the ceiling ~5× (0.10 → ~0.02 mN/m) and WIDENS the generation gap** vs the real MCF7 datum
γ ≈ 0.27 mN/m (Hosseini 2020, *Adv Sci*, AFM, interphase suspended; DOI 10.1002/advs.202001276, PMID
33042748 — **the γ-protocol claim should be re-confirmed from full text before relying on it; the paper
is real and from the Fischer-Friedrich group but its headline subject is EMT/mitotic rounding**, see
§11). Even a *maximal* foci→minifilament multiplier (~3) reaches only ~1.9/µm² → ceiling ~0.063 mN/m,
still ~4× under the datum. **Correct anchoring makes the gap WORSE — which is precisely why density
must NOT be chosen to make the ceiling reach the datum.** The conversion factor does not rescue the
gap; it sharpens the conclusion that the misattributed 3.0/µm² is unsupported and that a
literature-anchored cortical minifilament density is sub-1/µm².

> **This is a MAGIC-NUMBER CONTRACT change → PI sign-off required.** Changing
> `areal_density_per_um2` (or dropping the `# Salbreux 2012` comment) alters a sanity-gate input. Per
> the no-magic-number / no-gate-loosening rules, surface to PI rather than editing inline. This note is
> the surfacing artifact, not an edit.

---

## 9. What an experiment would measure (no MCF7 density exists)

To fill the genuine MCF7 gap (for PI / future overlay):

- **Super-resolution count of NMII minifilaments per µm² of the MCF7 cortex** — SIM + two-color dSTORM
  on MRLC (or NMIIA tail + head, the Beach 2014 / Truong Quang 2021 two-color scheme) to resolve
  individual ~300 nm bipolar filaments, then count per unit cortical area. Report **both** a raw
  foci/µm² **and** a stack-resolved minifilaments/µm² so the foci→minifilament factor is *measured for
  MCF7*, not assumed.
- **At the physiological adherent operating point** (resting turgor, real ~65 Pa·s cytoplasm viscosity,
  engaged FA), and ideally at the **apical/equatorial cortex** that actually sets the tension the gate
  cares about — not just the basal/medial cortex that confocal foci-counting reaches most easily.
- **Cross-calibrate to tension on the same cell** (paired AFM-confinement or cortex-tether γ) so the
  density→tension link the model encodes (`g_ceiling ∝ density`) is validated against a measured
  density–tension pair, closing the loop the current misattributed number papered over.

---

## 10. Search provenance (for reproducibility / audit)

- **Local full-text corpus (highest fidelity, used first):** `grep -inE` over `/tmp/refs_txt/*.txt`
  (109 PDFs). Confirmed in-corpus: AFINES `main.pdf` (0.2 motors/µm²), Chugh-Salbreux `ncb3525.pdf` ≡
  `emss-72183.pdf` (N_m/W, no /µm²), Bächer-Salbreux `fphy-09-753230.pdf` (continuum ζ), MEDYAN
  `pcbi.1004877.pdf`, Stam `pnas.201708625.pdf`, Nedelec `msb.20177796.pdf`, Miyazaki, Ennomani
  `nihms803347.pdf`, Gowrishankar supplement `mmc5.pdf`. **Mis-identifications corrected:**
  `emss-72183.pdf` = Chugh 2017 (NOT Svitkina EM); `nihms803347.pdf` = Ennomani 2016 (NOT Svitkina EM);
  `2021.06.28.450262v2.full.pdf` = Serwas 2021 CME cryo-ET (off-target). The interim Nie 2015 substitute
  is NOT in the local corpus (only Chugh's ref 6) — sourced directly from PMC.
- **PubMed:** targeted "MCF7 cortical myosin minifilament density" / "breast cortical myosin areal
  density" → **0 areal-density hits**; metadata + full-text + `convert_article_ids` used to verify every
  DOI/PMID/PMC and to **catch wrong-PMID defects** (§11). Nie 2015 full text fetched from PMC4361371 and
  density quotes verbatim-verified.
- **Web cross-check:** PMC / journal pages / DOI resolver for Truong Quang 2021, Hu 2017, Beach 2017,
  and the Salbreux review candidates (paywalled bodies not exhaustively word-read — flagged where so).
- **Closest non-MCF7 "epithelial":** zebrafish periderm (van Loon/Sonal 2021, *MBoC*, PMC8351741) —
  reports microridge density + minifilament *spacing* (281 nm doublet), NOT a /µm² count; not usable.

---

## 11. Citation-integrity confirmation

**All DOIs/PMIDs/PMCs below resolve to the stated real papers; none is on the project hallucination
list (`Yao2011_NatCommun`, `YapKovacs_JCS`, `NanoConvergence2021_Glioma`).** Verified by **local
full-text** where marked, otherwise by **PubMed metadata / full-text / ID-conversion**.

**Confirmed-correct, independently re-tied (DOI + identifier + paper-states-value):** Nie 2015 (PMID
25641802 / PMC4361371 — density quotes verbatim from PMC full text); Niederman & Pollard 1975 (240861);
Beach 2014 (24814144 / PMC4108432); Fenix 2016 (26960797 / PMC4850034); Hu 2017 (28114270); Beach 2017
(28114272); Shutova/Svitkina **2014** (25131674 / PMC4160463) — correctly disambiguated from Shutova
**2017** *J Cell Biol* (28701425 / PMC5584186, a different stress-fiber paper with no conversion
number); Morone 2006 (16954349 / PMC2064339 — mesh medians verbatim); Henson 2017 (28057763 /
PMC5328620); Chugh-Salbreux 2017 (28530659 / PMC5536221); Truong Quang 2021 (34764258 / PMC8586027);
Salbreux-Charras-Paluch 2012 (22871642); Clark-Salbreux 2014 (24845681).

**⚠ Defects in the raw research findings — CORRECTED and the affected rows DOWNGRADED (never presented
as fact):**

1. **AFINES PMID 28746850 is WRONG** → correct **28746855** / PMC5529201 (28746850 = an unrelated UNG2
   DNA-glycosylase paper). DOI 10.1016/j.bpj.2017.06.003 is correct. Row kept but marked **UNCONFIRMED /
   low** and labeled illustrative-in-silico (not a cortex datum).
2. **Stam 2017 PNAS PMID 29114043 is WRONG** → correct **29114058** / PMC**5703288** (not 5703302)
   (29114043 = a subplate-neurons neuroscience paper). DOI 10.1073/pnas.1708625114 is correct. No areal
   density anyway → **excluded from the value table.**
3. **Gowrishankar 2012 PMID 22863013 is WRONG** → correct **22682254** (22863013 = an NMDA-receptor
   paper). DOI correct. Actin-mesh context only.
4. **Ennomani 2016 PMID off-by-one** → correct **26898468** (26898467 = a plant-mitochondria paper).
   In-vitro free-parameter myosin count → not a density anchor; **excluded from the value table.**
5. **Miyazaki 2015 author list FABRICATED** in a raw finding ("Eto, Cao, Iwasa, Suzuki, Yumura,
   Miyamoto") → real authors **Miyazaki, Chiba, Eguchi, Ohki, Ishiwata**; DOI 10.1038/ncb3142 / PMID
   25799060 correct. In-vitro, no areal density → **excluded from the value table.**
6. **Hosseini 2020 (PMID 33042748) γ-protocol UNCONFIRMED:** the paper is real (Fischer-Friedrich group)
   but its headline subject is **EMT / mitotic rounding**; the specific "interphase-suspended γ ≈ 0.27
   mN/m" should be re-confirmed from full text before it is leaned on as the ceiling comparator. Used
   here only as the gap reference, with this caveat. (The sibling tension-datum note reads γ ≈ 0.27 off
   Fig. 2g directly — cross-check there.)
7. **Hu 2017 "~5–20 minifilaments/stack" is an AGENT GEOMETRIC ESTIMATE,** not a Hu-stated value
   (PMID 28114270 states only "several micrometres" stacks) → labeled as an estimate, down-weighted.

**Honest negatives:** no MCF7/breast cortical areal density was found (medium confidence — a negative
cannot be search-proven); no platinum-replica EM minifilaments-per-focus integer is retrievable; the
Salbreux 2012/2014 review *bodies* are paywalled (verdict rests on abstract + structural framing +
prior audit, strong but not literal full-body exhaustion).

---

## 12. SE-registration candidates (NOT auto-registered — PI-gated)

| Source (corrected identifiers) | Proposed KnowledgeClaim |
|---|---|
| **`Nie2015_Cytoskeleton`** (10.1002/cm.21207 / PMID 25641802 / PMC4361371) | "Cortical NMII minifilament areal density ≈ 0.6/µm² (0.3–0.9), HeLa interphase adhered medial cortex; intensity-calibrated MRLC-GFP foci ≡ minifilaments (40 ± 20 per 8×8 µm)." → **proposed `areal_density_per_um2` anchor (non-MCF7 proxy, LOW–MEDIUM confidence).** |
| **`NiedermanPollard1975_JCellBiol`** (10.1083/jcb.67.1.72 / PMID 240861) | "One NMII bipolar minifilament ≈ 28 myosin molecules, ≈ 320 nm — single minifilament ≈ optical diffraction limit (foci→minifilament lower bound = 1 only where resolved)." |
| **`Beach2014_CurrBiol`** (10.1016/j.cub.2014.03.071 / PMID 24814144) | "TIRF-SIM resolves a single ~300 nm bipolar filament as two ~300-nm-spaced puncta; dense regions unresolvable — anchors the foci→minifilament conversion lower bound (1) and its dense-region failure." |
| **`Fenix2016_MBoC`** (10.1091/mbc.E15-10-0725 / PMID 26960797) | "NMIIA filaments expand (dominant) + concatenate (~10%) into stacks; a multi-motor-group filament is ONE expanding filament, not stacked — motor-group count ≠ stack count." |
| **`Hu2017_NatCellBiol`** (10.1038/ncb3466 / PMID 28114270) | "Registered NMII filament stacks span up to several µm — qualitative upper bound for foci→minifilament (>1 in dense regions); per-focus integer NOT measured." |
| **`Morone2006_JCellBiol`** (10.1083/jcb.200606007 / PMID 16954349) | "Cortical ACTIN membrane-skeleton mesh side-length median 200 nm (NRK)/52 nm (FRSK), range 20–500 nm (platinum-replica EM tomography) — geometric context, NOT a myosin density." |
| **`ChughSalbreux2017_NatCellBiol`** (10.1038/ncb3525 / PMID 28530659) | "Cortex tension model parameterizes myosin as a 2D line-concentration N_m/W (T₀ ≈ 230 pN/µm), NEVER a per-µm² areal density — confirms 'Salbreux 3.0/µm²' is misattributed." → **relabel the runtime's misattributed density provenance.** |

**No new SourceEvidence/KnowledgeClaim was registered.** Registration is a separate PI-gated step;
this doc is the candidate record and the surfacing artifact for the magic-number / citation-integrity
fix.

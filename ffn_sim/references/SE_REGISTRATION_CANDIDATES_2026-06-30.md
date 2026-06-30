# SourceEvidence registration candidates — 2026-06-30 (γ-floor closeout)

Per CLAUDE.md, new papers stay BM25-only until a Notion SourceEvidence row is PI-authored; this file
is the candidate list (NOT an auto-registration). PI approved promotion (Decision 2, 2026-06-30) — the
rows below are ready for SE creation + the linked param/contract corrections.

## 1. Nie 2015 — the cortical NMII minifilament areal-density datum (replaces the Salbreux misattribution)

- **citation_key:** `Nie2015_Cytoskeleton`
- **Authors/Title:** Nie, Y. et al. (verify full author list at registration) — cortical NMII
  minifilament density via intensity-calibrated MRLC-GFP super-res.
- **Journal:** Cytoskeleton 72(1):29–46 (2015)
- **DOI / PMID / PMC:** 10.1002/cm.21207 / 25641802 / PMC4361371
- **Value to anchor:** cortical NMII minifilament **areal density ≈ 0.625/µm²** (40 ± 20 foci per
  8×8 µm HeLa medial cortex; range 0.31–0.94; 1 focus ≈ 1 minifilament, intensity-calibrated).
- **Confidence:** LOW–MEDIUM (HeLa, non-MCF7; single study).
- **Links:** KU-3.5; supersedes the runtime's "3.0/µm² = Salbreux 2012" (CONFIRMED misattribution,
  H7_CORTICAL_MYOSIN_DENSITY_DATUM_2026-06-07; FF use in `gamma_floor.py:NIE2015_DENSITY_UM2`).
- **Verdict expectation:** OK (real paper, verifiable DOI+PMID+PMC).

## 2. Blebbistatin myosin-fraction anchors (band is ~70% myosin) — for the KU-3.5 verdict

These support the 2026-06-30 authoritative record (band is myosin-dominated; genuine passive floor
~0.04 mN/m). Register if not already present:

- `FischerFriedrich2016_BiophysJ` — mitotic HeLa parallel-plate AFM, blebbistatin −68%, Y-27632 −73%
  (1.67 → 0.54 / 0.45 mN/m). PMID (verify at registration).
- `Tinevez2009_PNAS` — blebbistatin reduces cortical tension >50% (L929). DOI 10.1073/pnas.0903353106.
- `Warmt2021_NewJPhys` — MCF-10A suspended, blebbistatin ~91% (0.013 → 0.0012 mN/m); already SE307.

## 3. Param/contract CORRECTIONS to apply with the SE rows (not new papers — fixes)

- **`turgor_dP0` = 133 Pa is BAND-IMPLIED/TUNED** (no sourced row; PARAM_AUDIT_SIMUCELL3D_2026-06-25).
  → Decision 1: register in `params_manifest.yaml` as `tuned-controlled-variable`, OR replace with a
  state-dependent osmotic law (pending the 2026-06-30 turgor workflow). NOT to be cited as sourced.
- **`g_rigid` ≈ 0.57 mN/m "myosin-independent passive"** → reclassify as a turgor pressure-PARTNER
  double-book (Young-Laplace ΔP·R/2 is the pressure's partner, not an independent passive channel).
  The genuine blebbistatin-insensitive passive floor is ~0.04 mN/m.

## 5. FF↔Kim shear-modulus (Stage 6c, task c) — crosslinker + actin axial stiffness (PI-GATED)

The FF↔Kim absolute-shear-modulus work (FF_KIM_NETWORK_VALIDATION §mechanics; viz
`outputs/ff/figs/kim_shear_modulus.png`) needs these LIT-ANCHORED stiffnesses. They are currently
**unsourced code constants / placeholders** — register SE rows BEFORE correcting the live runtime
values (a magic-number/contract trigger). NONE are in the KB yet.

- **`Ferrer2008_PNAS` — crosslinker JUNCTION stiffness** ✅ **PI-APPROVED + APPLIED 2026-06-30** (the
  paper already cited in `hand_kmc.py` for the α-actinin off-rate k_off0=0.066/s; its STIFFNESS datum
  was skipped). Ferrer, Lee, Chen, Pelz, Nakamura, Kamm & Lang 2008, PNAS 105(27):9221–9226, DOI
  10.1073/pnas.0706124105. **Values applied:** `ALPHA_ACTININ.link_k` = **4.6e5 pN/µm** (455 pN/nm),
  `FILAMIN.link_k` = **8.2e5 pN/µm** (820 pN/nm) (κm = bond-well curvature from Dudko–Hummer–Szabo, the
  load-bearing junction stiffness = exactly what `link_k` is). Was 0.1 = ~4.5e6× too soft (1000×
  pN/µm-vs-pN/nm slip + AFINES soft-surrogate). The stiff value is SAFE in production via the
  crosslink-turnover resting baseline (Stage 6N, `equilibrate(crosslink_turnover=True)`): the active γ
  is identical (floored 1.46e-4 mN/m), γ_xl clean (force-free crosslinks). **STILL TODO:** create the
  Notion SourceEvidence row for the Ferrer stiffness datum + (optionally) a `params_manifest` entry so
  kb-check tracks it (currently link_k is a code constant, not a tracked YAML param → no drift flagged).
  PDF NOT in references/ — fetch + DOI-verify at SE creation. Note: filamin has a SECOND regime
  (entropic WLC ~2 pN/µm, Broedersz–Storm–MacKintosh 2009); the Ferrer κm is the load-bearing one for
  `link_k`.
- **`Kojima1994_PNAS` — actin axial stretching modulus EA.** Kojima, Ishijima & Yanagida 1994, PNAS
  91(26):12962, DOI 10.1073/pnas.91.26.12962. Direct glass-needle stretch 43.7 ± 4.6 pN/nm over 1 µm ⇒
  **EA = 4.4e-8 N = 4.4e4 pN** (⚠️ 1000× trap: 4.4e4 pN, not 4.4e7). Used in `kim_network.EA_ACTIN_PN`
  to set k_axial = EA/L_seg (measurement only, not yet a live config constant). No Kojima SE row exists.
- **Gardel/MacKintosh absolute cross-linked-actin G' band 0.1–1000 Pa** — MEDIUM confidence, canonical
  but NOT KB-anchored; Gardel 2004 Science PDF is cited-only (not in references/). Obtain PDF + DOI
  before any SE. Used only as an on-demand comparison band, never a live constant.

Also fix (clean 1000× slip): `configs/phase1_h3.yaml` actin intra-filament `k_xl/k_intra = 1.0e-7 N/m
(= 0.1 pN/µm)` → **0.1 pN/nm = 1e-4 N/m = 100 pN/µm** (lands at the registered PARAM-k_xl floor); its
in-code "KU-3.19 (Furuike 2001)" attribution is wrong (KB-3.19 is a Bell off-rate; Furuike has 0 KB
hits) — fix the comment.

## 6. Lamellipodium / Arp2/3 (Stage 6L increment 2, task B) — PI-GATED cites

The dendritic-array builder (`ff/architecture_spec.LAMELLIPODIUM`, `weave._build_lamellipodium_patch`)
uses these. The branch GEOMETRY constants are already lit-anchored + Magic-Number-Blocked in
`configs/phase1_h5.yaml` (Fäßler 2020); the others need SE rows before a deliverable cite.

- **`Fassler2020_EMBOJ`** — Arp2/3 in-cell cryo-ET branch angle **68 ± 9°** (→ θ₀=70°, k_angle=0.173
  pN·µm/rad² via kT/Var(θ)). EMBO J 39:e104254. ALREADY in configs/phase1_h5.yaml (audit C5); the SD=9°
  should be a KnowledgeClaim (the KB has the 70° MEAN as a bare scalar, KB-3.7/3.18, but NOT the SD —
  the SD is the load-bearing datum for the angle-harmonic kernel). PI-gated KB change.
- **`Mueller2017_Cell`** — lamellipodium **±35° two-mode** filament orientation (protrusion-axis-relative;
  reproduces Maly-Borisy 2001). In `references/downloaded/Mueller2017_Cell.pdf`. SE row to cite.
- **`Vinzenz2012_JCS`** — Arp2/3 branch **density 1.25/µm** (1 branch / 0.80 µm contour), inter-branch
  36.7/71.2 nm helical-repeat quantization. JCS 125:2775, DOI 10.1242/jcs.107623. NOT in references/ —
  fetch + SE before cite.
- Cross-check oracle (do NOT adopt): Cytosim fork:angular_stiffness 0.076 pN·µm/rad² (Akamatsu/Berro);
  our 0.173 is from the in-situ Fäßler σ=9° equipartition (both in the 0.05–0.17 decade).

## 7. Stress-fiber / microvillus / FA-traction (Stage 6L increment 3) — PI-GATED

- **SF sarcomere periodicity ~1.0 µm** (band 0.5–1.4): `Hotulainen2006_JCB` (J Cell Biol 173:383) /
  `Tojkander2012_JCS` / `Peterson2004_MBoC` — NOT in references/ or KB (0 paper_chunks). The
  alternation PATTERN (α-actinin Z-bodies @ barbed ends ↔ central NMIIA bands) IS KB-grounded
  (`Murrell2015_NRMCB`, nrm4012.pdf). Register one primary as SE before the SF periodicity is a fixed
  datum; `STRESS_FIBER.sarcomere_um=1.0` is a SWEPT/provisional value.
- **SF N_filaments/cross-section 7–30** (categorical RANGE, no point value): Cramer 1997 JCB 136:1287 +
  Tojkander 2012 (range provenance). Register as RANGE, no point count. Default 20 = swept.
- **SF active single-fiber tension 5–6 nN**: `Kassianidou2017_PNAS` (PNAS 114:2622) — already in the
  prior SE_REGISTRATION_CANDIDATES_2026-06-09_sf-mechanics.md (PI sign-off pending). Network total
  10–30 nN (Kumar 2006) EMERGES, not back-solved.
- **Microvillus constants (genuine KB gap, 0 rows)** — `DeRosier_Tilney` (brush-border paracrystal,
  count ~20–30, ~12 nm lateral c2c), `Bartles_espin`, `Loomis2003_JCB`, `Lange` (length 1–2 µm /
  diameter 80–100 nm / rootlet). espin/fimbrin(I-plastin)/villin bundler KINETICS un-sourced → the
  builder uses anchored FILAMIN as a stand-in (flagged). ⚠️ resolve the 12 nm (lateral c2c, what the
  metric measures) vs 33 nm (actin helical repeat, a different axial period) before a passing gate.
- **FA integrin clutch**: `Kong2009_Nature` (Nature 185:1275… verify; α5β1-FN catch-slip, KB-2.5,
  audit OK) — already registered. ⚠️ KB-CONSISTENCY FLAG: the recorded KB-2.5 params (k_catch=0.4/s
  Fc=7pN, k_slip=0.5/s Fs=30pN) give a lifetime peak F* = ln[(k_catch·Fs)/(k_slip·Fc)]/(1/Fc+1/Fs)
  ≈ **7 pN**, NOT the "F*≈30 pN" the claim text states — the catch-bond IS present but the peak-force
  claim is inconsistent with its own params; surfaced to PI (used as-recorded in `INTEGRIN_A5B1`, no
  tuning). FA clutch k_int=1 pN/nm + k_on=1/s are Phase-defaults (KB-2.4/2.18 "relation pending"),
  PI-gate. ⚠️ Use Kong (integrin-ECM), NOT Rakshit 2012 / KU-4.2 (that is cell-cell cadherin).
- **Gil-Redondo 2023** per-cell traction 102 nN — literature overlay only, NOT a model-validation gate
  (H.7 SF-array HALTED).

## 8. Cortical-tension DEFINITION / framing (FF_STAGE6P, 2026-07-01) — PI-gated

The cortical-tension framing audit (4 workflows) found the band provenance + cell-type were mis-labeled.
Register/correct:
- **`Chugh2017_NatCellBiol`** (Nat Cell Biol 19:689, DOI 10.1038/ncb3525) — the REAL anchor of the
  working 0.35–0.65 band (HeLa, interphase, AFM, T₀=230 pN/µm). "Salbreux band" is a MISNOMER (Salbreux
  2012 is a review, broad 0.1–1). Likely already in KB paper_chunks; ensure an SE row + fix the label.
- **`Hosseini2020_AdvSci`** (Adv Sci 7:2001276, DOI 10.1002/advs.202001276) — MCF7 interphase SUSPENDED
  cortical tension 0.27 mN/m (IQR 0.18–0.40). Already KB (audit OK). The MCF7-faithful overlay; adherent
  MCF7 γ = documented ABSENCE.
- **`Bohec2025_Jour`** — interphase basal T₀=0.47 mN/m (the genuine interphase mN/m anchor for the
  decomposition). Verify citation.
- **`Wagner2011`** — MCF7 R_cell = 7.5 µm (the FF cortex uses 10 µm, mislabeled "MCF7 Wagner 2011" — PI
  decision to correct R; γ ∝ 1/R).
- Band is TOTAL (~70% active); the active-only FF γ must compare to ACTIVE_FRACTION·band (~0.245–0.455
  mN/m), now in `gamma_estimator.active_band_pn_um()`. Turgor ΔP·R/2 stays OUT (double-book, confirmed).

## 4. Harvest

`OPS_HARVEST_CANDIDATES_2026-06-30.md` lists 27 un-harvested RunResults (mostly prior-session H3/H7/
Layer2/DCM, all corroborating the floor: gamma_soft ~1e-5–1e-6 mN/m) + 2 code modules. The FF 6d–6h
runs are research prototypes (not production-REPORT format) so are not in the harvest set. The broad
`harvest_ops --apply` is idempotent but cross-cuts DCM territory while the /loop dcm session is live →
recommend a coordinated apply, not a unilateral one.

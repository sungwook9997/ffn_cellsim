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

## 4. Harvest

`OPS_HARVEST_CANDIDATES_2026-06-30.md` lists 27 un-harvested RunResults (mostly prior-session H3/H7/
Layer2/DCM, all corroborating the floor: gamma_soft ~1e-5–1e-6 mN/m) + 2 code modules. The FF 6d–6h
runs are research prototypes (not production-REPORT format) so are not in the harvest set. The broad
`harvest_ops --apply` is idempotent but cross-cuts DCM territory while the /loop dcm session is live →
recommend a coordinated apply, not a unilateral one.

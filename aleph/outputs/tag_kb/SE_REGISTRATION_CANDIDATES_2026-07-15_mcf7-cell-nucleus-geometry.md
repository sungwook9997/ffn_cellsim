# SourceEvidence / KnowledgeClaim candidates — MCF7 cell & nucleus GEOMETRY (distribution)

**Date:** 2026-07-15 · **Origin:** PI challenge to R_cell / R_nuc accuracy (mech-hierarchy S3 review).
**Status:** DRAFT — not yet Notion rows. PI sign-off required before `harvest_ops`/registration.

## Why this exists

The FF code carried `R_nuc = 0.70·R_cell` with an in-code citation **"Moore 2016, N:C 1.9"** that a
TAG query confirmed is **NOT in the KB** (`se.citation_key LIKE 'Moore2016%'` → zero rows) — a phantom
citation. Separately, cell/nucleus size was a hard-coded **point value**, but real single-cell data is a
**distribution** (~15–30 % CV). This registers the real primary sources so geometry is sourced +
sampled, not asserted. Consumed by `aleph/common/cell_geometry.py` (2026-07-15).

## SourceEvidence candidates

### SE-G1 — MCF7 cell volume (the volume anchor)
- **Citation:** Gamcsik MP, Millis KK, Colvin OM. *Cancer Res* 1995;55(10):2012–6. PMID 7743493.
- **Via:** BioNumbers **BNID 115154**; cited in Wagner BA, Venkataraman S, Buettner GR. *Free Radic
  Biol Med* 2011;51(3):700. **doi:10.1016/j.freeradbiomed.2011.05.024**.
- **Value:** MCF7 cell volume **1760 µm³ (1.76 pL)** → sphere-equivalent radius **7.49 µm** (³¹P/¹³C NMR;
  cells grown as spheroids / on collagen sponge). Volume is the conserved quantity under turgor → the
  default central radius. NOTE: this already backs the existing `R_cell = 7.5 µm` (so R_cell was correct).

### SE-G2 — MCF7 cell + nucleus size distribution, N:C ratio (imaging flow cytometry, n=2164)
- **Citation:** comparison study, **PMC7000884** (10.1371/journal.pone…, "Determination of cell
  nucleus-to-cytoplasmic ratio using imaging flow cytometry and a combined ultrasound+photoacoustic
  technique"). Suspended cells in PBS.
- **Values (IFC, n = 2164):** cell diameter **18.88 ± 2.86 µm** (r 9.44 ± 1.43), nucleus diameter
  **12.68 ± 1.94 µm** (r 6.34 ± 0.97), **N:C radius ratio 0.68 ± 0.08**.
- **Caveat (from the authors):** IFC masking systematically OVER-estimates absolute size (→ the
  volume/photoacoustic anchor is ~1.25× smaller); take the SHAPE (CV, N:C) from IFC, the MEAN from SE-G1.

### SE-G3 — MCF7 cell + nucleus size, N:C ratio (UHF ultrasound / photoacoustic, n=37)
- **Citation:** *Int J Thermophys* 2016;37:118. **doi:10.1007/s10765-016-2129-y** (ADS 2016IJT....37..118M).
  Suspended-in-agarose (immobilised spherical).
- **Values (n = 37):** cell diameter **15.2 ± 3.5 µm** (r 7.6 — agrees with SE-G1), nucleus diameter
  **10.2 ± 3.5 µm** (r 5.1), **N:C 0.68 ± 0.19**. Independent confirmation of BOTH the 0.68 N:C ratio and
  the ~7.5 µm radius.

## KnowledgeClaim candidate

**KB-3.Bx — MCF7 single-cell geometry (suspended).** R_cell (sphere-equiv) central **7.5 µm** (volume
anchor SE-G1; IFC method-upper 9.44 µm SE-G2), radius CV **~0.15**; nucleus:cell RADIUS ratio **N:C = 0.68
± 0.08** (SE-G2, corroborated SE-G3). Nucleus co-varies with cell (conserved N:C). Cross-refs: E_nuc 399 Pa
MCF7 (Fischer 2020, already KB-3.B2.1); this claim is GEOMETRY only, orthogonal to the modulus.

## Corrections this triggers (surfaced to PI)

1. **Remove the phantom citation.** `R_NUC_FRAC = 0.70 "(Moore 2016, N:C 1.9)"` in
   `scripts/ff_s3_nucleus_compression.py`, `ff_stiffness_sensing.py`, `_gbook_fullcomp_prod.py` → replace
   with the sourced N:C **0.68** (SE-G2/G3). The 0.70 value was ~right (0.68 measured); only the citation
   was fake.
2. **Resolve the R_nuc conflict.** `R_nuc = 0.25·R` in `_gbook_fullcomp_prod.py` (line 39) and
   `ff_cell_on_substrate.py` (line 43, already ⚠-flagged "audit#15, MCF7 nuc larger") gives nucleus dia
   3.75 µm = **1/3 of the measured 10–12.6 µm → wrong**; adopt N:C 0.68 (or sample).
3. **Resolve the E_nuc conflict.** `E_nuc = 4700 Pa` in `_gbook_fullcomp_prod.py` (line 5) vs the ratified
   MCF7 in-situ **399 Pa** (Fischer 2020, KB-3.B2.1) — 12× off; the 4700 is the superseded isolated-nucleus
   value (audit#18/19). Adopt 399.
4. **Geometry is now a DISTRIBUTION** (`common/cell_geometry.py`): quenched ensemble samples (R_cell, R_nuc)
   per realization; central run uses the mean for back-compat. S3 to be re-run sampled (deferred).

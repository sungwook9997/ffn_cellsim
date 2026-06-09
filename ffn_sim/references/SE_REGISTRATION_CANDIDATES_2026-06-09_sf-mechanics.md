# SourceEvidence registration candidates — STRESS-FIBER MECHANICS (2026-06-09)

From the 2026-06-09 deep-research audit (102-agent fan-out, 3-vote adversarial verification)
of the SF backbone-stiffness provenance, run because the PI flagged the N_filaments / EA / SF
data as old and possibly unreliable. Outcome: EA is a CONFIRMED foundational constant (corrected
attribution), N_filaments is genuinely a categorical range (no single measured value, even in
modern cryo-ET) → swept + Kumar-validated. Feeds `cell/stress_fibers.py` (citation corrected
this commit) and the B5 active-SF gate.

**Pipeline (BATCH, full-rebuild — NOT run now; run at the next clean KB batch in the MAIN repo,
the platform worktree has no kb.duckdb):** `references_ingest.py` → `verify_sources.py`
(CrossRef/web hard-rule) → `references_to_se.py` (Notion SourceEvidence) → `link_se_claims_*.py`
→ `refresh.sh`.

| # | citation_key | title | DOI/PMID | suggested KnowledgeClaim | verdict |
|---|---|---|---|---|---|
| 1 | Kojima1994_PNAS | Direct measurement of stiffness of single actin filaments with and without tropomyosin | 10.1073/pnas.91.26.12962 / PMID 7809155 | ⭐ **EA_single** — F-actin AXIAL rigidity: 43.7±4.6 pN/nm·1µm (bare) ⇒ EA=4.37e-8 N; 65.3 pN/nm (+tropomyosin) ⇒ 6.5e-8. THE axial datum for μ_SF=N·EA_single | OK (3-0) |
| 2 | LiuPollack2002_BiophysJ | Mechanics of F-Actin Characterized with Microfabricated Cantilevers | Biophys J 83(5):2705 | EA independent bracket: 34.5±3.5 pN/nm ⇒ 3.45e-8 N (phalloidin-F-actin). With #1 brackets EA≈1.1-4.4e-8 (chem-environment spread) | OK (3-0) |
| 3 | Kumar2006_BiophysJ | Viscoelastic retraction of single living stress fibers … (laser nanoscissor) | 10.1529/biophysj.105.071506 / PMID 16500961 | **SF tension validation band ~10-30 nN** (bovine capillary endothelial, in-situ; actomyosin-origin via ROCK/MLCK inhibition). The B5 acceptance band | OK (3-0) |
| 4 | Deguchi2006_JBiomech | Tensile properties of single stress fibers isolated from cultured vascular smooth muscle cells | 10.1016/j.jbiomech.2005.08.026 / PMID 16216252 | isolated-SF: re-stretch ~10 nN, breaking 377 nN, Young's 1.45 MPa (low-strain) → ~104 MPa strain-stiffening @200%. Corroborates the tension band + bundle modulus | OK (3-0) |
| 5 | Deguchi2005_MolCellBiomech | Flow-induced hardening of endothelial actin cytoskeleton (isolated SF preexisting tension) | PMID 16708474 | isolated endothelial SF: preexisting tension ~4 nN, strain 0.24, Young's ~300 kPa | OK |

## ⚠️ CITATION CORRECTION (citation-integrity hard rule) — already applied to code this commit
- **Gittes, Mickey, Nettleton & Howard 1993 J Cell Biol 120:923** was MIS-ATTRIBUTED as a source of
  the F-actin AXIAL rigidity EA in `cell/stress_fibers.py`. Gittes 1993 measures **FLEXURAL**
  rigidity (EI = 7.3e-26 N·m², persistence length ~17.7 µm) and its 1.2 GPa Young's modulus is for
  **MICROTUBULES**, not actin axial stiffness. EA is properly **Kojima 1994** (#1). Gittes 1993 may
  still be registered as the source for actin FLEXURAL rigidity / persistence length (a DIFFERENT
  claim), but must NOT back the EA_single / μ_SF axial datum. (If Gittes1993 is already a
  SourceEvidence row linked to an axial-EA claim, re-link it to flexural-rigidity only.)

## N_filaments landscape (NO single measured value — do not register as a point datum)
- Tojkander 2012 JCS 125:1855 + Kassianidou & Kumar 2015 BBA: "~10-30 actin filaments" is
  CATEGORICAL, uncited on the definitional sentence, both tracing to Cramer/Siebert/Mitchison 1997
  JCB 136:1287 (serial-section EM, locomoting heart fibroblasts) WITHOUT re-deriving a count.
- 2024 Frontiers review variant: "7-20 filaments under significant tensile stress."
- Modern cryo-ET (2020-2023; Structure 2021 actin-polarity, Struwwel Tracer 2023) images SF actin
  but publishes NO direct per-cross-section count. ⇒ N_filaments stays a SWEPT range (7-30),
  Kumar-validated; register Cramer1997 / Tojkander2012 as the RANGE provenance, not a value.

## CrossRef-verified at fetch (deep-research 3-vote)
- Kojima 1994: PNAS 91(26):12962, PMID 7809155 — abstract quote verified ("65.3±6.3 and
  43.7±4.6 pN/nm"). Liu&Pollack 2002: Biophys J 83(5):2705 — "34.5±3.5 pN/nm" verified.
  Gittes 1993: JCB 120(4):923, PMID 8432732 — "7.3×10⁻²⁶ Nm²" + "~17.7 µm" verified.

# SourceEvidence registration candidates — CONTRACTION-MECHANISM (2026-06-05)

Fetched 2026-06-05 (gbook KAIST) after PI flagged the in-vitro reconstituted-actomyosin line was
under-read. This is the contraction-MECHANISM paper that pinpoints the KU-3.5 γ-floor (see
`docs/ACTIN_ARCHITECTURE_NOTES.md` §16 + ★ MECHANISM SYNTHESIS). Stages for the KB pipeline.

**Pipeline (BATCH, full-rebuild tools — NOT run now to avoid colliding with the in-flight
`outputs/tag_kb/` working-tree changes; run at the next clean KB batch):**
`references_ingest.py` (tag_corpus + BM25) → `verify_sources.py` (CrossRef/web hallucination
hard-rule) → `references_to_se.py` (Notion SourceEvidence) → `link_se_claims_*.py` → `refresh.sh`.

| # | citation_key | title | DOI | file | suggested KnowledgeClaim | in_corpus |
|---|---|---|---|---|---|---|
| 1 | Miyazaki2015_NatCellBiol | Cell-sized spherical confinement induces the spontaneous formation of contractile actomyosin rings in vitro | 10.1038/ncb3142 | Miyazaki2015_NatCellBiol_ContractileRing_confinement.pdf | ⭐ actomyosin-contraction-mechanism (THRESHOLD in effective myosin density/oligomer + active REMODELLING/CONDENSATION, volume-conserving; confinement→cortex-like shell, equatorial ring=min bending energy) — DIRECT root-cause of the KU-3.5 γ-floor (our rigid M-SHAKE backbone forbids condensation → r/r0=1.000) | BM25-pending |

## CrossRef-verified at fetch
- DOI 10.1038/ncb3142 → Nature Cell Biology 17(4):480 (2015-03-23). Authors: Makito Miyazaki,
  Masataka Chiba, Hiroki Eguchi, Takashi Ohki, Shin'ichi Ishiwata (Waseda Univ). Real PDF (28 pp,
  2.1 MB) — not a 403 stub. Verified by reading abstract/results/discussion (not hallucinated).

## KnowledgeClaim links to set when registering (append-not-replace)
- → KB-3.5 (cortical tension / contractility) — the contraction-mechanism anchor
- → KB-3.23 (cortex fiber-network model) — buckling/remodelling requirement for the runtime
- (cross-link) Layer-2 KB-5.13 (spheroid surface tension) — single-cell γ source the σ-bridge consumes

## Related works CITED by Miyazaki / our corpus (already-registered or separate) — same mechanism family
- Lenz, Thoresen, Gardel, Dinner 2012 PRL 108:238107 "Contractile units … arise from F-actin buckling"
  (`Lenz2012_PhysRevLett`) — buckling symmetry-breaking theory + intermediate-density window.
- Murrell & Gardel 2012 PNAS 109:20820 "F-actin buckling coordinates contractility and severing in a
  biomimetic actomyosin cortex" (`Murrell2012_ProcNatlAcadSciUSA` / `murrell-gardel-2012-...pdf`).
- Stam, Freedman, Banerjee, Weirich, Dinner, Gardel 2017 PNAS "Filament rigidity and connectivity tune
  the deformation modes …" (`Stam2017_ProcNatlAcadSciUSA` / `pnas.201708625.pdf`).
- Ennomani … Nédélec, Théry, Blanchoin 2016 Curr Biol 26:616 "Architecture and connectivity govern
  actin network contractility" (`Ennomani2016_CurrentBiology`).
- Sakamoto, Miyazaki, Maeda 2023 Phys Rev Research 5:013208 "State transitions of a confined actomyosin
  system …" (10.1103/PhysRevResearch.5.013208, OPEN ACCESS — APS HTML-blocked bare curl; fetch via
  browser/Unpaywall next batch) — confined-actomyosin phase behavior, directly relevant; NOT yet in
  references.

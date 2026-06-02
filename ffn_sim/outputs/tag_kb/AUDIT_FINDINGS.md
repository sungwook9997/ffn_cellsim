# SourceEvidence hallucination audit — consolidated findings

Date: 2026-06-02. Scope: all **243** SourceEvidence rows in the Notion Contract-Graph
(built by prior LLM sessions WITHOUT the actual PDFs → hallucination risk).
Method: CrossRef existence/metadata check (pass 1) → rule-based adjudication of
flagged rows (pass 2) → web/Scholar verification of the residual unconfirmed
(`verify_sources.py`, `verify_sources_pass2.py`, background web-verify subagent).

## Bottom line

**3 confirmed hallucinations / 243 (~1.2%).** The KB is overwhelmingly real —
prior sessions did **not** fabricate papers wholesale. The dominant issue is
**metadata quality** (missing/wrong DOIs, year/journal drift, inconsistent keys),
not fabrication. But the 3 fabrications each prop up a real KnowledgeClaim, and
one carries a parameter value → flagged for PI.

| outcome | n | |
|---|---|---|
| Verified real (CrossRef author+year match) | 109 | OK |
| Real, DOI missing (CrossRef found it → backfillable) | 102 | suggest DOI |
| Real, metadata drift (wrong DOI/year/journal) | 18 | fix metadata |
| **Confirmed hallucination / fabricated-as-described** | **3** | **PI decision** |
| (pass-1 false alarms cleared in pass-2/web) | 11 | — |

## The 3 confirmed hallucinations (PI action required)

Each is a fabricated/mis-cited source sitting under a *real* claim. The claim's
underlying physics is generally fine; the **source** is fake. Real substitutes found.

### 1. `Yao2011_NatCommun` → KB-3.18, **KB-3.19**  ⚠️ parameter at risk
- Claimed: "Yao et al 2011, Nat Commun (alpha-actinin)", α-actinin Bell-Evans
  `k_off0 ≈ 0.4 s⁻¹` catch-bond.
- Reality: **no such Nat Commun paper exists.** (Yao 2011 = J Mol Biol viscoelasticity;
  Yao 2013 = PRL gelation.) The α-actinin/actin bond-lifetime physics traces to
  **Ferrer 2008 PNAS** (10.1073/pnas.0706124105) and **Miyata 1996 BBA**
  (10.1016/0304-4165(96)00003-7).
- **Risk**: KB-3.19 feeds the cortex crosslinker Bell-Evans unbinding; the
  `k_off ≈ 0.4 s⁻¹` was attributed to a fabricated source → **re-anchor the number
  to Ferrer 2008 and re-verify it’s the same value** (magic-number hard rule).

### 2. `YapKovacs_JCS` → KB-4.1
- Claimed: "Yap & Kovacs, adherens junction review (JCS)".
- Reality: **no Yap-&-Kovacs JCS review exists.** Real Yap junction reviews:
  Yap/Gomez/Parton "Adherens Junctions Revisualized" Dev Cell 2015
  (10.1016/j.devcel.2015.09.012); Yap/Duszyc/Viasnoff CSH Perspect 2018. KB-4.1
  (E-cadherin adherens-junction structure) is sound; just re-anchor the source.

### 3. `NanoConvergence2021_Glioma` → KB-6.2.3
- Claimed: glioma whole-cell **AFM, T98G vs U87 MG**; PMC8253861, Nano Convergence 2021;8:20.
- Reality: PMC8253861 is **real but a different paper** — Ketebo et al., filamin-A in
  U87 via **soft pillars** (not AFM, U87-only, no T98G). The cited content is fabricated
  (real PMC ID mis-applied). Re-source KB-6.2.3’s glioma AFM data to a genuine T98G/U87
  AFM study, or relabel to the pillar paper’s actual finding.

## Metadata fixes (safe — real papers, correct refs in hand)

- **18 metadata-drift rows** with corrections, incl. wrong DOIs:
  `Hosseini2021_BiophysJ`→10.1016/j.bpj.2021.05.006; `Liew2024_CMBE`→10.1007/s12195-024-00811-4;
  `Smelser2015_BMMB`→10.1007/s10237-015-0677-x; `Funk2021_eLife`→**Funk 2019** (year wrong, DOI was dead);
  `Bi2015_PRX`→year 2016 (PRX 6:021011 correct); `Mui2016_PNAS`→J Cell Sci 10.1242/jcs.183699;
  `MegeIshiyama`→CSH Perspect 10.1101/cshperspect.a028738; plus year/journal fixes for review keys
  (`Geiger_NRMCB`, `Parsons_FAReview`, `RocaCusachs_IntegrinReviews`, `SensPlastino_NRMCB`,
  `TrepatSahai_NRC`, `PolacheckChen_ARBE`, `KhalilFriedl_TCB`, `AllenTildesley`, `deGennes`).
- **102 missing-DOI rows** with a CrossRef-suggested DOI (see `source_audit` table,
  `suggested_doi`) → backfillable.
- **Malformed keys** to normalize: `K562_BBRC2019` (real: 10.1016/j.bbrc.2019.06.054),
  `Yang_NRMCB` (real: Jing Yang EMT guidelines, 10.1038/s41580-020-0237-9 — add year).

## Reproduce
```
conda activate ffn_sim
cd ffn_sim/outputs/tag_kb
python verify_sources.py           # pass 1 -> source_audit_report.md + source_audit table
python verify_sources_pass2.py     # pass 2 -> source_audit_pass2.md (+ cached dossiers)
# web-verify of residual unconfirmed: launched as a subagent (WebSearch/WebFetch)
```
Verdicts live in `kb.duckdb` table `source_audit` (uid, citation_key, verdict,
suggested_doi, audit2, audit2_reason) — queryable via `tag_query.py`.

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

# FF Stage 6T — proteomics density reconciliation: the cortex is NOT myosin-limited

**Date:** 2026-07-01  **Engine:** FF (Warp, A5000)  **Branch:** dcm/main
**Prompted by PI:** resolve the imaging(Nie 0.6/µm²) vs proteomics(3–24/µm²) >10× cortical-myosin-density
discrepancy by fetching the ACTUAL MYH9 copy number — and the PI's cross-engine observation that in DCM,
turning cortical tension up to realistic values prevented the spheroid from even assembling.

## The retrieved data (grounded in deposited supplementary tables)

| quantity | value | source | grounding |
|---|---|---|---|
| MYH9 (NMIIA) heavy chains / cell | **7.17e6** | Itzhak 2016 eLife (doi 10.7554/eLife.16950), HeLa | RETRIEVED_SUPPL (supp1 row P35579) |
| MYH9 (NMIIA) heavy chains / cell | **1.88e6** | Hein 2015 Cell (doi 10.1016/j.cell.2015.09.053), HeLa | RETRIEVED_SUPPL (mmc4) |
| MYH10 (NMIIB) | 2.9e5 (Itzhak) / 1.0e4 (Hein) — poorly constrained | same | RETRIEVED_SUPPL |
| hexamers / minifilament | ~14–30 (bracket) | Billington 2013 / Descovich 2018 / Weißenbruch 2021 | RETRIEVED |
| monomer (soluble, non-filamentous) fraction | up to ~50% | Shutova 2014 (REF52 fibroblast) | RETRIEVED_FULLTEXT |

**Factor-of-2 resolved (adversarial verify):** MS proteomics counts HEAVY CHAINS; 1 hexamer (the force unit)
= 2 heavy chains, so hexamers = copies ÷ 2. NMIIA dominates (IIB/IIC minor). Cross-dataset spread ~4×.

## The conversion → cortical areal density (my arithmetic, rounded cell R=7.5µm, area 707µm²)

- Total NMII hexamers/cell ≈ (1.88e6…7.17e6)/2 ≈ **0.94e6 – 3.6e6**.
- Total minifilaments/cell = hexamers ÷ (14–30) ≈ **3e4 – 2.6e5**.
- × filamentous fraction (~50%) × cortical fraction (rounded cell, cortex-dominant, ~0.2–0.5) ÷ 707µm²
  → **cortical density ≈ 10 – 90 minifil/µm²** (conservative end ~6–11; central ~40–90).

**⭐ This is 15–150× ABOVE Nie's imaging (0.6/µm²) and it SPANS/EXCEEDS the tension band (16–21/µm²).**
For Nie (0.6) to be right, the cortical fraction would have to be **~0.3–1.3%** of total NMII — implausibly
low for a rounded cell where the cortex is the dominant NMII structure. **Weight of evidence: Nie's
calibrated-fluorescence imaging UNDERCOUNTS the dense cortical minifilaments** (resolution limit; ~30% "too
close to analyze" at <110nm; medial-cortex optical section). The cortex is NOT myosin-protein-limited.

## Does proteomics density close the FF floor? (A5000 sweep, native actin, interphase engagement 0.65)

| ρ_myo [/µm²] | n_myo | γ_myo [mN/m] | floor vs MCF7-active | label |
|---|---|---|---|---|
| 0.625 | 442 | 6.45e-5 | **1953×** | Nie imaging |
| 6 | 4241 | 6.40e-4 | 197× | proteomics LOW |
| 16 | 11310 | 1.71e-3 | 74× | generation-bound density |
| 21 | 14844 | 2.28e-3 | 55× | generation-bound hi |
| 60 | 42412 | 6.47e-3 | 19.5× | proteomics MID |
| 90 | 63617 | 9.63e-3 | **13×** | proteomics HIGH |

**Using the proteomics density shrinks the floor ~150× (1953×→13×) but does NOT close it.** γ_myo is exactly
linear in ρ_myo (confirmed). The residual ~13–20× at proteomics density is the **FF network SCREENING**:
at the generation-bound density (16–21/µm², where the theoretical-max stress ½·n·f·ℓ = the band) the FF
network gives floor 55–74×, i.e. the connected inextensible network screens contraction by **~55–74×**
(σ/σ_dipole < 1, exactly the FF_STAGE6H buckling-network finding — screening, not Ronceray amplification).

## Decomposition of the ~1950× interphase floor (superseding "missing datum")

```
floor(Nie) ~1950×  =  density undercount (Nie 0.6 vs proteomics ~90; ~150×)  ×  FF network screening (~55-74×)  [× engagement ~1.5×, folded]
```

- **Density (~150×): RESOLVED direction** — the myosin protein is present (proteomics); Nie imaging undercounts.
  Cortex is NOT myosin-limited. (Caveat: cortical fraction is a bracket, not pinned — no clean cortex-vs-SF-vs-
  soluble split exists for a rounded epithelial cell; all fraction data is adherent fibroblast.)
- **Network screening (~55–74×): THE remaining real question** — is the FF connected-inextensible-network
  screening physical (real cortex screens too) or an FF-model transmission artifact (real cortex transmits
  better)? FF_STAGE6H measured it (σ/σ_dipole 0.2–0.6). This — NOT the density, NOT a missing datum — is now
  the open lever.

## Cross-engine convergence (the PI's DCM observation)

Independently, DCM (`DCM_GAMMA_CONTROLLED_SWEEP_2026-06-29`) found: **rising γ to realistic values DEFLATES
cells (Young-Laplace) faster than cohesion holds them → a free aggregate does NOT assemble; only a confluent
(pre-packed) aggregate facets.** And the FF→DCM bridge: FF's turgor-γ ≈ faceting onset (K=2.5e3); FF's active
γ is ~4400× too weak to facet. So the isolated-cell BAND is incompatible with the aggregating/tissue state
from BOTH sides: FF can't generate it, and DCM can't assemble at it. **The band (Chugh/Salbreux, from rounded/
mitotic isolated cells = free-surface, maximally-contractile) is likely the WRONG reference for the tissue-
contact operating point** — where cortical tension is down-regulated at cadherin contacts (Maître/Lecuit-Lenne;
tissue surface tension ≠ single-cell γ). The FF-generated low γ may be the CORRECT contact-state tension.

## Net updated verdict

1. The cortex is **NOT myosin-protein-limited** — proteomics (grounded) gives ~10–90 minifil/µm², band-level;
   Nie imaging (0.6) undercounts ~15–150×. The "missing engaged-density datum" framing is **superseded**.
2. At proteomics density the FF floor is ~13–20×, and that residual is the **FF network screening** (~55–74×,
   FF_STAGE6H) — the physical-vs-artifact status of that screening is the new open question.
3. The isolated-cell **band is likely the wrong reference state** for the tissue/aggregating operating point
   (FF-can't-generate + DCM-can't-assemble converge; contact tension is down-regulated).

## Open items / next
- **Resolve the network-screening physical-vs-artifact question** — validate FF transmission against a
  CONTROLLED reconstituted actomyosin network of KNOWN density (Kim/Gardel/Murrell; Linsmeier R≈0.085 datum):
  does a known-density network generate the FF-predicted (screened) stress, or more? Non-circular.
- **Pin the cortical fraction** — needs a cortex-vs-SF-vs-soluble split for a rounded epithelial/MCF7 cell
  (currently only adherent-fibroblast monomer-fraction ~50% exists).
- SE candidates (Notion, PI-gated): Itzhak 2016, Hein 2015 (MYH9 copy#); Shutova 2014 (monomer fraction);
  Truong-Quang 2021 (engagement). NOT auto-registered.

Sources (PubMed/supp): Itzhak 2016 doi 10.7554/eLife.16950; Hein 2015 doi 10.1016/j.cell.2015.09.053;
Shutova 2014; Billington 2013; Descovich 2018 PMC6004588; Truong-Quang 2021 PMC8586027; Nie 2015 doi
10.1002/cm.21207. Related: [[project-gamma-floor-likely-deficit]], FF_STAGE6S, FF_STAGE6H, DCM_GAMMA_CONTROLLED_SWEEP.

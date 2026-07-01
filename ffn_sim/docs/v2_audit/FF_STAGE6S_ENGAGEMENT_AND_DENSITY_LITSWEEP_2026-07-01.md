# FF Stage 6S — engaged-overlap mechanism (Truong-Quang) + cross-cell-type density lit-sweep

**Date:** 2026-07-01  **Engine:** FF (Warp)  **Hardware:** gbook A5000  **Branch:** dcm/main
**Prompted by PI:** "does the cortical NMII density cluster across cell types, or is it unmeasurable / a
sweep-to-data spec?" → a 6-angle literature sweep + making myosin ENGAGEMENT an explicit mechanistic
variable in the γ-floor model.

---

## Part A — cross-cell-type / cross-method density lit-sweep (6 angles, PubMed + bioRxiv + Consensus)

**Question:** do independent measurements of cortical NMII minifilament areal density cluster near Nie
2015 (~0.6/µm², HeLa), or reach the ~16–21/µm² the tension band would require?

**Finding: the values do NOT cluster — there is exactly ONE clean whole-cortex areal-density number
(Nie), and the indirect routes sit 1–2 orders ABOVE it and disagree with each other.**

| route | value [minifil/µm²] | cell / method | grounding |
|---|---|---|---|
| direct imaging — **Nie 2015** | **0.6** (0.31–0.94) | HeLa, calibrated MRLC-GFP (40±20 / 64µm²) | RETRIEVED full-text (PMID 25641802, doi 10.1002/cm.21207) |
| direct imaging — Beach17/Fenix16/Weißenbruch21 | **no areal value** | MEF/U2OS/COS-7, super-res | composition/spacing/relative only (numbers live in figures) |
| proteomics (copy# ÷ ~29 hex ÷ area × cortical frac) | **3–24** | copy# TRAINING-unverified, frac 0.3–0.5 assumed | soft (Kulak/Itzhak suppl. not fetched) |
| EM close-packing (spacing→density) | **39–80 (LOCAL)** | stack-interior / leading edge (160nm, <110nm) | derived upper bound, NOT cortex-average |
| geometric packing ceiling | **~65** | Truong-Quang footprint, fully parallel | derived ceiling |
| model systems (fission-yeast node counts) | **~1–10** (bridged) | S. pombe ring→areal bridge | order-of-magnitude sanity |
| **tension band requires** | **16–21** | γ ~ n·f·ℓ | — |

Reads:
- **Nie (0.6/µm²) is the only clean whole-cortex imaging datum in the literature** — every other imaging
  paper reports per-filament composition, ROI counts, relative density, or nearest-neighbour spacing, but
  no minifil/µm². Confirmed across 2022–2026 (no new areal density for any cell type; **none for MCF7/breast**).
- **The band (16–21/µm²) is ~26–34× above Nie**, at the HIGH end of the (soft) proteomics bracket, ~1/3 of
  the geometric packing ceiling, and achievable LOCALLY inside EM stacks — i.e. **near-maximal parallel
  packing**, physically possible but not a whole-cortex mean by imaging.
- **⚠️ The real open discrepancy is imaging (Nie 0.6) vs proteomics (3–24) — >10×.** Either Nie undercounts
  unresolved dense minifilaments (SIM can't resolve <110nm; ~30% "too close to analyze" at the edge) OR
  most NMII is non-cortical (stress fibers + folded-monomer pool; cf. Liu/Billington ATP-driven
  depolymerization). This is now the sharpest lever, not "a missing MCF7 datum."
- **Hexamers/minifilament corrected: ~28–30** (Weißenbruch, Beach, Fenix all "up to 30"; Descovich 2018
  full-text: 16–56 motor subunits = 8–28 hexamers) — higher than the earlier "15–30" bracket used.

---

## Part B — ⭐ the "engaged fraction" is MEASURED (Truong-Quang 2021), not missing

The single most relevant datum for the γ-floor's open lever: **Truong-Quang et al. 2021 (PMC8586027, HeLa)
shows cortical tension is set by the myosin–actin OVERLAP (engagement), NOT the total minifilament count.**
In low-tension INTERPHASE ~35% of cortical NMII lies OUTSIDE the actin cortex (bound by one end → does not
transmit contractile stress → non-engaged); this overhang → ~0 in high-tension MITOSIS. The
interphase→mitosis ≥3× tension rise (Fischer-Friedrich: 0.2 → 1.6 mN/m) is driven by **increased engagement
at ~constant myosin amount**.

**This corrects the project's prior "engaged density = unmeasurable missing datum" framing:** the engaged
fraction IS measured, and it is a **~1.5× lever (0.65 interphase → 1.0 mitotic)** — NOT the ~30× (density)
or the floor-closing factor.

### Implemented as a mechanistic variable (`gamma_floor.py`)
`ENGAGED_FRACTION_INTERPHASE = 0.65`, `ENGAGED_FRACTION_MITOTIC = 1.0`; `measure_gamma(engaged_fraction=…)`
and `gamma_floor_run/ensemble/production` thread it. The FF 2D cortical shell does not resolve the RADIAL
overhang geometrically, so engagement enters as the measured scalar on the force-bearing count — and because
**γ_myo is exactly linear in the contributing-dipole count** (verified FF_STAGE6Q), the scalar is exact, not
a fudge. Production DEFAULT = the physiological interphase 0.65 (physiological-baseline rule). Test:
`test_engaged_fraction_scales_gamma_myo_truong_quang`.

### A5000 result — engagement DEEPENS the interphase floor (does NOT close it)

| state | engaged_fraction | γ_myo [mN/m] | floor vs MCF7-active |
|---|---|---|---|
| **interphase (physiological default)** | 0.65 | **4.82e-5** | **~2613×** |
| mitotic (fully engaged) | 1.0 | 7.42e-5 | ~1699× |

The ratio is exactly 1/0.65 = 1.54× (7.42/4.82). At the physiological resting engagement the floor is
DEEPER (~2613× vs the prior implicitly-mitotic ~1699×). **Making the model mechanistically faithful to the
best engagement datum CONFIRMS the floor — engagement is a real but ~1.5× lever.**

---

## Net verdict (updated)

- The cortical-tension DEFINITION + framework are resolved; the γ-floor is CONFIRMED and now more faithful
  (engagement is an explicit, lit-anchored, mechanistically-correct variable — tension ∝ overlap, per
  Truong-Quang).
- The floor is **~1700× (mitotic/full) to ~2600× (interphase) under the MCF7-active band** on the best
  whole-cortex imaging density (Nie) at physiological engagement.
- **Two honest open items** (neither closable by tuning):
  1. the **imaging↔proteomics >10× density discrepancy** (Nie 0.6 vs proteomics 3–24/µm²) — resolvable by
     fetching the MYH9 copy-number (Kulak 2014 / Itzhak 2016 suppl.) + a measured cortical NMII fraction;
     only the proteomics HIGH end (≈15–24/µm²) reaches the band, and even then only with ~50% cortical fraction.
  2. no MCF7/breast cortical minifilament-density datum exists (documented absence, re-confirmed 2022–2026).
- Engagement is now OFF the "unmeasured lever" list — it is measured (~1.5×) and does not close the gap.

## Sanity gate
- **Dimensional/linearity:** γ_myo = engaged_fraction × (raw γ_myo); raw γ_myo linear in f_myo → scalar exact. ✓
- **No magic number:** engaged_fraction set to the MEASURED Truong-Quang values (0.65/1.0), swept as a
  controlled variable; production default = physiological interphase; it DEEPENS (not shallows) the floor. ✓
- **Physiological baseline:** production defaults to the interphase (resting/adherent) engagement. ✓

## Files
- `ff/gamma_floor.py` — ENGAGED_FRACTION_* constants; engaged_fraction threaded through
  measure_gamma/run/ensemble/production (default 1.0 in measure, 0.65 in production).
- `tests/ff/test_gamma_floor.py` — `test_engaged_fraction_scales_gamma_myo_truong_quang`.

Related: [[project-gamma-floor-likely-deficit]], FF_STAGE6Q_LITFAITHFUL_DENSITY, FF_STAGE6P_CORTICAL_TENSION_DEFINITION.
Sources (PubMed): Nie 2015 doi 10.1002/cm.21207; Truong-Quang 2021 PMC8586027; Weißenbruch 2021 doi
10.7554/eLife.71888; Beach 2017 doi 10.1038/ncb3463; Fenix 2016 doi 10.1091/mbc.E15-10-0725; Descovich 2018 PMC6004588.

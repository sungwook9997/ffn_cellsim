---
kb_record:
  doc_id: FF_STAGE6P_CORTICAL_TENSION_DEFINITION_2026-07-01
  title: What IS the reference cortical tension? — γ-floor framing audit (active-vs-total, turgor, MCF7) — floor PERSISTS
  authoritative_as_of: 2026-07-01
  supersedes: []
  updates: [CORTICAL_TENSION_RECORD_2026-06-30, FF_STAGE6O_NATIVE_GAMMA_2026-07-01]
  status: current
---

# FF Stage 6P — what exactly is the reference cortical tension? (PI framing audit)

**Date:** 2026-07-01 · **Engine:** `ffn_sim/ff/` · branch `dcm/main` · **PI questions** (turgor; spheroid inner/outer;
"is cortical tension myosin-only?"; "what EXACTLY is the reference tension?"; "are the densities MCF7?")

Four parallel ultracode workflows (engaged-density · density-sweep · MCF7-anchoring · cortical-tension-
definition; ~3.2 M agent-tokens) audited whether the γ-floor compares the RIGHT quantity. **Verdict: the
γ-floor was mis-NAMED but is NOT a category error in its conclusion — it PERSISTS (~395× under) under
the corrected, MCF7-faithful, turgor-correct framing.**

## 1. What the measured cortical tension IS (precise definition)

The experimentally-measured cortical tension = **TOTAL effective surface tension** — the single in-plane γ
an instrument (AFM parallel-plate / micropipette) inverts from the Laplace balance ΔP=2γ/R on a rounded
cell. It never returns a decomposed component; "active" is a downstream blebbistatin-subtraction inference.
Decomposition (on Bohec 2025 interphase basal T₀ = 0.47 mN/m):

| channel | fraction | mN/m | inside the single apparent γ? |
|---|---|---|---|
| ACTIVE (actomyosin) | ~50–90% (central ~70%) | ~0.25–0.45 | YES — dominant |
| PASSIVE (cortex elastic/architecture) | ~9–12% | **~0.04** (blebbistatin-insensitive floor) | YES — thin floor |
| TURGOR-balanced (ΔP·R/2) | — | 0.20 (dP0=40 Pa, R=10µm) | **SEPARATE — Laplace PARTNER, not additive (double-book)** |
| MEMBRANE (lipid bilayer) | <10% | ~0.003–0.03 | separate, minor |

So cortical tension IS mostly myosin (~70%, blebbistatin −68% Fischer-Friedrich 2016 / >50% Tinevez 2009 /
~70–80% Warmt 2021 [the earlier "91%" was a citation error]). The "myosin-only is suspicious" worry resolves:
it is ~70% myosin + a thin ~0.04 passive floor; turgor and membrane are separate bookkeeping.

## 2. The band's TRUE source (a misnomer fixed)

The working band 0.35–0.65 mN/m is **TOTAL** tension, **interphase**, **HeLa**, AFM — and its real empirical
anchor is **Chugh et al. 2017 Nat Cell Biol 19:689** (T₀=230 pN/µm, coherent-sim peak T/T₀≈1.6 ⇒ ~0.37 mN/m),
NOT Salbreux. **"Salbreux band" is a MISNOMER** — Salbreux 2012 is a *review* carrying the broad 0.1–1 mN/m
range (Salbreux is a Chugh coauthor → the conflation). Chugh's own mechanistic sim REACHES 0.37 mN/m via a
myosin-stall-anchored T₀ — so the band IS achievable by a correctly-loaded cortex.

## 3. Is the γ-floor a category error? — mis-NAMED, conclusion stands

FF γ_active (method-of-planes over the actomyosin load path, `turgor=False`) is the ACTIVE in-cortex stress
only; it was compared to a band that is TOTAL. That is an apples-to-oranges **label**. The CORRECT comparison
is **FF(active) vs the ~70% active fraction of the band (~0.245–0.455 mN/m)** — NOT FF(active+turgor) vs total
(forbidden double-book). Native FF γ_active ≈ 6.2e-4 mN/m:
- vs full band lo 0.35 → ~565×
- **vs 70%-active target 0.245–0.455 → ~395×** ← the corrected, honest floor factor
- FF_total = active + genuine passive (6.2e-4 + 0.04) = 0.0406 → still ~9× under band lo

The framing correction tightens the gap ~1.4× and **does not close it**.

## 4. Turgor — double-book, osmotic (NOT myosin-generated), spheroid-interior inapplicable

- ΔP·R/2 is the Young-Laplace **PARTNER** of the apparent tension (ΔP=2γ/R; Stewart 2011, Taubenberger 2020),
  not an independent passive channel — adding it double-counts. The old "passive carries the band via
  ΔP·R/2≈0.665 mN/m" was **circular** (dP0=133 Pa was band-implied); re-anchored dP0=40 Pa (Fischer-Friedrich
  2014) gives only 0.20 mN/m, and even that is the partner. CONFIRMED double-book.
- **Myosin does NOT generate the turgor.** The bulk pressure is set by OSMOTIC pump-leak (Π_in≈5e5 pN/µm²,
  net excess ~40 Pa; orders of magnitude above cortex contractility). Myosin modulates the balance; it does
  not source the pressure. So cortical tension is NOT "turgor-mediated-via-myosin" — the active actomyosin γ
  and the osmotic ΔP are separate physics; the active channel is independently force-magnitude-floored.
- **Spheroid interior cells have NO free Laplace surface** (turgor balanced by neighbor contact, Newton-3) →
  ΔP·R/2 is doubly inapplicable there; a 3D spheroid has a TISSUE-level surface tension (Steinberg DAH,
  Manning-Foty) distinct from single-cell cortical tension — that is the DCM/tissue domain, not this single-cell γ.

## 5. MCF7 anchoring — ZERO floor inputs are MCF7 (softens ~1.5×, floor persists)

| input | FF value | true cell type | MCF7 datum |
|---|---|---|---|
| target band | 0.35–0.65 (Chugh) | **HeLa** interphase | MCF7 = **Hosseini 2020 0.27 mN/m** (suspended, IQR 0.18–0.40); adherent = **ABSENT** |
| myosin density | Nie 0.625/µm² | **HeLa** | none (no breast/any-line areal datum) |
| actin density | 38000→30/µm² | generic | none (also vs KB's own ~100/µm²) |
| R_cell | 10 µm | generic, **mislabeled "MCF7 Wagner 2011"** | **7.5 µm** (Wagner 2011) |
| η cytoplasm | 65.9 Pa·s | **MCF7 ✓** (Dessard 2024) | the only true-MCF7 input — doesn't enter the floor |

MCF7-correcting BOTH sides: target ×0.5–0.7 (Hosseini), R 10→7.5 (γ ×1.3–1.8), density (no rescue — HeLa is
the conservative-favorable choice; the band needs ρ≈16–21/µm² = 26–35× Nie, sourced nowhere). Floor goes from
~565× to still **~hundreds× under**. ⭐ **State mismatch:** all MCF7 γ is suspended/rounded; the FF/H.7 operating
point is spread-adherent → even the "right" MCF7 number is for the wrong state (adherent MCF7 γ = documented
absence). FF is a suspended sphere, so the suspended Hosseini 0.27 is the closest faithful target.

## 6. Density / engagement lever — cannot close it (independent confirmation)

- **No sourced load-engaged density above Nie** (any cell line). Truong-Quang 2021 is an overlap FRACTION
  (~35% overhang interphase → ~0% mitosis), not a density, and points DOWN ("one-end-bound … would not
  strongly contribute to tension"). Duty ratio LOWERS effective density (engaged = duty×present, duty<1).
- **A5000 density sweep:** γ ∝ ρ, but band needs ρ≈1792/µm² (≈2867× Nie, geometrically impossible); even
  40/µm² (64× Nie) is ~50× under — the inextensible network SCREENS the motor force ~100× below the
  active-gel envelope ½·n2D·f·ℓ. **FF already runs 5× above Nie (3/µm²) yet floored ~530×** — decisive.

## Bottom line

The γ-floor **survives every reframing**: corrected active-vs-active comparison (~395×), MCF7-faithful band
(Hosseini, still hundreds×), turgor (double-book, not additive; osmotic not myosin), density/engagement
(no datum reaches band; even impossible densities fall short via screening). It is a **real force-magnitude /
load-engaged-motor-density limit**, robust to all network mechanisms. The ONE irreducible open lever is the
**experimentally-measured load-engaged NMII density at the MCF7 adherent operating point — a documented
EXPERIMENTAL ABSENCE**, the same missing datum as the SF/NMII twin. NOT a tuning target.

## Actionable corrections (lit-anchoring/framing — NOT tuning; PI to ratify)

1. **Rename the band's provenance**: "Salbreux band" → Chugh 2017 (HeLa interphase total, AFM). Keep the
   value; fix the misnomer comment.
2. **Compare γ_active to the ~70% active fraction (~0.245–0.455 mN/m)**, not the full total band — document
   the ~395× corrected factor in the gate/record.
3. **Add the MCF7-faithful overlay**: Hosseini 2020 0.27 mN/m (suspended); flag the spread-adherent absence.
4. **R_cell**: the production cortex is mislabeled "MCF7" at 10 µm; MCF7 = 7.5 µm (Wagner 2011) — PI decision
   (foundational geometry; γ ∝ 1/R so it raises the floor target slightly).
5. **Turgor**: keep ΔP·R/2 OUT of the actomyosin γ (already is) — confirmed double-book.
6. Per-motor stall config divergence (gamma_floor 5 pN/side vs prod 56 pN/side) — reconcile.
SE candidates: Chugh2017, Hosseini2020, Bohec2025, Wagner2011 (`SE_REGISTRATION_CANDIDATES`).

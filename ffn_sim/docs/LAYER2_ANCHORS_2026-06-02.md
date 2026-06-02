# Layer-2 CBM — literature anchor provenance (dossier wf_6369ec0f-1bb, 2026-06-02)

> Distilled output of the MCF7-specific deep-research dossier (run `wf_6369ec0f-1bb`;
> 6 angles, 18 sources, 78 claims → 25 verified → 6 findings) plus the prior
> parameter-extraction pass. This is the **provenance + ratification record** behind the
> values in `configs/layer2_cbm.yaml`. Companion: `LAYER2_MULTICELL_DESIGN.md` (brief),
> `PI_EXP_VALIDATION_MAP.md` (Layer-2 framework, already-anchored non-MCF7 bands).
> HARD rule: literature-first; PI A/A₀ is overlay-only, never a fitting target.

## Status of the two MCF7 gaps the PI gated the build on

| Gap | Status | Value | Source / note |
|---|---|---|---|
| (1) MCF7 cell size → `r₀` | ✅ **CLOSED (measured)** | diameter **14.8 µm** (vol 1.70 pL) → use **15 µm** | Wagner 2011, *Free Radic Biol Med* 51:700, **PMC3147247** (Coulter Counter, MEASURED, open-access). 15 µm = low end of measured 15–20 µm band. |
| (2) MCF7 cell–cell cohesion → `D_e` | ✅ **CLOSED (measured)** | force-anchor: detachment **6–7 nN** @120s ⇒ `D_e = 2·F·contact_zone` | Iturri 2020 *Cells* PMC7227807 (MCF7-MCF7 SCFS, OA); Omidvar 2014/2016 corroborate. Direct-search "absence" (dossier) resolved via OA equivalents + source-authority audit. |

## (1) Cell size — resolved

- **Cite Wagner 2011 primary 14.8 µm / 1.70 pL**, NOT the BioNumbers-rounded 15.0 µm /
  1.76 pL (~3.5 % volume inflation). Conclusion (~15 µm) unchanged.
- REFUTED (do not cite): EMD Millipore Scepter brochure 15–17 µm (BNID 108926, 0-3);
  Gamcsik 1995 as the "primary" volume source (1-2, value incidental to an NMR study).
- Broader corroboration (all ≥ 15 µm): impedance cytometry 17.3 / 18.2 µm; imaging up to
  19.7 µm. So 15 µm is conservative. **Band-position within 15–20 µm is a documented
  modeling choice** (we take the low end); revisit if cell-volume sensitivity matters.
- Config: `cell.diameter: 1.5e-5` → `R_cell = 7.5 µm` (resolved).

## (2) Cell–cell cohesion `D_e` — RESOLVED via source-authority audit (measured)

The deep-research dossier reported a direct-search "absence," but a follow-up
**source-authority audit** (which is *whose* number is correct, not which value to pick)
recovered the measured anchor from open-access equivalents and corrected a misattribution.

**Authoritative anchor (MEASURED, OA):** **Iturri et al. 2020 *Cells* 9(4):935 (PMC7227807)**
— MCF7-MCF7 symmetric single-cell force spectroscopy: **de-adhesion force ~6–7 nN at 120 s
(mature), ~1–2 nN at 5 s (nascent); work ~180–200 fJ @120 s.** Corroborated by Omidvar 2014
*J Biomech* (PMID 25169659; ranking MCF7>T47D>MDA, nN-scale, figures paywalled — Cloudflare,
not gbook-recoverable) and Omidvar 2016 *Mol Cell Biochem* (MCF7 1.42–2.85 nN).

**How D_e is set (force-anchored):** the Morse MAX attractive force is pinned to the
measured detachment force, so cohesion is **nN-scale and balances nN traction** (the physics
constraint):
```
D_e = 2 · F_detach · contact_zone  = 2 · 6.5 nN · 1.5 µm ≈ 1.95e-14 J   (Morse max |F| = 6.5 nN)
```
Maturation range available: 1–2 nN (5 s) → 6–7 nN (120 s). The measured *work* (~200 fJ)
is a scale cross-check; it exceeds `2·F·contact_zone` because SCFS work integrates membrane-
tether pulling beyond the CBM contact range.

**⚠️ RETIRED — what was WRONG (the user's instinct, confirmed):**
- The cadherin `D_e = N_cad·⟨F⟩·Δx* = 1.2e-17 J` seed was **~4 orders of magnitude too
  small** (a pN single-bond product; nN traction would shred it).
- **Buckley 2014 = MISATTRIBUTION.** It is the αE-catenin–**actin** catch bond (intracellular),
  NOT the E-cadherin homophilic extracellular bond. So "Δx*=4 nm citing Buckley 2014"
  (KU-4.17 / `phase1_unit4_1.yaml`, propagated into `PI_EXP_VALIDATION_MAP.md`) is wrong on
  the molecule. The 4 nm value is not an E-cadherin parameter.

**E-cadherin catch-bond primary for the future L2.5 (KU-4.2) upgrade:** **Rakshit/Sivasankar
2012 *PNAS* 109:18815 (PMC3503169)** — biphasic catch (X-dimer, peak ~30 pN; the true source
of "30 pN") + slip (strand-swap, zero-force lifetime 0.63 s ⇒ k_off0 ≈ 1.6 s⁻¹), **sub-nm
x_β** (from the supplement / Manibog follow-ups). Implement as a two-pathway catch bond, not
a Bell slip-only proxy; strike the 4 nm.

## NEXT scale issue identified (analogous to D_e): per-cell MIGRATION drag
`γ = 6πηR` (water-Stokes, 9.77e-8 N·s/m) is correct for *thermal relaxation* (G1) but WRONG
for *cell migration* (nN traction / water-γ ⇒ unphysical ~mm/s). L2.2 quantitative spreading
needs a cell-migration drag (cytoplasm/substrate friction; Chen&Zou μ_cell/μ_s scale) —
another scale-bridge parameter to anchor before the A/A₀ sweep is quantitative.

## PI ratification queue
1. ~~D_e derivation~~ → **RESOLVED**: measured Iturri 2020 nN force-anchor (above). FYI.
2. **Fix `PI_EXP_VALIDATION_MAP.md`** (shared doc, awaiting your relay/approval): (a) strike
   the Buckley-2014 Δx*=4 nm E-cadherin attribution (it's αE-catenin–actin); (b) retire the
   pN-scale `N_cad·⟨F⟩·Δx*` cohesion seed for the Iturri nN force-anchor; (c) re-base the
   catch-bond on Rakshit 2012. (These are hallucination-class errors, same as the prior
   SourceEvidence audit — I can fix them on your OK.)
3. **Surface-tension validation target** — no MCF7 value; emergent-only vs non-MCF7 proxy
   (MCF10DCIS ~21 mN/m, Nagle 2022) under proxy-flag.
4. **Cell-size band position** — 15 µm (low end) vs mid-band 17–18 µm.

## Cheap follow-ups (not blocking)

- **gbook-recover Omidvar 2014** (PMID 25169659) → measured MCF7-MCF7 de-adhesion force +
  work → converts `D_e` from derived to measured, and gives the G5 ranking magnitudes.
- **MCF7 cortical tension** still unfound (Moazzeni 2021 pools 4 lines, refuted 1-2) —
  needed for the surface-tension = f(adhesion, cortical tension) cross-check.

## Other reconciliations

- **Whole-cell modulus** stays at **0.26 kPa** (Yubero 2020). The 2025 MCF7 AFM paper
  (Pørtner, *J R Soc Interface*) is **nuclear** microrheology (|E*| 70–170 Pa) — a
  different measurand, not a substitute. (A separate reconfirmation attempt was refuted
  0-3, so treat 0.26 kPa as single-method, not consensus.)
- **A/A₀ = a + b/R + c/R² has no literature analog** (candidate wetting laws refuted 0-3):
  genuinely novel. Platform job = reproduce the traction–cohesion partition the a/b/c terms
  encode, not match a literature curve. (Settled in `PI_EXP_VALIDATION_MAP.md`.)

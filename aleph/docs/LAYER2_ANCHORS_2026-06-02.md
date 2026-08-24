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

## L2.5 cadherin catch-bond anchor (Rakshit 2012, confirmed full-text 2026-06-03)

Source: **Rakshit, Zhang, Manibog, Shafraz, Sivasankar 2012, "Ideal, catch, and slip bonds
in cadherin adhesion," PNAS 109(46):18815 (PMC3503169, DOI 10.1073/pnas.1208349109)** —
single-molecule AFM force spectroscopy of E-cadherin. Confirmed from the open-access full
text (via PubMed):

| Conformation | Bond type | Measured | → kinetic anchor |
|---|---|---|---|
| **X-dimer** | **catch** | lifetime ↑ with force, **peak ~30 pN**, then ↓ | catch-peak **F* = 30 pN** |
| strand-swap | slip | lifetime ↓ with force; intrinsic τ₀ **0.63 s** | k_off0 ≈ **1.6 s⁻¹** |
| WT (0.3 s contact) | ideal | force-insensitive, τ **0.04 s** | (not used for cohesion) |

**⚠️ OPEN PARAMETER QUESTION (PI ratification — hard-rule "no unanchored magic number"):**
Rakshit fit a **sliding-rebinding model**; they state the **two-pathway (Pereverzev) form did
NOT fit** their data. But the project's *ratified* catch-slip MECHANISM is Pereverzev
two-pathway (KU-2.5 / PHASE_0_3 §D2). Two honest options for the L2.5 cohesion bond:
  (A) Use the Pereverzev form (sanctioned mechanism) parameterised to reproduce the two
      MEASURED observables — catch peak F*=30 pN and k_off(0)≈1.6 s⁻¹. This pins 2 of the 4
      (k_s,F_s,k_c,F_c); the remaining 2 are constrained by a documented choice (e.g. set
      F_c from F* via the F* identity + fix the slip x_β to the integrin KU-2.5 scale),
      flagged as a two-pathway REPRESENTATION of the Rakshit catch-slip, NOT a fit to their
      data.
  (B) Implement the Rakshit sliding-rebinding closed form directly (most faithful to the
      measurement, but a NEW oracle outside the ratified Pereverzev vocabulary).
Recommendation: (A) for the first L2.5 pass (reuses the validated `pereverzev_k_off` oracle +
`integrin_bonds.py` updater pattern), with F*=30 pN as the anchored physical observable; defer
(B) as a fidelity upgrade. **Needs PI sign-off before the cohesion magnitude enters results.**

The KEY physical point for fragmentation: the catch regime means cohesion **strengthens up to
~30 pN tension** — exactly the regime proliferation drives — so a catch-bond cohesion should
resist the L2.4 fragmentation that the static Morse well could not.

## Sliding-rebinding model (Lou & Zhu 2007) — equations extracted 2026-06-03

Source: **Lou JZ, Zhu C (2007) "A structure-based sliding-rebinding mechanism for catch
bonds," Biophys J 92(5):1471-1485** (PMID 17142266, **PMC1796828**, DOI
10.1529/biophysj.106.097048). [via PubMed]. (⚠️ NOT PMC1796818 — that's a PGK paper, the
`1709.pdf` mix-up; the references PDF still needs swapping.)

Two-pseudoatom catch-bond model. States: (1,1) both pairs bound; (1,0) one pair bound at the
original site; (0,1) a new interaction formed after sliding; (0,0) dissociated. Master eqs:
```
dp11/dt = -2 k-1(f) p11 + k+2 p01
dp10/dt =  2 k-1(f) p11 - k-1(f/2) p10 - k+1 p10 + pn k-1(f/2) p10
dp01/dt =  pn k-1(f/2) p10 - k-1(f/2) p01 - k+2 p01
dp00/dt = (1-pn) k-1(f/2) p10 + k-1(f/2) p01
```
- Bell dissociation of a single interaction: `k-1(f) = k0 · exp(f·a / kB T)` (a = length scale).
- Force shared equally between two interactions ⇒ each feels f/2 (the `k-1(f/2)` terms).
- New-interaction (rebinding-after-slide) probability: `pn = (f/f0)^2 / (1 + (f/f0)^2)` —
  rises 0→1 over force scale f0 (the hinge-opening threshold). This force-INCREASING rebinding
  is what creates the catch (longer lifetime under force) before Bell slip takes over → biphasic.
- Mean lifetime τ(F) = ∫ S(t)dt from the survival probability S(t)=Σ p_bound(t) (solve ODEs or MC).

Tabulated SELECTIN params (reference, NOT cadherin): P-selectin k0=30/s, k+1=30/s, k+2=1300/s,
f0=12 pN; L-selectin k0=70/s, k+1=30/s, k+2=500/s, f0=59.5 pN.

**E-cadherin X-dimer params = Rakshit 2012 SI Appendix Table S9 (STILL NEEDED).** Main text
anchors: catch peak F* ≈ 30 pN; strand-swap slip τ0 = 0.63 s; WT ideal τ = 0.04 s.
Implementation plan (option B, PI-chosen): a `sliding_rebinding` validation oracle computing
τ(F) from these ODEs, parameterised by the E-cadherin Table-S9 fit; the runtime Layer-2
`CadherinBondUpdater` samples bond lifetimes from it (mirrors `bridge/integrin_bonds.py`).

## E-cadherin sliding-rebinding FIT PARAMETERS — Rakshit 2012 SI Table S1 (extracted 2026-06-03)

From `references/1208349109_sapp.pdf` (SI p.2 model, p.4 Table S1) + `references/1471.pdf`
(Lou & Zhu 2007 source model). Faithful (PI-chosen option B). [via PubMed; DOIs
10.1073/pnas.1208349109, 10.1529/biophysj.106.097048.]

**Model (3-state, analytically solvable; W2A = X-dimer catch-slip):** states P11 (two pairs
bound), P10/P01 (one pair dissociated, pre-slide), P10' (new interaction after sliding), P00
(dissociated, absorbing). Survival S(t)=P11+P10+P10'.
```
dP11/dt  = 2 k+1 P10 + k+2 P10' − k-2 P11
dP10/dt  = k-2 P11 − 2(k+1 + k-1) P10
dP10'/dt = 2 Pn k+1 P10 − (k+2 + k-1) P10'
k-1(f) = k-1^0 · exp(−f·x / kB T)          # single-pair off-rate
k-2(f) = 2 · k-1(f/2)                        # two pairs share force
Pn(f)  = 0                       (f<0)
       = {0.5[1+sin(πf/f0 − π/2)]}^n  = sin^{2n}(πf/2f0)   (0≤f≤f0)
       = 1                       (f>f0)
```
Mean bond lifetime τ(f) = ∫S dt = −1ᵀ M⁻¹ [1,0,0]ᵀ (M = transient generator); effective
single-bond off-rate k_off(f) = 1/τ(f) — what the runtime CadherinBondUpdater samples.

**Table S1 (W2A E-cadherin X-dimer fit):**
| k-1^0 (s⁻¹) | x (nm) | k+1^0 (s⁻¹) | k+2^0 (s⁻¹) | f0 (pN) | n |
|---|---|---|---|---|---|
| 30.4 | 0.34 | 5.3 | 1985.9 | **29.2** | 4.8 |

Catch peak ≈ f0 = 29.2 pN (↔ measured ~30 pN). These are MEASURED-FIT constants (acceptance
values, NOT tuning knobs). The catch regime (lifetime ↑ to ~29 pN) is exactly the tension
proliferation generates ⇒ should resist the L2.4 fragmentation the static Morse could not.

# H.7 #7 — Hosseini-2020 γ=0.27 datum: independent primary-literature cross-check

**Date:** 2026-06-07
**Task:** Handoff §4 #7 — cross-check the MCF7 cortical-tension datum γ ≈ 0.27 mN/m
(flagged PROVISIONAL: "read off Fig 2g, γ-protocol unconfirmed").
**Method:** independent verification via PubMed metadata + PMC full-text read (not the
prior session's figure read). Companion to
[`H7_MCF7_CORTICAL_TENSION_DATUM_2026-06-07.md`](H7_MCF7_CORTICAL_TENSION_DATUM_2026-06-07.md)
(the original datum-search doc) — this doc VERIFIES it against the primary source and
states what is now resolved vs. intrinsically figure-bound.

> Source: Hosseini K, Taubenberger A, Werner C, Fischer-Friedrich E (2020).
> "EMT-Induced Cell-Mechanical Changes Enhance Mitotic Rounding Strength."
> *Adv Sci* 7(19):2001276. [DOI 10.1002/advs.202001276](https://doi.org/10.1002/advs.202001276)
> · PMID 33042748 · PMC7539203. (Verified via PubMed + PMC full text, 2026-06-07.)

---

## What the cross-check RESOLVED

1. **Paper identity — CONFIRMED.** PubMed metadata matches the cited record exactly
   (authors, journal, vol 7(19):2001276, 2020, DOI, PMID, PMC). The senior author is
   **Elisabeth Fischer-Friedrich** — the originator of the AFM parallel-plate
   confinement cortical-tension method — so the datum comes from the method's own lab,
   not a second-hand application.

2. **γ-protocol — CONFIRMED (flag resolved).** The PMC full text states the cortex
   mechanics are measured by "a previously established cell confinement setup based on
   AFM… initially round cells are dynamically confined between two parallel plates
   using a wedged cantilever." Cantilevers: tipless, 200–350 µm, force constants
   0.2–0.8 N/m, thermal-noise calibrated. Cortical tension is defined as **"the
   time-averaged value of [the effective cortical surface tension] during an oscillation
   period"** — i.e. a Laplace cortical surface tension γ, **the same observable the H.7
   model computes** (γ = F/[A(1/R₁+1/R₂)]), not a stiffness/modulus or a Pa-stress. The
   "protocol unconfirmed" flag is therefore resolved.

3. **Cell state — CONFIRMED suspended/rounded (NOT spread-adherent).** The text:
   "interphase cells are sampled in a state of suspension which also enhances
   comparability to measurements of rounded mitotic cells." So both the interphase
   (0.27) and mitotic (~0.75) MCF-7 values are **rounded/suspended** — the structural
   reason a spread-adherent MCF7 γ cannot be measured by this assay (no free spherical
   surface for Laplace's law). This re-confirms the datum doc's "spread-adherent MCF7 γ
   is an open experimental gap" and the GATE-B finding that the right observable for a
   spread adherent MCF7 may be **traction**, not cortical γ.

## What remains INTRINSICALLY figure-bound (not a fixable flag)

4. **The 0.27 central value is reported GRAPHICALLY (Fig 2g), not as an inline numeric.**
   The full text reports cortical tension only via the figure (and SI plots): there is no
   main-text "γ = 0.27 mN/m" sentence to quote. So the central value is necessarily a
   **figure read** — the prior session's read (median ≈ 0.27 mN/m, IQR ≈ 0.18–0.40,
   whiskers ≈ 0.05–0.62, from the Fig 2g top panel, y-axis 0–0.8 mN/m) stands as the best
   available value, with ~±0.03 read precision. This is a property of how the paper
   reports the data, not a gap in our verification. The "read off a figure" caveat should
   stay attached to the value, but it is no longer a *verification* gap — only a
   precision note.

## Net status of the datum

| Aspect | Before (handoff flag) | After this cross-check |
|---|---|---|
| Paper identity / DOI / PMID | assumed | **CONFIRMED** (PubMed) |
| Measurement protocol | "unconfirmed" | **CONFIRMED**: AFM parallel-plate confinement, Laplace γ (Fischer-Friedrich, senior author) |
| Observable matches model γ | assumed | **CONFIRMED** (time-averaged effective cortical surface tension) |
| Cell state | implied rounded | **CONFIRMED** suspended/rounded interphase (+ mitotic ~0.75) |
| Central value 0.27 mN/m | "read off Fig 2g" | **figure-read (intrinsic)**; median 0.27, IQR 0.18–0.40 — best available, ±0.03 |
| Spread-adherent MCF7 γ | open | **CONFIRMED non-existent** (assay needs a spherical free surface) |

**Conclusion:** the γ ≈ 0.27 mN/m anchor is **verification-strengthened**, not
overturned. It is a *rounded/suspended* MCF-7 Laplace cortical tension from the method's
own lab. It is the correct *family* of anchor for a steady-state cortical-tension gate
(vs. the Moazzeni ~10 mN/m short-pulse prestressed read), but it is a **rounded-cell**
reference — the H.7 spread-adherent operating point has no literature γ, by construction
of the assay. Keep the value tagged "rounded; figure-read ±0.03"; do **not** treat the
[0.35, 0.65] band as a spread-adherent MCF7 requirement.

## Citation-integrity

Verified against the live PubMed/PMC record (not memory). Not on the project's
hallucination list. No SourceEvidence/KnowledgeClaim auto-registered (PI-gated). The
SE-registration candidate stands and is now verification-strengthened:

- `Hosseini2020_AdvSci` → KC: "MCF-7 cortical tension ≈ 0.27 mN/m (interphase,
  suspended/rounded; AFM parallel-plate confinement; figure-read, IQR 0.18–0.40)" —
  proposed KU-3.5 *rounded-cell* anchor; protocol + identity verified 2026-06-07.

*Per PubMed terms: data retrieved from PubMed; primary source*
*[DOI 10.1002/advs.202001276](https://doi.org/10.1002/advs.202001276).*

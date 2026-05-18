# Questions for Sungwook — Phase 1 Unit 2.1 (Worker B)

Two open issues from Unit 2.1 implementation that need PI guidance before
Unit 2.2. Both are CLAUDE.md "gate-contract" surface-ups, not code bugs —
the bridge module itself is faithful to KU-2.4 / KU-2.5 / KU-2.8 / KU-2.18.

---

## Q1 — KU-2.5 catch-peak F\* with KU-2.18 defaults gives 7 pN, not 30 pN

**Issue.** With the KU-2.18 illustrative Pereverzev parameters

    k_off^slip  = 0.5 s⁻¹,  F_s = 30 pN
    k_off^catch = 0.4 s⁻¹,  F_c =  7 pN

the closed-form lifetime peak

    F* = F_s F_c / (F_s + F_c) · ln(k_c F_s / (k_s F_c))
       = 5.68 pN · ln(3.43) = 6.99 pN.

KU-2.5 separately quotes the *experimental* F\* ≈ 30 pN for α5β1–fibronectin
(Kong 2009 Nature, Elosegui-Artola 2016 Nat Mater). The two are **not
simultaneously satisfiable** with the KU-2.18 parameter set — recovering
F\* ≈ 30 pN with the same Pereverzev form needs either a much larger
`k_c` (≈ 20 s⁻¹) or substantially different F_s/F_c ratios.

**Brief acceptance.** "peak τ_max at F\* ≈ 30 pN within ±20 %." Cannot be
met with KU-2.18 defaults.

**Proposed resolution.** Keep KU-2.18 defaults; let the Unit 2.1 Sanity
Gate validate F\* against the *analytic* value from the chosen parameters
(rel err < 1e-3 ✓), and log the experimental-vs-analytic gap as info.
Refit the Pereverzev parameters to Kong 2009 force-clamp data in Phase 2.

**Decision needed.** Confirm OR (a) update KU-2.18 to a Kong 2009-fit
parameter set, or (b) keep illustrative defaults but rephrase the brief
acceptance to "within ±20 % of the closed-form F\* from the chosen
parameters".

---

## Q2 — KU-2.12 per-clutch force 5–20 pN band requires Unit 2.2 maturation

**Issue.** The brief's Task 7 asks for ≥ 90 % of per-clutch forces in
the KU-2.12 5–20 pN band on E = 5 kPa substrate. The simulation gives
mean 3.6 pN, median 2.4 pN, **27 % in band** → Sanity Gate FAILs.

KU-2.12's 5–20 pN range is measured on **mature FAs** (Plotnikov 2012
Cell, Trichet 2012 PNAS), where the integrin spring is reinforced by
vinculin binding to force-unfolded talin domains:

    k_int^eff = k_int^bare (1 + α N_vin)    (KU-2.7, Han 2021 eLife)

The Phase 1 Unit 2.1 brief explicitly excludes FA growth, talin, and
vinculin ("Unit 2.2의 일"). With bare `k_int = 1 pN/nm` and the
force-balance constraint at stall, per-clutch force is naturally
2–4 pN — at the lower end of the catch peak, *below* the literature
mature-FA window.

**Brief acceptance.** "90 % of values in 5–20 pN range." Inconsistent
with the brief's own scope exclusion.

**Proposed resolution.** Either

  (i) move the KU-2.12 acceptance into the Unit 2.2 contract (where
      vinculin enters and `k_int^eff` rises), or
  (ii) widen the Phase 1 Unit 2.1 acceptance to the nascent-FA range
       (1–8 pN), with the mature 5–20 pN band logged as info.

**Decision needed.** Pick (i), (ii), or alternative.

---

## Status note

The Unit 2.1 bridge module is **complete and correct** per KU-2.4 / 2.5 /
2.8 / 2.18; biphasic peak validation passes against the
Bangasser/Alonso-Matilla matched-stiffness analytic. The two gate FAILs
above are contract-scope issues, not implementation defects, and need
PI direction before the Unit 2.2 brief is finalised.

# Questions for Sungwook — Phase 1 Unit 2.1 (Worker B)

Status after PI conditional-acceptance review (2026-05-18):

- **Q1 (Pereverzev refit)** — still open, awaiting decision below.
- **Q2 (KU-2.12 scope)** — **RESOLVED** by PI: nascent 1–8 pN band for Unit 2.1, mature 5–20 pN band inherited by Unit 2.2.
- **Q3 (biphasic claim)** — **RESOLVED** by PI: rephrased as saturating, evidenced by replicated-seed validation (n=5 seeds, argmax spread 2.27 decades, prominence 0.34 %).
- **Q4 (`contact_radius` vs FA area)** — **RESOLVED** by PI: kept as independent model parameters, distinction documented in `substrate_stub.py` and YAML.

---

## Q1 (still open) — KU-2.5 catch-peak F\* with KU-2.18 defaults gives 7 pN, not 30 pN

**Issue.** With the KU-2.18 illustrative Pereverzev parameters

    k_off^slip  = 0.5 s⁻¹,  F_s = 30 pN
    k_off^catch = 0.4 s⁻¹,  F_c =  7 pN

the closed-form lifetime peak

    F* = F_s F_c / (F_s + F_c) · ln(k_c F_s / (k_s F_c))
       = 5.68 pN · ln(3.43) = 6.99 pN.

KU-2.5 separately quotes the *experimental* F\* ≈ 30 pN for α5β1–fibronectin
(Kong 2009 Nature, Elosegui-Artola 2016 Nat Mater). These are not
simultaneously satisfiable with the KU-2.18 parameter set — recovering
F\* ≈ 30 pN with the same Pereverzev form needs either `k_c ≈ 20 s⁻¹` or
substantially different F_s/F_c ratios.

**Phase 2 consequence.** Until the parameters are refit, the bridge
module cannot validate against the canonical Kong 2009 catch peak. The
saturating biphasic shape (Q3) is also a direct consequence of the
illustrative defaults: with F* ≈ 7 pN and N_m·F_stall/N_eng ≈ 4 pN, no
clutch reaches the slip regime, so the system never gets the
high-stiffness force dump that would produce a true peak.

**Proposed resolution.** Three options:

1. **Refit immediately in Unit 2.2.** Add a Kong 2009-fit Pereverzev
   parameter set to KU-2.18, alongside the illustrative one. The Unit
   2.2 brief explicitly mentions "Refit Pereverzev parameters" as a
   Phase 2 candidate but does not commit; promote it to first-class.

2. **Defer to Phase 2.** Keep KU-2.18 illustrative through Phase 1;
   refit in Phase 2 once Unit 2.2 has confirmed the vinculin
   reinforcement mechanism produces a true biphasic peak at the
   refitted parameters.

3. **Update KU-2.18 in-place now.** Change the canonical defaults to
   the Kong 2009 fit (estimated: `k_c ≈ 20 s⁻¹, F_s ≈ 25 pN, F_c ≈
   3 pN`, TBD). Re-run all Unit 2.1 validation; the saturating verdict
   may flip to "peaked" if F* moves into the per-clutch operating range.

**Recommendation.** Option 2 (defer to Phase 2) keeps Unit 2.1 honest
and Unit 2.2's contract clear: vinculin reinforcement is the mechanism
that should produce the peak, and a refit *plus* reinforcement is what
should reproduce Kong 2009 / Bangasser 2017 force-stiffness curves.

**Decision needed before Unit 2.2 brief is finalised.**

# KB-Registration Candidates — compartment missing-datum lit-ingest (2026-07-07)

**Origin:** PI approved "C" (build unblocked units + resolve blocks in parallel). This is the `compartment-lit-ingest`
workflow (12 agents, PubMed/bioRxiv/Consensus + BM25 corpus) researching the missing constants that block
compartment units. **All values PI-gated before registration** (SourceEvidence→KnowledgeClaim in Notion; confirm
source_audit.verdict=OK per the citation-integrity hard rule). No number invented; ABSENT verdicts are real HALTs.

## ⭐ Plan-changing findings (surfaced to PI)
1. **Myosin is LINEAR, not Hill.** No non-muscle-myosin Hill hyperbola exists; a/F₀≈0.25 is muscle-only (Hill 1938).
   The project's existing LINEAR force-velocity v=v₀(1−F/Fs), Fs≈4pN (Freedman/Cytosim/Tam) is the SOURCED law —
   so the CLAUDE.md "linear-stall = Wrong" entry is muscle-based and does NOT apply to non-muscle myosin IIA. The
   Myosin-Hill unit should be a Myosin-LINEAR unit. **PI: confirm.**
2. **E_nuc = 399 Pa (MCF7-specific, in-situ) = Fischer/Hayn/Mierke 2020** (10.3389/fcell.2020.00393) — NOT "Kim2020"
   (misattribution). R_nuc/R_cell = 0.68–0.77 (Moore 2016). Physiological-baseline rule favors in-situ 399 Pa over
   the isolated-nucleus 1–5 kPa the model currently sits at (audit#18/19). **PI: ratify in-situ 399 + R_nuc≈0.7.**
3. **MCF7 has NO caveolae** (Cav-1-deficient; Lavie 1998 / Fiucci 2002). The caveolae unit is biologically wrong
   for MCF7 → drop it; source RESERVOIR_STRAIN from folds/microvilli/ruffles instead (the magic 0.60 must go).

# KB-Registration Candidate Report — FF/MCF7 Literature Research Batch

Prepared for PI sign-off. All values below come only from the supplied research results; no new numbers were introduced. Low-confidence rows are flagged inline.

## (1) FOUND — registrable SourceEvidence → KnowledgeClaim candidates

| Datum | Value + unit | Primary source + DOI | MCF7 / proxy / generic | Conf. |
|---|---|---|---|---|
| NM IIA duty ratio (Myosin-Hill unit) | ~0.1 (NM IIA low end); non-muscle band 0.1–0.3 | Kovács, Wang, Hu, Zhang & Sellers 2003, JBC 278(40):38132 — 10.1074/jbc.M305453200 | Generic molecular constant (purified human NM IIA; not cell-line-specific) | medium |
| MCF7 nucleus in-situ Young modulus E_nuc | 399.01 ± 117.16 Pa (in-situ, adherent, n=55) | Fischer, Hayn & Mierke 2020, Front Cell Dev Biol — 10.3389/fcell.2020.00393 | **MCF7-specific** | high |
| Nucleus/cell radius ratio R_nuc/R_cell | 0.77 (Moore2016) / 0.68 (larger-n) → band **0.68–0.77** | Moore, Strohm & Kolios 2016, Int J Thermophys 37:118 — 10.1007/s10765-016-2129-y; corrob. Moore 2019 — 10.1117/1.JBO.24.10.106502 | **MCF7-specific** | high |
| Keratin single-IF force-extension (IF-cage bond, epithelial-appropriate) | Digitizable strain 0→~2.5, force 0→~1.5 nN | Lorenz, Forsting … Köster 2019, PRL 123:188102 — 10.1103/PhysRevLett.123.188102 | Proxy (recombinant single keratin; molecular property, legitimate for MCF7 epithelial IF) | high |
| Vimentin single-IF force-extension (model-oracle, two-state unfolding) | Strain 0→~2.5, plateau ~0.1–0.4 nN | Block … Janshoff & Köster 2017, PRL 118:048101 — 10.1103/PhysRevLett.118.048101 | Proxy (highest-res digitizable reference) | high |
| IF fold-stretch rupture bound | 2.6-fold avg (max 3.6-fold) stretch | Kreplak, Bär, Leterrier, Herrmann & Aebi 2005, JMB 354:569 — (verify DOI at registration) | Proxy | high |
| Spectrin membrane-skeleton 2D shear modulus μ₀ (Spectrin unit) | 6.3 µN/m (range 5.3–6.6) | Dao, Li & Suresh 2006, Mater Sci Eng C — 10.1016/j.msec.2005.08.020; corrob. Dao, Lim & Suresh 2003 — 10.1016/j.jmps.2003.09.019 | Proxy — **RBC/erythrocyte**, upper-bound-uncertain for non-erythroid | high (value) — but see HALT note |
| Nucleolus interfacial tension σ_no | ~1e-6 N/m (~1 µN/m) | Caragine, Haley & Zidovska 2018, PRL 121:148101 — 10.1103/PhysRevLett.121.148101 | Proxy (HeLa in vivo) | medium |
| Nucleoplasm viscosity η (governs nucleolar coalescence) | ~1e3 Pa·s (surrounding nucleoplasm); droplet interior ~1–10 Pa·s (Feric) | Caragine 2018 (as above); Feric et al. 2016, Cell 165:1686 — 10.1016/j.cell.2016.04.047 | Proxy (HeLa / Xenopus) | medium |
| Nesprin-SUN LINC resting tension (LINC unit) | ~8 pN (band 2–10) | Déjardin, Borghi et al. 2020, JCB — 10.1083/jcb.201908036; corrob. Arsenovic et al. 2016, Biophys J — 10.1016/j.bpj.2015.11.014 (relative FRET only) | Proxy (NIH3T3 fibroblast + MDCK epithelial) | medium |
| Glycocalyx brush thickness L_g (Glycocalyx unit) | ~30–200 nm (Muc1 Rg 32 nm; native ectodomain up to 200 nm) | Paszek et al. 2014, Nature — 10.1038/nature13535; grafting model Shurer et al. 2019, Cell — 10.1016/j.cell.2019.04.017 | Proxy (MCF10A mimetics; **not MCF7-native**) | **low — FLAG** |
| Fascin crosslink off-rate k_off (Filopodial-fascin unit, partial) | 0.12 s⁻¹ (in-vivo FRAP) to ~9 s⁻¹ (in-vitro); register band with 0.12 anchor | Aratyn, Schaus, Taylor & Borisy 2007, MBoC — 10.1091/mbc.e07-04-0346 | Proxy (B16 melanoma / in-vitro rabbit actin) | medium |

**Registration notes / flags**
- **Low-confidence FLAG — Glycocalyx L_g:** value is a proxy from engineered mammary-epithelial mimetics; register the *range* only, tagged non-MCF7. The companion grafting spacing *s* for native MCF7 is ABSENT (see §2).
- **Provenance tags mandatory** on every proxy row (Spectrin=RBC; Nucleolus=HeLa/Xenopus; LINC=fibroblast/epithelial; Keratin/Vimentin=recombinant single-filament; Fascin off-rate=melanoma/in-vitro). None may be labeled MCF7.
- **E_nuc citation correction (carry into KB):** the in-situ 399 Pa is **Fischer, Hayn & Mierke 2020**, *not* "Kim2020" — a task-lead misattribution. Adopt R_nuc≈0.7 as the unambiguous default fix; **E_nuc adoption stays PI-gated** (in-situ apparent whole-cell-over-nucleus 399 Pa vs isolated-nucleus material 1–5 kPa where the model currently sits). Do not silently swap 4700→399.
- Confirm `source_audit.verdict = OK` for each DOI before registering (citation-integrity hard rule). Kovács 2003, Fischer 2020, and the glycocalyx/LINC/nucleolus/spectrin primaries are **not yet in the BM25 corpus** → new SourceEvidence rows needed.

## (2) ABSENT / HALT — datums that gate a unit

| Datum (missing) | Unit / magic-number it blocks | Why absent | Recommended action |
|---|---|---|---|
| NMIIA Hill force-velocity curvature a/F_stall | Myosin-Hill force-velocity law (Phase 1) | Misattributed to Kovács 2003 (a solution-kinetics paper, no F-V curve); no non-muscle Hill hyperbola exists; a/F₀≈0.25 is muscle-only (Hill 1938). Project lineage is **linear**, not hyperbolic. | **Do not register a Hill curvature.** Register the sourced LINEAR law v=v₀(1−F/Fs), Fs≈4 pN (Freedman 2017 AFINES / Cytosim / Tam 2021 eq A.8). See §3. |
| Fascin single-molecule stiffness link_k | Filopodial-fascin crosslink bond | No AFM/single-molecule force-extension for fascin exists (only network/bundle rheology). Do NOT inherit filamin 8.2e5. | **HALT → PI.** Either an analytic/structural estimate flagged non-literature (magic-number block applies) or run without a sourced link_k. Off-rate half (§1) can proceed as proxy. |
| Septin cortical areal density (filaments·µm⁻²) | Septin unit + σ_c derivation | Genuinely absent everywhere (all data are yeast/Drosophila/HeLa). σ_c becomes fitted, not derivable → violates no-magic-number rule. | **HALT → PI.** Register c* curvature proxy (positive µm-scale, optimum radius ~1–3 µm; Bridges 2016 10.1083/jcb.201512029) as a named proxy at most; keep Septin unit gated OFF pending a mammalian density measurement or explicit PI proxy decision. Do not back-solve density from a target σ_c. |
| Septin bending rigidity κ_sep / L_p | Septin unit | Literature only qualitative ("highly flexible", composition-dependent, yeast); no clean MCF7-relevant digit extracted. | Do not hard-register a number; part of the Septin HALT above. |
| Caveolar area fraction, σ_flat, flattening kinetics | Membrane Reservoir/caveolae unit + RESERVOIR_STRAIN magic number | **Stronger than unmeasured:** native MCF-7 are Cav-1-deficient and devoid of caveolae (Lavie 1998; Fiucci 2002 10.1038/sj.onc.1205294). Caveolar area fraction ≈ 0. | **HALT → PI.** Either (a) set caveolar reservoir ≈0 for MCF7 and source RESERVOIR_STRAIN from a Cav1-negative-appropriate mechanism (membrane folds/microvilli/ruffles — itself a new lookup), or (b) drop the caveolae unit for the MCF7 target. |
| Glycocalyx native grafting spacing s (MCF7) | ADG-brush prefactor in Glycocalyx unit | Paszek model is a repulsive potential (no s); Shurer s is MCF10A/T47D/ZR-75-1, not MCF7. The load-bearing s is absent. | Register Shurer transition ~250 mucins/µm² → s≈63 nm as an explicitly-labeled PROXY band if the unit must proceed; surface that no MCF7-native s exists. Do not back-solve s from an assumed L_g. |
| LINC molecular stiffness k_linc (N/m) | LINC unit stiffness term | Nesprin-2-giant (~56 spectrin repeats) is strongly nonlinear (sequential unfolding ~25–35 pN plateaus); no single Hookean k for any cell line. | **HALT → PI.** Mark k_linc ABSENT. Proxy path = spectrin-repeat WLC nonlinear force-extension, NOT a linear k_linc. The yaml k_linc=1e-2 is a numerical control, not a physical constant. |
| Non-erythroid / MCF7 spectrin 2D shear modulus | Spectrin unit (μ₀ adoption) | No non-erythroid or breast-cancer 2D spectrin μ₀ exists; the only value (6.3 µN/m, §1) is RBC. | Register RBC μ₀=6.3 µN/m as **PI-gated proxy** with erythrocyte provenance tag (non-erythroid tetramers longer/less-dense, actin-integrated → upper-bound-uncertain). Do not treat as an MCF7 measurement. |

## (3) Focused sub-note — Myosin-Hill reconciliation (blocks Phase 1)

Two independent motor parameters were conflated in the leads; they do **not** bridge each other.

**A. Duty ratio — REGISTRABLE (candidate).**
NM IIA is a low-duty motor, ~0.1 (non-muscle umbrella band 0.1–0.3, IIA at the low end, IIB high ~0.2–0.4). Mechanistic basis (Kovács 2003): ADP-release rate ≫ steady-state ATPase ⇒ IIA spends only a small fraction of the cycle strongly actin-bound. This is a **molecular constant, not cell-line-specific** — MCF7 expresses NM IIA (MYH9, dominant) so the generic value is the correct physiological setpoint, not a back-solve. **Caveat to record:** tropomyosin raises NM II duty ratio 3–5× (Pathan-Chhatbar 2017, 10.1074/jbc.M117.796649) — the value is context-modulated. Confidence medium.

**B. Per-head stall force — CONFLICT, PI-gated, do NOT resolve here.**
- Single-molecule NM IIA unitary force: ~0.5–1 pN (3-bead assay).
- KB-3.18 currently carries 2 pN; AFINES/Cytosim canonical Fs≈4 pN (network/ensemble, stiffness × working-stroke basis).
- The ~4× gap is a **measurement-basis difference** (single-molecule unitary force vs stiffness×stroke), not an error the duty ratio can reconcile. NM IIA being weak/low-duty is *consistent* with the smaller single-head value.

**C. Force-velocity law.**
No non-muscle Hill hyperbola exists (the a/F_stall search returned ABSENT — see §2). The project's own lineage is **linear** v=v₀(1−F/Fs). Recommendation: for FF NMIIA, register the sourced **linear** relation with Fs≈4 pN (Freedman 2017 AFINES / Cytosim / Tam 2021 eq A.8) and load-free velocity v₀; do not invent a Hill curvature.

**HALT item for PI:** which per-head stall to adopt — IIA-specific 0.5–1 pN (single-molecule) vs generic 2 pN (KB-3.18) vs ensemble 4 pN (AFINES) — is a PI decision. The **duty-ratio datum is a clean registrable candidate independent of this**; only the stall-force sub-question is gated.

---
**Bottom line for sign-off:** 12 FOUND rows ready to register (2 MCF7-specific: E_nuc, R_nuc/R_cell; the rest proxies requiring provenance tags — glycocalyx L_g flagged low-confidence). 7 ABSENT/HALT items gate the Fascin-stiffness, Septin, Caveolae/Reservoir, LINC-stiffness, Spectrin-adoption, and Myosin-Hill-stall sub-questions — each with a recommended proxy-with-tag or gate-OFF action, no back-solving.
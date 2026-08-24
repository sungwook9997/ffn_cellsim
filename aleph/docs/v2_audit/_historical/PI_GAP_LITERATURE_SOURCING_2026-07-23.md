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

# PI-GAP Literature Sourcing — 2026-07-23

**Scope:** literature sourcing (KB `tag_query.py` + PubMed/bioRxiv) for the open PI-GAPs blocking the
`ac/engine` evidence ladder (per `AC_ENGINE_COMPLETION_ROADMAP_2026-07-23.md`). **No code changed, no
constants edited.** Every value carries a real citation; anything unsourced is marked **NOT FOUND — surface
to PI**. Nothing was invented. Values here are *candidates for PI ratification*, not auto-adopted.

Provenance: primary sources retrieved via **PubMed** (DOIs inline); KB rows cross-checked against
`outputs/tag_kb/kb.duckdb` + `docs/`. Where a DOI could not be independently re-verified this pass, it is
flagged "verify before KB-registering."

---

## 0. Executive triage

| Gap | Status | One-line |
|---|---|---|
| ⭐ resting bound-myosin (duty + per-head force) | **PARTIAL** | duty ratio sourced by isoform; the "0.5 vs 2 pN" conflict is a units artifact; a *direct single-molecule NM2 stall force* is NOT FOUND → PI must pick f_head + isoform |
| MT dynamic-instability rates | **SOURCED (proxy)** | full in-vivo mammalian-epithelial set (Rusan 2001); MCF7-specific NOT FOUND → ratify proxy |
| IF WLC EA / x_max | **PARTIAL** | vimentin well-sourced (+anisotropy caveat); keratin (primary MCF7 IF) precise EA/x_max NOT FOUND (paywalled) |
| filopodium k_fascin | **NOT FOUND → surface** | no single-molecule fascin stiffness exists anywhere; bundle-level derivation path is sourced |
| MCF7 ERM areal density | **SOURCED (estimate)** | ρ₀≈100 µm⁻² (range 100–1000) modeling baseline; MCF7 direct value NOT FOUND |
| nucleus I0-B2 | **PARTIAL** | E_nuc / two-regime knee / lamina modulus / LINC 8 pN tension sourced; k_linc *stiffness*, LINC density, clean rupture threshold, η_nuc DOIs, MCF7 aspect all NOT FOUND / verify |

---

## 1. ⭐ Resting bound-myosin fraction + per-head force  (GATE-A master blocker)

### (a) Duty ratio — SOURCED by isoform (no MCF7 value)
- **NM2A** low unloaded duty ~**0.05–0.1** — Kovács, Wang, Hu, Zhang, Sellers 2003, *JBC* 278(40):38132-40,
  [10.1074/jbc.M305453200](https://doi.org/10.1074/jbc.M305453200).
- **NM2B** higher duty ~**0.2–0.3** (slow ADP release, high actomyosin ADP affinity) — Wang, Kovács et al.
  2003, *JBC* 278(30):27439-48, [10.1074/jbc.M302510200](https://doi.org/10.1074/jbc.M302510200).
- **Duty ratio RISES under resistive load** (tension maintenance; ADP release slowed 5×/12× for 2A/2B) —
  Kovács, Thirumurugan, Knight, Sellers 2007, *PNAS* 104(24):9994-9,
  [10.1073/pnas.0701181104](https://doi.org/10.1073/pnas.0701181104).
- Filament-scale confirmation (NM2B processive, NM2A not, until load raises duty) — Melli et al. 2018,
  *eLife* 7:e32871, [10.7554/eLife.32871](https://doi.org/10.7554/eLife.32871).
- **MCF7/epithelial-specific duty: NOT FOUND** — use isoform kinetics + load-dependence.

### (b) Per-head isometric/stall force — the "0.5 vs 2 pN" conflict RESOLVED (mostly spurious)
- **"0.5" is NOT a force** — it is the *dimensionless* Hill curvature `a/F₀ = 0.5` (Kovács 2003),
  carried in `ac/motor/params_i0b3.yaml` as `kappa_hill`; a separate `F_stall_head = 0.5 pN` is an
  **AFINES placeholder** (Freedman 2017), not a Kovács single-molecule force. Likely a units conflation.
- **"2 pN"** is a **cortex-composition review figure** (Murrell, Oakes, Lenz, Gardel 2015, *NRMCB* 16(8):486-98,
  [10.1038/nrm4012](https://doi.org/10.1038/nrm4012)) = KB-3.18 lineage — not a single-molecule assay.
- **Direct single-molecule NM2A/NM2B isometric stall force (pN): NOT FOUND — surface to PI.** Closest primary
  mechanics is *step size*, not force — Norstrom, Smithback, Rock 2010, *JBC* 285(34):26326-34,
  [10.1074/jbc.M110.123851](https://doi.org/10.1074/jbc.M110.123851) (NM2B three-bead trap: +5.4/−5.9 nm steps,
  force-dependent detachment; no unitary stall force). Muscle-myosin trap history (~1–5 pN) makes a per-head
  NM2 value in **~0.5–2 pN physically plausible** but not NM2-sourced. **PI must select f_head, NM2A vs NM2B.**

### (c) Heads per minifilament + force budget
- **28–30 molecules per bipolar minifilament** is the structural literature count (Billington 2013;
  Niederman & Pollard 1975); **10** is the AFINES model simplification. Fine-grained mandate → use ≈28–30
  (N_side ≈ 14–15 per half).
- Budget (from `RESTING_SETPOINT_SOURCING_2026-07-23.md §5`): F_hoop ≈ **6,600 pN** to hold 40 Pa turgor at
  R=7.5 µm (γ≈140 pN/µm; anchored to Fischer-Friedrich 2014 HeLa-interphase 0.2 mN/m,
  [10.1038/srep06213](https://doi.org/10.1038/srep06213)). Placeholder content (N_side 10, f 0.5 pN, duty 0.10)
  is ~10–33× deficient (the recurring "~500× floor"); physiological content (≈28–30 heads, f≈2–4 pN, duty≈0.3
  rising under load) makes the *total* budget 1.2–4.8× F_hoop — whether it converts to hoop prestress is a
  native-transmission question only the full-native gbook run arbitrates.

### ⚠️ Fidelity flag (load-bearing for GATE-A, sourced not asserted)
The head–actin bond in `bell_kinetics_analytic.py` is a **pure Bell slip** (duty falls under load) — the
**opposite sign** of Kovács 2007 (resistive load *raises* duty for tension maintenance). If resting cortical
tension under-sustains, this sign is a prime suspect. Recommend PI review.

---

## 2. MT dynamic-instability rates — SOURCED (in-vivo mammalian-epithelial proxy)

**Recommended single internally-consistent dataset:** Rusan, Fagerstrom, Yvon & Wadsworth 2001,
*Mol Biol Cell* 12(4):971-980, PMID 11294900, [10.1091/mbc.12.4.971](https://doi.org/10.1091/mbc.12.4.971).
LLCPK-1α (pig-kidney epithelial), interphase, live-cell GFP-α-tubulin, **plus-end** (minus ends
centrosome-anchored, non-dynamic).

| Parameter | Value | Notes |
|---|---|---|
| v_grow | **11.50 ± 7.40 µm/min** (≈0.192 µm/s) | plus-end |
| v_shrink | **13.1 ± 8.43 µm/min** (≈0.218 µm/s) | plus-end |
| f_catastrophe | **0.026 ± 0.024 s⁻¹** | |
| f_rescue | **0.175 ± 0.104 s⁻¹** | |

- **KB reconciliation:** existing KB-DRAFT-3-07 holds `v_grow 2.0, v_shrink 17.0, f_cat 0.005, f_res 0.044`
  — these are **in-vitro tubulin** (Walker 1988, [10.1083/jcb.107.4.1437](https://doi.org/10.1083/jcb.107.4.1437))
  and the manifest already flags 2/4 mismatch Walker's raw data. In-cell growth is ~6× faster (MAPs/+TIPs).
  Per the physiological-baseline rule, **recommend PI replace the in-vitro draft with Rusan 2001**.
- Corroborating in-vivo mammalian: Vasquez 1997 (BSC-1 monkey-kidney epithelial, untreated baselines),
  [10.1091/mbc.8.6.973](https://doi.org/10.1091/mbc.8.6.973); Stepanova 2003 EB3 growth-only,
  [10.1523/JNEUROSCI.23-07-02655.2003](https://doi.org/10.1523/JNEUROSCI.23-07-02655.2003).
- **NOT FOUND:** MCF7/breast-specific DI (→ ratify LLCPK-1α proxy, as HeLa is used for cortical tension);
  in-cell minus-end rates (N/A, anchored). **Pause state:** interphase MTs spend ~73.5% in pause — the planned
  2-state grow/shrink switch omits it; PI decision whether to model a 3rd state.

---

## 3. IF WLC EA / x_max — PARTIAL (vimentin sourced; keratin precise values NOT FOUND)

**Critical caveat:** IFs are strongly anisotropic — *bending*-derived E (300–900 MPa) is ~100× the *axial
tensile* small-strain E (~7–15 MPa). **For a tensile WLC cable use the tensile numbers, not bending.** Response
is nonlinear/tri-phasic (entropic toe → α-helix→β-sheet unfolding plateau → stiffening).

Vimentin (well-sourced):
- Axial tensile tangent E ~**3–12 MPa** (low strain) → ~**185 MPa** (>150% strain, ~20× stiffening) — Qin,
  Kreplak & Buehler 2009, *PLoS ONE* 4:e7294, [10.1371/journal.pone.0007294](https://doi.org/10.1371/journal.pone.0007294).
- Bending E 300–400 MPa (native) / ≥900 MPa (fixed); diameter ~10 nm outer/3 nm inner (hollow) — Guzmán 2006,
  *JMB* 360:623, [10.1016/j.jmb.2006.05.030](https://doi.org/10.1016/j.jmb.2006.05.030).
- x_max ~**200–300% strain** (MD) / rupture ~1.0–1.5 strain (single-filament); >70% energy dissipated per
  cycle — Block 2018, *Sci Adv*, [10.1126/sciadv.aat1161](https://doi.org/10.1126/sciadv.aat1161); Forsting 2019,
  [10.1021/acs.nanolett.9b02972](https://doi.org/10.1021/acs.nanolett.9b02972). L_p ~1–2 µm (Mücke 2004;
  Nöding & Köster 2012).
- Desmin (type-III analog) x_max **3.4× (≈240% strain)**, tensile strength ≥240 MPa — Kreplak, Herrmann & Aebi
  2008, *Biophys J* 94:2790, [10.1529/biophysj.107.119826](https://doi.org/10.1529/biophysj.107.119826).

EA conversion (hollow A = 71.5 nm² = 7.15e-17 m²): toe EA ≈ **0.5–0.9 nN**, post-unfolding ≈ **13 nN**
(bending-derived 64 nN = anisotropy artifact, do not use for axial). **A (hollow vs solid) is a PI modeling choice.**

- **NOT FOUND (keratin — the primary MCF7-epithelial IF):** precise K8/K18 single-filament EA/Young's modulus
  and exact x_max — the definitive OT papers (Lorenz 2019, *PRL* 123:188102,
  [10.1103/PhysRevLett.123.188102](https://doi.org/10.1103/PhysRevLett.123.188102); Lorenz 2023, *Matter*
  6:2019, [10.1016/j.matt.2023.04.014](https://doi.org/10.1016/j.matt.2023.04.014)) are paywalled — only
  qualitative ">2.5-fold extension, 133% in-cell without rupture." **Must read full text/SI or surface.**
- Also NOT independently DOI-verified this pass: Kreplak 2005, Mücke 2004, Nöding & Köster 2012, Block 2017.
- EMT relevance: keratin elongates but retains stiffness (viscous sliding); vimentin softens but retains
  length (α-helix unfolding) — Lorenz 2019/2023.

---

## 4. filopodium k_fascin — NOT FOUND → surface to PI (bundle-level derivation is sourced)

- **No published single-molecule fascin crosslink stiffness (pN/nm) exists.** The KB already flags this
  (KB-DRAFT-3-23). **Do NOT inherit filamin 8.2e5.** Honest path = derive from bundle bending stiffness.
- Sourced bundle-level anchors:
  - Bundle L_p rises 10 µm (bare F-actin) → **150 µm** at fascin:actin 1:2; κ_bundle = L_p·k_BT — Takatsuki,
    Bengtsson & Månsson 2014, *BBA Gen Subj*, [10.1016/j.bbagen.2014.01.012](https://doi.org/10.1016/j.bbagen.2014.01.012).
  - Geometry: inter-filament ~8 nm, repeat ~36–38 nm, strictly parallel — Courson & Rock 2010, *JBC*,
    [10.1074/jbc.M110.123117](https://doi.org/10.1074/jbc.M110.123117).
  - Slip-like/avidity kinetics (no numeric k_off0) — Courson & Rock 2010; Maier/Lieleg/Bausch 2015, *EPJE*,
    [10.1140/epje/i2015-15050-3](https://doi.org/10.1140/epje/i2015-15050-3).
  - Euler buckling ~7.7 pN N²-tight bundle context — KB-3.8 (Mogilner & Rubinstein 2005; Pronk 2008).
- **Derivation caveat:** the coupled-bundle model papers needed (Claessens 2006 *Nat Mater*; Bathe 2008
  *Biophys J*; Vignjevic 2006 *JCB*) are **absent from the corpus** — ingest before the derivation is defensible.
- **NOT FOUND:** numeric fascin k_off0 (s⁻¹); per-tip formin count. Surface both.

---

## 5. MCF7 ERM areal density — SOURCED as a modeling estimate (MCF7 direct value NOT FOUND)

- **Recommend ρ₀ ≈ 100 µm⁻² (10¹⁴ m⁻²) baseline, sensitivity sweep 100–1000 µm⁻²** — Alert, Casademunt,
  Brugués, Sens 2015, *Biophys J* 108:1878, [10.1016/j.bpj.2015.02.027](https://doi.org/10.1016/j.bpj.2015.02.027).
  Generic animal-cell cortex, stated as a modeling estimate (not MCF7 measurement).
- **Independent cross-check:** adhesion energy J ≈ 6e-6 J/m² (Charras 2007, *Biophys J* 94:1836,
  [10.1529/biophysj.107.113605](https://doi.org/10.1529/biophysj.107.113605)) ÷ E_bond 3–25 k_BT → ρ ≈
  55–1400 µm⁻²; cortex mesh ξ≈30 nm route agrees. **Order 10²–10³ µm⁻² robust across 3 independent methods.**
- Consistency with KB: per-link critical force σ*/ρ₀ ≈ 16–18 pN matches KB single-ERM 11.4 pN / F_unbind ~50 pN;
  per-ERM stiffness 4600 pN/µm (Braunger) already sourced.
- Active-fraction caveat (physiological-baseline rule): ρ₀ = *available* linkers; the load-bearing
  (PIP₂+Thr-phospho) fraction is smaller (Bosk et al. 2011, *Biophys J* 100:1708,
  [10.1016/j.bpj.2011.02.039](https://doi.org/10.1016/j.bpj.2011.02.039)) → effective density may sit toward
  ~10²/µm². PI to decide whether to wire ρ₀ or an active-scaled value.
- **Wiring:** do NOT hard-code a per-node integer — set N_linker/node = ρ₀ × A_node (grid-invariant). "1/node"
  is physiological only if node area ≈ 0.01 µm².
- **NOT FOUND:** direct MCF7/epithelial ERM µm⁻² immuno-count.

---

## 6. Nucleus I0-B2 — PARTIAL

| Quantity | Value | Source |
|---|---|---|
| E_nuc (bulk) | in-cell ~5 kPa, isolated ~8 kPa (origin of 1–10 kPa band) | Caille 2002, *J Biomech* 35:177, PMID 11784536, [10.1016/s0021-9290(01)00201-9](https://doi.org/10.1016/s0021-9290(01)00201-9) |
| E_nuc MCF7-specific | project uses 399 Pa (MCF7)/157.7 Pa (MDA-MB-231), whole-cell-over-nucleus AFM | Fischer, Hayn & Mierke 2020, *Front Cell Dev Biol* 8:393, [10.3389/fcell.2020.00393](https://doi.org/10.3389/fcell.2020.00393) — **exact 399/157.7 not in abstract; verify body** |
| Chromatin (small-strain) spring k | ~0.70 nN/µm (HeLa; MEF 0.60; low-lamin 0.22); **governs <3 µm ext (<~30%)** | Stephens, Banigan, Adam, Goldman & Marko 2017, *MBoC* 28:1984, PMID 28057760, [10.1091/mbc.E16-09-0653](https://doi.org/10.1091/mbc.E16-09-0653) |
| Lamin-A/C strain-stiffening + **knee ≈ 3 µm ext** | large-strain regime, 1.5–2.5× | Stephens 2017 (same); mechanism Banigan, Stephens & Marko 2017, *Biophys J* 113:1654, [10.1016/j.bpj.2017.08.034](https://doi.org/10.1016/j.bpj.2017.08.034) |
| Nuclear lamina in-plane area modulus | **25 mN/m** | Dahl, Kahn, Wilson & Discher 2004, *J Cell Sci* 117:4779, PMID 15331638, [10.1242/jcs.01357](https://doi.org/10.1242/jcs.01357) |
| **LINC / nesprin resting tension** | **~8 pN (a TENSION, not stiffness)** — Déjardin anchor CONFIRMED | Déjardin et al. 2020, *JCB* 219(10):e201908036, PMID 32790861, [10.1083/jcb.201908036](https://doi.org/10.1083/jcb.201908036); corrob. Arsenovic 2016 [10.1016/j.bpj.2015.11.014](https://doi.org/10.1016/j.bpj.2015.11.014) — **8 pN is in body/FRET calibration, not abstract; verify** |
| R_nuc (MCF7) | ~5.1–6.3 µm (KB recommends 6 µm over provisional 5) | KB-3.33 (Medium). Pin N:C as diameter (0.68) vs volume (~0.30) ratio to avoid 2× mass error |

**NOT FOUND — surface to PI:**
1. **k_linc as a STIFFNESS (pN/µm)** — no literature value; Déjardin gives a *tension* only. Zero LINC/nesprin
   claim rows in KB. Spectrin-repeat unfolding 25–35 pN (Rief 1999) is a rupture scale, not resting stiffness.
   **PI-authored value + LINC areal density (also NOT FOUND) required.**
2. **Clean NE rupture threshold** (tension/strain/curvature) — WEAK. KB-3.32 asserts ~1% projected-area strain
   (Zhang-Lele 2018, MCF-10A) but that primary DOI/PMID could NOT be confirmed this pass → KB-internal until
   resolved. Denais/Lammerding 2016 (*Science* 352:353, [10.1126/science.aad7297](https://doi.org/10.1126/science.aad7297))
   confirms the phenomenon but gives NO numeric threshold.
3. **η_nuc (nucleoplasm viscosity)** value ~52 Pa·s bulk / G′ 18 Pa is KB-carried (Tseng 2004, *JCS* 117:2159,
   [10.1242/jcs.01073](https://doi.org/10.1242/jcs.01073); de Vries 2007) but **not PubMed-confirmed this pass**
   — resolve DOIs before hard-coding.
4. **MCF7 adherent nuclear oblateness/aspect ratio** — NOT FOUND (explicit KB-3.33 gap).

---

## 7. Consolidated NOT-FOUND → PI decision list

1. **Single-molecule NM2A/NM2B isometric stall force (pN)** — pick f_head + isoform (0.5/2 pN both non-primary).
2. **NMII head–actin bond sign** — Bell-slip (current) vs load-strengthening (Kovács 2007) — GATE-A relevant.
3. **MCF7 MT DI rates** — ratify LLCPK-1α (Rusan 2001) proxy; decide 2-state vs 3-state (pause).
4. **Keratin K8/K18 single-filament EA + x_max** — read Lorenz 2019/2023 full text or surface.
5. **Fascin single-molecule k_fascin + k_off0** — none exists; approve bundle-derivation path + ingest
   Claessens 2006 / Bathe 2008 / Vignjevic 2006.
6. **MCF7 ERM areal density** — adopt ρ₀≈100 µm⁻² (Alert 2015) estimate + wire as ρ₀×A_node; decide
   available vs active fraction.
7. **k_linc stiffness + LINC areal density**, **NE rupture threshold primary source (Zhang-Lele)**,
   **η_nuc DOIs**, **MCF7 nuclear aspect** — all NOT FOUND / verify.

No code or constants changed. This is a sourcing dossier for PI ratification only.

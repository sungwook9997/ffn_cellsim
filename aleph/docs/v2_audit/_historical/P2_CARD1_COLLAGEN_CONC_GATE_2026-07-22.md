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

# P2 · Card 1 — Collagen modulus = concentration-resolved conditional gate (PROPOSAL, PI sign-off required)

**Lane:** P2 evidence/SoT · **Date:** 2026-07-22 · **Status:** PROPOSED gate-contract change — surfaced to PI,
NOT applied. No gate value edited to pass. `make kb-check` unaffected (no committed KB artifact touched).
**SoT rule:** register in Notion Contract-Graph only; never hand-edit `kb.duckdb` / the Obsidian vault.

Executes `AC_DECISION_CARDS_2026-07-22.md` Card 1. The collagen conflict is a **missing-condition defect**,
not a value disagreement: a single Pa band was pinned without stating concentration / temperature / frequency /
strain amplitude / geometry, so two "bands" that describe *different concentrations* looked contradictory.

---

## 1. Evidence assembled (KB, verified this session)

| KB claim | Content (verbatim value_range_si) | Cites |
|---|---|---|
| **KB-1.32** | Reconstituted collagen-I gel **G′(2 mg/mL, 37 °C) ≈ 15 Pa**, **G′ ∝ c^2.1**, ξ ≈ 2 µm @2 mg/mL scaling ~c^−0.5. **1 mg/mL ~1–5 Pa; 3 mg/mL ~30–100 Pa.** | Yang-Leone-Kaufman 2009 BJ 97:2051 (10.1016/j.bpj.2009.07.035); Licup 2015 PNAS; Mickel 2008 BJ |
| **KB-1.V.2.1** | G′ ~0.5–350 Pa over c=0.5–4 mg/mL, **n ≈ 2.0–2.1** (37 °C, bending-dominated). Anchors **1.0→5, 3.0→55, 7.0→342 Pa**. Absolute G′ varies 2–10× by source/T. | Licup 2015 PNAS; Yang/Kaufman 2009 BJ PMC2756380 |
| **KB-1.V.2.2 / KB-1.7** | mesh **ξ ~ c^−0.5** (3D); ξ≈2 µm @1.5 mg/mL | Yang/Kaufman 2009; Licup 2015; Wolf 2013 |
| **KB-1.V.2.5** | strain-stiffening: differential modulus K rises >10–100× above **γ_c ~0.1–0.5** (γ_c **concentration-independent**, geometry-set); K~σ. **Must EMERGE, do NOT tune.** | Storm 2005 Nature; Licup 2015; Motte 2013 |
| **KB-1.30** | generic validation benchmark: G0 ~36 Pa (ref **30–100**), G(γ) stiffening, σ(r)~r⁻¹. Accept within 2×. | Jansen 2018; Licup 2015; Storm 2005; Han 2017 |

**Key resolution:** KB-1.32 states the 30–100 Pa band is the **3 mg/mL** value and 1 mg/mL is **1–5 Pa**. The two
code "bands" are therefore two *different concentrations*, exactly as Card 1 diagnosed.

---

## 2. Proposed conditional gate split (replaces the single-band gate)

All four gates specify **T, frequency, strain amplitude, geometry** (the dropped conditions). Small-strain shear
rheology unless noted. These are **acceptance gates the emergent Mikado modulus must satisfy** — the network is
never tuned to them (hard rule).

### `VG-ECM-G0(c)` — small-strain storage/shear modulus G′, per-concentration band
- **Conditions:** T = 37 °C; ω ≈ 0.1–1 rad/s (linear plateau); strain amplitude γ ≤ 0.02 (LVE regime); **3D bulk
  cone-plate/parallel-plate rheology** (NOT AFM indentation — see AFM caveat below).
- **Per-concentration bands (source anchors, KB-1.32 / KB-1.V.2.1):**
  | c (mg/mL) | G′ band (Pa) | anchor | disposition |
  |---|---|---|---|
  | 1.0 | 1–5 | KB-1.32 / KB-1.V.2.1 anchor 5 | reference anchor |
  | **1.5** | **[figure-extract] (interp ≈ 5–15, central ~11–12)** | Yang-Kaufman fig (see §4) | **PRODUCTION GATE — blocking on figure extract** |
  | 2.0 | ~15 (± source 2–10×) | KB-1.32 | cross-check |
  | 3.0 | **30–100** | KB-1.32 / KB-1.V.2.1 anchor 55 | **3 mg/mL validation** |
  | 7.0 | ~342 | KB-1.V.2.1 | high-c cross-check |
- **AFM caveat (KB-1.14):** AFM indentation with a ~1 µm tip reads 10–100× stiffer / softer than bulk G0 (E_AFM
  0.1–10 kPa vs bulk G0 30–200 Pa). The gate is **bulk rheology G′**; do not grade a bulk band against an AFM read.

### `VG-ECM-c-scaling` — G′ ∝ c^n, n ≈ 2.0–2.1
- **Conditions:** same as G0(c); fit over c = 0.5–4 mg/mL, 37 °C.
- **Band:** n ∈ [1.9, 2.3] (KB-1.V.2.1 n≈2.0–2.1; KB-1.32 c^2.1; absolute prefactor 2–10× source-variable).
- **⚠️ Distinct from the mesh exponent** — see §4 code note (this is the *emergent modulus* exponent, not the
  code's `conc_exponent`).

### `VG-ECM-stiffening` — differential modulus increase past a critical strain (emergent)
- **Conditions:** T = 37 °C; strain sweep to γ ~1; measure K = dσ/dγ; γ_c must be **concentration-independent**.
- **Band (KB-1.V.2.5):** K rises **>10× above γ_c ~0.1–0.5**; K~σ in the nonlinear regime; γ_c ~c-independent.
  **Emergent-only** (never tuned). Already a live acceptance target.

### `VG-ECM-normal-stress` — N1 sign and magnitude (NEW; currently unspecified)
- **Conditions:** T = 37 °C; shear at finite strain; measure first normal-stress difference N1.
- **Band:** collagen/biopolymer networks show **negative N1** (they pull inward / tend to contract under shear) —
  opposite to polymer melts. Sign is the primary acceptance criterion; magnitude source-dependent.
- **Evidence gap:** no KB claim yet carries N1 for collagen. **Blocking:** register a SourceEvidence row (candidate:
  Janmey et al. 2007 Nat Mater 6:48, "Negative normal stress in semiflexible biopolymer gels") before this gate is
  gradeable. Flagged to PI as a **new-source registration**, not auto-created.

---

## 3. Disposition of the historical bands / results

- `5–100 Pa` ([ff_ecm_validate.py:37](../../scripts/ff_ecm_validate.py)) → **exploration envelope (demoted)**: it is
  the 1–3 mg/mL span, useful as a coarse smoke bound, not a production gate.
- `30–100 Pa` ([ecm_library.py:111](../../ff/ecm_library.py)) → **relabeled 3 mg/mL validation** (KB-1.32).
- **1.5 mg/mL band → production gate**, value to be **extracted directly from the Yang-Kaufman source figure** (§4),
  NOT the interpolated 11 Pa. **Do NOT pin 11 Pa as a single-number gate** (Card 1 explicit).
- **Historical 11–15 Pa results → RE-EVALUATION PENDING** (neither PASS nor FAIL) until this contract is PI-approved.

---

## 4. Code inconsistencies — exact fix specification (surfaced to PI; NOT applied)

`aleph/laws/ecm_library.py`, `collagen_I` ECMSpec (lines 109–123) + `ff_ecm_validate.py` LIT (lines 36–38):

1. **Band↔concentration mislabel (gate value — PI to ratify).** `modulus_band_Pa=(30.0, 100.0)` is carried AT
   `ref_conc=1.5` but 30–100 Pa is the **3 mg/mL** band (KB-1.32). Fix: either tie the (30,100) band to a 3 mg/mL
   validation entry, or replace the reference band with the **1.5 mg/mL figure-extracted** value. Preferred: move to
   the per-concentration `VG-ECM-G0(c)` table above rather than one scalar band pinned to one concentration.

2. **`conc_exponent=0.5` is NOT an exponent conflict — Card's "second conflict" is a FALSE ALARM (finding).**
   `conc_exponent` governs the **mesh** law `ξ ∝ c^(−0.5)` ([ecm_library.py:72](../../ff/ecm_library.py) docstring;
   confirmed KB-1.7 / KB-1.V.2.2: ξ~c^−0.5). The **modulus** exponent n≈2.0–2.1 is a *different quantity* that
   **emerges** from the microstructure (ρ_L = ξ^−2 ∝ c → bending-dominated G′ ∝ c^~2; the code's own note line 123
   and the `collagen_concentration()` harness both confirm the emergent n≈2). **They are not in conflict.**
   Recommendation: **do not change the 0.5 value**; add a one-line comment/rename hint (`mesh_conc_exponent`) so it
   is never re-confused with the modulus exponent. Surfaced so PI does not "fix" a non-bug.

3. **Stale KB reference (metadata).** `sources=(... "KB-1.30" ...)` — KB-1.30 is the *generic* IF>15 benchmark claim;
   the concentration-resolved anchors are **KB-1.32** (per-conc G′ + ξ) and **KB-1.V.2.1** (G′(c) power law). Fix:
   add `KB-1.32` (and keep `KB-1.V.2.1`, already present); demote `KB-1.30` to a benchmark cross-reference.

4. **Yang-Kaufman DOI — TWO real papers, disambiguate (citation-integrity).** The code cites
   `Yang-Kaufman 2009 BiophysJ 10.1016/j.bpj.2008.10.063`. The KB carries **two rows under one citation_key**
   `YangKaufman2009_BiophysJ`:
   - `10.1016/j.bpj.2008.10.063` → Yang & Kaufman 2009 **BJ 96:1566** "Rheology + confocal reflectance of collagen
     (+HA) self-assembly" (anchor_status OK — the DOI resolves and the paper is real; **not a hallucination**).
   - `10.1016/j.bpj.2009.07.035` → Yang, Leone & Kaufman 2009 **BJ 97:2051** "Elastic moduli of collagen gels
     predicted from 2D confocal" — the paper KB-1.32 / KB-1.V.2.1 anchor to for **G′ ∝ c^2.1** and the c-anchors.
   Both are legitimate Yang-Kaufman 2009 collagen-rheology papers. **Finding:** the G′(c) anchors belong to the
   **97:2051** paper; the code's `bpj.2008.10.063` points to the **sibling 96:1566** paper. Fix: (a) in Notion,
   split the duplicated `YangKaufman2009_BiophysJ` into two keyed rows (`YangKaufman2009a_BJ96` /
   `YangKaufmanLeone2009_BJ97`); (b) pin each concentration anchor to the exact paper+figure it is read from;
   (c) update the code citation to `10.1016/j.bpj.2009.07.035` for the G′(c) scaling.

---

## 5. Blocking step for the production 1.5 mg/mL gate

The 1.5 mg/mL production band must be **digitized from the Yang-Kaufman 97:2051 G′(c) figure** (or KB-1.V.2.1's
anchor plot), not interpolated. The paper PDF is **not currently in `aleph/references/`** — surfaced to PI:
either (a) supply the PDF for figure extraction, or (b) approve the log-log interpolation ≈ 5–15 Pa (central ~11–12,
from anchors 1→5 / 3→55 with n≈2.18) as a **provisional** 1.5 mg/mL band tagged "interpolated, figure-extract
pending". Until then, `VG-ECM-G0(c)` at 1.5 mg/mL is registered as **[BLOCKED-ON-FIGURE]**.

---

## 6. What PI is asked to ratify

1. Split the single collagen modulus gate into `VG-ECM-{G0(c), c-scaling, stiffening, normal-stress}` with the
   conditions in §2.
2. Adopt the per-concentration `VG-ECM-G0(c)` table; relabel 30–100 Pa as 3 mg/mL; demote 5–100 Pa to envelope.
3. Approve either the figure PDF supply or the provisional interpolated 1.5 mg/mL band (§5).
4. Ratify the code-fix spec §4 items 1, 3, 4 (item 2 is a non-bug — confirm no change).
5. Register the Janmey 2007 N1 SourceEvidence for `VG-ECM-normal-stress`.
6. Historical 11–15 Pa collagen results marked **re-evaluation-pending** until (1)–(4) land.

---

## APPLICATION LOG — 2026-07-23 (PI pre-authorized package application)

Applied per PI pre-authorization of `P2_RATIFICATION_PACKAGE_2026-07-22.md`. The package FINALIZED the 1.5 mg/mL
production value as **G′N = 13.14 ± 0.04 Pa @ 37 °C** (Yang-Kaufman 2009 BJ96:1566 Table 1, a DIRECT table value
read from `paper_chunks` — supersedes the earlier §5 "[BLOCKED-ON-FIGURE]" / interpolated ≈11 Pa).

**Notion Contract-Graph (SoT):**
- Updated **KB-1.32** Value Range SI to carry the Table-1 37 °C anchors (0.5→1.39, 1.0→7.60, 1.5→13.14 Pa) +
  Table-2 32 °C; c^2.1 kept on BJ97:2051; ξ~c^−0.5.
- Split the duplicated **YangKaufman2009_BiophysJ** SourceEvidence (two DOI rows under one key) into two keyed
  rows: **YangKaufman2009a_BJ96** (10.1016/j.bpj.2008.10.063 → direct G′N(c) rheology) and
  **YangKaufmanLeone2009_BJ97** (10.1016/j.bpj.2009.07.035 → G′∝c^2.1 scaling).
- Registered N1 SourceEvidence **Janmey2007_NatMater** (Nat Mater 6(1):48-51, DOI 10.1038/nmat1810, PMID
  17187066 — PubMed-verified this session; abstract confirms negative N1).
- Created 4 ValidationGates: **VG-ECM-G0(c)** (per-c band, 1.5 mg/mL = 13.14 Pa PRODUCTION GATE, accept 10–16),
  **VG-ECM-c-scaling** (n∈[1.9,2.3]), **VG-ECM-stiffening** (K>10× above γ_c~0.1–0.5), **VG-ECM-normal-stress**
  (negative N1, anchored to Janmey2007_NatMater). All Type=Validation, Status=draft.

**Code/KB (this repo):**
- `aleph/laws/ecm_library.py` `collagen_I`: `modulus_band_Pa (30,100)→(10,16)`; sources `KB-1.30→KB-1.32`;
  Yang-Kaufman DOI disambiguated (BJ96 direct + BJ97 scaling added); notes carry the Table-1 anchors. **`conc_exponent=0.5`
  UNCHANGED** (mesh exponent; the "second exponent conflict" is a FALSE ALARM — comment added to prevent a non-fix).
- `aleph/scripts/ff_ecm_validate.py` `LIT["collagen_I"]`: `band (5,100)→(10,16)`; note re-anchored to the 13.14 Pa
  BJ96 Table-1 value (5–100 = exploration envelope, 30–100 = 3 mg/mL).
- No value was changed to make a gate pass; (10,16) is the literature-sourced 13.14±0.04 Pa (a tightening toward the
  sourced datum, not a loosening). `make kb-check`: [runs] gate OK · [params] gate OK · no drift (neither file is a
  manifest-tracked constant; ff_ecm_validate is a standalone script, not part of kb-check).

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

# Session status — 2026-07-01 (FF validation sprint + DCM /goal completion)

**Branch:** dcm/main  **Scope:** the overnight window ~2026-06-30 23:00 → 2026-07-01 13:42, in which two
engines advanced **disjointly on the same branch** (no file collisions; ff/main is an ancestor of dcm/main).
Every quantitative claim below was **independently recomputed + adversarially verified** by a 7-agent
read-only workflow (npz re-computation, git audit, gate re-runs) — verdicts are stated as
CONFIRMED / PARTIAL, not taken from commit messages.

---

## A. FF engine (going-forward fine-grained; this session's Lead work) — 12 commits, verdict ✅ CONFIRMED

### A-1. "a, b, c" — lit-faithful MCF7 densities + band-reporting (`94377c9`)
- **(a)** per-motor NMII stall documented as a literature RANGE (Stam-Hocky 5 / Stachowiak ~8 / Billington
  structural ~15 pN); `f_myo` is a swept controlled variable, `NMIIA_MINIFIL_STALL_PN=5` the conservative
  anchor — not tuned.
- **(b)** `gamma_floor_production` reports the floor vs the **MCF7-active band** (Hosseini 2020 180–400 pN/µm
  × 70% active) as PRIMARY; "Salbreux" relabeled Chugh 2017 (HeLa) misnomer.
- **(c)** every cortical density → direct measured areal value at MCF7 R=7.5 µm: **actin 100/µm² (KB-3.18)**
  → N=70686; **NMIIA 0.625/µm² (Nie 2015)** → N=442, DECOUPLED from the old n_fil/10 ratio; crosslinkers 1:1.

### A-2. ⭐ Substantive finding — the clean γ_myo floor metric
- At native density the `γ_active` SUM is **passive-contaminated**: γ_active(f_myo=0) ≈ 0.84 pN/µm, plateaus
  300→8000 steps (structural, not under-relaxation), carrier = γ_actin. **Adversarially verified physical**:
  scales linearly with density (0.0085 pN/µm per fil/µm², constant over 235×); **crosslinks-OFF → exactly 0**
  (bare great-circle arcs are geodesics) → it is crosslink-network prestress/frustration, not an arc/discretization
  artifact or a measure_gamma bug.
- Floor now reported on the clean **γ_myo** channel (linear, zero-intercept, γ_myo = 0.0199·f_myo). This is a
  pre-existing structural channel (Stage 6d), not invented to change the number; the switch makes the floor
  **DEEPER (1.87×) = honest**, not shallower.
- **Headline: γ_myo ≈ 7.4e-5 mN/m → floor ~1700× vs MCF7-active** (3303× vs Salbreux-active, 4719× vs total).
  Deeper than the prior ~416× because Nie myosin is 4.8× sparser + γ_myo removes the γ_active inflation.
  **CONFIRMS the force-magnitude γ-floor; does not rescue it.**

### A-3. Validation hardening (analytic ground-truth = PRIMARY always-on; oracle = cross-check)
- **κL/2R²** bending-energy magnitude anchor (always-on, binary-independent; corrected <0.5%, raw under-counts
  exactly (n−2)/(n−1)).
- **Turgor Guo-closure + DERIVED no-magic K_vol** = Π_in0/(1−vmin_frac) — supports the no-magic-number rule.
- **GPU↔CPU bit-parity** across every kernel + integrator at native scale (reshape 2.2e-16 → loaded-shell
  4.5e-14; link_spring 2.2e-11 = float64 atomic ordering); codified `test_gpu_parity.py` (skip-if-offline).
- **Production-envelope robustness (A5000):** finite/stable across N 1k–100k, k_xl 0.1→Ferrer 4.6e5 (4.6e6×),
  f_myo 0–40; native loaded-shell V/V0=1.00001. Stiff-crosslink relaxer (crosslink_turnover) confirmed working.
- **Unified architecture:** all 5 weave() architectures relax-stable under a 30 nm perturbation.
- **Open-lever lit re-check** (PubMed): no cortical/SF NMII areal density above Nie 2015 — documented gap stands.

### A-4. Artifacts
- Docs: `FF_STAGE6Q_LITFAITHFUL_DENSITY`, `FF_STAGE6R_GPU_PARITY`, `FF_VALIDATION_STATUS`, `ENGINE.md` (→ 6r).
- Figures: `outputs/ff/figs/gamma_floor_litfaithful_6q.png`, `validation_sprint_6qr.png` (4-panel).
- Tests: **FF suite 98 passed / 4 skipped** (102 collected / 20 files). Memory ×3, Notion milestones ×2.

**Independent audit verdict:** ✅ **CONFIRMED** — no magic-number tuning; densities are measured values,
K_vol derived, the γ_active→γ_myo switch physics-justified and floor-deepening.

---

## B. DCM /loop (separate autonomous session) — `/goal` phases completed

Phase 1→6 auto-advanced (plan: `DCM_FACETING_TO_PRODUCTION_PLAN_2026-07-01.md`), HEAD `9a8435d` @ 13:42.

| Phase | result | verdict |
|---|---|---|
| 1–4 aggregation | N=400 clean faceted foam (watertight χ=2, Q149, V/V0 1.000, pen-ratio 0.963) | ✅ gates PASS |
| 5 proliferation | N=200→280 (80 divisions, cfl 0 stable); **T47D interior-deforms-more gradient** | ⚠️ **PARTIAL** |
| 6 spreading | A/A0 top-down = 1.000 = does NOT spread = **traction-limited structural limit** | ⚠️ **PARTIAL** |

### Phase 5 (T47D gradient) — real but the commit over-states
- ✅ CORE-deforms-more-than-outer is **real + statistically significant** (recomputed from PHASE5d_prolif_n200.npz:
  CORE 0.021 > MID 0.013 > RIM 0.010; Welch p=5e-7; Spearman ρ=−0.38, p=3.6e-11). **Growth flips the sign**
  (growth-OFF equilibrium foam has ρ=+0.72) → a controlled before/after contrast, not tuned.
- ⚠️ Caveats: **single seed** (seed=23, no ensemble error bar); **MID-vs-RIM not significant** (defensible claim
  is "core > outer", not a clean CORE>MID>RIM staircase); the commit's **slope −0.00072 does not reproduce**
  (no analysis code committed — numbers live only in commit msg / viewer title); T47D source measures cell
  **area (µm²)**, sim measures **asphericity** → qualitative signature match, not like-for-like; **div-rate
  raised 0.01→0.08** (ref 0.04) to drive the demo — PI glance advised.

### Phase 6 (spreading limit) — honest, one caveat dropped
- ✅ A/A0=1.000 is real (recomputed from PHASE6_spread_bare.npz; correct top-down silhouette; maxZ & V/V0 held;
  cadherin stayed engaged). **No derived parameter lowered** (cad 29.2 pN, ecm 30 pN full-strength) → not
  gate-loosening. Consistent with 3 prior lines (Layer-2, lamellipodium-graft, Phase C): cohesive spheroids
  compact, only single cells spread.
- ⚠️ Caveats: **pen_frac ~3.1 contact-confound not disclosed** (the project previously flagged this value as
  disqualifying a "clean" verdict); **"traction-limited" label imprecise** (traction is at full strength; the
  real limiter is collective non-wetting, Douezan S<0); **no REPORT.md** (commit msg + viewer title only).

---

## C. Established scientific conclusions

1. **Active cortical-tension γ-floor** = a force-magnitude/engaged-density limit, ~1700× under the MCF7-active
   band, robust to every ablatable mechanism + framing. The ONE open lever = the missing experimental
   **MCF7-adherent load-engaged NMII density** datum (wet-lab). Model is sound; the gap is experimental.
2. **Aggregation** reproduced (clean faceted foam).
3. **Spreading**: cohesive spheroids compact / don't spread (structural limit); collective magnitude
   (A/A0 7–10) belongs to a coarser model — its FORM (a+b/R+c/R²) already reproduced by the Layer-2 line.
4. **Proliferation**: growth-driven "interior cells deform more" (T47D) signature reproduced (qualitative).

---

## D. Repo state

- ✅ Working tree clean (no uncommitted source); KB gates GREEN (runs 31 / params 42, none drifted); no
  conflicts / merge-in-progress; FF↔DCM files disjoint on dcm/main.
- ⚠️ **`aleph/outputs/tag_kb/kb.duckdb.tmp/` = 75 GB** — DuckDB temp spill from an interrupted TAG-KB refresh
  overnight. Live `kb.duckdb` (41 MB) intact, gates GREEN → KB integrity unaffected. Safe to delete:
  `rm -rf aleph/outputs/tag_kb/kb.duckdb.tmp`.
- ⚠️ Stray untracked `force.txt` + `messages.cmo` at repo root (not gitignored — accidental-add risk).

---

## E. Open items / PI decisions

1. Clean the 75 GB `kb.duckdb.tmp/` (safe) + gitignore/remove the stray root files.
2. **DCM Phase 5/6 documentation hardening**: multi-seed ensemble for the gradient, commit the slope/analysis
   code, write a Phase 5/6 REPORT.md, disclose the pen_frac 3.1 confound, correct the "traction-limited" label
   → "collective non-wetting". Confirm div-rate 0.08.
3. **CLAUDE.md ↔ KB actin-density conflict** (~38000 native = 30/µm² vs KB-3.18 100/µm²) — FF resolved toward
   KB-3.18 as a citation fix; the reconciliation is a PI call.
4. **Next FF unit** (WLC constitutive / loaded-shell production / Layer-2 coupling) — PI-roadmap-gated.
5. **The one measurement that could close the γ-floor** = wet-lab MCF7-adherent load-engaged NMII density.

---

## Verification provenance

7-agent read-only workflow (`overnight-report`, run wf_741f9651-395): 4 gather agents (DCM proliferation, DCM
spreading/goal, FF sprint audit, repo health) + 3 adversarial verifiers (T47D gradient → PARTIAL, spreading
limit → PARTIAL, FF γ-floor → CONFIRMED). 185 tool calls; every DCM number recomputed from committed npz.

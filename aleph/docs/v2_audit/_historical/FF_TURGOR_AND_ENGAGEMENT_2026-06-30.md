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

# Turgor (state-dependent osmotic) + motor engagement fraction — workflow verdicts (2026-06-30)

Two ultracode workflows (28 agents, cross-checked vs project files) resolving PI questions #1 (turgor
should be state-dependent) and #4 (force-bearing/engaged density). Both decisive.

## 1. Turgor — PI was right: state-dependent osmotic, NOT a fixed band-implied 133 Pa

- **Resting interphase turgor ≈ 40 Pa** (Fischer-Friedrich 2014 Sci Rep 4:6213, HeLa interphase
  ΔP=40±30 Pa, γ=0.17 mN/m; companion to Stewart 2011 Nature). MCF7-specific datum does NOT exist
  (HeLa proxy → flag to PI). The project's 133 Pa is ~3× this (partway-to-mitotic; mitotic ~320–400 Pa)
  and was band-implied (ΔP=2γ/R) — circular.
- **Must be state-dependent.** Animal cells (no wall) hold no static turgor: pump-leak (Kay & Blaustein
  2019) cancels the Donnan gradient → intracellular≈extracellular osmolarity, leaving only a sub-mM net
  excess (tens–hundreds Pa) balanced by/coupled to cortical tension via Laplace. The total osmotic each
  side is huge (~725 kPa, Wennerström 2022) but MATCHED; turgor is the net Δc·RT only.
- **First-principles + state law (implemented in FF):** Guo 2017 (PNAS 114:E8618) entropic
  excluded-volume closure Π(V)=N·kB·T/(V−Vmin), Vmin≈0.30·V0 (Venkova 2022 eLife / Adar 2025 Ponder
  fit), N from c≈200 mM cytoplasmic osmolytes (cross-validated vs cytoplasmic salt → not circular).
  This **DERIVES the bulk modulus** K_vol=Π_in0/(1−vmin)≈7.4e5 pN/µm² — which ≈ the DCM production
  value 7.73e5 (independent cross-validation) — so K_vol is no longer a free magic number. Young-Laplace
  ΔP=2γ/R demotes from the INPUT defining ΔP0 to an OUTPUT/consistency check (correct causal direction:
  Stewart 2011 osmotic generates, cortex maintains).
- **Applied (FF):** `gamma_floor.py` — `turgor_pressure` now the Guo state-dependent law, dP0=40 Pa,
  self-consistent V0 (`R0_mean`); rest ΔP=40.00 exactly, rises steeply on compression (Venkova
  volume-regulation regime). Kills TURGOR_DP0=133 + TURGOR_K_VOL=1e3 magic numbers. Passive turgor γ
  now 0.2 mN/m (sub-band) — consistent with "band is myosin, not passive".
- **DCM-side (PI / coordinate):** the same fix in `geometry.py` ResolvedDCM (turgor_dP0 133→40,
  K_vol→Guo closure) + the `dcm_warp` osmotic kernels (which already have the right structure:
  `_dp_from_vol_osm` + `osmotic_relax_kernel`). DCM territory + /loop dcm live → PI_DECISIONS Decision 1.

## 2. Motor engagement fraction — DECISIVE NEGATIVE: cannot close the floor (confirms density limit)

- **Engagement is sub-unity at every level:** head-level duty ratio NMIIA ~0.05–0.10 / NMIIB ~0.23–0.36
  (Stam 2015, Kovács 2003), up to ~0.52 under load (catch; Erdmann-Schwarz, Kovács 2007); filament-level
  engagement ~0.65 in interphase cortex (Truong-Quang 2021, the cleanest present-vs-engaged measurement),
  ~1.0 in mitosis; the project's own model realizes ~0.11.
- **So effective force-bearing density = engagement × 0.6/µm² (Nie 2015) is LOWER, not higher:**
  ~0.03–0.39/µm² vs the ~16–21/µm² needed → the gap WIDENS from ~30× to ~40–700×.
- **Verdict:** engagement/duty/clustering/stacking CANNOT manufacture force-bearers the cortex doesn't
  physically contain — they reorganize the same ~0.6/µm². The shortfall is a PRESENT-DENSITY /
  per-minifilament-force limit. The only legitimate lever is a genuinely higher *present* density
  (MCF7-specific, or a true SF cross-section count), which would need to be ~26–35× above two agreeing
  HeLa measurements = biophysically implausible. Confirms the γ-floor verdict (CORTICAL_TENSION_RECORD
  _2026-06-30) and the SF generation limit.

## Net

PI question #1 resolved (turgor is state-dependent osmotic; FF fixed, DCM PI-pending). PI question #4
resolved as a decisive NEGATIVE (engagement is a debit, not a rescue). Neither changes the γ-floor
verdict; both make the model more correct and kill magic numbers. No tuning.

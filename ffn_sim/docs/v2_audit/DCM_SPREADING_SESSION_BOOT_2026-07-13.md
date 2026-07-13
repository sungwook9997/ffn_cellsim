# DCM SPREADING SESSION — boot prompt (2026-07-13)

> Hand this file to a fresh DCM session. It is the bootable "프롬프트"; the full **계획서** is
> [`DCM_SPREADING_ACTIVE_AREA_PLAN_2026-07-13.md`](./DCM_SPREADING_ACTIVE_AREA_PLAN_2026-07-13.md)
> (read it in full at minute 0). Authored by the FF session as a cross-session handoff — the
> module map below is a discovery starting point, not a mandate on DCM internals (you own the engine).

## Your task (one line)
Make the DCM cell **spread** (single cell A/A0 → 2–8× stable pancake; then aggregate wetting) by adding the
**one missing physics — ACTIVE AREA GENERATION** — at physiological strength, so it beats the shell's
area-minimising surface tension + turgor. Prior attempts failed by being *under-powered*, not by a structural limit.

## Boot protocol (minute 0)
1. Read `CLAUDE.md` fully (esp. the physiological-baseline HARD rule + native-full-compartment HARD rule).
2. Read the plan `DCM_SPREADING_ACTIVE_AREA_PLAN_2026-07-13.md` (diagnosis, 3 terms, G1–G6 gates, increments).
3. Read the Dev Logs status board + Open items.
4. `git log --oneline -5`; confirm branch `dcm/main` + latest commit. `conda activate ffn_sim`.
5. Restate the task; start **Increment 1** (below). Ask PI only if genuinely ambiguous.

## The diagnosis (why this, and why prior attempts failed)
The DCM cell is a **turgid closed shell whose surface tension MINIMISES area** → it structurally OPPOSES spreading.
The area-*generating* machinery is absent/too weak. Real spreading = the cell does WORK to grow its footprint
(lamellipodial actin polymerisation + membrane unfolding), held by clutch traction, past the passive wetting angle.
Governing law = the active-wetting / Young–Dupré balance: `spread ⇔ W_adh ≳ γ  [+ active protrusion beyond passive]`.
Prior: `project-dcm-spreading-foundation` (single A/A0→3.8 but spheroid COMPACTS), `project-dcm-lamellipodium-graft`
(single→1.05, 400-cell RE-COMPACTS) — the drivers were **too weak to beat physiological γ + turgor**. The fix is to
ground every term physiologically AND add the missing area generator, and **prove the balance favours spreading
(G2) BEFORE any production run**.

## Implementation map (FF-session discovery — verify, then own)
The engine ALREADY has two of the three terms. **T-2 is the genuinely new build.**

- **T-1 — substrate wetting adhesion: EXISTS.** `ffn_sim/dcm/dcm_substrate_warp.py` →
  `DcmSubstrateForceGPU` (z-well `k_well=2·W_cs/rng²` + rigid dish floor) and `DcmSubstrateWettingGPU`
  (the in-plane xy-gradient of `U_adh = −W_cs·A_contact` over basal contact triangles — the conservative
  wetting drive, NOT a body-force proxy). **Job:** ground `W_cs` at the physiological integrin areal
  adhesion energy (~0.1–1 mJ/m²; pull exact from KB) and confirm the wetting drive is ON in the production config.
- **T-2 — rest-area growth `dA0/dt` (the AREA GENERATOR): NEW.** This is the core missing term. Locate how the
  DCM represents preferred/surface-tension area (the interfacial-tension term, `dcm_interfacial_tension_warp.py` /
  the surface-energy per-face law — the explicit per-cell `A0` may not exist yet; **confirming its representation
  is your first engine task**). Add: at basal-adhered / leading nodes, actively grow the preferred area while the
  edge is adhered and below a membrane-reservoir cap. Ground `k_prot` from the lamellipodium protrusion rate
  (µm/min edge advance × edge length) and the ~2–4× membrane-area reservoir. NO magic numbers → derive from KB.
- **T-3 — clutch traction + tension anisotropy: INFRASTRUCTURE EXISTS.** The active stack's `ActiveRimTraction`
  (coarse lamellipodium + contraction belt + FA clutch) is driven by `scripts/dcm_active_spheroid.py`
  (two-phase: aggregation → active spreading; flags `--f-act --spread-blocks --migrate-factor`), plus
  `dcm_ecm_clutch_host.py` and interfacial tension `dcm_interfacial_tension_warp.py`. ⚠️ **First verify the
  Warp-era location** — the legacy `cell/dcm_active.py` referenced by that driver is now under
  `archive/hoomd_legacy/cell/`; find/confirm the current Warp active-traction entry point before wiring T-3.
  **Job:** ground per-clutch traction (~5–25 pN, KB) and add the basal-vs-apical γ reduction so the pancake is stable.
- **Driver / readouts / biology-time:** `scripts/dcm_active_spheroid.py` (spreading driver, single-cell first
  via its single-cell path); A/A0 = **top-down xy silhouette** (`feedback-aa0-topdown-area`, NOT basal contact —
  that inflates 10–20×); metrics `dcm_jamming_metrics.py`, stress `dcm_virial_stress.py`; biology-time couples to
  the existing large-dt BDF2 + T1-rate `dcm_t1_rate.py` (`--t1-rate`) at Increment 4.

## First concrete step — Increment 1 (make-or-break, single cell)
1. **Ground the three constants** (`W_adh/W_cs`, `γ`, `k_prot`) at physiological values — pull from the KB
   (`tag_query.py`) + `project-mcf7-parameter-collection` (γ_MCF7 ≈ 1e-2 N/m, Moazzeni). Write down each with its
   KU/source. If a value is unknown, surface to PI — do NOT default to a convenient null.
2. **G2 balance check (the explicit fix for prior under-powered runs):** with the grounded values, verify
   `W_adh/γ + protrusion` predicts a spreading-favourable equilibrium (contact angle < 90°, A/A0 target > 2)
   BEFORE running. If it doesn't, the values or the mechanism are wrong — diagnose, don't run.
3. **T-1 + T-2 on a SINGLE native cell:** does it spread to A/A0 2–8× as a stable pancake, at constant volume?
   Gates G3 (V/V0=1), G4 (single-cell physiological spread), G5 (REAL flattening — basal area↑, maxZ↓, NOT
   node-ejection; **browser-verify the HTML viewer**, the metric has lied before). Only add T-3 (Increment 2)
   once a single cell spreads.

## Hard rules (from CLAUDE.md — do not violate)
- **Physiological baseline (HARD):** run FROM the resting physiological state — real cytoplasm viscosity, resting
  turgor, γ_MCF7 — never a convenient zero. Measure spreading as a perturbation FROM that baseline.
- **Native + full compartments (HARD, repeatedly stated):** validate/debug/visualize ONLY at native scale with
  every compartment ON. A coarse cell collapses on its own (artifact ≠ mechanism). A coarse/CPU run is a
  non-authoritative smoke test that MUST be reconfirmed on the native full cell before any conclusion. GPU work →
  gbook A5000 (`~/ff_scratch`, rsync edited files, full `~/miniconda3/envs/ffn_sim/bin/python` path, monitor by log file).
- **No magic numbers:** every constant derivable from KB/lit, grid-invariant, not chosen to pass a gate. `W_cs`,
  `k_prot`, γ, traction all get a source. If it can't be grounded, halt → PI.
- **No gate-loosening:** G1–G6 are contracts. If a gate is wrong, surface to PI for a contract change.
- **No production run until G2 confirms the balance favours spreading.**
- **Visualize at closeout:** interactive HTML cell-morphology viewer (`dcm_mesh_viewer_html`, peel/slab/cut),
  browser-verified (`scripts/browser_check.py`), full-res (no downsampling) + per-cell stress/strain. Not just PNG.

## File ownership / concurrent sessions (READ)
`dcm/main` is **shared by concurrent DCM(T1-KMC) + FF/M6 sessions**. Therefore:
- **Commit specific files only — NEVER `git add -A`.** Stage the exact paths you changed.
- You own the spreading modules (`dcm_substrate_warp.py`, the new T-2 area term, the active-traction wiring) +
  `scripts/dcm_active_spheroid.py` + the spreading outputs/figs. Do not touch FF (`ff/`, `outputs/ff/**`) or the
  other DCM session's in-flight files (`h_dcm_two_stage/**` T1-KMC viewers, etc.).
- Native full-res HTML viewers are gitignored under `outputs/ff/figs/` already; for DCM outputs, keep oversized
  (>100 MB) viewers LOCAL / release-uploaded (GitHub 100 MB limit) — commit PNG + curves + JSON as the record.

## Sanity gates (before first production — from the plan)
G1 dimensional · G2 balance (spreading-favourable at physiological values) · G3 volume conservation (V/V0=1) ·
G4 single-cell physiological spread first · G5 real flattening not node-ejection (browser-verify) · G6 aggregate
WETS not re-compacts (cadherin cohesion sets wet-vs-cohesive).

## Closeout (MUST — end of session / freeze-point / PI "wrap")
1. Commit on `dcm/main` (specific files). 2. Dev Logs status board + Open items + a `DCM spreading — Increment N`
milestone day-log (start/end commit, gate PASS/FAIL, KU refs). 3. Refresh figures (browser-verified) at any
milestone/finding. 4. `make kb-check` green (halt→PI on DRIFT). 5. If Notion unavailable → halt, surface to PI.
6. Final message ends with the literal line `Notion 업데이트 완료`.

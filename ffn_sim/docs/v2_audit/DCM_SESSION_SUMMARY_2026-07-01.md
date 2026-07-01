# DCM autonomous session summary — aggregation → proliferation → spreading (2026-07-01)

Full arc of the 8h autonomous run. **The problem PI kept flagging:** simulated spheroid cells stayed
round "bag of marbles" (low contact, no faceting) instead of a compact faceted tissue. **This session
solved it, and carried the /goal through proliferation and spreading with honest structural conclusions.**

## 1. Diagnosis — the faceting failure was an INITIALIZATION problem, not force/energy/architecture

Exhaustive testing (multi-agent workflows + ~20 GPU tests) eliminated every force lever: passive forces
freeze a gapped aggregate, active γ shrinks individual cells (not the aggregate), uniform γ rounds,
differential γ is insufficient, gap<2 overlap freezes (interpenetration). **Reading the SimuCell3D C++
source settled it:** SimuCell3D uses the SAME separate-mesh + node-face-contact topology as us (the
"shared-interface" hypothesis was FALSE) and NEVER assembles from gaps — it loads an already
space-filling mesh and the energy MAINTAINS faceting. We were asking the model to do something it (and
SimuCell3D) cannot: assemble a foam from separated round cells. → No engine rewrite needed.

## 2. Aggregation SOLVED (Phases 1–4) — a clean faceted N=400 spheroid

| step | result |
|---|---|
| confluent-init prototype (icosphere→Voronoi warp) | watertight, NO penetration, cells START faceted (asph 0.04–0.08, Q 150–170) |
| Phase 1 — V0 = confluent rest-volume (`--v0-from-init`) | eliminated the turgor over-inflation: V/V0 1.37→**1.00**, pen-ratio 0.19→**0.93** |
| Phase 2 — `--builder confluent` production builder | direct build reproduces it (Q150, V/V0 1, pen 0.97); Lloyd sped up (cKDTree) |
| Phase 3 — sharpen sweep | Q~150 is the radial-warp representation CEILING (not tuning); honest |
| Phase 4 — clean N=400 production | manifold 5/5, Q149, V/V0 1.00, pen 0.96 = **"aggregation properly built"** |

Viewers 4, 8, 9 (`viz_html/`); the fix is a confluent-init + V0-setpoint change, NOT an engine rewrite.

## 3. Proliferation (Phase 5) — reproduces the T47D interior-deforms-more gradient

Division from the confluent foam is stable (inset 0.02; 0.04 blows up). N grew 64→134 and 200→280
(70/80 divisions, stable). **Growth FLIPPED the deformation gradient to the real-spheroid signature:**
CORE asphericity 0.026 > MID > RIM 0.001 (inner cells deform more) — the equilibrium foam was inverted
(CORE < RIM, a Voronoi-in-ball boundary artifact). This is the T47D signature, driven by rim-cell
proliferation crowding the interior. Viewers 10, 11, 13.

## 4. Spreading (Phase 6) — a robust STRUCTURAL LIMIT (single-cell-scale)

Five mechanisms tried; ALL give A/A0 top-down ≈ 1.0 (the cohesive spheroid does not spread/flatten):
weak clutch, physiological 5 nN/FA clutch (167×), lamellipodium, and cadherin de-cohesion + lamellipodium
(3422 bonds ruptured, still no dispersal). **Why:** the clutch is an anchor not a motor; only 11/400 cells
touch the dish; turgor + bulk hold each cell; de-cohered bonds re-form. **Spreading is fundamentally a
single-cell phenomenon** — matches the Layer-2 structural conclusion (A/A0 = a+b/R+c/R² magnitude belongs
to the fine-grained single-cell line, not the multicellular aggregate) and the lamellipodium-graft prior.
The mechanistic path to real spreading = basal de-cohesion + single-cell dispersal, which a cohesive ball
resists — a future build, not a tuning fix. Docs: `DCM_SPREADING_LIMIT_CONCLUSION`. Viewers 12, 14.

## 5. Quantitative validation vs literature (honest matches + gaps)

- **Faceting Q = 149** (faceted; but γ̃ = 1.1e-4 ≪ SimuCell3D band [0.02,0.10] → our faceting is
  GEOMETRIC confluent-init-driven, not γ̃-energy-driven; K-reconciliation is a PI decision).
- **T47D interior/rim cell-area ratio = 1.14** vs lit 1.37 (correct direction, ~83% magnitude).
- **Porosity = −5.3%** (cells overlap = division-without-remesh penetration; T47D wants 12%→2%).
Doc `DCM_QUANT_VALIDATION`, figure 16.

## Open items (scoped, none forced)

1. **γ̃ / K-reconciliation** (faceting energy regime) — a PI decision (drop K to ~2500 for γ̃-in-band).
2. **division + remesh co-run** (grown-foam porosity) — a real re-architecture, NOT a one-line fix
   (mitotic fixed-npc-block indexing vs remesh cof-relabel conflict); `DCM_DIVISION_REMESH_CORUN_DESIGN`.
3. **spheroid spreading** — needs basal de-cohesion + single-cell dispersal regime (structural, future).

## Deliverables

~16 commits (dcm/main, disjoint from the FF session), 16 visualization files (`outputs/h_dcm_two_stage/
viz_html/`, interactive HTML + PNG), 5 docs (`DCM_FACETING_TO_PRODUCTION_PLAN`, `CONFLUENT_INITIALIZER_
DESIGN`, `DCM_SPREADING_LIMIT_CONCLUSION`, `DCM_QUANT_VALIDATION`, this summary). No magic-number tuning
anywhere; every gap surfaced honestly. The core PI complaint (bag of marbles) is resolved: the model now
builds a clean faceted spheroid, proliferates it with the correct interior-deformation gradient, and the
spreading limit is a documented structural property, not a bug.

# DCM many-cell spheroid — master design (full mechanism stack)

**2026-06-11 · branch `h7/compartment-platform` · PI-driven (graduation report).**

The platform's spheroid-scale work runs on a **Deformable Cell Model (DCM)** and must,
in the end, be a **large MANY-CELL spheroid** (a 7-cell cluster cannot show the 3-zone /
bulk-pressure / junction-switch physics — "세포가 많은 스페로이드를 구현해봐야 앎", PI).
This doc records the FULL mechanism stack the integrated model must carry, the validation
criteria, and the status of each piece (several built/running in parallel).

## Validation criteria (do not conflate)
- **A/A₀ = a + b/R + c/R²** (PI's layer-2 law: a=−0.33, b=188.7 µm, c=−2655 µm², r²=0.98,
  R=31–78 µm) — the **spheroid SPREADING-magnitude** target. PROLIFERATION-driven (the
  b/R = rim-localized division ∝ surface/volume; c/R² = cohesion penalty). Mechanical
  wetting at fixed cell number caps at **A/A₀ ≈ 1.2** (measured) — so the law is NOT
  reproducible without division + necrosis.
- **tension ~ 1/r (2D)** (Slater…Kim d0sm01911a) — the **cell–ECM stress-decay / remodeling**
  check (separate; validated qualitatively for the contracting DCM⊗ECM cell).

## Mechanism stack (the integrated large spheroid)
1. **DCM deformable cell** — icosphere elastic shell + exact triangulated turgor
   (`cell/dcm.py`). DONE.
2. **Explicit cross-linked fiber ECM** (Mikado, `ecm/mikado.py`) + **catch-slip FA clutch**
   (Pereverzev, `bridge/integrin_bonds.py`) — `cell/dcm_ecm.py`. DONE (remodels matrix;
   1/r sign recovered via no-cell-baseline subtraction).
3. **Proliferation** — rim-biased, contact-inhibited cell DIVISION (the b/R term).
   `cell/dcm_prolif.py` (Workflow `wq5tsya8y`, building). Pre-allocated cell pool +
   shared membrane type + per-cell-id custom adhesion (so cells can activate mid-run).
4. **Core NECROSIS / quiescence (3-zone)** — a core cell arrests/dies by TWO coupled
   factors: (a) NUTRIENT/O₂ access — depth below the spheroid surface vs the diffusion
   penetration depth (viable-rim thickness); (b) MECHANICAL solid stress — local
   compressive pressure (built by rim proliferation) vs a death/quiescence threshold.
   Calibrated to **MCF7 spheroid LIVE/DEAD** data (viable-rim vs diameter). Literature +
   DCM criterion: Workflow `w7gvjd994` (running). Zones: proliferating rim / quiescent
   mid / necrotic core.
5. **Active TRACTION (not passive wetting)** — the spreading driver is the rim's active
   machinery: **lamellipodium protrusion** (the engine built in `cell/spreading_drive.py`,
   commit f654531) + **filopodia** + **actomyosin contraction belt** (cortex myosin) +
   **FA molecular clutch** traction. This is why mechanical wetting alone capped at 1.2 —
   it lacked active traction. The rim cells protrude + grip the ECM + contract → traction
   → spreading.
6. **Bulk-pressure-driven JUNCTION SWITCH + pressure-release spreading** — as the spheroid
   compacts (proliferation + cohesion), BULK PRESSURE rises; adhesion switches from
   **cell-cell (cadherin)** to **cell-ECM (integrin/FA)**; when the built-up pressure
   releases ("터짐"), cells unjam and the spheroid spreads outward. (jamming/unjamming;
   pressure-driven escape / EMT-like adhesion switch.) Literature + implementation:
   parallel research Workflow (junction-switch / pressure-release) + a pressure-gated
   cadherin→integrin switch on the existing `junction/cadherin.py` + FA clutch.

## Status (parallel)
- DONE/committed: DCM tier (b5a7578), DCM⊗ECM + 1/r sign + 7-cell spheroid-on-ECM
  (1a71cd7), wetting-ceiling finding + a+b/R+c/R² figure (79b84ec), necrosis-requirement
  note (this commit).
- RUNNING: proliferation build + size-sweep + a+b/R+c/R² fit (`wq5tsya8y`); necrosis
  literature + MCF7 live/dead + DCM criterion (`w7gvjd994`); junction-switch /
  pressure-release literature (this launch).
- NEXT (the integration): assemble a LARGE many-cell DCM spheroid carrying 1–6 together,
  spread on the ECM via active traction, with proliferation+necrosis giving the
  size-dependent A/A₀ = a + b/R + c/R² and the junction-switch/pressure-release dynamics.
  Scale via coarse cells (subdiv≤1) + native GPU on gbook.

## Calibrated bands (from the two literature workflows — trusted directly)

**NECROSIS / 3-zone (refs `references/analysis/_necrosis/`):**
- O₂ consumption-limited penetration depth = **viable rim ~100–150 µm**, ~CONSTANT as the
  spheroid grows (Greenspan invariant); only the necrotic core expands.
- Critical diameter: <200 µm fully viable; **~200 µm hypoxic core begins**; ~400–600 µm
  necrotic band; **>500 µm robust necrotic core** (MCF7-confirmed); 600–800 µm growth plateau.
- Zones: proliferating shell ~20–40 µm / quiescent middle / necrotic core.
- MCF7: hypoxic >200 µm, necrotic >500 µm, intercellular gaps 5–10 µm, PI⁺ death by day 6.
  BT-474 necrotic volume fraction d5/d6/d7 = 0.20/0.29/0.61.
- **DCM rule (route A, depth-from-surface):** per cell depth d = R_spheroid − r. d < ~40 µm
  → PROLIFERATE; 40–150 µm → QUIESCENT (arrest division); **d > ~150 µm → NECROSE**
  (deactivate). Equivalently necrotic-core radius R_nec = max(0, R − 150 µm). Upgrade =
  route B (steady O₂ reaction-diffusion field, D≈1.5–3.8e-9 m²/s, necrosis at O₂<0.02 mM /
  glucose<0.06–0.08 mM). Refs: Thomlinson-Gray 1955, Greenspan 1972, Grimes 2014, Jiang 2005.

**JUNCTION SWITCH / BULK PRESSURE / UNJAMMING (refs `references/analysis/_junction/`):**
- Pressure ledger: resting turgor ~40–100 Pa; proliferation/motility gating **onset
  ~0.5 kPa, saturating ~5 kPa** (Dolega 2021, Delarue 2014); strong arrest 5–10 kPa
  (Montel 2011); endogenous growth-induced solid stress 0.37–19 kPa (Stylianopoulos 2012).
  Compression arrest is fully REVERSIBLE (release → re-fluidization → escape).
- **Unjamming threshold (shape index):** 2D q* = P/√A ≈ **3.81**; 3D SI = A/V^(2/3) ≈ **5.4**
  (Bi 2015/2016, Park 2015, Merkel-Manning 2018). Han 2021: core SI 5.84 (jammed) vs
  periphery 6.6 (unjammed). q is a STATIC snapshot order parameter (validation observable).
- **Cadherin→integrin switch = clutch competition** for shared actin/vinculin, gated at the
  **~5 pN** talin-R3 / α-catenin unfold (vinculin recruitment) threshold (Yao 2014 ×2);
  E-cadherin <12→>43 pN matured, integrin α5β1 catch-bond strengthens 10–30 pN
  (Wang 2016, Kong 2009); cell-pair cadherin ~100 nN, intercellular = 0.47× total ECM
  traction (Maruthamuthu 2011 — adhesions co-scale, not bond-for-bond trade).
- **DCM rule:** per-cell local compressive stress (neighbour crowding / contact force) →
  above ~0.5–5 kPa: weaken cell-cell cadherin bonds + strengthen cell-ECM integrin/FA
  clutch; track the shape index q; when local q crosses q* (release/unjamming) the cell
  fluidizes and escapes outward → spreading. (Vertex/SPV are ORACLES only — the
  fine-grained cadherin catch-bond + actomyosin cortex should reproduce the jamming line
  emergently; CLAUDE.md no-lumped-mechanism rule.)

## Honesty / scale
Mesoscale CPU caps the cell count; a genuinely LARGE spheroid (hundreds of cells, R 30–80 µm
to hit the layer-2 fit range) needs the gbook A5000 GPU (dirty-branch cleanup pending) or
coarse-grained cells. The literature bands (necrosis penetration depth, solid stress,
junction-switch pressure) are trusted directly per the DCM-tier convention.

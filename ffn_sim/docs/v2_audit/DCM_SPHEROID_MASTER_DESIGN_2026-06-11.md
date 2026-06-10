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

## Honesty / scale
Mesoscale CPU caps the cell count; a genuinely LARGE spheroid (hundreds of cells, R 30–80 µm
to hit the layer-2 fit range) needs the gbook A5000 GPU (dirty-branch cleanup pending) or
coarse-grained cells. The literature bands (necrosis penetration depth, solid stress,
junction-switch pressure) are trusted directly per the DCM-tier convention.

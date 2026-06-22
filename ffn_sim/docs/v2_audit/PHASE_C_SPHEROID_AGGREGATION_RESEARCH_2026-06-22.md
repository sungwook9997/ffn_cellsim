# Spheroid aggregation — research + decision (2026-06-22, ultracode workflow)

PI redirected: stop ad-hoc builder tinkering; research how (a) sim papers and (b) real biology
make spheroids, investigate our old GPU two-step aggregate (looked most spheroid-like), then
design the best method. 8-agent research workflow (web/PubMed/Consensus + local + SimuCell3D
deep-read + adversarial verify). HARD constraints: no magic-number tuning (rounding must EMERGE
from γ+cohesion), no one-by-one dynamic cell addition, result must be round + all-touching.

## ⭐ ROOT CAUSE — it is NOT the builder, NOT the drag.
The builder work (fcc faceted / voronoi loose / sphere RCP) was second-order. Measured from the
OLD aggregate the PI remembers (`outputs/h_dcm_two_stage/spheroids/spheroid_nf_n{400..1000}.npy`,
2026-06-13 "aggregation RESOLVED"), the real causes:

1. **Today's cells never DEFORM.** Every `dcm_warp_decohesion` run today: **V/V0 = 1.000 exactly,
   per-cell volume 1767 ± 0** = rigid icospheres translating as a pile of marbles. The OLD spheroid:
   **V/V0 = 1.11–1.13 with real ±50 spread** = turgor pressurised each cell and contact flattened it
   into polygonal junctions (Ψ 0.87). A pile of rigid spheres cannot round its envelope or close
   gaps; only deformable, pressurised cells flowing into junctions can. **The cortical surface-
   tension γ term is OFF by default** — the one force that flattens interfaces + rounds the envelope.
2. **The settle is orders of magnitude too short.** OLD ran ~30k+ steps from a GAPPED start; today's
   agg_N48 ran ~3.5k–14k steps = barely-started near-initial snapshots frozen near the placement.

**Exonerated:** (a) placement — OLD started gapped/non-touching (spacing 2.3) and rounded anyway;
voronoi starts near-ideal and still fails → builder is second-order. (b) drag — adversarial
verifier REFUTED "drag causes non-roundness": in overdamped dynamics high drag rescales relaxation
TIME, not the equilibrium; OLD ran at the SAME physiological η=65.9 Pa·s and rounded. Drag is a
wall-clock cost, not a cause.

## OLD spheroid vs TODAY (same metrics)
| metric | OLD n400 (spheroid) | TODAY fcc/voronoi/sphere |
|---|---|---|
| anisotropy λ1/λ3 | **1.02** | 1.10 / 1.03 / 1.13 |
| asphericity | **0.013** | 0.058 / 0.012 / 0.056 |
| contact @2µm | **1.000** | 1.000 / 0.896 / 1.000 |
| isolated cells | **0** | 0 / **5** / 0 |
| **V/V0 (turgor)** | **1.12** | **1.000** / 1.000 / 1.000 ← the key difference |
| per-cell Ψ | **0.87** | 0.95 / 0.99 / 0.99 (round = undeformed) |

OLD recipe: `dcm_two_stage_production --rep-strength 4e7 --adh-strength 5e7 --agg-spacing 2.3
--node-face-contact`, 30k+ steps, N=400–1000 (`scripts/agg_sweep_gbook.sh`).

## RECOMMENDED METHOD — replicate the OLD two-step PHYSICS (a relaxation-protocol + missing-γ fix, NOT an init change)
Respects all constraints (rounding emerges from turgor-distension + cohesion + γ; all N seeded
together = what real biology does too, prong B; round + all-touching):
- **A. Seed whole population gapped + soft:** N≥400 on a gapped spherical cluster (spacing ≈2.3),
  **soft lit contact `--rep-strength 4e7 --adh-strength 5e7`** (MCF7 ξ̄/ω̄, NOT the stiff 2e8 default
  → today's 2e8 causes ~11% interpenetration), `--node-face-contact`. (Biological "settled pile" IC.)
- **B. Make turgor ENGAGE** (audit why today's runs sit at V/V0=1.000) → cells distend to V/V0≈1.12,
  flatten into junctions. Deepest difference.
- **C. Turn ON cortical surface-tension γ during settle** (`surface_tension_kernel`) at the lit MCF7
  value (SimuCell3D Table-2 0.5–2.5e-3 N/m; ours ~1e-4 = low end). The emergent Steinberg rounding
  force — NOT a magic number (it's the measured cortical tension in the energy contract).
- **D. Run the FULL settle (~30k+ steps)** on the implicit IMEX integrator (buys larger dt to reach
  the ~5 s capillary time); do NOT lower the physiological drag.
- **E. Enable a REARRANGEMENT channel** (`--remesh-period` SWAP/SPLIT/COLLAPSE + cadherin break/form)
  so cells flow into gaps and relax facets — the difference between rounding and **arrested
  coalescence** (a stuck non-round pile). Rounding is rate-limited by active rearrangement /
  fluidization (verifier's strongest-supported physics), not passive relaxation.
- Co-set γ (C) and adhesion ω (A) to the Steinberg regime-II ratio (SimuCell3D Fig 3e/4c: ω̄=ωl/K
  vs γ̄=γ/Kl) to land round-AND-compact rather than loose.

## Adversarial-verification caveats (don't over-claim)
- SUPPORTED: rounding via surface-tension minimisation + active cell rearrangement/fluidization;
  real spheroids form by a population settling+compacting+rounding TOGETHER (not one-by-one).
- REFUTED: "RCP is the gold-standard sim init" — no single gold standard exists; init is 2nd-order.
- REFUTED: "free-space/high-drag prevents equilibrium" — drag rescales time only.

## First experiment
Re-run the OLD recipe on the WARP engine: N=400, gapped (2.3), rep 4e7 / adh 5e7, γ ON (lit),
remesh ON, implicit, full settle → measure V/V0 (must →~1.12), asphericity (→~0.01), contact
(→1.0), Ψ (→~0.87). If it reproduces the OLD spheroid, the diagnosis holds and the builder is moot.

## Honest open Qs
Why is turgor at V/V0=1.000 today (off, or not given time)? Is γ-on stable with the implicit step
at large dt? Does N≥400 implicit settle fit a tractable wall-clock?

Full research output: workflow wf_7d8cb163-5b7 (sim-methods, real-biology, simucell3d prongs).

# DCM cortical-tension (γ) controlled-variable sweep — 2026-06-29

Autonomous DCM-loop session, worktree `dcm/aggregation` (isolated from the FF session that
is committing Cytosim stages onto `dcm/main`). Non-gated lane: γ as a **controlled
independent variable** — vary γ, *report* the morphology that emerges, never tune γ to a
target (the magic-number rule). PI directive 2026-06-29: "코르티컬 텐션이 문제라면 …
스윕으로 텐션 여러개 동시에 돌려보자."

## What was built (all in `ffn_sim/dcm/`, my owned domain)

1. **`ffn_sim/dcm/gamma_sweep.py`** — a concurrent γ-sweep harness. Each γ point is one run of
   the parity-gated engine (`dcm_warp_decohesion`) as an isolated **free aggregate** (substrate
   well + wetting OFF, no division, no spread drivers) so γ is the only thing that changes.
   Workers launch as concurrent subprocesses (`--jobs`); each saves its per-frame mesh; the
   harness reads the final frame, computes morphology, and emits `gamma_sweep_summary.json` +
   `gamma_sweep_faceting.png` (faceting band overlaid). gbook-ready (`--device cuda:0`,
   `--n-cells`, `--subdiv`, `--jobs`); CPU-validated at small N.

2. **`--k-vol` exposed on `dcm_warp_decohesion`** (was hard-wired at `run_decohesion`'s
   `k_vol=7.73e5`). K is now a controllable variable and recorded in the run's out-dict. This
   was necessary to reach the faceting regime (below) and to make γ̃ = γ/(K·ℓ) auditable.

## Faceting dimensional finding (the headline)

The SimuCell3D faceting transition is governed by the dimensionless surface tension
**γ̃ = γ / (K·ℓ)**, ℓ = V₀^(1/3), with the faceting band **γ̃ ∈ [0.02, 0.10]** (Runser, Vetter
& Iber 2024). γ̃ depends critically on **K**:

| K [Pa] | source | γ̃ at γ = 1e-3 N/m | faceting reachable in measured γ range? |
|---|---|---|---|
| **7.73e5** | driver default (single-cell fried-egg spread tuning) | ~1e-4 | **NO** — γ̃ ≪ band for all γ ∈ [5e-4, 2.5e-3] |
| **2.5e3** | SimuCell3D / Fischer-Friedrich 2014 (HeLa cortex) | ~0.033 | **YES** — band crossed at γ ∈ [1e-3, 3e-3] |

→ **Cross-consistency issue surfaced for PI (not auto-fixed):** the K used for single-cell
spreading (7.73e5) and the K SimuCell3D faceting physics expects (~2.5e3) differ by ~300×. At
the driver default K, *no* cortical tension in the measured MCF7/HeLa range can facet the
aggregate (γ̃ stuck ~1e-4). This is the same class of question as the cortical-γ floor — which
K is the physiological osmotic/bulk modulus of the MCF7 cytoplasm, and is the single-cell-spread
K=7.73e5 a tuned value rather than a measured one? Surfacing, not resolving (PI call).

## Controlled sweep result (K = 2.5e3, N=12, subdiv 2, free aggregate)

`outputs/h_dcm_two_stage/gamma_sweep_facet_K2500/` — γ walked through and across the band:

| γ [N/m] | γ̃ | in band | mean V/V₀ | per-cell asph |
|---|---|---|---|---|
| 0 | 0.000 | – | 1.086 | 0.0020 |
| 5e-4 | 0.0167 | – | 1.033 | 0.0019 |
| 1e-3 | 0.0335 | ✓ | 0.980 | 0.0016 |
| 2e-3 | 0.0669 | ✓ | 0.868 | 0.0007 |
| 3e-3 | 0.1004 | – | 0.751 | 0.0002 |
| 5e-3 | 0.1673 | – | 0.507 | 0.0001 |
| 8e-3 | 0.2677 | – | 0.198 | 0.0001 |

**Monotone, physically-sensible turgor↔cortex balance as a function of γ̃:**
- γ̃ < band → **rounded marbles** (V/V₀ ≈ 1, turgor dominates; per-cell asphericity highest).
- γ̃ in band → cells **held but compressing** (V/V₀ 0.87–0.98) and rounding (asph ↓) — the
  tissue-faceting regime where cortical tension first competes with turgor.
- γ̃ > band → **cortex crushes the aggregate** (V/V₀ 0.75 → 0.20; runs stay finite, pen ≈ 0).

Caveat: per-cell asphericity (node-cloud gyration) measures whole-cell *roundness*, which
*decreases* under γ; true foam-like **contact-face flatness** is a distinct metric not captured
here (needs a contact-face planarity read-out — backlog). cluster asphericity is effectively
flat (~0.227): 12 free centroids barely rearrange over 8k steps — needs N≥400 + longer settle
on gbook for collective rounding.

## Division / proliferation baseline (CPU, `outputs/.../smoke/`)

Mitotic-round division (`--division`) works on CPU: parked icospheres activate, daughters split
along the Hertwig axis, V/V₀ regrows V0/2→V0; **mesh stays manifold and finite** (no break).
Interpenetration (pen ≈ 1.1–1.3) at subdiv 1 is a **contact-resolution artifact** (coarse mesh +
gap 2.05), not a division defect — faithful non-penetration needs subdiv 2 + IPC/remesh at
production N on gbook. Visualized: `div_smoke_n8_montage.png` / `_surface.mp4`.

## Iteration 2 additions (2026-06-29)

- **Contact-area-fraction faceting metric** added to the sweep (the direct SimuCell3D foam-like
  signal: fraction of surface on cell–cell contact faces). Re-run (K=2500, v2):

  | γ̃ | V/V₀ | contact-area frac |
  |---|---|---|
  | 0.000 | 1.086 | 0.078 |
  | 0.017 | 1.033 | 0.060 |
  | 0.033 (band) | 0.980 | 0.037 |
  | 0.067 (band) | 0.868 | 0.002 |
  | 0.100 | 0.751 | 0.000 |
  | 0.268 | 0.198 | 0.000 |

  **Honest finding:** in a *free aggregate at fixed centroid spacing*, rising cortical tension
  makes cells **shrink (V/V₀↓) and LOSE contact** rather than facet — contact area *falls* to 0.
  Faceting (contact-face flattening) needs the cells **held in apposition** (stronger cohesion or
  confinement) so the cortex flattens the junction instead of pulling cells apart. At the engine
  defaults the cohesion is overwhelmed by the soft-K cortex crush. → faceting emergence is a
  gbook-scale test with the "old-physics" recipe (soft rep + gap 2.3 + turgor + γ ON + remesh +
  long settle + N≥400; cf. [[project-phase-c-warp-migration]]), not a CPU free-aggregate.

- **Parity hardening:** the 10 DCM parity tests (`tests/dcm/*_parity.py`) were **silently
  SKIPPING** — they resolved fixtures via the pre-restructure path `warp_port/fixtures/` (the dir
  was renamed to `dcm/fixtures/`). Fixed the path component in all 10; `make warp-parity` now
  **runs 27 tests, all PASS, 0 skipped** (previously ~0 executed). This also confirms the
  `--k-vol` + comment edits did not perturb the engine.

## Iteration 4 — two-cell doublet faceting study (the decisive faceting test)

`outputs/h_dcm_two_stage/doublet_facet_adh{5e7,1e7}/` — N=2 free doublet, K=2500, γ ladder,
cohesion as a controlled variable (lit MCF7 5e7 vs default 1e7). Montage:
`doublet_facet_adh5e7/doublet_faceting_montage.png`.

| γ̃ | V/V₀ | contact-area (adh 5e7) | contact-area (adh 1e7) |
|---|---|---|---|
| 0.000 | 1.086 | 0.064 | 0.052 |
| 0.033 (band) | 0.981 | 0.051 | 0.039 |
| 0.067 (band) | 0.870 | 0.035 | — |
| 0.100 | 0.752 | 0.000 | 0.000 |
| 0.167 | 0.507 | 0.000 | 0.000 |

**Decisive finding (consistent with the iter-2 free aggregate):** a FREE doublet does **NOT
facet** even at lit MCF7 cohesion. Rising γ **deflates** the cells (Young–Laplace: equilibrium
V/V₀ = exp(−2γ/(R·K)) — predicts 0.81 at γ=2e-3, matches the measured 0.87) faster than node-node
cohesion can hold them apposed, so the contact area *falls to zero* (cells separate by γ̃≈0.1).
Higher cohesion (5e7 vs 1e7) only slightly delays separation. **Conclusion: faceting is a
CONFINEMENT phenomenon, not a cohesion one** — flat polygonal junctions need the cells held in
volume (tissue-interior confinement / enclosed-volume) so cortical tension flattens the interface
instead of freely deflating the cell. → faceting will emerge only in the **N≥400 confluent**
aggregate (interior cells confined by neighbours), the gbook production target; it is correctly
*absent* in small free clusters at any cohesion. The γ-sweep + nondim bridge are the ready tooling
for that run. (No magic-number: cohesion + γ are both controlled variables; the map is reported.)

## Iteration 3 — shared nondimensionalization (FF bridge)

`ffn_sim/dcm/nondim.py` (+ `tests/dcm/test_nondim.py`, 7 pass): the shared dimensionless language
so an FF-measured γ maps onto the DCM faceting regime without re-deriving scales. `describe_gamma(γ,
scales, K)` returns γ̃ = γ/(K·ℓ) + regime + companion balances (elastocapillary length, Young–Laplace
turgor/capillary ratio, Douezan s, viscous-capillary time). Measured scales only (R=7.5µm Wagner
2011, ℓ=12.09µm, η=65.9 Pa·s Dessard 2024, w_cs=2.85e-3, dP₀=133 Pa); γ stays a measured/controlled
INPUT, never tuned. Reinforces the K finding quantitatively: **at the driver K=7.73e5, faceting
needs γ ∈ [0.19, 0.94] N/m (~100× the MCF7 cortical tension) — unphysical**; at SimuCell3D K=2.5e3 an
FF-style γ=1e-3 → γ̃=0.033 (in band). Prototyped in `dcm/`; promote to `common/` after FF review.

## Co-run (division+remesh) — deferred with a design note

See `DCM_DIVISION_REMESH_CORUN_DESIGN_2026-06-29.md`: the "−2 sentinel one-liner" is insufficient
(parked-cell nodes double as cleave's dormant supply; mitotic uses fixed npc-blocks incompatible
with remesh relabelling). Two approaches recorded (index-range vs unified pool manager) +
the known cleave+remesh quality caveat; validate on gbook. Deferred from the CPU loop.

## Next (queued; gbook OFFLINE blocks N≥400 production)

- **cof-sentinel disambiguation** (parked cell = −2, remesh pool = −1; `split_edge` selects only
  −1): unblocks division + remesh co-run (currently force-disabled, `dcm_warp_decohesion.py:280`).
  Only the topology-agnostic *cleave* path is remesh-compatible (mitotic uses fixed npc-blocks).
- gbook production: faceting sweep at N≥400 subdiv 2 (collective rounding); division+remesh at
  the soft-K faceting regime; substrate hand-off (settle→spread) stability.
- PI decision on the **K reconciliation** (single-cell-spread 7.73e5 vs SimuCell3D ~2.5e3).

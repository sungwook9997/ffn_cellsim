# Dev Log — DCM spheroid: cohesion fix → lamellipodium → TWO-STAGE biological aggregation (2026-06-12)

> Session writeup for Notion Dev Logs (the Notion MCP connector dropped with the API
> reconnect, so this is the durable local copy — paste into the Dev Logs page
> `37c120da-ec5d-817d-86f7-db1ce13b8edf`, where ~8 progress comments were already posted
> live during the session). Branch `h7/compartment-platform`. PI prefers 한글; this is the
> consolidated record.

---

## Part A — cohesion fix + lamellipodium + arrest (h_dcm_gpu_lod)

**⭐ ROOT CAUSE found — cell-cell adhesion was DEAD.** The SimuCell3D bilinear tent only
adheres on node pairs in `[r_contact, c_adh)`. The GPU build had `c_adh=0.5µm` but
`r_contact = mean_edge = 4.37µm` (subdiv 1) → the band is EMPTY → cohesion silently dead.
The "spheroid" was a gapped lattice of mutually-repelling balls that dispersed (the old
hull-A/A₀ blow-up was scatter, not over-spread). Verified `scripts/dcm_cohesion_check.py`:
contact-node frac **0.0** at t=0 and after relaxation.

**FIX** — validated confluent bands into `ResolvedGpuDCM` (kept the stiff cortex + turgor):
`c_adh 0.5→5µm`, `adh 1e8→8e8`, `rep 1e8→2e8`, `spacing 2.3→2.0R`. Verified: **73% node-node
contact from t=0, stays bound** (Rg flat, no dispersal, no explosion).

**Spreading-arrest → default OFF (later fully removed by the parallel session).** Diagnostic
`probe3_on.json`: with arrest at full strength (gain→0, traction off) A/A₀ STILL blows up →
arrest is not the lever; it is a lumped, band-tuned over-spread cap (no-magic-number
violation), so it was disabled. This also fixed the failing parity test
`test_switched_integrin_gain_bit_parity` (GPU twin pinned arrest-OFF = the legacy law).

**Lamellipodium visualisation (PI request — "보이게"):** `scripts/dcm_lamellipodium_viz.py`
renders the rim lamellipodium on the cohesive spheroid as real cell SURFACES (not dots) —
`figs/lamellipodium_surface.{png,mp4}` (red basal caps = lamellipodia wetting the dish,
translucent bodies) + an analytic anatomy figure. Rim/basal/outward geometry recomputed
exactly as `ActiveRimTraction` and cross-checked against the live force object.

**VOLUME check (PI asked):** per-cell V/V0 stays ~1.06 flat through spreading — cells keep
their volume (turgor) while the spheroid SHAPE flattens. Not deflation.

Tests: fast DCM parity **14/14 PASS** (incl. the previously-failing one); fused force
`_footprint_radius` borrow-list bug fixed. Commits `007fd07`, `f7dd0d3`.

---

## Part B — TWO-STAGE v1 (REJECTED by PI)

First attempt at aggregation→spreading used an `AggregationDrive` central-pull force. **PI
correctly rejected it:** it is a forced central potential, not biology ("강제로 가까워지게 한
거"), and the cells barely moved (~13% squeeze). Withdrawn. Commit `c171bf8` (superseded).

---

## Part C — TWO-STAGE v2: BIOLOGICAL aggregation (active-matter coalescence)

**PI directive:** spreading must start from a spheroid that went through REAL aggregation;
reference SimuCell3D organoid aggregation; ALL cells full physics (no LOD/outer-only);
verify the aggregation is biologically correct BEFORE spreading.

**Mechanism (4-agent ultracode research → `AGGREGATION_REDESIGN_RESEARCH_2026-06-12.json`):**
active-matter COALESCENCE (search-and-capture), NOT an external pull:
- `DcmActiveMotilitySPP` — per-cell self-propulsion along a polarity `p_c`, NET `f_active·p_c`
  over the cell's nodes. ALL live cells (no rim gate, no LOD = the all-cells rule). Raises
  the collision rate so motile cells encounter + cadherin-adhere + coalesce.
- `DcmPolarityUpdater` — Ornstein-Uhlenbeck rotational diffusion (persistence τ_p), off the
  per-step force path (BAOAB-safe), |p|=1.
- `DcmDropConfinement` — soft one-sided spherical wall = hanging-drop meniscus. **Force-free
  interior**, inward only at r>R_drop. Bounds the search; does NOT aggregate/pre-strain/
  manufacture surface tension (explicitly not a central pull).
- `test_dcm_motility_gpu_parity.py` (4 PASS): net force per cell = `f_active·p_c`, every live
  cell propelled (none frozen), polarity unit-norm.

**✅ EMERGENCE VERIFIED (the PI's "생물학적으로 맞는지 확인" gate) — negative control:**
seeding cells JUST OUT of adhesion range (2.8R) so motility is required:

| | motility ON | motility OFF (negative control) |
|---|---|---|
| contact fraction | 0 → 0.05 → 0.11 → 0.16 → 0.19 (monotone rise) | 0 → 0 → … → 0 (never aggregates) |

Motility off → cells NEVER aggregate. Motility on → emergent coalescence + adhesion. The
aggregation is driven by the cells' OWN active motility, not a force or pre-placement.

**GPU production sweep (gbook A5000) — emergent, dense, round at scale:**

| N | aggregation contact (emergent) | asphericity | Rg compaction | spreading footprint | converged |
|---|---|---|---|---|---|
| 100 | 0.00 → 0.69 | 0.024→0.014 | 41.8→37.0µm | 365→7 443µm² | ✅ |
| 200 | 0.00 → 0.80 | 0.051→0.022 | 52.6→45.4µm | 389→10 430µm² | ✅ |
| 400 | 0.00 → 0.90 | 0.020→0.006 | 66.3→56.0µm | 1 179→13 666µm² | ✅ |
| 600 | 0.00 → 0.93 | →0.002 | 76→63µm | 515→15 890µm² | ✅ |
| 800 | 0.00 → 0.95 | →0.002 | 83→68µm | 1 588→17 369µm² | ✅ |

Denser AND rounder with N (more 3-D neighbours): a near-fully-confluent (contact 0.95)
near-perfect sphere (asph 0.002) by N=800. Volume held throughout (V/V0 ≈ 1.05–1.09).
Hand-off exact (Stage-2 frame 0 == Stage-1 final). N=800 (33 600 particles, ~15 GB) = the
fine-cell GPU-memory ceiling on the 16 GB A5000. Aggregation uses dt=1e-8 (cells loose →
stiff contact inactive) to cover the slow coalescence in feasible step counts.

Figures `figs/two_stage_n{100,200,400,600,800}.{png,_agg.mp4,_spread.mp4}` (4-panel:
loose seed → aggregated ball → spread with lamellipodia → V/V0·contact·Rg·footprint
diagnostics). Full report: `outputs/h_dcm_two_stage/REPORT.md`.

**⛔ Honest limit:** a DRAMATIC fully-dispersed suspension rounding into a dense spheroid is
a multi-day (12–48 h real) ≈ 50 h GPU process at the BAOAB per-step rate — the SAME wall
documented for spreading-magnitude + necrosis. So we seed near-spherical-but-out-of-range
and let the EMERGENT, negative-control-proven mechanism drive a real converged dense
aggregate in feasible wall-time, and state the timescale ceiling rather than fake it.

---

## Commits (branch h7/compartment-platform)
- `007fd07` cohesion fix (dead adhesion) + lamellipodium surface viz + arrest default-off
- `f7dd0d3` DCM tests under cohesion (fused `_footprint_radius` + arrest tests on legacy bands)
- `c171bf8` two-stage v1 (central-pull) — SUPERSEDED
- `17cc4b4` aggregation REDESIGN (active-matter coalescence, parity PASS)
- `1906502` emergence verified (negative control) + rounded-cell convergence
- `4d6e187`/`5f4b766`/`6a61e8a`/`5610fa0`/`3b21a89` two-stage N=100/200/400/600/800 figures

## Parallel session (after my API drop — for the record, NOT mine)
A parallel session is running a SPREADING dt-stability sweep (`dcm_spread_dt_sweep.py`,
reusing my N=100 aggregate): dt=1e-7 blows up (V/V0→1410×), dt=1e-8 deflates/collapses →
spreading needs dt=1e-9 (the spreading-magnitude timescale wall). It also fully removed the
arrest mechanism (deleted the script+test) — consistent with my "arrest is not real physics"
finding. These are that session's commits to make, not mine.

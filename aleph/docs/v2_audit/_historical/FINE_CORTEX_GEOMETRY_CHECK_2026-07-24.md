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

# Fine cortex config — CPU-weave geometry / percolation check (2026-07-24)

**Gate question.** Before spending a native GPU build on the derived FINE cortex config
(`cortex_seg_um=0.075`, `cortex_density_per_fil=40`), verify on the Mac CPU weave that it is a **valid,
single-spanning, ~75 nm physiological network** — not an under-connected one. Companion to
`CORTEX_MESH_FIDELITY_2026-07-24.md` §7.

**VERDICT — ✅ FINE CONFIG IS VALID. Build it.** At the full native population (70,686 filaments) the fine
network is **single-spanning** (one connected component, LCC = 100 % of nodes and 100 % of filaments, **zero**
isolated filaments), delivers the intended **~75 nm mesh** (contour-spacing median 75 nm), and its node/crosslink
counts match the memo exactly. **No density correction is needed** — percolation is held by `density_per_fil`
(filament-level degree ≈ 80 crosslinks/filament) and the mesoscale capture reach, **both independent of the mesh
resolution `seg_um`**, so refining the mesh does not threaten spanning. Robust across 3 RNG seeds.

## Method (read-only, no CUDA, no sim-code change)

Reused the viz weave-construction path — `ac.weave.woven_cell.weave_cell([cortex_region(...)],
overlap_free=False)` → `WovenCell.to_crosslinked_cortex()` (pure host NumPy; Warp imports only as a CPU dep and
is never launched). Script: `scratchpad/cortex_geom_check.py`. Ran **full native cell** (R_cortex = 7.40 µm,
`n_filaments = 70686` — the production FIXED value from `assemble.py:173`) for both configs, plus an R = 3 µm
full-areal-density proxy for cross-check and seed sweeps.

- **Connectivity:** node graph = actin backbone segments (consecutive nodes within a fiber) **+** α-actinin/
  filamin crosslinks (`xl_i–xl_j`); `scipy.sparse.csgraph.connected_components` (undirected). Filament-level
  component read at each fiber's first node (backbone makes all of a fiber's nodes share one label).
- **Mesh:** crosslink→crosslink **contour spacing** = gap (in segments) between consecutive crosslinked nodes on
  the same filament × `seg_um`. (Also a 3-D nearest-crosslinked-node distance, reported but secondary — see note.)
- **Degree:** crosslink endpoints per node and per filament; dangling = zero-crosslink nodes / filaments.

## Numbers — full native cell (n_filaments = 70,686)

| Metric | COARSE (seg 0.5, ρ20) | FINE (seg 0.075, ρ40) |
|---|---|---|
| actin nodes | **494,802** | **2,898,126** |
| crosslinks (`n_xl`) | **1,413,720** | **2,827,440** |
| beads / filament | 7 | 41 |
| **components (node graph)** | **1** | **1** |
| **LCC node fraction** | **1.000** | **1.000** |
| **LCC filament fraction** | **1.000** | **1.000** |
| isolated (singleton) filaments | **0** | **0** |
| **mesh contour spacing — median (IQR)** | **500 nm (500–500)** | **75 nm (75–75)** |
| mesh contour spacing — mean | 501 nm | 87 nm |
| crosslinks per **filament** (mean ± sd) | 40 ± 6.5 | 80 ± 9.7 |
| crosslinks per **node** (mean ± sd) | 5.71 ± 2.40 | 1.95 ± 1.41 |
| dangling **nodes** (0 crosslinks) | 1,683 (0.34 %) | 420,308 (14.5 %) |
| **dangling filaments** | **0** | **0** |
| NaN / self-loop / same-fiber / zero-rest bonds | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |

Sanity #4 ✅ — counts reproduce the memo (`494,802 / 1,413,720`; `2,898,126 / 2,827,440`) bit-exactly, no
degenerate bonds. Peak RSS for the fine full weave + analysis = **9.75 GB** (matches the ~10 GB estimate).

## Findings

1. **Percolation holds with huge headroom — the "under-connect" worry does not materialise.** Both configs are
   a single spanning component covering every filament, with **0** isolated filaments. The reason the finer mesh
   does *not* under-connect: crosslink count is set **per filament** (`n_xl = n_fil · density_per_fil`) and the
   capture reach is `mesoscale_reach = √(A_shell/n_fil) ≈ 0.10 µm` — **both independent of `seg_um`**. So every
   filament gets ~80 crosslink attachments (fine) / ~40 (coarse) regardless of mesh resolution, ~2 orders above
   the random-network spanning threshold (~1 link/filament). Confirmed single-spanning across seeds 0/1/2.

2. **What the 2× density actually buys is per-node connectivity, not spanning.** Crosslinks-per-node **drops**
   5.71 → 1.95 (fine) because nodes proliferate 5.86× while crosslinks only 2×. This is expected and harmless:
   the backbone rigidly chains all of a filament's nodes into one component, so filament-level spanning is
   insensitive to per-node degree. The 14.5 % "dangling" nodes in the fine mesh are interior actin beads that
   simply carry no crosslink — they are **not** disconnected (their filament is), and there are **zero** dangling
   filaments.

3. **The mesh is set by `seg_um`, confirming Issue B.** Contour spacing median = `seg_um` exactly in both
   (500 nm / 75 nm), because `density_per_fil` (20, 40) ≥ interior nodes (7, 41), so nearly every node carries a
   crosslink and adjacent crosslinked-node spacing collapses to one segment. The fine config **delivers the
   ~75 nm target**; the coarse ~500 nm is the discretization cap, exactly as the fidelity memo argued. The
   derivation `spacing = length/density` (150 nm coarse / 75 nm fine) is **not** what governs — the `seg_um` cap
   is; for fine they coincide only because ρ40 ≈ 41 nodes.

4. **Does holding percolation at seg 0.075 need ρ > 40? — No.** ρ40 spans with ~40× degree headroom; it would
   span even at ρ20 on the fine mesh. The ρ40 value is fixed by the **crosslink-spacing = mesh** requirement
   (Chugh 2018), *not* by percolation. There is **no density correction to derive** — the coarse ρ20 was a
   percolation proxy, but the fine ρ40 is already a sourced-spacing value that also percolates comfortably.

## Minor issue flagged (not blocking; no code changed)

The **production count uses a FIXED `n_filaments = 70686`** (`assemble.py:173`, = 100 µm⁻²·4π·7.5²), placed on
the R_cortex = **7.40 µm** shell. But the viz script's `--full-cell` path
(`scripts/ac_cortex_crosslink_viz.py:176`) *derives* `n_fil = round(100 · 4π · 7.40²) ≈ 68,806` while its inline
comment claims "= 70,686". That derivation (from R = 7.40, not 7.50) gives **68,806 filaments / 481,691 nodes**,
not the memo's 70,686 / 494,802 — a ~2.7 % undercount in the `--full-cell` viz only. It does **not** affect this
verdict (local statistics are identical; I ran the authoritative check at the production 70,686). Recommend the
Lead align the viz `--full-cell` path to the fixed `n_filaments=70686` (or read `CellConfig`) so its printed
counts match production.

*(Secondary metric note: the 3-D nearest-crosslinked-node distance (median ~17 nm coarse / ~8 nm fine) measures
the distance to a node's crosslink **partner** — i.e. the crosslink bond length / filament closest-approach —
not the pore size, so it is reported for completeness only; the contour spacing is the mesh metric.)*

## Bottom line

The fine config is a physically valid, single-spanning ~75 nm cortex network. **The native GPU build is worth
running.** Percolation is not at risk at ρ40 (nor would it be at lower ρ); the mesh target is met; counts and
bond integrity are clean.

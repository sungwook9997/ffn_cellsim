# The arena port plan — one 10M-node cell, on the device, measured

**Status:** PLAN. PI-requested 2026-08-20, session `d380041d`. Decides no constant, sets no threshold,
and selects no layout. Every number below is measured on this repository or stated by the PI in session;
where a figure is an extrapolation it says so on the same line.

**Supersedes** the staged port order proposed earlier the same day (cortex → membrane → turgor → ERM).
That order ranked pieces by porting difficulty, which three PI corrections retired.

---

## 0. What decided this plan

| PI statement, in session | What it retired |
|---|---|
| *"그 모든 런들을 아레나로 밀어서 사용하고 … 아레나에 이식하는 용도"* | the incumbent as a runtime |
| *"gpu only로 warp나 다른 방식으로 포팅"* | host-NumPy construction as the arena's build path |
| *"무조건 처음부터 네이티브로"* | develop-on-a-slice for this line |
| *"10M 먼저지 그렇게 풀 셀을 구현해야지 축소가 가능하니"* | multi-cell-first |
| *"측정한 것은 … 10m 노드 완벽하게 구현 안한 상황에서잖아"* | the 359x throughput figure |
| indentation ~3 min/cell, sweep needs ≈ real time | throughput as a nice-to-have |

---

## 1. The target, and what it is made of

`STATE.md` (b) native today is **551,434 nodes**.

⚠ **WHERE 10,000,000 CAME FROM, established 2026-08-20 after session A2 questioned it: it is a
CAPACITY number, not a derivation.** `arena_prestudy/w2_arena_capacity.json` carries
`"run_label": "w2_arena_capacity"` and `capacity: {nodes: 10000000, ...}` — it asked *"does a 10M-node
arena fit in 24 GiB"* and answered yes. **It never asked how many nodes a cell has.** This plan took
that figure as the target, computed §1's missing-compartment row as the subtraction residual
`10M − 5.13M`, and extrapolated §3's 6,514x throughput on top. Three layers rested on a capacity number.

**The budget derived so far, from A1 and A2 reporting densities x geometry:**

| population | nodes | basis |
|---|---:|---|
| cortex actin (A1) | 4,311,846 | 100/µm² x 4πR², seg 0.05 → 61 per filament |
| plasma membrane (A1) | 163,842 | subdiv 7; triangle edge 70.6 nm inside the sourced 50–100 nm band |
| nuclear envelope (A1) | 40,962 | ⚠ subdiv 6 **inferred** — lamin meshwork size is unsourced |
| microtubule (A2) | 74,001 | ⚠ per-cell count **unsourced** (proxy 250–600) |
| intermediate filament (A2) | 2,700 | ⚠ `n_fil` unsourced |
| filopodium (A2) | 61,000 | ⚠ per-cell count unsourced |
| lamellipodium (A2) | 4,200 | ⚠ areal density unsourced |
| stress fiber (A2) | 201,000 | ⚠ per-cell fibre count unsourced |
| **total** | **4,859,551** | **49% of the capacity figure** |

⚠ **Five of the eight rows are PI-GAPs, so 4.86M is not the answer either** — it is what the incumbent's
provisional counts give when discretised at the sourced segment length. Raising MT to the proxy ceiling
of 600, or putting filopodia and stress fibres at real per-cell numbers, moves it a lot. **The node
budget cannot be closed without one PI decision: the per-cell structure count for those five
populations.** That is now decision-queue item 1.

⚠ And the throughput arithmetic moves with it: at 4.86M the extrapolated step is 158 s and the factor
3,165x rather than 6,514x — which changes `k_xl`'s margin from 1.3x to 2.6x. Both are extrapolations
from a 551,434-node measurement and neither is a result. Measured from the
committed job-69 ledger, and the gap is not padding:

| | now | change | nodes |
|---|---:|---|---:|
| cortex actin | 494,802 | `seg_um` 0.5 → **0.05** (sourced; audit calls 0.5 *"5–10x TOO COARSE"*, provenance CONVENIENCE) | **4,311,846** |
| membrane | 642 (subdiv 3) | subdiv 7 | 163,842 |
| nucleus + rest | 16,312 | — | 16,312 |
| **subtotal** | | | **≈ 4.53 M** |
| MT, IF, filopodium, lamellipodium, SF | **0** | populations that do not exist yet | **the other ≈ 5.5 M** |

### 1.1 ⚠ The table above has no cytosol row, and the missing row is the largest allocation

Added 2026-08-20 after session E established that the cytosol is a **field, not a node population**
(`CYTOSOL_ARENA_REPRESENTATION_2026-08-20.md`). Every row above counts NODES. The fluid grid holds no
nodes, so it was never subtracted from the capacity figure and never entered §3's throughput
arithmetic — and it is plausibly the single biggest thing in the cell.

Reproduce with `python -m aleph.world.build.cytosol`. Box derived as `R_cell` + one OUTSIDE layer, so
these are cells the build would actually claim:

| dx [µm] | shape | cells | GB | CFL dt [s] | subcycles/step | relative fluid work |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 33³ | 35,937 | 0.002 | 7.50e-04 | 67 | 1× |
| 0.25 | 63³ | 250,047 | 0.015 | 1.88e-04 | 267 | 28× |
| 0.1 | 153³ | 3,581,577 | 0.215 | 3.00e-05 | 1,667 | **2,480×** |
| 0.05 | 303³ | 27,818,127 | 1.669 | 7.50e-06 | 6,667 | **77,027×** |

⚠ **CORRECTED against E's first draft**, which fixed a 16 µm cube and gave 33.1M cells / 1.985 GB /
91,586× at dx = 0.05. The conclusion is unchanged and the numbers above are the ones to quote. The lead
session quoted the draft figures once before this correction landed.

**Cost scales as dx⁻⁵** — cell count dx⁻³ times CFL subcycles dx⁻², from
`field_grid.py:126-128`: `safety · S · dx² / (2 · d · mobility)`.

**Two things follow, and neither had been priced:**

1. At the cortex's sourced 0.05 µm the grid alone is **27.8 M cells — 5.7× the entire 4.86 M solid
   budget**, and ~3.9× its bytes (solid ≈ 88 B/node ≈ 0.43 GB). The capacity question this plan
   inherited asked whether 10 M *nodes* fit. It never asked about a field.
2. It runs **against §3's only lever.** `k_xl` corrected into band buys ~6,500×; moving the grid from
   0.5 to 0.05 spends ~77,000×. **`dx_um` and `k_xl` have never appeared on the same page**, and the
   throughput target (AFM indentation ≈ 3 min/cell, sweep ≈ real time) is decided by their product.

**`dx_um` is a declared axis with no default; the builder refuses without it.** Band is sourced —
lower bound pore size ξ ≈ 14–40 nm (`KB-DRAFT-3.B-27`, `-3.B-18`, both draft), upper bound
√(D_p·dt_phys) = 1.58 µm (`KB-3.B3.2`, verified). The incumbent's 0.5 is inside the band but was not
chosen from it: `assemble.py` labels it *"numerical sizing; not literature"* — the same class of value
as the cortex's `seg_um = 0.5`. **Which value inside the band is a PI decision, not a derivation.**

**Nucleoplasm is ONE grid, not two** — settled 2026-08-20 by sessions F and E together, so `dx_um` goes
to the PI as a single number. The quantity that would split them, `ξ_nuc`, is absent from the corpus,
and `Moeendarbary2013_NatMater` p9 states its D_p was measured *"avoiding the nucleus"*. Nucleoplasm
becomes mask-indexed material properties on the shared grid — but not for free: `biot_pmass_update_kernel`
holds every non-FLUID cell, so the nucleus today is a **no-flux hole, not a pressure DOF**, and must be
promoted first. Reopen if `ξ_nuc` or `D_p,nuc` arrives; ⚠ note the sign is the opposite of the intuition
that raised it — lamin 0.4 µm and chromatin 0.6 µm against cortex 0.05 µm make the nucleus **8–12×
coarser**, so a separate grid would be justified to run the nucleus *coarser*, not finer.

⚠ **CORRECTED 2026-08-20 by session A1, and the correction is the point.** This table first read
4,948,020 for the cortex, computed as 494,802 x 10 — which assumes node count is inversely proportional
to segment length. It is not: **N = L/seg + 1**, and the +1 is per filament. Checked against the
committed ledger: 494,802 / 70,686 = 7.0 nodes per filament, and 3.0/0.5 + 1 = 7 exactly. At seg 0.05
the figure is 3.0/0.05 + 1 = 61 per filament, so **4,311,846**. The row is the derived number now.
A1 derived it from density x geometry as the plan asks; I had scaled a ratio.

⚠ **Half the target is compartments that have no population at all.** Reaching 10M is therefore not
"multiply the cortex"; it is building the compartments whose connectors are declared and whose runtime
objects do not exist. Sizing them from physiological density is work item **W2** below, not an
assumption here.

⚠ And the segment-length fix moves `kmax` as well as the node count: at the sourced length the backbone
is ~5.9e5 pN/µm and the crosslink ratio falls from 9.3x to 1.4x. **Node count and timestep move in
opposite directions**, which is why no throughput target can be settled before the cell is built.

---

## 2. The connector census, cell-only

From the committed device-run census (`connector_devicerun/native_record.json`, 4090-1):

| | count |
|---|---:|
| declared connectors | **38** |
| external (ECM / substrate / medium) — **out of scope for this plan** | 8 |
| **cell-internal** | **30** |
| of those, holding a runtime object | **2** (`nmii_cortex_motor`, `dorsal_arc_crosslink`) |
| `NOT_BOUND` — weaker than a stub, no object exists | **28** |

The 8 excluded: `ecm_crosslink`, `ecm_far_field_anchor`, `integrin_collagen_clutch`,
`membrane_ecm_contact`, `membrane_medium_traction`, `fa_actin_anchor`,
`filopodium_nascent_clutch`, `lamellipodium_nascent_clutch`.

**So 28 of 30 are not a port — they are a first implementation.** That is the honest shape of the work,
and it is also why the arena's family model matters: the declaration list contains duplicates that are
one physics each.

| declared as | count | arena family (estimate, to be confirmed one by one) |
|---|---:|---:|
| NMII motor (cortex / sf / filopodium / lamellipodium) | 4 | 1 |
| cytosol transfer (6 compartments) | 6 | 1 |
| LINC (mt / if / actin_cap) | 3 | 1 |
| membrane–cortex (contact + erm) | 2 | 1–2 |
| the rest, individually | 15 | 15 |
| **total** | **30** | **≈ 19–20** |

⚠ That folding is read off connector NAMES, not off their physics. Confirming it is work item **W1**.

---

## 3. The throughput requirement, and the only lever

Measured (job 69, 4090-1, 551,434 nodes, 4,000 inner, `dt_phys` 0.05 s): one outer step = **17.96 s**,
real-time factor **2.784e-03**. Linearly extrapolated to 10M nodes (the dominant kernels are
bandwidth-bound): **~326 s/step**, factor **1.535e-04 — 6,514x slower than real time.**

Against the PI's protocol — indentation ≈ 180 s simulated, sweep needs ≈ real time:

| | one point | 100-point sweep |
|---|---:|---:|
| today, extrapolated to 10M | **13.6 days** | **1,357 days** |
| real time | 3 min | 5.0 h |

**Required: ~6,500x.**

`k_xl` is the only known lever of that size: `dt_mu = 0.1/kmax` is crosslink-dominated at 8.2e5 pN/µm,
a ~2 Å binding-barrier curvature standing in for a 35–160 nm protein, and moving into the field band
(0.1–100 pN/µm) opens `dt` by 8,200x–8.2e6x. **The low end clears 6,514x by 1.3x** — and only if the
crosslink still dominates `kmax` in the band, which the built cell will decide. Solver optimisation is
not a lever: measured ceiling ~9%, recorded *"do not retry"*.

---

## 4. The phases

⚠ **The first version of this section was incoherent and the PI caught it.** It had W0 as *"a native
population world stands in device memory"* while putting the populations themselves in W2 — and §1 of
this same document says half the 10M is compartments that do not exist yet. **W0 could not have reached
native.** The corrected shape is the PI's: *"10m 세팅 구현되면 거기서 각각 의미 부여 하면서 연결하는
방향"* — stand the whole 10M setting up first, then give each part its meaning, then connect them.

That ordering is also the arena's own: `arena.py` allocates and claims, `strand`/`surface` give
structure, `bond` gives connection, and `laws/` gives meaning. Phases follow the layers.

---

### PHASE 1 — the 10M setting stands. Geometry only, no meaning.

**Every population, at native density x geometry, on the device.** Cortex at the sourced segment length,
membrane and envelope at their subdivisions, and the compartments that today have zero nodes — MT, IF,
filopodium, lamellipodium, SF — all claimed and built. Nodes, segments, angles, faces, hinges. No force,
no law, no parameter that is not a geometric one.

* **Proves:** ~10M nodes stand, `assert_partitioned` passes, and the count comes from **densities x
  geometry** rather than from the target. If the densities give 6M or 14M, that is the answer and 10M
  was the estimate.
* **Measures, and this is why it is first:** `exact_peak_gpu_bytes` — **a number no session has ever
  taken** — build wall time, and **bytes per node**. `arena.py` defers its solver-workspace layout
  because *"guessing would bake an answer into the layout"* and the measurement has not been taken.
  This is that measurement, and it cannot be taken on a partial cell.
* **May NOT claim:** any physics, any throughput of a STEP (nothing steps yet), and that the geometry is
  right because it built.
* **Design questions it must ANSWER rather than assume:** claims per population or per strand
  (`build_strand` takes 3 per strand today, so native cortex is 212,058 claims); and whether
  construction is a Warp kernel or a single host build uploaded once and never mirrored again.
* ⚠ **Where the PI-GAPs surface.** Every density is a declared axis under standing ruling #1, and the
  build refuses rather than defaulting. Expect refusals; they are the mechanism working, not a blocker.

### PHASE 2 — meaning. One population at a time.

Each population is told **what it is**: which `laws/` kernel reads it, what parameters that law carries,
with provenance, and its own render. `laws/` is bound, not ported — proven at `78942fc4`, where two
kernels ran over arena-addressed arrays unchanged at 1.108e-16 with exactly 0.0 pN outside the claim.
Same A/B per law.

* **Proves, per population:** its force channel is force-identical to the incumbent's on the same
  geometry, or the difference is explained.
* **Measures:** the marginal step cost of each law — which turns §3's single extrapolated number into a
  per-law budget, and makes the layout arithmetic instead of taste.
* ⚠ **Not everything is a bind.** Turgor traction and the ERM force live in `components/incumbent/`, and
  `∂V/∂x` exists in **three copies** — `components/nucleus/geometry.py:115` (host NumPy),
  `components/incumbent/membrane_pressure.py:48` (the only device version) and
  `resting_balance_oracle.py:94` (an oracle-local copy). Those are a rewrite and a **fold into `laws/`**.

### PHASE 3 — connection. The 30 become families.

Only now do populations couple. Each of the 30 cell-internal declarations becomes an arena `BondFamily`
or is shown to be a duplicate of one — **by physics, not by name**, confirming or refuting §2's estimate
of ~19–20.

* **Proves:** every family answers `BondCount`'s four questions — *how many, against what, for which
  cell, on whose authority* — or it is not built.
* **Why this is not bookkeeping:** the incumbent's ERM was ONE declared connector standing for a
  POPULATION, so **nobody had to answer how many**, and the answer turned out to be the icosphere vertex
  count. All three of this session's L0 runs measured that number and not a density.
* **28 of the 30 have no runtime object at all**, so this phase is mostly first implementation. The plan
  says so rather than calling it a port.

### PHASE 4 — the step, and what accepts it.

Forces accumulate, the solver runs, and a step is accepted or rejected.

* ⚠ **The one thing the arena must not inherit.** `ENGINE_FORWARD_ACCEPTANCE_2026-08-20.md` records that
  the incumbent accepts a step when two force channels cancel to within float64 roundoff — a predicate
  that cannot fail — while `inner_converged` is computed and read by no gate. The acceptance predicate
  here is written against the STATE, and it ships with a demonstrated negative control before it gates
  anything.
* **Only here** does the throughput requirement of §3 become testable, and only here does `kmax` get
  measured on the real cell — which is what decides whether `k_xl`'s 1.3x margin is real.

## 5. What this plan does not do

It selects no constant, sets no threshold, and picks no layout — W0's measurement is what makes the
layout decidable. It moves no tier-(a) row: `engine/` and `components/` own every native number until
the arena re-earns them. It does not touch the 8 external connectors. And it does not schedule the
`k_xl` decision, which is the PI's and now the highest-priority one, because W0's `kmax` measurement is
what tells them whether the 1.3x margin is real.

## 6. Carried forward, unresolved

1. **`k_xl`** — priority 1. The only lever large enough, with a 1.3x margin that is a measurement.
2. **Inter-cell ECM** — deferred with multi-cell, but it decides the claim model when it arrives.
3. **Multi-GPU** — deferred: one 10M cell is ~0.7 GB. It returns at ~29 cells, the one-card ceiling.
4. **`STRUCTURE.md`** — PI-authored, read-only, and now silent on `aleph/world/` entirely.
5. **Branch name** — `engine/main` while `aleph/engine/` is frozen port source.

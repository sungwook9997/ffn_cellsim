# The arena's real constraint is throughput, not memory — and the target is now a number

**Status:** MEASUREMENT + REQUIREMENT. Decides no layout. Written 2026-08-20, session `d380041d`, from
the PI's corrections in session: the target is not one cell, and *"동적인 시뮬레이션인 만큼 시뮬레이션이
빨리 돌아가야지 실제 시간과 시뮬레이션 타임 사이 차이가 적은데, 그걸 빼면 안되지"*.

Everything below is measured on this repository or stated by the PI. Nothing is extrapolated from
literature and no layout decision is taken.

---

## 1. What was measured, and on what

Full-native single cell, RTX 4090-1, Slurm job 69, 2026-08-20 — the L0 control point:

| | |
|---|---|
| population | 70,686 filaments / 494,802 actin / **551,434 total nodes** |
| inner budget | 4,000 iterations (matching the committed record at `4c165576`) |
| `dt_phys` | 0.05 s |
| outer step wall | **17.96 s** |
| build wall | 23.7 s (once, not per step) |

    real-time factor = 0.05 / 17.96 = 2.784e-03      i.e. 359x SLOWER than real time

⚠ Scope: this is ONE outer step of a rejected configuration on the PORT SOURCE engine, not a converged
trajectory. It bounds throughput; it says nothing about physics. The `STATE.md` (b) row measuring
*"6/11 components step in 243 ms; real-time ceiling ~0.64"* is a different, partial configuration and
is not the same number.

---

## 2. The requirement, from the PI's protocol

PI, in session: cortical tension is measured by **indenting one cell for about 3 minutes**, and a
parameter sweep needs throughput *"적어도 실제 시간과 비슷한 정도"* — at least near real time.

One AFM protocol = **180 s of simulated time**. Against the measured factor:

| throughput | one point | 10-point sweep | 100-point sweep |
|---|---:|---:|---:|
| **today** (2.78e-03) | **18.0 h** | **7 days** | **75 days** |
| x10 | 1.8 h | 18.0 h | 7 days |
| x100 | 0.2 h | 1.8 h | 18.0 h |
| **x359 = real time** | **3 min** | **0.5 h** | **5.0 h** |

**At the population measured, the requirement is ~359x.**

⚠ **BUT THAT POPULATION IS NOT THE CELL.** PI, correcting this document in session: *"측정한 것은 결국에는
10m 노드 완벽하게 구현 안한 상황에서잖아"* — and it is right. 551,434 nodes is **5.5% of the 10M-node
target**. The dominant kernels are bandwidth-bound, so step time is close to linear in nodes:

| | nodes | step | real-time factor | slower than real time |
|---|---:|---:|---:|---:|
| **measured**, job 69 | 551,434 | 17.96 s | 2.784e-03 | **359x** |
| **10M target**, linear extrapolation | 10,000,000 | **~326 s** | 1.535e-04 | **6,514x** |

Against the 3-minute protocol on a real 10M cell:

| throughput | one point | 100-point sweep |
|---|---:|---:|
| today, extrapolated | **13.6 days** | **1,357 days** |
| x100 | 3.3 h | 14 days |
| x1000 | 20 min | 1.4 days |
| real time | 3 min | 5.0 h |

**So the requirement is ~6,500x, not 359x**, and the earlier figure described a cell that does not exist.

⚠ And the missing 18x is not padding — it is physiology the model currently omits. The provenance audit
records `cortex_seg_um_deviation = "5–10x TOO COARSE (0.5 µm vs 0.05–0.10 µm)"`, provenance CONVENIENCE.
**Segment length alone accounts for most of the gap**, and the same document notes it also makes the
backbone artificially soft relative to the crosslink (ratio 1.4x instead of 9.3x at the sourced length).
Building the full cell therefore changes `kmax` as well as the node count, and the two move in opposite
directions — which is why the throughput target cannot be settled before the cell is built.

**This is also the PI's argument for building 10M FIRST:** *"그렇게 풀 셀을 구현해야지 축소가 가능하니"* —
a reduction is only defensible once there is a full thing to reduce FROM, and every number above shows
what happens when a partial configuration is quoted as if it were the target.

---

## 3. Memory says something different, and it is the less binding of the two

At ~80 B/node (position + snapshot + force + mobility, float64):

| configuration | nodes | GPU bytes |
|---|---:|---:|
| one cell, today's cortex | 551,434 | 0.04 GB |
| one cell, arena target | 10,000,000 | 0.7 GB |
| 10-cell spheroid | 100,000,000 | 7.5 GB |
| **29 cells — one 24 GiB card** | 290,000,000 | ~22 GB |
| 100 cells | 1,000,000,000 | 74.5 GB |
| 400 cells (DCM aggregation scale) | 4,000,000,000 | 298.0 GB |

Hardware on the box: 2 x RTX 4090 (24 GiB) + 1 x RTX 3090 (24 GiB) = 72 GiB, so ~88 cells of 10M nodes
before ECM. ⚠ **Inter-cell ECM is NOT in this table.** The only recorded ECM population is 200–220
collagen fibers with **zero crosslink bonds**; the matrix filling a spheroid is not sized anywhere.

⚠ **The two ceilings are nowhere near each other.** Memory allows 29 cells on one card. Those 29 cells
need **2 days of wall clock per simulated second** at today's factor. Spheroid aggregation is a
minutes-to-hours phenomenon, so the scale memory permits is one that throughput refused long before.
**Layout must therefore be chosen against step time, not against capacity** — and the first version of
this analysis presented only the memory table, which is the less binding axis.

---

## 4. Where the 359x could come from

**`k_xl` — the single largest known lever, and it is a PENDING PI DECISION.**

    dt_mu = 0.1 / kmax = 3.54e-08 s,  and kmax is crosslink-dominated at k_xl = 8.2e5 pN/um

`CROSSLINK_STIFFNESS_AXIS_2026-08-16.md` establishes that 8.2e5 is a **binding-barrier curvature over
~2 Å** used as the spring constant of a 35–160 nm protein, and that the field band is 0.1–100 pN/µm.
Moving into that band opens `dt` by **8,200x to 8.2e6x**. Against the corrected requirement of 6,514x
the **low end clears it by only 1.3x** — where against the retired 359x figure it looked like a 23x
margin. Two caveats now matter that did not before:

* it holds only **if the crosslink still dominates `kmax` in the band**, which is not yet known;
* building the cell to 10M **also changes `kmax`**, because the missing nodes come mostly from a segment
  length the audit calls 5–10x too coarse, and the sourced length raises the backbone stiffness.

**So `k_xl` alone is no longer comfortably sufficient**, and the margin is a measurement rather than an
assumption.

**What is NOT a lever:** solver optimisation. Measured ceiling ~9%, recorded as *"do not retry"*.

**What is unmeasured and layout-decides:** the arena's dominant kernels are **bandwidth-bound**
(`arena.py`'s own rationale for the live-prefix discipline). Throughput is then close to linear in bytes
touched per node per step, which makes these arithmetic rather than taste once measured:

* fields per node — today's 80 B assumes position + snapshot + force + mobility are all touched;
* float64 vs float32 — a factor of 2, and this engine is float64 throughout;
* whether the snapshot array is written every step or only on an accepted one;
* the inter-cell neighbour structure, which does not exist yet and may dominate in a spheroid.

⚠ `arena.py` deliberately does not allocate solver workspace, on the stated grounds that *"whether the
inner solve needs Krylov vectors, per-node clock state, or velocities depends on a measurement that has
not been taken yet, and guessing would bake an answer into the layout."* That is the right posture, and
**this document is the argument that the measurement is now the blocking one**: bytes per node per step
IS throughput, and throughput is whether the sweep exists.

---

## 5. What this changes about the port plan

The plan as of this morning ordered the port by porting difficulty — cortex first, because two of its
kernels already had arena parity. Three PI corrections in sequence retired that ordering:

1. **The world implementation must be GPU-only.** Four of the five arena modules are `CPU-importable`
   host NumPy today; only `arena.py` touches a device, and only to allocate.
2. **Start at native, never subdiv 3.** Vindicated by this session's own failure: all three L0 runs used
   membrane subdiv 3, where the ERM count (642) is the icosphere vertex count and not a density —
   `bond.py` already names that class of defect as *"a mesh number wearing a physiological label"*.
3. **The target is many cells plus inter-cell ECM**, which is what makes host-side assembly impossible
   rather than merely slow: at 10 cells there is no room to hold a host copy.

So the first piece is not a force port. It is: **a native-scale world standing on the device, with its
step time and its bytes-per-node measured.** Force porting is what that measurement makes decidable.

---

## 6. Open to the PI, in the order they now block

1. ~~First target scale~~ — **DECIDED 2026-08-20 by the PI: one 10M-node cell first**, then N.
   *"10M 먼저지 그렇게 풀 셀을 구현해야지 축소가 가능하니."* A reduction is only defensible from a full
   thing, and this document's own retracted 359x is the worked example of quoting a partial one.
2. **Inter-cell ECM** — a population in the same arena, or its own? Cells move, so the ECM bonds rewire,
   and that decides the claim model.
3. **Multi-GPU in the design now?** 29 cells is the one-card memory ceiling and 100–400 is inside the
   stated target, so partitioning cannot be added later without doing exactly what `arena.py` warns of.
4. **`k_xl`** — unchanged from 2026-08-16, and this document raises its priority: it is the only known
   lever large enough to reach the required 359x.

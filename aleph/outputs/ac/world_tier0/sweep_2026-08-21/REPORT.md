# The NMII shell has two constraints, and they had never been on one axis

The motors are placed out of reach of the cortical actin they exist to pull. This run measures that
at full native instead of inferring it, and finds that the recommended repair does not work either —
for the opposite reason, which no single run could have shown.

## ⚠ What is NOT claimed here, and why the headline carries no number

TIER 0 is existence and sign. Every distance below is a geometric property of the BUILT cell with no
force evaluated, in the same class as *"every node lies inside the built membrane"*. The capture
radius drawn as a reference band is the `hand_kmc` NMIIA proxy and is a **PI-GAP**;
`n_heads_per_side` is **UNRATIFIED** (PI queue item 2); the support area is PI queue item 19's
membrane sphere, held FIXED at its wrong value across every point so the ladder has one cause.

⚠ **Nothing was moved and no radius was chosen.** Picking a value off these curves is choosing a
constant so a gate passes unless it is derived, and the derivation is the PI's ruling on PI queue
item 14.

## What was run

| | |
|---|---|
| driver | `aleph/scripts/world_tier0_stations.py`, `aleph/scripts/world_phase1_native.py` |
| device | RTX 4090-1 `cuda:0`, inside Slurm job 73 via `~/bin/ffn-run` |
| build | `13233392`, hand-deployed; closure 22/22 hash-matched against the host **before** launch |
| population | 442 minifilaments / 32,708 NMII nodes, in the full native cell (11 populations + cytosol field) |
| axis | NMII shell radius, ladder `{6.95 … 7.40}` step 0.05 µm, **declared before the run and not narrowed after** |
| held fixed | count, support 706.86 µm², `n_bb`, head offset, `n_heads_per_side` = 30, shell thickness |

## The measurement

| shell radius [µm] | nearest cortex node, min / median / p95 [µm] | stations @ reach 0.20 µm | outermost NMII node [µm] | `assert_inside_membrane` |
|---|---|---|---|---|
| 6.95 (as built) | 0.2501 / 0.3532 / 0.4427 | 0 / 442 | 7.0536 | PASS |
| 7.00 | 0.2004 / 0.3034 / 0.3927 | 0 / 442 | 7.1036 | PASS |
| 7.05 | 0.1509 / 0.2540 / 0.3432 | 80 / 442 | 7.1535 | PASS |
| 7.10 | 0.1012 / 0.2046 / 0.2939 | 191 / 442 | 7.2035 | PASS |
| 7.15 | 0.0514 / 0.1553 / 0.2447 | 311 / 442 | 7.2535 | PASS |
| 7.20 | 0.0058 / 0.1061 / 0.1951 | 422 / 442 | 7.3035 | PASS (run) |
| 7.25 | 0.0040 / 0.0589 / 0.1455 | 442 / 442 | 7.3534 | not run |
| 7.30 | 0.0027 / 0.0250 / 0.0970 | 442 / 442 | 7.4034 | not run |
| 7.35 | 0.0017 / 0.0195 / 0.0498 | 442 / 442 | 7.4534 | PASS (run) |
| 7.40 (queue's option a) | 0.0013 / 0.0181 / 0.0294 | 442 / 442 | 7.5033 | **REFUSED — 414/32,708 nodes out** |

Figure: [`figs/reach_vs_containment.png`](figs/reach_vs_containment.png), rendered by
[`fig_reach_vs_containment.py`](fig_reach_vs_containment.py) from the records beside it.

## Three findings, in the order they matter

**1. Item 14's separation is confirmed by measurement, not by arithmetic.** The prior record
(`world_tier0/stations.json`, build `90a3b453`) predates the census gaining `nearest_cortex_node_um`
and carries it as `null`, so the separation had only ever been inferred from the station curve. It is
now read directly: at the as-built shell the CLOSEST minifilament in the whole population is farther
from actin than the capture proxy. Not "most cannot bind" — none can, and no rate constant reaches
across a gap.

**2. ⚠ The queue's option (a) is wrong in SIGN, and this is the correction.** Item 14 reasoned that
NMII at the cortical shell sits with its heads a few nm INSIDE the membrane. At that radius
`assert_inside_membrane` refuses the build with 414 nodes a few nm OUTSIDE it. The arithmetic and the
build disagree, and the build wins.

**3. ⚠ The excursion guard cannot fail for the thing it appears to protect.** `build/nmii.py:254`
raises when `excursion > thickness/2` — a FIT test, *"does a minifilament fit inside the shell"*.
But `nmii.py:156` draws the per-filament radius over the FULL thickness,
`r = radius + thickness * (randf - 0.5)`, with no inset for that excursion. A filament drawn at the
outer edge therefore puts its beads that excursion PAST the shell whatever the guard returns. The
measured overshoot is CONSTANT across the whole ladder and is SMALLER than the guard's own predicted
excursion — so the guard's number was never the problem, its question was. This is the family
`geometry.py`'s own self-check names when it points at `balance_ok` and `descent_ratio`.

## What this does not settle, and who settles it

The two constraints leave a window rather than a value, and the window's endpoints are themselves
resolved only to the ladder's step. **Choosing inside it is the PI's ruling on item 14**, and the
options are not equivalent:

* derive the shell from the cortex it must reach and subtract a clearance that is itself derived —
  the same shape as `tip_clearance_um`, which exists because a structure placed ON a boundary has no
  margin for the spread around it;
* or fix the placement draw so the shell means what the guard assumes, and re-ask the question;
* or declare this family non-contractile in the resting archetype, which is a physics claim owed its
  own justification.

⚠ The third finding is upstream of all three: while the guard tests fit and the draw ignores it, any
radius chosen from this table is chosen against a shell whose outer edge is not where the build puts
its nodes.

⚠ **`world/geometry.py` and `world/build/nmii.py` still hold two independent copies of this
placement**, and a third now reads it in each driver's `--nmii-radius-um`. That is deliberate and
temporary: one owner for the shell is the fix, and creating it is a change to the cell.

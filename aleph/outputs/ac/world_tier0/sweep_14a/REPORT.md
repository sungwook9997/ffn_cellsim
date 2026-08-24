# After the inset: the guard promises something true, and the window moved

The excursion guard tested FIT while the kernel drew a filament radius from the full shell thickness,
so a filament at the outer edge sat that excursion past the shell whatever the guard returned. The PI
ruled the ordering: fix the guard first, then the radius. This is the re-measure after the fix, on the
same declared ladder, extended one rung so the new refusal has a point on it.

## ⚠ What is NOT claimed here

TIER 0 is existence and sign. Every distance is a geometric property of the BUILT cell with no force
evaluated. The capture radius drawn as a reference band is the `hand_kmc` NMIIA proxy and is a
**PI-GAP**; `n_heads_per_side` is **UNRATIFIED**; the support area is queue 19's membrane sphere, held
FIXED across every point so the ladder has one cause. **No radius is chosen here** — that half of the
ruling is the PI's.

## What was run

| | |
|---|---|
| driver | `aleph/scripts/world_tier0_stations.py` |
| device | RTX 4090-1 `cuda:0`, Slurm job 74, via `~/bin/ffn-run` |
| build | `cffbfcf0`, hand-deployed; closure 22/22 hash-matched against the host **before** launch |
| population | 442 minifilaments / 32,708 NMII nodes, full native cell |
| axis | shell radius from 6.95 to 7.45, step 0.05, declared before the run |

## The measurement

| radius | nominal shell outer | measured outermost node | overshoot [nm] | nodes outside R | nearest actin min / p95 | stations @ 0.20 |
|---|---|---|---|---|---|---|
| 6.95 | 7.050 | 7.0491 | **−0.86** | 0 | 0.2546 / 0.4385 | 0 / 442 |
| 7.00 | 7.100 | 7.0991 | −0.86 | 0 | 0.2048 / 0.3886 | 0 / 442 |
| 7.05 | 7.150 | 7.1491 | −0.86 | 0 | 0.1553 / 0.3392 | 78 / 442 |
| 7.10 | 7.200 | 7.1991 | −0.86 | 0 | 0.1056 / 0.2900 | 191 / 442 |
| 7.15 | 7.250 | 7.2491 | −0.85 | 0 | 0.0557 / 0.2407 | 318 / 442 |
| 7.20 | 7.300 | 7.2991 | −0.85 | 0 | 0.0078 / 0.1911 | 428 / 442 |
| 7.25 | 7.350 | 7.3491 | −0.85 | 0 | 0.0025 / 0.1420 | **442 / 442** |
| 7.30 | 7.400 | 7.3991 | −0.85 | 0 | 0.0019 / 0.0931 | 442 / 442 |
| 7.35 | 7.450 | 7.4491 | −0.85 | 0 | 0.0030 / 0.0461 | 442 / 442 |
| **7.40** | **7.500** | **7.4991** | **−0.85** | **0** | 0.0028 / 0.0289 | 442 / 442 |
| 7.45 | 7.550 | 7.5491 | −0.85 | **7,678** | 0.0028 / 0.0544 | 442 / 442 |

Figure: [`figs/reach_vs_containment.png`](figs/reach_vs_containment.png).

## Three readings

**1. The guard's promise is now true, with margin.** The overshoot was a constant **+3.3 nm** at every
radius before the fix and is a constant **−0.85 nm** after — the outermost node now sits *inside* the
nominal shell. The sign flipped and the residue is small, which is what insetting by a conservative
formula should produce: the guard predicts 4.2–4.5 nm and the real excursion is less, so the inset
over-reserves by about 0.85 nm.

**2. ⚠ The admissible window MOVED because of the fix, not because of a new choice.** At 7.40
`assert_inside_membrane` refused the build with 414 of 32,708 nodes outside; it now passes with zero.
The refusal moved out to 7.45, where 7.45 + 0.10 = 7.55 exceeds the membrane and 7,678 nodes sit
49.1 nm past it — a refusal that is arithmetic rather than an artifact.

**3. ⚠ And the radius stopped being a choice inside a window.** With the fix, 7.40 is exactly the
cortical shell:

```
cortex   radius 7.4, thickness 0.2  ->  shell [7.300, 7.500]
NMII  at radius 7.4, thickness 0.2  ->  shell [7.300, 7.500]      identical
```

That is a **derivation** — *NMII lives in the cortical shell* — and not a value picked off a table to
make a gate pass. It is also the value `build/nmii.py`'s own self-check has always used. Before the
guard was fixed that derived answer was refused; the ordering the PI ruled is what made it available.

⚠ **What this does NOT settle.** `geometry.py` still derives the NMII radius as
`r_cell − tip_clearance − 0.30`, giving 6.95, where the 0.30 charges a tangential head arm as a radial
one. Making the cortical shell the built default is an edit to that line and is a change to the cell.
**The measurement is here; the ruling is not.**

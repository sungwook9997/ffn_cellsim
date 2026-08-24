# The τ run III figures — written beside their records, 2026-08-21

Produced by `aleph/scripts/world_tau_figures.py`. **Regenerate, do not hand-edit**, and re-run it when
seed 3 lands so all three seeds are on the same axes:

```
python aleph/scripts/world_tau_figures.py aleph/outputs/ac/world_phase4/tau3_seed{1,2,3}.json
```

⚠ **These plot; they do not judge.** Every stationarity verdict belongs to
`aleph/observe/stationarity.py` and to the records' own fields. **A refused series is drawn exactly
like an accepted one** — seed 2 does not reach stationarity inside its run and its trace is here in
full, because hiding a refused series is how a refusal stops being visible.

## `tau3_seed{1,2}_gamma_trace.png`

Every one of the 30,000 samples, no decimation (`feedback-viewer-no-downsample`). The abscissa is the
**sample index and not seconds** — the driver's `dt_s` is 1.0, so a sample is a step and labelling the
axis "s" would put a physical time on a clock that has none. That warning is on the axis rather than in
a caption, because a cropped panel loses its caption.

y runs from **0**, so nothing is truncated. The dashed line is the mean of the last 2,000 samples — a
description of where the series ended up, **not a definition of equilibrium and not in any contract**.
The dotted verticals are the three cuts characterised across seeds.

## `tau3_three_gates.png`

Same cuts and same columns for every seed, so a difference between lines is a difference between seeds
and not between treatments.

* **τ over the kept tail.** Seed 1 is flat (102 → 116 → 112); seed 2 falls monotonically without
  converging (971 → 556 → 138), which a settled series cannot do.
* **n_eff = kept / 2τ — the gate that decides.** The red band is the module's bar of 25, drawn as a
  reference band and **never recomputed here**: if it moves in `stationarity.py`, this figure is wrong
  and must be regenerated. ⚠ Seed 1 crosses **into** the band between 22,000 and 26,000, so its
  4,000-sample window is below the bar too — the only cut where the two seeds' `opens` verdicts agree
  is a cut where **neither has measured anything**.
* **Tail slope.** The axis is signed, not |·|, because the sign is the entire seed-1/seed-2 difference:
  seed 1 wanders across zero (+1.3e-04, −2.7e-04, +1.0e-04) while seed 2 is uniformly positive and
  **increasing as the window shortens** (+1.0e-03 → +1.4e-03).

## What was fixed by looking at them

The first render put five five-digit ticks on top of each other — `190002000021000220002300024000`,
reading as one number. **Found by opening the png**, which is the only way that class of defect is ever
found, and the same lesson as the app's status bar printing nothing for as long as it ran.

## What these figures are not

Nothing here is a criterion. The plateau, the percentile timings and the slope do not appear in any
contract and may not be used to accept or reject a run. The declared read order for contract III lives
in `aleph/scripts/world_tau_read.py` and in the contract itself.

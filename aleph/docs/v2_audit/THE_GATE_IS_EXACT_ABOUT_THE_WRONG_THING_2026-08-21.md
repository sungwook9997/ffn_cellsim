# The tolerance is exact to 0.000% and the residual it bounds swings by 1000% — and my first table compared two different cells

**2026-08-21, ~05:20 KST.** Supersedes §8 of
[`BALANCE_GATE_MEASURED_2026-08-21.md`](BALANCE_GATE_MEASURED_2026-08-21.md), written ninety minutes
earlier from **one run per configuration**.

## 1. The correction first

§8 tabulated *"membrane only / + cortex / all seven"* and drew a 654× tolerance growth and a "+9.8%
step cost" from it. Two things were wrong with it:

* ⚠ **The baseline was a different cell.** `phase4_first.json` and `phase4_cortex.json` stand
  **4,361,496 nodes and 0 populations**; the seven-population run stands **4,558,554 nodes and 11
  populations**. The rows were not a series — they were three builds, and the difference between
  their tolerances contains the difference between their membranes.
* ⚠ **One run cannot support a wall-clock claim.** "+9.8%" was a difference between two single
  timings, quoted as if it were a measurement.

This is the same defect class the rest of the night is about, committed by me, ninety minutes after
committing a document about it. It was found by re-running rather than by re-reading.

## 2. The measurement that replaces it

Five runs per configuration, **all on one cell** — 4,558,554 nodes, 11 populations, 20 steps —
differing only in which strand populations are bound. Records: `world_phase4/spread_*.json`.

| | A: membrane faces only | B: + cortex | C: all seven |
|---|---:|---:|---:|
| strand force terms | 0 | 8,188,866 | **8,489,866** |
| **tolerance [pN]** | 2.80399e-06 | 2.40232e-04 | **1.75492e-03** |
| **tolerance spread, 5 runs** | **0.000%** | **0.000%** | **0.000%** |
| residual range [pN] | 1.18e-13 – 2.95e-13 | 2.40e-14 – 2.73e-13 | 4.63e-12 – 5.45e-12 |
| **residual spread, 5 runs** | **149%** | **1037%** | **18%** |
| s/step range | 0.1423 – 0.1547 | 0.1513 – 0.1677 | 0.1595 – 0.1813 |
| s/step median | 0.1485 | 0.1590 | 0.1642 |

## 3. What this says that the single-run table could not

**The tolerance is deterministic to every digit printed. The quantity it bounds is not.** Across five
identical runs the tolerance does not move at all — 0.000%, not "small" — while the residual it is
supposed to bound moves by a factor of **eleven** in configuration B.

That is a sharper statement of the defect than "the tolerance is too loose". The gate is **exact about
the wrong thing**: it is a deterministic function of the configuration — how many terms, of what
magnitude — and the run's actual force imbalance is a stochastic quantity it never consults. Two runs
of the identical configuration, one with a residual an order of magnitude larger than the other, are
scored against a bound that is bit-identical between them.

**Tolerance C / A = 626×, exactly**, with no uncertainty to report, because there is none.

## 4. And what does NOT survive

⚠ **The step-cost claim survives only at the ends.** A and C are **disjoint** — A tops out at 0.1547
and C starts at 0.1595 — so *8.49 M strand force terms cost a resolved increase, median +10.6%*. But
**A vs B and B vs C each OVERLAP**, so the intermediate attribution is not resolved at n = 5 and
neither "+7.1%" nor "+3.3%" may be quoted. §8's "+9.8%" was a difference of two single draws from
distributions that are ±7% wide.

The qualitative conclusion is unchanged and is the one worth keeping: **adding twenty-one times the
force terms is nearly free in wall clock, and multiplies the acceptance tolerance by 626.** The model
gets richer at almost no cost; its gate gets three orders of magnitude weaker for nothing.

⚠ **Still no magnitude.** Seven declared test points, none sourced for this cell; zero of twenty steps
accepted in every run above.

## 5. Method note

The spread was measured because a re-run on a re-synced host gave 0.1773 s/step where the first gave
0.1620, and the difference had to be either the code or the noise. Five repeats put both inside a
0.1595–0.1813 range and settled it as noise. **That question is only askable because the re-run was
done at all** — and it was done only because the provenance work forced a second execution, not
because anyone doubted the first number.

#!/usr/bin/env python
"""Re-run the equilibration search over the tau run III traces, INSIDE and OUTSIDE its search range.

⚠ **This decides nothing and claims nothing.** It re-executes `aleph/observe/stationarity.py`'s own
`equilibration_start` and `integrated_autocorrelation_time` over the gamma traces the three seeds
already wrote. Every tau, n_eff and reliability below is NOT QUOTABLE — `STATE.md` (c) 3 and (c) 17
retire every gamma this engine has produced, and none of these numbers is a gamma. They are about
WHICH WINDOW IS A MEASUREMENT, which is what PI queue item 18 asks.

⚠ **`stationarity.py:607` was NOT patched.** A shared threshold must not be moved by whoever holds
the result it would rescue; the PI rules item 18. Note also that `engine/observe/stationarity.py` is
a second copy bound to this one by `test_stationarity_modules_agree.py`, so widening the range can
change past verdicts — another reason it is not a local edit.

**Reproduction check, and it is why the rest is worth reading:** seed 3's reliability at the chosen
cut comes out 10.95, exactly the `min_tau_reliability` the run recorded.

⚠ **THE OUT-OF-GRID COLUMN EXISTS BECAUSE THE FIRST VERSION OF THIS SCRIPT GOT IT WRONG.** That
version scanned only the candidate grid, found reliability pinned around 9-16 on every seed, and
concluded the refusal was a LENGTH refusal that no cut could fix. That is a statement about the
region OUTSIDE the range, drawn from data INSIDE it — and repairing line 607 means widening exactly
that range. Lane C refuted it from its own scan and the refutation reproduces here to the digit.
The error is left described rather than deleted, because item 18 is precisely the question "what
changes if the range is widened", and answering it from inside the range closes it against data that
cannot speak to it.

What the scan actually shows:

1. The chosen cut is an ENDPOINT on all three seeds — the LAST grid candidate on seeds 1 and 3,
   index 0 on seed 2. The search's own range bound is what stops it: a censored optimum, and a real
   defect.
2. n_eff is nearly flat across the grid (about one effective sample between best and worst), so the
   endpoint is picked by a difference the criterion cannot resolve. An argmax on a flat objective
   carries no information.
3. Past the range bound, tau falls by about an order of magnitude and seeds 1 and 3 cross the bar,
   with an INTERIOR optimum out there — reliability rises, peaks, then falls again as the remaining
   sample count shrinks. So the optimum the search wanted is outside the range it was allowed.
4. ⚠ Crossing the bar is not the same as being a measurement. Seed 3's tau falls monotonically
   across the far cuts, so its reliability rises because tau is still dropping, not because it
   settled. Seed 2 does not cross at any cut tried. **One seed of three produces a tau, and there is
   no tau scatter** — which is an input to D-3, not to D-2, and is what the run already established.

Run:  python tau3_cut_scan.py    (host-side postprocessing of stored output; no device needed)
"""
from __future__ import annotations

import inspect
import json
from pathlib import Path

import numpy as np

from aleph.observe.stationarity import (
    DEFAULT_MIN_TAU_WINDOWS,
    equilibration_start,
    integrated_autocorrelation_time,
)

HERE = Path(__file__).parent
#: The bar, IMPORTED rather than copied. Lane C caught the copy: the source is exactly the knob the
#: PI is ruling on right now (queue item 1 — this 50 against the engine copy's min_windows = 5.0), so
#: a literal here would go on printing `*` after the ruling moved it.
#: ⚠ NOT `/ 2`. That is the bar on n_eff (`stationarity.py:821`, `implied_min_n_effective`); the
#: column below is `reliability = T/tau`, which is compared against min_tau_windows itself
#: (`:951`, `required_s = min_tau_windows * tau_int_s`). The table proves it: where reliability
#: reads 108.8, n_eff reads 54.4. Halving it would have quietly halved the bar — the same trap as
#: the copy, wearing the other face.
BAR = DEFAULT_MIN_TAU_WINDOWS

#: Cuts to report. The first three are ON the candidate grid; the rest are PAST its bound, which is
#: the region item 18 is actually about. Declared here rather than chosen from the output.
CUTS = (0, 7_500, 14_625, 16_000, 18_900, 20_000, 22_000, 24_000, 26_000)


def main() -> int:
    n_cand = inspect.signature(equilibration_start).parameters["n_candidates"].default
    traces = {}
    for seed in (1, 2, 3):
        rec = json.loads((HERE / f"tau3_seed{seed}.json").read_text())
        traces[seed] = (np.asarray(rec["gamma"]["trace_pn_per_um"], float),
                        float(rec["coefficients"]["dt_phys_s"]))

    for seed in (1, 2, 3):
        x, dt = traces[seed]
        limit = max(1, x.size // 2)                       # stationarity.py:607, verbatim
        stride = max(1, limit // int(n_cand))
        grid = list(range(0, limit, stride))
        res = equilibration_start(x, dt)
        where = ("LAST candidate — the search hit its own range bound" if res.start_index == grid[-1]
                 else "index 0 — the other endpoint" if res.start_index == 0 else "interior")
        print(f"seed {seed}: N={x.size} limit=N//2={limit} stride={stride} grid=[0 .. {grid[-1]}]"
              f"   chosen={res.start_index} -> {where}"
              f"   opens_on_transient={res.opens_on_transient}")

    bound = max(1, traces[1][0].size // 2)
    print(f"\n⚠ cuts above {bound} are OUTSIDE the search range — widening it is what item 18 asks about")
    print("    cut  in-grid |   s1 tau  reliab   n_eff |   s2 tau  reliab   n_eff |   s3 tau  reliab"
          "   n_eff     (NOT QUOTABLE — (c) 3 / (c) 17)")
    for cut in CUTS:
        cells = []
        for seed in (1, 2, 3):
            x, dt = traces[seed]
            e = integrated_autocorrelation_time(x[cut:], dt)
            mark = "*" if e.reliability >= BAR else " "
            cells.append(f"{e.tau_samples:9.1f} {e.reliability:7.1f}{mark}{e.n_effective:7.1f}")
        print(f" {cut:6d}   {'YES' if cut < bound else 'no':>5}  |" + "|".join(cells))
    print(f"\n  * = reliability >= the bar of {BAR:.0f}. ⚠ Crossing it is NECESSARY, not sufficient: "
          "seed 3's tau\n    is still falling monotonically out there, so its reliability rises "
          "because tau drops, not\n    because it settled. Seed 2 never crosses. One seed of three "
          "produces a tau; there is no\n    tau scatter, and that is an input to D-3, not D-2.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

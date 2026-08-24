#!/usr/bin/env python
r"""The transient table: where each seed's rise ends, and what its tail is still doing.

**This describes; it does not judge.** The declared read order for contract III lives in
`world_tau_read.py` and in the contract itself. This exists because seed 2 refused where seed 1 did
not, and the difference turned out not to be "a longer transient" — the tail slopes differ by an
order of magnitude — so the three seeds need a description before any of them needs a verdict.

⚠ **Nothing here is a criterion.** The plateau is the mean of the last 2,000 samples, which is a
convenient description of where the series ended up and **not** a definition of equilibrium. The
percentile columns say when the series first reached a fraction of it. The slope is an ordinary
linear fit. **None of these appears in any contract**, and none may be used to accept or reject a run.

⚠ **The columns do NOT disagree — they lose power together, and a third gate is what decides.**
This module first said seed 2 at cut 26,000 was a case of two tests giving different answers:
`opens_on_transient` False at 1.43 sigma against a fit of +1.43e-03 per sample. **That framing was
wrong and session 21 corrected it.** Both detectors weaken as the window shortens — at 4,000 samples
seed 2 sits just under `opens_on_transient` (1.43) AND just under the module's own correlation-aware
`drift_sigma` (2.75), while the same series fails both decisively at 8,000 and beyond. Neither is
what rejects it. `n_eff` is: 4,000 samples at tau 137.9 is **14.5 effective**, against a bar of 25.

**A short window can pass every single test. The conjunction is the instrument, not any one column,**
and that is why `n_eff` is printed here beside the two that can be talked past.

Sanity Gate:
    * dimensions — γ is pN/µm and slopes are pN/µm per SAMPLE, never per second: the driver's dt is
      its own clock and is not a physical time.
    * boundary cases — a window shorter than 100 samples is skipped rather than fitted; a fit over a
      handful of points is a number without being a measurement.
    * conservation — not applicable; this reads a stored series and computes no physical quantity.
    * sign sense — the slope's sign is the whole content of the seed-1/seed-2 difference and is
      printed, never taken as an absolute value.
    * measurement protocol — every seed gets the same cuts and the same columns, so a difference
      between rows is a difference between seeds and not between treatments.

engine units: pN/µm, samples. Runtime: host; no device.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

#: Fractions of the end-state level at which the rise is timed. Descriptive, not a criterion.
FRACTIONS = (0.90, 0.95, 0.98, 0.99, 0.999)

#: Cuts at which the tail is characterised. Fixed across seeds so rows are comparable.
CUTS = (18_900, 22_000, 26_000)


def main(argv: list[str] | None = None) -> int:
    """Print the transient table for the records given."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("records", type=Path, nargs="+")
    args = ap.parse_args(argv)

    import numpy as np

    from aleph.observe.stationarity import (
        integrated_autocorrelation_time,
        window_opens_on_a_transient,
    )

    def tau(x: np.ndarray) -> float:
        t = integrated_autocorrelation_time(x, 1.0, c=5.0)
        return float(getattr(t, "tau_int_s", getattr(t, "tau", t)))

    series: dict[str, np.ndarray] = {}
    for f in args.records:
        d = json.loads(f.read_text())
        g = (d.get("gamma") or {}).get("trace_pn_per_um")
        if not g:
            print(f"refused: {f} carries no gamma trace — a record that keeps a verdict and not "
                  "its series cannot be described, only quoted.")
            return 2
        series[f.stem] = np.asarray(g, np.float64)

    print("⚠ DESCRIPTIVE. No column here is a criterion and none appears in any contract.\n")
    print("WHERE THE RISE ENDS — first index at a fraction of the end-state level")
    print(f"{'record':>16} " + "".join(f"{f'{p*100:g}%':>9}" for p in FRACTIONS) + f"{'end-state':>12}")
    for name, g in series.items():
        plat = float(g[-2000:].mean()) if g.size >= 2000 else float(g.mean())
        row = "".join(f"{int(np.argmax(g >= p * plat)):>9,}" for p in FRACTIONS)
        print(f"{name:>16} {row}{plat:>12.2f}")

    print("\nWHAT THE TAIL IS DOING — same cuts for every seed")
    print(f"{'record':>16} {'cut':>8} {'kept':>8} {'tau':>9} {'n_eff':>7} {'opens':>6} {'sigma':>7} "
          f"{'slope/sample':>14} {'rise over win':>14} {'sd':>7}")
    for name, g in series.items():
        for c in CUTS:
            k = g[c:]
            if k.size < 100:                      # a fit over a handful of points is not a measurement
                print(f"{name:>16} {c:>8,} {k.size:>8,}  (skipped: fewer than 100 samples)")
                continue
            flag, sig = window_opens_on_a_transient(g, c, 1.0, c=5.0)
            slope = float(np.polyfit(np.arange(k.size, dtype=float), k, 1)[0])
            t = tau(k)
            n_eff = k.size / (2.0 * t) if t > 0 else float("nan")     # the column that actually decides
            print(f"{name:>16} {c:>8,} {k.size:>8,} {t:>9.1f} {n_eff:>7.1f} {str(flag):>6} {sig:>7.2f} "
                  f"{slope:>+14.3e} {slope * k.size:>+14.3f} {k.std():>7.3f}")
    print("\n⚠ `opens` and the slope do NOT disagree — both weaken as the window shortens. Read n_eff "
          "first: a window can sit under every sigma bar and still be too short to have measured "
          "anything. n_eff = kept/(2*tau) against a bar of 25.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

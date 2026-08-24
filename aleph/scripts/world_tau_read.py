#!/usr/bin/env python
r"""Read the τ verdict out of one or more PHASE 4 records, in the declared order and no other.

**This script chooses nothing.** The read order, the bars and the resolution criterion are
`TAU_RUN_CONTRACT_III_2026-08-21.md` §4, written before the run; this only fetches the fields and
prints them in that order so the person reading cannot start from the field they hoped for.

⚠ **`dt_s = 1.0`, and every axis is reported in SAMPLES.** The driver's `dt_phys_s` is 0.05 s and its
own help says *"This driver's clock, NOT a physiological quantity"*; the mobility is unsourced too, so
two independent unsourced quantities stand between this series and a τ in seconds. There are three
candidate "seconds per step" — wall-clock throughput, the integrator's dt, and the step count — and
**feeding the first produces a complete, plausible report in which every `*_s` field is a measurement
of how fast the machine ran, wearing a physics label.** The module cannot detect that: both are
positive finite floats. So this passes 1.0 and reads the `*_samples` twins.

Sanity Gate:
    * dimensions — none. Every field printed is a count, a ratio, or a dimensionless verdict.
    * boundary cases — fewer replicates than the module's minimum is REFUSED by the module, not
      worked around here; a refusal is printed as the result it is.
    * conservation — not applicable: this reads, it does not compute a physical quantity.
    * sign sense — not applicable.
    * measurement protocol — the order below IS the protocol, and it is not this script's to change.

engine units: none (samples). Runtime: host; no device.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# ⚠ As `run_provenance.py` does. A script in `scripts/` that is invoked by path rather than as a
# module does not have the repository root on sys.path, and the failure is a ModuleNotFoundError at
# the first first-party import — after the argument parsing has already succeeded, so it looks like a
# runtime fault rather than an invocation one.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def main(argv: list[str] | None = None) -> int:
    """Print the declared read order for the records given."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("records", type=Path, nargs="+", help="PHASE 4 .json records, one per seed")
    args = ap.parse_args(argv)

    import numpy as np

    from aleph.world.observe_gamma import observe_gamma_seed_scatter

    traces: dict[int, np.ndarray] = {}
    for f in args.records:
        d = json.loads(f.read_text())
        g = (d.get("gamma") or {}).get("trace_pn_per_um")
        if not g:
            print(f"refused: {f} carries no gamma trace.")
            return 2
        seed = int((d.get("provenance", {}).get("argv", "") .split("--seed ")[-1] or "0").split()[0]
                   if "--seed " in d.get("provenance", {}).get("argv", "") else len(traces))
        traces[seed] = np.asarray(g, np.float64)
        print(f"[read] {f.name}: seed {seed}, {len(g):,} samples", flush=True)

    print(f"\n[read] dt_s = 1.0 — every axis below is in SAMPLES, not seconds. See the module "
          f"docstring for why.\n", flush=True)
    try:
        rep = observe_gamma_seed_scatter(traces, dt_s=1.0)
    except Exception as exc:  # noqa: BLE001 — a refusal IS the result
        print(f"[read] REFUSED by the module: {type(exc).__name__}: {exc}")
        return 0

    row = rep.as_dict() if hasattr(rep, "as_dict") else dict(rep)

    # ⚠ THESE PRINT FIRST AND ALWAYS, WHATEVER THEIR LENGTH. Trimming long fields for legibility put
    # `step_acceptance` — the module's own statement that this series was NOT sampled from judged
    # steps — into a "not printed" list, one minute after I wrote that it belongs in front of any tau
    # reading. A caveat that is suppressed for being long is a caveat that is suppressed.
    for pinned in ("step_acceptance", "magnitude_claim"):
        if pinned in row:
            print(f"[read] ⚠ {pinned}: {row[pinned]}\n", flush=True)

    order = ["min_tau_reliability", "max_tau_relative_error", "max_tau_samples",
             "equilibration_samples", "window_samples", "required_window_samples",
             "min_n_effective", "drift_directions_agree", "scatter_inflation"]
    print("[read] DECLARED ORDER (contract III §4) — nothing later may be quoted if something "
          "earlier refuses:")
    for k in order:
        if k in row:
            print(f"   {k:26} {row[k]}")
    # ⚠ Long values are NAMED, not printed. `per_seed` and `contract` run to thousands of characters
    # each, and a wall of JSON after the declared order buries the declared order — which defeats the
    # one property this script exists to have. They are in the record; this says where.
    extra = [k for k in sorted(row) if k not in order and not k.endswith("_s")]
    pinned_already = {"step_acceptance", "magnitude_claim"}
    short = [k for k in extra if len(repr(row[k])) <= 120 and k not in pinned_already]
    long_ = [k for k in extra if len(repr(row[k])) > 120 and k not in pinned_already]
    if short:
        print("[read] other dimensionless fields:")
        for k in short:
            print(f"   {k:26} {row[k]}")
    if long_:
        print(f"[read] {len(long_)} long field(s) NOT printed, read them in the record itself: "
              + ", ".join(f"{k} ({len(repr(row[k])):,} chars)" for k in long_))
    dropped = [k for k in sorted(row) if k.endswith("_s")]
    if dropped:
        print(f"[read] ⚠ NOT PRINTED — {len(dropped)} `*_s` field(s) carry a seconds label over a "
              f"dt_s of 1.0 and would read as physical times: {', '.join(dropped)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

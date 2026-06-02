"""Reproducible single-spheroid proliferation-growth smoke run + G4 gate verdict (L2.4).

Grows one MCF7 CBM spheroid by contact-inhibited proliferation and checks the G4 mechanism
bands (oracle config ``layer2_cbm.yaml`` g4): rim-localised division + sub-exponential
(contact-inhibited bulk) growth. Prints the metrics the gate is checked against; the band
comparison lives here (the runtime module stays oracle-free).

Usage: python -m ffn_sim.scripts.layer2_growth_smoke
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.spheroid.observables import effective_radius
from ffn_sim.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.spheroid.proliferation import run_growth

_ROOT = Path(__file__).resolve().parents[1]
_CFG = _ROOT / "configs" / "layer2_cbm.yaml"
_OCFG = _ROOT / "validation" / "oracles" / "configs" / "layer2_cbm.yaml"


def main(argv: list[str] | None = None) -> int:
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    g4 = yaml.safe_load(_OCFG.read_text())["spheroid"]["acceptance"]["g4"]

    n0 = 250  # a spheroid large enough to have a genuine (inhibited) bulk + a rim
    total_time = 2.0 * prolif.cycle_time_mean
    res = run_growth(
        resolved, prolif, n_cells_init=n0, total_time=total_time,
        epoch_steps=1200, settle_steps=1000, seed=42, max_cells=4000,
    )

    exp_ceiling = 2.0 ** (total_time / prolif.cycle_time_mean)
    print(f"[growth smoke] N0={n0}  R0={effective_radius(res['a0'])*1e6:.1f} µm  "
          f"cycle={prolif.cycle_time_mean/3600:.0f} h  total={total_time/3600:.0f} h")
    print(f"  N: {res['n_cells'][0]} → {res['n_cells'][-1]}   growth factor = {res['growth_factor']:.2f}")
    print(f"  A/A0 (final) = {res['area_over_a0'][-1]:.3f}")
    print(f"  rim-localised division fraction = {res['rim_fraction_mean']:.2f}  "
          f"({res['n_division_epochs']} division epochs)")

    ok_rim = res["rim_fraction_mean"] >= g4["rim_division_fraction_min"]
    ok_sub = res["growth_factor"] < exp_ceiling
    print(f"\n[G4 rim ≥ {g4['rim_division_fraction_min']}] {'PASS' if ok_rim else 'FAIL'}  "
          f"({res['rim_fraction_mean']:.2f})")
    print(f"[G4 sub-exponential] {'PASS' if ok_sub else 'FAIL'}  "
          f"(growth {res['growth_factor']:.2f} < ceiling {exp_ceiling:.1f})")
    if res["capped_at_max_cells"]:
        print("  [note] run hit max_cells cap (growth truncated — not a full doubling horizon)")
    return 0 if (ok_rim and ok_sub) else 1


if __name__ == "__main__":
    sys.exit(main())

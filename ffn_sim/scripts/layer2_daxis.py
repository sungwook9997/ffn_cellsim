"""Layer-2 D-axis run: does the substrate-reacted COLLECTIVE BASAL CRAWL close the magnitude gap?

Bracket the faithful substrate-crawl (substrate_crawl.py, crawl_mode='basal') at the anchored
magnitudes on the substrate+catch CBM at a PI-overlapping R0, matched to the §C decisive bracket
so the rim-only(eject) vs basal-collective(this) contrast is direct:
  f=0       (substrate only, no crawl)
  f=1.6 nN  (whole-cell anchor)
  f=9.4 nN  (protrusion anchor; the one that EJECTED as a rim-only 3D force)

READ: if basal-crawl @9.4 nN spreads coherently (raw A/A0 climbs toward the PI 7-10, no eject) ->
the magnitude gap was the rim-only force-ROUTING, and the CBM DOES reach the magnitude with the
faithful substrate-reacted collective crawl. If it stays ~1-2 / ejects -> deeper structural limit.

Usage:  python -m ffn_sim.scripts.layer2_daxis <N0> <seed> [cpu|gpu]
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.archive.hoomd_legacy.spheroid.motility_bridge import resolve_active_traction
from ffn_sim.archive.hoomd_legacy.spheroid.observables import effective_radius
from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.archive.hoomd_legacy.spheroid.proliferation import run_growth_pooled
from ffn_sim.archive.hoomd_legacy.spheroid.substrate import resolve_substrate

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


def _device(kind):
    import hoomd
    return hoomd.device.GPU(notice_level=0) if kind == "gpu" else hoomd.device.CPU(notice_level=0)


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    n, seed = int(args[0]), int(args[1])
    kind = args[2] if len(args) > 2 else "cpu"
    cohesion_mode = args[3] if len(args) > 3 else "brittle"  # 'brittle' | 'yield'
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved, yield_remodel=(cohesion_mode == "yield"))
    sub = resolve_substrate(resolved)
    at = resolve_active_traction(resolved)
    total_time = 2.0 * prolif.cycle_time_mean
    device = _device(kind)
    mc = int(2.5 * n) + 1000

    for label, f in [("f0", 0.0), ("f_wholecell", at.f_wholecell), ("f_protrusion", at.f_protrusion)]:
        t0 = time.perf_counter()
        res = run_growth_pooled(
            resolved, prolif, n_cells_init=n, total_time=total_time, seed=seed,
            cohesion="catch", cad=cad, substrate=sub, f_traction=f, crawl_mode="basal",
            max_cells=mc,
        )
        rec = {
            "arm": label, "crawl": "basal", "cohesion": cohesion_mode,
            "n": n, "seed": seed, "device": kind,
            "f_active_nN": round(f * 1e9, 2),
            "R0_um": round(float(effective_radius(res["a0"]) * 1e6), 1),
            "aa0_core": round(float(res["area_core_over_a0"][-1]), 3),
            "aa0_raw": round(float(res["area_raw_over_a0"][-1]), 3),
            "n_final": int(res["n_cells"][-1]),
            "ejected": bool(res.get("ejected", False)),
            "capped": bool(res.get("capped_at_max_cells", False)),
            "wall_s": round(time.perf_counter() - t0, 1),
        }
        print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

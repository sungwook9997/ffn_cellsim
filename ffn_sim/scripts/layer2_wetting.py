"""Layer-2 wetting axis: does stronger cell-SUBSTRATE adhesion flatten the 3D cap (raise A/A0)?

§E localised the residual ~5–9x magnitude gap to a 3D→2D WETTING/FLATTENING scope boundary:
geometrically a spheroid (vol ∝R³) flattening to a film of thickness h gives A/A0 ≈ 4R/(3h) ≈ 7–13
for R≈150 µm, h≈15–30 µm — i.e. the PI 7–10 magnitude is plausibly the spheroid melting into a
quasi-2D film, a cell-substrate-adhesion / wetting transition, NOT lateral collective migration.
But the current substrate adhesion is a PLACEHOLDER: `resolve_substrate` sets D_sub = adhesion_ratio·D_e
with adhesion_ratio=1.0 (= the cell-CELL cohesion), so the aggregate sits at balanced wetting (3D cap)
by construction (the L2.6 "MCF7 = 3D cap" finding).

This is the SENSITIVITY DIAGNOSTIC (PASSIVE wetting): sweep adhesion_ratio = D_sub/D_e and ask whether
the projected footprint A/A0 climbs toward the geometric flattening limit as substrate adhesion exceeds
cohesion (Young-Dupré / Brochard-Wyart aggregate-wetting: the spreading coefficient S = W_cs − 2γ flips
a rounded cap to a spread film as W_cs/γ grows). It does NOT anchor the ratio — the physiological
adhesion_ratio for MCF7-on-colI/laminin is being sourced by deep-research (→ SE rows) before any
production sweep. This run only tests whether the mechanism is RESPONSIVE in the right direction;
A/A0 vs the PI data stays overlay-only (never fitted).

Two modes:
  * PASSIVE (default, f_active=0): adhesion-only — tests whether stronger W_cs alone flattens.
  * ACTIVE (f_active>0, crawl_mode): the literature-defined complete-wetting test — Douezan/Beaune/
    Aslemarz-2024 show aggregate spreading is ACTIVE (motile cells pull a precursor MONOLAYER film
    from the edge) and the wetting state is set by S = W_cs − 2γ. So the decisive test couples the
    §E active polarity drive (`crawl_mode='plithotaxis'`) to the adhesion sweep, asking whether the
    S>0 regime WITH active motility produces a spreading monolayer (A/A0 jump) the passive run can't.

Usage:  python -m ffn_sim.scripts.layer2_wetting <N0> <seed> [cpu|gpu] [ratios e.g. 1,4,16]
                                                  [f_active_nN=0] [crawl_mode=edge]
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import yaml

from ffn_sim.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.spheroid.observables import effective_radius
from ffn_sim.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.spheroid.proliferation import run_growth_pooled
from ffn_sim.spheroid.substrate import resolve_substrate

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"
_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2" / "wetting"


def _device(kind):
    import hoomd
    return hoomd.device.GPU(notice_level=0) if kind == "gpu" else hoomd.device.CPU(notice_level=0)


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    n, seed = int(args[0]), int(args[1])
    kind = args[2] if len(args) > 2 else "cpu"
    ratios = [float(x) for x in args[3].split(",")] if len(args) > 3 else [1.0, 4.0, 16.0]
    f_active = (float(args[4]) * 1e-9) if len(args) > 4 else 0.0   # nN → N
    crawl_mode = args[5] if len(args) > 5 else "edge"
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved, yield_remodel=True)
    plith = None
    if crawl_mode == "plithotaxis":
        from ffn_sim.spheroid.plithotaxis import resolve_plithotaxis
        plith = resolve_plithotaxis(resolved)
    total_time = 2.0 * prolif.cycle_time_mean
    device = _device(kind)
    mc = int(2.5 * n) + 1000
    active = f_active > 0.0

    _OUT.mkdir(parents=True, exist_ok=True)
    tag = f"active{round(f_active*1e9)}_{crawl_mode}" if active else "passive"
    out_path = _OUT / f"wetting_N{n}_s{seed}_{tag}.jsonl"
    with out_path.open("w") as fh:
        for ratio in ratios:
            sub = resolve_substrate(resolved, adhesion_ratio=ratio)
            t0 = time.perf_counter()
            res = run_growth_pooled(
                resolved, prolif, n_cells_init=n, total_time=total_time, seed=seed,
                cohesion="catch", cad=cad, substrate=sub, f_traction=f_active,
                crawl_mode=(crawl_mode if active else "edge"), plith=plith,
                max_cells=mc, device=device,
            )
            rec = {
                "axis": ("wetting_active" if active else "wetting_passive"),
                "adhesion_ratio": ratio, "f_active_nN": round(f_active * 1e9, 2),
                "crawl_mode": (crawl_mode if active else "none"),
                "n": n, "seed": seed, "device": kind,
                "R0_um": round(float(effective_radius(res["a0"]) * 1e6), 1),
                "aa0_core": round(float(res["area_core_over_a0"][-1]), 3),
                "aa0_raw": round(float(res["area_raw_over_a0"][-1]), 3),
                "n_final": int(res["n_cells"][-1]),
                "ejected": bool(res.get("ejected", False)),
                "capped": bool(res.get("capped_at_max_cells", False)),
                "wall_s": round(time.perf_counter() - t0, 1),
            }
            print(json.dumps(rec), flush=True)
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""(c) Ligand-condition PRODUCTION: Bare / Pre / Lam4 A/A0(R0) at native N on GPU (FORM-level).

The deferred-but-tracked production axis (PI 2026-06-04). The §C–§F investigation established the
MAGNITUDE is the center-based structural limit (overlay-only); this run delivers the FORM-level
result the CBM CAN make: do the three PI ligand conditions SEPARATE and ORDER correctly across the
PI R0 range, from the mechanistic clutch-kinetics traction (A1) + the A4′ Lam4 partial-β1-uniformity?

Mechanism (all anchored, none fitted to PI):
  * A1 active edge-traction per condition from per-species integrin-clutch kinetics
    (`ligand_traction.resolve_ligand_traction`): col-I clutch (Bare/Pre) > laminin (Lam4, 0.61×);
    Bare<Pre = the flagged pV4D4 col-I density axis. → Bare 1.50 / Pre 2.50 / Lam4 1.53 nN (≤3 nN
    B1-stable band).
  * A4′ Lam4 = PARTIAL β1 uniformity (traction screening length Lp≈40 µm vs the edge 11 µm for
    Bare/Pre) — the collective resolution of the single-cell↔collective laminin split (REPORT §A4′):
    uniform-ish β1 engages the interior so Lam4's small-size penalty is removed/reversed, lifting it
    above Pre at large R0 (the PI Lam4>Pre>Bare collective ordering).
  * Common: catch/yield cohesion, substrate adhesion_ratio=1 (held common to isolate the ligand axis).

The PI medians (Bare 7.2 / Pre 7.5 / Lam4 10.0) are OVERLAY-ONLY (never fitted); the magnitude is the
known structural-limit scope boundary. The deliverable is the ORDERING + SEPARATION of the three
emergent A/A0 = a + b/R + c/R² curves.

Usage:  python -m ffn_sim.scripts.layer2_ligand_production [cpu|gpu] [n0s e.g. 1250,4000,8500,16000] [seeds e.g. 2]
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.spheroid.ligand_traction import resolve_ligand_traction
from ffn_sim.spheroid.observables import effective_radius
from ffn_sim.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.spheroid.proliferation import run_growth_pooled
from ffn_sim.spheroid.substrate import resolve_substrate

_ROOT = Path(__file__).resolve().parents[1]
_CFG = _ROOT / "configs" / "layer2_cbm.yaml"
_OUT = _ROOT / "outputs" / "layer2" / "ligand_prod"

_CONDITIONS = ("Bare", "Pre", "Lam4")
_LP_EDGE = 11.0e-6      # A1 edge screening length (Bare/Pre)
_LP_LAM4 = 40.0e-6      # A4′ partial-β1-uniformity sweet spot (Lam4)


def _device(kind):
    import hoomd
    return hoomd.device.GPU(notice_level=0) if kind == "gpu" else hoomd.device.CPU(notice_level=0)


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    kind = args[0] if len(args) > 0 else "gpu"
    n0s = [int(x) for x in args[1].split(",")] if len(args) > 1 else [1250, 4000, 8500, 16000]
    n_seeds = int(args[2]) if len(args) > 2 else 2

    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved, yield_remodel=True)
    sub = resolve_substrate(resolved, adhesion_ratio=1.0)
    ligs = {c: resolve_ligand_traction(c) for c in _CONDITIONS}
    lp = {"Bare": _LP_EDGE, "Pre": _LP_EDGE, "Lam4": _LP_LAM4}
    total_time = 2.0 * prolif.cycle_time_mean
    device = _device(kind)

    _OUT.mkdir(parents=True, exist_ok=True)
    out_path = _OUT / "ligand_prod.jsonl"
    done = set()
    if out_path.exists():  # resumable
        for line in out_path.read_text().splitlines():
            if line.strip().startswith("{"):
                d = json.loads(line); done.add((d["condition"], d["n"], d["seed"]))

    print("[c] ligand → active edge-traction (clutch kinetics; A4′ Lam4 Lp=40µm):")
    for c in _CONDITIONS:
        r = ligs[c]
        print(f"   {c:5s} {r.ligand:12s} f_traction={r.f_traction*1e9:.2f} nN  Lp={lp[c]*1e6:.0f}µm")

    with out_path.open("a") as fh:
        for c in _CONDITIONS:
            f_tr = ligs[c].f_traction
            for n in n0s:
                mc = int(2.5 * n) + 1000
                for s in range(n_seeds):
                    if (c, n, s) in done:
                        continue
                    t0 = time.perf_counter()
                    res = run_growth_pooled(
                        resolved, prolif, n_cells_init=n, total_time=total_time, seed=4000 + s,
                        cohesion="catch", cad=cad, substrate=sub, f_traction=f_tr, Lp=lp[c],
                        crawl_mode="edge", max_cells=mc, device=device,
                    )
                    rec = {
                        "axis": "ligand_production", "condition": c,
                        "f_traction_nN": round(f_tr * 1e9, 2), "Lp_um": round(lp[c] * 1e6, 1),
                        "n": n, "seed": s, "device": kind,
                        "R0_um": round(float(effective_radius(res["a0_core"]) * 1e6), 1),
                        "aa0_core": round(float(res["area_core_over_a0"][-1]), 3),
                        "aa0_raw": round(float(res["area_raw_over_a0"][-1]), 3),
                        "n_final": int(res["n_cells"][-1]),
                        "ejected": bool(res.get("ejected", False)),
                        "capped": bool(res.get("capped_at_max_cells", False)),
                        "wall_s": round(time.perf_counter() - t0, 1),
                    }
                    fh.write(json.dumps(rec) + "\n"); fh.flush()
                    print(f"  {c:5s} N0={n:5d} R0={rec['R0_um']:5.1f}µm "
                          f"A/A0 core={rec['aa0_core']:.2f} raw={rec['aa0_raw']:.2f} "
                          f"{'EJECT' if rec['ejected'] else ''} ({rec['wall_s']:.0f}s)", flush=True)

    try:
        subprocess.run([sys.executable, "-m", "ffn_sim.scripts.layer2_ligand_production_vis"],
                       check=False, timeout=300)
    except Exception as exc:  # noqa: BLE001
        print(f"[ligand_production] auto-viz skipped: {exc}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

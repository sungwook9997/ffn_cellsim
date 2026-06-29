"""Layer-2 §E run: does COHERENT COLLECTIVE TRACTION (SPP plithotaxis) close the magnitude gap?

§D eliminated cohesion as the ~5-9x A/A0 magnitude-gap suspect and localised it to DRIVING/
COORDINATION: a radially-symmetric basal crawl does net ~0 (isotropic pressure ⇒ static
equilibrium). This driver brackets the literature-faithful fix — the Smeets-2016 CIL-SPP
plithotaxis polarity field (`plithotaxis.py`, `crawl_mode='plithotaxis'`): persistent per-cell
polarity (D_r=0.05/min, MCF10A) + free-edge CIL repolarisation (f_cil=0.1/min) + the emergent
correlation from the catch+turnover cohesion (NO imposed Vicsek alignment — Garcia 2015/Henkes
2020). On the substrate+catch(yield) CBM at a PI-overlapping R0, matched to the §C/§D brackets:

  f=0       (substrate only, no active drive)                     — baseline
  f=1.6 nN  (MCF7 whole-cell net speed, motility_bridge)          — lower bound
  f=5.0 nN  (MCF10A self-propulsion v_m=1 µm/min · γ_cell, Smeets) — the PRINCIPLED breast-epi anchor
  f=9.4 nN  (Bieling protrusion, motility_bridge)                 — upper bound

READ: if the persistent+CIL polarity field spreads coherently (A/A0 climbs toward the PI 7-10,
no eject) the DRIVING/COORDINATION diagnosis is confirmed and the platform reaches the magnitude
with the literature-anchored collective mechanism. If A/A0 stays ~1-2 the gap is deeper still
(the center-based 1-particle structural limit; the fine-grained single-cell line is then required).

Anchors + citation-integrity: outputs/tag_kb/SE_REGISTRATION_CANDIDATES_2026-06-04_collective-
migration.md. Overlay-only: the PI A/A0 is NEVER fitted; the MCF10A/MDCK speeds are validation
anchors. Auto-viz at run end (production-driver-auto-viz rule).

Usage:  python -m ffn_sim.scripts.layer2_plithotaxis <N0> <seed> [cpu|gpu] [yield|brittle]
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.archive.hoomd_legacy.spheroid.motility_bridge import resolve_active_traction
from ffn_sim.archive.hoomd_legacy.spheroid.observables import effective_radius
from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.archive.hoomd_legacy.spheroid.plithotaxis import resolve_plithotaxis
from ffn_sim.archive.hoomd_legacy.spheroid.proliferation import run_growth_pooled
from ffn_sim.archive.hoomd_legacy.spheroid.substrate import resolve_substrate

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"
_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2" / "plithotaxis"

# MCF10A self-propulsion speed v_m ≈ 1 µm/min (Smeets et al. 2016, PNAS 113:14621; MSD fit) — the
# closest breast-epithelial anchor for the SPP self-propulsion. f_active = v_m·γ_cell (overdamped
# map, same as motility_bridge). NOT fitted to PI A/A0.
V0_MCF10A_SPP: float = 1.0e-6 / 60.0  # m/s


def _device(kind):
    import hoomd
    return hoomd.device.GPU(notice_level=0) if kind == "gpu" else hoomd.device.CPU(notice_level=0)


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    n, seed = int(args[0]), int(args[1])
    kind = args[2] if len(args) > 2 else "cpu"
    cohesion_mode = args[3] if len(args) > 3 else "yield"  # 'yield' (ductile) | 'brittle'
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved, yield_remodel=(cohesion_mode == "yield"))
    sub = resolve_substrate(resolved)
    plith = resolve_plithotaxis(resolved)        # Vicsek J=0 (default; correlations emergent)
    at = resolve_active_traction(resolved)
    f_mcf10a = V0_MCF10A_SPP * resolved.gamma_cell
    total_time = 2.0 * prolif.cycle_time_mean
    device = _device(kind)
    mc = int(2.5 * n) + 1000

    _OUT.mkdir(parents=True, exist_ok=True)
    out_path = _OUT / f"plithotaxis_N{n}_s{seed}_{cohesion_mode}.jsonl"
    arms = [
        ("f0", 0.0),
        ("f_wholecell", at.f_wholecell),
        ("f_mcf10a", f_mcf10a),
        ("f_protrusion", at.f_protrusion),
    ]
    with out_path.open("w") as fh:
        for label, f in arms:
            t0 = time.perf_counter()
            res = run_growth_pooled(
                resolved, prolif, n_cells_init=n, total_time=total_time, seed=seed,
                cohesion="catch", cad=cad, substrate=sub, f_traction=f,
                crawl_mode="plithotaxis", plith=plith, max_cells=mc, device=device,
            )
            rec = {
                "arm": label, "crawl": "plithotaxis", "cohesion": cohesion_mode,
                "n": n, "seed": seed, "device": kind,
                "f_active_nN": round(f * 1e9, 2),
                "D_r_per_min": round(plith.D_r * 60.0, 4),
                "f_cil_per_min": round(plith.f_cil * 60.0, 4),
                "psi": round(plith.psi, 3),
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

    # auto-viz at run end (best-effort; never fails the run) — production-driver-auto-viz rule.
    try:
        subprocess.run(
            [sys.executable, "-m", "ffn_sim.scripts.layer2_plithotaxis_vis", str(out_path)],
            check=False, timeout=300,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"[layer2_plithotaxis] auto-viz skipped: {exc}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

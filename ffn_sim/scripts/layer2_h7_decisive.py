"""DECISIVE re-run (2026-06-10): resolve the Layer-2 magnitude FORK with H.7's MECHANISTIC
single-cell traction instead of the heuristic anchors.

BACKGROUND
----------
The Layer-2 center-based CBM line is CLOSED: it reproduced the PI law A/A0 = a + b/R + c/R² in
FORM (r²=0.998, zero-calib) but its MAGNITUDE is bounded ~5–9× under the PI medians (7–10). The
``b/R`` term is the single-cell active traction. The original decisive experiment
(``layer2_decisive_traction.py``) bracketed it with HEURISTIC anchors (whole-cell 1.6 nN /
protrusion 9.4 nN) and found: small → A/A0 stays ~1–2; large → boundary cells EJECT. The fork's
two horns were "missing active magnitude" vs "1-particle CBM abstraction limit", but the bracket
was heuristic.

H.7 has now DERIVED the single-cell traction mechanistically (sarcomeric ventral-SF array,
``spheroid/h7_traction_seam.resolve_h7_traction``):
  * f_cell_test   ≈ 2.8 nN  — raw H.7 mechanism at the validated TEST-fiber density;
  * f_cell_physio = 102 nN  — same mechanism scaled to physiological density (Gil-Redondo 2023).

THE FORK, RE-RUN WITH MECHANISTIC ANCHORS
-----------------------------------------
Sweep ``f_active`` across the H.7 mechanistic markers + the old heuristics, at a matched R0, and
map A/A0(f_active) against the PI band [7–10]. READ IT:
  * if the curve ENTERS [7–10] at the H.7 mechanistic value → the gap was MISSING ACTIVE MAGNITUDE
    (H.7 supplies it);
  * if the curve SATURATES ~1–2 and then EJECTS before reaching [7–10] → the gap is the 1-particle
    CBM ABSTRACTION LIMIT — and H.7's mechanistic traction CONFIRMS route a (fine-grained
    single-cell) is where the magnitude lives, not a parameter the CBM was missing.

dt-substep (for f_active > the B1 ≤3 nN ceiling) is the SAME PI-authorized Layer-2-only pattern as
the original decisive script: a copy of the resolved params with dt_cfl/factor + epoch_steps×factor
keeps the biological time identical, only finer — the frozen integrator is NOT touched.

Usage:  python -m ffn_sim.scripts.layer2_h7_decisive <N0> <seed> [cpu|gpu] [--physio] [--smoke]
        -> one JSON line per f_active arm (f_active, A/A0 core+raw, ejected, wall) + a header
           line with the resolved H.7 seam values.
"""

from __future__ import annotations

import dataclasses
import json
import sys
import time
from pathlib import Path

import yaml

from ffn_sim.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.spheroid.h7_traction_seam import resolve_h7_traction
from ffn_sim.spheroid.observables import effective_radius
from ffn_sim.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.spheroid.proliferation import run_growth_pooled
from ffn_sim.spheroid.substrate import resolve_substrate

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"
_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2" / "h7_decisive"

# B1 stable-traction ceiling (ligand_traction.STABLE_TRACTION_CEILING). Above it the overdamped
# CBM step flings a marginally-detaching edge cell (v=F/γ·dt) → sub-step to resolve it.
_CEILING_N = 3.0e-9


def _device(kind: str):
    import hoomd
    return hoomd.device.GPU(notice_level=0) if kind == "gpu" else hoomd.device.CPU(notice_level=0)


def _dt_factor_for(f_active: float) -> int:
    """Sub-step factor: keep per-step displacement bounded for f_active above the ceiling.

    Scales so the >ceiling arms resolve edge-cell peeling as finely as the original decisive
    armB did at 9.4 nN (dt/5): factor ≈ ceil(5 · f_active / 9.4 nN), clamped to [1, 64].
    """
    if f_active <= _CEILING_N:
        return 1
    import math
    return max(1, min(64, math.ceil(5.0 * f_active / 9.4e-9)))


def _one_arm(resolved, prolif, cad, sub, *, label, n, seed, f_active, device, total_time,
             epoch_steps=1200, max_cells=None):
    r = resolved
    dtf = _dt_factor_for(f_active)
    es = epoch_steps
    if dtf != 1:
        r = dataclasses.replace(resolved, dt_cfl=resolved.dt_cfl / dtf)
        es = epoch_steps * dtf
    mc = max_cells if max_cells is not None else int(2.5 * n) + 1000
    t0 = time.perf_counter()
    res = run_growth_pooled(
        r, prolif, n_cells_init=n, total_time=total_time, seed=seed,
        cohesion="catch", cad=cad, substrate=sub, f_traction=f_active, max_cells=mc,
        epoch_steps=es,
    )
    return {
        "arm": label,
        "f_active_nN": round(f_active * 1e9, 3),
        "dt_factor": dtf,
        "R0_um": round(float(effective_radius(res["a0"]) * 1e6), 1),
        "aa0_core": round(float(res["area_core_over_a0"][-1]), 3),
        "aa0_raw": round(float(res["area_raw_over_a0"][-1]), 3),
        "n_final": int(res["n_cells"][-1]),
        "ejected": bool(res.get("ejected", False)),
        "capped": bool(res.get("capped_at_max_cells", False)),
        "wall_s": round(time.perf_counter() - t0, 1),
    }


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    n, seed = int(args[0]), int(args[1])
    kind = args[2] if len(args) > 2 else "cpu"
    physio = "--physio" in args
    smoke = "--smoke" in args

    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved)
    sub = resolve_substrate(resolved)              # common substrate (adhesion_ratio=1.0)
    h7 = resolve_h7_traction()
    total_time = 2.0 * prolif.cycle_time_mean
    device = _device(kind)

    # f_active sweep (N): 0 baseline; old heuristics 1.6/9.4 for comparison; H.7 mechanistic
    # test-density (~2.8 nN, ⭐); curve-fill 3.0(ceiling)/6.0; H.7 physiological 102 nN (⭐, --physio).
    f_test = h7.f_cell_test
    f_physio = h7.f_cell_physio
    arms = [
        ("baseline", 0.0),
        ("heur_wholecell_1p6", 1.6e-9),
        ("H7_mech_testdensity", f_test),     # ⭐ H.7 mechanistic, validated test density
        ("ceiling_3p0", 3.0e-9),
        ("mid_6p0", 6.0e-9),
        ("heur_protrusion_9p4", 9.4e-9),
    ]
    if smoke:
        arms = [("baseline", 0.0), ("H7_mech_testdensity", f_test), ("heur_protrusion_9p4", 9.4e-9)]
    if physio:
        arms.append(("H7_mech_physiological", f_physio))  # ⭐ Gil-Redondo physiological density

    header = {
        "header": True, "n": n, "seed": seed, "device": kind,
        "h7_source": h7.source, "h7_provisional": h7.provisional,
        "h7_per_SF_pN": round(h7.per_sf_traction * 1e12, 2),
        "h7_n_SF": h7.n_sf,
        "h7_f_cell_test_nN": round(f_test * 1e9, 3),
        "h7_f_cell_physio_nN": round(f_physio * 1e9, 1),
        "h7_density_gap_x": round(h7.density_gap_factor, 1),
        "pi_aa0_band": [7.0, 10.0],
    }
    print(json.dumps(header), flush=True)

    _OUT.mkdir(parents=True, exist_ok=True)
    recs = [header]
    for label, f in arms:
        rec = _one_arm(resolved, prolif, cad, sub, label=label, n=n, seed=seed,
                       f_active=f, device=device, total_time=total_time)
        rec = {"n": n, "seed": seed, "device": kind, **rec}
        print(json.dumps(rec), flush=True)
        recs.append(rec)
    (_OUT / f"h7_decisive_n{n}_s{seed}.json").write_text(json.dumps(recs, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

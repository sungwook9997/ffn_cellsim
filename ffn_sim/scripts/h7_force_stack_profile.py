"""H.7 force-stack profiler (Phase 0a) — decompose the dt=0 native ceiling.

WHY (NATIVE_HOT_LOOP_MIGRATION analysis): at the ratified production scale
``n_fil=1000`` (~17.7k particles) the ``native-dt0`` ceiling (~1,864 steps/s) is
only ~1.2x above the fast Python path (~1,543), so the HOOMD built-in force stack
(bond Harmonic + angle Harmonic + LJ excluded-volume) is ~80%+ of per-step time —
NOT Python re-entry. Removing Python (native BAOAB / forces / binder) can buy at
most ~1.2-1.35x at production scale; the 27.8k-55.6k steps/s target requires
cutting the force-computation cost itself. This script answers the single
decisive question that gates the whole migration plan:

    *Which built-in force dominates the dt=0 ceiling — bonds, angles, or the
     LJ/neighbor stack — and how sensitive is it to the neighbor-list buffer?*

It is a TIMING HARNESS, not a validation gate (dt=0, kT=0, no movement, no
physics claim). It records failures rather than hiding them.

Method: build the full production cortex once per (n_fil, subset), prune the
integrator's force list to a cumulative subset BEFORE the first ``run(0)`` (so
no force is ever attached-then-detached), set ``dt=0`` + a ``Brownian(kT=0)``
method, and time ``run(steps)``. Per-force cost = successive differences:

    baseline = forces:none ; bond = bond - none ;
    angle = (bond+angle) - bond ; LJ/pair = full - (bond+angle)

Run on gbook (GPU is the production target):

    PYTHONPATH=$HOME/ffn_cellsim python -m ffn_sim.scripts.h7_force_stack_profile \
        --device gpu --n-fil 1000 --steps 4000 --rbuff-sweep

CPU is a logic smoke only (``--device cpu --n-fil 20 --steps 200``).
"""

from __future__ import annotations

import argparse
import json
import time
from copy import deepcopy
from pathlib import Path
from typing import Any

import hoomd
import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation
from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"

# Cumulative force subsets; per-force cost is the successive difference.
SUBSETS: list[tuple[str, set[str]]] = [
    ("none", set()),
    ("bond", {"bond"}),
    ("bond+angle", {"bond", "angle"}),
    ("full", {"bond", "angle", "pair", "dihedral", "improper", "other"}),
]


def _device(name: str):
    if name == "gpu":
        return hoomd.device.GPU(notice_level=0)
    return hoomd.device.CPU(notice_level=0)


def _resolved(n_fil: int, batch_steps: int = 20000):
    """Resolve the same production params the hot-loop benchmark uses."""
    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["n_filaments"] = int(n_fil)
    cfg["cortex"]["demo_mode"] = True
    cfg.setdefault("cortex", {}).setdefault("myosin", {})
    cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    cfg["cortex"]["myosin"]["mesoscale_force_scaling"] = True
    cfg["cortex"]["myosin"]["batch_steps"] = int(batch_steps)
    p_cortex = resolve_h3_derived(cfg)
    p_myo = resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl, R_cell=p_cortex.R_cell)
    p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)
    return p_cortex, p_myo, p_xl


def _force_class(f) -> str:
    m = type(f).__module__
    if ".bond" in m:
        return "bond"
    if ".angle" in m:
        return "angle"
    if ".pair" in m:
        return "pair"
    if ".dihedral" in m:
        return "dihedral"
    if ".improper" in m:
        return "improper"
    return "other"


def _count_topology(sim: hoomd.Simulation) -> dict[str, int]:
    snap = sim.state.get_snapshot()  # setup-time only (NOT a per-step path)
    if snap.communicator.rank != 0:
        return {}
    return {
        "particles": int(snap.particles.N),
        "bonds": int(snap.bonds.N),
        "angles": int(snap.angles.N),
        "particle_types": len(list(snap.particles.types)),
        "bond_types": len(list(snap.bonds.types)),
    }


def _nlist_info(forces) -> dict[str, Any]:
    """Report the neighbor-list parameters of the first pair force (the LJ stack)."""
    for f in forces:
        nl = getattr(f, "nlist", None)
        if nl is not None:
            info: dict[str, Any] = {
                "pair_force": type(f).__name__,
                "nlist_type": type(nl).__name__,
            }
            for attr in ("buffer", "rebuild_check_delay", "default_r_cut"):
                try:
                    info[attr] = getattr(nl, attr)
                except Exception:  # noqa: BLE001
                    pass
            return info
    return {"pair_force": None}


def _build_dt0(*, n_fil: int, keep: set[str], device_name: str, seed: int,
               rbuff: float | None) -> tuple[hoomd.Simulation, dict[str, Any]]:
    p_cortex, p_myo, p_xl = _resolved(n_fil)
    hw = build_cortex_full_simulation(
        p_cortex, p_xlinks=p_xl, p_myosin=p_myo,
        device=_device(device_name), with_baoab=False, connected_mesh=True,
        rng=np.random.default_rng(seed),
    )
    sim = hw["sim"]
    # strip every Python per-step updater (myosin/xlink/etc.) — dt=0 ceiling.
    while len(sim.operations.updaters) > 0:
        sim.operations.updaters.remove(sim.operations.updaters[0])
    integ = sim.operations.integrator
    all_forces = list(integ.forces)
    classes = [_force_class(f) for f in all_forces]
    # prune forces NOT in `keep`, BEFORE the first run(0) (never attach→detach).
    for f, c in zip(all_forces, classes):
        if c not in keep:
            integ.forces.remove(f)
    kept = [type(f).__name__ for f in integ.forces]
    kept_classes = [_force_class(f) for f in integ.forces]
    nlist = _nlist_info(integ.forces)
    if rbuff is not None and nlist.get("pair_force"):
        for f in integ.forces:
            nl = getattr(f, "nlist", None)
            if nl is not None:
                try:
                    nl.buffer = float(rbuff)
                except Exception:  # noqa: BLE001
                    pass
    integ.dt = 0.0
    integ.methods.append(
        hoomd.md.methods.Brownian(
            filter=hoomd.filter.All(), kT=0.0, default_gamma=p_cortex.gamma_b
        )
    )
    meta = {
        "all_force_types": [type(f).__name__ for f in all_forces],
        "all_force_classes": classes,
        "kept_force_types": kept,
        "kept_force_classes": kept_classes,
        "nlist": nlist,
    }
    return sim, meta


def _time_subset(*, n_fil: int, label: str, keep: set[str], steps: int,
                 device_name: str, seed: int, rbuff: float | None) -> dict[str, Any]:
    try:
        sim, meta = _build_dt0(
            n_fil=n_fil, keep=keep, device_name=device_name, seed=seed, rbuff=rbuff
        )
        sim.run(0)  # force-eval-only warm; builds nlist once (no movement after)
        topo = _count_topology(sim)
        t0 = time.time()
        sim.run(int(steps))
        wall = time.time() - t0
        sps = float(steps / wall) if wall > 0 else float("nan")
        return {
            "label": label, "n_fil": int(n_fil), "steps": int(steps),
            "wall_s": wall, "steps_per_s": sps,
            "us_per_step": float(1e6 * wall / steps) if steps else float("nan"),
            "rbuff": rbuff, "topology": topo, **meta, "device": device_name,
        }
    except Exception as exc:  # noqa: BLE001 - record failures
        return {
            "label": label, "n_fil": int(n_fil), "status": "failed",
            "error_type": type(exc).__name__, "error": str(exc),
            "rbuff": rbuff, "device": device_name,
        }


def _breakdown(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Successive-difference per-force cost from the cumulative subset rows."""
    by = {r["label"]: r for r in rows if r.get("steps_per_s")}
    us = {k: by[k]["us_per_step"] for k in by if "us_per_step" in by[k]}
    out: dict[str, Any] = {"us_per_step": us}
    if {"none", "bond", "bond+angle", "full"} <= set(us):
        full = us["full"]
        comp = {
            "baseline(run-loop)": us["none"],
            "bond_harmonic": us["bond"] - us["none"],
            "angle_harmonic": us["bond+angle"] - us["bond"],
            "lj+nlist": us["full"] - us["bond+angle"],
        }
        out["per_force_us"] = comp
        out["per_force_pct_of_full"] = {
            k: (100.0 * v / full if full else float("nan")) for k, v in comp.items()
        }
        out["full_steps_per_s"] = by["full"]["steps_per_s"]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, nargs="+", default=[1000])
    ap.add_argument("--steps", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument(
        "--rbuff-sweep", action="store_true",
        help="also sweep the LJ neighbor-list buffer on the full stack",
    )
    ap.add_argument(
        "--rbuff-values", type=float, nargs="+",
        default=[0.0, 0.1, 0.25, 0.5, 1.0],
        help="nlist.buffer values (HOOMD length units) for the sweep",
    )
    ap.add_argument(
        "--json",
        default=str(PKG / "outputs" / "h7" / "production" / "force_stack_profile.json"),
    )
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    results: list[dict[str, Any]] = []
    for n_fil in args.n_fil:
        rows = []
        for label, keep in SUBSETS:
            print(f"RUN n_fil={n_fil} subset={label}", flush=True)
            row = _time_subset(
                n_fil=n_fil, label=label, keep=keep, steps=args.steps,
                device_name=args.device, seed=args.seed, rbuff=None,
            )
            rows.append(row)
            if row.get("status") == "failed":
                print(f"  FAIL {row['error_type']}: {row['error']}", flush=True)
            else:
                print(
                    f"  steps/s={row['steps_per_s']:.1f} "
                    f"us/step={row['us_per_step']:.1f}",
                    flush=True,
                )
        bd = _breakdown(rows)
        sweep = []
        if args.rbuff_sweep:
            for rb in args.rbuff_values:
                print(f"RUN n_fil={n_fil} rbuff={rb} (full stack)", flush=True)
                srow = _time_subset(
                    n_fil=n_fil, label=f"full@rbuff={rb}", keep=SUBSETS[-1][1],
                    steps=args.steps, device_name=args.device, seed=args.seed,
                    rbuff=rb,
                )
                sweep.append(srow)
                if srow.get("steps_per_s"):
                    print(f"  steps/s={srow['steps_per_s']:.1f}", flush=True)
        results.append({
            "n_fil": int(n_fil), "subset_rows": rows, "breakdown": bd,
            "rbuff_sweep": sweep,
        })
        if "per_force_us" in bd:
            print(f"\n=== n_fil={n_fil} per-force dt=0 breakdown (us/step) ===", flush=True)
            for k, v in bd["per_force_us"].items():
                pct = bd["per_force_pct_of_full"][k]
                print(f"  {k:24s} {v:9.1f} us  ({pct:5.1f}% of full)", flush=True)
            print(f"  full steps/s = {bd['full_steps_per_s']:.1f}\n", flush=True)

    out = Path(args.json)
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(
        {"results": results, "steps": int(args.steps), "device": args.device,
         "seed": int(args.seed)},
        open(out, "w"), indent=2,
    )
    print(f"WROTE {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

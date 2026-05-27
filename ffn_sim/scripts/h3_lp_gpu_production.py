#!/usr/bin/env python
"""Thin GPU production driver for the H.3 cortex L_p gate (A5000, Gbook).

Mirrors the ``H3_PRODUCTION_MEDIUM`` / ``H3_PRODUCTION_FULL`` fixture in
``ffn_sim/tests/test_cortex.py`` (equilibrate → snapshot loop → per-snapshot
checkpoint) but builds the simulation on ``hoomd.device.GPU()`` instead of the
CPU default, and writes trajectory + an L_p summary into
``ffn_sim/outputs/h3/production/``.

This is a *driver*, not a gate: it reports the ensemble L_p against the KU-1.1
first-principles band for PI review. The authoritative ✅ ratification stays
with the pytest gate (no-gate-loosening). Created 2026-05-27 per PI directive
(2-machine split: Gbook = A5000 GPU production).

Usage (on Gbook, env ffn_sim):
    python ffn_sim/scripts/h3_lp_gpu_production.py --scale medium --device gpu
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import yaml
import hoomd

from ffn_sim.cortex.cortex import build_cortex_simulation, resolve_h3_derived
from ffn_sim.common.filament_math import fit_persistence_length

PKG = Path(__file__).resolve().parents[1]          # ffn_sim/
CFG = PKG / "configs" / "phase1_h3.yaml"

# Scales copied verbatim from the test_cortex.py production_run fixture.
SCALES = {
    "smoke":  dict(n_fil=300,  n_eq=50_000,  n_snap=30,  interval=5_000),
    "medium": dict(n_fil=500,  n_eq=75_000,  n_snap=50,  interval=25_000),
    "full":   dict(n_fil=1000, n_eq=100_000, n_snap=100, interval=50_000),
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scale", choices=SCALES, default="medium")
    ap.add_argument("--device", choices=["gpu", "cpu"], default="gpu")
    args = ap.parse_args()
    s = SCALES[args.scale]

    with open(CFG) as f:
        cfg = yaml.safe_load(f)
    cfg["cortex"]["n_filaments"] = s["n_fil"]
    cfg["cortex"]["demo_mode"] = True          # relax cost ceiling at scale-down
    p = resolve_h3_derived(cfg)

    device = hoomd.device.GPU() if args.device == "gpu" else hoomd.device.CPU()
    sim, *_ = build_cortex_simulation(
        p, device=device, with_baoab=True, with_crosslinkers=False
    )
    n_part = p.n_filaments * p.beads_per_filament

    outdir = PKG / "outputs" / "h3" / "production"
    ckptdir = outdir / f"lp_{args.scale}_{args.device}"
    ckptdir.mkdir(parents=True, exist_ok=True)

    print(
        f"[start] scale={args.scale} device={args.device} "
        f"n_filaments={p.n_filaments} n_part={n_part} "
        f"hoomd={hoomd.version.version} gpu_build={hoomd.version.gpu_enabled}",
        flush=True,
    )

    t0 = time.time()
    sim.run(s["n_eq"])
    print(f"[equilibrate] {s['n_eq']} steps in {time.time() - t0:.1f}s", flush=True)

    frames = np.empty(
        (s["n_snap"], p.n_filaments, p.beads_per_filament, 3), dtype=np.float64
    )
    for k in range(s["n_snap"]):
        sim.run(s["interval"])
        with sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag)
            pos = np.asarray(snap.particles.position)
            inv = np.empty_like(tag)
            inv[tag] = np.arange(tag.size, dtype=tag.dtype)
            frames[k] = pos[inv].reshape(p.n_filaments, p.beads_per_filament, 3)
        np.savez_compressed(
            ckptdir / f"snapshot_{k:03d}.npz",
            frame=frames[k], k=k, n_snapshots=s["n_snap"],
            scale=args.scale, sim_timestep=int(sim.timestep),
        )
        el = time.time() - t0
        print(
            f"[snap {k + 1}/{s['n_snap']}] ts={int(sim.timestep)} "
            f"wall={el:.1f}s est_total={el * s['n_snap'] / (k + 1):.0f}s",
            flush=True,
        )

    # Ensemble L_p (same estimator as the pytest gate). Informational only.
    L_p_vals = [
        fit_persistence_length(frames[k], rest_length=p.rest_length).L_p_m
        for k in range(s["n_snap"])
    ]
    L_p_vals = [x for x in L_p_vals if np.isfinite(x)]
    L_p_mean = float(np.mean(L_p_vals)) if L_p_vals else float("nan")
    lo, hi = p.L_p_band_m
    in_band = bool(lo <= L_p_mean <= hi)

    result = {
        "scale": args.scale,
        "device": args.device,
        "n_filaments": int(p.n_filaments),
        "n_part": int(n_part),
        "n_snapshots": int(s["n_snap"]),
        "L_p_mean_m": L_p_mean,
        "L_p_mean_um": L_p_mean * 1e6,
        "KU11_band_um": [lo * 1e6, hi * 1e6],
        "in_band": in_band,
        "wall_s": time.time() - t0,
        "hoomd": hoomd.version.version,
        "gpu_build": bool(hoomd.version.gpu_enabled),
    }
    np.savez_compressed(outdir / f"lp_{args.scale}_{args.device}.npz", frames=frames)
    (outdir / f"lp_{args.scale}_{args.device}_result.json").write_text(
        json.dumps(result, indent=2)
    )
    print("RESULT " + json.dumps(result), flush=True)

    # Visualize-at-closeout (CLAUDE.md hard rule): every production run
    # auto-refreshes its figures so the physics is always inspectable —
    # baked into the driver so it can never be forgotten (PI 2026-05-28).
    # Best-effort: a viz failure must NOT discard the production result.
    try:
        import subprocess
        subprocess.run(
            [sys.executable, str(PKG / "scripts" / "h3_lp_vis.py"),
             "--scale", args.scale, "--device", args.device],
            check=True, cwd=str(PKG.parent),
        )
        print("FIGS auto-generated (visualize-at-closeout)", flush=True)
    except Exception as e:  # noqa: BLE001  (viz is non-critical to the gate)
        print(f"WARN visualize-at-closeout failed (result still valid): {e}",
              flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""KU-5.1 dendritic density driver — H.5 lamellipodium production gate.

H.5 brief §Validation acceptance — Gate KU-5.1:
    ≈ 100 barbed ends per μm² of WAVE plane at steady state.

Builds a Cell with `with_lamellipodium=True` (H.5 단계 2 wiring landed
2026-05-29, commit `3306eb1`), runs warmup + production, samples the
per-snapshot dendritic density (live barbed-end count ÷ WAVE plane area).
Mirrors the KU-3.5 driver structure (`h3_ku35_tension.py`) so the
operational pattern is consistent: per-sample ``PROGRESS`` print + diag
list, auto-viz subprocess at the end (deferred — see §Auto-viz below),
``--out`` for explicit output paths in multi-seed sweeps.

Scope (scaffold, 2026-05-29 autonomous /loop):
    * Driver structure + diag schema established.
    * Single-seed CLI working; multi-seed parallel via per-seed ``--out``.
    * Density measurement reads `lamellipodium_state.barbed_end_tags`
      live; capping ablation visible via the `n_capped` field.
    * Auto-viz subprocess wired to `h5_ku51_density_vis.py` (TO BE WRITTEN
      by a future iter or PI follow-up — best-effort, warns on failure).
    * dt-factor default `0.001` matches the H.5 단계 1 demo config; the
      KU-5.1 production-grade dt + n_warmup + n_sample sweep parameters
      need PI calibration once a Mac smoke (n_fil=60, n_sample=4) confirms
      the assembly is stable end-to-end.

NOT in scope (deferred to PI / next session):
    * Multi-seed sweep analysis script (h5_ku51_sweep_analysis.py) —
      can be modelled after h3_ku35_sweep_analysis.py once we have
      multiple seed JSONs to drive the design.
    * KU-5.2 Bieling F-V curve driver (separate gate, separate driver).
    * KU-5.3 Funk abortive branching driver (separate gate).
    * Validation acceptance band tightening: H.5 brief reports ≈ 100/μm²
      as the target. Tolerance band [50, 200]/μm² is a placeholder pending
      Bieling/Funk-anchored derivation (no empirical magic numbers rule).

Usage:
    python ffn_sim/scripts/h5_ku51_density.py --n-fil 120 --dt-factor 0.001
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

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cell.lamellipodium import resolve_h5_lamellipodium
from ffn_sim.cell.cell import build_cortex_full_simulation
from ffn_sim.common.production_policy import (
    add_production_device_args,
    require_production_device,
    validate_production_device_args,
)

PKG = Path(__file__).resolve().parents[1]
CFG_H3 = PKG / "configs" / "phase1_h3.yaml"
CFG_H5 = PKG / "configs" / "phase1_h5.yaml"


def _density_per_um2(lamel_state, wave_area_m2: float) -> float:
    """Per-μm² dendritic density from live barbed-end count + WAVE area.

    Mirrors :func:`ffn_sim.tests.validation.test_ku5x_lamellipodium.dendritic_density`
    so the driver, the analysis pipeline, and the gate test all share one
    definition. ``wave_area_m2`` comes from ``ResolvedH5.wave_area``.
    """
    n_be = len(lamel_state.barbed_end_tags)
    if wave_area_m2 <= 0:
        return float("nan")
    return n_be / (wave_area_m2 * 1.0e12)


def run(
    *,
    n_fil: int,
    seed: int,
    dt_factor: float = 0.001,
    n_warmup: int = 40_000,
    n_sample: int = 80,
    interval: int = 5_000,
    device: str = "gpu",
    allow_cpu_dev: bool = False,
    out: str | None = None,
) -> dict:
    require_production_device(
        device, allow_cpu_dev=allow_cpu_dev, hoomd_module=hoomd
    )

    cfg_h3 = yaml.safe_load(open(CFG_H3))
    cfg_h5 = yaml.safe_load(open(CFG_H5))
    cfg_h3["cortex"]["n_filaments"] = n_fil
    cfg_h3["cortex"]["demo_mode"] = True
    p_h3 = resolve_h3_derived(cfg_h3)
    nca = p_h3.n_filaments * p_h3.beads_per_filament

    # Lamellipodium dt + box geometry inherited from cortex (D1 spec).
    tau_bend = p_h3.gamma_b * p_h3.rest_length ** 3 / p_h3.bending_modulus
    dtc = dt_factor * tau_bend
    p_h5 = resolve_h5_lamellipodium(
        cfg_h5, L_box=p_h3.L_box, dt=dtc, kT=p_h3.kT,
    )

    dev = hoomd.device.GPU() if device == "gpu" else hoomd.device.CPU(notice_level=0)

    # KU-5.1 isolates the lamellipodium branched-network steady state:
    # cortex backbone runs but xlinks / myosin / ERM are off so the
    # density measurement is unconfounded by cortex contractility.
    hc = build_cortex_full_simulation(
        p_h3,
        p_xlinks=None,
        p_myosin=None,
        p_lamellipodium=p_h5,
        device=dev,
        with_baoab=True,
        constrained=False,
        rng=np.random.default_rng(seed),
    )
    sim = hc["sim"]
    lamel_state = hc["lamellipodium_state"]
    lamel_layout = hc["lamellipodium_layout"]
    elong = hc["elong_action"]
    branch = hc["branch_action"]
    cap = hc["cap_action"]

    sim.run(0)
    sim.run(n_warmup)

    diag = []
    t0 = time.time()
    print(
        f"[start] device={device} n_fil={n_fil} n_part_cortex={nca} "
        f"n_WAVE={p_h5.n_WAVE} wave_area_um2={p_h5.wave_area * 1e12:.2f} "
        f"dt_s={dtc:.3e} dt_factor={dt_factor} "
        f"n_warmup={n_warmup} n_sample={n_sample} interval={interval} "
        f"hoomd={hoomd.version.version} gpu_build={hoomd.version.gpu_enabled}",
        flush=True,
    )
    for k in range(n_sample):
        sim.run(interval)
        n_be = len(lamel_state.barbed_end_tags)
        n_capped = len(lamel_state.capped_tags)
        density = _density_per_um2(lamel_state, p_h5.wave_area)
        wall = time.time() - t0
        eta_min = wall * (n_sample - (k + 1)) / max(k + 1, 1) / 60.0
        diag.append(
            dict(
                step=int(sim.timestep),
                n_barbed_ends=n_be,
                n_capped=n_capped,
                density_per_um2=density,
            )
        )
        print(
            f"PROGRESS sample={k + 1}/{n_sample} step={int(sim.timestep)} "
            f"wall={wall:.1f}s eta={eta_min:.1f}min "
            f"barbed={n_be} capped={n_capped} "
            f"density={density:.2f}/um2 (target ≈ 100)",
            flush=True,
        )

    # Plateau (last third) density mean — KU-5.1 gate quantity.
    last_third = diag[max(0, len(diag) * 2 // 3):]
    plateau_density = float(np.mean([d["density_per_um2"] for d in last_third]))

    outdir = PKG / "outputs" / "h5" / "production" / "ku51"
    outdir.mkdir(parents=True, exist_ok=True)
    if out is None:
        out_json = outdir / f"ku51_n{n_fil}_s{seed}.json"
    else:
        out_json = Path(out)
        out_json.parent.mkdir(parents=True, exist_ok=True)

    result = dict(
        n_fil=int(n_fil),
        seed=int(seed),
        dt_s=float(dtc),
        dt_factor=float(dt_factor),
        wave_area_m2=float(p_h5.wave_area),
        n_WAVE=int(p_h5.n_WAVE),
        n_warmup=int(n_warmup),
        n_sample=int(n_sample),
        interval=int(interval),
        plateau_density_per_um2=plateau_density,
        target_density_per_um2=100.0,  # H.5 brief KU-5.1 target
        final_density_per_um2=diag[-1]["density_per_um2"],
        n_barbed_ends_final=diag[-1]["n_barbed_ends"],
        n_capped_final=diag[-1]["n_capped"],
        note=(
            "KU-5.1 dendritic density at WAVE plane. Driver scaffold "
            "(2026-05-29 autonomous /loop) — params + plateau band PI calibration pending."
        ),
        diag=diag,
    )
    out_json.write_text(json.dumps(result, indent=2))
    print("RESULT " + json.dumps({k: v for k, v in result.items() if k != "diag"}), flush=True)
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=120)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dt-factor", type=float, default=0.001)
    ap.add_argument("--n-warmup", type=int, default=40_000)
    ap.add_argument("--n-sample", type=int, default=80)
    ap.add_argument("--interval", type=int, default=5_000)
    add_production_device_args(ap, default="gpu")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    validate_production_device_args(ap, args)
    run(
        n_fil=args.n_fil,
        seed=args.seed,
        dt_factor=args.dt_factor,
        n_warmup=args.n_warmup,
        n_sample=args.n_sample,
        interval=args.interval,
        device=args.device,
        allow_cpu_dev=args.allow_cpu_dev,
        out=args.out,
    )

    # Visualize-at-closeout (CLAUDE.md hard rule + memory
    # feedback_production_driver_auto_viz): subprocess to a per-seed viz
    # script. h5_ku51_density_vis.py is NOT YET WRITTEN — scaffold stage
    # left this as the next-PR follow-up. Best-effort: a viz failure or a
    # missing vis script must NOT mask a valid seed JSON.
    try:
        import subprocess
        vis = PKG / "scripts" / "h5_ku51_density_vis.py"
        if not vis.exists():
            print(
                "WARN visualize-at-closeout skipped: "
                "h5_ku51_density_vis.py not yet written (seed result still valid).",
                flush=True,
            )
        else:
            indir = (Path(args.out).resolve().parent if args.out
                     else PKG / "outputs" / "h5" / "production" / "ku51")
            subprocess.run(
                [sys.executable, str(vis), "--indir", str(indir)],
                check=True, cwd=str(PKG.parent),
            )
            print("FIGS auto-generated (visualize-at-closeout)", flush=True)
    except Exception as e:  # noqa: BLE001
        print(
            f"WARN visualize-at-closeout failed (seed result still valid): {e}",
            flush=True,
        )


if __name__ == "__main__":
    main()

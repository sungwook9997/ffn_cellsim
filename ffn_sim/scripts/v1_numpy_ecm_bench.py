"""v1 numpy ECM wall-time benchmark, apples-to-apples vs HOOMD H.1.

Drives the same canonical Mikado topology (``ffn_sim/configs/phase1_h1.yaml``)
through ``N_steps`` of **Euler-Maruyama** updates with the v1-equivalent
numpy force kernel (`ffn_sim.validation.oracles.ecm.fiber_mechanics.compute_forces`
+ `ffn_sim.validation.oracles.ecm.cross_links.compute_xl_energy_and_forces`).
These are the same oracle modules that the HOOMD KU-1.24 energy gate
validates against to 1e-6 relative — so this run is the v1 reference
implementation on the H.1 topology.

Caveat: v1 ECM did not have D7 LJ (added in v2). HOOMD bench includes
LJ. The direct ratio is therefore an *upper bound* on the HOOMD/numpy
wall-time regression — bonds + angles + xl dominates the per-step cost
at this density, so the LJ contribution is small but present.

Usage::

    PYTHONPATH=. python ffn_sim/scripts/v1_numpy_ecm_bench.py \
        --n-steps 1000

Writes ``ffn_sim/outputs/h1/v1_numpy_bench.json``.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.ecm.mikado import resolve_derived
from ffn_sim.validation.oracles.ecm.cross_links import generate_cross_links
from ffn_sim.validation.oracles.ecm.fiber_mechanics import compute_forces
from ffn_sim.validation.oracles.ecm.fiber_network import generate_2d_fiber_network


CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_h1.yaml"
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h1"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-steps", type=int, default=1000)
    args = ap.parse_args()

    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    p = resolve_derived(cfg)
    print(
        f"[v1_numpy_ecm_bench] n_fibers={p.n_fibers}, "
        f"N_beads={p.n_fibers*p.beads_per_fiber}, dt={p.dt_cfl:.3e} s"
    )

    t_build = time.time()
    net = generate_2d_fiber_network(
        L_box=p.L_box, n_fibers=p.n_fibers, L_fiber=p.L_fiber,
        beads_per_fiber=p.beads_per_fiber, S_order=0.0, theta0=0.0,
        seed=p.seed,
        params={"resolved_from": "v1_numpy_ecm_bench"},
    )
    xls = generate_cross_links(net, stiffness=p.xl_stiffness)
    t_build_done = time.time()
    print(
        f"[v1_numpy_ecm_bench] topology + xl: {t_build_done - t_build:.2f} s "
        f"(n_xl={len(xls)})"
    )

    pos = net.bead_positions.copy()  # (F, N, 2) in [0, L)
    rng = np.random.default_rng(p.seed)
    noise_amp = float(np.sqrt(2.0 * p.kT * p.dt_cfl / p.gamma_b))
    drift_pref = p.dt_cfl / p.gamma_b

    # Warm-up: 1 force eval so first-step timing isn't biased by import-time
    # JIT / cache misses.
    _ = compute_forces(
        pos,
        rest_length=p.rest_length,
        stretching_modulus=p.stretching_modulus,
        bending_modulus=p.bending_modulus,
        box_size=p.L_box,
        cross_links=xls,
    )

    t_prod = time.time()
    for _ in range(args.n_steps):
        F = compute_forces(
            pos,
            rest_length=p.rest_length,
            stretching_modulus=p.stretching_modulus,
            bending_modulus=p.bending_modulus,
            box_size=p.L_box,
            cross_links=xls,
        )
        drift = drift_pref * F
        noise = rng.standard_normal(pos.shape) * noise_amp
        pos = pos + drift + noise
        pos = np.mod(pos, p.L_box)
    t_prod_done = time.time()

    prod_s = t_prod_done - t_prod
    steps_per_s = args.n_steps / prod_s
    wall_per_sim_s = prod_s / (args.n_steps * p.dt_cfl)

    print(
        f"[v1_numpy_ecm_bench] EulerMaruyama production: {prod_s:.2f} s "
        f"for {args.n_steps} steps -> {steps_per_s:.1f} steps/s"
    )
    print(f"[v1_numpy_ecm_bench] wall_per_simulated_s = {wall_per_sim_s:.3e}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out = {
        "integrator": "EulerMaruyama (v1 reference)",
        "force_kernel": "ffn_sim.validation.oracles.ecm (numpy)",
        "includes_lj": False,
        "n_fibers": p.n_fibers,
        "n_beads": p.n_fibers * p.beads_per_fiber,
        "n_xl_bonds": len(xls),
        "dt_s": p.dt_cfl,
        "build_s": t_build_done - t_build,
        "production_s": prod_s,
        "production_steps": args.n_steps,
        "production_steps_per_s": steps_per_s,
        "simulated_time_s": args.n_steps * p.dt_cfl,
        "wall_per_simulated_second": wall_per_sim_s,
    }
    out_path = OUTPUT_DIR / "v1_numpy_bench.json"
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"[v1_numpy_ecm_bench] wrote {out_path}")


if __name__ == "__main__":
    main()

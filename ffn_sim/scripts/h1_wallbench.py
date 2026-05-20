"""H.1 ECM Mikado wall-time benchmark (M2-rest).

Drives the canonical N≈66 k Mikado config (``ffn_sim/configs/phase1_h1.yaml``)
through the standard equilibration prelude + ``--n-steps`` (default
10 000) BAOAB no-shear steps, reports steps/s, and writes a small JSON
artefact under ``ffn_sim/outputs/h1/``.

Reference: H.1 brief §Validation acceptance — "M1 Max CPU wall-time per
simulated second ≤ 5× the v1 numpy time." Phase 0.2 baseline polymer
(100-bead, ``hoomd_polymer_sanity.py``) hits ~91 k steps/s with HOOMD-
native Brownian; the L-M BAOAB Updater pays a ~9× per-step penalty for
the Python ``cpu_local_snapshot`` re-entry. ECM (much larger N, sparse
LJ pair list) is expected to be slower per step still — this bench
measures by how much.

Usage::

    python ffn_sim/scripts/h1_wallbench.py [--n-steps 10000] [--n-baoab 1000]

CLI flags:
    --n-steps : Production BAOAB steps after the prelude (default 10000).
    --n-softstart : Soft-start clipped-Brownian steps (default 100).
    --n-baoab : BAOAB no-shear prelude steps (default 1000).
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.ecm.equilibrate import equilibrate_no_shear
from ffn_sim.ecm.mikado import build_mikado_simulation, resolve_derived


CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_h1.yaml"
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h1"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-steps", type=int, default=10000)
    ap.add_argument("--n-softstart", type=int, default=100)
    ap.add_argument("--n-baoab", type=int, default=1000)
    args = ap.parse_args()

    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    p = resolve_derived(cfg)
    print(f"[h1_wallbench] n_fibers={p.n_fibers}, N_beads={p.n_fibers*p.beads_per_fiber}, dt={p.dt_cfl:.3e} s")

    t_build = time.time()
    sim, updater, action = build_mikado_simulation(p, with_cross_links=True)
    t_build_done = time.time()
    print(f"[h1_wallbench] build_mikado_simulation: {t_build_done - t_build:.2f} s")

    t_prelude = time.time()
    diag = equilibrate_no_shear(
        sim, action, updater,
        n_softstart=args.n_softstart, n_baoab=args.n_baoab,
        rest_length=p.rest_length, gamma_b=p.gamma_b,
    )
    t_prelude_done = time.time()
    prelude_s = t_prelude_done - t_prelude
    prelude_steps = args.n_softstart + args.n_baoab
    print(
        f"[h1_wallbench] equilibrate_no_shear: {prelude_s:.2f} s for "
        f"{prelude_steps} steps -> {prelude_steps/prelude_s:.1f} steps/s"
    )
    print(f"[h1_wallbench]   max|F| before={diag['max_force_before']:.3e}, "
          f"after soft={diag['max_force_after_soft']:.3e}, "
          f"after BAOAB={diag['max_force_after_baoab']:.3e}")

    t_prod = time.time()
    sim.run(args.n_steps)
    t_prod_done = time.time()
    prod_s = t_prod_done - t_prod
    print(
        f"[h1_wallbench] BAOAB production: {prod_s:.2f} s for {args.n_steps} steps "
        f"-> {args.n_steps/prod_s:.1f} steps/s"
    )

    with sim.state.cpu_local_snapshot as s:
        F = np.asarray(s.particles.net_force)
        post_force_max = float(np.abs(F).max())
    print(f"[h1_wallbench] post-production max|F|: {post_force_max:.3e} N")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out = {
        "n_fibers": p.n_fibers,
        "n_beads": p.n_fibers * p.beads_per_fiber,
        "n_xl_bonds": sim.state.N_bonds - p.n_fibers * (p.beads_per_fiber - 1),
        "dt_s": p.dt_cfl,
        "build_s": t_build_done - t_build,
        "prelude_s": prelude_s,
        "prelude_steps_total": prelude_steps,
        "prelude_steps_per_s": prelude_steps / prelude_s,
        "prelude_max_force_before_N": diag["max_force_before"],
        "prelude_max_force_after_soft_N": diag["max_force_after_soft"],
        "prelude_max_force_after_baoab_N": diag["max_force_after_baoab"],
        "production_s": prod_s,
        "production_steps": args.n_steps,
        "production_steps_per_s": args.n_steps / prod_s,
        "post_production_max_force_N": post_force_max,
        "simulated_time_s": args.n_steps * p.dt_cfl,
        "wall_per_simulated_second": prod_s / (args.n_steps * p.dt_cfl),
    }
    with open(OUTPUT_DIR / "h1_wallbench.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"[h1_wallbench] wrote {OUTPUT_DIR / 'h1_wallbench.json'}")


if __name__ == "__main__":
    main()

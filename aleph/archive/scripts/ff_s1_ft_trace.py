"""S1 F(t) saturation study — plate force vs relaxation TIME at native.

Records the plate force through the whole relaxation (metrics['trace']) so one can watch
whether the force SATURATES or keeps changing over time — the poroelastic / relaxation
time-dependence. The saturated (plateau) force is the trustworthy one; if it never plateaus,
that is itself the finding (and no single modulus is well-defined without a defined time).
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from aleph.laws.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from aleph.laws.network_warp import simulate_whole_cell_compression_on_device
from aleph.scripts.ff_s1_sphere import LOAD_PHYSIO


def main() -> None:
    ap = argparse.ArgumentParser(description="S1 plate-force F(t) saturation study (native).")
    ap.add_argument("--nfil", type=int, default=70686)
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--trace-every", type=int, default=200)
    ap.add_argument("--strain", type=float, default=0.05)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--out", default="aleph/outputs/mech_hier/s1_sphere/ft_trace_native.json")
    a = ap.parse_args()

    cx = build_crosslinked_cortex(CortexParams(), n_filaments=a.nfil, n_xl=a.nfil,
                                  n_myo=max(1, a.nfil // 10), rng=np.random.default_rng(1))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    _, m = simulate_whole_cell_compression_on_device(
        cx, NMIIA_MINIFIL_STALL_PN, strain=a.strain, n_steps=a.steps,
        device=a.device, trace_every=a.trace_every, **LOAD_PHYSIO)
    dt = m.get("dt_s", 0.0)
    tr = m["trace"]
    for r in tr:
        r["t_s"] = r["step"] * dt
    with open(a.out, "w") as f:
        json.dump({"note": f"native F(t), NF={a.nfil}, strain={a.strain}, steps={a.steps}",
                   "dt_s": dt, "final_F_pN": m["F_plate_pN"], "strain": a.strain,
                   "R0_um": cx.R0_mean, "trace": tr}, f, indent=2)
    print(f"F(t): {len(tr)} pts, F {tr[0]['F_plate_pN']:.0f}->{tr[-1]['F_plate_pN']:.0f} pN, "
          f"final {m['F_plate_pN']:.0f} pN, dt={dt:.2e}s, T={a.steps*dt:.3f}s", flush=True)
    print(f"FT TRACE DONE -> {a.out}", flush=True)


if __name__ == "__main__":
    main()

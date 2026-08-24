"""S1 Step-2 numerical convergence — native rigid-plate E_fit vs relaxation steps.

Convergence-first protocol (PI 2026-07-14): before any physics claim, establish that the
native rigid-plate modulus is CONVERGED. Sweeps n_steps at fixed NF and reports E_fit(n_steps)
+ the constitutive fit; the plateau (<5% change) is the trustworthy native S1 modulus.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from aleph.scripts.ff_s1_sphere import s1_compression_sweep, gate_s1_hertz, LOAD_PHYSIO


def main() -> None:
    ap = argparse.ArgumentParser(description="S1 relaxation-step convergence (native rigid).")
    ap.add_argument("--nfil", type=int, default=70686)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--nsteps", default="2000,4000,8000,16000")
    ap.add_argument("--n-strain", type=int, default=3)
    ap.add_argument("--out", default="aleph/outputs/mech_hier/s1_sphere/convergence_nsteps_native.json")
    a = ap.parse_args()

    steps_list = [int(x) for x in a.nsteps.split(",")]
    strains = np.linspace(0.005, 0.03, a.n_strain)
    rows = []
    prev = None
    print("n_steps   E_fit[Pa]   shell_r2   verdict            d(E_fit)", flush=True)
    for ns in steps_list:
        sw = s1_compression_sweep(strains, n_filaments=a.nfil, n_steps=ns, device=a.device,
                                  load=LOAD_PHYSIO, with_stress=False, quiet=True)
        g = gate_s1_hertz(sw)
        e = g["E_fit_Pa"]
        d = "" if prev is None else f"{(e-prev)/prev*100:+.1f}%"
        rows.append({"n_steps": ns, "E_fit_Pa": e, "k_shell_pN_per_um": g.get("k_shell_pN_per_um"),
                     "shell_r2": g.get("shell_r_squared"), "verdict": g["verdict"]})
        print(f"{ns:6d}   {e:9.1f}   {g.get('shell_r_squared'):.3f}    {g['verdict']:22s} {d}", flush=True)
        prev = e

    with open(a.out, "w") as f:
        json.dump({"note": f"native rigid-plate relaxation-step convergence, NF={a.nfil}",
                   "rows": rows, "mcf7_band": [224.0, 279.0]}, f, indent=2)
    print(f"CONVERGENCE DONE -> {a.out}", flush=True)


if __name__ == "__main__":
    main()

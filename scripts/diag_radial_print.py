"""Read csf_penetration_diag.json and print radial profile of c, |grad c|."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


def main(json_path: Path) -> None:
    d = json.loads(json_path.read_text(encoding="utf-8"))
    print(f"Radial profile from {json_path}")
    print(f"Box centre at r/R0=0; spheroid surface at r/R0=1.0")
    print()
    for n_pass in [0, 2, 4]:
        rec = d["radial_profiles"][f"n_pass={n_pass}"]
        r_over_R0 = np.array(rec["r_over_R0"])
        c_r = np.array(rec["c_radial"])
        g_r = np.array(rec["grad_c_mag_radial"])
        print(f"--- n_passes = {n_pass} ---")
        print(f"{'r/R0':>8} {'c':>10} {'|grad c|':>10}")
        for i, r in enumerate(r_over_R0):
            if 0.0 <= r <= 1.5:
                print(f"{r:>8.3f} {c_r[i]:>10.4f} {g_r[i]:>10.4f}")
        print()


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1
              else "results/stage1a_pilot/csf_penetration_diag.json"))

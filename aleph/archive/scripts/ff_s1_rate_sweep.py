"""S1 poroelastic rate sweep — E_fit vs loading rate (cytoplasm drainage), saved to JSON.

Rate-dependence of the bare-cortex apparent modulus IS cell poroelasticity (Moeendarbary
2013): fast loading traps the interstitial fluid (stiff), slow loading drains it (soft).
Uses the validated ff_hertz_validation harness (drained biphasic + turnover + membrane).
"""

from __future__ import annotations

import argparse
import json

from aleph.scripts import ff_hertz_validation as fh


def main() -> None:
    ap = argparse.ArgumentParser(description="S1 poroelastic rate sweep (E_fit vs loading rate).")
    ap.add_argument("--nfil", type=int, default=2000)
    ap.add_argument("--steps", type=int, default=1600)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--rates", default="5.0,0.05,0.005", help="comma-separated loading rates [um/s]")
    ap.add_argument("--lp", type=float, default=1.0e-7)
    ap.add_argument("--kdrained", type=float, default=300.0)
    ap.add_argument("--koff", type=float, default=0.066)
    ap.add_argument("--fexcess", type=float, default=0.25)
    ap.add_argument("--out", default="aleph/outputs/mech_hier/s1_sphere/poroelastic_rate_native.json")
    a = ap.parse_args()

    rates = [float(x) for x in a.rates.split(",")]
    rows = []
    print("rate[um/s]  E_fit[Pa]  verdict", flush=True)
    for v in rates:
        r = fh.run(
            v_load_um_s=v, Lp=a.lp, K_drained_Pa=a.kdrained, xl_koff_per_s=a.koff,
            f_excess=a.fexcess, n_filaments=a.nfil, n_steps=a.steps, device=a.device, quiet=True,
        )
        rows.append({"rate_um_s": v, "E_fit_Pa": r["E_fit_Pa"], "verdict": r["verdict"]})
        print(f"{v:9.4g}  {r['E_fit_Pa']:8.1f}  {r['verdict']}", flush=True)

    payload = {
        "note": f"poroelastic rate sweep, NF={a.nfil}, {a.device}, "
                f"Lp={a.lp}, K_drained={a.kdrained}, koff={a.koff}, f_excess={a.fexcess}",
        "rows": rows,
        "mcf7_band": [224.0, 279.0],
    }
    with open(a.out, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"RATE SWEEP DONE -> {a.out}", flush=True)


if __name__ == "__main__":
    main()

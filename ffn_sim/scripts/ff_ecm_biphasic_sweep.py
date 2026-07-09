"""S6 lever-3 cross-check — clutch traction vs collagen stiffness (Chan-Odde / Bangasser biphasic).

Sweeps the collagen fiber persistence length ``--ecm-lp-um`` (κ = kBT·Lp, the ECM stiffness knob) and reads
the ENGAGED-clutch traction from each run's JSON (``ecm_metrics`` + ``traction_nN``). The molecular-clutch
model (Chan & Odde 2008; Bangasser et al. 2013) predicts traction is BIPHASIC in substrate stiffness — an
inverted-U with a peak at an intermediate optimum (~2–300 kPa for many cells). This is a LITERATURE CROSS-CHECK,
never a gate to pass: the honest caveat is that reconstituted collagen-I gels are soft, so at physiological Lp
the cell sits on the rising limb and the full peak only appears once the swept effective stiffness reaches ~kPa.

Runs the real driver (``ff_crawl_on_substrate.py``) as a subprocess per Lp so the physics + adhesion setup are
identical to a production run. Emits a traction-vs-stiffness figure (log-x) + a JSON of the swept points.

The effective per-node stiffness is estimated analytically as k_node ≈ κ/ℓ₀³ = kBT·Lp/ℓ₀³ (ℓ₀ = seg length),
and mapped to an equivalent Young's modulus via the Hertzian ``substrate.py`` bridge k_sub = 2·E·a/(1−ν²) so the
x-axis can be read in kPa alongside the Bangasser 2–300 kPa band. Both mappings are approximate and PI-gated —
they calibrate the axis for the cross-check, they are not registered constants.

Usage (run on the A5000, AFTER the main native run so the GPU is free):
    python ff_ecm_biphasic_sweep.py --cortex-fil 38000 --from-resting --microtubules --implicit \
        --dt-impl 0.05 --myosin-linear --ecm-fibers 3000 --steps 400 --kmc-every 50 \
        --lps 2,20,200,2000,20000,200000 --device cuda:0 --tag ff_ecm_biphasic
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import numpy as np

_SEG_UM = 0.5            # collagen discretisation segment ℓ₀ (ecm_mikado default); k_node ≈ κ/ℓ₀³
_KBT = 4.28e-3         # pN·µm (FF units) — matches ffn_sim.ff.units.KBT (kB·T_physio / PN_UM)
_NU = 0.45             # substrate Poisson ratio (substrate.py)
_A_UM = 0.05           # integrin contact radius a [µm] (substrate.py Hertzian bridge)


def _e_eff_kpa(lp_um: float) -> float:
    """Approximate equivalent Young's modulus [kPa] for a collagen Lp, via k_node=κ/ℓ₀³ and k_sub=2·E·a/(1−ν²).

    Both steps are order-of-magnitude axis calibrations (PI-gated), not registered constants — they let the
    biphasic x-axis be read in kPa next to the Bangasser 2–300 kPa optimum band.
    """
    kappa = _KBT * lp_um                       # pN·µm²
    k_node = kappa / _SEG_UM ** 3              # pN/µm  (bending stiffness per node scale)
    e_pa = k_node * (1.0 - _NU ** 2) / (2.0 * _A_UM)   # invert k_sub = 2·E·a/(1−ν²); E in pN/µm² = Pa (µm·pN units → Pa)
    return e_pa / 1e3                          # kPa


def _run_point(lp_um: float, args, tag: str) -> dict:
    """Run the driver once at collagen Lp=lp_um; return the parsed JSON (or {} on failure)."""
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    cmd = [sys.executable, os.path.join(os.path.dirname(__file__), "ff_crawl_on_substrate.py"),
           "--ecm", "--ecm-lp-um", str(lp_um), "--ecm-fibers", str(args.ecm_fibers),
           "--cortex-fil", str(args.cortex_fil), "--steps", str(args.steps),
           "--kmc-every", str(args.kmc_every), "--record-every", str(max(args.steps, 1)),
           "--dt-impl", str(args.dt_impl), "--device", args.device, "--tag", tag, "--out", args.out]
    if args.from_resting:
        cmd.append("--from-resting")
    if args.microtubules:
        cmd.append("--microtubules")
    if args.implicit:
        cmd.append("--implicit")
    if args.myosin_linear:
        cmd.append("--myosin-linear")
    env = dict(os.environ, PYTHONPATH=root, PYTHONUNBUFFERED="1")
    print(f"[sweep] Lp={lp_um:g} µm  (E_eff≈{_e_eff_kpa(lp_um):.2g} kPa) → running ...", flush=True)
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    jpath = f"{args.out}/figs/{tag}.json"
    if r.returncode != 0 or not os.path.exists(jpath):
        print(f"[sweep] Lp={lp_um:g} FAILED (rc={r.returncode})\n{r.stdout[-800:]}\n{r.stderr[-800:]}", flush=True)
        return {}
    d = json.load(open(jpath))
    m = d.get("ecm_metrics", {})
    print(f"[sweep] Lp={lp_um:g}: traction={d['clutch_on'].get('traction_nN'):.3f} nN  "
          f"recruit={m.get('recruit_nm', float('nan')):.1f} nm  bound={d['clutch_on'].get('bound_frac'):.2f}", flush=True)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lps", default="2,20,200,2000,20000,200000",
                    help="comma list of collagen Lp [µm] to sweep (thin fibril → thick bundle)")
    ap.add_argument("--cortex-fil", type=int, default=38000)
    ap.add_argument("--ecm-fibers", type=int, default=3000)
    ap.add_argument("--steps", type=int, default=400)
    ap.add_argument("--kmc-every", type=int, default=50)
    ap.add_argument("--dt-impl", type=float, default=0.05)
    ap.add_argument("--from-resting", action="store_true")
    ap.add_argument("--microtubules", action="store_true")
    ap.add_argument("--implicit", action="store_true")
    ap.add_argument("--myosin-linear", action="store_true")
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--out", default="ffn_sim/outputs/ff")
    ap.add_argument("--tag", default="ff_ecm_biphasic")
    args = ap.parse_args()

    lps = [float(x) for x in args.lps.split(",") if x.strip()]
    pts = []
    for lp in lps:
        d = _run_point(lp, args, f"_biph_lp{lp:g}")
        if not d:
            continue
        m = d.get("ecm_metrics", {})
        pts.append(dict(lp_um=lp, e_eff_kpa=_e_eff_kpa(lp), traction_nN=d["clutch_on"].get("traction_nN"),
                        recruit_nm=m.get("recruit_nm"), dens_nm=m.get("dens_nm"), coh=m.get("coh"),
                        bound_frac=d["clutch_on"].get("bound_frac"), n_grip=m.get("n_grip")))
    if not pts:
        print("[sweep] no successful points — aborting"); return
    json.dump(pts, open(f"{args.out}/figs/{args.tag}.json", "w"), indent=2)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        lp = np.array([p["lp_um"] for p in pts]); ekpa = np.array([p["e_eff_kpa"] for p in pts])
        tr = np.array([p["traction_nN"] for p in pts]); rc = np.array([p["recruit_nm"] for p in pts])
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
        ax[0].plot(lp, tr, "o-", color="#1f77b4", lw=2, ms=7)
        ax[0].set_xscale("log"); ax[0].set_xlabel("collagen persistence length Lp [µm]  (κ = kBT·Lp)")
        ax[0].set_ylabel("engaged-clutch traction [nN]"); ax[0].set_title("Traction vs collagen stiffness (biphasic cross-check)")
        ax[0].grid(True, which="both", alpha=0.3)
        ax[1].plot(ekpa, tr, "o-", color="#d62728", lw=2, ms=7)
        ax[1].axvspan(2, 300, color="#2ca02c", alpha=0.12, label="Bangasser 2013 optimum band (2–300 kPa)")
        ax[1].set_xscale("log"); ax[1].set_xlabel("equivalent Young's modulus E_eff [kPa]  (approx, PI-gated axis)")
        ax[1].set_ylabel("engaged-clutch traction [nN]"); ax[1].set_title("vs equivalent substrate modulus")
        ax[1].grid(True, which="both", alpha=0.3); ax[1].legend(loc="best", fontsize=8)
        fig.suptitle("S6 lever-3: molecular-clutch traction is biphasic in collagen stiffness (Chan-Odde / Bangasser) — LITERATURE CROSS-CHECK, not a gate", fontsize=10)
        fig.tight_layout()
        png = f"{args.out}/figs/{args.tag}.png"
        fig.savefig(png, dpi=130); print(f"[sweep] wrote {png} + {args.tag}.json ({len(pts)} points)")
    except Exception as e:
        print(f"[sweep] figure skipped ({e}); JSON written")


if __name__ == "__main__":
    main()

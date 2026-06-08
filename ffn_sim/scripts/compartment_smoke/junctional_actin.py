"""ENABLED-PATH smoke — junctional actin coupling (SCALAR ONLY; STUB).

junctional_actin is a STUB: the explicit HOOMD topology build path
(extend_snapshot_with_junctional_actin) raises NotImplementedError EVEN when all
constants are anchored (the enabled+anchored build is reserved for the
PI-sanctioned release), so NO HOOMD assembly smoke is possible. This harness
therefore (1) exercises the RUNNABLE pure scalar physics laws
(coupling_force_magnitude, catch_off_rate) with SMOKE-ONLY anchored constants,
(2) ASSERTS the HOOMD build path is blocked (documents the STUB barrier), and
(3) leaves a PLATFORM_PI_QUEUE entry for the build-path + constant ratification.

NOT a physics claim. The catch-bond constants are the DETAILED_PLANS candidates
(x_catch/x_slip are KB-4.17 Buckley-anchored; k_catch0/k_slip0 are ORDER_ESTIMATE;
the rest are DERIVED H.3 transfers) — used SMOKE-ONLY to evaluate the SHAPE. No
band pass/fail.

Run:  python ffn_sim/scripts/compartment_smoke/junctional_actin.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
for _p in (_HERE.parents[3], _HERE.parents[0]):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import _smoke_common as sc
from ffn_sim.junction.junctional_actin import (
    catch_off_rate,
    coupling_force_magnitude,
    extend_snapshot_with_junctional_actin,
    resolve_junctional_actin,
)

# --------------------------------------------------------------------------
# SMOKE-ONLY anchored constants (DETAILED_PLANS junctional_actin candidates)
# --------------------------------------------------------------------------
ANCHORED = dict(
    k_couple=1.0e-6,      # N/m  DERIVED (H.3 attach transfer ×10)
    k_anchor=1.0e-5,      # N/m  DERIVED
    anchor_r0=5.0e-9,     # m    DERIVED (catenin footprint)
    max_couple_dist=6.0e-8,  # m  DERIVED (H.3 acceptor reach)
    k_catch0=1.0,         # 1/s  ORDER_ESTIMATE
    x_catch=4.0e-9,       # m    SOLID (Buckley 2014 / KB-4.17)
    k_slip0=0.02,         # 1/s  ORDER_ESTIMATE
    x_slip=0.4e-9,        # m    SOLID (Buckley 2014 / KB-4.17)
    k_on=10.0,            # 1/s  DERIVED (H.3 crosslinker transfer)
)
N_CAD = 50


def run() -> dict:
    p = resolve_junctional_actin(
        {"junction": {"junctional_actin": {"enabled": True, **ANCHORED}}},
        kT=sc.KT, dt=6.98e-7, n_cadherin=N_CAD,
    )
    is_anchored = bool(getattr(p, "is_anchored", False))

    # (1) Scalar law: tensile-only coupling clutch (sign-sense).
    ext = np.linspace(-5e-9, 20e-9, 51)
    fcoup = np.array([coupling_force_magnitude(p, e) for e in ext])
    tensile_only = bool(np.all(fcoup[ext <= 0] == 0.0) and np.all(fcoup[ext > 0] > 0.0))

    # (2) Scalar law: biphasic catch-slip off-rate; locate F* numerically.
    F = np.linspace(0.0, 40e-12, 400)
    koff = np.array([catch_off_rate(p, f) for f in F])
    i_star = int(np.argmin(koff))
    F_star_num = float(F[i_star])
    F_star_analytic = (sc.KT / (ANCHORED["x_catch"] + ANCHORED["x_slip"])) * math.log(
        (ANCHORED["k_catch0"] * ANCHORED["x_catch"])
        / (ANCHORED["k_slip0"] * ANCHORED["x_slip"])
    )
    biphasic = bool(0.0 < F_star_num < 40e-12 and koff[i_star] < koff[0])
    engaged_fraction = ANCHORED["k_on"] / (
        ANCHORED["k_on"] + ANCHORED["k_catch0"] + ANCHORED["k_slip0"]
    )

    # (3) ASSERT the HOOMD build path is blocked (the STUB barrier).
    build_blocked = False
    build_msg = ""
    try:
        extend_snapshot_with_junctional_actin(object(), p)
    except NotImplementedError as exc:
        build_blocked = True
        build_msg = str(exc).split(".")[0]
    except Exception as exc:  # any other failure is informative too
        build_msg = f"{type(exc).__name__}: {exc}"

    ok = is_anchored and tensile_only and biphasic and build_blocked
    result = {
        "compartment": "junctional_actin",
        "verdict": "SCALAR_OK" if ok else "SCALAR_FAIL",
        "physics_claim": False,
        "hoomd_build": "BLOCKED (NotImplementedError, STUB reserved path)",
        "params": {**ANCHORED, "n_cadherin": N_CAD, "kT_J": sc.KT},
        "scalar_laws": {
            "is_anchored": is_anchored,
            "coupling_tensile_only": tensile_only,
            "catch_biphasic": biphasic,
            "F_star_numeric_pN": F_star_num * 1e12,
            "F_star_analytic_pN": F_star_analytic * 1e12,
            "F_star_match": bool(abs(F_star_num - F_star_analytic) < 1e-12),
            "koff_at_0_per_s": float(koff[0]),
            "koff_min_per_s": float(koff[i_star]),
            "engaged_fraction_at_F0": engaged_fraction,
        },
        "build_path": {
            "blocked": build_blocked,
            "message": build_msg,
        },
        "constants_smoke_only": [
            "k_catch0=1.0, k_slip0=0.02 (ORDER_ESTIMATE)",
            "k_couple/k_anchor/k_on/max_couple_dist (DERIVED H.3 transfers)",
        ],
        "note": (
            "SCALAR smoke ONLY — junctional_actin is a STUB: the HOOMD build path "
            "raises NotImplementedError even when anchored, so no assembly/step "
            "smoke exists. This exercises the runnable scalar laws + asserts the "
            "build barrier. PI queue: build-path implementation + catch-set "
            "ratification (depends on cadherin_junction GATE-J first)."
        ),
    }

    _figure(ext, fcoup, F, koff, F_star_num, result)
    return result


def _figure(ext, fcoup, F, koff, F_star, result) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.6))

    ax0.plot(ext * 1e9, fcoup * 1e12, "-", color="#01665e", lw=1.8)
    ax0.axvline(0, color="0.6", lw=1, ls="--")
    ax0.set_xlabel("bond extension Δr [nm]")
    ax0.set_ylabel("coupling force [pN]")
    ax0.set_title("Tensile-only coupling clutch (sign-sense)")

    ax1.plot(F * 1e12, koff, "-", color="#8c510a", lw=1.8)
    ax1.axvline(F_star * 1e12, color="crimson", lw=1.5,
                label=f"F* = {F_star*1e12:.1f} pN (catch min)")
    ax1.axvspan(5.0, 10.0, color="orange", alpha=0.12, label="5-10 pN target shape")
    ax1.set_xlabel("force F [pN]")
    ax1.set_ylabel("k_off(F) [1/s]")
    ax1.set_title("α-catenin catch-slip off-rate (biphasic)")
    ax1.legend(fontsize=8)

    sl = result["scalar_laws"]
    fig.suptitle(
        f"junctional_actin SCALAR smoke — {result['verdict']} | HOOMD build "
        f"BLOCKED (STUB) | engaged f={sl['engaged_fraction_at_F0']:.2f} (SMOKE-ONLY)",
        fontsize=9.5,
    )
    sc.save_fig(fig, "smoke_junctional_actin")
    plt.close(fig)


def main() -> int:
    res = run()
    path = sc.save_json("smoke_junctional_actin", res)
    sl = res["scalar_laws"]
    print(f"[junctional_actin smoke] {res['verdict']}  (SCALAR only; HOOMD build BLOCKED)")
    print(f"  F*={sl['F_star_numeric_pN']:.2f} pN (analytic {sl['F_star_analytic_pN']:.2f}), "
          f"engaged f={sl['engaged_fraction_at_F0']:.3f}, "
          f"biphasic={sl['catch_biphasic']}, build_blocked={res['build_path']['blocked']}")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "SCALAR_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())

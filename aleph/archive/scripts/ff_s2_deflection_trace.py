"""S2 deflection(t) saturation / convergence study.

Applies a fixed localized patch force and traces the patch deflection over time so we can see
whether the dimple reaches an elastic equilibrium (plateau) or keeps creeping under constant load
(the PI's "does the force saturate or keep changing over time" question, for the localized-load
stage). One run per force; writes deflection(t) JSON + a figure. NATIVE by default.

Usage: python -m aleph.scripts.ff_s2_deflection_trace --nfil 70686 --f -8000 --steps 20000 --device cuda:0
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from aleph.laws.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from aleph.laws.network_warp import simulate_whole_cell_compression_on_device
from aleph.scripts.ff_s1_sphere import LOAD_PHYSIO


def run_trace(nfil, f_node, n_steps, trace_every, device, seed=1):
    """One fixed-force run; return (dt_s, trace list of {step, patch_deflection_um}, final_defl)."""
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=nfil, n_xl=nfil,
                                  n_myo=max(1, nfil // 10), rng=np.random.default_rng(seed))
    _, m = simulate_whole_cell_compression_on_device(
        cx, NMIIA_MINIFIL_STALL_PN, n_steps=n_steps,
        turgor_every=(50 if device.startswith("cuda") else 20), device=device,
        trace_every=trace_every,
        local_load={"axis": (0.0, 0.0, 1.0), "cos_thresh": 0.92, "f_node_pN": float(f_node)},
        **LOAD_PHYSIO)
    return float(m["dt_s"]), m["trace"], float(m["patch_deflection_um"])


def main() -> None:
    ap = argparse.ArgumentParser(description="S2 deflection(t) saturation study.")
    ap.add_argument("--nfil", type=int, default=70686)
    ap.add_argument("--f", type=float, default=-8000.0, help="per-node patch force [pN] (<0 inward)")
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--trace-every", type=int, default=500)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--out", default="aleph/outputs/mech_hier/s1_sphere/s2_deflection_trace.json")
    ap.add_argument("--fig-out", default="aleph/outputs/mech_hier/figs/s2_deflection_saturation.png")
    a = ap.parse_args()

    dt_s, trace, final = run_trace(a.nfil, a.f, a.steps, a.trace_every, a.device)
    steps = np.array([t["step"] for t in trace], float)
    defl = np.array([t["patch_deflection_um"] for t in trace], float)
    t_s = steps * dt_s

    # saturation diagnosis: last-decade slope relative to the total excursion
    if len(defl) >= 4:
        tail = slice(len(defl) * 3 // 4, None)                       # final quarter of the trace
        d_tail = defl[tail]; t_tail = t_s[tail]
        slope = float(np.polyfit(t_tail, d_tail, 1)[0]) if np.ptp(t_tail) > 0 else 0.0   # µm/s
        frac_per_decade = abs(slope) * (t_s[-1]) / max(abs(defl[-1]), 1e-12)
        saturated = frac_per_decade < 0.05                           # <5% further change over the run's timespan
    else:
        slope = 0.0; frac_per_decade = float("nan"); saturated = False

    out = {"nfil": a.nfil, "f_node_pN": a.f, "n_steps": a.steps, "dt_s": dt_s,
           "sim_time_s": float(t_s[-1]) if len(t_s) else 0.0,
           "final_deflection_um": final, "tail_slope_um_per_s": slope,
           "tail_frac_change": frac_per_decade, "saturated": bool(saturated),
           "t_s": t_s.tolist(), "deflection_um": defl.tolist()}
    with open(a.out, "w") as fjson:
        json.dump(out, fjson, indent=2)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.6, 4.2))
        ax.plot(t_s * 1e3, defl, "-o", ms=3, color="#0072B2")
        ax.set_xlabel("sim time [ms]"); ax.set_ylabel("patch deflection Δ [µm]")
        ax.axhline(final, color="#999", ls="--", lw=0.8, label=f"final {final:+.3f} µm")
        ax.set_title(f"S2 deflection(t) — f={a.f:.0f} pN/node, NF={a.nfil}\n"
                     f"{'SATURATED' if saturated else 'STILL CREEPING'} "
                     f"(tail Δ/run={frac_per_decade:.2%}, reaches {t_s[-1]*1e3:.2f} ms)")
        ax.legend(fontsize=8); fig.tight_layout()
        fig.savefig(a.fig_out, dpi=130); plt.close(fig)
        print(f"wrote {a.fig_out}", flush=True)
    except Exception as e:
        print(f"[fig skip] {e}", flush=True)

    print(f"# S2 DEFLECTION(t): f={a.f:.0f}pN final={final:+.3f}um  reaches {t_s[-1]*1e3:.3f}ms  "
          f"tailΔ/run={frac_per_decade:.2%} -> {'SATURATED' if saturated else 'STILL CREEPING'}", flush=True)
    print(f"wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()

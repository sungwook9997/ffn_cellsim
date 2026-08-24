"""Visualize the on-device (GPU) loaded-shell FF cortex run (Stage 6j dynamic complement).

Runs ``network_warp.simulate_loaded_shell_on_device`` — the fully GPU-resident dynamic cortex
(bending + crosslinker springs + myosin contraction + STATE-DEPENDENT osmotic turgor, Guo 2017 +
NF2007 inextensibility reshape) — and renders:

  (A) the shell-state trajectory V/V0, ΔP(V), R_mean vs step — showing the turgor-pressurised shell
      is VOLUME-STABLE (the tiny expansion self-relieves the resting osmotic excess to ΔP≈0), the
      numerically-honest result once the CFL includes the stiff turgor breathing-mode (omitting it
      gives a pure CFL-artifact runaway V/V0→1e3, which is NOT physics).
  (B) the settled 3D actin shell.

The physics is device-agnostic (GPU↔CPU bit-parity, ``tests/ff/test_network_warp``); this figure is
generated on CPU but is identical to the A5000 run. Writes ``outputs/_archive/ff/figs/loaded_shell_gpu.png``.

Run: ``python -m aleph.viz.viz_loaded_shell``  (``--n 1000`` fibers, ``--device cuda:0`` on gbook).
"""

from __future__ import annotations

import os

import numpy as np

from aleph.laws.gamma_floor import (
    NMIIA_MINIFIL_STALL_PN,
    CortexParams,
    build_crosslinked_cortex,
)
from aleph.laws.network_warp import simulate_loaded_shell_on_device

OUTDIR = os.path.join(os.path.dirname(__file__), "..", "outputs", "_archive", "ff", "figs")


def render(n_filaments: int = 1000, n_steps: int = 40000, seed: int = 0,
           device: str = "cpu", outdir: str = OUTDIR) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    params = CortexParams()
    cx = build_crosslinked_cortex(params, n_filaments=n_filaments, n_xl=n_filaments,
                                  n_myo=max(1, n_filaments // 10), rng=np.random.default_rng(seed))
    pos, traj = simulate_loaded_shell_on_device(
        cx, NMIIA_MINIFIL_STALL_PN, n_steps=n_steps, turgor_every=20,
        record_every=max(1, n_steps // 80), device=device)
    steps = np.array([t["step"] for t in traj])
    vv0 = np.array([t["V_over_V0"] for t in traj])
    dP = np.array([t["dP"] for t in traj])
    R = np.array([t["R_mean"] for t in traj])

    os.makedirs(outdir, exist_ok=True)
    fig = plt.figure(figsize=(16, 5.5))

    axA = fig.add_subplot(1, 3, 1)
    axA.plot(steps, vv0, "-", color="navy", lw=2)
    axA.axhline(1.0, color="grey", ls=":", lw=1)
    axA.set_xlabel("step"); axA.set_ylabel("V / V₀")
    axA.set_ylim(0.9, 1.1)
    axA.set_title(f"Volume stability — loaded shell (N={n_filaments})\n"
                  f"V/V₀ holds ≈1 (incompressible osmotic shell)")

    axB = fig.add_subplot(1, 3, 2)
    axB.plot(steps, dP, "-", color="crimson", lw=2, label="ΔP(V) [pN/µm²]")
    axB.set_xlabel("step"); axB.set_ylabel("ΔP(V)  [pN/µm²]", color="crimson")
    axB.tick_params(axis="y", labelcolor="crimson")
    axB.axhline(40.0, color="crimson", ls=":", lw=1, alpha=0.6)
    axR = axB.twinx()
    axR.plot(steps, R, "-", color="seagreen", lw=2, label="R_mean [µm]")
    axR.set_ylabel("R_mean  [µm]", color="seagreen")
    axR.tick_params(axis="y", labelcolor="seagreen")
    axB.set_title("State-dependent turgor (Guo 2017) + shell radius\n"
                  "ΔP self-relieves 40→0 via the tiny expansion")

    axC = fig.add_subplot(1, 3, 3, projection="3d")
    off = cx.net.fiber_offsets
    segs = [pos[int(off[f]):int(off[f + 1])] for f in range(cx.net.n_fibers)]
    axC.add_collection3d(Line3DCollection(segs, colors="steelblue", linewidths=0.4, alpha=0.5))
    Rax = params.R_um * 1.2
    axC.set_xlim(-Rax, Rax); axC.set_ylim(-Rax, Rax); axC.set_zlim(-Rax, Rax)
    for a in (axC.set_xlabel, axC.set_ylabel, axC.set_zlabel):
        a("µm")
    axC.set_title(f"settled actin shell (GPU-resident, device={device})\n"
                  f"R={R[-1]:.1f}µm, V/V₀={vv0[-1]:.3f}")

    fig.suptitle("FF loaded cortex shell — fully GPU-resident dynamic run "
                 "(bending + crosslink + myosin + state-dependent turgor + reshape, Warp)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = os.path.join(outdir, "loaded_shell_gpu.png")
    fig.savefig(path, dpi=135)
    plt.close(fig)
    print(f"wrote {path}  (V/V0 {vv0[0]:.4f}->{vv0[-1]:.4f}, dP {dP[0]:.0f}->{dP[-1]:.0f}, "
          f"R {R[0]:.2f}->{R[-1]:.2f}µm)")
    return path


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--steps", type=int, default=40000)
    ap.add_argument("--device", type=str, default="cpu")
    a = ap.parse_args()
    render(n_filaments=a.n, n_steps=a.steps, device=a.device)

"""Cytosim parity oracle — validate FF kernels against the independent C++ Cytosim (Stage 6f).

ENGINE.md's design rests on Cytosim (Nédélec & Foethke's own C++ code) being the INDEPENDENT parity
oracle — "the active layer is never validated against only itself". This module drives a built
``sim`` binary and compares FF quantities to Cytosim's on identical inputs. It is NOT a repo
dependency: every entry point degrades gracefully (returns/raises a clear skip) when no binary is
found, so the committed test suite stays green on machines without Cytosim.

Build (macOS, one-time, user-authorized): clone https://gitlab.com/f-nedelec/cytosim, then
``cmake -B build -DCMAKE_BUILD_TYPE=Release && make -C build sim`` (needs only BLAS/LAPACK — Apple
Accelerate; no GUI). Point ``CYTOSIM_SIM`` at the resulting ``build/bin/sim``.

Cytosim's config units are **seconds · micrometres · pico-Newtons** — identical to the FF unit
system (`ff.units`), so quantities compare directly with no conversion.

FINDING (2026-06-30, `docs/v2_audit/FF_CYTOSIM_PARITY_2026-06-30.md`): on the bending-energy of a
fixed circular arc, **Cytosim matches the continuum κL/2R² at every resolution, while FF's NF2007
interior-triple discrete energy under-counts by the end-factor (n−2)/(n−1)** (exact: FF·(n−1)/(n−2)
recovers the continuum; converges as n→∞; ~17 % low at the cortex's n=7). Both are valid
discretisations of the same operator (Cytosim end-corrects; FF is faithful to the paper's published
interior sum) — surfaced to PI as a quantitative-fidelity choice, not silently changed (bending is
sub-dominant to the myosin/turgor terms in γ, so the γ-floor conclusion is unaffected).
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile

import numpy as np

from aleph.laws.fiber_network import build_fiber_network
from aleph.laws.forces_warp import bending_energy


def find_cytosim_sim() -> str | None:
    """Locate the Cytosim ``sim`` binary (env ``CYTOSIM_SIM`` or PATH). None if unavailable."""
    env = os.environ.get("CYTOSIM_SIM")
    if env and os.path.exists(env):
        return env
    return shutil.which("sim")


def _arc_points(R_arc: float, L: float, n: int) -> np.ndarray:
    """``n`` equispaced points on a centred circular arc of radius ``R_arc``, contour length ``L``."""
    s = np.linspace(-L / 2, L / 2, n)
    ang = s / R_arc
    pts = np.stack([R_arc * np.sin(ang), R_arc * np.cos(ang), np.zeros(n)], axis=1)
    return pts - pts.mean(axis=0)


def cytosim_fiber_energy(pts: np.ndarray, kappa: float, *, sim_bin: str,
                         viscosity: float = 0.05) -> float:
    """Run Cytosim on a fiber placed at the exact ``pts`` and return its reported bending energy.

    Zero time-step run (``run 0``) at kT=0 → reports the static elastic energy of the given shape.
    """
    n = pts.shape[0]
    L = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
    seg = L / (n - 1)
    shape = ", ".join(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}" for p in pts)
    with tempfile.TemporaryDirectory() as d:
        cfg = (f"set simul system {{ time_step=0.001; viscosity={viscosity}; kT=0; }}\n"
               f"set space cell {{ shape=sphere; }}\nnew cell {{ radius={max(20.0, 3 * L)}; }}\n"
               f"set fiber fil {{ rigidity={kappa}; segmentation={seg}; }}\n"
               f"new 1 fil {{ shape={shape}; }}\n"
               f"run 0 system {{ nb_frames=1; }}\n"
               f"report fiber:energy energy.txt {{}}\n")
        open(os.path.join(d, "p.cym"), "w").write(cfg)
        r = subprocess.run([sim_bin, os.path.join(d, "p.cym")], cwd=d,
                           capture_output=True, text=True, timeout=120)
        if r.returncode != 0:
            raise RuntimeError(f"cytosim sim failed: {r.stderr[-500:]}")
        lines = [ln for ln in open(os.path.join(d, "energy.txt"))
                 if ln.strip() and not ln.startswith("%")]
        if not lines:
            raise RuntimeError("no energy reported by cytosim")
        return float(lines[0].split()[2])


def bending_energy_parity(R_arc: float = 5.0, L: float = 2.0, kappa: float = 20.0,
                          n_list=(5, 9, 17, 33), *, sim_bin: str | None = None) -> list[dict]:
    """Compare FF vs Cytosim vs analytic bending energy of a circular arc across resolutions.

    Returns a list of per-n dicts (ff, cytosim, analytic, ratios). Raises if no ``sim`` binary.
    The analytic continuum value is κ·L/(2 R_arc²) (constant-curvature arc).
    """
    sim_bin = sim_bin or find_cytosim_sim()
    if sim_bin is None:
        raise RuntimeError("no Cytosim 'sim' binary (set CYTOSIM_SIM); parity oracle unavailable")
    E_an = 0.5 * kappa * (1.0 / R_arc**2) * L
    rows = []
    for n in n_list:
        pts = _arc_points(R_arc, L, n)
        net = build_fiber_network([pts], kappa=kappa)
        E_ff = bending_energy(net)                                   # default: end-corrected
        E_ff_raw = bending_energy(net, end_correction=False)         # paper-literal interior sum
        E_cy = cytosim_fiber_energy(pts, kappa, sim_bin=sim_bin)
        rows.append({"n": n, "ff": E_ff, "ff_raw": E_ff_raw, "cytosim": E_cy, "analytic": E_an,
                     "ff_over_analytic": E_ff / E_an, "ff_raw_over_analytic": E_ff_raw / E_an,
                     "cyto_over_analytic": E_cy / E_an, "ff_over_cytosim": E_ff / E_cy})
    return rows


if __name__ == "__main__":
    sb = find_cytosim_sim()
    if sb is None:
        print("Cytosim 'sim' not found — set CYTOSIM_SIM=/path/to/build/bin/sim")
        raise SystemExit
    print(f"Cytosim parity (bending energy of a circular arc), sim={sb}")
    print(f"{'n':>4} {'FF(corr)':>10} {'FF(raw)':>10} {'Cytosim':>10} {'analytic':>10} "
          f"{'FFc/cyto':>9} {'FFraw/an':>9}")
    for r in bending_energy_parity():
        print(f"{r['n']:>4} {r['ff']:>10.5f} {r['ff_raw']:>10.5f} {r['cytosim']:>10.5f} "
              f"{r['analytic']:>10.5f} {r['ff_over_cytosim']:>9.4f} {r['ff_raw_over_analytic']:>9.4f}")

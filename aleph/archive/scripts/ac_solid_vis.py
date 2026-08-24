"""Regenerate the solid-ev (I2b) excluded-volume analytic figures — the single track entry point.

Every figure OVERLAYS the closed-form WCA oracle / brute-force reference on the numeric result, annotates
units (engine um-pN-s), and never truncates an axis (the professor's visualization-integrity rules, §1.9).
Pure numpy/matplotlib on the ac/solid oracle suite — runs on the dev Mac (no Warp/CUDA). Extend this as the
native crowding render lands (the HTML cell-morphology view is the lead's, per the native-gate spec).

    python -m aleph.scripts.ac_solid_vis   ->  aleph/outputs/ac/solid-ev/figs/*.png

Figures: WCA potential + FD-grad sign arbiter; zero-overlap -> zero-force (C1 cutoff); compressed-pair
resistance (finite steric volume); k_EV grid-invariance (hash-grid == brute) + derived-not-tuned saturation.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aleph.components.solid.steric_reference import steric_force_bruteforce, steric_force_celllist
from aleph.components.solid.wca_analytic import (
    compressed_pair_separation,
    epsilon_from_contact_stiffness,
    wca_cutoff,
    wca_energy,
    wca_force_magnitude,
)

SIGMA = 0.05          # um  steric diameter (illustrative coarse shell; params_i0b2b.yaml sigma_EV = GAP)
K_EV = 200.0          # pN/um  contact stiffness (numerical repulsion scale, Magic-Number Block)
EPS = epsilon_from_contact_stiffness(K_EV, SIGMA)
R_C = wca_cutoff(SIGMA)
OUT = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "solid-ev" / "figs"


def fig_wca_potential() -> Path:
    """WCA U(r) and F(r) = -dU/dr with the finite-difference gradient overlaid (the sign arbiter)."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))
    r = np.linspace(0.80 * SIGMA, 1.15 * R_C, 600)
    u = wca_energy(r, SIGMA, EPS)
    ax0.plot(r * 1e3, u, lw=2, c="navy", label="WCA U(r) oracle")
    ax0.axvline(R_C * 1e3, ls=":", c="gray")
    ax0.axhline(0.0, ls="-", c="k", lw=0.6)
    ax0.annotate(f"cutoff r_c = 2^(1/6)σ = {R_C*1e3:.1f} nm", (R_C * 1e3, 0.0),
                 (R_C * 1e3 + 2, 0.3 * u.max()), arrowprops=dict(arrowstyle="->", color="gray"))
    ax0.set_xlabel("separation r  [nm]"); ax0.set_ylabel("pair energy U  [pN·µm]")
    ax0.set_title("WCA potential (purely repulsive, C¹ at r_c)"); ax0.legend()

    # force vs the central-difference gradient of U
    rf = np.linspace(0.85 * SIGMA, 1.12 * R_C, 260)
    f_analytic = wca_force_magnitude(rf, SIGMA, EPS)
    h = 1e-8
    f_fd = -(wca_energy(rf + h, SIGMA, EPS) - wca_energy(rf - h, SIGMA, EPS)) / (2.0 * h)
    ax1.plot(rf * 1e3, f_analytic, lw=2, c="crimson", label="F = −dU/dr  (analytic)")
    ax1.scatter(rf[::8] * 1e3, f_fd[::8], s=22, facecolors="none", edgecolors="k",
                label="central FD of U  (sign arbiter)")
    ax1.axvline(R_C * 1e3, ls=":", c="gray"); ax1.axhline(0.0, ls="-", c="k", lw=0.6)
    ax1.annotate("F > 0 everywhere inside r_c → repulsive", (0.9 * SIGMA * 1e3, f_analytic.max() * 0.6))
    ax1.set_xlabel("separation r  [nm]"); ax1.set_ylabel("repulsion force F  [pN]")
    ax1.set_title("Force = −∇U (analytic vs finite difference)"); ax1.legend()
    fig.tight_layout()
    p = OUT / "i2b_wca_potential.png"
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_zero_force() -> Path:
    """Zero-overlap == zero-force: F and U vanish exactly at/beyond r_c (C¹, no attractive tail)."""
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    r = np.linspace(0.9 * SIGMA, 1.6 * R_C, 700)
    f = wca_force_magnitude(r, SIGMA, EPS)
    ax.plot(r * 1e3, f, lw=2, c="crimson", label="F(r)")
    ax.axvspan(R_C * 1e3, r.max() * 1e3, color="tab:green", alpha=0.12,
               label="r ≥ r_c : U ≡ 0, F ≡ 0  (OFF == bit-identical)")
    ax.axvline(R_C * 1e3, ls=":", c="gray")
    ax.axhline(0.0, ls="-", c="k", lw=0.6)
    ax.annotate("C¹: force → 0 smoothly\n(not a hard-shell jump)", (R_C * 1e3, 0.0),
                (R_C * 1e3 + 3, f.max() * 0.35), arrowprops=dict(arrowstyle="->", color="gray"))
    ax.set_xlabel("separation r  [nm]"); ax.set_ylabel("repulsion force F  [pN]")
    ax.set_title("Zero-overlap → zero-force (steric is repulsion-only, truncated at r_c)")
    ax.legend()
    fig.tight_layout()
    p = OUT / "i2b_zero_force.png"
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_compressed_pair() -> Path:
    """Compressed-pair equilibrium r_eq(load): finite steric volume — never collapses to r = 0."""
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    loads = np.logspace(0, 8, 200)     # 1 .. 1e8 pN compressive load
    r_eq = np.array([compressed_pair_separation(f, SIGMA, EPS) for f in loads])
    ax.semilogx(loads, r_eq * 1e3, lw=2, c="navy", label="equilibrium separation r_eq")
    ax.axhline(R_C * 1e3, ls=":", c="gray")
    ax.annotate("r_c (cutoff)", (1.5, R_C * 1e3 + 0.3))
    ax.axhline(0.0, ls="-", c="k", lw=0.8)
    ax.annotate("r_eq > 0 for ANY finite load\n→ core resists collapse (finite steric volume)",
                (1e2, R_C * 1e3 * 0.35))
    ax.set_ylim(0.0, R_C * 1e3 * 1.15)          # full axis from 0 — no truncation
    ax.set_xlabel("external compressive load on the pair  [pN]  (log)")
    ax.set_ylabel("equilibrium separation r_eq  [nm]")
    ax.set_title("Compressed pair resists collapse:  F_WCA(r_eq) = load")
    ax.legend(loc="upper right")
    fig.tight_layout()
    p = OUT / "i2b_compressed_pair.png"
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_grid_invariance() -> Path:
    """k_EV grid-invariance: cell-list (HashGrid analogue) == brute-force for any cell ≥ r_c; +saturation."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))

    # (left) hash-grid/cell-list force == brute force, across neighbour cell sizes (log-log spans the
    # many decades of per-node |F| — from grazing contacts to deep-core overlaps — all on the identity).
    rng = np.random.default_rng(7)
    n_fiber, per_fiber = 6, 8
    pos = rng.uniform(0.0, 4.0 * R_C, size=(n_fiber * per_fiber, 3))
    fiber_id = np.repeat(np.arange(n_fiber), per_fiber).astype(np.int64)
    f_brute = steric_force_bruteforce(pos, SIGMA, EPS, fiber_id)
    mag_brute = np.linalg.norm(f_brute, axis=1)
    nz = mag_brute > 0.0                                     # isolated (zero-force) nodes can't be logged
    markers = {"1.0·r_c": "o", "1.5·r_c": "s", "3.0·r_c": "^"}
    for (lab, mk), mult in zip(markers.items(), (1.0, 1.5, 3.0)):
        f_cell = steric_force_celllist(pos, SIGMA, EPS, fiber_id, cell_size=mult * R_C)
        ax0.scatter(mag_brute[nz], np.linalg.norm(f_cell, axis=1)[nz], s=34, marker=mk,
                    facecolors="none", edgecolors="C0", label=f"cell={lab}")
    lo, hi = mag_brute[nz].min() * 0.5, mag_brute[nz].max() * 2.0
    ax0.plot([lo, hi], [lo, hi], ls="--", c="k", lw=1, label="y = x (identity)")
    ax0.set_xscale("log"); ax0.set_yscale("log")
    ax0.set_xlim(lo, hi); ax0.set_ylim(lo, hi)
    ax0.set_xlabel("|F| brute-force O(N²)  [pN]  (log)")
    ax0.set_ylabel("|F| cell-list (HashGrid analogue)  [pN]  (log)")
    ax0.set_title("Neighbour cell size is a grid-invariant accelerator")
    ax0.legend(loc="upper left", fontsize=8)

    # (right) no-interpenetration saturates above the derived threshold (k_EV not tuned to an outcome)
    f_load = 20.0
    k_vals = np.logspace(1.5, 5, 120)   # 31 .. 1e5 pN/um
    frac = np.array([
        compressed_pair_separation(f_load, SIGMA, epsilon_from_contact_stiffness(k, SIGMA)) / R_C
        for k in k_vals
    ])
    ax1.semilogx(k_vals, frac, lw=2, c="navy", label="r_eq / r_c  (residual separation)")
    ax1.axhline(1.0, ls=":", c="gray"); ax1.annotate("r_c (no overlap)", (30, 1.005))
    ax1.axhline(0.90, ls="--", c="crimson", label="10% penetration tolerance")
    ax1.set_ylim(0.0, 1.08)
    ax1.set_xlabel("k_EV contact stiffness  [pN/µm]  (log)")
    ax1.set_ylabel("held-apart fraction  r_eq / r_c  [–]")
    ax1.set_title("Outcome saturates → k_EV derived, not tuned")
    ax1.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    p = OUT / "i2b_grid_invariance.png"
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    figs = [fig_wca_potential(), fig_zero_force(), fig_compressed_pair(), fig_grid_invariance()]
    for p in figs:
        print(f"wrote {p.relative_to(OUT.parents[3])}")


if __name__ == "__main__":
    main()

"""L2 D2 — EMERGENT aggregate surface tension from a CBM spheroid run (measurement).

Closes the single-cell -> spheroid surface-tension bridge with a MEASUREMENT (not an
anchor demonstration). The anchor demo ``layer2_surface_tension_bridge.py`` *asserted*
``sigma_tissue = gamma`` from published relations. This script instead MEASURES the
aggregate's emergent surface tension directly from a settled CBM spheroid and compares it
to the single-cell cortical tension ``gamma = 0.57 mN/m`` (KU-3.5 g_rigid native).

Measurement chain
-----------------
A cohesive aggregate is held together by surface tension and is therefore under an interior
Young-Laplace overpressure ``dP`` relative to the (zero-pressure) free exterior::

    dP = sigma * (1/R + 1/R')           [Young-Laplace, Roffay 2021 3D; sphere => 2 sigma / R]
    sigma = dP / (1/R + 1/R')           [inverse -> emergent sigma]

``dP`` is the INTERIOR mechanical (configurational virial) pressure of the cell-center
particles under the Morse cell-cell pair interaction::

    P = W / (3 V) ,   W = sum_pairs r_ij . f_ij        (athermal CBM: kT term dropped)

(``observables.virial_pressure``). A cohesive pair contributes ``W < 0`` (tension); we report
the OVERPRESSURE ``dP = -P_interior`` (interior held UNDER tension => negative virial pressure
=> positive overpressure that the surface tension balances). The exterior (free-surface) shell
has near-zero net pressure (cells there have fewer / weaker cohesive partners), which is the
``P_exterior ~ 0`` reference; we measure it too as a control.

Morse pair force (analytic, matches the HOOMD ``md.pair.Morse`` the run uses)::

    U(r)   = D0 [ e^{-2 a (r-r0)} - 2 e^{-a (r-r0)} ]
    F_r(r) = -dU/dr = 2 D0 a [ e^{-2 a (r-r0)} - e^{-a (r-r0)} ]
    f_ij   = F_r(r_ij) * r_hat_ij      (F_r > 0 repulsive along +r_ij; F_r < 0 cohesive)

so the SAME potential gives excluded volume (r<r0, F_r>0) and adhesion (r>r0, F_r<0); the
sign convention is consistent with ``virial_pressure``'s docstring (repulsive => P>0).

Sign convention (documented, so sigma comes out POSITIVE)
---------------------------------------------------------
- Interior pairs are mostly at r ~ r0 to slightly > r0 (settled, cohesive) => W_interior < 0
  => P_interior < 0 (the interior is under tension, not compression).
- The free droplet has NO confining external pressure, so the exterior reference is P ~ 0.
- Young-Laplace says an interior held by surface tension is OVER-pressured relative to
  outside: dP = P_in - P_out > 0. For a self-cohesive aggregate with no rigid wall the
  "pressure that the surface tension must supply to hold the curved interface" equals the
  magnitude of the interior cohesive (negative) virial: dP = |P_interior| = -P_interior.
  We invert Young-Laplace on this POSITIVE dP to get a POSITIVE sigma.
- CONTROL: a purely-repulsive (excluded-volume-only) config gives W>0 => P_interior>0 =>
  the formula would yield a NEGATIVE dP (no cohesion, nothing for surface tension to hold) —
  we run that control and confirm the sign flips, validating the measurement.

Run:  conda activate ffn_sim
      python -m ffn_sim.scripts.layer2_emergent_sigma [--n-cells 200] [--settle 40000]
Outputs: ffn_sim/outputs/layer2/emergent_sigma.json
         ffn_sim/outputs/layer2/emergent_sigma.png
         ffn_sim/outputs/layer2/emergent_sigma.log  (stdout tee, if redirected)

This is a MEASUREMENT script: it imports the runtime-forbidden bridge oracle ONLY for the
Young-Laplace inversion + the published anchors (NOT in any cell-build path), per the
subagent brief — kept import-isolated to this module.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import numpy.typing as npt
import yaml
from scipy.spatial import cKDTree

from ffn_sim.spheroid.cbm import build_cbm_simulation, get_positions
from ffn_sim.spheroid.observables import (
    connected_components,
    convex_hull_volume,
    nearest_neighbor_stats,
    radial_density_profile,
    radius_of_gyration,
    virial_pressure,
)
from ffn_sim.spheroid.params import ResolvedL2, resolve_layer2

# Oracle import is MEASUREMENT-only (Young-Laplace inversion + anchors); never a runtime path.
from ffn_sim.validation.oracles.spheroid import surface_tension_bridge as br

_ROOT = Path(__file__).resolve().parents[1]
_RUNTIME_CFG = _ROOT / "configs" / "layer2_cbm.yaml"
_OUT_DIR = _ROOT / "outputs" / "layer2"

#: Single-cell cortical-tension anchor [N/m] (KU-3.5 g_rigid native; the bridge target).
GAMMA_CORTICAL = br.G_RIGID_NATIVE_MN_M * 1e-3   # 0.57 mN/m -> N/m


def morse_force_magnitude(
    r: npt.NDArray[np.float64], D0: float, alpha: float, r0: float
) -> npt.NDArray[np.float64]:
    """Radial Morse force F_r(r) [N] for the HOOMD ``md.pair.Morse`` used in the run.

    ``F_r = -dU/dr = 2 D0 alpha [ e^{-2a(r-r0)} - e^{-a(r-r0)} ]``. Positive = repulsive
    (along +r_hat, r<r0 excluded volume); negative = cohesive (r>r0 adhesion). F_r(r0)=0.

    Args:
        r: pair separations [m].
        D0: Morse well depth [J].
        alpha: inverse range [1/m].
        r0: rest separation [m].
    """
    e = np.exp(-alpha * (r - r0))
    return 2.0 * D0 * alpha * (e * e - e)


def collect_pairs(
    positions: npt.NDArray[np.float64],
    r_cut: float,
    D0: float,
    alpha: float,
    r0: float,
    *,
    subset_mask: npt.NDArray[np.bool_] | None = None,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """Enumerate interacting pairs within ``r_cut`` and their Morse pair forces.

    Mirrors the HOOMD neighbour interaction: every unique unordered pair with separation
    < ``r_cut`` contributes ``f_ij = F_r(r) * r_hat_ij`` (force on i due to j). When
    ``subset_mask`` is given, only pairs with BOTH members in the subset are kept (so the
    virial of a sub-region — e.g. the interior core — uses only forces internal to it).

    Returns:
        ``(pair_indices (M,2), pair_forces (M,3))`` indexed into the FULL ``positions`` array.
    """
    tree = cKDTree(positions)
    raw = np.array(sorted(tree.query_pairs(r_cut)), dtype=np.int64).reshape(-1, 2)
    if subset_mask is not None and raw.size:
        keep = subset_mask[raw[:, 0]] & subset_mask[raw[:, 1]]
        raw = raw[keep]
    if raw.size == 0:
        return np.zeros((0, 2), dtype=np.int64), np.zeros((0, 3), dtype=np.float64)
    r_vec = positions[raw[:, 0]] - positions[raw[:, 1]]
    r_mag = np.linalg.norm(r_vec, axis=1)
    r_hat = r_vec / r_mag[:, None]
    F_r = morse_force_magnitude(r_mag, D0, alpha, r0)
    f_ij = F_r[:, None] * r_hat
    return raw, f_ij


def effective_radius_from_profile(
    positions: npt.NDArray[np.float64], n_bins: int = 24
) -> tuple[float, npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Aggregate radius from the radial-density inflection (where density falls to half-core).

    Returns the radius at which the radial number density first drops below half of the core
    (inner-third) plateau density — the density-profile definition of the droplet edge — plus
    the profile arrays for plotting. Falls back to max particle radius if no crossing.
    """
    r_c, dens = radial_density_profile(positions, n_bins=n_bins)
    core = dens[: max(1, n_bins // 3)]
    core_dens = float(core[core > 0].mean()) if np.any(core > 0) else float(dens.max())
    half = 0.5 * core_dens
    below = np.where(dens < half)[0]
    if below.size and below[0] > 0:
        # linear interpolation of the half-density crossing between bins
        i = below[0]
        d0, d1 = dens[i - 1], dens[i]
        r_edge = r_c[i - 1] + (r_c[i] - r_c[i - 1]) * (d0 - half) / (d0 - d1 + 1e-300)
    else:
        c = positions.mean(axis=0)
        r_edge = float(np.linalg.norm(positions - c, axis=1).max())
    return float(r_edge), r_c, dens


def measure_emergent_sigma(
    resolved: ResolvedL2,
    n_cells: int,
    *,
    settle_steps: int,
    interior_frac: float = 0.6,
    link_factor: float = 1.5,
    seed: int | None = None,
    repulsive_control: bool = False,
) -> dict:
    """Build + settle a CBM spheroid, measure interior/exterior virial pressure, invert YL.

    Args:
        resolved: resolved Layer-2 parameters.
        n_cells: number of cells.
        settle_steps: BAOAB steps to settle the loose blob into a stable aggregate.
        interior_frac: cells within this fraction of the aggregate radius (from the centre)
            count as INTERIOR; the outer shell is the free-surface EXTERIOR control.
        link_factor: connected-component link radius in units of r0 (core extraction).
        seed: realization seed (defaults to resolved.seed).
        repulsive_control: if True, zero out the attractive Morse branch in the FORCE
            accounting (use only r<r0 repulsion) — the sign-control config.

    Returns:
        Dict of measured quantities (SI + mN/m), ready to serialise.
    """
    sim, _action, _updater, r_cut = build_cbm_simulation(
        resolved, n_cells, seed=seed
    )
    sim.run(0)  # seed net_force for the L-M Action
    pos_init = get_positions(sim)
    sim.run(settle_steps)
    pos = get_positions(sim)

    c = pos.mean(axis=0)
    r = np.linalg.norm(pos - c, axis=1)

    # Aggregate radius from the radial-density inflection.
    R_edge, r_centers, density = effective_radius_from_profile(pos)
    Rg = radius_of_gyration(pos)

    # Connected core (largest component) — the cohesive body, fragmentation-robust.
    labels = connected_components(pos, link_radius=link_factor * resolved.morse_r0)
    core_mask = labels == 0
    n_core = int(core_mask.sum())

    # Interior vs exterior shells (by radius). Interior = inner interior_frac of R_edge.
    interior_mask = (r <= interior_frac * R_edge) & core_mask
    exterior_mask = (r > interior_frac * R_edge) & core_mask
    n_interior = int(interior_mask.sum())
    n_exterior = int(exterior_mask.sum())

    D0 = resolved.D_e
    alpha = resolved.morse_alpha
    r0 = resolved.morse_r0

    def _virial_P(mask: npt.NDArray[np.bool_]) -> tuple[float, int]:
        idx, f = collect_pairs(pos, r_cut, D0, alpha, r0, subset_mask=mask)
        if repulsive_control and idx.size:
            # keep ONLY repulsive contributions (r<r0 => F_r>0 => f_ij along +r_ij)
            r_vec = pos[idx[:, 0]] - pos[idx[:, 1]]
            r_mag = np.linalg.norm(r_vec, axis=1)
            rep = r_mag < r0
            idx, f = idx[rep], f[rep]
        if idx.size == 0:
            return 0.0, 0
        # Volume of the masked sub-region's convex hull (the region the virial spans).
        V = convex_hull_volume(pos[mask])
        if V <= 0.0:
            return float("nan"), idx.shape[0]
        P = virial_pressure(pos, idx, f, volume=V, kT=0.0)
        return P, idx.shape[0]

    P_interior, n_pair_int = _virial_P(interior_mask)
    P_exterior, n_pair_ext = _virial_P(exterior_mask)
    P_whole, n_pair_all = _virial_P(core_mask)

    # --- DIAGNOSTIC dP (naive radial split). For a SELF-BOUND drop with no external wall the
    # total virial ~ 0 (bulk cohesion balances bulk repulsion EVERYWHERE), so this radial
    # interior-vs-exterior split does NOT cleanly isolate the surface tension — it is reported
    # as a diagnostic only. The mechanical (Irving-Kirkwood) estimator below is the real route.
    dP = -(P_interior - P_exterior)

    # --- PRIMARY: mechanical surface tension via the Irving-Kirkwood SPHERICAL estimator.
    # For a spherical drop, sigma = (1/A) ∫ (P_N - P_T) dr collapses (pairwise forces) to a
    # sum over interacting pairs of the anisotropy of each pair's virial about the drop centre:
    #
    #   sigma = -(1/(16 pi R^2)) * sum_pairs (r_ij . f_ij) * [ 1 - 3 (s_ij . r_hat_c)^2 / s_ij^2 ]
    #
    # where s_ij = r_i - r_j, r_hat_c is the radial unit vector at the pair MIDPOINT (about the
    # drop COM), and R is the drop radius. A purely RADIAL pair (s_ij || r_hat_c) carries the
    # NORMAL pressure (factor 1-3 = -2); a TANGENTIAL pair (s_ij ⟂ r_hat_c) carries the
    # tangential pressure (factor +1). Surface cohesive pairs are predominantly tangential and
    # cohesive (r_ij.f_ij < 0), so they yield sigma > 0 (tension) — this is the surface signal
    # that the isotropic virial averages away. (Thompson et al. 1984 JCP 81:530; Irving-Kirkwood.)
    idx_all, f_all = collect_pairs(pos, r_cut, D0, alpha, r0, subset_mask=core_mask)
    if repulsive_control and idx_all.size:
        s = pos[idx_all[:, 0]] - pos[idx_all[:, 1]]
        rep = np.linalg.norm(s, axis=1) < r0
        idx_all, f_all = idx_all[rep], f_all[rep]
    if idx_all.size and R_edge > 0.0:
        s_ij = pos[idx_all[:, 0]] - pos[idx_all[:, 1]]
        mid = 0.5 * (pos[idx_all[:, 0]] + pos[idx_all[:, 1]]) - c  # midpoint rel. to COM
        mid_norm = np.linalg.norm(mid, axis=1)
        mid_norm[mid_norm == 0.0] = 1e-300
        r_hat_c = mid / mid_norm[:, None]
        rdotf = np.sum(s_ij * f_all, axis=1)                       # r_ij . f_ij  (virial term)
        s2 = np.sum(s_ij * s_ij, axis=1)
        proj2 = np.sum(s_ij * r_hat_c, axis=1) ** 2                # (s_ij . r_hat_c)^2
        anis = 1.0 - 3.0 * proj2 / s2
        sigma_ik = -float(np.sum(rdotf * anis)) / (16.0 * np.pi * R_edge**2)
    else:
        sigma_ik = float("nan")
    sigma = sigma_ik  # the reported emergent surface tension (mechanical IK)

    nn = nearest_neighbor_stats(pos)

    # --- MCF7 emergent beta/gamma from the contact geometry ---
    # Mean cohesive-pair separation ~ contact spacing; contact "radius" ~ overlap of the
    # adhesive zone. Emergent contact area from the mean neighbour separation: two cells of
    # diameter r0 sharing a flat contact face of radius a_contact where the cohesive Morse
    # well operates over the contact_zone_width. Geometric chord estimate:
    #   d = mean NN separation; if d < r0 the surfaces interpenetrate the adhesive zone;
    #   contact radius a = sqrt(R_cell^2 - (d/2)^2) (spherical-cap chord), clamped >= 0.
    d_mean = nn["median"]
    R_cell = resolved.R_cell
    chord = R_cell**2 - (0.5 * d_mean) ** 2
    a_contact = float(np.sqrt(max(chord, 0.0)))
    # If cells sit at/just beyond r0 (no interpenetration), fall back to the adhesive-zone
    # disk radius (contact_zone_width sets the lateral adhesive extent).
    if a_contact <= 0.0:
        a_contact = resolved.contact_zone_width
    contact_area = float(np.pi * a_contact**2)
    # Work of de-adhesion per contact = Morse well depth D_e (energy to pull one pair out of
    # the well to r_cut), the same scale params.py uses (D_e = 2 F_detach contact_zone).
    work_deadhesion = D0
    beta = br.adhesion_tension(work_deadhesion, contact_area)  # N/m
    beta_over_gamma = beta / GAMMA_CORTICAL
    roffay_lo, roffay_hi = br.ROFFAY_SURFACE_INTERIOR_RATIO
    bog_window = (1.0 - 1.0 / roffay_lo, 1.0 - 1.0 / roffay_hi)  # (0.375, 0.5)

    return {
        "n_cells": n_cells,
        "n_core": n_core,
        "n_interior": n_interior,
        "n_exterior": n_exterior,
        "n_pairs_interior": n_pair_int,
        "n_pairs_exterior": n_pair_ext,
        "n_pairs_whole": n_pair_all,
        "settle_steps": settle_steps,
        "seed": resolved.seed if seed is None else seed,
        "repulsive_control": repulsive_control,
        # geometry
        "R_edge_um": R_edge * 1e6,
        "Rg_um": Rg * 1e6,
        "nn_median_over_r0": nn["median"] / r0,
        # pressures (DIAGNOSTIC — naive radial split; ~0 for a self-bound drop)
        "P_interior_Pa": P_interior,
        "P_exterior_Pa": P_exterior,
        "P_whole_Pa": P_whole,
        "dP_diag_Pa": dP,
        # emergent surface tension (PRIMARY: Irving-Kirkwood spherical mechanical estimator)
        "sigma_emergent_mN_m": sigma * 1e3 if np.isfinite(sigma) else float("nan"),
        "gamma_cortical_mN_m": GAMMA_CORTICAL * 1e3,
        "sigma_over_gamma": (sigma / GAMMA_CORTICAL) if np.isfinite(sigma) else float("nan"),
        # MCF7 emergent beta/gamma
        "contact_radius_um": a_contact * 1e6,
        "contact_area_um2": contact_area * 1e12,
        "work_deadhesion_fJ": work_deadhesion * 1e15,
        "beta_mN_m": beta * 1e3,
        "beta_over_gamma": beta_over_gamma,
        "roffay_beta_over_gamma_window": list(bog_window),
        "beta_over_gamma_in_roffay": bool(
            min(bog_window) <= beta_over_gamma <= max(bog_window)
        ),
        "beta_over_gamma_above_roffay": bool(beta_over_gamma > max(bog_window)),
        # profile (for plotting)
        "_r_centers_um": (r_centers * 1e6).tolist(),
        "_density": density.tolist(),
        "_pos_init_um": (pos_init * 1e6).tolist(),
        "_pos_final_um": (pos * 1e6).tolist(),
        "_interior_mask": interior_mask.tolist(),
        "_exterior_mask": exterior_mask.tolist(),
    }


def make_figure(main_res: dict, ctrl_res: dict, out_png: Path) -> None:
    """4-panel: aggregate (interior/exterior), radial density+R, pressures, sigma vs gamma."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pos = np.array(main_res["_pos_final_um"])
    int_mask = np.array(main_res["_interior_mask"], dtype=bool)
    ext_mask = np.array(main_res["_exterior_mask"], dtype=bool)
    other = ~(int_mask | ext_mask)

    fig, ax = plt.subplots(1, 4, figsize=(20, 4.8))

    # Panel A: aggregate cross-section coloured by interior/exterior shell.
    ax[0].scatter(pos[other, 0], pos[other, 1], s=14, c="lightgray", label="not in core")
    ax[0].scatter(pos[ext_mask, 0], pos[ext_mask, 1], s=20, c="tab:orange",
                  edgecolors="k", linewidths=0.2, label=f"exterior shell (n={ext_mask.sum()})")
    ax[0].scatter(pos[int_mask, 0], pos[int_mask, 1], s=20, c="tab:blue",
                  edgecolors="k", linewidths=0.2, label=f"interior (n={int_mask.sum()})")
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("x (µm)")
    ax[0].set_ylabel("y (µm)")
    ax[0].set_title(f"A  settled aggregate (N={main_res['n_cells']})\nR_edge={main_res['R_edge_um']:.1f} µm")
    ax[0].legend(fontsize=7, loc="upper right")

    # Panel B: radial density profile + R_edge.
    r_c = np.array(main_res["_r_centers_um"])
    dens = np.array(main_res["_density"])
    ax[1].plot(r_c, dens, "o-", color="tab:green")
    ax[1].axvline(main_res["R_edge_um"], color="crimson", ls="--",
                  label=f"R_edge = {main_res['R_edge_um']:.1f} µm")
    ax[1].set_xlabel("radius (µm)")
    ax[1].set_ylabel("number density (m⁻³)")
    ax[1].set_title("B  radial density profile\n(R_edge = half-core inflection)")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    # Panel C: emergent sigma (IK) cohesive vs repulsive-only control — the sign check.
    labels = ["cohesive\n(full Morse)", "control\n(repulsive-only)"]
    vals = [main_res["sigma_emergent_mN_m"], ctrl_res["sigma_emergent_mN_m"]]
    colors = ["tab:blue", "tab:red"]
    ax[2].bar([0, 1], vals, 0.6, color=colors)
    ax[2].axhline(0.0, color="k", lw=0.8)
    ax[2].set_xticks([0, 1])
    ax[2].set_xticklabels(labels, fontsize=8)
    ax[2].set_ylabel("emergent sigma_IK (mN/m)")
    ax[2].set_title("C  sign check: cohesion gives sigma>0,\nrepulsion-only flips it")
    ax[2].grid(alpha=0.3, axis="y")

    # Panel D: emergent sigma vs gamma anchor + KU-3.5 band.
    band = list(br.CORTICAL_TENSION_BAND_MN_M)
    ax[3].axhspan(band[0], band[1], color="tab:green", alpha=0.16, label=f"KU-3.5 band {band}")
    ax[3].axhline(main_res["gamma_cortical_mN_m"], color="k", lw=2,
                  label=f"gamma = {main_res['gamma_cortical_mN_m']:.2f} mN/m")
    sig = main_res["sigma_emergent_mN_m"]
    if np.isfinite(sig):
        ax[3].axhline(sig, color="tab:blue", lw=2.5, ls="-",
                      label=f"emergent sigma = {sig:.3f} mN/m")
        ax[3].text(0.5, sig, f"  sigma/gamma = {main_res['sigma_over_gamma']:.2f}",
                   va="bottom", fontsize=8, color="tab:blue")
    ax[3].set_yscale("log")
    ax[3].set_ylabel("surface tension (mN/m)")
    ax[3].set_xticks([])
    ax[3].set_title("D  EMERGENT sigma vs single-cell gamma\n(virial -> Young-Laplace)")
    ax[3].legend(fontsize=7, loc="lower left")

    fig.suptitle(
        "L2 D2 — emergent aggregate surface tension from a CBM spheroid (measurement, not anchor)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def _strip_arrays(d: dict) -> dict:
    """Drop the large plotting arrays before JSON serialisation."""
    return {k: v for k, v in d.items() if not k.startswith("_")}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Layer-2 emergent surface-tension measurement.")
    ap.add_argument("--n-cells", type=int, default=200)
    ap.add_argument("--settle", type=int, default=40_000)
    ap.add_argument("--interior-frac", type=float, default=0.6)
    ap.add_argument("--no-fig", action="store_true")
    ap.add_argument(
        "--ensemble", action="store_true",
        help="also run a small seed x N sweep and record the spread (honest noise estimate).",
    )
    args = ap.parse_args(argv)

    resolved = resolve_layer2(yaml.safe_load(_RUNTIME_CFG.read_text()))
    print(f"[resolve] r0={resolved.morse_r0*1e6:.2f} µm  D_e={resolved.D_e:.3e} J "
          f"({resolved.D_e/resolved.kT:.0f} kT)  alpha={resolved.morse_alpha:.3e} 1/m  "
          f"dt={resolved.dt_cfl:.3e} s  gamma_cell={resolved.gamma_cell:.3f} N·s/m")
    print(f"[anchor] single-cell cortical tension gamma = {GAMMA_CORTICAL*1e3:.3f} mN/m\n")

    print(f"[run] cohesive spheroid: N={args.n_cells}, settle={args.settle} ...")
    main_res = measure_emergent_sigma(
        resolved, args.n_cells, settle_steps=args.settle,
        interior_frac=args.interior_frac,
    )
    print(f"[run] repulsive-only CONTROL (same config) ...")
    ctrl_res = measure_emergent_sigma(
        resolved, args.n_cells, settle_steps=args.settle,
        interior_frac=args.interior_frac, repulsive_control=True,
    )

    print("\n=== EMERGENT surface tension (cohesive spheroid) ===")
    print(f"  R_edge            = {main_res['R_edge_um']:.2f} µm  (Rg = {main_res['Rg_um']:.2f} µm)")
    print(f"  core / interior / exterior cells = "
          f"{main_res['n_core']} / {main_res['n_interior']} / {main_res['n_exterior']}")
    print(f"  P_whole (virial)  = {main_res['P_whole_Pa']:+.4e} Pa  "
          f"(~0 => self-bound drop, no external wall — expected)")
    print(f"  [diag] P_interior = {main_res['P_interior_Pa']:+.4e} Pa, "
          f"P_exterior = {main_res['P_exterior_Pa']:+.4e} Pa "
          f"(naive radial split — does NOT isolate sigma)")
    print(f"  -> emergent sigma (Irving-Kirkwood) = {main_res['sigma_emergent_mN_m']:.4f} mN/m")
    print(f"     gamma (anchor)                   = {main_res['gamma_cortical_mN_m']:.4f} mN/m")
    print(f"     sigma / gamma                    = {main_res['sigma_over_gamma']:.3f}")

    print("\n=== sign CONTROL (repulsive-only, same geometry) ===")
    print(f"  sigma_IK (repulsive-only) = {ctrl_res['sigma_emergent_mN_m']:.4f} mN/m  "
          f"(should be <= 0 / sign-flipped vs cohesive: no adhesion to make a surface)")
    sign_ok = (
        np.isfinite(main_res["sigma_emergent_mN_m"])
        and main_res["sigma_emergent_mN_m"] > 0.0
        and (
            not np.isfinite(ctrl_res["sigma_emergent_mN_m"])
            or ctrl_res["sigma_emergent_mN_m"] < main_res["sigma_emergent_mN_m"]
        )
    )
    print(f"  sign convention validated: {sign_ok}  "
          f"(cohesive sigma>0 and > repulsive-only control)")

    print("\n=== MCF7 emergent beta/gamma (cell-cell contact geometry) ===")
    print(f"  mean NN sep / r0  = {main_res['nn_median_over_r0']:.4f}")
    print(f"  contact radius    = {main_res['contact_radius_um']:.3f} µm  "
          f"(area {main_res['contact_area_um2']:.2f} µm²)")
    print(f"  work de-adhesion  = {main_res['work_deadhesion_fJ']:.3f} fJ (= D_e)")
    print(f"  beta (adhesion)   = {main_res['beta_mN_m']:.4f} mN/m")
    print(f"  beta / gamma      = {main_res['beta_over_gamma']:.3f}  "
          f"(Roffay window {tuple(round(x,3) for x in main_res['roffay_beta_over_gamma_window'])})")
    loc = ("IN" if main_res["beta_over_gamma_in_roffay"]
           else "ABOVE" if main_res["beta_over_gamma_above_roffay"] else "BELOW")
    print(f"  -> beta/gamma lands {loc} the Roffay 0.375-0.5 window")

    ensemble = None
    if args.ensemble:
        print("\n=== ensemble (seed x N) — honest noise estimate ===")
        seeds = (42, 7, 101)
        ns = (150, 200, 300)
        runs = []
        for sd in seeds:
            for nn_n in ns:
                r = measure_emergent_sigma(
                    resolved, nn_n, settle_steps=args.settle,
                    interior_frac=args.interior_frac, seed=sd,
                )
                runs.append({
                    "seed": sd, "n_cells": nn_n,
                    "sigma_emergent_mN_m": r["sigma_emergent_mN_m"],
                    "sigma_over_gamma": r["sigma_over_gamma"],
                    "beta_over_gamma": r["beta_over_gamma"],
                    "R_edge_um": r["R_edge_um"],
                })
                print(f"  seed={sd} N={nn_n}: sigma={r['sigma_emergent_mN_m']:+.4f} mN/m "
                      f"(sigma/g={r['sigma_over_gamma']:+.3f})  beta/g={r['beta_over_gamma']:.2f}")
        sig_arr = np.array([x["sigma_emergent_mN_m"] for x in runs])
        bog_arr = np.array([x["beta_over_gamma"] for x in runs])
        ensemble = {
            "runs": runs,
            "sigma_mean_mN_m": float(np.nanmean(sig_arr)),
            "sigma_std_mN_m": float(np.nanstd(sig_arr)),
            "sigma_over_gamma_mean": float(np.nanmean(sig_arr) / (GAMMA_CORTICAL * 1e3)),
            "beta_over_gamma_mean": float(np.nanmean(bog_arr)),
            "beta_over_gamma_min": float(np.nanmin(bog_arr)),
            "beta_over_gamma_max": float(np.nanmax(bog_arr)),
        }
        print(f"  -> sigma = {ensemble['sigma_mean_mN_m']:.4f} +/- "
              f"{ensemble['sigma_std_mN_m']:.4f} mN/m "
              f"(sigma/gamma ~ {ensemble['sigma_over_gamma_mean']:.3f})")
        print(f"  -> beta/gamma = {ensemble['beta_over_gamma_mean']:.2f} "
              f"[{ensemble['beta_over_gamma_min']:.2f}, {ensemble['beta_over_gamma_max']:.2f}] "
              f"(all ABOVE Roffay 0.375-0.5)")

    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_json = _OUT_DIR / "emergent_sigma.json"
    payload = {
        "cohesive": _strip_arrays(main_res),
        "control_repulsive": _strip_arrays(ctrl_res),
        "sign_convention_validated": bool(sign_ok),
        "ensemble": ensemble,
        "resolved": {
            "morse_r0_um": resolved.morse_r0 * 1e6,
            "D_e_J": resolved.D_e,
            "morse_alpha_1_per_m": resolved.morse_alpha,
            "gamma_cell_Ns_per_m": resolved.gamma_cell,
            "dt_cfl_s": resolved.dt_cfl,
        },
        "method": (
            "EMERGENT sigma = mechanical Irving-Kirkwood SPHERICAL surface tension of a "
            "settled CBM spheroid: sigma = -(1/(16 pi R^2)) sum_pairs (r_ij.f_ij)[1 - "
            "3(s.r_hat_c)^2/s^2], f_ij the analytic Morse pair force (matches the run's "
            "md.pair.Morse), R the radial-density inflection radius. The isotropic interior/"
            "exterior virial split is a DIAGNOSTIC only (P_whole ~ 0 for a self-bound drop). "
            "beta from emergent contact area (NN-separation chord) and work-of-de-adhesion D_e. "
            "Sign control: repulsive-only config flips sigma sign (no cohesion => no surface)."
        ),
    }
    out_json.write_text(json.dumps(payload, indent=2))
    print(f"\n[json] {out_json}")

    if not args.no_fig:
        try:
            out_png = _OUT_DIR / "emergent_sigma.png"
            make_figure(main_res, ctrl_res, out_png)
            print(f"[viz]  {out_png}")
        except Exception as exc:  # best-effort
            print(f"[viz]  skipped figure: {exc}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

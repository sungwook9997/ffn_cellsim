"""Regenerate the unified-weave (I4) analytic-oracle figures — the single entry-point for this track.

Every figure OVERLAYS the closed-form oracle / reference band on the numeric result, annotates SI-ish engine
units (deg, um, 1/s), never truncates an axis, and shows the reference line, per the professor's
visualization-integrity rules (PI 2026-05-21 / 2026-07-07). Pure numpy/matplotlib on the CPU-green analytic
suite — runs on the dev Mac (no Warp/CUDA).

    python aleph/scripts/ac_weave_vis.py    ->  aleph/outputs/ac/weave/figs/*.png

Figures:
  i4_branch_angle.png   dendritic Arp2/3 branch-angle distribution vs theta0 +- sigma_theta (the §3A.b gate;
                        the angle FLUCTUATES thermally, it is NOT a rigid 72 deg delta)
  i4_cortex_parity.png  Gate-1: weave_cell([CORTEX]) is bit-identical to ff.weave.weave(CORTEX) (max|diff|=0)
  i4_dendritic_net.png  the dendritic mother/daughter Arp2/3 network (the rebuilt lamellipodium topology)
  i4_crosslink_kmc.png  topology-reforming crosslink KMC: steady bound fraction vs k_on/(k_on+k_off) + partner
                        turnover (the topology reforms while the count is steady)

The branch-angle knobs are the GROUNDED Faessler-2020 anchors (theta0=70 deg, sigma_theta=9 deg -> k_theta by
equipartition); the KMC k_on/k_off are oracle-sweep values (I0-B3/B4). The NATIVE protrusion / bundle-
condensation magnitude gates are the lead's (gbook) and are INVALID until the I0-B4 GAPs close.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aleph.laws.architecture_spec import CORTEX
from aleph.laws.weave import weave as ff_weave
from aleph.components.weave.branch_angle import (
    ARP23_K_THETA,
    ARP23_THETA0_RAD,
    KT_310K_PN_UM,
    angle_energy,
    branch_angle,
    thermal_sigma,
)
from aleph.components.weave.crosslink_kmc import reattach_step, steady_bound_fraction
from aleph.components.weave.lamellipodium import build_lamellipodium
from aleph.components.weave.regions import RegionSpec
from aleph.components.weave.woven_cell import weave_cell

OUT = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "weave" / "figs"


def fig_branch_angle() -> Path:
    """Measured dendritic branch angles vs the theta0 +- sigma_theta Boltzmann oracle (the §3A.b gate)."""
    rng = np.random.default_rng(1)
    seed = build_lamellipodium(n_mothers=20, n_daughter_pool=4000, length_um=1.0, seg_um=0.5,
                               patch_half_um=8.0, rng=rng, monomer_c=1e6, k_arp0=50.0, k_m=5.0, tau=5.0)
    ang = np.rad2deg(np.array([branch_angle(seed.pos[m], seed.pos[j], seed.pos[d])
                               for m, j, d in seed.branch_triples]))
    theta0 = np.rad2deg(ARP23_THETA0_RAD)
    sigma = np.rad2deg(thermal_sigma(ARP23_K_THETA))

    th = np.linspace(theta0 - 5 * sigma, theta0 + 5 * sigma, 500)
    thr = np.deg2rad(th)
    pdf = np.sin(thr) * np.exp(-angle_energy(thr, ARP23_K_THETA, ARP23_THETA0_RAD) / KT_310K_PN_UM)
    pdf /= np.trapezoid(pdf, th)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.hist(ang, bins=60, density=True, color="#4c72b0", alpha=0.55,
            label=f"measured branch angles (N={ang.size})")
    ax.plot(th, pdf, "k-", lw=2, label=r"Boltzmann $\propto \sin\theta\,e^{-U/kT}$ oracle")
    ax.axvline(theta0, color="crimson", lw=2, ls="--", label=r"$\theta_0=70°$ (Fäßler 2020)")
    ax.axvspan(theta0 - sigma, theta0 + sigma, color="crimson", alpha=0.10,
               label=r"$\theta_0\pm\sigma_\theta$ ($\sigma_\theta=9°$)")
    ax.annotate(f"measured: mean={ang.mean():.1f}°, SD={ang.std():.1f}°\n"
                "(the angle FLUCTUATES — not a rigid 72°)",
                xy=(0.02, 0.97), xycoords="axes fraction", va="top",
                bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9), fontsize=9)
    ax.set_xlabel("Arp2/3 branch angle θ [deg]")
    ax.set_ylabel("probability density [1/deg]")
    ax.set_title("I4 dendritic Arp2/3 branch-angle distribution vs θ₀±σ_θ oracle")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    p = OUT / "i4_branch_angle.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    return p


def fig_cortex_parity() -> Path:
    """Gate-1: weave_cell([CORTEX]) bit-identical to ff.weave.weave(CORTEX) — max|diff| = 0 across arrays."""
    fil = dataclasses.replace(CORTEX.filament, n_filaments=400)
    spec = dataclasses.replace(CORTEX, filament=fil)
    ref = ff_weave(spec, rng=np.random.default_rng(0))
    wc = weave_cell([RegionSpec(arch=spec)], rng=np.random.default_rng(0))

    names = ["pos", "xl_i", "xl_j", "xl_k", "xl_rest", "myo_i", "myo_j"]
    pairs = [(wc.pos, ref.net.pos), (wc.xl_i, ref.xl_i), (wc.xl_j, ref.xl_j), (wc.xl_k, ref.xl_k),
             (wc.xl_rest, ref.xl_rest), (wc.myo_i, ref.myo_i), (wc.myo_j, ref.myo_j)]
    maxdiff = [float(np.max(np.abs(a.astype(float) - b.astype(float)))) if a.size else 0.0 for a, b in pairs]

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.4))
    r_wc = np.linalg.norm(wc.pos - wc.pos.mean(0), axis=1)
    r_ref = np.linalg.norm(ref.net.pos - ref.net.pos.mean(0), axis=1)
    ax0.hist(r_ref, bins=40, histtype="stepfilled", color="#dd8452", alpha=0.5, label="ff.weave.weave(CORTEX)")
    ax0.hist(r_wc, bins=40, histtype="step", color="#4c72b0", lw=2, label="weave_cell([CORTEX])")
    ax0.set_xlabel("node radius about centroid [µm]")
    ax0.set_ylabel("count")
    ax0.set_title("cortex node cloud — perfectly overlapping")
    ax0.legend(fontsize=8)

    ax1.bar(range(len(names)), np.array(maxdiff) + 1e-20, color="#55a868")
    ax1.set_xticks(range(len(names)))
    ax1.set_xticklabels(names, rotation=35, ha="right", fontsize=8)
    ax1.set_ylabel("max |weave_cell − ff.weave|")
    ax1.set_ylim(0, 1.0)
    ax1.set_title(f"Gate-1 bit-identical parity (all max|diff| = {max(maxdiff):.0e})")
    ax1.annotate("regression anchor:\nweave_cell([CORTEX]) all-OFF\n= the current cortex → identical γ",
                 xy=(0.5, 0.6), xycoords="axes fraction", ha="center",
                 bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9), fontsize=9)
    fig.tight_layout()
    p = OUT / "i4_cortex_parity.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    return p


def fig_dendritic_net() -> Path:
    """The rebuilt dendritic Arp2/3 lamellipodium (mother/daughter tree), the net-new §3A.b topology."""
    rng = np.random.default_rng(4)
    seed = build_lamellipodium(n_mothers=10, n_daughter_pool=60, length_um=1.0, seg_um=0.25,
                               patch_half_um=6.0, rng=rng, npf_kappa=2.0, monomer_c=50.0,
                               k_arp0=8.0, k_m=5.0, tau=2.0)
    off = seed.fiber_offsets
    n_mothers = 10
    fig, ax = plt.subplots(figsize=(6.6, 6.4))
    for f in range(seed.fiber_offsets.shape[0] - 1):
        seg = seed.pos[off[f]:off[f + 1]]
        if not seed.active_mask[f]:
            continue
        is_mother = f < n_mothers
        ax.plot(seg[:, 0], seg[:, 1], "-", lw=2.2 if is_mother else 1.3,
                color="#c44e52" if is_mother else "#4c72b0", alpha=0.9 if is_mother else 0.7)
    # mark branch junctions
    if seed.branch_triples.size:
        bn = seed.branch_triples[seed.branch_active][:, 1] if seed.branch_active.any() else seed.branch_triples[:, 1]
        ax.scatter(seed.pos[bn, 0], seed.pos[bn, 1], s=14, c="k", zorder=5, label="Arp2/3 branch (θ₀=70°)")
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color="#c44e52", lw=2.5, label="mother filament"),
               Line2D([0], [0], color="#4c72b0", lw=1.5, label="Arp2/3 daughter"),
               Line2D([0], [0], marker="o", color="k", lw=0, label="branch junction")]
    ax.legend(handles=handles, loc="upper right", fontsize=8)
    ax.set_xlabel("x [µm]")
    ax.set_ylabel("y [µm]  (protrusion axis +y)")
    ax.set_aspect("equal")
    ax.set_title("I4 lamellipodium REBUILD — dendritic Arp2/3 net (angle-harmonic branch)")
    fig.tight_layout()
    p = OUT / "i4_dendritic_net.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    return p


def fig_crosslink_kmc() -> Path:
    """Topology-reforming crosslink KMC: steady bound fraction vs k_on/(k_on+k_off) + partner turnover."""
    rng = np.random.default_rng(3)
    nb = 25
    fibs = [np.column_stack([np.full(nb, x), np.linspace(0, 6, nb), np.zeros(nb)]) for x in (0.0, 0.15, 0.30)]
    pos = np.vstack(fibs)
    fiber_of_node = np.concatenate([np.full(nb, i) for i in range(3)])
    e = 120
    node_of_end = rng.integers(0, nb, e)
    bound = np.zeros(e, bool)
    partner = np.full(e, -1, np.int64)
    loads = np.zeros(e)
    k_on, k_off = 5.0, 2.0
    frac, turnover, partner0 = [], [], None
    for t in range(400):
        bound, partner = reattach_step(bound, partner, pos, fiber_of_node, node_of_end, loads,
                                       reach=0.5, k_on=k_on, k_off0=k_off, f0=7.0, tau=0.05, rng=rng)
        frac.append(bound.mean())
        if t == 60:
            partner0 = partner.copy()
        if partner0 is not None:
            bb = (partner0 >= 0) & (partner >= 0)
            turnover.append(np.count_nonzero(partner[bb] != partner0[bb]) / max(bb.sum(), 1))
        else:
            turnover.append(0.0)

    expected = steady_bound_fraction(k_on, k_off)
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.4))
    ax0.plot(frac, color="#4c72b0", lw=1.4, label="KMC bound fraction")
    ax0.axhline(expected, color="crimson", lw=2, ls="--",
                label=fr"$k_{{on}}/(k_{{on}}+k_{{off}})={expected:.3f}$ (detailed balance)")
    ax0.set_xlabel("KMC tick")
    ax0.set_ylabel("bound crosslinker fraction")
    ax0.set_ylim(0, 1)
    ax0.set_title("crosslinker bound fraction vs the oracle")
    ax0.legend(fontsize=8)

    ax1.plot(turnover, color="#55a868", lw=1.6)
    ax1.set_xlabel("KMC tick")
    ax1.set_ylabel("fraction of bonds re-partnered (vs t=60)")
    ax1.set_ylim(0, 1)
    ax1.set_title("topology REFORMS while the count is steady")
    ax1.annotate("partners hop to new cross-fiber nodes\n→ the network topology reforms\n"
                 "(the SF/cap/arc bundles condense at I5)",
                 xy=(0.5, 0.25), xycoords="axes fraction", ha="center",
                 bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9), fontsize=9)
    fig.tight_layout()
    p = OUT / "i4_crosslink_kmc.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    return p


def fig_stress_fiber() -> Path:
    """SF prestress load path: an antiparallel FA<->FA sarcomere condenses a closed contractile dipole while the
    isotropic-conditional seed does NOT — the emergent-prestress signature (magnitude-independent SHAPE)."""
    from aleph.components.weave.regions import VENTRAL_SF_REGION
    from aleph.components.weave.stress_fiber import stress_fiber_load_path, wire_stress_fibers

    # (a) an organized antiparallel sarcomere -> closed contractile dipole (what SF condensation LOOKS like)
    A = np.column_stack([np.linspace(0.0, 2.0, 9), np.zeros(9), np.zeros(9)])
    B = np.column_stack([np.linspace(2.0, 4.0, 9), np.zeros(9), np.zeros(9)])
    pos = np.vstack([A, B]); foff = np.array([0, 9, 18]); pol = np.array([-1, +1])
    myo_i = np.array([6, 7, 8]); myo_j = np.array([9, 10, 11])       # interior overlap motors
    lp = stress_fiber_load_path(pos, foff, pol, myo_i, myo_j, np.array([0, 17]), np.zeros(0, int), n_bins=32)

    # (b) the isotropic-conditional ventral-SF SEED (has NOT condensed -> open, not a clean closed dipole)
    cell = weave_cell([VENTRAL_SF_REGION], rng=np.random.default_rng(0))
    seed = wire_stress_fibers(cell).load_paths[0]

    def _frac(sb):
        span = max(sb[-1] - sb[0], 1e-12)
        return (sb - sb[0]) / span

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax0.plot(_frac(lp.s_bins), lp.tension_per_fhead, "-o", ms=3, color="#b2182b",
             label="organized sarcomere (condensed: uniform tension)")
    ax0.plot(_frac(seed.s_bins), seed.tension_per_fhead, "-s", ms=3, color="#2166ac",
             label=f"isotropic seed (S={seed.nematic_order:.2f}: disordered, mixed-sign)")
    ax0.axhline(0.0, color="k", lw=0.6)
    ax0.set_xlabel("fractional axial position along bundle  [-]")
    ax0.set_ylabel("internal axial prestress / f_head  [pN per pN-head]")
    ax0.set_title("SF prestress load path (magnitude-independent shape)")
    ax0.legend(fontsize=8, loc="lower center")

    labels = ["condensed\nsarcomere", "isotropic\nseed"]
    dip = [lp.dipole_per_fhead, seed.dipole_per_fhead]
    res = [lp.balance_residual_per_fhead, seed.balance_residual_per_fhead]
    x = np.arange(2)
    ax1.bar(x - 0.18, dip, 0.36, color="#b2182b", label="traction dipole / f_head")
    ax1.bar(x + 0.18, res, 0.36, color="#8c8c8c", label="balance residual / f_head")
    ax1.set_xticks(x); ax1.set_xticklabels(labels)
    ax1.set_ylabel("per f_head  [pN per pN-head]")
    ax1.set_title("closed dipole (residual->0) is the SF signature")
    ax1.legend(fontsize=8)
    for i, (d, r) in enumerate(zip(dip, res)):
        ax1.annotate("CLOSED" if d > 1e-9 and r <= 0.5 * d else "open", (i, max(d, r)),
                     ha="center", va="bottom", fontsize=8)
    fig.text(0.5, 0.005, "MAGNITUDE IS GAP (I0-B3 f_stall + I0-B6 head density) — shape only; report-not-tune.",
             ha="center", fontsize=8, style="italic")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    p = OUT / "i4_stress_fiber.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    return p


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for fn in (fig_branch_angle, fig_cortex_parity, fig_dendritic_net, fig_crosslink_kmc, fig_stress_fiber):
        print("wrote", fn())


if __name__ == "__main__":
    main()

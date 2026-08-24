"""Regenerate the emergence-detector (I5 detector-core) analytic-oracle figures.

The single regeneration entry-point for the emergence track (extend it as the deferred native emergence
proof lands). Every figure OVERLAYS the closed-form null band / oracle / pre-registered bar on the
measured detector output — never the measurement alone — annotates axes (S, z, fractions are
dimensionless [-]), and never truncates an axis (the professor's visualization-integrity rules). Pure
NumPy/SciPy/matplotlib on the synthetic oracle configs — runs on the dev Mac (no Warp/CUDA).

    python aleph/scripts/ac_emergence_vis.py    ->  aleph/outputs/ac/emergence/figs/*.png

Figures (per the parallel-sessions per-track spec):
  * i5_null_vs_aligned   — S~0 finite-N isotropic null band vs S=1 aligned, + the 3/(2N) analytic anchor;
  * i5_localization      — a planted bundle CONDENSING out of an isotropic background (local-S map + recovery);
  * i5_invariance        — rotational + permutation/label invariance of S and of the recovered centroid;
  * i5_effect_size       — the pre-registered effect-size z rule: dose-response + specificity vs Z_CRIT.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aleph.components.emergence import synthetic as syn
from aleph.components.emergence.condensation import compute_local_fields, localize_bundles
from aleph.components.emergence.nematic import fiber_axes, nematic_order
from aleph.components.emergence.null_model import (
    EFFECT_SIZE_Z_CRIT,
    effect_size,
    nematic_null_band,
    nematic_null_scale_sq,
)

OUT = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "emergence" / "figs"


def fig_null_vs_aligned() -> Path:
    """S~0 finite-N isotropic null band vs S=1 aligned + the analytic 3/(2N) sum-of-squares anchor."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ns = np.array([25, 50, 100, 200, 400, 800, 1600, 3200])
    means, stds = [], []
    for n in ns:
        band = nematic_null_band(int(n), n_mc=2500, seed=17)
        means.append(band.mean); stds.append(band.std)
    means, stds = np.array(means), np.array(stds)

    ax0.fill_between(ns, means - stds, means + stds, alpha=0.25, color="steelblue",
                     label="isotropic null band  μ ± σ")
    ax0.plot(ns, means, "o-", color="steelblue", lw=2, label="null mean  μ_null(N)")
    # per-realisation isotropic draws (thin markers) land in the band
    rng = np.random.default_rng(0)
    for n in ns:
        s_draws = [nematic_order(fiber_axes(*syn.isotropic_network(int(n), seed=int(rng.integers(1e6)))))
                   for _ in range(8)]
        ax0.scatter(np.full(len(s_draws), n), s_draws, s=9, color="gray", alpha=0.5, zorder=1)
    # 1/sqrt(N) scaling guide anchored at the largest-N null mean
    guide = means[-1] * np.sqrt(ns[-1] / ns)
    ax0.plot(ns, guide, "--", color="navy", lw=1.2, label="1/√N scaling (analytic bias law)")
    # aligned corner: S = 1 for every N
    s_aligned = [nematic_order(fiber_axes(*syn.aligned_network(int(n), seed=1))) for n in ns]
    ax0.plot(ns, s_aligned, "s-", color="crimson", lw=2, label="aligned bundle  S = 1")
    ax0.set_xscale("log"); ax0.set_ylim(0, 1.05)
    ax0.set_xlabel("fiber count  N  [–]"); ax0.set_ylabel("nematic order  S = λ_max(Q)  [–]")
    ax0.set_title("Isotropic null (S~0) vs aligned (S=1)"); ax0.legend(loc="center right", fontsize=8)

    # analytic cross-check: MC E[Σλ²] overlaid on the closed-form 3/(2N)
    mc = [nematic_null_band(int(n), n_mc=2500, seed=17).mc_sumsq for n in ns]
    analytic_sumsq = np.array([nematic_null_scale_sq(int(n)) for n in ns])  # = 3/(2N)
    ax1.plot(ns, analytic_sumsq, "-", color="black", lw=2, label="analytic  E[Σλ²] = 3/(2N)")
    ax1.scatter(ns, mc, color="crimson", zorder=5, label="Monte-Carlo  ⟨Σλ²⟩")
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("fiber count  N  [–]"); ax1.set_ylabel("Σλ²(Q)  [–]")
    ax1.set_title("Null cross-check: MC vs closed form 3/(2N)"); ax1.legend(loc="upper right")
    fig.suptitle("I5 emergence — nematic order corners + analytic null anchor")
    fig.tight_layout()
    pth = OUT / "i5_null_vs_aligned.png"; fig.savefig(pth, dpi=130); plt.close(fig)
    return pth


def fig_localization() -> Path:
    """A planted bundle condensing out of an isotropic background: local-S map + label-blind recovery."""
    pos, off, mask, center = syn.planted_bundle(400, 40, seed=3, director=(0.0, 0.0, 1.0))
    fields = compute_local_fields(pos, off)
    bundles = localize_bundles(fields)
    top = max(bundles, key=lambda b: b.size)
    rec = set(top.members.tolist()); planted = set(np.flatnonzero(mask).tolist())
    tp = len(rec & planted)
    precision, recall = tp / len(rec), tp / len(planted)
    cen_err = float(np.linalg.norm(top.center - center))
    xz = fields.centroids[:, [0, 2]]  # project onto the (x, z=director) plane

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11.5, 5.0))
    sc = ax0.scatter(xz[:, 0], xz[:, 1], c=fields.s_local, s=16, cmap="viridis", vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax0, label="local nematic order  S_local  [–]")
    ax0.set_title("Local order field (label-blind): the bundle patch lights up")
    ax0.set_xlabel("x  [µm]"); ax0.set_ylabel("z (director)  [µm]"); ax0.set_aspect("equal")

    ax1.scatter(xz[~mask, 0], xz[~mask, 1], s=12, color="lightgray", label="isotropic background (truth)")
    ax1.scatter(xz[mask, 0], xz[mask, 1], s=18, color="steelblue", label="planted bundle (truth)")
    mem = top.members
    ax1.scatter(xz[mem, 0], xz[mem, 1], s=60, facecolors="none", edgecolors="crimson", lw=1.2,
                label="detector-recovered members")
    ax1.scatter([center[0]], [center[2]], marker="x", s=140, color="black", label="planted centre")
    ax1.scatter([top.center[0]], [top.center[2]], marker="*", s=180, color="crimson",
                label="recovered centre")
    ax1.set_title(f"Localization: recall={recall:.2f} precision={precision:.2f} |Δc|={cen_err:.2f} µm")
    ax1.set_xlabel("x  [µm]"); ax1.set_ylabel("z (director)  [µm]"); ax1.set_aspect("equal")
    ax1.legend(loc="upper left", fontsize=7)
    fig.suptitle("I5 emergence — planted-bundle condensation & label-blind localization")
    fig.tight_layout()
    pth = OUT / "i5_localization.png"; fig.savefig(pth, dpi=130); plt.close(fig)
    return pth


def fig_invariance() -> Path:
    """Rotational + permutation/label invariance of S; rotational equivariance of the recovered centroid."""
    pos, off = syn.partially_aligned_network(400, aligned_frac=0.4, seed=6)
    s_ref = nematic_order(fiber_axes(pos, off))
    rng = np.random.default_rng(0)
    n_trials = 40
    s_rot = np.array([nematic_order(fiber_axes(syn.rotate_config(pos, syn.random_rotation(rng)), off))
                      for _ in range(n_trials)])
    s_perm = []
    for _ in range(n_trials):
        perm = rng.permutation(off.shape[0] - 1)
        pp, oo = syn.permute_fibers(pos, off, perm)
        s_perm.append(nematic_order(fiber_axes(pp, oo)))
    s_perm = np.array(s_perm)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax0.axhline(s_ref, color="black", lw=2, label=f"reference S = {s_ref:.4f}")
    ax0.plot(np.abs(s_rot - s_ref), "-", color="steelblue", alpha=0.4, lw=0.8)
    ax0.scatter(range(n_trials), np.abs(s_rot - s_ref), s=14, color="steelblue", label="rotations |ΔS|")
    ax0.plot(np.abs(s_perm - s_ref), "-", color="crimson", alpha=0.4, lw=0.8)
    ax0.scatter(range(n_trials), np.abs(s_perm - s_ref), s=14, color="crimson", label="permutations |ΔS|")
    ax0.set_yscale("log"); ax0.set_ylim(1e-16, 1e-1)
    ax0.axhline(1e-12, ls=":", color="gray", label="machine-precision floor")
    ax0.set_xlabel("trial  [–]"); ax0.set_ylabel("|S − S_ref|  [–]  (log)")
    ax0.set_title("Rotational + label invariance of S (|ΔS| ~ 1e-13)")
    ax0.legend(loc="upper right", fontsize=8)

    # localization equivariance: recovered centre rotates by the same R
    pos_b, off_b, mask_b, center_b = syn.planted_bundle(400, 40, seed=3)
    top0 = max(localize_bundles(compute_local_fields(pos_b, off_b)), key=lambda b: b.size)
    errs = []
    for _ in range(8):
        r = syn.random_rotation(rng)
        pr = syn.rotate_config(pos_b, r, about=center_b)
        topr = max(localize_bundles(compute_local_fields(pr, off_b)), key=lambda b: b.size)
        expected = center_b + (top0.center - center_b) @ r.T
        errs.append(float(np.linalg.norm(topr.center - expected)))
    ax1.scatter(range(len(errs)), errs, s=30, color="seagreen")
    ax1.plot(errs, "-", color="seagreen", alpha=0.4)
    ax1.axhline(0.0, color="black", lw=1)
    ax1.set_ylim(-1e-6, max(1e-6, 2 * (max(errs) if errs else 1e-6)))
    ax1.set_xlabel("rotation trial  [–]")
    ax1.set_ylabel("|recovered centre − R·reference|  [µm]")
    ax1.set_title("Localization is rotation-equivariant (centre error ~ 0)")
    fig.suptitle("I5 emergence — invariance / equivariance (anti-coupling firewall property)")
    fig.tight_layout()
    pth = OUT / "i5_invariance.png"; fig.savefig(pth, dpi=130); plt.close(fig)
    return pth


def fig_effect_size() -> Path:
    """The pre-registered effect-size z rule: dose-response crossing Z_CRIT + isotropic specificity."""
    band = nematic_null_band(2000, n_mc=2500, seed=21)
    fracs = np.linspace(0.0, 0.6, 13)
    z = np.array([effect_size(nematic_order(fiber_axes(
        *syn.partially_aligned_network(2000, aligned_frac=p, seed=22))), band) for p in fracs])

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax0.axhspan(-3, EFFECT_SIZE_Z_CRIT, color="lightgray", alpha=0.5, label="not ordered (z < Z_crit)")
    ax0.axhline(EFFECT_SIZE_Z_CRIT, color="crimson", lw=2, ls="--",
                label=f"pre-registered Z_crit = {EFFECT_SIZE_Z_CRIT:g} (5σ)")
    ax0.plot(fracs, z, "o-", color="steelblue", lw=2, label="effect size  z(aligned fraction)")
    ax0.set_xlabel("aligned (condensation) fraction  p  [–]"); ax0.set_ylabel("effect size  z  [–]")
    ax0.set_title("Dose-response: order clears the pre-registered bar"); ax0.legend(loc="upper left", fontsize=8)

    # specificity vs sensitivity: z distributions
    iso_z = np.array([effect_size(nematic_order(fiber_axes(*syn.isotropic_network(400, seed=1000 + k))),
                                  nematic_null_band(400, n_mc=2000, seed=13)) for k in range(60)])
    ali_z = np.array([effect_size(nematic_order(fiber_axes(*syn.aligned_network(400, seed=200 + k, jitter_deg=8.0))),
                                  nematic_null_band(400, n_mc=2000, seed=13)) for k in range(12)])
    ax1.axvline(EFFECT_SIZE_Z_CRIT, color="crimson", lw=2, ls="--", label=f"Z_crit = {EFFECT_SIZE_Z_CRIT:g}")
    ax1.hist(iso_z, bins=20, color="steelblue", alpha=0.7, label="isotropic draws (specificity)")
    # aligned z are enormous; show them as a rug at the right edge with an annotation
    ax1.scatter(np.clip(ali_z, None, ax1.get_xlim()[1]), np.full(ali_z.size, 1.0), marker="^",
                color="crimson", zorder=5, label=f"aligned draws  z≈{ali_z.mean():.0f} (off-scale →)")
    ax1.set_xlabel("effect size  z  [–]"); ax1.set_ylabel("count  [–]")
    ax1.set_title("Specificity: isotropic z stays below Z_crit"); ax1.legend(loc="upper right", fontsize=8)
    fig.suptitle("I5 emergence — pre-registered effect-size / falsifiability rule")
    fig.tight_layout()
    pth = OUT / "i5_effect_size.png"; fig.savefig(pth, dpi=130); plt.close(fig)
    return pth


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    figs = [fig_null_vs_aligned(), fig_localization(), fig_invariance(), fig_effect_size()]
    for f in figs:
        print(f"  wrote {f.relative_to(Path(__file__).resolve().parents[2])}")


if __name__ == "__main__":
    main()

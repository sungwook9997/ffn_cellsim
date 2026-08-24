"""Regenerate the motor-nmii (I3) analytic-oracle figures — the single entry-point for this track.

Every figure OVERLAYS the closed-form oracle / reference band on the numeric result, annotates SI-ish engine
units (pN, µm/s, 1/s), never truncates an axis, and shows stochastic spread + mean for the ensemble plot
(the professor's visualization-integrity rules, PI 2026-05-21 / 2026-07-07). Pure numpy/matplotlib on the
CPU-green analytic suite — runs on the dev Mac (no Warp/CUDA).

    python aleph/scripts/ac_motor_vis.py    ->  aleph/outputs/ac/motor-nmii/figs/*.png

⚠ Every numeric constant used below is an I0-B3 GAP or provisional value (params_i0b3.yaml), swept here as an
ORACLE variable — NOT a chosen value. The figures show the SHAPE gates (magnitude-independent); the native
gamma-floor magnitude gate is the lead's, and is INVALID until the GAPs close.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aleph.components.motor.bell_kinetics_analytic import (
    bell_f0_from_x_beta,
    bell_off_rate,
    engaged_fraction_steady,
)
from aleph.components.motor.ensemble_stall_analytic import (
    ensemble_stall_meanfield,
    ensemble_stall_stochastic,
)
from aleph.components.motor.hill_fv_analytic import (
    hill_force,
    hill_velocity,
    linear_velocity,
    power_per_head,
)
from aleph.components.motor.minifilament_topology import head_newton_residual, working_stroke_strain

# ── I0-B3 GAP/provisional oracle-sweep values (NOT chosen; see params_i0b3.yaml) ─────────────────
V0 = 0.2          # um/s   GAP: 0.12 vs 0.2
FS_HEAD = 2.0     # pN     GAP: 0.5 vs 2.0
K_OFF0 = 0.35     # 1/s    provisional (Stam-Hocky/Tam)
X_BETA = 0.6e-3   # um     provisional (Veigel 2002) -> f0 = kBT/x_beta
K_ON = 50.0       # 1/s    GAP (config-chosen)
N_SIDE = 28       # heads  GAP: 10 (AFINES) vs 28-30 (Billington)
F0_BELL = bell_f0_from_x_beta(X_BETA)
KB318_BAND = (50.0, 100.0)   # KB-3.18 per-minifilament stall band — REFERENCE overlay, NOT a target
GAP_TAG = "I0-B3 GAP — oracle sweep, not a chosen value"

OUT = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "motor-nmii" / "figs"


def fig_hill_fv() -> Path:
    """Per-head Hill force-velocity for several curvatures + the Hill-1938 point; power sign/single-max."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.4))
    f = np.linspace(0.0, FS_HEAD, 400)
    kappas = [(0.25, "Hill-1938 muscle κ=0.25"), (0.5, "Kovács/archived NMII κ=0.5"), (1.0, "κ=1.0")]
    for kap, lab in kappas:
        ax0.plot(f, hill_velocity(f, V0, FS_HEAD, kap), lw=2, label=lab)
    ax0.plot(f, linear_velocity(f, V0, FS_HEAD), "k--", lw=2, label="linear (κ→∞, PI 2026-07-07)")
    # the distinct Hill-1938 reference point v/v0 = 1/6 at F=Fs/2 for κ=0.25
    ax0.scatter([FS_HEAD / 2], [V0 / 6.0], c="crimson", zorder=6, label="Hill-1938 ref: v/v₀=1/6 @ F=F_s/2")
    ax0.scatter([0.0, FS_HEAD], [V0, 0.0], c="k", marker="o", zorder=6)
    ax0.annotate("force-free\nv(0)=v₀", (0, V0), (0.4, V0 * 0.78),
                 arrowprops=dict(arrowstyle="->", color="gray"))
    ax0.annotate("stall\nv(F_s)=0", (FS_HEAD, 0), (FS_HEAD * 0.62, V0 * 0.18),
                 arrowprops=dict(arrowstyle="->", color="gray"))
    ax0.set_xlabel("resisting load  F  [pN]"); ax0.set_ylabel("shortening velocity  v  [µm/s]")
    ax0.set_title("Single-head Hill FV vs Hill-1938"); ax0.set_xlim(0, FS_HEAD); ax0.set_ylim(0, V0 * 1.05)
    ax0.legend(fontsize=7.5, loc="upper right")

    for kap, _lab in kappas:
        p = power_per_head(f, V0, FS_HEAD, kap)
        ax1.plot(f, p, lw=2, label=f"κ={kap}")
    ax1.axhline(0.0, ls=":", c="gray")
    ax1.scatter([0.0, FS_HEAD], [0.0, 0.0], c="k", zorder=6, label="P=0 at v₀ and at stall")
    ax1.set_xlabel("resisting load  F  [pN]"); ax1.set_ylabel("mechanical power  P=F·v  [pN·µm/s]")
    ax1.set_title("Work sign: P≥0, single interior max"); ax1.set_xlim(0, FS_HEAD); ax1.legend(fontsize=8)
    fig.suptitle(f"I3 Hill force-velocity oracle   (v₀={V0} µm/s, F_s={FS_HEAD} pN — {GAP_TAG})")
    fig.tight_layout()
    pth = OUT / "i3_hill_fv.png"; fig.savefig(pth, dpi=130); plt.close(fig)
    return pth


def fig_bell() -> Path:
    """Bell slip off-rate vs load (monotone slip) + the emergent, self-limiting engaged fraction."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.4))
    f = np.linspace(0.0, 30.0, 400)
    ax0.semilogy(f, bell_off_rate(f, K_OFF0, F0_BELL), lw=2, label="p_off = k_off0·exp(|f|/f₀)")
    ax0.scatter([0.0], [K_OFF0], c="crimson", zorder=6, label=f"p_off(0)=k_off0={K_OFF0}/s")
    ax0.axvline(F0_BELL, ls=":", c="gray")
    ax0.annotate(f"f₀=kBT/x_β≈{F0_BELL:.1f} pN", (F0_BELL, K_OFF0 * 3),
                 (F0_BELL + 2, K_OFF0 * 6), arrowprops=dict(arrowstyle="->", color="gray"))
    ax0.set_xlabel("head load  f  [pN]"); ax0.set_ylabel("off-rate  p_off  [1/s]  (log)")
    ax0.set_title("Bell slip off-rate (monotone in load)"); ax0.set_xlim(0, 30); ax0.legend(fontsize=8)

    phi = engaged_fraction_steady(K_ON, f, K_OFF0, F0_BELL)
    ax1.plot(f, phi, lw=2, label="φ_b = k_on/(k_on+p_off(f))")
    phi0 = K_ON / (K_ON + K_OFF0)
    ax1.scatter([0.0], [phi0], c="crimson", zorder=6, label=f"φ_b(0)={phi0:.3f} (unloaded)")
    ax1.axvline(FS_HEAD, ls="--", c="darkgreen")
    ax1.annotate(f"per-head stall F_s={FS_HEAD} pN\nφ_b still ≈{engaged_fraction_steady(K_ON, FS_HEAD, K_OFF0, F0_BELL):.3f}",
                 (FS_HEAD, engaged_fraction_steady(K_ON, FS_HEAD, K_OFF0, F0_BELL)),
                 (FS_HEAD + 3, 0.6), arrowprops=dict(arrowstyle="->", color="darkgreen"), fontsize=8)
    ax1.set_xlabel("head load  f  [pN]"); ax1.set_ylabel("engaged fraction  φ_b  [–]")
    ax1.set_title("Emergent engaged fraction (self-limiting)"); ax1.set_xlim(0, 30); ax1.set_ylim(0, 1.02)
    ax1.legend(fontsize=8, loc="lower left")
    fig.suptitle(f"I3 Bell slip kinetics oracle   (k_off0={K_OFF0}/s, x_β=0.6 nm, k_on={K_ON}/s — {GAP_TAG})")
    fig.tight_layout()
    pth = OUT / "i3_bell_kinetics.png"; fig.savefig(pth, dpi=130); plt.close(fig)
    return pth


def fig_ensemble_stall() -> Path:
    """Ensemble stall EMERGES from bound heads (stochastic spread + mean) — NOT the imposed N_side·F_head."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.4))
    rng = np.random.default_rng(44)
    mf = ensemble_stall_meanfield(N_SIDE, FS_HEAD, K_ON, K_OFF0, F0_BELL)
    draws = ensemble_stall_stochastic(N_SIDE, FS_HEAD, K_ON, K_OFF0, F0_BELL, n_realizations=40_000, rng=rng)
    ax0.hist(draws, bins=np.arange(draws.min() - FS_HEAD, draws.max() + 2 * FS_HEAD, FS_HEAD),
             density=True, alpha=0.5, color="steelblue", label="per-realisation ensemble (stochastic)")
    ax0.axvline(mf["f_ensemble"], color="crimson", lw=2, label=f"mean-field E[F]={mf['f_ensemble']:.1f} pN")
    ax0.axvline(mf["f_naive"], color="k", ls="--", lw=2, label=f"naive N_side·F_head={mf['f_naive']:.0f} pN (imposed)")
    ax0.axvspan(*KB318_BAND, color="green", alpha=0.10, label="KB-3.18 band 50–100 pN (reference)")
    ax0.set_xlabel("per-minifilament stall force  F_ensemble  [pN]"); ax0.set_ylabel("density  [1/pN]")
    ax0.set_title("Ensemble stall EMERGES (not imposed)"); ax0.legend(fontsize=7.5, loc="upper left")

    koff = np.logspace(-1, 2, 60)
    fe = np.array([ensemble_stall_meanfield(N_SIDE, FS_HEAD, K_ON, k, F0_BELL)["f_ensemble"] for k in koff])
    ax1.semilogx(koff, fe, lw=2, label="emergent E[F] vs zero-force off-rate")
    ax1.axhline(mf["f_naive"], color="k", ls="--", label="naive N_side·F_head (all heads bound)")
    ax1.axhspan(*KB318_BAND, color="green", alpha=0.10, label="KB-3.18 band (reference, NOT a target)")
    ax1.axvline(K_OFF0, ls=":", c="gray"); ax1.annotate(f"k_off0={K_OFF0}", (K_OFF0, fe[0] * 0.5))
    ax1.set_xlabel("zero-force off-rate  k_off0  [1/s]  (log)"); ax1.set_ylabel("emergent stall  E[F]  [pN]")
    ax1.set_title("Self-limiting: more shedding → less force"); ax1.legend(fontsize=7.5, loc="lower left")
    fig.suptitle(f"I3 ensemble-stall emergence   (N_side={N_SIDE}, F_head={FS_HEAD} pN — {GAP_TAG})")
    fig.tight_layout()
    pth = OUT / "i3_ensemble_stall.png"; fig.savefig(pth, dpi=130); plt.close(fig)
    return pth


def fig_newton_and_kxb() -> Path:
    """Per-head Newton closure residual → machine precision + the k_xb master-knob strain arbiter."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.4))
    res = head_newton_residual(x_init=0.5, k_xb=500.0, k_anchor=1000.0, r0_head=0.0, x0_a=0.02, n_iter=4)
    resid = np.array(res["residuals"])
    resid_floor = np.where(resid <= 0, 1e-18, resid)   # show machine-zero on the log axis without a flip
    ax0.semilogy(range(len(resid_floor)), resid_floor, "o-", lw=2, label="|dE/dx| after Newton step")
    ax0.axhline(1e-12, ls=":", c="gray", label="closure tol 1e-12")
    ax0.set_xlabel("Newton iteration  [–]"); ax0.set_ylabel("force residual  |dE/dx|  [pN]  (log)")
    ax0.set_title("Per-head Newton closure (1-step, linear)"); ax0.legend(fontsize=8)

    kxb = np.logspace(-0.5, 3.5, 200)   # 0.3 .. 3000 pN/um
    strain_nm = np.array([working_stroke_strain(FS_HEAD, k) for k in kxb]) * 1e3  # um -> nm
    ax1.loglog(kxb, strain_nm, lw=2, label="strain = F_head/k_xb")
    ax1.axhspan(5.0, 20.0, color="green", alpha=0.12, label="physical working stroke 5–20 nm")
    ax1.axvspan(100.0, 1000.0, color="green", alpha=0.08, label="physical k_xb 100–1000 pN/µm")
    ax1.axhline(301.0, color="purple", ls="--", label="minifilament length ≈301 nm")
    k_broken = 1.0
    ax1.scatter([k_broken], [working_stroke_strain(FS_HEAD, k_broken) * 1e3], c="crimson", zorder=6,
                label="archived 1 pN/µm → 2000 nm (BREAKS stall)")
    ax1.set_xlabel("crossbridge stiffness  k_xb  [pN/µm]  (log)"); ax1.set_ylabel("working-stroke strain  [nm]  (log)")
    ax1.set_title("k_xb MASTER knob: '1 pN/µm breaks stall'"); ax1.legend(fontsize=7, loc="upper right")
    fig.suptitle(f"I3 per-head Newton closure + k_xb strain arbiter   (F_head={FS_HEAD} pN — {GAP_TAG})")
    fig.tight_layout()
    pth = OUT / "i3_newton_kxb.png"; fig.savefig(pth, dpi=130); plt.close(fig)
    return pth


def fig_weave_joint() -> Path:
    """I4-weave: the actin↔motor JOINT — per-segment barbed directions + polarity-correctness on a woven cortex.

    LEFT: a small woven cortex patch (weave_cell([CORTEX])) with each actin segment's barbed-end unit direction
    (the walk_dir a bound head is overwritten with) quivered at the segment midpoint. RIGHT: the polarity gate —
    every segment's barbed dir projects POSITIVELY on the node→barbed-end direction (correct sign), and a
    straight fiber reproduces the node-anchored walk_dir exactly (the two polarity paths agree). Pure host
    geometry (no magnitude); the on-device overwrite is the native gate.
    """
    import dataclasses

    from aleph.laws.architecture_spec import CORTEX
    from aleph.components.weave.regions import RegionSpec
    from aleph.components.weave.walk_dir import barbed_end_node, walk_dir_from_polarity
    from aleph.components.weave.woven_cell import weave_cell

    fil = dataclasses.replace(CORTEX.filament, n_filaments=90)
    wc = weave_cell([RegionSpec(arch=dataclasses.replace(CORTEX, filament=fil))], rng=np.random.default_rng(0))
    seg_a, seg_b, seg_pol = wc.actin_segment_topology()
    seg_dir = wc.segment_barbed_directions()
    mid = 0.5 * (wc.pos[seg_a] + wc.pos[seg_b])

    fig = plt.figure(figsize=(12, 5.0))
    ax0 = fig.add_subplot(1, 2, 1, projection="3d")
    # thin the quiver for legibility (draw a representative sample; report the total in the title)
    step = max(seg_a.shape[0] // 400, 1)
    sl = slice(None, None, step)
    ax0.quiver(mid[sl, 0], mid[sl, 1], mid[sl, 2],
               seg_dir[sl, 0], seg_dir[sl, 1], seg_dir[sl, 2],
               length=0.6, normalize=True, color="crimson", linewidth=0.6, arrow_length_ratio=0.4)
    ax0.scatter(wc.pos[:, 0], wc.pos[:, 1], wc.pos[:, 2], s=1.5, c="steelblue", alpha=0.35)
    ax0.set_title(f"Woven cortex: per-segment barbed dir (walk_dir)\n{seg_a.shape[0]} segments "
                  f"(every {step}th drawn)", fontsize=9)
    ax0.set_xlabel("x [µm]"); ax0.set_ylabel("y [µm]"); ax0.set_zlabel("z [µm]")

    # polarity-correctness: projection of segment barbed dir on the node->barbed-end direction (must be > 0)
    barbed = barbed_end_node(wc.fiber_offsets, wc.polarity)
    fib = wc.node_fiber[seg_a]
    proj = []
    for s in range(seg_a.shape[0]):
        nd = walk_dir_from_polarity(int(seg_a[s]), int(barbed[fib[s]]), wc.pos)
        if np.linalg.norm(nd) < 1e-9:
            continue
        proj.append(float(np.dot(seg_dir[s], nd)))
    proj = np.array(proj)
    ax1 = fig.add_subplot(1, 2, 2)
    ax1.hist(proj, bins=40, color="seagreen", alpha=0.8)
    ax1.axvline(0.0, color="crimson", ls="--", lw=2, label="sign boundary (all must be > 0)")
    ax1.set_xlim(-1.05, 1.05)
    ax1.set_xlabel("segment barbed dir · (node→barbed-end) direction  [–]")
    ax1.set_ylabel("segment count")
    ax1.set_title(f"Polarity gate: {int(np.sum(proj > 0))}/{proj.shape[0]} segments point toward the barbed end\n"
                  f"(min projection = {proj.min():.3f} > 0)", fontsize=9)
    ax1.legend(fontsize=8, loc="upper left")
    fig.suptitle("I4-weave actin↔motor joint — barbed polarity is correct across the woven network "
                 "(geometry, not a magnitude)")
    fig.tight_layout()
    pth = OUT / "i4_weave_joint.png"; fig.savefig(pth, dpi=130); plt.close(fig)
    return pth


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    figs = [fig_hill_fv(), fig_bell(), fig_ensemble_stall(), fig_newton_and_kxb(), fig_weave_joint()]
    root = Path(__file__).resolve().parents[2]
    for f in figs:
        print(f"  wrote {f.relative_to(root)}")


if __name__ == "__main__":
    main()

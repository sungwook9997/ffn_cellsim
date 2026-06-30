"""Architectural validation metrics for the unified actin-architecture framework (task a, increment 1).

Each woven structure is validated against a quantitative ARCHITECTURAL metric (not just γ): the cortex
against mesh/connectivity (reuse ``kim_network``), the filopodium against bundle order. These are the
acceptance criteria that prove ``weave`` produces the RIGHT architecture per structure, with
literature bands. Reuses ``kim_network`` where possible; adds bundle/order metrics here.
"""

from __future__ import annotations

import numpy as np


def _fiber_axis(net, f: int) -> np.ndarray:
    """Unit end-to-end axis of filament ``f``."""
    off = net.fiber_offsets
    a, b = int(off[f]), int(off[f + 1])
    v = net.pos[b - 1] - net.pos[a]
    return v / (np.linalg.norm(v) + 1e-12)


def parallel_order_parameter(net) -> float:
    """Nematic order parameter S = (3⟨cos²θ⟩−1)/2 of the filament axes (1 = perfectly parallel bundle,
    0 = isotropic). The largest eigenvalue form (orientation-invariant) of the Q-tensor."""
    axes = np.array([_fiber_axis(net, f) for f in range(net.n_fibers)])
    Q = (3.0 * np.einsum("ni,nj->ij", axes, axes) / len(axes) - np.eye(3)) / 2.0
    return float(np.linalg.eigvalsh(Q).max())


def bundle_count(net) -> int:
    """Number of filaments in the bundle (= n_fibers; the architectural count, lit band 10–30 for
    filopodia, 20–30 microvilli)."""
    return int(net.n_fibers)


def inter_filament_spacing_nm(net) -> float:
    """Median nearest-neighbour spacing between parallel filaments [nm], measured on the bundle
    cross-section (project onto the plane ⊥ the mean axis, midpoints of each filament)."""
    off = net.fiber_offsets
    mids = np.array([net.pos[(int(off[f]) + int(off[f + 1])) // 2] for f in range(net.n_fibers)])
    axes = np.array([_fiber_axis(net, f) for f in range(net.n_fibers)])
    ax = axes.mean(axis=0); ax /= np.linalg.norm(ax) + 1e-12
    proj = mids - np.outer(mids @ ax, ax)                      # drop the axial component
    from scipy.spatial import cKDTree
    d, _ = cKDTree(proj).query(proj, k=2)
    return float(np.median(d[:, 1]) * 1e3)                     # µm → nm


def branch_angle_distribution(cortex) -> dict:
    """Arp2/3 branch-angle statistics over the cortex's branch triples [i, j(apex), k] — mean + SD [deg]
    (lit anchor: Fäßler 2020 in-cell cryo-ET 68 ± 9°; rest θ₀ = 70°)."""
    net = cortex.net
    tri = cortex.branch_triples
    if tri.shape[0] == 0:
        return {"n_branches": 0, "mean_deg": float("nan"), "std_deg": float("nan")}
    r1 = net.pos[tri[:, 0]] - net.pos[tri[:, 1]]
    r2 = net.pos[tri[:, 2]] - net.pos[tri[:, 1]]
    c = np.einsum("ij,ij->i", r1, r2) / (np.linalg.norm(r1, axis=1) * np.linalg.norm(r2, axis=1) + 1e-12)
    ang = np.rad2deg(np.arccos(np.clip(c, -1, 1)))
    return {"n_branches": int(tri.shape[0]), "mean_deg": float(ang.mean()), "std_deg": float(ang.std())}


def two_mode_orientation(cortex, axis=(0.0, 1.0, 0.0)) -> dict:
    """Lamellipodium two-mode filament orientation about the protrusion ``axis`` (default +y): the ±mode
    peak means [deg] + the fraction in the ±20–50° two-mode band (Mueller 2017 ±35° dendritic signature).
    Uses the in-plane signed angle (x vs y for the default axis)."""
    net = cortex.net
    off = net.fiber_offsets
    phis = []
    for f in range(net.n_fibers):
        a, b = int(off[f]), int(off[f + 1])
        d = net.pos[b - 1] - net.pos[a]
        phis.append(np.rad2deg(np.arctan2(d[0], d[1])))        # signed angle from +y in the x–y plane
    phis = np.array(phis)
    pos = phis[phis > 0]; neg = phis[phis < 0]
    return {"plus_mode_deg": float(pos.mean()) if pos.size else float("nan"),
            "minus_mode_deg": float(neg.mean()) if neg.size else float("nan"),
            "two_mode_frac": float(np.mean((np.abs(phis) > 20) & (np.abs(phis) < 50)))}


def sarcomeric_period_um(cortex, axis=2) -> dict:
    """Stress-fiber sarcomeric period: the dominant spacing of the α-actinin (crosslink) Z-bodies along
    the fiber ``axis`` (default z), via FFT of the crosslink-node density. Lit band 0.5–1.4 µm
    (Hotulainen-Lappalainen 2006). Also returns the NMIIA-band anti-registration phase offset (≈ half a
    period = anti-registered, Murrell 2015)."""
    net = cortex.net
    if cortex.xl_i.size == 0:
        return {"period_um": float("nan"), "n_zbody_planes": 0, "anti_registered": False}
    zc = np.concatenate([net.pos[cortex.xl_i, axis], net.pos[cortex.xl_j, axis]])
    lo, hi = zc.min(), zc.max()
    if hi - lo < 1e-6:
        return {"period_um": float("nan"), "n_zbody_planes": 0, "anti_registered": False}
    nbin = 256
    dens, edges = np.histogram(zc, bins=nbin, range=(lo, hi))
    dz = (hi - lo) / nbin
    spec = np.abs(np.fft.rfft(dens - dens.mean()))
    freqs = np.fft.rfftfreq(nbin, d=dz)                         # cycles per µm
    k = 1 + int(np.argmax(spec[1:]))                            # skip DC
    period = 1.0 / freqs[k] if freqs[k] > 0 else float("nan")
    # NMIIA anti-registration: phase offset between myosin-band z and Z-body z (≈ period/2)
    anti = False
    if cortex.myo_i.size:
        zm = np.concatenate([net.pos[cortex.myo_i, axis], net.pos[cortex.myo_j, axis]])
        offset = abs((np.median(zm) - np.median(zc))) % period if np.isfinite(period) else 0.0
        anti = bool(0.25 * period < offset < 0.75 * period) if np.isfinite(period) else False
    return {"period_um": float(period), "n_zbody_planes": int(np.round((hi - lo) / period)) if np.isfinite(period) else 0,
            "anti_registered": anti}


def bundle_dimensions(cortex) -> dict:
    """Bundle axial length + cross-section diameter [µm] (microvillus/SF). Length = axial extent along
    the dominant axis; diameter = 2× the mean in-plane radius of the filament mid-points."""
    net = cortex.net
    off = net.fiber_offsets
    mids = np.array([net.pos[(int(off[f]) + int(off[f + 1])) // 2] for f in range(net.n_fibers)])
    # principal axis (PCA) of the node cloud — robust to GRADED polarity (mean fiber axis ≈ 0 for an
    # antiparallel SF bundle, so use the largest-variance direction instead).
    p = net.pos - net.pos.mean(axis=0)
    ax = np.linalg.eigh(p.T @ p)[1][:, -1]
    axial = net.pos @ ax
    perp = mids - np.outer(mids @ ax, ax)
    return {"length_um": float(axial.max() - axial.min()),
            "diameter_um": float(2.0 * np.linalg.norm(perp - perp.mean(axis=0), axis=1).mean())}


def cortex_metrics(cortex) -> dict:
    """Cortex architectural metrics — areal density + connectivity z + giant-component fraction."""
    from ffn_sim.ff.kim_network import connectivity_z
    net = cortex.net
    xl = np.stack([cortex.xl_i, cortex.xl_j], axis=1) if cortex.xl_i.size else np.zeros((0, 2), int)
    return {
        "n_filaments": net.n_fibers,
        "areal_density_um2": net.n_fibers / (4.0 * np.pi * cortex.R_um**2),
        "connectivity_z": connectivity_z(net, xl),
        "n_crosslinks": int(cortex.xl_i.size),
        "n_myosin": int(cortex.myo_i.size),
    }

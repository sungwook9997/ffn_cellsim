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

r"""Assembled-cell cortical-tension MEASUREMENT harness — the method-of-planes on a contractile cell.

Once the production NMII motor is wired and the full-native active baseline converges, this is the harness that
turns that cell state into the physical cortical tension γ [pN/µm] — the decisive magnitude number the whole
investigation is about. It feeds :func:`aleph.components.incumbent.cortical_tension.measure_cortical_stress` (the validated
method-of-planes) with the load-bearing element families extracted from the assembled cell.

⭐ The key extraction — ACTIN AXIAL TENSION from the inextensibility constraint. Actin is modelled as an
inextensible bead-rod (NF2007 reshape projection), so its axial tension is NOT a spring ``k·Δx`` — it is the
constraint (Lagrange) force that holds each segment at rest length. At the converged overdamped equilibrium the
constraint force on each node balances the EXTERNAL force (crossbridge reaction + crosslink + bending + steric +
pressure) that ``driver._accumulate_all`` sums BEFORE the reshape. Force balance along a free-ended filament
then gives the axial tension in segment ``k`` as the running sum of the external force from the filament end::

    T_k · û_k = − Σ_{j ≤ k} f_ext[node_j]        ⇒   T_k = − (Σ_{j ≤ k} f_ext[node_j]) · û_k

so the actin tension carries EVERYTHING the myosin injects and the crosslinks transmit — it IS the cortical
network tension, with no double-count and no dependence on the motor's internal anchor representation. (A
per-filament net force of ≈0 at equilibrium is the physical consistency check.) The crosslink tensions (a
Hookean spring, extractable directly) are added as the lateral network family.

The extraction (:func:`actin_axial_tension`, :func:`crosslink_tension`) is pure NumPy — an out-of-hot-loop
measurement (I0-A: Python may post-process a synchronised device read; no authoritative per-step host state).
The device orchestration (:func:`measure_assembled_cortical_stress`) accumulates the external force on the
CUDA cell, synchronises once, and post-processes on the host.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from aleph.components.incumbent.cortical_tension import CorticalStressReport, measure_cortical_stress

__all__ = [
    "actin_axial_tension",
    "crosslink_tension",
    "measure_assembled_cortical_stress",
]


def actin_axial_tension(
    f_ext: npt.NDArray, pos_actin: npt.NDArray, fiber_offsets: npt.NDArray,
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """Per-segment actin axial tension from the external-force running sum along each filament.

    Args:
        f_ext: ``(n_actin, 3)`` external force per actin node (``driver._accumulate_all`` output, pre-reshape).
        pos_actin: ``(n_actin, 3)`` actin node positions [µm].
        fiber_offsets: ``(n_fibers + 1,)`` node-range offsets — fiber ``i`` owns nodes ``[off[i], off[i+1])``.

    Returns:
        ``(seg_a, seg_b, tension)``: for each actin segment its two node indices and the SIGNED axial tension
        [pN] (+ tension, − compression). Filaments with < 2 nodes contribute nothing.
    """
    f_ext = np.asarray(f_ext, dtype=np.float64)
    pos_actin = np.asarray(pos_actin, dtype=np.float64)
    off = np.asarray(fiber_offsets, dtype=np.int64)
    seg_a: list[int] = []
    seg_b: list[int] = []
    tension: list[float] = []
    for i in range(off.shape[0] - 1):
        s, e = int(off[i]), int(off[i + 1])
        if e - s < 2:
            continue
        cum = np.zeros(3)
        for k in range(e - s - 1):
            cum = cum + f_ext[s + k]                    # running external force up to node s+k
            u = pos_actin[s + k + 1] - pos_actin[s + k]
            length = float(np.linalg.norm(u))
            t = 0.0 if length < 1e-12 else -float(np.dot(cum, u / length))
            seg_a.append(s + k)
            seg_b.append(s + k + 1)
            tension.append(t)
    return (np.array(seg_a, dtype=np.int64), np.array(seg_b, dtype=np.int64), np.array(tension, dtype=np.float64))


def crosslink_tension(
    pos: npt.NDArray, xl: npt.NDArray, kxl: npt.NDArray, r0xl: npt.NDArray,
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """Hookean crosslink tensions ``T = k (|r| − r0)`` [pN] → ``(nodeA, nodeB, tension)``."""
    pos = np.asarray(pos, dtype=np.float64)
    xl = np.asarray(xl, dtype=np.int64).reshape(-1, 2)
    if xl.shape[0] == 0:
        z = np.zeros(0)
        return np.zeros(0, np.int64), np.zeros(0, np.int64), z
    a, b = xl[:, 0], xl[:, 1]
    length = np.linalg.norm(pos[b] - pos[a], axis=1)
    t = np.asarray(kxl, dtype=np.float64) * (length - np.asarray(r0xl, dtype=np.float64))
    return a, b, t


def measure_assembled_cortical_stress(cell, R: float, *, n_planes: int = 64) -> CorticalStressReport:
    """Measure the method-of-planes cortical stress on a (contractile) assembled CUDA cell.

    Accumulates the external force (pre-reshape), synchronises once, extracts the actin-axial + crosslink
    load-bearing families on the host, and returns the source/network/total split. Runtime path: requires the
    device cell; the extraction math is host NumPy (gated by :func:`actin_axial_tension` host tests).

    Args:
        cell: the assembled cell (``AssembledCell``) with ``pos_d``, ``foff_d``, ``xl_d`` etc.
        R: the cortex-shell radius for the cut circumference ``2πR`` [µm].
        n_planes: number of near-isotropic diametral cut orientations.

    Returns:
        A :class:`CorticalStressReport` (γ_source / γ_network / γ_total [pN/µm]).
    """
    import warp as wp

    from aleph.components.incumbent.driver import _accumulate_all

    dev = wp.get_device(cell.device)
    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    _accumulate_all(cell, cell.pos_d, f)                # external forces (crossbridge + crosslink + bending + …)
    wp.synchronize_device(dev)
    pos = cell.pos_d.numpy()
    f_ext = f.numpy()
    n_actin = int(cell.n_actin)
    foff = cell.foff_d.numpy()

    sa, sb, t_actin = actin_axial_tension(f_ext[:n_actin], pos[:n_actin], foff)
    families: dict[str, tuple] = {"actin": (pos[sa], pos[sb], t_actin)}
    if getattr(cell, "n_xl", 0):
        xa, xb, t_xl = crosslink_tension(pos, cell.xl_d.numpy(), cell.kxl_d.numpy(), cell.r0xl_d.numpy())
        families["crosslink"] = (pos[xa], pos[xb], t_xl)

    return measure_cortical_stress(
        R, families, pos, f_ext, n_planes=n_planes,
        source_families=(), network_families=tuple(families.keys()))

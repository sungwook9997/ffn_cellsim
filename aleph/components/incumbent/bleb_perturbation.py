r"""Local membrane–cortex de-adhesion — the Oracle-B blebbing perturbation (Warp-CUDA).

Charras et al. (2005/2008) nucleate a bleb by locally detaching the plasma membrane from the actin cortex:
inside the ablated patch the ERM linkers no longer restrain the bilayer, so the resting pore pressure drives
it outward and a bleb grows *emergently* from the existing membrane (bending + tension) and pressure physics.
This module supplies ONLY that perturbation — an angular-cap de-adhesion of the ERM tethers — never the bleb
itself; the bulge, its expansion, and the pressure transient are outputs of the assembled cell, not scripted.

The de-adhesion is applied to the SAME ``erm_bound_d`` state the commit-safe rupture kernel owns, and with the
same discipline: it is IRREVERSIBLE biology, so it may fire only once at an ACCEPTED outer physical-time
boundary (``accepted[0]==1``) — never inside the quasi-static inner iteration, where a rejected mechanical
candidate must not leave a bleb behind. A tether is de-adhered iff its membrane endpoint lies within the cap
``dot((x_mem − centroid)/|·|, patch_axis) ≥ cos(half_angle)``. Unbinding is monotone (``1→0`` only); it never
re-binds a tether, so replaying the kernel is idempotent.

Sanity Gate:
    * Geometry: a tether exactly on the cap axis (cos=1) always de-adheres for any half-angle < 90°; one on the
      opposite pole (cos=−1) never does; the boundary ``cos == cos(half_angle)`` de-adheres (``≥``).
    * Time: the kernel writes ``bound`` only when ``accepted[0]==1``; solver iterations are not biological time.
    * Monotonicity: ``bound`` only ever goes 1→0; an already-unbound tether is untouched (idempotent replay).
    * Conservation: de-adhesion removes tethers from the force sum only; it injects no force and moves no node.
    * Residency: positions, indices, bound state and the accept flag are CUDA arrays; no host readback.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import warp as wp

__all__ = [
    "erm_patch_deadhesion_kernel",
    "erm_patch_deadhesion_reference",
    "cap_cos_half_angle",
    "count_tethers_in_cap",
]


@wp.kernel
def erm_patch_deadhesion_kernel(
    pos: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    centroid: wp.vec3d,
    patch_axis: wp.vec3d,
    cos_half_angle: wp.float64,
    accepted: wp.array(dtype=wp.int32),
) -> None:
    """Irreversibly de-adhere the ERM tethers whose membrane endpoint lies in the ablation cap.

    Commit-only: like ``erm_tether_rupture_kernel`` this may run once per ACCEPTED outer step. ``patch_axis``
    must be a unit vector; ``cos_half_angle`` is the cosine of the cap half-angle (larger ⇒ smaller patch).
    """
    k = wp.tid()
    if accepted[0] == 0 or bound[k] == 0:
        return
    r = pos[membrane_idx[k]] - centroid
    d = wp.length(r)
    if d < wp.float64(1.0e-12):
        return
    if wp.dot(r / d, patch_axis) >= cos_half_angle:
        bound[k] = 0


def erm_patch_deadhesion_reference(
    membrane_pos: npt.NDArray[np.float64],
    bound: npt.NDArray[np.integer],
    centroid: npt.NDArray[np.float64],
    patch_axis: npt.NDArray[np.float64],
    cos_half_angle: float,
) -> npt.NDArray[np.int32]:
    """Pure-NumPy acceptance oracle: the ``bound`` state after one accepted de-adhesion.

    ``membrane_pos`` is the membrane endpoint of each tether (``pos[membrane_idx]``), shape (Ne, 3).
    """
    x = np.asarray(membrane_pos, dtype=np.float64)
    b = np.asarray(bound, dtype=np.int32).copy()
    axis = np.asarray(patch_axis, dtype=np.float64)
    r = x - np.asarray(centroid, dtype=np.float64)
    d = np.linalg.norm(r, axis=1)
    safe = d > 1.0e-12
    cos = np.full(d.shape, -2.0)
    cos[safe] = (r[safe] @ axis) / d[safe]
    unbind = safe & (b == 1) & (cos >= float(cos_half_angle))
    b[unbind] = 0
    return b


def cap_cos_half_angle(half_angle_deg: float) -> float:
    """Cosine of a cap half-angle in degrees (the kernel's ``cos_half_angle`` input)."""
    if not 0.0 < half_angle_deg < 90.0:
        raise ValueError("half_angle_deg must be in (0, 90)")
    return float(np.cos(np.deg2rad(half_angle_deg)))


def count_tethers_in_cap(
    membrane_pos: npt.NDArray[np.float64],
    bound: npt.NDArray[np.integer],
    centroid: npt.NDArray[np.float64],
    patch_axis: npt.NDArray[np.float64],
    cos_half_angle: float,
) -> int:
    """How many currently-bound tethers the cap would de-adhere — for sizing the ablation patch."""
    before = np.asarray(bound, dtype=np.int32)
    after = erm_patch_deadhesion_reference(membrane_pos, before, centroid, patch_axis, cos_half_angle)
    return int(np.count_nonzero((before == 1) & (after == 0)))

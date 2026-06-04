"""Layer-2 D-axis: explicit SUBSTRATE-reacted collective crawl (the magnitude-gap test).

The §C decisive bracket showed that the lamellipodial-protrusion-anchored traction (9.4 nN),
applied as the existing EDGE-localized, 3D-radial, cell-cell-transmitted force
(``spreading.edge_outward_forces``), EJECTS rim cells: a single rim cell pulls away from the
stationary bulk and its lone catch contact (≤6.5 nN) ruptures. That is the artifact of *how*
the force was applied, not of the magnitude.

This module applies the SAME anchored force the FAITHFUL way — the way a real lamellipodium
spreads a colony:

  * **Substrate-reacted, not cohesion-competing.** A basal cell's lamellipodium grips the dish
    (integrin clutch) and the protrusion force is reacted by the SUBSTRATE (Newton's 3rd), not
    by its neighbours. In the overdamped CBM the substrate coupling is the clutch drag γ_cell
    (KU-2.18), so the crawl velocity is f_crawl/γ; cohesion is a SEPARATE parallel constraint
    that holds the sheet, it does not oppose the propulsion.
  * **In-plane (xy), along the dish.** The crawl is tangent to the z=0 substrate plane (the
    cell migrates ACROSS the dish), so the force has NO z-component — unlike the 3D-radial
    edge force, which also pushed cap cells up/out.
  * **Collective (every basal cell), not rim-only.** All cells touching the substrate
    polarise outward (toward free dish / away from the in-plane centroid) and crawl together,
    so the sheet expands COHERENTLY — neighbouring cells move at nearly the same velocity, the
    inter-cell catch bonds stretch only at the (small) velocity *difference*, and no single
    bond bears the full protrusion load. This is the load-distribution the rim-only forcing
    lacked.

So this is the literature's "collective active wetting / basal crawling" (Pérez-González 2019;
the lamellipodium-as-2D-substrate-motility role, Mattila 2008) rendered mechanistically on the
center-based CBM — a new D-axis force, NOT a tuned parameter (the magnitude is the anchored
``motility_bridge`` protrusion value; the direction and basal selection are geometry).

Sanity Gate
-----------
- Dimensional: positions, z_substrate, basal_band [m]; f_crawl [N]; returns force [N].
- Boundary: f_crawl=0 ⇒ zero force (reduces to the substrate-confined catch CBM). A cell far
  above the dish (z ≫ z_substrate+basal_band) feels ~0 (only basal cells crawl).
- Sign/sense: force points OUTWARD from the in-plane centroid (spreading), z-component ≡ 0.
- Conservation: the net in-plane force is an outward "expansion" (no single direction biased);
  COM stays put in xy to within the basal asymmetry; z-momentum untouched (substrate handles z).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["substrate_crawl_forces"]


def substrate_crawl_forces(
    positions: npt.NDArray[np.float64],
    *,
    f_crawl: float,
    z_substrate: float,
    basal_band: float,
) -> npt.NDArray[np.float64]:
    """In-plane, substrate-reacted, collective outward crawl force per cell (N).

    Every cell within ``basal_band`` of the substrate plane (``z ≤ z_substrate + basal_band``)
    gets an outward IN-PLANE (xy) force of magnitude ``f_crawl`` directed away from the basal
    in-plane centroid (the spreading direction); the z-component is zero (the cell migrates
    along the dish, the substrate wall reacts z). Non-basal (apical) cells feel no crawl. A
    smooth half-cosine taper over the band edge avoids a discontinuity as a cell lifts off.

    Args:
        positions: (N, 3) cell centers (m).
        f_crawl: per-basal-cell outward crawl magnitude (N) — the anchored protrusion force.
        z_substrate: substrate plane height (m); cells rest near here (= substrate r0_sub).
        basal_band: thickness above the plane within which a cell counts as basal (m), ~r0.

    Returns:
        (N, 3) crawl force (N), z-component identically 0.
    """
    p = np.asarray(positions, dtype=np.float64)
    n = p.shape[0]
    out = np.zeros((n, 3), dtype=np.float64)
    if f_crawl == 0.0 or n == 0:
        return out
    if basal_band <= 0.0:
        raise ValueError("basal_band must be strictly positive [m].")

    # Basal membership: full weight at the substrate, smooth half-cosine taper to 0 over the
    # band, 0 above. A geometry weight (which cells touch the dish), not a physics constant.
    dz = p[:, 2] - z_substrate
    w = np.zeros(n, dtype=np.float64)
    in_band = (dz <= basal_band) & (dz >= -basal_band)
    w[in_band] = 0.5 * (1.0 + np.cos(np.pi * np.clip(np.abs(dz[in_band]) / basal_band, 0.0, 1.0)))

    # In-plane outward direction from the BASAL in-plane centroid (the spreading front).
    basal = w > 0.0
    if not basal.any():
        return out
    com_xy = p[basal, :2].mean(axis=0)
    d_xy = p[:, :2] - com_xy
    r_xy = np.linalg.norm(d_xy, axis=1)
    rhat = np.zeros((n, 2), dtype=np.float64)
    nz = r_xy > 0.0
    rhat[nz] = d_xy[nz] / r_xy[nz, None]

    out[:, 0] = f_crawl * w * rhat[:, 0]
    out[:, 1] = f_crawl * w * rhat[:, 1]
    # out[:, 2] stays 0 — crawl is tangent to the dish.
    return out

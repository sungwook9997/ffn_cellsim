"""Young–Dupré cell-contact-angle oracle (doublet + triplet).

ACCEPTANCE ORACLE — closed-form only, runtime-import-forbidden (CLAUDE.md). Used to
validate that the DCM engine's emergent cell–cell contact mechanics obey the capillary
force balance at a three-phase contact line, i.e. that an aggregate is *mechanically*
correctly assembled (not just geometrically touching). This is the gate SimuCell3D uses
to validate its own contact model (Runser, Vetter & Iber, Nat. Comput. Sci. 4, 299–309,
2024, Supplementary Fig. 5); the analytic curves there are from the Lagrangian doublet
framework of **Maître et al., Nature 536, 344–348 (2016)** — the literature anchor here.

Relations (transcribed exactly from SimuCell3D Supp Fig 5, no Young's-law factor of 2):

  Doublet  cos θ = (γ_l,1 + γ_l,2) / (2·γ_a)          [symmetric: cos θ = γ_l / γ_a]
  Triplet  cos(φ/2) = η/2,  η = (γ_l,2 + γ_l,3)/(γ_l,1 + γ_l,2)

  γ_a   = apical / free surface tension (cell–medium interface) [N/m]
  γ_l,i = lateral surface tension contributed by cell i's membrane to the shared
          cell–cell interface [N/m]; adhesion LOWERS γ_l relative to γ_a.
  θ     = per-cell half-angle at the contact line, between the flat cell–cell interface
          and the curved free (medium-facing) surface. θ→0° as γ_l/γ_a→1 (no adhesion,
          tangential touch); θ→90° as γ_l/γ_a→0 (strong adhesion, full flattening).

DCM mapping. Our free faces carry cortical tension γ_surf; the cell–cell interface
membrane tension is reduced by the cadherin adhesion energy density w_adh [J/m² = N/m]:
γ_l = γ_surf − w_adh. With a symmetric doublet (γ_l,1 = γ_l,2):

  cos θ = γ_l / γ_a = (γ_surf − w_adh) / γ_surf = 1 − w_adh/γ_surf

so w_adh ∈ [0, γ_surf] maps θ ∈ [0°, 90°]; w_adh > γ_surf is unphysical (the two
membranes would gain energy by unbounded interface growth — the spreading instability).

Sanity Gate
- Dimensional: γ_a, γ_l, w_adh all [N/m]; θ dimensionless (rad). cos θ a pure ratio. OK.
- Boundary: w_adh=0 → cos θ=1 → θ=0 (touching spheres, no flattening). w_adh=γ_surf →
  cos θ=0 → θ=90° (hemispheres, fully wetted). Both verified in unit tests below.
- Sign-sense: increasing adhesion (w_adh↑) lowers cos θ, raises θ (more flattening). OK.
- Measurement-protocol: θ extracted from a relaxed doublet by fitting each cell's FREE
  surface to a sphere (centre c, radius r) and the interface plane (⊥ to the line of
  centres at the contact); then cos θ = d/r with d = |c − plane|. Limits: d=r → θ=0
  (sphere tangent to plane, a_contact=0); d=0 → θ=90° (centre on the plane, hemisphere).
"""

from __future__ import annotations

import numpy as np


def doublet_angle(gamma_a: float, gamma_l1: float, gamma_l2: float | None = None) -> float:
    """Young–Dupré doublet contact angle θ [rad].

    Args:
        gamma_a: apical / free-surface tension γ_a [N/m].
        gamma_l1: cell-1 lateral (interface) membrane tension γ_l,1 [N/m].
        gamma_l2: cell-2 lateral tension γ_l,2 [N/m]; defaults to γ_l,1 (symmetric).

    Returns:
        Contact half-angle θ [rad] in [0, π/2]. Raises ValueError if the tension ratio
        falls outside [0, 1] (cos θ must be a physical cosine).
    """
    gl2 = gamma_l1 if gamma_l2 is None else gamma_l2
    if gamma_a <= 0:
        raise ValueError(f"gamma_a must be > 0, got {gamma_a}")
    cos_t = (gamma_l1 + gl2) / (2.0 * gamma_a)
    if not (0.0 <= cos_t <= 1.0):
        raise ValueError(f"cos θ = {cos_t:.3f} ∉ [0,1] — unphysical tension ratio "
                         f"(γ̄_l={0.5 * (gamma_l1 + gl2):.2e}, γ_a={gamma_a:.2e})")
    return float(np.arccos(cos_t))


def doublet_angle_from_adhesion(gamma_surf: float, w_adh: float) -> float:
    """DCM-mapped doublet angle θ [rad] from cortical tension + adhesion energy density.

    cos θ = 1 − w_adh/γ_surf  (symmetric doublet, γ_l = γ_surf − w_adh).

    Args:
        gamma_surf: cortical / free-surface tension γ_surf [N/m].
        w_adh: effective cell–cell adhesion energy density [J/m² = N/m], in [0, γ_surf].

    Returns:
        Contact half-angle θ [rad].
    """
    gamma_l = gamma_surf - w_adh
    return doublet_angle(gamma_surf, gamma_l, gamma_l)


def triplet_angle(gamma_l1: float, gamma_l2: float, gamma_l3: float) -> float:
    """Young–Dupré tricellular-junction angle φ [rad]: cos(φ/2) = η/2,
    η = (γ_l,2 + γ_l,3)/(γ_l,1 + γ_l,2)."""
    eta = (gamma_l2 + gamma_l3) / (gamma_l1 + gamma_l2)
    c = eta / 2.0
    if not (0.0 <= c <= 1.0):
        raise ValueError(f"cos(φ/2) = {c:.3f} ∉ [0,1] — unphysical η={eta:.3f}")
    return float(2.0 * np.arccos(c))


def angle_from_doublet_geometry(centre: np.ndarray, radius: float,
                                plane_point: np.ndarray, plane_normal: np.ndarray) -> float:
    """Contact angle θ [rad] measured from a relaxed doublet's geometry.

    cos θ = d/r, d = signed distance from the free-surface sphere centre to the interface
    plane, r = free-surface sphere radius. (d=r → θ=0 tangent; d=0 → θ=90° hemisphere.)

    Args:
        centre: (3,) fitted centre of the cell's FREE-surface sphere [m].
        radius: fitted free-surface sphere radius r [m].
        plane_point: (3,) any point on the cell–cell interface plane [m].
        plane_normal: (3,) unit normal of the interface plane.
    """
    n = np.asarray(plane_normal, float)
    n = n / max(np.linalg.norm(n), 1e-30)
    d = abs(float(np.dot(np.asarray(centre, float) - np.asarray(plane_point, float), n)))
    cos_t = min(1.0, d / max(radius, 1e-30))
    return float(np.arccos(cos_t))


def fit_sphere(points: np.ndarray) -> tuple[np.ndarray, float]:
    """Least-squares sphere fit (algebraic, Coope 1993) to (N,3) points → (centre, radius)."""
    P = np.asarray(points, float)
    A = np.hstack([2.0 * P, np.ones((P.shape[0], 1))])
    b = (P ** 2).sum(axis=1)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    centre = sol[:3]
    radius = float(np.sqrt(max(sol[3] + centre @ centre, 0.0)))
    return centre, radius


# --- Sanity-gate unit checks (run as a module: python -m ...young_dupre) -------------
def _sanity() -> None:
    import math
    # boundary cases
    assert abs(doublet_angle_from_adhesion(1e-4, 0.0) - 0.0) < 1e-9, "no adhesion → θ=0"
    assert abs(doublet_angle_from_adhesion(1e-4, 1e-4) - math.pi / 2) < 1e-9, "full → θ=90°"
    # monotone: more adhesion → larger angle
    ts = [doublet_angle_from_adhesion(1e-4, w) for w in (0.0, 2e-5, 5e-5, 8e-5)]
    assert all(ts[i] < ts[i + 1] for i in range(len(ts) - 1)), "θ must rise with adhesion"
    # geometry measurement matches analytic for a known cap (d/r = cos θ)
    for cos_t in (0.2, 0.5, 0.8):
        r = 7.5e-6
        c = np.array([cos_t * r, 0.0, 0.0])             # centre at distance d=cos_t·r
        th = angle_from_doublet_geometry(c, r, np.zeros(3), np.array([1.0, 0, 0]))
        assert abs(math.cos(th) - cos_t) < 1e-9, f"geom cos θ mismatch {cos_t}"
    # triplet symmetric (all equal) → η=1 → cos(φ/2)=0.5 → φ=120°
    assert abs(triplet_angle(1.0, 1.0, 1.0) - math.radians(120.0)) < 1e-9, "equal triplet → 120°"
    print("young_dupre sanity gate: PASS (boundaries, monotonicity, geometry, triplet-120°)")


if __name__ == "__main__":
    _sanity()

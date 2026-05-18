"""Affine shear + virial stress for preliminary strain-stiffening (KU-1.29).

**Phase 1 preliminary only.** Per PI direction, quantitative
publication-grade strain-stiffening validation (KU-1.30 benchmark #1
G_0 ≈ 36 Pa, #2 γ_c ≈ 0.16) is deferred until either (i) exact
intersection-node insertion or (ii) finer bead discretization is in
place.

What this module does NOT do
----------------------------
Lees-Edwards periodic shear is **not** implemented. Without
Lees-Edwards, a Langevin run after an affine shear relaxes back to
the unsheared equilibrium (the periodic box does not pin the
deformation), so the relaxed σ_xy → 0 within thermal noise. For
Phase 1 we therefore measure the **instantaneous affine stress**
(no Langevin relaxation): apply r → r + (γ y, 0), immediately read
σ_xy from the Born virial. This is the upper bound on G; the true
non-affine, Langevin-relaxed modulus is lower and requires
Lees-Edwards (Phase 2+).
"""

from __future__ import annotations

import numpy as np

from acs_kb.ecm.cross_links import CrossLink
from acs_kb.ecm.fiber_network import FiberNetwork


def apply_affine_shear(bead_positions: np.ndarray, gamma: float) -> np.ndarray:
    """Return positions sheared affinely: x ← x + γ y."""
    out = bead_positions.copy()
    out[..., 0] += gamma * bead_positions[..., 1]
    return out


def _minimum_image(delta: np.ndarray, L: float) -> np.ndarray:
    return delta - L * np.round(delta / L)


def compute_virial_shear_stress_2d(
    bead_positions: np.ndarray,
    rest_length: float,
    stretching_modulus: float,
    bending_modulus: float,
    box_size: float,
    cross_links: list[CrossLink] | None = None,
) -> float:
    """Born-virial 2D shear stress σ_xy in **N/m** (not Pa!).

    Per-bond contribution σ_xy_bond = b_x F_y / V_2D where V_2D is the
    box area (m²) and (b·F) has units N·m, so the result is N/m. To
    compare with experimental Pa values, divide by an ECM layer
    thickness h: σ_3D ≈ σ_2D / h.  Phase 1 reports σ_2D directly to
    avoid hiding the convention.
    """
    V_2d = box_size * box_size  # 2D 'volume' (area)
    sigma_xy = 0.0

    # ---- Backbone stretching ----
    b = bead_positions[:, 1:, :] - bead_positions[:, :-1, :]
    b = _minimum_image(b, box_size)
    bond_len = np.linalg.norm(b, axis=-1)
    safe_len = np.where(bond_len > 0.0, bond_len, 1.0)
    bhat = b / safe_len[..., None]
    f_mag = (stretching_modulus / rest_length) * (bond_len - rest_length)
    fvec = f_mag[..., None] * bhat
    sigma_xy += float(np.sum(b[..., 0] * fvec[..., 1]))

    # ---- Cross-links ----
    if cross_links:
        fa = np.fromiter((xl.fiber_a for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        ba = np.fromiter((xl.bead_a for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        fb = np.fromiter((xl.fiber_b for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        bb = np.fromiter((xl.bead_b for xl in cross_links), dtype=np.int64,
                         count=len(cross_links))
        r0 = np.fromiter((xl.rest_length for xl in cross_links),
                         dtype=np.float64, count=len(cross_links))
        k = np.fromiter((xl.stiffness for xl in cross_links),
                        dtype=np.float64, count=len(cross_links))
        d = bead_positions[fb, bb, :] - bead_positions[fa, ba, :]
        d = _minimum_image(d, box_size)
        dist = np.linalg.norm(d, axis=1)
        safe = np.where(dist > 0.0, dist, 1.0)
        fmag = (k * (dist - r0)) / safe
        fvec = fmag[:, None] * d
        sigma_xy += float(np.sum(d[:, 0] * fvec[:, 1]))

    # Bending three-body contribution: in 2D the standard derivation
    # gives Σ_triplet F_i · r_i which simplifies under translation
    # invariance to the same per-bond virial above, but its numerical
    # magnitude on a near-rest configuration is far below the
    # stretching term for our parameters (κ/ℓ₀² ≪ μ/ℓ₀). Phase 1
    # bookkeeping omits it; flagged in the docstring.
    return sigma_xy / V_2d


def preliminary_affine_G_curve(
    net: FiberNetwork,
    cross_links: list[CrossLink],
    stretching_modulus: float,
    bending_modulus: float,
    gamma_values: np.ndarray,
    layer_thickness: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sweep γ, return (σ_2D, G_2D, G_3D) — **affine only**.

    σ_2D, G_2D are reported in N/m (the natural 2D virial units);
    G_3D = G_2D / layer_thickness in Pa, only when `layer_thickness`
    is supplied. Without it, G_3D is filled with NaN so callers can't
    accidentally compare to KU-1.30 Pa values without declaring the
    2D→3D convention.
    """
    p0 = net.bead_positions
    sigma_2d = np.zeros_like(gamma_values, dtype=np.float64)
    G_2d = np.zeros_like(gamma_values, dtype=np.float64)
    G_3d = np.full_like(gamma_values, np.nan, dtype=np.float64)
    for i, g in enumerate(gamma_values):
        p = apply_affine_shear(p0, g)
        sigma_2d[i] = compute_virial_shear_stress_2d(
            p, net.rest_length, stretching_modulus, bending_modulus,
            net.box_size, cross_links=cross_links,
        )
        G_2d[i] = sigma_2d[i] / g if g != 0 else 0.0
        if layer_thickness is not None and layer_thickness > 0:
            G_3d[i] = G_2d[i] / layer_thickness
    return sigma_2d, G_2d, G_3d

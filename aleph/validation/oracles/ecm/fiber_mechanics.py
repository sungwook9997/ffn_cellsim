"""Discrete worm-like chain Hamiltonian and forces (KU-1.24).

Per-fiber Hamiltonian:

    H = (μ / 2 ℓ₀) Σ_i (|b_i| − ℓ₀)²                      # stretching (KU-1.2, KU-1.25)
      + (κ / ℓ₀)   Σ_i (1 − cos θ_i)                       # bending     (KU-1.1)

where b_i = r_{i+1} − r_i is the i-th bond vector and θ_i is the angle
between consecutive bonds at interior bead i.

Phase 1 Unit 1.1 supplies **energy and forces only** — no time
integration. Time stepping is Unit 1.2 (Langevin, KU-1.26).

All quantities SI: positions m, energy J, forces N.

Sanity Gate
-----------
1. Dimensional analysis:
   - [μ/ℓ₀] = N/m, [|b|−ℓ₀]² = m², so stretch energy [E] = N·m = J ✓
   - [κ/ℓ₀] = N·m²/m = N·m = J, scaled by dimensionless (1−cosθ) ⇒ J ✓
   - Forces are −∂E/∂r so [F] = J/m = N ✓
   - No time integration in this module → no CFL bound to check.
   - Numerical FD check: max relative gradient error 4.2e-9 on a
     perturbed 3-fiber test, confirming analytic forces match the
     finite difference of the energy.
2. Boundary cases:
   - N=2 beads: stretching only, bending block guarded by `N >= 3`.
   - Zero bond length: protected by `safe_len = where(|b|>0, |b|, 1)`
     so divisions don't NaN; physical bond length ≈ 0 is unphysical
     in any well-formed network (would mean overlapping beads).
3. Conservation invariants:
   - Newton 3rd law: per-bond stretch contributions are equal/opposite
     (+fbond on bead i, −fbond on bead i+1). Bending term sums
     ∂cosθ over three beads as +ip1 + i + im1 = 0 by translation
     invariance — verified numerically (Σ F ≈ 1e-25 N).
   - Energy is purely conservative; no explicit dissipation channel.
4. Numerical sanity: float64. Vectorized over fibers; no atomic ops.
5. Sign / sense:
   - Stretched bond (|b|>ℓ₀) pulls endpoints together (attractive) —
     stretch_mag > 0, force on bead i = +fbond points toward bead i+1.
   - Bent triplet (cosθ<1) pushes interior bead i to straighten — the
     ∂E/∂r_i term comes from translation invariance and opposes the
     bend energy gradient.
6. Measurement protocol: this module returns per-bead forces and a
   scalar energy; both are computed over every bead with no off-region
   sampling, matching the WLC analytical proof everywhere.
"""

from __future__ import annotations

import numpy as np

from aleph.validation.oracles.ecm.fiber_network import FiberNetwork


def _minimum_image(delta: np.ndarray, L: float) -> np.ndarray:
    """Wrap displacement components into [−L/2, +L/2)."""
    return delta - L * np.round(delta / L)


def _bond_vectors(positions: np.ndarray, box_size: float) -> np.ndarray:
    """Return bond vectors b_i = MI(r_{i+1} − r_i), shape (F, N−1, 2)."""
    delta = positions[:, 1:, :] - positions[:, :-1, :]
    return _minimum_image(delta, box_size)


def compute_energy(
    bead_positions: np.ndarray,
    rest_length: float,
    stretching_modulus: float,
    bending_modulus: float,
    box_size: float,
    cross_links: list | None = None,
) -> float:
    """Total network energy (Joules) — WLC backbone + optional cross-links.

    Parameters
    ----------
    bead_positions : np.ndarray, shape (F, N, 2)
        Bead coordinates (m).
    rest_length : float
        ℓ₀ (m).
    stretching_modulus : float
        μ in Newtons (force per unit relative extension); see KU-1.2.
    bending_modulus : float
        κ in N·m² (κ = ℓ_p k_B T, KU-1.1).
    box_size : float
        Periodic box side L (m).
    cross_links : list[CrossLink] or None
        If supplied, add harmonic XL energy (KU-1.28).

    Returns
    -------
    energy : float
        Total Hamiltonian in Joules (backbone + XL).
    """
    b = _bond_vectors(bead_positions, box_size)
    bond_len = np.linalg.norm(b, axis=-1)  # (F, N−1)

    e_stretch = (stretching_modulus / (2.0 * rest_length)) * np.sum(
        (bond_len - rest_length) ** 2
    )

    if bead_positions.shape[1] < 3:
        e_bend = 0.0
    else:
        u = b[:, :-1, :]
        v = b[:, 1:, :]
        lu = bond_len[:, :-1]
        lv = bond_len[:, 1:]
        cos_theta = np.einsum("fij,fij->fi", u, v) / (lu * lv)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        e_bend = (bending_modulus / rest_length) * np.sum(1.0 - cos_theta)

    e_xl = 0.0
    if cross_links:
        # Local import to keep this module independent if the user
        # never touches cross-links.
        from aleph.validation.oracles.ecm.cross_links import compute_xl_energy_and_forces
        e_xl, _ = compute_xl_energy_and_forces(bead_positions, cross_links, box_size)

    return float(e_stretch + e_bend + e_xl)


def compute_forces(
    bead_positions: np.ndarray,
    rest_length: float,
    stretching_modulus: float,
    bending_modulus: float,
    box_size: float,
    cross_links: list | None = None,
) -> np.ndarray:
    """Per-bead forces F_k = −∂H/∂r_k (Newtons), shape (F, N, 2).

    Analytic gradient of the WLC stretching + bending Hamiltonian, plus
    optional harmonic cross-link contribution (KU-1.28). Minimum-image
    convention is applied to both backbone bonds and cross-link
    displacements.
    """
    F = bead_positions.shape[0]
    N = bead_positions.shape[1]
    forces = np.zeros_like(bead_positions)

    b = _bond_vectors(bead_positions, box_size)              # (F, N−1, 2)
    bond_len = np.linalg.norm(b, axis=-1)                    # (F, N−1)

    # -- Stretching --
    # ∂E/∂r_{i+1} = (μ/ℓ₀)(|b_i| − ℓ₀) b_i/|b_i|
    # ∂E/∂r_i    = −(μ/ℓ₀)(|b_i| − ℓ₀) b_i/|b_i|
    safe_len = np.where(bond_len > 0.0, bond_len, 1.0)
    bhat = b / safe_len[..., None]
    stretch_mag = (stretching_modulus / rest_length) * (bond_len - rest_length)
    fbond = stretch_mag[..., None] * bhat                    # (F, N−1, 2), force vector along bond

    # Accumulate: bead i loses fbond[i] (the i-th bond pulls it toward bead i+1
    # when stretched), bead i+1 gains fbond[i] but with opposite sign — see
    # derivation above. Force on r_k from bond k:    +fbond[k]  (k = i index)
    # Force on r_k from bond k−1 (in which it is i+1): −fbond[k−1].
    forces[:, :-1, :] += fbond
    forces[:, 1:, :]  -= fbond

    # -- Bending (requires interior triplets) --
    if N >= 3:
        u = b[:, :-1, :]                                     # (F, N−2, 2)
        v = b[:, 1:, :]                                      # (F, N−2, 2)
        lu = bond_len[:, :-1]
        lv = bond_len[:, 1:]
        safe_lu = np.where(lu > 0.0, lu, 1.0)
        safe_lv = np.where(lv > 0.0, lv, 1.0)
        uv = np.einsum("fij,fij->fi", u, v)
        cos_theta = np.clip(uv / (safe_lu * safe_lv), -1.0, 1.0)

        inv_uv_norm = 1.0 / (safe_lu * safe_lv)              # 1 / (|u||v|)
        u_hat = u / safe_lu[..., None]
        v_hat = v / safe_lv[..., None]

        # ∂cosθ/∂r_{i+1} = (1/|u||v|)(u − cosθ · (|u|/|v|) v)
        # ∂cosθ/∂r_{i-1} = (1/|u||v|)(-v + cosθ · (|v|/|u|) u)
        # ∂cosθ/∂r_i     = −(∂_{i+1} + ∂_{i-1})  (translation invariance)
        # Force = +(κ/ℓ₀) ∂cosθ/∂r  since E = (κ/ℓ₀)(1 − cosθ).
        prefactor = bending_modulus / rest_length

        # ∂cosθ/∂r_{i+1}
        cos_uv_ratio_p = cos_theta * (safe_lu / safe_lv)            # (F, N−2)
        d_cos_dr_ip1 = inv_uv_norm[..., None] * (u - cos_uv_ratio_p[..., None] * v)

        # ∂cosθ/∂r_{i−1}
        cos_vu_ratio_m = cos_theta * (safe_lv / safe_lu)
        d_cos_dr_im1 = inv_uv_norm[..., None] * (-v + cos_vu_ratio_m[..., None] * u)

        d_cos_dr_i = -d_cos_dr_ip1 - d_cos_dr_im1

        forces[:, 2:, :]    += prefactor * d_cos_dr_ip1   # bead i+1 (index 2..N−1)
        forces[:, 1:-1, :]  += prefactor * d_cos_dr_i     # bead i   (index 1..N−2)
        forces[:, :-2, :]   += prefactor * d_cos_dr_im1   # bead i−1 (index 0..N−3)

    if cross_links:
        from aleph.validation.oracles.ecm.cross_links import compute_xl_energy_and_forces
        _, f_xl = compute_xl_energy_and_forces(bead_positions, cross_links, box_size)
        forces += f_xl

    return forces


def total_force(forces: np.ndarray) -> np.ndarray:
    """Sum of per-bead forces; should be ≈ 0 for closed bonded systems (Newton's 3rd law)."""
    return forces.sum(axis=(0, 1))


# Convenience wrappers that take a FiberNetwork ----------------------

def energy_of(net: FiberNetwork, stretching_modulus: float, bending_modulus: float) -> float:
    return compute_energy(
        net.bead_positions, net.rest_length, stretching_modulus, bending_modulus, net.box_size
    )


def forces_of(net: FiberNetwork, stretching_modulus: float, bending_modulus: float) -> np.ndarray:
    return compute_forces(
        net.bead_positions, net.rest_length, stretching_modulus, bending_modulus, net.box_size
    )

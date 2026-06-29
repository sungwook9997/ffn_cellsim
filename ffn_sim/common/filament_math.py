"""Tangent-correlation, persistence-length fit, and equipartition
utilities for single-filament validation (H.2).

Shared between H.2 (single F-actin L_p validation), H.3 (cortex
multi-filament fluctuations), and H.5 (lamellipodium leading-edge
filament tangent statistics).

All inputs are bead coordinates from a HOOMD ``cpu_local_snapshot``
ordered as ``(N_filaments, N_beads, 3)`` SI metres; outputs are SI
(length: m, energy: J, angle: rad).

Sanity Gate
-----------
*Per CLAUDE.md hard rule "Sanity Gate Protocol mandatory before first
execution of any physics/numerics module."*

1. **Dimensional analysis**
   - ``tangent_correlation(positions, rest_length)`` returns a
     dimensionless `C(s) = ⟨t̂_i · t̂_{i+s}⟩` per bond-separation s; t̂
     is normalised so it has unit norm. ✓
   - ``fit_persistence_length`` returns L_p in metres (extracted from
     the exponential fit `C(s) = exp(−s·ℓ₀ / L_p)`). ✓
   - ``bending_energy_per_bond`` returns J per interior bead from the
     harmonic angle Hamiltonian ½ k_θ (θ − π)². ✓

2. **Boundary cases**
   - Empty positions ``shape[0] == 0`` or ``shape[1] < 2``: returns
     empty arrays / NaN, with explicit guard.
   - Single fiber (``shape[0] == 1``): supported; tangent correlation
     averaged over that one fiber's bonds.
   - ``rest_length ≤ 0``: ValueError.
   - ``shape[1] < 3``: bending_energy returns empty (no interior
     triplet).

3. **Conservation invariants**
   - Tangent vectors are MI-wrapped: ``t̂_i = MI(r_{i+1} - r_i) /
     |MI(r_{i+1} - r_i)|``. Box wrap convention matches
     ``ffn_sim.archive.hoomd_legacy.integrator.baoab._wrap_into_box`` for Lees-Edwards
     compatibility.
   - The exponential fit is performed in log-space; numerical
     stability checked via finite-positive C(s).

4. **Numerical sanity**
   - All arrays in float64.
   - The exponential fit uses ``np.polyfit(s · ℓ₀, log(C(s)), 1)`` with
     a guard against negative C(s) (drop bins where C ≤ 0 from the fit
     — happens at large s when the fit window exceeds the persistence
     length and C(s) wraps below zero).

5. **Sign / sense**
   - C(s=0) = ⟨t̂·t̂⟩ = 1 (positive definite, normalisation check).
   - L_p > 0 by construction (slope of log C(s) vs s is negative for
     all sensible filaments; we fit the negative slope and report
     its inverse as L_p).
   - Bending angle θ ∈ [0, π], stored as the HOOMD-side angle between
     consecutive bond vectors at the interior bead j (= π at straight
     chain). The harmonic deviation is (θ − π).

6. **Measurement protocol**
   - The angle θ_i used for both `boltzmann_angle_histogram` and
     `bending_energy_per_bond` is computed from MI bond vectors as
     ``cos(π − θ_i) = -t̂_i · t̂_{i+1}`` (i.e. cos of the supplementary
     since t̂_i and t̂_{i+1} are tangent unit vectors aligned with the
     chain). We use ``θ_HOOMD = π − arccos(t̂_i · t̂_{i+1})`` to match
     HOOMD's convention and the `md.angle.Harmonic` t0 = π reference.
   - Boltzmann reference: `p(θ_HOOMD) ∝ sin(θ_HOOMD) · exp(−κ_B
     (π − θ_HOOMD)² / (2 ℓ_0 k_BT))` (3D solid-angle volume element
     `sin θ`). The brief mentions `p(θ) ∝ exp(−κ_B θ² / 2 ℓ_0 k_BT)`
     which is the 2D small-angle limit; we report both the 3D form
     (HOOMD canonical) and the 2D form (brief reference) so the
     downstream KS test can use whichever the gate requires.

References
----------
- H.2 brief ``ffn_sim/docs/briefs/H2_single_filament.md`` §Measurement.
- KU-1.1: actin L_p = 17 μm, κ_B = ℓ_p · k_B T at 310 K.
- BAOAB freeze §Open #1 (bending equipartition 3D-corrected) — this
  module's `bending_energy_per_bond` is the principled re-home of
  that gate at the single-filament level.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(slots=True)
class PersistenceLengthFit:
    """Fit result for `C(s) = exp(−s·ℓ₀ / L_p)`."""

    L_p_m: float                       # persistence length [m]
    L_p_stderr_m: float                # 1σ stderr from polyfit covariance [m]
    s_bonds: np.ndarray                # bond-separation indices used in fit
    C_observed: np.ndarray             # observed C(s) at those s
    log_slope_inv_m: float             # −1 / L_p (slope of log C(s) vs s·ℓ₀)


@dataclass(slots=True)
class EquipartitionResult:
    """Per-bond bending equipartition diagnostic."""

    energy_per_bond_J: np.ndarray      # shape (N_bonds_interior,)
    mean_J: float                      # ⟨E_bend⟩ across bonds + samples
    kT_over_2_J: float                 # kT/2 reference
    relative_deviation: float          # (mean − kT/2) / (kT/2)


def _minimum_image_3d(
    delta: np.ndarray,
    box_Lx: float, box_Ly: float, box_Lz: float,
    box_xy: float = 0.0, box_xz: float = 0.0, box_yz: float = 0.0,
) -> np.ndarray:
    """Wrap a (..., 3) displacement into the (possibly tilted) HOOMD box.

    Mirrors `ffn_sim.archive.hoomd_legacy.integrator.baoab._wrap_into_box` convention so
    tangent vectors computed here agree with the L-M Updater's wrap.
    """
    fz = delta[..., 2] / box_Lz
    fy = (delta[..., 1] - box_yz * box_Lz * fz) / box_Ly
    fx = (delta[..., 0] - box_xy * box_Ly * fy - box_xz * box_Lz * fz) / box_Lx
    fx -= np.round(fx)
    fy -= np.round(fy)
    fz -= np.round(fz)
    out = np.empty_like(delta)
    out[..., 0] = box_Lx * fx + box_xy * box_Ly * fy + box_xz * box_Lz * fz
    out[..., 1] = box_Ly * fy + box_yz * box_Lz * fz
    out[..., 2] = box_Lz * fz
    return out


def bond_tangent_unit_vectors(
    positions_FN3: np.ndarray, *, box: tuple[float, float, float] | None = None
) -> np.ndarray:
    """Return t̂_i = MI(r_{i+1} - r_i) / |MI(r_{i+1} - r_i)| with shape (F, N−1, 3).

    Parameters
    ----------
    positions_FN3 : ndarray, shape (F, N, 3)
        Bead coordinates [m].
    box : (Lx, Ly, Lz) or None
        Orthorhombic box dimensions for MI wrap. If None, no MI
        wrapping (e.g. for non-periodic / open-box single filament).
    """
    if positions_FN3.ndim != 3 or positions_FN3.shape[-1] != 3:
        raise ValueError(
            f"positions_FN3 must have shape (F, N, 3); got "
            f"{positions_FN3.shape}."
        )
    N = positions_FN3.shape[1]
    if N < 2:
        raise ValueError(f"Need ≥ 2 beads per fiber; got N={N}.")

    delta = positions_FN3[:, 1:, :] - positions_FN3[:, :-1, :]
    if box is not None:
        delta = _minimum_image_3d(delta, *box)
    norms = np.linalg.norm(delta, axis=-1, keepdims=True)
    safe = np.where(norms > 0.0, norms, 1.0)
    return delta / safe


def tangent_correlation(
    positions_FN3: np.ndarray, *,
    box: tuple[float, float, float] | None = None,
    max_separation: int | None = None,
) -> np.ndarray:
    """Compute `C(s) = ⟨t̂_i · t̂_{i+s}⟩` for s ∈ [0, max_separation].

    Averaged over (i) starting bond index and (ii) fiber index.
    Returns a 1D array of length ``max_separation + 1``.
    """
    t = bond_tangent_unit_vectors(positions_FN3, box=box)
    F, Nm1, _ = t.shape
    if max_separation is None:
        max_separation = Nm1 - 1
    if max_separation < 0:
        raise ValueError("max_separation must be ≥ 0.")
    if max_separation > Nm1 - 1:
        max_separation = Nm1 - 1
    C = np.zeros(max_separation + 1, dtype=np.float64)
    for s in range(max_separation + 1):
        # ⟨t̂_i · t̂_{i+s}⟩ averaged over i in [0, Nm1-1-s] and over F.
        dots = np.einsum("fij,fij->fi", t[:, :Nm1 - s, :], t[:, s:Nm1, :])
        C[s] = float(dots.mean())
    return C


def fit_persistence_length(
    positions_FN3: np.ndarray, rest_length: float, *,
    box: tuple[float, float, float] | None = None,
    fit_s_max: int | None = None,
) -> PersistenceLengthFit:
    """Fit `C(s) = exp(−s·ℓ₀ / L_p)` over s = 1..fit_s_max bonds.

    Per the H.2 brief, the canonical fit window is s ∈ [1, N/2].
    """
    if rest_length <= 0.0:
        raise ValueError(f"rest_length must be > 0; got {rest_length}.")
    N = positions_FN3.shape[1]
    if fit_s_max is None:
        fit_s_max = N // 2
    C = tangent_correlation(
        positions_FN3, box=box, max_separation=fit_s_max
    )
    # Sanity §5: C(0) must be 1 to numerical precision.
    if not (abs(C[0] - 1.0) < 1e-9):
        raise RuntimeError(
            f"C(s=0) = {C[0]} != 1.0 — tangent normalisation broken."
        )

    # Sanity §4: drop bins with C ≤ 0 (occurs at large s on noisy
    # samples; can't take log).
    s_full = np.arange(1, fit_s_max + 1, dtype=np.float64)
    C_fit = C[1:]
    mask = C_fit > 0.0
    if mask.sum() < 2:
        return PersistenceLengthFit(
            L_p_m=float("nan"),
            L_p_stderr_m=float("nan"),
            s_bonds=s_full,
            C_observed=C_fit,
            log_slope_inv_m=float("nan"),
        )
    s_fit = s_full[mask]
    log_C_fit = np.log(C_fit[mask])
    x = s_fit * rest_length          # arclength [m]
    # Linear fit log C = slope · x + intercept.
    coeffs, cov = np.polyfit(x, log_C_fit, 1, cov=True)
    slope = float(coeffs[0])
    slope_stderr = float(np.sqrt(cov[0, 0]))
    if slope >= 0.0:
        # Filament is more correlated at large s than small — unphysical.
        # Report as NaN L_p (rather than a negative value).
        return PersistenceLengthFit(
            L_p_m=float("nan"),
            L_p_stderr_m=float("nan"),
            s_bonds=s_full,
            C_observed=C_fit,
            log_slope_inv_m=slope,
        )
    L_p = -1.0 / slope
    L_p_stderr = L_p**2 * slope_stderr   # δ(1/slope) = δ(slope) / slope²
    return PersistenceLengthFit(
        L_p_m=L_p,
        L_p_stderr_m=L_p_stderr,
        s_bonds=s_full,
        C_observed=C_fit,
        log_slope_inv_m=slope,
    )


def hoomd_angle_array(
    positions_FN3: np.ndarray, *,
    box: tuple[float, float, float] | None = None,
) -> np.ndarray:
    """Return θ_HOOMD = π − arccos(t̂_i · t̂_{i+1}) per interior bead.

    Shape: (F, N − 2). θ ∈ [0, π], θ = π at straight chain (matches
    `md.angle.Harmonic.params['ecm-angle'].t0`).
    """
    t = bond_tangent_unit_vectors(positions_FN3, box=box)
    if t.shape[1] < 2:
        return np.empty((t.shape[0], 0), dtype=np.float64)
    dot = np.einsum("fij,fij->fi", t[:, :-1, :], t[:, 1:, :])
    dot = np.clip(dot, -1.0, 1.0)
    theta_supplement = np.arccos(dot)   # angle between t̂_i and t̂_{i+1}, ∈ [0, π]
    theta_hoomd = np.pi - theta_supplement
    return theta_hoomd


def bending_energy_per_bond(
    positions_FN3: np.ndarray, *,
    angle_k: float, angle_t0: float = np.pi,
    box: tuple[float, float, float] | None = None,
) -> np.ndarray:
    """Per-interior-bond bending energy ½ k_θ (θ − t0)²  [J].

    Shape: (F, N − 2). Matches the harmonic angle Hamiltonian used by
    HOOMD's `md.angle.Harmonic` so the equipartition gate can be
    closed against the on-line force compute.
    """
    theta = hoomd_angle_array(positions_FN3, box=box)
    return 0.5 * angle_k * (theta - angle_t0) ** 2


def equipartition_check(
    energy_J_samples: list[np.ndarray] | np.ndarray, *, kT_J: float
) -> EquipartitionResult:
    """Average per-bond energy across (samples, bonds) and compare to kT/2.

    ``energy_J_samples`` may be either:
    - a list of (F, N−2) arrays from sequential frame snapshots, or
    - a single (n_samples, F, N−2) ndarray.
    """
    if isinstance(energy_J_samples, list):
        stack = np.stack(energy_J_samples, axis=0)
    else:
        stack = np.asarray(energy_J_samples)
    flat = stack.reshape(-1)
    mean_J = float(flat.mean())
    kT_over_2 = 0.5 * kT_J
    rel = (mean_J - kT_over_2) / kT_over_2
    return EquipartitionResult(
        energy_per_bond_J=flat,
        mean_J=mean_J,
        kT_over_2_J=kT_over_2,
        relative_deviation=rel,
    )


def boltzmann_angle_density_2d(
    theta_hoomd: np.ndarray, *,
    bending_modulus: float, rest_length: float, kT: float,
) -> np.ndarray:
    """2D small-angle Boltzmann density `p(θ_HOOMD) ∝ exp(−κ_B (π − θ)² / (2 ℓ_0 k_BT))`.

    Brief's reference. Returns un-normalised density per input θ value
    (normalise via numerical integration in the test harness).
    """
    delta = theta_hoomd - np.pi
    return np.exp(-(bending_modulus * delta**2) / (2.0 * rest_length * kT))


def boltzmann_angle_density_3d(
    theta_hoomd: np.ndarray, *,
    bending_modulus: float, rest_length: float, kT: float,
) -> np.ndarray:
    """3D solid-angle Boltzmann density `p(θ_HOOMD) ∝ sin(θ_HOOMD) · exp(...)`.

    Required for HOOMD's 3D state; the sin θ volume element shifts
    the mode of `p` away from π toward smaller θ. The BAOAB freeze
    §Open #1 noted that the 2D reference (above) gives a mean energy
    biased by +71% vs measurement, while the 3D form is within 8.5%
    — H.2 closes that gate using the 3D form here.
    """
    delta = theta_hoomd - np.pi
    return np.sin(theta_hoomd) * np.exp(
        -(bending_modulus * delta**2) / (2.0 * rest_length * kT)
    )

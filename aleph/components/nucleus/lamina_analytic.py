"""Closed-form lamina references — the analytic ground truth for the I2 deformable-nucleus solver.

Host-side acceptance oracles (pure NumPy, NO Warp, NO simulation runtime): the Warp-CUDA envelope
kernels are gated AGAINST these (founding rule — analytic ground truth is primary; peer engines are
on-demand cross-checks only). Two physics families:

  (A) Seung–Nelson dihedral BENDING — reuses the membrane pattern: per interior edge shared by two
      outward-wound triangles E_e = κ̃(1 − n̂₁·n̂₂); Σ over edges = discrete Helfrich/Willmore
      (∮(κ/2)(2H)²dA → 8πκ on a sphere). κ̃ is CALIBRATED to the resting mesh (κ̃ = 8π·κ_NE/Σ_ref) —
      the icosphere dihedral sum Σ_ref ≈ 7.50 (NOT the ideal-lattice 4√3π; measured, resolution- and
      R-independent). Force = −∂E/∂x, FD-verified (the sign arbiter).

  (B) framework-#6 LAMIN SPLIT areal tension — a strain-stiffening bilinear that resolves the lumped
      chromatin→lamin into three physical contributions:
          ε ≤ knee :  soft  = CHROMATIN (internal polymer net) + LAMIN-B (uniform baseline meshwork)
          ε > knee :  stiff = LAMIN-A/C (strain-stiffening lamina engages)
          ε > rupt :  torn  = 0  (EMERGENT envelope rupture; NOT scripted)
      The tangent areal modulus jumps ×(K_laminAC / K_soft) at the knee — the measured nuclear
      strain-stiffening (KB-3.B2.1/.2). Young-Laplace ΔP = 2σ/R closes the pressurised-shell balance.

Magnitude gates (E_nuc-dependent σ, rupture strain, knee) are INVALID until I0-B2 closes; the
STRUCTURAL gates here — 8πκ topology, FD-gradient sign, continuity at the knee, tangent-ratio jump
UP (stiffening not softening), rupture on/off, Young-Laplace — are magnitude-independent and pass now.

References:
  Seung & Nelson 1988 PRA 38:1005 (discrete bending); Helfrich 1973; Willmore energy 8πκ.
  Nuclear lamina strain-stiffening + rupture: Stephens 2017/2019, Denais 2016, Lammerding —
  KB-3.B2.1 (chromatin↔lamin ratio 1.4–5×), KB-3.B2.2 (knee ~10%), KB-3.B2.3/.4 (rupture).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

__all__ = [
    "LaminaParams",
    "dihedral_bending_energy",
    "dihedral_bending_forces",
    "calibrate_kappa_tilde",
    "sphere_willmore_energy",
    "lamina_tension",
    "lamina_tangent_modulus",
    "young_laplace_pressure",
    "is_ruptured",
]


# ---------------------------------------------------------------------------------------------------
# (A) dihedral bending — the exact host reference the Warp kernel force must match
# ---------------------------------------------------------------------------------------------------
def dihedral_bending_energy(
    verts: npt.NDArray[np.float64], hinges: npt.NDArray[np.int64], kappa_tilde: float
) -> float:
    """Discrete Helfrich energy ``E = Σ_hinges κ̃·(1 − n̂₁·n̂₂)`` [pN·µm]. Independent re-derivation
    (not a re-import of the kernel under test) — the sphere→8πκ and FD-gradient ground truth."""
    if hinges.shape[0] == 0:
        return 0.0
    xi = verts[hinges[:, 0]]
    xj = verts[hinges[:, 1]]
    xk = verts[hinges[:, 2]]
    xl = verts[hinges[:, 3]]
    n1 = np.cross(xj - xi, xk - xi)
    n2 = np.cross(xi - xj, xl - xj)
    n1 /= np.linalg.norm(n1, axis=1, keepdims=True) + 1e-300
    n2 /= np.linalg.norm(n2, axis=1, keepdims=True) + 1e-300
    c = np.sum(n1 * n2, axis=1)
    return float(kappa_tilde * np.sum(1.0 - c))


def dihedral_bending_forces(
    verts: npt.NDArray[np.float64], hinges: npt.NDArray[np.int64], kappa_tilde: float
) -> npt.NDArray[np.float64]:
    """Analytic bending force ``f = −∂E/∂x`` (Nv,3) [pN] — the EXACT formula the Warp
    ``helfrich_bending_kernel`` implements (i,j on the shared edge feel both triangles; k,l only
    their own). FD-verified against ``dihedral_bending_energy`` (the sign arbiter). Σf = 0."""
    f = np.zeros_like(verts)
    if hinges.shape[0] == 0:
        return f
    xi = verts[hinges[:, 0]]
    xj = verts[hinges[:, 1]]
    xk = verts[hinges[:, 2]]
    xl = verts[hinges[:, 3]]
    n1 = np.cross(xj - xi, xk - xi)
    n2 = np.cross(xi - xj, xl - xj)
    a1 = np.linalg.norm(n1, axis=1, keepdims=True)
    a2 = np.linalg.norm(n2, axis=1, keepdims=True)
    n1 = n1 / (a1 + 1e-300)
    n2 = n2 / (a2 + 1e-300)
    c = np.sum(n1 * n2, axis=1, keepdims=True)
    q1 = n2 - c * n1
    q2 = n1 - c * n2
    inv1 = kappa_tilde / (a1 + 1e-300)                    # κ̃/(2A₁)
    inv2 = kappa_tilde / (a2 + 1e-300)                    # κ̃/(2A₂)
    fi = -inv1 * np.cross(xk - xj, q1) - inv2 * np.cross(xj - xl, q2)
    fj = -inv1 * np.cross(xi - xk, q1) - inv2 * np.cross(xl - xi, q2)
    fk = -inv1 * np.cross(xj - xi, q1)
    fl = -inv2 * np.cross(xi - xj, q2)
    np.add.at(f, hinges[:, 0], fi)
    np.add.at(f, hinges[:, 1], fj)
    np.add.at(f, hinges[:, 2], fk)
    np.add.at(f, hinges[:, 3], fl)
    return f


def calibrate_kappa_tilde(
    kappa_ne: float, verts: npt.NDArray[np.float64], hinges: npt.NDArray[np.int64]
) -> float:
    """κ̃ = 8π·κ_NE / Σ_ref [pN·µm], with Σ_ref the raw dihedral sum of the RESTING mesh, so the
    discrete sphere reproduces the continuum Willmore energy 8πκ_NE. Grid-invariant (Σ_ref is a
    topological constant of the icosphere ≈ 7.50, not a fit) ⇒ the Magic-Number Block is satisfied."""
    sigma = dihedral_bending_energy(verts, hinges, 1.0)   # raw Σ(1 − n̂₁·n̂₂)
    if sigma <= 1e-9:
        raise ValueError("degenerate/flat mesh: Σ ≈ 0, cannot calibrate κ̃")
    return 8.0 * np.pi * float(kappa_ne) / float(sigma)


def sphere_willmore_energy(kappa_ne: float) -> float:
    """Continuum Willmore/Helfrich bending energy of a sphere ``8πκ`` [pN·µm] — the bending gate
    target (R-independent; catches winding/normal bugs)."""
    return float(8.0 * np.pi * kappa_ne)


# ---------------------------------------------------------------------------------------------------
# (B) framework-#6 lamin-split areal tension
# ---------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class LaminaParams:
    """Areal-elasticity parameters of the nuclear lamina (framework #6 split). All moduli are 2-D
    (areal) tensions [pN/µm] = E·h. Values close through ``params_i0b2.yaml`` (I0-B2) — do NOT tune.

    Attributes:
        k_chrom: chromatin internal-polymer-net areal modulus (soft, small-strain) [pN/µm].
        k_lamin_b: lamin-B uniform-meshwork baseline areal modulus (soft) [pN/µm].
        k_lamin_ac: lamin-A/C strain-stiffening areal modulus (stiff, past the knee) [pN/µm].
        knee_strain: areal strain at lamin-A/C engagement (KB-3.B2.2 ~0.10).
        eps_rupture: areal strain at EMERGENT envelope rupture (KB-3.B2.3/.4; GAP — do NOT tune).
        kappa_ne: envelope bending rigidity [pN·µm] (~20 kBT lamina+bilayer).
    """
    k_chrom: float
    k_lamin_b: float
    k_lamin_ac: float
    knee_strain: float
    eps_rupture: float
    kappa_ne: float

    @property
    def k_soft(self) -> float:
        """Sub-knee tangent = chromatin + lamin-B baseline [pN/µm]."""
        return self.k_chrom + self.k_lamin_b

    @property
    def stiffening_ratio(self) -> float:
        """Tangent-modulus jump at the knee K_laminAC / K_soft (strain-stiffening ⇒ > 1)."""
        return self.k_lamin_ac / self.k_soft


def lamina_tension(eps: npt.ArrayLike, p: LaminaParams) -> npt.NDArray[np.float64]:
    """Strain-stiffening areal tension σ(ε) [pN/µm], continuous at the knee, 0 past rupture.

    ``σ = K_soft·ε``                            (ε ≤ knee: chromatin + lamin-B)
    ``σ = K_soft·knee + K_laminAC·(ε−knee)``    (knee < ε ≤ rupture: lamin-A/C engaged)
    ``σ = 0``                                   (ε > rupture: torn)
    """
    e = np.asarray(eps, dtype=np.float64)
    soft = p.k_soft * e
    stiff = p.k_soft * p.knee_strain + p.k_lamin_ac * (e - p.knee_strain)
    sigma = np.where(e <= p.knee_strain, soft, stiff)
    sigma = np.where(e > p.eps_rupture, 0.0, sigma)
    return sigma


def lamina_tangent_modulus(eps: npt.ArrayLike, p: LaminaParams) -> npt.NDArray[np.float64]:
    """Tangent areal modulus dσ/dε [pN/µm]: K_soft below the knee, K_laminAC above, 0 if ruptured."""
    e = np.asarray(eps, dtype=np.float64)
    k = np.where(e <= p.knee_strain, p.k_soft, p.k_lamin_ac)
    k = np.where(e > p.eps_rupture, 0.0, k)
    return k


def young_laplace_pressure(sigma: float, radius: float) -> float:
    """Young-Laplace pressure jump across a spherical shell at areal tension σ: ``ΔP = 2σ/R`` [pN/µm²
    = Pa]. The pressurised-nucleus balance the nucleoplasm volume constraint must reproduce."""
    if radius <= 0.0:
        raise ValueError("radius must be positive")
    return float(2.0 * sigma / radius)


def is_ruptured(eps: npt.ArrayLike, eps_rupture: float) -> npt.NDArray[np.bool_]:
    """EMERGENT per-face rupture mask ``ε > eps_rupture`` — rupture appears past a sourced threshold
    and is absent below it (the falsifiable gate; the threshold is GAP, not tuned to a gate)."""
    return np.asarray(eps, dtype=np.float64) > float(eps_rupture)

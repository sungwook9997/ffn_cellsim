r"""Method-of-planes cortical-tension estimator for the ac/ engine — source / network / total split.

The assembled-cell activation path (:mod:`aleph.components.incumbent.myosin_activate`) reports actomyosin "cortical
tension" through :func:`~aleph.components.incumbent.myosin_activate._radial_contractile_force`, which sums only the myosin
force on actin nodes and projects it onto the radial direction: ``Σ_i F_i·r̂_i``. For a closed, isotropic
contractile shell that scalar **near-cancels** (the internal force vectors obey Newton's 3rd law, so their sum
is identically zero, and the *radial projection* averages to ~0 as ``1/√N`` by symmetry) even while the true
surface **tension** — a cut-force / stress, not a net vector force — is finite and positive. Reporting that
scalar as "tension" is the ``~1015×`` self-cancellation artifact flagged in the magnitude audit.

This module measures the physically-correct quantity by **method of planes**, reusing the validated FF
estimator :func:`aleph.laws.gamma_estimator.method_of_planes_gamma` (the same instrument the BAOAB-MD and FF
γ-floor experiments used — no re-implementation of the plane math here). It sums the tensile elements that
cross each diametral plane, ``γ = mean_n̂ |(1/2πR) Σ_{crossing} T·|û·n̂||``, and separates:

* ``gamma_source``  — the active crossbridge dipoles alone (what myosin injects);
* ``gamma_network`` — the actin inextensibility constraint tensions + crosslink tensions (what the network
  actually carries onward);
* ``gamma_total``   — every load-bearing element (the observable comparable to an AFM/Laplace band).

The net-force vector and the legacy radial scalar are still reported (``net_force_vector_pn``,
``summed_inward_radial_force_pn``) for **regression only** — they are diagnostics of internal-force closure,
NOT a tension, and are labelled as such so no downstream code mistakes them for γ again.

Units: engine convention µm-pN-s, so γ is in ``pN/µm`` (``1 pN/µm = 1e-6 N/m = 1e-3 mN/m``; see the band
provenance note in :mod:`aleph.laws.gamma_estimator`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import numpy.typing as npt

from aleph.laws.gamma_estimator import method_of_planes_gamma

__all__ = ["CorticalStressReport", "method_of_planes_from_elements", "measure_cortical_stress"]

# element family = (endpoint A positions (M,3), endpoint B positions (M,3), signed scalar tension (M,))
_ElementFamily = tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]


@dataclass(frozen=True)
class CorticalStressReport:
    """Separated cortical-stress observables for one cell state (all γ in pN/µm).

    Attributes:
        gamma_source_pn_per_um: method-of-planes γ over the active crossbridge dipoles alone.
        gamma_network_pn_per_um: method-of-planes γ over the actin-constraint + crosslink tensions.
        gamma_total_pn_per_um: method-of-planes γ over every load-bearing element.
        net_force_vector_pn: Σ of the per-node force vectors [pN] — a Newton-closure diagnostic (≈0 for
            internal forces), NOT a tension.
        summed_inward_radial_force_pn: LEGACY ``-Σ_i F_i·r̂_i`` [pN] (the old ``_radial_contractile_force``
            observable) — retained for regression comparison ONLY; it is NOT a surface tension.
        gamma_by_family_pn_per_um: per-family γ (``crossbridge`` / ``actin`` / ``crosslink`` / …).
        n_planes: number of diametral cut orientations averaged.
        R_um: great-circle radius used for the ``2πR`` cut circumference.
    """

    gamma_source_pn_per_um: float
    gamma_network_pn_per_um: float
    gamma_total_pn_per_um: float
    net_force_vector_pn: npt.NDArray[np.float64]
    summed_inward_radial_force_pn: float
    gamma_by_family_pn_per_um: Mapping[str, float]
    n_planes: int
    R_um: float


def method_of_planes_from_elements(
    rA: npt.NDArray[np.float64],
    rB: npt.NDArray[np.float64],
    tension: npt.NDArray[np.float64],
    R: float,
    *,
    n_planes: int = 64,
    centre: npt.NDArray[np.float64] | None = None,
) -> float:
    """Method-of-planes γ [pN/µm] over one element list (thin adapter over the validated FF estimator).

    Args:
        rA: (M, 3) element endpoint-A positions [µm].
        rB: (M, 3) element endpoint-B positions [µm].
        tension: (M,) signed scalar element tensions [pN] (+ tension, − compression).
        R: cell radius for the cut circumference ``2πR`` [µm].
        n_planes: number of near-isotropic diametral cut orientations.
        centre: cell centre; defaults inside the FF estimator to the endpoint mean.

    Returns:
        Mean over orientations of ``|γ|`` [pN/µm]; ``0.0`` if there are no elements.
    """
    rA = np.asarray(rA, dtype=np.float64)
    if rA.size == 0:
        return 0.0
    return float(
        method_of_planes_gamma(
            rA, np.asarray(rB, dtype=np.float64), np.asarray(tension, dtype=np.float64),
            R, n_planes=n_planes, centre=centre,
        )
    )


def _concat_families(
    families: Mapping[str, _ElementFamily], names: tuple[str, ...]
) -> _ElementFamily:
    """Concatenate the (rA, rB, tension) arrays of the named families (missing/empty families skipped)."""
    rAs, rBs, Ts = [], [], []
    for name in names:
        fam = families.get(name)
        if fam is None:
            continue
        rA, rB, T = fam
        rA = np.asarray(rA, dtype=np.float64)
        if rA.size == 0:
            continue
        rAs.append(rA)
        rBs.append(np.asarray(rB, dtype=np.float64))
        Ts.append(np.asarray(T, dtype=np.float64))
    if not rAs:
        empty = np.zeros((0, 3), dtype=np.float64)
        return empty, empty, np.zeros((0,), dtype=np.float64)
    return np.vstack(rAs), np.vstack(rBs), np.concatenate(Ts)


def measure_cortical_stress(
    R: float,
    families: Mapping[str, _ElementFamily],
    node_pos: npt.NDArray[np.float64],
    node_force: npt.NDArray[np.float64],
    *,
    centre: npt.NDArray[np.float64] | None = None,
    n_planes: int = 64,
    source_families: tuple[str, ...] = ("crossbridge",),
    network_families: tuple[str, ...] = ("actin", "crosslink"),
) -> CorticalStressReport:
    """Measure the separated cortical stress for one cell state.

    Args:
        R: cell radius [µm] for the cut circumference.
        families: element lists keyed by family name, each ``(rA, rB, tension)`` — e.g. ``crossbridge``
            (active myosin dipoles), ``actin`` (segment inextensibility tensions), ``crosslink``.
        node_pos: (N, 3) node positions [µm] (for the legacy net-radial / net-vector diagnostics).
        node_force: (N, 3) net force on each node [pN] (e.g. the myosin force on actin nodes).
        centre: cell centre; defaults to ``node_pos`` mean.
        n_planes: diametral cut orientations to average.
        source_families / network_families: which families roll into γ_source / γ_network.

    Returns:
        A :class:`CorticalStressReport`. ``gamma_*`` are the physical method-of-planes tensions;
        ``net_force_vector_pn`` and ``summed_inward_radial_force_pn`` are Newton-closure diagnostics, NOT γ.
    """
    node_pos = np.asarray(node_pos, dtype=np.float64)
    node_force = np.asarray(node_force, dtype=np.float64)
    if centre is None:
        centre = node_pos.mean(axis=0) if node_pos.size else np.zeros(3)
    centre = np.asarray(centre, dtype=np.float64)

    gamma_by_family = {
        name: method_of_planes_from_elements(*fam, R, n_planes=n_planes, centre=centre)
        for name, fam in families.items()
    }
    gamma_source = method_of_planes_from_elements(
        *_concat_families(families, source_families), R, n_planes=n_planes, centre=centre
    )
    gamma_network = method_of_planes_from_elements(
        *_concat_families(families, network_families), R, n_planes=n_planes, centre=centre
    )
    gamma_total = method_of_planes_from_elements(
        *_concat_families(families, tuple(families.keys())), R, n_planes=n_planes, centre=centre
    )

    # Legacy diagnostics (NOT a tension): net internal-force vector + its inward radial projection.
    if node_pos.size:
        rhat = node_pos - centre
        rhat = rhat / (np.linalg.norm(rhat, axis=1, keepdims=True) + 1e-30)
        summed_inward_radial = float(-np.sum(node_force * rhat))
        net_vec = node_force.sum(axis=0)
    else:
        summed_inward_radial = 0.0
        net_vec = np.zeros(3, dtype=np.float64)

    return CorticalStressReport(
        gamma_source_pn_per_um=gamma_source,
        gamma_network_pn_per_um=gamma_network,
        gamma_total_pn_per_um=gamma_total,
        net_force_vector_pn=net_vec,
        summed_inward_radial_force_pn=summed_inward_radial,
        gamma_by_family_pn_per_um=gamma_by_family,
        n_planes=n_planes,
        R_um=float(R),
    )

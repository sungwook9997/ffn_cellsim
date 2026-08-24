"""Static bending-energy magnitude — ANALYTIC ground-truth, always-on (no Cytosim binary).

The runnable-Cytosim parity oracle (`cytosim_parity.py`) validates the same anchor but is SKIPPED
offline (no `sim` binary). Per the project rule that FF is validated against analytic ground truth as
the PRIMARY always-on test (oracle = on-demand cross-check, not truth), this pins the bending kernel's
ABSOLUTE energy magnitude — complementing test_cytosim_bending (force=−∇E shape) and
test_bending_dispersion (κq⁴ force) and test_relax (κ/ζ·q⁴ rate), none of which check the magnitude.

Ground truth: a filament held as a constant-curvature circular arc of radius R and contour length L
has continuum bending energy

    E_bend = (κ / 2) ∫ (1/R)² ds = κ L / (2 R²)            [Euler–Bernoulli, constant curvature]

The FF interior-triple operator (NF2007 p9) sums only interior vertices, so the RAW discrete energy
under-counts by exactly (n−2)/(n−1) (two fewer curvature samples than the continuum integral); the
default end-correction g = (n−1)/(n−2) restores it. Both facts are checked here analytically.
"""

import numpy as np
import pytest

from aleph.validation.cytosim_parity import _arc_points
from aleph.laws.fiber_network import build_fiber_network
from aleph.laws.forces_warp import bending_energy

KAPPA = 20.0      # pN·µm² (the NF2007 microtubule value; magnitude check is κ-linear)
R_ARC = 5.0       # µm radius of curvature
L = 2.0           # µm contour length


def _arc_energy(R: float, n: int, kappa: float = KAPPA, *, end_correction: bool = True) -> float:
    net = build_fiber_network([_arc_points(R, L, n)], kappa=kappa)
    return bending_energy(net, end_correction=end_correction)


def test_corrected_energy_matches_kappaL_over_2Rsq():
    """End-corrected FF bending energy = κL/2R² to <0.5% at every resolution (incl. coarse n=5)."""
    E_an = 0.5 * KAPPA * L / R_ARC**2
    for n in (5, 9, 17, 33):
        assert _arc_energy(R_ARC, n) == pytest.approx(E_an, rel=5e-3)


def test_raw_interior_sum_undercounts_by_end_correction():
    """The RAW (paper-literal interior-triple) energy under-counts by exactly (n−2)/(n−1) — the
    documented end-correction, confirmed against the analytic κL/2R² (not against Cytosim)."""
    E_an = 0.5 * KAPPA * L / R_ARC**2
    for n in (5, 9, 17, 33, 65):
        raw_ratio = _arc_energy(R_ARC, n, end_correction=False) / E_an
        assert raw_ratio == pytest.approx((n - 2) / (n - 1), rel=2e-3)


def test_energy_converges_to_continuum_with_resolution():
    """Refining the discretization drives the corrected energy monotonically toward κL/2R²
    (continuum limit); the residual error shrinks."""
    E_an = 0.5 * KAPPA * L / R_ARC**2
    errs = [abs(_arc_energy(R_ARC, n) - E_an) / E_an for n in (5, 17, 65)]
    assert errs[0] > errs[1] > errs[2]        # monotone refinement
    assert errs[-1] < 1e-4                     # fine grid ≈ continuum


def test_one_over_Rsq_curvature_scaling():
    """Bending energy ∝ 1/R² — doubling the radius (halving the curvature) QUARTERS the energy."""
    e_R = _arc_energy(R_ARC, 33)
    e_2R = _arc_energy(2 * R_ARC, 33)
    assert e_2R / e_R == pytest.approx(0.25, rel=5e-3)


def test_energy_is_kappa_linear():
    """E_bend is linear in κ (same geometry, doubled rigidity → doubled energy)."""
    e1 = _arc_energy(R_ARC, 17, kappa=KAPPA)
    e2 = _arc_energy(R_ARC, 17, kappa=2 * KAPPA)
    assert e2 / e1 == pytest.approx(2.0, rel=1e-6)

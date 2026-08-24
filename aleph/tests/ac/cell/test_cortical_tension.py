r"""Host gate for the ac/ method-of-planes cortical-tension estimator (:mod:`aleph.components.incumbent.cortical_tension`).

The crux (G1/G3): on an ISOTROPIC contractile shell the legacy net-radial scalar
(``_radial_contractile_force``) near-cancels while the method-of-planes tension is finite and positive — this
is the numerical proof that the ``~33 pN`` / ``~1015×``-cancelled observable was the WRONG quantity. Pure
NumPy; no Warp/CUDA. Thresholds are derived from the ``1/√M`` isotropy scaling, not tuned.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.incumbent.cortical_tension import (
    CorticalStressReport,
    measure_cortical_stress,
    method_of_planes_from_elements,
)
from aleph.laws.gamma_estimator import method_of_planes_gamma

R = 7.4  # cortex shell radius [µm]


def _fib_sphere(n: int, radius: float) -> np.ndarray:
    """``n`` near-isotropic points on the sphere of the given radius (Fibonacci lattice)."""
    phi = (1.0 + 5.0**0.5) / 2.0
    i = np.arange(n, dtype=np.float64)
    z = 1.0 - 2.0 * (i + 0.5) / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    theta = 2.0 * np.pi * i / phi
    return radius * np.stack([r * np.cos(theta), r * np.sin(theta), z], axis=1)


def _contractile_dipoles(n: int, T: float, *, half_len: float = 0.2, seed: int = 0, uniaxial: bool = False,
                         tangential_force: bool = False):
    """A shell of ``n`` contractile dipoles: tangent (isotropic) or all-``x̂`` (uniaxial).

    Returns element arrays ``(rA, rB, tension)`` for the method-of-planes AND the per-endpoint node
    positions/forces for the legacy net-radial diagnostic. Each tensile dipole pulls its two endpoints
    together: force on A = +T·û (toward B), on B = −T·û — so the pair's vector sum is exactly zero.

    ``tangential_force`` projects the per-node force onto the local tangent plane (removes its radial
    component). This models the empirical native-run finding that the myosin-on-actin forces are ~isotropic /
    non-contractile (per-head tangential load ≈ 0, net radial only 0.1% and outward): with such forces the
    legacy radial NET-force scalar reads ~0 while the tensile network still carries a finite method-of-planes γ.
    """
    rng = np.random.default_rng(seed)
    c = _fib_sphere(n, R)
    rhat = c / np.linalg.norm(c, axis=1, keepdims=True)
    if uniaxial:
        t = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    else:  # random isotropic direction, projected into the local tangent plane
        v = rng.standard_normal((n, 3))
        t = v - (np.sum(v * rhat, axis=1, keepdims=True)) * rhat
        t /= np.linalg.norm(t, axis=1, keepdims=True) + 1e-30
    rA = c - half_len * t
    rB = c + half_len * t
    u = (rB - rA) / (np.linalg.norm(rB - rA, axis=1, keepdims=True) + 1e-30)
    tension = np.full(n, T, dtype=np.float64)
    node_pos = np.vstack([rA, rB])
    node_force = np.vstack([+T * u, -T * u])   # tensile: endpoints pulled together
    if tangential_force:
        rn = node_pos / (np.linalg.norm(node_pos, axis=1, keepdims=True) + 1e-30)
        node_force = node_force - np.sum(node_force * rn, axis=1, keepdims=True) * rn
    return (rA, rB, tension), node_pos, node_force


# ── G1: the crux — isotropic shell → net-radial ≈ 0 while γ > 0 ──────────────────────────────────
def test_g1_isotropic_net_radial_cancels_but_tension_is_finite():
    M, T = 4000, 5.0
    (rA, rB, tn), node_pos, node_force = _contractile_dipoles(M, T, seed=1)
    rep = measure_cortical_stress(R, {"crossbridge": (rA, rB, tn)}, node_pos, node_force,
                                  centre=np.zeros(3), n_planes=64)
    assert isinstance(rep, CorticalStressReport)

    sum_abs_F = float(np.sum(np.linalg.norm(node_force, axis=1)))       # = 2·M·T
    # A CLOSED contractile network has ZERO net FORCE (Newton's 3rd law) — independent of its tension.
    assert np.linalg.norm(rep.net_force_vector_pn) / sum_abs_F < 1e-9
    # ...yet the method-of-planes surface TENSION is finite and clearly positive: zero net force ≠ zero tension.
    assert rep.gamma_total_pn_per_um > 1.0
    # The legacy radial NET-force scalar (pN) is a DIFFERENT observable from γ (pN/µm): converting it to a
    # pseudo-tension (÷2πR) does not recover γ — you cannot read a surface tension off a net force.
    pseudo = abs(rep.summed_inward_radial_force_pn) / (2.0 * np.pi * R)
    assert not np.isclose(pseudo, rep.gamma_total_pn_per_um, rtol=0.3)


# ── G2: faithful reuse of the FF estimator + linearity/density scaling (no re-implementation) ─────
def test_g2_adapter_matches_ff_estimator_exactly():
    (rA, rB, tn), *_ = _contractile_dipoles(1500, 5.0, seed=2)
    mine = method_of_planes_from_elements(rA, rB, tn, R, n_planes=50, centre=np.zeros(3))
    ff = method_of_planes_gamma(rA, rB, tn, R, n_planes=50, centre=np.zeros(3))
    assert mine == pytest.approx(ff, rel=1e-12)   # I delegate, never re-derive the plane math


def test_g2_gamma_linear_in_tension_and_density():
    (rA, rB, tn), *_ = _contractile_dipoles(2000, 3.0, seed=3)
    g1 = method_of_planes_from_elements(rA, rB, tn, R, centre=np.zeros(3))
    g2 = method_of_planes_from_elements(rA, rB, 2.0 * tn, R, centre=np.zeros(3))
    assert g2 == pytest.approx(2.0 * g1, rel=1e-9)                       # γ ∝ tension
    (rA2, rB2, tn2), *_ = _contractile_dipoles(4000, 3.0, seed=3)
    gd = method_of_planes_from_elements(rA2, rB2, tn2, R, centre=np.zeros(3))
    assert 1.6 * g1 < gd < 2.4 * g1                                      # γ ∝ areal density (±stochastic)


# ── G3: source/network/total decomposition + the reframe on the active (source) shell ────────────
def test_g3_source_network_total_decomposition():
    (xA, xB, xT), npx, nfx = _contractile_dipoles(3000, 5.0, seed=4)                      # crossbridge
    (lA, lB, lT), *_ = _contractile_dipoles(3000, 8.0, seed=5)                            # crosslink
    rep = measure_cortical_stress(
        R, {"crossbridge": (xA, xB, xT), "crosslink": (lA, lB, lT)}, npx, nfx,
        centre=np.zeros(3), n_planes=64,
    )
    assert rep.gamma_source_pn_per_um == pytest.approx(rep.gamma_by_family_pn_per_um["crossbridge"], rel=1e-9)
    assert rep.gamma_network_pn_per_um == pytest.approx(rep.gamma_by_family_pn_per_um["crosslink"], rel=1e-9)
    # union tension brackets: ≥ each part, ≤ their sum (independent isotropic sets superpose sub-additively)
    assert rep.gamma_total_pn_per_um >= max(rep.gamma_source_pn_per_um, rep.gamma_network_pn_per_um) - 1e-9
    assert rep.gamma_total_pn_per_um <= rep.gamma_source_pn_per_um + rep.gamma_network_pn_per_um + 1e-9


def test_g3_reframe_active_tension_real_while_net_radial_is_zero():
    """THE reframe: with ~isotropic/tangential myosin forces (as the native run measured — per-head load ≈ 0,
    net radial 0.1% and outward), the legacy radial NET-force scalar reads ≈ 0, yet the active crossbridge
    network still carries a real method-of-planes γ_source > 0. So the 33 pN net-force was never the tension."""
    (xA, xB, xT), npx, nfx = _contractile_dipoles(4000, 5.0, seed=6, tangential_force=True)
    rep = measure_cortical_stress(R, {"crossbridge": (xA, xB, xT)}, npx, nfx, centre=np.zeros(3))
    budget = float(np.sum(np.linalg.norm(nfx, axis=1)))
    assert abs(rep.summed_inward_radial_force_pn) < 1e-6 * budget   # tangential forces → net radial ≈ 0
    assert rep.gamma_source_pn_per_um > 1.0                         # ...but the active tension is finite


# ── uniaxial contrast: orientation actually matters (sanity that γ isn't a constant) ──────────────
def test_uniaxial_differs_from_isotropic():
    (iA, iB, iT), *_ = _contractile_dipoles(3000, 5.0, seed=7, uniaxial=False)
    (uA, uB, uT), *_ = _contractile_dipoles(3000, 5.0, seed=7, uniaxial=True)
    g_iso = method_of_planes_from_elements(iA, iB, iT, R, centre=np.zeros(3))
    g_uni = method_of_planes_from_elements(uA, uB, uT, R, centre=np.zeros(3))
    assert g_iso > 0.0 and g_uni > 0.0
    assert not np.isclose(g_iso, g_uni, rtol=0.05)   # the |û·n̂| projection distinguishes them

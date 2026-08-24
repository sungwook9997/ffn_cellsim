"""FF unit system (ff/units) — pN·µm·s nondimensionalization (Stage 6c-a).

Anchors: SI↔FF converters round-trip; the grounded κ_actin matches its SI value (7.3e-26 N·m²);
the NF2007 mobility law gives a sane positive drag; and — the POINT of nondimensionalizing — the
reused implicit overdamped solver relaxes a bent actin fiber at FF scale using its DEFAULT DCM
tolerances (newton_tol=1e-10, cg_tol=1e-8), whereas the same solver false-converges (Δx≈0) on the
SI-scale fiber with those same default tolerances.
"""

import numpy as np
import pytest

from aleph.laws import units as U
from aleph.laws.constraints import segment_lengths
from aleph.laws.fiber_network import build_fiber_network
from aleph.laws.forces_warp import bending_energy, make_bending_force_fn


def test_converters_roundtrip():
    x = np.array([1.0, 2.5, -3.0])
    assert np.allclose(U.length_to_si(U.length_to_ff(x * U.UM)), x * U.UM)
    assert U.force_to_ff(5e-12) == pytest.approx(5.0)
    assert U.force_to_si(5.0) == pytest.approx(5e-12)


def test_kappa_actin_matches_si():
    """κ_actin in FF units, converted back to SI, equals the SI value used by the bending tests."""
    assert U.kappa_to_si(U.KAPPA_ACTIN) == pytest.approx(7.3e-26, rel=2e-2)
    # κ = k_B·T·L_p identity holds in FF units too
    assert U.KAPPA_ACTIN == pytest.approx(U.KBT * U.LP_ACTIN_UM)


def test_mobility_and_drag_sane():
    mu = U.fiber_mobility(0.7)                 # 0.7 µm actin fiber in MCF7 cytoplasm
    assert mu > 0.0
    gamma = U.fiber_point_drag(0.7, 7)
    assert gamma == pytest.approx(1.0 / (8 * mu))
    assert gamma > 0.0
    with pytest.raises(ValueError):           # L_h must exceed δ
        U.fiber_mobility(0.005)


def _bent_arc_ff(n=8):
    """A quarter-arc actin fiber built in FF units (µm), seg ≈ 0.1 µm."""
    th = np.linspace(0.0, np.pi / 2, n)
    R = 0.5  # µm
    nodes = np.stack([R * np.cos(th), R * np.sin(th), np.zeros_like(th)], axis=1)
    return build_fiber_network([nodes], kappa=U.KAPPA_ACTIN)


def test_implicit_relaxes_at_ff_scale_with_default_tolerances():
    """The payoff: at FF (pN·µm) scale the DCM implicit solver relaxes a bent fiber using its
    DEFAULT tolerances — no FF-specific tolerance tuning needed (cf. the SI false-convergence)."""
    from aleph.dcm.dcm_warp_implicit import implicit_overdamped_step

    net = _bent_arc_ff()
    L = float(np.sum(net.seg_rest))           # fiber length [µm]
    p = net.n_nodes - 1
    gamma = U.fiber_point_drag(L, p)          # grounded per-point drag [pN·s/µm]
    ffn = make_bending_force_fn(net)
    e0 = bending_energy(net)
    x = net.pos.reshape(-1).copy()
    for _ in range(5):                         # a few large implicit steps
        x, info = implicit_overdamped_step(x, ffn, gamma=gamma, dt=1.0, n_newton=2)
        net.pos = x.reshape(-1, 3)
    assert np.all(np.isfinite(x))
    assert bending_energy(net) < 0.1 * e0     # genuinely relaxed (DEFAULT tols, FF scale)


def test_si_scale_false_converges_with_default_tolerances():
    """Control: the IDENTICAL bent fiber built in SI metres, with the solver's default tolerances,
    barely moves — the false convergence the FF unit system is designed to avoid."""
    from aleph.dcm.dcm_warp_implicit import implicit_overdamped_step

    th = np.linspace(0.0, np.pi / 2, 8)
    R = 0.5 * U.UM
    nodes = np.stack([R * np.cos(th), R * np.sin(th), np.zeros_like(th)], axis=1)
    net = build_fiber_network([nodes], kappa=U.kappa_to_si(U.KAPPA_ACTIN))   # SI κ (7.3e-26)
    ffn = make_bending_force_fn(net)
    e0 = bending_energy(net)
    x = net.pos.reshape(-1).copy()
    gamma_si = U.fiber_point_drag(0.5, 7) * (U.PN * U.S / U.UM)              # γ in SI N·s/m
    for _ in range(5):
        x, info = implicit_overdamped_step(x, ffn, gamma=gamma_si, dt=1.0, n_newton=2)
        net.pos = x.reshape(-1, 3)
    assert bending_energy(net) > 0.9 * e0      # barely relaxed → false convergence at SI scale

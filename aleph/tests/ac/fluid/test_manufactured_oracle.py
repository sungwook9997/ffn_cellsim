"""Self-test of the manufactured-solution oracles (pure NumPy — no Warp/CUDA).

Uses the I0-B1 closure values so the eigenmode/consistency checks exercise the real ledger numbers:
c_v = mobility / S  must hold, i.e. 50 um^2/s = 5e-3 um^2/(Pa*s) / 1e-4 /Pa.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.fluid.manufactured import (
    cosine_mode,
    cosine_mode_decay_rate,
    mms_source,
    moving_interval_content,
    moving_interval_content_rate,
)

# I0-B1 ledger closure (um-pN-s-Pa), consistent by construction.
C_V = 50.0        # um^2/s
MOBILITY = 5.0e-3  # k/mu, um^2/(Pa*s)   (k=5e-6 um^2, mu=1e-3 Pa*s)
STORAGE_S = 1.0e-4  # 1/M, 1/Pa


def test_ledger_values_are_consistent() -> None:
    """c_v = mobility / S must hold or the eigenmode source will not vanish."""
    assert C_V == pytest.approx(MOBILITY / STORAGE_S, rel=1e-12)


def test_decay_rate_is_laplacian_eigenvalue() -> None:
    """lam = c_v |k|^2 for a 3-D cosine mode."""
    k = np.array([1.0, 2.0, 0.5])
    assert cosine_mode_decay_rate(k, C_V) == pytest.approx(C_V * (1.0 + 4.0 + 0.25))


def test_eigenmode_has_zero_mms_source_1d() -> None:
    """A true eigenmode is source-free: S p_t - mobility laplacian(p) = 0 to round-off."""
    k = 1.0
    lam = cosine_mode_decay_rate([k], C_V)
    x = np.linspace(0.0, np.pi / k, 501)  # zero-flux at both ends (sin=0)
    t = 0.3
    p = cosine_mode(x, t, amplitude=40.0, wavevector=[k], c_v=C_V)
    p_t = -lam * p
    p_lap = -(k**2) * p
    s = mms_source(STORAGE_S, MOBILITY, p_t, p_lap)
    assert np.max(np.abs(s)) < 1e-12 * np.max(np.abs(p))


def test_eigenmode_shape_3d() -> None:
    """3-D cosine mode returns one value per point and decays in time at lam."""
    k = np.array([1.0, 1.0, 1.0])
    pts = np.random.default_rng(0).uniform(0.0, 1.0, size=(50, 3))
    p0 = cosine_mode(pts, 0.0, 1.0, k, C_V)
    p1 = cosine_mode(pts, 0.01, 1.0, k, C_V)
    assert p0.shape == (50,)
    ratio = p1 / p0
    assert np.allclose(ratio, np.exp(-cosine_mode_decay_rate(k, C_V) * 0.01))


def test_moving_interval_conservation_identity() -> None:
    """dPhi/dt (analytic moving-face + flux + source) matches a finite-difference of Phi(t).

    Manufactured field p*(x,t) = (a + b x + c x^2) exp(-w t) on a breathing interval L(t).
    """
    a, b, c, w = 30.0, 1.5, -0.4, 0.7
    l0, eps, omega = 6.0, 0.2, 0.9

    def p(x, t):
        return (a + b * x + c * x**2) * np.exp(-w * t)

    def px(x, t):
        return (b + 2.0 * c * x) * np.exp(-w * t)

    def length(t):
        return l0 * (1.0 + eps * np.sin(omega * t))

    def length_rate(t):
        return l0 * eps * omega * np.cos(omega * t)

    def source(x, t):
        # s_water = S p_t - mobility p_xx = S(-w p) - mobility(2c e^{-wt})
        return STORAGE_S * (-w * p(x, t)) - MOBILITY * (2.0 * c * np.exp(-w * t))

    t = 0.5
    rate = moving_interval_content_rate(
        p, px, source, length, length_rate, t, STORAGE_S, MOBILITY
    )
    # direct central difference of the content integral
    dt = 1e-5
    phi_p = moving_interval_content(p, length, t + dt, STORAGE_S)
    phi_m = moving_interval_content(p, length, t - dt, STORAGE_S)
    dphi_fd = (phi_p - phi_m) / (2.0 * dt)

    assert rate["total"] == pytest.approx(dphi_fd, rel=1e-6)
    # all three channels are individually non-trivial (the gate exercises each)
    assert abs(rate["moving_face"]) > 1e-6
    assert abs(rate["face_flux"]) > 1e-6
    assert abs(rate["source"]) > 1e-6

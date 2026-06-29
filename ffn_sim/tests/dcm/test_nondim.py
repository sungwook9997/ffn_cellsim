"""Unit checks for the shared DCM/FF nondimensionalization (``ffn_sim.dcm.nondim``).

Pure arithmetic on measured scales — guards the faceting-group bridge against regression so an
FF-measured γ keeps mapping onto the same DCM regime.
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.dcm.nondim import (
    CellScales, FACETING_BAND, describe_gamma, douezan_spreading, faceting_regime,
    gamma_from_tilde, gamma_tilde,
)


def test_reference_length_scale():
    """ℓ = V₀^(1/3) for the MCF7 reference radius ≈ 12.09 µm (PARAM_AUDIT value)."""
    s = CellScales()
    assert s.ell == pytest.approx(12.09e-6, rel=1e-3)


def test_simucell3d_K_reaches_faceting_band():
    """At the SimuCell3D bulk modulus, a ~mN/m γ lands in the faceting band."""
    s = CellScales()
    gt = gamma_tilde(1.0e-3, s.K_simucell3d, s.ell)
    assert FACETING_BAND[0] <= gt <= FACETING_BAND[1]
    assert "FACETING" in faceting_regime(gt)


def test_driver_K_cannot_facet_in_measured_range():
    """At the driver default K, no γ in the measured MCF7 range (≤10 mN/m) can facet."""
    s = CellScales()
    for g in (1e-4, 1e-3, 1e-2):
        assert gamma_tilde(g, s.K_vol, s.ell) < FACETING_BAND[0]


def test_gamma_from_tilde_roundtrips():
    s = CellScales()
    for gt_target in (0.02, 0.05, 0.10):
        g = gamma_from_tilde(gt_target, s.K_simucell3d, s.ell)
        assert gamma_tilde(g, s.K_simucell3d, s.ell) == pytest.approx(gt_target, rel=1e-9)


def test_regime_boundaries():
    assert faceting_regime(0.01).startswith("rounded")
    assert "FACETING" in faceting_regime(0.05)
    assert faceting_regime(0.5).startswith("cortex-crushed")


def test_douezan_sign():
    """w_cs/(2γ) − 1: small γ ⇒ wetting (s>0); large γ ⇒ non-wetting (s<0)."""
    assert douezan_spreading(2.85e-3, 1e-3) > 0
    assert douezan_spreading(2.85e-3, 1e-2) < 0


def test_describe_gamma_keys_and_finiteness():
    s = CellScales()
    d = describe_gamma(1.0e-3, s, K=s.K_simucell3d)
    assert d["in_faceting_band"] is True
    for k in ("gamma_tilde", "elastocapillary_length_um", "turgor_capillary_ratio",
              "douezan_s", "viscous_capillary_time_s"):
        assert np.isfinite(d[k])

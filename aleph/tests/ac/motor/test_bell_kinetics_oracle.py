"""Self-test of the Bell slip + Pereverzev catch-slip kinetics oracle (pure NumPy — no Warp/CUDA)."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.motor.bell_kinetics_analytic import (
    attach_probability,
    bell_f0_from_x_beta,
    bell_off_rate,
    catch_slip_off_rate,
    catch_slip_peak_force,
    detach_probability,
    engaged_fraction_catch_slip,
    engaged_fraction_steady,
)

# illustrative (provisional/GAP) params — SHAPE gates are magnitude-independent (params_i0b3.yaml)
K_OFF0 = 0.35     # 1/s  (Stam-Hocky/Tam, provisional)
X_BETA = 0.6e-3   # um   (Veigel 2002, provisional)
K_ON = 50.0       # 1/s  (config-chosen, GAP)

KBT = 4.28e-3     # pN*um (37 C engine value)
# Illustrative NM2B-like catch-slip params (provisional/GAP — PI; magnitudes constrained by Kovacs 2007
# 5x/12x ADP-release slowdown + NM2B unloaded duty 0.2-0.3, NOT a single-molecule fit). Shape gates below
# are magnitude-independent; these values only exercise the biphasic form and the Kovacs load-raises-duty sign.
CS_K_CATCH0 = 3.0     # 1/s  catch-pathway (load-slowed ADP release) prefactor
CS_X_CATCH = 1.0e-3   # um   catch bond length
CS_K_SLIP0 = 0.5      # 1/s  slip-pathway (forced unbinding) prefactor
CS_X_SLIP = 0.6e-3    # um   slip (Bell) bond length
CS_K_ON = 1.2         # 1/s  attach rate giving unloaded duty ~0.25 (NM2B band)


def test_off_rate_zero_force_and_slip_monotonicity() -> None:
    """p_off(0) = k_off0 and p_off strictly increases with |load| (pure slip bond)."""
    f0 = bell_f0_from_x_beta(X_BETA)
    assert float(bell_off_rate(0.0, K_OFF0, f0)) == pytest.approx(K_OFF0, rel=1e-12)
    f = np.linspace(0.0, 30.0, 100)
    p = bell_off_rate(f, K_OFF0, f0)
    assert np.all(np.diff(p) > 0.0)                          # slip: rises with load
    # symmetric in the sign of the load (magnitude only)
    assert np.allclose(bell_off_rate(f, K_OFF0, f0), bell_off_rate(-f, K_OFF0, f0))


def test_bell_f0_from_x_beta_matches_ff_value() -> None:
    """f0 = kBT/x_beta reproduces the ff/hand_kmc NMIIA value ~7.1 pN at x_beta = 0.6 nm."""
    f0 = bell_f0_from_x_beta(X_BETA)
    assert f0 == pytest.approx(7.13, abs=0.1)


def test_probabilities_bounded_and_poisson_limits() -> None:
    """Per-tick attach/detach probabilities in [0,1); ~ rate*tau for small tau; -> 1 as tau -> inf."""
    assert 0.0 <= attach_probability(1e-4, K_ON) < 1.0
    assert attach_probability(1e-6, K_ON) == pytest.approx(1e-6 * K_ON, rel=1e-3)  # small-tau linear
    assert attach_probability(1e6, K_ON) == pytest.approx(1.0, abs=1e-9)           # saturates
    p = detach_probability(1e-4, np.array([0.35, 1.0, 5.0]))
    assert np.all((p >= 0.0) & (p < 1.0))


def test_engaged_fraction_bounds_and_unloaded_value() -> None:
    """phi_b in (0,1]; phi_b(0) = k_on/(k_on+k_off0); the unloaded engaged fraction."""
    f0 = bell_f0_from_x_beta(X_BETA)
    phi0 = float(engaged_fraction_steady(K_ON, 0.0, K_OFF0, f0))
    assert phi0 == pytest.approx(K_ON / (K_ON + K_OFF0), rel=1e-12)
    assert 0.0 < phi0 <= 1.0


def test_engaged_fraction_self_limiting_under_load() -> None:
    """phi_b strictly decreases with load (heads shed under force — self-limiting)."""
    f0 = bell_f0_from_x_beta(X_BETA)
    f = np.linspace(0.0, 30.0, 100)
    phi = engaged_fraction_steady(K_ON, f, K_OFF0, f0)
    assert np.all(np.diff(phi) < 0.0)
    assert np.all((phi > 0.0) & (phi <= 1.0))


def test_engaged_fraction_binding_limits() -> None:
    """k_on -> inf => phi_b -> 1 (always bound); k_on -> 0 => phi_b -> 0 (never bound)."""
    f0 = bell_f0_from_x_beta(X_BETA)
    assert float(engaged_fraction_steady(1e9, 5.0, K_OFF0, f0)) == pytest.approx(1.0, abs=1e-6)
    assert float(engaged_fraction_steady(0.0, 5.0, K_OFF0, f0)) == pytest.approx(0.0, abs=1e-15)


# ── Pereverzev catch-slip (NMII head-actin; Kovacs 2007 correction) ────────────────────────────────
def _cs(force):
    return catch_slip_off_rate(force, CS_K_CATCH0, CS_X_CATCH, CS_K_SLIP0, CS_X_SLIP, KBT)


def test_catch_slip_biphasic_with_minimum_at_peak_force() -> None:
    """Catch-slip off-rate FALLS with load (catch) then RISES (slip); the minimum is at the analytic F*."""
    fstar = catch_slip_peak_force(CS_K_CATCH0, CS_X_CATCH, CS_K_SLIP0, CS_X_SLIP, KBT)
    assert fstar > 0.0
    # off-rate drops below its unloaded value on the catch branch, then exceeds it on the slip branch
    assert float(_cs(fstar)) < float(_cs(0.0))
    assert float(_cs(3.0 * fstar)) > float(_cs(0.0))
    # F* is the global minimum of a dense sweep (argmin coincides with the closed form within the grid)
    f = np.linspace(0.0, 4.0 * fstar, 4001)
    p = _cs(f)
    assert f[int(np.argmin(p))] == pytest.approx(fstar, abs=(f[1] - f[0]) * 2)
    # biphasic shape: strictly decreasing before F*, strictly increasing after
    lo = p[f < fstar]
    hi = p[f > fstar]
    assert np.all(np.diff(lo) < 0.0)
    assert np.all(np.diff(hi) > 0.0)


def test_catch_slip_symmetric_in_load_sign() -> None:
    """Only the load magnitude matters (|f|), matching bell_off_rate and the device wp.abs(f)."""
    f = np.linspace(0.0, 20.0, 50)
    assert np.allclose(_cs(f), _cs(-f))


def test_catch_slip_reduces_to_bell_slip_when_catch_vanishes() -> None:
    """k_catch0 -> 0 leaves a pure slip term k_slip0*exp(|f|*x_slip/kT) == bell_off_rate(f, k_slip0, kT/x_slip)."""
    f = np.linspace(0.0, 25.0, 60)
    slip_only = catch_slip_off_rate(f, 1e-12, CS_X_CATCH, CS_K_SLIP0, CS_X_SLIP, KBT)
    bell = bell_off_rate(f, CS_K_SLIP0, KBT / CS_X_SLIP)
    assert np.allclose(slip_only, bell, rtol=1e-9, atol=1e-12)


def test_catch_slip_matches_pereverzev_single_source_of_truth() -> None:
    """The oracle delegates to ff.hand_kmc.pereverzev_off_rate (project SoT) bit-for-formula on |f|."""
    from aleph.laws.hand_kmc import pereverzev_off_rate

    f = np.linspace(0.0, 15.0, 40)
    sot = pereverzev_off_rate(np.abs(f), CS_K_CATCH0, CS_X_CATCH, CS_K_SLIP0, CS_X_SLIP, KBT)
    assert np.allclose(_cs(f), sot, rtol=1e-12, atol=0.0)


def test_catch_slip_peak_force_requires_catch_dominant_low_load() -> None:
    """No interior catch peak when k_catch0*x_catch <= k_slip0*x_slip (slip-only) -> ValueError."""
    with pytest.raises(ValueError):
        catch_slip_peak_force(0.01, CS_X_CATCH, CS_K_SLIP0, CS_X_SLIP, KBT)  # catch too weak


def test_engaged_fraction_catch_slip_rises_under_load_then_falls() -> None:
    """THE Kovacs-2007 FIX: duty RISES with resistive load up to F*, then falls — opposite of pure slip.

    phi_b(F*) > phi_b(0) (load raises the bound duty for tension maintenance); the peak of phi_b coincides
    with the off-rate minimum F*; and unlike engaged_fraction_steady (slip) it is NOT monotone-decreasing.
    """
    fstar = catch_slip_peak_force(CS_K_CATCH0, CS_X_CATCH, CS_K_SLIP0, CS_X_SLIP, KBT)
    phi0 = float(engaged_fraction_catch_slip(CS_K_ON, 0.0, CS_K_CATCH0, CS_X_CATCH, CS_K_SLIP0, CS_X_SLIP, KBT))
    phistar = float(
        engaged_fraction_catch_slip(CS_K_ON, fstar, CS_K_CATCH0, CS_X_CATCH, CS_K_SLIP0, CS_X_SLIP, KBT)
    )
    # load raises duty (catch-bond tension maintenance)
    assert phistar > phi0
    # phi_b peaks exactly at the off-rate minimum F*
    f = np.linspace(0.0, 4.0 * fstar, 4001)
    phi = engaged_fraction_catch_slip(CS_K_ON, f, CS_K_CATCH0, CS_X_CATCH, CS_K_SLIP0, CS_X_SLIP, KBT)
    assert f[int(np.argmax(phi))] == pytest.approx(fstar, abs=(f[1] - f[0]) * 2)
    # explicitly NON-monotone under load (would be strictly decreasing for a pure slip bond)
    assert not np.all(np.diff(phi) < 0.0)
    assert np.all((phi > 0.0) & (phi <= 1.0))


def test_unloaded_nm2b_duty_in_band() -> None:
    """Unloaded engaged fraction phi_b(0) = k_on/(k_on + k_catch0 + k_slip0) lands in the NM2B 0.2-0.3 band."""
    phi0 = float(engaged_fraction_catch_slip(CS_K_ON, 0.0, CS_K_CATCH0, CS_X_CATCH, CS_K_SLIP0, CS_X_SLIP, KBT))
    assert phi0 == pytest.approx(CS_K_ON / (CS_K_ON + CS_K_CATCH0 + CS_K_SLIP0), rel=1e-12)
    assert 0.2 <= phi0 <= 0.3


def test_catch_slip_off_rate_rejects_nonpositive_params() -> None:
    """Every catch-slip constant must be finite and positive (no silent slip-only fallback)."""
    for bad in ("k_catch0", "x_catch", "k_slip0", "x_slip", "kT"):
        kwargs = dict(k_catch0=CS_K_CATCH0, x_catch=CS_X_CATCH, k_slip0=CS_K_SLIP0, x_slip=CS_X_SLIP, kT=KBT)
        kwargs[bad] = 0.0
        with pytest.raises(ValueError):
            catch_slip_off_rate(1.0, **kwargs)

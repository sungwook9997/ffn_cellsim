"""Active-gel polarization (ff/polarization_activegel) — piece 2/5 of active movement.

Validates the spectral solver against its ANALYTIC linear dispersion λ(k) (oracle-is-crosscheck rule), the
marginal-stability threshold ζ_c, and the nonlinear spontaneous single-cap symmetry-breaking (Bois-Jülicher-Grill
2011 actomyosin contractile-flow instability). Constants (D, k_off, ℓ) are lit-anchored (Bois 2011 / Mayer 2010);
ζ is the swept control variable; γ is a normalization (absolute flow speed is NOT asserted — PI-gated).
"""

import numpy as np
import pytest

from ffn_sim.ff.polarization_activegel import (
    dispersion,
    make_grid,
    measure_growth_rate,
    resolve_activegel,
    step,
    zeta_critical,
    ELL_HYDRO_UM,
)


def test_resolve_and_dispersion_form():
    """Ring length L=2πR; λ(k) has the Bois form (active − diffusion − turnover)."""
    p = resolve_activegel(R_um=7.5, zeta=100.0)
    assert p.L == pytest.approx(2 * np.pi * 7.5)
    assert p.ell == ELL_HYDRO_UM
    # at k→0 the active term ~k² and diffusion ~k² vanish → λ(0) = -k_off (pure decay)
    assert dispersion(1e-6, p) == pytest.approx(-p.k_off, abs=1e-6)
    with pytest.raises(ValueError):
        resolve_activegel(R_um=-1.0)


def test_solver_matches_analytic_dispersion():
    """The spectral solver's numerical growth rate per Fourier mode == analytic λ(k) (the ground truth)."""
    p0 = resolve_activegel(R_um=7.5, zeta=1.0)
    zc1 = zeta_critical(p0, 2 * np.pi / p0.L)
    p = resolve_activegel(R_um=7.5, zeta=3.0 * zc1)      # several modes unstable
    for mode in (1, 2, 3, 6, 8):
        lam_num, lam_an = measure_growth_rate(p, mode, N=256, eps=1e-5, n_steps=150)
        assert lam_num == pytest.approx(lam_an, rel=0.05, abs=1e-3)
    # band structure: low modes grow, high modes decay
    assert dispersion(2 * np.pi / p.L, p) > 0 and dispersion(8 * 2 * np.pi / p.L, p) < 0


def test_threshold_brackets_zeta_c():
    """Below ζ_c the homogeneous state is stable (λ<0); above, it is unstable (λ>0) — onset at ζ_c."""
    p0 = resolve_activegel(R_um=7.5, zeta=1.0)
    zc1 = zeta_critical(p0, 2 * np.pi / p0.L)
    lo, _ = measure_growth_rate(resolve_activegel(R_um=7.5, zeta=0.9 * zc1), 1, n_steps=150)
    hi, _ = measure_growth_rate(resolve_activegel(R_um=7.5, zeta=1.1 * zc1), 1, n_steps=150)
    assert lo < 0 < hi


def test_nonlinear_single_cap_polarization():
    """Above threshold, random IC → a SINGLE high-myosin cap (spontaneous front-rear symmetry-breaking)."""
    p0 = resolve_activegel(R_um=7.5, zeta=1.0)
    zc1 = zeta_critical(p0, 2 * np.pi / p0.L)
    p = resolve_activegel(R_um=7.5, zeta=4.0 * zc1)
    x, dx, kg = make_grid(p, N=256)
    rng = np.random.default_rng(0)
    c = p.c0 + 1e-3 * rng.standard_normal(x.size)
    dt = 0.05 * min(dx ** 2 / p.D, 1.0 / p.k_off)
    for _ in range(60000):
        c, _ = step(c, p, dx, dt, kg)
    assert np.isfinite(c).all()
    dom = int(np.argmax(np.abs(np.fft.rfft(c - c.mean()))[1:]) + 1)
    contrast = (c.max() - c.min()) / c.mean()
    assert dom == 1 and contrast > 0.1                    # one cap, strongly polarized

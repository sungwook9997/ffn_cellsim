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
    evolve,
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


def test_single_cap_robust_near_threshold():
    """In the POLARIZATION regime — ζ just above ζ_c1 where ONLY the cell-perimeter mode k1 is unstable — a
    single myosin cap forms for EVERY random seed (seed-independent spontaneous front-rear symmetry-breaking).

    (This is the honest robust regime. NOT ζ≫ζ_c: far above threshold shorter modes are unstable too and the
    system forms transient MULTI-cap states whose coarsening to one cap is slow + seed-dependent — see
    test_multicap_above_threshold_is_not_robust. The earlier single-seed dom==1 at ζ=4·ζ_c1 was a
    seed-cherry-pick, corrected 2026-07-02 audit.)"""
    p0 = resolve_activegel(R_um=7.5, zeta=1.0)
    k1 = 2 * np.pi / p0.L
    zc1, zc2 = zeta_critical(p0, k1), zeta_critical(p0, 2 * k1)
    zeta = 1.2 * zc1
    assert zc1 < zeta < zc2                               # ONLY mode 1 linearly unstable
    p = resolve_activegel(R_um=7.5, zeta=zeta)
    assert dispersion(k1, p) > 0 and dispersion(2 * k1, p) < 0
    for seed in range(4):
        _, dom, contrast = evolve(p, seed=seed, N=256, n_steps=90000)
        assert dom == 1 and contrast > 0.0               # single cap, polarized, EVERY seed


def test_multicap_above_threshold_is_not_robust():
    """Honesty guard: far above threshold (ζ=4·ζ_c1, several modes unstable) the single-cap outcome is NOT
    seed-robust — some seeds settle on 2–3 caps. We must therefore NOT claim a robust single cap there."""
    p0 = resolve_activegel(R_um=7.5, zeta=1.0)
    zc1 = zeta_critical(p0, 2 * np.pi / p0.L)
    p = resolve_activegel(R_um=7.5, zeta=4.0 * zc1)
    doms = [evolve(p, seed=s, N=256, n_steps=60000)[1] for s in range(6)]
    assert any(d != 1 for d in doms)                     # multi-cap for at least one seed (seed-dependent)

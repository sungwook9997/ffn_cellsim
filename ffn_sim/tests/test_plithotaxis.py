"""Sanity-Gate tests for the §E SPP-plithotaxis polarity field (plithotaxis.py).

Pins the anchored constants (Smeets-2016 MCF10A: D_r, f_cil, ψ≈1), the force contract
(in-plane / z≡0 / basal-gated / f_active=0 ⇒ 0), the persistence law (⟨cosΔθ⟩=e^{−D_r t}),
the CIL free-edge outward repolarisation, and the pure-noise null (no coherent front without
CIL — the §D "radial does nothing" recovered as a limit).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.spheroid.params import resolve_layer2
from ffn_sim.spheroid.plithotaxis import (
    PolarizationField,
    ResolvedPlithotaxis,
    resolve_plithotaxis,
)

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


def _resolved():
    return resolve_layer2(yaml.safe_load(_CFG.read_text()))


# --- anchored constants -----------------------------------------------------------------------

def test_anchored_constants_si_and_psi():
    """D_r=0.05/min, f_cil=0.1/min in SI; τ=1/D_r=20 min; ψ=f_cil/(2D_r)=1 (Smeets MCF10A)."""
    p = resolve_plithotaxis(_resolved())
    assert p.D_r == pytest.approx(0.05 / 60.0, rel=1e-12)
    assert p.f_cil == pytest.approx(0.1 / 60.0, rel=1e-12)
    assert p.tau == pytest.approx(20.0 * 60.0, rel=1e-9)        # 20 min in seconds
    assert p.psi == pytest.approx(1.0, rel=1e-9)
    assert p.J == 0.0                                          # Vicsek OFF by default
    assert p.interaction_radius == pytest.approx(1.3 * _resolved().morse_r0, rel=1e-12)


def test_vicsek_J_negative_rejected():
    with pytest.raises(ValueError):
        resolve_plithotaxis(_resolved(), J=-1.0)


# --- force contract ---------------------------------------------------------------------------

def _field(n, plith=None, seed=0):
    return PolarizationField(n, plith or resolve_plithotaxis(_resolved()), seed=seed)


def test_zero_factive_is_zero_force():
    f = _field(10)
    pos = np.random.default_rng(0).normal(size=(10, 3))
    act = np.ones(10, dtype=bool)
    out = f.forces(pos, act, f_active=0.0, z_substrate=0.0, basal_band=1.0)
    assert np.allclose(out, 0.0)


def test_force_in_plane_zero_z():
    f = _field(10)
    pos = np.random.default_rng(1).normal(size=(10, 3)) * 0.1   # near substrate plane
    act = np.ones(10, dtype=bool)
    out = f.forces(pos, act, f_active=9.4e-9, z_substrate=0.0, basal_band=1.0)
    assert np.allclose(out[:, 2], 0.0)                          # crawl tangent to dish


def test_force_magnitude_equals_factive_at_substrate():
    """A cell exactly on the substrate plane (w_basal=1) feels |F|=f_active along p̂."""
    f = _field(1)
    pos = np.array([[2.0, 1.0, 0.0]])                           # dz=0 ⇒ w_basal=1
    out = f.forces(pos, np.array([True]), f_active=5.0e-9, z_substrate=0.0, basal_band=1.0)
    assert np.linalg.norm(out[0]) == pytest.approx(5.0e-9, rel=1e-9)


def test_apical_cell_no_force():
    f = _field(2)
    pos = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 10.0]])         # one basal, one apical (z=10≫band)
    out = f.forces(pos, np.ones(2, dtype=bool), f_active=1e-9, z_substrate=0.0, basal_band=1.0)
    assert np.linalg.norm(out[1]) == pytest.approx(0.0, abs=1e-18)
    assert np.linalg.norm(out[0]) > 0.0


def test_empty_active_returns_empty():
    f = _field(5)
    out = f.forces(np.zeros((5, 3)), np.zeros(5, dtype=bool), f_active=1e-9, z_substrate=0.0, basal_band=1.0)
    assert out.shape == (0, 3)


# --- persistence (rotational diffusion): ⟨cosΔθ⟩ = exp(−D_r t) -------------------------------

def test_persistence_autocorrelation_matches_Dr():
    """Isolated cells (no neighbours ⇒ no CIL): polarity autocorrelation decays as exp(−D_r t)."""
    D_r = 0.01                                                  # s⁻¹ (pick clean numbers)
    plith = ResolvedPlithotaxis(
        D_r=D_r, f_cil=0.0, J=0.0, interaction_radius=0.1, n_bulk=6, tau=1.0 / D_r, psi=0.0
    )
    n = 4000
    f = PolarizationField(n, plith, seed=7)
    theta0 = f.theta.copy()
    # cells placed far apart (spacing ≫ interaction_radius) ⇒ no neighbours ⇒ pure diffusion
    pos = np.zeros((n, 3))
    pos[:, 0] = np.arange(n) * 10.0
    act = np.ones(n, dtype=bool)
    dt = 5.0
    T = 0.0
    for _ in range(20):
        f.update(pos, act, dt, z_substrate=0.0, basal_band=1.0)
        T += dt
    corr = np.cos(f.theta - theta0).mean()
    assert corr == pytest.approx(np.exp(-D_r * T), abs=0.05)    # D_r·T = 1.0 ⇒ ~0.368


# --- CIL free-edge outward repolarisation ----------------------------------------------------

def test_cil_repolarises_edge_cell_outward():
    """A free-edge cell initially pointing INWARD rotates toward the outward (free-space) dir."""
    # a connected line of cells along +x (spacing 1.0 < interaction_radius); the rightmost cell's
    # only neighbour is inward (−x), so its free-space direction θ^free = +x (outward).
    pos = np.array([[float(i), 0.0, 0.0] for i in range(5)])     # x = 0,1,2,3,4
    n = pos.shape[0]
    plith = ResolvedPlithotaxis(
        D_r=0.0, f_cil=0.05, J=0.0, interaction_radius=1.5, n_bulk=6, tau=np.inf, psi=0.0
    )
    f = PolarizationField(n, plith, seed=0)
    f._theta[:] = 2.5            # inward-ish (cos<0), NOT exactly antipodal (avoids the sin=0 saddle)
    act = np.ones(n, dtype=bool)
    for _ in range(60):
        f.update(pos, act, 5.0, z_substrate=0.0, basal_band=5.0)
    rim_angle = f.theta[4]       # rightmost cell
    # rim should now point outward (+x): cos(angle) > 0, rotated from 2.5 rad toward 0
    assert np.cos(rim_angle) > 0.5


def test_bulk_cell_no_cil_drive():
    """A fully-surrounded cell (w_edge≈0) gets ~no CIL torque (only diffusion)."""
    # dense ring of neighbours around a central cell ⇒ central has ≥ n_bulk neighbours ⇒ w_edge=0
    ang = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    ring = np.column_stack([np.cos(ang), np.sin(ang), np.zeros(8)])
    pos = np.vstack([[[0.0, 0.0, 0.0]], ring])
    n = pos.shape[0]
    plith = ResolvedPlithotaxis(
        D_r=0.0, f_cil=0.05, J=0.0, interaction_radius=1.5, n_bulk=6, tau=np.inf, psi=0.0
    )
    f = PolarizationField(n, plith, seed=0)
    f._theta[:] = 0.5
    act = np.ones(n, dtype=bool)
    for _ in range(20):
        f.update(pos, act, 5.0, z_substrate=0.0, basal_band=5.0)
    # central cell barely moved (no noise, w_edge≈0 ⇒ no CIL): angle stays near 0.5
    assert f.theta[0] == pytest.approx(0.5, abs=1e-6)


# --- pure-noise null: no coherent front without CIL (the §D "radial does nothing" limit) ------

def test_pure_noise_no_coherent_polarity():
    """f_cil=0 ⇒ random polarities ⇒ ensemble-mean polarity ~0 (no net front); CIL breaks it."""
    rng = np.random.default_rng(5)
    pos = rng.normal(scale=2.0, size=(200, 3)) * np.array([1.0, 1.0, 0.0])
    n = pos.shape[0]
    act = np.ones(n, dtype=bool)
    # NO CIL, only diffusion
    plith0 = ResolvedPlithotaxis(D_r=0.02, f_cil=0.0, J=0.0, interaction_radius=1.5, n_bulk=6, tau=50.0, psi=0.0)
    f0 = PolarizationField(n, plith0, seed=1)
    for _ in range(30):
        f0.update(pos, act, 5.0, z_substrate=0.0, basal_band=5.0)
    p0 = f0.polarization(act)
    mean_mag0 = np.linalg.norm(p0.mean(axis=0))
    # WITH CIL: edge cells coherently point radially-outward ⇒ outward-radial order > the null
    plithC = ResolvedPlithotaxis(D_r=0.02, f_cil=0.1, J=0.0, interaction_radius=1.5, n_bulk=6, tau=50.0, psi=2.5)
    fC = PolarizationField(n, plithC, seed=1)
    for _ in range(30):
        fC.update(pos, act, 5.0, z_substrate=0.0, basal_band=5.0)
    pC = fC.polarization(act)
    rhat = pos[:, :2] / np.linalg.norm(pos[:, :2], axis=1, keepdims=True).clip(1e-12)
    outward_C = (pC * rhat).sum(axis=1).mean()       # mean outward-radial projection with CIL
    outward_0 = (p0 * rhat).sum(axis=1).mean()        # without CIL
    assert mean_mag0 < 0.3                            # no coherent global polarity from noise alone
    assert outward_C > 0.08                           # CIL induces a clear net OUTWARD polarity
    assert outward_C > 2.0 * abs(outward_0)           # ... far above the no-CIL isotropic null


# --- numerical guard: large dt sub-steps without blow-up --------------------------------------

def test_large_dt_substeps_stay_bounded():
    plith = resolve_plithotaxis(_resolved())
    f = PolarizationField(50, plith, seed=0)
    pos = np.random.default_rng(0).normal(scale=1e-5, size=(50, 3))
    pos[:, 2] = np.abs(pos[:, 2])
    act = np.ones(50, dtype=bool)
    # a deliberately huge epoch dt (would give Δθ≫1 in a single step) must sub-step, not blow up
    f.update(pos, act, 1e5, z_substrate=0.0, basal_band=1e-5)
    assert np.all(np.isfinite(f.theta))
    assert np.all(np.abs(f.theta) <= np.pi + 1e-9)   # wrapped into (−π, π]

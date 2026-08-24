r"""Compliant elastic substrate (ff/substrate) — E→k_sub bridge + Bangasser-Odde series compliance.

FF-native validation of the compliant-substrate machinery that carries the PI's Bare/Pre/Lam4 PAA-gel axis.
Everything here is an ANALYTIC ground-truth check (oracle-is-crosscheck), not a fit:

  1. FLAT-PUNCH bridge — ``resolve_substrate`` maps the substrate Young's modulus E to the Winkler anchor
     stiffness by the Boussinesq/Sneddon rigid-circular-punch relation ``k_sub = 2·E·a/(1−ν²)`` exactly, and
     k_sub ∝ E exactly (the E-dependence that drives the biphasic/durotaxis physics is independent of the
     surfaced patch radius ``a``).
  2. SERIES COMPLIANCE — the equilibrium anchor kernel places the movable anchor so the actin end feels the
     clutch↔substrate SERIES stiffness ``k_series = k_int·k_sub/(k_int+k_sub)`` (Bangasser-Odde). Derived
     from first principles (not from the kernel) and checked against the kernel's anchor placement.
  3. RIGID / SOFT LIMITS — k_sub→∞ ⇒ anchor→rest (the fixed-pin path recovered byte-for-byte, the off-path
     invariance); k_sub→0 ⇒ anchor→loaded actin (no traction transmitted); unbound ⇒ anchor==rest.
  4. CFL + input validation guards.

Magnitude note: the absolute ``k_sub`` carries the SURFACED patch radius ``a`` (KB-blocked, same missing-datum
family as ρ_L); only the E-DEPENDENCE and the series shape are asserted as gates — no magnitude match is claimed.
"""

import numpy as np
import pytest
import warp as wp

from aleph.laws.substrate import (
    A_ADHESION_UM_DEFAULT,
    E_SUB_DEFAULT_PA,
    NU_SUB_DEFAULT,
    SubstrateParams,
    clutch_anchor_reaction_kernel,
    resolve_substrate,
    substrate_anchor_equilibrium_kernel,
    substrate_cfl_dt,
)


# ----------------------------------------------------------------------------- E → k_sub flat-punch bridge
def test_flat_punch_k_sub_is_exact():
    """k_sub = 2·E·a/(1−ν²) (Boussinesq/Sneddon flat-punch) — exact, for a spread of (E, ν, a)."""
    for E, nu, a in [(5.0e3, 0.45, 0.05), (1.0e2, 0.0, 0.1), (5.0e4, 0.49, 0.2), (2.0e3, 0.3, 0.03)]:
        s = resolve_substrate(E_pa=E, nu=nu, a_um=a)
        assert s.k_sub == pytest.approx(2.0 * E * a / (1.0 - nu**2), rel=1e-12)
        assert isinstance(s, SubstrateParams) and s.E_pn_um2 == E and s.a_um == a


def test_k_sub_is_linear_in_E():
    """The stiffness axis (biphasic / durotaxis) rides k_sub ∝ E exactly — independent of the surfaced ``a``."""
    a, nu = 0.05, 0.45
    k1 = resolve_substrate(E_pa=1.0e3, nu=nu, a_um=a).k_sub
    k2 = resolve_substrate(E_pa=2.0e3, nu=nu, a_um=a).k_sub
    k10 = resolve_substrate(E_pa=1.0e4, nu=nu, a_um=a).k_sub
    assert k2 == pytest.approx(2.0 * k1, rel=1e-12)
    assert k10 == pytest.approx(10.0 * k1, rel=1e-12)
    # ``a`` only rescales the magnitude, never the E-slope
    ka = resolve_substrate(E_pa=1.0e3, nu=nu, a_um=2 * a).k_sub
    assert ka == pytest.approx(2.0 * k1, rel=1e-12)


def test_defaults_are_kb_anchored():
    """The Phase-1 defaults are the KB-1.5 PAA-gel values, not tuned constants."""
    assert E_SUB_DEFAULT_PA == 5.0e3 and NU_SUB_DEFAULT == 0.45
    assert A_ADHESION_UM_DEFAULT == 0.05                       # nascent-adhesion ~50 nm (surfaced, order-of-magnitude)
    s = resolve_substrate()
    assert s.k_sub == pytest.approx(2.0 * 5.0e3 * 0.05 / (1.0 - 0.45**2), rel=1e-12)


def test_input_validation_raises():
    for bad in [dict(E_pa=0.0), dict(E_pa=-1.0), dict(nu=0.5), dict(nu=0.6), dict(a_um=0.0), dict(a_um=-0.1)]:
        with pytest.raises(ValueError):
            resolve_substrate(**bad)


# ----------------------------------------------------------------------------- series compliance (Bangasser-Odde)
def _run_equilibrium(p_actin, p_rest, bound, k_int, k_sub):
    """Launch substrate_anchor_equilibrium_kernel on one anchor and return the resulting anchor position."""
    d = "cpu"
    pos = wp.array(np.array([p_actin], float), dtype=wp.vec3d, device=d)
    ac = wp.array([0], dtype=wp.int32, device=d)
    anchor = wp.zeros(1, dtype=wp.vec3d, device=d)
    rest = wp.array(np.array([p_rest], float), dtype=wp.vec3d, device=d)
    bd = wp.array([int(bound)], dtype=wp.int32, device=d)
    wp.launch(substrate_anchor_equilibrium_kernel, dim=1,
              inputs=[pos, ac, anchor, rest, bd, wp.float64(k_int), wp.float64(k_sub)], device=d)
    return anchor.numpy()[0]


def test_equilibrium_anchor_is_series_balance():
    """anchor = (k_int·p_actin + k_sub·p_rest)/(k_int+k_sub) — the clutch↔substrate force balance (kernel parity)."""
    p_actin = np.array([0.30, 0.0, 0.0]); p_rest = np.array([0.0, 0.0, 0.0])
    k_int, k_sub = 1000.0, 400.0
    got = _run_equilibrium(p_actin, p_rest, 1, k_int, k_sub)
    want = (k_int * p_actin + k_sub * p_rest) / (k_int + k_sub)
    assert np.allclose(got, want, atol=1e-12)


def test_actin_feels_series_stiffness():
    """FIRST-PRINCIPLES: with the anchor at series balance, the actin end feels a restoring force
    ``k_series·|Δ|`` with ``k_series = k_int·k_sub/(k_int+k_sub)`` (Bangasser-Odde compliant-clutch)."""
    p_actin = np.array([0.20, 0.0, 0.0]); p_rest = np.array([0.0, 0.0, 0.0])
    for k_int, k_sub in [(1000.0, 100.0), (1000.0, 1000.0), (500.0, 5000.0)]:
        anchor = _run_equilibrium(p_actin, p_rest, 1, k_int, k_sub)
        f_on_actin = k_int * (anchor - p_actin)               # linear clutch spring pulls actin toward anchor
        k_series = k_int * k_sub / (k_int + k_sub)
        assert np.linalg.norm(f_on_actin) == pytest.approx(k_series * np.linalg.norm(p_actin - p_rest), rel=1e-10)
        # and it is restoring (points from loaded actin back toward rest)
        assert np.dot(f_on_actin, p_rest - p_actin) > 0


def test_rigid_limit_recovers_fixed_pin():
    """k_sub→∞ ⇒ anchor→rest: the compliant path collapses onto the fixed-pin path (off-path invariance)."""
    p_actin = np.array([0.5, 0.1, -0.2]); p_rest = np.array([0.0, 0.0, 0.0])
    anchor = _run_equilibrium(p_actin, p_rest, 1, 1000.0, 1.0e18)
    assert np.allclose(anchor, p_rest, atol=1e-6)


def test_soft_limit_transmits_no_traction():
    """k_sub→0 ⇒ anchor→loaded actin ⇒ clutch length→rest ⇒ zero series stiffness (no force transmitted)."""
    p_actin = np.array([0.5, 0.1, -0.2]); p_rest = np.array([0.0, 0.0, 0.0])
    anchor = _run_equilibrium(p_actin, p_rest, 1, 1000.0, 1.0e-18)
    assert np.allclose(anchor, p_actin, atol=1e-6)


def test_unbound_anchor_relaxes_to_rest():
    """A detached clutch (bound=0) parks its anchor at the dish rest point regardless of actin load."""
    p_actin = np.array([0.7, -0.3, 0.4]); p_rest = np.array([0.1, 0.1, 0.1])
    anchor = _run_equilibrium(p_actin, p_rest, 0, 1000.0, 400.0)
    assert np.allclose(anchor, p_rest, atol=1e-12)


def test_reaction_is_newton_pair_of_clutch_spring():
    """clutch_anchor_reaction_kernel returns the clutch-spring reaction ON the anchor = k_int·(L−rest)·d̂,
    the Newton pair of the force the clutch exerts on the actin. Unbound → 0."""
    d = "cpu"
    L, rest = 0.20, 0.05
    pos = wp.array(np.array([[L, 0.0, 0.0], [L, 0.0, 0.0]]), dtype=wp.vec3d, device=d)
    ac = wp.array([0, 1], dtype=wp.int32, device=d)
    anchor = wp.array(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]), dtype=wp.vec3d, device=d)
    bd = wp.array([1, 0], dtype=wp.int32, device=d)            # anchor 0 bound, 1 detached
    out = wp.zeros(2, dtype=wp.vec3d, device=d)
    k_int = 1000.0
    wp.launch(clutch_anchor_reaction_kernel, dim=2,
              inputs=[pos, ac, anchor, bd, wp.float64(k_int), wp.float64(rest), out], device=d)
    r = out.numpy()
    assert r[0] == pytest.approx(np.array([k_int * (L - rest), 0.0, 0.0]))   # magnitude k_int·(L−rest), toward actin
    assert np.allclose(r[1], 0.0)                              # detached anchor feels nothing


# ----------------------------------------------------------------------------- CFL guard
def test_cfl_dt_matches_definition():
    """Explicit anchor CFL dt = safety·γ/k_sub, and larger k_sub ⇒ smaller stable dt."""
    assert substrate_cfl_dt(500.0, 1.0e-3, safety=0.1) == pytest.approx(0.1 * 1.0e-3 / 500.0, rel=1e-12)
    assert substrate_cfl_dt(5000.0, 1.0e-3) < substrate_cfl_dt(500.0, 1.0e-3)

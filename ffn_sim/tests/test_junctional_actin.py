"""Sanity-Gate tests for ``ffn_sim.junction.junctional_actin`` (STUB).

Resolver-level + OFF-identity + sign-sense + no-γ-contamination checks. These
run cheaply (no full Cell, no HOOMD sim). The STUB's enabled-but-un-anchored
physics is asserted to raise ``NotImplementedError`` (Hard Rule 3: do not
invent constants to make it run).

Maps to the module Sanity Gate sections:
  §2 boundary/OFF-identity, §1 dimensional/derived, §5 sign-sense (catch +
  tensile clutch), Hard-Rule-5 denylist prefix.
"""

from __future__ import annotations

import math
import re

import numpy as np
import pytest

from ffn_sim.junction import junctional_actin as ja


KT = 4.28e-21    # J  (KU thermal energy at ~37 °C, the project default)
DT = 1.0e-8      # s  (representative host dt)


# ---------------------------------------------------------------------------
# §2 OFF-identity — disabled resolver + no-op builder
# ---------------------------------------------------------------------------
def test_resolve_default_off_when_missing():
    """A config with no ``enabled`` key resolves DEFAULT-OFF (Hard Rule 2)."""
    p = ja.resolve_junctional_actin({}, kT=KT, dt=DT)
    assert p.enabled is False
    assert p.n_couplings_max == 0
    # Every un-anchored physics field is None (no invented constants).
    for fld in (
        "max_couple_dist", "anchor_r0", "k_anchor", "k_couple",
        "k_catch0", "x_catch", "k_slip0", "x_slip", "k_on",
    ):
        assert getattr(p, fld) is None
    assert p.is_anchored is False


def test_resolve_explicit_false():
    """``enabled: false`` resolves DEFAULT-OFF, even nested under junction."""
    p = ja.resolve_junctional_actin(
        {"junction": {"junctional_actin": {"enabled": False}}}, kT=KT, dt=DT
    )
    assert p.enabled is False
    assert p.n_couplings_max == 0


def test_builder_off_identity_returns_input_unchanged():
    """Disabled builder is a bit-for-bit no-op: returns its input snapshot."""
    p = ja.resolve_junctional_actin({"enabled": False}, kT=KT, dt=DT)
    sentinel = object()                 # any snapshot-like object
    out = ja.extend_snapshot_with_junctional_actin(sentinel, p)
    assert out is sentinel              # same object, unchanged (zero particles/bonds)


def test_builder_off_identity_adds_zero_particles_and_bonds():
    """Disabled builder adds zero particles + zero bonds to a tiny snapshot."""
    p = ja.resolve_junctional_actin({"enabled": False}, kT=KT, dt=DT)

    class _TinySnap:
        def __init__(self):
            self.n_particles = 3
            self.n_bonds = 0

    snap = _TinySnap()
    out = ja.extend_snapshot_with_junctional_actin(snap, p)
    assert out is snap
    assert out.n_particles == 3        # unchanged
    assert out.n_bonds == 0            # unchanged


# ---------------------------------------------------------------------------
# §1 dimensional / parameter sanity — resolver validation
# ---------------------------------------------------------------------------
def test_resolve_rejects_nonpositive_kT_dt():
    with pytest.raises(ValueError):
        ja.resolve_junctional_actin({}, kT=0.0, dt=DT)
    with pytest.raises(ValueError):
        ja.resolve_junctional_actin({}, kT=KT, dt=-1.0e-9)


def test_resolve_rejects_bad_bins_batch():
    with pytest.raises(ValueError):
        ja.resolve_junctional_actin({"n_bins": 0}, kT=KT, dt=DT)
    with pytest.raises(ValueError):
        ja.resolve_junctional_actin({"batch_steps": 0}, kT=KT, dt=DT)


def test_resolve_rejects_negative_n_cadherin():
    with pytest.raises(ValueError):
        ja.resolve_junctional_actin({}, kT=KT, dt=DT, n_cadherin=-5)


def test_resolve_enabled_plans_topology_but_keeps_physics_none():
    """enabled=True plans the coupling-head count but invents NO constants."""
    p = ja.resolve_junctional_actin(
        {"enabled": True}, kT=KT, dt=DT, n_cadherin=100
    )
    assert p.enabled is True
    assert p.n_couplings_max == 100    # one head per cadherin (topology bound)
    assert p.is_anchored is False      # constants still un-anchored (STUB)
    assert p.k_couple is None
    # batch_dt derived correctly (dimensional check §1).
    assert math.isclose(p.batch_dt, p.batch_steps * DT, rel_tol=1e-12)


def test_resolve_enabled_rejects_explicit_nonpositive_constant():
    """A PI-supplied physics constant is range-checked (no silent override)."""
    with pytest.raises(ValueError):
        ja.resolve_junctional_actin(
            {"enabled": True, "k_couple": -1.0}, kT=KT, dt=DT
        )
    with pytest.raises(ValueError):
        ja.resolve_junctional_actin(
            {"enabled": True, "max_couple_dist": 0.0}, kT=KT, dt=DT
        )


def test_bin_rest_lengths_dimensional():
    """Bin centers are strictly inside (0, max] and monotone (§1)."""
    n_bins = 8
    max_d = 60.0e-9
    r0 = ja.junc_actin_couple_bin_rest_lengths(n_bins, max_d)
    assert r0.shape == (n_bins,)
    assert (r0 > 0.0).all()                 # first center > 0 (no zero-r0 degeneracy)
    assert (r0 <= max_d).all()
    assert np.all(np.diff(r0) > 0.0)        # monotone increasing
    with pytest.raises(ValueError):
        ja.junc_actin_couple_bin_rest_lengths(0, max_d)
    with pytest.raises(ValueError):
        ja.junc_actin_couple_bin_rest_lengths(n_bins, 0.0)


# ---------------------------------------------------------------------------
# STUB enforcement — Hard Rule 3 (no invented constants; un-anchored = raise)
# ---------------------------------------------------------------------------
def test_enabled_unanchored_build_raises_notimplemented():
    """Enabled but un-anchored BUILD raises NotImplementedError (STUB)."""
    p = ja.resolve_junctional_actin(
        {"enabled": True}, kT=KT, dt=DT, n_cadherin=10
    )
    with pytest.raises(NotImplementedError):
        ja.extend_snapshot_with_junctional_actin(object(), p)


def test_pure_laws_raise_while_unanchored():
    """Scalar coupling / catch laws refuse to run with None constants."""
    p = ja.resolve_junctional_actin(
        {"enabled": True}, kT=KT, dt=DT, n_cadherin=1
    )
    with pytest.raises(NotImplementedError):
        ja.coupling_force_magnitude(p, 5.0e-9)
    with pytest.raises(NotImplementedError):
        ja.catch_off_rate(p, 1.0e-12)


def test_updater_refuses_unanchored():
    """The coupling Updater scaffold refuses to construct un-anchored."""
    p = ja.resolve_junctional_actin(
        {"enabled": True}, kT=KT, dt=DT, n_cadherin=4
    )
    with pytest.raises(NotImplementedError):
        ja.JunctionalActinCouplingUpdater(
            p=p, cadherin_tags=np.arange(4), n_same_cell_cortex=50
        )


def test_pi_decisions_nonempty_and_name_the_constants():
    """PI_DECISIONS lists the un-anchored constants (Hard Rule 3 honesty)."""
    assert isinstance(ja.PI_DECISIONS, list)
    assert len(ja.PI_DECISIONS) >= 1
    joined = " ".join(ja.PI_DECISIONS)
    assert "k_couple" in joined
    assert "catch" in joined.lower()


# ---------------------------------------------------------------------------
# §5 sign-sense — once anchored, the scalar laws have the right sign
# ---------------------------------------------------------------------------
def _anchored_params() -> ja.ResolvedJunctionalActin:
    """A hand-anchored params object for SIGN-SENSE testing ONLY.

    These values are NOT a sanctioned production parameter set — they are an
    arbitrary positive, dimensionally-correct catch-slip set used purely to
    exercise the sign of the scalar laws (the production constants are
    un-anchored, see PI_DECISIONS). Passed in via cfg so the resolver's own
    range checks apply; never used to build a real sim.
    """
    return ja.resolve_junctional_actin(
        {
            "enabled": True,
            "max_couple_dist": 60.0e-9,   # m  (placeholder, sign-test only)
            "anchor_r0": 30.0e-9,         # m
            "k_anchor": 1.0e-4,           # N/m
            "k_couple": 1.0e-4,           # N/m
            "k_catch0": 1.0,              # 1/s  catch dominates at low F
            "x_catch": 1.0e-9,            # m
            "k_slip0": 1.0e-3,            # 1/s
            "x_slip": 0.3e-9,             # m
            "k_on": 1.0,                  # 1/s
        },
        kT=KT,
        dt=DT,
        n_cadherin=10,
    )


def test_anchored_resolver_marks_anchored():
    p = _anchored_params()
    assert p.is_anchored is True


def test_coupling_clutch_is_tensile_only():
    """§5: stretched coupling bond pulls inward (F ≥ 0); compression = 0 N."""
    p = _anchored_params()
    # Stretched (extension > 0): positive tensile force pulling endpoints together.
    f_stretch = ja.coupling_force_magnitude(p, +5.0e-9)
    assert f_stretch > 0.0
    assert math.isclose(f_stretch, p.k_couple * 5.0e-9, rel_tol=1e-12)
    # Compressed (extension < 0): no compression push (tensile-only clutch).
    f_compress = ja.coupling_force_magnitude(p, -5.0e-9)
    assert f_compress == 0.0
    # Zero extension: zero force.
    assert ja.coupling_force_magnitude(p, 0.0) == 0.0


def test_catch_signature_off_rate_falls_then_rises():
    """§5: α-catenin/actin CATCH — k_off FALLS with force at low F (Buckley 2014).

    The defining catch signature: there is a peak force F* where k_off is
    minimal, and k_off(F*) < k_off(0). The off-rate must DECREASE just above
    zero force (the slope-sign that distinguishes a catch from a slip bond).
    """
    p = _anchored_params()
    k0 = ja.catch_off_rate(p, 0.0)

    # Analytic peak force F* of the two-pathway catch-slip law.
    f_star = (KT / (p.x_catch + p.x_slip)) * math.log(
        (p.k_catch0 * p.x_catch) / (p.k_slip0 * p.x_slip)
    )
    assert f_star > 0.0                      # a real catch regime exists
    k_star = ja.catch_off_rate(p, f_star)

    # Catch: the minimum off-rate sits BELOW the zero-force rate.
    assert k_star < k0
    # Slope just above zero force is negative (force stabilises the bond).
    f_small = 0.05 * f_star
    assert ja.catch_off_rate(p, f_small) < k0
    # Slip recovery: well beyond F* the off-rate climbs back above k0.
    assert ja.catch_off_rate(p, 5.0 * f_star) > k_star


# ---------------------------------------------------------------------------
# Hard Rule 5 — no γ contamination (denylist prefix on every bond type)
# ---------------------------------------------------------------------------
def test_gamma_denylist_prefix_nonempty():
    assert isinstance(ja.GAMMA_DENYLIST_PREFIX, str)
    assert ja.GAMMA_DENYLIST_PREFIX == "junc_actin"
    assert len(ja.GAMMA_DENYLIST_PREFIX) > 0


def test_every_bond_type_carries_denylist_prefix():
    """Every bond type the builder would create starts with the prefix."""
    names = ja.junc_actin_bond_type_names(n_bins=10)
    assert len(names) == 11             # 1 anchor + 10 coupling bins
    for name in names:
        assert name.startswith(ja.GAMMA_DENYLIST_PREFIX), name
    # Coupling-bin names specifically.
    bin_names = ja.junc_actin_couple_bin_names(7)
    assert len(bin_names) == 7
    pat = re.compile(rf"^{re.escape(ja.GAMMA_DENYLIST_PREFIX)}_couple_b\d+$")
    for name in bin_names:
        assert pat.match(name), name


if __name__ == "__main__":  # pragma: no cover
    import sys
    sys.exit(pytest.main([__file__, "-q"]))

"""Unit tests for Phase 1 Unit 4.1 cadherin Bell-Evans slip-only kinetics.

These tests exercise the off-rate formula and the stochastic bond
update in isolation, with no Cell or Cortex dependency, so any failure
here points at the cadherin module itself (and not at the integration
with Worker C's Cell or with the Young equation).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from acs_kb.cell.cell import Cell
from acs_kb.cell.cortex import generate_cortex
from acs_kb.junction.cadherin import (
    KU417_DEFAULTS,
    bond_off_rate,
    update_bonds,
)
from acs_kb.junction.types import make_ecadherin_junction


def _make_pair_junction(*, n_bonds_total: int = 100):
    ctx_a = generate_cortex(
        R_cell=10e-6, n_cortex_fibers=20, L_cortex_fiber=3e-6, beads_per_fiber=5,
        cell_center=np.array([-9e-6, 0.0]), seed=1, params={"gamma_cortex": 5e-4},
    )
    ctx_b = generate_cortex(
        R_cell=10e-6, n_cortex_fibers=20, L_cortex_fiber=3e-6, beads_per_fiber=5,
        cell_center=np.array([+9e-6, 0.0]), seed=2, params={"gamma_cortex": 5e-4},
    )
    cell_a = Cell.from_cortex(0, ctx_a)
    cell_b = Cell.from_cortex(1, ctx_b)
    return cell_a, cell_b, make_ecadherin_junction(cell_a, cell_b, n_bonds_total=n_bonds_total)


def test_bond_off_rate_zero_force_matches_kof0() -> None:
    """At F = 0 the slip rate equals the zero-force off-rate k_off0."""
    assert bond_off_rate(0.0) == pytest.approx(KU417_DEFAULTS["k_off0"])


def test_bond_off_rate_monotonic_in_force() -> None:
    """k_off(F) is strictly increasing in F (slip bond sign-sense)."""
    F = np.array([0.0, 0.5e-12, 1.0e-12, 5.0e-12, 1.0e-11])
    k = bond_off_rate(F, dx_star=0.1e-9)
    # use float() because k is a 0-d array for scalars but a 1-d for vectors here
    diffs = np.diff(k)
    assert np.all(diffs > 0), f"k_off should be monotonically increasing; got {k}"


def test_bond_off_rate_exponential_form_phase1_dx_star() -> None:
    """k_off(F) = k_off0 · exp(F·Δx*/kT) verified to float64 tolerance."""
    F = 5.0e-12
    dx_star = 0.2e-9
    expected = KU417_DEFAULTS["k_off0"] * math.exp(
        F * dx_star / KU417_DEFAULTS["kT"]
    )
    got = bond_off_rate(F, dx_star=dx_star)
    assert got == pytest.approx(expected, rel=1e-12)


def test_bond_off_rate_extreme_force_does_not_overflow() -> None:
    """``exp`` is guarded by an exponent clamp; no inf, no NaN at huge F."""
    F = 1.0e-3  # 1 mN — absurd, but must not crash.
    result = bond_off_rate(F)
    assert np.isfinite(result), f"k_off should stay finite; got {result}"
    assert result > 0


def test_update_bonds_rejects_negative_force_and_zero_dt() -> None:
    _, _, j = _make_pair_junction()
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="F_total must be non-negative"):
        update_bonds(j, F_total=-1.0e-12, dt=0.1, rng=rng)
    with pytest.raises(ValueError, match="dt must be positive"):
        update_bonds(j, F_total=0.0, dt=0.0, rng=rng)


def test_update_bonds_population_bounds() -> None:
    """``n_bonds_engaged`` stays in ``[0, n_bonds_total]`` over many steps."""
    _, _, j = _make_pair_junction(n_bonds_total=50)
    rng = np.random.default_rng(7)
    for _ in range(2000):
        update_bonds(j, F_total=0.0, dt=0.1, rng=rng)
        assert 0 <= j.n_bonds_engaged <= 50
        assert j.bond_forces.shape == (j.n_bonds_engaged,)


def test_update_bonds_zero_force_steady_state() -> None:
    """At F = 0 the steady state matches k_on / (k_on + k_off0) within MC noise."""
    _, _, j = _make_pair_junction(n_bonds_total=200)
    rng = np.random.default_rng(42)
    # Warm up to steady state, then average over a long tail.
    for _ in range(500):
        update_bonds(j, F_total=0.0, dt=0.1, rng=rng)
    samples = []
    for _ in range(1500):
        update_bonds(j, F_total=0.0, dt=0.1, rng=rng)
        samples.append(j.n_bonds_engaged)
    mean = float(np.mean(samples))
    expected = (
        KU417_DEFAULTS["k_on"]
        / (KU417_DEFAULTS["k_on"] + KU417_DEFAULTS["k_off0"])
        * 200
    )
    # ±5% tolerance — well above one-sigma binomial MC fluctuation at N = 200.
    assert mean == pytest.approx(expected, rel=0.05), (
        f"steady state {mean:.2f} ≠ analytic {expected:.2f}"
    )


def test_update_bonds_force_drives_population_lower() -> None:
    """Higher F_total reduces steady-state ``n_bonds_engaged`` (slip-bond sense)."""
    rng = np.random.default_rng(13)
    # Two parallel junctions with identical N, driven at low vs high F.
    runs = {}
    for label, F_per_bond in (("low", 1.0e-12), ("high", 2.0e-11)):
        _, _, j = _make_pair_junction(n_bonds_total=100)
        for _ in range(500):
            F_total = j.n_bonds_engaged * F_per_bond
            update_bonds(j, F_total=F_total, dt=0.1, rng=rng, dx_star=0.1e-9)
        # Average the tail.
        tail = []
        for _ in range(1000):
            F_total = j.n_bonds_engaged * F_per_bond
            update_bonds(j, F_total=F_total, dt=0.1, rng=rng, dx_star=0.1e-9)
            tail.append(j.n_bonds_engaged)
        runs[label] = float(np.mean(tail))
    assert runs["high"] < runs["low"], (
        f"slip-bond steady state should be lower under higher load; "
        f"got low={runs['low']:.1f}, high={runs['high']:.1f}"
    )


def test_update_bonds_age_accumulates_with_dt() -> None:
    _, _, j = _make_pair_junction()
    rng = np.random.default_rng(0)
    for _ in range(123):
        update_bonds(j, F_total=0.0, dt=0.05, rng=rng)
    assert j.age == pytest.approx(123 * 0.05, rel=1e-12)


def test_update_bonds_bond_forces_share_load_equally() -> None:
    """When engaged, every entry of ``bond_forces`` equals F_total / n_engaged."""
    _, _, j = _make_pair_junction()
    rng = np.random.default_rng(3)
    # Wait until at least one bond engages.
    for _ in range(50):
        update_bonds(j, F_total=2.0e-10, dt=0.1, rng=rng, dx_star=0.1e-9)
        if j.n_bonds_engaged > 0:
            break
    assert j.n_bonds_engaged > 0
    expected = 2.0e-10 / j.n_bonds_engaged
    assert np.allclose(j.bond_forces, expected, rtol=1e-12)

"""Tests for the L2.5 adiabatic catch-bond cohesion (spheroid/cadherin_bonds.py).

Pure-function checks of the effective force law + resolver derivations, plus a small HOOMD
G1-with-catch-cohesion stability run. The single-cadherin kinetics are validated separately
in test_cadherin_catch_bond.py (the oracle).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.spheroid.cadherin_bonds import (
    catch_cohesion_force,
    occupancy,
    resolve_cadherin,
    run_g1_catch,
)
from ffn_sim.spheroid.params import resolve_layer2

_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


@pytest.fixture(scope="module")
def resolved():
    return resolve_layer2(yaml.safe_load(_CONFIG.read_text()))


@pytest.fixture(scope="module")
def cad(resolved):
    return resolve_cadherin(resolved)


def test_resolve_anchors_from_measured(resolved, cad):
    """N_cad = F_detach/f0; k_bond puts per-cad force at f0 at one contact-zone extension."""
    assert cad.n_cad == pytest.approx(resolved.deadhesion_force_mature / cad.catch.f0)
    assert cad.n_cad == pytest.approx(223.0, rel=0.02)
    # per-cadherin force at ext = contact_zone equals f0 (by construction)
    f_at_zone = cad.k_bond * resolved.contact_zone_width / cad.n_cad
    assert f_at_zone == pytest.approx(cad.catch.f0, rel=1e-9)
    assert cad.wca_sigma * 2.0 ** (1.0 / 6.0) == pytest.approx(cad.r0, rel=1e-9)  # WCA min at r0


def test_occupancy_bounds_and_catch(cad):
    """φ ∈ (0,1); φ is higher at the catch force f0 than deep in the slip regime."""
    phi = occupancy(np.array([0.0, cad.catch.f0, 3.0 * cad.catch.f0]), cad)
    assert np.all((phi > 0.0) & (phi < 1.0))
    assert phi[1] > phi[2]                       # catch peak occupancy > deep-slip occupancy


def test_cohesion_force_zero_outside_window(cad):
    """No cohesion at/below r0 (WCA handles overlap) or beyond r_cut (slip-ruptured)."""
    assert catch_cohesion_force(np.array([cad.r0]), cad)[0] == 0.0
    assert catch_cohesion_force(np.array([0.9 * cad.r0]), cad)[0] == 0.0
    assert catch_cohesion_force(np.array([cad.r_cut * 1.01]), cad)[0] == 0.0


def test_cohesion_peak_at_f0(resolved, cad):
    """Force-strengthening: the holding force peaks where per-cadherin force ≈ f0 (catch)."""
    d = np.linspace(cad.r0, cad.r_cut, 600)
    F = catch_cohesion_force(d, cad)
    d_peak = d[int(np.argmax(F))]
    f_at_peak = cad.k_bond * (d_peak - cad.r0) / cad.n_cad
    assert f_at_peak == pytest.approx(cad.catch.f0, rel=0.1)   # peak at the catch force
    # peak holding force is a substantial fraction of the measured de-adhesion force
    assert 0.3 < F.max() / resolved.deadhesion_force_mature < 1.0


def test_g1_catch_is_stable_aggregate(resolved, cad):
    """A loose blob settles into a stable cohesive aggregate under the catch cohesion (G1)."""
    res = run_g1_catch(resolved, cad, n_cells=120, settle_steps=6000, measure_steps=3000, seed=2)
    assert 0.90 <= res["nn_median_over_r0"] <= 1.20      # settled to contact
    assert res["detached_fraction"] <= 0.02              # cohesive, no gas
    assert res["rg_growth_factor"] <= 1.50               # not dispersing

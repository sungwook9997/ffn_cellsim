"""Tests for the L2 D2 catch-bond emergent surface-tension measurement.

Covers the new ``--catch-bond`` path of ``scripts/layer2_emergent_sigma.py``:
the tabulated catch force sign convention, the Irving-Kirkwood spherical estimator,
and the stretch-activation finding (rest σ≈0, σ_eff(ε) engages positive under strain).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.scripts.layer2_emergent_sigma import (
    catch_pair_force_magnitude,
    ik_spherical_sigma,
    measure_emergent_sigma_catch,
)
from ffn_sim.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.spheroid.params import resolve_layer2

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


@pytest.fixture(scope="module")
def cad():
    resolved = resolve_layer2(yaml.safe_load(_CFG.read_text()))
    return resolve_cadherin(resolved), resolved


def test_catch_force_sign_convention(cad):
    """F_r>0 (repulsive) for d<r0, F_r<0 (cohesive) for r0<d<r_cut, 0 for d>=r_cut."""
    c, _ = cad
    r0, r_cut = c.r0, c.r_cut
    d = np.array([0.9 * r0, 0.5 * (r0 + r_cut), 1.001 * r_cut])
    F = catch_pair_force_magnitude(d, c)
    assert F[0] > 0.0          # excluded-volume repulsion below rest
    assert F[1] < 0.0          # catch cohesion (attractive) in the band
    assert F[2] == 0.0         # slip-ruptured beyond the cutoff


def test_catch_cohesion_dormant_at_rest(cad):
    """At d=r0 exactly, both WCA (truncated at r0) and catch cohesion (ext<=0) vanish."""
    c, _ = cad
    assert abs(float(catch_pair_force_magnitude(np.array([c.r0]), c)[0])) < 1e-18


def test_ik_empty_and_tangential(cad):
    """IK returns NaN with no pairs; a tangential cohesive pair yields σ>0."""
    c, _ = cad
    com = np.zeros(3)
    assert np.isnan(ik_spherical_sigma(np.zeros((2, 3)),
                                       np.zeros((0, 2), dtype=np.int64),
                                       np.zeros((0, 3)), 1.0, com))
    # two cells straddling the +x surface, separated tangentially (along y), pulled together
    R = c.r0
    pos = np.array([[R, 0.5 * R, 0.0], [R, -0.5 * R, 0.0]])
    idx = np.array([[0, 1]], dtype=np.int64)
    # cohesive force on 0 due to 1 points -y (toward 1): r_ij.f_ij < 0 (tension) -> σ>0
    f = np.array([[0.0, -1e-9, 0.0]])
    assert ik_spherical_sigma(pos, idx, f, R, com) > 0.0


def test_measure_catch_rest_zero_and_strain_engages(cad):
    """Tiny run: rest σ≈0 (cohesion dormant) and σ_eff(ε) rises positive under strain."""
    c, resolved = cad
    res = measure_emergent_sigma_catch(
        resolved, c, 120, settle_steps=8_000, n_snapshots=3, snapshot_interval=1_000, seed=1
    )
    # rest IK surface tension is ~0 (|σ| well under the γ band floor 0.35 mN/m)
    assert abs(res["sigma_emergent_mN_m"]) < 0.05
    strain = np.array(res["sigma_strain_mN_m"])
    grid = np.array(res["strain_grid"])
    assert grid[0] == 0.0
    # engages positive and is non-decreasing with strain (the catch holds under tension)
    assert strain[-1] > 0.0
    assert np.all(np.diff(strain) > -1e-4)

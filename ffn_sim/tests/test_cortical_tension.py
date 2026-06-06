"""Tests for the unified 3-channel cortical-tension (γ) estimator (B4 / GATE-B).

Builds a small physiological-baseline MCF7 cell (cortex + enclosed-volume
turgor + nucleus + membrane-surface compartments) and checks the three γ
channels of :func:`ffn_sim.cortex.cortical_tension.measure_cortical_tension`:

* all three channels finite;
* ``gamma_passive`` recovers the Young-Laplace turgor analytically
  (``turgor_dP0 · R_cell / 2``);
* ``gamma_soft``, ``gamma_rigid`` ≥ 0;
* ``gamma_structural == gamma_soft + gamma_rigid`` (turgor NOT folded in);
* ``gamma_passive`` is reported SEPARATELY, not inside ``gamma_structural``.

Kept fast: a small cortex (``n_filaments=120``, demo_mode) + small nucleus, no
production run — the estimator reads the constructed mechanical state directly.
"""

from __future__ import annotations

import math
from copy import deepcopy

import numpy as np
import pytest

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.cortex.cortical_tension import (
    CORTICAL_TENSION_BAND_N_PER_M,
    measure_cortical_tension,
)


@pytest.fixture(scope="module")
def small_cell():
    """A small physiological-baseline MCF7 cell for fast γ measurement."""
    m = deepcopy(load_manifest("mcf7_baseline.yaml"))
    m["cortex_overrides"] = {"cortex": {"n_filaments": 120, "demo_mode": True}}
    m["compartments"]["nucleus"]["n_beads"] = 300
    return build_baseline_cell(manifest=m, device=None)


def _measure(cell):
    handles = cell.extras["handles"]
    return measure_cortical_tension(
        cell.simulation,
        R_cell=cell.p_cortex.R_cell,
        p_enclosed_volume=cell.p_enclosed_volume,
        lambda_accumulator=handles.get("baoab_action"),
        dt=cell.simulation.operations.integrator.dt,
    )


def test_all_three_channels_finite(small_cell):
    res = _measure(small_cell)
    assert math.isfinite(res["gamma_soft"]), res["gamma_soft"]
    assert math.isfinite(res["gamma_rigid"]), res["gamma_rigid"]
    assert math.isfinite(res["gamma_passive"]), res["gamma_passive"]
    assert math.isfinite(res["gamma_structural"]), res["gamma_structural"]


def test_passive_gamma_is_young_laplace_turgor(small_cell):
    """gamma_passive == turgor_dP0 · R_cell / 2 (analytic, tight tolerance)."""
    res = _measure(small_cell)
    dP = float(small_cell.p_enclosed_volume.turgor_dP0)
    R = float(small_cell.p_cortex.R_cell)
    expected = dP * R / 2.0
    assert res["gamma_passive"] == pytest.approx(expected, rel=1e-12)
    # Manifest turgor is 133 Pa, R = 7.5e-6 m -> ~4.99e-4 N/m.
    assert res["gamma_passive"] == pytest.approx(4.9875e-4, rel=1e-6)
    # The passive channel must also report the pressure it used.
    assert res["channels"]["passive"]["dP"] == pytest.approx(dP, rel=1e-12)


def test_soft_and_rigid_nonnegative(small_cell):
    res = _measure(small_cell)
    assert res["gamma_soft"] >= 0.0, res["gamma_soft"]
    assert res["gamma_rigid"] >= 0.0, res["gamma_rigid"]


def test_structural_is_soft_plus_rigid(small_cell):
    res = _measure(small_cell)
    assert res["gamma_structural"] == pytest.approx(
        res["gamma_soft"] + res["gamma_rigid"], rel=1e-12, abs=1e-18
    )


def test_passive_reported_separately_not_in_structural(small_cell):
    """GATE-B B3 rule: turgor (passive) is NEVER folded into structural."""
    res = _measure(small_cell)
    # Structural must NOT include the passive turgor term.
    assert res["gamma_structural"] != pytest.approx(
        res["gamma_soft"] + res["gamma_rigid"] + res["gamma_passive"]
    ) or res["gamma_passive"] == 0.0
    # And the passive channel is a distinct, non-zero key (turgor is ON here).
    assert res["gamma_passive"] > 0.0
    assert "gamma_passive" in res
    assert res["channels"]["passive"]["gamma"] == res["gamma_passive"]


def test_band_constant_overlay_present(small_cell):
    """The literature band is surfaced for overlay (NOT used for tuning)."""
    res = _measure(small_cell)
    assert res["band_N_per_m"] == CORTICAL_TENSION_BAND_N_PER_M
    lo, hi = res["band_N_per_m"]
    assert 0.0 < lo < hi
    # Salbreux/Charras/Paluch 2012: [0.35, 0.65] mN/m.
    assert lo == pytest.approx(0.35e-3) and hi == pytest.approx(0.65e-3)


def test_rigid_channel_zero_when_unconstrained(small_cell):
    """Unconstrained baseline -> no Lagrange buffer -> rigid channel absent."""
    res = _measure(small_cell)
    # The baseline build is unconstrained (no M-SHAKE), so the rigid channel
    # has no λ buffer: it must report 0.0 / unavailable, not crash or NaN.
    assert res["channels"]["rigid"]["available"] is False
    assert res["gamma_rigid"] == 0.0

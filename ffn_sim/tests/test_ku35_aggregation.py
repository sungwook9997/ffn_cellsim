"""Sanity-gate tests for the KU-3.5-ACTIVE force-AGGREGATION estimator.

The estimator (``scripts/h3_ku35_aggregation.py``) is a MEASUREMENT module: it
quantifies *why* per-head myosin force does not aggregate into shell tension
(the open KU-3.5-active gate, ``docs/KU_3_5_GATE_STRUCTURE_2026-06-04.md``). A
measurement module must itself be validated against KNOWN inputs before it can
diagnose unknowns — these tests are that validation.

Per CLAUDE.md Sanity-Gate Protocol:
* DIMENSIONAL — γ in N/m, F in N, efficiencies dimensionless ∈ [0,1].
* BOUNDARY — empty bond set → 0; n_engaged=0 → NO-MOTORS verdict.
* CONSERVATION/CALIBRATION — corrected method-of-planes ≈ Irving-Kirkwood hoop
  on a synthetic isotropic shell (same 8πR² convention as the estimator audit).
* SIGN-SENSE — tension (+) builds γ, compression (−) cancels it.
* MEASUREMENT-PROTOCOL — known force → known aggregation VERDICT:
    tangential+coherent → OK; radial → η_agg≈0; compressive medium → η_medium≪1;
    tiny per-head force → GENERATION-LIMITED.
"""
from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.scripts.h3_ku35_aggregation import (
    CH_ATTACH,
    CH_XLINK,
    aggregation_ledger_from_bonds,
    bond_geometry,
    classify_channel,
    ik_ceiling,
    ik_hoop,
    mop_tension,
    self_test,
    _fib_sphere,
    _tangential_shell,
)

R = 5.0e-6  # synthetic shell radius [m]


# --------------------------------------------------------------------------- #
# Channel classification
# --------------------------------------------------------------------------- #
def test_channel_classification():
    assert classify_channel("cortex_myosin_attach_b3") == CH_ATTACH
    assert classify_channel("cortex_myosin_backbone") == "myosin_internal"
    assert classify_channel("cortex_myosin_head_backbone") == "myosin_internal"
    assert classify_channel("xlink_attach_b0") == CH_XLINK
    assert classify_channel("xlink_intra") == CH_XLINK
    assert classify_channel("cortex-bond") == "other"


# --------------------------------------------------------------------------- #
# Calibration: corrected MOP ≈ IK hoop on an isotropic tensile shell
# --------------------------------------------------------------------------- #
def test_mop_ik_calibration_isotropic_shell():
    pos, bg = _tangential_shell(1200, R)
    u, L, sin2 = bond_geometry(pos, bg)
    # k-NN tangential mesh: bonds lie ~in the tangent plane.
    assert sin2.mean() > 0.9
    T = np.full(bg.shape[0], 1e-13)
    g_mop = mop_tension(pos, bg, T, R)
    g_ik = ik_hoop(L, sin2, T, R)
    assert g_mop > 0 and g_ik > 0
    # The whole point of the 8πR² convention: the two estimators agree.
    assert 0.9 <= g_mop / g_ik <= 1.1


def test_boundary_empty_bondset():
    pos = _fib_sphere(50, R)
    empty = np.empty((0, 2), dtype=np.int64)
    assert mop_tension(pos, empty, np.empty(0), R) == 0.0
    assert ik_hoop(np.empty(0), np.empty(0), np.empty(0), R) == 0.0
    assert ik_ceiling(np.empty(0), np.empty(0), R) == 0.0


# --------------------------------------------------------------------------- #
# Measurement-protocol: known force → known aggregation verdict
# --------------------------------------------------------------------------- #
def _tangential_attach_ledger(T_val, F_stall):
    pos, bg = _tangential_shell(1200, R)
    chan = np.full(bg.shape[0], CH_ATTACH, dtype=object)
    T = np.full(bg.shape[0], T_val)
    return aggregation_ledger_from_bonds(
        pos=pos, bg=bg, channel=chan, T=T, R_cell=R, F_stall=F_stall)


def test_case_tangential_coherent_efficiency_one():
    """Tangential tensile attach with adequate generation → η_agg & η_medium ≈ 1."""
    led = _tangential_attach_ledger(2e-10, 2e-10)
    eff = led["efficiency"]
    assert eff["eta_agg"] > 0.85
    assert abs(eff["eta_medium"] - 1.0) < 0.05
    assert eff["ceiling_over_band"] > 1.0           # generation adequate
    assert "RADIAL" not in led["verdict"]
    assert "MEDIUM" not in led["verdict"]


def test_case_radial_kills_eta_agg():
    """Radial spokes (force toward centre) → η_agg ≈ 0 (no hoop tension)."""
    n_b = 800
    outer = _fib_sphere(n_b, R)
    inner = outer * 0.7
    pos = np.vstack([outer, inner])
    bg = np.stack([np.arange(n_b), np.arange(n_b) + n_b], axis=1)
    chan = np.full(n_b, CH_ATTACH, dtype=object)
    T = np.full(n_b, 5e-12)
    led = aggregation_ledger_from_bonds(
        pos=pos, bg=bg, channel=chan, T=T, R_cell=R, F_stall=5e-12)
    assert led["efficiency"]["eta_agg"] < 0.15
    assert led["tension"]["mean_tangential_fraction"] < 0.05
    # realized MOP tension ≈ 0 despite finite generation
    assert led["tension"]["g_mop_attach_Npm"] < 1e-6


def test_case_compressive_medium_cancels():
    """Tangential attach + 90% compressive xlink medium → η_medium ≪ 1, η_agg ≈ 1."""
    pos, bg = _tangential_shell(1200, R)
    bg_c = np.vstack([bg, bg])
    chan = np.array([CH_ATTACH] * bg.shape[0] + [CH_XLINK] * bg.shape[0],
                    dtype=object)
    T = np.concatenate([np.full(bg.shape[0], 2e-10),
                        np.full(bg.shape[0], -1.8e-10)])
    led = aggregation_ledger_from_bonds(
        pos=pos, bg=bg_c, channel=chan, T=T, R_cell=R, F_stall=2e-10)
    eff = led["efficiency"]
    assert eff["eta_agg"] > 0.85          # generators still tangential
    assert eff["eta_medium"] < 0.5        # medium cancels
    assert "MEDIUM" in led["verdict"]


def test_case_under_generation_is_generation_limited():
    """Tiny per-head force → coherent ceiling below band → GENERATION-LIMITED."""
    led = _tangential_attach_ledger(1e-15, 5e-12)
    assert led["efficiency"]["ceiling_over_band"] < 1.0
    assert "GENERATION" in led["verdict"]


def test_no_motors_baseline():
    """Empty attach channel → NO-MOTORS verdict (passive baseline)."""
    pos, bg = _tangential_shell(400, R)
    chan = np.full(bg.shape[0], CH_XLINK, dtype=object)  # only xlinks, no motors
    T = np.full(bg.shape[0], 1e-13)
    led = aggregation_ledger_from_bonds(
        pos=pos, bg=bg, channel=chan, T=T, R_cell=R, F_stall=5e-12)
    assert led["generation"]["n_engaged"] == 0
    assert "NO-MOTORS" in led["verdict"]


# --------------------------------------------------------------------------- #
# Sign-sense: tension builds γ, compression reduces it
# --------------------------------------------------------------------------- #
def test_sign_sense_tension_vs_compression():
    pos, bg = _tangential_shell(800, R)
    u, L, sin2 = bond_geometry(pos, bg)
    g_tension = ik_hoop(L, sin2, np.full(bg.shape[0], 1e-12), R)
    g_compress = ik_hoop(L, sin2, np.full(bg.shape[0], -1e-12), R)
    assert g_tension > 0
    assert g_compress < 0
    assert g_tension == pytest.approx(-g_compress, rel=1e-9)


# --------------------------------------------------------------------------- #
# The bundled self-test passes end to end
# --------------------------------------------------------------------------- #
def test_self_test_passes():
    assert self_test(verbose=False) is True

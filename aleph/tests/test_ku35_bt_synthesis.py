"""Tests for the KU-3.5-active binding-throughput SYNTHESIS assay.

Guards the deciding arithmetic: the active g_soft floor is GENERATION-LIMITED
(per-head force x N), NOT binding-throughput, so no literature-anchored
binding-kinetics lever crosses the band. These assertions are the auditable
"floor is real at physiological kinetics" contract.
"""

from __future__ import annotations

import math

from aleph.scripts import h3_ku35_bt_synthesis as syn


def test_self_test_passes():
    """The assay's own sanity gate must pass."""
    syn._self_test()  # raises on failure


def test_off_rate_is_overset_vs_NMIIA_literature():
    """Model k_off0=10/s is ~29x the NMIIA reference 0.35/s (Tam 2021 p37)."""
    audit = syn.synthesize()
    assert abs(audit["anchors"]["off_rate_overset_factor"] - 10.0 / 0.35) < 1e-9


def test_reanchoring_off_rate_raises_duty_not_lowers():
    """TRACK-1 falsification: lit off-rate RAISES duty toward saturation."""
    a = syn.synthesize()["anchors"]
    assert a["equilibrium_duty_lit"] > a["equilibrium_duty_model"]
    assert a["equilibrium_duty_lit"] > 0.99


def test_ideal_coherent_ceiling_is_below_band():
    """THE DECIDER: even 100% bound / fully stalled / coherent is < band."""
    ceiling = syn.absolute_coherent_ceiling()
    f_band = syn.band_required_cut_force(syn.R_CELL_MCF7_M)
    assert ceiling < f_band
    assert 0.01 < ceiling / f_band < 0.1  # ~0.03x


def test_band_required_force_matches_aggregation_assay_formula():
    """F_band = gamma_lo * 2*pi*R must match h3_ku35_aggregation.py:250."""
    R = syn.R_CELL_MCF7_M
    expected = syn.BAND_LO_NPM * 2.0 * math.pi * R
    assert abs(syn.band_required_cut_force(R) - expected) < 1e-30


def test_best_principled_lever_still_under_band():
    """The lit-off-rate + CFL-max + equilibrated lever stays far under band."""
    v = {x["scenario"]: x for x in syn.synthesize()["verdicts"]}
    lit = v["lit_off-rate + CFL-max + equilibrated"]
    assert lit["realistic_over_band"] < 1.0
    assert lit["under_band_factor"] > 10.0


def test_gsoft_does_not_scale_with_completion():
    """Completion ~doubles but g_soft moves <1 decade (generation-limited)."""
    v = {x["scenario"]: x for x in syn.synthesize()["verdicts"]}
    model = v["model_GATE-B (10/s, 40-tick)"]
    lit = v["lit_off-rate + CFL-max + equilibrated"]
    g_ratio = lit["g_soft_realistic_Npm"] / model["g_soft_realistic_Npm"]
    assert g_ratio < 10.0


def test_larger_radius_is_not_easier():
    """Band-required force grows with R, so a bigger cell stays generation-limited."""
    v = {x["scenario"]: x for x in syn.synthesize()["verdicts"]}
    assert v["IDEAL @ generic R=10um"]["is_generation_limited"]


def test_floor_real_at_physiological_kinetics():
    """The headline decision flag."""
    d = syn.synthesize()["decision"]
    assert d["floor_real_at_physiological_kinetics"] is True
    assert d["generation_limited"] is True

"""Mesoscale myosin force-scaling — active-stress CONSERVATION test (PI 2026-06-07).

The ×40 mesoscopic coarse-graining reduces the myosin minifilament count
(native 3/µm²·4πR² → n_motors_per_cell). `mesoscale_force_scaling` (Route B,
PI-ratified 2026-05-31) compensates by scaling each effective minifilament's
force ×factor so the AGGREGATE active stress is conserved:

    N_eff · F_eff  ==  N_native · F_native           (conservation)

PI (2026-06-07) flagged force-scaling=OFF as a platform-wide bug (it undercounts
the active stress by `factor` ≈ 21×). These tests pin the conservation law, the
Bell-Evans per-native-head invariance, and that the physiological PRODUCTION
manifest (mcf7_baseline.yaml) resolves Route B by default.
"""

from __future__ import annotations

import math
from copy import deepcopy

import pytest

from ffn_sim.archive.hoomd_legacy.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin


def _myosin_cfg(**over):
    """Real phase1_h3.yaml myosin block (grip_walk) + per-test overrides."""
    base = load_manifest("phase1_h3.yaml")["cortex"]["myosin"]
    cfg = deepcopy(base)
    cfg["stepping_mode"] = "grip_walk"
    cfg.update(over)
    return {"cortex": {"myosin": cfg}}


R_CELL = 7.5e-6


def test_force_scaling_conserves_active_stress():
    """N_eff·F_eff == N_native·F_native (the aggregate active stress is invariant)."""
    base = resolve_cortex_myosin(_myosin_cfg(mesoscale_force_scaling=False), dt=1e-9, R_cell=R_CELL)
    scaled = resolve_cortex_myosin(_myosin_cfg(mesoscale_force_scaling=True), dt=1e-9, R_cell=R_CELL)
    assert scaled.extras["mesoscale_force_scaling"] is True
    factor = scaled.extras["mesoscale_force_factor"]
    native_n = scaled.extras["native_n_motors"]
    # native density 0.6/µm² · 4πR² minifilaments (Nie 2015; re-anchored from the misattributed
    # 3.0 "Salbreux 2012" per PI 2026-06-07 — see H7_CORTICAL_MYOSIN_DENSITY_DATUM_2026-06-07.md).
    assert native_n == pytest.approx(0.6e12 * 4 * math.pi * R_CELL ** 2, rel=1e-9)
    assert factor == pytest.approx(native_n / 100.0, rel=1e-9)
    # CONSERVATION: effective (n_motors · F_scaled) == native (native_n · F_native).
    n_eff, n_native = base.n_motors_per_cell, native_n
    lhs = n_eff * scaled.F_stall_per_head
    rhs = n_native * base.F_stall_per_head
    assert lhs == pytest.approx(rhs, rel=1e-9), "active stress not conserved under ×40"
    # springs also scale ×factor (parallel native bundle).
    assert scaled.k_head_spring == pytest.approx(base.k_head_spring * factor, rel=1e-9)
    assert scaled.k_head_actin == pytest.approx(base.k_head_actin * factor, rel=1e-9)


def test_force_scaling_preserves_bell_evans_per_native_head():
    """x_β ÷ factor so k_off(F_eff; x_β/factor) == k_off(F_native; x_β)."""
    base = resolve_cortex_myosin(_myosin_cfg(mesoscale_force_scaling=False), dt=1e-9, R_cell=R_CELL)
    scaled = resolve_cortex_myosin(_myosin_cfg(mesoscale_force_scaling=True), dt=1e-9, R_cell=R_CELL)
    factor = scaled.extras["mesoscale_force_factor"]
    assert scaled.head_actin_x_beta == pytest.approx(base.head_actin_x_beta / factor, rel=1e-9)
    # the load-sensitivity product is invariant: F_eff·(x_β/factor) == F_native·x_β.
    assert scaled.F_stall_per_head * scaled.head_actin_x_beta == pytest.approx(
        base.F_stall_per_head * base.head_actin_x_beta, rel=1e-9)


def test_s_grip_max_invariant_under_scaling():
    """s_grip_max = F_stall/k stays physical (both ×factor → ratio invariant)."""
    base = resolve_cortex_myosin(_myosin_cfg(mesoscale_force_scaling=False), dt=1e-9, R_cell=R_CELL)
    scaled = resolve_cortex_myosin(_myosin_cfg(mesoscale_force_scaling=True), dt=1e-9, R_cell=R_CELL)
    assert (scaled.F_stall_per_head / scaled.k_head_spring) == pytest.approx(
        base.F_stall_per_head / base.k_head_spring, rel=1e-9)


def test_production_manifest_resolves_route_b():
    """The physiological production manifest must run Route B (grip_walk + scaling)."""
    m = deepcopy(load_manifest("mcf7_baseline.yaml"))
    m["cortex_overrides"]["cortex"]["n_filaments"] = 120
    m["cortex_overrides"]["cortex"]["demo_mode"] = True
    cell = build_baseline_cell(manifest=m, device=None)
    assert cell.p_myosin.stepping_mode == "grip_walk"
    assert cell.p_myosin.extras["mesoscale_force_scaling"] is True
    assert cell.p_myosin.extras["mesoscale_force_factor"] > 1.0


def test_binned_r0_stays_unscaled_legacy_default():
    """Legacy binned_r0 ignores force-scaling (byte-identical legacy proxy)."""
    cfg = _myosin_cfg(stepping_mode="binned_r0", mesoscale_force_scaling=True)
    p = resolve_cortex_myosin(cfg, dt=1e-9, R_cell=R_CELL)
    assert p.extras["mesoscale_force_scaling"] is False
    assert p.F_stall_per_head == pytest.approx(0.5e-12, rel=1e-9)

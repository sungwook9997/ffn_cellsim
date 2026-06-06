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
    ADHESION_BOND_TYPES,
    CORTICAL_TENSION_BAND_N_PER_M,
    cortical_bond_typeid_mask,
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


# --- Cortical bond-type filter (B4 refinement: exclude the adhesion load path) -


def test_cortical_bond_typeid_mask_denylist_vs_allowlist():
    """The bond-type mask: denylist drops only adhesion types; allowlist keeps
    only the named ones. Pure unit test (no cell build)."""
    bond_types = [
        "cortex-bond",
        "xlink_intra",
        "cortex_myosin_attach_b3",
        "integrin_ligand",
        "fa_actin_clutch",
        "fa_actin_clutch_b7",
        "lamel_actin_bond",  # a future actin structure (must NOT be excluded)
    ]
    # Denylist default: every type cortical EXCEPT the 3 adhesion ones.
    mask = cortical_bond_typeid_mask(bond_types, None)
    assert mask.tolist() == [True, True, True, False, False, False, True]
    # Allowlist override: only the explicitly named types are cortical.
    allow = {"cortex-bond", "xlink_intra"}
    mask2 = cortical_bond_typeid_mask(bond_types, allow)
    assert mask2.tolist() == [True, True, False, False, False, False, False]
    # Both canonical adhesion constants are recognised.
    assert "integrin_ligand" in ADHESION_BOND_TYPES
    assert "fa_actin_clutch" in ADHESION_BOND_TYPES


@pytest.fixture(scope="module")
def fa_cell():
    """A small FA-adhered MCF7 cell, settled onto the substrate.

    FA enabled, small cortex (n_filaments=120, demo_mode), nucleus n_beads at
    the manifest default 3000 (CFL), built constrained + equilibrated exactly
    like the GATE-B smoke so the focal-adhesion clutch bonds are present and
    loaded.
    """
    m = deepcopy(load_manifest("mcf7_baseline.yaml"))
    m["optional_subsystems"]["fa"]["enabled"] = True
    m["cortex_overrides"] = {"cortex": {"n_filaments": 120, "demo_mode": True}}
    # Nucleus n_beads kept at manifest default (3000) for CFL — do NOT shrink.
    return build_baseline_cell(
        manifest=m,
        device=None,
        seed=1,
        constrained=True,
        equilibrate=True,
        equilibrate_steps=120,
        equilibrate_softstart_steps=100,
    )


def test_fa_clutch_bonds_present(fa_cell):
    """Sanity: the FA build actually created adhesion (clutch) bonds, so the
    filter has something to exclude (otherwise the test is vacuous)."""
    bond_types = list(fa_cell.simulation.state.bond_types)
    assert "integrin_ligand" in bond_types
    assert any(t.startswith("fa_actin_clutch") for t in bond_types)


def test_fa_soft_channel_excludes_clutch_bonds(fa_cell):
    """gamma_soft is NOT inflated by the focal-adhesion clutch bonds.

    With the cortical bond-type filter (denylist default) the soft channel
    counts cortical/actomyosin bonds only and drops the adhesion load path
    (integrin_ligand, fa_actin_clutch[_b*]). The filtered active tension must be
    a finite, physical value far below the clutch-inflated ~49 mN/m the
    unfiltered sum produced on this FA-adhered cell (GATE-B smoke).
    """
    R = float(fa_cell.p_cortex.R_cell)
    handles = fa_cell.extras["handles"]
    act = handles.get("baoab_action")
    dt = fa_cell.simulation.operations.integrator.dt

    res = measure_cortical_tension(
        fa_cell.simulation, R_cell=R,
        p_enclosed_volume=fa_cell.p_enclosed_volume,
        lambda_accumulator=act, dt=dt,
        cortical_bond_types=None,  # robust adhesion-denylist default
    )
    soft = res["channels"]["soft"]
    # Adhesion bonds were demonstrably excluded.
    assert soft["n_bonds_excluded"] > 0, soft
    assert soft["n_bonds"] == soft["n_bonds_total"] - soft["n_bonds_excluded"]
    # The filtered active tension is finite and far below the ~49 mN/m
    # clutch-inflated GATE-B smoke number (use a generous 10 mN/m ceiling —
    # this is an exclusion check, NOT a band gate; the band stays an overlay).
    assert math.isfinite(res["gamma_soft"])
    assert res["gamma_soft"] >= 0.0
    assert res["gamma_soft"] < 10.0e-3, res["gamma_soft"]


def test_fa_unfiltered_soft_far_exceeds_filtered(fa_cell):
    """The clutch bonds are the inflators: counting them (allowlist that
    INCLUDES the adhesion types) gives a soft channel >> the cortical-only one.

    This isolates the adhesion contribution directly: build an allowlist of
    (all cortical types + the adhesion types) vs (the same set minus the
    adhesion types) and show the former is much larger.
    """
    R = float(fa_cell.p_cortex.R_cell)
    handles = fa_cell.extras["handles"]
    act = handles.get("baoab_action")
    dt = fa_cell.simulation.operations.integrator.dt
    all_types = set(fa_cell.simulation.state.bond_types)
    adhesion = {t for t in all_types
                if t in ADHESION_BOND_TYPES or t.startswith("fa_actin_clutch")}
    assert adhesion, "expected adhesion bond types on an FA-adhered cell"
    cortical_only = all_types - adhesion

    res_with = measure_cortical_tension(
        fa_cell.simulation, R_cell=R,
        p_enclosed_volume=fa_cell.p_enclosed_volume,
        lambda_accumulator=act, dt=dt,
        cortical_bond_types=all_types,  # INCLUDE adhesion (the bug behaviour)
    )
    res_without = measure_cortical_tension(
        fa_cell.simulation, R_cell=R,
        p_enclosed_volume=fa_cell.p_enclosed_volume,
        lambda_accumulator=act, dt=dt,
        cortical_bond_types=cortical_only,  # EXCLUDE adhesion (the fix)
    )
    g_with = res_with["gamma_soft"]
    g_without = res_without["gamma_soft"]
    # Adhesion-included soft channel is much larger (clutch bonds dominate the
    # FA load path); the cortical-only value is far smaller. >5× is the same
    # threshold the GATE-B inflation guard used.
    assert g_with > 5.0 * max(g_without, 1e-12), (g_with, g_without)
    # And the denylist default matches the explicit cortical-only allowlist
    # (both exclude exactly the adhesion load path).
    res_default = measure_cortical_tension(
        fa_cell.simulation, R_cell=R,
        p_enclosed_volume=fa_cell.p_enclosed_volume,
        lambda_accumulator=act, dt=dt,
        cortical_bond_types=None,
    )
    assert res_default["gamma_soft"] == pytest.approx(g_without, rel=1e-9)

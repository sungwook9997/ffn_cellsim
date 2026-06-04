"""A1 tests: ligand identity → active edge-traction (clutch-kinetics anchored ordering).

These check the ANCHORED part of A1 — that the col-I vs laminin ordering falls out of the
measured registry kinetics (not a guess), that the resolved tractions sit in the B1 stable
band, and that the resolver's boundary/sign-sense contracts hold. The density axis (Bare<Pre)
is a flagged modeling knob and is checked only for the qualitative ordering it asserts.
"""

from __future__ import annotations

import math

import pytest

from ffn_sim.bridge.ligand_species import DEFAULT_LIGAND_FOR_CONDITION, K_BT, LIGAND_REGISTRY
from ffn_sim.spheroid.ligand_traction import (
    KU_2_4_K_ON,
    LP_EDGE,
    LP_UNIFORM,
    STABLE_TRACTION_CEILING,
    T_REF_DEFAULT,
    clutch_strength,
    resolve_ligand_traction,
)

_CONDITIONS = ("Bare", "Pre", "Lam4")


def test_clutch_strength_colI_reference_is_unity():
    """The reference ligand (col-I) has clutch strength 1.0 by definition."""
    assert clutch_strength("collagen_I") == pytest.approx(1.0)


def test_laminin_clutch_weaker_than_colI_from_measured_kinetics():
    """ANCHORED ordering: laminin's φ·F_s < col-I's (higher k_off0 AND smaller F_s) → < 1.

    Matches the measured direction (breast-epithelial traction lower on LN-111). The value
    falls out of the registry kinetics, recomputed here independently from k_off0/x_β."""
    s = clutch_strength("laminin_111")
    assert 0.0 < s < 1.0
    # independent recompute: φ·F_s ratio
    def phi_Fs(lig):
        L = LIGAND_REGISTRY[lig]
        phi = KU_2_4_K_ON / (KU_2_4_K_ON + L.k_off0)
        return phi * (K_BT / L.x_beta)
    assert s == pytest.approx(phi_Fs("laminin_111") / phi_Fs("collagen_I"), rel=1e-12)
    assert s == pytest.approx(0.61, abs=0.02)  # ~0.611 with the registry anchors


def test_resolved_tractions_in_stable_band_and_finite():
    """Every condition resolves to 0 < f_traction ≤ the B1 stable ceiling (no detachment)."""
    for c in _CONDITIONS:
        r = resolve_ligand_traction(c)
        assert 0.0 < r.f_traction <= STABLE_TRACTION_CEILING
        assert math.isfinite(r.engaged_fraction) and 0.0 < r.engaged_fraction < 1.0


def test_condition_ligand_mapping_matches_registry():
    """Bare/Pre → col-I (not proxy); Lam4 → laminin-111 (proxy-flagged)."""
    assert resolve_ligand_traction("Bare").ligand == "collagen_I"
    assert resolve_ligand_traction("Pre").ligand == "collagen_I"
    assert resolve_ligand_traction("Lam4").ligand == "laminin_111"
    assert resolve_ligand_traction("Bare").proxy is False
    assert resolve_ligand_traction("Lam4").proxy is True
    for c in _CONDITIONS:
        assert resolve_ligand_traction(c).ligand == DEFAULT_LIGAND_FOR_CONDITION[c]


def test_bare_pre_same_ligand_differ_only_by_density():
    """Bare and Pre share col-I kinetics (same strength) and differ ONLY in density factor."""
    b, p = resolve_ligand_traction("Bare"), resolve_ligand_traction("Pre")
    assert b.clutch_strength == pytest.approx(p.clutch_strength)   # same ligand
    assert b.density_factor < p.density_factor                     # Bare lower availability
    assert b.f_traction < p.f_traction                            # → less traction (flagged axis)


def test_sign_sense_lower_koff_gives_more_traction():
    """Sign-sense: a lower off-rate (longer-lived bond) → higher occupancy → stronger clutch."""
    # col-I (k_off0 1.3) vs laminin (1.85): col-I lower k_off → higher strength
    assert clutch_strength("collagen_I") > clutch_strength("laminin_111")


def test_density_factor_scales_traction_linearly():
    """Boundary/scaling: f_traction ∝ density (density→0 ⇒ traction→0)."""
    full = resolve_ligand_traction("Pre", density_factors={"Pre": 1.0})
    half = resolve_ligand_traction("Pre", density_factors={"Pre": 0.5})
    assert half.f_traction == pytest.approx(0.5 * full.f_traction)


def test_a4_beta1_distribution_sets_traction_localization():
    """A4 axis: Lam4 'uniform β1' → uniform (large-Lp) traction; col-I edge-localized (LP_EDGE).

    The MAGNITUDE is unchanged by the distribution axis (still the A1 clutch-kinetics value);
    only Lp (where the traction engages: rim vs whole footprint) differs."""
    bare, pre, lam = (resolve_ligand_traction(c) for c in _CONDITIONS)
    assert bare.beta1_distribution == "diffuse" and pre.beta1_distribution == "peripheral"
    assert lam.beta1_distribution == "uniform"
    # col-I conditions are edge-localized; Lam4 is uniform (Lp far larger, engages the interior)
    assert bare.Lp == LP_EDGE and pre.Lp == LP_EDGE
    assert lam.Lp == LP_UNIFORM and lam.uniform_beta1 is True
    assert lam.Lp > 10.0 * pre.Lp                       # genuinely uniform (>> spheroid scale)
    assert bare.uniform_beta1 is False and pre.uniform_beta1 is False


def test_a4_uniform_flag_forces_uniform_without_changing_magnitude():
    """The uniform_beta1=True override forces uniform Lp on any condition; f_traction unchanged."""
    edge = resolve_ligand_traction("Pre")
    forced = resolve_ligand_traction("Pre", uniform_beta1=True)
    assert edge.Lp == LP_EDGE and forced.Lp == LP_UNIFORM and forced.uniform_beta1 is True
    assert forced.f_traction == pytest.approx(edge.f_traction)   # distribution ≠ magnitude


def test_resolver_rejects_bad_inputs():
    """Unknown condition, non-positive t_ref, and over-ceiling traction all raise (surfaced)."""
    with pytest.raises(ValueError):
        resolve_ligand_traction("FN")
    with pytest.raises(ValueError):
        resolve_ligand_traction("Pre", t_ref=0.0)
    with pytest.raises(ValueError):  # density 10× pushes col-I past the 3 nN ceiling
        resolve_ligand_traction("Pre", density_factors={"Pre": 10.0})


def test_active_traction_bridge_anchors_and_ceiling():
    """Lamellipodium→CBM bridge: protrusion (Bieling) vs whole-cell anchors + ceiling flag."""
    import yaml
    from pathlib import Path
    from ffn_sim.spheroid.params import resolve_layer2
    from ffn_sim.spheroid.ligand_traction import (
        resolve_active_traction, STABLE_TRACTION_CEILING,
    )
    cfg = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"
    r = resolve_layer2(yaml.safe_load(cfg.read_text()))
    a = resolve_active_traction(r)
    # f = v0 * gamma_cell (overdamped); Bieling 31.3 nm/s * 0.30 N*s/m ≈ 9.4 nN
    assert a.f_protrusion == pytest.approx(a.v0_protrusion * r.gamma_cell, rel=1e-9)
    assert a.f_wholecell == pytest.approx(a.v0_wholecell * r.gamma_cell, rel=1e-9)
    assert 8e-9 < a.f_protrusion < 11e-9          # ~9.4 nN
    assert 1e-9 < a.f_wholecell < 2.5e-9          # ~1.6 nN
    # the literature-first protrusion anchor exceeds the B1 overdamped ejection ceiling
    assert a.protrusion_exceeds_ceiling is True
    assert a.ceiling == STABLE_TRACTION_CEILING
    assert a.f_protrusion / a.f_wholecell == pytest.approx(11.6 * 2.7e-9 / (0.32e-6 / 60.0), rel=1e-6)

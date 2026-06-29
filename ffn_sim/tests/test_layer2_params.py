"""Tests for ffn_sim.archive.hoomd_legacy.spheroid.params.resolve_layer2 — the Magic-Number-Block derivations.

Loads the real runtime config (configs/layer2_cbm.yaml) and verifies every DERIVED SI
quantity is computed correctly from the literature-anchored primaries. Pure python, no
HOOMD. This is the "derivations verified before the physics run" discipline.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2

_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


@pytest.fixture(scope="module")
def resolved():
    cfg = yaml.safe_load(_CONFIG.read_text())
    return resolve_layer2(cfg)


def test_cell_radius_is_half_diameter(resolved):
    assert resolved.R_cell == pytest.approx(resolved.diameter / 2.0)
    assert resolved.R_cell == pytest.approx(7.5e-6)  # MCF7 15 µm diameter


def test_adhesion_well_depth_from_measured_deadhesion_force(resolved):
    """D_e set so the Morse max force (D_e·alpha/2) == measured MCF7-MCF7 de-adhesion force.

    Anchor: Iturri 2020 (MCF7-MCF7 SCFS, measured) ~6-7 nN mature. The Morse cohesion is
    therefore nN-scale (balances nN traction) — NOT the retired pN single-bond seed.
    """
    expected = 2.0 * resolved.deadhesion_force_mature * resolved.contact_zone_width
    assert resolved.D_e == pytest.approx(expected)
    morse_max_force = resolved.D_e * resolved.morse_alpha / 2.0
    assert morse_max_force == pytest.approx(resolved.deadhesion_force_mature, rel=1e-9)
    assert resolved.deadhesion_force_mature == pytest.approx(6.5e-9, rel=1e-9)  # nN-scale


def test_morse_rest_separation_is_diameter(resolved):
    """Two R_cell spheres touch (surfaces) at centre-to-centre 2R = diameter."""
    assert resolved.morse_r0 == pytest.approx(resolved.diameter)
    assert resolved.morse_r0 == pytest.approx(2.0 * resolved.R_cell)


def test_cortical_tension_is_config_input_not_cortex_import(resolved):
    """Layer-2 consumes γ as a physiological config value; Track B produces it elsewhere."""
    assert resolved.cortical_tension == pytest.approx(5.7e-4)
    lo, hi = resolved.cortical_tension_band
    assert lo <= resolved.cortical_tension <= hi

    spheroid_dir = Path(__file__).resolve().parents[1] / "spheroid"
    offenders = []
    for path in spheroid_dir.glob("*.py"):
        text = path.read_text()
        if "ffn_sim.archive.hoomd_legacy.cortex" in text:
            offenders.append(path.name)
    assert offenders == []


def test_morse_alpha_inverse_contact_zone(resolved):
    assert resolved.morse_alpha == pytest.approx(1.0 / resolved.contact_zone_width)


def test_morse_spring_constant_is_curvature_at_minimum(resolved):
    """k_spring = U''(r0) = 2 D_e alpha² for the Morse potential."""
    expected = 2.0 * resolved.D_e * resolved.morse_alpha**2
    assert resolved.morse_k_spring == pytest.approx(expected)


def test_per_cell_migration_drag_from_clutch(resolved):
    """gamma = n_eng·kappa/k_off (clutch-ensemble MIGRATION drag, KU-2.18) ~ 0.3 N·s/m.

    RETIRES water-Stokes (~1e-7, ~6.5 OOM too small for crawling). Sanity: a 1 nN net
    traction => v = F/gamma ~ µm/min (slow epithelial MCF7), not mm/s.
    """
    assert resolved.gamma_cell == pytest.approx(0.30, rel=0.05)
    v_at_1nn = 1e-9 / resolved.gamma_cell  # m/s
    assert 1e-9 < v_at_1nn < 1e-7  # ~0.06-6 µm/min, physical (not mm/s)


def test_cfl_timestep_positive_and_overdamped(resolved):
    """dt = safety · gamma / k_spring; finite, positive, and sub-relaxation-time."""
    tau = resolved.gamma_cell / resolved.morse_k_spring
    assert resolved.dt_cfl == pytest.approx(resolved.cfl_safety_factor * tau)
    assert 0.0 < resolved.dt_cfl < tau  # CFL-safe (safety < 1)


def test_resolve_rejects_nonpositive_primary():
    bad = {
        "spheroid": {
            "temperature": 310.0, "kT": 4.28e-21, "water_viscosity": 6.913e-4, "seed": 1,
            "cell": {"diameter": -1.0},  # invalid
            "mechanics": {
                "cortical_tension": 5.7e-4,
                "cortical_tension_band": [3.5e-4, 6.5e-4],
            },
            "adhesion": {
                "deadhesion_force_mature": 6.5e-9, "deadhesion_force_nascent": 1.5e-9,
                "contact_zone_width": 1.5e-6,
            },
            "dynamics": {"cfl_safety_factor": 0.1},
        }
    }
    with pytest.raises(ValueError):
        resolve_layer2(bad)

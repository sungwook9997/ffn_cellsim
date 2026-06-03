"""Static tests for the single-cell-gamma -> spheroid surface-tension BRIDGE.

Validates the closed-form bridge oracle
``ffn_sim.validation.oracles.spheroid.surface_tension_bridge`` (Young-Laplace, DITH
Gamma=cortical-adhesion, Fastabend R=lambda/sigma, Okuda 3D-cap, Roffay outer/interior ratio)
and the runtime virial-pressure observable ``ffn_sim.spheroid.observables.virial_pressure``
against synthetic configurations with KNOWN closed-form answers. Pure numpy/scipy, no HOOMD.

Provenance for the relations under test: ``docs/CORTICAL_TENSION_TRIAGE_2026-06-03.md``.
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.spheroid import observables as obs
from ffn_sim.validation.oracles.spheroid import surface_tension_bridge as br

RNG = np.random.default_rng(20260603)


# --------------------------------------------------------------------------- #
# Oracle: Gamma = cortical - adhesion (DITH; Okuda)
# --------------------------------------------------------------------------- #
def test_adhesion_tension_is_energy_density():
    # 200 fJ over a 7 um^2 contact -> beta = W / A.
    beta = br.adhesion_tension(2.0e-13, 7.0e-12)
    assert beta == pytest.approx(2.0e-13 / 7.0e-12)


def test_adhesion_tension_guards():
    with pytest.raises(ValueError):
        br.adhesion_tension(-1.0, 1.0)
    with pytest.raises(ValueError):
        br.adhesion_tension(1.0, 0.0)


def test_interfacial_tension_reduces_cortical_by_adhesion():
    assert br.interfacial_tension(0.57e-3, 0.20e-3) == pytest.approx(0.37e-3)


def test_interfacial_tension_clamps_at_zero_on_wetting():
    # adhesion > cortical => wetting, no cohesive surface => clamp to 0 (not negative).
    assert br.interfacial_tension(0.3e-3, 0.5e-3) == 0.0


def test_aggregate_surface_tension_equals_cortical():
    assert br.aggregate_surface_tension(0.57e-3) == pytest.approx(0.57e-3)


def test_surface_interior_ratio_hits_roffay_band():
    # Roffay outer/interior ~1.6-2.0 is reproduced at beta/gamma ~ 0.375-0.5.
    gamma = 0.57e-3
    lo, hi = br.ROFFAY_SURFACE_INTERIOR_RATIO
    r_lo = br.surface_interior_ratio(gamma, 0.375 * gamma)  # -> 1/(1-0.375)=1.6
    r_hi = br.surface_interior_ratio(gamma, 0.50 * gamma)   # -> 1/(1-0.5)  =2.0
    assert r_lo == pytest.approx(lo, rel=1e-6)
    assert r_hi == pytest.approx(hi, rel=1e-6)


def test_surface_interior_ratio_infinite_on_full_wetting():
    assert br.surface_interior_ratio(0.3e-3, 0.5e-3) == float("inf")


# --------------------------------------------------------------------------- #
# Oracle: Young-Laplace (Roffay 3D) + inverse + Fastabend radius balance
# --------------------------------------------------------------------------- #
def test_young_laplace_sphere_is_two_sigma_over_R():
    sigma, R = 0.57e-3, 50e-6
    assert br.young_laplace_pressure(sigma, R) == pytest.approx(2.0 * sigma / R)


def test_young_laplace_two_radii_mean_curvature():
    sigma, R, R2 = 1.0e-3, 40e-6, 80e-6
    assert br.young_laplace_pressure(sigma, R, R2) == pytest.approx(
        sigma * (1.0 / R + 1.0 / R2)
    )


def test_pressure_inversion_roundtrips():
    sigma, R, R2 = 0.57e-3, 33e-6, 71e-6
    dP = br.young_laplace_pressure(sigma, R, R2)
    assert br.surface_tension_from_pressure(dP, R, R2) == pytest.approx(sigma)


def test_young_laplace_vanishes_at_infinite_radius():
    # flat interface -> no Laplace overpressure (boundary case).
    assert br.young_laplace_pressure(0.57e-3, 1e30) == pytest.approx(0.0, abs=1e-20)


def test_laplace_radius_balance():
    # R = lambda / sigma (Fastabend 2D analog) -> dimensional [N]/[N/m] = [m].
    assert br.laplace_radius(1.14e-8, 0.57e-3) == pytest.approx(1.14e-8 / 0.57e-3)


def test_oracle_guards():
    for fn in (br.young_laplace_pressure, br.surface_tension_from_pressure):
        with pytest.raises(ValueError):
            fn(1.0, -1.0)
        with pytest.raises(ValueError):
            fn(-1.0, 1.0)
    with pytest.raises(ValueError):
        br.laplace_radius(0.0, 1.0)


# --------------------------------------------------------------------------- #
# Oracle: Okuda 3D-cap criterion (underwrites L2.6)
# --------------------------------------------------------------------------- #
def test_three_d_cap_when_free_surface_exceeds_threshold():
    # free-surface tension > 0.2 * cell-cell tension => 3D cap (Okuda).
    assert br.is_three_d_cap(0.57e-3, 0.37e-3) is True       # 0.57 > 0.2*0.37
    assert br.is_three_d_cap(0.05e-3, 0.57e-3) is False      # 0.05 < 0.2*0.57


# --------------------------------------------------------------------------- #
# Observable: virial pressure (the emergent-tension measurement route)
# --------------------------------------------------------------------------- #
def test_convex_hull_volume_of_cube():
    # 8 corners of a cube of side L -> hull volume L^3.
    L = 3.0
    corners = np.array(
        [[x, y, z] for x in (0, L) for y in (0, L) for z in (0, L)], dtype=float
    )
    assert obs.convex_hull_volume(corners) == pytest.approx(L**3)


def test_convex_hull_volume_degenerate():
    assert obs.convex_hull_volume(np.zeros((3, 3))) == 0.0  # < 4 points


def test_virial_pressure_repulsive_is_positive():
    # two cells pushed apart: f on i points along +r_ij => W>0 => P>0.
    pos = np.array([[0.0, 0, 0], [1e-6, 0, 0], [0, 1e-6, 0], [0, 0, 1e-6]])
    pairs = np.array([[1, 0]])
    forces = np.array([[1e-9, 0, 0]])  # force on cell 1 (at +x) points further +x
    P = obs.virial_pressure(pos, pairs, forces, volume=1e-18)
    assert P > 0.0


def test_virial_pressure_cohesive_is_negative():
    # cohesive pair: force on i points toward j (opposite r_ij) => W<0 => P<0 (tension).
    pos = np.array([[0.0, 0, 0], [1e-6, 0, 0], [0, 1e-6, 0], [0, 0, 1e-6]])
    pairs = np.array([[1, 0]])
    forces = np.array([[-1e-9, 0, 0]])  # cell 1 pulled back toward cell 0
    P = obs.virial_pressure(pos, pairs, forces, volume=1e-18)
    assert P < 0.0


def test_virial_pressure_matches_closed_form():
    # P = (N*kT + W/3)/V with a hand-computed W.
    pos = np.array([[0.0, 0, 0], [2e-6, 0, 0], [0, 2e-6, 0], [0, 0, 2e-6]])
    pairs = np.array([[1, 0], [2, 0]])
    forces = np.array([[3e-9, 0, 0], [0, 5e-9, 0]])
    r10 = pos[1] - pos[0]
    r20 = pos[2] - pos[0]
    W = float(r10 @ forces[0] + r20 @ forces[1])
    V = 4.0e-18
    assert obs.virial_pressure(pos, pairs, forces, volume=V) == pytest.approx(W / 3.0 / V)


def test_virial_pressure_ideal_gas_kinetic_term():
    # no pairs, only the kinetic term: P = N*kT/V.
    pos = RNG.uniform(0, 1e-6, size=(10, 3))
    P = obs.virial_pressure(pos, np.empty((0, 2)), np.empty((0, 3)), volume=1e-18, kT=4e-21)
    assert P == pytest.approx(10 * 4e-21 / 1e-18)


def test_virial_pressure_guards():
    pos = np.array([[0.0, 0, 0], [1e-6, 0, 0], [0, 1e-6, 0], [0, 0, 1e-6]])
    with pytest.raises(ValueError):  # bad pair shape
        obs.virial_pressure(pos, np.array([[0]]), np.zeros((1, 3)), volume=1.0)
    with pytest.raises(ValueError):  # force/pair mismatch
        obs.virial_pressure(pos, np.array([[0, 1]]), np.zeros((2, 3)), volume=1.0)
    with pytest.raises(ValueError):  # index out of range
        obs.virial_pressure(pos, np.array([[0, 9]]), np.zeros((1, 3)), volume=1.0)
    with pytest.raises(ValueError):  # non-positive volume
        obs.virial_pressure(pos, np.array([[0, 1]]), np.zeros((1, 3)), volume=0.0)


# --------------------------------------------------------------------------- #
# End-to-end: a synthetic cohesive droplet recovers its input surface tension
# --------------------------------------------------------------------------- #
def test_youngLaplace_roundtrip_on_synthetic_droplet():
    """A droplet with imposed interior overpressure dP and radius R recovers sigma = dP*R/2.

    This is the measurement chain the L2 D2 gate uses: observable measures dP (virial), oracle
    inverts Young-Laplace to sigma. Here we feed a known dP and confirm the oracle returns the
    sigma that produced it.
    """
    sigma_true, R = 0.57e-3, 50e-6
    dP = br.young_laplace_pressure(sigma_true, R)        # forward (what a droplet would show)
    sigma_meas = br.surface_tension_from_pressure(dP, R)  # inverse (the measurement)
    assert sigma_meas == pytest.approx(sigma_true, rel=1e-9)

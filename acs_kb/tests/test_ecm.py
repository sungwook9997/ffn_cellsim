"""Phase 1 Unit 1.1 validation tests (post P1 corrections).

Six KU-driven checks:

1. Generation contract — KU-1.27 / KU-1.22.
2. Nematic order — KU-1.9.
3. Emergent ⟨z⟩ matches Mikado theory — KU-1.27 (theory ± periodic loss).
4. Emergent ℓ_c matches the inverted target — KU-1.27 self-consistency.
5. Backbone + XL force analytic vs finite difference — KU-1.24, KU-1.28.
6. Full Phase 1 ECM Sanity Gate passes (and reports the biology gap).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from acs_kb.common.derived_params import load_config
from acs_kb.common.sanity_gate import gate_phase1_ecm
from acs_kb.ecm.cross_links import (
    compute_xl_energy_and_forces,
    generate_cross_links,
    measure_coordination,
    measure_xl_per_fiber,
)
from acs_kb.ecm.fiber_mechanics import compute_energy, compute_forces
from acs_kb.ecm.fiber_network import (
    generate_2d_fiber_network,
    measure_nematic_order,
)

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_unit1.yaml"


@pytest.fixture(scope="module")
def cfg():
    return load_config(CONFIG_PATH)


@pytest.fixture(scope="module")
def network(cfg):
    ecm = cfg["ecm"]
    d = ecm["derived"]
    return generate_2d_fiber_network(
        L_box=ecm["L_box"],
        n_fibers=d["n_fibers"],
        L_fiber=ecm["L_fiber"],
        beads_per_fiber=ecm["beads_per_fiber"],
        S_order=0.0,
        seed=ecm["seed"],
    )


@pytest.fixture(scope="module")
def links(cfg, network):
    return generate_cross_links(network, stiffness=cfg["ecm"]["xl_stiffness"])


def test_generation_contract(cfg, network):
    """KU-1.27 / KU-1.22: shape, wrap, rest length, endpoints exist."""
    ecm = cfg["ecm"]
    d = ecm["derived"]
    n_f, n_b, dim = network.bead_positions.shape
    assert n_f == d["n_fibers"]
    assert n_b == ecm["beads_per_fiber"]
    assert dim == 2
    assert 0.0 <= network.bead_positions.min()
    assert network.bead_positions.max() < ecm["L_box"]
    assert network.rest_length == pytest.approx(d["rest_length"], rel=1e-12)
    assert network.fiber_endpoints.shape == (d["n_fibers"], 2, 2)


def test_nematic_order_aligned(cfg):
    """KU-1.9: aligned variant (S=0.5) within statistical 1/√(2N) band."""
    ecm = cfg["ecm"]
    d = ecm["derived"]
    net = generate_2d_fiber_network(
        L_box=ecm["L_box"], n_fibers=d["n_fibers"], L_fiber=ecm["L_fiber"],
        beads_per_fiber=ecm["beads_per_fiber"], S_order=0.5,
        seed=ecm["seed"] + 2,
    )
    S = measure_nematic_order(net.fiber_orientations)
    sigma = 1.0 / np.sqrt(2.0 * d["n_fibers"])
    assert abs(S - 0.5) <= 3 * sigma, f"S={S:.4f}, target 0.5, σ={sigma:.4f}"


def test_emergent_z_matches_mikado_theory(cfg, network, links):
    """KU-1.27 / KU-1.3: measured ⟨z⟩ within 10% of theoretical prediction.

    Discrepancy budget = periodic-boundary loss (~5% boundary fibers in
    Phase 1 default), absorbed by the 10% tolerance.
    """
    ecm = cfg["ecm"]
    d = ecm["derived"]
    n_beads_total = network.bead_positions.shape[0] * network.bead_positions.shape[1]
    z = measure_coordination(links, n_beads_total, ecm["beads_per_fiber"])
    expected = d["expected_total_z"]
    assert abs(z - expected) / expected < 0.10, \
        f"emergent ⟨z⟩={z:.4f}, predicted {expected:.4f}"


def test_emergent_segment_length_matches_target(cfg, network, links):
    """KU-1.27 self-consistency: measured ℓ_c matches target ξ within 10%."""
    ecm = cfg["ecm"]
    target = ecm["target_segment_length"]
    # Empirical ℓ_c from realized links: ℓ_c = L_fiber / (n_int per fiber)
    xl_per_f = measure_xl_per_fiber(links, ecm["derived"]["n_fibers"])
    measured_ell_c = ecm["L_fiber"] / max(xl_per_f, 1e-9)
    assert abs(measured_ell_c - target) / target < 0.10, \
        f"emergent ℓ_c={measured_ell_c:.3e} m, target {target:.3e}"


def test_analytic_forces_match_finite_difference(cfg):
    """KU-1.24 + KU-1.28: backbone + XL gradient agrees with FD (<1e-6)."""
    ecm = cfg["ecm"]
    rng = np.random.default_rng(0)
    net = generate_2d_fiber_network(
        L_box=ecm["L_box"], n_fibers=80, L_fiber=ecm["L_fiber"],
        beads_per_fiber=ecm["beads_per_fiber"], S_order=0.0, seed=1,
    )
    links = generate_cross_links(net, stiffness=ecm["xl_stiffness"])
    p = net.bead_positions + rng.normal(scale=2e-7, size=net.bead_positions.shape)
    p = np.mod(p, net.box_size)
    mu = ecm["stretching_modulus"]; kappa = ecm["bending_modulus"]; L = net.box_size
    F_an = compute_forces(p, net.rest_length, mu, kappa, L, cross_links=links)
    # Newton 3rd law (network + links is a closed bonded system):
    assert np.allclose(F_an.sum(axis=(0, 1)), 0.0, atol=1e-20)
    # Sample 100 components rather than all 800 for speed.
    h = 1e-11
    errs = []
    coords = [(f, b, c) for f in range(p.shape[0]) for b in range(p.shape[1])
              for c in range(2)]
    for f, b, c in coords[:100]:
        pp = p.copy(); pp[f, b, c] += h
        pm = p.copy(); pm[f, b, c] -= h
        Ep = compute_energy(pp, net.rest_length, mu, kappa, L, cross_links=links)
        Em = compute_energy(pm, net.rest_length, mu, kappa, L, cross_links=links)
        F_fd = -(Ep - Em) / (2 * h)
        F_a = F_an[f, b, c]
        errs.append(abs(F_fd - F_a) / (abs(F_a) + 1e-30))
    assert max(errs) < 1e-6, f"max rel err {max(errs):.2e}"


def test_xl_forces_pairwise_newton_third_law(cfg):
    """KU-1.28: each cross-link's two-bead force is equal and opposite."""
    ecm = cfg["ecm"]
    rng = np.random.default_rng(0)
    net = generate_2d_fiber_network(
        L_box=ecm["L_box"], n_fibers=80, L_fiber=ecm["L_fiber"],
        beads_per_fiber=ecm["beads_per_fiber"], S_order=0.0, seed=1,
    )
    links = generate_cross_links(net, stiffness=ecm["xl_stiffness"])
    p = net.bead_positions + rng.normal(scale=2e-7, size=net.bead_positions.shape)
    p = np.mod(p, net.box_size)
    _, F_xl = compute_xl_energy_and_forces(p, links, net.box_size)
    assert np.allclose(F_xl.sum(axis=(0, 1)), 0.0, atol=1e-20)


def test_sanity_gate_passes_on_resolved_defaults(cfg, network, links):
    """CLAUDE.md hard rule: Phase 1 ECM gate passes; biology gap logged."""
    ecm = cfg["ecm"]
    d = ecm["derived"]
    n_beads_total = network.bead_positions.shape[0] * network.bead_positions.shape[1]
    z = measure_coordination(links, n_beads_total, ecm["beads_per_fiber"])
    xl_per_f = measure_xl_per_fiber(links, d["n_fibers"])
    measured_ell_c = ecm["L_fiber"] / max(xl_per_f, 1e-9)
    rep = gate_phase1_ecm(
        biological_mesh_target=ecm["biological_mesh"],
        measured_segment_length=measured_ell_c,
        expected_segment_length=d["mikado_ell_c_predicted"],
        measured_z=z,
        expected_z=d["expected_total_z"],
        ku13_reference_z=ecm["target_coordination_ref"],
        n_fibers=d["n_fibers"], n_beads_total=n_beads_total, n_links=len(links),
        acceptance=ecm["acceptance"], demo_mode=ecm["demo_mode"],
    )
    assert rep.passed, rep.summary()

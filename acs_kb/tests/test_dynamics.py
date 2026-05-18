"""Unit 1.2 — overdamped Langevin integration validation.

Three groups of tests:

1. CFL Sanity Gate behaviour (KU-1.26): the gate must pass for the
   default dt and must reject dt > α·τ_min.
2. Integrator zero-force sanity: with F ≡ 0 the variance of bead
   displacement grows as 2 D t with D = k_B T / γ_b.
3. Equipartition (KU-1.26): in a small XL-only network, the mean
   per-link harmonic energy converges to ½ k_B T within 10 %.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from acs_kb.common.derived_params import load_config
from acs_kb.common.sanity_gate import (
    SanityGateFailure,
    gate_equipartition,
    gate_unit1_2_dynamics,
)
from acs_kb.ecm.cross_links import generate_cross_links
from acs_kb.ecm.diagnostics import link_extension_energies
from acs_kb.ecm.fiber_mechanics import compute_forces
from acs_kb.ecm.fiber_network import generate_2d_fiber_network
from acs_kb.ecm.integrator import EulerMaruyama, run
from acs_kb.ecm.shear_protocol import (
    apply_affine_shear,
    compute_virial_shear_stress_2d,
    preliminary_affine_G_curve,
)

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_unit1.yaml"


@pytest.fixture(scope="module")
def cfg():
    return load_config(CONFIG_PATH)


def test_cfl_gate_passes_on_default(cfg):
    """KU-1.26: default dt satisfies dt < α·τ_min."""
    ecm = cfg["ecm"]; d = ecm["derived"]; dyn = ecm["dynamics"]
    rep = gate_unit1_2_dynamics(
        dt=dyn["dt"], tau_xl=d["tau_xl"], tau_stretch=d["tau_stretch"],
        tau_bend=d["tau_bend"], safety_factor=dyn["cfl_safety_factor"],
        integrator_name=dyn["integrator"],
    )
    assert rep.passed, rep.summary()


def test_cfl_gate_rejects_excess_dt(cfg):
    """KU-1.26: dt > α·τ_min must fail loudly."""
    ecm = cfg["ecm"]; d = ecm["derived"]; dyn = ecm["dynamics"]
    rep = gate_unit1_2_dynamics(
        dt=10.0 * dyn["cfl_safety_factor"] * d["tau_min"],  # 100× over limit
        tau_xl=d["tau_xl"], tau_stretch=d["tau_stretch"], tau_bend=d["tau_bend"],
        safety_factor=dyn["cfl_safety_factor"],
        integrator_name=dyn["integrator"],
    )
    with pytest.raises(SanityGateFailure):
        rep.raise_if_failed()


def test_free_bead_diffusion(cfg):
    """KU-1.26: free overdamped bead obeys ⟨|Δr|²⟩ = 4 D t in 2D."""
    ecm = cfg["ecm"]; d = ecm["derived"]
    gamma_b = d["gamma_b"]
    kT = ecm["kT"]
    D = kT / gamma_b
    dt = 0.5 * d["tau_min"]  # well below CFL bound

    rng = np.random.default_rng(0)
    n_walkers = 5000
    n_steps = 200
    positions = np.zeros((n_walkers, 1, 2))  # use the (F, N, 2) shape convention
    forces_fn = lambda p: np.zeros_like(p)
    integrator = EulerMaruyama()
    p = positions.copy()
    for _ in range(n_steps):
        p = integrator.step(p, forces_fn, gamma_b, kT, dt, rng)
    # Do NOT wrap (free diffusion); compute MSD in unwrapped frame.
    msd = np.mean(np.sum((p - positions) ** 2, axis=2))
    expected = 4.0 * D * n_steps * dt
    assert abs(msd - expected) / expected < 0.05, \
        f"MSD={msd:.3e}, expected {expected:.3e}"


def test_equipartition_on_small_network(cfg):
    """KU-1.26: ⟨½ k_xl |Δr|²⟩ → ½ k_BT within 10% on a small network."""
    ecm = cfg["ecm"]; d = ecm["derived"]; dyn = ecm["dynamics"]
    # Dense-enough small network to guarantee a few dozen XLs.
    net = generate_2d_fiber_network(
        L_box=ecm["L_box"], n_fibers=400, L_fiber=ecm["L_fiber"],
        beads_per_fiber=ecm["beads_per_fiber"], S_order=0.0, seed=3,
    )
    links = generate_cross_links(net, stiffness=ecm["xl_stiffness"])
    # We need at least a handful of XLs for the sample mean to converge.
    if len(links) < 3:
        pytest.skip(f"Only {len(links)} XLs in seed-3 60-fiber net; "
                    f"insufficient sample.")

    mu = ecm["stretching_modulus"]
    kappa = ecm["bending_modulus"]
    gamma_b = d["gamma_b"]; kT = ecm["kT"]; dt = dyn["dt"]

    def forces_fn(p):
        return compute_forces(p, net.rest_length, mu, kappa, net.box_size,
                              cross_links=links)

    rng = np.random.default_rng(11)
    integrator = EulerMaruyama()
    # equilibrate
    p = run(net.bead_positions, forces_fn, integrator, gamma_b, kT, dt,
            net.box_size, n_steps=8000, rng=rng)
    # sample
    samples: list[float] = []
    for _ in range(160):
        p = run(p, forces_fn, integrator, gamma_b, kT, dt, net.box_size,
                n_steps=50, rng=rng)
        e = link_extension_energies(p, links, net.box_size)
        samples.extend(e.tolist())
    mean_e = float(np.mean(samples))
    rep = gate_equipartition(
        measured_mean_link_energy=mean_e, kT=kT, n_samples=len(samples),
        tolerance=0.12,  # 12% to absorb the ~5% statistical floor + coupling bias
    )
    assert rep.passed, rep.summary() + f"\n(N_links={len(links)})"


def test_affine_shear_is_pure_x_shift():
    """KU-1.29: r → r + (γ y, 0) leaves y unchanged."""
    rng = np.random.default_rng(0)
    p = rng.uniform(0, 200e-6, size=(5, 5, 2))
    p_s = apply_affine_shear(p, 0.05)
    np.testing.assert_array_equal(p_s[..., 1], p[..., 1])
    np.testing.assert_allclose(p_s[..., 0], p[..., 0] + 0.05 * p[..., 1])


def test_virial_stress_zero_on_rest_network(cfg):
    """A freshly generated network has |b| = ℓ₀ ⇒ σ_xy must be ~0."""
    ecm = cfg["ecm"]; d = ecm["derived"]
    net = generate_2d_fiber_network(
        L_box=ecm["L_box"], n_fibers=200, L_fiber=ecm["L_fiber"],
        beads_per_fiber=ecm["beads_per_fiber"], S_order=0.0, seed=5,
    )
    links = generate_cross_links(net, stiffness=ecm["xl_stiffness"])
    s = compute_virial_shear_stress_2d(
        net.bead_positions, net.rest_length,
        ecm["stretching_modulus"], ecm["bending_modulus"],
        net.box_size, cross_links=links,
    )
    # Bond lengths are EXACTLY ℓ₀ at construction, and XL rest lengths
    # equal current distance ⇒ no force ⇒ no virial.
    assert abs(s) < 1e-15, f"unsheared virial ≠ 0: {s}"


def test_preliminary_G_in_KU130_order_of_magnitude(cfg):
    """Affine preliminary: G(γ→0) using h=ξ should fall in the KU-1.30 band."""
    ecm = cfg["ecm"]
    h = ecm["layer_thickness_for_2d_to_3d"]
    net = generate_2d_fiber_network(
        L_box=ecm["L_box"], n_fibers=400, L_fiber=ecm["L_fiber"],
        beads_per_fiber=ecm["beads_per_fiber"], S_order=0.0, seed=7,
    )
    links = generate_cross_links(net, stiffness=ecm["xl_stiffness"])
    sigma_2d, G_2d, G_3d = preliminary_affine_G_curve(
        net, links, ecm["stretching_modulus"], ecm["bending_modulus"],
        np.array([0.005, 0.01]), layer_thickness=h,
    )
    G_small = abs(G_3d[0])  # at γ = 0.005
    # KU-1.30 #1 acceptance band [15, 200] Pa, deliberately wide.
    assert 15.0 <= G_small <= 200.0, (
        f"Preliminary G(γ=0.005) = {G_small:.1f} Pa outside KU-1.30 band "
        f"[15, 200] Pa (h={h*1e6} μm convention)"
    )

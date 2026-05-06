"""Stage 1a physics — fast conservation checks (small N, few steps).

These run on whatever Taichi backend is available (cuda preferred, cpu fallback).
They are deliberately tiny so they fit in the existing pytest smoke suite.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

ti = pytest.importorskip("taichi")


@pytest.fixture(scope="module")
def solver():
    from acs.gpu import init_taichi
    from acs.physics.mlsmpm import MLSMPMSolver, SolverConfig

    backend = os.environ.get("ACS_GPU_BACKEND") or "auto"
    init_taichi(backend=backend)

    cfg = SolverConfig(
        n_particles=200,
        grid_n=32,
        domain_star=4.0,
        radius_star=1.0,
        dt_star=0.01,
        K_star=1.0,
        mu_star=0.3,
        tau_star=1.0,
        capillary_number=0.01,
        drag_xi_star=1.0,
        density_star=1.0,
        free_surface_threshold=0.6,
        seed=7,
    )
    s = MLSMPMSolver(cfg)
    s.initialize_sphere(np.full(3, cfg.domain_star * 0.5, dtype=np.float32))
    yield s
    ti.reset()


def test_initial_state_is_finite(solver):
    inv = solver.invariants()
    assert inv["nan_count"] == 0
    assert inv["mass_star"] > 0
    assert np.all(np.isfinite(inv["momentum_star"]))


def test_mass_exactly_conserved(solver):
    inv0 = solver.invariants()
    for _ in range(20):
        solver.step()
    inv1 = solver.invariants()
    drift = abs(inv1["mass_star"] - inv0["mass_star"]) / inv0["mass_star"]
    # Mass is set once at construction; drift here detects writes-to-mass bugs.
    assert drift < 1e-12, f"mass drift {drift:.2e} should be at round-off"


def test_no_nan_after_many_steps(solver):
    for _ in range(100):
        solver.step()
    inv = solver.invariants()
    assert inv["nan_count"] == 0
    assert np.all(np.isfinite(inv["momentum_star"]))
    assert np.isfinite(inv["kinetic_energy_star"])
    assert np.isfinite(inv["strain_energy_star"])


def test_solver_rejects_too_few_particles():
    from acs.physics.mlsmpm import MLSMPMSolver, SolverConfig

    cfg = SolverConfig(
        n_particles=4, grid_n=32, domain_star=4.0, radius_star=1.0, dt_star=0.01,
        K_star=1.0, mu_star=0.3, tau_star=1.0, capillary_number=0.01,
        drag_xi_star=1.0, density_star=1.0, free_surface_threshold=0.6,
    )
    with pytest.raises(ValueError):
        MLSMPMSolver(cfg)

"""KU-2.4 / KU-2.8 biphasic regression test.

The Chan-Odde model predicts that mean traction force is a non-monotonic
function of substrate stiffness, peaked at::

    k_sub* ≈ N_m · F_stall / (v_unloaded · τ_max)

where ``τ_max`` is the maximum bond lifetime from the Pereverzev catch
peak (KU-2.5). Converting to Young's modulus through the substrate stub::

    E* = k_sub* / (π a) · (1 − ν²)

This test sweeps ``E`` across 5 decades and confirms the location of the
biphasic peak agrees with the closed-form ``E*`` within a factor of 2,
per the Brief Task 7 acceptance criterion.

The test is marked ``slow`` so the default ``pytest`` invocation does
not run it. Use ``pytest -m slow`` or call the module directly.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from acs_kb.bridge.motor_clutch import (
    MotorClutchFA,
    MotorClutchParams,
    run_steady_state,
)
from acs_kb.bridge.substrate_stub import LinearElasticSubstrate
from acs_kb.bridge.types import make_focal_adhesion
from acs_kb.common.derived_params import load_bridge_config

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_unit2_1.yaml"


def _mean_traction_at_E(
    E: float,
    *,
    nu: float,
    a: float,
    dt: float,
    n_steps: int,
    rng: np.random.Generator,
) -> float:
    """Run a single FA for ``n_steps`` at substrate Young's modulus ``E``."""
    sub = LinearElasticSubstrate(young_modulus=E, poisson_ratio=nu,
                                 contact_radius=a)
    fa = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=50)
    mc = MotorClutchFA(fa, sub, MotorClutchParams(n_clutches=50))
    result = run_steady_state(mc, dt=dt, n_steps=n_steps, rng=rng,
                              burn_in_fraction=0.5)
    return result["mean_force_total"]


@pytest.mark.slow
def test_KU_2_4_biphasic_peak_within_factor_2():
    """KU-2.4 / KU-2.8: argmax E from a 5-decade sweep ≈ closed-form E*."""
    cfg = load_bridge_config(CONFIG_PATH)
    b = cfg["bridge"]
    sw = b["biphasic_sweep"]
    factor = float(b["acceptance"]["biphasic_peak_factor"])
    E_star_analytic = float(b["derived"]["biphasic_EStar"])

    E_grid = np.logspace(sw["log10_E_min"], sw["log10_E_max"], int(sw["n_points"]))
    mean_F = np.empty_like(E_grid)
    seed_root = int(sw["seed"])
    for i, E in enumerate(E_grid):
        rng = np.random.default_rng(seed_root + i)
        mean_F[i] = _mean_traction_at_E(
            float(E), nu=b["substrate"]["poisson_ratio"],
            a=b["substrate"]["contact_radius"],
            dt=b["dynamics"]["dt"], n_steps=int(sw["n_steps_per_E"]),
            rng=rng,
        )

    E_argmax = float(E_grid[int(np.argmax(mean_F))])
    ratio = E_argmax / E_star_analytic
    assert (1.0 / factor) <= ratio <= factor, (
        f"E_argmax={E_argmax:.3e} Pa, E*_analytic={E_star_analytic:.3e} Pa, "
        f"ratio={ratio:.3f}, factor={factor}. "
        f"Sweep: E_grid={E_grid.tolist()}, mean_F={mean_F.tolist()}"
    )
    # Biphasic shape sanity: peak should NOT be at the boundary of the sweep.
    assert 0 < int(np.argmax(mean_F)) < len(E_grid) - 1, (
        "peak landed at sweep boundary — widen the log10_E range"
    )

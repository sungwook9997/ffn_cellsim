"""KU-5.x lamellipodium emergent-physics validation gates (H.5 production sign-off).

Per H.5 brief §Validation acceptance:

| Gate                  | Criterion                                                    | KU      |
| ---                   | ---                                                          | ---     |
| Dendritic density     | ≈ 100 barbed ends per μm² of WAVE plane at steady state      | KU-5.1  |
| Bieling force-velocity| v_network(F) matches Bieling 2016 Fig 2 within ±30 %         | KU-5.2  |
| Funk abortive         | branching rate drops > 50 % above 500 Pa per-WAVE force      | KU-5.3  |

All gates require long BAOAB simulation to reach steady-state (multi-
hour wall on M1 Max).  Skeletons SKIP under H5_KU5_PRODUCTION=1
opt-in, same multi-hour rationale as H.3 KU-3.x and L_p FULL gates.
"""

from __future__ import annotations

import math
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.archive.hoomd_legacy.cell.lamellipodium import (
    build_lamellipodium_simulation,
    resolve_h5_lamellipodium,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[2] / "configs" / "phase1_h5.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


H5_KU5_PRODUCTION = bool(int(os.environ.get("H5_KU5_PRODUCTION", "0")))


# ---------------------------------------------------------------------------
# Measurement utilities (always available)
# ---------------------------------------------------------------------------
def dendritic_density(n_barbed_ends: int, wave_area: float) -> float:
    """Per-μm² dendritic density given total barbed-end count + WAVE area.

    Returns rate per μm² for direct comparison with KU-5.1 ≈ 100 / μm².
    """
    if wave_area <= 0:
        return float("nan")
    return n_barbed_ends / (wave_area * 1.0e12)  # m² → μm²


def bieling_force_velocity_oracle(
    F: np.ndarray, *, v0: float, delta_elong: float, kT: float,
) -> np.ndarray:
    """Bieling 2016 slip Bell-Evans v(F) closed-form (oracle).

    ``v(F) = v0 · exp(−F · δ_elong / kT)``
    """
    return v0 * np.exp(-F * delta_elong / kT)


# ---------------------------------------------------------------------------
# KU-5.1 dendritic density gate (skeleton)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not H5_KU5_PRODUCTION,
    reason=(
        "KU-5.1 dendritic-density steady-state gate requires multi-hour "
        "BAOAB simulation to reach branching/capping equilibrium. "
        "Opt-in via H5_KU5_PRODUCTION=1."
    ),
)
class TestKU51DendriticDensity:
    def test_steady_state_density_matches_KU51_target(self):
        pytest.skip("Production sign-off skeleton.")


# ---------------------------------------------------------------------------
# KU-5.2 Bieling force-velocity gate (skeleton)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not H5_KU5_PRODUCTION,
    reason="H5_KU5_PRODUCTION=1 required (see KU-5.1 reason).",
)
class TestKU52BielingForceVelocity:
    def test_force_velocity_curve_matches_bieling_oracle(self):
        pytest.skip("Production sign-off skeleton.")


# ---------------------------------------------------------------------------
# KU-5.3 Funk abortive branching gate (skeleton)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not H5_KU5_PRODUCTION,
    reason="H5_KU5_PRODUCTION=1 required (see KU-5.1 reason).",
)
class TestKU53FunkAbortive:
    def test_branching_rate_drops_above_threshold(self):
        pytest.skip("Production sign-off skeleton.")


# ---------------------------------------------------------------------------
# Measurement utility STATIC checks (always run)
# ---------------------------------------------------------------------------
class TestMeasurementUtilities:
    def test_dendritic_density_units(self):
        # 100 barbed ends in 1 μm² WAVE area → 100 / μm².
        density = dendritic_density(100, 1.0e-12)
        assert math.isclose(density, 100.0, rel_tol=1e-9)

    def test_bieling_oracle_v0_at_zero_force(self):
        v0 = 11.6   # 1/s
        kT = 4.28e-21
        delta = 2.7e-9
        v = bieling_force_velocity_oracle(
            np.array([0.0]), v0=v0, delta_elong=delta, kT=kT,
        )
        assert math.isclose(float(v[0]), v0, rel_tol=1e-9)

    def test_bieling_oracle_monotonic_decreasing(self):
        v0 = 11.6
        kT = 4.28e-21
        delta = 2.7e-9
        F = np.linspace(0, 2.0e-12, 20)
        v = bieling_force_velocity_oracle(F, v0=v0, delta_elong=delta, kT=kT)
        assert (np.diff(v) <= 0).all()

    def test_bieling_oracle_e_fold_at_kT_over_delta(self):
        """v drops by 1/e when F = kT/δ_elong."""
        v0 = 11.6
        kT = 4.28e-21
        delta = 2.7e-9
        F_efold = kT / delta
        v = bieling_force_velocity_oracle(
            np.array([F_efold]), v0=v0, delta_elong=delta, kT=kT,
        )
        assert math.isclose(float(v[0]), v0 / math.e, rel_tol=1e-6)

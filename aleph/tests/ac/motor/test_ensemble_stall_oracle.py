"""Self-test of the ensemble-stall EMERGENCE oracle (pure NumPy — no Warp/CUDA).

The gate the build plan names explicitly: ensemble stall EMERGES from bound heads, it is NOT the imposed
product N_side * F_head (I3, §3 line 243).
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.motor.bell_kinetics_analytic import bell_f0_from_x_beta
from aleph.components.motor.ensemble_stall_analytic import (
    ensemble_stall_meanfield,
    ensemble_stall_stochastic,
)

# illustrative (GAP/provisional) params — structural gates are magnitude-independent (params_i0b3.yaml)
N_SIDE = 28       # GAP: 10 (AFINES) vs 28-30 (Billington)
F_HEAD = 2.0      # pN, GAP: 0.5 vs 2.0
K_ON = 50.0
K_OFF0 = 0.35
F0 = bell_f0_from_x_beta(0.6e-3)


def test_meanfield_is_emergent_product_not_naive() -> None:
    """F_ensemble = N_side * phi_b * F_head with phi_b computed from kinetics, and F_ensemble <= N_side*F_head."""
    r = ensemble_stall_meanfield(N_SIDE, F_HEAD, K_ON, K_OFF0, F0)
    assert 0.0 < r["phi_b"] <= 1.0
    assert r["f_ensemble"] == pytest.approx(N_SIDE * r["phi_b"] * F_HEAD, rel=1e-12)
    assert r["f_naive"] == pytest.approx(N_SIDE * F_HEAD, rel=1e-12)
    assert r["f_ensemble"] <= r["f_naive"] + 1e-12       # strict bound (equality only if phi_b == 1)
    assert r["f_ensemble"] < r["f_naive"]                # phi_b < 1 here => strictly below the naive product


def test_equality_only_when_all_heads_bound() -> None:
    """As k_off0 -> 0 (or k_on -> inf) phi_b -> 1 and the emergent stall approaches the naive product."""
    r = ensemble_stall_meanfield(N_SIDE, F_HEAD, K_ON, k_off0=1e-9, f0=F0)
    assert r["phi_b"] == pytest.approx(1.0, abs=1e-6)
    assert r["f_ensemble"] == pytest.approx(r["f_naive"], rel=1e-5)


def test_stochastic_mean_converges_to_meanfield() -> None:
    """The per-realisation stochastic ensemble mean converges to the mean-field value (heads are summed)."""
    rng = np.random.default_rng(44)
    r = ensemble_stall_meanfield(N_SIDE, F_HEAD, K_ON, K_OFF0, F0)
    draws = ensemble_stall_stochastic(N_SIDE, F_HEAD, K_ON, K_OFF0, F0, n_realizations=200_000, rng=rng)
    sem = draws.std() / np.sqrt(draws.size)
    assert abs(draws.mean() - r["f_ensemble"]) < 5.0 * sem
    # per-realisation spread is real (binomial), not a delta at the mean
    assert draws.std() > 0.0


def test_self_limiting_monotone_in_off_rate() -> None:
    """Higher zero-force off-rate (or lower Bell force f0) => lower engaged fraction => lower ensemble stall."""
    r_lo = ensemble_stall_meanfield(N_SIDE, F_HEAD, K_ON, k_off0=0.35, f0=F0)
    r_hi = ensemble_stall_meanfield(N_SIDE, F_HEAD, K_ON, k_off0=20.0, f0=F0)
    assert r_hi["f_ensemble"] < r_lo["f_ensemble"]


def test_not_imposed_arbiter_kinetics_move_the_force() -> None:
    """Changing the binding kinetics changes F_ensemble — a hard-coded N_side*F_head would not move."""
    r_a = ensemble_stall_meanfield(N_SIDE, F_HEAD, k_on=50.0, k_off0=K_OFF0, f0=F0)
    r_b = ensemble_stall_meanfield(N_SIDE, F_HEAD, k_on=2.0, k_off0=K_OFF0, f0=F0)
    assert r_a["f_ensemble"] != pytest.approx(r_b["f_ensemble"])
    assert r_b["f_ensemble"] < r_a["f_ensemble"]         # weaker binding => fewer engaged => less force

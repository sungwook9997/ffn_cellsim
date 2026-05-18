"""Phase 1 Unit 4.1 two-cell pair integration test (brief Task 4).

Acceptance (from the Worker D Unit 4.1 brief):

1. Contact angle ≈ π/2 ± 0.2 rad (KU-4.4 nominal symmetric pair).
2. Junction force ``sum(bond_forces) ∈ [1, 10] nN``.
3. ``n_bonds_engaged ≥ 30`` (mature regime).
4. 60 s simulated time runs in < 30 s wall-clock.

Parameters are loaded from ``acs_kb/configs/phase1_unit4_1.yaml`` and
follow the documented Phase 1 calibration (KU-4.17 with the Bell 1978
lower-bound Δx* override; see ``acs_kb/configs/phase1_unit4_1.yaml`` and
``acs_kb/outputs/phase1/unit4_1/REPORT.md`` for the full Sanity-Gate
rationale).
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pytest
import yaml

from acs_kb.cell.cell import Cell
from acs_kb.cell.cortex import generate_cortex
from acs_kb.junction.cadherin import update_bonds
from acs_kb.junction.contact_angle import compute_contact_angle
from acs_kb.junction.types import make_ecadherin_junction


_CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_unit4_1.yaml"
)


@pytest.fixture(scope="module")
def cfg() -> dict:
    with _CONFIG_PATH.open("r") as fh:
        return yaml.safe_load(fh)["junction"]


def _build_pair(cfg: dict):
    pt = cfg["pair_test"]
    R = float(pt["R_cell"])
    sep = float(pt["cell_separation"])
    params_cortex = {"gamma_cortex": float(pt["gamma_cortex"])}
    ctx_a = generate_cortex(
        R_cell=R,
        n_cortex_fibers=int(pt["n_cortex_fibers"]),
        L_cortex_fiber=float(pt["L_cortex_fiber"]),
        beads_per_fiber=int(pt["beads_per_fiber"]),
        cell_center=np.array([-0.5 * sep, 0.0]),
        seed=int(pt["seed_cell_a"]),
        params=params_cortex,
    )
    ctx_b = generate_cortex(
        R_cell=R,
        n_cortex_fibers=int(pt["n_cortex_fibers"]),
        L_cortex_fiber=float(pt["L_cortex_fiber"]),
        beads_per_fiber=int(pt["beads_per_fiber"]),
        cell_center=np.array([+0.5 * sep, 0.0]),
        seed=int(pt["seed_cell_b"]),
        params=params_cortex,
    )
    cell_a = Cell.from_cortex(0, ctx_a)
    cell_b = Cell.from_cortex(1, ctx_b)
    junction = make_ecadherin_junction(
        cell_a, cell_b, n_bonds_total=int(cfg["n_bonds_total"]),
    )
    return cell_a, cell_b, junction


def test_two_cell_pair_steady_state_matches_acceptance(cfg: dict) -> None:
    """End-to-end check of the four brief Task-4 acceptance criteria."""
    pt = cfg["pair_test"]
    acc = cfg["acceptance"]

    cell_a, cell_b, junction = _build_pair(cfg)
    rng = np.random.default_rng(int(pt["seed_rng"]))

    dt = float(pt["dt"])
    n_steps = int(pt["n_steps"])
    F_target = float(pt["F_total_target"])
    dx_star = float(cfg["dx_star"])
    k_off0 = float(cfg["k_off0"])
    k_on = float(cfg["k_on"])
    kT = float(cfg["kT"])
    gamma_J_override = float(pt["gamma_J_override"])

    n_history = []
    sum_force_history = []

    wall_start = time.perf_counter()
    # Drive the junction at a constant external load F_target. When
    # n_bonds_engaged > 0, this load is shared equally across engaged
    # bonds (mean field); when n_bonds_engaged = 0, the force is
    # "applied" but transmits nothing — we still call update_bonds so
    # the on-rate arm can re-engage bonds.
    for _ in range(n_steps):
        update_bonds(
            junction,
            F_total=F_target,
            dt=dt,
            rng=rng,
            k_off0=k_off0, dx_star=dx_star, k_on=k_on, kT=kT,
        )
        n_history.append(junction.n_bonds_engaged)
        sum_force_history.append(float(junction.bond_forces.sum()))
    wall_elapsed = time.perf_counter() - wall_start

    # Compute the contact angle once at the end of the simulation
    # using the γ_J override (the Young equation is tested as a
    # geometric balance, decoupled from the bond population — see
    # contact_angle.py docstring).
    theta = compute_contact_angle(
        cell_a, cell_b, junction,
        params={
            "gamma_c_a": float(pt["gamma_cortex"]),
            "gamma_c_b": float(pt["gamma_cortex"]),
            "gamma_J": gamma_J_override,
        },
    )
    junction.contact_angle = theta

    # ---- Acceptance band checks ----

    # 1. Contact angle ≈ π/2 ± 0.2 rad.
    theta_target = float(acc["contact_angle_target"])
    theta_tol = float(acc["contact_angle_tol"])
    assert abs(theta - theta_target) <= theta_tol, (
        f"contact angle θ = {theta:.4f} rad ({np.degrees(theta):.1f}°) "
        f"out of {theta_target:.4f} ± {theta_tol}"
    )

    # 2. Junction force ∈ [1, 10] nN — averaged over the last 20 s of
    #    the simulation so transient fluctuations do not dominate.
    tail = max(1, int(20.0 / dt))
    mean_force = float(np.mean(sum_force_history[-tail:]))
    F_min = float(acc["F_total_min"])
    F_max = float(acc["F_total_max"])
    assert F_min <= mean_force <= F_max, (
        f"mean junction force in last 20 s = {mean_force * 1e9:.3f} nN "
        f"out of [{F_min*1e9:.1f}, {F_max*1e9:.1f}] nN"
    )

    # 3. n_engaged ≥ 30 in the last 20 s.
    mean_n = float(np.mean(n_history[-tail:]))
    n_min = int(acc["n_engaged_min"])
    assert mean_n >= n_min, (
        f"mean n_engaged in last 20 s = {mean_n:.1f} below {n_min}"
    )

    # 4. Wall-clock < 30 s for 60 s of simulated time.
    wall_max = float(acc["pair_sim_wall_s_max"])
    assert wall_elapsed <= wall_max, (
        f"60 s of simulated time took {wall_elapsed:.2f} s wall "
        f"(budget {wall_max:.1f} s)"
    )


def test_two_cell_pair_contact_points_lie_inward(cfg: dict) -> None:
    """The contact points of the freshly built junction face each other.

    A symmetric pair with cell A on the left and cell B on the right
    should have ``contact_position_a.x > cell_a.center_position.x`` and
    ``contact_position_b.x < cell_b.center_position.x`` — i.e. each
    contact point lies on the inward-facing cortex arc. This is a
    geometric sanity check for the Cell.compute_cortex_boundary_position
    call inside ``make_ecadherin_junction``.
    """
    cell_a, cell_b, junction = _build_pair(cfg)
    assert junction.contact_position_a[0] > cell_a.center_position[0]
    assert junction.contact_position_b[0] < cell_b.center_position[0]


def test_two_cell_pair_canonical_id_order(cfg: dict) -> None:
    """Junction always orders ``cell_a_id < cell_b_id`` regardless of input order."""
    cell_a, cell_b, _ = _build_pair(cfg)
    # Build the junction in reversed order; the factory should swap.
    j_rev = make_ecadherin_junction(cell_b, cell_a, n_bonds_total=int(cfg["n_bonds_total"]))
    assert j_rev.cell_a_id == min(cell_a.id, cell_b.id)
    assert j_rev.cell_b_id == max(cell_a.id, cell_b.id)

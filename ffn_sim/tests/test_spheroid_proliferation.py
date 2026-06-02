"""Tests for ffn_sim.spheroid.proliferation — L2.4 contact-inhibited division.

Pure-geometry tests (no HOOMD): the contact-shell count, free-direction bud, the division
gate (timer + contact inhibition), rim localisation, and the resolver derivations. Plus two
small HOOMD growth runs: the dilute-limit doubling (exponential) boundary and the
proliferation-OFF reduction to the G1 stable aggregate (count conserved).

These exercise the G4 acceptance bands (oracle config ``layer2_cbm.yaml`` g4): rim-localised
division, dilute doubling, sub-exponential populated growth.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.spheroid.params import (
    ResolvedProliferation,
    resolve_layer2,
    resolve_proliferation,
)
from ffn_sim.spheroid.proliferation import (
    apply_divisions,
    best_bud_direction,
    build_pool_simulation,
    candidate_directions,
    first_shell_counts,
    run_growth,
    run_growth_pooled,
    sample_cycle_targets,
)

_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


@pytest.fixture(scope="module")
def resolved():
    return resolve_layer2(yaml.safe_load(_CONFIG.read_text()))


@pytest.fixture(scope="module")
def prolif(resolved):
    return resolve_proliferation(yaml.safe_load(_CONFIG.read_text()), resolved)


# --------------------------------------------------------------------------- resolver

def test_resolve_proliferation_derivations(resolved, prolif):
    """Length scales derive from r0; kissing number is the geometric Z=12."""
    assert prolif.kissing_number == 12
    assert prolif.shell_cutoff == pytest.approx(prolif.shell_factor * resolved.morse_r0)
    assert prolif.split_distance == pytest.approx(prolif.split_factor * resolved.morse_r0)
    assert prolif.split_distance == pytest.approx(resolved.morse_r0)  # split_factor = 1
    assert prolif.cycle_time_mean == pytest.approx(1.08e5)  # 30 h MCF7 doubling


def test_resolve_proliferation_rejects_bad_cv(resolved):
    cfg = yaml.safe_load(_CONFIG.read_text())
    cfg["spheroid"]["proliferation"]["cycle_time_cv"] = 1.5  # out of [0,1)
    with pytest.raises(ValueError):
        resolve_proliferation(cfg, resolved)


def test_resolve_proliferation_rejects_shell_below_rest(resolved):
    cfg = yaml.safe_load(_CONFIG.read_text())
    cfg["spheroid"]["proliferation"]["shell_factor"] = 0.5  # below the rest separation
    with pytest.raises(ValueError):
        resolve_proliferation(cfg, resolved)


# --------------------------------------------------------------------------- shell counts

def test_first_shell_counts_octahedron():
    """A central cell with 6 neighbours at the rest separation has count 6 (the rest farther)."""
    r0 = 1.0
    pos = np.array(
        [[0, 0, 0], [r0, 0, 0], [-r0, 0, 0], [0, r0, 0], [0, -r0, 0], [0, 0, r0], [0, 0, -r0]],
        dtype=float,
    )
    counts = first_shell_counts(pos, shell_cutoff=1.25 * r0)
    assert counts[0] == 6  # central cell sees all 6
    assert counts[1] == 1  # an arm sees only the centre (others are >1.25 r0 away)


def test_first_shell_counts_single_cell():
    assert first_shell_counts(np.zeros((1, 3)), shell_cutoff=1.0).tolist() == [0]


# --------------------------------------------------------------------------- bud direction

def test_candidate_directions_unit_norm():
    d = candidate_directions(12)
    assert d.shape == (12, 3)
    assert np.allclose(np.linalg.norm(d, axis=1), 1.0)


def test_best_bud_direction_points_into_free_space():
    """Neighbours clustered in -x ⇒ the daughter buds toward +x (away from the crowd)."""
    r0 = 1.0
    pos = np.array([[0, 0, 0], [-r0, 0, 0], [-r0, r0, 0], [-r0, -r0, 0]], dtype=float)
    nbr = np.array([1, 2, 3])
    nhat, gap = best_bud_direction(pos, 0, nbr, split=r0)
    assert nhat[0] > 0.0  # x-component positive => budding into the empty +x half-space
    assert gap > 0.0


def test_best_bud_direction_no_neighbours_is_valid_unit():
    nhat, gap = best_bud_direction(np.zeros((1, 3)), 0, np.array([], dtype=int), split=1.0)
    assert np.linalg.norm(nhat) == pytest.approx(1.0)
    assert gap == float("inf")


# --------------------------------------------------------------------------- division gate

def _prolif(cycle_time=1.0e5, cv=0.0, r0=1.0):
    return ResolvedProliferation(
        cycle_time_mean=cycle_time, cycle_time_cv=cv, kissing_number=12,
        shell_factor=1.25, split_factor=1.0, shell_cutoff=1.25 * r0, split_distance=r0,
        min_gap=0.9 * r0,  # Morse repulsive-core onset (r0 − 0.1·r0 contact zone)
    )


def test_due_isolated_cell_divides():
    """An isolated, timer-elapsed cell divides → count +1, daughters at the rest separation."""
    p = _prolif()
    rng = np.random.default_rng(0)
    pos = np.zeros((1, 3))
    out = apply_divisions(pos, np.array([2.0e5]), np.array([1.0e5]), p, rng)
    assert out["positions"].shape[0] == 2
    sep = np.linalg.norm(out["positions"][0] - out["positions"][1])
    assert sep == pytest.approx(p.split_distance, rel=1e-9)  # born at rest separation r0
    assert out["divided_idx"].tolist() == [0]


def test_not_due_cell_does_not_divide():
    p = _prolif()
    rng = np.random.default_rng(0)
    out = apply_divisions(np.zeros((1, 3)), np.array([0.5e5]), np.array([1.0e5]), p, rng)
    assert out["positions"].shape[0] == 1  # age < target → no division
    assert out["divided_idx"].size == 0


def test_contact_inhibition_blocks_buried_cell():
    """A cell with a FULL first shell (>= kissing number) is quiescent even when due."""
    p = _prolif()
    rng = np.random.default_rng(0)
    # central cell surrounded by 12 neighbours on a shell at the rest separation
    dirs = candidate_directions(12)
    pos = np.vstack([np.zeros(3), dirs * 1.0])
    ages = np.full(pos.shape[0], 2.0e5)
    targets = np.full(pos.shape[0], 1.0e5)
    out = apply_divisions(pos, ages, targets, p, rng)
    # the central cell (index 0, 12 contacts) must NOT be among the divided
    assert 0 not in out["divided_idx"].tolist()


def test_division_is_rim_localized(resolved):
    """In a REAL settled aggregate, the cells that divide are biased to the rim (G4).

    A fixed neighbour-count threshold mis-classifies here (the settled liquid-like bulk
    coordination is only ≈11); the free-space gate is what produces rim localisation. We
    therefore test on a genuinely-settled blob, not a sparse synthetic cloud.
    """
    from ffn_sim.spheroid.cbm import build_cbm_simulation, get_positions

    p = resolve_proliferation(yaml.safe_load(_CONFIG.read_text()), resolved)
    sim, _a, _u, _rc = build_cbm_simulation(resolved, 250, seed=1)
    sim.run(0)
    sim.run(20_000)
    pos = get_positions(sim)
    ages = np.full(pos.shape[0], 2.0 * p.cycle_time_mean)  # all due
    targets = np.full(pos.shape[0], p.cycle_time_mean)
    out = apply_divisions(pos, ages, targets, p, np.random.default_rng(1))
    assert out["divided_idx"].size > 0
    rim_frac = float(np.mean(out["divided_radial"] >= out["median_radial"]))
    assert rim_frac >= 0.70  # G4 rim_division_fraction_min


def test_sample_cycle_targets_positive_and_mean():
    p = _prolif(cycle_time=1.0e5, cv=0.15)
    rng = np.random.default_rng(0)
    t = sample_cycle_targets(5000, p, rng)
    assert (t > 0).all()
    assert t.mean() == pytest.approx(1.0e5, rel=0.02)


# --------------------------------------------------------------------------- HOOMD growth

@pytest.mark.parametrize("seed", [7])
def test_dilute_limit_doubles_exponentially(resolved, seed):
    """A single founder cell (no crowding) doubles repeatedly → exponential, not linear."""
    p = _prolif(cycle_time=3.0e3, cv=0.0, r0=resolved.morse_r0)
    res = run_growth(
        resolved, p, n_cells_init=1, total_time=1.2e4, epoch_steps=300,
        settle_steps=100, seed=seed,
    )
    # ~4 cycles of head-room; from a random phase expect >= 3 doublings (>=8 cells), well
    # above any linear/contact-limited growth, and not above the exponential ceiling.
    t_final = float(res["t"][-1])
    ceiling = 2 ** (t_final / p.cycle_time_mean + 1.0)  # +1 cycle slack
    assert res["n_cells"][-1] >= 6
    assert res["n_cells"][-1] <= ceiling


def test_proliferation_off_reduces_to_stable_aggregate(resolved):
    """Effectively-infinite cycle time ⇒ no division ⇒ count conserved (the G1 limit)."""
    p = _prolif(cycle_time=1.0e15, cv=0.0, r0=resolved.morse_r0)
    res = run_growth(
        resolved, p, n_cells_init=60, total_time=6.0e3, epoch_steps=500,
        settle_steps=200, seed=3,
    )
    assert res["growth_factor"] == pytest.approx(1.0)
    assert res["n_cells"][-1] == 60
    assert res["n_division_epochs"] == 0


# ----------------------------------------------------------------- pooled (leak-free, L2.4b)

def test_pooled_build_pool_types_and_active_mask(resolved):
    """The pool has n_max particles: n_active 'cell' + the rest parked 'void'."""
    import hoomd

    dev = hoomd.device.CPU(notice_level=0)
    sim, active, _rc = build_pool_simulation(resolved, n_active_init=30, n_max=200, device=dev, seed=1)
    assert sim.state.N_particles == 200       # fixed-N pool (no rebuild ever)
    assert active.sum() == 30                  # 30 active cells
    snap = sim.state.get_snapshot()
    typeid = np.asarray(snap.particles.typeid)
    assert (typeid == 1).sum() == 170          # 170 parked voids
    del sim


def test_pooled_growth_grows_and_is_bounded(resolved):
    """Pooled growth increases N (proliferation), never exceeds the pool, stays rim-localised."""
    p = _prolif(cycle_time=4.0e3, cv=0.0, r0=resolved.morse_r0)
    res = run_growth_pooled(
        resolved, p, n_cells_init=120, total_time=1.2e4, epoch_steps=400,
        settle_steps=300, max_cells=600, seed=7,
    )
    assert res["n_cells"][-1] > res["n_cells"][0]          # it grew
    assert res["n_cells"][-1] <= 600                        # never exceeds the pool
    assert res["growth_factor"] < 2.0 ** (1.2e4 / p.cycle_time_mean) + 0.5  # sub-exponential


def test_pooled_no_growth_when_cycle_infinite(resolved):
    """Effectively-infinite cycle ⇒ no division ⇒ active count conserved (pooled G1 limit)."""
    p = _prolif(cycle_time=1.0e15, cv=0.0, r0=resolved.morse_r0)
    res = run_growth_pooled(
        resolved, p, n_cells_init=80, total_time=4.0e3, epoch_steps=500,
        settle_steps=200, max_cells=400, seed=3,
    )
    assert res["n_cells"][-1] == 80
    assert res["growth_factor"] == pytest.approx(1.0)

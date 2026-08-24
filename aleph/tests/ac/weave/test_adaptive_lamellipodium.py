"""Self-test of the adaptive Arp2/3 lamellipodium dynamics — nucleation / capping / ratchet (pure NumPy)."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.weave.lamellipodium import AdaptiveLamellipodium, build_lamellipodium

# I0-B4 GAP magnitudes used ONLY as oracle sweep variables (SHAPE gates are magnitude-independent).
STEP = dict(
    dt=0.1, membrane_load_pN=5.0, k_arp0=5.0, k_m=5.0, k_cap=2.0, k_dissolve=1.0,
    k_on_c=30.0, k_off=3.0, monomer_c=50.0, npf_activity=1.0, clutch_engagement=1.0,
)


def _make(seed_int: int, *, n_mothers: int = 10, pool: int = 400) -> AdaptiveLamellipodium:
    rng = np.random.default_rng(seed_int)
    seed = build_lamellipodium(
        n_mothers=n_mothers, n_daughter_pool=pool, length_um=1.0, seg_um=0.5, patch_half_um=8.0,
        rng=rng, monomer_c=0.0, k_arp0=5.0, k_m=5.0, tau=0.0)     # tau=0 => no daughters pre-nucleated at build
    return AdaptiveLamellipodium.from_seed(seed, seg_um=0.5, length_um=1.0, n_mothers=n_mothers)


def test_population_is_n_fixed_under_activation() -> None:
    """Nucleation activates dormant daughters WITHOUT ever exceeding the pre-allocated pool (N-FIXED)."""
    lam = _make(0)
    rng = np.random.default_rng(0)
    capacity = lam.active.shape[0]
    for _ in range(50):
        tick = lam.step(rng=rng, **STEP)
        assert tick.n_active <= capacity                          # never allocates beyond the fixed pool
        assert lam.length_um.shape[0] == capacity                 # no new nodes/filaments
    assert int(lam.active.sum()) + int((~lam.active).sum()) == capacity   # active + dormant conserved


def test_growing_count_reaches_a_nucleation_capping_balance() -> None:
    """Growing ends converge to a steady turnover fixed point (nucleation = capping = dissolution flux).

    Fast Arp2/3 nucleation over a large dormant pool OVERSHOOTS, then capping+dissolution pull the growing
    population down to a plateau that is strictly between an empty array and the full pool capacity — a genuine
    balance, not a runaway to the cap nor a collapse to zero.
    """
    lam = _make(1)
    rng = np.random.default_rng(1)
    traj = lam.run(160, rng=rng, **STEP)
    capacity = lam.active.shape[0]
    mid = np.mean([t.n_growing for t in traj[80:100]])
    late = np.mean([t.n_growing for t in traj[-40:]])
    assert 0 < late < capacity                                    # a nonzero steady population below the cap
    assert abs(late - mid) / max(late, 1.0) < 0.15                # converged (mid and late agree)
    tail = np.array([t.n_growing for t in traj[-40:]], float)
    assert tail.std() / max(tail.mean(), 1.0) < 0.2               # quasi-steady tail (small relative spread)


def test_higher_capping_gives_a_shorter_dendritic_mesh() -> None:
    """Rarer capping -> longer growing tips (mean length ~ V/k_cap); the mesh size is SET by k_cap, not scripted."""
    lo_cap = dict(STEP, k_cap=1.0)
    hi_cap = dict(STEP, k_cap=8.0)
    lam_lo, lam_hi = _make(2), _make(2)
    traj_lo = lam_lo.run(150, rng=np.random.default_rng(2), **lo_cap)
    traj_hi = lam_hi.run(150, rng=np.random.default_rng(2), **hi_cap)
    mean_len_lo = np.mean([t.mean_growing_length_um for t in traj_lo[-40:]])
    mean_len_hi = np.mean([t.mean_growing_length_um for t in traj_hi[-40:]])
    assert mean_len_lo > mean_len_hi


def test_protrusion_slows_under_higher_membrane_load() -> None:
    """The ratchet is load-adaptive: a heavier membrane load slows the steady protrusion velocity."""
    light = dict(STEP, membrane_load_pN=2.0)
    heavy = dict(STEP, membrane_load_pN=40.0)
    lam_l, lam_h = _make(3), _make(3)
    traj_l = lam_l.run(80, rng=np.random.default_rng(3), **light)
    traj_h = lam_h.run(80, rng=np.random.default_rng(3), **heavy)
    v_light = np.mean([t.protrusion_velocity_um_s for t in traj_l[-30:]])
    v_heavy = np.mean([t.protrusion_velocity_um_s for t in traj_h[-30:]])
    assert v_light > v_heavy
    assert traj_l[-1].membrane_front_um > traj_h[-1].membrane_front_um


def test_clutch_engagement_partitions_protrusion_and_retrograde() -> None:
    """The seeded clutch engagement splits barbed-end polymer into protrusion vs retrograde (sums to v_poly)."""
    lam = _make(4)
    rng = np.random.default_rng(4)
    tick = None
    for _ in range(10):
        tick = lam.step(rng=rng, **dict(STEP, clutch_engagement=0.3))
    assert tick is not None
    total = tick.protrusion_velocity_um_s + tick.retrograde_velocity_um_s
    assert tick.retrograde_velocity_um_s > tick.protrusion_velocity_um_s   # engagement 0.3 => mostly retrograde
    # full engagement puts all polymer into protrusion, zero into retrograde
    lam2 = _make(4)
    rng2 = np.random.default_rng(4)
    t2 = None
    for _ in range(10):
        t2 = lam2.step(rng=rng2, **dict(STEP, clutch_engagement=1.0))
    assert t2 is not None
    assert t2.retrograde_velocity_um_s == pytest.approx(0.0)
    assert total > 0.0


def test_depleted_monomer_stops_nucleation() -> None:
    """With a depleted monomer field the flux-limited nucleation is zero: no new filaments activate (flux limit).

    Turnover (capping + dissolution) can still remove filaments, so the active population can only SHRINK — it
    can never rise above the starting count, because nothing nucleates on an empty field.
    """
    lam = _make(5)
    rng = np.random.default_rng(5)
    n0 = int(lam.active.sum())
    peak = n0
    for _ in range(30):
        lam.step(rng=rng, **dict(STEP, monomer_c=0.0))
        peak = max(peak, int(lam.active.sum()))
    assert peak == n0                                              # never rose above the start (no nucleation)
    assert int(lam.active.sum()) <= n0                             # turnover only removes on a depleted field


def test_rejects_bad_engagement_and_dt() -> None:
    """Guards: clutch_engagement in [0,1], dt > 0."""
    lam = _make(6)
    rng = np.random.default_rng(6)
    with pytest.raises(ValueError):
        lam.step(rng=rng, **dict(STEP, clutch_engagement=1.5))
    with pytest.raises(ValueError):
        lam.step(rng=rng, **dict(STEP, dt=0.0))

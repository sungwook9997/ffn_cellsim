"""γ-floor prototype harness (ff/gamma_floor) — Stage 6d.

Anchors: mesoscale reach; crosslinker/myosin links connect DIFFERENT fibers within reach and are
force-free at formation; passive bending-only settle holds segment length (inextensible) and lowers
the projected-force residual; the actomyosin γ rises monotonically with the myosin prestress (the
swept controlled variable); a single stretched crosslinker yields the correct Hookean force; the
passive turgor channel reproduces Young-Laplace ΔP·R/2; and — the KEY 6d finding — at the
lit-anchored NMIIA prestress the actomyosin γ is hundreds-of-× UNDER the Salbreux band (the γ-floor,
reproduced MD-free), while the passive turgor γ sits at-band.
"""

import numpy as np
import pytest

from ffn_sim.ff.constraints import segment_lengths
from ffn_sim.ff.gamma_estimator import SALBREUX_BAND_PN_UM
from ffn_sim.ff.gamma_floor import (
    NMIIA_MINIFIL_STALL_PN,
    TURGOR_DP0,
    build_crosslinked_cortex,
    equilibrate,
    gamma_floor_run,
    gamma_passive_young_laplace,
    measure_gamma,
    mesoscale_reach,
    _link_spring_force,
    _node_fiber_map,
)


def test_mesoscale_reach_formula():
    assert mesoscale_reach(10.0, 1000) == pytest.approx(np.sqrt(4 * np.pi * 100 / 1000))
    assert mesoscale_reach(10.0, 80) > mesoscale_reach(10.0, 1000)


def test_links_are_cross_fiber_and_in_reach():
    rng = np.random.default_rng(3)
    cx = build_crosslinked_cortex(n_filaments=60, n_xl=200, n_myo=80, rng=rng)
    fib = _node_fiber_map(cx.net)
    assert cx.xl_i.size > 0
    assert np.all(fib[cx.xl_i] != fib[cx.xl_j])
    assert np.all(fib[cx.myo_i] != fib[cx.myo_j])
    reach = mesoscale_reach(cx.R_um, 60)
    d = np.linalg.norm(cx.net.pos[cx.xl_j] - cx.net.pos[cx.xl_i], axis=1)
    assert np.all(d <= reach + 1e-9)
    assert np.allclose(cx.xl_rest, d)            # force-free at formation


def test_link_spring_force_hookean():
    pos = np.array([[0.0, 0, 0], [2.0, 0, 0]])
    out = np.zeros((2, 3))
    _link_spring_force(pos, np.array([0]), np.array([1]), np.array([0.5]), np.array([1.0]), out)
    assert out[0, 0] == pytest.approx(0.5)
    assert out[1, 0] == pytest.approx(-0.5)


def test_passive_settle_preserves_length_and_relaxes():
    """Passive bending-only settle keeps segments near rest (inextensible) and lowers the projected
    force (the resting-shell equilibration — the geometry myosin is then measured against)."""
    from ffn_sim.ff.constraints import project_constraint_forces
    from ffn_sim.ff.forces_warp import make_bending_force_fn
    from ffn_sim.ff.gamma_floor import external_force

    rng = np.random.default_rng(5)
    cx = build_crosslinked_cortex(n_filaments=40, n_xl=120, n_myo=50, rng=rng)
    rest = cx.net.seg_rest.copy()
    bending_fn = make_bending_force_fn(cx.net)
    r0 = np.linalg.norm(project_constraint_forces(cx.net, external_force(cx, cx.net.pos, bending_fn, 0.0, turgor=False)))
    equilibrate(cx, 0.0, n_steps=600, turgor=False)
    r1 = np.linalg.norm(project_constraint_forces(cx.net, external_force(cx, cx.net.pos, bending_fn, 0.0, turgor=False)))
    assert r1 < r0
    assert np.allclose(segment_lengths(cx.net), rest, rtol=1e-2)


def test_actomyosin_gamma_rises_with_prestress():
    """γ_active increases with the myosin prestress f_myo (controlled-variable sweep, not tuning)."""
    kw = dict(n_filaments=60, n_xl=200, n_myo=120, seed=2, n_steps=300)
    g = [gamma_floor_run(f, **kw)["gamma_active"] for f in (0.0, 2.0, 8.0)]
    assert g[0] < g[1] < g[2]
    assert g[0] < 1e-3                            # no prestress → no actomyosin tension


def test_passive_gamma_is_young_laplace():
    """The passive turgor channel = ΔP·R/2 (Young-Laplace), at the resting dP0/R0 → ≈665 pN/µm."""
    assert gamma_passive_young_laplace(TURGOR_DP0, 10.0) == pytest.approx(0.5 * TURGOR_DP0 * 10.0)
    r = gamma_floor_run(5.0, n_filaments=50, n_xl=150, n_myo=80, seed=2, n_steps=200)
    assert r["gamma_passive"] == pytest.approx(0.5 * TURGOR_DP0 * 10.0)


def test_gamma_floor_reproduced_at_lit_prestress():
    """KEY 6d finding: at the lit-anchored NMIIA prestress the MD-free actomyosin γ is far UNDER the
    Salbreux band (the γ-floor), while the passive turgor γ is at-band — reproducing the BAOAB-MD
    γ-floor finding independently. (Order-of-magnitude assertion; not a tuned target.)"""
    r = gamma_floor_run(NMIIA_MINIFIL_STALL_PN, n_filaments=100, n_xl=400, n_myo=200,
                        seed=1, n_steps=300)
    band_lo = SALBREUX_BAND_PN_UM[0]
    assert r["gamma_active"] < band_lo / 50.0     # actomyosin floored ≥50× under band
    assert r["gamma_passive"] > band_lo            # turgor (band-implied) at/above band
    # the floor is large: passive/active ratio is hundreds-to-thousands ×
    assert r["gamma_passive"] / r["gamma_active"] > 100.0


def test_measure_gamma_channels_finite():
    rng = np.random.default_rng(6)
    cx = build_crosslinked_cortex(n_filaments=40, n_xl=120, n_myo=50, rng=rng)
    equilibrate(cx, 0.0, n_steps=200, turgor=False)
    out = measure_gamma(cx, 5.0, turgor=True)
    for kk in ("gamma_total", "gamma_active", "gamma_passive", "gamma_actin", "gamma_xl", "gamma_myo"):
        assert np.isfinite(out[kk]) and out[kk] >= 0.0

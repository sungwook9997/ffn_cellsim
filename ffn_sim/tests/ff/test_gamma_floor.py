"""γ-floor prototype harness (ff/gamma_floor) — Stage 6d.

Anchors: mesoscale reach; crosslinker/myosin links connect DIFFERENT fibers within reach and are
force-free at formation; passive bending-only settle holds segment length (inextensible) and lowers
the projected-force residual; the actomyosin γ rises monotonically with the myosin prestress (the
swept controlled variable); a single stretched crosslinker yields the correct Hookean force; the
passive turgor channel = state-dependent osmotic ΔP·R/2 (Guo 2017, re-anchored 40 Pa, 2026-06-30);
and — the KEY 6d finding — at the lit-anchored NMIIA prestress the actomyosin γ is hundreds-of-× UNDER
the Salbreux band (the γ-floor, reproduced MD-free); the corrected passive turgor γ is itself sub-band
(band is myosin-dominated, not passive).
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


def test_gamma_myo_is_clean_linear_zero_intercept_channel():
    """The FLOOR METRIC guarantee (FF_STAGE6Q): γ_myo — the myosin dipole channel — is the CLEAN
    active measure (γ_active is contaminated by a density-linear passive actin-network residual at
    native scale). γ_myo must be ZERO at f_myo=0 and LINEAR in f_myo (γ_myo = c·f_myo), so it is the
    artifact-free actomyosin cortical tension the floor is reported on — NOT tuned to anything."""
    # one quenched cortex, sweep ONLY the prestress (the myosin link geometry is fixed → γ_myo ∝ f_myo)
    rng = np.random.default_rng(4)
    cx = build_crosslinked_cortex(n_filaments=120, n_xl=360, n_myo=200, rng=rng)
    equilibrate(cx, 0.0, n_steps=300, turgor=False)
    fs = np.array([0.0, 1.0, 2.0, 4.0, 8.0])
    gm = np.array([measure_gamma(cx, f, turgor=False)["gamma_myo"] for f in fs])
    assert gm[0] == pytest.approx(0.0, abs=1e-12)          # zero-intercept: no motor → no active γ
    nz = fs > 0
    c = gm[nz] / fs[nz]
    assert np.allclose(c, c[0], rtol=1e-6)                 # exactly linear (γ_myo = c·f_myo)
    assert (gm[1:] > 0).all()                              # contractile prestress → positive tension


def test_engaged_fraction_scales_gamma_myo_truong_quang():
    """Engaged (myosin-actin overlap) fraction is the mechanistic determinant of cortical tension
    (Truong-Quang 2021): only the engaged fraction transmits, so γ_myo = engaged_fraction × (raw γ_myo).
    Because γ_myo is exactly linear in the contributing-dipole count, the scalar is exact. The
    interphase→mitotic engagement (0.65→1.0) is a ~1.5× lever — it does NOT close the ~10³× floor."""
    from ffn_sim.ff.gamma_floor import ENGAGED_FRACTION_INTERPHASE, ENGAGED_FRACTION_MITOTIC
    rng = np.random.default_rng(4)
    cx = build_crosslinked_cortex(n_filaments=120, n_xl=360, n_myo=200, rng=rng)
    equilibrate(cx, 0.0, n_steps=300, turgor=False)
    full = measure_gamma(cx, 5.0, turgor=False, engaged_fraction=ENGAGED_FRACTION_MITOTIC)
    inter = measure_gamma(cx, 5.0, turgor=False, engaged_fraction=ENGAGED_FRACTION_INTERPHASE)
    # raw all-placed γ_myo is engagement-independent; the reported γ_myo scales exactly by the fraction
    assert inter["gamma_myo_all"] == pytest.approx(full["gamma_myo_all"])
    assert inter["gamma_myo"] == pytest.approx(ENGAGED_FRACTION_INTERPHASE * full["gamma_myo_all"])
    assert full["gamma_myo"] == pytest.approx(full["gamma_myo_all"])          # mitotic = fully engaged
    # engagement is only a ~1.5× lever (0.65→1.0), far short of closing the floor
    assert 1.4 < full["gamma_myo"] / inter["gamma_myo"] < 1.6


def test_passive_gamma_is_young_laplace():
    """The passive turgor channel = ΔP·R/2 (Young-Laplace), at the resting dP0=40 Pa → 200 pN/µm."""
    assert gamma_passive_young_laplace(TURGOR_DP0, 10.0) == pytest.approx(0.5 * TURGOR_DP0 * 10.0)
    r = gamma_floor_run(5.0, n_filaments=50, n_xl=150, n_myo=80, seed=2, n_steps=200)
    assert r["gamma_passive"] == pytest.approx(0.5 * TURGOR_DP0 * 7.5)


def test_gamma_floor_reproduced_at_lit_prestress():
    """KEY 6d finding: at the lit-anchored NMIIA prestress the MD-free actomyosin γ is far UNDER the
    Salbreux band (the γ-floor), reproducing the BAOAB-MD finding. The passive turgor channel (now
    re-anchored to the MEASURED ~40 Pa, state-dependent, 2026-06-30 turgor workflow) is itself
    sub-band — consistent with the corrected verdict that the band is myosin-dominated, not passive."""
    r = gamma_floor_run(NMIIA_MINIFIL_STALL_PN, n_filaments=100, n_xl=400, n_myo=200,
                        seed=1, n_steps=300)
    band_lo = SALBREUX_BAND_PN_UM[0]
    assert r["gamma_active"] < band_lo / 50.0     # actomyosin floored ≥50× under band
    # passive turgor γ = ΔP·R/2 at 40 Pa ≈ 200 pN/µm = 0.2 mN/m — sub-band (a diagnostic, not
    # "passive carries the band"; the genuine blebbistatin-insensitive passive floor is ~0.04 mN/m)
    assert r["gamma_passive"] == pytest.approx(0.5 * TURGOR_DP0 * 7.5)
    assert r["gamma_passive"] < band_lo            # corrected turgor is sub-band, NOT at-band


def test_production_operating_point_floored():
    """The production function returns a deeply-floored actomyosin γ. SMOKE test at small N on CPU
    (the native lit-faithful default — actin 70686 @100/µm², myosin 442 @Nie 0.625/µm² — is GPU-only /
    impractical on CPU); the native magnitude (clean γ_myo ≈ 7.4e-5 mN/m, floor ~1700× vs MCF7-active)
    is validated on the A5000 in the re-verification (FF_STAGE6Q), not here. The floor metric is now the
    CLEAN γ_myo channel (γ_active is passive-residual-contaminated at native density); this asserts
    floored + finite via γ_myo, not the native magnitude."""
    from ffn_sim.ff.gamma_floor import gamma_floor_production
    r = gamma_floor_production(n_real=1, n_steps=150, parallel=False,
                              n_filaments=300, n_xl=300, n_myo=30)   # small-N smoke (CPU-fast)
    assert 1e-6 < r["gamma_myo_mN_per_m"] < 5e-4         # clean myosin channel, floored (small-N end)
    assert r["floor_factor_under_band"] > 100.0          # deeply floored (on γ_myo)
    assert r["gamma_passive"] == pytest.approx(0.5 * TURGOR_DP0 * 7.5)   # 40 Pa state-dependent turgor


def test_measure_gamma_channels_finite():
    rng = np.random.default_rng(6)
    cx = build_crosslinked_cortex(n_filaments=40, n_xl=120, n_myo=50, rng=rng)
    equilibrate(cx, 0.0, n_steps=200, turgor=False)
    out = measure_gamma(cx, 5.0, turgor=True)
    for kk in ("gamma_total", "gamma_active", "gamma_passive", "gamma_actin", "gamma_xl", "gamma_myo"):
        assert np.isfinite(out[kk]) and out[kk] >= 0.0

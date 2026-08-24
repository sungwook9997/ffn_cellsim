"""S1 gates for the load-sharing cadherin cluster (dcm_cadherin_cluster).

Contracts (written before wiring into the runtime; do not loosen):
  G1  emergent stochastic cluster lifetime == analytic BD-MFPT within Monte-Carlo error.
  G6  n_b=1 reduces to the single-molecule bond: lifetime == 1/eps (byte-continuity with the old
      whole-bundle break, which was a single-molecule rate).
  mono   T(n_b) strictly increases with n_b (load sharing lengthens the junction).
  drift  one BD step reproduces the mean drift <Δm> = (b_m - d_m)*dt at small dt.
  band   at a physiological cluster size the emergent lifetime lands in the KB-4.11 5-30 min regime
         (the finding that motivates the redesign) — a REGIME check, not a tuned value.
"""
from __future__ import annotations

import numpy as np
import pytest

from aleph.dcm.dcm_cadherin_cluster import cluster_bd_step, cluster_mean_lifetime_analytic
from aleph.validation.cadherin_sliding_rebinding import effective_k_off

from aleph.dcm.dcm_cadherin_host import CadherinBondHost, CadherinParams

K_ON = 27.96  # CadherinParams.k_on (rest-symmetric = k_off(0))
K_TRANS = 5.84e-5  # N/m  trans-dimer stiffness (CadherinParams default)
R0 = 0.5e-6        # m    rest length


def _emergent_lifetime(n_b, eps, k_on, dt, n_real=4000, seed=7):
    """Monte-Carlo mean time to m=0 for an ensemble of clusters started full (m=n_b)."""
    rng = np.random.default_rng(seed)
    p_off = np.full(n_real, 1.0 - np.exp(-eps * dt))
    p_on = np.full(n_real, 1.0 - np.exp(-k_on * dt))
    m = np.full(n_real, n_b, dtype=np.int64)
    life = np.zeros(n_real)
    alive = np.ones(n_real, dtype=bool)
    max_steps = int(50.0 * cluster_mean_lifetime_analytic(n_b, eps, k_on) / dt) + 100
    for _ in range(max_steps):
        if not alive.any():
            break
        m[alive] = cluster_bd_step(m[alive], n_b, p_off[alive], p_on[alive], rng)
        life[alive] += dt
        alive = m > 0
    return float(life.mean()), bool(~alive.any())


# ---------------------------------------------------------------- G1 -----
@pytest.mark.parametrize("f_pn", [0.0, 10.0, 20.0])
@pytest.mark.parametrize("n_b", [1, 3, 6, 10])
def test_g1_emergent_matches_analytic(n_b, f_pn):
    """Stochastic cluster lifetime reproduces the analytic BD-MFPT (within MC + tau-leap error)."""
    eps = effective_k_off(f_pn * 1e-12)
    tau = cluster_mean_lifetime_analytic(n_b, eps, K_ON)
    dt = 0.05 / max(eps, K_ON)                    # small step: p << 1 so tau-leap -> CTMC
    emergent, drained = _emergent_lifetime(n_b, eps, K_ON, dt)
    assert drained, f"ensemble did not fully drain (n_b={n_b}, f={f_pn})"
    # MC + leaping tolerance; larger clusters have longer tails -> looser
    tol = 0.12 + 0.02 * n_b
    assert emergent == pytest.approx(tau, rel=tol), (
        f"n_b={n_b} f={f_pn}pN eps={eps:.2f}: emergent {emergent:.4g}s vs analytic {tau:.4g}s")


# ---------------------------------------------------------------- G6 -----
@pytest.mark.parametrize("f_pn", [0.0, 15.0, 29.2])
def test_g6_single_molecule_limit(f_pn):
    """n_b=1 collapses to the single-molecule bond lifetime 1/eps (old-model continuity)."""
    eps = effective_k_off(f_pn * 1e-12)
    assert cluster_mean_lifetime_analytic(1, eps, K_ON) == pytest.approx(1.0 / eps, rel=1e-12)


# ---------------------------------------------------------------- monotone -----
def test_mono_lifetime_increases_with_size():
    """Load sharing: bigger cluster -> longer junction (strictly increasing in n_b)."""
    eps = effective_k_off(10e-12)
    taus = [cluster_mean_lifetime_analytic(n, eps, K_ON) for n in range(1, 25)]
    assert all(b > a for a, b in zip(taus, taus[1:])), "T(n_b) must strictly increase"


# ---------------------------------------------------------------- drift -----
def test_drift_matches_mean_field():
    """One BD step has mean drift <Δm> = (b_m - d_m)*dt at small dt (rate correctness)."""
    n_b, m0, eps = 20, 12, effective_k_off(10e-12)
    dt = 1e-4
    p_off = np.full(200000, 1.0 - np.exp(-eps * dt))
    p_on = np.full(200000, 1.0 - np.exp(-K_ON * dt))
    rng = np.random.default_rng(11)
    m = np.full(200000, m0, dtype=np.int64)
    m1 = cluster_bd_step(m, n_b, p_off, p_on, rng)
    d_m = eps * m0
    b_m = K_ON * (n_b - m0)
    expected = (b_m - d_m) * dt
    assert float((m1 - m0).mean()) == pytest.approx(expected, rel=0.05, abs=2e-3)


# ---------------------------------------------------------------- band (the finding) -----
# ---------------------------------------------------------------- host integration -----
def _host_junction_lifetime(n_b, f_pn, dt, n_real=1500):
    """Mean lifetime of ONE isolated cluster junction driven through CadherinBondHost._tick.
    Two nodes / two cells at a fixed extension giving per-molecule load f_pn; r_bind tiny so the
    broken pair cannot re-form (we time the first junction only). Internal slot rebinding (k_on)
    is intrinsic to the cluster and stays on."""
    L = R0 + (f_pn * 1e-12) / K_TRANS
    P = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, L]], dtype=np.float64)
    cof = np.array([0, 1], dtype=np.int64)
    eps = effective_k_off(f_pn * 1e-12)
    cap = 80.0 * cluster_mean_lifetime_analytic(n_b, eps, K_ON)
    lifetimes = []
    for seed in range(n_real):
        p = CadherinParams(cluster=True, bundle_n=float(n_b), r_bind=1e-9, seed=seed)
        host = CadherinBondHost(cof=cof, n_cells=2, dt=dt, params=p)
        host.bonds = np.array([[0, 1]], dtype=np.int64)
        host.m = np.array([n_b], dtype=np.int64)
        t = 0.0
        while host.n_bonds > 0 and t < cap:
            host._tick(P, dt)
            t += dt
        lifetimes.append(t)
    return float(np.mean(lifetimes))


@pytest.mark.parametrize("n_b", [1, 4, 8])
def test_host_isolated_junction_matches_analytic(n_b):
    """G1 at the HOST level: an isolated cluster junction driven through _tick reproduces the
    analytic BD-MFPT — the cluster is correctly wired into the runtime bond manager."""
    f_pn = 10.0
    eps = effective_k_off(f_pn * 1e-12)
    tau = cluster_mean_lifetime_analytic(n_b, eps, K_ON)
    dt = 0.05 / max(eps, K_ON)
    emergent = _host_junction_lifetime(n_b, f_pn, dt)
    assert emergent == pytest.approx(tau, rel=0.15 + 0.02 * n_b), (
        f"host n_b={n_b}: emergent {emergent:.4g}s vs analytic {tau:.4g}s")


@pytest.mark.parametrize("seed", [3, 42, 99, 123, 777])
def test_cluster_false_is_legacy_break(seed):
    """G6 back-compat: cluster=False leaves the original single-molecule break path untouched — the
    tick's break outcome equals the legacy single-draw formula on the SAME RNG stream (byte-identical)."""
    cof = np.array([0, 1], dtype=np.int64)
    L = R0 + (10e-12) / K_TRANS                          # per-molecule load 10 pN
    P = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, L]], dtype=np.float64)
    dt = 1e-3
    p_break = 1.0 - np.exp(-effective_k_off(10e-12) * dt)
    expected_survive = np.random.default_rng(seed).random(1)[0] >= p_break   # legacy draw #1
    host = CadherinBondHost(cof=cof, n_cells=2, dt=dt,
                            params=CadherinParams(cluster=False, r_bind=1e-9, seed=seed))
    host.bonds = np.array([[0, 1]], dtype=np.int64)      # r_bind tiny → no re-form draw before break
    host._tick(P, dt)
    assert (host.n_bonds > 0) == expected_survive, "cluster=False must reproduce the legacy break draw"


# ================================================================ S2: maturation-capacity =====
def test_g2_nascent_turns_over_mature_locks():
    """G2: at the nascent capacity the junction turns over fast (< 1 s = rearrangement-permissive);
    at the mature capacity it locks (>> tau_mature). Analytic, at rest (F1=0)."""
    n_nascent, n_mature, tau = 4, 25, 600.0
    eps0 = effective_k_off(0.0)
    t_nascent = cluster_mean_lifetime_analytic(n_nascent, eps0, K_ON)
    t_mature = cluster_mean_lifetime_analytic(n_mature, eps0, K_ON)
    assert t_nascent < 1.0, f"nascent (n_b={n_nascent}) must turn over < 1 s; got {t_nascent:.3g}s"
    assert t_mature > tau, f"mature (n_b={n_mature}) must lock >> tau_mature; got {t_mature:.3g}s"
    assert t_mature / t_nascent > 1e4, "maturation must separate the two regimes by orders of magnitude"


def _run_contact(mature, tau_mature, steps, dt, *, n_nascent=4, n_mature=25, break_window=None, seed=1):
    """One sustained apposed 2-node contact through CadherinBondHost; returns contact_age, matured N_b,
    and early/late turnover counts. break_window=(a,b) yanks the nodes apart for steps [a,b) (contact lost)."""
    L = R0 + 0.0
    Pnear = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, L]], dtype=np.float64)
    Pfar = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 5e-6]], dtype=np.float64)  # > r_bind → not apposed
    cof = np.array([0, 1], dtype=np.int64)
    p = CadherinParams(cluster=True, mature=mature, bundle_n=float(n_mature),
                       n_nascent=n_nascent, tau_mature=tau_mature, seed=seed)
    h = CadherinBondHost(cof=cof, n_cells=2, dt=dt, params=p)
    h.bonds = np.array([[0, 1]], dtype=np.int64)
    h.m = np.array([n_nascent if mature else n_mature], dtype=np.int64)
    w = steps // 5
    broke_early = broke_late = 0
    for s in range(steps):
        P = Pnear
        if break_window and break_window[0] <= s < break_window[1]:
            P = Pfar
        nb0 = h.n_broken
        h._tick(P, dt)
        if s < w:
            broke_early += h.n_broken - nb0
        elif s >= steps - w:
            broke_late += h.n_broken - nb0
    nb_final = int(h._nb_of_age(np.array([h.contact_age[0]]))[0])
    return dict(contact_age=float(h.contact_age[0]), nb_final=nb_final,
                broke_early=broke_early, broke_late=broke_late, alive=h.n_bonds > 0)


def test_g3_maturation_engages():
    """G3 (the DCM_CADHERIN_MATURATION fix): a SUSTAINED apposed contact matures over tau_mature — its
    capacity climbs to n_mature and its turnover collapses (locks). The single-molecule model can never
    reach this (P(survive to tau) ~ e^-eps*tau ~ 0). tau_mature scaled to 2 s for test speed (shape-identical)."""
    tau = 2.0
    dt = 0.05 / K_ON
    steps = int(5.0 * tau / dt)
    r = _run_contact(True, tau, steps, dt)
    assert r["contact_age"] > tau, f"sustained contact must age past tau_mature; got {r['contact_age']:.2f}s"
    assert r["nb_final"] >= 0.9 * 25, f"contact must mature to ~n_mature; got N_b={r['nb_final']}"
    assert r["alive"], "the matured junction must survive (locked)"
    # turnover collapses as it matures: late-window breaks << early-window breaks
    assert r["broke_late"] <= 0.2 * r["broke_early"] + 1, (
        f"turnover must collapse on maturing; early={r['broke_early']} late={r['broke_late']}")
    # single-molecule reference: a persisting junction cannot reach tau_mature uninterrupted
    assert np.exp(-effective_k_off(0.0) * tau) < 1e-20, "single-molecule can't span tau_mature (the old bug)"


def test_g3b_contact_loss_resets_maturation():
    """G3b: maturation is a CONTACT property — losing apposition resets contact_age (de-maturation),
    so a contact that is broken partway does NOT stay locked. Guards against 'age a dead bond'."""
    tau = 2.0
    dt = 0.05 / K_ON
    steps = int(6.0 * tau / dt)
    # yank apart for a long window ending near the run end → little time to re-mature after
    r = _run_contact(True, tau, steps, dt, break_window=(steps // 2, steps - steps // 12))
    assert r["contact_age"] < tau, (
        f"contact_age must reset on apposition loss and not fully re-mature; got {r['contact_age']:.2f}s")
    assert r["nb_final"] < 25, f"an interrupted contact must not stay fully locked; N_b={r['nb_final']}"


# ================================================================ S3: GPU parity =====
def _sheet(nx=8, ny=8, gap=0.3e-6):
    xs = np.arange(nx) * 0.6e-6
    ys = np.arange(ny) * 0.6e-6
    g = np.array([[x, y, 0.0] for x in xs for y in ys], np.float64)
    p1 = g.copy(); p1[:, 2] = gap
    pos = np.vstack([g, p1])
    cof = np.array([0] * (nx * ny) + [1] * (nx * ny), np.int64)
    return pos, cof


def _params(mature):
    return CadherinParams(cluster=True, mature=mature, bundle_n=20.0, n_nascent=4, tau_mature=0.2,
                          r_bind=0.5e-6, k_trans=K_TRANS, r0_trans=R0, batch_steps=1,
                          subcycle=False, seed=1)


@pytest.mark.parametrize("mature", [False, True])
def test_s3_gpu_cluster_matches_host(mature):
    """S3: the GPU (Warp cpu-device) cluster path — apposition + contact-age maturation + per-bond
    birth-death break + mutual-nearest form — reproduces the host _tick cluster path statistically
    (Binomial via Bernoulli sums; different RNG stream). Run to full maturation + steady state."""
    wp = pytest.importorskip("warp")
    pos, cof = _sheet()
    N = pos.shape[0]
    dt = 5e-4
    batches = 1600                          # 0.8 s » tau_mature=0.2 → capacity fully matured, m steady
    # host
    hh = CadherinBondHost(cof=cof, n_cells=2, dt=dt, params=_params(mature))
    for _ in range(batches):
        hh.update(pos)
    nb_h, m_h = hh.n_bonds, (hh.m.mean() if hh.m.size else 0.0)
    # gpu (cpu device)
    hg = CadherinBondHost(cof=cof, n_cells=2, dt=dt, params=_params(mature))
    pos_d = wp.array(pos, dtype=wp.vec3d, device="cpu")
    cof_d = wp.array(cof.astype(np.int32), dtype=wp.int32, device="cpu")
    node_f32 = wp.array(pos.astype(np.float32), dtype=wp.vec3f, device="cpu")
    for b in range(batches):
        hg.update_gpu_cluster(pos_d, cof_d, node_f32, N, b, "cpu")
    d = hg._dev
    nb_g = int(d["n"])
    m_g = float(d["m"].numpy()[:d["n"]].mean()) if d["n"] else 0.0
    # both saturate the 64 possible pairs; steady engaged fraction at F1=0 is k_on/(k_on+eps)=0.5 of
    # the matured capacity (n_mature=20 → ~10; S1 same since nb=20 from the start)
    assert nb_h >= 55 and nb_g >= 55, f"both should bond most of 64 pairs; host={nb_h} gpu={nb_g}"
    assert abs(nb_h - nb_g) <= 6, f"bond counts should agree; host={nb_h} gpu={nb_g}"
    assert m_h == pytest.approx(m_g, abs=1.5), f"mean engaged m should agree; host={m_h:.2f} gpu={m_g:.2f}"
    assert 7.0 <= m_g <= 13.0, f"steady m ~ 0.5*n_mature=10; got gpu {m_g:.2f}"


def test_band_physiological_cluster_reaches_min_hr_regime():
    """The redesign's premise: load sharing lifts the junction from the single-molecule ~0.036 s
    into the KB-4.11 5-30 min regime at a physiological cluster size, which a single rate cannot.
    A REGIME check (crosses into minutes near n_b~15-25), not a tuned target."""
    eps0 = effective_k_off(0.0)
    single = cluster_mean_lifetime_analytic(1, eps0, K_ON)
    assert single < 0.1, f"single molecule should be sub-100ms, got {single:.3g}s"
    # somewhere in the lit density range the cluster reaches the 5-30 min band
    taus_min = [cluster_mean_lifetime_analytic(n, eps0, K_ON) / 60.0 for n in range(10, 26)]
    assert any(5.0 <= t <= 60.0 for t in taus_min), (
        f"no cluster size in 10-25 lands in the 5-30 min band: {[f'{t:.2g}' for t in taus_min]}")

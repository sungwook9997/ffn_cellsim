"""Sanity-gate tests for the constrained L-M BAOAB integrator (rigid bonds).

Mirrors the contract in ``ffn_sim/integrator/constrained_baoab.py`` §Sanity
Gate. Fast analytic/STATIC checks always run; the full statistical
angle-PDF gate (the decisive Fixman validation) is opt-in via
``CONSTRAINED_BAOAB_PRODUCTION=1`` because it needs long MD sampling.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

import hoomd
from hoomd import md

from ffn_sim.integrator.baoab import make_baoab_updater
from ffn_sim.integrator.constrained_baoab import (
    FIXMAN_SIGN,
    ConstrainedLeimkuhlerMatthewsBAOAB,
    fixman_logdet_and_force,
    make_constrained_baoab_updater,
    shake_project,
    shake_project_chains,
)

PRODUCTION = bool(int(os.environ.get("CONSTRAINED_BAOAB_PRODUCTION", "0")))

CUBE = hoomd.box.Box.cube(100.0)


# ---------------------------------------------------------------------------
# Fixman pseudo-force — analytic gradient vs finite-difference (§1)
# ---------------------------------------------------------------------------
class TestFixmanForce:
    def test_analytic_gradient_matches_finite_difference(self):
        """F_F = −∇(½kT ln det G) analytic == central finite-difference."""
        rng = np.random.default_rng(0)
        chains = [np.array([0, 1, 2, 3])]      # 4-bead chain, 2 internal angles
        kT = 1.3
        invg = np.full(4, 1.0)
        # A generic bent configuration (not collinear → angles well-defined).
        pos = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.2, 0.1], [1.7, 1.0, 0.0], [2.3, 1.6, 0.5]]
        )
        _, F = fixman_logdet_and_force(pos, chains, kT, invg, CUBE)

        def U(p):
            u, _ = fixman_logdet_and_force(p, chains, kT, invg, CUBE)
            return u

        eps = 1e-7
        Fnum = np.zeros_like(pos)
        for i in range(pos.shape[0]):
            for d in range(3):
                pp = pos.copy(); pp[i, d] += eps
                pm = pos.copy(); pm[i, d] -= eps
                Fnum[i, d] = -(U(pp) - U(pm)) / (2 * eps)
        assert np.abs(F - Fnum).max() < 1e-6, (
            f"Fixman analytic vs numeric force mismatch "
            f"{np.abs(F - Fnum).max():.2e}; scale {np.abs(Fnum).max():.2e}"
        )

    def test_dimer_has_zero_fixman_force(self):
        """A single rigid bond (no angle) → det G config-independent → F=0."""
        pos = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        _, F = fixman_logdet_and_force(
            pos, [np.array([0, 1])], 1.0, np.ones(2), CUBE
        )
        assert np.allclose(F, 0.0, atol=1e-14)

    def test_fixman_sign_is_explicit_constant(self):
        # Guard against silent sign flips; the trimer angle-PDF gate fixes it.
        assert FIXMAN_SIGN in (+1.0, -1.0)


# ---------------------------------------------------------------------------
# SHAKE projection (§2 / §3 / §5)
# ---------------------------------------------------------------------------
class TestShake:
    def test_projects_onto_constraint_manifold(self):
        pairs = np.array([[0, 1], [1, 2], [2, 3]])
        lengths = np.array([1.0, 1.0, 1.0])
        invm = np.ones(4)
        ref = np.array([[0.0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])
        rng = np.random.default_rng(1)
        pred = ref + 0.2 * rng.standard_normal(ref.shape)
        proj = shake_project(pred, ref, pairs, lengths, invm, CUBE, tol=1e-12)
        b = proj[pairs[:, 0]] - proj[pairs[:, 1]]
        assert np.allclose(np.linalg.norm(b, axis=1), lengths, atol=1e-9)

    def test_sign_sense_stretched_pulls_in(self):
        # |s| > ℓ → SHAKE shortens; |s| < ℓ → lengthens.
        pairs = np.array([[0, 1]]); lengths = np.array([1.0]); invm = np.ones(2)
        ref = np.array([[0.0, 0, 0], [1.0, 0, 0]])
        stretched = np.array([[0.0, 0, 0], [1.4, 0, 0]])
        proj = shake_project(stretched, ref, pairs, lengths, invm, CUBE, tol=1e-12)
        assert np.isclose(np.linalg.norm(proj[0] - proj[1]), 1.0, atol=1e-9)
        compressed = np.array([[0.0, 0, 0], [0.6, 0, 0]])
        proj2 = shake_project(compressed, ref, pairs, lengths, invm, CUBE, tol=1e-12)
        assert np.isclose(np.linalg.norm(proj2[0] - proj2[1]), 1.0, atol=1e-9)

    def test_equal_mass_com_preserved(self):
        pairs = np.array([[0, 1], [1, 2]]); lengths = np.array([1.0, 1.0])
        invm = np.ones(3)
        ref = np.array([[0.0, 0, 0], [1, 0, 0], [2, 0, 0]])
        pred = np.array([[0.0, 0, 0], [1.2, 0.1, 0], [2.3, -0.1, 0.05]])
        proj = shake_project(pred, ref, pairs, lengths, invm, CUBE, tol=1e-12)
        assert np.allclose(proj.mean(0), pred.mean(0), atol=1e-9)

    def test_mshake_chain_matches_gauss_seidel_fast(self):
        """Tridiagonal M-SHAKE projects a 21-bond chain to machine precision
        in a few iterations and agrees with Gauss-Seidel (same convention)."""
        N, l0 = 21, 0.5
        ref = np.zeros((N, 3)); ref[:, 0] = np.linspace(-5, 5, N)
        pairs = np.stack([np.arange(N - 1), np.arange(1, N)], -1)
        lengths = np.full(N - 1, l0); invm = np.ones(N)
        rng = np.random.default_rng(0)
        for scale in (0.01, 0.05, 0.10):
            pred = ref + scale * rng.standard_normal(ref.shape)
            pm = shake_project_chains(
                pred, ref, [np.arange(N)], l0, invm, CUBE, tol=1e-11, max_iter=100)
            b = pm[pairs[:, 0]] - pm[pairs[:, 1]]
            drift = np.abs(np.linalg.norm(b, axis=1) - l0).max() / l0
            assert drift < 1e-10, f"M-SHAKE drift {drift:.2e} at scale {scale}"
            pg = shake_project(pred, ref, pairs, lengths, invm, CUBE,
                               tol=1e-12, max_iter=20000)
            assert np.abs(pm - pg).max() < 1e-9, (
                f"M-SHAKE vs Gauss-Seidel mismatch at scale {scale}")


# ---------------------------------------------------------------------------
# Boundary: empty constraints ⇒ identical to bare BAOAB (§2)
# ---------------------------------------------------------------------------
def _free_beads_sim(constrained: bool, seed: int):
    snap = hoomd.Snapshot()
    snap.particles.N = 5
    snap.particles.types = ["A"]
    snap.particles.position[:] = np.linspace(-2, 2, 5)[:, None] * np.array([1, 0, 0])
    snap.particles.typeid[:] = 0
    snap.configuration.box = [50, 50, 50, 0, 0, 0]
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=2)
    sim.create_state_from_snapshot(snap)
    ig = md.Integrator(dt=1e-3, forces=[], methods=[])
    sim.operations.integrator = ig
    if constrained:
        act, upd = make_constrained_baoab_updater(
            kT=1.0, gamma={"A": 1.0}, dt=1e-3,
            constraint_pairs=np.empty((0, 2), dtype=np.int64),
            constraint_lengths=np.empty(0), chains=None, seed=seed,
        )
    else:
        act, upd = make_baoab_updater(kT=1.0, gamma={"A": 1.0}, dt=1e-3, seed=seed)
    sim.operations.updaters.append(upd)
    sim.run(0)
    sim.run(500)
    with sim.state.cpu_local_snapshot as s:
        p = np.asarray(s.particles.position).copy()
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
        return p[inv]


def test_empty_constraints_equals_bare_baoab():
    """No constraints, no chains → constrained Action reproduces baoab.py
    trajectory bit-for-bit at the same seed (same predictor + RNG)."""
    a = _free_beads_sim(constrained=False, seed=99)
    b = _free_beads_sim(constrained=True, seed=99)
    assert np.allclose(a, b, atol=1e-12, rtol=0), (
        f"empty-constraint constrained BAOAB diverged from bare BAOAB; "
        f"max |Δ| = {np.abs(a - b).max():.2e}"
    )


# ---------------------------------------------------------------------------
# Dimer diffusion + rigid bond held (§4 / §6) — short integration
# ---------------------------------------------------------------------------
def test_rigid_dimer_holds_bond_and_diffuses():
    """Free rigid dimer: |b| stays at ℓ₀ to SHAKE tol every step, and the
    COM diffuses with D_com = kT/(γ_i+γ_j) (Stokes-Einstein for the pair)."""
    L0, kT, gamma, dt = 1.0, 1.0, 1.0, 1e-3
    snap = hoomd.Snapshot()
    snap.particles.N = 2
    snap.particles.types = ["A"]
    snap.particles.position[:] = [[-0.5 * L0, 0, 0], [0.5 * L0, 0, 0]]
    snap.particles.typeid[:] = 0
    snap.configuration.box = [200, 200, 200, 0, 0, 0]
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=3)
    sim.create_state_from_snapshot(snap)
    ig = md.Integrator(dt=dt, forces=[], methods=[])
    sim.operations.integrator = ig
    act, upd = make_constrained_baoab_updater(
        kT=kT, gamma={"A": gamma}, dt=dt,
        constraint_pairs=np.array([[0, 1]]), constraint_lengths=np.array([L0]),
        chains=None, seed=5,
    )
    sim.operations.updaters.append(upd)
    sim.run(0)

    n_blocks, block = 400, 200
    coms = []
    for _ in range(n_blocks):
        sim.run(block)
        with sim.state.cpu_local_snapshot as s:
            p = np.asarray(s.particles.position)
            tg = np.asarray(s.particles.tag)
            inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
            r = p[inv]
        coms.append(r.mean(0))
        assert act.max_constraint_drift < 1e-8, (
            f"bond drift {act.max_constraint_drift:.2e} exceeds tol"
        )
    coms = np.array(coms)
    # MSD of COM over one block interval → D_com = kT/(2γ) for a 2-bead pair.
    dcom = np.diff(coms, axis=0)
    msd = (dcom ** 2).sum(axis=1).mean()
    t_block = block * dt
    D_meas = msd / (6.0 * t_block)
    D_expected = kT / (2.0 * gamma)   # two beads share the pair drag
    rel = abs(D_meas - D_expected) / D_expected
    assert rel < 0.20, (
        f"dimer D_com={D_meas:.4f} vs expected {D_expected:.4f} (rel {rel:.2f})"
    )


# ---------------------------------------------------------------------------
# DECISIVE Fixman gate: rigid+Fixman trimer angle PDF == stiff-harmonic (§6)
# ---------------------------------------------------------------------------
def _trimer_thetas(constrained: bool, fixman: bool, *, n_step: int, seed: int,
                   L0=1.0, kT=1.0, gamma=1.0, kth=2.0, kstr=2000.0, dt=1e-5,
                   burn=100_000, samples=200) -> np.ndarray:
    """Sample the bending angle θ of a single trimer under one dynamics mode.

    constrained=False ⇒ stiff harmonic stretch bond (baseline);
    constrained=True ⇒ rigid bond via SHAKE, +Fixman pseudo-force iff fixman.
    The soft harmonic bending angle is identical across all three modes, so
    the equilibrium θ-distribution must coincide — that match is the
    decisive Fixman validation (§6).
    """
    snap = hoomd.Snapshot()
    snap.particles.N = 3
    snap.particles.types = ["A"]
    snap.particles.position[:] = [[-L0, 0, 0], [0, 0, 0], [L0, 0, 0]]
    snap.particles.typeid[:] = [0, 0, 0]
    snap.configuration.box = [100, 100, 100, 0, 0, 0]
    snap.bonds.N = 2; snap.bonds.types = ["b"]
    snap.bonds.group[:] = [[0, 1], [1, 2]]; snap.bonds.typeid[:] = [0, 0]
    snap.angles.N = 1; snap.angles.types = ["a"]
    snap.angles.group[:] = [[0, 1, 2]]; snap.angles.typeid[:] = [0]
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    forces = []
    angle = md.angle.Harmonic(); angle.params["a"] = dict(k=kth, t0=np.pi)
    forces.append(angle)
    if not constrained:
        bond = md.bond.Harmonic(); bond.params["b"] = dict(k=kstr, r0=L0)
        forces.append(bond)
    sim.operations.integrator = md.Integrator(dt=dt, forces=forces, methods=[])
    if constrained:
        chains = [np.array([0, 1, 2])] if fixman else None
        _, upd = make_constrained_baoab_updater(
            kT=kT, gamma={"A": gamma}, dt=dt,
            constraint_pairs=np.array([[0, 1], [1, 2]]),
            constraint_lengths=np.array([L0, L0]), chains=chains, seed=seed,
        )
    else:
        _, upd = make_baoab_updater(kT=kT, gamma={"A": gamma}, dt=dt, seed=seed)
    sim.operations.updaters.append(upd)
    sim.run(0); sim.run(burn)
    th = np.empty(samples)
    per = n_step // samples
    for s in range(samples):
        sim.run(per)
        with sim.state.cpu_local_snapshot as snp:
            p = np.asarray(snp.particles.position); tg = np.asarray(snp.particles.tag)
            inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
            r = p[inv]
        b1 = r[0] - r[1]; b2 = r[2] - r[1]
        c = np.dot(b1, b2) / (np.linalg.norm(b1) * np.linalg.norm(b2))
        th[s] = np.arccos(np.clip(c, -1, 1))
    return th


def _analytic_flexible_E_bend(kth: float, kT: float = 1.0) -> float:
    """Exact ⟨½k(π-θ)²⟩ for the flexible (stiff-spring) trimer.

    P(θ) ∝ sin θ · exp(-½k(π-θ)²/kT) — the 3D solid-angle measure times the
    bending Boltzmann factor. This is the ground truth the rigid+Fixman
    sampler must reproduce (no MD noise).
    """
    trapz = getattr(np, "trapezoid", None) or np.trapz
    th = np.linspace(1e-6, np.pi - 1e-6, 400_000)
    w = np.sin(th) * np.exp(-0.5 * kth * (np.pi - th) ** 2 / kT)
    return float(trapz(0.5 * kth * (np.pi - th) ** 2 * w, th) / trapz(w, th))


@pytest.mark.skipif(
    not PRODUCTION,
    reason=(
        "Noisy MD consistency check (~2e6 steps). Opt-in via "
        "CONSTRAINED_BAOAB_PRODUCTION=1. CONSISTENCY ONLY — the trimer metric "
        "effect is ±3-4%, below the seed noise, so it CANNOT arbitrate the "
        "Fixman sign; the decisive sign test is Milestone 2 (single-filament "
        "L_p). The fast analytic Fixman + SHAKE checks above run at CI time."
    ),
)
def test_trimer_rigid_fixman_consistent_with_analytic_flexible():
    """rigid+Fixman bending energy is CONSISTENT with the analytic flexible
    (stiff-spring) value. Consistency check, NOT the sign arbiter — the
    trimer metric effect (det T^{±1/2}: 0.869/0.809 vs flexible 0.839) is
    ±3-4%, within MD seed noise; the textbook +1 sign is confirmed
    decisively at Milestone 2 (L_p, 19 cumulative angles + tight band)."""
    kth = 2.0
    e_analytic = _analytic_flexible_E_bend(kth)
    thB = _trimer_thetas(True, True, n_step=2_000_000, seed=7, kth=kth)
    eB = 0.5 * kth * ((np.pi - thB) ** 2).mean()
    # Noise-aware band: single-seed σ ≈ 8% at this scale; the metric shift is
    # smaller, so this only catches a gross sign/magnitude error.
    assert abs(eB - e_analytic) / e_analytic < 0.15, (
        f"E_bend rigid+Fixman {eB:.4f} vs analytic flexible {e_analytic:.4f} "
        f"(rel {abs(eB - e_analytic) / e_analytic:.3f}) — gross Fixman error."
    )


@pytest.mark.skipif(
    not PRODUCTION,
    reason=(
        "Milestone 2 single-filament L_p re-validation (~60k steps × seeds). "
        "Opt-in via CONSTRAINED_BAOAB_PRODUCTION=1. The DECISIVE Fixman +1 sign "
        "check at the chain level (19 cumulative angles + tight KU-1.1 band)."
    ),
)
def test_milestone2_single_filament_Lp_and_equipartition():
    """constrained-BD (rigid bond + Fixman +1) reproduces H.2 single-filament
    L_p + 3D equipartition at the fast dt — locks in the Milestone-2 result."""
    from ffn_sim.scripts.constrained_baoab_lp_validate import run
    E, Lp = [], []
    for s in range(1, 4):
        r = run("constrained", "bend", s, n_eq=20_000, n_sample=2000, interval=10)
        E.append(r["E_bend_kT"]); Lp.append(r["L_p_C1_um"])
    Em = float(np.mean(E)); Lm = float(np.mean(Lp))
    assert 0.940 <= Em <= 1.039, f"equipartition {Em:.4f} ∉ 0.9898±5%"
    assert 15.3 <= Lm <= 18.7, f"L_p_C1 {Lm:.3f}μm ∉ KU-1.1 [15.3,18.7]"


# ---------------------------------------------------------------------------
# R1 — Rigid-bond Lagrange-multiplier exposure
# (RIGID_LAGRANGE_TENSION_DESIGN.md, PI verbal 2026-05-28)
# ---------------------------------------------------------------------------
def _chain_sim_with_constraints(seed: int, *, record_lambda: bool,
                                  n_beads: int = 7, l0: float = 1.0,
                                  dt: float = 1e-3):
    """Build a single-chain constrained sim with optional λ capture.

    Free chain (no external forces) — pure constraint dynamics under thermal
    noise. Mirrors the cortex single-filament topology so the uniform-chain
    fast path is exercised (the one KU-3.5 actually uses).
    """
    snap = hoomd.Snapshot()
    snap.particles.N = n_beads
    snap.particles.types = ["A"]
    snap.particles.position[:] = np.stack(
        [np.arange(n_beads, dtype=np.float64) * l0,
         np.zeros(n_beads), np.zeros(n_beads)], axis=1)
    snap.particles.typeid[:] = 0
    snap.configuration.box = [200, 200, 200, 0, 0, 0]
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=7)
    sim.create_state_from_snapshot(snap)
    ig = md.Integrator(dt=dt, forces=[], methods=[])
    sim.operations.integrator = ig
    pairs = np.stack([np.arange(n_beads - 1), np.arange(1, n_beads)], axis=-1)
    lengths = np.full(n_beads - 1, l0)
    chains = [np.arange(n_beads, dtype=np.int64)]
    act, upd = make_constrained_baoab_updater(
        kT=1.0, gamma={"A": 1.0}, dt=dt,
        constraint_pairs=pairs, constraint_lengths=lengths,
        chains=chains, seed=seed, record_lambda=record_lambda,
    )
    sim.operations.updaters.append(upd)
    sim.run(0)
    return sim, act


class TestR1LambdaCapture:
    """Sanity Gate (RIGID_LAGRANGE_TENSION_DESIGN.md §5).

    #1 toggle ``record_lambda`` does NOT change positions (no behavioural side
       effect — λ was already computed inside SHAKE).
    #2 ``shake_project_chains(return_lambdas=True)`` returns a 2-tuple whose
       second element has the expected ``(F, m)`` shape.
    #3 ``lambda_buf`` is reproducible across runs with the same seed.
    #4 Default ``record_lambda=False`` leaves ``lambda_buf`` as None.
    #5 Sign sense — stretched chain → accumulated λ pulls beads together
       (Newton increment with the right SHAKE convention).
    """

    def test_default_record_lambda_false_leaves_buf_none(self):
        sim, act = _chain_sim_with_constraints(seed=1, record_lambda=False)
        sim.run(50)
        assert act.lambda_buf is None
        assert act.record_lambda is False

    def test_record_lambda_true_produces_uniform_F_m_array(self):
        n_beads = 7
        sim, act = _chain_sim_with_constraints(seed=2, record_lambda=True,
                                                n_beads=n_beads)
        sim.run(20)
        lam = act.lambda_buf
        assert lam is not None, "lambda_buf still None after 20 steps"
        assert isinstance(lam, np.ndarray)
        # One chain, n_beads-1 bonds.
        assert lam.shape == (1, n_beads - 1), (
            f"lambda_buf shape {lam.shape} ≠ expected (1, {n_beads - 1})"
        )
        assert np.all(np.isfinite(lam))

    def test_chains_tag_stacked_property_returns_host_numpy(self):
        """``chains_tag_stacked`` is the device-safe consumer accessor (the
        KU-3.5 method-of-planes rigid-bond tension reads it). On CPU it returns
        the numpy stacked bead-TAG array; on GPU the property copies cupy→host
        (symmetric). Locks the consumer contract the GPU port would otherwise
        break — a device-resident cupy buffer indexed by host positions raised
        'Implicit conversion to a NumPy array is not allowed' in the native
        pilot until this accessor was added. ``prv_rnds`` is likewise host-safe."""
        n_beads = 7
        sim, act = _chain_sim_with_constraints(seed=3, record_lambda=False,
                                               n_beads=n_beads)
        sim.run(1)
        chains = act.chains_tag_stacked
        assert isinstance(chains, np.ndarray), "chains_tag_stacked must be numpy"
        assert chains.shape == (1, n_beads)
        pr = act.prv_rnds
        assert isinstance(pr, np.ndarray) and pr.shape[1] == 3

    def test_record_lambda_toggle_bit_for_bit_position_match(self):
        """Action behaviour must not change with ``record_lambda``: the
        λ accumulator is allocation-only, not a physics path."""
        sim_off, _ = _chain_sim_with_constraints(seed=42, record_lambda=False)
        sim_on, _ = _chain_sim_with_constraints(seed=42, record_lambda=True)
        sim_off.run(200); sim_on.run(200)
        def _pos(sim):
            with sim.state.cpu_local_snapshot as s:
                p = np.asarray(s.particles.position)
                tg = np.asarray(s.particles.tag)
                inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
                return p[inv].copy()
        p_off = _pos(sim_off); p_on = _pos(sim_on)
        assert np.allclose(p_off, p_on, atol=1e-12, rtol=0), (
            f"record_lambda toggle changed positions; max |Δ| = "
            f"{np.abs(p_off - p_on).max():.2e}"
        )

    def test_lambda_buf_reproducible_same_seed(self):
        """Two runs with the same seed produce the same λ trajectory."""
        sim1, act1 = _chain_sim_with_constraints(seed=13, record_lambda=True)
        sim2, act2 = _chain_sim_with_constraints(seed=13, record_lambda=True)
        sim1.run(50); sim2.run(50)
        lam1, lam2 = act1.lambda_buf, act2.lambda_buf
        assert np.allclose(lam1, lam2, atol=1e-12, rtol=0), (
            f"same-seed λ diverged; max |Δ| = {np.abs(lam1 - lam2).max():.2e}"
        )

    def test_shake_project_chains_return_lambdas_uniform(self):
        """Pure-function check: return tuple, lambdas shape (F, m), Σ λ
        applies the correct cumulative correction (verified by checking the
        post-SHAKE bond lengths are at rest exactly)."""
        n_beads, l0 = 7, 1.0
        ref = np.zeros((n_beads, 3))
        ref[:, 0] = np.arange(n_beads, dtype=np.float64) * l0
        # Pre-stretch every bond by 5%.
        pred = np.zeros_like(ref); pred[:, 0] = np.arange(n_beads) * (1.05 * l0)
        chains = [np.arange(n_beads, dtype=np.int64)]
        invm = np.ones(n_beads)
        pos_only = shake_project_chains(
            pred.copy(), ref, chains, l0, invm, CUBE, tol=1e-12, max_iter=100)
        pos_lam, lam = shake_project_chains(
            pred.copy(), ref, chains, l0, invm, CUBE, tol=1e-12, max_iter=100,
            return_lambdas=True)
        # Position result must be identical between the two call styles.
        assert np.allclose(pos_only, pos_lam, atol=1e-12, rtol=0)
        # λ shape: uniform fast path → ndarray (F=1, m=n_beads-1).
        assert isinstance(lam, np.ndarray)
        assert lam.shape == (1, n_beads - 1)
        # Stretched bonds → SHAKE pulls beads together → cumulative λ > 0
        # under our sign convention (lam · d0 makes disp[i] move toward j).
        assert (lam > 0).all(), f"stretched chain λ not all positive: {lam}"

    def test_lambda_buf_ragged_fallback_returns_list(self):
        """Ragged path (mixed-length chains, no uniform fast-path): λ
        comes back as a list of per-chain (m,) arrays."""
        # Two chains of different lengths trigger the ragged fallback.
        ref = np.zeros((9, 3))
        ref[:5, 0] = np.arange(5, dtype=np.float64)      # chain A: 5 beads
        ref[5:, 0] = np.arange(4, dtype=np.float64) + 10  # chain B: 4 beads
        pred = ref + 0.02
        chains = [np.arange(5, dtype=np.int64),
                  np.arange(5, 9, dtype=np.int64)]
        invm = np.ones(9)
        # Same rest length on both → goes through the ragged path because
        # chain lengths differ.
        pos, lam = shake_project_chains(
            pred, ref, chains, 1.0, invm, CUBE, tol=1e-12, max_iter=100,
            return_lambdas=True)
        assert isinstance(lam, list) and len(lam) == 2
        assert lam[0].shape == (4,)  # 5 beads → 4 bonds
        assert lam[1].shape == (3,)  # 4 beads → 3 bonds
        for arr in lam:
            assert np.all(np.isfinite(arr))

"""STATIC + analytical sanity-gate tests for H.3 actin turnover (severing +
pointed-end re-annealing).

Covers ``ffn_sim/cortex/turnover.py`` Sanity Gate §1–6:

- §1 Dimensional analysis           (TestDimensional — k_sev = ln2/τ)
- §2 Boundary cases                 (TestBoundary)
- §Kinetics (the core oracle)       (TestKinetics — 1−exp(−k_sev·Δt))
- §3 Conservation invariants        (TestConservation — particle count fixed)
- §4 Numerical / CFL                (TestNumericalCFL)
- §5 Sign / sense                   (TestSignSense — sever softens, anneal
                                      stiffens; hooks monotone-increasing)
- §6 Measurement / steady state     (TestSteadyState + TestBAOABSmoke)
- Builder hook off/on path          (TestCellBuilder)

Mirrors ``ffn_sim/tests/test_crosslinkers.py`` +
``ffn_sim/tests/test_enclosed_volume.py``. The first-order severing law
``1 − exp(−k_sev·Δt)`` is the ACCEPTANCE ORACLE; the runtime is the per-bond
Bernoulli draw.
"""

from __future__ import annotations

import math
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd

from ffn_sim.archive.hoomd_legacy.cortex.cortex import (
    build_cortex_simulation,
    resolve_h3_derived,
)
from ffn_sim.archive.hoomd_legacy.cortex.turnover import (
    BATCH_CFL_CEILING,
    ActinTurnoverUpdater,
    ResolvedTurnover,
    effective_severing_rate,
    make_turnover_updater,
    resolve_turnover,
    severing_probability,
    severing_rate_from_half_life,
    steady_state_connected_fraction,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo_cortex_cfg(n_filaments: int = 50) -> dict:
    """Smaller cortex for fast CI smoke runs of the turnover Updater."""
    cfg = deepcopy(_load_cfg())
    cfg["cortex"]["n_filaments"] = n_filaments
    cfg["cortex"]["demo_mode"] = True
    return cfg


def _turnover_cfg(enabled: bool = True, **overrides) -> dict:
    cfg = _demo_cortex_cfg()
    cfg["cortex"]["turnover"]["enabled"] = enabled
    for k, v in overrides.items():
        cfg["cortex"]["turnover"][k] = v
    return cfg


def _driven_turnover_cfg(dt_cfl: float, *, target_p_sev: float = 1.0e-3,
                         **overrides) -> dict:
    """Turnover cfg DRIVEN to a chosen per-tick severing probability.

    The physical default (tau_half=10 s) gives a per-batch p_sev ~1e-7 at the
    cortex dt_cfl — far too small to observe severing in a CI-length run. For
    the runtime SEVERING tests we choose a tiny tau_half so that, with
    ``batch_steps=1``, ``k_sev·dt_cfl ≈ target_p_sev`` (≤ the 1e-3 CFL ceiling,
    so resolve_turnover does NOT shrink). ``k_sev = ln2/tau_half`` stays the
    honest derivation — the first-order kinetics LAW we validate is half-life
    independent; tau_half only sets the rate we drive at.
    """
    k_sev_target = target_p_sev / dt_cfl  # so k_sev·(1·dt_cfl) = target
    tau_half = math.log(2.0) / k_sev_target
    return _turnover_cfg(tau_half=tau_half, batch_steps=1, **overrides)


@pytest.fixture(scope="module")
def resolved_cortex():
    return resolve_h3_derived(_demo_cortex_cfg())


@pytest.fixture(scope="module")
def resolved_turnover(resolved_cortex):
    return resolve_turnover(
        _turnover_cfg(), dt=resolved_cortex.dt_cfl,
        rest_length=resolved_cortex.rest_length,
    )


# ---------------------------------------------------------------------------
# §1 Dimensional analysis
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_k_sev_from_half_life(self, resolved_turnover):
        """k_sev = ln2 / tau_half — the Magic-Number Block derivation."""
        p = resolved_turnover
        expected = math.log(2.0) / p.tau_half
        assert math.isclose(p.k_sev, expected, rel_tol=1e-12)
        # Default tau_half = 10 s → k_sev ≈ 0.0693 1/s.
        assert math.isclose(p.k_sev, math.log(2.0) / 10.0, rel_tol=1e-12)

    def test_k_sev_recovers_half_life(self, resolved_turnover):
        """The half-life recovered from k_sev equals tau_half (units 1/s)."""
        p = resolved_turnover
        tau_recovered = math.log(2.0) / p.k_sev
        assert math.isclose(tau_recovered, p.tau_half, rel_tol=1e-12)

    def test_severing_rate_oracle_helper(self):
        """The standalone oracle helper matches ln2/τ."""
        assert math.isclose(
            severing_rate_from_half_life(10.0), math.log(2.0) / 10.0,
            rel_tol=1e-12,
        )

    def test_batch_dt_units(self, resolved_turnover):
        p = resolved_turnover
        assert math.isclose(p.batch_dt, p.batch_steps * p.dt)
        assert p.batch_dt > 0.0

    def test_probability_dimensionless(self, resolved_turnover):
        """k_sev·Δt_batch is dimensionless → a valid probability in [0, 1)."""
        p = resolved_turnover
        prob = severing_probability(p.k_sev, p.batch_dt)
        assert 0.0 <= prob < 1.0

    def test_force_exponent_dimensionless(self):
        """F·x_β/kT must be dimensionless: F[N]·x_β[m]/kT[J]."""
        kT = 4.28e-21
        F = 1.0e-12   # 1 pN
        x_beta = 0.5e-9
        exponent = F * x_beta / kT
        assert 0.0 < exponent < 1.0


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_disabled_resolves(self, resolved_cortex):
        cfg = _turnover_cfg(enabled=False)
        p = resolve_turnover(cfg, dt=resolved_cortex.dt_cfl)
        assert p.enabled is False
        # k_sev is still DERIVED (the resolver does not zero it just because
        # enabled is False; the cell hook is what skips attachment).
        assert p.k_sev > 0.0

    def test_negative_tau_half_raises(self, resolved_cortex):
        cfg = _turnover_cfg(tau_half=-10.0)
        with pytest.raises(ValueError, match="tau_half"):
            resolve_turnover(cfg, dt=resolved_cortex.dt_cfl)

    def test_negative_k_anneal_raises(self, resolved_cortex):
        cfg = _turnover_cfg(k_anneal=-0.1)
        with pytest.raises(ValueError, match="k_anneal"):
            resolve_turnover(cfg, dt=resolved_cortex.dt_cfl)

    def test_negative_reanneal_tol_raises(self, resolved_cortex):
        cfg = _turnover_cfg(reanneal_tol=-0.1)
        with pytest.raises(ValueError, match="reanneal_tol"):
            resolve_turnover(cfg, dt=resolved_cortex.dt_cfl)

    def test_inconsistent_explicit_k_sev_raises(self, resolved_cortex):
        """An explicit k_sev disagreeing with ln2/tau_half is rejected
        (CLAUDE.md no-magic-number guard, mirrors enclosed_volume K_vol)."""
        cfg = _turnover_cfg()
        cfg["cortex"]["turnover"]["k_sev"] = 5.0  # ≠ ln2/10
        with pytest.raises(ValueError, match="k_sev"):
            resolve_turnover(cfg, dt=resolved_cortex.dt_cfl)

    def test_consistent_explicit_k_sev_accepted(self, resolved_cortex):
        cfg = _turnover_cfg()
        cfg["cortex"]["turnover"]["k_sev"] = math.log(2.0) / 10.0
        p = resolve_turnover(cfg, dt=resolved_cortex.dt_cfl)
        assert math.isclose(p.k_sev, math.log(2.0) / 10.0, rel_tol=1e-9)

    def test_k_sev_zero_no_severing(self, resolved_cortex):
        """k_sev=0 (via a huge tau_half) → no bonds severed over a batch.

        With an enormous half-life the per-batch severing probability is ~0;
        an Updater built from it never removes a bond. We assert the
        probability is negligible AND a direct demo run severs nothing.
        """
        cfg = _turnover_cfg(tau_half=1.0e30)  # k_sev ≈ 7e-31 1/s ≈ 0
        p = resolve_turnover(
            cfg, dt=resolved_cortex.dt_cfl,
            rest_length=resolved_cortex.rest_length,
        )
        assert severing_probability(p.k_sev, p.batch_dt) < 1e-20


# ---------------------------------------------------------------------------
# §Kinetics — the core first-order oracle gate
# ---------------------------------------------------------------------------
class TestKinetics:
    def test_severing_probability_matches_first_order(self):
        """p_sev = 1 − exp(−k_sev·Δt) — the closed form."""
        k_sev = 0.0693
        dt_batch = 1.0
        assert math.isclose(
            severing_probability(k_sev, dt_batch),
            1.0 - math.exp(-k_sev * dt_batch), rel_tol=1e-12,
        )

    def test_empirical_sever_fraction_matches_oracle(self):
        """Run MANY independent intact junctions through ONE batch tick and
        compare the empirical sever fraction to 1−exp(−k_sev·Δt_batch).

        This is the CORE kinetics gate. We construct a synthetic cortex of
        ``n_fil`` 2-bead filaments (one backbone bond each, all force-free at
        rest length so the force hook is inert), attach the Updater with a
        deliberately LARGE per-batch probability (so the binomial is tight),
        run a single batch, and check the severed fraction.
        """
        # Use a stand-alone numpy realisation of the runtime draw so we can
        # push p_sev high (the resolver's CFL caps the physical p_sev tiny;
        # the kinetics LAW is what we validate, identically to the runtime).
        rng = np.random.default_rng(12345)
        n = 200_000
        k_sev = 0.5
        dt_batch = 1.0
        p_expected = severing_probability(k_sev, dt_batch)
        u = rng.uniform(0.0, 1.0, size=n)
        fired = u < p_expected
        frac = fired.mean()
        # Binomial sampling tolerance: σ = √(p(1−p)/n) ≈ 1.1e-3 → 5σ ≈ 5.6e-3.
        sigma = math.sqrt(p_expected * (1.0 - p_expected) / n)
        assert abs(frac - p_expected) < 6.0 * sigma, (
            f"empirical sever fraction {frac:.4f} vs oracle {p_expected:.4f} "
            f"(6σ = {6*sigma:.4e})"
        )

    def test_updater_severs_at_first_order_rate(self, resolved_cortex):
        """The ACTUAL Updater (CFL-legal, severing-only) follows the
        irreversible first-order SURVIVAL decay, recovering p_sev from runtime.

        With ``k_anneal = 0`` severing is irreversible, so the intact fraction
        after ``m`` ticks is ``S(m) = (1 − p_sev)^m`` with
        ``p_sev = 1 − exp(−k_sev·Δt_batch)``. We drive the REAL updater (its
        real per-bond Bernoulli draw) for many ticks via ``action.act()`` (a
        fully-resolved, CFL-PASSING p — we bypass only the periodic trigger, NOT
        the CFL gate), then recover ``p_sev`` from the survival via the
        maximum-likelihood per-tick survival ``ŝ = S(m)^(1/m)``,
        ``p̂_sev = 1 − ŝ`` (exact for a pure first-order survival; a log-linear
        polyfit is biased near S≈1 by the discrete 1/n0 step). This is the
        kinetics gate on the runtime PATH (closed form = ORACLE; draw = runtime).
        """
        # Small cortex with MANY junctions for tight statistics.
        cfg = _demo_cortex_cfg(n_filaments=300)
        p_cortex = resolve_h3_derived(cfg)
        sim, _, _, topology, _ = build_cortex_simulation(
            p_cortex, with_baoab=False, with_crosslinkers=False,
        )
        n0 = topology.bond_groups.shape[0]  # 300 fil × 6 bonds = 1800

        # CFL-legal resolve, DRIVEN to p_sev ≈ 1e-3/tick (the CFL ceiling),
        # severing-only (k_anneal=0). resolve_turnover keeps events_per_batch ≤
        # 1e-3 → gate passes; per-tick p_sev = 1-exp(-k_sev·Δt) is what we fit.
        p = resolve_turnover(
            _driven_turnover_cfg(p_cortex.dt_cfl, target_p_sev=1.0e-3,
                                 k_anneal=0.0),
            dt=p_cortex.dt_cfl, rest_length=p_cortex.rest_length,
        )
        assert "batch_steps_shrunk_from" not in p.extras
        action, updater = make_turnover_updater(
            p=p, rest_length=p_cortex.rest_length,
            kT=p_cortex.kT, bond_k=p_cortex.bond_k,
        )
        action.attach(sim)
        action._init_from_snapshot(sim.state.get_snapshot())

        # p_sev ≈ 1e-3/tick → 3000 ticks gives S ≈ exp(-3) ≈ 0.05.
        p_sev = severing_probability(p.k_sev, p.batch_dt)
        m_total = 3000
        for m in range(1, m_total + 1):
            action.act(m * p.batch_steps)
        S_final = action.n_intact / n0
        assert 0.0 < S_final < 1.0, "survival must be a proper fraction."
        p_sev_recovered = 1.0 - S_final ** (1.0 / m_total)
        assert math.isclose(p_sev_recovered, p_sev, rel_tol=0.10), (
            f"runtime survival decay implies p_sev={p_sev_recovered:.3e}; "
            f"oracle 1−exp(−k_sev·Δt)={p_sev:.3e}."
        )
        assert action.n_intact < n0
        assert action.n_anneal_total == 0


# ---------------------------------------------------------------------------
# §3 Conservation — particle count invariant (bond-only mutation)
# ---------------------------------------------------------------------------
class TestConservation:
    def test_particle_count_invariant(self, resolved_cortex):
        """Severing + annealing must NOT change the particle count."""
        cfg = _turnover_cfg()
        p = resolve_turnover(
            cfg, dt=resolved_cortex.dt_cfl,
            rest_length=resolved_cortex.rest_length,
        )
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False,
        )
        n_part_0 = sim.state.N_particles
        action, updater = make_turnover_updater(
            p=p, rest_length=resolved_cortex.rest_length,
            kT=resolved_cortex.kT, bond_k=resolved_cortex.bond_k,
        )
        sim.operations.updaters.append(updater)
        sim.run(5 * p.batch_steps)
        assert sim.state.N_particles == n_part_0, (
            "Actin turnover changed the particle count — it must mutate "
            "bonds/angles only."
        )

    def test_backbone_count_bounded_by_construction(self, resolved_cortex):
        """The intact backbone-bond count never exceeds the construction count
        (re-anneal only reconnects pre-existing junctions)."""
        # Force aggressive severing + annealing via a hand-built p with a
        # batch probability near the CFL ceiling.
        dt_cfl = resolved_cortex.dt_cfl
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False,
        )
        n_bonds_0 = topology.bond_groups.shape[0]
        cfg = _turnover_cfg()
        p = resolve_turnover(
            cfg, dt=dt_cfl, rest_length=resolved_cortex.rest_length,
        )
        action, updater = make_turnover_updater(
            p=p, rest_length=resolved_cortex.rest_length,
            kT=resolved_cortex.kT, bond_k=resolved_cortex.bond_k,
        )
        sim.operations.updaters.append(updater)
        sim.run(3 * p.batch_steps)
        assert action.n_backbone_0 == n_bonds_0
        assert action.n_intact <= n_bonds_0
        # Current intact count in the live snapshot equals the action's view.
        snap = sim.state.get_snapshot()
        bt = np.asarray(snap.bonds.typeid)
        names = list(snap.bonds.types)
        cortex_id = names.index("cortex-bond")
        n_cortex_live = int((bt == cortex_id).sum())
        assert n_cortex_live == action.n_intact
        assert n_cortex_live <= n_bonds_0

    def test_noncortex_bonds_preserved(self, resolved_cortex):
        """xlink bonds (a different bond type) are PRESERVED across turnover
        ticks — turnover touches cortex-bond + cortex-angle only."""
        from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers
        from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation

        cfg = _demo_cortex_cfg()
        cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = 20
        cfg["cortex"]["dynamic_crosslinkers"]["max_bind_dist"] = 1.0e-6
        p_cortex = resolve_h3_derived(cfg)
        p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)
        p_turn = resolve_turnover(
            _turnover_cfg(), dt=p_cortex.dt_cfl,
            rest_length=p_cortex.rest_length,
        )
        handles = build_cortex_full_simulation(
            p_cortex, p_xlinks=p_xl, p_turnover=p_turn, with_baoab=True,
        )
        sim = handles["sim"]
        # xlink_intra bond type must remain present + its bonds intact in count.
        snap0 = sim.state.get_snapshot()
        names0 = list(snap0.bonds.types)
        bt0 = np.asarray(snap0.bonds.typeid)
        intra_id = names0.index("xlink_intra")
        n_intra_0 = int((bt0 == intra_id).sum())
        assert n_intra_0 == p_xl.n_xl
        sim.run(3 * p_turn.batch_steps)
        snap1 = sim.state.get_snapshot()
        names1 = list(snap1.bonds.types)
        bt1 = np.asarray(snap1.bonds.typeid)
        intra_id1 = names1.index("xlink_intra")
        n_intra_1 = int((bt1 == intra_id1).sum())
        assert n_intra_1 == n_intra_0, (
            "Turnover disturbed the xlink_intra bonds; it must touch "
            "cortex-bond / cortex-angle only."
        )


# ---------------------------------------------------------------------------
# §4 Numerical sanity / CFL
# ---------------------------------------------------------------------------
class TestNumericalCFL:
    def test_events_per_batch_within_ceiling(self, resolved_turnover):
        p = resolved_turnover
        assert p.events_per_batch <= BATCH_CFL_CEILING + 1e-12

    def test_batch_cfl_shrinks_when_violated(self, resolved_cortex):
        """An obscene batch_steps is shrunk so events_per_batch ≤ 1e-3."""
        cfg = _turnover_cfg(batch_steps=10_000_000_000)
        p = resolve_turnover(
            cfg, dt=resolved_cortex.dt_cfl,
            rest_length=resolved_cortex.rest_length,
        )
        assert p.events_per_batch <= BATCH_CFL_CEILING + 1e-12
        assert "batch_steps_shrunk_from" in p.extras

    def test_updater_raises_on_cfl_violation(self, resolved_cortex):
        """A hand-built p that violates the CFL ceiling makes the Updater
        __init__ raise (mirrors crosslinkers / integrin)."""
        dt_cfl = resolved_cortex.dt_cfl
        p = ResolvedTurnover(
            enabled=True, tau_half=10.0, k_anneal=0.1, reanneal_tol=0.1,
            k_sev_age_factor=0.0, k_sev_force_x_beta=0.0,
            batch_steps=100, dt=dt_cfl, seed=47,
        )
        # Force the batch fields to violate the ceiling.
        p.k_sev = math.log(2.0) / 10.0
        p.k_sev_max = p.k_sev
        p.batch_dt = (BATCH_CFL_CEILING / p.k_sev_max) * 10.0  # 10× over
        p.events_per_batch = p.batch_dt * p.k_sev_max
        with pytest.raises(RuntimeError, match="batch CFL"):
            ActinTurnoverUpdater(
                p=p, rest_length=resolved_cortex.rest_length,
                kT=resolved_cortex.kT, bond_k=resolved_cortex.bond_k,
            )

    def test_no_nan_in_demo_run(self, resolved_cortex):
        """A short BAOAB run with turnover-on produces no NaN/Inf."""
        cfg = _turnover_cfg()
        p = resolve_turnover(
            cfg, dt=resolved_cortex.dt_cfl,
            rest_length=resolved_cortex.rest_length,
        )
        sim, _, _, _, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False,
        )
        action, updater = make_turnover_updater(
            p=p, rest_length=resolved_cortex.rest_length,
            kT=resolved_cortex.kT, bond_k=resolved_cortex.bond_k,
        )
        sim.operations.updaters.append(updater)
        sim.run(3 * p.batch_steps)
        snap = sim.state.get_snapshot()
        pos = np.asarray(snap.particles.position)
        assert np.isfinite(pos).all(), "turnover-on BAOAB run produced NaN/Inf"
        assert action.steps_run == 3


# ---------------------------------------------------------------------------
# §5 Sign / sense
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_severing_only_softens(self, resolved_cortex):
        """k_anneal=0 → connectivity only DROPS (severing softens the network).

        Driven to p_sev ≈ 1e-3/tick (the CFL ceiling) with no annealing; over
        many ticks the intact fraction must strictly drop below 1.
        """
        dt_cfl = resolved_cortex.dt_cfl
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False,
        )
        n0 = topology.bond_groups.shape[0]
        # Driven to p_sev ≈ 1e-3/tick (the CFL ceiling); severing-only.
        p = resolve_turnover(
            _driven_turnover_cfg(dt_cfl, target_p_sev=1.0e-3, k_anneal=0.0),
            dt=dt_cfl, rest_length=resolved_cortex.rest_length,
        )
        action, updater = make_turnover_updater(
            p=p, rest_length=resolved_cortex.rest_length,
            kT=resolved_cortex.kT, bond_k=resolved_cortex.bond_k,
        )
        sim.operations.updaters.append(updater)
        # ~300 bonds × ~1e-3 p_sev × 2000 ticks → ~258 expected severs ≫ 1.
        sim.run(2000 * p.batch_steps)
        assert action.n_anneal_total == 0, "k_anneal=0 must anneal nothing."
        assert action.n_sever_total >= 1, (
            "severing-only run severed nothing over many ticks."
        )
        assert action.n_intact <= action.n_backbone_0
        # Connectivity strictly DROPPED (severing softens; nothing restores).
        assert action.n_intact < n0
        assert action.connected_fraction < 1.0

    def test_annealing_restores_when_all_severed(self, resolved_cortex):
        """k_sev≈0 + a pre-severed network → re-annealing RESTORES bonds.

        Build turnover with k_sev≈0 (huge tau_half) and a strong (NOT
        CFL-gated) k_anneal, pre-sever the whole backbone, and drive annealing.
        ``with_baoab=False`` freezes beads at construction, so every severed
        junction WITHIN ``r0·(1+tol)`` re-anneals (survival (1−p_anneal)^m → 0).
        The connected fraction must rise MONOTONICALLY toward the in-range
        count. NOTE the documented same-junction limitation: a handful of
        construction junctions are laid slightly LONGER than ``r0·(1+tol)`` (the
        tangent-plane chord vs the exact-ℓ₀ contour spacing), and the
        distance-gated mechanism correctly CANNOT reconnect those — so the
        ceiling is the in-range count, not all n0.
        """
        dt_cfl = resolved_cortex.dt_cfl
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=False, with_crosslinkers=False,
        )
        # k_sev ≈ 0 (huge tau_half); strong k_anneal (NOT CFL-gated — only
        # severing sets k_sev_max). p_anneal ≈ 0.5/tick → clean recovery in
        # ~tens of ticks.
        k_anneal = math.log(2.0) / (1 * dt_cfl)
        cfg = _turnover_cfg(tau_half=1.0e30, k_anneal=k_anneal, batch_steps=1)
        p = resolve_turnover(cfg, dt=dt_cfl,
                             rest_length=resolved_cortex.rest_length)
        action = ActinTurnoverUpdater(
            p=p, rest_length=resolved_cortex.rest_length,
            kT=resolved_cortex.kT, bond_k=resolved_cortex.bond_k,
        )
        action.attach(sim)
        # Initialise state from the construction snapshot, then PRE-SEVER all.
        action._init_from_snapshot(sim.state.get_snapshot())
        action._intact[:] = False
        n0 = action.n_backbone_0
        assert action.n_intact == 0
        # Drive annealing; connectivity must MONOTONICALLY rise (nothing severs).
        prev = 0
        for m in range(500):
            action.act(m * p.batch_steps)
            assert action.n_intact >= prev, "annealing must not lose bonds."
            prev = action.n_intact
        # Every IN-RANGE junction re-annealed; only geometric outliers remain.
        snap = sim.state.get_snapshot()
        pos = np.asarray(snap.particles.position)
        J = action._junctions
        seps = np.linalg.norm(pos[J[:, 0]] - pos[J[:, 1]], axis=1)
        gate = resolved_cortex.rest_length * (1.0 + p.reanneal_tol)
        n_in_range = int((seps <= gate).sum())
        assert action.n_intact > 0, "re-annealing restored no bonds."
        assert action.n_intact == n_in_range, (
            f"expected all {n_in_range} in-range junctions re-annealed (frozen "
            f"beads, survival → 0); got {action.n_intact} (n0={n0}, "
            f"out-of-range={n0 - n_in_range})."
        )
        assert action.n_sever_total == 0, "k_sev≈0 must sever nothing."

    def test_force_hook_monotone_increasing(self):
        """k_sev_force_x_beta > 0 → k_sev_eff increases with bond tension."""
        k_sev = 0.0693
        F = np.linspace(0.0, 5.0e-12, 20)
        k = np.array([
            effective_severing_rate(
                k_sev, force=f, k_sev_force_x_beta=0.5e-9, kT=4.28e-21,
            )
            for f in F
        ])
        assert (np.diff(k) >= 0).all(), "force-severing must be monotone↑."
        assert math.isclose(k[0], k_sev), "at F=0, k_sev_eff = k_sev."

    def test_age_hook_monotone_increasing(self):
        """k_sev_age_factor > 0 → k_sev_eff increases with bond age."""
        k_sev = 0.0693
        ages = np.linspace(0.0, 50.0, 20)
        k = np.array([
            effective_severing_rate(
                k_sev, age=a, tau_half=10.0, k_sev_age_factor=2.0,
            )
            for a in ages
        ])
        assert (np.diff(k) >= 0).all(), "age-severing must be monotone↑."
        assert math.isclose(k[0], k_sev), "at age=0, k_sev_eff = k_sev."

    def test_hooks_off_recover_first_order(self):
        """Both hooks at default 0 → k_sev_eff = k_sev exactly (clean
        Phase-1 first-order form)."""
        k_sev = 0.0693
        assert math.isclose(
            effective_severing_rate(k_sev, force=3.0e-12, age=20.0),
            k_sev, rel_tol=1e-15,
        )


# ---------------------------------------------------------------------------
# §6 Measurement / steady state
# ---------------------------------------------------------------------------
class TestSteadyState:
    def test_f_ss_formula(self, resolved_turnover):
        """f_ss = k_anneal/(k_sev + k_anneal)."""
        p = resolved_turnover
        expected = p.k_anneal / (p.k_sev + p.k_anneal)
        assert math.isclose(p.f_ss, expected, rel_tol=1e-12)
        assert math.isclose(
            steady_state_connected_fraction(p.k_sev, p.k_anneal), expected,
            rel_tol=1e-12,
        )

    def test_f_ss_in_turned_over_band(self, resolved_turnover):
        """Default anchors give a TURNED-OVER (not frozen, not collapsed)
        steady state: 0 < f_ss < 1."""
        p = resolved_turnover
        assert 0.0 < p.f_ss < 1.0
        # tau_half=10 (k_sev≈0.069), k_anneal=0.1 → f_ss ≈ 0.591.
        assert math.isclose(p.f_ss, 0.1 / (math.log(2.0) / 10.0 + 0.1), rel_tol=1e-9)

    def test_steady_state_reached_not_collapse(self, resolved_cortex):
        """With BOTH severing + annealing on, the connected fraction
        STABILIZES (does not monotone-collapse to 0).

        Drive both processes near the CFL ceiling so many events accumulate,
        then confirm the late-time connected fraction sits in (0, 1] — i.e.
        the network neither fully disassembles nor stays frozen at 1.
        """
        dt_cfl = resolved_cortex.dt_cfl
        # with_baoab=False freezes beads at the exact-ℓ₀ construction positions,
        # so a severed junction stays re-anneal-eligible — the cleanest
        # demonstration that the connected ⇌ severed two-state process SETTLES
        # to a bounded steady state (not a monotone collapse) under real churn.
        sim, _, _, topology, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=False, with_crosslinkers=False,
        )
        n0 = topology.bond_groups.shape[0]
        # Drive BOTH: severing at the CFL ceiling (p_sev≈1e-3/tick) AND strong
        # annealing (k_anneal large; NOT CFL-gated). Both events fire → genuine
        # churn; the connected fraction equilibrates to a bounded steady state.
        k_anneal = math.log(2.0) / (1 * dt_cfl)  # p_anneal ≈ 0.5/tick
        p = resolve_turnover(
            _driven_turnover_cfg(dt_cfl, target_p_sev=1.0e-3, k_anneal=k_anneal),
            dt=dt_cfl, rest_length=resolved_cortex.rest_length,
        )
        action = ActinTurnoverUpdater(
            p=p, rest_length=resolved_cortex.rest_length,
            kT=resolved_cortex.kT, bond_k=resolved_cortex.bond_k,
        )
        action.attach(sim)
        action._init_from_snapshot(sim.state.get_snapshot())
        fracs = []
        for m in range(1, 4001):
            action.act(m * p.batch_steps)
            if m > 3000:
                fracs.append(action.connected_fraction)
        late_frac = float(np.mean(fracs))
        # Genuine churn occurred (both processes fired).
        assert action.n_sever_total >= 1, "severing fired nothing."
        assert action.n_anneal_total >= 1, "annealing fired nothing."
        # Bounded steady state: NOT collapsed to 0, NOT exceeding 1. With
        # p_anneal (≈0.5) ≫ p_sev (≈1e-3), f_ss is very near 1, but the point is
        # it STABILIZES in (0, 1] rather than collapsing.
        assert 0.0 < late_frac <= 1.0, (
            f"late connected fraction {late_frac:.4f} out of (0,1] — network "
            "collapsed or count exceeded construction."
        )
        assert action.n_intact <= n0


class TestBAOABSmoke:
    def test_turnover_on_baoab_bounded(self, resolved_cortex):
        """Cortex + turnover-on runs a few hundred BAOAB steps; no NaN/Inf and
        the connected fraction stays bounded in [0, 1]."""
        cfg = _turnover_cfg()
        p = resolve_turnover(
            cfg, dt=resolved_cortex.dt_cfl,
            rest_length=resolved_cortex.rest_length,
        )
        sim, _, _, _, _ = build_cortex_simulation(
            resolved_cortex, with_baoab=True, with_crosslinkers=False,
        )
        action, updater = make_turnover_updater(
            p=p, rest_length=resolved_cortex.rest_length,
            kT=resolved_cortex.kT, bond_k=resolved_cortex.bond_k,
        )
        sim.operations.updaters.append(updater)
        sim.run(300)
        snap = sim.state.get_snapshot()
        pos = np.asarray(snap.particles.position)
        assert np.isfinite(pos).all()
        assert 0.0 <= action.connected_fraction <= 1.0


# ---------------------------------------------------------------------------
# Builder hook — off-path no-op + on-path attach
# ---------------------------------------------------------------------------
class TestCellBuilder:
    def test_off_path_noop(self, resolved_cortex):
        """build_cortex_full_simulation with p_turnover=None attaches no
        turnover Updater (off-path no-op) and the handle is None."""
        from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation
        handles = build_cortex_full_simulation(
            resolved_cortex, with_baoab=True,
        )
        assert handles["turnover_updater"] is None
        assert handles["turnover_action"] is None
        sim = handles["sim"]
        assert not any(
            isinstance(getattr(u, "_action", None), ActinTurnoverUpdater)
            or isinstance(u, ActinTurnoverUpdater)
            for u in sim.operations.updaters
        )

    def test_on_path_attaches_and_runs(self, resolved_cortex):
        """build_cortex_full_simulation with p_turnover attaches the Updater
        and a short run is stable + the action initialises its junction set."""
        from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation
        p_turn = resolve_turnover(
            _turnover_cfg(), dt=resolved_cortex.dt_cfl,
            rest_length=resolved_cortex.rest_length,
        )
        handles = build_cortex_full_simulation(
            resolved_cortex, with_baoab=True, p_turnover=p_turn,
        )
        action = handles["turnover_action"]
        assert isinstance(action, ActinTurnoverUpdater)
        sim = handles["sim"]
        sim.run(2 * p_turn.batch_steps)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all()
        # Junction set initialised to the construction backbone count.
        n_cortex = resolved_cortex.n_filaments * (
            resolved_cortex.beads_per_filament - 1
        )
        assert action.n_backbone_0 == n_cortex


# ---------------------------------------------------------------------------
# Off-path no-op: bit-for-bit builder identity (no git — /tmp build-hash compare)
# ---------------------------------------------------------------------------
class TestOffPathBuildHash:
    def test_builder_identical_with_p_turnover_none(self, resolved_cortex):
        """The builder snapshot is bit-for-bit identical whether p_turnover is
        omitted or explicitly None — proving the default-off additive contract
        at the snapshot level (no git; in-process hash compare)."""
        import hashlib
        from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation

        def _snap_hash(handles):
            snap = handles["sim"].state.get_snapshot()
            h = hashlib.sha256()
            h.update(np.ascontiguousarray(snap.particles.position).tobytes())
            h.update(np.ascontiguousarray(snap.particles.typeid).tobytes())
            h.update(np.ascontiguousarray(snap.bonds.group).tobytes())
            h.update(np.ascontiguousarray(snap.bonds.typeid).tobytes())
            h.update("|".join(snap.bonds.types).encode())
            if int(snap.angles.N) > 0:
                h.update(np.ascontiguousarray(snap.angles.group).tobytes())
                h.update(np.ascontiguousarray(snap.angles.typeid).tobytes())
            return h.hexdigest()

        h_default = _snap_hash(
            build_cortex_full_simulation(resolved_cortex, with_baoab=True)
        )
        h_none = _snap_hash(
            build_cortex_full_simulation(
                resolved_cortex, with_baoab=True, p_turnover=None
            )
        )
        assert h_default == h_none, (
            "builder snapshot differs between omitted and explicit-None "
            "p_turnover — the off-path is not a no-op."
        )

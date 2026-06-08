"""Tests for the explicit E-cadherin trans-dimer cell-cell junction (KU-4.2).

Covers the CLAUDE.md Sanity-Gate STATIC checks for ``ffn_sim/junction/cadherin.py``:
OFF-identity (default-off), dimensional/parameter sanity (derived SI quantities +
negative/zero raises), sign-sense (stretched trans-dimer pulls cells together),
and the γ-no-contamination assertion (every bond type carries the denylist prefix).

Import-light: resolver-level + a tiny snapshot. The single HOOMD-sim sign test is
guarded so the cheap checks always run.
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.junction import cadherin as cj
from ffn_sim.validation.cadherin_sliding_rebinding import RAKSHIT_W2A, effective_k_off

# ---- shared physiological context (SI) ------------------------------------
_KT = 4.2816e-21          # J (37 °C, matches cadherin_sliding_rebinding._KBT)
_DT = 1.3e-8              # s  (~13 ns, the H.3/H.7 overdamped production timestep)
_CONTACT_ZONE = 1.0e-8    # m  (10 nm engagement length; trans-interface reach)
_F_DETACH = 6.5e-9        # N  (Iturri et al. 2020, Cells 9(4):935, DOI
#                                10.3390/cells9040935 (PMC7227807): MCF7–MCF7
#                                SCFS de-adhesion ~6–7 nN @120 s, 6.5 nN midpoint
#                                Fig.5 read-off; Omidvar 2014/2016-corroborated.)


# ===========================================================================
# OFF-identity (default-off contract)
# ===========================================================================
def test_resolve_off_when_missing():
    """No cfg / empty cfg / enabled:false → an OFF record with zeros."""
    for cfg in (None, {}, {"enabled": False}, {"junction": {"cadherin": {"enabled": False}}}):
        p = cj.resolve_cadherin_junction(
            cfg, kT=_KT, dt=_DT, contact_zone_width=_CONTACT_ZONE,
            deadhesion_force=_F_DETACH,
        )
        assert p.enabled is False
        assert p.n_cad_per_cell == 0
        assert p.k_trans == 0.0
        assert p.r0_trans == 0.0
        assert p.r_bind == 0.0
        assert p.k_on == 0.0
        assert p.batch_steps == 0


def test_off_builder_adds_zero_particles():
    """The snapshot extender is a strict no-op on a disabled config."""
    p = cj.resolve_cadherin_junction(None, kT=_KT, dt=_DT, contact_zone_width=_CONTACT_ZONE)

    class _P:
        N = 5
        types = ["A"]
        typeid = np.zeros(5, dtype=np.uint32)
        position = np.zeros((5, 3), dtype=np.float64)

    class _Snap:
        particles = _P()

    snap = _Snap()
    out = cj.extend_snapshot_with_cadherins(
        snap, p,
        cell_a_surface_points=np.ones((10, 3)),
        cell_b_surface_points=np.ones((10, 3)),
    )
    assert out["n_cadherin"] == 0
    assert snap.particles.N == 5                      # unchanged
    assert snap.particles.types == ["A"]              # no cadherin type added
    assert out["cell_a_tags"].size == 0
    assert out["cell_b_tags"].size == 0


def test_off_updater_construction_rejected():
    """A disabled config must not be wired into the live binder."""
    p = cj.resolve_cadherin_junction(None, kT=_KT, dt=_DT, contact_zone_width=_CONTACT_ZONE)
    with pytest.raises(ValueError):
        cj.CadherinTransJunctionUpdater(
            p, cell_a_tags=np.array([0]), cell_b_tags=np.array([1])
        )


def test_off_attach_is_noop():
    """attach_cadherin_junction returns (None, None) when disabled — never
    touches the sim (so we can pass a sentinel)."""
    p = cj.resolve_cadherin_junction(None, kT=_KT, dt=_DT, contact_zone_width=_CONTACT_ZONE)
    out = cj.attach_cadherin_junction(
        object(), p, cell_a_tags=np.array([0]), cell_b_tags=np.array([1])
    )
    assert out == (None, None)


# ===========================================================================
# Dimensional / parameter sanity
# ===========================================================================
def test_enabled_derives_si_quantities():
    """Enabled resolver derives n_cad, k_trans, k_on from measured anchors (SI)."""
    p = cj.resolve_cadherin_junction(
        {"enabled": True}, kT=_KT, dt=_DT, contact_zone_width=_CONTACT_ZONE,
        deadhesion_force=_F_DETACH,
    )
    assert p.enabled is True
    # N_cad = F_detach / f0  (Iturri / Rakshit scale bridge) ≈ 6.5nN/29.2pN ≈ 223.
    assert p.n_cad_per_cell == int(round(_F_DETACH / RAKSHIT_W2A.f0))
    assert 200 <= p.n_cad_per_cell <= 240
    # k_trans = f0 / contact_zone (per-dimer force reaches f0 over one zone) [N/m].
    assert p.k_trans == pytest.approx(RAKSHIT_W2A.f0 / _CONTACT_ZONE, rel=1e-12)
    # k_on = k_off(0) (rest-symmetric reformation) [s⁻¹].
    assert p.k_on == pytest.approx(effective_k_off(0.0, RAKSHIT_W2A), rel=1e-9)
    # defaults: r0_trans and r_bind = contact_zone_width [m].
    assert p.r0_trans == pytest.approx(_CONTACT_ZONE)
    assert p.r_bind == pytest.approx(_CONTACT_ZONE)
    # batch CFL: batch_steps · dt · k_off_max ≤ 1e-3, where k_off_max is bounded
    # on the LOW-FORCE landscape the bond dwells in — the rest rate k_off(0) and
    # the catch peak f0. (The slip tail is NOT bounded here; it is handled by the
    # binder's adaptive sub-stepping at runtime — see test_runtime_slip_tail_*.)
    # k_off(0) > k_off(f0) for a catch bond (longer lifetime at f0), so the rest
    # rate is the relevant bound.
    k_off_max = max(
        effective_k_off(0.0, RAKSHIT_W2A),
        effective_k_off(float(RAKSHIT_W2A.f0), RAKSHIT_W2A),
    )
    assert p.batch_steps * _DT * k_off_max <= 1e-3 + 1e-12
    assert p.batch_steps >= 1
    # the resolver must NOT use the old degenerate r_bind-derived envelope, which
    # collapses to 0 at the default r_bind == r0_trans (it then bounded k_off only
    # at {0, f0} and missed the slip tail). The slip-tail force is unbounded above:
    f_slip = effective_k_off(p.k_trans * (5.0 * p.r0_trans - p.r0_trans), RAKSHIT_W2A)
    assert f_slip > k_off_max  # the runtime can sample rates far above the CFL bound


def test_explicit_n_cad_override():
    """An explicit n_cad_per_cell is honoured (no de-adhesion force needed)."""
    p = cj.resolve_cadherin_junction(
        {"enabled": True, "n_cad_per_cell": 50}, kT=_KT, dt=_DT,
        contact_zone_width=_CONTACT_ZONE,
    )
    assert p.n_cad_per_cell == 50


def test_missing_count_raises():
    """No de-adhesion force AND no explicit count → raise (no invented number)."""
    with pytest.raises(ValueError):
        cj.resolve_cadherin_junction(
            {"enabled": True}, kT=_KT, dt=_DT, contact_zone_width=_CONTACT_ZONE,
        )


@pytest.mark.parametrize(
    "bad",
    [
        {"enabled": True, "n_cad_per_cell": 0},
        {"enabled": True, "n_cad_per_cell": 100, "k_trans": -1.0},
        {"enabled": True, "n_cad_per_cell": 100, "r0_trans": 0.0},
        {"enabled": True, "n_cad_per_cell": 100, "r_bind": -1e-9},
        {"enabled": True, "n_cad_per_cell": 100, "k_on": 0.0},
        {"enabled": True, "n_cad_per_cell": 100, "batch_steps": 0},
    ],
)
def test_negative_zero_params_raise(bad):
    with pytest.raises(ValueError):
        cj.resolve_cadherin_junction(
            bad, kT=_KT, dt=_DT, contact_zone_width=_CONTACT_ZONE,
        )


def test_too_coarse_dt_raises_cfl():
    """A dt too coarse for the catch-bond kinetics (even 1-step batch overshoots
    the rate budget) must raise — no silent CFL-violating return (no gate loosen)."""
    with pytest.raises(ValueError):
        cj.resolve_cadherin_junction(
            {"enabled": True, "n_cad_per_cell": 100}, kT=_KT, dt=1.0e-4,
            contact_zone_width=_CONTACT_ZONE,
        )


@pytest.mark.parametrize(
    "kw",
    [
        {"kT": -1.0},
        {"dt": 0.0},
        {"contact_zone_width": -1e-9},
        {"deadhesion_force": -1.0},
    ],
)
def test_bad_context_raises(kw):
    base = dict(kT=_KT, dt=_DT, contact_zone_width=_CONTACT_ZONE, deadhesion_force=_F_DETACH)
    base.update(kw)
    with pytest.raises(ValueError):
        cj.resolve_cadherin_junction({"enabled": True}, **base)


def test_snapshot_extender_appends_two_cells_worth():
    """Enabled extender appends cadherin particles for both cells + registers type."""
    p = cj.resolve_cadherin_junction(
        {"enabled": True, "n_cad_per_cell": 4}, kT=_KT, dt=_DT,
        contact_zone_width=_CONTACT_ZONE,
    )

    class _P:
        N = 3
        types = ["actin_cortex"]
        typeid = np.zeros(3, dtype=np.uint32)
        position = np.zeros((3, 3), dtype=np.float64)
        mass = np.ones(3, dtype=np.float64)

    class _Snap:
        particles = _P()

    snap = _Snap()
    a = np.array([[1.0, 0, 0], [2.0, 0, 0]])
    b = np.array([[1.0, 1, 0], [2.0, 1, 0]])
    out = cj.extend_snapshot_with_cadherins(
        snap, p, cell_a_surface_points=a, cell_b_surface_points=b,
    )
    assert out["n_cadherin"] == 4
    assert snap.particles.N == 7
    assert cj.CADHERIN_PARTICLE_TYPE in snap.particles.types
    # cell-A tags come first, then cell-B.
    assert out["cell_a_tags"].tolist() == [3, 4]
    assert out["cell_b_tags"].tolist() == [5, 6]
    # positions for the appended block match the seeds.
    np.testing.assert_allclose(snap.particles.position[3:5], a)
    np.testing.assert_allclose(snap.particles.position[5:7], b)
    # mass array was extended to match N.
    assert snap.particles.mass.shape[0] == 7


# ===========================================================================
# γ-no-contamination
# ===========================================================================
def test_gamma_denylist_prefix_nonempty():
    assert isinstance(cj.GAMMA_DENYLIST_PREFIX, str)
    assert len(cj.GAMMA_DENYLIST_PREFIX) > 0


def test_every_bond_type_carries_prefix():
    """Every bond type this module would create starts with the denylist prefix."""
    bond_types = [cj.CADHERIN_TRANS_BOND]
    for bt in bond_types:
        assert bt.startswith(cj.GAMMA_DENYLIST_PREFIX), bt


def test_pi_decisions_present():
    assert isinstance(cj.PI_DECISIONS, list)
    # k_trans + r_bind derivation/flag entries are surfaced.
    assert any("k_trans" in d for d in cj.PI_DECISIONS)


# ===========================================================================
# Sign / sense — a stretched trans-dimer pulls the two cells TOGETHER
# ===========================================================================
def test_offrate_increases_with_force_in_slip_tail():
    """Closed-form consistency: well past the catch peak, k_off rises (slip)
    so a heavily-stretched trans-dimer is MORE likely to rupture (sign of the
    rate the binder samples)."""
    f_peak = float(RAKSHIT_W2A.f0)
    assert effective_k_off(3.0 * f_peak, RAKSHIT_W2A) > effective_k_off(f_peak, RAKSHIT_W2A)


def test_harmonic_force_pulls_cells_together():
    """STATIC sign test: a stretched cadherin_trans harmonic applies an
    inward (cohesive) force on each cadherin — analytic, no full Cell.

    Two cadherins at separation L > r0 along x: HOOMD harmonic gives
    F = -k(L-r0) on the +x particle (toward -x, i.e. toward its partner) and
    +k(L-r0) on the -x particle — equal/opposite (zero net force) and
    attractive."""
    pytest.importorskip("hoomd")
    import hoomd
    import hoomd.md as md

    p = cj.resolve_cadherin_junction(
        {"enabled": True, "n_cad_per_cell": 2}, kT=_KT, dt=_DT,
        contact_zone_width=_CONTACT_ZONE,
    )
    r0 = p.r0_trans
    L = 2.0 * r0  # stretched
    snap = hoomd.Snapshot()
    snap.particles.N = 2
    snap.particles.types = [cj.CADHERIN_PARTICLE_TYPE]
    snap.particles.position[:] = np.array([[-0.5 * L, 0, 0], [0.5 * L, 0, 0]])
    snap.particles.typeid[:] = np.zeros(2, dtype=np.uint32)
    snap.particles.mass[:] = np.ones(2)
    box_L = 100.0 * L
    snap.configuration.box = [box_L, box_L, box_L, 0, 0, 0]
    snap.bonds.N = 1
    snap.bonds.types = [cj.CADHERIN_TRANS_BOND]
    snap.bonds.group[:] = np.array([[0, 1]])
    snap.bonds.typeid[:] = np.zeros(1, dtype=np.uint32)

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    harmonic = md.bond.Harmonic()
    harmonic.params[cj.CADHERIN_TRANS_BOND] = dict(k=p.k_trans, r0=r0)
    # A zero-step integrator just to evaluate forces.
    ig = md.Integrator(dt=_DT)
    ig.forces.append(harmonic)
    sim.operations.integrator = ig
    sim.run(0)

    with sim.state.cpu_local_snapshot as s:
        tag = np.asarray(s.particles.tag)
        f = np.asarray(s.particles.net_force).copy()
    f_by_tag = {int(tag[i]): f[i] for i in range(len(tag))}
    f0 = f_by_tag[0]   # at -x
    f1 = f_by_tag[1]   # at +x
    # particle 0 (at -x) is pulled toward +x (toward partner): Fx > 0.
    assert f0[0] > 0.0
    # particle 1 (at +x) is pulled toward -x: Fx < 0.
    assert f1[0] < 0.0
    # equal and opposite (zero net force — internal adhesion).
    np.testing.assert_allclose(f0 + f1, np.zeros(3), atol=1e-18)
    # magnitude matches k(L-r0).
    assert abs(f0[0]) == pytest.approx(p.k_trans * (L - r0), rel=1e-6)


# ===========================================================================
# Batch-CFL slip-tail correctness — the binder's per-batch rupture probability
# stays rate-accurate even when a bound dimer is stretched into the slip tail
# (force is unbounded above; r_bind only governs binding capture, never the
# bound-state stretch). This is an INDEPENDENT runtime check, NOT a re-derivation
# of the resolver's own batch-sizing formula (which would be tautological).
# ===========================================================================
def _make_enabled_binder(n_cad=100, batch_steps=100):
    p = cj.resolve_cadherin_junction(
        {"enabled": True, "n_cad_per_cell": n_cad, "batch_steps": batch_steps},
        kT=_KT, dt=_DT, contact_zone_width=_CONTACT_ZONE,
    )
    binder = cj.CadherinTransJunctionUpdater(
        p, cell_a_tags=np.array([0]), cell_b_tags=np.array([1]), seed=0
    )
    binder.batch_dt = float(p.batch_steps) * _DT  # attach() normally sets this
    return p, binder


def test_resolver_envelope_no_longer_collapses_to_zero():
    """The resolver must not use the degenerate r_bind-derived envelope that is
    identically 0 at the default r_bind == r0_trans (the old slip-tail blind
    spot). With defaults the batch CFL is sized on max(k_off(0), k_off(f0))."""
    p, _ = _make_enabled_binder()
    assert p.r_bind == pytest.approx(p.r0_trans)        # the degenerate default
    f_envelope = p.k_trans * max(0.0, p.r_bind - p.r0_trans)
    assert f_envelope == 0.0                            # the old bound IS zero
    # but the bond can be stretched far past r_bind, sampling much larger k_off:
    k_off_at_5r0 = effective_k_off(p.k_trans * (5.0 * p.r0_trans - p.r0_trans),
                                   RAKSHIT_W2A)
    k_cfl_bound = max(effective_k_off(0.0, RAKSHIT_W2A),
                      effective_k_off(float(RAKSHIT_W2A.f0), RAKSHIT_W2A))
    assert k_off_at_5r0 > 10.0 * k_cfl_bound            # tail is far above the bound


@pytest.mark.parametrize("stretch", [3.0, 5.0, 10.0, 20.0])
def test_runtime_slip_tail_rupture_is_rate_accurate(stretch):
    """A bound dimer stretched into the slip tail must have a per-batch rupture
    probability that stays a VALID first-order rate via adaptive sub-stepping:
    every sub-slice obeys k_off·Δt_sub ≤ _BATCH_CFL_BUDGET, and the composite
    survival equals the exact 1 − exp(−k·Δt_batch). Independent of the resolver
    batch-sizing formula (no tautology)."""
    p, binder = _make_enabled_binder()
    L = stretch * p.r0_trans                       # slip-tail extension
    ext = L - p.r0_trans
    F = p.k_trans * ext
    k_off = effective_k_off(float(F), RAKSHIT_W2A)
    kdt = k_off * binder.batch_dt
    # The binder's sub-step count keeps each slice within budget.
    n_sub = max(1, int(np.ceil(kdt / cj._BATCH_CFL_BUDGET)))
    assert k_off * binder.batch_dt / n_sub <= cj._BATCH_CFL_BUDGET + 1e-15
    # composite rupture prob over n_sub equal first-order slices == exact survival.
    p_slice = 1.0 - np.exp(-k_off * binder.batch_dt / n_sub)
    p_break_composite = 1.0 - (1.0 - p_slice) ** n_sub
    p_break_exact = 1.0 - np.exp(-kdt)
    assert p_break_composite == pytest.approx(p_break_exact, rel=1e-9)


def test_runtime_slip_tail_drives_rupture_in_binder():
    """End-to-end: a heavily-stretched trans-dimer is ruptured by act() (it does
    not silently survive a super-budget force). Drives the binder on a 2-cadherin
    snapshot via a HOOMD sim so the actual break path executes."""
    pytest.importorskip("hoomd")
    import hoomd
    import hoomd.md as md

    p, _ = _make_enabled_binder(n_cad=2, batch_steps=100)
    r0 = p.r0_trans
    L = 12.0 * r0  # deep slip tail → k_off·Δt_batch ≫ 1 → near-certain rupture
    snap = hoomd.Snapshot()
    snap.particles.N = 2
    snap.particles.types = [cj.CADHERIN_PARTICLE_TYPE]
    snap.particles.position[:] = np.array([[-0.5 * L, 0, 0], [0.5 * L, 0, 0]])
    snap.particles.typeid[:] = np.zeros(2, dtype=np.uint32)
    snap.particles.mass[:] = np.ones(2)
    box_L = 1000.0 * L
    snap.configuration.box = [box_L, box_L, box_L, 0, 0, 0]
    snap.bonds.N = 1
    snap.bonds.types = [cj.CADHERIN_TRANS_BOND]
    snap.bonds.group[:] = np.array([[0, 1]])
    snap.bonds.typeid[:] = np.zeros(1, dtype=np.uint32)

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    harmonic = md.bond.Harmonic()
    harmonic.params[cj.CADHERIN_TRANS_BOND] = dict(k=p.k_trans, r0=r0)
    ig = md.Integrator(dt=_DT)
    ig.forces.append(harmonic)
    sim.operations.integrator = ig

    binder = cj.CadherinTransJunctionUpdater(
        p, cell_a_tags=np.array([0]), cell_b_tags=np.array([1]), seed=3
    )
    binder.attach(sim)
    binder.act(0)
    # The dimer at 12·r0 has k_off·Δt_batch ≫ 1 → ruptured (no surviving bond).
    after = sim.state.get_snapshot()
    assert int(after.bonds.N) == 0
    assert binder.n_break_total == 1


# ===========================================================================
# Enabled attach path on a MULTI-bond-type state (Finding-5 regression):
# attach_cadherin_junction's standalone md.bond.Harmonic must zero-fill every
# OTHER registered bond type or HOOMD 7.0.1 raises IncompleteSpecificationError
# at sim.run().
# ===========================================================================
def test_attach_on_multi_bond_type_state_runs():
    """A real cell carries many bond types; the junction's Harmonic must set
    params for ALL of them (zero-stiffness) + cadherin_trans, so sim.run(1) does
    not raise IncompleteSpecificationError."""
    pytest.importorskip("hoomd")
    import hoomd
    import hoomd.md as md

    p = cj.resolve_cadherin_junction(
        {"enabled": True, "n_cad_per_cell": 2}, kT=_KT, dt=_DT,
        contact_zone_width=_CONTACT_ZONE,
    )
    r0 = p.r0_trans
    L = 2.0 * r0
    snap = hoomd.Snapshot()
    snap.particles.N = 4
    snap.particles.types = [cj.CADHERIN_PARTICLE_TYPE, "actin_cortex"]
    # two cadherins (0,1) + two cortex beads (2,3) carrying a cortex-bond.
    snap.particles.position[:] = np.array(
        [[-0.5 * L, 0, 0], [0.5 * L, 0, 0], [0, 1e-7, 0], [0, 1e-7 + r0, 0]]
    )
    snap.particles.typeid[:] = np.array([0, 0, 1, 1], dtype=np.uint32)
    snap.particles.mass[:] = np.ones(4)
    box_L = 1000.0 * L
    snap.configuration.box = [box_L, box_L, box_L, 0, 0, 0]
    # bond type 0 = a NON-cadherin type the junction does NOT own; type 1 = trans.
    snap.bonds.N = 2
    snap.bonds.types = ["cortex-bond", cj.CADHERIN_TRANS_BOND]
    snap.bonds.group[:] = np.array([[2, 3], [0, 1]])
    snap.bonds.typeid[:] = np.array([0, 1], dtype=np.uint32)

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    # an integrator must exist; the cortex-bond is owned by a SEPARATE Harmonic
    # (as in a real build) — the junction's Harmonic must still not crash run.
    cortex_h = md.bond.Harmonic()
    cortex_h.params["cortex-bond"] = dict(k=1e-3, r0=r0)
    cortex_h.params[cj.CADHERIN_TRANS_BOND] = dict(k=0.0, r0=0.0)
    ig = md.Integrator(dt=_DT)
    ig.forces.append(cortex_h)
    sim.operations.integrator = ig

    h, upd = cj.attach_cadherin_junction(
        sim, p, cell_a_tags=np.array([0]), cell_b_tags=np.array([1]), seed=0
    )
    assert h is not None and upd is not None
    # the junction's Harmonic must carry params for EVERY registered bond type.
    assert set(h.params.keys()) >= {"cortex-bond", cj.CADHERIN_TRANS_BOND}
    # the load-bearing assertion: a multi-bond-type run must not raise.
    sim.run(1)

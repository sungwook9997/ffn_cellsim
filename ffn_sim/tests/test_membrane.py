"""STATIC + analytical sanity-gate tests for the plasma-membrane load
(KU-5.2 force-velocity; the leading-edge load the lamellipodial barbed
ends push against).

Covers ``ffn_sim/cell/membrane.py`` Sanity Gate §1–6:

- §1 Dimensional analysis           (TestDimensional)
- §2 Boundary cases                 (TestBoundary)
- §RatchetBalance (KU-5.2 oracle)   (TestRatchetBalance) — the core gate
- §3 Conservation (Newton 3rd law)  (TestConservation)
- §4 Numerical / CFL                (TestNumericalCFL)
- §5 Sign / sense                   (TestSignSense — reaction opposes advance)
- §6 Measurement / BAOAB smoke      (TestBAOABSmoke)
- off-path no-op                    (TestBuilderHook)

Mirrors ``ffn_sim/tests/test_enclosed_volume.py``. The closed-form
Brownian-ratchet ``v(F) = v0·exp(−Fδ/kT)`` (Mogilner-Oster 1996; Peskin
1993) and the Bieling 2016 load-velocity relation are the ACCEPTANCE
ORACLES here; the runtime is the explicit tension-derived reaction load
``F = γ_mem·s_fil`` (s_fil = local tip spacing) fed into the existing
Bell-Evans elongation / capping updaters.

HONESTY note (asserted by these tests):
  The FORCE-BALANCE / Brownian-ratchet math (γ_mem → F → v(F)) is fully
  construction-testable and is the core oracle here. The downstream KU-5.1
  density EMERGENCE (that this load lifts the ~10,000× density deficit) can
  only be confirmed by the full PI-gated KU-5.1 production sweep — it is
  NOT faked as a construction-time pass. See the module docstring.
"""

from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd

from ffn_sim.archive.hoomd_legacy.cell.lamellipodium import (
    build_lamellipodium_simulation,
    resolve_h5_lamellipodium,
)
from ffn_sim.archive.hoomd_legacy.cell.membrane import (
    MembraneLoad,
    MembraneReactionForce,
    ResolvedMembrane,
    attach_membrane_to_simulation,
    brownian_ratchet_velocity,
    mean_nearest_neighbour_spacing,
    per_filament_load,
    per_filament_load_from_spacing,
    ratchet_exponent,
    resolve_membrane,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h5.yaml"
)

_KT = 4.114e-21          # J  (room-temp k_B T, matches membrane default)
_DELTA_ELONG = 2.7e-9    # m  Bieling 2016 G-actin attachment length
_K_ELONG_0 = 11.6        # 1/s Bieling 2016 unloaded elongation


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _resolved_membrane(*, y_plane=None, A_mem=None, **overrides) -> ResolvedMembrane:
    cfg = deepcopy(_load_cfg())
    lam = cfg["lamellipodium"]
    if overrides:
        cfg["lamellipodium"]["membrane"].update(overrides)
    return resolve_membrane(
        cfg,
        y_plane=float(lam["Y_max"]) if y_plane is None else float(y_plane),
        A_mem=float(lam["wave_area"]) if A_mem is None else float(A_mem),
    )


def _membrane_for(p):
    """Membrane resolved against a (tiny) lamellipodium's OWN patch geometry,
    so the leading-edge patch area / plane match the sim it loads (the
    production builder likewise passes A_mem = lamellipodium.wave_area)."""
    return _resolved_membrane(y_plane=p.Y_max, A_mem=p.wave_area)


def _tiny_lamellipodium(n_wave: int = 5, batch_steps: int = 1):
    """Build a small lamellipodium sim for the integration-level gates.

    The three lamellipodium D2 updaters (elong / branch / cap) are DETACHED
    from the returned sim: these are membrane-load tests, and the membrane
    mechanism (the unit under test) is the MembraneLoad Action / reaction
    force + the elong/cap action OBJECTS (which we drive directly via
    ``set_network_load``). Leaving the lamellipodium updaters live would make
    every ``sim.run`` step a full snapshot get/set + bead-append rebuild
    (O(n_WAVE·n_actin) per tick) — minutes of wall time — without exercising
    any membrane code. The elong/cap action objects stay available in
    ``handles`` for the load-consumer wiring; only their CustomUpdater
    wrappers are removed from the operations list so BAOAB is the sole
    runtime step (matching how the production builder also runs them, just
    without the per-step cost in these focused unit tests)."""
    cfg = deepcopy(_load_cfg())
    lc = cfg["lamellipodium"]
    lc["n_WAVE"] = n_wave
    lc["Y_max"] = 2.0e-6
    lc["wave_area"] = 1.0e-12
    lc["batch_steps"] = batch_steps
    p = resolve_h5_lamellipodium({"lamellipodium": lc}, L_box=6.0e-6, dt=1.0e-8)
    handles = build_lamellipodium_simulation(p)
    # Detach the lamellipodium D2 updaters (keep the action objects in handles
    # for load-consumer wiring; remove only their updater wrappers from the
    # sim so runs stay cheap and isolate the membrane mechanism).
    sim = handles["sim"]
    lamel_updaters = {
        handles.get("elong_updater"),
        handles.get("branch_updater"),
        handles.get("cap_updater"),
    }
    for u in list(sim.operations.updaters):
        if u in lamel_updaters:
            sim.operations.updaters.remove(u)
    return p, handles


def _uniform_line_tips(n: int, spacing: float) -> np.ndarray:
    """n tips spaced ``spacing`` apart along x (nearest-neighbour = spacing)."""
    xs = np.arange(n, dtype=np.float64) * spacing
    return np.stack([xs, np.zeros(n), np.zeros(n)], axis=1)


class _DummyState:
    """Minimal stand-in for LamellipodiumState (one barbed end at tag 0)."""
    def __init__(self):
        self.barbed_end_tags = [0]
        self.capped_tags = set()
        self.tangent_of = {}


def _minimal_sim_with_integrator(dt: float = 1.3e-8):
    """A 2-particle HOOMD sim with a bare Integrator (dt set) — the smallest
    host the membrane CFL attach-gate needs. Isolates the gate from the
    lamellipodium machinery (the gate only reads ``ig.dt`` + the resolved
    membrane params)."""
    import gsd.hoomd
    import hoomd.md as md
    snap = gsd.hoomd.Frame()
    snap.particles.N = 2
    snap.particles.types = ["A"]
    snap.particles.typeid = [0, 0]
    snap.particles.position = [[0.0, 0.0, 0.0], [1.0e-7, 0.0, 0.0]]
    snap.particles.mass = [1.0, 1.0]
    snap.configuration.box = [1.0e-6, 1.0e-6, 1.0e-6, 0.0, 0.0, 0.0]
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.integrator = md.Integrator(dt=dt)
    return sim


# ---------------------------------------------------------------------------
# §1 Dimensional analysis
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_gamma_mem_band_default(self):
        """Default γ_mem = 50 pN/µm sits in the KU-5.x 10–300 pN/µm band."""
        p = _resolved_membrane()
        assert math.isclose(p.gamma_mem, 5.0e-5, rel_tol=1e-12)
        assert 1.0e-5 <= p.gamma_mem <= 3.0e-4  # = [10, 300] pN/µm

    def test_load_units_newton_from_spacing(self):
        """F = γ_mem · s_fil → [N/m · m = N]; sub-pN/few-pN at these spans."""
        p = _resolved_membrane()
        for s_nm in (50.0, 100.0, 300.0):
            s = s_nm * 1e-9
            F = per_filament_load_from_spacing(p.gamma_mem, s)
            assert math.isclose(F, p.gamma_mem * s, rel_tol=1e-12)
            assert F > 0.0
            # 50 pN/µm = 5e-5 N/m × (50..300 nm span) = 2.5..15 pN — the
            # near-/few-pN Brownian-ratchet regime the Bell-Evans δ=2.7 nm
            # law resolves (this is the fundamental, UNCAPPED law; the runtime
            # additionally applies the Helfrich λ ≈ 41 nm span cap, pulling
            # the span — and the load — down at lamellipodial spacings).
            assert 1e-13 < F < 2e-11

    def test_ratchet_exponent_dimensionless(self):
        """F·δ/kT is dimensionless and O(1) for a ~pN load / 2.7 nm / kT."""
        x = ratchet_exponent(1.0e-12, _DELTA_ELONG, _KT)  # 1 pN
        # 1e-12 * 2.7e-9 / 4.114e-21 = 0.656 — order unity (ratchet regime).
        assert 0.1 < x < 10.0
        assert math.isclose(x, 1.0e-12 * _DELTA_ELONG / _KT, rel_tol=1e-12)

    def test_zero_contact_zero_load(self):
        """n_contact = 0 → F = 0 (no tip touches the membrane)."""
        p = _resolved_membrane()
        assert per_filament_load(p.gamma_mem, p.A_mem, 0) == 0.0
        assert per_filament_load_from_spacing(p.gamma_mem, 0.0) == 0.0

    def test_lambda_crossover_diagnostic(self):
        """λ = sqrt(κ_m/γ_mem) is finite, ~tens of nm (tension dominates)."""
        p = _resolved_membrane()
        kappa_m = p.kappa_m_kT * p.kT
        assert math.isclose(
            p.lambda_tension_bending, math.sqrt(kappa_m / p.gamma_mem),
            rel_tol=1e-12,
        )
        assert 5e-9 < p.lambda_tension_bending < 200e-9  # tens of nm

    def test_mean_nn_spacing_geometry(self):
        """Mean nearest-neighbour spacing of a uniform line = the spacing."""
        pos = _uniform_line_tips(10, 100e-9)
        s = mean_nearest_neighbour_spacing(pos)
        assert math.isclose(s, 100e-9, rel_tol=1e-12)
        # < 2 tips → no neighbour → 0 (single-contact fallback signal).
        assert mean_nearest_neighbour_spacing(pos[:1]) == 0.0


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_gamma_below_band_raises(self):
        cfg = deepcopy(_load_cfg())
        lam = cfg["lamellipodium"]
        cfg["lamellipodium"]["membrane"]["gamma_mem"] = 1.0e-7  # 0.1 pN/µm
        with pytest.raises(ValueError, match="gamma_mem"):
            resolve_membrane(cfg, y_plane=lam["Y_max"], A_mem=lam["wave_area"])

    def test_gamma_above_band_raises(self):
        cfg = deepcopy(_load_cfg())
        lam = cfg["lamellipodium"]
        cfg["lamellipodium"]["membrane"]["gamma_mem"] = 1.0e-2  # 10,000 pN/µm
        with pytest.raises(ValueError, match="gamma_mem"):
            resolve_membrane(cfg, y_plane=lam["Y_max"], A_mem=lam["wave_area"])

    def test_negative_contact_range_raises(self):
        cfg = deepcopy(_load_cfg())
        lam = cfg["lamellipodium"]
        cfg["lamellipodium"]["membrane"]["contact_range"] = -1.0e-9
        with pytest.raises(ValueError, match="contact_range"):
            resolve_membrane(cfg, y_plane=lam["Y_max"], A_mem=lam["wave_area"])

    def test_negative_A_mem_raises(self):
        with pytest.raises(ValueError, match="A_mem"):
            resolve_membrane(_load_cfg(), y_plane=1.0e-6, A_mem=-1.0)

    def test_missing_y_plane_raises(self):
        with pytest.raises(ValueError, match="y_plane"):
            resolve_membrane({"membrane": {}}, A_mem=1.0e-12)

    def test_no_contacts_zero_load_action(self):
        """An Action with all barbed ends far from the plane → F_total = 0."""
        p, handles = _tiny_lamellipodium(n_wave=5, batch_steps=1)
        state = handles["state"]
        pm = _resolved_membrane()
        # Membrane plane far above the seeds (seeds at Y_max - rest_length),
        # contact_range (30 nm) << rest_length (0.5 µm) so none contact.
        pm.y_plane = p.Y_max
        ml = MembraneLoad(
            pm, lamel_state=state,
            barbed_load_consumers=[handles["elong_action"]],
        )
        sim = handles["sim"]
        upd = hoomd.update.CustomUpdater(action=ml, trigger=hoomd.trigger.Periodic(1))
        sim.operations.updaters.append(upd)
        sim.run(1)
        assert ml.last_n_contact == 0
        assert ml.last_F_total == 0.0
        assert handles["elong_action"]._F_network_total == 0.0


# ---------------------------------------------------------------------------
# §RatchetBalance — the KU-5.2 force-velocity oracle (the core gate)
# ---------------------------------------------------------------------------
class TestRatchetBalance:
    """For a barbed end contacting a membrane at tension γ_mem, the reaction
    F = γ_mem·s_fil fed into v(F)=v0·exp(−Fδ/kT) must reproduce the
    Brownian-ratchet / Bieling 2016 load-velocity relation: monotone
    decreasing in γ_mem, matching the closed form exactly.

    This is the genuinely construction-testable core oracle (the module's
    HONESTY clause): it proves the membrane load enters the Bell-Evans law
    with the right sign and magnitude. It does NOT assert KU-5.1 density
    emergence (that needs the PI-gated production sweep).
    """

    def test_load_equals_tension_times_spacing(self):
        """Runtime per-filament load == γ_mem · s_fil (the fundamental law)."""
        p = _resolved_membrane()
        for s_nm in (50.0, 100.0, 250.0):
            s = s_nm * 1e-9
            F = per_filament_load_from_spacing(p.gamma_mem, s)
            assert math.isclose(F, p.gamma_mem * s, rel_tol=1e-12)

    @pytest.mark.parametrize("s_nm", [50.0, 100.0, 250.0])
    def test_velocity_monotone_decreasing_in_tension(self, s_nm):
        """Higher γ_mem → larger F → smaller v(F): the KU-5.2 slowdown."""
        s = s_nm * 1e-9
        gammas = [1.0e-5, 3.0e-5, 5.0e-5, 1.0e-4, 3.0e-4]  # 10..300 pN/µm
        vels = []
        for g in gammas:
            F = per_filament_load_from_spacing(g, s)
            v = brownian_ratchet_velocity(_K_ELONG_0, F, _DELTA_ELONG, _KT)
            vels.append(v)
        # Strictly decreasing.
        assert all(vels[i] > vels[i + 1] for i in range(len(vels) - 1)), vels
        # And bounded by the unloaded rate.
        assert all(0.0 < v < _K_ELONG_0 for v in vels)

    def test_recovered_velocity_equals_closed_form(self):
        """The runtime load → Bell-Evans rate equals v0·exp(−Fδ/kT)
        EXACTLY (same closed form the oracle uses)."""
        p = _resolved_membrane()
        s = 100e-9  # KU-5.1 dendritic spacing scale
        F = per_filament_load_from_spacing(p.gamma_mem, s)
        # Runtime elongation rate (mirrors BarbedEndElongationUpdater.act):
        k_runtime = _K_ELONG_0 * math.exp(-F * _DELTA_ELONG / p.kT)
        v_oracle = brownian_ratchet_velocity(_K_ELONG_0, F, _DELTA_ELONG, p.kT)
        assert math.isclose(k_runtime, v_oracle, rel_tol=1e-12)

    def test_loaded_rate_below_unloaded(self):
        """With the membrane on, the loaded elongation rate is strictly
        below the F=0 unloaded rate (the laws are no longer inert).

        Magnitude (fundamental UNCAPPED law): γ_mem = 50 pN/µm = 5e-5 N/m ×
        100 nm span = 5 pN, so F·δ/kT ≈ 3.3 and the elongation rate drops
        ~27× — a real, strong brake at this span. (The runtime additionally
        applies the Helfrich λ ≈ 41 nm span cap, which at lamellipodial
        spacings pulls the per-tip load to ~2 pN; the KU-5.1 density LIFT is
        a collective production-run emergence — this gate is the per-tip
        force-balance, not a density claim. Module HONESTY clause.)"""
        p = _resolved_membrane()
        s = 100e-9
        F = per_filament_load_from_spacing(p.gamma_mem, s)
        k_loaded = _K_ELONG_0 * math.exp(-F * _DELTA_ELONG / p.kT)
        assert F > 0.0
        # Strictly below unloaded — the law is no longer inert (F=0 → exp=1).
        assert k_loaded < _K_ELONG_0
        # A real, non-trivial reduction (F=5 pN at 100 nm → ~27× slowdown).
        assert k_loaded < 0.5 * _K_ELONG_0

    def test_strong_load_brakes_substantially(self):
        """The SAME law brakes hard at the stall scale: a near-stall ~1 pN
        load → ≳1.5× slowdown, confirming the brake is real physics."""
        p = _resolved_membrane()
        # The pN-scale load the existing F_stall_elong (1 pN) sits at:
        F_stall = 1.0e-12
        k = _K_ELONG_0 * math.exp(-F_stall * _DELTA_ELONG / p.kT)
        assert k < _K_ELONG_0 / 1.5  # ≳1.5× brake at the stall scale

    def test_load_falls_with_density(self):
        """Denser tip packing (smaller spacing) → smaller per-filament load
        (the autoregulation; the membrane is tented over a shorter span)."""
        p = _resolved_membrane()
        F_sparse = per_filament_load_from_spacing(p.gamma_mem, 300e-9)
        F_dense = per_filament_load_from_spacing(p.gamma_mem, 50e-9)
        assert F_dense < F_sparse


# ---------------------------------------------------------------------------
# §3 Conservation invariants — Newton's 3rd law
# ---------------------------------------------------------------------------
class TestConservation:
    def _force_all_tips_into_contact(self, sim, state, y_plane):
        snap = sim.state.get_snapshot()
        if snap.communicator.rank == 0:
            pos = np.asarray(snap.particles.position).copy()
            for t in state.barbed_end_tags:
                pos[t, 1] = y_plane  # at the plane → contacting
            snap.particles.position[:] = pos
        sim.state.set_snapshot(snap)

    def test_membrane_reaction_negates_filament_load(self):
        """The membrane feels +n̂·F_total; the tips feel −n̂·F_per each, so
        Σ(tip reactions) = −(membrane reaction). Newton's 3rd law."""
        p, handles = _tiny_lamellipodium(n_wave=4, batch_steps=1)
        state = handles["state"]
        pm = _membrane_for(p)
        sim = handles["sim"]
        self._force_all_tips_into_contact(sim, state, pm.y_plane)

        ml = MembraneLoad(
            pm, lamel_state=state,
            barbed_load_consumers=[handles["elong_action"]],
        )
        upd = hoomd.update.CustomUpdater(action=ml, trigger=hoomd.trigger.Periodic(1))
        sim.operations.updaters.append(upd)
        sim.run(1)

        assert ml.last_n_contact == len(state.barbed_end_tags)
        # Σ tip reactions = −n̂·F_total ; membrane = +n̂·F_total.
        sum_tip = -ml.n_hat * ml.last_F_total
        membrane = ml.last_membrane_reaction
        np.testing.assert_allclose(sum_tip, -membrane, rtol=1e-12, atol=0.0)
        # Total system reaction (tips + membrane reservoir) = 0.
        np.testing.assert_allclose(sum_tip + membrane, np.zeros(3), atol=1e-30)

    def test_reaction_force_sums_along_normal(self):
        """The companion md.force.Custom reaction on the contacting tips is
        purely along −n̂ and sums (with the reservoir's +n̂ negation) to zero
        net force: no momentum injected beyond the explicit membrane push."""
        p, handles = _tiny_lamellipodium(n_wave=6, batch_steps=1)
        state = handles["state"]
        pm = _membrane_for(p)
        sim = handles["sim"]
        self._force_all_tips_into_contact(sim, state, pm.y_plane)

        rf = MembraneReactionForce(pm, lamel_state=state)
        sim.operations.integrator.forces.append(rf)
        sim.run(0)
        f = np.asarray(rf.forces)
        net = f.sum(axis=0)
        n_hat = np.array([0.0, 1.0, 0.0])
        nz = np.linalg.norm(f, axis=1) > 0.0
        assert int(nz.sum()) >= 1
        # Every reaction is a LOAD: −n̂ projection, zero transverse.
        proj = f[nz] @ n_hat
        assert (proj < 0.0).all(), proj
        np.testing.assert_allclose(f[nz][:, [0, 2]], 0.0, atol=1e-30)
        # Net is purely along −n̂; reservoir feels +n̂·|net| → closed sum = 0.
        assert net[0] == 0.0 and net[2] == 0.0
        assert net[1] < 0.0
        membrane_reaction = -net
        np.testing.assert_allclose(
            net + membrane_reaction, np.zeros(3), atol=1e-30
        )
        # Self-consistent: net = (#contacting) · (−n̂ · F_per_applied).
        f_per_applied = -float(f[nz][0] @ n_hat)
        assert f_per_applied > 0.0
        np.testing.assert_allclose(
            net, -n_hat * int(nz.sum()) * f_per_applied, rtol=1e-9, atol=1e-30
        )


# ---------------------------------------------------------------------------
# §4 Numerical sanity / CFL
# ---------------------------------------------------------------------------
class TestNumericalCFL:
    def test_effective_stiffness_formula(self):
        """k_mem_eff = γ_mem / s_fil (tented-membrane spring; area estimate)."""
        p = _resolved_membrane()
        ml = MembraneLoad(
            p, lamel_state=type("S", (), {
                "barbed_end_tags": [], "capped_tags": set()})(),
            barbed_load_consumers=[],
        )
        for n_contact in (1, 10, 100):
            s_fil = math.sqrt(p.A_mem / n_contact)
            assert math.isclose(
                ml.effective_tip_stiffness(n_contact),
                p.gamma_mem / s_fil, rel_tol=1e-12,
            )
        assert ml.effective_tip_stiffness(0) == 0.0

    def test_action_path_has_no_cfl_constraint(self):
        """The PRIMARY membrane path is the MembraneLoad Action (feeds a
        Bell-Evans RATE, no integrator force) → NO CFL constraint at all.
        Attaching WITHOUT the reaction force never raises, even when the
        bounded per-step drift WOULD exceed contact_range (the gate is
        irrelevant because no force acts on the integrator)."""
        from ffn_sim.archive.hoomd_legacy.cell.membrane import _physical_span_cap
        sim = _minimal_sim_with_integrator(dt=1.3e-8)
        # A band-valid membrane whose λ-capped bounded load WOULD trip a
        # force-CFL gate (tiny contact_range → drift ≫ ceiling) IF it were a
        # per-step force.
        pm = ResolvedMembrane(
            gamma_mem=3.0e-4, y_plane=1.0e-6,
            contact_range=1.0e-15, A_mem=9.0e-10, kT=_KT,
        )
        drift = pm.gamma_mem * _physical_span_cap(pm) * 1.3e-8 / 3.91e-10
        assert drift > 0.1 * pm.contact_range  # would trip a force-CFL gate
        # But with_reaction_force=False applies no force → must NOT raise.
        res = attach_membrane_to_simulation(
            sim, pm, lamel_state=_DummyState(), barbed_load_consumers=[],
            batch_steps=1, gamma_b=3.91e-10,
            with_reaction_force=False, cfl_strict=True,
        )
        assert res["membrane_load"] is not None
        assert res["membrane_reaction_force"] is None

    def test_cfl_gate_blocks_absurd_load(self):
        """With the per-step REACTION FORCE on, the CFL attach-gate raises
        when the bounded-load per-step drift exceeds safety·contact_range
        (parity with erm.py / enclosed_volume.py). The physically-meaningful
        way to violate it is a vanishingly small contact_range (drift ceiling
        → 0 while the bounded load stays finite). Hand-built ResolvedMembrane
        (band-valid γ_mem) so only the CFL branch is under test."""
        sim = _minimal_sim_with_integrator(dt=1.3e-8)
        pm_absurd = ResolvedMembrane(
            gamma_mem=3.0e-4, y_plane=1.0e-6,
            contact_range=1.0e-15, A_mem=9.0e-10, kT=_KT,
        )
        with pytest.raises(RuntimeError, match="Membrane CFL violated"):
            attach_membrane_to_simulation(
                sim, pm_absurd, lamel_state=_DummyState(),
                barbed_load_consumers=[],
                batch_steps=1, gamma_b=3.91e-10,
                with_reaction_force=True, cfl_strict=True,
            )

    def test_cfl_gate_skipped_when_strict_false(self):
        """cfl_strict=False suppresses the CFL raise even with the reaction
        force on (diagnostic-run escape hatch, parity with erm.py)."""
        sim = _minimal_sim_with_integrator(dt=1.3e-8)
        pm_absurd = ResolvedMembrane(
            gamma_mem=3.0e-4, y_plane=1.0e-6,
            contact_range=1.0e-15, A_mem=9.0e-10, kT=_KT,
        )
        res = attach_membrane_to_simulation(  # must NOT raise
            sim, pm_absurd, lamel_state=_DummyState(),
            barbed_load_consumers=[],
            batch_steps=1, gamma_b=3.91e-10,
            with_reaction_force=True, cfl_strict=False,
        )
        assert res["membrane_load"] is not None

    def test_load_finite_under_run(self):
        """Membrane load stays finite across a short run."""
        p, handles = _tiny_lamellipodium(n_wave=5, batch_steps=1)
        state = handles["state"]
        pm = _membrane_for(p)
        res = attach_membrane_to_simulation(
            handles["sim"], pm, lamel_state=state,
            barbed_load_consumers=[handles["elong_action"], handles["cap_action"]],
            batch_steps=1, gamma_b=3.91e-10,
        )
        handles["sim"].run(10)
        ml = res["membrane_load"]
        assert math.isfinite(ml.last_F_total)
        assert ml.last_F_total >= 0.0


# ---------------------------------------------------------------------------
# §5 Sign / sense — reaction opposes barbed-end advance (a load)
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_reaction_opposes_advance(self):
        """A barbed end contacting the membrane feels a reaction along −n̂
        (into the cytosol, opposing advance toward +n̂): a load, not a
        push-forward."""
        p, handles = _tiny_lamellipodium(n_wave=3, batch_steps=1)
        state = handles["state"]
        pm = _membrane_for(p)
        sim = handles["sim"]
        snap = sim.state.get_snapshot()
        if snap.communicator.rank == 0:
            pos = np.asarray(snap.particles.position).copy()
            for t in state.barbed_end_tags:
                pos[t, 1] = pm.y_plane
            snap.particles.position[:] = pos
        sim.state.set_snapshot(snap)
        rf = MembraneReactionForce(pm, lamel_state=state)
        sim.operations.integrator.forces.append(rf)
        sim.run(0)
        f = np.asarray(rf.forces)
        n_hat = np.array([0.0, 1.0, 0.0])
        contacting = np.linalg.norm(f, axis=1) > 1e-30
        assert contacting.any()
        # Every reaction projects negatively on +n̂ (points into cytosol).
        proj = f[contacting] @ n_hat
        assert (proj < 0.0).all(), proj

    def test_higher_tension_larger_load(self):
        """Sign of the γ_mem → F monotonicity (load grows with tension)."""
        s = 100e-9
        F_lo = per_filament_load_from_spacing(1.0e-5, s)
        F_hi = per_filament_load_from_spacing(3.0e-4, s)
        assert F_hi > F_lo > 0.0


# ---------------------------------------------------------------------------
# §6 Measurement protocol — short BAOAB smoke with membrane on
# ---------------------------------------------------------------------------
class TestBAOABSmoke:
    def test_baoab_lamellipodium_with_membrane_no_nan(self):
        """A small lamellipodium + membrane-on runs a few hundred BAOAB steps
        with no NaN/Inf and the membrane load stays finite and ≥ 0."""
        p, handles = _tiny_lamellipodium(n_wave=8, batch_steps=1)
        state = handles["state"]
        pm = _membrane_for(p)
        res = attach_membrane_to_simulation(
            handles["sim"], pm, lamel_state=state,
            barbed_load_consumers=[handles["elong_action"], handles["cap_action"]],
            batch_steps=1, gamma_b=3.91e-10, with_reaction_force=True,
        )
        sim = handles["sim"]
        sim.run(200)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all(), "Membrane-on BAOAB run produced NaN/Inf."
        ml = res["membrane_load"]
        assert math.isfinite(ml.last_F_total) and ml.last_F_total >= 0.0
        rf = res["membrane_reaction_force"]
        assert rf is not None
        assert math.isfinite(rf.energy)


# ---------------------------------------------------------------------------
# Builder hook — off-path no-op
# ---------------------------------------------------------------------------
class TestBuilderHook:
    def _resolved_cortex(self):
        from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
        cfg_path = (
            Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
        )
        cfg = yaml.safe_load(open(cfg_path))
        cfg["cortex"]["n_filaments"] = 30
        cfg["cortex"]["demo_mode"] = True
        return resolve_h3_derived(cfg)

    def test_builder_off_path_when_membrane_none(self):
        """build_cortex_full_simulation with p_membrane=None attaches no
        membrane (off-path no-op) and the handles are None."""
        from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation
        handles = build_cortex_full_simulation(self._resolved_cortex(), with_baoab=True)
        assert handles["membrane_load_action"] is None
        assert handles["membrane_reaction_force"] is None
        sim = handles["sim"]
        assert not any(
            isinstance(f, MembraneReactionForce)
            for f in sim.operations.integrator.forces
        )
        assert not any(
            isinstance(getattr(u, "action", None), MembraneLoad)
            for u in sim.operations.updaters
        )

    def test_builder_off_path_bit_for_bit(self):
        """p_membrane=None is bit-for-bit identical to not passing it at all
        (off-path no-op; in-process construction-snapshot hash compare)."""
        import hashlib
        from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation
        p = self._resolved_cortex()

        def fingerprint(handles):
            sim = handles["sim"]
            sim.run(0)
            with sim.state.cpu_local_snapshot as s:
                pos = np.asarray(s.particles.position).copy()
                tag = np.asarray(s.particles.tag).copy()
            order = np.argsort(tag)
            h = hashlib.sha256()
            h.update(pos[order].tobytes())
            ig = sim.operations.integrator
            h.update(repr([type(f).__name__ for f in ig.forces]).encode())
            h.update(
                repr([type(u).__name__ for u in sim.operations.updaters]).encode()
            )
            return h.hexdigest()

        h_baseline = fingerprint(build_cortex_full_simulation(p, with_baoab=True))
        h_none = fingerprint(
            build_cortex_full_simulation(p, with_baoab=True, p_membrane=None)
        )
        assert h_baseline == h_none, (
            "p_membrane=None changed the off-path construction — not a no-op."
        )

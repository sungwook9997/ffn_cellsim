"""Tests for cortex/osmotic_regulation.py — dynamic osmotic / volume
regulation (RVD / RVI water-flux setpoint Updater on the enclosed-volume force).

Import-light + fast: resolver-level, OFF-identity, dimensional/sign sanity, and
a bare-Action ``act`` step that needs NO full Cell and NO long sim. The single
test that wants a built simulation is skipped (heavy/HOOMD-state).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ffn_sim.cortex.enclosed_volume import (
    EnclosedVolumePressure,
    resolve_enclosed_volume,
)
from ffn_sim.cortex.osmotic_regulation import (
    GAMMA_DENYLIST_PREFIX,
    PI_DECISIONS,
    R_GAS,
    OsmoticRegulationUpdater,
    ResolvedOsmoticRegulation,
    attach_osmotic_regulation_to_simulation,
    resolve_osmotic_regulation,
    vant_hoff_pressure,
    water_flux_volume_step,
)

R_CELL = 7.5e-6  # m   MCF7 radius (Wagner 2011)
DT = 1.0e-8      # s   representative BAOAB timestep


def _ev(turgor: float = 40.0):
    """Resolve a host enclosed-volume force at the MCF7 operating point."""
    return resolve_enclosed_volume(
        {"turgor_dP0": turgor}, R_cell=R_CELL
    )


# ---------------------------------------------------------------------------
# OFF-identity
# ---------------------------------------------------------------------------
class TestOffIdentity:
    def test_resolve_disabled_by_default(self):
        """Missing block => enabled=False, zeroed transport coefficients."""
        p = resolve_osmotic_regulation({}, R_cell=R_CELL, dt=DT, p_enclosed_volume=_ev())
        assert p.enabled is False
        assert p.Lp == 0.0
        assert p.A_mem == 0.0
        assert p.batch_steps == 0
        assert p.batch_dt == 0.0
        assert p.tau_RVD == 0.0

    def test_resolve_explicit_false(self):
        p = resolve_osmotic_regulation(
            {"enabled": False, "Lp": 1e-13}, R_cell=R_CELL, dt=DT,
            p_enclosed_volume=_ev(),
        )
        assert p.enabled is False
        assert p.Lp == 0.0

    def test_builder_disabled_is_noop(self):
        """Disabled config => builder attaches nothing, returns (None, None)."""
        p_off = resolve_osmotic_regulation({}, R_cell=R_CELL, dt=DT, p_enclosed_volume=_ev())
        ev = EnclosedVolumePressure(_ev(), (0, 100))
        # No real sim needed: the disabled path must early-return before
        # touching sim at all. Pass a deliberately-unusable sentinel to prove it.
        action, updater = attach_osmotic_regulation_to_simulation(
            sim=object(), p=p_off, ev_force=ev  # type: ignore[arg-type]
        )
        assert action is None
        assert updater is None
        # The enclosed-volume setpoint is untouched by an OFF build.
        assert ev.p.V0 == _ev().V0

    def test_builder_noop_when_ev_force_none(self):
        p_on = resolve_osmotic_regulation(
            {"enabled": True}, R_cell=R_CELL, dt=DT, p_enclosed_volume=_ev()
        )
        action, updater = attach_osmotic_regulation_to_simulation(
            sim=object(), p=p_on, ev_force=None  # type: ignore[arg-type]
        )
        assert action is None
        assert updater is None


# ---------------------------------------------------------------------------
# No-gamma-contamination (this module adds NO bonds → empty prefix)
# ---------------------------------------------------------------------------
class TestGammaDenylist:
    def test_prefix_is_empty_no_bonds(self):
        """This compartment modulates an existing setpoint; it adds no bonds.

        The contract for a bond-adding compartment is a non-empty prefix; for a
        no-bond compartment the prefix is '' and the builder must add zero bond
        types. Both are satisfied here (the Updater touches no bonds array).
        """
        assert GAMMA_DENYLIST_PREFIX == ""

    def test_pi_decisions_present_and_typed(self):
        assert isinstance(PI_DECISIONS, list)
        # Lp has no MCF7-specific datum → an open decision must be recorded.
        assert any("Lp" in s for s in PI_DECISIONS)


# ---------------------------------------------------------------------------
# Dimensional / parameter sanity
# ---------------------------------------------------------------------------
class TestResolverSanity:
    def test_vant_hoff_dimensional(self):
        """Π = R_gas·T·Δc [Pa]; check against a hand value."""
        # 1 mOsm/L = 1 mol/m³; at 310.15 K, Π = 8.314·310.15·1 ≈ 2579 Pa.
        Pi = vant_hoff_pressure(1.0, 310.15)
        assert math.isclose(Pi, R_GAS * 310.15 * 1.0, rel_tol=1e-12)
        assert 2500.0 < Pi < 2650.0

    def test_tau_rvd_is_reference_osmotic_timescale(self):
        """τ_RVD = V0/(Lp·A·Π_osm) is the REFERENCE osmotic-osmometer timescale.

        Π_osm is NOT the stiffness of the integrated trajectory (that is K_vol —
        see test_tau_kvol_is_integrated_mode). τ_RVD is retained only as a
        labelled reference / Hoffmann-band overlay.
        """
        ev = _ev()
        p = resolve_osmotic_regulation(
            {"enabled": True, "Lp": 1e-13, "batch_steps": 1000},
            R_cell=R_CELL, dt=DT, p_enclosed_volume=ev,
        )
        assert math.isclose(p.Pi_osm, vant_hoff_pressure(p.c_phys, p.temperature_K),
                            rel_tol=1e-12)
        expected = ev.V0 / (p.Lp * p.A_mem * p.Pi_osm)
        assert math.isclose(p.tau_RVD, expected, rel_tol=1e-12)
        assert p.tau_RVD > 0.0
        # Π_osm (~7.7e5 Pa) is far larger than K_vol (~1.3e3 Pa); the two
        # timescales therefore differ by that ratio (~580×).
        assert p.Pi_osm > 100.0 * ev.K_vol

    def test_tau_kvol_is_integrated_mode_and_gates_cfl(self):
        """τ_Kvol = V0/(Lp·A·K_vol) is the mode the updater actually integrates.

        The water law drives dV0/dt = -Lp·A·(dP_mech - dP_target) with
        dP_mech = Π₀ - K_vol·(V-V0)/V0, so the only V0-dependent restoring
        stiffness is the mechanical K_vol. τ_Kvol must equal V0/(Lp·A·K_vol)
        and must be the (much larger) integrated-mode timescale, ~580× the
        reference τ_RVD here.
        """
        ev = _ev()
        p = resolve_osmotic_regulation(
            {"enabled": True, "Lp": 1e-13, "batch_steps": 1000},
            R_cell=R_CELL, dt=DT, p_enclosed_volume=ev,
        )
        expected = ev.V0 / (p.Lp * p.A_mem * p.K_vol)
        assert math.isclose(p.tau_Kvol, expected, rel_tol=1e-12)
        assert p.tau_Kvol > 0.0
        # Integrated-mode timescale is K_vol-governed, ~580× the osmotic ref.
        assert math.isclose(p.tau_Kvol / p.tau_RVD, p.Pi_osm / p.K_vol,
                            rel_tol=1e-9)

    def test_tau_rvd_in_seconds_to_minutes_band(self):
        """Hoffmann 2009 RVD/RVI relaxation is seconds-to-minutes.

        With Lp in the Olbrich band and the MCF7 geometry, the REFERENCE
        osmotic τ_RVD must land in that physiological window (Hoffmann-band
        overlay consistency check, not a fit).
        """
        ev = _ev()
        p = resolve_osmotic_regulation(
            {"enabled": True, "Lp": 1e-13, "batch_steps": 1000},
            R_cell=R_CELL, dt=DT, p_enclosed_volume=ev,
        )
        assert 1.0 <= p.tau_RVD <= 600.0  # 1 s .. 10 min

    def test_default_target_is_resting_turgor(self):
        """No Δc / dP_target => relax toward the resting turgor (no perturb)."""
        ev = _ev(turgor=40.0)
        p = resolve_osmotic_regulation(
            {"enabled": True}, R_cell=R_CELL, dt=DT, p_enclosed_volume=ev
        )
        assert math.isclose(p.dP_target, 40.0, rel_tol=1e-12)

    def test_delta_c_sets_target_via_vant_hoff(self):
        """A hypotonic Δc raises the osmotic target above resting turgor."""
        ev = _ev(turgor=40.0)
        p = resolve_osmotic_regulation(
            {"enabled": True, "delta_c": 0.5}, R_cell=R_CELL, dt=DT,
            p_enclosed_volume=ev,
        )
        expected = 40.0 + vant_hoff_pressure(0.5, p.temperature_K)
        assert math.isclose(p.dP_target, expected, rel_tol=1e-12)
        assert p.dP_target > 40.0

    @pytest.mark.parametrize(
        "bad",
        [
            {"enabled": True, "Lp": 0.0},
            {"enabled": True, "Lp": -1e-13},
            {"enabled": True, "A_mem": 0.0},
            {"enabled": True, "A_mem": -1.0},
            {"enabled": True, "batch_steps": 0},
            {"enabled": True, "batch_steps": -5},
            {"enabled": True, "temperature_K": 0.0},
            {"enabled": True, "V0_min": -1.0},
        ],
    )
    def test_negative_or_zero_params_raise(self, bad):
        with pytest.raises(ValueError):
            resolve_osmotic_regulation(
                bad, R_cell=R_CELL, dt=DT, p_enclosed_volume=_ev()
            )

    def test_v0min_must_be_below_v0ref(self):
        ev = _ev()
        with pytest.raises(ValueError, match="V0_min"):
            resolve_osmotic_regulation(
                {"enabled": True, "V0_min": 2.0 * ev.V0},
                R_cell=R_CELL, dt=DT, p_enclosed_volume=ev,
            )

    def test_slow_mode_cfl_raises_when_batch_too_coarse(self):
        """batch_dt ≥ τ_Kvol must raise (integrated-mode explicit-Euler bound)."""
        ev = _ev()
        # Force τ_Kvol tiny by a huge Lp so any sane batch overshoots it. The
        # gate now anchors on the integrated (K_vol) mode, not the osmotic ref.
        with pytest.raises(ValueError, match="slow-mode CFL.*τ_Kvol"):
            resolve_osmotic_regulation(
                {"enabled": True, "Lp": 1.0, "batch_steps": 10**9},
                R_cell=R_CELL, dt=DT, p_enclosed_volume=ev,
            )

    def test_slow_mode_cfl_inert_without_host_force(self):
        """No enclosed-volume host => K_vol=0 => τ_Kvol=0 => CFL gate inert.

        Without a restoring stiffness there is no integrated mode to
        destabilise, so even an absurd batch must NOT raise (finding #7: the
        gate is the integrated-mode anchor, not the osmotic τ_RVD)."""
        p = resolve_osmotic_regulation(
            {"enabled": True, "Lp": 1.0, "batch_steps": 10**9},
            R_cell=R_CELL, dt=DT, p_enclosed_volume=None,
        )
        assert p.K_vol == 0.0
        assert p.tau_Kvol == 0.0
        # τ_RVD is still a non-zero, finite reference timescale.
        assert p.tau_RVD > 0.0
        assert math.isfinite(p.tau_RVD)


# ---------------------------------------------------------------------------
# Sign / sense of the water-flux setpoint step (the force-of-record motion)
# ---------------------------------------------------------------------------
class TestSignSense:
    def _resolved(self, target_dP: float = 40.0):
        ev = _ev(turgor=target_dP)
        return resolve_osmotic_regulation(
            {"enabled": True, "Lp": 1e-13, "batch_steps": 1000,
             "dP_target": target_dP},
            R_cell=R_CELL, dt=DT, p_enclosed_volume=ev,
        ), ev

    def test_cortex_pushes_harder_water_leaves_v0_shrinks(self):
        """ΔP_mech > ΔP_target => water leaves => V0 decreases (RVD)."""
        p, ev = self._resolved(target_dP=40.0)
        V0_new = water_flux_volume_step(p, dP_mech=80.0, V0=ev.V0)
        assert V0_new < ev.V0

    def test_cortex_pushes_less_water_enters_v0_grows(self):
        """ΔP_mech < ΔP_target => water enters => V0 increases (RVI)."""
        p, ev = self._resolved(target_dP=40.0)
        V0_new = water_flux_volume_step(p, dP_mech=10.0, V0=ev.V0)
        assert V0_new > ev.V0

    def test_at_target_is_fixed_point(self):
        """ΔP_mech == ΔP_target => dV/dt = 0 => V0 unchanged."""
        p, ev = self._resolved(target_dP=40.0)
        V0_new = water_flux_volume_step(p, dP_mech=40.0, V0=ev.V0)
        assert math.isclose(V0_new, ev.V0, rel_tol=1e-15)

    def test_volume_floor_clamps(self):
        """A huge outward water loss cannot drive V0 to/through zero."""
        p, ev = self._resolved(target_dP=40.0)
        # Step with an absurd mechanical overshoot AND a large batch_dt proxy
        # by calling many times; ensure it never crosses the floor.
        V0 = ev.V0
        for _ in range(10_000):
            V0 = water_flux_volume_step(p, dP_mech=1.0e9, V0=V0)
        assert V0 >= p.V0_min
        assert V0 > 0.0

    def test_euler_step_matches_closed_form(self):
        """One explicit-Euler step reproduces the closed-form ΔV exactly."""
        p, ev = self._resolved(target_dP=40.0)
        dP_mech = 70.0
        V0_new = water_flux_volume_step(p, dP_mech=dP_mech, V0=ev.V0)
        dVdt = -p.Lp * p.A_mem * (dP_mech - p.dP_target)
        expected = ev.V0 + dVdt * p.batch_dt
        assert math.isclose(V0_new, expected, rel_tol=1e-15)


# ---------------------------------------------------------------------------
# Updater act() — bare Action, no sim required
# ---------------------------------------------------------------------------
class TestUpdaterAct:
    def test_act_steps_setpoint_no_particle_motion(self):
        """act() mutates V0 on the live force and moves NO particle/bond."""
        ev_cfg = _ev(turgor=40.0)
        ev_force = EnclosedVolumePressure(ev_cfg, (0, 100))
        p = resolve_osmotic_regulation(
            {"enabled": True, "Lp": 1e-13, "batch_steps": 1000,
             "dP_target": 40.0},
            R_cell=R_CELL, dt=DT, p_enclosed_volume=ev_cfg,
        )
        action = OsmoticRegulationUpdater(p, ev_force)
        # Simulate that set_forces ran and the cortex is over-pressured.
        ev_force.last_pressure = 90.0
        V0_before = ev_force.p.V0
        action.act(timestep=1000)
        # Cortex over target => water leaves => V0 shrinks.
        assert ev_force.p.V0 < V0_before
        assert action.n_ticks == 1
        assert math.isclose(action.last_V0, ev_force.p.V0, rel_tol=1e-15)
        assert math.isclose(action.last_dP_mech, 90.0, rel_tol=1e-15)

    def test_act_firsttick_nan_pressure_falls_back_to_turgor(self):
        """Before any set_forces (last_pressure=NaN) act() uses resting turgor."""
        ev_cfg = _ev(turgor=40.0)
        ev_force = EnclosedVolumePressure(ev_cfg, (0, 100))
        # last_pressure defaults to NaN until set_forces runs.
        assert math.isnan(ev_force.last_pressure)
        p = resolve_osmotic_regulation(
            {"enabled": True, "Lp": 1e-13, "batch_steps": 1000,
             "dP_target": 40.0},
            R_cell=R_CELL, dt=DT, p_enclosed_volume=ev_cfg,
        )
        action = OsmoticRegulationUpdater(p, ev_force)
        action.act(timestep=1000)
        # turgor (40) == target (40) => fixed point => V0 essentially unchanged.
        assert math.isclose(ev_force.p.V0, ev_cfg.V0, rel_tol=1e-12)
        assert math.isclose(action.last_dP_mech, 40.0, rel_tol=1e-15)

    def test_disabled_config_cannot_build_updater(self):
        ev_cfg = _ev()
        ev_force = EnclosedVolumePressure(ev_cfg, (0, 100))
        p_off = ResolvedOsmoticRegulation(
            enabled=False, Lp=0.0, A_mem=0.0, batch_steps=0, dt=DT,
            batch_dt=0.0, dP_target=0.0, temperature_K=310.15, V0_min=0.0,
        )
        with pytest.raises(ValueError, match="disabled"):
            OsmoticRegulationUpdater(p_off, ev_force)

    def test_updater_reasserts_cfl_on_hand_built_dataclass(self):
        """Finding #3: a hand-built / post-mutated enabled resolved dataclass
        that bypasses the resolver's gate (batch_dt ≥ τ_Kvol) must still be
        rejected at Updater construction."""
        ev_cfg = _ev()
        ev_force = EnclosedVolumePressure(ev_cfg, (0, 100))
        # Hand-build an enabled config with a finite τ_Kvol but batch_dt above
        # it (the resolver never saw this — it is constructed directly).
        # batch_steps chosen so batch_dt = steps·dt genuinely exceeds τ_Kvol
        # (τ_Kvol = V0/(Lp·A·K_vol) ≈ 1.9e4 s here; dt = 1e-8 s ⇒ need
        # steps ≳ 1.9e12). 3e12 ⇒ batch_dt = 3e4 s > τ_Kvol.
        bad_steps = 3 * 10**12
        p_bad = ResolvedOsmoticRegulation(
            enabled=True, Lp=1e-13, A_mem=4.0 * math.pi * R_CELL**2,
            batch_steps=bad_steps, dt=DT, batch_dt=bad_steps * DT,
            dP_target=40.0, temperature_K=310.15, V0_min=0.01 * ev_cfg.V0,
            K_vol=ev_cfg.K_vol, V0_ref=ev_cfg.V0,
            tau_Kvol=ev_cfg.V0 / (1e-13 * 4.0 * math.pi * R_CELL**2 * ev_cfg.K_vol),
        )
        assert p_bad.batch_dt >= p_bad.tau_Kvol
        with pytest.raises(ValueError, match="slow-mode CFL.*τ_Kvol"):
            OsmoticRegulationUpdater(p_bad, ev_force)


# ---------------------------------------------------------------------------
# Heavy: full-sim integration (skipped — needs a built cortex + integrator)
# ---------------------------------------------------------------------------
@pytest.mark.skip(reason="heavy: needs a built cortex sim + integrator")
def test_attach_to_built_simulation_relaxes_setpoint():
    """Smoke: over many ticks V0(t) relaxes toward the target volume."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# LIVE activation wiring (2026-06-09; PI 소유권 허용): manifest + post-build
# attach + registry LIVE. Full-build tests (~0.5 s each).
# ---------------------------------------------------------------------------
def test_osmotic_off_build_is_bit_identity():
    """Baseline (osmotic absent) => p_osmotic_regulation None; build unchanged."""
    from ffn_sim.cell.manifest import (
        build_baseline_cell, load_manifest, resolve_baseline,
    )
    rb = resolve_baseline(load_manifest("mcf7_baseline.yaml"))
    assert rb.p_osmotic_regulation is None
    cell = build_baseline_cell("mcf7_baseline.yaml", seed=1)
    assert cell.simulation.state.N_particles > 0


def test_osmotic_on_attaches_updater_without_contamination():
    """osmotic_rvd => +1 updater, tau_RVD in band, IDENTICAL bond inventory."""
    from ffn_sim.cell.compartment_registry import REGISTRY, load_recipe
    from ffn_sim.cell.manifest import (
        build_baseline_cell, load_manifest, resolve_baseline,
    )
    base = load_manifest("mcf7_baseline.yaml")
    manifest, deferred = REGISTRY.compose_manifest(
        load_recipe("osmotic_rvd"), base_manifest=base, strict=True
    )
    assert deferred == []
    rb = resolve_baseline(manifest)
    assert rb.p_osmotic_regulation is not None and rb.p_osmotic_regulation.enabled
    # PRIMARY gate (analytic): tau_RVD in the Hoffmann 2009 seconds-minutes band.
    assert 3.0 <= rb.p_osmotic_regulation.tau_RVD <= 600.0

    cell_off = build_baseline_cell("mcf7_baseline.yaml", seed=1)
    cell_on = build_baseline_cell("mcf7_baseline.yaml", manifest=manifest, seed=1)
    n_off = len(cell_off.simulation.operations.updaters)
    n_on = len(cell_on.simulation.operations.updaters)
    assert n_on == n_off + 1, (n_off, n_on)   # exactly the osmotic updater added

    s_off = cell_off.simulation.state.get_snapshot()
    s_on = cell_on.simulation.state.get_snapshot()
    # NO gamma contamination: osmotic adds ZERO bonds => identical bond inventory.
    assert int(s_on.bonds.N) == int(s_off.bonds.N)
    assert sorted(s_on.bonds.types) == sorted(s_off.bonds.types)
    assert int(s_on.particles.N) == int(s_off.particles.N)

"""STATIC + demo sanity-gate tests for H.3 + H.5 Cell composition (v2).

Covers ``ffn_sim/cell/cell.py`` Sanity Gate §1–6:

- §2 Boundary cases               (TestBoundary)
- §3 Conservation invariants      (TestCountsAndTagRanges)
- §4 Numerical sanity             (TestERMCFLPropagation)
- §6 Measurement protocol         (TestDiagnostics + demo run)

H.5 단계 2 lamellipodium wiring tests live in
:class:`TestLamellipodiumWiring` — mirrors the with_myosin / with_erm
patterns.  Boundary case (missing ``p_lamellipodium`` raises), additive
particle / bond / type wiring, and the full-assembly compatibility check
(``with_lamellipodium=True`` × ``with_myosin=True`` × ``with_xlinks=True``
× ``with_erm=True``).
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.archive.hoomd_legacy.cell import Cell, CellBuildOptions, resolve_h5_lamellipodium
from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.archive.hoomd_legacy.cortex.erm import ResolvedERM, resolve_erm


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo_cfg(n_filaments: int = 50) -> dict:
    cfg = deepcopy(_load_cfg())
    cfg["cortex"]["n_filaments"] = n_filaments
    cfg["cortex"]["demo_mode"] = True
    return cfg


@pytest.fixture(scope="module")
def p_cortex():
    return resolve_h3_derived(_demo_cfg())


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_bare_cortex_no_addons(self, p_cortex):
        """Cell with no add-ons builds equivalent to bare cortex."""
        cell = Cell.build(p_cortex)
        assert cell.n_particles == p_cortex.n_filaments * p_cortex.beads_per_filament
        assert cell.n_xlink_heads == 0
        assert cell.erm_force is None
        assert cell.xlink_action is None

    def test_with_erm_requires_p_erm(self, p_cortex):
        opts = CellBuildOptions(with_erm=True)
        with pytest.raises(ValueError, match="requires p_erm"):
            Cell.build(p_cortex, options=opts)

    def test_with_erm_propagates_cfl_violation(self, p_cortex):
        """Cell.build propagates ERM CFL violation when caller passes
        an old-style k_ERM that violates CFL at cortex dt_cfl.

        KU-3.18 RE-RATIFIED 2026-05-26 to k_ERM=1e-4 N/m (CFL-safe).
        The yaml default no longer raises — but the gate is still
        guarded against caller-supplied old values (e.g. legacy 0.1 N/m
        from external configs).
        """
        import math
        # Construct an ERM with the OLD 0.1 N/m value (NOT from yaml).
        p_erm_old = ResolvedERM(
            k_ERM=0.1, R_cell=p_cortex.R_cell, cell_center=(0.0, 0.0, 0.0),
        )
        p_erm_old.sigma_radial_thermal = math.sqrt(p_cortex.kT / p_erm_old.k_ERM)
        opts = CellBuildOptions(with_erm=True)
        with pytest.raises(RuntimeError, match="ERM CFL violated"):
            Cell.build(p_cortex, p_erm=p_erm_old, options=opts)

    def test_with_zero_xlinks_idle(self, p_cortex):
        """with_crosslinkers=True + n_xl=0 builds a bare cortex (no xlink
        layer attached). The Cell's xlink slots stay None."""
        cfg = _demo_cfg()
        cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = 0
        p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)
        opts = CellBuildOptions(with_crosslinkers=True)
        cell = Cell.build(p_cortex, p_xlinks=p_xl, options=opts)
        # n_xl=0 → xlink layer SKIPPED (matches the Cell.build branch).
        assert cell.xlink_action is None
        assert cell.n_xlink_heads == 0


# ---------------------------------------------------------------------------
# §3 Conservation — counts + tag ranges
# ---------------------------------------------------------------------------
class TestCountsAndTagRanges:
    def test_bead_count_summary_sums_to_particle_total(self, p_cortex):
        cell = Cell.build(p_cortex)
        summary = cell.bead_count_summary()
        assert sum(summary.values()) == cell.n_particles

    def test_tag_ranges_contiguous_non_overlapping(self, p_cortex):
        cell = Cell.build(p_cortex)
        ranges = cell.tag_ranges()
        prev_end = 0
        for name, (lo, hi) in ranges.items():
            assert lo == prev_end, (
                f"Tag range gap or overlap at {name}: {lo} != {prev_end}"
            )
            assert hi >= lo
            prev_end = hi
        assert prev_end == cell.n_particles

    def test_with_xlinks_particle_count(self, p_cortex):
        cfg = _demo_cfg()
        # Widen max_bind_dist for demo cortex (sparse)
        cfg["cortex"]["dynamic_crosslinkers"]["max_bind_dist"] = 1.0e-6
        cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = 20
        p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)
        opts = CellBuildOptions(with_crosslinkers=True)
        cell = Cell.build(p_cortex, p_xlinks=p_xl, options=opts)
        expected = p_cortex.n_filaments * p_cortex.beads_per_filament + 2 * 20
        assert cell.n_particles == expected
        assert cell.n_xlink_heads == 40


# ---------------------------------------------------------------------------
# §4 Numerical sanity — ERM CFL gate propagation
# ---------------------------------------------------------------------------
class TestERMCFLPropagation:
    def test_soft_ERM_builds_and_runs(self, p_cortex):
        """With a CFL-safe soft k_ERM (matches test_erm.py demo), Cell.build
        succeeds and the simulation runs without stiff-spring runaway."""
        import math
        soft_k_ERM = 1.0e-4
        p_erm_soft = ResolvedERM(
            k_ERM=soft_k_ERM, R_cell=p_cortex.R_cell,
        )
        p_erm_soft.sigma_radial_thermal = math.sqrt(p_cortex.kT / soft_k_ERM)
        opts = CellBuildOptions(with_erm=True)
        cell = Cell.build(p_cortex, p_erm=p_erm_soft, options=opts)
        cell.run(500)
        with cell.simulation.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all()


# ---------------------------------------------------------------------------
# §6 Measurement protocol — diagnostics + demo run
# ---------------------------------------------------------------------------
class TestDiagnostics:
    def test_diagnostics_dict_has_expected_keys(self, p_cortex):
        cell = Cell.build(p_cortex)
        d = cell.diagnostics()
        for key in (
            "particle_total", "bead_counts", "tag_ranges",
            "cortex_R_cell_m", "cortex_L_box_m", "dt_cfl_s",
            "with_erm", "with_crosslinkers", "with_baoab",
        ):
            assert key in d, f"diagnostics() missing key {key}"
        assert d["particle_total"] == cell.n_particles
        assert d["cortex_R_cell_m"] == p_cortex.R_cell

    def test_demo_bare_cortex_500_steps_no_NaN(self, p_cortex):
        """Bare cortex (no ERM / xl) BAOAB smoke."""
        cell = Cell.build(p_cortex)
        cell.run(500)
        with cell.simulation.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all()


# ---------------------------------------------------------------------------
# H.5 단계 2 — Lamellipodium composition wiring
# ---------------------------------------------------------------------------
H5_CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h5.yaml"
)


def _load_h5_cfg() -> dict:
    with open(H5_CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo_h5_cfg(n_WAVE: int = 10) -> dict:
    """Small WAVE count keeps the snapshot-extension overhead modest in CI."""
    cfg = deepcopy(_load_h5_cfg())
    cfg["lamellipodium"]["n_WAVE"] = n_WAVE
    return cfg


@pytest.fixture(scope="module")
def p_lamel(p_cortex):
    """Lamellipodium params resolved against the cortex dt_cfl + L_box so
    the shared BAOAB integrator dt matches."""
    cfg = _demo_h5_cfg()
    return resolve_h5_lamellipodium(
        cfg, L_box=p_cortex.L_box, dt=p_cortex.dt_cfl, kT=p_cortex.kT,
    )


class TestLamellipodiumWiring:
    """H.5 단계 2 — :meth:`Cell.build(with_lamellipodium=True)` adds WAVE
    + mother-actin beads, the WAVE membrane plane pin, and the three
    D1/D2 Bieling/Funk Updaters on top of an existing cortex sim."""

    def test_with_lamellipodium_requires_p_lamellipodium(self, p_cortex):
        opts = CellBuildOptions(with_lamellipodium=True)
        with pytest.raises(ValueError, match="requires p_lamellipodium"):
            Cell.build(p_cortex, options=opts)

    def test_cell_with_lamellipodium_alone_dispatch(self, p_cortex, p_lamel):
        """``with_lamellipodium=True`` alone (no other add-ons) wires the
        full lamellipodium stack on a bare cortex."""
        opts = CellBuildOptions(with_lamellipodium=True)
        cell = Cell.build(
            p_cortex, p_lamellipodium=p_lamel, options=opts,
        )

        # Lamellipodium handles populated.
        assert cell.lamellipodium_layout is not None
        assert cell.lamellipodium_state is not None
        assert cell.wave_pin_force is not None
        assert cell.elong_action is not None
        assert cell.elong_updater is not None
        assert cell.branch_action is not None
        assert cell.branch_updater is not None
        assert cell.cap_action is not None
        assert cell.cap_updater is not None

        # Per-WAVE bead bookkeeping: n_wave + n_mother seeds = 2·n_WAVE.
        assert cell.n_wave_particles == p_lamel.n_WAVE
        assert cell.n_lamellipodium_actin == p_lamel.n_WAVE
        expected = (
            p_cortex.n_filaments * p_cortex.beads_per_filament
            + 2 * p_lamel.n_WAVE
        )
        assert cell.n_particles == expected

        # New particle types appeared in the simulation state.
        ptypes = list(cell.simulation.state.particle_types)
        assert "wave_particle" in ptypes
        assert "actin_lamel" in ptypes

        # bead_count_summary sums to n_particles.
        summary = cell.bead_count_summary()
        assert sum(summary.values()) == cell.n_particles
        assert summary["wave_particle"] == p_lamel.n_WAVE
        assert summary["lamellipodium_actin"] == p_lamel.n_WAVE

        # Tag ranges contiguous and non-overlapping.
        ranges = cell.tag_ranges()
        prev_end = 0
        for name, (lo, hi) in ranges.items():
            assert lo == prev_end, f"tag range gap at {name}"
            assert hi >= lo
            prev_end = hi
        assert prev_end == cell.n_particles

        # Diagnostics flag propagated.
        d = cell.diagnostics()
        assert d["with_lamellipodium"] is True
        assert d["n_wave_particles"] == p_lamel.n_WAVE

    def test_cell_with_lamellipodium_short_baoab_no_NaN(self, p_cortex, p_lamel):
        """Short BAOAB window on cortex + lamellipodium — no NaN/Inf,
        WAVE pin keeps WAVE particles near the membrane plane, the three
        D2-batched Updaters tick at ``p_lamel.batch_steps`` cadence."""
        opts = CellBuildOptions(with_lamellipodium=True)
        cell = Cell.build(
            p_cortex, p_lamellipodium=p_lamel, options=opts,
        )
        n_steps = max(p_lamel.batch_steps * 5, 500)
        cell.run(n_steps)

        with cell.simulation.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
            tag = np.asarray(s.particles.tag).copy()
        assert np.isfinite(pos).all()

        # WAVE tags sit at [N_cortex, N_cortex + n_WAVE); their y-coordinate
        # should stay close to Y_max (within ~σ_radial_thermal of the pin).
        wave_lo, wave_hi = cell.tag_ranges()["wave_particle"]
        wave_rows = np.where((tag >= wave_lo) & (tag < wave_hi))[0]
        if wave_rows.size:
            ys = pos[wave_rows, 1]
            # σ_pin = √(kT / k_wave_pin); allow 10σ envelope.
            sigma_pin = float(np.sqrt(p_cortex.kT / p_lamel.k_wave_pin))
            assert np.all(np.abs(ys - p_lamel.Y_max) < 10.0 * sigma_pin + 1e-7)

        # Updater tick count = n_steps // batch_steps.
        expected_ticks = n_steps // p_lamel.batch_steps
        assert cell.elong_action.steps_run == expected_ticks
        assert cell.branch_action.steps_run == expected_ticks
        assert cell.cap_action.steps_run == expected_ticks

    def test_cell_full_assembly_with_all_subsystems(self, p_cortex, p_lamel):
        """Combined ``with_lamellipodium`` × ``with_myosin`` × ``with_xlinks``
        × ``with_erm`` builds the full Phase 1 mechanobiology stack on a
        single Cell.  Counts + tag ranges remain consistent."""
        # Local imports of the resolvers used only in this assembly test
        # to keep the rest of the file's import surface minimal.
        from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin

        cfg = _demo_cfg()
        cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = 10
        cfg["cortex"]["dynamic_crosslinkers"]["max_bind_dist"] = 1.0e-6
        cfg["cortex"]["myosin"]["n_motors_per_cell"] = 3

        p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)
        p_myo = resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl)
        p_erm = resolve_erm(cfg, R_cell=p_cortex.R_cell, kT=p_cortex.kT)

        opts = CellBuildOptions(
            with_crosslinkers=True,
            with_erm=True,
            with_myosin=True,
            with_lamellipodium=True,
        )
        cell = Cell.build(
            p_cortex,
            p_xlinks=p_xl,
            p_erm=p_erm,
            p_myosin=p_myo,
            p_lamellipodium=p_lamel,
            options=opts,
        )

        # All subsystem handles populated.
        assert cell.xlink_action is not None
        assert cell.erm_force is not None
        assert cell.myosin_action is not None
        assert cell.elong_action is not None

        # Particle bookkeeping: cortex + 2·n_xl + n_myosin·n_per_motor +
        # 2·n_WAVE.
        expected = (
            p_cortex.n_filaments * p_cortex.beads_per_filament
            + 2 * p_xl.n_xl
            + p_myo.n_motors_per_cell * p_myo.n_particles_per_motor
            + 2 * p_lamel.n_WAVE
        )
        assert cell.n_particles == expected

        # bead_count_summary sums to n_particles and tag_ranges
        # contiguous-and-non-overlapping (Sanity Gate §3 invariants).
        summary = cell.bead_count_summary()
        assert sum(summary.values()) == cell.n_particles
        ranges = cell.tag_ranges()
        prev_end = 0
        for name, (lo, hi) in ranges.items():
            assert lo == prev_end, f"tag range gap at {name}"
            assert hi >= lo
            prev_end = hi
        assert prev_end == cell.n_particles

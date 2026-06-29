"""Integration tests for build_cortex_full_simulation + Cell.build(with_myosin=True).

Covers the H.3 단계 5 unified cortex+xlinks+myosin composition path:

- §3 Conservation: particle / bond counts match the sum of subsystems.
- §4 Numerical sanity: demo BAOAB run for the 3-way sim, no NaN/Inf,
  no force runaway (image wrap ≤ 1).
- §6 Measurement: MyosinStepUpdater fires + xlink Updater fires on the
  3-way sim during the demo run.

Production-scale 3-way integration gates (60 s simulated, KU-3.x bands)
are skip-marked under ``H3_INTEGRATION_PRODUCTION=1`` per the same
multi-hour rationale as the KU-3.x skeletons.
"""

from __future__ import annotations

import math
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.archive.hoomd_legacy.cell import Cell, CellBuildOptions, build_cortex_full_simulation
from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo_3way_cfg() -> dict:
    """Small cortex + few xlinks + few motors — CI-friendly smoke."""
    cfg = deepcopy(_load_cfg())
    cfg["cortex"]["n_filaments"] = 40
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = 15
    cfg["cortex"]["dynamic_crosslinkers"]["max_bind_dist"] = 1.0e-6
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = 5
    return cfg


@pytest.fixture(scope="module")
def cfg_3way():
    return _demo_3way_cfg()


@pytest.fixture(scope="module")
def p_cortex(cfg_3way):
    return resolve_h3_derived(cfg_3way)


@pytest.fixture(scope="module")
def p_xl(p_cortex, cfg_3way):
    return resolve_crosslinkers(cfg_3way, dt=p_cortex.dt_cfl)


@pytest.fixture(scope="module")
def p_myo(p_cortex, cfg_3way):
    return resolve_cortex_myosin(cfg_3way, dt=p_cortex.dt_cfl)


# ---------------------------------------------------------------------------
# §3 Conservation — particle / bond counts
# ---------------------------------------------------------------------------
class TestThreeWayCounts:
    def test_helper_builds_with_all_three(self, p_cortex, p_xl, p_myo):
        h = build_cortex_full_simulation(
            p_cortex, p_xlinks=p_xl, p_myosin=p_myo, with_baoab=False,
        )
        sim = h["sim"]
        expected_particles = (
            h["n_cortex_actin"] + h["n_xlink_heads"] + h["n_myosin_particles"]
        )
        assert sim.state.N_particles == expected_particles

    def test_helper_builds_with_xlinks_only(self, p_cortex, p_xl):
        h = build_cortex_full_simulation(
            p_cortex, p_xlinks=p_xl, p_myosin=None, with_baoab=False,
        )
        assert h["myosin_layout"] is None
        assert h["myosin_action"] is None
        assert h["n_myosin_particles"] == 0
        assert h["n_xlink_heads"] == 2 * p_xl.n_xl

    def test_helper_builds_with_myosin_only(self, p_cortex, p_myo):
        h = build_cortex_full_simulation(
            p_cortex, p_xlinks=None, p_myosin=p_myo, with_baoab=False,
        )
        assert h["xlink_layout"] is None
        assert h["xlink_action"] is None
        assert h["n_xlink_heads"] == 0
        assert h["n_myosin_particles"] == (
            p_myo.n_motors_per_cell * p_myo.n_particles_per_motor
        )

    def test_helper_builds_bare_cortex(self, p_cortex):
        h = build_cortex_full_simulation(
            p_cortex, p_xlinks=None, p_myosin=None, with_baoab=False,
        )
        assert h["n_xlink_heads"] == 0
        assert h["n_myosin_particles"] == 0
        assert h["sim"].state.N_particles == h["n_cortex_actin"]


# ---------------------------------------------------------------------------
# §4 Numerical sanity — demo BAOAB no NaN/Inf, no runaway
# ---------------------------------------------------------------------------
class TestThreeWayBAOAB:
    def test_500_step_no_NaN(self, p_cortex, p_xl, p_myo):
        h = build_cortex_full_simulation(
            p_cortex, p_xlinks=p_xl, p_myosin=p_myo, with_baoab=True,
        )
        sim = h["sim"]
        sim.run(500)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
            image = np.asarray(s.particles.image).copy()
        assert np.isfinite(pos).all(), "NaN/Inf in positions after 3-way demo"
        # No particle should have wrapped images more than ±1 in 500 steps
        # on a 30 μm cube box.
        assert (np.abs(image) <= 1).all(), (
            f"Image wrap > 1 on 3-way demo (max |image| = "
            f"{int(np.max(np.abs(image)))}); possible runaway."
        )

    def test_myosin_updater_ticks_on_demo_run(self, p_cortex, p_xl, p_myo):
        """500 steps / batch_steps=100 → expect 5 Updater ticks for each
        of XlinkBondUpdater and MyosinStepUpdater (D2 batch CFL guarantees
        batch_steps is preserved by resolve_*)."""
        h = build_cortex_full_simulation(
            p_cortex, p_xlinks=p_xl, p_myosin=p_myo, with_baoab=True,
        )
        sim = h["sim"]
        sim.run(500)
        assert h["xlink_action"].steps_run == 5
        assert h["myosin_action"].steps_run == 5


# ---------------------------------------------------------------------------
# §6 Measurement — Cell.build dispatch path
# ---------------------------------------------------------------------------
class TestCellBuildWithMyosin:
    def test_cell_build_with_myosin_dispatch(self, p_cortex, p_xl, p_myo):
        """Cell.build(with_myosin=True) routes through the unified
        build_cortex_full_simulation and exposes myosin handles."""
        opts = CellBuildOptions(with_crosslinkers=True, with_myosin=True)
        cell = Cell.build(
            p_cortex, p_xlinks=p_xl, p_myosin=p_myo, options=opts,
        )
        assert cell.myosin_action is not None
        assert cell.myosin_updater is not None
        assert cell.myosin_layout is not None
        assert cell.xlink_action is not None
        assert cell.n_myosin_particles == (
            p_myo.n_motors_per_cell * p_myo.n_particles_per_motor
        )

    def test_cell_with_myosin_requires_p_myosin(self, p_cortex):
        opts = CellBuildOptions(with_myosin=True)
        with pytest.raises(ValueError, match="requires p_myosin"):
            Cell.build(p_cortex, options=opts)

    def test_cell_bead_count_summary_includes_myosin(self, p_cortex, p_xl, p_myo):
        opts = CellBuildOptions(with_crosslinkers=True, with_myosin=True)
        cell = Cell.build(
            p_cortex, p_xlinks=p_xl, p_myosin=p_myo, options=opts,
        )
        summary = cell.bead_count_summary()
        # Sums to particle total
        assert sum(summary.values()) == cell.n_particles
        # Myosin sub-counts
        assert summary["myosin_backbone"] == (
            p_myo.n_motors_per_cell * p_myo.n_backbone
        )
        assert summary["myosin_head"] == (
            p_myo.n_motors_per_cell * 2 * p_myo.n_heads_per_side
        )

    def test_cell_tag_ranges_contiguous_with_myosin(self, p_cortex, p_xl, p_myo):
        opts = CellBuildOptions(with_crosslinkers=True, with_myosin=True)
        cell = Cell.build(
            p_cortex, p_xlinks=p_xl, p_myosin=p_myo, options=opts,
        )
        ranges = cell.tag_ranges()
        prev_end = 0
        for name, (lo, hi) in ranges.items():
            assert lo == prev_end
            assert hi >= lo
            prev_end = hi
        # Ends at total particle count (cortex + xlink + myosin).
        # lamellipodium / fa zero-width at end → prev_end == cell.n_particles.
        assert prev_end == cell.n_particles

    def test_cell_diagnostics_reports_n_realised(self, p_cortex, p_xl, p_myo):
        opts = CellBuildOptions(with_crosslinkers=True, with_myosin=True)
        cell = Cell.build(
            p_cortex, p_xlinks=p_xl, p_myosin=p_myo, options=opts,
        )
        d = cell.diagnostics()
        assert d["with_myosin"] is True
        assert d["n_myosin_motors_realised"] == p_myo.n_motors_per_cell
        assert d["n_xlinks_realised"] == p_xl.n_xl


# ---------------------------------------------------------------------------
# Opt-in production-scale 3-way gate (60 s simulated, multi-hour wall)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not bool(int(os.environ.get("H3_INTEGRATION_PRODUCTION", "0"))),
    reason=(
        "Full 3-way 60 s production gate; opt-in via "
        "H3_INTEGRATION_PRODUCTION=1. Multi-hour wall (same rationale "
        "as KU-3.x production opt-ins). Also requires ERM CFL resolution."
    ),
)
class TestThreeWayProduction:
    def test_60s_simulation_completes_without_failure(self):
        pytest.skip("Production sign-off skeleton.")

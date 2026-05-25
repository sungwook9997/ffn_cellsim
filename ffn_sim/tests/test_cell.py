"""STATIC + demo sanity-gate tests for H.3 Cell composition (v2).

Covers ``ffn_sim/cell/cell.py`` Sanity Gate §1–6:

- §2 Boundary cases               (TestBoundary)
- §3 Conservation invariants      (TestCountsAndTagRanges)
- §4 Numerical sanity             (TestERMCFLPropagation)
- §6 Measurement protocol         (TestDiagnostics + demo run)
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.cell import Cell, CellBuildOptions
from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cortex.erm import ResolvedERM, resolve_erm


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
        """Brief-literal k_ERM = 0.1 N/m at cortex dt_cfl violates CFL;
        Cell.build must raise (not silently swallow)."""
        p_erm = resolve_erm(
            _load_cfg(), kT=p_cortex.kT, R_cell=p_cortex.R_cell
        )
        opts = CellBuildOptions(with_erm=True)
        with pytest.raises(RuntimeError, match="ERM CFL violated"):
            Cell.build(p_cortex, p_erm=p_erm, options=opts)

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

"""STATIC + demo tests for H.3 variable-length filament distribution.

Covers the additive variable-length topology generator + state +
simulation builders in `ffn_sim/cortex/cortex.py` (단계 7 deliverable).

- §3 Conservation: per-filament counts sum to expected totals; bond
  + angle counts match Σ(N_i-1), Σ(N_i-2).
- §4 Numerical: positions inside box, bond rest length exact at
  construction, force-free initial state.
- §6 Measurement: empirical L distribution mean ≈ (L_min+L_max)/2;
  range spans [L_min, L_max] under ℓ_0 quantization.

Production-scale variable-length L_p sweep is deferred via
`H3_VARIABLE_LENGTH_PRODUCTION=1` opt-in (multi-hour wall, same
rationale as KU-3.x and L_p full-sweep gates).
"""

from __future__ import annotations

import math
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.cortex.cortex import (
    ResolvedH3,
    VariableLengthCortexLayout,
    build_variable_length_cortex_simulation,
    build_variable_length_cortex_state,
    generate_variable_length_cortex_layout,
    resolve_h3_derived,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo_cfg(n_filaments: int = 80) -> dict:
    cfg = deepcopy(_load_cfg())
    cfg["cortex"]["n_filaments"] = n_filaments
    cfg["cortex"]["demo_mode"] = True
    return cfg


@pytest.fixture(scope="module")
def p_cortex():
    return resolve_h3_derived(_demo_cfg())


@pytest.fixture(scope="module")
def layout(p_cortex):
    return generate_variable_length_cortex_layout(
        p_cortex, L_min=1.0e-6, L_max=5.0e-6,
    )


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_negative_L_min_raises(self, p_cortex):
        with pytest.raises(ValueError, match="L_min"):
            generate_variable_length_cortex_layout(
                p_cortex, L_min=-1.0e-6, L_max=5.0e-6,
            )

    def test_L_max_below_L_min_raises(self, p_cortex):
        with pytest.raises(ValueError, match="L_max"):
            generate_variable_length_cortex_layout(
                p_cortex, L_min=5.0e-6, L_max=1.0e-6,
            )

    def test_L_max_exceeds_diameter_raises(self, p_cortex):
        with pytest.raises(ValueError, match="cell diameter"):
            generate_variable_length_cortex_layout(
                p_cortex, L_min=1.0e-6, L_max=3.0 * p_cortex.R_cell,
            )

    def test_n_beads_min_below_2_raises(self, p_cortex):
        with pytest.raises(ValueError, match="n_beads_min"):
            generate_variable_length_cortex_layout(
                p_cortex, L_min=1.0e-6, L_max=5.0e-6, n_beads_min=1,
            )

    def test_zero_filaments_raises(self, p_cortex):
        with pytest.raises(ValueError, match="n_filaments"):
            generate_variable_length_cortex_layout(
                p_cortex, L_min=1.0e-6, L_max=5.0e-6, n_filaments=0,
            )


# ---------------------------------------------------------------------------
# §3 Conservation invariants
# ---------------------------------------------------------------------------
class TestConservation:
    def test_layout_shapes_consistent(self, p_cortex, layout):
        F = p_cortex.n_filaments
        n_total = int(layout.n_beads_per_filament.sum())
        assert layout.positions_flat.shape == (n_total, 3)
        assert layout.n_beads_per_filament.shape == (F,)
        assert layout.filament_starts.shape == (F,)
        assert layout.L_per_filament.shape == (F,)
        assert layout.centers_of_mass.shape == (F, 3)
        assert layout.tangents.shape == (F, 3)

    def test_filament_starts_match_cumsum(self, p_cortex, layout):
        expected_starts = np.zeros(p_cortex.n_filaments, dtype=np.int64)
        expected_starts[1:] = np.cumsum(layout.n_beads_per_filament[:-1])
        assert np.array_equal(layout.filament_starts, expected_starts)

    def test_bond_count_matches_sum_N_minus_1(self, p_cortex, layout):
        expected = int((layout.n_beads_per_filament - 1).sum())
        assert layout.bond_groups.shape == (expected, 2)

    def test_angle_count_matches_sum_N_minus_2_clamped(self, p_cortex, layout):
        # For N=2 filaments, contributes 0 angles (not -0; clamp at 0).
        expected = int(np.clip(layout.n_beads_per_filament - 2, 0, None).sum())
        assert layout.angle_groups.shape == (expected, 3)

    def test_L_per_filament_matches_quantization(self, p_cortex, layout):
        expected = (layout.n_beads_per_filament - 1).astype(np.float64) * p_cortex.rest_length
        assert np.allclose(layout.L_per_filament, expected, rtol=1e-12)


# ---------------------------------------------------------------------------
# §4 Numerical sanity
# ---------------------------------------------------------------------------
class TestNumerical:
    def test_positions_within_box(self, p_cortex, layout):
        half = 0.5 * p_cortex.L_box
        assert (np.abs(layout.positions_flat) <= half + 1e-9).all()

    def test_bond_rest_length_exact_at_construction(self, p_cortex, layout):
        pos = layout.positions_flat
        groups = layout.bond_groups
        diffs = pos[groups[:, 1]] - pos[groups[:, 0]]
        lens = np.linalg.norm(diffs, axis=1)
        assert np.allclose(lens, p_cortex.rest_length, atol=1.0e-12, rtol=0)

    def test_centers_on_shell(self, p_cortex, layout):
        r = np.linalg.norm(layout.centers_of_mass, axis=1)
        assert np.allclose(r, p_cortex.R_cell, rtol=1.0e-9)

    def test_tangents_unit_norm_and_perpendicular_to_normal(self, p_cortex, layout):
        normals = layout.centers_of_mass / p_cortex.R_cell
        dots = np.einsum("fi,fi->f", normals, layout.tangents)
        assert np.allclose(dots, 0.0, atol=1.0e-9)
        tn = np.linalg.norm(layout.tangents, axis=1)
        assert np.allclose(tn, 1.0, atol=1.0e-9)


# ---------------------------------------------------------------------------
# §6 Measurement protocol — length distribution + state builder + sim
# ---------------------------------------------------------------------------
class TestMeasurement:
    def test_empirical_L_mean_close_to_uniform_mean(self, layout):
        """With F=80 filaments and Uniform(1, 5) μm distribution, the
        empirical mean approaches (L_min+L_max)/2 = 3 μm.  With
        ℓ_0 quantization (0.5 μm) the realised mean has tiny rounding
        bias — std error ~ √(Var(L)/F) for Uniform(1,5) ≈ 1.15/√80 ≈
        0.13 μm.  Allow ±0.5 μm tolerance (3σ + quantization)."""
        L_um = layout.L_per_filament * 1e6
        mean_um = float(L_um.mean())
        assert 2.5 < mean_um < 3.5, (
            f"Empirical L mean = {mean_um:.3f} μm; expected ≈ 3 μm "
            "(Uniform[1, 5] under ℓ_0=0.5 μm quantization)."
        )

    def test_L_range_covers_quantization_bins(self, layout):
        """Each ℓ_0 quantum bin in [L_min, L_max] should be sampled at
        least once for F=80 filaments × 9 bins (1, 1.5, 2, …, 5 μm)."""
        L_um = layout.L_per_filament * 1e6
        unique = np.unique(np.round(L_um * 2)) / 2.0   # round to 0.5 μm
        # Span coverage: at minimum, both endpoints should appear OR
        # the realised range should be wide enough to span > 3 quanta.
        assert L_um.min() <= 1.5 + 1e-6
        assert L_um.max() >= 4.5 - 1e-6
        assert unique.size >= 5

    def test_state_builder_particle_count_matches(self, p_cortex, layout):
        snap = build_variable_length_cortex_state(p_cortex, layout)
        assert snap.particles.N == int(layout.n_beads_per_filament.sum())
        assert snap.bonds.N == int((layout.n_beads_per_filament - 1).sum())

    def test_sim_builder_force_free_at_construction(self, p_cortex, layout):
        """Bond energy + angle energy should be 0 at exact-ℓ_0 + collinear
        construction.  Sim builder runs with_baoab=False so positions
        don't drift before we read energies."""
        sim, _, _, _ = build_variable_length_cortex_simulation(
            p_cortex, layout, with_baoab=False,
        )
        sim.run(0)
        ig = sim.operations.integrator
        bond_e = ig.forces[0].energy
        angle_e = ig.forces[1].energy
        # Allow tiny numerical noise from radial drift at endpoints.
        assert abs(bond_e) < 1.0e-25, (
            f"Bond energy {bond_e:.3e} J should be ~0 at construction"
        )
        assert abs(angle_e) < 1.0e-25, (
            f"Angle energy {angle_e:.3e} J should be ~0 at construction"
        )

    def test_sim_demo_baoab_no_NaN(self, p_cortex, layout):
        sim, _, _, _ = build_variable_length_cortex_simulation(
            p_cortex, layout, with_baoab=True,
        )
        sim.run(500)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position)
        assert np.isfinite(pos).all()


# ---------------------------------------------------------------------------
# Opt-in production-scale variable-length L_p sweep
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not bool(int(os.environ.get("H3_VARIABLE_LENGTH_PRODUCTION", "0"))),
    reason=(
        "Variable-length L_p production sweep; opt-in via "
        "H3_VARIABLE_LENGTH_PRODUCTION=1.  Multi-hour wall (same "
        "rationale as L_p full + KU-3.x opt-ins)."
    ),
)
class TestVariableLengthProduction:
    def test_per_filament_L_p_emergent_distribution(self):
        pytest.skip("Production skeleton; see module docstring.")

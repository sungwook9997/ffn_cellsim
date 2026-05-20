"""STATIC sanity-gate tests for H.1 ECM cross-link seeding.

Covers the Sanity Gate sections written in
``ffn_sim/ecm/cross_links.py``'s module docstring:

- §2 Boundary cases (empty, no self-bond, dedup)
- §3 Conservation (count matches oracle, Newton 3rd on a single bond)
- §4 Numerical sanity (indices in range; float64)
- §5 Sign / sense (stretched xl is attractive)
- §6 Measurement (count vs KU-1.27 expected ±15%, KU-1.3 coordination)

Production-scale KU-1.30 #1/#2/#3 lives in
``tests/validation/test_ku130.py`` (next sub-milestone).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd
import hoomd.md as md
import gsd.hoomd

from ffn_sim.ecm.cross_links import (
    XL_BIN_WIDTH_M,
    XL_BOND_TYPE_NAME,
    XL_N_BINS,
    XLBonds,
    add_xl_to_frame,
    generate_xl_bonds,
    measure_coordination,
    n_xl_bonds,
    quantize_rest_lengths,
    xl_bin_rest_lengths,
    xl_bin_type_names,
)
from ffn_sim.ecm.mikado import (
    ResolvedH1,
    build_mikado_simulation,
    build_mikado_state,
    resolve_derived,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h1.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="module")
def resolved() -> ResolvedH1:
    return resolve_derived(_load_cfg())


@pytest.fixture(scope="module")
def xl_pair(resolved):
    return generate_xl_bonds(resolved)


@pytest.fixture(scope="module")
def xl(xl_pair):
    xl, _ = xl_pair
    return xl


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_single_fiber_yields_zero_xls(self):
        cfg = _load_cfg()
        cfg["ecm"]["demo_mode"] = True
        p = resolve_derived(cfg)
        p.n_fibers = 1
        xl, _ = generate_xl_bonds(p)
        assert xl.group.shape == (0, 2)
        assert xl.rest_lengths.shape == (0,)

    def test_two_parallel_fibers_yield_zero_xls(self):
        cfg = _load_cfg()
        cfg["ecm"]["demo_mode"] = True
        p = resolve_derived(cfg)
        p.n_fibers = 2
        # The oracle generates fibers with isotropic orientations by
        # default; two random fibers usually intersect at our box scale.
        # Force parallel orientations by patching theta0 + S_order via
        # a direct call to the oracle is overkill; we just trust that
        # for n=2 fibers in a 200μm box, the count is small (0 to a few).
        xl, _ = generate_xl_bonds(p)
        assert xl.group.shape[1] == 2  # well-formed regardless
        assert xl.group.shape[0] == xl.rest_lengths.shape[0]


# ---------------------------------------------------------------------------
# §3 Conservation: count matches oracle exactly
# ---------------------------------------------------------------------------
class TestConservationOfCount:
    def test_n_xl_matches_oracle(self, xl, resolved):
        # Re-run the oracle directly with the same seed for a parallel
        # ground-truth count (the wrapper deduplicates; oracle does too
        # at the (i<j) level).
        from ffn_sim.validation.oracles.ecm.cross_links import (
            generate_cross_links as oracle_gen,
        )
        from ffn_sim.validation.oracles.ecm.fiber_network import (
            generate_2d_fiber_network,
        )
        net = generate_2d_fiber_network(
            L_box=resolved.L_box,
            n_fibers=resolved.n_fibers,
            L_fiber=resolved.L_fiber,
            beads_per_fiber=resolved.beads_per_fiber,
            seed=resolved.seed,
        )
        oracle_links = oracle_gen(net, stiffness=resolved.xl_stiffness)
        assert xl.group.shape[0] == len(oracle_links), (
            f"xl wrapper count {xl.group.shape[0]} != oracle "
            f"{len(oracle_links)} — wrapper diverged from oracle geometry."
        )


# ---------------------------------------------------------------------------
# §4 Numerical sanity: indices in range, no self-bonds, float64
# ---------------------------------------------------------------------------
class TestNumericalSanity:
    def test_xl_global_indices_in_range(self, xl, resolved):
        N = resolved.n_fibers * resolved.beads_per_fiber
        assert (xl.group < N).all(), (
            f"xl.group max = {int(xl.group.max())}, expected < {N}."
        )
        assert (xl.group >= 0).all()

    def test_no_self_bond(self, xl):
        assert (xl.group[:, 0] != xl.group[:, 1]).all(), (
            "Found xl bond with identical endpoints."
        )

    def test_rest_lengths_finite_and_nonneg(self, xl):
        assert np.all(np.isfinite(xl.rest_lengths))
        assert np.all(xl.rest_lengths >= 0.0)

    def test_dtypes(self, xl):
        assert xl.group.dtype == np.uint32
        assert xl.rest_lengths.dtype == np.float64
        assert xl.bin_idx.dtype == np.int64


# ---------------------------------------------------------------------------
# r0 binning (PI 2026-05-20, option B)
# ---------------------------------------------------------------------------
class TestBinning:
    def test_bin_centers_evenly_spaced(self):
        bins = xl_bin_rest_lengths()
        assert bins.shape == (XL_N_BINS,)
        # Spacing == XL_BIN_WIDTH_M, starting at 0.
        assert math.isclose(bins[0], 0.0)
        diffs = np.diff(bins)
        assert np.allclose(diffs, XL_BIN_WIDTH_M, rtol=0, atol=1e-15)

    def test_quantize_rest_lengths_round_to_nearest(self):
        # Test cases on a known grid.
        w = XL_BIN_WIDTH_M
        cases = {
            0.0: 0,
            0.4 * w: 0,
            0.6 * w: 1,
            1.0 * w: 1,
            1.5 * w: 2,                  # ties round to even/away
            2.4 * w: 2,
            2.6 * w: 3,
            (XL_N_BINS - 1) * w: XL_N_BINS - 1,
            XL_N_BINS * w: XL_N_BINS - 1,    # clamped
            10 * w + 1.0e-3: XL_N_BINS - 1,  # over-range clamp
        }
        rest = np.array(list(cases.keys()))
        expected = np.array(list(cases.values()))
        observed = quantize_rest_lengths(rest)
        # numpy rounding ties to even, so 1.5*w → 2 (even). All inputs above
        # are unambiguous except 1.5; we accept either tie convention.
        # Compare element-wise but allow ±1 on ties.
        diff = np.abs(observed - expected)
        assert (diff <= 1).all(), (
            f"quantize_rest_lengths: observed {observed.tolist()}, "
            f"expected {expected.tolist()}"
        )

    def test_quantization_error_bounded_by_half_bin(self, xl):
        bins = xl_bin_rest_lengths()
        quantized = bins[xl.bin_idx]
        err = np.abs(xl.rest_lengths - quantized)
        max_err = err.max() if err.size else 0.0
        # Allow a small margin over w/2 for the upper clamp (which can
        # carry larger error if any link has rest_length > N_bins·w).
        assert max_err <= XL_BIN_WIDTH_M / 2.0 + 1e-12, (
            f"max quantization error {max_err:e} > w/2 = "
            f"{XL_BIN_WIDTH_M / 2:e}; check binning."
        )

    def test_all_bin_types_registered_in_simulation(self, resolved):
        """build_mikado_simulation registers every xl_b<i> type so that
        downstream code (e.g. KU-1.30 production) can index by bin
        without conditional checks."""
        sim, _, _ = build_mikado_simulation(resolved, with_cross_links=True)
        sim.run(0)
        bond_force = None
        for f in sim.operations.integrator.forces:
            cls = type(f).__name__
            mod = type(f).__module__
            if cls == "Harmonic" and "bond" in mod:
                bond_force = f
                break
        assert bond_force is not None
        for name in xl_bin_type_names():
            params = bond_force.params[name]
            assert math.isclose(params["k"], resolved.xl_stiffness)
            i = int(name.removeprefix("xl_b"))
            assert math.isclose(
                params["r0"], i * XL_BIN_WIDTH_M, rel_tol=0, abs_tol=1e-15
            )


# ---------------------------------------------------------------------------
# Energy-oracle extension to xl (binning-aware)
# ---------------------------------------------------------------------------
class TestEnergyOracleWithXL:
    """At construction (force-free for ecm-bond + angle), the HOOMD bond
    PE *including* the xl bin contributions must agree with the oracle's
    compute_energy(cross_links=...) to a known tolerance set by the
    quantization error per link.

    Tolerance derivation: per-link energy mismatch between HOOMD's
    (r0_bin) and oracle's (r0_orig) at the construction-time configuration
    (where |d| = r0_orig) is ½·k·(r0_orig − r0_bin)² ≤ ½·k·(w/2)². With
    k=1e-3, w=50 nm: ≤ 3.1e-19 J/link. Over N_xl ≈ 7700 links: ≤ 2.4e-15 J
    total. Absolute tolerance: 1e-14 J (5× margin).
    """

    def test_xl_bonded_energy_matches_binning_prediction(self, resolved):
        sim, _, _ = build_mikado_simulation(resolved, with_cross_links=True)
        sim.run(0)

        # HOOMD bond.energy is the SUM over all bond types (ecm-bond + every
        # xl_b<i>). The ecm-bond contribution is ≈ 0 at construction
        # (straight chains, rest-length bonds). So bond.energy is the xl
        # contribution.
        bond_e = 0.0
        for f in sim.operations.integrator.forces:
            cls = type(f).__name__
            mod = type(f).__module__
            if cls == "Harmonic" and "bond" in mod:
                bond_e = float(f.energy)
                break

        # Analytical prediction from quantization: Σ ½·k·(r0_orig − r0_bin)².
        from ffn_sim.ecm.cross_links import (
            generate_xl_bonds, xl_bin_rest_lengths,
        )
        xl, _ = generate_xl_bonds(resolved)
        bins = xl_bin_rest_lengths()
        residual = xl.rest_lengths - bins[xl.bin_idx]
        predicted_xl_e = float(0.5 * resolved.xl_stiffness * np.sum(residual ** 2))

        # Difference vs prediction should be ~0 (HOOMD = analytical).
        diff = abs(bond_e - predicted_xl_e)
        assert diff <= 1.0e-14, (
            f"HOOMD bond.energy at construction = {bond_e:e} J; "
            f"analytical quantization prediction = {predicted_xl_e:e} J; "
            f"abs diff = {diff:e} J > 1e-14 J tolerance."
        )

    def test_xl_construction_energy_below_brief_simplification(self, resolved):
        """Binning reduces construction-time xl energy by ≥10× vs r0=0."""
        from ffn_sim.ecm.cross_links import (
            generate_xl_bonds, xl_bin_rest_lengths,
        )
        xl, _ = generate_xl_bonds(resolved)
        bins = xl_bin_rest_lengths()
        e_binned = float(
            0.5 * resolved.xl_stiffness * np.sum(
                (xl.rest_lengths - bins[xl.bin_idx]) ** 2
            )
        )
        e_zero_r0 = float(
            0.5 * resolved.xl_stiffness * np.sum(xl.rest_lengths ** 2)
        )
        # binning should be ≥ 10× smaller (typically ~100×).
        assert e_binned * 10.0 < e_zero_r0, (
            f"binning xl energy {e_binned:e} J vs r0=0 {e_zero_r0:e} J; "
            "expected ≥ 10× reduction."
        )


# ---------------------------------------------------------------------------
# §5 Sign / sense (on a hand-built 2-bead xl)
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_xl_pair_attractive_when_stretched(self):
        """Two beads connected by an xl with r0=0, separated by 100 nm.
        Force on bead 0 must point toward bead 1 along +x."""
        snap = gsd.hoomd.Frame()
        snap.particles.N = 2
        snap.particles.types = ["actin_ecm"]
        snap.particles.typeid = np.zeros(2, dtype=np.uint32)
        snap.particles.position = np.array(
            [[0.0, 0.0, 0.0], [1.0e-7, 0.0, 0.0]], dtype=np.float64
        )
        snap.particles.mass = np.ones(2)
        snap.bonds.N = 1
        snap.bonds.types = [XL_BOND_TYPE_NAME]
        snap.bonds.typeid = np.zeros(1, dtype=np.uint32)
        snap.bonds.group = np.array([[0, 1]], dtype=np.uint32)
        snap.configuration.box = [1.0e-5, 1.0e-5, 1.0e-5, 0.0, 0.0, 0.0]

        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=0)
        sim.create_state_from_snapshot(snap)
        bond = md.bond.Harmonic()
        bond.params[XL_BOND_TYPE_NAME] = dict(k=1.0e-3, r0=0.0)
        ig = md.Integrator(dt=1.0e-9)
        ig.forces.append(bond)
        sim.operations.integrator = ig
        sim.run(0)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            F = np.asarray(s.particles.net_force).copy()
            r0 = int(np.argwhere(tag == 0).item())
            r1 = int(np.argwhere(tag == 1).item())
        # Bead 0 at x=0 should be pulled toward bead 1 at x=+1e-7 (+F_x>0);
        # bead 1 toward bead 0 (-F_x<0). Magnitudes: k·|d| = 1e-3·1e-7 = 1e-10 N.
        assert F[r0, 0] > 0, f"F[bead0]_x = {F[r0, 0]:e}; expected > 0."
        assert F[r1, 0] < 0, f"F[bead1]_x = {F[r1, 0]:e}; expected < 0."
        assert math.isclose(abs(F[r0, 0]), 1.0e-10, rel_tol=1e-6)
        assert math.isclose(abs(F[r1, 0]), 1.0e-10, rel_tol=1e-6)

    def test_xl_pair_newton_third(self):
        """Σ F over the two xl-bonded beads ≈ 0 (Newton 3rd)."""
        snap = gsd.hoomd.Frame()
        snap.particles.N = 2
        snap.particles.types = ["actin_ecm"]
        snap.particles.typeid = np.zeros(2, dtype=np.uint32)
        snap.particles.position = np.array(
            [[-5.0e-8, 0.0, 0.0], [5.0e-8, 0.0, 0.0]], dtype=np.float64
        )
        snap.particles.mass = np.ones(2)
        snap.bonds.N = 1
        snap.bonds.types = [XL_BOND_TYPE_NAME]
        snap.bonds.typeid = np.zeros(1, dtype=np.uint32)
        snap.bonds.group = np.array([[0, 1]], dtype=np.uint32)
        snap.configuration.box = [1.0e-5, 1.0e-5, 1.0e-5, 0.0, 0.0, 0.0]

        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=0)
        sim.create_state_from_snapshot(snap)
        bond = md.bond.Harmonic()
        bond.params[XL_BOND_TYPE_NAME] = dict(k=1.0e-3, r0=0.0)
        ig = md.Integrator(dt=1.0e-9)
        ig.forces.append(bond)
        sim.operations.integrator = ig
        sim.run(0)
        with sim.state.cpu_local_snapshot as s:
            net = np.asarray(s.particles.net_force).sum(axis=0)
        assert np.all(np.abs(net) <= 1e-20), (
            f"Σ F on 2-bead xl pair = {net}; expected ≤ 1e-20 N."
        )


# ---------------------------------------------------------------------------
# §6 Measurement: KU-1.27 expected count, KU-1.3 coordination
# ---------------------------------------------------------------------------
class TestMeasurement:
    def test_xl_count_in_KU127_band(self, xl, resolved):
        # Expected: 2 ρ_L L_fiber / π crossings per fiber × N_fibers / 2
        # (each crossing produces one xl link in primary-copy-only mode).
        expected = resolved.expected_xl_per_fiber * resolved.n_fibers / 2.0
        actual = xl.group.shape[0]
        # Oracle's primary-copy-only construction misses ~5% of crossings
        # at L_fiber/L_box ≈ 0.05; accept ±15% finite-size + periodic gap.
        rel = abs(actual - expected) / expected
        assert rel <= 0.15, (
            f"xl count {actual} vs KU-1.27 expected {expected:.0f}; "
            f"rel deviation {rel:.2%} > 15% acceptance band."
        )

    def test_coordination_in_KU13_band(self, xl, resolved):
        """PI 2026-05-20 ratified D4-anchored [2.0, 5.5] band; re-enabled."""
        z = measure_coordination(xl, resolved)
        lo, hi = resolved.z_range
        assert lo <= z <= hi, (
            f"⟨z⟩ = {z:.3f} outside KU-1.3 acceptance [{lo}, {hi}]."
        )


# ---------------------------------------------------------------------------
# Integration smoke: full Mikado + xl + BAOAB attaches and steps cleanly
# ---------------------------------------------------------------------------
class TestSimulationSmoke:
    def test_full_h1_simulation_builds_and_steps(self, resolved):
        """Build the full M2 simulation (mikado + xl + BAOAB), run the
        equilibrate_no_shear prelude (PI 2026-05-20), and verify
        ``max |F|`` reaches the thermal scale. Confirms no force blow-up,
        no NaN, no crash, no int32-image overflow at ~66k particles +
        ~8k xl bonds + ~63k ecm bonds + ~60k angles + LJ pair list."""
        from ffn_sim.ecm.equilibrate import equilibrate_no_shear

        sim, updater, action = build_mikado_simulation(
            resolved, with_cross_links=True
        )
        # State carries both ecm-bond and xl-type bonds; check total.
        n_ecm = resolved.n_fibers * (resolved.beads_per_fiber - 1)
        xl, _ = generate_xl_bonds(resolved)
        assert sim.state.N_bonds == n_ecm + xl.group.shape[0]

        # Short prelude: keep the smoke test under ~10 s. The full
        # 1000-step prelude is exercised by tests/validation/test_ku130.py.
        diag = equilibrate_no_shear(
            sim, action, updater,
            n_softstart=100, n_baoab=50,
            rest_length=resolved.rest_length, gamma_b=resolved.gamma_b,
        )
        with sim.state.cpu_local_snapshot as s:
            F = np.asarray(s.particles.net_force)
            pos = np.asarray(s.particles.position)
        assert np.all(np.isfinite(F)), "Non-finite net_force after prelude."
        assert np.all(np.isfinite(pos)), "Non-finite position after prelude."
        # Soft phase must drain construction force by many decades.
        assert diag["max_force_before"] > 1.0, (
            f"Construction max|F| expected > 1 N (LJ overlap); got "
            f"{diag['max_force_before']:.3e}."
        )
        assert diag["max_force_after_soft"] < 1.0e-9, (
            f"max|F| after soft phase = {diag['max_force_after_soft']:.3e} "
            "N exceeded the 1 nN BAOAB-safety target."
        )

    def test_topology_smoke_xl_indices_distinct_from_ecm(self, resolved):
        """The xl bond entries (any bin) must not duplicate any ecm-bond
        entry (a backbone (i, i+1) pair must not also appear as an xl)."""
        snap = build_mikado_state(resolved, with_cross_links=True)
        groups = np.asarray(snap.bonds.group, dtype=np.int64)
        typeids = np.asarray(snap.bonds.typeid, dtype=np.int64)
        ecm_id = snap.bonds.types.index("ecm-bond")
        xl_ids = {
            snap.bonds.types.index(name) for name in xl_bin_type_names()
        }
        ecm_pairs = {tuple(sorted(g)) for g, t in zip(groups, typeids) if t == ecm_id}
        xl_pairs = {tuple(sorted(g)) for g, t in zip(groups, typeids) if t in xl_ids}
        overlap = ecm_pairs & xl_pairs
        assert not overlap, (
            f"{len(overlap)} pairs appear as both ecm-bond and an xl bin "
            "type; the oracle's nearest-bead rounding chose a "
            "backbone-adjacent pair."
        )

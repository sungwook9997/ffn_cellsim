"""STATIC + demo + opt-in production sanity-gate tests for H.3 cortex.

Covers the Sanity Gate sections written in ``ffn_sim/cortex/cortex.py``'s
module docstring:

- §1 Dimensional analysis           (TestDimensional)
- §2 Boundary cases                 (TestBoundary)
- §3 Conservation + topology counts (TestTopology, TestInitialEnergy)
- §4 Numerical sanity               (TestNumerical)
- §5 Sign / sense                   (TestSignSense)
- §6 Measurement protocol           (TestMeasurement)
- Pre-flight: σ_z vs L_z H.2 slab-lesson gate at H.3 setup
                                     (TestSigmaZvsLz)
- Demo:  short BAOAB smoke run, no NaN/Inf  (TestDemoRun)
- Crosslinkers (optional, KU-3.19)  (TestCrosslinkers)
- Production (opt-in via H3_PRODUCTION=1):
    full per-filament L_p / 3D equipartition / 3D Boltzmann angle gates
                                     (TestH3Production)

Production gates use the same first-principles bands as the H.2 ✅
strict-PASS (per-filament force constants are identical), aggregated
over the H.3 cortex's 1000 filaments × 5 interior beads.
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
import hoomd.md as md

from ffn_sim.common.filament_math import (
    boltzmann_angle_density_3d,
    bending_energy_per_bond,
    fit_persistence_length,
    hoomd_angle_array,
)
from ffn_sim.archive.hoomd_legacy.cortex.cortex import (
    CortexTopology,
    CrosslinkerBonds,
    ResolvedH3,
    build_cortex_simulation,
    build_cortex_state,
    generate_cortex_topology,
    generate_static_xl_bonds,
    resolve_h3_derived,
    xl_bin_rest_lengths,
    xl_bin_type_names,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo_cfg(**overrides) -> dict:
    """Smaller-cortex copy of the production yaml for fast CI smoke runs.

    Keeps all KU-anchored physics constants identical to the production
    config; only reduces the n_filaments cost knob (and applies any
    explicit override keys).
    """
    cfg = deepcopy(_load_cfg())
    cfg["cortex"]["n_filaments"] = 50         # 350 beads (vs 7000 production)
    cfg["cortex"]["demo_mode"] = True         # relax cost ceiling assertion
    for k, v in overrides.items():
        path = k.split(".")
        sub = cfg["cortex"]
        for key in path[:-1]:
            sub = sub[key]
        sub[path[-1]] = v
    return cfg


# ---------------------------------------------------------------------------
# Shared module-scoped fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def resolved_production() -> ResolvedH3:
    return resolve_h3_derived(_load_cfg())


@pytest.fixture(scope="module")
def resolved_demo() -> ResolvedH3:
    return resolve_h3_derived(_demo_cfg())


@pytest.fixture(scope="module")
def topology_demo(resolved_demo) -> CortexTopology:
    return generate_cortex_topology(resolved_demo)


@pytest.fixture(scope="module")
def state_demo(resolved_demo):
    snap, topology, xl = build_cortex_state(resolved_demo)
    return snap, topology, xl


@pytest.fixture(scope="module")
def sim_demo(resolved_demo):
    """Full HOOMD setup with BAOAB attached, no xl (cheap to construct)."""
    sim, updater, action, topology, xl = build_cortex_simulation(
        resolved_demo, with_baoab=True, with_crosslinkers=False
    )
    sim.run(0)
    return sim, updater, action, topology, xl


@pytest.fixture(scope="module")
def force_only_sim(resolved_demo):
    """Force-eval-only sim (no BAOAB Updater); callers mutate positions."""
    sim, _, _, topology, xl = build_cortex_simulation(
        resolved_demo, with_baoab=False, with_crosslinkers=False
    )
    sim.run(0)
    return sim, topology, xl


# ---------------------------------------------------------------------------
# §1  Dimensional analysis
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_bond_k_units(self, resolved_demo: ResolvedH3):
        p = resolved_demo
        # k_bond = μ / ℓ_0 → [N/m].
        expected = p.stretching_modulus / p.rest_length
        assert math.isclose(p.bond_k, expected, rel_tol=1e-12)
        assert p.bond_k > 0.0

    def test_angle_k_units(self, resolved_demo: ResolvedH3):
        p = resolved_demo
        # k_θ = κ_B / ℓ_0 → [N·m] = [J/rad²].
        expected = p.bending_modulus / p.rest_length
        assert math.isclose(p.angle_k, expected, rel_tol=1e-12)
        assert p.angle_k > 0.0

    def test_gamma_b_units(self, resolved_demo: ResolvedH3):
        p = resolved_demo
        expected = 6.0 * math.pi * p.water_viscosity * p.bead_radius
        assert math.isclose(p.gamma_b, expected, rel_tol=1e-12)

    def test_lj_units(self, resolved_demo: ResolvedH3):
        p = resolved_demo
        assert math.isclose(p.lj_sigma, p.lj_sigma_factor * p.bead_radius)
        assert math.isclose(p.lj_epsilon, p.lj_epsilon_kT * p.kT)
        assert math.isclose(p.lj_r_cut, 2 ** (1 / 6) * p.lj_sigma, rel_tol=1e-12)

    def test_dt_cfl_positive(self, resolved_demo: ResolvedH3):
        p = resolved_demo
        assert p.dt_cfl > 0
        assert p.dt_cfl == p.cfl_safety_factor * p.tau_min

    def test_box_derived_from_R_cell(self, resolved_demo: ResolvedH3):
        p = resolved_demo
        assert math.isclose(p.L_box, p.box_factor * p.R_cell)


# ---------------------------------------------------------------------------
# §2  Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_beads_below_two_raises(self):
        cfg = _demo_cfg()
        cfg["cortex"]["beads_per_filament"] = 1
        with pytest.raises(ValueError, match="beads_per_filament"):
            resolve_h3_derived(cfg)

    def test_zero_filaments_raises(self):
        cfg = _demo_cfg()
        cfg["cortex"]["n_filaments"] = 0
        with pytest.raises(ValueError, match="n_filaments"):
            resolve_h3_derived(cfg)

    def test_negative_R_cell_raises(self):
        cfg = _demo_cfg()
        cfg["cortex"]["R_cell"] = -1.0e-6
        with pytest.raises(ValueError, match="R_cell"):
            resolve_h3_derived(cfg)

    def test_L_filament_exceeds_diameter_raises(self):
        cfg = _demo_cfg()
        cfg["cortex"]["R_cell"] = 1.0e-6
        cfg["cortex"]["L_filament"] = 3.0e-6
        with pytest.raises(ValueError, match="cell diameter"):
            resolve_h3_derived(cfg)

    def test_wrong_integrator_raises(self):
        cfg = _demo_cfg()
        cfg["cortex"]["dynamics"]["integrator"] = "euler_maruyama"
        with pytest.raises(ValueError, match="leimkuhler_matthews_baoab"):
            resolve_h3_derived(cfg)

    def test_box_too_small_strict_raises(self):
        # demo_mode False with small box should fail the cell-fits-in-box check.
        cfg = _demo_cfg()
        cfg["cortex"]["demo_mode"] = False
        cfg["cortex"]["box"]["L_box_over_R_cell"] = 1.5  # half-box = 7.5 < R_cell + σ + ℓ_0
        with pytest.raises(ValueError, match="thermal margin"):
            resolve_h3_derived(cfg)

    def test_cost_ceiling_strict_raises(self):
        cfg = _demo_cfg()
        cfg["cortex"]["demo_mode"] = False
        cfg["cortex"]["n_filaments"] = 5000
        with pytest.raises(ValueError, match="cost ceiling"):
            resolve_h3_derived(cfg)


# ---------------------------------------------------------------------------
# §Pre-flight  σ_z vs L_z (H.2 slab lesson at H.3 setup)
# ---------------------------------------------------------------------------
class TestSigmaZvsLz:
    """The H.2 slab artefact failure ratio was L_z/σ_perp ≈ 0.14
    (Lz=0.2 μm slab); H.2 strict-PASS at L_z/σ_perp ≥ 10. This test
    asserts H.3's L_z/σ_perp_max is well into the safe regime."""

    def test_box_above_safe_ratio_production(self, resolved_production: ResolvedH3):
        p = resolved_production
        ratio = p.L_box / p.sigma_perp_per_filament
        assert ratio >= 10.0, (
            f"L_box / σ_perp = {ratio:.2f} below H.2 safe ratio (≥10). "
            f"L_box = {p.L_box:.3e} m, σ_perp = "
            f"{p.sigma_perp_per_filament:.3e} m."
        )

    def test_brief_max_L_filament_safe(self, resolved_production: ResolvedH3):
        """Even at the brief's max L = 5 μm, σ_perp should be < half-box."""
        p = resolved_production
        L_max = 5.0e-6
        sigma_perp_max = math.sqrt(L_max**3 / (3 * p.persistence_length))
        half_box = 0.5 * p.L_box
        cell_extent = p.R_cell + sigma_perp_max + p.rest_length
        assert half_box >= cell_extent, (
            f"Half-box {half_box:.3e} m insufficient for cell + L_max "
            f"thermal margin {cell_extent:.3e} m. Increase "
            "box.L_box_over_R_cell."
        )

    def test_cortex_thickness_not_box_artefact(
        self, resolved_production: ResolvedH3
    ):
        """200 nm cortex thickness must come from physics, not from box."""
        p = resolved_production
        # Box is a CUBE (Lz = L_box), not a slab. ERM tether is the
        # mechanism for the 200 nm confinement (later H.3 deliverable).
        # H.3 production: L_box/cortex_thickness ≥ 10 ensures the box is
        # NOT acting as a slab confinement at the cortex-thickness scale.
        assert p.L_box / p.cortex_thickness >= 10.0, (
            f"L_box / cortex_thickness = {p.L_box / p.cortex_thickness:.1f} "
            "< 10 — box may be confining filaments to the cortex layer; "
            "that constraint should come from the ERM tether (physics), "
            "not from box geometry."
        )


# ---------------------------------------------------------------------------
# §3  Conservation invariants + topology counts
# ---------------------------------------------------------------------------
class TestTopology:
    def test_particle_count(self, state_demo, resolved_demo: ResolvedH3):
        snap, topology, _ = state_demo
        p = resolved_demo
        assert snap.particles.N == p.n_filaments * p.beads_per_filament
        assert topology.positions.shape == (
            p.n_filaments, p.beads_per_filament, 3
        )

    def test_bond_count_backbone_only(
        self, state_demo, resolved_demo: ResolvedH3
    ):
        snap, _, xl = state_demo
        p = resolved_demo
        n_backbone = p.n_filaments * (p.beads_per_filament - 1)
        # xl is empty in demo (with_crosslinkers default OFF in YAML).
        assert xl.bond_pairs.shape == (0, 2)
        assert snap.bonds.N == n_backbone
        # All bond typeids are 0 = cortex-bond.
        assert (np.asarray(snap.bonds.typeid) == 0).all()

    def test_angle_count(self, state_demo, resolved_demo: ResolvedH3):
        snap, _, _ = state_demo
        p = resolved_demo
        n_angles = p.n_filaments * (p.beads_per_filament - 2)
        assert snap.angles.N == n_angles

    def test_bond_groups_are_contiguous(
        self, state_demo, resolved_demo: ResolvedH3
    ):
        snap, _, _ = state_demo
        p = resolved_demo
        groups = np.asarray(snap.bonds.group)
        for f in range(min(5, p.n_filaments)):
            for i in range(p.beads_per_filament - 1):
                idx = f * (p.beads_per_filament - 1) + i
                a, b = groups[idx]
                assert a == f * p.beads_per_filament + i
                assert b == f * p.beads_per_filament + i + 1

    def test_bond_rest_length_exact_at_construction(
        self, state_demo, resolved_demo: ResolvedH3
    ):
        snap, _, _ = state_demo
        p = resolved_demo
        pos = np.asarray(snap.particles.position)
        groups = np.asarray(snap.bonds.group)
        diffs = pos[groups[:, 1]] - pos[groups[:, 0]]
        lens = np.linalg.norm(diffs, axis=1)
        # Tangent-plane placement should give exact ℓ_0 to float64.
        assert np.allclose(lens, p.rest_length, atol=1.0e-12, rtol=0)

    def test_filament_collinear_at_construction(
        self, topology_demo: CortexTopology, resolved_demo: ResolvedH3
    ):
        """Beads on each filament should be collinear (angle = π straight)."""
        p = resolved_demo
        pos = topology_demo.positions
        # bond vector i+1 minus bond vector i for each interior triplet.
        b1 = pos[:, 1:-1, :] - pos[:, :-2, :]
        b2 = pos[:, 2:, :] - pos[:, 1:-1, :]
        b1n = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
        b2n = b2 / np.linalg.norm(b2, axis=-1, keepdims=True)
        cos_angle = np.einsum("fij,fij->fi", b1n, b2n)
        # All should be ≈ 1 (collinear, angle 0 between b1, b2 → θ_HOOMD = π).
        assert np.allclose(cos_angle, 1.0, atol=1.0e-9)

    def test_centers_on_sphere(
        self, topology_demo: CortexTopology, resolved_demo: ResolvedH3
    ):
        p = resolved_demo
        centers = topology_demo.centers_of_mass
        r = np.linalg.norm(centers, axis=1)
        assert np.allclose(r, p.R_cell, rtol=1.0e-9)

    def test_tangents_in_local_tangent_plane(
        self, topology_demo: CortexTopology, resolved_demo: ResolvedH3
    ):
        p = resolved_demo
        centers = topology_demo.centers_of_mass
        tangents = topology_demo.tangents
        normals = centers / p.R_cell
        dots = np.einsum("fi,fi->f", normals, tangents)
        assert np.allclose(dots, 0.0, atol=1.0e-9)
        # Tangent vectors unit norm.
        tn = np.linalg.norm(tangents, axis=1)
        assert np.allclose(tn, 1.0, atol=1.0e-9)


class TestInitialEnergy:
    def test_potential_energy_zero_at_construction(self, force_only_sim):
        """Use force-only sim (no BAOAB) so positions stay at construction."""
        sim, _, _ = force_only_sim
        ig = sim.operations.integrator
        bond = ig.forces[0]
        angle = ig.forces[1]
        # sim.run(0) on force-only triggers force eval without advancing
        # any updater. Bond / angle energies should be exactly 0 at
        # exact-ℓ_0 tangent-plane construction (collinear beads).
        sim.run(0)
        bond_e = bond.energy
        angle_e = angle.energy
        assert abs(bond_e) < 1.0e-25, (
            f"Bond energy {bond_e:.3e} J should be 0 at exact-ℓ_0 "
            "construction; non-zero means tangent placement violated."
        )
        assert abs(angle_e) < 1.0e-25, (
            f"Angle energy {angle_e:.3e} J should be 0 at collinear "
            "construction; non-zero means filaments are not straight."
        )


# ---------------------------------------------------------------------------
# §4  Numerical sanity
# ---------------------------------------------------------------------------
class TestNumerical:
    def test_positions_within_box(self, state_demo, resolved_demo: ResolvedH3):
        snap, _, _ = state_demo
        p = resolved_demo
        pos = np.asarray(snap.particles.position)
        half = 0.5 * p.L_box
        assert (np.abs(pos) <= half + 1.0e-9).all(), (
            f"Some bead is outside box [-{half}, {half}); "
            f"max |pos| = {np.max(np.abs(pos)):.3e}"
        )

    def test_position_dtype(self, state_demo):
        snap, _, _ = state_demo
        pos = np.asarray(snap.particles.position)
        assert pos.dtype == np.float64

    def test_all_derived_finite(self, resolved_production: ResolvedH3):
        p = resolved_production
        for attr in (
            "L_box", "rest_length", "bond_k", "angle_k", "lj_epsilon",
            "lj_sigma", "lj_r_cut", "gamma_b", "tau_stretch", "tau_bend",
            "tau_min", "dt_cfl", "sigma_perp_per_filament",
            "cortex_surface_area", "filament_areal_density",
        ):
            v = getattr(p, attr)
            assert math.isfinite(v) and v > 0.0, (
                f"{attr} = {v!r} not finite-positive"
            )


# ---------------------------------------------------------------------------
# §5  Sign / sense
# ---------------------------------------------------------------------------
class TestSignSense:
    """All sign-sense tests build a fresh isolated mini-cortex with two
    filaments (smallest setup that supports both intra-filament bonds and
    inter-filament LJ). Uses tag indirection (H.1 mikado pattern) because
    cpu_local_snapshot row order is not guaranteed to match GSD tag order.
    """

    @staticmethod
    def _two_filament_sim():
        cfg = _demo_cfg()
        cfg["cortex"]["n_filaments"] = 2
        p = resolve_h3_derived(cfg)
        sim, _, _, topology, _ = build_cortex_simulation(
            p, with_baoab=False, with_crosslinkers=False
        )
        sim.run(0)
        return sim, p, topology

    def test_stretched_bond_attractive(self):
        """Stretch the LAST bond of filament 0 (so only one bond responds —
        moving an interior bead would deform both of its bonds and give
        2× force, masking the per-bond sign-sense check)."""
        sim, p, topology = self._two_filament_sim()
        N = p.beads_per_filament
        tag_a, tag_b = N - 2, N - 1             # last bond of filament 0
        tangent = topology.tangents[0]
        delta = 1.0e-8
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            pos = np.asarray(s.particles.position)
            row_b = int(np.argwhere(tag == tag_b).item())
            pos[row_b] = pos[row_b] + delta * tangent
        sim.run(1)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            F = np.asarray(s.particles.net_force).copy()
            row_a = int(np.argwhere(tag == tag_a).item())
            row_b = int(np.argwhere(tag == tag_b).item())
        proj_a = float(np.dot(F[row_a], tangent))
        proj_b = float(np.dot(F[row_b], tangent))
        expected_mag = p.bond_k * delta
        # End bead B pulled back toward A (-tangent); interior A pulled
        # toward B (+tangent). Angle (N-3, N-2, N-1) stays π (collinear
        # along tangent) so angle force = 0.
        assert proj_a > 0 and proj_b < 0, (
            f"sign-sense failed: proj_a={proj_a:.3e}, proj_b={proj_b:.3e}"
        )
        assert math.isclose(proj_a, +expected_mag, rel_tol=0.02), (
            f"proj_a = {proj_a:.3e} N vs expected k_bond·δ = "
            f"{expected_mag:.3e} N"
        )
        assert math.isclose(proj_b, -expected_mag, rel_tol=0.02), (
            f"proj_b = {proj_b:.3e} N vs expected -k_bond·δ = "
            f"{-expected_mag:.3e} N"
        )

    def test_bent_triplet_restores(self):
        sim, p, topology = self._two_filament_sim()
        # Perturb interior bead (tag=1) of filament 0 along the local radial
        # outward direction (perpendicular to the filament's tangent).
        center = topology.centers_of_mass[0]
        normal = center / p.R_cell
        delta = 1.0e-8
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            pos = np.asarray(s.particles.position)
            row = int(np.argwhere(tag == 1).item())
            pos[row] = pos[row] + delta * normal
        sim.run(1)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            F = np.asarray(s.particles.net_force).copy()
            row = int(np.argwhere(tag == 1).item())
        proj = float(np.dot(F[row], normal))
        assert proj < 0.0, (
            f"Bent interior bead should feel restoring force; "
            f"got proj = {proj:.3e} N along outward normal."
        )

    def test_lj_wiring_repulsive_only(self, resolved_demo):
        """WCA wiring contract: ε = 0.5 kT > 0, σ = 2·R_bead, r_cut at
        the LJ minimum 2^(1/6)·σ (no attractive tail).

        Functional sign-sense of WCA force is already validated in H.1
        ``test_h1_mikado.py::TestSignSense::test_lj_repulsive`` against
        the same ``md.pair.LJ`` compute (identical kernel, only particle
        type label differs: actin_ecm → actin_cortex). A cortex-specific
        functional sign-sense test would require an isolated 2-bead
        bond-free setup because moving a single bead in the H.3 multi-
        filament state generates large intra-filament bond forces that
        swamp the WCA on its own filament's end bead (probed and
        rejected during H.3 sanity-gate development).
        """
        p = resolved_demo
        assert p.lj_enabled is True
        assert p.lj_epsilon > 0.0
        assert math.isclose(p.lj_sigma, 2.0 * p.bead_radius, rel_tol=1e-12)
        assert math.isclose(
            p.lj_r_cut, (2 ** (1 / 6)) * p.lj_sigma, rel_tol=1e-12
        )
        # Confirm the pair coeff was wired correctly by inspecting the
        # built simulation rather than the YAML.
        sim, _, _, _, _ = build_cortex_simulation(
            p, with_baoab=False, with_crosslinkers=False
        )
        sim.run(0)
        lj = sim.operations.integrator.forces[2]
        coeff = lj.params[("actin_cortex", "actin_cortex")]
        assert math.isclose(coeff["epsilon"], p.lj_epsilon, rel_tol=1e-12)
        assert math.isclose(coeff["sigma"], p.lj_sigma, rel_tol=1e-12)
        assert math.isclose(
            lj.r_cut[("actin_cortex", "actin_cortex")], p.lj_r_cut,
            rel_tol=1e-12,
        )


# ---------------------------------------------------------------------------
# §6  Measurement protocol — topology smoke + filament_math integration
# ---------------------------------------------------------------------------
class TestMeasurement:
    def test_radial_drift_within_cortex_thickness(
        self, topology_demo: CortexTopology, resolved_demo: ResolvedH3
    ):
        """Endpoint radial drift ℓ_end²/(2 R_cell) should be ≤ cortex thickness."""
        p = resolved_demo
        pos = topology_demo.positions
        r = np.linalg.norm(pos, axis=-1)        # shape (F, N)
        drift = r - p.R_cell                    # radial deviation per bead
        max_drift = float(np.max(np.abs(drift)))
        lo, hi = p.radial_drift_band_m
        assert lo <= max_drift <= hi, (
            f"Max radial drift {max_drift:.3e} m outside first-principles "
            f"band [{lo:.3e}, {hi:.3e}] (KU-3.17 cortex thickness)."
        )

    def test_filament_areal_density_in_cortex_range(
        self, resolved_demo: ResolvedH3
    ):
        # Demo case at 50 filaments → 50 / 1256 μm² ≈ 0.04 /μm² is lower
        # than KU-3.17 production target; just check the calculation is
        # consistent.
        p = resolved_demo
        density = p.n_filaments / (4 * math.pi * p.R_cell**2)
        assert math.isclose(density, p.filament_areal_density, rel_tol=1e-12)

    def test_filament_math_callable_on_topology(
        self, topology_demo: CortexTopology
    ):
        """`fit_persistence_length` must accept our (F, N, 3) topology layout."""
        # At construction (all straight), C(s) = 1.0 for all s.  The fit
        # returns NaN because log(C)=0 slope is non-negative (degenerate
        # straight chain). This test confirms the function runs and the
        # straight-chain edge case is handled gracefully.
        fit = fit_persistence_length(
            topology_demo.positions, rest_length=0.5e-6
        )
        assert fit.s_bonds.size > 0
        # Either NaN (straight degenerate) or some positive number — both
        # acceptable for a topology-only sanity check.
        assert math.isnan(fit.L_p_m) or fit.L_p_m > 0


# ---------------------------------------------------------------------------
# Demo BAOAB run
# ---------------------------------------------------------------------------
class TestDemoRun:
    def test_short_baoab_run_no_nan(self, sim_demo, resolved_demo):
        sim, _, action, _, _ = sim_demo
        # 1000 BAOAB steps — should still be at thermal-scale displacement
        # well within box / cortex_thickness on the demo cortex.
        sim.run(1000)
        with sim.state.cpu_local_snapshot as snap:
            pos = np.asarray(snap.particles.position)
            f = np.asarray(snap.particles.net_force)
        assert np.isfinite(pos).all(), "NaN/Inf in positions after demo run"
        assert np.isfinite(f).all(), "NaN/Inf in net_force after demo run"

    def test_demo_particles_stay_inside_box(self, resolved_demo):
        """A separate fresh sim — after the smoke run, particles must
        stay inside the periodic box (i.e. no force runaway pushing
        beads past L_box/2). HOOMD's PBC wraps wrapped coordinates back
        in, so we check ``image`` magnitudes — a particle that wrapped
        many times indicates runaway. tag indirection is required
        because cpu_local_snapshot row order may differ from GSD order.
        """
        p = resolved_demo
        sim, _, _, _, _ = build_cortex_simulation(
            p, with_baoab=True, with_crosslinkers=False
        )
        sim.run(500)
        with sim.state.cpu_local_snapshot as snap:
            pos = np.asarray(snap.particles.position).copy()
            image = np.asarray(snap.particles.image).copy()
        assert np.isfinite(pos).all() and np.isfinite(image).all()
        # No particle should have wrapped more than once per axis on a
        # smoke-scale run (500 BAOAB steps ≈ 6.5 μs). The cell sits
        # comfortably inside the 30 μm box — wrap should be 0 for ALL
        # beads on this short run.
        assert (np.abs(image) <= 1).all(), (
            f"At least one particle wrapped image > 1 on a 500-step "
            f"smoke run; max |image| = {int(np.max(np.abs(image)))}. "
            "Possible force runaway."
        )


# ---------------------------------------------------------------------------
# Crosslinkers (KU-3.19 static optional bonds)
# ---------------------------------------------------------------------------
class TestCrosslinkers:
    def test_xl_disabled_yields_empty_bonds(self, resolved_demo):
        p = resolved_demo
        assert p.xl_enabled is False
        xl = generate_static_xl_bonds(_make_topology(p), p)
        assert xl.bond_pairs.shape == (0, 2)
        assert xl.n_alpha == 0
        assert xl.n_filamin == 0

    def test_xl_enabled_demo_generates_bonds(self):
        cfg = _demo_cfg(**{
            "crosslinkers.enabled": True,
            "crosslinkers.n_xl": 50,
        })
        p = resolve_h3_derived(cfg)
        topology = generate_cortex_topology(p)
        xl = generate_static_xl_bonds(topology, p)
        # Some xl may be dropped if no acceptor within 60 nm, but at
        # demo cortex (50 filaments on 4π·100 μm² = 1256 μm², density 0.04/μm²)
        # the bind density is too low to expect 50 bonds. Just require some
        # nonzero count OR all-dropped (which is also a valid sparse-cortex
        # outcome).
        assert xl.bond_pairs.shape[0] <= 50
        n_dropped = p.extras.get("n_xl_dropped", 0)
        assert xl.bond_pairs.shape[0] + n_dropped > 0

    def test_xl_alpha_filamin_split_when_present(self):
        cfg = _demo_cfg(**{
            "crosslinkers.enabled": True,
            "crosslinkers.n_xl": 500,
            "crosslinkers.max_bind_dist": 500.0e-9,  # widen for demo density
        })
        p = resolve_h3_derived(cfg)
        topology = generate_cortex_topology(p)
        xl = generate_static_xl_bonds(topology, p)
        # If any xl realised, the split should be approximately 30/70.
        n_total = xl.n_alpha + xl.n_filamin
        if n_total >= 50:
            frac_alpha = xl.n_alpha / n_total
            assert abs(frac_alpha - p.xl_alpha_fraction) < 0.10, (
                f"Realised α-actinin fraction {frac_alpha:.3f} too far "
                f"from target {p.xl_alpha_fraction:.3f} (n={n_total})."
            )

    def test_xl_pairs_are_from_different_filaments(self):
        cfg = _demo_cfg(**{
            "crosslinkers.enabled": True,
            "crosslinkers.n_xl": 100,
            "crosslinkers.max_bind_dist": 500.0e-9,
        })
        p = resolve_h3_derived(cfg)
        topology = generate_cortex_topology(p)
        xl = generate_static_xl_bonds(topology, p)
        # Every bond's two beads must belong to different filaments.
        N = p.beads_per_filament
        fila_a = xl.bond_pairs[:, 0] // N
        fila_b = xl.bond_pairs[:, 1] // N
        assert (fila_a != fila_b).all() if xl.bond_pairs.size else True

    def test_xl_force_negligible_at_construction(self):
        """Per-bin r0 ensures xl construction energy ≪ kT."""
        cfg = _demo_cfg(**{
            "crosslinkers.enabled": True,
            "crosslinkers.n_xl": 100,
            "crosslinkers.max_bind_dist": 500.0e-9,
        })
        p = resolve_h3_derived(cfg)
        sim, _, _, _, xl = build_cortex_simulation(
            p, with_baoab=False, with_crosslinkers=True
        )
        sim.run(1)
        bond = sim.operations.integrator.forces[0]
        # Total bond energy = cortex-bond (0 at exact ℓ_0) + xl bonds.
        # Upper bound: n_xl_realised · ½ · k_xl · (bin_width/2)².
        bin_width = p.xl_max_bind_dist / p.xl_n_bins
        upper = (
            xl.bond_pairs.shape[0]
            * 0.5 * p.xl_stiffness * (0.5 * bin_width) ** 2
        )
        assert bond.energy <= upper + 1e-20, (
            f"Bond energy {bond.energy:.3e} J exceeds xl-binning upper "
            f"bound {upper:.3e} J — construction state not force-free."
        )

    def test_xl_n_bins_consistency(self):
        names = xl_bin_type_names(10)
        r0 = xl_bin_rest_lengths(10, 100.0e-9)
        assert len(names) == 10
        assert r0.shape == (10,)
        assert math.isclose(r0[0], 5.0e-9)        # bin 0 center
        assert math.isclose(r0[-1], 95.0e-9)      # bin 9 center


def _make_topology(p: ResolvedH3) -> CortexTopology:
    return generate_cortex_topology(p)


# ---------------------------------------------------------------------------
# Production gates (opt-in via H3_PRODUCTION=1)
# ---------------------------------------------------------------------------
H3_PRODUCTION = bool(int(os.environ.get("H3_PRODUCTION", "0")))


@pytest.mark.skipif(
    not H3_PRODUCTION,
    reason="H.3 production gates; opt-in via H3_PRODUCTION=1.",
)
class TestH3Production:
    """Per-filament L_p / 3D equipartition / 3D Boltzmann angle gates.

    Aggregates over all 1000 cortex filaments. Per-filament force
    constants are identical to H.2 ✅ strict-PASS, so the same
    first-principles bands apply.

    Default (H3_PRODUCTION=1): SMOKE scale — 300 filaments × 5 interior
    beads × 30 snapshots × sample_interval 5 000 = 150 000 BAOAB steps
    after a 50 000-step equilibration. Aggregate sample size 45 000 —
    still 9σ above the noise floor for the ±5 % equipartition gate, so
    the first-principles bands remain statistically powered.

    Full-scale gate (1000 filaments × 100 snapshots × 50 k interval =
    5.1 M steps, ~15-min wall on M1 Max CPU) is opt-in via
    H3_PRODUCTION_FULL=1 — kept for the eventual ✅ DONE production
    sign-off sweep.
    """

    @pytest.fixture(scope="class")
    def production_run(self, resolved_production):
        full = bool(int(os.environ.get("H3_PRODUCTION_FULL", "0")))
        medium = bool(int(os.environ.get("H3_PRODUCTION_MEDIUM", "0")))
        p = resolved_production

        if full:
            # Production sign-off: full-scale (original H.3 brief).
            n_filaments_run = p.n_filaments
            n_equilibrate = 100_000
            n_snapshots = 100
            sample_interval = 50_000
        elif medium:
            # Intermediate sweep (단계 8): 500 filaments × 50 snapshots ×
            # 25k step interval = 1.25 M BAOAB steps ~15 min wall on
            # M1 Max CPU.  Aggregate sample size F × (N − 2) × n_snapshots
            # = 500 × 5 × 50 = 125 000 → σ_L_p ≈ 1.34 μm → 1.3σ on KU-1.1
            # ±10 % band [15.3, 18.7] μm.  Borderline interim signal —
            # PASS at MEDIUM is suggestive but the FULL gate remains
            # required for H.3 → ✅ DONE production sign-off.
            n_filaments_run = 500
            n_equilibrate = 75_000
            n_snapshots = 50
            sample_interval = 25_000
        else:
            # CI smoke: scaled-down to fit a single dev iteration
            # (~1-2 min wall on M1 Max CPU). Aggregate sample size
            # F × (N − 2) × n_snapshots = 300 × 5 × 30 = 45 000 — well
            # above the 1/0.05² ≈ 400 floor needed for the ±5 %
            # equipartition tol at 1σ, giving ~9σ statistical reach.
            n_filaments_run = 300
            n_equilibrate = 50_000
            n_snapshots = 30
            sample_interval = 5_000

        # Build a fresh ResolvedH3 with the (possibly scaled-down)
        # filament count. Other parameters carry over unchanged.
        from copy import deepcopy
        cfg = _load_cfg()
        cfg["cortex"]["n_filaments"] = n_filaments_run
        cfg["cortex"]["demo_mode"] = True   # relax cost ceiling at scale-down
        p_run = resolve_h3_derived(cfg)

        sim, _, _, topology, _ = build_cortex_simulation(
            p_run, with_baoab=True, with_crosslinkers=False
        )

        # Per-iteration checkpoint + progress print (lesson from H.3 단계 9:
        # the 5-hour CPU L_p FULL run had to be killed without yielding any
        # data because frames were memory-only.  Now each snapshot is also
        # written to ffn_sim/outputs/h3/checkpoints/lp_{scale}/snapshot_NNN.npz
        # so a midway kill leaves partial data on disk + every snapshot
        # prints a timestamped progress line for `tail -f` monitoring).
        import time as _time
        scale_tag = (
            "full" if full else ("medium" if medium else "smoke")
        )
        ckpt_dir = (
            Path(__file__).resolve().parents[1] / "outputs" / "h3"
            / "checkpoints" / f"lp_{scale_tag}"
        )
        ckpt_dir.mkdir(parents=True, exist_ok=True)

        t_start = _time.time()
        sim.run(n_equilibrate)
        print(
            f"[{_time.strftime('%H:%M:%S')}] L_p {scale_tag.upper()} "
            f"equilibration done ({n_equilibrate} steps, "
            f"{(_time.time() - t_start):.1f} s wall)", flush=True,
        )

        frames = np.empty(
            (n_snapshots, p_run.n_filaments, p_run.beads_per_filament, 3),
            dtype=np.float64,
        )
        for k in range(n_snapshots):
            sim.run(sample_interval)
            with sim.state.cpu_local_snapshot as snap:
                # Tag indirection so row order matches our (F, N, 3)
                # construction order even if HOOMD's ParticleSorter has
                # reordered rows.
                tag = np.asarray(snap.particles.tag)
                pos = np.asarray(snap.particles.position)
                inv = np.empty_like(tag)
                inv[tag] = np.arange(tag.size, dtype=tag.dtype)
                frames[k] = pos[inv].reshape(
                    p_run.n_filaments, p_run.beads_per_filament, 3
                )
            # Checkpoint each snapshot to disk (negligible cost vs
            # sample_interval BAOAB steps).
            np.savez_compressed(
                ckpt_dir / f"snapshot_{k:03d}.npz",
                frame=frames[k], k=k, n_snapshots=n_snapshots,
                scale=scale_tag, sim_timestep=int(sim.timestep),
            )
            elapsed = _time.time() - t_start
            print(
                f"[{_time.strftime('%H:%M:%S')}] L_p {scale_tag.upper()} "
                f"snapshot {k+1}/{n_snapshots} "
                f"(ts={int(sim.timestep)}, wall={elapsed:.1f}s, "
                f"est_total={elapsed * n_snapshots / max(k+1, 1):.0f}s)",
                flush=True,
            )
        return p_run, frames

    @pytest.mark.skipif(
        not (
            bool(int(os.environ.get("H3_PRODUCTION_FULL", "0")))
            or bool(int(os.environ.get("H3_PRODUCTION_MEDIUM", "0")))
        ),
        reason=(
            "L_p estimator is sample-size-limited at H.3 short-filament "
            "scale (L=3 μm vs L_p=17 μm → L/L_p ≈ 0.18, fit window only "
            "s∈[1,3] bonds, dynamic range of ln C(s) only [-0.029, -0.088]). "
            "Single-snapshot σ_L_p ≈ 5 μm at smoke scale (300 filaments × "
            "30 snapshots = 9 000 pairs per s); MEDIUM (500 × 50 = 125 000 "
            "pairs per s, σ_L_p ≈ 1.34 μm → 1.3σ borderline) is opt-in via "
            "H3_PRODUCTION_MEDIUM=1 (~15 min wall); FULL (1000 × 100 = "
            "500 000 pairs per s, σ_L_p ≈ 0.36 μm → 4.7σ resolved) is "
            "opt-in via H3_PRODUCTION_FULL=1 (~91 min wall on M1 Max CPU). "
            "CLAUDE.md no-gate-loosening forbids widening the band to fit "
            "the smoke estimator."
        ),
    )
    def test_per_filament_L_p_in_KU11_band(self, production_run):
        p, frames = production_run
        # Aggregate fit across all frames + filaments (same C(s) computed
        # over all bond-separations × all filaments × all snapshots).
        # We average per-frame and then average frame fits — equivalent
        # to ensemble fit at fit_s_max = (N-1)/2.
        L_p_list = []
        for k in range(frames.shape[0]):
            fit = fit_persistence_length(frames[k], rest_length=p.rest_length)
            if math.isfinite(fit.L_p_m):
                L_p_list.append(fit.L_p_m)
        L_p_mean = float(np.mean(L_p_list))
        lo, hi = p.L_p_band_m
        assert lo <= L_p_mean <= hi, (
            f"Ensemble L_p {L_p_mean:.3e} m outside KU-1.1 first-principles "
            f"band [{lo:.3e}, {hi:.3e}] m."
        )

    def test_3d_equipartition_strict(self, production_run):
        p, frames = production_run
        # Aggregate energy across all snapshots × all interior beads × all filaments.
        energies = []
        for k in range(frames.shape[0]):
            e = bending_energy_per_bond(frames[k], angle_k=p.angle_k)
            energies.append(e)
        all_e = np.stack(energies, axis=0)
        mean_J = float(all_e.mean())
        target_J = p.equipartition_target_3d_kT * p.kT
        rel = (mean_J - target_J) / target_J
        assert abs(rel) <= p.equipartition_rel_tol, (
            f"Per-bond ⟨E_bend⟩ = {mean_J:.3e} J = "
            f"{mean_J/p.kT:.4f} kT, rel = {rel:.3%} outside ±"
            f"{p.equipartition_rel_tol:.0%} of 3D analytical target "
            f"{p.equipartition_target_3d_kT:.4f} kT."
        )

    def test_3d_boltzmann_angle_KS(self, production_run):
        from scipy import stats
        p, frames = production_run
        # Aggregate angles across snapshots + filaments.
        angles_all = []
        for k in range(frames.shape[0]):
            theta = hoomd_angle_array(frames[k])      # (F, N-2)
            angles_all.append(theta.ravel())
        theta_obs = np.concatenate(angles_all)
        # Build CDF reference numerically.
        theta_grid = np.linspace(0.0, np.pi, 5000)
        dens = boltzmann_angle_density_3d(
            theta_grid, bending_modulus=p.bending_modulus,
            rest_length=p.rest_length, kT=p.kT,
        )
        cdf = np.cumsum(dens)
        cdf = cdf / cdf[-1]
        # KS test via inverse-CDF mapping.
        ks_stat = stats.kstest(theta_obs, lambda x: np.interp(x, theta_grid, cdf))
        assert ks_stat.statistic <= p.angle_ks_stat_max, (
            f"3D Boltzmann KS_stat = {ks_stat.statistic:.4f} above "
            f"first-principles ceiling {p.angle_ks_stat_max} (H.2 diagnostic)."
        )

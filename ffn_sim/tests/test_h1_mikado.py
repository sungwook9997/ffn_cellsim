"""STATIC sanity-gate tests for H.1 ECM Mikado (PHASE_0_3_DECISIONS D3+D4+D7).

Covers the Sanity Gate sections written in
``ffn_sim/ecm/mikado.py``'s module docstring:

- §1 Dimensional analysis  (TestDimensional)
- §2 Boundary cases        (TestBoundary)
- §3 Conservation invariants + topology counts (TestTopology, TestInitialEnergy)
- §5 Sign / sense          (TestSignSense)
- §6 Measurement protocol  (TestEnergyOracle — the H.1 brief's ≤ 1e-6 gate)

Empirical KU-1.30 (#1 G_0, #2 strain stiffening, #3 point dipole) lives
separately in ``tests/validation/test_ku130.py`` (Milestone 2). This
file is the **topology + energy-oracle freeze** test set: it must PASS
before any KU-1.30 production run.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd
import hoomd.md as md

from ffn_sim.ecm.mikado import (
    ResolvedH1,
    build_mikado_simulation,
    build_mikado_state,
    resolve_derived,
)
from ffn_sim.validation.oracles.ecm.fiber_mechanics import compute_energy
from ffn_sim.validation.oracles.ecm.fiber_network import (
    generate_2d_fiber_network,
)


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
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
def sim_triplet(resolved):
    """One full HOOMD setup with BAOAB attached; expensive — cache module-wide.

    M1 fixture: ``with_cross_links=False`` so bond / angle counts match the
    M1 topology contract. Cross-link-aware tests live in
    ``test_h1_cross_links.py`` with their own fixtures.
    """
    sim, updater, action = build_mikado_simulation(
        resolved, with_cross_links=False
    )
    sim.run(0)  # trigger force evaluation
    return sim, updater, action


@pytest.fixture(scope="module")
def force_eval_sim(resolved):
    """Force-evaluation-only Simulation (no BAOAB Updater attached).

    Used by the M1 energy-oracle and sign-sense tests so callers can
    write positions via ``cpu_local_snapshot`` and re-evaluate forces
    with ``sim.run(1)`` without an Updater advancing the state. HOOMD
    7's second ``run(0)`` does *not* re-evaluate forces (it caches), so
    ``run(1)`` is required to refresh ``force.energy`` after a write.
    M1 contract: cross-links disabled.
    """
    sim, _, _ = build_mikado_simulation(
        resolved, with_baoab=False, with_cross_links=False
    )
    sim.run(0)
    return sim


def _bonded_pe(sim: hoomd.Simulation) -> float:
    """Return HOOMD bond + angle potential energy (Joules), excluding LJ.

    The energy-oracle gate compares this against the oracle's
    ``compute_energy`` (which is WLC stretch + bend only — no LJ). The
    full ``ThermodynamicQuantities.potential_energy`` includes the LJ
    pair contribution, which is dominated by the inevitable inter-fiber
    bead overlaps in a random Mikado at construction (KU-1.30
    production runs use an equilibration phase to push these apart;
    the energy-oracle gate is BEFORE that equilibration).
    """
    bond_e = 0.0
    angle_e = 0.0
    for f in sim.operations.integrator.forces:
        cls = type(f).__name__
        if cls == "Harmonic":
            # md.bond.Harmonic and md.angle.Harmonic both share the
            # class name "Harmonic". Disambiguate via module path.
            mod = type(f).__module__
            if "bond" in mod:
                bond_e += float(f.energy)
            elif "angle" in mod:
                angle_e += float(f.energy)
    return bond_e + angle_e


def _read_pos_by_tag(sim) -> np.ndarray:
    """Read bead positions ordered by particle tag (stable identity).

    HOOMD's ParticleSorter reorders rows between calls; we always read
    in tag order so per-fiber bookkeeping (which the bond/angle indices
    rely on) stays valid.
    """
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tag = np.asarray(s.particles.tag)
        out = np.empty_like(pos)
        out[tag] = pos
        return out.copy()


# ---------------------------------------------------------------------------
# §1 Dimensional analysis
# ---------------------------------------------------------------------------
class TestDimensional:
    """The derivations in resolve_derived must reproduce SI-balanced values."""

    def test_rest_length_matches_L_fiber_over_N_minus_one(self, resolved):
        expected = resolved.L_fiber / (resolved.beads_per_fiber - 1)
        assert math.isclose(resolved.rest_length, expected, rel_tol=0, abs_tol=0)

    def test_bond_k_is_mu_over_ell0(self, resolved):
        expected = resolved.stretching_modulus / resolved.rest_length
        assert math.isclose(resolved.bond_k, expected, rel_tol=0, abs_tol=0)

    def test_angle_k_is_kappa_over_ell0(self, resolved):
        expected = resolved.bending_modulus / resolved.rest_length
        assert math.isclose(resolved.angle_k, expected, rel_tol=0, abs_tol=0)

    def test_gamma_b_is_stokes_6pi_eta_R(self, resolved):
        expected = 6.0 * math.pi * resolved.water_viscosity * resolved.bead_radius
        assert math.isclose(resolved.gamma_b, expected, rel_tol=0, abs_tol=0)

    def test_lj_sigma_is_2R(self, resolved):
        assert math.isclose(resolved.lj_sigma, 2.0 * resolved.bead_radius)

    def test_lj_r_cut_is_2pow1_6_sigma(self, resolved):
        assert math.isclose(
            resolved.lj_r_cut, 2.0 ** (1.0 / 6.0) * resolved.lj_sigma,
            rel_tol=1e-12, abs_tol=0,
        )

    def test_dt_cfl_is_alpha_times_tau_min(self, resolved):
        expected = resolved.cfl_safety_factor * resolved.tau_min
        assert math.isclose(resolved.dt_cfl, expected, rel_tol=0, abs_tol=0)

    def test_tau_min_is_minimum_of_three(self, resolved):
        assert resolved.tau_min == min(
            resolved.tau_xl, resolved.tau_stretch, resolved.tau_bend
        )


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_beads_per_fiber_below_2_raises(self):
        cfg = _load_cfg()
        cfg["ecm"]["beads_per_fiber"] = 1
        with pytest.raises(ValueError, match="beads_per_fiber"):
            resolve_derived(cfg)

    def test_dimensions_not_3_raises(self):
        cfg = _load_cfg()
        cfg["ecm"]["topology"]["dimensions"] = 2
        with pytest.raises(ValueError, match="dimensions must be 3"):
            resolve_derived(cfg)

    def test_non_baoab_integrator_raises(self):
        cfg = _load_cfg()
        cfg["ecm"]["dynamics"]["integrator"] = "euler_maruyama"
        with pytest.raises(ValueError, match="D3 canonical"):
            resolve_derived(cfg)

    def test_zero_Lz_raises(self):
        cfg = _load_cfg()
        cfg["ecm"]["topology"]["L_z_over_L_box"] = 0.0
        with pytest.raises(ValueError, match="L_z_over_L_box"):
            resolve_derived(cfg)


# ---------------------------------------------------------------------------
# §3 Topology / counts + §6 measurement: bond / angle / particle counts
# ---------------------------------------------------------------------------
class TestTopology:
    def test_particle_count(self, sim_triplet, resolved):
        sim, _, _ = sim_triplet
        n_expected = resolved.n_fibers * resolved.beads_per_fiber
        assert sim.state.N_particles == n_expected

    def test_bond_count(self, sim_triplet, resolved):
        sim, _, _ = sim_triplet
        n_expected = resolved.n_fibers * (resolved.beads_per_fiber - 1)
        assert sim.state.N_bonds == n_expected

    def test_angle_count(self, sim_triplet, resolved):
        sim, _, _ = sim_triplet
        n_expected = resolved.n_fibers * (resolved.beads_per_fiber - 2)
        assert sim.state.N_angles == n_expected

    def test_box_is_3d_cubic(self, sim_triplet, resolved):
        sim, _, _ = sim_triplet
        box = sim.state.box
        assert math.isclose(box.Lx, resolved.L_box)
        assert math.isclose(box.Ly, resolved.L_box)
        assert math.isclose(box.Lz, resolved.L_z)
        assert box.xy == 0 and box.xz == 0 and box.yz == 0

    def test_mikado_line_density_KU_1_27(self, resolved):
        # ρ_L = N L_f / L_box²
        rho_L_expected = resolved.n_fibers * resolved.L_fiber / resolved.L_box**2
        assert math.isclose(resolved.line_density, rho_L_expected, rel_tol=1e-15)
        # ℓ_c_predicted = π / (2 ρ_L)
        ell_c_expected = math.pi / (2.0 * rho_L_expected)
        assert math.isclose(
            resolved.mikado_ell_c_predicted, ell_c_expected, rel_tol=1e-15
        )

    def test_positions_inside_hoomd_box(self, sim_triplet, resolved):
        sim, _, _ = sim_triplet
        pos = _read_pos_by_tag(sim)
        # HOOMD box is [-L/2, L/2). The build shifts oracle's [0, L) by -L/2,
        # so values must satisfy |x| ≤ L/2 + 1e-9.
        L_half = 0.5 * resolved.L_box
        assert (np.abs(pos[:, 0]) <= L_half + 1e-9).all()
        assert (np.abs(pos[:, 1]) <= L_half + 1e-9).all()
        # z must be exactly 0 in the initial 2D-extruded slab.
        assert np.all(pos[:, 2] == 0.0)

    def test_positions_are_float64(self, sim_triplet):
        sim, _, _ = sim_triplet
        with sim.state.cpu_local_snapshot as s:
            assert np.asarray(s.particles.position).dtype == np.float64


# ---------------------------------------------------------------------------
# §3 + §6 Initial energy / Newton 3rd law
# ---------------------------------------------------------------------------
class TestInitialEnergy:
    """At the construction-time configuration (all bonds at ℓ₀, all chains
    straight), HOOMD potential_energy must be 0 to numerical precision,
    and the oracle compute_energy on the same 2D-projected positions
    must agree.
    """

    def test_hoomd_bonded_energy_is_zero_at_init(self, force_eval_sim, resolved):
        sim = force_eval_sim
        pe_bonded = _bonded_pe(sim)
        # All bonds at rest length and all chains straight → bond + angle
        # PE vanish to float64 round-off. LJ is NOT included here (see
        # _bonded_pe docstring) because at construction the random
        # Mikado has inter-fiber bead overlaps that produce a large LJ
        # contribution; that's a real physical state to be relaxed by
        # the production equilibration, not a topology error.
        # Acceptance: |bonded PE| ≤ kT · 1e-8 (well below the energy
        # scale kT = 4.28e-21 J — see Mikado §6 measurement gate).
        assert abs(pe_bonded) <= 1e-8 * resolved.kT, (
            f"HOOMD bond+angle PE at initial Mikado = {pe_bonded:e} J; "
            "expected ≤ 1e-8 kT for the all-straight construction."
        )

    def test_oracle_potential_energy_is_zero_at_init(self, resolved):
        net = generate_2d_fiber_network(
            L_box=resolved.L_box,
            n_fibers=resolved.n_fibers,
            L_fiber=resolved.L_fiber,
            beads_per_fiber=resolved.beads_per_fiber,
            seed=resolved.seed,
        )
        E = compute_energy(
            bead_positions=net.bead_positions,
            rest_length=resolved.rest_length,
            stretching_modulus=resolved.stretching_modulus,
            bending_modulus=resolved.bending_modulus,
            box_size=resolved.L_box,
        )
        assert abs(E) <= 1e-8 * resolved.kT, (
            f"Oracle compute_energy at initial Mikado = {E:e} J; "
            "expected ≤ 1e-8 kT for the all-straight construction."
        )


# ---------------------------------------------------------------------------
# §6 Energy oracle agreement (the H.1 brief's ≤ 1e-6 rel error gate)
# ---------------------------------------------------------------------------
class TestEnergyOracle:
    """At a small-amplitude in-plane perturbation of the Mikado, the HOOMD
    harmonic bond + harmonic angle potential energy must agree with the
    oracle's WLC stretch + (1-cos) bend potential energy to relative
    error ≤ 1e-6.

    Tolerance derivation: the angle forms differ by φ⁴/24 (from the
    Taylor expansion of (1-cos φ) − ½φ²). For per-angle deflection
    |φ| ≲ 2 mrad the per-angle relative error is φ²/12 ≲ 3·10⁻⁷, well
    inside the 1e-6 acceptance band. We perturb bead positions by a
    Gaussian with σ = 0.5 nm in xy and exactly 0 in z, keeping every
    bond stretch and every angle deflection in the small-amplitude
    regime.
    """

    def test_perturbed_hoomd_matches_oracle_1e_minus_6(
        self, resolved, force_eval_sim
    ):
        sim = force_eval_sim
        rng = np.random.default_rng(123)
        sigma_pos = 0.5e-9  # 0.5 nm — see tolerance derivation above

        F, N = resolved.n_fibers, resolved.beads_per_fiber
        # Read current (initial, straight) positions in tag order.
        pos_flat = _read_pos_by_tag(sim)
        delta = np.zeros_like(pos_flat)
        delta[:, :2] = rng.normal(scale=sigma_pos, size=(pos_flat.shape[0], 2))
        new_flat = pos_flat + delta

        # Write the perturbed positions back, respecting tag → row mapping.
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            np.asarray(s.particles.position)[:] = new_flat[tag]

        # HOOMD bond+angle energy at perturbed state (LJ excluded — see
        # _bonded_pe docstring). Use run(1) to force re-eval (run(0)
        # after the initial fixture run(0) does not refresh — HOOMD 7
        # caches forces between identical run(0) calls).
        sim.run(1)
        pe_hoomd = _bonded_pe(sim)

        # Oracle energy on the 2D projection of the same positions.
        pos2d = new_flat[:, :2].reshape(F, N, 2)
        # Shift back to [0, L) for tidiness; the oracle uses minimum-image
        # so this is mathematically a no-op, but matches the generator's
        # output convention.
        pos2d_oracle = (pos2d + 0.5 * resolved.L_box) % resolved.L_box
        pe_oracle = compute_energy(
            bead_positions=pos2d_oracle,
            rest_length=resolved.rest_length,
            stretching_modulus=resolved.stretching_modulus,
            bending_modulus=resolved.bending_modulus,
            box_size=resolved.L_box,
        )

        # Reset positions back to the initial state for any later test
        # in this module (force_eval_sim is module-scoped).
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            np.asarray(s.particles.position)[:] = pos_flat[tag]
        sim.run(1)

        assert pe_oracle > 0, (
            "Oracle energy at perturbed state should be positive; "
            f"got {pe_oracle:e} — perturbation may be too small to register."
        )
        rel = abs(pe_hoomd - pe_oracle) / abs(pe_oracle)
        assert rel <= 1e-6, (
            f"HOOMD vs oracle PE rel error = {rel:.3e} at σ_pos={sigma_pos:e}m; "
            f"HOOMD={pe_hoomd:e} J, oracle={pe_oracle:e} J. "
            "Acceptance ≤ 1e-6 per H.1 brief §Validation acceptance."
        )


# ---------------------------------------------------------------------------
# §5 Sign / sense
# ---------------------------------------------------------------------------
class TestSignSense:
    """Sign / sense tests on an LJ-free, single-fiber clean config.

    The full Mikado has random inter-fiber bead overlaps that produce
    large LJ forces at construction — those are physically real (KU-1.30
    production equilibrates them) but they swamp the per-bond / per-angle
    sign tests. We isolate sign behaviour on a single straight fiber,
    LJ disabled, no BAOAB Updater, force-eval only.
    """

    def _isolated_fiber(self, resolved):
        cfg = _load_cfg()
        cfg["ecm"]["excluded_volume"]["enabled"] = False
        cfg["ecm"]["demo_mode"] = True
        p = resolve_derived(cfg)
        p.n_fibers = 1
        # Single fiber → no inter-fiber cross-links possible; disable
        # explicitly so the test isolates bond+angle physics.
        sim, _, _ = build_mikado_simulation(
            p, with_baoab=False, with_cross_links=False
        )
        sim.run(0)
        return sim, p

    def test_stretched_bond_produces_attractive_force(self, resolved):
        """Stretch bead 1 along +x by δ; bead 0 must be pulled toward +x
        (attracted to bead 1) and bead 1 must be pulled toward -x."""
        sim, _ = self._isolated_fiber(resolved)

        # Move bead with tag=1 along +x by 10 nm (clearly above any
        # rounding noise; small enough to keep linear regime).
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            pos = np.asarray(s.particles.position)
            row1 = int(np.argwhere(tag == 1).item())
            pos[row1, 0] += 1.0e-8

        sim.run(1)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            F = np.asarray(s.particles.net_force).copy()
            r0 = int(np.argwhere(tag == 0).item())
            r1 = int(np.argwhere(tag == 1).item())
        assert F[r0, 0] > 0, (
            f"Bead 0 net_force_x = {F[r0, 0]:e}; expected > 0 (attracted "
            "toward stretched bead 1)."
        )
        assert F[r1, 0] < 0, (
            f"Bead 1 net_force_x = {F[r1, 0]:e}; expected < 0 (attracted "
            "back toward bead 0)."
        )

    def test_newton_third_law_on_isolated_fiber(self, resolved):
        """Σ F over a single straight fiber's beads ≈ 0 by Newton 3rd law."""
        sim, p = self._isolated_fiber(resolved)
        with sim.state.cpu_local_snapshot as s:
            F = np.asarray(s.particles.net_force).sum(axis=0)
        # On a clean isolated fiber at construction (bonds at rest length,
        # straight, no LJ), forces should be exactly zero to round-off.
        scale = p.kT / p.bead_radius  # natural force scale
        assert np.all(np.abs(F) <= 1e-6 * scale), (
            f"Isolated-fiber Σ F = {F}; expected ≤ 1e-6 · kT/R "
            f"(scale={scale:e} N)."
        )

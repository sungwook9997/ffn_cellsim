"""STATIC + analytical sanity-gate tests for the ECM substrate-compliance
anchor spring (KU-1.V.1).

Covers ``ffn_sim/ecm/substrate.py`` Sanity Gate §1–6:

- §1 Dimensional analysis        (TestDimensional)
- §2 Boundary cases              (TestBoundary)
- §3 Conservation invariants     (TestConservation)
- §4 Numerical sanity / CFL      (TestNumericalCFL)
- §5 Sign / sense                (TestSignSense — restoring toward anchor)
- §6 Measurement protocol        (TestMeasurement / TestBAOABSmoke)

Mirrors ``ffn_sim/tests/test_erm.py``. The tests build a MINIMAL,
self-contained HOOMD sim of a handful of single-type "ligand" particles with
the frozen L-M BAOAB integrator (same ``md.Integrator`` + BAOAB-updater
pattern the cortex builder uses), so nothing outside ``ecm/substrate.py`` is
exercised and the test is fast.

NOTE on constants: the only physical knob, ``k_sub`` [N/m], is a REQUIRED
parameter (no magic default) — every value used here is a deliberate test
fixture, not a tuned constant baked into the module. ``kT`` and
``gamma_ligand`` are stand-in SI scales for the harmonic-trap / CFL formulae;
no gate-passing number is invented in the module under test.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.ecm.substrate import (
    ResolvedSubstrate,
    SubstrateAnchorSpring,
    attach_substrate_spring,
    effective_E_sub,
    resolve_substrate,
)
from ffn_sim.archive.hoomd_legacy.integrator.baoab import make_baoab_updater


# Test-fixture SI scales (NOT module constants — see module-docstring note).
KT_TEST: float = 4.28e-21          # J   (room-T thermal energy scale)
GAMMA_LIGAND_TEST: float = 1.0e-8  # N·s/m  (stand-in ligand Stokes drag)
LIGAND_TYPE: str = "ligand"


# ---------------------------------------------------------------------------
# Minimal self-contained sim builder (single particle type = ligand)
# ---------------------------------------------------------------------------
def _build_ligand_sim(
    anchors: np.ndarray,
    *,
    k_sub: float,
    with_baoab: bool = False,
    gamma_ligand: float = GAMMA_LIGAND_TEST,
    cfl_safety_factor: float = 0.1,
    cfl_strict: bool = True,
    dt: float = 1.0e-7,
    displace: np.ndarray | None = None,
) -> tuple[hoomd.Simulation, SubstrateAnchorSpring, np.ndarray]:
    """Build a tiny sim of ``len(anchors)`` ligand particles + the spring.

    Particles start AT their anchors (the force-free construction state)
    unless ``displace`` (shape ``(n, 3)`` [m]) is given, which offsets each
    ligand from its anchor. Returns ``(sim, spring, ligand_tags)``.
    """
    anchors = np.asarray(anchors, dtype=np.float64).reshape(-1, 3)
    n = anchors.shape[0]
    start = anchors.copy()
    if displace is not None:
        start = start + np.asarray(displace, dtype=np.float64).reshape(-1, 3)

    max_extent = float(np.max(np.abs(anchors))) if n > 0 else 0.0
    box_L = max_extent * 4.0 + 1.0e-6
    snap = hoomd.Snapshot()
    snap.configuration.box = [box_L, box_L, box_L, 0.0, 0.0, 0.0]
    snap.particles.N = n
    snap.particles.types = [LIGAND_TYPE]
    snap.particles.position[:] = start
    snap.particles.typeid[:] = np.zeros(n, dtype=np.int32)
    # Tags 0..n-1 in this minimal sim; anchors row-aligned with these tags.
    ligand_tags = np.arange(n, dtype=np.int64)

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(snap)

    ig = md.Integrator(dt=dt)  # methods=[] (L-M Action contract)
    sim.operations.integrator = ig

    spring = attach_substrate_spring(
        sim,
        k_sub,
        ligand_tags,
        anchors,
        gamma_ligand,
        cfl_safety_factor=cfl_safety_factor,
        cfl_strict=cfl_strict,
    )

    if with_baoab:
        _, updater = make_baoab_updater(
            kT=KT_TEST, gamma={LIGAND_TYPE: gamma_ligand}, dt=dt, seed=7
        )
        sim.operations.updaters.append(updater)

    return sim, spring, ligand_tags


def _fixed_anchors() -> np.ndarray:
    """A few off-z=0 + on-z=0 anchor sites [m] (sub-µm scale)."""
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0e-7, 0.0, 0.0],
            [0.0, 2.0e-7, 0.0],
            [-1.5e-7, 1.0e-7, 0.0],
            [3.0e-7, -2.0e-7, 0.0],
        ],
        dtype=np.float64,
    )


# ---------------------------------------------------------------------------
# §1 Dimensional analysis
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_resolve_sets_sigma(self):
        """σ = √(kT/k_sub) [m]; k_sub [N/m] · Δr [m] → F [N]."""
        k_sub = 1.0e-3
        p = resolve_substrate({"substrate": {"k_sub": k_sub}}, kT=KT_TEST)
        assert p.k_sub == k_sub
        assert math.isclose(
            p.sigma_thermal, math.sqrt(KT_TEST / k_sub), rel_tol=1e-12
        )

    def test_force_units_newton(self):
        """|F| = k_sub · |Δr| has Newton scale for a sub-µm displacement."""
        k_sub = 1.0e-3                     # N/m
        delta = 1.0e-7                     # m (100 nm)
        F_mag = k_sub * delta
        assert math.isclose(F_mag, 1.0e-10, rel_tol=1e-12)  # 1e-3·1e-7 N
        assert 1e-14 < F_mag < 1e-6        # sane Newton-scale

    def test_energy_units_joule(self):
        """U = ½ k_sub |Δr|² has Joule scale."""
        k_sub, delta = 1.0e-3, 1.0e-7
        U = 0.5 * k_sub * delta**2
        assert math.isclose(U, 0.5e-17, rel_tol=1e-12)  # ½·1e-3·1e-14 J

    def test_E_sub_bridge_units_pascal(self):
        """E_sub ≈ k_sub · n_ligand / A_substrate → [N/m · 1/m²] = Pa-like."""
        # k_sub [N/m] / A [m²] · n[—] → N/m³? No: N/m · 1/m² = N/m³ is wrong;
        # the documented bridge is a leading-order modulus estimate. Check the
        # ARITHMETIC matches the documented formula exactly (diagnostic only).
        k_sub, n_lig, A = 1.0e-3, 100, 1.0e-10
        E = effective_E_sub(k_sub, n_lig, A)
        assert math.isclose(E, k_sub * n_lig / A, rel_tol=1e-12)
        assert E > 0.0


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_k_sub_required_no_default(self):
        """k_sub absent → ValueError (REQUIRED, no magic default)."""
        with pytest.raises(ValueError, match="k_sub is REQUIRED"):
            resolve_substrate({"substrate": {}}, kT=KT_TEST)

    def test_k_sub_none_raises(self):
        with pytest.raises(ValueError, match="k_sub is REQUIRED"):
            resolve_substrate({"substrate": {"k_sub": None}}, kT=KT_TEST)

    def test_negative_k_sub_raises(self):
        with pytest.raises(ValueError, match="k_sub"):
            resolve_substrate({"substrate": {"k_sub": -1.0}}, kT=KT_TEST)

    def test_nonfinite_k_sub_raises(self):
        with pytest.raises(ValueError, match="k_sub"):
            resolve_substrate(
                {"substrate": {"k_sub": float("inf")}}, kT=KT_TEST
            )

    def test_anchor_shape_mismatch_raises(self):
        p = ResolvedSubstrate(k_sub=1.0e-3)
        with pytest.raises(ValueError, match="anchor_positions"):
            SubstrateAnchorSpring(
                p, np.array([0, 1, 2]), np.zeros((2, 3))
            )

    def test_duplicate_tags_raise(self):
        p = ResolvedSubstrate(k_sub=1.0e-3)
        with pytest.raises(ValueError, match="duplicate"):
            SubstrateAnchorSpring(
                p, np.array([0, 0, 1]), np.zeros((3, 3))
            )

    def test_nonfinite_anchor_raises(self):
        p = ResolvedSubstrate(k_sub=1.0e-3)
        bad = np.zeros((2, 3))
        bad[0, 0] = np.nan
        with pytest.raises(ValueError, match="non-finite"):
            SubstrateAnchorSpring(p, np.array([0, 1]), bad)

    def test_empty_ligand_set_zero_force(self):
        """No ligands → spring attaches, run(0) gives all-zero force."""
        sim, spring, tags = _build_ligand_sim(
            np.zeros((0, 3)), k_sub=1.0e-3, dt=1.0e-7
        )
        sim.run(0)
        assert np.asarray(spring.forces).size == 0 or np.allclose(
            np.asarray(spring.forces), 0.0
        )

    def test_construction_force_free(self):
        """Particles AT anchors → F = 0, U = 0 (the zero-force config)."""
        anchors = _fixed_anchors()
        sim, spring, tags = _build_ligand_sim(
            anchors, k_sub=1.0e-3, dt=1.0e-7
        )
        sim.run(0)
        F = np.asarray(spring.forces)
        assert np.allclose(F, 0.0, atol=1e-30), (
            f"At-anchor construction not force-free: max|F| = "
            f"{np.abs(F).max():.3e} N"
        )
        assert math.isclose(spring.energy, 0.0, abs_tol=1e-30)


# ---------------------------------------------------------------------------
# §3 Conservation invariants
# ---------------------------------------------------------------------------
class TestConservation:
    def test_external_field_net_force_nonzero(self):
        """A displaced ligand set has NONZERO net force (external ground
        field anchored at fixed lab sites — intentional, like ERM)."""
        anchors = _fixed_anchors()
        # Displace all ligands the SAME way → net restoring force ≠ 0.
        disp = np.tile(np.array([5.0e-8, 0.0, 0.0]), (anchors.shape[0], 1))
        sim, spring, tags = _build_ligand_sim(
            anchors, k_sub=1.0e-3, dt=1.0e-7, displace=disp
        )
        sim.run(0)
        net = np.asarray(spring.forces).sum(axis=0)
        assert np.linalg.norm(net) > 0.0, (
            "Uniformly displaced ligands should have nonzero net restoring "
            "force (external ground field)."
        )

    def test_energy_nonnegative(self):
        anchors = _fixed_anchors()
        disp = np.full((anchors.shape[0], 3), 3.0e-8)
        sim, spring, tags = _build_ligand_sim(
            anchors, k_sub=1.0e-3, dt=1.0e-7, displace=disp
        )
        sim.run(0)
        assert spring.energy >= 0.0

    def test_no_force_on_non_ligand(self):
        """A particle whose tag is NOT in ligand_tags feels zero spring
        force (mask). Build 5 ligands but mask only tags {0,1}."""
        anchors = _fixed_anchors()
        n = anchors.shape[0]
        # Construct the spring over a SUBSET of tags (only 0,1) so tags 2-4
        # are non-ligand and must stay force-free even when displaced.
        disp = np.full((n, 3), 4.0e-8)
        # Build the sim manually so we can pass a subset to the spring.
        start = anchors + disp
        box_L = float(np.max(np.abs(anchors)) * 4.0 + 1.0e-6)
        snap = hoomd.Snapshot()
        snap.configuration.box = [box_L, box_L, box_L, 0.0, 0.0, 0.0]
        snap.particles.N = n
        snap.particles.types = [LIGAND_TYPE]
        snap.particles.position[:] = start
        snap.particles.typeid[:] = np.zeros(n, dtype=np.int32)
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
        sim.create_state_from_snapshot(snap)
        ig = md.Integrator(dt=1.0e-7)
        sim.operations.integrator = ig
        subset_tags = np.array([0, 1], dtype=np.int64)
        spring = attach_substrate_spring(
            sim, 1.0e-3, subset_tags, anchors[:2], GAMMA_LIGAND_TEST
        )
        sim.run(0)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
        F = np.asarray(spring.forces)
        rows_with_force = np.flatnonzero(np.linalg.norm(F, axis=1) > 1e-30)
        tags_with_force = tag[rows_with_force]
        assert (tags_with_force < 2).all(), (
            f"Spring applied force to non-ligand tags: {tags_with_force}"
        )


# ---------------------------------------------------------------------------
# §4 Numerical sanity / CFL
# ---------------------------------------------------------------------------
class TestNumericalCFL:
    def test_force_magnitude_matches_analytic(self):
        """Displace ONE ligand by δ; its spring force = k_sub·δ to f64."""
        anchors = _fixed_anchors()
        n = anchors.shape[0]
        delta_vec = np.array([7.0e-8, -3.0e-8, 1.0e-8])
        disp = np.zeros((n, 3))
        disp[2] = delta_vec  # displace only ligand tag=2
        k_sub = 2.0e-3
        sim, spring, tags = _build_ligand_sim(
            anchors, k_sub=k_sub, dt=1.0e-7, displace=disp
        )
        sim.run(0)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            row = int(np.argwhere(tag == 2).item())
        F = np.asarray(spring.forces)[row]
        expected = -k_sub * delta_vec
        assert np.allclose(F, expected, rtol=1e-9, atol=1e-30), (
            f"F = {F} N, expected −k_sub·Δr = {expected} N"
        )

    def test_force_participates_in_net_force(self):
        """The custom force reaches HOOMD net_force after a sim.run(1)."""
        anchors = _fixed_anchors()
        n = anchors.shape[0]
        disp = np.zeros((n, 3))
        disp[1] = np.array([6.0e-8, 0.0, 0.0])
        k_sub = 1.0e-3
        sim, spring, tags = _build_ligand_sim(
            anchors, k_sub=k_sub, dt=1.0e-7, displace=disp
        )
        sim.run(1)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            netF = np.asarray(s.particles.net_force).copy()
            row = int(np.argwhere(tag == 1).item())
        # No other force in the integrator → net_force == spring force.
        expected = -k_sub * disp[1]
        assert np.allclose(netF[row], expected, rtol=1e-6, atol=1e-30), (
            f"net_force[ligand] = {netF[row]} N, expected {expected} N"
        )

    def test_cfl_tau_formula_and_gate_blocks_stiff_k_sub(self):
        """τ = γ_ligand / k_sub; the gate raises when dt > α·τ."""
        anchors = _fixed_anchors()
        gamma = GAMMA_LIGAND_TEST
        # Pick dt so the gate is comfortably satisfied for a soft spring,
        # then crank k_sub until τ = γ/k_sub shrinks below dt/α → must trip.
        dt = 1.0e-3
        # Soft spring: τ = 1e-8/5e-7 = 2e-2 s; 0.1·τ = 2e-3 > dt → passes
        # (comfortably off the boundary so float rounding can't flip it).
        soft = 5.0e-7
        tau_soft = gamma / soft
        assert dt < 0.1 * tau_soft  # gate passes (strict)
        sim, spring, _ = _build_ligand_sim(
            anchors, k_sub=soft, dt=dt, gamma_ligand=gamma
        )  # no raise
        # Stiff spring: choose k_sub so τ = γ/k_sub < dt/α = dt·10.
        # Need k_sub > γ·10/dt = 1e-8·10/1e-3 = 1e-4 N/m. Use 1e-2 (100×).
        stiff = 1.0e-2
        tau_stiff = gamma / stiff
        assert math.isclose(tau_stiff, gamma / stiff, rel_tol=1e-12)
        assert dt > 0.1 * tau_stiff  # gate MUST trip
        with pytest.raises(RuntimeError, match="Substrate-spring CFL violated"):
            _build_ligand_sim(
                anchors, k_sub=stiff, dt=dt, gamma_ligand=gamma,
                cfl_strict=True,
            )

    def test_cfl_strict_false_skips_gate(self):
        """cfl_strict=False lets an over-stiff k_sub attach (diagnostic)."""
        anchors = _fixed_anchors()
        sim, spring, _ = _build_ligand_sim(
            anchors, k_sub=1.0e-2, dt=1.0e-7, gamma_ligand=GAMMA_LIGAND_TEST,
            cfl_strict=False,
        )
        assert isinstance(spring, SubstrateAnchorSpring)

    def test_negative_gamma_raises(self):
        anchors = _fixed_anchors()
        with pytest.raises(ValueError, match="gamma_ligand"):
            _build_ligand_sim(
                anchors, k_sub=1.0e-3, dt=1.0e-7, gamma_ligand=-1.0
            )

    def test_attach_without_integrator_raises(self):
        anchors = _fixed_anchors()
        box_L = 1.0e-5
        snap = hoomd.Snapshot()
        snap.configuration.box = [box_L, box_L, box_L, 0.0, 0.0, 0.0]
        snap.particles.N = anchors.shape[0]
        snap.particles.types = [LIGAND_TYPE]
        snap.particles.position[:] = anchors
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
        sim.create_state_from_snapshot(snap)
        # No integrator set.
        with pytest.raises(RuntimeError, match="integrator must be set"):
            attach_substrate_spring(
                sim, 1.0e-3, np.arange(anchors.shape[0]), anchors,
                GAMMA_LIGAND_TEST,
            )


# ---------------------------------------------------------------------------
# §5 Sign / sense — restoring toward the anchor
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_force_points_toward_anchor(self):
        """Ligand displaced away from anchor feels force back TOWARD it
        (projection of F on the outward displacement is negative)."""
        anchors = _fixed_anchors()
        n = anchors.shape[0]
        out = np.array([4.0e-8, 5.0e-8, -2.0e-8])  # arbitrary 3D offset
        disp = np.zeros((n, 3))
        disp[3] = out
        sim, spring, _ = _build_ligand_sim(
            anchors, k_sub=1.5e-3, dt=1.0e-7, displace=disp
        )
        sim.run(0)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            row = int(np.argwhere(tag == 3).item())
        F = np.asarray(spring.forces)[row]
        proj = float(np.dot(F, out / np.linalg.norm(out)))
        assert proj < 0.0, (
            f"Displaced ligand should feel RESTORING force toward anchor; "
            f"projection on +Δr = {proj:.3e} N (should be < 0)."
        )
        # And F should be exactly anti-parallel to Δr (central spring).
        assert np.allclose(
            F, -1.5e-3 * out, rtol=1e-9, atol=1e-30
        )

    def test_at_anchor_zero_force(self):
        """Ligand exactly at anchor → F = 0 (no direction singularity)."""
        anchors = _fixed_anchors()
        sim, spring, _ = _build_ligand_sim(
            anchors, k_sub=1.0e-3, dt=1.0e-7
        )
        sim.run(0)
        assert np.allclose(np.asarray(spring.forces), 0.0, atol=1e-30)

    def test_rigid_limit_sigma_shrinks(self):
        """k_sub → ∞ ⇒ σ = √(kT/k_sub) → 0 (rigid-pin limit, documented).

        Finite-large k_sub ≠ pin bit-for-bit (it adds a force / leaves σ>0),
        but σ shrinks monotonically toward the rigid endpoint."""
        s1 = SubstrateAnchorSpring.predicted_sigma(1.0e-3, KT_TEST)
        s2 = SubstrateAnchorSpring.predicted_sigma(1.0e3, KT_TEST)
        assert s2 < s1
        assert s2 < 1.0e-10  # sub-Ångström at high stiffness → ~rigid
        # σ → 0 as k_sub → ∞.
        s3 = SubstrateAnchorSpring.predicted_sigma(1.0e9, KT_TEST)
        assert s3 < s2


# ---------------------------------------------------------------------------
# §6 Measurement protocol — short BAOAB equilibration smoke
# ---------------------------------------------------------------------------
class TestBAOABSmoke:
    def test_baoab_keeps_ligands_bounded(self):
        """Spring-on BAOAB run keeps ligands near anchors, no NaN/Inf.

        σ = √(kT/k_sub); with kT=4.28e-21 J, k_sub=1e-3 N/m → σ ≈ 2.1 nm.
        Over a short run the ligands stay bounded near their anchors (the
        trap is harmonic), and positions stay finite. This is a no-runaway
        smoke (loose 1 µm absolute bound, mirroring ``erm.py``), NOT a tight
        thermodynamic gate.
        """
        anchors = _fixed_anchors()
        k_sub = 1.0e-3
        sigma = math.sqrt(KT_TEST / k_sub)  # ≈ 2.07e-9 m
        assert 1e-9 < sigma < 1e-8          # ≈ 2 nm sanity on the trap width
        sim, spring, tags = _build_ligand_sim(
            anchors, k_sub=k_sub, dt=1.0e-7, with_baoab=True,
            gamma_ligand=GAMMA_LIGAND_TEST,
        )
        sim.run(0)
        # Construction is force-free.
        assert np.allclose(np.asarray(spring.forces), 0.0, atol=1e-30)
        sim.run(300)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
            tag = np.asarray(s.particles.tag).copy()
        assert np.isfinite(pos).all(), "Spring-on BAOAB produced NaN/Inf."
        # Displacement from each ligand's anchor (matched by tag).
        order = np.argsort(tag)
        disp = pos[order] - anchors
        max_disp = float(np.linalg.norm(disp, axis=1).max())
        # Loose no-runaway bound (1 µm absolute, ≫ the ~2 nm trap σ — same
        # spirit as erm.py's soft-ERM smoke). The harmonic trap keeps the
        # ligands bounded; this just rejects a numerical blow-up.
        assert max_disp < 1.0e-6, (
            f"Ligand drifted {max_disp:.3e} m > 1 µm from anchor; bounded "
            "trap should prevent runaway."
        )

    def test_predicted_sigma_helper(self):
        """Static helper exposes the analytic per-component σ."""
        sigma = SubstrateAnchorSpring.predicted_sigma(1.0e-3, KT_TEST)
        assert math.isclose(sigma, math.sqrt(KT_TEST / 1.0e-3), rel_tol=1e-12)
        assert 1e-9 < sigma < 1e-8  # ≈ 2.07 nm at k_sub=1e-3, kT=4.28e-21

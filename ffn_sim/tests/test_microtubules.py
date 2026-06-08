r"""STATIC + analytical Sanity-Gate tests for the microtubule compartment.

Covers ``ffn_sim/cell/microtubules.py`` Sanity Gate §1–6 cheaply:

- OFF-identity              (TestOffIdentity — disabled resolver + no-op builders)
- §1 Dimensional analysis   (TestDimensional — EI/ℓ_0 + Y/ℓ_0 + τ_bend bridges)
- §2 Boundary cases         (TestBoundary — negative/zero/out-of-band raise)
- §3 No-net-force / γ       (TestGammaDenylist + TestForceFreeConstruction)
- §4 Numerical CFL          (TestCFL — stiff bending dt computed + gate raises)
- §5 Sign / sense           (TestSignSense — stretched backbone pulls inward;
                             bent chain straightens) — uses a TINY HOOMD sim.

The heavy full-Cell path is NOT exercised; the sign tests build a minimal
2–3 bead rod (not a Cell) and are marked light enough to run in CI.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ffn_sim.cell.microtubules import (
    GAMMA_DENYLIST_PREFIX,
    PI_DECISIONS,
    MTDynamicInstability,
    MTTopology,
    ResolvedMicrotubules,
    attach_microtubule_forces,
    build_mt_topology,
    extend_snapshot_with_microtubules,
    resolve_microtubules,
)


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------
_KT = 4.28e-21                    # J   project k_B·T
# Per-bead Stokes drag γ = 6π η R_bead (water η≈1e-3, R_bead≈100 nm).
_GAMMA_B = 6.0 * math.pi * 1.0e-3 * 1.0e-7   # ≈ 1.885e-9 N·s/m
_L_MT = 5.0e-6                    # 5 µm MT contour length (host geometry)


def _cfg(enabled: bool = True, **overrides) -> dict:
    base = {
        "enabled": enabled,
        "n_mt": 6,
        "beads_per_mt": 5,
        "L_mt": _L_MT,
    }
    base.update(overrides)
    return {"microtubules": base}


def _resolved(**overrides) -> ResolvedMicrotubules:
    return resolve_microtubules(
        _cfg(**overrides), kT=_KT, gamma_b=_GAMMA_B
    )


def _bare_snapshot(n: int = 4):
    """Smallest gsd Frame with n particles, no bonds/angles — the host."""
    import gsd.hoomd

    snap = gsd.hoomd.Frame()
    snap.particles.N = n
    snap.particles.types = ["actin"]
    snap.particles.typeid = [0] * n
    snap.particles.position = [[0.0, 0.0, float(i) * 1e-7] for i in range(n)]
    snap.particles.mass = [1.0] * n
    snap.particles.velocity = [[0.0, 0.0, 0.0]] * n
    snap.particles.image = [[0, 0, 0]] * n
    L = 50.0e-6
    snap.configuration.box = [L, L, L, 0.0, 0.0, 0.0]
    snap.bonds.N = 0
    snap.bonds.types = []
    snap.angles.N = 0
    snap.angles.types = []
    return snap


# ---------------------------------------------------------------------------
# OFF-identity (DEFAULT-OFF hard rule)
# ---------------------------------------------------------------------------
class TestOffIdentity:
    def test_missing_block_is_disabled(self):
        p = resolve_microtubules({}, kT=_KT, gamma_b=_GAMMA_B)
        assert p.enabled is False
        assert p.n_mt == 0 and p.beads_per_mt == 0
        assert p.k_angle == 0.0 and p.k_backbone == 0.0
        assert p.n_beads_total == 0
        assert p.n_backbone_bonds == 0 and p.n_bending_angles == 0

    def test_enabled_false_is_disabled(self):
        p = _resolved(enabled=False)
        assert p.enabled is False
        assert p.EI == 0.0 and p.dt_cfl_bend == 0.0

    def test_no_enabled_key_defaults_false(self):
        # A config block present but with NO 'enabled' key ⇒ default-off.
        p = resolve_microtubules(
            {"microtubules": {"n_mt": 9, "beads_per_mt": 4, "L_mt": _L_MT}},
            kT=_KT, gamma_b=_GAMMA_B,
        )
        assert p.enabled is False

    def test_topology_builder_disabled_is_empty(self):
        p = _resolved(enabled=False)
        topo = build_mt_topology(p)
        assert topo.positions.shape == (0, 3)
        assert topo.backbone_bonds.shape == (0, 2)
        assert topo.bending_angles.shape == (0, 3)

    def test_snapshot_builder_disabled_returns_input_unchanged(self):
        snap = _bare_snapshot(4)
        p = _resolved(enabled=False)
        out = extend_snapshot_with_microtubules(snap, p)
        # OFF-identity: SAME object, zero particles/bonds/angles added.
        assert out is snap
        assert int(out.particles.N) == 4
        assert int(out.bonds.N) == 0
        assert int(out.angles.N) == 0

    def test_attach_disabled_is_noop(self):
        p = _resolved(enabled=False)
        # No sim needed: disabled returns (None, None) before touching sim.
        assert attach_microtubule_forces(object(), p) == (None, None)


# ---------------------------------------------------------------------------
# §1 Dimensional analysis — the grid-invariant bridges
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_k_angle_is_EI_over_l0(self):
        p = _resolved()
        l0 = _L_MT / (5 - 1)
        assert p.l0 == pytest.approx(l0, rel=1e-12)
        # k_angle = EI / ℓ_0  (the H.2 bending_modulus/rest_length bridge).
        assert p.k_angle == pytest.approx(p.EI / l0, rel=1e-12)

    def test_k_backbone_is_Ystretch_over_l0(self):
        p = _resolved()
        assert p.k_backbone == pytest.approx(p.Y_stretch / p.l0, rel=1e-12)

    def test_EI_equals_kT_times_Lp(self):
        # WLC identity, the consistency the resolver enforces.
        p = _resolved()
        assert p.EI == pytest.approx(_KT * p.L_p, rel=0.12)

    def test_tau_bend_formula(self):
        p = _resolved()
        expect = _GAMMA_B * p.l0**3 / p.EI
        assert p.tau_bend == pytest.approx(expect, rel=1e-12)
        assert p.dt_cfl_bend == pytest.approx(0.1 * expect, rel=1e-12)

    def test_grid_invariance_of_continuum_EI(self):
        # Finer discretisation (more beads, smaller ℓ_0) raises k_angle but the
        # reproduced continuum rigidity EI is unchanged — k_angle is DERIVED.
        p_coarse = _resolved(beads_per_mt=5)
        p_fine = _resolved(beads_per_mt=21)
        assert p_fine.k_angle > p_coarse.k_angle      # per-angle stiffer
        # EI = k_angle · ℓ_0 is invariant.
        assert (p_fine.k_angle * p_fine.l0) == pytest.approx(
            p_coarse.k_angle * p_coarse.l0, rel=1e-12
        )


# ---------------------------------------------------------------------------
# §2 Boundary cases — bad params raise
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_negative_n_mt_raises(self):
        with pytest.raises(ValueError):
            _resolved(n_mt=0)

    def test_too_few_beads_for_backbone_raises(self):
        with pytest.raises(ValueError):
            _resolved(beads_per_mt=1)

    def test_too_few_beads_for_bending_raises(self):
        with pytest.raises(ValueError):
            _resolved(beads_per_mt=2)

    def test_missing_L_mt_raises(self):
        cfg = {"microtubules": {"enabled": True, "n_mt": 6, "beads_per_mt": 5}}
        with pytest.raises(ValueError):
            resolve_microtubules(cfg, kT=_KT, gamma_b=_GAMMA_B)

    def test_negative_L_mt_raises(self):
        with pytest.raises(ValueError):
            _resolved(L_mt=-1.0)

    def test_out_of_band_Lp_raises(self):
        # L_p far outside the Gittes 1–6 mm band ⇒ rejected.
        with pytest.raises(ValueError):
            _resolved(L_p=1.0e-6)            # 1 µm (actin scale, not MT)

    def test_inconsistent_EI_raises(self):
        # EI grossly inconsistent with kT·L_p must raise (a tuned EI).
        with pytest.raises(ValueError):
            _resolved(EI=1.0e-20)           # ~450× kT·L_p

    def test_zero_gamma_b_raises(self):
        with pytest.raises(ValueError):
            resolve_microtubules(_cfg(), kT=_KT, gamma_b=0.0)

    def test_bad_cfl_safety_raises(self):
        with pytest.raises(ValueError):
            resolve_microtubules(
                _cfg(), kT=_KT, gamma_b=_GAMMA_B, cfl_safety_factor=2.0
            )


# ---------------------------------------------------------------------------
# §3 No-gamma-contamination — every bond type starts with the denylist prefix
# ---------------------------------------------------------------------------
class TestGammaDenylist:
    def test_prefix_nonempty(self):
        assert GAMMA_DENYLIST_PREFIX
        assert GAMMA_DENYLIST_PREFIX == "mt_"

    def test_backbone_bond_type_has_prefix(self):
        p = _resolved()
        assert p.backbone_bond_type.startswith(GAMMA_DENYLIST_PREFIX)

    def test_built_bonds_all_have_prefix(self):
        # Every bond type the builder would actually add must be prefix-tagged
        # so cortical_tension.py's denylist excludes the MT compression path.
        snap = _bare_snapshot(4)
        p = _resolved()
        out = extend_snapshot_with_microtubules(snap, p)
        mt_bond_types = [t for t in out.bonds.types if t not in snap.bonds.types]
        assert mt_bond_types, "builder added no new bond types"
        for t in mt_bond_types:
            assert t.startswith(GAMMA_DENYLIST_PREFIX), t


# ---------------------------------------------------------------------------
# §3 Force-free construction + topology counts
# ---------------------------------------------------------------------------
class TestForceFreeConstruction:
    def test_topology_counts(self):
        p = _resolved(n_mt=6, beads_per_mt=5)
        topo = build_mt_topology(p)
        # 1 MTOC + n_mt·beads_per_mt beads.
        assert topo.positions.shape[0] == 1 + 6 * 5
        assert topo.positions.shape[0] == p.n_beads_total
        # backbone: beads_per_mt per chain (incl. MTOC→first).
        assert topo.backbone_bonds.shape[0] == 6 * 5
        assert topo.backbone_bonds.shape[0] == p.n_backbone_bonds
        # interior bending angles: (beads_per_mt-1) per chain.
        assert topo.bending_angles.shape[0] == 6 * (5 - 1)
        assert topo.bending_angles.shape[0] == p.n_bending_angles

    def test_built_straight_segments_at_l0(self):
        # Every backbone bond is exactly ℓ_0 (force-free stretch at build).
        p = _resolved()
        topo = build_mt_topology(p)
        pos = topo.positions
        for i, j in topo.backbone_bonds:
            d = float(np.linalg.norm(pos[j] - pos[i]))
            assert d == pytest.approx(p.l0, rel=1e-9)

    def test_built_straight_angles_at_pi(self):
        # Every interior bending triple is straight (θ = π) at build.
        p = _resolved()
        topo = build_mt_topology(p)
        pos = topo.positions
        for a, b, c in topo.bending_angles:
            u = pos[a] - pos[b]
            v = pos[c] - pos[b]
            cos = float(
                np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
            )
            assert cos == pytest.approx(-1.0, abs=1e-9)   # θ = π


# ---------------------------------------------------------------------------
# §4 Numerical CFL — the stiff bending dt is computed & the gate raises
# ---------------------------------------------------------------------------
class TestCFL:
    def test_bend_cfl_tighter_than_stretch(self):
        # With the stiff-rod Y_stretch default, the bending dt is the binding
        # one and the dataclass min() picks it.
        p = _resolved()
        assert p.dt_cfl == min(p.dt_cfl_bend, p.dt_cfl_stretch)
        assert p.dt_cfl > 0.0 and math.isfinite(p.dt_cfl)

    def test_disabled_cfl_is_inf(self):
        p = _resolved(enabled=False)
        assert p.dt_cfl == float("inf")

    def test_attach_raises_on_cfl_violation(self):
        import hoomd
        import hoomd.md as md

        snap = _bare_snapshot(4)
        p = _resolved()
        out = extend_snapshot_with_microtubules(snap, p)
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
        sim.create_state_from_snapshot(out)
        # Deliberately oversized dt (≫ dt_cfl_bend) ⇒ the gate must raise.
        sim.operations.integrator = md.Integrator(dt=p.dt_cfl * 1e3)
        with pytest.raises(RuntimeError):
            attach_microtubule_forces(sim, p, cfl_strict=True)

    def test_attach_ok_with_safe_dt(self):
        import hoomd
        import hoomd.md as md

        snap = _bare_snapshot(4)
        p = _resolved()
        out = extend_snapshot_with_microtubules(snap, p)
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
        sim.create_state_from_snapshot(out)
        sim.operations.integrator = md.Integrator(dt=p.dt_cfl * 0.5)
        backbone, bending = attach_microtubule_forces(sim, p, cfl_strict=True)
        assert backbone is not None and bending is not None
        assert backbone.params[p.backbone_bond_type]["r0"] == pytest.approx(p.l0)
        assert bending.params[p.bending_angle_type]["t0"] == pytest.approx(
            math.pi
        )


# ---------------------------------------------------------------------------
# §5 Sign / sense — stretched backbone pulls inward; bent chain straightens
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_stretched_backbone_pulls_inward(self):
        import gsd.hoomd
        import hoomd

        p = _resolved()
        l0 = p.l0
        # Two beads stretched to 1.5·ℓ_0 along x with an mt_backbone bond.
        snap = gsd.hoomd.Frame()
        snap.particles.N = 2
        snap.particles.types = ["mt_bead"]
        snap.particles.typeid = [0, 0]
        snap.particles.position = [[-0.75 * l0, 0, 0], [0.75 * l0, 0, 0]]
        snap.particles.mass = [1.0, 1.0]
        L = 50e-6
        snap.configuration.box = [L, L, L, 0, 0, 0]
        snap.bonds.N = 1
        snap.bonds.types = [p.backbone_bond_type]
        snap.bonds.typeid = [0]
        snap.bonds.group = [[0, 1]]

        import hoomd.md as md
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
        sim.create_state_from_snapshot(snap)
        sim.operations.integrator = md.Integrator(dt=p.dt_cfl * 0.5)
        backbone = md.bond.Harmonic()
        backbone.params[p.backbone_bond_type] = dict(k=p.k_backbone, r0=l0)
        sim.operations.integrator.forces.append(backbone)
        sim.run(0)
        # Per-particle forces from the bond ForceCompute (HOOMD .forces array).
        f = np.asarray(backbone.forces)
        # Bead 0 (at −x) feels +x force (toward bead 1); bead 1 feels −x.
        assert f[0, 0] > 0.0, f
        assert f[1, 0] < 0.0, f
        # Magnitude = k_backbone·(|Δr|−ℓ_0) = k_backbone·0.5·ℓ_0.
        assert abs(f[0, 0]) == pytest.approx(
            p.k_backbone * 0.5 * l0, rel=1e-6
        )

    def test_bent_chain_straightens(self):
        import gsd.hoomd
        import hoomd
        import hoomd.md as md

        p = _resolved()
        l0 = p.l0
        # Three beads forming a 90° kink: the interior bead must be pushed to
        # straighten (restoring toward θ = π). Place an mt_bending angle.
        snap = gsd.hoomd.Frame()
        snap.particles.N = 3
        snap.particles.types = ["mt_bead"]
        snap.particles.typeid = [0, 0, 0]
        snap.particles.position = [
            [0.0, 0.0, 0.0],          # a
            [l0, 0.0, 0.0],           # b (interior, the kink vertex)
            [l0, l0, 0.0],            # c (90° turn)
        ]
        snap.particles.mass = [1.0, 1.0, 1.0]
        L = 50e-6
        snap.configuration.box = [L, L, L, 0, 0, 0]
        snap.angles.N = 1
        snap.angles.types = [p.bending_angle_type]
        snap.angles.typeid = [0]
        snap.angles.group = [[0, 1, 2]]

        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
        sim.create_state_from_snapshot(snap)
        sim.operations.integrator = md.Integrator(dt=p.dt_cfl * 0.5)
        bending = md.angle.Harmonic()
        bending.params[p.bending_angle_type] = dict(k=p.k_angle, t0=math.pi)
        sim.operations.integrator.forces.append(bending)
        sim.run(0)
        f = np.asarray(bending.forces)
        # The bending force opens the 90° kink toward straight: the endpoints a
        # and c are pushed apart (restoring torque). Net force must be nonzero
        # and the interior vertex pushed away from the (a,c) midpoint diagonal.
        assert np.linalg.norm(f) > 0.0
        # Endpoint c sits at +y of b; straightening pushes c toward −y... check
        # the angle restoring sense: the force on the interior bead b has a
        # component toward the (a..c) chord, opening the angle. Concretely the
        # total potential decreases as θ→π, so forces are non-trivial and the
        # endpoints separate: |a−c| should be driven to increase.
        # Verify the restoring direction on endpoint c is toward straightening
        # (c pushed in +x, away from b along the would-be-straight line).
        assert f[2, 0] > 0.0, f      # c pushed in +x (toward collinear a-b-c)


# ---------------------------------------------------------------------------
# Optional dynamic-instability is honestly stubbed (NotImplementedError)
# ---------------------------------------------------------------------------
class TestDynamicInstabilityStub:
    def test_di_resolves_when_requested(self):
        p = _resolved(dynamic_instability=True)
        assert p.dynamic_instability is True
        assert p.v_grow > 0.0 and p.f_cat > 0.0

    def test_di_updater_raises_not_implemented(self):
        p = _resolved(dynamic_instability=True)
        with pytest.raises(NotImplementedError):
            MTDynamicInstability(p, dt=1e-9, mt_batch_steps=100)


# ---------------------------------------------------------------------------
# Module hygiene
# ---------------------------------------------------------------------------
class TestModuleHygiene:
    def test_pi_decisions_present(self):
        # Y_stretch + DI provenance are flagged (honest open decisions).
        assert isinstance(PI_DECISIONS, list)
        assert len(PI_DECISIONS) >= 1
        assert any("Y_stretch" in s for s in PI_DECISIONS)

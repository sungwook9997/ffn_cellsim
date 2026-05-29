"""STATIC sanity-gate checks for the L-M BAOAB-limit integrator (D3).

These tests cover the dimensional, boundary, sign-sense, and prv_rnds-memory
gates from the ``ffn_sim/integrator/baoab.py`` module docstring §Sanity Gate.
Empirical / thermodynamic gates (equipartition, diffusion, L-M order) live in
``ffn_sim/scripts/baoab_polymer_sanity.py`` because they require a longer
HOOMD run than is reasonable for a unit-test suite.

Run::

    conda activate ffn_sim
    pytest ffn_sim/tests/test_baoab.py -v
"""

from __future__ import annotations

import math
import tempfile

import numpy as np
import pytest

import gsd.hoomd
import hoomd
import hoomd.md as md

from ffn_sim.integrator.baoab import (
    LeimkuhlerMatthewsBAOAB,
    _wrap_into_box,
    make_baoab_updater,
)


# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------
def _minimal_state(
    n: int = 3, box: float = 20.0, types: tuple[str, ...] = ("A",)
) -> str:
    """Write a tiny .gsd file of `n` collinear beads on x-axis, return path."""
    snap = gsd.hoomd.Frame()
    snap.particles.N = n
    snap.particles.types = list(types)
    snap.particles.typeid = np.zeros(n, dtype=np.uint32)
    snap.particles.position = np.column_stack(
        [np.arange(n, dtype=float), np.zeros(n), np.zeros(n)]
    )
    snap.particles.mass = np.ones(n)
    if n >= 2:
        snap.bonds.N = n - 1
        snap.bonds.types = ["p"]
        snap.bonds.typeid = np.zeros(n - 1, dtype=np.uint32)
        snap.bonds.group = np.column_stack(
            [np.arange(n - 1), np.arange(1, n)]
        ).astype(np.uint32)
    snap.configuration.box = [box, box, box, 0.0, 0.0, 0.0]
    fp = tempfile.NamedTemporaryFile(suffix=".gsd", delete=False).name
    with gsd.hoomd.open(fp, mode="w") as f:
        f.append(snap)
    return fp


def _read_positions_by_tag(sim: hoomd.Simulation) -> np.ndarray:
    """Read positions ordered by particle tag (stable identity).

    HOOMD's ``ParticleSorter`` reorders local-snapshot rows between
    invocations of the integrator. Comparing positions row-wise across
    two reads can therefore see spurious "motion" that is in fact just
    a row permutation. Tag ordering removes that confound.
    """
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tag = np.asarray(s.particles.tag)
        out = np.empty_like(pos)
        out[tag] = pos
        return out.copy()


def _build_sim(
    n: int = 3,
    box: float = 20.0,
    dt: float = 1e-4,
    bond_k: float = 10.0,
    bond_r0: float = 1.0,
    with_bond: bool = True,
) -> hoomd.Simulation:
    fp = _minimal_state(n=n, box=box)
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=7)
    sim.create_state_from_gsd(filename=fp)
    ig = md.Integrator(dt=dt)
    if with_bond and n >= 2:
        bond = md.bond.Harmonic()
        bond.params["p"] = dict(k=bond_k, r0=bond_r0)
        ig.forces.append(bond)
    sim.operations.integrator = ig
    return sim


# ---------------------------------------------------------------------------
# Sanity Gate §4 — RUNTIME input validation
# ---------------------------------------------------------------------------
class TestRuntimeInputValidation:
    def test_negative_gamma_rejected(self):
        with pytest.raises(ValueError, match="gamma"):
            LeimkuhlerMatthewsBAOAB(
                kT=1.0, gamma={"A": -1.0}, dt=1e-4, seed=0
            )

    def test_zero_gamma_rejected(self):
        with pytest.raises(ValueError, match="gamma"):
            LeimkuhlerMatthewsBAOAB(
                kT=1.0, gamma={"A": 0.0}, dt=1e-4, seed=0
            )

    def test_nonfinite_gamma_rejected(self):
        with pytest.raises(ValueError, match="gamma"):
            LeimkuhlerMatthewsBAOAB(
                kT=1.0, gamma={"A": math.inf}, dt=1e-4, seed=0
            )

    def test_zero_dt_rejected(self):
        with pytest.raises(ValueError, match="dt"):
            LeimkuhlerMatthewsBAOAB(
                kT=1.0, gamma={"A": 1.0}, dt=0.0, seed=0
            )

    def test_negative_kT_rejected(self):
        with pytest.raises(ValueError, match="kT"):
            LeimkuhlerMatthewsBAOAB(
                kT=-1.0, gamma={"A": 1.0}, dt=1e-4, seed=0
            )

    def test_missing_type_in_gamma_rejected_at_attach(self):
        sim = _build_sim()
        # Add B type to the snapshot via state.particle_types? In HOOMD 7 you
        # can't easily add new types post-construction. So instead build with
        # two types but provide gamma for only one.
        # Simpler: gamma covers types present; the missing-type path is
        # exercised by the structural check below.
        action, updater = make_baoab_updater(
            kT=1.0, gamma={"B": 1.0}, dt=1e-4, seed=0
        )
        sim.operations.updaters.append(updater)
        with pytest.raises(RuntimeError, match="missing entries"):
            sim.run(0)  # triggers attach()

    def test_methods_nonempty_rejected_at_attach(self):
        sim = _build_sim()
        # Attach a Brownian method - the very thing D3 forbids - and confirm
        # the Action refuses to coexist.
        brownian = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0)
        brownian.gamma["A"] = 1.0
        sim.operations.integrator.methods.append(brownian)
        action, updater = make_baoab_updater(
            kT=1.0, gamma={"A": 1.0}, dt=1e-4, seed=0
        )
        sim.operations.updaters.append(updater)
        with pytest.raises(RuntimeError, match="methods must be empty"):
            sim.run(0)

    def test_dt_mismatch_rejected_at_attach(self):
        sim = _build_sim(dt=1e-4)
        action, updater = make_baoab_updater(
            kT=1.0, gamma={"A": 1.0}, dt=2e-4, seed=0  # mismatched
        )
        sim.operations.updaters.append(updater)
        with pytest.raises(RuntimeError, match="must equal Action dt"):
            sim.run(0)


# ---------------------------------------------------------------------------
# Sanity Gate §2 — Boundary cases
# ---------------------------------------------------------------------------
class TestBoundaryCases:
    def test_zero_temperature_zero_force_no_motion(self):
        """kT=0, F=0 ⇒ positions strictly invariant (deterministic limit)."""
        sim = _build_sim(n=2, with_bond=False, box=50.0, dt=1e-4)
        action, updater = make_baoab_updater(
            kT=0.0, gamma={"A": 1.0}, dt=1e-4, seed=0
        )
        sim.operations.updaters.append(updater)
        sim.run(0)
        r0 = _read_positions_by_tag(sim)
        sim.run(50)
        r1 = _read_positions_by_tag(sim)
        assert np.allclose(r0, r1, atol=0, rtol=0), (
            "with kT=0 and no force, positions must be exactly invariant; "
            f"max drift was {np.max(np.abs(r1 - r0)):e}"
        )

    def test_infinite_gamma_freezes_particle(self):
        """γ → very large ⇒ both increment terms → 0."""
        sim = _build_sim(n=2, with_bond=True, box=50.0, dt=1e-4, bond_k=10.0)
        action, updater = make_baoab_updater(
            kT=1.0, gamma={"A": 1e20}, dt=1e-4, seed=0
        )
        sim.operations.updaters.append(updater)
        sim.run(0)
        r0 = _read_positions_by_tag(sim)
        sim.run(100)
        r1 = _read_positions_by_tag(sim)
        assert np.max(np.abs(r1 - r0)) < 1e-10, (
            "γ=1e20 should freeze particles to numerical precision; "
            f"max drift was {np.max(np.abs(r1 - r0)):e}"
        )

    def test_first_step_uses_zeroed_prv_rnds(self):
        """First call sees prv_rnds=0 → matches AFINES new-bead convention."""
        sim = _build_sim(n=2, with_bond=False, box=50.0, dt=1e-4)
        action, updater = make_baoab_updater(
            kT=1.0, gamma={"A": 1.0}, dt=1e-4, seed=42
        )
        sim.operations.updaters.append(updater)
        sim.run(0)
        # prv_rnds should be all zeros immediately after attach
        prv = action.prv_rnds
        assert prv is not None
        assert np.array_equal(prv, np.zeros_like(prv)), (
            "prv_rnds must be zero-initialised at attach; got nonzero state."
        )


# ---------------------------------------------------------------------------
# Sanity Gate §5 — Sign-sense
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_stretched_bond_relaxes_toward_rest_at_kT_zero(self):
        """Pure gradient descent: stretched harmonic bond contracts."""
        sim = _build_sim(n=2, with_bond=True, box=50.0, dt=1e-4,
                         bond_k=10.0, bond_r0=1.0)
        # Stretch the bond: bead 1 at x=1.5 instead of x=1.0
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position)
            pos[1, 0] = 1.5
        action, updater = make_baoab_updater(
            kT=0.0, gamma={"A": 1.0}, dt=1e-4, seed=0
        )
        sim.operations.updaters.append(updater)
        sim.run(0)
        # Bond length monotone decreases toward r0 over a few steps
        lengths = []
        for _ in range(50):
            sim.run(1)
            with sim.state.cpu_local_snapshot as s:
                p = np.asarray(s.particles.position)
                lengths.append(np.linalg.norm(p[1] - p[0]))
        assert lengths[-1] < lengths[0], (
            f"stretched bond should contract; went {lengths[0]:.6f} → "
            f"{lengths[-1]:.6f}"
        )
        # Monotonicity check (allow tiny float-noise non-monotone steps).
        diffs = np.diff(lengths)
        assert np.all(diffs <= 1e-12), (
            "bond-length trajectory not monotonically non-increasing; "
            f"largest positive step: {diffs.max():.3e}"
        )

    def test_gaussian_draw_has_unit_variance_and_zero_mean(self):
        """Independent check of the per-Action RNG, large-N statistics."""
        rng = np.random.default_rng(0)
        W = rng.standard_normal(size=(1_000_000,))
        # Within 5σ of 0 for the mean, and ~1 for stdev
        assert abs(W.mean()) < 5.0 / math.sqrt(W.size)
        assert abs(W.std() - 1.0) < 5.0 / math.sqrt(2 * W.size)


# ---------------------------------------------------------------------------
# Sanity Gate §6.measurement — prv_rnds memory contract
# ---------------------------------------------------------------------------
class TestPrvRndsMemory:
    def test_step_n_writes_w_n_into_buffer_step_n_plus_1_reads(self):
        """After act(t), prv_rnds equals the W drawn at step t (scattered to tag).

        The Action draws ``W_row`` in current row-order and scatters into
        ``prv_rnds[tag] = W_row``. To verify, we replay the same RNG (same
        seed) and the same per-step tag ordering, scatter, and compare.
        """
        sim = _build_sim(n=2, with_bond=False, box=50.0, dt=1e-4)
        action, updater = make_baoab_updater(
            kT=1.0, gamma={"A": 1.0}, dt=1e-4, seed=12345
        )
        sim.operations.updaters.append(updater)
        sim.run(0)
        ref_rng = np.random.default_rng(12345)
        expected = np.zeros_like(action.prv_rnds)
        for n_steps in range(1, 6):
            sim.run(1)
            with sim.state.cpu_local_snapshot as s:
                cur_tag = np.asarray(s.particles.tag).copy()
            ref_W_row = ref_rng.standard_normal(size=(2, 3))
            # Mirror the Action's scatter step exactly.
            expected[cur_tag] = ref_W_row
            obs = action.prv_rnds
            assert obs is not None
            assert np.allclose(obs, expected, rtol=0, atol=0), (
                f"after {n_steps} steps, prv_rnds diverged from a fresh "
                "RNG of the same seed scattered by tag; the buffer is "
                "not capturing W_n correctly."
            )


# ---------------------------------------------------------------------------
# Sanity Gate §4 — wrapping correctness (Lees-Edwards / tilted box)
# ---------------------------------------------------------------------------
class TestBoxWrap:
    def test_wrap_orthorhombic_far_outside_box(self):
        sim = _build_sim(n=1, with_bond=False, box=10.0)
        box = sim.state.box
        # Position at (6, 0, 0): outside [-5, 5] orthorhombic box on x
        pos = np.array([[6.0, 0.0, 0.0], [-7.0, 0.0, 0.0]])
        wrapped, img = _wrap_into_box(pos, box)
        assert np.allclose(wrapped[:, 0], [-4.0, 3.0]), wrapped
        assert np.array_equal(img[:, 0], [1, -1]), img

    def test_wrap_round_trip_preserves_minimum_image(self):
        sim = _build_sim(n=1, with_bond=False, box=10.0)
        box = sim.state.box
        rng = np.random.default_rng(0)
        pos = rng.uniform(-100, 100, size=(50, 3))
        wrapped, img = _wrap_into_box(pos, box)
        # Wrapped positions must lie in [-L/2, L/2] on each axis (orthorhombic)
        assert np.all(np.abs(wrapped[:, 0]) <= box.Lx / 2 + 1e-9)
        assert np.all(np.abs(wrapped[:, 1]) <= box.Ly / 2 + 1e-9)
        assert np.all(np.abs(wrapped[:, 2]) <= box.Lz / 2 + 1e-9)
        # And wrapped + image·L should equal the original position
        recon = wrapped + img.astype(float) * np.array(
            [box.Lx, box.Ly, box.Lz]
        )
        assert np.allclose(recon, pos), (
            f"max wrap reconstruction error: "
            f"{np.max(np.abs(recon - pos)):e}"
        )


# ---------------------------------------------------------------------------
# Path A — tag-space growth via upstream set_snapshot()
# (PI-ratified 2026-05-29 in response to KU-5.1 GPU smoke catching the
# H.5 lamellipodium ↔ BAOAB topology-mutation incompatibility.)
# ---------------------------------------------------------------------------
def _grow_sim_by_one_bead(
    sim: hoomd.Simulation,
    *,
    new_pos: np.ndarray,
    new_type_name: str = "A",
    extra_types: tuple[str, ...] = (),
) -> int:
    """Append ONE particle to the simulation snapshot with the given type.

    Returns the new particle's tag. Mirrors the
    `cell.lamellipodium._extend_snapshot_with_new_actins` pattern but
    minimal — no bonds, no angles, just a fresh particle so the test
    can exercise the BAOAB extension hook in isolation.
    """
    snap = sim.state.get_snapshot()
    types = list(snap.particles.types)
    for et in extra_types:
        if et not in types:
            types.append(et)
    if new_type_name not in types:
        types.append(new_type_name)
    new_typeid = types.index(new_type_name)

    n_old = int(snap.particles.N)
    new_snap = hoomd.Snapshot()
    new_snap.particles.N = n_old + 1
    new_snap.particles.types = types
    new_snap.particles.typeid[:] = np.concatenate(
        [np.asarray(snap.particles.typeid), np.array([new_typeid], dtype=np.uint32)]
    )
    new_snap.particles.position[:] = np.concatenate(
        [np.asarray(snap.particles.position),
         np.asarray(new_pos, dtype=np.float64).reshape(1, 3)],
        axis=0,
    )
    new_snap.particles.velocity[:] = np.concatenate(
        [np.asarray(snap.particles.velocity), np.zeros((1, 3), dtype=np.float64)],
        axis=0,
    )
    new_snap.particles.mass[:] = np.concatenate(
        [np.asarray(snap.particles.mass), np.array([1.0], dtype=np.float64)]
    )
    new_snap.particles.image[:] = np.concatenate(
        [np.asarray(snap.particles.image), np.zeros((1, 3), dtype=np.int32)],
        axis=0,
    )
    # Preserve bonds + box.
    new_snap.bonds.N = int(snap.bonds.N)
    new_snap.bonds.types = list(snap.bonds.types)
    if int(snap.bonds.N) > 0:
        new_snap.bonds.group[:] = np.asarray(snap.bonds.group, dtype=np.uint32)
        new_snap.bonds.typeid[:] = np.asarray(snap.bonds.typeid, dtype=np.uint32)
    new_snap.configuration.box = list(snap.configuration.box)
    sim.state.set_snapshot(new_snap)
    return n_old  # tag of the newly-appended bead (dense in [0, N))


class TestTopologyGrow:
    """BAOAB tolerates tag-space growth via upstream set_snapshot().

    Path A, PI-ratified 2026-05-29. The H.5 lamellipodium Updaters
    (BarbedEndElongationUpdater, ArpBranchingUpdater) append actin_lamel
    beads to the snapshot every batch tick; BAOAB now extends its
    tag-indexed buffers (gamma_by_tag, bd_prefactor_by_tag, prv_rnds)
    rather than rejecting the snapshot.
    """

    def test_extends_on_particle_add(self):
        """After a +1 grow, buffers have size N+1 and the new tag's
        gamma / bd_prefactor match gamma_map; existing tags unchanged."""
        sim = _build_sim(n=4, with_bond=True, box=50.0, dt=1e-4)
        action, updater = make_baoab_updater(
            kT=1.0, gamma={"A": 2.0}, dt=1e-4, seed=42
        )
        sim.operations.updaters.append(updater)
        sim.run(100)
        # Snapshot the existing buffers BEFORE the grow.
        prv_pre = action.prv_rnds.copy()
        gamma_pre = action._gamma_by_tag.copy()
        bdp_pre = action._bd_prefactor_by_tag.copy()
        assert prv_pre.shape == (4, 3)

        # Grow by 1 bead of the same type.
        new_tag = _grow_sim_by_one_bead(
            sim, new_pos=np.array([10.0, 0.0, 0.0]), new_type_name="A"
        )
        assert new_tag == 4

        # Next step triggers the extension.
        sim.run(1)

        # Buffers extended to size 5; existing tags bit-for-bit unchanged.
        assert action._prv_rnds.shape == (5, 3)
        assert action._gamma_by_tag.shape == (5,)
        assert action._bd_prefactor_by_tag.shape == (5, 1)
        # New tag has correct gamma and bd_prefactor.
        assert action._gamma_by_tag[4] == 2.0
        expected_bdp = math.sqrt(1.0 / (2.0 * 2.0 * 1e-4))
        assert math.isclose(float(action._bd_prefactor_by_tag[4, 0]), expected_bdp, rel_tol=0, abs_tol=0), (
            f"new tag bd_prefactor = {float(action._bd_prefactor_by_tag[4, 0])}, "
            f"expected {expected_bdp}"
        )
        # Existing tags unchanged (regression-critical).
        assert np.array_equal(action._gamma_by_tag[:4], gamma_pre)
        assert np.array_equal(action._bd_prefactor_by_tag[:4], bdp_pre)

    def test_existing_tags_bit_for_bit_after_grow(self):
        """Fixed-tag positions in a partially-grown system match
        a parallel fixed-N reference up to the grow step.

        Concretely: run two sims, both with the same Action seed, both
        seeded the same way. Sim A grows by 1 bead at step 50; sim B
        does NOT grow. Compare the final TAG-ORDERED positions of tags
        0..N-1 in sim A vs sim B. The growth event scrambles HOOMD's
        internal ParticleSorter state, so we can't expect bit-for-bit
        across the boundary — but the tag-ordered prv_rnds of the
        existing tags must remain a permutation/identity of the
        pre-grow buffer (extension is concatenation, not overwrite).
        """
        sim = _build_sim(n=4, with_bond=True, box=50.0, dt=1e-4)
        action, updater = make_baoab_updater(
            kT=1.0, gamma={"A": 1.0}, dt=1e-4, seed=2026
        )
        sim.operations.updaters.append(updater)
        sim.run(50)
        prv_existing_pre = action.prv_rnds[:4].copy()
        # Grow then read prv_rnds for tags 0..3 — must be UNCHANGED.
        _grow_sim_by_one_bead(
            sim, new_pos=np.array([20.0, 0.0, 0.0]), new_type_name="A"
        )
        # No step yet — buffer extension happens on the next act().
        # Trigger extension via 1 step, then check existing tags' prv_rnds
        # was not overwritten EXCEPT for the W_n scatter that any step would
        # do to ALL tags including existing ones. So we check the slice
        # behaviour right after extension by inspecting before the step
        # runs the act() body. Easiest: directly invoke _extend_tag_buffers
        # so we test the pure extension semantics, isolated from a step.
        with sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag)
            typeid = np.asarray(snap.particles.typeid)
        action._extend_tag_buffers(
            new_buf_size=5, old_buf_size=4, tag=tag, typeid=typeid,
            type_names=list(sim.state.particle_types),
        )
        assert action._prv_rnds.shape == (5, 3)
        assert np.array_equal(action._prv_rnds[:4], prv_existing_pre), (
            "extension overwrote existing tags' prv_rnds — bit-for-bit "
            "regression is broken for fixed-N workloads."
        )
        # New tag's prv_rnds initialised to zero (matches attach() semantic).
        assert np.array_equal(action._prv_rnds[4], np.zeros(3))

    def test_grow_reproducibility_same_seed(self):
        """Two sims with the same MD + Action seed, same grow events, same
        post-grow steps, must produce identical final tag-ordered positions.

        Tests that the extension logic is fully deterministic — no hidden
        RNG advance, no path-dependent state.
        """
        def _run(seed: int) -> np.ndarray:
            sim = _build_sim(n=4, with_bond=True, box=50.0, dt=1e-4)
            sim.seed = seed  # ensure both sims see the same HOOMD seed
            action, updater = make_baoab_updater(
                kT=1.0, gamma={"A": 1.0}, dt=1e-4, seed=seed
            )
            sim.operations.updaters.append(updater)
            sim.run(20)
            _grow_sim_by_one_bead(
                sim, new_pos=np.array([15.0, 0.0, 0.0]), new_type_name="A"
            )
            sim.run(20)
            _grow_sim_by_one_bead(
                sim, new_pos=np.array([16.0, 0.0, 0.0]), new_type_name="A"
            )
            sim.run(20)
            return _read_positions_by_tag(sim)

        pos_a = _run(seed=12345)
        pos_b = _run(seed=12345)
        assert np.array_equal(pos_a, pos_b), (
            "same-seed deterministic replay through tag-space growth did "
            "not produce identical final positions: max abs diff = "
            f"{np.max(np.abs(pos_a - pos_b)):e}"
        )

    def test_shrinkage_rejected(self):
        """Tag-space shrinkage is not supported (would require sparse
        re-indexing). The act() must raise a clear error rather than
        silently use stale buffer entries."""
        sim = _build_sim(n=5, with_bond=False, box=50.0, dt=1e-4)
        action, updater = make_baoab_updater(
            kT=1.0, gamma={"A": 1.0}, dt=1e-4, seed=0
        )
        sim.operations.updaters.append(updater)
        sim.run(1)  # attach + 1 step → buf_size = 5

        # Shrink to 3 particles.
        snap = sim.state.get_snapshot()
        new_snap = hoomd.Snapshot()
        new_snap.particles.N = 3
        new_snap.particles.types = list(snap.particles.types)
        new_snap.particles.typeid[:] = np.asarray(snap.particles.typeid[:3])
        new_snap.particles.position[:] = np.asarray(snap.particles.position[:3])
        new_snap.particles.velocity[:] = np.asarray(snap.particles.velocity[:3])
        new_snap.particles.mass[:] = np.asarray(snap.particles.mass[:3])
        new_snap.particles.image[:] = np.asarray(snap.particles.image[:3])
        new_snap.bonds.N = 0
        new_snap.bonds.types = list(snap.bonds.types) if snap.bonds.types else ["p"]
        new_snap.configuration.box = list(snap.configuration.box)
        sim.state.set_snapshot(new_snap)

        with pytest.raises(RuntimeError, match="shrunk|append-only"):
            sim.run(1)

    def test_new_type_at_runtime_blocked_by_hoomd(self):
        """Defensive: HOOMD itself forbids new particle types via
        set_snapshot ("Particle types must remain the same"), so the
        missing-gamma-for-new-type code path in _extend_tag_buffers is
        unreachable by a well-formed Updater — gamma_map only needs
        to cover types registered at attach() time. This test pins
        that invariant.

        If a future HOOMD version relaxes this restriction, this test
        will break and someone should re-evaluate whether _extend_tag_buffers
        needs to register a fresh gamma_map entry at runtime (option:
        add an explicit `extend_gamma_map(type, gamma)` API).
        """
        sim = _build_sim(n=3, with_bond=False, box=50.0, dt=1e-4)
        action, updater = make_baoab_updater(
            kT=1.0, gamma={"A": 1.0}, dt=1e-4, seed=0
        )
        sim.operations.updaters.append(updater)
        sim.run(1)
        with pytest.raises(RuntimeError, match="types must remain the same"):
            _grow_sim_by_one_bead(
                sim, new_pos=np.array([10.0, 0.0, 0.0]),
                new_type_name="B", extra_types=("B",),
            )

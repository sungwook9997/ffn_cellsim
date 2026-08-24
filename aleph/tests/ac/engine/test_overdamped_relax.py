"""The component integrator — its sign, its derived bound, and what it refuses.

The dangerous failures here are silent ones: a sign error grows every mode, and a `dt` above the
stability bound produces a trajectory that looks like physics for a while. Both are pinned.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.engine.overdamped_relax import (
    OverdampedStep,
    build_overdamped_step,
    derive_gershgorin_lambda_max,
    derive_point_drag,
)
from aleph.laws import units as U
from aleph.laws.hand_kmc import ALPHA_ACTININ

wp.init()
_CUDA_DEVICE = next((str(d) for d in wp.get_devices() if d.is_cuda), None)


class _Topology:
    """The shape `SFMechanicsTopology` already has, reduced to what the builder reads."""

    def __init__(self, n_nodes, links, link_k, link_r0, bend_triples=None, bend_alpha=None):
        self.n_nodes = n_nodes
        self.links = np.asarray(links, dtype=np.int64)
        self.link_k = np.asarray(link_k, dtype=np.float64)
        self.link_r0 = np.asarray(link_r0, dtype=np.float64)
        self.bend_triples = None if bend_triples is None else np.asarray(bend_triples, dtype=np.int64)
        self.bend_alpha = None if bend_alpha is None else np.asarray(bend_alpha, dtype=np.float64)


def _null(*args, **kwargs) -> None:
    """A recorder-free launcher: these gates are host arithmetic, and CPU execution is forbidden."""


def _chain(n_nodes=5, k=100.0, seg=0.1, alpha=None):
    links = [[i, i + 1] for i in range(n_nodes - 1)]
    triples = [[i, i + 1, i + 2] for i in range(n_nodes - 2)] if alpha is not None else None
    return _Topology(
        n_nodes, links, [k] * len(links), [seg] * len(links),
        triples, None if alpha is None else [alpha] * len(triples),
    )


def test_gamma_is_the_sourced_nf2007_law_not_a_local_formula() -> None:
    """It must agree with `units.fiber_point_drag` exactly — same law, read off a topology."""
    seg, n_seg = 0.1, 4
    gamma = derive_point_drag(
        segment_rest_um=np.full(n_seg, seg), n_nodes=n_seg + 1, n_filaments=1,
    )
    assert gamma == pytest.approx(U.fiber_point_drag(seg * n_seg, n_seg))


def test_an_empty_topology_refuses_rather_than_defaulting() -> None:
    """A component with no segments has no fiber mobility; a default would be a chosen constant."""
    with pytest.raises(ValueError, match="chosen constant"):
        derive_point_drag(segment_rest_um=np.array([]), n_nodes=0, n_filaments=0)


def test_gershgorin_counts_two_k_per_link_endpoint() -> None:
    """A Hookean link puts `2k` in each endpoint's row: its diagonal plus the off-diagonal."""
    lam = derive_gershgorin_lambda_max(
        links=np.array([[0, 1]]), link_k=np.array([100.0]), n_nodes=2,
    )
    assert lam == pytest.approx(200.0)

    interior = derive_gershgorin_lambda_max(
        links=np.array([[0, 1], [1, 2]]), link_k=np.array([100.0, 100.0]), n_nodes=3,
    )
    assert interior == pytest.approx(400.0), "the shared node carries both links"


def test_bending_contributes_sixteen_alpha() -> None:
    """λ_max = 16κ/seg³ and α = κ/seg³ — the identity laws/relax.py's CFL already uses."""
    lam = derive_gershgorin_lambda_max(
        links=np.zeros((0, 2), dtype=np.int64), link_k=np.array([]), n_nodes=3,
        bend_triples=np.array([[0, 1, 2]]), bend_alpha=np.array([2.0]),
    )
    assert lam == pytest.approx(32.0)


def test_the_arc_crosslink_stiffness_is_IN_the_bound() -> None:
    """α-actinin is ~460× the axial k. Omitting it does not loosen the bound — it deletes it.

    Measured when this was found: a node carrying one axial pair and one `arc_joint` bounded at
    4,000 pN/µm without the joint and 924,000 with it, so the admitted `dt` was **231× too large**.
    Nothing raises; the trajectory just looks like physics for a while.
    """
    links, link_k = np.array([[0, 1], [1, 2]]), np.array([1000.0, 1000.0])
    joints, arc_k = np.array([[1, 3]]), np.array([ALPHA_ACTININ.link_k])

    without = derive_gershgorin_lambda_max(links=links, link_k=link_k, n_nodes=4)
    with_arc = derive_gershgorin_lambda_max(
        links=links, link_k=link_k, n_nodes=4, extra_pairs=((joints, arc_k),),
    )
    assert without == pytest.approx(4000.0)
    assert with_arc == pytest.approx(4000.0 + 2.0 * ALPHA_ACTININ.link_k)
    assert with_arc / without > 100.0, "the omission is orders of magnitude, not a rounding"


def test_the_builder_reads_arc_joints_off_a_real_topology_shape() -> None:
    """The builder must find them itself — passing them by hand is what was already missed once."""
    plain = _chain()
    with_joints = _chain()
    with_joints.arc_joints = np.array([[1, 3]], dtype=np.int64)
    with_joints.arc_k = np.array([ALPHA_ACTININ.link_k])

    bare = build_overdamped_step(plain, n_filaments=1, borrowed_pairs=(), launch=_null)
    joined = build_overdamped_step(with_joints, n_filaments=1, borrowed_pairs=(), launch=_null)
    assert joined.lambda_max_pn_per_um > bare.lambda_max_pn_per_um
    assert joined.dt_max_s < bare.dt_max_s, "a stiffer graph must admit a SMALLER step"


def test_a_stiffness_array_the_builder_cannot_sum_is_REFUSED() -> None:
    """The next `arc_k` must break the build, not the run.

    `arc_k` was omitted silently because nothing required the builder to account for every stiffness
    it was handed. A named allowlist alone would repeat that; this refuses the unrecognised one.
    """
    topology = _chain()
    topology.membrane_k = np.array([500.0])
    with pytest.raises(ValueError, match="does not sum into its Gershgorin"):
        build_overdamped_step(topology, n_filaments=1, borrowed_pairs=(), launch=_null)


def test_an_empty_extra_stiffness_array_is_not_refused() -> None:
    """A component that declares the field but has no such elements is not a defect."""
    topology = _chain()
    topology.arc_joints = np.zeros((0, 2), dtype=np.int64)
    topology.arc_k = np.array([])
    step = build_overdamped_step(topology, n_filaments=1, borrowed_pairs=(), launch=_null)
    assert step.lambda_max_pn_per_um == pytest.approx(400.0), "chain of k=100, interior node"


def test_a_connector_borrowing_the_arrays_raises_the_bound() -> None:
    """`sf_cortex_transient` scatters α-actinin into `sf_owner.force_d` while owning no array here.

    The topology scan cannot see it — that stiffness is on the CONNECTOR, not the component — so
    `borrowed_pairs` is the only channel, and it is why that argument has no default.
    """
    joints, k = np.array([[1, 3]]), np.array([ALPHA_ACTININ.link_k])
    alone = build_overdamped_step(_chain(), n_filaments=1, borrowed_pairs=(), launch=_null)
    borrowed = build_overdamped_step(
        _chain(), n_filaments=1, borrowed_pairs=((joints, k),), launch=_null,
    )
    assert borrowed.lambda_max_pn_per_um == pytest.approx(
        alone.lambda_max_pn_per_um + 2.0 * ALPHA_ACTININ.link_k
    )
    assert borrowed.dt_max_s < alone.dt_max_s / 100.0


def test_borrowed_pairs_has_no_default_so_a_borrower_cannot_be_forgotten() -> None:
    """A default of `()` is what let `arc_k` sit outside the bound and admit a dt 459× too large."""
    with pytest.raises(TypeError, match="borrowed_pairs"):
        build_overdamped_step(_chain(), n_filaments=1, launch=_null)  # type: ignore[call-arg]


def test_the_row_sum_bounds_the_true_nonlinear_tangent() -> None:
    """The system is NOT linear: a stretched bond carries `k_t = k(L−r0)/L` transverse to itself.

    The scalar row sum omits that term, so "Gershgorin bounds λ_max" is a CLAIM about the real 3N×3N
    tangent, not a theorem we inherit. Assemble the exact tangent and check the ratio never reaches 1
    — including the stiff, high-coordination, heavily stretched shape `arc_joints` actually creates.
    """
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(200):
        n = int(rng.integers(3, 10))
        pos = rng.normal(0.0, 1.0, (n, 3))
        bonds = [(i, i + 1) for i in range(n - 1)]
        for _ in range(int(rng.integers(0, 2 * n))):          # extra pairs: arc-joint coordination
            i, j = rng.choice(n, 2, replace=False)
            bonds.append((int(i), int(j)))
        k = rng.uniform(1.0, ALPHA_ACTININ.link_k, len(bonds))
        stretch = float(rng.uniform(1.0, 50.0))
        r0 = [float(np.linalg.norm(pos[i] - pos[j])) / stretch for i, j in bonds]

        tangent = np.zeros((3 * n, 3 * n))
        for (i, j), k_b, r0_b in zip(bonds, k, r0):
            d = pos[i] - pos[j]
            length = float(np.linalg.norm(d))
            unit = d / length
            k_t = k_b * (length - r0_b) / length
            block = k_b * np.outer(unit, unit) + k_t * (np.eye(3) - np.outer(unit, unit))
            for a, sa in ((i, 1.0), (j, -1.0)):
                for b, sb in ((i, 1.0), (j, -1.0)):
                    tangent[3 * a:3 * a + 3, 3 * b:3 * b + 3] += sa * sb * block

        bound = derive_gershgorin_lambda_max(
            links=np.zeros((0, 2), dtype=np.int64), link_k=np.array([]), n_nodes=n,
            extra_pairs=((np.asarray(bonds), k),),
        )
        worst = max(worst, float(np.linalg.eigvalsh(tangent).max()) / bound)

    assert worst < 1.0, f"the row sum is NOT a bound on the nonlinear tangent (ratio {worst:.4f})"
    assert worst > 0.5, "a bound this loose would cost subcycles for nothing; it should stay tight"


def test_the_stability_bound_is_exact_with_no_safety_factor() -> None:
    """`2γ/λ_max`. A safety factor inside a stability statement reads as physics and is not."""
    step = OverdampedStep(gamma_pn_s_per_um=4.0, lambda_max_pn_per_um=8.0, n_nodes=2)
    assert step.dt_max_s == pytest.approx(1.0)

    free = OverdampedStep(gamma_pn_s_per_um=4.0, lambda_max_pn_per_um=0.0, n_nodes=2)
    assert free.dt_max_s == float("inf"), "no stiffness means no stability limit"


def test_an_unstable_dt_is_refused_and_the_message_carries_the_bound() -> None:
    """Integrating past the bound produces a trajectory that looks like physics for a while."""
    step = OverdampedStep(4.0, 8.0, 2, launch=_null)
    with pytest.raises(ValueError, match="Subcycle 2 times"):
        step.advance(_fake(2), _fake(2), dt_phys=1.5)


def test_zero_dt_moves_nothing_and_negative_dt_is_refused() -> None:
    launched = []
    step = OverdampedStep(4.0, 8.0, 2, launch=lambda *a, **k: launched.append(a))
    step.advance(_fake(2), _fake(2), dt_phys=0.0)
    assert launched == []
    with pytest.raises(ValueError, match="nonnegative"):
        step.advance(_fake(2), _fake(2), dt_phys=-0.1)


def test_a_mismatched_array_is_refused() -> None:
    step = OverdampedStep(4.0, 8.0, 5, launch=_null)
    with pytest.raises(ValueError, match="derived for 5 nodes"):
        step.advance(_fake(2), _fake(2), dt_phys=0.1)


class _FakeArray:
    def __init__(self, n):
        self.shape = (n,)
        self.device = "none"


def _fake(n):
    return _FakeArray(n)


def test_the_builder_derives_both_numbers_from_the_topology_alone() -> None:
    step = build_overdamped_step(_chain(alpha=0.5), n_filaments=1, borrowed_pairs=(), launch=_null)
    assert step.gamma_pn_s_per_um > 0.0
    assert step.lambda_max_pn_per_um > 0.0
    assert step.dt_max_s > 0.0
    assert step.n_nodes == 5


# ── device: the sign, and that a stable step actually descends ────────────────────────────────────

@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: the integrator kernel requires a CUDA GPU")
def test_a_node_moves_ALONG_its_force() -> None:
    """The opposite sign grows every mode, and a trajectory would not obviously show it."""
    pos = wp.array(np.zeros((2, 3)), dtype=wp.vec3d, device=_CUDA_DEVICE)
    force = wp.array(np.array([[1.0, 0.0, 0.0], [0.0, -2.0, 0.0]]), dtype=wp.vec3d,
                     device=_CUDA_DEVICE)
    step = OverdampedStep(gamma_pn_s_per_um=2.0, lambda_max_pn_per_um=0.0, n_nodes=2)
    step.advance(pos, force, dt_phys=0.5)
    wp.synchronize_device(wp.get_device(_CUDA_DEVICE))
    moved = pos.numpy()
    assert moved[0][0] == pytest.approx(0.25), "dt/gamma = 0.25, force +1 -> +0.25"
    assert moved[1][1] == pytest.approx(-0.5)


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: the integrator kernel requires a CUDA GPU")
def test_a_stretched_hookean_pair_relaxes_toward_its_rest_length() -> None:
    """End to end at a stable dt: the step must REDUCE the extension, not merely change it."""
    k, rest = 100.0, 1.0
    step = build_overdamped_step(
        _Topology(2, [[0, 1]], [k], [rest]), n_filaments=1, borrowed_pairs=(),
    )
    pos_np = np.array([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]])
    pos = wp.array(pos_np, dtype=wp.vec3d, device=_CUDA_DEVICE)

    extension_before = float(np.linalg.norm(pos_np[1] - pos_np[0]) - rest)
    for _ in range(20):
        current = pos.numpy()
        delta = current[1] - current[0]
        length = float(np.linalg.norm(delta))
        pair = k * (length - rest) * delta / length
        force = wp.array(np.array([pair, -pair]), dtype=wp.vec3d, device=_CUDA_DEVICE)
        step.advance(pos, force, dt_phys=step.dt_max_s * 0.25)
    wp.synchronize_device(wp.get_device(_CUDA_DEVICE))

    final = pos.numpy()
    extension_after = float(np.linalg.norm(final[1] - final[0]) - rest)
    assert abs(extension_after) < abs(extension_before), "a stable step must relax the spring"
    assert abs(extension_after) < 0.5 * abs(extension_before)

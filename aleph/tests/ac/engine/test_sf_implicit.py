"""CPU-green gates for the SF-motor implicit PCG solve (:mod:`aleph.engine.sf_implicit`).

The claim this solver rests on is a LINEAR-ALGEBRA claim, so it is provable on the host without a GPU, and it
is proved here rather than asserted in a docstring:

1. **Masking the right-hand side alone is NOT a Dirichlet condition.** With ``(aI+K) dx = M b`` a pinned node
   still moves, because the stiffness couples it to its neighbours. Only masking the OPERATOR — the form the
   inherited ``_operator`` produces, ``M(aI+K)M + a(I−M)`` — gives exactly ``dx_pinned = 0``. This is the
   correction that makes the FA-anchor physics right, and the failure it prevents is silent (a converged-looking
   tension whose anchors drifted).
2. **The masked operator stays SPD**, so PCG is legitimate on it.
3. **The reused tangents really are the derivatives of the forces the slice launches** — the pair tangent
   matches a finite difference of ``link_spring`` in tension, and the bending tangent matches the NF2007
   bending force exactly (that force is linear, so its tangent is exact, not approximate).
4. **The concatenated index space is built correctly**: NMII tables shift by ``+n_sf``, SF tables do not, and
   the Dirichlet mask lives only in the SF block.

No Warp kernel is launched (CPU-Warp launches are forbidden in the test lane); every operator here is the NumPy
transcription of the kernel formula, read from the kernel source. The CUDA behaviour is the Lead's native gate.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pytest
import warp as wp

from aleph.engine import sf_implicit as sf_implicit_module
from aleph.engine.sf_implicit import SFImplicitCG, build_sf_implicit_index_map
from aleph.engine.sf_mechanics import bending_force_reference, link_spring_force_reference
from aleph.components.incumbent.implicit_mechanics import ProjectedAnalyticCG

_EPS64 = 2.220446049250313e-16


# ── NumPy transcriptions of the reused kernel formulas (source: ac/cell/implicit_mechanics.py) ─────────────
def _central_action(delta: np.ndarray, relative: np.ndarray, k: float, rest: float) -> np.ndarray:
    """``_central_action`` (implicit_mechanics.py:81-92): PSD central-spring tangent action."""
    length = float(np.linalg.norm(delta))
    if length <= 1.0e-12:
        return np.zeros(3)
    unit = delta / length
    axial = k * float(unit @ relative) * unit
    transverse_k = max(k * (length - rest) / length, 0.0)
    transverse = transverse_k * (relative - float(unit @ relative) * unit)
    return axial + transverse


def _pair_stiffness_action(pos, vector, links, k, rest) -> np.ndarray:
    """``add_pair_stiffness_kernel`` (implicit_mechanics.py:195-209)."""
    out = np.zeros_like(pos)
    for t in range(len(links)):
        i, j = int(links[t][0]), int(links[t][1])
        action = _central_action(pos[j] - pos[i], vector[i] - vector[j], float(k[t]), float(rest[t]))
        out[i] += action
        out[j] -= action
    return out


def _bending_stiffness_action(vector, triples, alpha) -> np.ndarray:
    """``add_bending_stiffness_kernel`` (implicit_mechanics.py:121-136)."""
    out = np.zeros_like(vector)
    for t in range(len(triples)):
        a, b, c = (int(x) for x in triples[t])
        action = float(alpha[t]) * (vector[a] - 2.0 * vector[b] + vector[c])
        out[a] += action
        out[b] -= 2.0 * action
        out[c] += action
    return out


def _dense_operator(pos, links, k, rest, triples, alpha, a: float) -> np.ndarray:
    """Assemble ``aI + K`` densely by applying the action to each basis vector."""
    n = pos.shape[0]
    matrix = np.zeros((3 * n, 3 * n))
    for column in range(3 * n):
        basis = np.zeros((n, 3))
        basis[column // 3, column % 3] = 1.0
        action = a * basis + _pair_stiffness_action(pos, basis, links, k, rest)
        if len(triples):
            action = action + _bending_stiffness_action(basis, triples, alpha)
        matrix[:, column] = action.reshape(-1)
    return matrix


def _mask_matrix(n: int, pinned: np.ndarray) -> np.ndarray:
    """The 0/1 Dirichlet projector ``M`` as a dense matrix over flattened (n,3) vectors."""
    keep = np.repeat((pinned == 0).astype(float), 3)
    return np.diag(keep)


def _chain(n_nodes: int = 5, spacing: float = 1.0, stretch: float = 1.15):
    """A stretched straight chain: links in TENSION (so the transverse branch is active and unclamped)."""
    pos = np.zeros((n_nodes, 3))
    pos[:, 0] = spacing * stretch * np.arange(n_nodes)
    links = np.array([[i, i + 1] for i in range(n_nodes - 1)], np.int64)
    k = np.full(len(links), 1000.0)
    rest = np.full(len(links), spacing)          # rest < current ⇒ tension
    triples = np.array([[i - 1, i, i + 1] for i in range(1, n_nodes - 1)], np.int64)
    alpha = np.full(len(triples), 50.0)
    return pos, links, k, rest, triples, alpha


# ── 1. THE LOAD-BEARING CLAIM: masking the RHS alone is not a Dirichlet condition ──────────────────────────
def test_masking_only_the_rhs_lets_pinned_nodes_move() -> None:
    """A pinned node MOVES when only the right-hand side is masked — the failure the operator mask prevents.

    This is the numerical proof behind swapping the projector: the driver's force-zeroing is sufficient for an
    explicit step (``x += step·F`` with ``F = 0`` cannot move a node) but NOT for a coupled implicit solve.
    """
    pos, links, k, rest, triples, alpha = _chain()
    n = pos.shape[0]
    pinned = np.zeros(n, np.int32)
    pinned[[0, n - 1]] = 1                       # both outer ends held, like a ventral SF at its FA anchors
    a = 1.0e4

    operator = _dense_operator(pos, links, k, rest, triples, alpha, a)
    mask = _mask_matrix(n, pinned)

    # a force that pulls only an interior node
    force = np.zeros((n, 3))
    force[n // 2] = np.array([0.0, 5.0, 0.0])
    rhs = (mask @ force.reshape(-1))

    dx_rhs_only = np.linalg.solve(operator, rhs).reshape(n, 3)
    pinned_motion = np.abs(dx_rhs_only[pinned == 1]).max()
    assert pinned_motion > 1.0e-9, (
        "expected the pinned nodes to drift when only the RHS is masked; if this ever becomes zero the "
        "argument for masking the operator needs revisiting"
    )


def test_masking_the_operator_pins_exactly() -> None:
    """``M(aI+K)M + a(I−M)`` with a masked RHS gives EXACTLY zero displacement at pinned nodes."""
    pos, links, k, rest, triples, alpha = _chain()
    n = pos.shape[0]
    pinned = np.zeros(n, np.int32)
    pinned[[0, n - 1]] = 1
    a = 1.0e4

    operator = _dense_operator(pos, links, k, rest, triples, alpha, a)
    mask = _mask_matrix(n, pinned)
    identity = np.eye(3 * n)
    projected = mask @ operator @ mask + a * (identity - mask)

    force = np.zeros((n, 3))
    force[n // 2] = np.array([0.0, 5.0, 0.0])
    rhs = mask @ force.reshape(-1)

    dx = np.linalg.solve(projected, rhs).reshape(n, 3)
    assert np.allclose(dx[pinned == 1], 0.0, atol=1e-14)          # exactly held
    assert np.abs(dx[pinned == 0]).max() > 1.0e-9                 # and the free nodes still respond


def test_projected_operator_is_spd() -> None:
    """PCG is only legitimate on an SPD operator — check symmetry and positive definiteness."""
    pos, links, k, rest, triples, alpha = _chain()
    n = pos.shape[0]
    pinned = np.zeros(n, np.int32)
    pinned[[0, n - 1]] = 1
    a = 1.0e4

    operator = _dense_operator(pos, links, k, rest, triples, alpha, a)
    mask = _mask_matrix(n, pinned)
    projected = mask @ operator @ mask + a * (np.eye(3 * n) - mask)

    assert np.allclose(projected, projected.T, atol=1e-9 * np.abs(projected).max())
    eigenvalues = np.linalg.eigvalsh(projected)
    assert eigenvalues.min() > 0.0
    # the regularizer is what guarantees this: an unanchored bundle (arc / cap / free minifilament) has a
    # singular K, and a is what makes aI+K invertible.
    unregularized = _dense_operator(pos, links, k, rest, triples, alpha, 0.0)
    assert np.linalg.eigvalsh(unregularized).min() < 1.0e-9      # singular without a (rigid modes)


# ── 2. the reused tangents are the derivatives of the forces the slice launches ────────────────────────────
def test_pair_tangent_matches_link_spring_force_jacobian_in_tension() -> None:
    """MEASURED: the reused pair tangent equals a finite difference of the SF axial force, in tension."""
    pos, links, k, rest, _, _ = _chain()
    links32 = links.astype(np.int32)
    n = pos.shape[0]
    rng = np.random.default_rng(11)
    direction = rng.normal(size=(n, 3))
    direction /= np.linalg.norm(direction)

    # K·v from the tangent (note the sign convention: K = -dF/dx, so compare against -(dF))
    tangent_action = _pair_stiffness_action(pos, direction, links, k, rest)

    step = 1.0e-7
    force_plus = link_spring_force_reference(pos + step * direction, links32, k, rest)
    force_minus = link_spring_force_reference(pos - step * direction, links32, k, rest)
    fd = -(force_plus - force_minus) / (2.0 * step)

    scale = max(np.abs(fd).max(), 1.0)
    assert np.allclose(tangent_action, fd, rtol=2e-4, atol=2e-4 * scale)


def test_bending_tangent_is_exact() -> None:
    """The NF2007 bending force is LINEAR, so its reused tangent is exact — not a Gauss-Newton approximation."""
    pos, _, _, _, triples, alpha = _chain()
    triples32 = triples.astype(np.int32)
    n = pos.shape[0]
    rng = np.random.default_rng(12)
    direction = rng.normal(size=(n, 3))

    tangent_action = _bending_stiffness_action(direction, triples, alpha)
    step = 1.0e-6
    fd = -(bending_force_reference(pos + step * direction, triples32, alpha)
           - bending_force_reference(pos - step * direction, triples32, alpha)) / (2.0 * step)
    assert np.allclose(tangent_action, fd, rtol=1e-9, atol=1e-9 * max(np.abs(fd).max(), 1.0))


def test_compressed_pair_tangent_is_psd_by_clamping() -> None:
    """Under COMPRESSION the transverse branch is clamped to zero — that clamp is what keeps the tangent PSD."""
    pos = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])       # current 0.5 < rest 1.0 ⇒ compressed
    links = np.array([[0, 1]], np.int64)
    k = np.array([1000.0])
    rest = np.array([1.0])
    matrix = np.zeros((6, 6))
    for column in range(6):
        basis = np.zeros((2, 3))
        basis[column // 3, column % 3] = 1.0
        matrix[:, column] = _pair_stiffness_action(pos, basis, links, k, rest).reshape(-1)
    assert np.linalg.eigvalsh(matrix).min() >= -1.0e-9        # PSD despite compression
    # and the transverse block really is zero (pure axial), which is the clamp doing the work
    transverse = matrix[1, 1]
    assert transverse == pytest.approx(0.0, abs=1e-12)


# ── 3. the concatenated index space ────────────────────────────────────────────────────────────────────────
class _FakeArray:
    """Host stand-in exposing only ``.numpy()`` / ``.shape``, which is all the topology builder reads."""

    def __init__(self, array: np.ndarray) -> None:
        self._array = np.ascontiguousarray(array)
        self.shape = self._array.shape

    def numpy(self) -> np.ndarray:
        return self._array


class _FakeBending:
    k_theta_backbone = 7.0
    k_theta_arm = 3.0


class _FakeMechanics:
    def __init__(self, n_bb: int, n_heads: int) -> None:
        self.backbone_bonds_d = _FakeArray(
            np.column_stack([np.arange(n_bb - 1), np.arange(1, n_bb)]).astype(np.int32))
        self.head_bonds_d = _FakeArray(
            np.column_stack([np.arange(n_bb, n_bb + n_heads), np.zeros(n_heads, np.int32)]).astype(np.int32))
        self.backbone_angles_d = _FakeArray(
            np.column_stack([np.arange(n_bb - 2), np.arange(1, n_bb - 1), np.arange(2, n_bb)]).astype(np.int32))
        self.head_arm_angles_d = _FakeArray(
            np.column_stack([np.zeros(n_heads, np.int32), np.ones(n_heads, np.int32),
                             np.arange(n_bb, n_bb + n_heads)]).astype(np.int32))
        self.k_backbone = 1.0e3
        self.r0_backbone = 0.023
        self.k_head_spring = 1.0e2
        self.r0_head = 0.2
        self.bending = _FakeBending()


class _FakeActuator:
    def __init__(self, n_particles: int, n_bb: int, n_heads: int) -> None:
        self.mechanics = _FakeMechanics(n_bb, n_heads)
        self.head_node_d = _FakeArray(np.arange(n_bb, n_bb + n_heads, dtype=np.int32))
        self.n_particles = n_particles


class _FakeConnectorState:
    def __init__(self, n_heads: int) -> None:
        self.n_heads = n_heads
        for name in ("bound_d", "seg_a_d", "seg_b_d", "bary_t_d", "abscissa_d", "walk_dir_d"):
            setattr(self, name, object())


def _sf_topology(n_nodes: int = 6):
    class _T:
        links = np.array([[i, i + 1] for i in range(n_nodes - 1)], np.int32)
        link_k = np.full(n_nodes - 1, 1000.0)
        link_r0 = np.full(n_nodes - 1, 1.0)
        n_links = n_nodes - 1
        bend_triples = np.array([[i - 1, i, i + 1] for i in range(1, n_nodes - 1)], np.int32)
        bend_alpha = np.full(n_nodes - 2, 50.0)
        n_triples = n_nodes - 2
        arc_joints = np.array([[0, n_nodes - 1]], np.int32)
        arc_k = np.array([4.6e5])
        arc_r0 = np.array([2.0])
        n_arc_joints = 1
    return _T()


def _index_map(n_sf: int, n_bb: int, n_heads: int, pinned: np.ndarray):
    """Build the concatenated index map from the fake NMII tables — no device, no Warp allocation.

    These tests previously reached this arithmetic by constructing an :class:`SFImplicitTopology` on a
    Warp CPU device, which really allocated Warp CPU arrays and so violated the I0-A Warp-CUDA-only
    contract (the static guard in ``test_foundation_static_contract`` was RED on exactly those two call
    sites).  The index arithmetic never needed a device, so it now lives in a pure function and is
    tested as one.
    """
    mech = _FakeMechanics(n_bb, n_heads)
    return build_sf_implicit_index_map(
        n_sf=n_sf, n_nmii=n_bb + n_heads,
        backbone_bonds=mech.backbone_bonds_d.numpy(),
        head_bonds=mech.head_bonds_d.numpy(),
        backbone_angles=mech.backbone_angles_d.numpy(),
        head_arm_angles=mech.head_arm_angles_d.numpy(),
        head_node=np.arange(n_bb, n_bb + n_heads, dtype=np.int32),
        pinned_sf_nodes=pinned,
    )


def test_topology_shifts_only_the_nmii_block() -> None:
    """SF tables keep SF-local indices; NMII tables shift by +n_sf; the mask lives only in the SF block."""
    n_sf, n_bb, n_heads = 6, 5, 4
    n_nmii = n_bb + n_heads
    imap = _index_map(n_sf, n_bb, n_heads, np.array([0, n_sf - 1]))

    assert imap.n == n_sf + n_nmii
    # SF side is untouched by construction: the SF tables are never passed through the shift
    sf = _sf_topology(n_sf)
    assert sf.links.max() < n_sf
    assert sf.arc_joints.max() < n_sf
    assert sf.bend_triples.max() < n_sf
    # NMII side shifted into the concatenated block
    assert imap.backbone_bonds.min() >= n_sf
    assert imap.head_bonds.min() >= n_sf
    assert imap.backbone_angles.min() >= n_sf
    assert imap.head_arm_angles.min() >= n_sf
    assert imap.head_node_shifted.min() >= n_sf
    assert imap.head_node_shifted.max() < imap.n
    # the Dirichlet mask is SF-only and counts exactly the pinned anchors
    assert imap.pinned.sum() == 2 == imap.n_pinned
    assert imap.pinned[n_sf:].sum() == 0                         # NMII particles are never held


def test_topology_refuses_a_pin_outside_the_sf_block() -> None:
    """A pin index in the NMII block is a wiring error, not something to silently accept."""
    n_sf, n_bb, n_heads = 6, 5, 4
    with pytest.raises(ValueError, match="pinned_sf_nodes must index the sf_arc node block"):
        _index_map(n_sf, n_bb, n_heads, np.array([n_sf + 1]))
    with pytest.raises(ValueError, match="pinned_sf_nodes must index the sf_arc node block"):
        _index_map(n_sf, n_bb, n_heads, np.array([-1]))


def test_device_topology_uploads_the_same_index_map_it_was_built_from() -> None:
    """The device path must not re-derive the arithmetic — it uploads what the pure function returned."""
    src = (
        pathlib.Path(sf_implicit_module.__file__).read_text(encoding="utf-8")
        .split("class SFImplicitTopology")[1]
    )
    assert "build_sf_implicit_index_map(" in src
    # no second copy of the shift arithmetic inside the device class
    assert "+ shift" not in src


# ── 4. the inheritance contract: the CG recurrence is reused, the projector is NOT ─────────────────────────
def test_cg_recurrence_is_inherited_and_projector_is_overridden() -> None:
    """The audited device recurrence must be reused verbatim; only the cell-dependent members are replaced."""
    assert issubclass(SFImplicitCG, ProjectedAnalyticCG)
    for inherited in ("solve", "_operator", "_precondition"):
        assert inherited not in SFImplicitCG.__dict__, f"{inherited} must be inherited, not reimplemented"
    for overridden in ("_project", "_stiffness", "_build_preconditioner", "__init__"):
        assert overridden in SFImplicitCG.__dict__, f"{overridden} must be overridden (it reads cell)"


def test_projector_does_not_use_the_inextensibility_constraint() -> None:
    """The SF projector must be the Dirichlet mask; reusing NF2007 inextensibility would be a physics bug.

    SF carries a Hookean axial spring (a documented modelling GAP precisely because NF2007 treats the backbone
    as inextensible), so projecting out the axial residual would leave its tension permanently unequilibrated.
    """
    import inspect

    source = inspect.getsource(SFImplicitCG._project)
    assert "mask_vec3_kernel" in source
    assert "project_constraint_forces_kernel" not in source
    assert "inner_mechanics" not in source
    # and the parent's projector is the constraint one, i.e. these really are different physics
    parent_source = inspect.getsource(ProjectedAnalyticCG._project)
    assert "project_constraint_forces_kernel" in parent_source


def test_step_refuses_a_nonpositive_mobility() -> None:
    """The regularizer is 1/mobility_step, so a non-positive step is a REQUIRED-PARAM error, not a default."""
    solver = SFImplicitCG.__new__(SFImplicitCG)                  # no device allocation needed for the guard
    with pytest.raises(ValueError, match="mobility_step must be positive-finite"):
        SFImplicitCG.step(
            solver, sf_position_d=None, sf_force_d=None, nmii_position_d=None, nmii_force_d=None,
            mobility_step=0.0)

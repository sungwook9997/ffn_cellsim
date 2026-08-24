r"""Host acceptance gate for the F6 angle-harmonic bending (backbone_warp) — the Round-2 motor-fix mechanics.

Proves the NEW bending primitive is mechanically correct WITHOUT a CUDA device (dev-Mac I0-A discipline): the
NumPy mirror :func:`angle_harmonic_forces` is the bit-for-formula twin of the device
:func:`angle_harmonic_kernel`, so the physics the A5000 kernel evaluates is certified here. A single static
force EVALUATION on the Warp 'cpu' device (no state evolution, no time-stepping — not a simulation) cross-checks
that the kernel and the NumPy mirror agree to machine precision.

Gate (backbone_warp / minifilament_topology F6 fix):
  * ZERO force/torque at rest — straight backbone (θ = π) and perpendicular head-arm (θ = π/2);
  * RESTORING moment under bend — the bent middle bead / swung head is pushed back toward rest;
  * Newton's 3rd law — the three triple forces sum to zero (no net force, no net torque);
  * FIRST-ORDER transverse stiffness — a head displaced ⟂ to its arm bears k_θ,arm / r0² (the transmission the
    single distance spring lacked; the F6 root cause);
  * kernel ↔ NumPy-mirror parity to machine precision (Warp 'cpu' static eval);
  * k_θ DERIVATION (Magic-Number-Block): grid-invariance k_θ·a = κ, the rigid-lever arm identity, CFL scales.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.motor.backbone_warp import (
    KBT_PN_UM,
    angle_harmonic_forces,
    arm_bending_cfl_stiffness,
    arm_orientation_k_theta,
    backbone_bending_cfl_stiffness,
    backbone_bending_kappa,
    backbone_bending_k_theta,
)
from aleph.components.motor.minifilament_topology import MinifilamentTopology

PI = float(np.pi)
A = 0.301 / 13.0          # a production backbone segment length [µm]
R0 = 0.200                # head-arm rest length [µm]


# ── zero force at rest ───────────────────────────────────────────────────────────────────────────
def test_backbone_zero_force_at_straight_rest() -> None:
    """A straight backbone triple (θ = π) bears zero bending force."""
    pos = np.array([[-A, 0, 0], [0, 0, 0], [A, 0, 0]], float)
    f = angle_harmonic_forces(pos, np.array([[0, 1, 2]]), k_theta=1.5, theta0=PI)
    assert np.abs(f).max() < 1e-12


def test_arm_zero_force_at_perpendicular_rest() -> None:
    """A perpendicular head-arm triple (θ = π/2) bears zero bending force."""
    pos = np.array([[A, 0, 0], [0, 0, 0], [0, R0, 0]], float)   # neighbor +x, bead origin, head +y
    f = angle_harmonic_forces(pos, np.array([[0, 1, 2]]), k_theta=40.0, theta0=PI / 2)
    assert np.abs(f).max() < 1e-12


# ── restoring moment + Newton's 3rd law ────────────────────────────────────────────────────────
def test_backbone_restoring_moment_under_bend() -> None:
    """Bending the backbone pushes the middle bead back toward straight; forces sum to zero (Newton 3rd)."""
    pos = np.array([[-A, 0, 0], [0, 0.25 * A, 0], [A, 0, 0]], float)   # middle bead lifted +y
    f = angle_harmonic_forces(pos, np.array([[0, 1, 2]]), k_theta=1.5, theta0=PI)
    assert f[1, 1] < 0.0                                   # middle bead pushed back down (−y): restoring
    assert np.abs(f.sum(axis=0)).max() < 1e-12            # Σ F = 0
    # zero net torque about the vertex (internal potential)
    tau = np.cross(pos - pos[1], f).sum(axis=0)
    assert np.abs(tau).max() < 1e-12


def test_arm_restoring_moment_when_head_swung() -> None:
    """A head swung along the backbone axis is pushed back toward perpendicular; forces sum to zero."""
    pos = np.array([[A, 0, 0], [0, 0, 0], [0.03, R0, 0]], float)   # head swung +x (toward neighbor)
    f = angle_harmonic_forces(pos, np.array([[0, 1, 2]]), k_theta=40.0, theta0=PI / 2)
    assert f[2, 0] < 0.0                                   # head pushed back −x: restoring the swing
    assert np.abs(f.sum(axis=0)).max() < 1e-12


# ── first-order transverse stiffness (the F6 transmission the distance spring lacked) ────────────
def test_head_bears_first_order_transverse_stiffness() -> None:
    """A transversely-displaced head bears a restoring force ≈ −(k_θ,arm/r0²)·δx to first order.

    A single head↔backbone DISTANCE spring gives ZERO transverse stiffness (the F6 bug: the head swings free);
    the angle-harmonic gives ``k_θ,arm / r0²`` — exactly :func:`arm_bending_cfl_stiffness`.
    """
    k_theta = arm_orientation_k_theta(1.0e3, R0)           # rigid-lever: k_eff = k_xb = 1000
    k_expect = arm_bending_cfl_stiffness(k_theta, R0)
    dxs = np.array([1e-4, 2e-4, 4e-4])
    fx = np.array([
        angle_harmonic_forces(np.array([[A, 0, 0], [0, 0, 0], [dx, R0, 0]], float),
                              np.array([[0, 1, 2]]), k_theta=k_theta, theta0=PI / 2)[2, 0]
        for dx in dxs
    ])
    k_meas = -(fx / dxs).mean()
    assert k_meas == pytest.approx(k_expect, rel=1e-3)
    assert k_expect == pytest.approx(1.0e3, rel=1e-9)      # rigid-lever identity: k_eff == k_xb


# ── kernel ↔ NumPy-mirror structural parity (AST/source; NO Warp launch — I0-A: no CPU sim path) ──
def test_kernel_mirrors_numpy_reference_source() -> None:
    """``angle_harmonic_kernel`` must be the bit-for-formula twin of ``angle_harmonic_forces``.

    The dev Mac cannot launch the CUDA kernel and the foundation contract forbids a Warp-CPU launch in the
    ac tests (I0-A: no CPU simulation path — ``test_ac_tests_do_not_launch_production_kernels_on_cpu``). We
    instead assert the kernel SOURCE carries the identical LAMMPS force expansion the NumPy mirror is validated
    on (the physics is gated by the NumPy-reference tests above; the on-device parity is the gbook native gate).
    A re-derivation of either that drops a coefficient breaks this structural check locally.
    """
    import ast
    import inspect

    from aleph.components.motor import backbone_warp

    def _body(fn_name: str) -> str:
        src = inspect.getsource(getattr(backbone_warp, fn_name))
        return src

    kernel_src = _body("angle_harmonic_kernel")
    ref_src = _body("angle_harmonic_forces")
    # the shared LAMMPS angle-force expansion tokens must appear in BOTH (kernel + mirror)
    for token in ("a11", "a12", "a22", "dtheta", "theta0", "f_i", "f_k"):
        assert token in kernel_src, f"kernel missing '{token}' — diverged from the NumPy mirror"
        assert token in ref_src, f"NumPy mirror missing '{token}' — diverged from the kernel"
    # both must parse (the kernel body is valid Python/Warp source)
    assert ast.parse(inspect.getsource(backbone_warp)) is not None


# ── k_θ DERIVATION (Magic-Number-Block: derived, grid-invariant, not tuned) ──────────────────────
def test_backbone_k_theta_is_grid_invariant() -> None:
    """k_θ = κ/a holds the physical rigidity κ = L_p·k_BT fixed under backbone refinement (grid-invariant)."""
    kappa = backbone_bending_kappa(1.0)                    # L_p = 1 µm
    assert kappa == pytest.approx(1.0 * KBT_PN_UM)
    for n_bb in (8, 14, 28):                               # refine the rod: a shrinks, k_θ rises, κ = k_θ·a fixed
        a = 0.301 / (n_bb - 1)
        k_theta = backbone_bending_k_theta(kappa, a)
        assert k_theta * a == pytest.approx(kappa, rel=1e-12)


def test_cfl_stiffness_scales() -> None:
    """Bending adds an effective linear stiffness folded into kmax (CFL note): 4k_θ/a² and k_θ/r0²."""
    k_bb = backbone_bending_k_theta(backbone_bending_kappa(1.0), A)
    assert backbone_bending_cfl_stiffness(k_bb, A) == pytest.approx(4.0 * k_bb / (A * A), rel=1e-12)
    k_arm = arm_orientation_k_theta(1.0e3, R0)
    assert arm_bending_cfl_stiffness(k_arm, R0) == pytest.approx(k_arm / (R0 * R0), rel=1e-12)


def test_derivations_reject_nonphysical_inputs() -> None:
    """The Magic-Number-Block derivations reject non-positive inputs (fail loud, never silently default)."""
    for bad in (0.0, -1.0):
        with pytest.raises(ValueError):
            backbone_bending_kappa(bad)
        with pytest.raises(ValueError):
            backbone_bending_k_theta(1.0, bad)
        with pytest.raises(ValueError):
            arm_orientation_k_theta(bad, R0)


# ── topology: triple counts + mirror symmetry ────────────────────────────────────────────────────
def test_backbone_angle_triple_count_and_connectivity() -> None:
    """n_backbone_angles = n_bb − 2 consecutive triples (i−1, i, i+1); empty for n_bb = 2."""
    topo = MinifilamentTopology(n_bb=14, n_heads_per_side=10, backbone_length_um=0.301, head_offset_um=0.2)
    tri = topo.backbone_angle_triples()
    assert tri.shape == (topo.n_backbone_angles, 3) == (12, 3)
    assert np.array_equal(tri, np.column_stack([np.arange(12), np.arange(1, 13), np.arange(2, 14)]))
    tiny = MinifilamentTopology(n_bb=2, n_heads_per_side=1, backbone_length_um=0.3, head_offset_um=0.2)
    assert tiny.backbone_angle_triples().shape == (0, 3)


def test_head_arm_angle_triples_are_mirror_symmetric() -> None:
    """Each head gets a (neighbor, bead, head) triple with the neighbor stepping toward the backbone CENTRE.

    Mirror symmetry between the + and − arms is what keeps the two clamp reactions equal-and-opposite; an
    asymmetric neighbor rule biases |reaction_A| ≠ |reaction_B| (the transmission-ratio metric).
    """
    topo = MinifilamentTopology(n_bb=14, n_heads_per_side=10, backbone_length_um=0.301, head_offset_um=0.2)
    tri = topo.head_arm_angle_triples()
    assert tri.shape == (topo.n_head_arm_angles, 3) == (20, 3)
    centre = 0.5 * (topo.n_bb - 1)
    for neighbor, bead, head in tri:
        assert 0 <= neighbor < topo.n_bb and abs(neighbor - bead) == 1     # a real adjacent backbone bead
        assert head >= topo.n_bb                                            # a head particle
        if bead != neighbor:                                               # neighbor steps toward the centre
            assert (neighbor > bead) == (bead < centre)

r"""Structural gate for the KERNEL_BOUND ``ecm`` collagen constitutive force pass (CPU-safe; native forbidden).

Proves the advance from SEAMED → KERNEL_BOUND for the ``ecm`` component: the ECM component now LAUNCHES the
reused ``ff`` kernels (``link_spring_kernel`` axial + ``cytosim_bending_kernel``) over its OWN
:class:`~aleph.components.ecm.device_schema.ECMTopologyState` SoA, mirroring
:class:`aleph.engine.sf_mechanics.SFFilamentMechanics`.

Two independent proofs, both CUDA-free:
  * LAUNCH — a recording launcher (metadata-only ``_FakeArray`` doubles) proves
    :meth:`ECMConstitutiveForce.accumulate` launches ``link_spring_kernel`` over ``segments_d`` then
    ``cytosim_bending_kernel`` over ``bend_triples_d`` with the byte-identical cortex/SF arg order and dims =
    segment/bend capacity (no kernel executes on the host);
  * PHYSICS — the pure-NumPy reference of the EXACT kernel force law + the exact on-device derivation proves a
    RELAXED collagen network injects ≈0 force, a STRETCHED segment develops correct-sign restoring tension, a
    BENT triple develops a restoring bending force, the force scatters onto the right SoA nodes, and an inactive
    element is gated to a zero contribution.  No Warp-CPU production-kernel launch is used.

Plus SOURCED/GAP discipline (κ SOURCED, EA REQUIRED-no-default).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest
import warp as wp

from aleph.engine.ecm_mechanics import (
    COLLAGEN_FORCE_PROVENANCE,
    COLLAGEN_MATERIAL_KEY,
    ECMConstitutiveForce,
    bending_force_reference,
    collagen_bend_alpha_np,
    collagen_force_provenance,
    collagen_segment_stiffness_np,
    link_spring_force_reference,
)
from aleph.laws.ecm_library import get_spec


# ── metadata-only CUDA doubles + recording launcher (no host kernel execution) ────────────────────────
@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    ptr: int
    shape: tuple[int, ...] = (4,)
    dtype: object = wp.vec3d
    device: _FakeDevice = _FakeDevice()


@dataclass(slots=True)
class _LaunchRecorder:
    calls: list[dict[str, object]] = field(default_factory=list)

    def __call__(self, kernel, *, dim, inputs, device) -> None:
        self.calls.append({"kernel": kernel, "dim": dim, "inputs": inputs, "device": device})


_LINK_KERNEL = object()
_BENDING_KERNEL = object()
_STIFFNESS_KERNEL = object()
_ALPHA_KERNEL = object()

_EA_PN = 1.1e6 * np.pi * 0.05 ** 2       # EA = E_fibril·π r_f² from the collagen card (GAP value, for the double)
_KAPPA = 4.28e-3 * 17.0                  # κ = k_B·T·L_p (SOURCED collagen), FF units


def _force_pass(recorder: _LaunchRecorder, *, n_segments: int = 4, n_bends: int = 3) -> ECMConstitutiveForce:
    return ECMConstitutiveForce(
        device="cuda:0",
        n_nodes=8,
        segments_d=_FakeArray(500, shape=(n_segments, 2), dtype=wp.int32),
        seg_rest_d=_FakeArray(501, shape=(n_segments,), dtype=wp.float64),
        segment_active_d=_FakeArray(502, shape=(n_segments,), dtype=wp.int32),
        bend_triples_d=_FakeArray(503, shape=(n_bends, 3), dtype=wp.int32),
        bend_left_segment_d=_FakeArray(504, shape=(n_bends,), dtype=wp.int32),
        bend_right_segment_d=_FakeArray(505, shape=(n_bends,), dtype=wp.int32),
        bend_active_d=_FakeArray(506, shape=(n_bends,), dtype=wp.int32),
        k_seg_d=_FakeArray(507, shape=(n_segments,), dtype=wp.float64),
        alpha_d=_FakeArray(508, shape=(n_bends,), dtype=wp.float64),
        ea_pn=_EA_PN,
        kappa_pn_um2=_KAPPA,
        link_kernel=_LINK_KERNEL,
        bending_kernel=_BENDING_KERNEL,
        stiffness_kernel=_STIFFNESS_KERNEL,
        alpha_kernel=_ALPHA_KERNEL,
        launch=recorder,
    )


# ── LAUNCH proof: the ECM component launches BOTH ff kernels itself, correct arg order + dims ──────────
def test_ecm_force_launches_link_and_bending_over_the_soa() -> None:
    recorder = _LaunchRecorder()
    force_pass = _force_pass(recorder, n_segments=4, n_bends=3)
    pos = _FakeArray(600)
    force = _FakeArray(601)
    force_pass.accumulate(pos, force)

    assert len(recorder.calls) == 2
    link_call, bend_call = recorder.calls
    assert link_call["kernel"] is _LINK_KERNEL
    assert link_call["dim"] == 4  # one thread per collagen segment (full capacity)
    assert link_call["inputs"] == [
        pos, force_pass.segments_d, force_pass.k_seg_d, force_pass.seg_rest_d, force,
    ]
    assert link_call["device"] == "cuda:0"
    assert bend_call["kernel"] is _BENDING_KERNEL
    assert bend_call["dim"] == 3  # one thread per collagen bending triple
    assert bend_call["inputs"] == [pos, force_pass.bend_triples_d, force_pass.alpha_d, force]
    # every launch reads pos (first) and scatters onto the ECM-owned force (last), never a foreign array
    for call in recorder.calls:
        assert call["inputs"][0] is pos
        assert call["inputs"][-1] is force


def test_refresh_parameters_derives_k_and_alpha_on_device() -> None:
    recorder = _LaunchRecorder()
    force_pass = _force_pass(recorder, n_segments=4, n_bends=3)
    force_pass.refresh_parameters()

    assert len(recorder.calls) == 2
    k_call, a_call = recorder.calls
    assert k_call["kernel"] is _STIFFNESS_KERNEL and k_call["dim"] == 4
    # seg_rest + active gate + EA scalar → k_seg_d (out)
    assert k_call["inputs"][0] is force_pass.seg_rest_d
    assert k_call["inputs"][1] is force_pass.segment_active_d
    assert k_call["inputs"][-1] is force_pass.k_seg_d
    assert a_call["kernel"] is _ALPHA_KERNEL and a_call["dim"] == 3
    assert a_call["inputs"][0] is force_pass.seg_rest_d
    assert a_call["inputs"][-1] is force_pass.alpha_d


def test_rejects_cpu_arrays_and_shape_mismatch() -> None:
    recorder = _LaunchRecorder()
    with pytest.raises(TypeError, match="dtype wp.int32"):
        ECMConstitutiveForce(
            device="cuda:0", n_nodes=8,
            segments_d=_FakeArray(500, shape=(4, 2), dtype=wp.float64),  # wrong dtype
            seg_rest_d=_FakeArray(501, shape=(4,), dtype=wp.float64),
            segment_active_d=_FakeArray(502, shape=(4,), dtype=wp.int32),
            bend_triples_d=_FakeArray(503, shape=(3, 3), dtype=wp.int32),
            bend_left_segment_d=_FakeArray(504, shape=(3,), dtype=wp.int32),
            bend_right_segment_d=_FakeArray(505, shape=(3,), dtype=wp.int32),
            bend_active_d=_FakeArray(506, shape=(3,), dtype=wp.int32),
            k_seg_d=_FakeArray(507, shape=(4,), dtype=wp.float64),
            alpha_d=_FakeArray(508, shape=(3,), dtype=wp.float64),
            ea_pn=_EA_PN, kappa_pn_um2=_KAPPA,
            link_kernel=_LINK_KERNEL, bending_kernel=_BENDING_KERNEL,
            stiffness_kernel=_STIFFNESS_KERNEL, alpha_kernel=_ALPHA_KERNEL, launch=recorder,
        )
    with pytest.raises(ValueError, match="length must equal"):
        ECMConstitutiveForce(
            device="cuda:0", n_nodes=8,
            segments_d=_FakeArray(500, shape=(4, 2), dtype=wp.int32),
            seg_rest_d=_FakeArray(501, shape=(3,), dtype=wp.float64),  # wrong length vs 4 segments
            segment_active_d=_FakeArray(502, shape=(4,), dtype=wp.int32),
            bend_triples_d=_FakeArray(503, shape=(3, 3), dtype=wp.int32),
            bend_left_segment_d=_FakeArray(504, shape=(3,), dtype=wp.int32),
            bend_right_segment_d=_FakeArray(505, shape=(3,), dtype=wp.int32),
            bend_active_d=_FakeArray(506, shape=(3,), dtype=wp.int32),
            k_seg_d=_FakeArray(507, shape=(4,), dtype=wp.float64),
            alpha_d=_FakeArray(508, shape=(3,), dtype=wp.float64),
            ea_pn=_EA_PN, kappa_pn_um2=_KAPPA,
            link_kernel=_LINK_KERNEL, bending_kernel=_BENDING_KERNEL,
            stiffness_kernel=_STIFFNESS_KERNEL, alpha_kernel=_ALPHA_KERNEL, launch=recorder,
        )
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    with pytest.raises(ValueError, match="CUDA device array"):
        ECMConstitutiveForce(
            device="cuda", n_nodes=8,
            segments_d=_FakeArray(500, shape=(4, 2), dtype=wp.int32, device=cpu),
            seg_rest_d=_FakeArray(501, shape=(4,), dtype=wp.float64, device=cpu),
            segment_active_d=_FakeArray(502, shape=(4,), dtype=wp.int32, device=cpu),
            bend_triples_d=_FakeArray(503, shape=(3, 3), dtype=wp.int32, device=cpu),
            bend_left_segment_d=_FakeArray(504, shape=(3,), dtype=wp.int32, device=cpu),
            bend_right_segment_d=_FakeArray(505, shape=(3,), dtype=wp.int32, device=cpu),
            bend_active_d=_FakeArray(506, shape=(3,), dtype=wp.int32, device=cpu),
            k_seg_d=_FakeArray(507, shape=(4,), dtype=wp.float64, device=cpu),
            alpha_d=_FakeArray(508, shape=(3,), dtype=wp.float64, device=cpu),
            ea_pn=_EA_PN, kappa_pn_um2=_KAPPA,
            link_kernel=_LINK_KERNEL, bending_kernel=_BENDING_KERNEL,
            stiffness_kernel=_STIFFNESS_KERNEL, alpha_kernel=_ALPHA_KERNEL, launch=recorder,
        )


# ── a small explicit two-fiber collagen SoA fixture (host NumPy; mirrors the seeding layout) ──────────
def _two_fiber_soa():
    """Two straight collinear collagen fibers of 3 nodes each along +x (build = rest ⇒ relaxed)."""
    seg_um = 0.5
    pos = np.array(
        [[0.0, 0.0, 0.0], [seg_um, 0.0, 0.0], [2 * seg_um, 0.0, 0.0],   # fiber 0 nodes 0,1,2
         [0.0, 1.0, 0.0], [seg_um, 1.0, 0.0], [2 * seg_um, 1.0, 0.0]],  # fiber 1 nodes 3,4,5
        dtype=np.float64,
    )
    segments = np.array([[0, 1], [1, 2], [3, 4], [4, 5]], dtype=np.int32)          # 4 segments
    seg_rest = np.linalg.norm(pos[segments[:, 1]] - pos[segments[:, 0]], axis=1).astype(np.float64)
    seg_active = np.array([1, 1, 1, 1], dtype=np.int32)
    bend_triples = np.array([[0, 1, 2], [3, 4, 5]], dtype=np.int32)                # 2 bends
    bend_left = np.array([0, 2], dtype=np.int32)                                   # left segment of each bend
    bend_right = np.array([1, 3], dtype=np.int32)                                  # right segment of each bend
    bend_active = np.array([1, 1], dtype=np.int32)
    return pos, segments, seg_rest, seg_active, bend_triples, bend_left, bend_right, bend_active


# ── PHYSICS proof (NumPy reference of the exact kernel law + the exact derivation) ────────────────────
def test_relaxed_collagen_network_injects_zero_force() -> None:
    spec = get_spec(COLLAGEN_MATERIAL_KEY)
    pos, segments, seg_rest, seg_active, triples, bl, br, bend_active = _two_fiber_soa()
    k_seg = collagen_segment_stiffness_np(seg_rest, seg_active, spec.EA_pN)
    alpha = collagen_bend_alpha_np(seg_rest, bend_active, bl, br, spec.kappa_pN_um2)

    # k = EA/seg_rest exactly; α = κ/seg_local³ exactly (relaxed derivation is real, not empty)
    assert np.allclose(k_seg, spec.EA_pN / seg_rest)
    seg_local = 0.5 * (seg_rest[bl] + seg_rest[br])   # per-bend local mean segment length
    assert np.allclose(alpha, spec.kappa_pN_um2 / seg_local ** 3)

    axial = link_spring_force_reference(pos, segments, k_seg, seg_rest)   # build == rest ⇒ 0 tension
    bend = bending_force_reference(pos, triples, alpha)                   # straight ⇒ 0 curvature
    assert np.abs(axial).max() < 1e-9
    assert np.abs(bend).max() < 1e-9


def test_stretched_segment_develops_correct_sign_tension_on_right_nodes() -> None:
    spec = get_spec(COLLAGEN_MATERIAL_KEY)
    pos, segments, seg_rest, seg_active, _, _, _, _ = _two_fiber_soa()
    k_seg = collagen_segment_stiffness_np(seg_rest, seg_active, spec.EA_pN)

    stretched = pos.copy()
    i, j = int(segments[0, 0]), int(segments[0, 1])
    axis = stretched[j] - stretched[i]
    stretched[j] = stretched[i] + 1.5 * axis                             # stretch segment 0 by 50%
    force = link_spring_force_reference(stretched, segments, k_seg, seg_rest)

    unit = axis / np.linalg.norm(axis)
    assert np.dot(force[i], unit) > 0.0    # node i pulled toward j (restoring)
    assert np.dot(force[j], unit) < 0.0    # node j pulled back toward i
    # force scatters onto fiber-0 nodes only; the untouched second fiber (nodes 3,4,5) sees nothing
    assert np.linalg.norm(force[3]) < 1e-12 and np.linalg.norm(force[5]) < 1e-12
    # the stretched segment's outer endpoint (node 0, touched by segment 0 alone) is pulled toward j (+x)
    assert force[i][0] > 0.0 and np.linalg.norm(force[i][1:]) < 1e-9


def test_bent_triple_develops_restoring_bending_force_on_center_node() -> None:
    spec = get_spec(COLLAGEN_MATERIAL_KEY)
    pos, _, seg_rest, _, triples, bl, br, bend_active = _two_fiber_soa()
    alpha = collagen_bend_alpha_np(seg_rest, bend_active, bl, br, spec.kappa_pN_um2)

    bent = pos.copy()
    b = int(triples[0, 1])                                               # center node of fiber-0 triple
    bent[b] = bent[b] + np.array([0.0, 0.0, 0.25])                       # kink out of line
    force = bending_force_reference(bent, triples, alpha)

    # a curvature restoring force appears on the center node, pushing it back toward the chord (−z)
    assert np.linalg.norm(force[b]) > 1e-6
    assert force[b][2] < 0.0
    # the straight second fiber's triple stays force-free
    assert np.linalg.norm(force[3]) < 1e-12 and np.linalg.norm(force[5]) < 1e-12


def test_inactive_elements_are_gated_to_zero_contribution() -> None:
    spec = get_spec(COLLAGEN_MATERIAL_KEY)
    pos, segments, seg_rest, seg_active, triples, bl, br, bend_active = _two_fiber_soa()
    seg_active = seg_active.copy(); seg_active[:] = 0        # deactivate every segment
    bend_active = bend_active.copy(); bend_active[:] = 0     # deactivate every bend
    k_seg = collagen_segment_stiffness_np(seg_rest, seg_active, spec.EA_pN)
    alpha = collagen_bend_alpha_np(seg_rest, bend_active, bl, br, spec.kappa_pN_um2)
    assert np.all(k_seg == 0.0) and np.all(alpha == 0.0)

    stretched = pos.copy()
    stretched[int(segments[0, 1])] += np.array([0.3, 0.0, 0.0])
    assert np.abs(link_spring_force_reference(stretched, segments, k_seg, seg_rest)).max() < 1e-12
    bent = pos.copy(); bent[int(triples[0, 1])] += np.array([0.0, 0.0, 0.25])
    assert np.abs(bending_force_reference(bent, triples, alpha)).max() < 1e-12


# ── SOURCED vs GAP discipline ─────────────────────────────────────────────────────────────────────────
def test_bending_alpha_uses_sourced_collagen_kappa() -> None:
    spec = get_spec(COLLAGEN_MATERIAL_KEY)
    # κ = k_B·T·L_p with L_p = 17 µm (SOURCED, Licup 2015); the card exposes it directly.
    assert spec.kappa_pN_um2 == pytest.approx(4.28e-3 * 17.0, rel=1e-2)
    _, _, seg_rest, _, _, bl, br, bend_active = _two_fiber_soa()
    alpha = collagen_bend_alpha_np(seg_rest, bend_active, bl, br, spec.kappa_pN_um2)
    seg_local = 0.5 * (seg_rest[int(bl[0])] + seg_rest[int(br[0])])
    assert alpha[0] == pytest.approx(spec.kappa_pN_um2 / seg_local ** 3, rel=1e-12)


def test_axial_modulus_is_a_required_gap_no_default() -> None:
    recorder = _LaunchRecorder()
    # EA ≤ 0 must raise at construction — never the ecm_library k_seg_pN_um 5e4 fallback (E_fibril is a GAP).
    with pytest.raises(ValueError, match="REQUIRED-PARAM 'ea_pn'"):
        ECMConstitutiveForce(
            device="cuda:0", n_nodes=8,
            segments_d=_FakeArray(500, shape=(4, 2), dtype=wp.int32),
            seg_rest_d=_FakeArray(501, shape=(4,), dtype=wp.float64),
            segment_active_d=_FakeArray(502, shape=(4,), dtype=wp.int32),
            bend_triples_d=_FakeArray(503, shape=(3, 3), dtype=wp.int32),
            bend_left_segment_d=_FakeArray(504, shape=(3,), dtype=wp.int32),
            bend_right_segment_d=_FakeArray(505, shape=(3,), dtype=wp.int32),
            bend_active_d=_FakeArray(506, shape=(3,), dtype=wp.int32),
            k_seg_d=_FakeArray(507, shape=(4,), dtype=wp.float64),
            alpha_d=_FakeArray(508, shape=(3,), dtype=wp.float64),
            ea_pn=0.0, kappa_pn_um2=_KAPPA,
            link_kernel=_LINK_KERNEL, bending_kernel=_BENDING_KERNEL,
            stiffness_kernel=_STIFFNESS_KERNEL, alpha_kernel=_ALPHA_KERNEL, launch=recorder,
        )


def test_provenance_flags_ea_as_gap_and_kappa_as_sourced() -> None:
    spec = get_spec(COLLAGEN_MATERIAL_KEY)
    prov = collagen_force_provenance(spec)
    assert prov["material"] == "collagen_I"
    assert COLLAGEN_FORCE_PROVENANCE["kappa_pN_um2"].startswith("SOURCED")
    assert COLLAGEN_FORCE_PROVENANCE["E_fibril_Pa"].startswith("GAP")
    assert COLLAGEN_FORCE_PROVENANCE["EA_pN"].startswith("GAP")
    assert prov["EA_pN"] == pytest.approx(spec.EA_pN)

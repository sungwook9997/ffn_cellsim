"""Oracle-B blebbing perturbation — local ERM de-adhesion cap selection + commit-safety.

CPU tests exercise the pure-NumPy acceptance oracle (cap geometry, monotonicity); a CUDA-gated test asserts
the Warp kernel reproduces the oracle bit-for-bit and honours the accept-gate.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.components.incumbent.bleb_perturbation import (
    cap_cos_half_angle,
    count_tethers_in_cap,
    erm_patch_deadhesion_kernel,
    erm_patch_deadhesion_reference,
)

_CUDA_DEVICE = next((d for d in wp.get_devices() if d.is_cuda), None)
_Z = np.array([0.0, 0.0, 1.0])
_ORIGIN = np.zeros(3)


# ── CPU: cap geometry + monotonicity (the acceptance oracle) ─────────────────────────────────────────

def test_axis_pole_deadheres_opposite_pole_survives() -> None:
    pos = np.array([[0, 0, 1.0], [0, 0, -1.0]])          # +z pole, −z pole
    bound = np.ones(2, np.int32)
    out = erm_patch_deadhesion_reference(pos, bound, _ORIGIN, _Z, cap_cos_half_angle(30.0))
    assert out[0] == 0 and out[1] == 1


def test_cap_boundary_is_inclusive() -> None:
    half = 30.0
    r = np.deg2rad(half)
    on_edge = np.array([[np.sin(r), 0.0, np.cos(r)]])    # exactly at the 30° rim
    just_out = np.array([[np.sin(r + 1e-3), 0.0, np.cos(r + 1e-3)]])
    cos_half = cap_cos_half_angle(half)
    assert erm_patch_deadhesion_reference(on_edge, np.ones(1, np.int32), _ORIGIN, _Z, cos_half)[0] == 0
    assert erm_patch_deadhesion_reference(just_out, np.ones(1, np.int32), _ORIGIN, _Z, cos_half)[0] == 1


def test_deadhesion_is_monotone_and_idempotent() -> None:
    rng = np.random.default_rng(0)
    pos = rng.standard_normal((200, 3))
    bound = (rng.random(200) > 0.3).astype(np.int32)     # some already unbound
    cos_half = cap_cos_half_angle(45.0)
    once = erm_patch_deadhesion_reference(pos, bound, _ORIGIN, _Z, cos_half)
    twice = erm_patch_deadhesion_reference(pos, once, _ORIGIN, _Z, cos_half)
    assert np.all(once <= bound), "bound may only go 1->0"
    assert np.array_equal(once, twice), "replay must be idempotent"
    # an already-unbound tether inside the cap is never reactivated
    assert np.all(once[bound == 0] == 0)


def test_count_scales_with_patch_angle() -> None:
    rng = np.random.default_rng(1)
    pos = rng.standard_normal((2000, 3))
    bound = np.ones(2000, np.int32)
    n20 = count_tethers_in_cap(pos, bound, _ORIGIN, _Z, cap_cos_half_angle(20.0))
    n60 = count_tethers_in_cap(pos, bound, _ORIGIN, _Z, cap_cos_half_angle(60.0))
    assert 0 < n20 < n60 < 2000, "a wider cap de-adheres strictly more tethers"


def test_cap_cos_half_angle_rejects_out_of_range() -> None:
    for bad in (0.0, 90.0, 120.0, -5.0):
        with pytest.raises(ValueError):
            cap_cos_half_angle(bad)


# ── CUDA: the kernel reproduces the oracle and honours the accept-gate ────────────────────────────────

@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_kernel_matches_reference_and_respects_accept() -> None:
    dev = str(_CUDA_DEVICE)
    rng = np.random.default_rng(7)
    ne = 500
    mem_pos = rng.standard_normal((ne, 3)) + np.array([2.0, 0.0, 0.0])   # off-origin centroid too
    centroid = mem_pos.mean(axis=0)
    axis = np.array([0.3, -0.4, 0.866]); axis /= np.linalg.norm(axis)
    bound0 = (rng.random(ne) > 0.2).astype(np.int32)
    cos_half = cap_cos_half_angle(35.0)
    ref = erm_patch_deadhesion_reference(mem_pos, bound0, centroid, axis, cos_half)

    with wp.ScopedDevice(dev):
        pos = wp.array(mem_pos, dtype=wp.vec3d, device=dev)
        midx = wp.array(np.arange(ne, dtype=np.int32), dtype=wp.int32, device=dev)
        bound = wp.array(bound0.copy(), dtype=wp.int32, device=dev)
        # accepted == 0 -> no-op (rejected mechanical candidate leaves no bleb)
        rej = wp.zeros(1, dtype=wp.int32, device=dev)
        wp.launch(erm_patch_deadhesion_kernel, dim=ne,
                  inputs=[pos, midx, bound, wp.vec3d(*centroid), wp.vec3d(*axis),
                          wp.float64(cos_half), rej], device=dev)
        assert np.array_equal(bound.numpy(), bound0), "accepted=0 must not touch bound"
        # accepted == 1 -> matches the oracle
        acc = wp.ones(1, dtype=wp.int32, device=dev)
        wp.launch(erm_patch_deadhesion_kernel, dim=ne,
                  inputs=[pos, midx, bound, wp.vec3d(*centroid), wp.vec3d(*axis),
                          wp.float64(cos_half), acc], device=dev)
        assert np.array_equal(bound.numpy(), ref), "kernel must match the NumPy oracle"

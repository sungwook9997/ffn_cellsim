"""CUDA parity gates for the NF2007 variable-length constraint-force projector."""

from __future__ import annotations

import numpy as np
import pytest

wp = pytest.importorskip("warp")

from aleph.components.incumbent.inner_mechanics import project_constraint_forces_kernel  # noqa: E402
from aleph.laws.constraints import _project_one_fiber  # noqa: E402

wp.init()
_CUDA_DEVICE = next((device for device in wp.get_devices() if device.is_cuda), None)
pytestmark = pytest.mark.skipif(
    _CUDA_DEVICE is None,
    reason="I0-A: constraint projector kernel gate requires a CUDA GPU",
)
DEV = str(_CUDA_DEVICE) if _CUDA_DEVICE is not None else ""


def test_constraint_projector_matches_independent_dense_oracle() -> None:
    """Variable 4/3-node fibers match the NumPy dense solve and preserve net force."""
    pos = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.1, 0.0],
        [1.0, 0.25, 0.05],
        [1.45, 0.5, 0.1],
        [-0.2, 0.3, 0.0],
        [0.1, 0.75, 0.1],
        [0.25, 1.2, 0.3],
    ], dtype=np.float64)
    force = np.array([
        [1.0, -0.2, 0.1],
        [-0.3, 0.7, -0.4],
        [0.2, -0.6, 0.8],
        [-0.9, 0.1, -0.5],
        [0.4, 0.3, -0.2],
        [-0.7, 0.5, 0.9],
        [0.1, -0.8, -0.7],
    ], dtype=np.float64)
    fiber_off = np.array([0, 4, 7], dtype=np.int32)
    seg_off = np.array([0, 3, 5], dtype=np.int32)

    expected = force.copy()
    expected[0:4] = _project_one_fiber(pos[0:4], force[0:4], np.ones(3))
    expected[4:7] = _project_one_fiber(pos[4:7], force[4:7], np.ones(2))

    pos_d = wp.array(pos, dtype=wp.vec3d, device=DEV)
    force_d = wp.array(force, dtype=wp.vec3d, device=DEV)
    projected_d = wp.empty_like(force_d)
    wp.copy(projected_d, force_d)
    finite_d = wp.ones(1, dtype=wp.int32, device=DEV)
    wp.launch(
        project_constraint_forces_kernel,
        dim=2,
        inputs=[
            pos_d,
            force_d,
            wp.array(fiber_off, dtype=wp.int32, device=DEV),
            wp.array(seg_off, dtype=wp.int32, device=DEV),
            projected_d,
            wp.zeros(5, dtype=wp.float64, device=DEV),
            wp.zeros(5, dtype=wp.float64, device=DEV),
            finite_d,
        ],
        device=DEV,
    )
    wp.synchronize_device(DEV)
    measured = projected_d.numpy()

    np.testing.assert_allclose(measured, expected, rtol=2.0e-13, atol=2.0e-13)
    np.testing.assert_allclose(measured.sum(axis=0), force.sum(axis=0), rtol=0.0, atol=5.0e-15)
    assert int(finite_d.numpy()[0]) == 1

    # Projected motion is tangent: J(Pf) = 0 independently on every segment.
    for start, stop in ((0, 4), (4, 7)):
        segment = pos[start + 1:stop] - pos[start:stop - 1]
        jpf = 2.0 * np.einsum(
            "ij,ij->i",
            segment,
            measured[start + 1:stop] - measured[start:stop - 1],
        )
        np.testing.assert_allclose(jpf, 0.0, rtol=0.0, atol=2.0e-14)

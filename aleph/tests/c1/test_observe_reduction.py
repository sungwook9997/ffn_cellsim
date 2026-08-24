"""The device observation reduction must equal the obvious host arithmetic, on a small field.

WHY THIS EXISTS.  The observers moved onto the GPU for a real reason — a 513^3 field is 1.08 GB and 120
samples per point would have moved 130 GB across PCIe to compute three scalars — but a hand-written
atomic reduction is exactly the kind of code that is off by a mean, a square root, or a cell count while
still producing a clean-looking exponential decay.  Nothing downstream would catch that: the scorer
would fit a crossing, the log-log fit would return a slope, and the number would be wrong.

So the reduction is checked against ``numpy`` on a field small enough that the host path is trivially
right, with a NON-UNIFORM field and a NON-ZERO mean so a bug in the mean subtraction cannot cancel.

This runs on whatever device Warp resolves, which on the dev host is CPU.  That is not a physics
measurement — no law is evaluated, no timestep is taken, and nothing here may be quoted as a result.
It is a check that two ways of adding the same numbers agree.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.scripts.c1_biot_timescale import _observe_kernel, _sum_kernel

N = 17


def _field() -> np.ndarray:
    """A lumpy, strictly non-uniform field with a non-zero mean."""
    ax = np.linspace(-1.0, 1.0, N)
    x, y, z = np.meshgrid(ax, ax, ax, indexing="ij")
    return (0.7 + np.exp(-(x ** 2 + y ** 2 + z ** 2) * 3.0) + 0.1 * np.sin(5.0 * x) * np.cos(3.0 * y))


def _patch() -> np.ndarray:
    ax = np.linspace(-1.0, 1.0, N)
    x, y, z = np.meshgrid(ax, ax, ax, indexing="ij")
    return ((x ** 2 + y ** 2 + z ** 2) <= 0.25).astype(np.int32)


def test_device_reduction_matches_numpy() -> None:
    p_np, patch_np = _field(), _patch()
    dev = str(wp.get_device())
    p = wp.array(np.ascontiguousarray(p_np), dtype=wp.float64, device=dev)
    patch = wp.array(np.ascontiguousarray(patch_np), dtype=wp.int32, device=dev)

    total = wp.zeros(1, dtype=wp.float64, device=dev)
    wp.launch(_sum_kernel, dim=(N, N, N), inputs=[p, total], device=dev)
    assert float(total.numpy()[0]) == pytest.approx(float(p_np.sum()), rel=1e-12)

    mean = wp.array(np.array([p_np.mean()], np.float64), dtype=wp.float64, device=dev)
    acc = wp.zeros(4, dtype=wp.float64, device=dev)
    ctr = N // 2
    wp.launch(_observe_kernel, dim=(N, N, N), inputs=[p, patch, mean, wp.int32(ctr), acc], device=dev)
    s2, patch_sum, patch_n, centre = (float(v) for v in acc.numpy())

    d = p_np - p_np.mean()
    assert np.sqrt(s2) == pytest.approx(float(np.linalg.norm(d)), rel=1e-11)
    assert patch_n == pytest.approx(float(patch_np.sum()), rel=1e-12)
    assert patch_sum / patch_n == pytest.approx(float(d[patch_np == 1].mean()), rel=1e-11)
    assert centre == pytest.approx(float(d[ctr, ctr, ctr]), rel=1e-12)


def test_the_mean_subtraction_is_not_a_no_op() -> None:
    """The control for the test above: with a non-zero mean, subtracting it must CHANGE the answer.

    Without this, a reduction that silently ignored `mean` would pass the comparison whenever the field
    happened to be centred — and the fields this run uses are not centred.
    """
    p_np = _field()
    assert abs(p_np.mean()) > 0.1
    assert np.linalg.norm(p_np) != pytest.approx(np.linalg.norm(p_np - p_np.mean()), rel=1e-3)

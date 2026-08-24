"""Oracle A (runtime) — the solid->fluid dilatation source is physically correct on the GPU.

CUDA-gated companion to ``test_fsi_on_audit.py`` (which certifies the wiring + setpoint statically on CPU).
This drives the real ``SolidDilatationCoupling`` on a Warp FieldGrid and checks the sign/magnitude sanity the
``fsi_coupling`` docstring claims but that had no dedicated test:

  * NULL — a rigid translation of the whole node cloud produces ``div(v_s) == 0`` (no spurious pore source).
  * CONTRACTION — a uniform contraction ``v = -k (x - c)`` has the exact analytic divergence ``-3k``; the
    reconstructed ``div_vs`` must reproduce it in the cloud interior, so ``-alpha*div(v_s) = +3k alpha`` is a
    POSITIVE pressure source (squeezing the skeleton pressurises the pore fluid).

Both use a linear velocity field, which the Peskin spread + central-difference divergence reproduce exactly
in the interior, so the acceptance bands are tight (not fudge factors). Runs on the gbook A5000.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.components.fluid.field_grid import FieldGrid
from aleph.components.incumbent.fsi_coupling import SolidDilatationCoupling

_CUDA_DEVICE = next((d for d in wp.get_devices() if d.is_cuda), None)
pytestmark = pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")

_N = 24          # cells per axis
_DX = 0.5        # um
_NODE_LO, _NODE_HI, _NODE_STEP = 3.0, 9.0, 0.5   # dense node lattice, well inside the grid


def _node_lattice() -> np.ndarray:
    ax = np.arange(_NODE_LO, _NODE_HI + 1e-9, _NODE_STEP)
    xs, ys, zs = np.meshgrid(ax, ax, ax, indexing="ij")
    return np.stack([xs.ravel(), ys.ravel(), zs.ravel()], axis=1).astype(np.float64)


def _run_div_vs(velocity_np: np.ndarray, pos_np: np.ndarray) -> np.ndarray:
    """Spread the given nodal solid velocity to the grid and return the reconstructed div_vs field."""
    dev = str(_CUDA_DEVICE)
    grid = FieldGrid((_N, _N, _N), _DX, origin=(0.0, 0.0, 0.0), device=dev)
    coupling = SolidDilatationCoupling(grid)
    n = pos_np.shape[0]
    with wp.ScopedDevice(dev):
        pos = wp.array(pos_np, dtype=wp.vec3d, device=dev)
        vs = wp.array(velocity_np, dtype=wp.vec3d, device=dev)
        vol = wp.array(np.full(n, _NODE_STEP**3, np.float64), dtype=wp.float64, device=dev)
        active = wp.zeros(n, dtype=wp.int32, device=dev)   # 0 => load-bearing (>=0)
        coupling.update(pos, vs, vol, active)
    return coupling.div_vs_host()


def _interior(field: np.ndarray) -> np.ndarray:
    """Cells strictly inside the node cloud (avoid the cloud-edge and grid-edge one-sided stencils)."""
    lo = int(np.ceil((_NODE_LO + 1.0) / _DX))
    hi = int(np.floor((_NODE_HI - 1.0) / _DX))
    return field[lo:hi, lo:hi, lo:hi]


def test_rigid_translation_makes_no_pore_source() -> None:
    """NULL: a uniform translation has zero divergence -> no spurious fluid source."""
    pos = _node_lattice()
    vel = np.tile(np.array([0.13, -0.07, 0.21]), (pos.shape[0], 1))
    div = _interior(_run_div_vs(vel, pos))
    assert np.max(np.abs(div)) < 1e-9, "rigid translation must not create a pore-pressure source"


def test_uniform_contraction_is_a_positive_pressure_source() -> None:
    """CONTRACTION: v = -k(x-c) has div = -3k exactly; div_vs must reproduce it (negative)."""
    k = 0.05
    pos = _node_lattice()
    centroid = pos.mean(axis=0)
    vel = -k * (pos - centroid)
    div = _interior(_run_div_vs(vel, pos))
    assert np.mean(div) == pytest.approx(-3.0 * k, rel=0.05), "interior div_vs must equal the analytic -3k"
    assert np.max(div) < 0.0, "a contracting skeleton must give div(v_s) < 0 -> -alpha*div > 0 (source)"

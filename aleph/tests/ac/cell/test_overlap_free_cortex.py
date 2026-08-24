"""ac/ overlap-free cortex integration — build_cell(overlap_free_cortex=True) starts with zero steric force.

Certifies the CellConfig.overlap_free_cortex flag, threaded down build_cell -> weave_cell -> ff.weave.weave
-> build_cortex_network, removes the ~64k build-interpenetration steric artifact (2634 pN) from the assembled
production cell. CUDA-gated (build_cell needs the Warp device); runs on the gbook A5000.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

_CUDA_DEVICE = next((d for d in wp.get_devices() if d.is_cuda), None)


def _steric_max(cell) -> tuple[float, int]:
    dev = cell.device
    with wp.ScopedDevice(dev):
        f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=dev)
        cell.steric.accumulate(cell.state, f)
        mag = np.linalg.norm(f.numpy(), axis=1)
    return float(mag.max()), int((mag > 1e-9).sum())


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_overlap_free_cortex_removes_steric_artifact() -> None:
    from aleph.components.incumbent.assemble import CellConfig, build_cell
    common = dict(n_filaments=20000, with_membrane=False, with_nucleus=False,
                  with_pressure=False, with_myosin=False)
    off_max, off_n = _steric_max(build_cell(CellConfig(overlap_free_cortex=False, **common)))
    on_max, on_n = _steric_max(build_cell(CellConfig(overlap_free_cortex=True, **common)))
    assert off_max > 100.0 and off_n > 0, "the zero-thickness shell must show the build interpenetration"
    assert on_max < 1e-6 and on_n == 0, f"overlap-free cortex must start with zero steric force (got {on_max})"

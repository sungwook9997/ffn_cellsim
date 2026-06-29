"""CPU-path bit-parity: DcmTurgorForceGPU vs the live DcmTurgorForce.

On a CUDA-less machine only the CPU dispatch path of
``cell/dcm_gpu_forces.py::DcmTurgorForceGPU`` is exercisable. This test asserts
that path is BIT-IDENTICAL to the original ``cell/dcm.py::DcmTurgorForce`` (the
deliverable's hard gate, max abs force diff < 1e-12 N).

Both turgor forces are constructed with the SAME parameters (faces / face_cell /
n_cells / V0 / turgor_dP0 / K_vol) on the same DCM spheroid, each computes its
force array over the same snapshot, and the two arrays are compared. The GPU
path itself (cupy / gpu_local) is structurally constructed but only
gbook-A5000-validated — there is no CUDA GPU here.

The CPU twin routes through ``gpu_opt.kernels_cpu.mesh_pressure_forces`` (the K1
kernel) while the reference computes the same divergence-theorem volume +
face-normal pressure inline; the two are the same FP ops in the same order, so
the result is bit-identical (not merely close).

Run:  ~/miniconda3/envs/ffn_sim/bin/python -m pytest \
          ffn_sim/tests/test_dcm_turgor_gpu_parity.py -q
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.archive.hoomd_legacy.cell.dcm import (
    DcmTurgorForce,
    ResolvedDCM,
    build_dcm_simulation,
)
from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_forces import DcmTurgorForceGPU, on_gpu


def _force_of(force) -> np.ndarray:
    """Force array (N,3) of a Custom force after its set_forces ran (CPU path)."""
    with force.cpu_local_force_arrays as arr:
        return np.asarray(arr.force).copy()


def _build(n_cells=3):
    """Build a small DCM spheroid, attach a twin DcmTurgorForce + DcmTurgorForceGPU.

    Returns (sim, turgor_ref, turgor_gpu, handles). The reference turgor is the
    builder-wired ``DcmTurgorForce``; we additionally append a ``DcmTurgorForceGPU``
    constructed from the SAME faces/face_cell/V0 so both run over the identical
    configuration each step.
    """
    p = ResolvedDCM(subdivisions=1, dt=1.0e-9, seed=7)
    h = build_dcm_simulation(p, n_cells)
    sim = h["sim"]
    ig = sim.operations.integrator
    ref = h["turgor"]
    assert isinstance(ref, DcmTurgorForce)

    gpu = DcmTurgorForceGPU(
        faces=ref.faces, face_cell=ref.face_cell, n_cells=ref.n_cells,
        V0=ref.V0, turgor_dP0=ref.turgor_dP0, K_vol=ref.K_vol)
    ig.forces.append(gpu)
    sim.run(0)
    return sim, ref, gpu, h, p


def test_cpu_path_is_cpu_not_gpu():
    """Guard: on this machine the dispatch must pick the CPU path (no CUDA)."""
    sim, _, gpu, _, _ = _build(n_cells=2)
    assert not on_gpu(sim)
    assert gpu._dispatch().gpu is False
    assert gpu._dispatch().xp is np


@pytest.mark.parametrize("n_cells", [1, 3, 5])
def test_turgor_gpu_cpu_path_bit_parity(n_cells):
    """DcmTurgorForceGPU (CPU path) == DcmTurgorForce, diff < 1e-12 N."""
    sim, ref, gpu, _, _ = _build(n_cells=n_cells)

    # Run a little so the shells deform off the perfect sphere (non-trivial volume
    # deviation → non-trivial pressure), then recompute both at the SAME snapshot.
    sim.run(200)
    ref.set_forces(sim.timestep)
    gpu.set_forces(sim.timestep)
    F_ref = _force_of(ref)
    F_gpu = _force_of(gpu)

    assert np.all(np.isfinite(F_ref))
    assert np.all(np.isfinite(F_gpu))
    # the turgor must actually be engaged (else a trivial all-zero pass)
    assert np.abs(F_ref).max() > 0.0, "no turgor force engaged — test is trivial"

    max_abs = float(np.abs(F_ref - F_gpu).max())
    assert max_abs < 1e-12, f"turgor CPU-path diff {max_abs:.3e} N > 1e-12"


def test_turgor_gpu_parity_at_t0():
    """Parity also holds at t=0 (perfect-sphere baseline pressure)."""
    sim, ref, gpu, _, _ = _build(n_cells=4)
    ref.set_forces(0)
    gpu.set_forces(0)
    F_ref = _force_of(ref)
    F_gpu = _force_of(gpu)
    max_abs = float(np.abs(F_ref - F_gpu).max())
    assert max_abs < 1e-12, f"turgor t=0 diff {max_abs:.3e} N > 1e-12"

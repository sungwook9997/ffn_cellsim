"""CPU-path bit-parity: DcmActiveRimTractionGPUVec vs DcmActiveRimTractionGPU.

The vectorized active rim traction (``cell/dcm_gpu_forces.py::
DcmActiveRimTractionGPUVec``) removes the per-cell Python loop that the gbook
A5000 profile showed is the dominant per-step cost (~22 ms/step, ~57% of total at
N=200). This test asserts the vectorized force is BIT-IDENTICAL (max abs diff
< 1e-10 N) to the parent ``DcmActiveRimTractionGPU`` on the CPU dispatch path, AND
that the rim diagnostics match — so the vectorization is a pure speedup, not a
physics change. (The parent is itself CPU-bit-parity to the live ActiveRimTraction
via ``test_dcm_active_gpu_parity.py``, so this chains parity to production.)

Built on the GPU-FRIENDLY build (uniform contiguous per-cell node blocks, which the
vectorized class requires), with the active rim traction engaged.

Run:  ~/miniconda3/envs/ffn_sim/bin/python -m pytest \
          ffn_sim/tests/test_dcm_active_vec_parity.py -q
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_build import ResolvedGpuDCM, build_gpu_dcm_simulation
from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_forces import (
    DcmActiveRimTractionGPU,
    DcmActiveRimTractionGPUVec,
)


def _force_of(force) -> np.ndarray:
    with force.cpu_local_force_arrays as arr:
        return np.asarray(arr.force).copy()


def _build(n_cells):
    """Active GPU-friendly spheroid; append a parent + a vec twin sharing state."""
    p = ResolvedGpuDCM(subdivisions=1, seed=7)
    h = build_gpu_dcm_simulation(p, n_cells, active=True)
    sim = h["sim"]
    ig = sim.operations.integrator
    parent = h["traction"]
    assert isinstance(parent, DcmActiveRimTractionGPU)

    common = dict(
        cell_of_node=h["cell_of_node"], ranges=h["ranges"], active=parent.active,
        int_mult=parent.int_mult, R_cell=p.R_cell, z0=p.z_substrate,
        f_act=parent.f_act, f_cap=parent.f_cap, ramp_steps=parent.ramp_steps,
        contact_band=parent.contact_band,
        neighbour_factor=parent.r_neigh / p.R_cell,
        max_neighbours=parent.max_neigh,
        integrin_switch_gain=parent.switch_gain, belt_factor=parent.belt)
    vec = DcmActiveRimTractionGPUVec(**common)
    ig.forces.append(vec)
    sim.run(0)
    return sim, parent, vec, h, p


@pytest.mark.parametrize("n_cells", [12, 40])
def test_active_vec_bit_parity(n_cells):
    """Vectorized active traction == parent (loop) traction, diff < 1e-10 N."""
    sim, parent, vec, _, _ = _build(n_cells)
    # past the ramp → full magnitude, engaged forces.
    sim.run(6000)
    parent.set_forces(sim.timestep)
    vec.set_forces(sim.timestep)
    F_par = _force_of(parent)
    F_vec = _force_of(vec)

    assert np.all(np.isfinite(F_par)) and np.all(np.isfinite(F_vec))
    assert np.abs(F_par).max() > 0.0, "no rim traction engaged — trivial test"

    max_abs = float(np.abs(F_par - F_vec).max())
    assert max_abs < 1e-10, f"vec active-traction diff {max_abs:.3e} N > 1e-10"

    # rim diagnostics must match (order-independent)
    assert np.array_equal(np.sort(parent.rim_cells), np.sort(vec.rim_cells))
    # f_per_cell keys (cells that contributed) must match
    assert set(parent.f_per_cell) == set(vec.f_per_cell)
    for c in parent.f_per_cell:
        assert abs(parent.f_per_cell[c] - vec.f_per_cell[c]) < 1e-18


def test_active_vec_switched_integrin_parity():
    """int_mult (junction-switch integrin gain) branch is also bit-identical."""
    sim, parent, vec, h, p = _build(40)
    parent.int_mult[0] = 3.0
    parent.int_mult[5] = 3.0
    sim.run(6000)
    parent.set_forces(sim.timestep)
    vec.set_forces(sim.timestep)
    F_par = _force_of(parent)
    F_vec = _force_of(vec)
    max_abs = float(np.abs(F_par - F_vec).max())
    assert max_abs < 1e-10, f"vec int_mult diff {max_abs:.3e} N > 1e-10"


def test_active_vec_belt_off_parity():
    """belt_factor=0 (no apical contraction belt) is also bit-identical."""
    p = ResolvedGpuDCM(subdivisions=1, seed=7)
    h = build_gpu_dcm_simulation(p, 24, active=True)
    sim = h["sim"]
    ig = sim.operations.integrator
    parent = h["traction"]
    common = dict(
        cell_of_node=h["cell_of_node"], ranges=h["ranges"], active=parent.active,
        int_mult=parent.int_mult, R_cell=p.R_cell, z0=p.z_substrate,
        f_act=parent.f_act, f_cap=parent.f_cap, ramp_steps=parent.ramp_steps,
        contact_band=parent.contact_band,
        neighbour_factor=parent.r_neigh / p.R_cell,
        max_neighbours=parent.max_neigh,
        integrin_switch_gain=parent.switch_gain, belt_factor=0.0)
    # rebuild parent with belt=0 too for a fair compare
    parent0 = DcmActiveRimTractionGPU(**common)
    vec0 = DcmActiveRimTractionGPUVec(**common)
    ig.forces.append(parent0)
    ig.forces.append(vec0)
    sim.run(6000)
    parent0.set_forces(sim.timestep)
    vec0.set_forces(sim.timestep)
    F_par = _force_of(parent0)
    F_vec = _force_of(vec0)
    max_abs = float(np.abs(F_par - F_vec).max())
    assert max_abs < 1e-10, f"vec belt-off diff {max_abs:.3e} N > 1e-10"

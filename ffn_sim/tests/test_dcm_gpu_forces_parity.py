"""CPU-path bit-parity: DcmTentContactGPU vs the live DcmTentContact.

On a CUDA-less machine only the CPU dispatch path of ``cell/dcm_gpu_forces.py``
is exercisable. This test asserts that path is BIT-IDENTICAL to the existing
``cell/dcm_contact.py::DcmTentContact`` (the deliverable's hard gate, < 1e-12 N).

Both forces are attached to the SAME small native-mesh + tent spheroid, each
computes its force array over the same snapshot, and the two arrays are compared.
The GPU path itself (cupy / gpu_local) is structurally constructed but only
gbook-A5000-validated — there is no CUDA GPU here.

Run:  ~/miniconda3/envs/ffn_sim/bin/python -m pytest ffn_sim/tests/test_dcm_gpu_forces_parity.py -q
"""

from __future__ import annotations

import numpy as np
import pytest

import hoomd

from ffn_sim.archive.hoomd_legacy.cell.dcm_contact import DcmTentContact
from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_forces import DcmTentContactGPU, on_gpu
from ffn_sim.archive.hoomd_legacy.cell.dcm_native_shell import (
    ResolvedNativeDCM,
    build_native_dcm_simulation,
)


def _force_of(force) -> np.ndarray:
    """Force array (N,3) of a Custom force after its set_forces ran (CPU path)."""
    with force.cpu_local_force_arrays as arr:
        return np.asarray(arr.force).copy()


def _build_pair_sim(n_cells=4, *, c_adh=2.0e-6, adh=1.0e8, rep=1.0e8):
    """Build a small native-mesh DCM sim and attach a twin DcmTentContactGPU.

    Returns (sim, tent_ref, tent_gpu). The reference tent is the one wired by
    build_native_dcm_simulation; the GPU twin is built with the SAME parameters
    and appended to the integrator so HOOMD calls its set_forces too.
    """
    p = ResolvedNativeDCM(subdivisions=2, cluster="3d", spacing_factor=2.15,
                          adh_strength=adh, rep_strength=rep, c_adh=c_adh,
                          dt=3e-10)
    h = build_native_dcm_simulation(p, n_cells, substrate=False, contact=True)
    sim, tent_ref = h["sim"], h["tent"]

    nv = h["nv"]
    patch_area = 4.0 * np.pi * p.R_cell ** 2 / nv
    tent_gpu = DcmTentContactGPU(
        cell_of_node=h["cell_of_tag"], r_contact=1.05 * h["mean_edge"],
        c_adh=p.c_adh, rep_strength=p.rep_strength,
        adh_strength=p.adh_strength, patch_area=patch_area,
        force_cap=p.force_cap)
    sim.operations.integrator.forces.append(tent_gpu)
    sim.run(0)  # attach + first force eval
    return sim, tent_ref, tent_gpu, p


def test_cpu_path_is_cpu_not_gpu():
    """Guard: on this machine the dispatch must pick the CPU path (no CUDA)."""
    sim, _, tent_gpu, _ = _build_pair_sim(n_cells=2)
    assert not on_gpu(sim)
    assert tent_gpu._dispatch().gpu is False
    # array module is numpy, kernel module is kernels_cpu on the CPU path
    import numpy as _np
    from ffn_sim.archive.hoomd_legacy.gpu_opt import kernels_cpu
    assert tent_gpu._dispatch().xp is _np
    assert tent_gpu._dispatch().kernels is kernels_cpu


@pytest.mark.parametrize("n_cells", [2, 4, 8])
def test_tent_contact_gpu_cpu_path_bit_parity(n_cells):
    """DcmTentContactGPU (CPU path) == DcmTentContact, max abs diff < 1e-12 N."""
    sim, tent_ref, tent_gpu, _ = _build_pair_sim(n_cells=n_cells)

    # Run a few steps so cells move into genuine contact (non-trivial forces),
    # then recompute both at the SAME snapshot.
    sim.run(800)
    tent_ref.set_forces(sim.timestep)
    tent_gpu.set_forces(sim.timestep)
    F_ref = _force_of(tent_ref)
    F_gpu = _force_of(tent_gpu)

    assert np.all(np.isfinite(F_ref))
    assert np.all(np.isfinite(F_gpu))
    # the contact must actually be engaged (else a trivial all-zero pass)
    assert np.abs(F_ref).max() > 0.0, "no contact engaged — test is trivial"

    max_abs = float(np.abs(F_ref - F_gpu).max())
    assert max_abs < 1e-12, f"CPU-path force diff {max_abs:.3e} N exceeds 1e-12"


def test_cad_mult_path_bit_parity():
    """The junction-switch cad_mult branch is also bit-identical on the CPU path."""
    p = ResolvedNativeDCM(subdivisions=2, cluster="2d", spacing_factor=2.05,
                          adh_strength=1.0e8, rep_strength=1.0e8, c_adh=2.0e-6,
                          dt=3e-10)
    h = build_native_dcm_simulation(p, 3, substrate=False, contact=True)
    sim = h["sim"]
    nv = h["nv"]
    patch_area = 4.0 * np.pi * p.R_cell ** 2 / nv
    cad = np.array([0.3, 1.0, 2.5])  # per-cell adhesion multipliers

    ref = DcmTentContact(
        cell_of_node=h["cell_of_tag"], r_contact=1.05 * h["mean_edge"],
        c_adh=p.c_adh, rep_strength=p.rep_strength, adh_strength=p.adh_strength,
        patch_area=patch_area, force_cap=p.force_cap, cad_mult=cad)
    gpu = DcmTentContactGPU(
        cell_of_node=h["cell_of_tag"], r_contact=1.05 * h["mean_edge"],
        c_adh=p.c_adh, rep_strength=p.rep_strength, adh_strength=p.adh_strength,
        patch_area=patch_area, force_cap=p.force_cap, cad_mult=cad)
    sim.operations.integrator.forces.append(ref)
    sim.operations.integrator.forces.append(gpu)
    sim.run(0)
    sim.run(800)
    ref.set_forces(sim.timestep)
    gpu.set_forces(sim.timestep)
    F_ref = _force_of(ref)
    F_gpu = _force_of(gpu)
    max_abs = float(np.abs(F_ref - F_gpu).max())
    assert max_abs < 1e-12, f"cad_mult CPU-path diff {max_abs:.3e} N exceeds 1e-12"

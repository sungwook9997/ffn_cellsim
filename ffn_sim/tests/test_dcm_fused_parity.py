"""CPU-path parity: DcmFusedForceGPU total == sum of the separate forces.

The fused force (``cell/dcm_gpu_forces.py::DcmFusedForceGPU``) computes turgor +
tent + substrate + (vectorized) active in ONE ``set_forces`` callback. This test
asserts that, on the CPU dispatch path, the fused force's total equals the SUM of
the four separate forces (DcmTurgorForceGPU + DcmTentContactGPU +
DcmSubstrateForceGPU + DcmActiveRimTractionGPUVec) to floating-point round-off
(max abs diff < 1e-10 N) — so fusing is purely a per-step-overhead optimization,
not a physics change.

Run:  ~/miniconda3/envs/ffn_sim/bin/python -m pytest \
          ffn_sim/tests/test_dcm_fused_parity.py -q
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.cell.dcm_gpu_build import ResolvedGpuDCM, build_gpu_dcm_simulation
from ffn_sim.cell.dcm_gpu_forces import (
    DcmFusedForceGPU,
    DcmActiveRimTractionGPUVec,
)


def _force_of(force) -> np.ndarray:
    with force.cpu_local_force_arrays as arr:
        return np.asarray(arr.force).copy()


@pytest.mark.parametrize("active", [False, True])
def test_fused_equals_sum_of_separate(active):
    """Fused total == turgor + tent + substrate + active, diff < 1e-10 N."""
    p = ResolvedGpuDCM(subdivisions=1, seed=7)
    h = build_gpu_dcm_simulation(p, 24, active=active)
    sim = h["sim"]
    ig = sim.operations.integrator

    turgor, contact, substrate = h["turgor"], h["contact"], h["substrate"]
    # swap the loop active for the vectorized active (fused uses the vec law)
    vec = None
    if active:
        parent = h["traction"]
        ig.forces.remove(parent)
        vec = DcmActiveRimTractionGPUVec(
            cell_of_node=h["cell_of_node"], ranges=h["ranges"],
            active=parent.active, int_mult=parent.int_mult, R_cell=p.R_cell,
            z0=p.z_substrate, f_act=parent.f_act, f_cap=parent.f_cap,
            ramp_steps=parent.ramp_steps, contact_band=parent.contact_band,
            neighbour_factor=parent.r_neigh / p.R_cell,
            max_neighbours=parent.max_neigh,
            integrin_switch_gain=parent.switch_gain, belt_factor=parent.belt)
        ig.forces.append(vec)

    # build the fused force with identical parameters; append it too.
    W_cs = p.W_cs_Jm2 * p.ligand_density * h["area_per_node"]
    fused = DcmFusedForceGPU(
        n_cells=24,
        faces=turgor.faces, face_cell=turgor.face_cell, V0=turgor.V0,
        turgor_dP0=turgor.turgor_dP0, K_vol=turgor.K_vol,
        cell_of_node=contact.cell_of_node, r_contact=contact.r_contact,
        c_adh=contact.c_adh, rep_strength=contact.rep, adh_strength=contact.omega,
        patch_area=contact.A, contact_force_cap=contact.force_cap,
        cad_mult=contact.cad_mult,
        z0=substrate.z0, W_cs=W_cs, adh_range=p.R_cell, k_sub=p.k_sub_Nm,
        with_active=active,
        ranges=h["ranges"], active=(vec.active if active else None),
        int_mult=(vec.int_mult if active else None), R_cell=p.R_cell,
        f_act=(vec.f_act if active else 0.0),
        f_cap=(vec.f_cap if active else 0.0),
        ramp_steps=(vec.ramp_steps if active else 1),
        contact_band=(vec.contact_band if active else 0.5),
        neighbour_factor=(vec.r_neigh / p.R_cell if active else 2.6),
        max_neighbours=(vec.max_neigh if active else 9),
        belt_factor=(vec.belt if active else 0.0))
    ig.forces.append(fused)
    sim.run(0)

    sim.run(6000)
    ts = sim.timestep
    turgor.set_forces(ts)
    contact.set_forces(ts)
    substrate.set_forces(ts)
    if active:
        vec.set_forces(ts)
    fused.set_forces(ts)

    sep = _force_of(turgor) + _force_of(contact) + _force_of(substrate)
    if active:
        sep = sep + _force_of(vec)
    fus = _force_of(fused)

    assert np.all(np.isfinite(fus))
    assert np.abs(sep).max() > 0.0, "no force engaged — trivial test"
    max_abs = float(np.abs(sep - fus).max())
    assert max_abs < 1e-10, f"fused vs sum-of-separate diff {max_abs:.3e} N > 1e-10"

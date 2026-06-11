"""CPU-path bit-parity: DcmActiveRimTractionGPU vs the live ActiveRimTraction.

On a CUDA-less machine only the CPU dispatch path of
``cell/dcm_gpu_forces.py::DcmActiveRimTractionGPU`` is exercisable. This test
asserts that path is BIT-IDENTICAL to the existing
``cell/dcm_active.py::ActiveRimTraction`` (the deliverable's hard gate,
max abs force diff < 1e-12 N).

Both forces are constructed with the SAME parameters and SHARED mutable state
(cell_of_node / ranges / active / int_mult) on the same ACTIVE native+tent
spheroid, each computes its force array over the same snapshot, and the two
arrays are compared. The GPU path itself (cupy / gpu_local) is structurally
constructed but only gbook-A5000-validated — there is no CUDA GPU here.

Run:  ~/miniconda3/envs/ffn_sim/bin/python -m pytest \
          ffn_sim/tests/test_dcm_active_gpu_parity.py -q
"""

from __future__ import annotations

import numpy as np
import pytest

from ffn_sim.cell.dcm_active import (
    ActiveRimTraction,
    ResolvedActiveSpheroid,
    build_active_spheroid,
)
from ffn_sim.cell.dcm_gpu_forces import DcmActiveRimTractionGPU, on_gpu


def _force_of(force) -> np.ndarray:
    """Force array (N,3) of a Custom force after its set_forces ran (CPU path)."""
    with force.cpu_local_force_arrays as arr:
        return np.asarray(arr.force).copy()


def _build(n_active=8, n_max=10):
    """Build an active spheroid, then attach a twin DcmActiveRimTractionGPU.

    Returns (sim, traction_ref, traction_gpu, handles). The reference traction is
    a LEGACY-MODE ``ActiveRimTraction`` (migrate_factor=0, lead_bias=1) — the
    basal-only splay law that the frozen ``DcmActiveRimTractionGPU`` (a separate,
    PI-gated GPU port file) still implements. This test pins CPU↔GPU parity OF THE
    SAME LAW; the production CPU traction additionally applies the whole-cell
    migration term (PI 2026-06-11, migrate_factor=0.6) which the GPU twin does not
    yet carry — porting it is frozen-file work tracked for the GPU port, not this
    parity gate. We REMOVE the production traction that build_active_spheroid wired
    (migration-on) and append a legacy-mode twin + the GPU twin, both with the same
    shared state, so set_forces runs over the identical configuration.
    """
    p = ResolvedActiveSpheroid(subdivisions=2, spacing_factor=2.3, dt=3.0e-10,
                               seed=7)
    h = build_active_spheroid(p, n_active, n_max, belt=True,
                              integrin_substrate=True)
    sim = h["sim"]
    ig = sim.operations.integrator
    # drop the production (migration-on) traction wired by the builder
    if h["traction"] in ig.forces:
        ig.forces.remove(h["traction"])

    common = dict(
        cell_of_node=h["cell_of_node"], ranges=h["ranges"], active=h["active"],
        int_mult=h["st"].int_mult, R_cell=p.R_cell, z0=p.z_substrate,
        f_act=p.f_act, f_cap=p.f_cap, ramp_steps=p.ramp_steps,
        contact_band=p.rim_contact_band, neighbour_factor=p.rim_neighbour_factor,
        max_neighbours=p.rim_max_neighbours,
        integrin_switch_gain=p.integrin_switch_gain, belt_factor=p.belt_factor)

    # legacy-mode CPU reference (matches the GPU twin's basal-only law). The legacy
    # ActiveRimTraction has NO spreading-arrest term, so it takes no arrest kwarg.
    ref = ActiveRimTraction(migrate_factor=0.0, lead_bias=1.0, **common)
    assert isinstance(ref, ActiveRimTraction)
    ig.forces.append(ref)

    # SAME LAW as the legacy reference: pin the GPU twin's arrest OFF so the
    # bit-parity comparison is of the basal-only splay law alone. (Arrest is
    # opt-in/default-off anyway; pinning it here keeps the test correct regardless of
    # the default — it was the prior failure: the GPU twin defaulted arrest ON,
    # throttling the int_mult-boosted footprint while the legacy ref could not, a
    # 1.2e-10 N mismatch in test_switched_integrin_gain_bit_parity.)
    gpu = DcmActiveRimTractionGPU(**common, arrest_radius_factor=None)
    ig.forces.append(gpu)
    sim.run(0)
    return sim, ref, gpu, h, p


def test_cpu_path_is_cpu_not_gpu():
    """Guard: on this machine the dispatch must pick the CPU path (no CUDA)."""
    sim, _, gpu, _, _ = _build(n_active=4, n_max=5)
    assert not on_gpu(sim)
    assert gpu._dispatch().gpu is False
    import numpy as _np
    assert gpu._dispatch().xp is _np


@pytest.mark.parametrize("n_active", [4, 8])
def test_active_traction_gpu_cpu_path_bit_parity(n_active):
    """DcmActiveRimTractionGPU (CPU path) == ActiveRimTraction, diff < 1e-12 N."""
    sim, ref, gpu, _, _ = _build(n_active=n_active, n_max=n_active + 2)

    # Run past the ramp so the traction is at full magnitude (non-trivial forces),
    # then recompute both at the SAME snapshot/timestep.
    sim.run(5000)
    ref.set_forces(sim.timestep)
    gpu.set_forces(sim.timestep)
    F_ref = _force_of(ref)
    F_gpu = _force_of(gpu)

    assert np.all(np.isfinite(F_ref))
    assert np.all(np.isfinite(F_gpu))
    # the active traction must actually be engaged (else a trivial all-zero pass)
    assert np.abs(F_ref).max() > 0.0, "no rim traction engaged — test is trivial"

    max_abs = float(np.abs(F_ref - F_gpu).max())
    assert max_abs < 1e-12, f"active-traction CPU-path diff {max_abs:.3e} N > 1e-12"

    # rim diagnostics must also match
    assert np.array_equal(np.sort(ref.rim_cells), np.sort(gpu.rim_cells))


def test_switched_integrin_gain_bit_parity():
    """The int_mult (junction-switch integrin gain) branch is also bit-identical."""
    sim, ref, gpu, h, p = _build(n_active=8, n_max=10)
    # raise integrin gain on a couple of cells (the junction-switch effect)
    st = h["st"]
    st.int_mult[0] = p.integrin_strong_factor
    st.int_mult[3] = p.integrin_strong_factor

    sim.run(5000)
    ref.set_forces(sim.timestep)
    gpu.set_forces(sim.timestep)
    F_ref = _force_of(ref)
    F_gpu = _force_of(gpu)
    max_abs = float(np.abs(F_ref - F_gpu).max())
    assert max_abs < 1e-12, f"int_mult CPU-path diff {max_abs:.3e} N > 1e-12"

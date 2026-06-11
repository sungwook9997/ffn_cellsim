"""Spreading-arrest law on the active rim traction — behaviour + parity.

The arrest law (``cell/dcm_gpu_forces.py``: ``arrest_radius_factor`` on
``DcmActiveRimTractionGPU`` / ``…GPUVec`` / ``DcmFusedForceGPU``) smoothly switches
the outward rim traction OFF as a rim cell's radial spread approaches a physiological
cap ``r_max = arrest_radius_factor · R0_cluster``. These tests assert:

  1. ``arrest_radius_factor=None`` (arrest OFF) reproduces the legacy traction EXACTLY
     (the no-arrest blow-up path is preserved bit-for-bit for the A/B comparison).
  2. With arrest ON, once cells have spread the per-rim traction magnitude is
     STRICTLY ≤ the arrest-OFF magnitude (and < it for the outer rim) — the arrest
     genuinely reduces the outward push.
  3. The vectorized + fused arrest paths are bit-parity to the loop arrest path.
  4. The arrest gain is the documented smooth clip law (analytic check of the helper).

Run:  ~/miniconda3/envs/ffn_sim/bin/python -m pytest \
          ffn_sim/tests/test_dcm_spreading_arrest.py -q
"""

from __future__ import annotations

import numpy as np

from ffn_sim.cell.dcm_gpu_build import ResolvedGpuDCM, build_gpu_dcm_simulation
from ffn_sim.cell.dcm_gpu_forces import (
    DcmActiveRimTractionGPU,
    DcmActiveRimTractionGPUVec,
    DcmFusedForceGPU,
)


def _force_of(force) -> np.ndarray:
    with force.cpu_local_force_arrays as arr:
        return np.asarray(arr.force).copy()


def _dispersal_p():
    """NON-cohesive (dispersal-prone) bands — the regime the spreading-arrest mechanism
    targets. With the cohesion fix now the default (c_adh=5µm etc.), a production spheroid
    does NOT over-spread (cells stay adhered), so the arrest's "the basal footprint grows
    past the cap" premise only holds in the LEGACY non-cohesive regime. Building the arrest
    tests on these legacy bands (a) restores that premise so the arrest actually engages and
    (b) keeps them fast — cohesive contact is ~4× heavier per step (cells actually touch).
    The arrest law itself is band-agnostic; this only controls whether the rim over-spreads.
    """
    return ResolvedGpuDCM(subdivisions=1, seed=7, c_adh=5.0e-7, adh_strength=1.0e8,
                          rep_strength=1.0e8, spacing_factor=2.3)


def _common(parent, p, h, *, arrest_radius_factor):
    return dict(
        cell_of_node=h["cell_of_node"], ranges=h["ranges"], active=parent.active,
        int_mult=parent.int_mult, R_cell=p.R_cell, z0=p.z_substrate,
        f_act=parent.f_act, f_cap=parent.f_cap, ramp_steps=parent.ramp_steps,
        contact_band=parent.contact_band,
        neighbour_factor=parent.r_neigh / p.R_cell,
        max_neighbours=parent.max_neigh,
        integrin_switch_gain=parent.switch_gain, belt_factor=parent.belt,
        arrest_radius_factor=arrest_radius_factor)


def test_arrest_off_is_legacy_bit_identical():
    """arrest_radius_factor=None == the legacy traction (no arrest), diff = 0."""
    p = _dispersal_p()
    # build WITHOUT arrest so the in-build traction is the legacy path
    h = build_gpu_dcm_simulation(p, 30, active=True, arrest=False)
    sim = h["sim"]
    ig = sim.operations.integrator
    legacy = h["traction"]
    assert legacy.arrest_factor is None

    # a second legacy force constructed independently must match the build's.
    twin = DcmActiveRimTractionGPU(**_common(legacy, p, h,
                                             arrest_radius_factor=None))
    ig.forces.append(twin)
    sim.run(6000)
    legacy.set_forces(sim.timestep)
    twin.set_forces(sim.timestep)
    F0 = _force_of(legacy)
    F1 = _force_of(twin)
    assert np.abs(F0).max() > 0.0
    assert float(np.abs(F0 - F1).max()) == 0.0


def test_arrest_reduces_traction():
    """With arrest ON, the per-rim traction is ≤ the arrest-OFF traction.

    Build a spreading spheroid, run it past the ramp so the rim has spread, then
    evaluate an arrest-ON and an arrest-OFF traction on the SAME state. The arrest-ON
    rim-force magnitude must be ≤ arrest-OFF everywhere and STRICTLY smaller for at
    least the outer rim cells (whose radial distance has grown toward r_max).
    """
    p = _dispersal_p()
    h = build_gpu_dcm_simulation(p, 30, active=True, arrest=True,
                                 arrest_radius_factor=1.2)  # tight cap → clear effect
    sim = h["sim"]
    ig = sim.operations.integrator
    on = h["traction"]
    off = DcmActiveRimTractionGPUVec(**_common(on, p, h,
                                               arrest_radius_factor=None))
    ig.forces.append(off)
    # spread well past the ramp so the rim's radial distance grows past (1-w)*r_max.
    sim.run(20000)
    on.set_forces(sim.timestep)
    off.set_forces(sim.timestep)
    F_on = np.linalg.norm(_force_of(on), axis=1)
    F_off = np.linalg.norm(_force_of(off), axis=1)

    assert on._R0_cluster is not None and on._R0_cluster > 0
    # arrest never INCREASES the force (per node)
    assert np.all(F_on <= F_off + 1e-18)
    # and the arrest mean gain is below 1 (the rim is being throttled)
    assert on.arrest_gain_mean < 1.0
    # net outward push is reduced
    assert F_on.sum() < F_off.sum()


def test_arrest_vec_loop_fused_parity():
    """Loop / vec / fused arrest paths are bit-parity (< 1e-10 N) with arrest ON."""
    p = _dispersal_p()
    h = build_gpu_dcm_simulation(p, 24, active=True, arrest=True,
                                 arrest_radius_factor=1.5)
    sim = h["sim"]
    ig = sim.operations.integrator
    vec = h["traction"]                      # build wires the vec by default
    loop = DcmActiveRimTractionGPU(**_common(vec, p, h, arrest_radius_factor=1.5))
    fused = DcmFusedForceGPU(
        n_cells=h["n_cells"], faces=h["faces"], face_cell=h["face_cell"],
        V0=h["V0"], turgor_dP0=p.turgor_dP0, K_vol=p.K_vol,
        cell_of_node=h["cell_of_node"], r_contact=p.r_contact_factor * h["mean_edge"],
        c_adh=p.c_adh, rep_strength=p.rep_strength, adh_strength=p.adh_strength,
        patch_area=h["area_per_node"], contact_force_cap=p.contact_force_cap,
        cad_mult=h["cad_mult"], z0=p.z_substrate, W_cs=0.0, adh_range=p.R_cell,
        with_active=True, ranges=h["ranges"], active=vec.active,
        int_mult=vec.int_mult, R_cell=p.R_cell, f_act=vec.f_act, f_cap=vec.f_cap,
        ramp_steps=vec.ramp_steps, contact_band=vec.contact_band,
        neighbour_factor=vec.r_neigh / p.R_cell, max_neighbours=vec.max_neigh,
        belt_factor=vec.belt, arrest_radius_factor=1.5)
    ig.forces.append(loop)
    sim.run(20000)
    vec.set_forces(sim.timestep)
    loop.set_forces(sim.timestep)
    F_vec = _force_of(vec)
    F_loop = _force_of(loop)
    assert np.abs(F_vec).max() > 0.0
    assert float(np.abs(F_vec - F_loop).max()) < 1e-10

    # fused active term equals the vec force (on the same captured R0). The fused
    # force captures its own R0 from its own first call → drive it on the same state.
    d = vec._dispatch()
    pos_g = None
    with d.snapshot() as snap:
        tag = np.asarray(snap.particles.tag)
        pos = np.asarray(snap.particles.position, dtype=np.float64)
        pos_g = pos[np.argsort(tag)]
    F_fused_act = np.asarray(fused._active_term(d, pos_g, sim.timestep))
    # vec F in global (tag) order for comparison
    with d.snapshot() as snap:
        tag = np.asarray(snap.particles.tag)
    perm = np.argsort(tag)
    F_vec_g = F_vec[perm]
    assert float(np.abs(F_vec_g - F_fused_act).max()) < 1e-10


def test_arrest_gain_helper_law():
    """The arrest gain helper is the documented smooth clip g=clip((r_max-r)/(w·r_max),0,1)."""
    p = _dispersal_p()
    h = build_gpu_dcm_simulation(p, 12, active=True, arrest=True,
                                 arrest_radius_factor=2.0)
    f = h["traction"]
    f._R0_cluster = 10.0e-6          # set a known R0
    f.arrest_factor = 2.0
    f.arrest_width = 0.4
    r_max = 2.0 * 10.0e-6
    rn = np.array([0.0, (1 - 0.4) * r_max, r_max - 0.2 * 0.4 * r_max,
                   r_max, 1.5 * r_max])
    g = np.asarray(f._arrest_gain(np, rn))
    expect = np.clip((r_max - rn) / (0.4 * r_max), 0.0, 1.0)
    assert np.allclose(g, expect)
    assert g[0] == 1.0          # deep inside → full traction
    assert g[1] == 1.0          # at the edge of the band → still full
    assert 0.0 < g[2] < 1.0     # inside the band → partial
    assert g[3] == 0.0          # at r_max → zero
    assert g[4] == 0.0          # beyond r_max → clipped zero

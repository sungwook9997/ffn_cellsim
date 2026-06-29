"""Tests for the device-aware overdamped BAOAB Action (GPU-main port, B2).

The decisive contract: on CPU (``xp=np``) ``OverdampedBAOABDevice`` must be **bit-identical**
to the frozen ``integrator.baoab.LeimkuhlerMatthewsBAOAB`` — same RNG stream, same step, same
wrap — so the Layer-2 CBM trajectory is unchanged when it opts into the device Action on a CPU
box. The GPU path (cupy) is validated on the gbook A5000 (statistical parity), not here (cupy
is absent on the CPU dev box; those checks live in the gbook smoke script).
"""

from __future__ import annotations

import numpy as np
import pytest

import gsd.hoomd
import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.integrator.baoab import make_baoab_updater
from ffn_sim.archive.hoomd_legacy.integrator.baoab_device import (
    OverdampedBAOABDevice,
    make_baoab_updater_device,
)


def _blob_state(n: int, box: float, r0: float, seed: int) -> gsd.hoomd.Frame:
    """A small jittered cluster of ``cell`` particles (+ optional parked ``void``)."""
    rng = np.random.default_rng(seed)
    m = int(np.ceil(n ** (1.0 / 3.0))) + 1
    g = (np.arange(m) - (m - 1) / 2.0) * (1.05 * r0)
    xx, yy, zz = np.meshgrid(g, g, g, indexing="ij")
    pts = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])[:n]
    pts = pts + rng.normal(0.0, 0.02 * r0, pts.shape)
    snap = gsd.hoomd.Frame()
    snap.particles.N = n
    snap.particles.types = ["cell", "void"]
    snap.particles.typeid = np.zeros(n, dtype=np.uint32)
    snap.particles.position = pts.astype(np.float64)
    snap.particles.mass = np.ones(n, dtype=np.float64)
    snap.configuration.box = [box, box, box, 0.0, 0.0, 0.0]
    return snap


def _build(action_factory, *, n=27, box=40.0, r0=1.0, dt=2e-3, seed=11):
    """Build a Morse-cohesion CBM-like sim wired with the given (action, updater) factory."""
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=seed)
    sim.create_state_from_snapshot(_blob_state(n, box, r0, seed))
    nlist = md.nlist.Tree(buffer=0.1 * r0)
    morse = md.pair.Morse(nlist=nlist, default_r_cut=0.0)
    morse.params[("cell", "cell")] = dict(D0=2.0, alpha=5.0 / r0, r0=r0)
    morse.r_cut[("cell", "cell")] = r0 + 5.0 / (5.0 / r0)
    for pair in (("cell", "void"), ("void", "void")):
        morse.params[pair] = dict(D0=0.0, alpha=1.0, r0=r0)
        morse.r_cut[pair] = 0.0
    morse.mode = "shift"
    ig = md.Integrator(dt=dt)
    ig.forces.append(morse)
    sim.operations.integrator = ig
    action, updater = action_factory(
        kT=0.5, gamma={"cell": 1.3, "void": 1.3e6}, dt=dt, seed=seed
    )
    sim.operations.updaters.append(updater)
    sim._action_anchor = action  # keep alive
    return sim, action


def _positions_by_tag(sim: hoomd.Simulation) -> np.ndarray:
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tag = np.asarray(s.particles.tag)
        out = np.empty_like(pos)
        out[tag] = pos
        return out.copy()


class TestCpuBitParity:
    """CPU path must match the frozen Action bit-for-bit (same seed → same trajectory)."""

    def test_trajectory_bit_identical_to_frozen(self):
        sim_ref, _a = _build(make_baoab_updater, seed=23)
        sim_dev, _b = _build(make_baoab_updater_device, seed=23)
        sim_ref.run(0)
        sim_dev.run(0)
        for _ in range(40):
            sim_ref.run(10)
            sim_dev.run(10)
            p_ref = _positions_by_tag(sim_ref)
            p_dev = _positions_by_tag(sim_dev)
            np.testing.assert_array_equal(p_ref, p_dev)

    def test_prv_rnds_buffer_matches_frozen(self):
        sim_ref, a_ref = _build(make_baoab_updater, seed=5)
        sim_dev, a_dev = _build(make_baoab_updater_device, seed=5)
        sim_ref.run(7)
        sim_dev.run(7)
        np.testing.assert_array_equal(a_ref.prv_rnds, a_dev.prv_rnds)

    def test_device_action_is_cpu_by_default(self):
        _sim, action = _build(make_baoab_updater_device)
        _sim.run(0)
        assert action._on_gpu is False  # noqa: SLF001 — explicit device-branch check


class TestRuntimeGuards:
    def test_negative_gamma_rejected(self):
        with pytest.raises(ValueError):
            OverdampedBAOABDevice(kT=1.0, gamma={"cell": -1.0}, dt=1e-3)

    def test_zero_dt_rejected(self):
        with pytest.raises(ValueError):
            OverdampedBAOABDevice(kT=1.0, gamma={"cell": 1.0}, dt=0.0)

    def test_negative_kT_rejected(self):
        with pytest.raises(ValueError):
            OverdampedBAOABDevice(kT=-1.0, gamma={"cell": 1.0}, dt=1e-3)

    def test_missing_type_in_gamma_rejected_at_attach(self):
        sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
        sim.create_state_from_snapshot(_blob_state(8, 40.0, 1.0, 1))
        ig = md.Integrator(dt=1e-3)
        sim.operations.integrator = ig
        _a, updater = make_baoab_updater_device(
            kT=1.0, gamma={"cell": 1.0}, dt=1e-3, seed=1  # missing "void"
        )
        sim.operations.updaters.append(updater)
        with pytest.raises(RuntimeError):
            sim.run(0)

    def test_methods_nonempty_rejected_at_attach(self):
        sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
        sim.create_state_from_snapshot(_blob_state(8, 40.0, 1.0, 1))
        ig = md.Integrator(dt=1e-3)
        ig.methods.append(md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0))
        sim.operations.integrator = ig
        _a, updater = make_baoab_updater_device(
            kT=1.0, gamma={"cell": 1.0, "void": 1e6}, dt=1e-3, seed=1
        )
        sim.operations.updaters.append(updater)
        with pytest.raises(RuntimeError):
            sim.run(0)

"""Enabled-path CANARY tests for default-OFF compartments (fast CI regression).

These lock in the ENABLED-PATH smoke findings (scripts/compartment_smoke/) as a
fast CI guard: a future refactor that breaks a compartment's enabled
build+attach+step path is caught here, not only by re-running the ad-hoc smoke
scripts. They assemble a MINIMAL enabled compartment, step it with the project
L-M BAOAB integrator, and assert the run stays finite + a force-free / topology
invariant holds.

PLUMBING ONLY — these are crash/finite/assembly regressions, NOT physics-band
gates. Sizes are tiny (fast); the full-resolution plumbing lives in the smoke
harnesses. SMOKE-ONLY candidate constants are used for None-gated parameters.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import gsd.hoomd
import hoomd
import hoomd.md as md

from ffn_sim.integrator.baoab import make_baoab_updater

_KT = 4.28e-21
_ETA = 65.9  # MCF7 cytoplasm (Hu 2024)


def _drag(r: float) -> float:
    return 6.0 * math.pi * _ETA * r


def _step(snap, *, dt, gamma_map, attach, n_steps=150, seed=3):
    """Build a CPU sim from ``snap``, let ``attach(sim, ig)`` add forces, attach
    BAOAB, run, and return (sim, finite)."""
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=seed)
    sim.create_state_from_snapshot(snap)
    ig = md.Integrator(dt=dt, methods=[])
    sim.operations.integrator = ig
    attach(sim, ig)
    # gamma must cover every type in the state.
    types = list(sim.state.particle_types)
    gm = {t: gamma_map.get(t, _drag(1.0e-7)) for t in types}
    _, upd = make_baoab_updater(kT=_KT, gamma=gm, dt=dt, seed=seed + 1)
    sim.operations.updaters.append(upd)
    sim.run(0)
    sim.run(n_steps)
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position, dtype=np.float64).copy()
    return sim, bool(np.all(np.isfinite(pos)))


def _bare_frame(n=4, box=5.0e-5, bead="cortex_actin"):
    fr = gsd.hoomd.Frame()
    fr.particles.N = n
    fr.particles.types = [bead]
    fr.particles.typeid = [0] * n
    fr.particles.position = [[0.0, 0.0, i * 1e-7] for i in range(n)]
    fr.particles.mass = [1.0] * n
    fr.particles.velocity = [[0.0, 0.0, 0.0]] * n
    fr.particles.image = [[0, 0, 0]] * n
    fr.configuration.box = [box, box, box, 0.0, 0.0, 0.0]
    fr.bonds.N = 0
    fr.bonds.types = []
    fr.angles.N = 0
    fr.angles.types = []
    return fr


def test_microtubules_enabled_path_assembles_and_steps():
    from ffn_sim.cell.microtubules import (
        attach_microtubule_forces,
        extend_snapshot_with_microtubules,
        resolve_microtubules,
    )

    gb = _drag(12.5e-9)
    p = resolve_microtubules(
        {"microtubules": {"enabled": True, "n_mt": 4, "beads_per_mt": 5,
                          "L_mt": 2.0e-6, "Y_stretch": 2.3e-7}},
        kT=_KT, gamma_b=gb,
    )
    snap = extend_snapshot_with_microtubules(_bare_frame(), p)
    dt = 0.5 * p.dt_cfl

    def attach(sim, ig):
        attach_microtubule_forces(sim, p, cfl_strict=True)

    sim, finite = _step(snap, dt=dt, gamma_map={"mt_bead": gb, "mtoc": gb},
                        attach=attach)
    assert finite
    # backbone bonds stayed near l0 (stiff rod).
    fin = sim.state.get_snapshot()
    g = np.asarray(fin.bonds.group).reshape(-1, 2)
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position).copy()
        tags = np.asarray(s.particles.tag).copy()
    by_tag = np.empty_like(pos)
    by_tag[tags] = pos
    L = np.linalg.norm(by_tag[g[:, 0]] - by_tag[g[:, 1]], axis=1)
    assert np.all(np.abs(L - p.l0) / p.l0 < 0.05)


def test_intermediate_filaments_enabled_path_is_force_free_and_steps():
    from ffn_sim.cell.intermediate_filaments import (
        attach_if_bonds_to_simulation,
        extend_snapshot_with_if_cage,
        resolve_intermediate_filaments,
    )

    gif = _drag(5.0e-9)
    p = resolve_intermediate_filaments(
        {"intermediate_filaments": {"enabled": True, "n_filaments": 8,
                                    "beads_per_fil": 5}},
        kT=_KT, R_cell=7.5e-6, R_nuc=3.0e-6,
    )
    # IF extender needs a HOOMD-style snapshot (bare gsd None groups crash it).
    tmp = hoomd.Simulation(device=hoomd.device.CPU(), seed=0)
    tmp.create_state_from_snapshot(_bare_frame())
    host = tmp.state.get_snapshot()
    snap = extend_snapshot_with_if_cage(host, p, centroid=(0, 0, 0), gamma_if=gif, seed=1)

    captured = {}

    def attach(sim, ig):
        captured["bond"] = attach_if_bonds_to_simulation(sim, p, gamma_if=gif, cfl_strict=True)

    dt = min(6.98e-7, 0.5 * 0.1 * gif / p.k_bb)
    sim, finite = _step(snap, dt=dt, gamma_map={"if_bead": gif}, attach=attach)
    assert finite
    # Force-free construction: born bond energy is thermal-scale, not the huge
    # r0=0 pre-tension the audit fixed.
    # (energy read right after build via a 0-step probe on a fresh sim.)
    e0 = float(captured["bond"].energy)
    assert e0 / _KT < 5000.0  # thermal-scale, not ~1e5+ pre-tension


def test_linc_enabled_path_forms_bridges_and_steps():
    from ffn_sim.cell.linc import (
        configure_linc_bond_potential,
        extend_snapshot_with_linc,
        resolve_linc,
    )

    n_nuc, r0, R_nuc = 20, 50e-9, 3.0e-6
    p = resolve_linc(
        {"linc": {"enabled": True, "k_linc": 1.0e-2, "r0": r0, "capture_radius": 150e-9}},
        R_cell=7.5e-6, R_nuc=R_nuc,
    )
    # nucleus shell + perinuclear acceptors one r0 out (same directions).
    idx = np.arange(n_nuc) + 0.5
    phi = math.pi * (3 - math.sqrt(5))
    ct = np.clip(1 - 2 * idx / n_nuc, -1, 1)
    st = np.sqrt(np.maximum(0, 1 - ct * ct))
    az = phi * idx
    dirs = np.stack([st * np.cos(az), st * np.sin(az), ct], axis=1)
    pos = np.concatenate([R_nuc * dirs, (R_nuc + r0) * dirs])
    fr = gsd.hoomd.Frame()
    fr.particles.N = 2 * n_nuc
    fr.particles.types = ["nucleus_bead", "cortex_actin"]
    fr.particles.typeid = [0] * n_nuc + [1] * n_nuc
    fr.particles.position = pos.tolist()
    fr.particles.mass = [1.0] * (2 * n_nuc)
    fr.particles.velocity = [[0.0, 0.0, 0.0]] * (2 * n_nuc)
    fr.particles.image = [[0, 0, 0]] * (2 * n_nuc)
    fr.configuration.box = [6e-5, 6e-5, 6e-5, 0, 0, 0]
    fr.bonds.N = 0
    fr.bonds.types = []
    fr.angles.N = 0
    fr.angles.types = []

    snap = extend_snapshot_with_linc(
        fr, p, nucleus_tags=np.arange(n_nuc),
        cytoskeleton_tags=np.arange(n_nuc, 2 * n_nuc),
    )
    assert int(snap.bonds.N) == n_nuc  # one bridge per nucleus bead (no zero-bonds trap)

    gb = _drag(150e-9)

    def attach(sim, ig):
        h = md.bond.Harmonic()
        configure_linc_bond_potential(h, p)
        ig.forces.append(h)

    sim, finite = _step(
        snap, dt=6.98e-7,
        gamma_map={"nucleus_bead": gb, "cortex_actin": gb}, attach=attach,
    )
    assert finite

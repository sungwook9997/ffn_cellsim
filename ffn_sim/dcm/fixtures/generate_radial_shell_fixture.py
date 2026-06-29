"""Generate + commit HOOMD reference fixtures for the B2 radial-shell parity test.

Runs the three REAL production compartment ``md.force.Custom`` forces
(``NucleusConfinement``, ``MembraneSurfaceTension``, ``EnclosedVolumePressure``)
inside a HOOMD CPU sim on the SAME configs the native parity script uses
(``scripts/h7_native_radial_force_parity.py``), and commits the per-bead force +
energy as the bit-parity TARGET for ``tests/warp_port/test_radial_shell_warp_parity.py``.
The Warp kernel is graded against THIS committed numpy-force output (guard-rail 2).

Run (writes three .npz fixtures next to this file)::

    python ffn_sim/warp_port/fixtures/generate_radial_shell_fixture.py
"""

from __future__ import annotations

import os

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.cell.nucleus import resolve_nucleus, NucleusConfinement
from ffn_sim.archive.hoomd_legacy.cell.membrane_surface import resolve_membrane_surface, MembraneSurfaceTension
from ffn_sim.archive.hoomd_legacy.cortex.enclosed_volume import resolve_enclosed_volume, EnclosedVolumePressure

HERE = os.path.dirname(os.path.abspath(__file__))

# Match the native parity script exactly (known-correct construction).
N = 2048
R_CELL = 7.5e-6
R_NUC = 3.0e-6
SEED = 3


def _cloud(n, R, spread, seed):
    rng = np.random.default_rng(seed)
    r = R + (rng.random(n) - 0.5) * 2.0 * spread
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    return (r[:, None] * u).astype(np.float64)


def _build_sim(pos):
    n = pos.shape[0]
    snap = hoomd.Snapshot()
    L = 60.0e-6
    snap.configuration.box = [L, L, L, 0, 0, 0]
    snap.particles.N = n
    snap.particles.types = ["shell"]
    snap.particles.position[:] = pos
    snap.particles.typeid[:] = np.zeros(n, dtype=np.int32)
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.tuners.clear()  # keep rows tag-ordered (row i == tag i)
    sim.operations.integrator = md.Integrator(dt=1e-6)
    return sim


def _force_of(sim, force):
    sim.operations.integrator.forces = [force]
    sim.run(0)
    return (np.asarray(force.forces, np.float64).copy(),
            np.asarray(force.energies, np.float64).copy())


def main() -> None:
    rng_range = (0, N)

    p_nuc = resolve_nucleus(
        {"nucleus": {"E_nuc": 5e3, "ratio_lamin": 3.0, "knee_strain": 0.10,
                     "critical_pore_area_um2": 7.0}},
        R_nuc=R_NUC, n_beads=N)
    p_mem = resolve_membrane_surface(
        {"membrane_surface": {"gamma_mem": 1.0e-4}}, R_cell=R_CELL)
    p_ev = resolve_enclosed_volume({}, R_cell=R_CELL)

    cl_nuc = _cloud(N, R_NUC, 3.0 * p_nuc.d_knee, SEED)
    cl_shell = _cloud(N, R_CELL, 0.05 * R_CELL, SEED)

    cases = [
        ("nucleus", 0, NucleusConfinement(p_nuc, rng_range), cl_nuc,
         dict(R0=p_nuc.R_nuc, pa=p_nuc.k_chrom, pb=p_nuc.k_lamin,
              pc=p_nuc.d_knee, pd=p_nuc.F_knee)),
        ("membrane", 1, MembraneSurfaceTension(p_mem, rng_range), cl_shell,
         dict(R0=0.0, pa=p_mem.gamma_mem, pb=p_mem.K_A, pc=p_mem.A0, pd=0.0)),
        ("turgor", 2, EnclosedVolumePressure(p_ev, rng_range), cl_shell,
         dict(R0=0.0, pa=p_ev.turgor_dP0, pb=p_ev.K_vol, pc=p_ev.V0, pd=0.0)),
    ]

    for name, law, force, pos, prm in cases:
        sim = _build_sim(pos)
        F, U = _force_of(sim, force)
        out = os.path.join(HERE, f"radial_ref_{name}.npz")
        np.savez(
            out,
            name=name, law=law, N=N, t0=rng_range[0], t1=rng_range[1],
            pos=pos, tag=np.arange(N, dtype=np.uint32),
            R0=float(prm["R0"]), pa=float(prm["pa"]), pb=float(prm["pb"]),
            pc=float(prm["pc"]), pd=float(prm["pd"]),
            ref_force=F, ref_energy=U,
        )
        print(f"[{name:8s}] law={law} wrote {out}  "
              f"max|F|={np.abs(F).max():.3e}  max|U|={np.abs(U).max():.3e}  "
              f"params(R0={prm['R0']:.3e},pa={prm['pa']:.3e},pb={prm['pb']:.3e},"
              f"pc={prm['pc']:.3e},pd={prm['pd']:.3e})")


if __name__ == "__main__":
    main()

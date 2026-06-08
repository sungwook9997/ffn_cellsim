"""ENABLED-PATH smoke harness — microtubule aster (PLUMBING ONLY).

Assembles a default-OFF microtubule aster (1 MTOC + n_mt stiff chains) onto a
tiny inert cortex stub, steps it with the project L-M BAOAB integrator at the
physiological cytoplasm drag, and confirms the enabled build+attach+step path
runs without crashing and stays rod-like. This is the FIRST real execution of
``microtubules`` enabled (the suite only does sim.run(0) static force checks).

NOT a physics claim. Constants that have no ratified value use the
``COMPARTMENT_ACTIVATION_DETAILED_PLANS_2026-06-09`` PI-candidate values, used
SMOKE-ONLY to make the plumbing run (the production enabled path still raises /
is PI-gated). No band pass/fail.

Run:  python ffn_sim/scripts/compartment_smoke/microtubules.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
for _p in (_HERE.parents[3], _HERE.parents[0]):   # repo root + this dir
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import _smoke_common as sc
from ffn_sim.cell.microtubules import (
    GAMMA_DENYLIST_PREFIX,
    attach_microtubule_forces,
    extend_snapshot_with_microtubules,
    resolve_microtubules,
)

# --------------------------------------------------------------------------
# SMOKE-ONLY candidate constants (COMPARTMENT_ACTIVATION_DETAILED_PLANS_2026-06-09)
# --------------------------------------------------------------------------
N_MT = 20            # PI-candidate ORDER_ESTIMATE (module default); SMOKE-ONLY
BEADS_PER_MT = 25    # discretisation knob (grid-invariant); SMOKE-ONLY
L_MT = 5.0e-6        # m  PI-candidate ORDER_ESTIMATE (3-5 um aster); SMOKE-ONLY
Y_STRETCH = 2.3e-7   # N  PI-candidate ORDER_ESTIMATE (Kis2002/Pampaloni2006); SMOKE-ONLY
R_MTBEAD = 12.5e-9   # m  MT tube radius (Gittes geometry; DERIVED, cited)
# EI / L_p use the module defaults (Gittes 1993 SOLID_LITERATURE — not pending).
N_STEPS = 3000


def run() -> dict:
    gamma_b = sc.stokes_drag(R_MTBEAD)        # cytoplasm drag at the MT bead radius
    p = resolve_microtubules(
        {"microtubules": {
            "enabled": True, "n_mt": N_MT, "beads_per_mt": BEADS_PER_MT,
            "L_mt": L_MT, "Y_stretch": Y_STRETCH,
        }},
        kT=sc.KT, gamma_b=gamma_b,
    )

    host = sc.cortex_stub_frame(n_beads=12, r_shell=1.0e-6, box_l=6.0e-5)
    n_host = int(host.particles.N)
    snap = extend_snapshot_with_microtubules(host, p)

    sim = sc.make_sim(snap, seed=7)
    dt = 0.5 * p.dt_cfl                        # safely below the (stretch) CFL
    sc.set_integrator(sim, dt)
    backbone, bending = attach_microtubule_forces(sim, p, cfl_strict=True)
    gamma_map = sc.auto_gamma_map(
        sim, radius_by_type={"mt_bead": R_MTBEAD, "mtoc": R_MTBEAD}
    )
    sc.attach_baoab(sim, dt=dt, gamma_map=gamma_map, seed=8)

    pos0 = sc.positions_by_tag(sim)
    sim.run(0)                                 # BAOAB attach + first force eval
    sim.run(N_STEPS)
    pos1 = sc.positions_by_tag(sim)
    finite = sc.all_finite(sim)

    # ---- observables (plumbing sanity, NOT a band) ----
    final = sim.state.get_snapshot()
    groups = np.asarray(final.bonds.group, dtype=np.int64).reshape(-1, 2)
    btypeid = np.asarray(final.bonds.typeid, dtype=np.int64).reshape(-1)
    btypes = list(final.bonds.types)
    mt_mask = np.array(
        [btypes[t].startswith(GAMMA_DENYLIST_PREFIX) for t in btypeid], dtype=bool
    )
    mt_groups = groups[mt_mask]
    lengths = np.linalg.norm(pos1[mt_groups[:, 0]] - pos1[mt_groups[:, 1]], axis=1)
    rel_dev = (lengths - p.l0) / p.l0

    mtoc_tag = n_host                          # MTOC is the first appended particle
    mt_tags = np.arange(n_host, pos1.shape[0])
    bead_tags = mt_tags[mt_tags != mtoc_tag]
    dist_from_mtoc = np.linalg.norm(pos1[bead_tags] - pos1[mtoc_tag], axis=1)
    max_disp = float(np.max(np.linalg.norm(pos1 - pos0, axis=1)))

    result = {
        "compartment": "microtubules",
        "verdict": "PLUMBING_OK" if finite else "PLUMBING_FAIL",
        "physics_claim": False,
        "n_steps": N_STEPS,
        "dt_s": dt,
        "params": {
            "n_mt": p.n_mt, "beads_per_mt": p.beads_per_mt, "l0_m": p.l0,
            "EI_Nm2": p.EI, "L_p_m": p.L_p, "Y_stretch_N": p.Y_stretch,
            "k_backbone_Npm": p.k_backbone, "k_angle_Nm": p.k_angle,
            "gamma_b_Nspm": gamma_b,
            "dt_cfl_stretch_s": p.dt_cfl_stretch, "dt_cfl_bend_s": p.dt_cfl_bend,
        },
        "topology": {
            "n_host_stub": n_host, "n_mt_particles": int(p.n_beads_total),
            "n_mt_backbone_bonds": int(mt_groups.shape[0]),
        },
        "observables": {
            "backbone_len_mean_m": float(np.mean(lengths)),
            "backbone_len_min_m": float(np.min(lengths)),
            "backbone_len_max_m": float(np.max(lengths)),
            "backbone_rel_dev_rms": float(np.sqrt(np.mean(rel_dev**2))),
            "frac_bonds_within_20pct_l0": float(np.mean(np.abs(rel_dev) < 0.20)),
            "max_bead_dist_from_mtoc_m": float(np.max(dist_from_mtoc)),
            "L_mt_m": L_MT,
            "max_particle_disp_m": max_disp,
            "all_finite": finite,
        },
        "constants_smoke_only": [
            "n_mt=20 (ORDER_ESTIMATE)", "L_mt=5e-6 (ORDER_ESTIMATE)",
            "Y_stretch=2.3e-7 (ORDER_ESTIMATE)",
        ],
        "note": (
            "PLUMBING smoke: confirms the microtubules enabled path "
            "assembles+steps without crashing and stays rod-like. NOT the "
            "PRIMARY L_p-recovery activation gate."
        ),
    }

    _figure(p, pos0, pos1, mtoc_tag, bead_tags, lengths, result)
    return result


def _figure(p, pos0, pos1, mtoc_tag, bead_tags, lengths, result) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.6))

    # (a) backbone bond-length distribution vs rest length l0
    ax0.hist(lengths * 1e9, bins=30, color="#4477aa", alpha=0.85)
    ax0.axvline(p.l0 * 1e9, color="crimson", lw=2,
                label=f"rest l0 = {p.l0*1e9:.1f} nm")
    ax0.set_xlabel("mt_backbone bond length [nm]")
    ax0.set_ylabel("count")
    ax0.set_title("MT backbone lengths (rod stability)")
    ax0.legend(fontsize=8)

    # (b) aster x-z projection, initial (grey) vs final (blue), MTOC marked
    ax1.scatter(pos0[bead_tags, 0] * 1e6, pos0[bead_tags, 2] * 1e6,
                s=6, color="0.7", label="initial")
    ax1.scatter(pos1[bead_tags, 0] * 1e6, pos1[bead_tags, 2] * 1e6,
                s=6, color="#1f77b4", label="final")
    ax1.scatter([pos1[mtoc_tag, 0] * 1e6], [pos1[mtoc_tag, 2] * 1e6],
                s=80, marker="*", color="crimson", label="MTOC", zorder=5)
    ax1.set_xlabel("x [um]")
    ax1.set_ylabel("z [um]")
    ax1.set_title(f"MT aster (n_mt={p.n_mt}, {result['n_steps']} steps)")
    ax1.set_aspect("equal", adjustable="datalim")
    ax1.legend(fontsize=8)

    fig.suptitle(
        f"microtubules ENABLED-PATH smoke — {result['verdict']} "
        f"(dt={result['dt_s']:.2e}s, SMOKE-ONLY constants)", fontsize=11
    )
    sc.save_fig(fig, "smoke_microtubules")
    plt.close(fig)


def main() -> int:
    res = run()
    path = sc.save_json("smoke_microtubules", res)
    obs = res["observables"]
    print(f"[microtubules smoke] {res['verdict']}  "
          f"backbone len mean={obs['backbone_len_mean_m']*1e9:.2f}nm "
          f"(l0={res['params']['l0_m']*1e9:.2f}nm, rms dev="
          f"{obs['backbone_rel_dev_rms']*100:.2f}%), "
          f"max_disp={obs['max_particle_disp_m']*1e9:.2f}nm, "
          f"dt={res['dt_s']:.2e}s")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PLUMBING_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())

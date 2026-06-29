"""MULTI-COMPARTMENT coexistence smoke — MT aster + IF cage (PLUMBING ONLY).

The per-compartment smokes each used their OWN md.bond.Harmonic. The real cell
assembles MANY compartments on ONE shared bonded force (HOOMD 7.0.1 demands params
for EVERY bond type on EVERY Harmonic — two competing Harmonics crash; this is the
class of the audit's #6/#12 standalone-Harmonic bugs). This harness builds two
internal compartments together — a microtubule aster (mt_backbone bond +
mt_bending angle) AND an intermediate-filament cage (if_backbone + per-r0-bin
if_crosslink) — registers ALL bond types on ONE shared md.bond.Harmonic + the
angle on ONE md.angle.Harmonic + ONE gamma_map covering every type, and steps with
BAOAB. Confirms the multi-compartment integration pattern assembles + steps
without crashing and both compartments stay intact.

NOT a physics claim. SMOKE-ONLY constants as in the per-compartment harnesses.

Run:  python ffn_sim/scripts/compartment_smoke/coexistence.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
for _p in (_HERE.parents[3], _HERE.parents[0]):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import hoomd.md as md

import _smoke_common as sc
from ffn_sim.archive.hoomd_legacy.cell.intermediate_filaments import (
    extend_snapshot_with_if_cage,
    register_if_bond_params,
    resolve_intermediate_filaments,
)
from ffn_sim.archive.hoomd_legacy.cell.microtubules import (
    extend_snapshot_with_microtubules,
    resolve_microtubules,
)

R_MTBEAD, R_IFBEAD = 12.5e-9, 5.0e-9
N_STEPS = 2500


def run() -> dict:
    gb_mt, gif = sc.stokes_drag(R_MTBEAD), sc.stokes_drag(R_IFBEAD)
    p_mt = resolve_microtubules(
        {"microtubules": {"enabled": True, "n_mt": 8, "beads_per_mt": 10,
                          "L_mt": 4.0e-6, "Y_stretch": 2.3e-7}},
        kT=sc.KT, gamma_b=gb_mt,
    )
    p_if = resolve_intermediate_filaments(
        {"intermediate_filaments": {"enabled": True, "n_filaments": 30,
                                    "beads_per_fil": 6}},
        kT=sc.KT, R_cell=7.5e-6, R_nuc=3.0e-6,
    )

    # Chain the two extenders onto one snapshot (MT first, then IF).
    host = sc.to_hoomd_snapshot(
        sc.cortex_stub_frame(n_beads=8, r_shell=1.0e-6, box_l=6.0e-5)
    )
    n_host = int(host.particles.N)
    snap = extend_snapshot_with_microtubules(host, p_mt)
    snap = extend_snapshot_with_if_cage(
        snap, p_if, centroid=(0.0, 0.0, 0.0), gamma_if=gif, seed=2
    )

    sim = sc.make_sim(snap, seed=71)
    dt = 0.5 * min(p_mt.dt_cfl, 0.1 * gif / p_if.k_bb)
    ig = sc.set_integrator(sim, dt)

    # ONE shared bond Harmonic covering BOTH families (the integration pattern).
    shared_bond = md.bond.Harmonic()
    register_if_bond_params(shared_bond, p_if)                 # if_backbone + if_crosslink_b*
    shared_bond.params["mt_backbone"] = dict(k=p_mt.k_backbone, r0=p_mt.l0)
    ig.forces.append(shared_bond)
    # ONE angle Harmonic for the MT bending term.
    shared_angle = md.angle.Harmonic()
    shared_angle.params["mt_bending"] = dict(k=p_mt.k_angle, t0=p_mt.angle_t0)
    ig.forces.append(shared_angle)

    gamma_map = sc.auto_gamma_map(
        sim, radius_by_type={"mtoc": R_MTBEAD, "mt_bead": R_MTBEAD,
                             "if_bead": R_IFBEAD},
    )
    sc.attach_baoab(sim, dt=dt, gamma_map=gamma_map, seed=72)

    pos0 = sc.positions_by_tag(sim)
    sim.run(0)
    sim.run(N_STEPS)
    pos1 = sc.positions_by_tag(sim)
    finite = sc.all_finite(sim)

    # ---- observables ----
    final = sim.state.get_snapshot()
    btypes = list(final.bonds.types)
    groups = np.asarray(final.bonds.group, dtype=np.int64).reshape(-1, 2)
    btypeid = np.asarray(final.bonds.typeid, dtype=np.int64).reshape(-1)
    ptypes = list(final.particles.types)

    def _len(name_pred):
        m = np.array([name_pred(btypes[t]) for t in btypeid], bool)
        g = groups[m]
        return np.linalg.norm(pos1[g[:, 0]] - pos1[g[:, 1]], axis=1) if g.size else np.array([])

    mt_len = _len(lambda s: s == "mt_backbone")
    if_len = _len(lambda s: s == "if_backbone")
    n_bond_force = sum(1 for f in ig.forces if isinstance(f, md.bond.Harmonic))

    both_present = all(t in ptypes for t in ("mt_bead", "mtoc", "if_bead"))
    ok = (finite and both_present and n_bond_force == 1 and mt_len.size > 0
          and if_len.size > 0)
    result = {
        "compartment": "coexistence(microtubules+intermediate_filaments)",
        "verdict": "COEXIST_OK" if ok else "COEXIST_FAIL",
        "physics_claim": False,
        "n_steps": N_STEPS,
        "dt_s": dt,
        "topology": {
            "n_host_stub": n_host,
            "n_bond_types_on_shared_force": len(btypes),
            "bond_types": btypes,
            "n_bond_force_computes": n_bond_force,
            "particle_types": ptypes,
        },
        "observables": {
            "mt_backbone_rel_dev_rms": float(np.sqrt(np.mean(((mt_len - p_mt.l0) / p_mt.l0) ** 2))),
            "if_backbone_rel_dev_rms": float(np.sqrt(np.mean(((if_len - p_if.l_seg) / p_if.l_seg) ** 2))),
            "both_compartments_present": both_present,
            "single_shared_bond_force": n_bond_force == 1,
            "max_particle_disp_m": float(np.max(np.linalg.norm(pos1 - pos0, axis=1))),
            "all_finite": finite,
        },
        "note": (
            "PLUMBING coexistence smoke: MT (mt_backbone + mt_bending) and IF "
            "(if_backbone + if_crosslink_b*) on ONE shared md.bond.Harmonic + ONE "
            "md.angle.Harmonic + ONE gamma_map, stepped together. Validates the "
            "single-shared-force integration pattern (NOT per-compartment "
            "standalone Harmonics, which would crash on each other's types)."
        ),
    }

    _figure(pos1, n_host, p_mt, p_if, mt_len, if_len, result)
    return result


def _figure(pos1, n_host, p_mt, p_if, mt_len, if_len, result) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.6))

    ax0.hist(mt_len * 1e9, bins=25, color="#1f77b4", alpha=0.7, label="mt_backbone")
    ax0.hist(if_len * 1e9, bins=25, color="#aa6644", alpha=0.7, label="if_backbone")
    ax0.axvline(p_mt.l0 * 1e9, color="#1f77b4", ls="--", lw=1.2)
    ax0.axvline(p_if.l_seg * 1e9, color="#aa6644", ls="--", lw=1.2)
    ax0.set_xlabel("backbone bond length [nm]")
    ax0.set_ylabel("count")
    ax0.set_title("Both families stable on ONE shared Harmonic")
    ax0.legend(fontsize=8)

    n_mt = 1 + p_mt.n_mt * p_mt.beads_per_mt
    mt0 = n_host
    mt1 = n_host + n_mt
    ax1.scatter(pos1[mt0:mt1, 0] * 1e6, pos1[mt0:mt1, 2] * 1e6, s=5,
                color="#1f77b4", label="MT aster")
    ax1.scatter(pos1[mt1:, 0] * 1e6, pos1[mt1:, 2] * 1e6, s=4,
                color="#aa6644", label="IF cage")
    ax1.set_xlabel("x [um]")
    ax1.set_ylabel("z [um]")
    ax1.set_title("Coexisting internal compartments")
    ax1.set_aspect("equal", adjustable="datalim")
    ax1.legend(fontsize=8)

    t = result["topology"]
    fig.suptitle(
        f"MT+IF coexistence smoke — {result['verdict']} "
        f"({t['n_bond_types_on_shared_force']} bond types on "
        f"{t['n_bond_force_computes']} shared Harmonic, SMOKE)", fontsize=10
    )
    sc.save_fig(fig, "smoke_coexistence_mt_if")
    plt.close(fig)


def main() -> int:
    res = run()
    path = sc.save_json("smoke_coexistence_mt_if", res)
    t, o = res["topology"], res["observables"]
    print(f"[coexistence MT+IF smoke] {res['verdict']}  "
          f"{t['n_bond_types_on_shared_force']} bond types on "
          f"{t['n_bond_force_computes']} shared Harmonic; "
          f"mt rms dev={o['mt_backbone_rel_dev_rms']*100:.2f}%, "
          f"if rms dev={o['if_backbone_rel_dev_rms']*100:.2f}%, dt={res['dt_s']:.2e}s")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "COEXIST_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())

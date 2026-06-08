"""CAPSTONE coexistence — confined_migration internal stack (PLUMBING ONLY).

The `confined_migration` recipe declares LINC + intermediate_filaments +
microtubules together (around the nucleus). This harness assembles that internal
stack — a nucleus envelope, a microtubule aster, an intermediate-filament
perinuclear cage, and LINC nesprin bridges coupling the nucleus to perinuclear
acceptors — registering ALL FOUR bond families (linc_nesprin + mt_backbone +
if_backbone + if_crosslink_b*) on ONE shared md.bond.Harmonic + the MT bending
angle on ONE md.angle.Harmonic + ONE BAOAB gamma_map covering every type, and
steps them together. Confirms the multi-compartment internal-stack integration
pattern assembles + steps without crashing.

NOT a physics claim. SMOKE-ONLY constants as in the per-compartment harnesses.

Run:  python ffn_sim/scripts/compartment_smoke/stack_confined_migration.py
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

_HERE = Path(__file__).resolve()
for _p in (_HERE.parents[3], _HERE.parents[0]):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import gsd.hoomd
import hoomd.md as md

import _smoke_common as sc
from ffn_sim.cell.intermediate_filaments import (
    extend_snapshot_with_if_cage,
    register_if_bond_params,
    resolve_intermediate_filaments,
)
from ffn_sim.cell.linc import (
    configure_linc_bond_potential,
    extend_snapshot_with_linc,
    resolve_linc,
)
from ffn_sim.cell.microtubules import (
    extend_snapshot_with_microtubules,
    resolve_microtubules,
)

R_NUC, R_CELL = 3.0e-6, 7.5e-6
N_NUC = 30
R0_LINC = 50.0e-9
R_MTBEAD, R_IFBEAD, R_NUCBEAD = 12.5e-9, 5.0e-9, 150e-9
N_STEPS = 2500


def _host_frame():
    """cortex stub + nucleus envelope (R_nuc) + perinuclear acceptors (R_nuc+r0)."""
    n_stub = 6
    stub = sc.cortex_stub_frame(n_beads=n_stub, r_shell=R_CELL, box_l=6e-5).particles.position
    idx = np.arange(N_NUC) + 0.5
    phi = math.pi * (3 - math.sqrt(5))
    ct = np.clip(1 - 2 * idx / N_NUC, -1, 1)
    st = np.sqrt(np.maximum(0, 1 - ct * ct))
    az = phi * idx
    dirs = np.stack([st * np.cos(az), st * np.sin(az), ct], axis=1)
    nuc = R_NUC * dirs
    acc = (R_NUC + R0_LINC) * dirs
    pos = np.concatenate([np.asarray(stub), nuc, acc], axis=0)
    n = pos.shape[0]
    fr = gsd.hoomd.Frame()
    fr.particles.N = n
    fr.particles.types = ["cortex_actin", "nucleus_bead"]
    fr.particles.typeid = [0] * n_stub + [1] * N_NUC + [0] * N_NUC  # acceptors=cortex_actin
    fr.particles.position = pos.tolist()
    fr.particles.mass = [1.0] * n
    fr.particles.velocity = [[0.0, 0.0, 0.0]] * n
    fr.particles.image = [[0, 0, 0]] * n
    fr.configuration.box = [6e-5, 6e-5, 6e-5, 0, 0, 0]
    fr.bonds.N = 0
    fr.bonds.types = []
    fr.angles.N = 0
    fr.angles.types = []
    return fr, n_stub


def run() -> dict:
    gb_mt, gif, gnuc = sc.stokes_drag(R_MTBEAD), sc.stokes_drag(R_IFBEAD), sc.stokes_drag(R_NUCBEAD)
    p_mt = resolve_microtubules(
        {"microtubules": {"enabled": True, "n_mt": 6, "beads_per_mt": 8,
                          "L_mt": 4.0e-6, "Y_stretch": 2.3e-7}},
        kT=sc.KT, gamma_b=gb_mt,
    )
    p_if = resolve_intermediate_filaments(
        {"intermediate_filaments": {"enabled": True, "n_filaments": 20,
                                    "beads_per_fil": 5}},
        kT=sc.KT, R_cell=R_CELL, R_nuc=R_NUC,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # perinuclear-cap capture warning (expected)
        p_linc = resolve_linc(
            {"linc": {"enabled": True, "k_linc": 1.0e-2, "r0": R0_LINC,
                      "capture_radius": 150e-9}},
            R_cell=R_CELL, R_nuc=R_NUC,
        )

    host, n_stub = _host_frame()
    nuc_tags = np.arange(n_stub, n_stub + N_NUC)
    acc_tags = np.arange(n_stub + N_NUC, n_stub + 2 * N_NUC)

    # LINC first (mutates the gsd frame in place), then MT, then IF (append-only,
    # so the LINC bond tags stay valid).
    host = extend_snapshot_with_linc(host, p_linc, nucleus_tags=nuc_tags,
                                     cytoskeleton_tags=acc_tags)
    n_linc = int(host.bonds.N)
    snap = extend_snapshot_with_microtubules(host, p_mt)
    snap = extend_snapshot_with_if_cage(snap, p_if, centroid=(0, 0, 0), gamma_if=gif, seed=2)

    sim = sc.make_sim(snap, seed=81)
    dt = 0.5 * min(p_mt.dt_cfl, 0.1 * gif / p_if.k_bb)
    ig = sc.set_integrator(sim, dt)

    # ONE shared bond Harmonic covering ALL FOUR families.
    shared_bond = md.bond.Harmonic()
    register_if_bond_params(shared_bond, p_if)
    shared_bond.params["mt_backbone"] = dict(k=p_mt.k_backbone, r0=p_mt.l0)
    configure_linc_bond_potential(shared_bond, p_linc)
    ig.forces.append(shared_bond)
    shared_angle = md.angle.Harmonic()
    shared_angle.params["mt_bending"] = dict(k=p_mt.k_angle, t0=p_mt.angle_t0)
    ig.forces.append(shared_angle)

    gamma_map = sc.auto_gamma_map(
        sim, radius_by_type={"mtoc": R_MTBEAD, "mt_bead": R_MTBEAD,
                             "if_bead": R_IFBEAD, "nucleus_bead": R_NUCBEAD,
                             "cortex_actin": R_NUCBEAD},
    )
    sc.attach_baoab(sim, dt=dt, gamma_map=gamma_map, seed=82)

    pos0 = sc.positions_by_tag(sim)
    sim.run(0)
    sim.run(N_STEPS)
    pos1 = sc.positions_by_tag(sim)
    finite = sc.all_finite(sim)

    final = sim.state.get_snapshot()
    btypes = set(final.bonds.types)
    ptypes = set(final.particles.types)
    n_bond_force = sum(1 for f in ig.forces if isinstance(f, md.bond.Harmonic))
    families = {
        "linc_nesprin": "linc_nesprin" in btypes,
        "mt_backbone": "mt_backbone" in btypes,
        "if_backbone": "if_backbone" in btypes,
        "if_crosslink": any(t.startswith("if_crosslink") for t in btypes),
    }
    ok = (finite and all(families.values()) and n_bond_force == 1
          and n_linc > 0 and {"nucleus_bead", "mt_bead", "if_bead"} <= ptypes)

    result = {
        "compartment": "stack_confined_migration(nucleus+MT+IF+LINC)",
        "verdict": "STACK_OK" if ok else "STACK_FAIL",
        "physics_claim": False,
        "n_steps": N_STEPS, "dt_s": dt,
        "topology": {
            "n_linc_bridges": n_linc,
            "n_bond_types_total": len(btypes),
            "bond_families_present": families,
            "n_bond_force_computes": n_bond_force,
            "particle_types": sorted(ptypes),
        },
        "observables": {
            "all_four_families_on_one_shared_force": bool(all(families.values()) and n_bond_force == 1),
            "max_particle_disp_m": float(np.max(np.linalg.norm(pos1 - pos0, axis=1))),
            "all_finite": finite,
        },
        "note": (
            "PLUMBING capstone: the confined_migration internal stack (nucleus "
            "envelope + MT aster + IF cage + LINC bridges) on ONE shared "
            "md.bond.Harmonic (4 bond families) + ONE angle force + ONE gamma_map, "
            "stepped together. Validates the recipe-declared internal-compartment "
            "integration. NOT the confined-migration physics gate."
        ),
    }

    _figure(pos1, n_stub, p_mt, result)
    return result


def _figure(pos1, n_stub, p_mt, result) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    types = result["topology"]["particle_types"]
    # crude per-type coloring by tag ranges is fiddly post-extend; color by radius.
    r = np.linalg.norm(pos1, axis=1) * 1e6
    sc_plot = ax.scatter(pos1[:, 0] * 1e6, pos1[:, 2] * 1e6, c=r, s=6,
                         cmap="viridis")
    ax.set_xlabel("x [um]")
    ax.set_ylabel("z [um]")
    ax.set_aspect("equal", adjustable="datalim")
    t = result["topology"]
    ax.set_title(
        f"confined_migration internal stack — {result['verdict']}\n"
        f"{t['n_bond_types_total']} bond types ({sum(t['bond_families_present'].values())}/4 "
        f"families) on {t['n_bond_force_computes']} shared Harmonic, "
        f"{t['n_linc_bridges']} LINC bridges", fontsize=9)
    fig.colorbar(sc_plot, ax=ax, label="radius [um]")
    sc.save_fig(fig, "smoke_stack_confined_migration")
    plt.close(fig)


def main() -> int:
    res = run()
    path = sc.save_json("smoke_stack_confined_migration", res)
    t = res["topology"]
    fam = sum(t["bond_families_present"].values())
    print(f"[stack confined_migration] {res['verdict']}  "
          f"{fam}/4 bond families ({t['n_bond_types_total']} types) on "
          f"{t['n_bond_force_computes']} shared Harmonic, "
          f"{t['n_linc_bridges']} LINC bridges, dt={res['dt_s']:.2e}s")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "STACK_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())

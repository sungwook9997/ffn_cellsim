"""ENABLED-PATH smoke harness — intermediate-filament cage (PLUMBING ONLY).

Assembles the default-OFF keratin/vimentin IF perinuclear cage (linear backbone
+ per-r0-bin force-free crosslinks) onto a tiny inert cortex stub, registers the
shared-bond params, steps with the project L-M BAOAB integrator, and confirms the
enabled build+attach+step path runs without crashing, that construction is
force-free (the audit's r0=0 pre-tension fix), and the cage stays bounded.

NOT a physics claim. The nonlinear strain-stiffening law is PI-pending and stays
NotImplementedError — this smoke runs ONLY the ratified linear small-strain path.
E_if/Lp/d_if use the module's cited defaults; n_filaments/beads_per_fil are
SMOKE-ONLY layout knobs. No band pass/fail.

Run:  python ffn_sim/scripts/compartment_smoke/intermediate_filaments.py
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

import _smoke_common as sc
from ffn_sim.archive.hoomd_legacy.cell.intermediate_filaments import (
    GAMMA_DENYLIST_PREFIX,
    attach_if_bonds_to_simulation,
    extend_snapshot_with_if_cage,
    resolve_intermediate_filaments,
)

# --------------------------------------------------------------------------
# SMOKE-ONLY / module-default constants
# --------------------------------------------------------------------------
N_FIL = 60            # IF count — SMOKE-ONLY layout knob (biology count)
BEADS_PER_FIL = 8     # SMOKE-ONLY layout knob (grid-invariant)
R_CELL = 7.5e-6       # m  MCF7 (Wagner 2011)
R_NUC = 3.0e-6        # m  ~0.4 R_cell (cage seeds just outside the nucleus)
R_IFBEAD = 5.0e-9     # m  IF radius (d_if/2 = 5 nm; canonical 10 nm assembled IF)
# E_if=6e6 (ORDER_ESTIMATE), Lp=0.5e-6 keratin, d_if=10e-9, max_stretch_ratio=3.0
# all use the module's cited defaults (in-band).
N_STEPS = 3000


def run() -> dict:
    p = resolve_intermediate_filaments(
        {"intermediate_filaments": {
            "enabled": True, "n_filaments": N_FIL, "beads_per_fil": BEADS_PER_FIL,
        }},
        kT=sc.KT, R_cell=R_CELL, R_nuc=R_NUC,
    )
    gamma_if = sc.stokes_drag(R_IFBEAD)

    host = sc.to_hoomd_snapshot(
        sc.cortex_stub_frame(n_beads=12, r_shell=1.0e-6, box_l=6.0e-5)
    )
    n_host = int(host.particles.N)
    snap = extend_snapshot_with_if_cage(
        host, p, centroid=(0.0, 0.0, 0.0), gamma_if=gamma_if, seed=3
    )

    sim = sc.make_sim(snap, seed=11)
    # IF backbone is soft (k_bb ~ 1e-3 N/m); CFL dt = 0.1*gamma_if/k_bb is huge.
    dt = min(6.98e-7, 0.5 * 0.1 * gamma_if / p.k_bb)
    sc.set_integrator(sim, dt)
    bond = attach_if_bonds_to_simulation(
        sim, p, gamma_if=gamma_if, cfl_strict=True
    )
    gamma_map = sc.auto_gamma_map(sim, radius_by_type={"if_bead": R_IFBEAD})
    sc.attach_baoab(sim, dt=dt, gamma_map=gamma_map, seed=12)

    pos0 = sc.positions_by_tag(sim)
    sim.run(0)
    e0 = float(bond.energy)                    # construction bond energy (force-free?)
    sim.run(N_STEPS)
    pos1 = sc.positions_by_tag(sim)
    finite = sc.all_finite(sim)
    e1 = float(bond.energy)

    # ---- observables (plumbing sanity) ----
    final = sim.state.get_snapshot()
    groups = np.asarray(final.bonds.group, dtype=np.int64).reshape(-1, 2)
    btypeid = np.asarray(final.bonds.typeid, dtype=np.int64).reshape(-1)
    btypes = list(final.bonds.types)
    is_backbone = np.array([btypes[t] == "if_backbone" for t in btypeid], bool)
    is_xl = np.array(
        [btypes[t].startswith("if_crosslink") for t in btypeid], bool
    )
    bb_groups = groups[is_backbone]
    bb_len = np.linalg.norm(pos1[bb_groups[:, 0]] - pos1[bb_groups[:, 1]], axis=1)
    bb_rel_dev = (bb_len - p.l_seg) / p.l_seg
    n_xl = int(np.sum(is_xl))
    n_xl_bin_types = len({btypes[t] for t in btypeid if btypes[t].startswith("if_crosslink")})

    cage_tags = np.arange(n_host, pos1.shape[0])
    r_cage = np.linalg.norm(pos1[cage_tags], axis=1)   # dist from origin centroid
    max_disp = float(np.max(np.linalg.norm(pos1 - pos0, axis=1)))
    e0_kT = e0 / sc.KT

    result = {
        "compartment": "intermediate_filaments",
        "verdict": "PLUMBING_OK" if finite else "PLUMBING_FAIL",
        "physics_claim": False,
        "n_steps": N_STEPS,
        "dt_s": dt,
        "params": {
            "n_filaments": p.n_filaments, "beads_per_fil": p.beads_per_fil,
            "l_seg_m": p.l_seg, "E_if_Pa": p.E_if, "Lp_m": p.Lp,
            "d_if_m": p.d_if, "k_bb_Npm": p.k_bb, "k_xl_Npm": p.k_xl,
            "ratio_xl": p.ratio_xl, "max_stretch_ratio": p.max_stretch_ratio,
            "R_cage_inner_m": p.R_cage_inner, "R_cage_outer_m": p.R_cage_outer,
            "gamma_if_Nspm": gamma_if,
        },
        "topology": {
            "n_host_stub": n_host,
            "n_if_beads": int(p.n_filaments * p.beads_per_fil),
            "n_backbone_bonds": int(bb_groups.shape[0]),
            "n_crosslink_bonds": n_xl,
            "n_crosslink_bin_types": n_xl_bin_types,
        },
        "observables": {
            "backbone_len_mean_m": float(np.mean(bb_len)),
            "backbone_rel_dev_rms": float(np.sqrt(np.mean(bb_rel_dev**2))),
            "construction_bond_energy_kT": e0_kT,
            "construction_force_free": bool(e0_kT < 10.0 * max(1, n_xl)),
            "post_run_bond_energy_kT": e1 / sc.KT,
            "cage_r_min_m": float(np.min(r_cage)),
            "cage_r_max_m": float(np.max(r_cage)),
            "max_particle_disp_m": max_disp,
            "all_finite": finite,
        },
        "constants_smoke_only": [
            "n_filaments=60 (layout)", "beads_per_fil=8 (layout)",
            "E_if=6e6 (ORDER_ESTIMATE module default)",
        ],
        "note": (
            "PLUMBING smoke: linear small-strain IF path only (nonlinear law is "
            "PI-pending / NotImplementedError). Confirms shared-bond registration "
            "of if_backbone + per-r0-bin if_crosslink_b{i}, force-free "
            "construction (r0=0 pre-tension audit fix), and BAOAB stepping. NOT "
            "the delta-gamma activation gate."
        ),
    }

    _figure(p, pos0, pos1, cage_tags, bb_len, result)
    return result


def _figure(p, pos0, pos1, cage_tags, bb_len, result) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.6))

    ax0.hist(bb_len * 1e9, bins=30, color="#aa6644", alpha=0.85)
    ax0.axvline(p.l_seg * 1e9, color="navy", lw=2,
                label=f"rest l_seg = {p.l_seg*1e9:.0f} nm")
    ax0.set_xlabel("if_backbone bond length [nm]")
    ax0.set_ylabel("count")
    ax0.set_title("IF backbone lengths (rod stability)")
    ax0.legend(fontsize=8)

    th = np.linspace(0, 2 * np.pi, 200)
    for R, c, lab in ((R_NUC, "0.5", "R_nuc"),
                      (p.R_cage_inner, "#88bb88", "cage inner"),
                      (p.R_cage_outer, "#226622", "cage outer")):
        ax1.plot(R * 1e6 * np.cos(th), R * 1e6 * np.sin(th), "--",
                 color=c, lw=1, label=lab)
    ax1.scatter(pos0[cage_tags, 0] * 1e6, pos0[cage_tags, 1] * 1e6,
                s=4, color="0.7", label="initial")
    ax1.scatter(pos1[cage_tags, 0] * 1e6, pos1[cage_tags, 1] * 1e6,
                s=4, color="#aa6644", label="final")
    ax1.set_xlabel("x [um]")
    ax1.set_ylabel("y [um]")
    ax1.set_title(f"IF cage (n_fil={p.n_filaments}, {result['n_steps']} steps)")
    ax1.set_aspect("equal", adjustable="datalim")
    ax1.legend(fontsize=7, loc="upper right")

    fig.suptitle(
        f"intermediate_filaments ENABLED-PATH smoke — {result['verdict']} "
        f"(construction {result['observables']['construction_bond_energy_kT']:.0f} kT, "
        f"SMOKE-ONLY)", fontsize=10
    )
    sc.save_fig(fig, "smoke_intermediate_filaments")
    plt.close(fig)


def main() -> int:
    res = run()
    path = sc.save_json("smoke_intermediate_filaments", res)
    obs, top = res["observables"], res["topology"]
    print(f"[intermediate_filaments smoke] {res['verdict']}  "
          f"backbone rms dev={obs['backbone_rel_dev_rms']*100:.2f}%, "
          f"construction={obs['construction_bond_energy_kT']:.1f} kT "
          f"(force-free={obs['construction_force_free']}), "
          f"{top['n_backbone_bonds']} backbone + {top['n_crosslink_bonds']} xl "
          f"({top['n_crosslink_bin_types']} bins), dt={res['dt_s']:.2e}s")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PLUMBING_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())

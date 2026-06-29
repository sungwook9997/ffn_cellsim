"""ENABLED-PATH smoke harness — cadherin cell-cell junction (PLUMBING ONLY).

Builds a minimal TWO-interface patch: cell-A cadherin tips and cell-B cadherin
tips facing each other across a ~10 nm contact zone, registers the cadherin_trans
bond type, attaches the Rakshit catch-bond binder + elastic force, and steps with
the project L-M BAOAB integrator. Confirms the enabled extend+attach+bind path
runs without crashing, the binder forms trans-dimers ONLY between A and B (no
intra-cell), and the dimer tension is finite.

NOT a physics claim. n_cad_per_cell=40 is a SMOKE-ONLY layout count (NOT the
Iturri-2020 N_cad=223 scale bridge); k_trans/r_bind/r0_trans are the DERIVED
defaults from the contact-zone bridge. The per-dimer 1-5 pN Sim-2015 band is the
GATE-J activation gate (settled equilibrated doublet), NOT this plumbing smoke.

Run:  python ffn_sim/scripts/compartment_smoke/cadherin_junction.py
"""

from __future__ import annotations

import math
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
from ffn_sim.archive.hoomd_legacy.junction.cadherin import (
    CADHERIN_TRANS_BOND,
    attach_cadherin_junction,
    extend_snapshot_with_cadherins,
    resolve_cadherin_junction,
)

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------
N_CAD = 24            # cadherins per cell — SMOKE-ONLY layout (NOT Iturri 223)
CONTACT_ZONE = 1.0e-8  # m  engagement length (spheroid scale bridge); 10 nm
GAP = 9.0e-9          # m  A-B vertical separation (< r_bind, ~ r0_trans)
SPACING = 30.0e-9     # m  lateral cadherin spacing (> r_bind → 1:1 pairing)
R_CAD = 1.0e-6        # m  EFFECTIVE membrane-anchored drag radius — SMOKE-ONLY.
#   Cadherins are cortex/membrane-tethered, NOT freely 3D-diffusing; at a bare
#   5 nm tip drag they random-walk ~10 nm over the run and the facing pairs drift
#   out of r_bind before binding. A high effective drag (low mobility) keeps the
#   interface registered so the bind→elastic-force path is exercised. Physical
#   f0/k_trans are unchanged (from the contact-zone bridge, not from R_CAD).
# Binding is a per-batch Bernoulli at rate k_on (~28/s); over a run of length T
# the cumulative bind prob is ~1-exp(-k_on*T). batch_steps is CFL-capped small,
# so we run long enough (T ~ 0.07 s) that most facing pairs engage — a plumbing
# choice to exercise the bind→elastic-force path, NOT a kinetics measurement.
N_STEPS = 100000


def _interface_points():
    """Two facing patches: A at z=-GAP/2, B at z=+GAP/2 on an aligned grid."""
    side = int(math.ceil(math.sqrt(N_CAD)))
    xs = (np.arange(side) - (side - 1) / 2.0) * SPACING
    gx, gy = np.meshgrid(xs, xs)
    grid = np.stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)], axis=1)[:N_CAD]
    pts_a = grid + np.array([0.0, 0.0, -GAP / 2.0])
    pts_b = grid + np.array([0.0, 0.0, +GAP / 2.0])
    return pts_a, pts_b


def run() -> dict:
    p = resolve_cadherin_junction(
        {"junction": {"cadherin": {"enabled": True, "n_cad_per_cell": N_CAD,
                                   "batch_steps": 100}}},
        kT=sc.KT, dt=6.98e-7, contact_zone_width=CONTACT_ZONE,
    )
    pts_a, pts_b = _interface_points()

    host = sc.cortex_stub_frame(
        n_beads=4, r_shell=0.5e-6, box_l=1.0e-5, bead_type="cortex_actin"
    )
    info = extend_snapshot_with_cadherins(
        host, p, cell_a_surface_points=pts_a, cell_b_surface_points=pts_b
    )
    # Register the cadherin_trans bond type (0 bonds) so the binder can add them.
    host.bonds.types = [CADHERIN_TRANS_BOND]
    host.bonds.N = 0

    sim = sc.make_sim(host, seed=41)
    dt = 6.98e-7
    sc.set_integrator(sim, dt)
    gamma_cad = sc.stokes_drag(R_CAD)
    gamma_map = sc.auto_gamma_map(sim, radius_by_type={"cadherin": R_CAD})
    sc.attach_baoab(sim, dt=dt, gamma_map=gamma_map, seed=42)

    cell_a = info["cell_a_tags"]
    cell_b = info["cell_b_tags"]
    harmonic, _ = attach_cadherin_junction(
        sim, p, cell_a_tags=cell_a, cell_b_tags=cell_b, seed=43,
        gamma_cad=gamma_cad, cfl_strict=True,
    )

    sim.run(0)
    sim.run(N_STEPS)
    pos1 = sc.positions_by_tag(sim)
    finite = sc.all_finite(sim)

    # ---- observables ----
    final = sim.state.get_snapshot()
    btypes = list(final.bonds.types)
    n_b = int(final.bonds.N)
    groups = (
        np.asarray(final.bonds.group, dtype=np.int64).reshape(-1, 2)
        if n_b else np.empty((0, 2), np.int64)
    )
    btypeid = (
        np.asarray(final.bonds.typeid, dtype=np.int64).reshape(-1)
        if n_b else np.empty((0,), np.int64)
    )
    trans_idx = (
        CADHERIN_TRANS_BOND in btypes and btypes.index(CADHERIN_TRANS_BOND)
    )
    is_trans = np.array(
        [btypes[t] == CADHERIN_TRANS_BOND for t in btypeid], dtype=bool
    )
    trans_groups = groups[is_trans]
    n_trans = int(trans_groups.shape[0])

    a_set, b_set = set(cell_a.tolist()), set(cell_b.tolist())
    cross_cell = sum(
        1 for i, j in trans_groups
        if ({int(i), int(j)} & a_set) and ({int(i), int(j)} & b_set)
    )
    intra_cell = n_trans - cross_cell

    if n_trans:
        L = np.linalg.norm(pos1[trans_groups[:, 0]] - pos1[trans_groups[:, 1]], axis=1)
        tension_pN = p.k_trans * np.maximum(0.0, L - p.r0_trans) * 1e12
    else:
        L = np.array([]); tension_pN = np.array([])

    ok = finite and n_trans > 0 and intra_cell == 0
    result = {
        "compartment": "cadherin_junction",
        "verdict": "PLUMBING_OK" if ok else "PLUMBING_FAIL",
        "physics_claim": False,
        "n_steps": N_STEPS,
        "dt_s": dt,
        "params": {
            "n_cad_per_cell": p.n_cad_per_cell, "k_trans_Npm": p.k_trans,
            "r0_trans_m": p.r0_trans, "r_bind_m": p.r_bind, "k_on_per_s": p.k_on,
            "batch_steps": p.batch_steps, "contact_zone_m": CONTACT_ZONE,
            "gamma_cad_Nspm": gamma_cad, "f0_pN": p.catch.f0 * 1e12,
        },
        "topology": {
            "n_cadherin_A": int(cell_a.size), "n_cadherin_B": int(cell_b.size),
            "n_trans_bonds_formed": n_trans,
            "cross_cell_bonds": cross_cell, "intra_cell_bonds": intra_cell,
            "engaged_fraction": (n_trans / N_CAD) if N_CAD else 0.0,
        },
        "observables": {
            "trans_len_mean_m": float(np.mean(L)) if n_trans else None,
            "per_dimer_tension_mean_pN": float(np.mean(tension_pN)) if n_trans else None,
            "per_dimer_tension_max_pN": float(np.max(tension_pN)) if n_trans else None,
            "all_A_B_no_intra": bool(intra_cell == 0),
            "all_finite": finite,
        },
        "constants_smoke_only": ["n_cad_per_cell=40 (layout; not Iturri 223)"],
        "note": (
            "PLUMBING smoke: confirms the catch-bond binder forms cadherin_trans "
            "dimers ONLY between A and B (no intra-cell) and the elastic force "
            "evaluates finite. NOT the GATE-J 1-5 pN equilibrated-doublet gate."
        ),
    }

    _figure(pts_a, pts_b, pos1, info, trans_groups, tension_pN, result)
    return result


def _figure(pts_a, pts_b, pos1, info, trans_groups, tension_pN, result) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.6))

    a = pos1[info["cell_a_tags"]]
    b = pos1[info["cell_b_tags"]]
    ax0.scatter(a[:, 0] * 1e9, a[:, 2] * 1e9, s=18, color="#2166ac", label="cell A")
    ax0.scatter(b[:, 0] * 1e9, b[:, 2] * 1e9, s=18, color="#b2182b", label="cell B")
    for i, j in trans_groups:
        ax0.plot([pos1[i, 0] * 1e9, pos1[j, 0] * 1e9],
                 [pos1[i, 2] * 1e9, pos1[j, 2] * 1e9],
                 color="0.4", lw=0.6, zorder=0)
    ax0.set_xlabel("x [nm]")
    ax0.set_ylabel("z [nm]")
    ax0.set_title(f"trans-dimers ({result['topology']['n_trans_bonds_formed']} "
                  f"A↔B, {result['topology']['intra_cell_bonds']} intra)")
    ax0.legend(fontsize=8)

    ax1.axvspan(1.0, 5.0, color="orange", alpha=0.15,
                label="Sim-2015 1-5 pN (GATE-J, NOT this smoke)")
    if tension_pN.size:
        ax1.hist(tension_pN, bins=20, color="#5aae61", alpha=0.85)
    ax1.set_xlabel("per-dimer tension [pN]")
    ax1.set_ylabel("count")
    ax1.set_title("Trans-dimer tension (context only)")
    ax1.legend(fontsize=7)

    fig.suptitle(
        f"cadherin_junction ENABLED-PATH smoke — {result['verdict']} "
        f"(engaged {result['topology']['engaged_fraction']*100:.0f}%, SMOKE-ONLY)",
        fontsize=10,
    )
    sc.save_fig(fig, "smoke_cadherin_junction")
    plt.close(fig)


def main() -> int:
    res = run()
    path = sc.save_json("smoke_cadherin_junction", res)
    top, obs = res["topology"], res["observables"]
    print(f"[cadherin_junction smoke] {res['verdict']}  "
          f"{top['n_trans_bonds_formed']} trans bonds "
          f"({top['cross_cell_bonds']} A↔B, {top['intra_cell_bonds']} intra), "
          f"engaged {top['engaged_fraction']*100:.0f}%, "
          f"tension mean={obs['per_dimer_tension_mean_pN']} pN, dt={res['dt_s']:.2e}s")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PLUMBING_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())

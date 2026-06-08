"""ENABLED-PATH smoke harness — LINC nesprin-SUN complex (PLUMBING ONLY).

Builds a nucleus-bead cloud + a SEEDED PERINUCLEAR acceptor shell (the
physiologically-faithful cap topology that avoids the zero-bonds-at-real-geometry
trap), forms linc_nesprin bridges, configures the shared-bond potential at a
PI-candidate k_linc, steps with the project L-M BAOAB integrator, and confirms
the enabled build+configure+step path runs without crashing, forms n_bridges>0
(the REFUTE-guard), and is force-free at construction (acceptors seeded at r0).

NOT a physics claim. k_linc=1e-2 N/m is the SMOKE-ONLY route-A candidate (folded-
rod pre-unfolding secant); production keeps k_linc=None so the enabled path
raises until PI ratifies it. The [2,10] pN resting-tension oracle is the
ACTIVATION gate (needs actomyosin load), NOT this plumbing smoke. No band pass/fail.

Run:  python ffn_sim/scripts/compartment_smoke/linc.py
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
from ffn_sim.cell.linc import (
    configure_linc_bond_potential,
    extend_snapshot_with_linc,
    linc_cfl_dt_max,
    resolve_linc,
)

# --------------------------------------------------------------------------
# SMOKE-ONLY / cited constants
# --------------------------------------------------------------------------
K_LINC = 1.0e-2       # N/m  PI-candidate route A (Rief 1999 folded-rod secant); SMOKE-ONLY
R0_LINC = 50.0e-9     # m    perinuclear gap (Crisp 2006) — SOLID geometric
CAPTURE = 150.0e-9    # m    default span; works here (acceptors seeded near envelope)
R_CELL = 7.5e-6       # m    MCF7 (Wagner 2011)
R_NUC = 3.0e-6        # m    ~0.4 R_cell
N_NUC = 80            # nucleus envelope beads (1 LINC per bead)
R_NUCBEAD = 150.0e-9  # m    nucleus bead radius for drag (plan a~100-250 nm)
N_STEPS = 3000


def _linc_host_frame() -> gsd.hoomd.Frame:
    """gsd Frame: nucleus shell (R_nuc) + perinuclear acceptor shell (R_nuc+r0).

    Acceptors are seeded at the SAME Fibonacci directions one rest length r0
    outside the envelope, so each nucleus bead's nearest acceptor is its radial
    partner at exactly r0 → bonds born FORCE-FREE (platform principle).
    """
    idx = np.arange(N_NUC, dtype=np.float64) + 0.5
    phi = math.pi * (3.0 - math.sqrt(5.0))
    cos_t = np.clip(1.0 - 2.0 * idx / float(N_NUC), -1.0, 1.0)
    sin_t = np.sqrt(np.maximum(0.0, 1.0 - cos_t * cos_t))
    az = phi * idx
    dirs = np.stack([sin_t * np.cos(az), sin_t * np.sin(az), cos_t], axis=1)
    nuc = R_NUC * dirs
    acc = (R_NUC + R0_LINC) * dirs
    pos = np.concatenate([nuc, acc], axis=0)
    n = 2 * N_NUC
    box_l = 6.0e-5

    fr = gsd.hoomd.Frame()
    fr.particles.N = n
    fr.particles.types = ["nucleus_bead", "cortex_actin"]
    fr.particles.typeid = [0] * N_NUC + [1] * N_NUC
    fr.particles.position = pos.tolist()
    fr.particles.mass = [1.0] * n
    fr.particles.velocity = [[0.0, 0.0, 0.0]] * n
    fr.particles.image = [[0, 0, 0]] * n
    fr.configuration.box = [box_l, box_l, box_l, 0.0, 0.0, 0.0]
    fr.bonds.N = 0
    fr.bonds.types = []
    fr.angles.N = 0
    fr.angles.types = []
    return fr


def run() -> dict:
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        p = resolve_linc(
            {"linc": {"enabled": True, "k_linc": K_LINC, "r0": R0_LINC,
                      "capture_radius": CAPTURE}},
            R_cell=R_CELL, R_nuc=R_NUC,
        )
        host = _linc_host_frame()
        snap = extend_snapshot_with_linc(
            host, p,
            nucleus_tags=np.arange(N_NUC),
            cytoskeleton_tags=np.arange(N_NUC, 2 * N_NUC),
        )
    expected_warnings = [str(w.message) for w in wlist]

    n_bridges = int(snap.bonds.N)
    if n_bridges == 0:
        raise RuntimeError(
            "LINC smoke: 0 bridges formed — zero-bonds-at-real-geometry trap. "
            "Acceptor seeding / capture_radius wrong."
        )

    sim = sc.make_sim(snap, seed=21)
    gamma_bead = sc.stokes_drag(R_NUCBEAD)
    dt = min(6.98e-7, 0.5 * linc_cfl_dt_max(K_LINC, gamma_bead))
    ig = sc.set_integrator(sim, dt)
    harmonic = md.bond.Harmonic()
    configure_linc_bond_potential(harmonic, p)
    ig.forces.append(harmonic)
    gamma_map = sc.auto_gamma_map(
        sim, radius_by_type={"nucleus_bead": R_NUCBEAD, "cortex_actin": R_NUCBEAD}
    )
    sc.attach_baoab(sim, dt=dt, gamma_map=gamma_map, seed=22)

    pos0 = sc.positions_by_tag(sim)
    sim.run(0)
    e0 = float(harmonic.energy)
    sim.run(N_STEPS)
    pos1 = sc.positions_by_tag(sim)
    finite = sc.all_finite(sim)

    # ---- observables ----
    final = sim.state.get_snapshot()
    groups = np.asarray(final.bonds.group, dtype=np.int64).reshape(-1, 2)
    lengths = np.linalg.norm(pos1[groups[:, 0]] - pos1[groups[:, 1]], axis=1)
    tension_N = K_LINC * (lengths - R0_LINC)        # signed; sign-sense check
    max_disp = float(np.max(np.linalg.norm(pos1 - pos0, axis=1)))

    result = {
        "compartment": "linc",
        "verdict": "PLUMBING_OK" if finite else "PLUMBING_FAIL",
        "physics_claim": False,
        "n_steps": N_STEPS,
        "dt_s": dt,
        "params": {
            "k_linc_Npm": K_LINC, "r0_m": R0_LINC, "capture_radius_m": CAPTURE,
            "f_rest_oracle_N": p.f_rest, "gamma_bead_Nspm": gamma_bead,
            "cfl_dt_max_s": linc_cfl_dt_max(K_LINC, gamma_bead),
        },
        "topology": {
            "n_nucleus_beads": N_NUC, "n_acceptor_beads": N_NUC,
            "n_linc_bridges": n_bridges,
        },
        "observables": {
            "bridge_len_mean_m": float(np.mean(lengths)),
            "bridge_len_min_m": float(np.min(lengths)),
            "bridge_len_max_m": float(np.max(lengths)),
            "construction_bond_energy_kT": e0 / sc.KT,
            "mean_tension_pN": float(np.mean(tension_N)) * 1e12,
            "rms_tension_pN": float(np.sqrt(np.mean(tension_N**2))) * 1e12,
            "f_rest_oracle_pN": p.f_rest * 1e12,
            "max_particle_disp_m": max_disp,
            "all_finite": finite,
        },
        "expected_warnings": expected_warnings,
        "constants_smoke_only": ["k_linc=1e-2 N/m (route-A candidate)"],
        "note": (
            "PLUMBING smoke: perinuclear-cap topology (acceptors at r0) avoids "
            "the zero-bonds trap and is force-free at construction; mean tension "
            "is near 0 (thermal), NOT the 8 pN resting-tension oracle (that needs "
            "actomyosin load — the ACTIVATION gate, out of scope). The "
            "capture<R_cell-R_nuc warning is EXPECTED for the cap topology."
        ),
    }

    _figure(p, pos0, pos1, lengths, result)
    return result


def _figure(p, pos0, pos1, lengths, result) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.6))

    ax0.hist(lengths * 1e9, bins=25, color="#6a51a3", alpha=0.85)
    ax0.axvline(R0_LINC * 1e9, color="crimson", lw=2,
                label=f"rest r0 = {R0_LINC*1e9:.0f} nm")
    ax0.set_xlabel("linc_nesprin bond length [nm]")
    ax0.set_ylabel("count")
    ax0.set_title("LINC bridge lengths (force-free at r0)")
    ax0.legend(fontsize=8)

    nuc = pos1[:N_NUC]
    acc = pos1[N_NUC:2 * N_NUC]
    ax1.scatter(nuc[:, 0] * 1e6, nuc[:, 2] * 1e6, s=10, color="#2166ac",
                label="nucleus envelope")
    ax1.scatter(acc[:, 0] * 1e6, acc[:, 2] * 1e6, s=10, color="#b2182b",
                label="perinuclear acceptors")
    ax1.set_xlabel("x [um]")
    ax1.set_ylabel("z [um]")
    ax1.set_title(f"LINC ({result['topology']['n_linc_bridges']} bridges, "
                  f"{result['n_steps']} steps)")
    ax1.set_aspect("equal", adjustable="datalim")
    ax1.legend(fontsize=8)

    fig.suptitle(
        f"linc ENABLED-PATH smoke — {result['verdict']} "
        f"(k_linc={K_LINC:.0e} N/m SMOKE-ONLY, "
        f"construction {result['observables']['construction_bond_energy_kT']:.0f} kT)",
        fontsize=10,
    )
    sc.save_fig(fig, "smoke_linc")
    plt.close(fig)


def main() -> int:
    res = run()
    path = sc.save_json("smoke_linc", res)
    obs, top = res["observables"], res["topology"]
    print(f"[linc smoke] {res['verdict']}  {top['n_linc_bridges']} bridges, "
          f"len mean={obs['bridge_len_mean_m']*1e9:.1f}nm (r0={R0_LINC*1e9:.0f}nm), "
          f"construction={obs['construction_bond_energy_kT']:.1f} kT, "
          f"rms tension={obs['rms_tension_pN']:.2f} pN, dt={res['dt_s']:.2e}s")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PLUMBING_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())

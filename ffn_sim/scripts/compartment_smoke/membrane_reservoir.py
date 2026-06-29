"""ENABLED-PATH smoke harness — membrane reservoir / tether mesh (PLUMBING ONLY).

Builds an OWN radially-offset mem_node layer outside the cortex shell (resolving
BLOCKER-1: with a shared shell the tethers would be degenerate self-pairs), seeds
the breakable mem_tether mesh (each membrane bead → nearest cortex bead within
max_tether_dist), attaches the builtin mem_tether harmonic, and steps with the
project L-M BAOAB integrator. Confirms the STATIC tether-mesh enabled path runs
without crashing, every tether is cross-layer (membrane↔cortex, no self-pair),
and construction is force-free.

NOT a physics claim. The bleb-rupture updater (MembraneTetherUpdater) is
DELIBERATELY NOT constructed: it raises NotImplementedError unconditionally
(σ_crit_bleb=None PI-pending AND the Bell-Evans act() loop is unimplemented), and
f_excess=None blocks the reservoir-release path. This smoke covers ONLY the static
tether mesh (W_MCA/k_tether/max_tether_dist are anchored, not pending). No band
pass/fail.

Run:  python ffn_sim/scripts/compartment_smoke/membrane_reservoir.py
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

import gsd.hoomd
import hoomd.md as md

import _smoke_common as sc
from ffn_sim.archive.hoomd_legacy.cell.membrane_reservoir import (
    GAMMA_DENYLIST_PREFIX,
    attach_membrane_tether_force,
    build_membrane_tethers,
    resolve_membrane_reservoir,
)

# --------------------------------------------------------------------------
# Constants (W_MCA/k_tether/max_tether_dist are anchored module defaults)
# --------------------------------------------------------------------------
N_SHELL = 120         # cortex (and membrane) beads — SMOKE-ONLY layout
R_CELL = 7.5e-6       # m  MCF7 (Wagner 2011)
OFFSET = 100.0e-9     # m  membrane layer radial offset (< max_tether_dist=200nm)
R_MEMBEAD = 50.0e-9   # m  bead radius (drag)
N_STEPS = 3000


def _two_layer_frame():
    """Cortex shell (radius R_cell) + own membrane layer (radius R_cell+OFFSET)
    on the SAME Fibonacci directions → each membrane bead's nearest cortex is its
    radial partner at OFFSET (cross-layer, non-degenerate)."""
    idx = np.arange(N_SHELL, dtype=np.float64) + 0.5
    phi = math.pi * (3.0 - math.sqrt(5.0))
    cos_t = np.clip(1.0 - 2.0 * idx / float(N_SHELL), -1.0, 1.0)
    sin_t = np.sqrt(np.maximum(0.0, 1.0 - cos_t * cos_t))
    az = phi * idx
    dirs = np.stack([sin_t * np.cos(az), sin_t * np.sin(az), cos_t], axis=1)
    cortex = R_CELL * dirs
    membrane = (R_CELL + OFFSET) * dirs
    pos = np.concatenate([cortex, membrane], axis=0)
    n = 2 * N_SHELL
    fr = gsd.hoomd.Frame()
    fr.particles.N = n
    fr.particles.types = ["cortex_actin", "mem_node"]
    fr.particles.typeid = [0] * N_SHELL + [1] * N_SHELL
    fr.particles.position = pos.tolist()
    fr.particles.mass = [1.0] * n
    fr.configuration.box = [6.0e-5, 6.0e-5, 6.0e-5, 0.0, 0.0, 0.0]
    fr.bonds.N = 0
    fr.bonds.types = []
    fr.angles.N = 0
    fr.angles.types = []
    return fr


def run() -> dict:
    p = resolve_membrane_reservoir(
        {"membrane_reservoir": {"enabled": True}}, R_cell=R_CELL
    )  # W_MCA/k_tether/max_tether_dist defaults; sigma_crit/f_excess stay None
    host = _two_layer_frame()
    snap, layout = build_membrane_tethers(
        host, p,
        membrane_tag_range=(N_SHELL, 2 * N_SHELL),
        cortex_tag_range=(0, N_SHELL),
    )
    if layout.n_tether == 0:
        raise RuntimeError(
            "membrane_reservoir smoke: 0 tethers formed — BLOCKER-1 geometry "
            "(degenerate shared shell?) or offset > max_tether_dist."
        )

    sim = sc.make_sim(snap, seed=61)
    gamma_b = sc.stokes_drag(R_MEMBEAD)
    dt = 6.98e-7
    sc.set_integrator(sim, dt)
    bond = attach_membrane_tether_force(
        sim, p, layout, gamma_b=gamma_b, cfl_strict=True
    )
    gamma_map = sc.auto_gamma_map(
        sim, radius_by_type={"cortex_actin": R_MEMBEAD, "mem_node": R_MEMBEAD}
    )
    sc.attach_baoab(sim, dt=dt, gamma_map=gamma_map, seed=62)

    pos0 = sc.positions_by_tag(sim)
    sim.run(0)
    e0 = float(bond.energy)
    sim.run(N_STEPS)
    pos1 = sc.positions_by_tag(sim)
    finite = sc.all_finite(sim)

    # ---- observables ----
    final = sim.state.get_snapshot()
    groups = np.asarray(final.bonds.group, dtype=np.int64).reshape(-1, 2)
    btypeid = np.asarray(final.bonds.typeid, dtype=np.int64).reshape(-1)
    btypes = list(final.bonds.types)
    is_tether = np.array([btypes[t].startswith("mem_") for t in btypeid], bool)
    tg = groups[is_tether]
    tlen = np.linalg.norm(pos1[tg[:, 0]] - pos1[tg[:, 1]], axis=1)
    # Cross-layer check: one endpoint in cortex [0,N), other in membrane [N,2N).
    cortex_set = set(range(N_SHELL))
    mem_set = set(range(N_SHELL, 2 * N_SHELL))
    cross = sum(
        1 for i, j in tg
        if ({int(i), int(j)} & cortex_set) and ({int(i), int(j)} & mem_set)
    )
    self_pairs = int(tg.shape[0]) - cross
    r0_med = float(np.median(layout.tether_r0))

    ok = finite and layout.n_tether > 0 and self_pairs == 0
    result = {
        "compartment": "membrane_reservoir",
        "verdict": "PLUMBING_OK" if ok else "PLUMBING_FAIL",
        "physics_claim": False,
        "n_steps": N_STEPS,
        "dt_s": dt,
        "params": {
            "W_MCA_Jm2": p.W_MCA, "k_tether_Npm": p.k_tether,
            "max_tether_dist_m": p.max_tether_dist, "offset_m": OFFSET,
            "rupture_force_N": layout.rupture_force, "r0_median_m": r0_med,
            "sigma_crit_bleb": p.sigma_crit_bleb, "f_excess": p.f_excess,
            "gamma_b_Nspm": gamma_b,
        },
        "topology": {
            "n_cortex_beads": N_SHELL, "n_mem_nodes": N_SHELL,
            "n_tethers": int(layout.n_tether),
            "cross_layer_tethers": cross, "self_pair_tethers": self_pairs,
        },
        "observables": {
            "tether_len_mean_m": float(np.mean(tlen)),
            "tether_len_vs_r0_rms_dev": float(
                np.sqrt(np.mean(((tlen - r0_med) / r0_med) ** 2))
            ),
            "construction_bond_energy_kT": e0 / sc.KT,
            "max_particle_disp_m": float(np.max(np.linalg.norm(pos1 - pos0, axis=1))),
            "all_finite": finite,
        },
        "blocked_paths": {
            "MembraneTetherUpdater": "NotImplementedError (unconditional: σ_crit_bleb None + act() unimplemented) — NOT constructed",
            "reservoir_release": "f_excess None — blocked",
        },
        "constants_smoke_only": [],
        "note": (
            "PLUMBING smoke: STATIC tether mesh only (W_MCA/k_tether/"
            "max_tether_dist are anchored, not pending). Bleb rupture + reservoir "
            "release are PI-blocked and intentionally NOT exercised. Confirms own "
            "mem_node layer (BLOCKER-1) → cross-layer tethers, force-free build, "
            "BAOAB stepping. NOT the bleb-nucleation / tension-buffering gate."
        ),
    }

    _figure(pos1, tg, tlen, r0_med, result)
    return result


def _figure(pos1, tg, tlen, r0_med, result) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.6))

    ax0.hist(tlen * 1e9, bins=25, color="#3690c0", alpha=0.85)
    ax0.axvline(r0_med * 1e9, color="crimson", lw=2,
                label=f"median r0 = {r0_med*1e9:.0f} nm")
    ax0.set_xlabel("mem_tether bond length [nm]")
    ax0.set_ylabel("count")
    ax0.set_title("Tether lengths (force-free at construction)")
    ax0.legend(fontsize=8)

    cortex = pos1[:result["topology"]["n_cortex_beads"]]
    mem = pos1[result["topology"]["n_cortex_beads"]:]
    ax1.scatter(cortex[:, 0] * 1e6, cortex[:, 2] * 1e6, s=8, color="#1b9e77",
                label="cortex")
    ax1.scatter(mem[:, 0] * 1e6, mem[:, 2] * 1e6, s=8, color="#d95f02",
                label="mem_node")
    ax1.set_xlabel("x [um]")
    ax1.set_ylabel("z [um]")
    ax1.set_title(f"{result['topology']['n_tethers']} cross-layer tethers")
    ax1.set_aspect("equal", adjustable="datalim")
    ax1.legend(fontsize=8)

    fig.suptitle(
        f"membrane_reservoir ENABLED-PATH smoke — {result['verdict']} "
        f"(static mesh; rupture updater PI-blocked, SMOKE)", fontsize=10
    )
    sc.save_fig(fig, "smoke_membrane_reservoir")
    plt.close(fig)


def main() -> int:
    res = run()
    path = sc.save_json("smoke_membrane_reservoir", res)
    top, obs = res["topology"], res["observables"]
    print(f"[membrane_reservoir smoke] {res['verdict']}  "
          f"{top['n_tethers']} tethers ({top['cross_layer_tethers']} cross-layer, "
          f"{top['self_pair_tethers']} self), len mean={obs['tether_len_mean_m']*1e9:.1f}nm, "
          f"construction={obs['construction_bond_energy_kT']:.1f} kT, dt={res['dt_s']:.2e}s "
          f"[rupture updater PI-blocked]")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PLUMBING_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())

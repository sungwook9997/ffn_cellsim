"""ENABLED-PATH smoke harness — ventral stress fibers (PLUMBING ONLY).

Fabricates FA clutch anchor particles, builds explicit FA→FA actomyosin bundles
(sf_actin backbone + α-actinin bands + FA anchors) at a PI-candidate N_filaments
(→ μ_SF), registers the four sf_ bond types on a shared md.bond.Harmonic, steps
with the project L-M BAOAB integrator, and measures the basal sf_actin tension.
Confirms the enabled layout+register+step+measure path runs without crashing and
the PASSIVE backbone is force-free (no NMII).

NOT a physics claim. N_filaments=20 → μ_SF is the SMOKE-ONLY candidate (production
keeps μ_SF=None → enabled build HALTS until PI ratifies). NMII (active contraction)
is the BLOCKER #2 (SF NMII reuses cortex_myosin_* → would contaminate cortical γ;
needs a distinct sf_myosin_* prefix) and is DELIBERATELY NOT attached here, so the
measured tension is ~0 (passive) — NOT compared to the Kumar 10-30 nN band (that
is the activation gate, with NMII on). No band pass/fail.

Run:  python ffn_sim/scripts/compartment_smoke/stress_fibers.py
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

import gsd.hoomd
import hoomd.md as md

import _smoke_common as sc
from ffn_sim.archive.hoomd_legacy.cell.stress_fibers import (
    extend_snapshot_with_stress_fibers,
    measure_sf_tension,
    register_stress_fiber_bond_params,
    resolve_stress_fibers,
    stress_fiber_dt_cfl,
)

# --------------------------------------------------------------------------
# SMOKE-ONLY / cited constants
# --------------------------------------------------------------------------
N_SF = 6              # bundles — SMOKE-ONLY layout
NB = 24              # beads per SF (grid-invariant)
N_FIL = 20           # vSF cross-section count (Cramer 1997 10-30); SMOKE-ONLY → μ_SF
EA_SINGLE = 4.3e-8   # N  F-actin axial rigidity (Gittes 1993) — SOLID, module default
SPAN = 8.0e-6        # m  FA→FA bundle length (ventral SF few-tens µm)
SEED = 47            # = resolver default seed (we replicate its FA pairing)
R_SFBEAD = 50.0e-9   # m  SF actin bead radius (drag)
R_FA = 1.0e-6        # m  FA anchor effective substrate-anchored drag — SMOKE-ONLY
N_STEPS = 3000


def _fa_host_frame():
    """12 FA anchors placed so the resolver's deterministic permutation pairs
    each bundle's two ends exactly SPAN apart → uniform ell0 → force-free."""
    m = 2 * N_SF
    # Replicate generate_stress_fiber_layout's pairing: order=permutation(M),
    # pairs=order[:2n].reshape(n,2). Same SEED → same pairing.
    order = np.random.default_rng(SEED).permutation(m)
    pairs = order[: 2 * N_SF].reshape(N_SF, 2)
    pos = np.zeros((m, 3), dtype=np.float64)
    for b in range(N_SF):
        ia, ib = int(pairs[b, 0]), int(pairs[b, 1])
        pos[ia] = (b * 3.0e-6, 0.0, 0.0)
        pos[ib] = (b * 3.0e-6, SPAN, 0.0)
    fr = gsd.hoomd.Frame()
    fr.particles.N = m
    fr.particles.types = ["fa_actin_clutch"]
    fr.particles.typeid = [0] * m
    fr.particles.position = pos.tolist()
    fr.particles.mass = [1.0] * m
    fr.configuration.box = [6.0e-5, 6.0e-5, 6.0e-5, 0.0, 0.0, 0.0]
    fr.bonds.N = 0
    fr.bonds.types = []
    fr.angles.N = 0
    fr.angles.types = []
    return fr, pos


def run() -> dict:
    ell0 = SPAN / (NB - 1)
    mu_SF = N_FIL * EA_SINGLE
    k_actin = mu_SF / ell0          # DERIVED (= k_bond); set so measure_sf_tension works
    p = resolve_stress_fibers(
        {"stress_fibers": {
            "enabled": True, "n_SF": N_SF, "n_beads_per_SF": NB,
            "N_filaments": N_FIL, "k_actin": k_actin, "seed": SEED,
        }},
    )
    host, fa_pos = _fa_host_frame()
    fa_tags = np.arange(2 * N_SF, dtype=np.int64)
    snap = extend_snapshot_with_stress_fibers(host, p, fa_pos, fa_tags)
    layout = snap.stress_fiber_layout
    ell0_mean = float(np.mean(layout.ell0_actin))

    sim = sc.make_sim(snap, seed=51)
    gamma_sf = sc.stokes_drag(R_SFBEAD)
    dt = 0.5 * stress_fiber_dt_cfl(p, gamma_b=gamma_sf, ell0=ell0_mean)
    ig = sc.set_integrator(sim, dt)
    bond = md.bond.Harmonic()
    register_stress_fiber_bond_params(bond, p, layout)
    ig.forces.append(bond)
    gamma_map = sc.auto_gamma_map(
        sim, radius_by_type={"sf_actin": R_SFBEAD, "sf_xlink_head": R_SFBEAD,
                             "fa_actin_clutch": R_FA},
    )
    sc.attach_baoab(sim, dt=dt, gamma_map=gamma_map, seed=52)

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
    sf_bond_types_present = sorted({btypes[t] for t in btypeid if btypes[t].startswith("sf_")})

    ten = measure_sf_tension(pos1, groups, btypeid, btypes, p, r0=ell0_mean)

    is_chain = np.array([btypes[t] == "sf_actin_bond" for t in btypeid], bool)
    chain_g = groups[is_chain]
    chain_len = np.linalg.norm(pos1[chain_g[:, 0]] - pos1[chain_g[:, 1]], axis=1)
    rel_dev = (chain_len - ell0_mean) / ell0_mean
    n_sf_beads = int(np.sum([1 for t in np.asarray(final.particles.typeid)
                             if final.particles.types[t] == "sf_actin"]))

    ok = finite and len(sf_bond_types_present) == 4 and chain_g.shape[0] > 0
    result = {
        "compartment": "ventral_stress_fibers",
        "verdict": "PLUMBING_OK" if ok else "PLUMBING_FAIL",
        "physics_claim": False,
        "n_steps": N_STEPS,
        "dt_s": dt,
        "params": {
            "n_SF": p.n_SF, "n_beads_per_SF": p.n_beads_per_SF,
            "N_filaments": p.N_filaments, "mu_SF_N": mu_SF, "EA_single_N": EA_SINGLE,
            "ell0_mean_m": ell0_mean, "k_actin_Npm": k_actin,
            "span_m": SPAN, "gamma_sf_Nspm": gamma_sf,
            "sf_tension_band_N": list(p.sf_tension_band_N),
        },
        "topology": {
            "n_fa_anchors": int(2 * N_SF), "n_sf_actin_beads": n_sf_beads,
            "n_backbone_bonds": int(chain_g.shape[0]),
            "sf_bond_types_present": sf_bond_types_present,
        },
        "observables": {
            "backbone_rel_dev_rms": float(np.sqrt(np.mean(rel_dev**2))),
            "construction_bond_energy_kT": e0 / sc.KT,
            "sf_tension_mean_nN": ten["mean_tension_N"] * 1e9,
            "sf_tension_max_nN": ten["max_tension_N"] * 1e9,
            "sf_tension_in_kumar_band": bool(ten["in_band"]),
            "r0_is_live": ten["r0_is_live"],
            "max_particle_disp_m": float(np.max(np.linalg.norm(pos1 - pos0, axis=1))),
            "all_finite": finite,
        },
        "constants_smoke_only": ["N_filaments=20 → μ_SF (ORDER_ESTIMATE)",
                                 "n_SF=6 layout", "k_actin=μ_SF/ell0 DERIVED"],
        "note": (
            "PLUMBING smoke: PASSIVE backbone, NMII NOT attached (BLOCKER #2: SF "
            "NMII reuses cortex_myosin_* → would contaminate cortical γ; needs a "
            "distinct sf_myosin_* prefix). Tension ~0 is EXPECTED (no active "
            "motor); NOT compared to the Kumar 10-30 nN band (the activation gate)."
        ),
    }

    _figure(pos1, layout, chain_len, ell0_mean, result)
    return result


def _figure(pos1, layout, chain_len, ell0_mean, result) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.6))

    ax0.hist(chain_len * 1e9, bins=30, color="#cc6677", alpha=0.85)
    ax0.axvline(ell0_mean * 1e9, color="navy", lw=2,
                label=f"rest ell0 = {ell0_mean*1e9:.0f} nm")
    ax0.set_xlabel("sf_actin_bond length [nm]")
    ax0.set_ylabel("count")
    ax0.set_title("SF backbone lengths (force-free passive)")
    ax0.legend(fontsize=8)

    for b in range(len(layout.actin_tags)):
        tags = layout.actin_tags[b]
        ax1.plot(pos1[tags, 0] * 1e6, pos1[tags, 1] * 1e6, "-", lw=1.2)
        ax1.scatter([pos1[layout.fa_endpoints[b, 0], 0] * 1e6,
                     pos1[layout.fa_endpoints[b, 1], 0] * 1e6],
                    [pos1[layout.fa_endpoints[b, 0], 1] * 1e6,
                     pos1[layout.fa_endpoints[b, 1], 1] * 1e6],
                    s=30, marker="s", color="k", zorder=5)
    ax1.set_xlabel("x [um]")
    ax1.set_ylabel("y [um]")
    ax1.set_title(f"{len(layout.actin_tags)} FA→FA bundles (□ = FA anchors)")
    ax1.set_aspect("equal", adjustable="datalim")

    obs = result["observables"]
    fig.suptitle(
        f"ventral_stress_fibers ENABLED-PATH smoke — {result['verdict']} "
        f"(passive T={obs['sf_tension_mean_nN']:.2e} nN, NMII off, SMOKE-ONLY)",
        fontsize=10,
    )
    sc.save_fig(fig, "smoke_stress_fibers")
    plt.close(fig)


def main() -> int:
    res = run()
    path = sc.save_json("smoke_stress_fibers", res)
    top, obs = res["topology"], res["observables"]
    print(f"[stress_fibers smoke] {res['verdict']}  "
          f"{res['params']['n_SF']} bundles, {top['n_sf_actin_beads']} sf_actin beads, "
          f"{len(top['sf_bond_types_present'])}/4 sf bond types, "
          f"backbone rms dev={obs['backbone_rel_dev_rms']*100:.2f}%, "
          f"passive T={obs['sf_tension_mean_nN']:.2e} nN (NMII off), dt={res['dt_s']:.2e}s")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PLUMBING_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())

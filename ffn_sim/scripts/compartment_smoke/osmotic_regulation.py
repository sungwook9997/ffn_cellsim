"""ENABLED-PATH smoke harness — dynamic osmotic regulation (PLUMBING ONLY).

Builds a cortex shell + the LIVE enclosed-volume turgor force at the
physiological resting turgor (133 Pa), attaches the osmotic-regulation batch
updater referencing the SAME live force object, imposes a hyperosmotic shock,
steps with the project L-M BAOAB integrator, and confirms the enabled
resolve+attach+tick path runs without crashing, the updater reads the live
ev_force.last_pressure and mutates ev_force.p.V0 each batch, and V0 moves in the
RVD sign (hyperosmotic → V0 down).

NOT a physics claim. The implemented setpoint law relaxes on the seconds-to-
minutes timescale (tau_RVD ~ seconds; tau_Kvol ~ minutes), reported ANALYTICALLY
from the resolver; a short HOOMD run (ms) shows only the INITIAL SLOPE / sign,
NOT full relaxation, and is NOT compared to the Hoffmann band. Lp=1e-12 is the
SMOKE-ONLY MCF7/AQP5 candidate (production keeps 1e-13 / PI-pending).

Run:  python ffn_sim/scripts/compartment_smoke/osmotic_regulation.py
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
from ffn_sim.cortex.enclosed_volume import (
    attach_enclosed_volume_to_simulation,
    resolve_enclosed_volume,
)
from ffn_sim.cortex.osmotic_regulation import (
    attach_osmotic_regulation_to_simulation,
    resolve_osmotic_regulation,
)

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------
N_SHELL = 200         # cortex shell beads (volume estimate)
R_CELL = 7.5e-6       # m  MCF7 (Wagner 2011)
TURGOR = 133.0        # Pa physiological resting turgor (registry; baseline ON)
LP = 1.0e-12          # m/(s·Pa) MCF7/AQP5 candidate (Jung 2011); SMOKE-ONLY
BATCH_STEPS = 50
DELTA_C = -100.0      # mol/m³ hyperosmotic step (Jung 2011 100 mM sorbitol) → RVD
N_STEPS = 3000
N_CHUNKS = 30
R_CORTEXBEAD = 1.0e-7


def run() -> dict:
    host = sc.cortex_stub_frame(
        n_beads=N_SHELL, r_shell=R_CELL, box_l=6.0e-5, bead_type="cortex_actin"
    )
    sim = sc.make_sim(host, seed=31)
    dt = 6.98e-7
    sc.set_integrator(sim, dt)
    gamma_b = sc.stokes_drag(R_CORTEXBEAD)

    p_ev = resolve_enclosed_volume(
        {"enclosed_volume": {"turgor_dP0": TURGOR}}, R_cell=R_CELL
    )
    ev = attach_enclosed_volume_to_simulation(
        sim, p_ev, shell_tag_range=(0, N_SHELL), gamma_b=gamma_b, cfl_strict=True
    )
    gamma_map = sc.auto_gamma_map(sim, radius_by_type={"cortex_actin": R_CORTEXBEAD})
    sc.attach_baoab(sim, dt=dt, gamma_map=gamma_map, seed=32)

    p_osmo = resolve_osmotic_regulation(
        {"osmotic_regulation": {
            "enabled": True, "Lp": LP, "batch_steps": BATCH_STEPS,
            "delta_c": DELTA_C,
        }},
        R_cell=R_CELL, dt=dt, p_enclosed_volume=p_ev,
    )
    action, _ = attach_osmotic_regulation_to_simulation(sim, p_osmo, ev)

    V0_ref = float(p_osmo.V0_ref)
    sim.run(0)
    traj_t, traj_V0, traj_P = [0.0], [float(ev.p.V0)], [float(ev.last_pressure)]
    per = N_STEPS // N_CHUNKS
    for c in range(N_CHUNKS):
        sim.run(per)
        traj_t.append((c + 1) * per * dt)
        traj_V0.append(float(ev.p.V0))
        traj_P.append(float(ev.last_pressure))

    V0_initial, V0_final = traj_V0[0], traj_V0[-1]
    finite = sc.all_finite(sim) and np.all(np.isfinite(traj_V0))
    # RVD: hyperosmotic (delta_c<0) → dP_target<turgor → water leaves → V0 down.
    sign_ok = V0_final < V0_initial
    n_ticks = int(action.n_ticks)

    result = {
        "compartment": "osmotic_regulation",
        "verdict": "PLUMBING_OK" if (finite and n_ticks > 0 and sign_ok) else "PLUMBING_FAIL",
        "physics_claim": False,
        "n_steps": N_STEPS,
        "dt_s": dt,
        "params": {
            "Lp_mspa": LP, "A_mem_m2": float(p_osmo.A_mem),
            "turgor_dP0_Pa": TURGOR, "delta_c_molm3": DELTA_C,
            "dP_target_Pa": float(p_osmo.dP_target),
            "K_vol_Pa": float(p_osmo.K_vol), "Pi_osm_Pa": float(p_osmo.Pi_osm),
            "batch_steps": BATCH_STEPS, "batch_dt_s": float(p_osmo.batch_dt),
        },
        "timescales": {
            "tau_RVD_s": float(p_osmo.tau_RVD),
            "tau_Kvol_s": float(p_osmo.tau_Kvol),
            "run_wallclock_s": float(traj_t[-1]),
            "note": "run time << tau_RVD: smoke shows the initial slope/sign, not full relaxation",
        },
        "observables": {
            "n_ticks": n_ticks,
            "V0_ref_m3": V0_ref,
            "V0_initial_m3": V0_initial,
            "V0_final_m3": V0_final,
            "V0_frac_change": (V0_final - V0_initial) / V0_ref,
            "RVD_sign_correct": bool(sign_ok),
            "last_pressure_initial_Pa": traj_P[0],
            "last_pressure_final_Pa": traj_P[-1],
            "V0_above_floor": bool(V0_final > p_osmo.V0_min),
            "all_finite": finite,
        },
        "constants_smoke_only": ["Lp=1e-12 m/(s·Pa) (MCF7/AQP5 candidate)"],
        "note": (
            "PLUMBING smoke: confirms the updater reads the LIVE "
            "ev_force.last_pressure and mutates ev_force.p.V0 each batch via "
            "object identity, V0 moves with the RVD sign, no crash. NOT the "
            "tau_RVD Hoffmann-band activation gate."
        ),
    }

    _figure(np.array(traj_t), np.array(traj_V0) / V0_ref, np.array(traj_P), result)
    return result


def _figure(t, V0_norm, P, result) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.6))

    ax0.plot(t * 1e3, V0_norm, "-o", ms=3, color="#1b7837")
    ax0.set_xlabel("time [ms]")
    ax0.set_ylabel("V0(t) / V0_ref")
    ax0.set_title("Reference-volume setpoint (RVD: should ramp DOWN)")
    ax0.axhline(1.0, color="0.6", ls="--", lw=1, label="V0_ref")
    ax0.legend(fontsize=8)

    ax1.plot(t * 1e3, P, "-o", ms=3, color="#762a83")
    ax1.axhline(result["params"]["dP_target_Pa"], color="crimson", ls="--", lw=1,
                label=f"dP_target = {result['params']['dP_target_Pa']:.0f} Pa")
    ax1.axhline(result["params"]["turgor_dP0_Pa"], color="0.5", ls=":", lw=1,
                label=f"turgor0 = {result['params']['turgor_dP0_Pa']:.0f} Pa")
    ax1.set_xlabel("time [ms]")
    ax1.set_ylabel("ev_force.last_pressure [Pa]")
    ax1.set_title("Live mechanical pressure read by updater")
    ax1.legend(fontsize=8)

    ts = result["timescales"]
    fig.suptitle(
        f"osmotic_regulation ENABLED-PATH smoke — {result['verdict']} "
        f"({result['observables']['n_ticks']} ticks, "
        f"tau_RVD={ts['tau_RVD_s']:.1f}s, Lp=1e-12 SMOKE-ONLY)", fontsize=10
    )
    sc.save_fig(fig, "smoke_osmotic_regulation")
    plt.close(fig)


def main() -> int:
    res = run()
    path = sc.save_json("smoke_osmotic_regulation", res)
    obs, ts = res["observables"], res["timescales"]
    print(f"[osmotic_regulation smoke] {res['verdict']}  {obs['n_ticks']} ticks, "
          f"V0 frac change={obs['V0_frac_change']*100:.2f}% "
          f"(RVD sign={obs['RVD_sign_correct']}), "
          f"tau_RVD={ts['tau_RVD_s']:.1f}s / tau_Kvol={ts['tau_Kvol_s']:.0f}s, "
          f"run={ts['run_wallclock_s']*1e3:.1f}ms")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PLUMBING_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())

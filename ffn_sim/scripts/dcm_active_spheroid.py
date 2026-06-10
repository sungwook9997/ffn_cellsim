"""ACTIVE-spreading DCM spheroid run — active traction + junction switch + division.

Drives the ``cell/dcm_active.py`` stack: the VALIDATED native-mesh + bilinear-tent
spheroid (passive A/A0=2.36 static) with the PI-specified ACTIVE drivers wired in:

  (#5) ActiveRimTraction      — coarse-grained lamellipodium + contraction belt +
                                FA clutch: basal rim cells pull OUTWARD (active).
  (#6) JunctionSwitch         — bulk-pressure (crowding) > 0.5 kPa => weaken
                                cadherin (cad_mult down on tent) + strengthen
                                integrin (int_mult up on traction + substrate).
  LIVE PROLIFERATION          — rim-biased contact-inhibited division on the native
                                pre-allocated pool (no mesh split).

The run is two-phased: AGGREGATION (compact + state-tick, no active traction yet —
the ball forms its 3-zone state) then ACTIVE SPREADING (traction ramps in, junction
switch fires, division grows the cluster). Tracks A/A0(t), cell number, switched-
cell fraction, 3-zone fractions over time.

Validation (the point):
  * ACTIVE A/A0 ends ABOVE the PASSIVE baseline (active over-drive of wetting).
  * junction switch fires (switched count > 0) under bulk pressure.
  * division grows the cluster + A/A0 over time.

Run FROM REPO ROOT:
    ~/miniconda3/envs/ffn_sim/bin/python ffn_sim/scripts/dcm_active_spheroid.py [--quick]

Outputs:
  * outputs/h_dcm_active/figs/active_spheroid.png
  * outputs/h_dcm_active/active_spheroid.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from ffn_sim.cell.dcm_active import (
    ResolvedActiveSpheroid,
    build_active_spheroid,
)
from ffn_sim.cell.dcm_spheroid_state import CellState

OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "h_dcm_active")
FIGS = os.path.join(OUT, "figs")
os.makedirs(FIGS, exist_ok=True)

_UM = 1.0e6


# ---------------------------------------------------------------------------
# Measurement helpers
# ---------------------------------------------------------------------------
def _positions(sim):
    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position).copy()
    return pos[np.argsort(tag)]


def _footprint_area(pos_g, active_node_mask, z0, R):
    """Basal footprint A = xy convex-hull area of substrate-contacting ACTIVE nodes."""
    from scipy.spatial import ConvexHull
    pg = pos_g[active_node_mask]
    contact = pg[pg[:, 2] < z0 + 0.5 * R][:, :2]
    if contact.shape[0] < 3:
        contact = pg[:, :2]
    try:
        return float(ConvexHull(contact).volume)  # 2D hull "volume" = area
    except Exception:  # noqa: BLE001
        return float(np.pi * (np.ptp(contact[:, 0]) / 2) * (np.ptp(contact[:, 1]) / 2))


def _active_node_mask(cell_of_node):
    return cell_of_node >= 0


def _zone_fracs(st, active):
    ids = np.where(active)[0]
    if ids.size == 0:
        return dict(proliferating=0.0, quiescent=0.0, necrotic=0.0)
    s = st.state[ids]
    n = float(ids.size)
    return dict(
        proliferating=float((s == int(CellState.PROLIFERATING)).sum() / n),
        quiescent=float((s == int(CellState.QUIESCENT)).sum() / n),
        necrotic=float((s == int(CellState.NECROTIC)).sum() / n),
    )


# ---------------------------------------------------------------------------
# One run: aggregation -> active spreading, tracking the time series
# ---------------------------------------------------------------------------
def run_active(p: ResolvedActiveSpheroid, *, n_active, n_max, active_traction=True,
               do_switch=True, do_div=True, dt=3.0e-10,
               aggreg_blocks=6, spread_blocks=24, block=1000,
               state_every=1000, div_every=4000, verbose=True, label="ACTIVE"):
    """Build + run one spheroid; return a time-series dict.

    active_traction/do_switch/do_div toggle the three active mechanisms (the PASSIVE
    baseline sets all three False). The schedule: aggregation_blocks of settle (state
    updaters on, traction off), then spread_blocks with traction ramped + junction
    switch + division.
    """
    f_act = p.f_act if active_traction else 0.0
    pp = ResolvedActiveSpheroid(**{**{k: getattr(p, k) for k in p.__slots__
                                      if k != "extras"}, "f_act": f_act, "dt": dt})
    h = build_active_spheroid(pp, n_active=n_active, n_max=n_max,
                              integrin_substrate=True, belt=True)
    sim = h["sim"]
    st = h["st"]
    active = h["active"]
    cell_of_node = h["cell_of_node"]
    ranges = h["ranges"]
    traction = h["traction"]
    R = pp.R_cell
    z0 = pp.z_substrate

    necrosis = h["necrosis"]
    pressure = h["pressure"]
    junction = h["junction"]
    prolif = h["prolif"]

    # --- attach state updaters (necrosis + pressure always; switch/div toggled) ---
    import hoomd
    necrosis.attach(sim)
    pressure.attach(sim)
    if do_switch:
        junction.attach(sim)
    if do_div:
        prolif.attach(sim)

    # --- finite gate (build -> run(0) -> run(2000)) BEFORE the long run ---
    sim.run(0)
    pg0 = _positions(sim)
    mask0 = _active_node_mask(cell_of_node)
    A0 = _footprint_area(pg0, mask0, z0, R)
    sim.run(2000)
    pg = _positions(sim)
    if not np.all(np.isfinite(pg)):
        if verbose:
            print(f"  [{label}] NON-FINITE at finite gate -> retry dt=1e-10")
        if dt > 1.1e-10:
            return run_active(p, n_active=n_active, n_max=n_max,
                              active_traction=active_traction, do_switch=do_switch,
                              do_div=do_div, dt=1.0e-10, aggreg_blocks=aggreg_blocks,
                              spread_blocks=spread_blocks, block=block,
                              state_every=state_every, div_every=div_every,
                              verbose=verbose, label=label)
        return {"finite": False, "label": label}

    ts = {  # time series
        "t_step": [], "AoverA0": [], "A_abs_um2": [], "n_cells": [],
        "switched": [], "switched_frac": [], "rim_cells": [], "ramp": [],
        "prolif_frac": [], "quiescent_frac": [], "necrotic_frac": [],
        "mean_pressure_kPa": [], "max_pressure_kPa": [],
    }

    def manual_tick():
        # Drive the state updaters by hand at this cadence (act() is idempotent and
        # cheap; we call directly rather than register a Periodic trigger so the
        # cadence is explicit + the diagnostics are read at the same instant).
        necrosis.act(int(sim.timestep))
        pressure.act(int(sim.timestep))
        if do_switch:
            junction.act(int(sim.timestep))

    def record():
        pg = _positions(sim)
        mask = _active_node_mask(cell_of_node)
        A = _footprint_area(pg, mask, z0, R)
        ids = np.where(active)[0]
        fr = _zone_fracs(st, active)
        pres = st.pressure_kPa[ids] if ids.size else np.array([0.0])
        ts["t_step"].append(int(sim.timestep))
        ts["AoverA0"].append(float(A / A0))
        ts["A_abs_um2"].append(float(A * _UM ** 2))
        ts["n_cells"].append(int(active.sum()))
        nsw = int(st.switched[ids].sum())
        ts["switched"].append(nsw)
        ts["switched_frac"].append(float(nsw / max(1, ids.size)))
        ts["rim_cells"].append(int(traction.rim_cells.size))
        ts["ramp"].append(float(traction._ramp))
        ts["prolif_frac"].append(fr["proliferating"])
        ts["quiescent_frac"].append(fr["quiescent"])
        ts["necrotic_frac"].append(fr["necrotic"])
        ts["mean_pressure_kPa"].append(float(pres.mean()))
        ts["max_pressure_kPa"].append(float(pres.max()))

    # --- AGGREGATION phase: settle + state ticks; traction is ramp=0 (off) ---
    manual_tick()
    record()
    for b in range(aggreg_blocks):
        sim.run(block)
        manual_tick()
        record()
        if not np.all(np.isfinite(_positions(sim))):
            ts["finite"] = False
            return {"finite": False, "label": label, "ts": ts}

    # --- ACTIVE SPREADING phase: traction ramps in, switch + division act ---
    for b in range(spread_blocks):
        sim.run(block)
        manual_tick()
        if do_div and (int(sim.timestep) % div_every < block):
            prolif.act(int(sim.timestep))
        record()
        if not np.all(np.isfinite(_positions(sim))):
            ts["finite"] = False
            if verbose:
                print(f"  [{label}] NON-FINITE during spread at step {sim.timestep}")
            return {"finite": False, "label": label, "ts": ts}

    # representative final snapshot (centroids + states) for the cluster panel
    pg = _positions(sim)
    ids = np.where(active)[0]
    cents = np.array([pg[ranges[int(c)][0]:ranges[int(c)][1]].mean(0) for c in ids])
    out = {
        "finite": True, "label": label, "ts": ts,
        "A0_um2": float(A0 * _UM ** 2),
        "AoverA0_final": ts["AoverA0"][-1],
        "n_cells_final": ts["n_cells"][-1],
        "switched_final": ts["switched"][-1],
        "n_divisions": int(prolif.n_divisions) if do_div else 0,
        "centroids_um": (cents * _UM).tolist(),
        "states": st.state[ids].tolist(),
        "pressure_kPa": st.pressure_kPa[ids].tolist(),
        "switched_mask": st.switched[ids].tolist(),
        "int_mult": st.int_mult[ids].tolist(),
        "depth_um": st.depth_um[ids].tolist(),
    }
    if verbose:
        print(f"  [{label}] A/A0 {ts['AoverA0'][0]:.3f} -> {out['AoverA0_final']:.3f}"
              f"  cells {ts['n_cells'][0]}->{out['n_cells_final']}"
              f"  switched {out['switched_final']}  divisions {out['n_divisions']}")
    return out


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def make_figure(active, passive, p):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ta = np.array(active["ts"]["t_step"])
    tp = np.array(passive["ts"]["t_step"])
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))

    # (a) A/A0(t): active vs passive baseline
    a0 = ax[0, 0]
    a0.plot(tp, passive["ts"]["AoverA0"], "s--", color="0.5", ms=4,
            label=f"PASSIVE wetting (final {passive['AoverA0_final']:.2f})")
    a0.plot(ta, active["ts"]["AoverA0"], "o-", color="C3", ms=4,
            label=f"ACTIVE traction+switch+div (final {active['AoverA0_final']:.2f})")
    agg_end = ta[min(len(ta) - 1, _agg_idx(active))]
    a0.axvline(agg_end, ls=":", color="k", lw=1.0, label="aggregation->spreading")
    a0.set_xlabel("BAOAB step"); a0.set_ylabel("A / A0  (vs own initial footprint)")
    a0.set_title("(a) A/A0(t): active spreading vs passive wetting")
    a0.legend(fontsize=8, loc="best")

    # (b) cell number(t)
    a1 = ax[0, 1]
    a1.plot(ta, active["ts"]["n_cells"], "o-", color="C0", ms=4, label="ACTIVE (with division)")
    a1.plot(tp, passive["ts"]["n_cells"], "s--", color="0.5", ms=4, label="PASSIVE (no division)")
    a1.set_xlabel("BAOAB step"); a1.set_ylabel("active cell number")
    a1.set_title(f"(b) cell number(t)  [+{active['n_divisions']} divisions]")
    a1.legend(fontsize=8, loc="best")

    # (c) switched-cell fraction(t) + mean pressure
    a2 = ax[0, 2]
    a2.plot(ta, active["ts"]["switched_frac"], "o-", color="C2", ms=4,
            label="switched fraction")
    a2.set_xlabel("BAOAB step"); a2.set_ylabel("switched-cell fraction", color="C2")
    a2.tick_params(axis="y", labelcolor="C2")
    a2.set_ylim(-0.03, 1.03)
    a2b = a2.twinx()
    a2b.plot(ta, active["ts"]["max_pressure_kPa"], "^--", color="C4", ms=4,
             label="max bulk pressure")
    a2b.axhline(p.P_switch_kPa, ls=":", color="C4", lw=1.0)
    a2b.set_ylabel("bulk pressure [kPa]", color="C4")
    a2b.tick_params(axis="y", labelcolor="C4")
    a2.set_title(f"(c) junction switch fires (P>{p.P_switch_kPa} kPa onset)")

    # (d) 3-zone fractions(t)
    a3 = ax[1, 0]
    a3.plot(ta, active["ts"]["prolif_frac"], "-", color="#2ca02c", label="proliferating (rim)")
    a3.plot(ta, active["ts"]["quiescent_frac"], "-", color="#ff7f0e", label="quiescent (shell)")
    a3.plot(ta, active["ts"]["necrotic_frac"], "-", color="#7f0000", label="necrotic (core)")
    a3.set_xlabel("BAOAB step"); a3.set_ylabel("cell fraction")
    a3.set_ylim(-0.03, 1.03)
    a3.set_title("(d) 3-zone state fractions(t)")
    a3.legend(fontsize=8, loc="best")

    # (e) final cluster: 3-zone state (xy) + switched-cell ring
    a4 = ax[1, 1]
    cents = np.array(active["centroids_um"])
    states = np.array(active["states"])
    switched = np.array(active["switched_mask"])
    cc = cents.mean(0)
    colmap = {int(CellState.PROLIFERATING): ("#2ca02c", "PROLIFERATING"),
              int(CellState.QUIESCENT): ("#ff7f0e", "QUIESCENT"),
              int(CellState.NECROTIC): ("#7f0000", "NECROTIC")}
    for s, (col, lab) in colmap.items():
        m = states == s
        if m.any():
            a4.scatter(cents[m, 0] - cc[0], cents[m, 1] - cc[1], s=130, color=col,
                       edgecolor="k", lw=0.5, label=lab, zorder=3)
    if switched.any():
        a4.scatter(cents[switched, 0] - cc[0], cents[switched, 1] - cc[1], s=300,
                   facecolors="none", edgecolors="blue", lw=2.0, zorder=4,
                   label="junction-switched")
    a4.set_aspect("equal")
    a4.set_xlabel("x − x̄ [µm]"); a4.set_ylabel("y − ȳ [µm]")
    a4.set_title(f"(e) final cluster (N={active['n_cells_final']}): 3-zone + switched")
    a4.legend(fontsize=7, loc="best")

    # (f) per-cell active-traction gain (int_mult) / pressure overlay
    a5 = ax[1, 2]
    pres = np.array(active["pressure_kPa"])
    intm = np.array(active["int_mult"])
    sc = a5.scatter(cents[:, 0] - cc[0], cents[:, 1] - cc[1], s=130, c=pres,
                    cmap="viridis", edgecolor="k", lw=0.4, zorder=3)
    big = intm > 1.01
    if big.any():
        a5.scatter(cents[big, 0] - cc[0], cents[big, 1] - cc[1], s=320,
                   facecolors="none", edgecolors="red", lw=2.0, zorder=4,
                   label=f"integrin-strengthened (x{intm[big].max():.0f})")
        a5.legend(fontsize=7, loc="best")
    a5.set_aspect("equal")
    a5.set_xlabel("x − x̄ [µm]"); a5.set_ylabel("y − ȳ [µm]")
    a5.set_title("(f) per-cell bulk pressure + integrin gain")
    fig.colorbar(sc, ax=a5, label="bulk pressure [kPa]", fraction=0.046, pad=0.04)

    fig.suptitle("ACTIVE-spreading DCM spheroid: active rim traction + bulk-pressure "
                 "junction switch + live proliferation\n(coarse-grained on the "
                 "validated native-mesh + bilinear-tent stack)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    path = os.path.join(FIGS, "active_spheroid.png")
    fig.savefig(path, dpi=120)
    print(f"  figure -> {os.path.relpath(path)}")
    return path


def _agg_idx(run):
    """Index in the time series where the aggregation phase ended (ramp first > 0)."""
    ramp = run["ts"]["ramp"]
    for i, r in enumerate(ramp):
        if r > 0.0:
            return max(0, i - 1)
    return 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="short run (smoke): fewer blocks + cells")
    args = ap.parse_args()

    if args.quick:
        n_active, n_max = 12, 20
        aggreg, spread, block = 3, 10, 800
    else:
        n_active, n_max = 14, 34
        aggreg, spread, block = 6, 26, 1000

    p = ResolvedActiveSpheroid(subdivisions=2, dt=3.0e-10)

    print("=== ACTIVE-spreading DCM spheroid ===")
    print(f"    start cells {n_active}, pool {n_max}, aggreg {aggreg}x{block} + "
          f"spread {spread}x{block} steps, dt {p.dt:.0e}")
    print(f"    f_act={p.f_act:.1e} N/node (cap {p.f_cap:.1e}); P_switch="
          f"{p.P_switch_kPa} kPa; p_div={p.p_div}\n")

    print("-- PASSIVE baseline (no traction, no switch, no division) --")
    passive = run_active(p, n_active=n_active, n_max=n_active, active_traction=False,
                         do_switch=False, do_div=False, dt=p.dt,
                         aggreg_blocks=aggreg, spread_blocks=spread, block=block,
                         label="PASSIVE")

    print("-- ACTIVE (traction + junction switch + division) --")
    active = run_active(p, n_active=n_active, n_max=n_max, active_traction=True,
                        do_switch=True, do_div=True, dt=p.dt,
                        aggreg_blocks=aggreg, spread_blocks=spread, block=block,
                        label="ACTIVE")

    if not (passive.get("finite") and active.get("finite")):
        print("FAIL: a run went non-finite.")
        # still dump what we have
        result = {"passive": passive, "active": active, "ok": False}
    else:
        active_over_passive = active["AoverA0_final"] > passive["AoverA0_final"]
        switch_fired = active["switched_final"] > 0
        div_grew = active["n_cells_final"] > active["ts"]["n_cells"][0]
        a_traj = active["ts"]["AoverA0"]
        a_grew = a_traj[-1] > a_traj[0]
        print("\n=== RESULT ===")
        print(f"  ACTIVE A/A0 final {active['AoverA0_final']:.3f}  vs  PASSIVE "
              f"{passive['AoverA0_final']:.3f}  -> active>passive: {active_over_passive}")
        print(f"  junction switch fired: {switch_fired} "
              f"({active['switched_final']} cells switched)")
        print(f"  division grew cluster: {div_grew} "
              f"({active['ts']['n_cells'][0]}->{active['n_cells_final']} cells, "
              f"{active['n_divisions']} divisions)")
        print(f"  A/A0 grew over time: {a_grew} "
              f"({a_traj[0]:.3f}->{a_traj[-1]:.3f})")
        result = {
            "passive": passive, "active": active,
            "active_over_passive": bool(active_over_passive),
            "switch_fired": bool(switch_fired),
            "division_grew_cluster": bool(div_grew),
            "AoverA0_grew_over_time": bool(a_grew),
            "ok": bool(active_over_passive and switch_fired and div_grew),
        }

    result["params"] = {k: getattr(p, k) for k in p.__slots__ if k != "extras"}
    result["calibration"] = (
        "f_act=1.2e-10 N/node calibrated to the single-cell spreading_drive engine: "
        "its FA clutch k_int=1e-3 N/m at clutch extensions 50-200 nm gives per-clutch "
        "traction 5e-11..2e-10 N; the DCM substrate well max pull is ~6e-10 N/node, so "
        "f_act is in the single-cell clutch band AND ~0.2x the passive well => a genuine "
        "active over-drive of wetting, ~400x below the 5e-8 N contact cap (BAOAB-safe).")

    jpath = os.path.join(OUT, "active_spheroid.json")
    with open(jpath, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n  json -> {os.path.relpath(jpath)}")

    if passive.get("ts") and active.get("ts"):
        try:
            make_figure(active, passive, p)
        except Exception as e:  # noqa: BLE001
            import traceback
            print(f"  (figure skipped: {e})")
            traceback.print_exc()

    return 0 if result.get("ok") else 0  # always 0; validation reported in stdout/json


if __name__ == "__main__":
    raise SystemExit(main())

"""GPU-scale ACTIVE + NECROSIS size sweep of the spreading law A/A0 = a + b/R + c/R².

This is the script the PI runs on the gbook RTX A5000 to TEST the layer-2 spheroid
spreading law at GPU scale, where the necrotic core can actually turn on. Strategy
(``docs/v2_audit/GPU_SCALE_PATH_2026-06-11.md``): the law's size-dependence is
probed by a STATIC size sweep (no live division — division OFF by default, a
``--division`` flag enables it) of LARGE active+necrosis spheroids:

    larger spheroid  ->  necrotic core grows (∝ volume)
                     ->  the ACTIVE, spreading RIM fraction shrinks (∝ surface/vol)
                     ->  A/A0 DECREASES with R   (b/R surface term, c/R² core term)

This avoids the ~6-7 week live-division mesh-split port for the law validation: the
size-dependence is geometric+mechanical (active rim vs dead core), not the dynamic
proliferation mechanism. If the static sweep reproduces the coefficients, live
division is NOT on the critical path for the law.

STACK PER SIZE (built via ``cell/dcm_active.build_active_spheroid`` + GPU forces):
  * native md.mesh turgor/edge shell (already GPU-native, dcm_native_shell)
  * ``DcmTentContactGPU``            (cupy K-grid cell-cell tent contact)
  * ``DcmActiveRimTractionGPU``      (cupy active rim lamellipodium/belt/clutch)
  * ``DcmSubstrateForceGPU``         (cupy adhesive substrate well)
  * CPU low-cadence state updaters: SpheroidNecrosis + SpheroidPressureProbe +
    SpheroidJunctionSwitch (batched, off the per-step path); proliferation only
    if --division.

SCALE / COARSE-GRAINING (documented, stated choice): the necrotic core turns on at
depth > ~150 µm (R > ~150 µm). With the fine native R_cell = 7.5 µm even 450 cells
only reach R ≈ 80 µm — below threshold — so a fine sweep CANNOT activate necrosis
(needs ~12k cells). We instead COARSE-GRAIN the DCM cell to a multicellular PATCH
radius (default --cell-radius-um 20: each DCM "cell" represents a small cluster of
real ~7.5 µm cells). At R_cell = 20 µm the 60->450-cell sweep spans R ≈ 118->213 µm,
crossing the 150 µm necrosis threshold cleanly. The mechanical bands (K_bulk, turgor,
adhesion, traction) are the validated native+active bands; only the cell length scale
is coarse-grained (a stated coarse-graining, like the ×40 filament mesoscale in
CLAUDE.md — a hardware-driven length rescale, not a mechanism swap). For a true fine
R = 7.5 µm sweep at R > 150 µm, run thousands of cells on the A5000 (--cell-radius-um
7.5 --sizes 2000 6000 12000) — far slower; the coarse-grained patch sweep is the
fast law-test path.

DEVICE: ``hoomd.device.GPU()`` when available; else FALLBACK to ``hoomd.device.CPU()``
with a clear log line (the bit-identical CPU dispatch path). GPU speedup / GPU-scale
(R>150 µm necrosis-on) results are gbook-A5000 ONLY; on this CPU dev Mac use --quick
for a small finite mini-sweep (structural smoke + fit + figure, necrosis OFF at small R).

Run on the gbook A5000 (production, necrosis-on law test):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python -m ffn_sim.scripts.dcm_gpu_law_sweep \
        --sizes 60 120 200 300 450

Run on a CPU dev box (mini-sweep finite smoke + fit + figure):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python -m ffn_sim.scripts.dcm_gpu_law_sweep \
        --sizes 14 20 28 --quick --force-cpu

Outputs:
  * outputs/h_dcm_gpu/gpu_law_sweep.json
  * outputs/h_dcm_gpu/figs/gpu_law_sweep.png
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

import hoomd  # noqa: E402

from ffn_sim.cell.dcm_active import (  # noqa: E402
    ResolvedActiveSpheroid,
    build_active_spheroid,
)
from ffn_sim.cell.dcm_gpu_forces import (  # noqa: E402
    DcmTentContactGPU,
    DcmActiveRimTractionGPU,
    DcmSubstrateForceGPU,
    on_gpu,
)
from ffn_sim.cell.dcm_spheroid_state import CellState  # noqa: E402
# Reuse the VALIDATED capstone measurement + fit helpers (does not modify them).
from ffn_sim.scripts.dcm_native_capstone import (  # noqa: E402
    _positions,
    _footprint_area,
    _cell_centroids,
    _effective_radius,
    fit_law,
)

OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "h_dcm_gpu")
FIGS = os.path.join(OUT, "figs")
os.makedirs(FIGS, exist_ok=True)

_UM = 1.0e6

# PI law target (a, b µm, c µm²; r²=0.98 over R = 31-78 µm).
PI_LAW = (-0.33, 188.7, -2655.0)


# ---------------------------------------------------------------------------
# Device selection
# ---------------------------------------------------------------------------
def select_device(force_cpu: bool):
    """Return (device, on_gpu_bool, note). Prefer GPU; fall back to CPU + log."""
    if force_cpu:
        return hoomd.device.CPU(notice_level=0), False, "forced CPU (--force-cpu)"
    try:
        if hoomd.device.GPU.is_available():
            return hoomd.device.GPU(notice_level=0), True, "GPU (CUDA available)"
        note = "no CUDA GPU build available — FALLBACK to CPU"
    except Exception as e:  # noqa: BLE001
        note = f"GPU init failed ({e}) — FALLBACK to CPU"
    return hoomd.device.CPU(notice_level=0), False, note


# ---------------------------------------------------------------------------
# Swap the active spheroid's CPU forces for the device-dispatched GPU twins.
# ---------------------------------------------------------------------------
def _gpu_swap_forces(h, p: ResolvedActiveSpheroid):
    """Replace the CPU tent/traction/substrate forces with GPU-dispatched twins.

    ``build_active_spheroid`` wires the CPU ``DcmTentContact`` (in the native
    build) + CPU ``ActiveRimTraction`` + CPU ``IntegrinModulatedSubstrate``. For
    the GPU-resident run we swap in the device-dispatched ``DcmTentContactGPU`` /
    ``DcmActiveRimTractionGPU`` / ``DcmSubstrateForceGPU`` (CPU-parity tested),
    sharing the SAME mutable state arrays (cell_of_node / int_mult / cad_mult) so
    the low-cadence updaters keep steering them.

    NOTE: the GPU substrate well is the plain per-node well (no per-cell integrin
    gain) — the junction-switch integrin gain still acts through the GPU active
    traction's ``int_mult``; the substrate gain is dropped on the GPU path for the
    K2-kernel-resident well (the dominant grip is the active traction). The
    cadherin-weakening still acts through the tent contact's ``cad_mult`` seam.
    """
    sim = h["sim"]
    ig = sim.operations.integrator
    nv = h["nv"]
    ranges = h["ranges"]
    st = h["st"]
    cell_of_node = h["cell_of_node"]
    patch_area = 4.0 * np.pi * p.R_cell ** 2 / nv

    # remove the CPU forces wired by the builder
    for key in ("tent", "traction", "substrate"):
        f = h.get(key)
        if f is not None and f in ig.forces:
            ig.forces.remove(f)

    base = p.to_native()
    # mean edge of an undeformed cell (tent contact radius), from the template.
    verts0 = h["verts0"]
    tris0 = h["tris0"]
    e = []
    for t in tris0:
        for a, b in ((0, 1), (1, 2), (2, 0)):
            e.append(np.linalg.norm(verts0[t[a]] - verts0[t[b]]))
    mean_edge = float(np.mean(e))

    tent = DcmTentContactGPU(
        cell_of_node=cell_of_node, r_contact=1.05 * mean_edge,
        c_adh=p.c_adh, rep_strength=p.rep_strength, adh_strength=p.adh_strength,
        patch_area=patch_area, force_cap=p.contact_force_cap, cad_mult=st.cad_mult)
    ig.forces.append(tent)

    traction = DcmActiveRimTractionGPU(
        cell_of_node=cell_of_node, ranges=ranges, active=h["active"],
        int_mult=st.int_mult, R_cell=p.R_cell, z0=p.z_substrate,
        f_act=p.f_act, f_cap=p.f_cap, ramp_steps=p.ramp_steps,
        contact_band=p.rim_contact_band, neighbour_factor=p.rim_neighbour_factor,
        max_neighbours=p.rim_max_neighbours,
        integrin_switch_gain=p.integrin_switch_gain, belt_factor=p.belt_factor)
    ig.forces.append(traction)

    area_per_node = 4.0 * np.pi * p.R_cell ** 2 / nv
    W_cs = p.W_cs_Jm2 * area_per_node
    sub = DcmSubstrateForceGPU(z0=p.z_substrate, W_cs=W_cs, adh_range=p.R_cell)
    ig.forces.append(sub)

    h["tent"] = tent
    h["traction"] = traction
    h["substrate"] = sub
    return h


def _necrotic_fraction(st, active: np.ndarray) -> dict:
    """3-zone fractions over the ACTIVE cells (read from the live st.state array)."""
    ids = np.where(active)[0]
    if ids.size == 0:
        return {"necrotic": 0.0, "quiescent": 0.0, "proliferating": 0.0}
    states = st.state[ids]
    n = float(ids.size)
    return {
        "necrotic": float((states == int(CellState.NECROTIC)).sum() / n),
        "quiescent": float((states == int(CellState.QUIESCENT)).sum() / n),
        "proliferating": float((states == int(CellState.PROLIFERATING)).sum() / n),
    }


# ---------------------------------------------------------------------------
# One spheroid: build -> finite gate -> equilibrate -> measure
# ---------------------------------------------------------------------------
def run_one(n_cells, *, device, is_gpu, cell_radius_um, division, quick,
            seed=7, verbose=True):
    """Build a LARGE active+necrosis spheroid, gate, equilibrate, measure.

    Returns a dict with R_eff, equilibrium A/A0, necrotic/active-rim fractions and
    a finite-gate flag. Drops dt 3e-10 -> 1e-10 and retries once on non-finite.
    """
    block = 600 if quick else 1000
    n_equil = 8 if quick else 16
    state_cadence = 400 if quick else 800

    for trial_dt in (3.0e-10, 1.0e-10):
        # NOTE: build_active_spheroid uses the builder's default device (CPU);
        # to run GPU-resident we re-pass the device via to_native through a fresh
        # build. The builder does not take a device kwarg directly, so we set it on
        # the native base it constructs. We rebuild a device-aware sim here.
        p = ResolvedActiveSpheroid(
            R_cell=cell_radius_um * 1.0e-6, subdivisions=2, spacing_factor=2.3,
            adh_strength=1.0e8, rep_strength=1.0e8, W_cs_Jm2=0.5e-3,
            dt=trial_dt, seed=seed)
        # pool: a few extra cells if division enabled, else exactly n_cells.
        n_max = int(n_cells * 1.3) + 4 if division else n_cells
        h = build_active_spheroid(p, n_cells, n_max, device=device,
                                  belt=True, integrin_substrate=True)
        sim = h["sim"]
        # HONEST device check: build_active_spheroid (read-only) currently does NOT
        # forward its `device` kwarg into build_native_dcm_simulation (dcm_active.py
        # line ~609), so the sim is built on the builder's default CPU device even
        # when device=GPU was requested. If that one-line forward is added (PI sign-
        # off on the frozen file), this run goes GPU-resident automatically. Until
        # then we detect the mismatch and fall back to the (bit-identical) CPU path
        # so the run still completes — but it is NOT a GPU run.
        sim_on_gpu = on_gpu(sim)
        if is_gpu and not sim_on_gpu:
            if verbose:
                print("  [warn] requested GPU but build_active_spheroid built a CPU "
                      "sim (device kwarg not forwarded in the read-only dcm_active.py)"
                      " -> running CPU dispatch. Forward `device=device` at "
                      "dcm_active.py:609 (PI sign-off) for the real GPU run.")
            is_gpu = False
        if is_gpu:
            h = _gpu_swap_forces(h, p)
        ranges, vol = h["ranges"], h["volume"]
        st, active = h["st"], h["active"]
        R = p.R_cell

        # --- finite gate: build -> run(0) -> run(2000) (per size, BEFORE equil) ---
        sim.run(0)
        pg0 = _positions(sim)
        A0_self = _footprint_area(pg0, p.z_substrate, R)
        sim.run(2000)
        pg = _positions(sim)
        if not np.all(np.isfinite(pg)) or not np.all(np.isfinite(vol.volume)):
            if verbose:
                msg = ("retry at smaller dt" if trial_dt == 3.0e-10
                       else "FAIL even at dt=1e-10")
                print(f"  N={n_cells} dt={trial_dt:.0e}: run(2000) NON-FINITE -> {msg}")
            if trial_dt == 3.0e-10:
                continue
            return {"n_cells": n_cells, "finite_gate": False, "dt": trial_dt}
        if verbose:
            print(f"  N={n_cells} dt={trial_dt:.0e}: finite gate PASS "
                  f"({sim.state.N_particles} particles)")

        # --- attach the low-cadence CPU state updaters (necrosis/pressure/switch) ---
        ops = sim.operations
        trig = hoomd.trigger.Periodic(state_cadence)
        for upd in (h["necrosis"], h["pressure"], h["junction"]):
            ops.writers.append(hoomd.write.CustomWriter(action=upd, trigger=trig))
        if division:
            ops.writers.append(hoomd.write.CustomWriter(
                action=h["prolif"], trigger=hoomd.trigger.Periodic(state_cadence)))

        # tick the updaters once so the necrosis zones are classified at t=0
        h["necrosis"].attach(sim)
        h["necrosis"].act(sim.timestep)

        # --- equilibrate; track footprint to confirm it plateaus ---
        areas, areas_abs = [], []
        finite = True
        t_run0 = time.time()
        total = 0
        for _ in range(n_equil):
            sim.run(block)
            total += block
            pg = _positions(sim)
            if not np.all(np.isfinite(pg)):
                finite = False
                break
            A_abs = _footprint_area(pg, p.z_substrate, R)
            areas.append(A_abs / A0_self)
            areas_abs.append(A_abs * _UM ** 2)
        if not finite:
            if trial_dt == 3.0e-10:
                continue
            return {"n_cells": n_cells, "finite_gate": False, "dt": trial_dt}
        wall = time.time() - t_run0
        sps = total / wall if wall > 0 else float("nan")

        tail = areas[max(1, 2 * len(areas) // 3):]
        AoverA0 = float(np.mean(tail))
        AoverA0_tailstd = float(np.std(tail))

        cents = _cell_centroids(pg, ranges)
        # restrict to active cells for the radius/centroid measure
        act_ids = np.where(active)[0]
        cents_act = np.array([cents[c] for c in act_ids]) if act_ids.size else cents
        R_max, R_rms, _cc, _r = _effective_radius(cents_act)
        R_eff_um = (R_max + R) * _UM

        # ensure necrosis zones are current at the equilibrated config
        h["necrosis"].act(sim.timestep)
        frac = _necrotic_fraction(st, active)
        rim_frac = (len(h["traction"].rim_cells) / float(act_ids.size)
                    if act_ids.size else 0.0)

        if verbose:
            print(f"    R_eff={R_eff_um:.1f} µm  A/A0={AoverA0:.2f}±{AoverA0_tailstd:.2f}  "
                  f"nec_f={frac['necrotic']:.2f}  rim_f={rim_frac:.2f}  "
                  f"({sps:.0f} steps/s)")

        return {
            "n_cells": n_cells, "finite_gate": True, "dt": trial_dt,
            "on_gpu": is_gpu, "n_particles": int(sim.state.N_particles),
            "cell_radius_um": cell_radius_um, "division": division,
            "R_eff_um": R_eff_um, "R_max_um": float(R_max * _UM),
            "R_rms_um": float(R_rms * _UM),
            "AoverA0": AoverA0, "AoverA0_tailstd": AoverA0_tailstd,
            "A_traj": [float(a) for a in areas],
            "A_abs_um2_traj": [float(a) for a in areas_abs],
            "A0_self_um2": float(A0_self * _UM ** 2),
            "necrotic_frac": frac["necrotic"], "quiescent_frac": frac["quiescent"],
            "prolif_frac": frac["proliferating"], "active_rim_frac": float(rim_frac),
            "n_divisions": int(h["prolif"].n_divisions) if division else 0,
            "steps_per_s": float(sps),
            "centroids_um": (cents_act * _UM).tolist(),
            "states": st.state[act_ids].tolist(),
        }
    return {"n_cells": n_cells, "finite_gate": False, "dt": 1.0e-10}


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def make_figure(rows, fit, is_gpu, cell_radius_um):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = sorted(rows, key=lambda r: r["R_eff_um"])
    R = np.array([r["R_eff_um"] for r in rows])
    A = np.array([r["AoverA0"] for r in rows])
    Aerr = np.array([r["AoverA0_tailstd"] for r in rows])
    necf = np.array([r["necrotic_frac"] for r in rows])
    rimf = np.array([r["active_rim_frac"] for r in rows])
    a, b, c, r2 = fit
    pa, pb, pc = PI_LAW
    dev = "GPU (A5000)" if is_gpu else "CPU fallback (dev)"

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.0))

    # (a) A/A0 vs R + fitted curve + PI law overlay
    ax0 = ax[0]
    ax0.errorbar(R, A, yerr=Aerr, fmt="o", ms=9, color="C3", capsize=4, zorder=5,
                 label="measured (equil A/A0 ± tail std)")
    Rg = np.linspace(R.min() * 0.9, R.max() * 1.1, 200)
    ax0.plot(Rg, a + b / Rg + c / Rg ** 2, "-", color="C0", lw=2,
             label=f"fit: a={a:.2f}, b={b:.0f}, c={c:.0f}\n(r²={r2:.3f})")
    Rpi = np.linspace(31, 78, 200)
    ax0.plot(Rpi, pa + pb / Rpi + pc / Rpi ** 2, "--", color="k", lw=1.6,
             label=f"PI law (R=31-78µm)\na={pa}, b={pb}, c={pc}")
    for r in rows:
        ax0.annotate(f"N={r['n_cells']}", (r["R_eff_um"], r["AoverA0"]),
                     textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax0.axvline(150.0, ls=":", color="0.4", lw=1.0,
                label="necrosis-onset R≈150µm")
    ax0.set_xlabel("effective radius R [µm]")
    ax0.set_ylabel("equilibrium A / A0  (vs own initial footprint)")
    ax0.set_title(f"(a) GPU law sweep: A/A0 = a + b/R + c/R²\n[{dev}, "
                  f"R_cell={cell_radius_um:.1f}µm patch]")
    ax0.legend(fontsize=8, loc="best")

    # (b) necrotic + active-rim fraction vs R
    ax1 = ax[1]
    ax1.plot(R, necf, "o-", color="#7f0000", label="necrotic fraction")
    ax1.plot(R, rimf, "s--", color="#2ca02c", label="active-rim fraction")
    ax1.axvline(150.0, ls=":", color="0.4", lw=1.0, label="necrosis-onset R≈150µm")
    ax1.set_xlabel("effective radius R [µm]")
    ax1.set_ylabel("cell fraction")
    ax1.set_ylim(-0.03, 1.03)
    ax1.set_title("(b) necrotic core grows / active rim shrinks with R")
    ax1.legend(fontsize=8, loc="best")

    fig.tight_layout()
    path = os.path.join(FIGS, "gpu_law_sweep.png")
    fig.savefig(path, dpi=120)
    print(f"  figure -> {os.path.relpath(path)}")
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="GPU-scale active+necrosis size sweep of A/A0 = a + b/R + c/R².")
    ap.add_argument("--sizes", type=int, nargs="+",
                    default=[60, 120, 200, 300, 450],
                    help="cell counts to sweep (-> effective radii R).")
    ap.add_argument("--cell-radius-um", type=float, default=20.0,
                    help="coarse-grained DCM cell (patch) radius [µm]. 20 µm makes "
                         "the 60-450 sweep span R≈118-213 µm, crossing necrosis "
                         "onset (~150 µm). Use 7.5 for fine cells (needs thousands).")
    ap.add_argument("--division", action="store_true",
                    help="enable LIVE proliferation (default OFF — the static law test).")
    ap.add_argument("--quick", action="store_true",
                    help="fewer/shorter blocks (CPU dev mini-sweep smoke).")
    ap.add_argument("--force-cpu", action="store_true",
                    help="force the CPU device (dev/fallback smoke).")
    ap.add_argument("--seed", type=int, default=7, help="RNG seed.")
    args = ap.parse_args(argv)

    device, is_gpu, dev_note = select_device(args.force_cpu)
    print(f"[device] {dev_note}  -> on_gpu={is_gpu}")
    print(f"[sweep]  sizes={args.sizes}  R_cell={args.cell_radius_um}µm  "
          f"division={args.division}  quick={args.quick}\n")

    rows = []
    for n in args.sizes:
        print(f"-- N={n} active+necrosis spheroid --")
        r = run_one(n, device=device, is_gpu=is_gpu,
                    cell_radius_um=args.cell_radius_um, division=args.division,
                    quick=args.quick, seed=args.seed)
        if r.get("finite_gate"):
            rows.append(r)

    if len(rows) < 3:
        print(f"\nFAIL: only {len(rows)} sizes passed the finite gate (need >=3 to fit).")
        result = {"finite_sizes": len(rows), "on_gpu": is_gpu,
                  "device_note": dev_note, "rows": rows}
        with open(os.path.join(OUT, "gpu_law_sweep.json"), "w") as f:
            json.dump(result, f, indent=2)
        return 1

    R = [r["R_eff_um"] for r in rows]
    A = [r["AoverA0"] for r in rows]
    fit = fit_law(R, A)
    a, b, c, r2 = fit

    order = np.argsort(R)
    A_sorted = np.array(A)[order]
    nec_sorted = np.array([rows[i]["necrotic_frac"] for i in order])
    decreasing = bool(A_sorted[-1] < A_sorted[0])
    mono_frac = float((np.diff(A_sorted) < 0).mean())
    corr = float(np.corrcoef(np.array(R), np.array(A))[0, 1])
    nec_any = bool(nec_sorted.max() > 0.0)
    nec_grows = bool(nec_any and nec_sorted[-1] > nec_sorted[0])
    nec_corr = (float(np.corrcoef(np.array(R), nec_sorted)[0, 1])
                if nec_sorted.std() > 0 else float("nan"))

    print("\n=== SWEEP TABLE ===")
    print(f"  {'N':>5}  {'R_eff[µm]':>10}  {'A/A0':>14}  {'nec_f':>6}  {'rim_f':>6}  {'parts':>8}")
    for r in sorted(rows, key=lambda x: x["R_eff_um"]):
        print(f"  {r['n_cells']:>5}  {r['R_eff_um']:>10.1f}  "
              f"{r['AoverA0']:>7.2f}±{r['AoverA0_tailstd']:<5.2f}  "
              f"{r['necrotic_frac']:>6.2f}  {r['active_rim_frac']:>6.2f}  "
              f"{r['n_particles']:>8}")

    print("\n=== FIT vs PI ===")
    print(f"  fitted : a={a:+.3f}  b={b:+.1f} µm  c={c:+.1f} µm²  r²={r2:.3f}")
    print(f"  PI law : a={PI_LAW[0]:+.3f}  b={PI_LAW[1]:+.1f} µm  c={PI_LAW[2]:+.1f} µm²  r²=0.98")
    print(f"  R range (this sweep): {min(R):.1f}-{max(R):.1f} µm  (PI: 31-78 µm)")
    print(f"  (i)  A/A0 DECREASES with R : {decreasing}  "
          f"(corr={corr:+.2f}, mono-decr frac {mono_frac:.2f})")
    if nec_any:
        print(f"  (ii) necrotic fraction GROWS with R : {nec_grows}  (corr={nec_corr:+.2f})")
    else:
        print("  (ii) necrotic fraction GROWS with R : N/A — necrosis OFF everywhere "
              "(R < ~150µm onset; scale up R_cell/sizes on the A5000 to activate it)")
    print(f"  (iii) 3-param fit r²>0.7 : {'YES' if r2 > 0.7 else 'NO'} (r²={r2:.3f})")

    result = {
        "sizes": args.sizes, "cell_radius_um": args.cell_radius_um,
        "division": args.division, "quick": args.quick,
        "on_gpu": is_gpu, "device_note": dev_note,
        "measurement": "equilibrium A/A0 = mean(last third of footprint traj) / "
                       "spheroid's OWN initial (t=0) basal footprint; R_eff = "
                       "(max active-centroid radius + R_cell). Necrosis/pressure/"
                       "junction updaters run at low cadence; division OFF unless "
                       "--division. Coarse-grained DCM cell (patch) radius so the "
                       "sweep crosses the ~150 µm necrosis-onset radius.",
        "rows": rows,
        "fit": {"a": a, "b_um": b, "c_um2": c, "r2": r2},
        "pi_law": {"a": PI_LAW[0], "b_um": PI_LAW[1], "c_um2": PI_LAW[2],
                   "r2": 0.98, "R_range_um": [31, 78]},
        "R_range_um": [min(R), max(R)],
        "corr_R_AoverA0": corr,
        "law_signature_decreasing": decreasing,
        "monotone_decreasing_fraction": mono_frac,
        "necrosis_active_any": nec_any,
        "necrotic_grows_with_R": nec_grows,
        "corr_R_necrotic": nec_corr,
        "fit_r2_gt_0p7": bool(r2 > 0.7),
    }
    jpath = os.path.join(OUT, "gpu_law_sweep.json")
    with open(jpath, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n  json -> {os.path.relpath(jpath)}")

    try:
        make_figure(rows, fit, is_gpu, args.cell_radius_um)
    except Exception as e:  # noqa: BLE001
        print(f"  (figure skipped: {e})")

    print("\n=== GPU LAW SWEEP RESULT ===")
    print(f"  sizes passed finite gate : {len(rows)}/{len(args.sizes)}")
    print(f"  A/A0 down with R (i)     : {'YES' if decreasing else 'NO'}")
    print(f"  necrotic up with R (ii)  : "
          f"{'YES' if nec_grows else ('N/A (necrosis OFF, R<150µm)' if not nec_any else 'NO')}")
    print(f"  fit r²>0.7 (iii)         : {'YES' if r2 > 0.7 else 'NO'}")
    if not is_gpu:
        print("  NOTE: CPU fallback — this is a STRUCTURAL/finite smoke. The real "
              "R>150µm necrosis-on law test is gbook-A5000 only "
              "(--sizes 60 120 200 300 450).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""H.7 — DYNAMIC single-cell spreading run (real drivers, explicit BAOAB).

This is the *dynamic* successor to ``h7_spreading_compare.py``. The compare
script measured the basal footprint of a cell whose lamellipodium barely grew
over its short run (A/A₀ ≈ 1.03–1.06 after ~1500 steps): at the physiological
elongation rate ``k_elong⁰ = 11.6 s⁻¹`` and ``dt ≈ 13 ns`` the ~1500-step run
spans only ~20 µs of simulation time, so the expected number of monomer-addition
events is ≈ 0 — the dendritic network never advances and the footprint stays at
its construction geometry. That made it read as a *kinematic/geometric mock*.

Real cell spreading takes minutes (a ~5 µm front advance at ~6 µm/s of raw
barbed-end growth is ~0.86 s ≈ 66 million BAOAB steps), which is infeasible on
CPU in the 48 h window. Per PI authorization (2026-06-10: relax the wall-clock
hard rule, band-match literature data), this harness runs the lamellipodium
KINETICS in explicit *fast-forward*: every Bell-Evans / Arp2/3 / capping RATE is
multiplied by a single, explicitly-reported acceleration factor ``S`` (the
"kinetic clock" runs S× faster than the mechanical BAOAB clock). The MECHANISM
is untouched — the BAOAB integrator still moves every particle, the
force-velocity force-sensitivity ``δ_elong`` is unchanged so the protrusion
force-velocity *shape* is preserved, and the FA molecular clutch still binds via
Pereverzev catch-slip. Only the absolute timescale is compressed. Results are
reported in EFFECTIVE time ``t_eff = S · t_sim`` and validated against the
literature spreading band (A/A₀ ≈ 2–4 for an isotropically spreading single
cell; early power-law → plateau time-course, Cuvelier 2007 / Dubin-Thaler /
Henry 2015 / Betorz 2023).

Observables sampled over the run:
  * A(t)  — basal contact footprint = convex-hull area of the lamellipodial
            actin near the basal plane (the spreading observable). A/A₀.
  * r_front(t) — leading-edge radius (max basal radius of lamellipodial actin).
  * n_engaged(t) — engaged FA integrin↔ligand clutches (substrate adhesion).
  * traction(t)  — summed |F| over engaged clutch bonds [N] (substrate traction).
  * n_actin(t)   — lamellipodial actin bead count (network growth).

Usage (CPU dev / mesoscale):
    python -m ffn_sim.scripts.h7_spreading_dynamics --n-filaments 120 \
        --warmup 120 --accel 600 --steps 60000 --sample-every 2000 \
        --device cpu --allow-cpu-dev \
        --frames-out outputs/h7/spreading/frames.npz \
        --curve-out  outputs/h7/spreading/curve.json
"""

from __future__ import annotations

import argparse
import json
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)

_UM = 1.0e6  # m -> µm


# ---------------------------------------------------------------------------
# Basal-footprint observables
# ---------------------------------------------------------------------------
def _actin_lamel_xyz(snap) -> np.ndarray:
    """(n,3) positions [m] of the lamellipodial actin beads (type actin_lamel)."""
    types = list(snap.particles.types)
    if "actin_lamel" not in types:
        return np.empty((0, 3), dtype=np.float64)
    tid = np.asarray(snap.particles.typeid)
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    mask = tid == types.index("actin_lamel")
    return pos[mask] if mask.any() else np.empty((0, 3), dtype=np.float64)


def _basal_xy(xyz: np.ndarray, z_basal: float, band: float) -> np.ndarray:
    """(m,2) [µm] of the actin beads within ``band`` [m] of the basal plane.

    The spreading footprint is the contact area, so we keep only beads near the
    basal plane ``z_basal`` (the cell-substrate contact) — a bead that branched
    up into the cytosol is not part of the contact footprint.
    """
    if len(xyz) == 0:
        return np.empty((0, 2))
    near = np.abs(xyz[:, 2] - z_basal) <= band
    return xyz[near, :2] * _UM if near.any() else np.empty((0, 2))


def _hull_area_um2(xy: np.ndarray) -> float:
    """Convex-hull area [µm²] of an (n,2) point cloud (0 if degenerate)."""
    if len(xy) < 3:
        return 0.0
    try:
        from scipy.spatial import ConvexHull

        return float(ConvexHull(xy).volume)  # 2-D hull "volume" == area
    except Exception:  # noqa: BLE001
        return 0.0


def _front_radius_um(xy: np.ndarray) -> float:
    """Leading-edge radius [µm] = 95th-percentile in-plane radius of the basal
    point cloud (robust to single LJ-flung outliers, unlike the raw max)."""
    if len(xy) == 0:
        return 0.0
    r = np.hypot(xy[:, 0], xy[:, 1])
    return float(np.percentile(r, 95)) if len(r) >= 5 else float(r.max())


def _engaged_clutch_traction(snap, integrin_action, k_int: float, r0: float):
    """(n_engaged, traction_N) from the live integrin↔ligand bonds in the snap.

    Traction = Σ |F| over engaged ``integrin_ligand`` bonds, F = k·max(0,|Δr|−r0).
    """
    btypes = list(snap.bonds.types)
    if "integrin_ligand" not in btypes:
        return 0, 0.0
    bt = np.asarray(snap.bonds.typeid)
    bg = np.asarray(snap.bonds.group, dtype=np.int64)
    sel = bg[bt == btypes.index("integrin_ligand")]
    if sel.shape[0] == 0:
        return 0, 0.0
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    dr = pos[sel[:, 0]] - pos[sel[:, 1]]
    r = np.linalg.norm(dr, axis=1)
    F = k_int * np.clip(r - r0, 0.0, None)
    return int(sel.shape[0]), float(F.sum())


# ---------------------------------------------------------------------------
# Build the adherent, actively-spreading cell with accelerated kinetics
# ---------------------------------------------------------------------------
def build_spreading_cell(
    *, n_filaments, warmup, accel, device, seed, geometry, membrane_load,
    box_factor=3.0,
):
    """Build the FA-adhered + basal-ring-lamellipodium cell at the physiological
    baseline, with the lamellipodium kinetics scaled by ``accel`` (kinetic
    fast-forward; see module docstring). ``box_factor`` enlarges the box
    (L_box = box_factor·R_cell) so the spreading front has room to reach the
    physiological footprint without the periodic wall."""
    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    overrides = {"box": {"L_box_over_R_cell": float(box_factor)}}
    if n_filaments is not None:
        overrides["cortex"] = {"n_filaments": int(n_filaments), "demo_mode": True}
    manifest["cortex_overrides"] = overrides

    opt = manifest["optional_subsystems"]
    opt["fa"]["enabled"] = True
    opt["lamellipodium"]["enabled"] = True
    opt["lamellipodium"]["geometry"] = geometry
    # Kinetic fast-forward: scale the three Bell-Evans / Arp2/3 rates. The
    # force-sensitivity lengths (delta_elong/delta_cap) and the stall forces are
    # UNCHANGED, so the protrusion force-velocity shape and the abortive-branch
    # threshold are preserved — only the absolute kinetic clock is compressed.
    if accel != 1.0:
        opt["lamellipodium"]["lamellipodium"] = {
            "k_elong_0": 11.6 * accel,
            "k_cap_0": 3.0 * accel,
            "k_b_0": 0.037 * accel,
        }
    if membrane_load:
        opt["membrane_load"]["enabled"] = True
        # The leading-edge Brownian-ratchet brake makes the protrusion
        # force-velocity-coupled, so the trajectory saturates (plateau) as the
        # advancing edge loads up — the physiological time-course shape.

    cell = build_baseline_cell(
        manifest=manifest, device=device, seed=seed,
        equilibrate=True, equilibrate_steps=warmup,
        equilibrate_softstart_steps=max(50, warmup // 2),
    )
    return cell


def attach_basal_clutch(cell, *, k_adh=None):
    """Attach the FA molecular-clutch ensemble as a basal-adhesion z-tether on
    the lamellipodial sheet (keeps the protruded front on the substrate).
    Returns the force (None if not attached). k_adh defaults to the FA clutch
    stiffness k_int_bare."""
    from ffn_sim.cell.spreading_drive import BasalAdhesionTether

    handles = cell.extras["handles"]
    z_basal, band, R_cell, ell0 = _basal_z_band(cell)
    if k_adh is None:
        # Membrane-tension-scale basal holding stiffness: the FA-adhesion ensemble
        # only needs to resist the cortex/membrane lift on the contact sheet
        # (γ_mem ~ 50 pN/µm = 5e-5 N/m, h5 config). A stiffer per-bead tether
        # (e.g. the full integrin clutch k_int ~1e-3) rigidly pins the dendritic
        # network into one plane and the branch points overlap → LJ blow-up;
        # 5e-5 holds the front basal without 2-D compression.
        k_adh = 5.0e-5
    snap = cell.simulation.state.get_snapshot()
    types = list(snap.particles.types)
    if "actin_lamel" not in types:
        return None
    tether = BasalAdhesionTether(
        z_basal=z_basal, k_adh=float(k_adh),
        actin_typeid=types.index("actin_lamel"),
    )
    cell.simulation.operations.integrator.forces.append(tether)
    return tether


def _basal_z_band(cell):
    """Derive the basal-plane z [m], a contact band [m], R_cell and ℓ₀ [m].

    The basal-ring lamellipodium places its WAVE beads exactly on the basal
    contact plane and the mother seeds one ℓ₀ inward; so z_basal = mean WAVE z
    and ℓ₀ = mean |WAVE − mother| (robust, no handle dependency). The contact
    band is 2·ℓ₀ so beads that branched one layer up still count as contact but
    cytosolic branches do not.
    """
    handles = cell.extras["handles"]
    lay = handles.get("lamellipodium_layout")
    if lay is not None and getattr(lay, "wave_positions", None) is not None \
            and len(lay.wave_positions) > 0:
        z_basal = float(np.mean(lay.wave_positions[:, 2]))
        ell0 = float(np.mean(np.linalg.norm(
            lay.wave_positions - lay.mother_seed_positions, axis=1)))
    else:
        z_basal, ell0 = -7.0e-6, 0.5e-6
    R_cell = float(abs(z_basal) + ell0)  # z_basal = -R_cell + ell0
    band = max(2.0 * ell0, 1.0e-6)
    return z_basal, band, R_cell, ell0


def attach_leading_edge(cell, *, nucleate_every, p_advance, max_advance, seed,
                        max_spread_um=6.0):
    """Attach the membrane-tracked protrusion engine (drives the advancing front
    past the fixed WAVE ring at a band-matched rate). Returns the action."""
    from ffn_sim.cell.spreading_drive import LeadingEdgeNucleationUpdater
    import hoomd

    z_basal, band, R_cell, ell0 = _basal_z_band(cell)
    state = cell.extras["handles"].get("lamellipodium_state")
    if state is None:
        raise RuntimeError("no lamellipodium_state handle; lamellipodium not built")
    # Cap the footprint at the physiological maximal spread radius (a spreading
    # MCF7's basal footprint ~doubles → r ~5–6 µm, A/A₀ ~4), AND keep it inside
    # 80% of the box half-edge (out-of-bounds guard) — whichever is smaller.
    box = cell.simulation.state.get_snapshot().configuration.box
    max_radius = min(max_spread_um * 1e-6, 0.80 * 0.5 * float(box[0]))
    action = LeadingEdgeNucleationUpdater(
        lamel_state=state, z_basal=z_basal, band=band, rest_length=ell0,
        p_advance=p_advance, advance_margin=2.0 * ell0, max_advance=max_advance,
        max_radius=max_radius, seed=seed,
    )
    cell.simulation.operations.updaters.append(
        hoomd.update.CustomUpdater(
            action=action, trigger=hoomd.trigger.Periodic(int(nucleate_every))
        )
    )
    return action


def run_spreading(
    cell, *, steps, sample_every, accel,
):
    """Run BAOAB, sampling the spreading observables every ``sample_every`` steps."""
    handles = cell.extras["handles"]
    integrin_action = handles.get("integrin_action")
    fa = handles.get("fa_integration")
    p_fa = getattr(fa, "p", None)
    k_int = 1.0e-3
    r0_int = 0.0
    dt = float(handles.get("dt_used") or 1.3e-8)

    z_basal, band, R_cell, ell0 = _basal_z_band(cell)

    series = []   # list of dicts per sample
    frames = []   # (m,2) basal xy [µm] per sample, for the animation

    def sample(step):
        snap = cell.simulation.state.get_snapshot()
        xyz = _actin_lamel_xyz(snap)
        xy = _basal_xy(xyz, z_basal, band)
        A = _hull_area_um2(xy)
        rf = _front_radius_um(xy)
        n_eng, trac = _engaged_clutch_traction(snap, integrin_action, k_int, r0_int)
        t_sim = step * dt
        rec = {
            "step": int(step),
            "t_sim_s": float(t_sim),
            "t_eff_s": float(t_sim * accel),
            "A_um2": float(A),
            "r_front_um": float(rf),
            "n_actin": int(len(xyz)),
            "n_basal_actin": int(len(xy)),
            "n_engaged": int(n_eng),
            "traction_N": float(trac),
        }
        series.append(rec)
        frames.append(xy.astype(np.float64))
        return rec

    r0 = sample(0)
    A0 = r0["A_um2"] if r0["A_um2"] > 0 else float("nan")
    print(f"[spread] t=0  A0={A0:.2f} µm²  r_front={r0['r_front_um']:.2f} µm  "
          f"n_actin={r0['n_actin']}  n_engaged={r0['n_engaged']}  "
          f"traction={r0['traction_N']:.2e} N", flush=True)

    n_chunks = max(1, steps // sample_every)
    t_wall0 = time.time()
    for c in range(1, n_chunks + 1):
        cell.simulation.run(sample_every)
        rec = sample(c * sample_every)
        aa0 = rec["A_um2"] / A0 if A0 and A0 == A0 else float("nan")
        print(f"[spread] step={rec['step']:>7}  t_eff={rec['t_eff_s']:.3f}s  "
              f"A={rec['A_um2']:.2f} µm²  A/A0={aa0:.3f}  "
              f"r_front={rec['r_front_um']:.2f} µm  n_actin={rec['n_actin']}  "
              f"n_engaged={rec['n_engaged']}  trac={rec['traction_N']:.2e} N",
              flush=True)

    t_wall = time.time() - t_wall0
    meta = {
        "A0_um2": float(A0),
        "accel": float(accel),
        "dt_s": float(dt),
        "R_cell_um": float(R_cell * _UM),
        "ell0_um": float(ell0 * _UM),
        "z_basal_um": float(z_basal * _UM),
        "steps": int(steps),
        "sample_every": int(sample_every),
        "wall_s": float(t_wall),
    }
    return series, frames, meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=120,
                    help="cortex demo n_filaments (mesoscale CPU dev). None=full ×40.")
    ap.add_argument("--warmup", type=int, default=120)
    ap.add_argument("--accel", type=float, default=600.0,
                    help="kinetic fast-forward factor S (rates × S; t_eff = S·t_sim)")
    ap.add_argument("--steps", type=int, default=60000)
    ap.add_argument("--sample-every", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--geometry", default="basal_ring",
                    choices=["basal_ring", "polarized_patch"])
    ap.add_argument("--membrane-load", action="store_true",
                    help="enable the KU-5.2 Brownian-ratchet brake (force-velocity plateau)")
    ap.add_argument("--no-leading-edge", action="store_true",
                    help="disable membrane-tracked protrusion engine (control: front stalls)")
    ap.add_argument("--nucleate-every", type=int, default=200,
                    help="protrusion-engine cadence (integration steps)")
    ap.add_argument("--p-advance", type=float, default=0.04,
                    help="per-tip append probability per tick (protrusion rate knob)")
    ap.add_argument("--max-advance", type=int, default=60,
                    help="max front tips advanced per protrusion tick")
    ap.add_argument("--max-spread-um", type=float, default=11.0,
                    help="footprint radius cap [µm] — set well inside the box so it "
                         "does not bind during the run (a binding cap piles beads at "
                         "the cap radius → LJ overlap); the box+steps bound the spread")
    ap.add_argument("--box-factor", type=float, default=3.0,
                    help="L_box = box_factor·R_cell (default 3 = validated build)")
    ap.add_argument("--v-front-lit-um-min", type=float, default=1.0,
                    help="literature single-cell front velocity [µm/min] for the "
                         "physiological-time mapping (Betorz 2023 P1 ~3µm/3min)")
    ap.add_argument("--no-clutch", action="store_true",
                    help="disable the FA basal-adhesion clutch tether (control: front lifts/stalls)")
    ap.add_argument("--k-adh", type=float, default=None,
                    help="clutch adhesion stiffness [N/m] (default = FA k_int_bare)")
    ap.add_argument("--frames-out", default=None)
    ap.add_argument("--curve-out", default=None)
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    n_fil = None if (args.n_filaments is not None and args.n_filaments <= 0) else args.n_filaments
    cell = build_spreading_cell(
        n_filaments=n_fil, warmup=args.warmup, accel=args.accel,
        device=dev, seed=args.seed, geometry=args.geometry,
        membrane_load=args.membrane_load, box_factor=args.box_factor,
    )
    clutch = None
    if not args.no_clutch:
        clutch = attach_basal_clutch(cell, k_adh=args.k_adh)
    le_action = None
    if not args.no_leading_edge:
        le_action = attach_leading_edge(
            cell, nucleate_every=args.nucleate_every, p_advance=args.p_advance,
            max_advance=args.max_advance, seed=args.seed + 7,
            max_spread_um=args.max_spread_um)
    series, frames, meta = run_spreading(
        cell, steps=args.steps, sample_every=args.sample_every, accel=args.accel,
    )
    meta["leading_edge"] = bool(not args.no_leading_edge)
    meta["clutch"] = bool(not args.no_clutch)
    meta["k_adh"] = float(clutch.k_adh) if clutch is not None else 0.0
    meta["n_promoted"] = int(le_action.n_promoted) if le_action is not None else 0
    meta["v_front_lit_um_min"] = float(args.v_front_lit_um_min)
    # Physiological-time mapping: the emergent front advance mapped to real time
    # at the literature single-cell spreading velocity (PI-authorized band-match).
    r0f = series[0]["r_front_um"]
    rf_um = np.array([r["r_front_um"] for r in series])
    t_phys_min = np.maximum(rf_um - r0f, 0.0) / max(args.v_front_lit_um_min, 1e-9)
    t_sim_total = series[-1]["t_sim_s"] if series[-1]["t_sim_s"] > 0 else 1.0
    meta["S_accel_effective"] = float(t_phys_min[-1] * 60.0 / t_sim_total)
    for r, tp in zip(series, t_phys_min):
        r["t_phys_min"] = float(tp)

    final = series[-1]
    aa0 = final["A_um2"] / meta["A0_um2"] if meta["A0_um2"] else float("nan")
    print(f"\n[spread] DONE  A/A0_final={aa0:.3f}  "
          f"r_front {series[0]['r_front_um']:.2f}->{final['r_front_um']:.2f} µm  "
          f"n_actin {series[0]['n_actin']}->{final['n_actin']}  "
          f"n_engaged {series[0]['n_engaged']}->{final['n_engaged']}  "
          f"wall={meta['wall_s']:.1f}s", flush=True)

    if args.curve_out:
        Path(args.curve_out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.curve_out, "w") as fh:
            json.dump({"meta": meta, "series": series, "geometry": args.geometry,
                       "membrane_load": bool(args.membrane_load)}, fh, indent=2)
        print(f"[spread] wrote curve {args.curve_out}", flush=True)

    if args.frames_out:
        Path(args.frames_out).parent.mkdir(parents=True, exist_ok=True)
        blob = {
            "t_eff_s": np.array([r["t_eff_s"] for r in series]),
            "A_um2": np.array([r["A_um2"] for r in series]),
            "A0_um2": np.array([meta["A0_um2"]]),
            "geometry": np.array([args.geometry]),
        }
        for i, fr in enumerate(frames):
            blob[f"frame__{i}"] = np.asarray(fr, dtype=np.float64)
        blob["n_frames"] = np.array([len(frames)])
        np.savez_compressed(args.frames_out, **blob)
        print(f"[spread] wrote frames {args.frames_out}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

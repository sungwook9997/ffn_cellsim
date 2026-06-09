"""H.7 active-γ ARCHITECTURE step-1 — load-path vs connectivity construction sweep.

active-γ saga closed: the floor is generation/TRANSMISSION-bound. The transmission
wall (WALL A in ``h7_active_force_budget.py``) is the actin-network channel carrying
only ~28 % of the myosin-dipole γ, because the contraction load-path between anchors
is ~1 backbone segment (investigation #2: 73 % of inter-anchor spans are
single-segment, mean 1.52 seg). Chugh 2017 / Truong Quang 2021: real cortical tension
is set by actin-filament-length / actin-myosin OVERLAP — i.e. by a LONG load-path —
largely independent of myosin number.

PI fork (2026-06-09): *"아키텍처 갈아엎어야 되는 거 아닌가?"* — is the cortex
CONSTRUCTION the lever? This probe answers the prerequisite **without any contraction
run (no binder cost)**: across the construction knobs that set the load-path
(``z_struct`` anchor density, ``L_long_mean`` formin backbone length, ``bundle_mult``
bundling, ``formin_fraction``), can the mesh produce a LONG load-path (mean inter-anchor
span ≫ 1 segment, low single-segment fraction) while KEEPING connectivity
(giant-component ≥ 0.9)?

  - If YES (a knob region gives long spans AND giant ≥ 0.9): the short load-path is a
    construction choice, not a topological necessity → architecture IS the lever →
    proceed to the loading-phase WALL-A confirmation, then the rebuild.
  - If NO (long span ALWAYS fragments the mesh): connectivity and load-path are
    coupled in the point-crosslink construction → the rebuild must move connectivity
    onto BUNDLE/overlap (Chugh-style) rather than point crosslinks, or it is a deeper
    representation limit → surface to PI before a big build.

This is a STATIC mesh-topology measurement (seed only, ``with_simulation=False``):
milliseconds-to-seconds per knob set, no HOOMD integrator, no myosin updater. It does
NOT change any production config — the production cortex stays at its ratified knobs;
this is mechanism scoping, not tuning (band LOCKED, no gate touched).

Usage:
    python -m ffn_sim.scripts.h7_loadpath_architecture_sweep --n-filaments 300
    python -m ffn_sim.scripts.h7_loadpath_architecture_sweep --n-filaments 300 \
        --z-list 1.6,2.0,2.4,2.8,3.2 --llong-list 3,5,8,12
"""
from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

import numpy as np

from ffn_sim.cell.manifest import load_manifest, _deep_merge
from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cortex.connected_mesh import build_connected_cortex

_UM = 1.0e6


def _resolve_params(n_filaments):
    """Resolve (p_cortex, p_xl) on the SAME config path resolve_baseline uses, so
    the swept mesh is the production MCF7 cortex topology (knobs aside)."""
    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    base_cfg = load_manifest(manifest["base_cortex_config"])
    cfg = _deep_merge(base_cfg, manifest.get("cortex_overrides"))
    cfg.setdefault("cortex", {})["R_cell"] = float(manifest["R_cell"])
    if n_filaments is not None:
        cfg["cortex"]["n_filaments"] = int(n_filaments)
        cfg["cortex"]["demo_mode"] = True
    p_cortex = resolve_h3_derived(cfg)
    p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)
    return p_cortex, p_xl


def _spans_from_anchors(filament_starts, n_beads_per_filament, anchored_set, ell0):
    """Inter-anchor free-span distribution along each filament contour.

    A filament's beads are contiguous & contour-ordered:
    ``[start, start+n)``. Anchored beads partition the chain into spans; a span
    between anchors at contour positions p0<p1 covers (p1-p0) backbone segments.
    Returns the per-span segment counts (load-path lengths) as an array.
    """
    seg_counts = []
    F = len(filament_starts)
    for i in range(F):
        s = int(filament_starts[i])
        n = int(n_beads_per_filament[i])
        anchors = [p for p in range(n) if (s + p) in anchored_set]
        for p0, p1 in zip(anchors[:-1], anchors[1:]):
            if p1 - p0 > 0:
                seg_counts.append(p1 - p0)
    return np.asarray(seg_counts, dtype=np.int64)


def _measure(handles, ell0):
    """Static topology metrics for one built seed."""
    seed = handles.seed
    layout = handles.layout
    fstart = np.asarray(layout.filament_starts)
    nbpf = np.asarray(layout.n_beads_per_filament)

    # crosslink-anchored actin beads (column 1 of seeded_attach = actin bead tag).
    sa = np.asarray(seed.seeded_attach)
    xl_anchor = set(int(b) for b in sa[:, 1]) if sa.size else set()
    # branch-junction beads (mother + daughter-base) also tie filaments → anchors.
    bb = np.asarray(layout.branch_bonds)
    branch_anchor = set()
    if bb.size:
        branch_anchor = set(int(x) for x in bb.reshape(-1))
    all_anchor = xl_anchor | branch_anchor

    spans_xl = _spans_from_anchors(fstart, nbpf, xl_anchor, ell0)
    spans_all = _spans_from_anchors(fstart, nbpf, all_anchor, ell0)

    def _stats(spans):
        if spans.size == 0:
            return dict(n_spans=0, mean_seg=None, frac_1seg=None,
                        frac_2plus=None, mean_path_um=None)
        return dict(
            n_spans=int(spans.size),
            mean_seg=float(spans.mean()),
            frac_1seg=float((spans == 1).mean()),
            frac_2plus=float((spans >= 2).mean()),
            mean_path_um=float(spans.mean() * ell0 * _UM),
        )

    return {
        "giant_fraction": float(seed.giant_fraction),
        "z_struct_realised": float(seed.z_struct_realised),
        "L_over_lc": float(seed.L_over_lc),
        "n_xl": int(seed.n_xl),
        "n_homeless": int(seed.n_homeless),
        "n_cortex_actin": int(handles.n_cortex_actin),
        "spans_crosslink_only": _stats(spans_xl),   # comparable to investigation #2
        "spans_all_anchors": _stats(spans_all),     # transmission load-path (xl+branch)
    }


def run_sweep(*, n_filaments, z_list, llong_list, bundle_list, formin_list, seed):
    p_cortex, p_xl = _resolve_params(n_filaments)
    ell0 = float(p_cortex.rest_length)
    F = int(n_filaments) if n_filaments else int(p_cortex.n_filaments)

    print(f"  ℓ0 = {ell0*1e9:.0f} nm   F = {F} filaments   R = {p_cortex.R_cell*_UM:.2f} µm",
          flush=True)
    print("=" * 100, flush=True)
    print(f"  {'z':>4} {'Llong':>6} {'bnd':>4} {'fmn':>5} | "
          f"{'giant':>6} {'z_real':>6} {'L/lc':>5} | "
          f"{'spans':>6} {'mean_seg':>8} {'1seg%':>6} {'≥2seg%':>6} {'path_µm':>7}",
          flush=True)
    print("-" * 100, flush=True)

    rows = []
    for z in z_list:
        for llong in llong_list:
            for bnd in bundle_list:
                for fmn in formin_list:
                    rng = np.random.default_rng(seed)
                    handles = build_connected_cortex(
                        p_cortex, p_xl,
                        formin_fraction=fmn,
                        L_long_mean=llong * 1e-6,
                        z_struct=z,
                        bundle_mult=int(bnd),
                        n_filaments=F,
                        with_baoab=False,
                        with_simulation=False,
                        rng=rng,
                    )
                    m = _measure(handles, ell0)
                    sx = m["spans_all_anchors"]
                    rec = dict(z_struct=z, L_long_um=llong, bundle_mult=int(bnd),
                               formin_fraction=fmn, **m)
                    rows.append(rec)
                    print(f"  {z:>4.1f} {llong:>6.1f} {bnd:>4d} {fmn:>5.2f} | "
                          f"{m['giant_fraction']*100:>5.1f}% {m['z_struct_realised']:>6.2f} "
                          f"{m['L_over_lc']:>5.1f} | "
                          f"{(sx['n_spans'] or 0):>6d} "
                          f"{(sx['mean_seg'] or 0):>8.2f} "
                          f"{(sx['frac_1seg'] or 0)*100:>5.1f}% "
                          f"{(sx['frac_2plus'] or 0)*100:>5.1f}% "
                          f"{(sx['mean_path_um'] or 0):>7.3f}", flush=True)
    print("=" * 100, flush=True)
    return rows, dict(ell0_nm=ell0 * 1e9, n_filaments=F, R_cell_um=p_cortex.R_cell * _UM)


def _verdict(rows):
    """Name the conclusion: is there a long-load-path + connected region?"""
    print("  READING — architecture-as-lever prerequisite:", flush=True)
    conn = [r for r in rows if r["giant_fraction"] >= 0.9]
    if not conn:
        print("   ⚠ NO knob set kept giant-component ≥ 0.9 at this scale "
              "(connectivity floor not met) — widen z or check scale.", flush=True)
        return
    # the baseline (production-ish) and the longest-load-path connected point
    best = max(conn, key=lambda r: (r["spans_all_anchors"]["mean_seg"] or 0))
    sb = best["spans_all_anchors"]
    base = min(rows, key=lambda r: abs(r["z_struct"] - 2.8) + abs(r["L_long_um"] - 5.0))
    sbase = base["spans_all_anchors"]
    print(f"   baseline  (z={base['z_struct']}, Llong={base['L_long_um']}µm): "
          f"giant={base['giant_fraction']*100:.0f}%  mean_seg={sbase['mean_seg']:.2f}  "
          f"1seg={sbase['frac_1seg']*100:.0f}%", flush=True)
    print(f"   best-conn (z={best['z_struct']}, Llong={best['L_long_um']}µm): "
          f"giant={best['giant_fraction']*100:.0f}%  mean_seg={sb['mean_seg']:.2f}  "
          f"1seg={sb['frac_1seg']*100:.0f}%  path={sb['mean_path_um']:.2f}µm", flush=True)
    ratio = (sb["mean_seg"] or 0) / (sbase["mean_seg"] or 1)
    if (sb["mean_seg"] or 0) >= 2.0 and best["giant_fraction"] >= 0.9:
        print(f"   → LONG load-path ({sb['mean_seg']:.1f} seg, {ratio:.1f}× baseline) IS "
              f"achievable WHILE connected → architecture is a candidate lever.", flush=True)
        print(f"     NEXT: loading-phase WALL-A force-budget at this knob set "
              f"(does g_actin/g_myo rise?).", flush=True)
    else:
        print(f"   → load-path stays short even at the connectivity edge "
              f"(best mean_seg={sb['mean_seg']:.1f}) → point-crosslink connectivity and "
              f"load-path are COUPLED; rebuild must move connectivity onto bundle/overlap. "
              f"Surface to PI.", flush=True)
    print("=" * 100, flush=True)


def _figure(rows, meta, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
    zs = sorted(set(r["z_struct"] for r in rows))
    llongs = sorted(set(r["L_long_um"] for r in rows))

    ax = axes[0]
    for ll in llongs:
        sub = sorted([r for r in rows if r["L_long_um"] == ll], key=lambda r: r["z_struct"])
        ax.plot([r["z_struct"] for r in sub],
                [r["giant_fraction"] * 100 for r in sub], "o-",
                label=f"Llong={ll}µm")
    ax.axhline(90, color="green", ls="--", lw=1.2, label="giant ≥ 90% (connected)")
    ax.set_xlabel("z_struct (anchor density)"); ax.set_ylabel("giant-component (%)")
    ax.set_title("connectivity vs anchor density"); ax.legend(fontsize=8)

    ax = axes[1]
    for ll in llongs:
        sub = sorted([r for r in rows if r["L_long_um"] == ll], key=lambda r: r["z_struct"])
        ax.plot([r["z_struct"] for r in sub],
                [(r["spans_all_anchors"]["mean_seg"] or 0) for r in sub], "s-",
                label=f"Llong={ll}µm")
    ax.axhline(2.0, color="purple", ls="--", lw=1.2, label="≥2 seg (multi-segment path)")
    ax.set_xlabel("z_struct (anchor density)")
    ax.set_ylabel("mean inter-anchor span (segments = load-path)")
    ax.set_title("load-path vs anchor density")
    ax.legend(fontsize=8)
    fig.suptitle(f"H.7 architecture step-1: load-path vs connectivity "
                 f"(F={meta['n_filaments']}, ℓ0={meta['ell0_nm']:.0f}nm)",
                 fontweight="bold")
    fig.savefig(out_png, dpi=130)
    print(f"  figure → {out_png}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=300)
    ap.add_argument("--z-list", type=str, default="1.6,2.0,2.4,2.8,3.2",
                    help="comma list of z_struct (anchor density) values")
    ap.add_argument("--llong-list", type=str, default="3,5,8,12",
                    help="comma list of L_long_mean values in µm")
    ap.add_argument("--bundle-list", type=str, default="3",
                    help="comma list of bundle_mult values")
    ap.add_argument("--formin-list", type=str, default="0.12",
                    help="comma list of formin_fraction values")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-json", type=str,
                    default="ffn_sim/outputs/h7/production/h7_loadpath_architecture_sweep.json")
    ap.add_argument("--out-png", type=str,
                    default="ffn_sim/outputs/h7/figs/h7_loadpath_architecture_sweep.png")
    args = ap.parse_args()

    z_list = [float(x) for x in args.z_list.split(",")]
    llong_list = [float(x) for x in args.llong_list.split(",")]
    bundle_list = [int(x) for x in args.bundle_list.split(",")]
    formin_list = [float(x) for x in args.formin_list.split(",")]

    rows, meta = run_sweep(
        n_filaments=args.n_filaments, z_list=z_list, llong_list=llong_list,
        bundle_list=bundle_list, formin_list=formin_list, seed=args.seed,
    )
    _verdict(rows)

    out = {"meta": meta, "note": "STATIC mesh-topology sweep (seed only); no production "
           "config change; mechanism scoping for architecture-as-lever (band LOCKED).",
           "rows": rows}
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2))
    print(f"  json → {args.out_json}", flush=True)
    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)
    _figure(rows, meta, args.out_png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""H.7 — compare the two single-cell lamellipodium spreading geometries.

Builds the full physiological MCF7 cell (FA-adhered, lamellipodium ON, membrane
OFF for step-1) with EACH leading-edge geometry and measures the basal spreading
footprint A(t)/A0 + the WAVE azimuthal distribution, so the two models can be
compared against the PI's experimental A/A0 data (the model must match data to
be usable, PI 2026-06-07):

  * basal_ring      — single-cell ISOTROPIC spreading (circumferential rim);
                      footprint grows ~radially-symmetric.
  * polarized_patch — single-cell MIGRATING (one leading edge about p̂);
                      footprint extends DIRECTIONALLY along p̂.

This is the comparison TOOL + a pipeline SMOKE. The authoritative A/A0 vs the PI
experimental overlay needs the full production scale on GPU (gbook) over a long
spreading run; a small --n-filaments run is labelled a smoke.

Usage (smoke):
    python -m ffn_sim.scripts.h7_spreading_compare --n-filaments 120 \
        --warmup 80 --chunks 4 --chunk-steps 150 --device cpu --allow-cpu-dev
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest  # noqa: E402
from ffn_sim.common.production_policy import (  # noqa: E402
    add_production_device_args,
    validate_production_device_args,
)

_GEOMS = ("basal_ring", "polarized_patch")
_UM = 1.0e6


def _basal_actin_xy(cell):
    """(x, y) [µm] of the lamellipodium actin beads (the spreading structure)."""
    snap = cell.simulation.state.get_snapshot()
    types = list(snap.particles.types)
    tid = np.asarray(snap.particles.typeid)
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    names = [t for t in types if "lamel" in t.lower() and "actin" in t.lower()]
    if not names:  # fall back to the WAVE mother actin type if naming differs
        names = [t for t in types if "actin_lamel" in t]
    mask = np.zeros(len(pos), dtype=bool)
    for nm in names:
        mask |= tid == types.index(nm)
    return pos[mask, :2] * _UM if mask.any() else np.empty((0, 2))


def _footprint_area(xy: np.ndarray) -> float:
    """Convex-hull area [µm²] of an (n,2) basal point cloud (0 if degenerate)."""
    if len(xy) < 3:
        return 0.0
    try:
        from scipy.spatial import ConvexHull

        return float(ConvexHull(xy).volume)  # 2-D hull "volume" == area
    except Exception:  # noqa: BLE001
        return 0.0


def _polarization_index(xy: np.ndarray) -> float:
    """Circular resultant length of the basal azimuths: 0 = isotropic ring
    (directions cancel), -> 1 = tightly polarized patch (branch-cut robust,
    unlike a max-min span)."""
    if len(xy) == 0:
        return 0.0
    phi = np.arctan2(xy[:, 1], xy[:, 0])
    return float(np.hypot(np.cos(phi).mean(), np.sin(phi).mean()))


def run_geometry(geom, *, n_filaments, warmup, chunks, chunk_steps, device, seed):
    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    if n_filaments is not None:
        manifest["cortex_overrides"] = {
            "cortex": {"n_filaments": int(n_filaments), "demo_mode": True}
        }
    manifest["optional_subsystems"]["fa"]["enabled"] = True
    manifest["optional_subsystems"]["lamellipodium"]["enabled"] = True
    manifest["optional_subsystems"]["lamellipodium"]["geometry"] = geom

    cell = build_baseline_cell(
        manifest=manifest, device=device, seed=seed,
        equilibrate=True, equilibrate_steps=warmup,
        equilibrate_softstart_steps=max(50, warmup // 2),
    )
    xy0 = _basal_actin_xy(cell)
    a0 = _footprint_area(xy0)
    series = [(0, a0, 1.0)]
    frames = [xy0]                       # per-chunk basal-actin xy for the animation
    for c in range(1, chunks + 1):
        cell.simulation.run(chunk_steps)
        xy = _basal_actin_xy(cell)
        a = _footprint_area(xy)
        series.append((c * chunk_steps, a, (a / a0 if a0 > 0 else float("nan"))))
        frames.append(xy)
    xy_final = frames[-1]
    return {
        "geom": geom,
        "n_wave": int(cell.extras["handles"].get("n_wave_particles") or 0),
        "n_clutch": int(cell.extras["handles"].get("n_fa_clutch_bonds") or 0),
        "polarization_index": _polarization_index(xy0),
        "A0_um2": a0,
        "series": series,        # [(step, A_um2, A/A0)]
        "xy0": xy0, "xy_final": xy_final,
        "frames": frames,        # list of (n,2) basal xy per chunk (cell-movement anim)
    }


def _figure(results, out_path: Path):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), constrained_layout=True)
    colors = {"basal_ring": "#2f6fb0", "polarized_patch": "#e23b3b"}
    # Panels 1-2: basal (x,y) footprint per geometry (initial grey + final colour)
    for ax, r in zip(axes[:2], results):
        ax.scatter(r["xy0"][:, 0], r["xy0"][:, 1], s=6, c="0.7", label="t=0")
        ax.scatter(r["xy_final"][:, 0], r["xy_final"][:, 1], s=6,
                   c=colors[r["geom"]], alpha=0.6, label="t=final")
        ax.set_aspect("equal"); ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
        ax.set_title(f"{r['geom']}\nbasal footprint (polarization "
                     f"{r['polarization_index']:.2f}; 0=isotropic→1=polar)", fontsize=10)
        ax.legend(fontsize=8)
    # Panel 3: A/A0 over the spreading run
    ax = axes[2]
    for r in results:
        steps = [s for s, _, _ in r["series"]]
        aa0 = [v for _, _, v in r["series"]]
        ax.plot(steps, aa0, "-o", color=colors[r["geom"]], label=r["geom"])
    ax.set_xlabel("sample step"); ax.set_ylabel("A / A0 (basal footprint)")
    ax.set_title("spreading footprint A/A0\n(SMOKE — full A/A0 vs PI data = gbook GPU)",
                 fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.suptitle("H.7 — single-cell lamellipodium spreading geometries "
                 "(isotropic rim vs polarized patch), membrane OFF step-1",
                 fontsize=12, fontweight="bold")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=None)
    ap.add_argument("--warmup", type=int, default=120)
    ap.add_argument("--chunks", type=int, default=4)
    ap.add_argument("--chunk-steps", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument(
        "--out",
        default=str(Path(__file__).resolve().parents[1]
                    / "outputs" / "h7" / "figs" / "h7_spreading_compare.png"),
    )
    ap.add_argument("--frames-out", default=None,
                    help="npz path to save per-chunk basal-xy frames (for the cell-movement animation)")
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    results = []
    for geom in _GEOMS:
        r = run_geometry(
            geom, n_filaments=args.n_filaments, warmup=args.warmup,
            chunks=args.chunks, chunk_steps=args.chunk_steps, device=dev, seed=args.seed,
        )
        results.append(r)
        scale = "FULL x40" if args.n_filaments is None else f"smoke n_fil={args.n_filaments}"
        print(f"[{geom}] ({scale}) WAVE={r['n_wave']} polarization={r['polarization_index']:.2f} "
              f"clutches={r['n_clutch']} A0={r['A0_um2']:.2f} µm² "
              f"A/A0_final={r['series'][-1][2]:.3f}", flush=True)

    _figure(results, Path(args.out))
    print(f"[h7-spread] wrote {args.out}", flush=True)

    if args.frames_out:
        blob = {"geoms": np.array([r["geom"] for r in results]),
                "n_frames": np.array([len(r["frames"]) for r in results], dtype=int),
                "chunk_steps": np.array([args.chunk_steps]),
                "R_cell_um": np.array([results[0].get("R_cell_um", 7.5)])}
        for r in results:
            for i, fr in enumerate(r["frames"]):
                blob[f"{r['geom']}__{i}"] = np.asarray(fr, dtype=np.float64)
        Path(args.frames_out).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.frames_out, **blob)
        print(f"[h7-spread] wrote frames {args.frames_out} "
              f"(render with scripts/h7_spreading_anim.py)", flush=True)

    print("PROVISIONAL smoke — authoritative A/A0 vs PI experimental data needs the "
          "full production scale on gbook GPU over a long spreading run.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

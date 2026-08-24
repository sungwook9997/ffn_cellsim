#!/usr/bin/env python3
r"""Activate the head-resolved NMII motor on the assembled cell — attach -> power-stroke -> Bell detach.

Drives the KMC that ``ac.motor`` provides (``MyosinForce.step_kinetics``: attach free heads in capture range,
advance the walked abscissa by tau*Hill(load), detach with the Bell slip rate) coupled to the inner
mechanical solve, so the bound-head fraction EMERGES from k_on/k_off (never an imposed duty) and a contractile
cortical tension develops. Reports the magnitude HONESTLY — the master §6.2 density-floor is expected and is
NOT closed here by adding heads/density.

The attach kernel needs, per head, the nearest actin node + its distance (``nearest_dist``, ``nearest_idx``):
a free head binds only within ``capture_radius``. We compute that on the host with a KD-tree over the actin
nodes (out-of-hot-loop, refreshed every ``--nn-every`` ticks as the cortex creeps). If the reference
minifilament geometry places the heads OUTSIDE capture of every actin node (the I4-weave per-head actin
hand-off is not wired in this first assembled cell), NO head can bind — that is reported as a finding, not
patched by loosening the capture radius.

STERIC is off by default so the motor's contractile signal is isolated. Runtime: Warp-CUDA only (I0-A).

    ~/miniconda3/envs/ffn_sim/bin/python -m aleph.components.incumbent.myosin_activate --n-filaments 70686 --ticks 40
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import _accumulate_all, step_myosin_kinetics
from aleph.laws.network_warp import axpy_kernel, reshape_kernel


def _nearest_actin(actin_pos: np.ndarray, head_pos: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nearest actin node index + distance for each head (KD-tree; numpy chunked fallback)."""
    try:
        from scipy.spatial import cKDTree
        d, idx = cKDTree(actin_pos).query(head_pos, k=1)
        return d.astype(np.float64), idx.astype(np.int32)
    except Exception:  # noqa: BLE001 — fallback if scipy missing
        d = np.empty(head_pos.shape[0], np.float64)
        idx = np.empty(head_pos.shape[0], np.int32)
        for s in range(0, head_pos.shape[0], 4096):
            hp = head_pos[s:s + 4096]
            dd = np.linalg.norm(actin_pos[None, :, :] - hp[:, None, :], axis=2)
            idx[s:s + 4096] = dd.argmin(1)
            d[s:s + 4096] = dd.min(1)
        return d, idx


def _myosin_only_force(cell, pos: wp.array) -> np.ndarray:
    """Per-actin-node |F| from the MyosinForce primitive ALONE (fresh scratch; no state mutation) [pN]."""
    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    cell.myosin.accumulate(pos, f)
    wp.synchronize_device(cell.device)
    return np.linalg.norm(f.numpy()[: cell.n_actin], axis=1)


def _radial_contractile_force(cell, pos: wp.array, com: np.ndarray) -> float:
    """Net INWARD (contractile) radial component of the myosin force on actin nodes [pN] (sum)."""
    f = wp.zeros(cell.n_total, dtype=wp.vec3d, device=cell.device)
    cell.myosin.accumulate(pos, f)
    wp.synchronize_device(cell.device)
    fa = f.numpy()[: cell.n_actin]
    pa = pos.numpy()[: cell.n_actin]
    rhat = pa - com
    rhat /= (np.linalg.norm(rhat, axis=1, keepdims=True) + 1e-30)
    fr = np.sum(fa * rhat, axis=1)          # + outward, - inward
    return float(-fr.sum())                 # net inward (contractile) magnitude


def run_myosin(cfg: CellConfig, *, ticks: int = 40, tau: float = 2.0e-3, inner_iters: int = 20,
               nn_every: int = 10, seed: int = 0, out_json: str | None = None) -> dict:
    """Build the cell, bind/walk/detach the NMII heads for ``ticks`` KMC steps; report bound fraction + tension."""
    t0 = time.time()
    cell = build_cell(cfg)
    if cell.myosin is None:
        raise RuntimeError("no myosin composed (do not pass --no-myosin).")
    dev = cell.device
    myo = cell.myosin
    n_heads = int(myo.state["bound"].shape[0])
    pos = cell.pos_d
    com = pos.numpy()[: cell.n_actin].mean(axis=0)

    head_node = myo.head_node.numpy().astype(np.int64)

    def _refresh_nn() -> dict:
        allp = pos.numpy()
        d, idx = _nearest_actin(allp[: cell.n_actin], allp[head_node])
        return {"d": d, "idx": idx,
                "nd": wp.array(d, dtype=wp.float64, device=dev),
                "ni": wp.array(idx, dtype=wp.int32, device=dev)}

    nn = _refresh_nn()
    cap = float(myo.params.capture_radius)
    within = int((nn["d"] <= cap).sum())
    diag = {
        "n_heads": n_heads, "capture_radius_um": cap,
        "head_actin_dist_min_um": float(nn["d"].min()), "head_actin_dist_p50_um": float(np.median(nn["d"])),
        "head_actin_dist_p05_um": float(np.percentile(nn["d"], 5)),
        "heads_within_capture": within, "heads_within_capture_frac": within / max(n_heads, 1),
        "heads_within_2x_capture": int((nn["d"] <= 2 * cap).sum()),
        "heads_within_0p2um": int((nn["d"] <= 0.2).sum()),
    }

    # ── KMC loop coupled to the inner mechanical solve ───────────────────────────────────────────────
    bound_frac_series: list[float] = []
    tension_series: list[float] = []
    for tk in range(ticks):
        with wp.ScopedDevice(dev):
            for it in range(inner_iters):                     # relax mechanics at the current binding state
                _accumulate_all(cell, pos, cell.f_d)
                wp.launch(axpy_kernel, dim=cell.n_total,
                          inputs=[pos, wp.float64(cell.dt_mu), cell.f_d], device=dev)
            if cell.n_fibers:
                wp.launch(reshape_kernel, dim=cell.n_fibers,
                          inputs=[pos, cell.foff_d, cell.soff_d, cell.srest_d, wp.int32(2)], device=dev)
            wp.synchronize_device(dev)
        myo.compute_loads(pos)                                 # per-head crossbridge tension for Hill/Bell
        # attach -> fill_walk_dir (I4 hand-off: walk the REAL actin polarity) -> step -> detach (device RNG)
        step_myosin_kinetics(cell, nn["nd"], nn["ni"], tau, seed + tk)
        bf = float(myo.state["bound"].numpy().mean())
        bound_frac_series.append(bf)
        tension_series.append(_radial_contractile_force(cell, pos, com))
        if nn_every and (tk + 1) % nn_every == 0:
            nn = _refresh_nn()

    fmyo = _myosin_only_force(cell, pos)
    bound_final = float(myo.state["bound"].numpy().mean())
    n_bound = int(myo.state["bound"].numpy().sum())
    report = {
        **diag,
        "n_filaments": cfg.n_filaments, "n_actin": cell.n_actin, "n_total": cell.n_total,
        "ticks": ticks, "tau_s": tau, "inner_iters": inner_iters,
        "k_on": float(myo.params.k_on), "k_off0": float(myo.params.k_off0),
        "bound_fraction_final": bound_final, "n_bound_heads": n_bound,
        "bound_fraction_series": bound_frac_series,
        "expected_bound_frac_if_in_range": float(myo.params.k_on) / (float(myo.params.k_on) + float(myo.params.k_off0)),
        "net_inward_contractile_force_pN": tension_series[-1] if tension_series else 0.0,
        "myosin_only_force_max_pN": float(fmyo.max()), "myosin_only_force_mean_pN": float(fmyo.mean()),
        "wall_s": time.time() - t0,
    }
    if out_json:
        with open(out_json, "w") as fh:
            json.dump(report, fh, indent=2, default=float)
    return report


def _print(r: dict) -> None:
    print("\n=== NMII activation (attach -> power-stroke -> Bell detach) ===", flush=True)
    print(f"  n_filaments={r['n_filaments']:,}  n_actin={r['n_actin']:,}  n_heads={r['n_heads']:,}", flush=True)
    print(f"  BINDING GEOMETRY: capture_radius={r['capture_radius_um']} um  |  "
          f"head->actin dist min={r['head_actin_dist_min_um']:.4f} p05={r['head_actin_dist_p05_um']:.4f} "
          f"p50={r['head_actin_dist_p50_um']:.4f} um", flush=True)
    print(f"    heads within capture={r['heads_within_capture']:,} "
          f"({100*r['heads_within_capture_frac']:.2f}%)  within 2x={r['heads_within_2x_capture']:,}  "
          f"within 0.2um={r['heads_within_0p2um']:,}", flush=True)
    print(f"  KINETICS: k_on={r['k_on']}/s  k_off0={r['k_off0']}/s  -> bound frac if in range ~"
          f"{r['expected_bound_frac_if_in_range']:.3f}", flush=True)
    print(f"  RESULT: bound fraction={r['bound_fraction_final']:.4g}  (n_bound={r['n_bound_heads']:,})", flush=True)
    print(f"    net inward contractile force={r['net_inward_contractile_force_pN']:.4g} pN   "
          f"myosin-only |F|: max={r['myosin_only_force_max_pN']:.4g} mean={r['myosin_only_force_mean_pN']:.4g} pN",
          flush=True)
    print(f"  wall={r['wall_s']:.1f}s", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description="Activate the NMII motor on the assembled cell.")
    p.add_argument("--n-filaments", type=int, default=70686)
    p.add_argument("--ticks", type=int, default=40)
    p.add_argument("--tau", type=float, default=2.0e-3, help="KMC tick [s]")
    p.add_argument("--inner-iters", type=int, default=20)
    p.add_argument("--nn-every", type=int, default=10)
    p.add_argument("--no-steric", action="store_true")
    p.add_argument("--device", type=str, default=None, help="Warp CUDA device alias (default: current CUDA)")
    p.add_argument("--out", type=str, default="")
    a = p.parse_args()
    cfg = CellConfig(n_filaments=a.n_filaments, with_myosin=True, with_steric=not a.no_steric,
                     with_pressure=True, device=a.device)
    print(f"[myo] building composed cell n_filaments={cfg.n_filaments} steric={cfg.with_steric}", flush=True)
    r = run_myosin(cfg, ticks=a.ticks, tau=a.tau, inner_iters=a.inner_iters, nn_every=a.nn_every,
                   out_json=a.out or None)
    _print(r)
    if a.out:
        print(f"[myo] wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()

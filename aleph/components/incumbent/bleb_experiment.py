#!/usr/bin/env python3
r"""Oracle-B blebbing experiment on the assembled Active Cell — apply local de-adhesion, read the nucleation criterion.

Wires :mod:`ac.cell.bleb_perturbation` into the REAL composed cell (``build_cell``), the way ``fsi_demo`` /
``myosin_activate`` reuse the frozen builder without touching the shared driver loop. It:

  1. builds the full cell (cortex + membrane + ERM tethers + Biot pressure at Π₀),
  2. de-adheres the membrane–cortex ERM tethers inside an angular cap (the Charras ablation), commit-safe,
  3. reads the EMERGENT nucleation criterion from sourced parameters — no tuning.

Nucleation physics (Charras 2005/2008). With the ERM linkers removed over a cap of geodesic radius
``a = R_cell·sin(θ½)``, the resting pore pressure ΔP pushes the bare bilayer outward; the membrane's own
tension γ_mem resists via Laplace ``ΔP_Laplace = 2 γ_mem / a``. A bleb nucleates iff the *measured* local
ΔP exceeds that resistance, i.e. iff the de-adhered patch is larger than the critical radius

    a_crit = 2 γ_mem / ΔP      →      nucleates ⇔ a > a_crit .

Both γ_mem (KB-3.B1.1) and the pore pressure ΔP (≈Π₀=40 Pa, measured here from the live field, not assumed)
are sourced, so ``a_crit`` is a prediction, not a fit. This is the c_v-independent nucleation observable of
the blebbing validation plan; the full multi-step emergent growth + ΔP non-equilibration (~10 µm / ~10 s) is
the next step and needs the mechanical outer-loop.

Runtime: Warp-CUDA only (I0-A). Runs on the gbook A5000.

    ~/miniconda3/envs/ffn_sim/bin/python -m aleph.components.incumbent.bleb_experiment \
        --n-filaments 70686 --patch-deg 25 --out aleph/outputs/ac/bleb/bleb.npz
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.bleb_perturbation import (
    cap_cos_half_angle,
    count_tethers_in_cap,
    erm_patch_deadhesion_kernel,
)
from aleph.components.fluid.field_grid import FLUID

__all__ = ["run_bleb_experiment"]


def _cap_pressure_jump(grid, centroid: np.ndarray, axis: np.ndarray, cos_half: float, p_ext: float) -> tuple:
    """Mean/median (p − p_ext) over the FLUID cells whose direction from ``centroid`` lies in the cap [Pa]."""
    p = grid.p.numpy()
    mask = grid.mask.numpy()
    nx, ny, nz = grid.shape
    ii, jj, kk = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    xyz = np.stack([ii, jj, kk], -1) * grid.dx + np.asarray(grid.origin)
    r = xyz - centroid
    d = np.linalg.norm(r, axis=-1)
    safe = d > 1e-12
    cos = np.full(d.shape, -2.0)
    cos[safe] = (r[safe] @ axis) / d[safe]
    sel = (mask == FLUID) & (cos >= cos_half)
    if not np.any(sel):
        return 0.0, 0.0, 0
    dp = p[sel] - p_ext
    return float(dp.mean()), float(np.median(dp)), int(sel.sum())


def run_bleb_experiment(cfg: CellConfig, *, patch_deg: float = 25.0, axis: np.ndarray | None = None,
                        out_npz: str | None = None) -> dict:
    """Apply a commit-safe local ERM de-adhesion to the composed cell and report the nucleation criterion."""
    t0 = time.time()
    cell = build_cell(cfg)
    dev = cell.device
    mem, grid = cell.membrane, cell.grid
    if mem is None or grid is None:
        raise RuntimeError("bleb experiment needs the membrane + Biot substrate (no --no-membrane/--no-pressure).")
    if mem.n_erm == 0:
        raise RuntimeError("no ERM tethers on the membrane — cannot de-adhere.")

    ax = np.array([0.0, 0.0, 1.0]) if axis is None else np.asarray(axis, np.float64)
    ax = ax / np.linalg.norm(ax)
    cos_half = cap_cos_half_angle(patch_deg)

    pos = cell.pos_d.numpy()
    centroid = pos[: cell.n_actin].mean(axis=0)
    mem_idx = mem.erm_m_d.numpy()
    mem_pos = pos[mem_idx]                                   # membrane endpoint of every tether
    bound0 = mem.erm_bound_d.numpy()
    R_cell = float(np.linalg.norm(mem_pos - centroid, axis=1).mean())

    n_target = count_tethers_in_cap(mem_pos, bound0, centroid, ax, cos_half)

    # ── apply the commit-safe de-adhesion on the REAL device state (accepted = 1) ─────────────────────
    with wp.ScopedDevice(dev):
        accepted = wp.ones(1, dtype=wp.int32, device=dev)
        wp.launch(erm_patch_deadhesion_kernel, dim=mem.n_erm,
                  inputs=[cell.pos_d, mem.erm_m_d, mem.erm_bound_d,
                          wp.vec3d(*centroid), wp.vec3d(*ax), wp.float64(cos_half), accepted], device=dev)
    bound1 = mem.erm_bound_d.numpy()
    n_deadhered = int(((bound0 == 1) & (bound1 == 0)).sum())

    # ── nucleation criterion from sourced γ_mem + the MEASURED local pore pressure ────────────────────
    p_ext = float(cell.substrate.p_ext) if hasattr(cell.substrate, "p_ext") else 0.0
    dp_mean, dp_med, n_cells = _cap_pressure_jump(grid, centroid, ax, cos_half, p_ext)
    a_patch = R_cell * np.sin(np.deg2rad(patch_deg))         # geodesic cap radius [µm]
    dp = dp_mean if dp_mean > 0 else float(getattr(cell.substrate, "p_ext", 0.0))
    a_crit = (2.0 * mem.gamma_mem / dp) if dp > 0 else float("inf")   # 2 γ_mem / ΔP [µm]
    dp_laplace = 2.0 * mem.gamma_mem / a_patch if a_patch > 0 else float("inf")

    report = {
        "n_filaments": cfg.n_filaments, "n_erm_total": mem.n_erm, "n_erm_bound_pre": int(bound0.sum()),
        "patch_half_angle_deg": patch_deg, "patch_axis": ax.tolist(),
        "n_tethers_deadhered": n_deadhered, "n_tethers_predicted": n_target,
        "deadhered_matches_prediction": bool(n_deadhered == n_target),
        "frac_membrane_deadhered": n_deadhered / max(int(bound0.sum()), 1),
        "R_cell_um": R_cell, "gamma_mem_pN_um": mem.gamma_mem, "f_rupt_pN": mem.f_rupt, "k_erm_pN_um": mem.k_erm,
        "adhesion_force_removed_pN": n_deadhered * mem.f_rupt,
        "dP_local_mean_Pa": dp_mean, "dP_local_median_Pa": dp_med, "n_fluid_cells_in_cap": n_cells,
        "a_patch_um": float(a_patch), "a_crit_um": float(a_crit), "dP_laplace_Pa": float(dp_laplace),
        "nucleates": bool(a_patch > a_crit),
        "wall_s": time.time() - t0,
    }
    if out_npz:
        np.savez_compressed(out_npz, bound_pre=bound0.astype(np.int8), bound_post=bound1.astype(np.int8),
                            mem_pos=mem_pos.astype(np.float32), centroid=centroid.astype(np.float64),
                            patch_axis=ax, report_json=np.array(json.dumps(report), dtype=object))
        with open(str(out_npz).rsplit(".npz", 1)[0] + ".json", "w") as fh:
            json.dump(report, fh, indent=2, default=float)
    return report


def _print(r: dict) -> None:
    print("\n=== Oracle-B bleb nucleation (local ERM de-adhesion on the composed cell) ===", flush=True)
    print(f"  n_filaments={r['n_filaments']:,}  ERM tethers={r['n_erm_total']:,} (bound {r['n_erm_bound_pre']:,})",
          flush=True)
    print(f"  patch: half-angle={r['patch_half_angle_deg']}°  axis={np.round(r['patch_axis'],3).tolist()}  "
          f"-> de-adhered {r['n_tethers_deadhered']:,} tethers "
          f"(predicted {r['n_tethers_predicted']:,}, match={r['deadhered_matches_prediction']})", flush=True)
    print(f"  removed adhesion = {r['adhesion_force_removed_pN']:.4g} pN  "
          f"(f_rupt={r['f_rupt_pN']:.3g} pN, k_erm={r['k_erm_pN_um']:.4g} pN/µm)", flush=True)
    print(f"  MEASURED local ΔP = {r['dP_local_mean_Pa']:.3g} Pa (median {r['dP_local_median_Pa']:.3g}) over "
          f"{r['n_fluid_cells_in_cap']} cap fluid cells;  γ_mem={r['gamma_mem_pN_um']:.3g} pN/µm", flush=True)
    print(f"  nucleation: a_patch={r['a_patch_um']:.3g} µm  vs  a_crit=2γ/ΔP={r['a_crit_um']:.3g} µm  "
          f"(ΔP_Laplace={r['dP_laplace_Pa']:.3g} Pa)", flush=True)
    print(f"  ==> BLEB NUCLEATES: {r['nucleates']}   (wall={r['wall_s']:.1f}s)", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description="Oracle-B bleb nucleation via local ERM de-adhesion.")
    p.add_argument("--n-filaments", type=int, default=70686)
    p.add_argument("--patch-deg", type=float, default=25.0, help="ablation cap half-angle [deg]")
    p.add_argument("--axis", type=float, nargs=3, default=None, help="cap axis (default +z)")
    p.add_argument("--no-myosin", action="store_true")
    p.add_argument("--no-steric", action="store_true")
    p.add_argument("--device", type=str, default=None)
    p.add_argument("--out", type=str, default="")
    a = p.parse_args()
    cfg = CellConfig(n_filaments=a.n_filaments, with_myosin=not a.no_myosin, with_steric=not a.no_steric,
                     with_membrane=True, with_pressure=True, device=a.device)
    print(f"[bleb] building composed cell n_filaments={cfg.n_filaments}", flush=True)
    r = run_bleb_experiment(cfg, patch_deg=a.patch_deg, axis=a.axis, out_npz=a.out or None)
    _print(r)
    if a.out:
        print(f"[bleb] wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()

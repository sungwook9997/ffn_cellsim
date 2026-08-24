#!/usr/bin/env python3
r"""Dump the COUPLED active+fluid run state (positions + pressure p + v_f + per-node force) for viz.

Runs ONE coupled outer step on the assembled cell with the FSI back-coupling LIVE (a controlled cortex
compression -> div(v_s) -> pore-pressure response) and, optionally, the NMII motor bound, then writes the
full state a follow-up dynamics visualiser needs:

  * actin node positions (post-step)                            [um]
  * the Eulerian pore-pressure field p                          [Pa]
  * the reconstructed pore-fluid velocity v_f = v_s + q/phi     [um/s]  (phi PROVISIONAL — I0-B1b GAP)
  * the Darcy discharge q                                       [um/s]
  * the solid dilatation-rate source div(v_s)                   [1/s]
  * the composed per-node force |F| (bending+xlink+myosin+steric+pressure) [pN]

STERIC is off by default so the coupled poroelastic + motor signal is isolated from the sigma_EV confound.
Reuses ``build_cell`` + ``SolidDilatationCoupling`` + ``DarcyVelocity`` read-only; edits no track / ff module.
Runtime: Warp-CUDA only (I0-A). Runs on the gbook A5000; the Mac renders the npz.
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import _accumulate_all, _residual_host, step_myosin_kinetics
from aleph.components.incumbent.fsi_coupling import SolidDilatationCoupling
from aleph.components.incumbent.fsi_demo import _sub_cycle_biot
from aleph.components.incumbent.myosin_activate import _nearest_actin
from aleph.components.fluid.field_grid import FLUID
from aleph.components.fluid.velocity import DarcyVelocity
from aleph.laws.network_warp import axpy_kernel, reshape_kernel

PHI_PROVISIONAL = 0.5  # I0-B1b GAP


def dump_coupled(cfg: CellConfig, out_npz: str, *, eps: float = 0.01, dt_phys: float = 0.05,
                 activate_myosin: bool = True, myo_ticks: int = 20, phi: float = PHI_PROVISIONAL) -> dict:
    """Advance one coupled active+fluid step and dump positions + p + v_f + q + div_vs + per-node force."""
    t0 = time.time()
    cell = build_cell(cfg)
    dev, g, sub = cell.device, cell.grid, cell.substrate
    n_actin, n_total = cell.n_actin, cell.n_total
    pos = cell.pos_d
    com = pos.numpy()[:n_actin].mean(axis=0)

    # ── optional: bind + walk the NMII heads for a few KMC ticks (native-density active motor) ────────
    n_bound = 0
    bound_frac = 0.0
    if activate_myosin and cell.myosin is not None:
        myo = cell.myosin
        head_node = myo.head_node.numpy().astype(np.int64)
        allp = pos.numpy()
        d, idx = _nearest_actin(allp[:n_actin], allp[head_node])
        nd = wp.array(d, dtype=wp.float64, device=dev)
        ni = wp.array(idx, dtype=wp.int32, device=dev)
        for tk in range(myo_ticks):
            with wp.ScopedDevice(dev):
                for _ in range(12):
                    _accumulate_all(cell, pos, cell.f_d)
                    wp.launch(axpy_kernel, dim=n_total, inputs=[pos, wp.float64(cell.dt_mu), cell.f_d], device=dev)
                if cell.n_fibers:
                    wp.launch(reshape_kernel, dim=cell.n_fibers,
                              inputs=[pos, cell.foff_d, cell.soff_d, cell.srest_d, wp.int32(2)], device=dev)
                wp.synchronize_device(dev)
            myo.compute_loads(pos)
            step_myosin_kinetics(cell, nd, ni, 2.0e-3, tk)   # attach → I4 walk_dir hand-off → step/detach
        n_bound = int(myo.state["bound"].numpy().sum())
        bound_frac = n_bound / max(myo.state["bound"].shape[0], 1)

    # ── FSI: compress the cortex over dt_phys, spread v_s -> div(v_s), step the Biot substrate ────────
    pos_old = pos.numpy().copy()
    pos_new = com + (1.0 - eps) * (pos_old - com)
    vs = (pos_new - pos_old) / dt_phys
    active = np.concatenate([np.zeros(n_actin, np.int32), np.full(n_total - n_actin, -1, np.int32)])
    pos.assign(np.ascontiguousarray(pos_new, np.float64))
    vs_d = wp.array(np.ascontiguousarray(vs, np.float64), dtype=wp.vec3d, device=dev)
    active_d = wp.array(active, dtype=wp.int32, device=dev)
    coupling = SolidDilatationCoupling(g)
    coupling.update(pos, vs_d, cell.node_volume_d, active_d)
    n_sub = _sub_cycle_biot(sub, cell.membrane_bc, dt_phys, with_membrane=True)

    # ── reconstruct the pore-fluid velocity v_f = v_s + q/phi on the grid ────────────────────────────
    dv = DarcyVelocity(g, mobility=sub.mobility, phi=phi)
    dv.discharge()
    with wp.ScopedDevice(dev):
        dv.pore_velocity(coupling._vs_grid)                 # v_f = v_s_grid + q/phi
    q = dv.q.numpy()
    v_f = dv.v_f.numpy()
    p_field = g.p.numpy()
    div_vs = g.div_vs.numpy()
    mask = g.mask.numpy()
    fluid = mask == FLUID

    # ── composed per-node force at the coupled state ─────────────────────────────────────────────────
    res, fmag_all = _residual_host(cell, pos, cell.f_d)
    fmag = fmag_all[: cell.n_actin]
    pos_post = pos.numpy()

    qmag = np.linalg.norm(q, axis=-1)[fluid]
    vfmag = np.linalg.norm(v_f, axis=-1)[fluid]
    report = {
        "n_filaments": cfg.n_filaments, "n_actin": n_actin, "n_total": n_total,
        "with_steric": cfg.with_steric, "with_myosin": cfg.with_myosin, "activate_myosin": activate_myosin,
        "n_bound_heads": n_bound, "bound_fraction": bound_frac,
        "eps_compression": eps, "dt_phys": dt_phys, "n_biot_subcycles": n_sub,
        "c_v_um2_s": sub.c_v, "phi_provisional_GAP": phi,
        "p_excess_max_Pa": float(np.abs(p_field[fluid] - 40.0).max()),
        "max_q_um_s": float(qmag.max()), "max_vf_um_s": float(vfmag.max()),
        "composed_force_max_pN": float(res), "composed_force_mean_pN": float(fmag.mean()),
        "wall_s": time.time() - t0,
    }
    np.savez_compressed(
        out_npz,
        actin_pos=pos_post[:n_actin].astype(np.float32),
        f_node_mag=fmag.astype(np.float32),
        p_field=p_field.astype(np.float32), div_vs=div_vs.astype(np.float32),
        q_field=q.astype(np.float32), v_f_field=v_f.astype(np.float32),
        mask=mask.astype(np.int32), grid_shape=np.array(g.shape, np.int32),
        grid_dx=np.float64(g.dx), grid_origin=np.array(g.origin, np.float64),
        fiber_offsets=cell.foff_d.numpy().astype(np.int32),
        report_json=np.array(json.dumps(report), dtype=object),
        ledger_json=np.array(json.dumps(cell.ledger, default=float), dtype=object),
    )
    with open(str(out_npz).rsplit(".npz", 1)[0] + ".json", "w") as fh:
        json.dump({"report": report, "ledger": cell.ledger}, fh, indent=2, default=float)
    print(f"[coupled] wrote {out_npz}  (p_excess_max={report['p_excess_max_Pa']:.4g} Pa  "
          f"max|v_f|={report['max_vf_um_s']:.4g} um/s  bound={n_bound}  |F|max={res:.4g} pN  "
          f"wall={report['wall_s']:.1f}s)", flush=True)
    return report


def main() -> None:
    p = argparse.ArgumentParser(description="Dump the coupled active+fluid run state for viz.")
    p.add_argument("--n-filaments", type=int, default=70686)
    p.add_argument("--eps", type=float, default=0.01)
    p.add_argument("--dt-phys", type=float, default=0.05)
    p.add_argument("--no-steric", action="store_true")
    p.add_argument("--no-activate-myosin", action="store_true")
    p.add_argument("--phi", type=float, default=PHI_PROVISIONAL)
    p.add_argument("--device", type=str, default=None, help="Warp CUDA device alias (default: current CUDA)")
    p.add_argument("--out", type=str, required=True)
    a = p.parse_args()
    cfg = CellConfig(n_filaments=a.n_filaments, with_myosin=True, with_steric=not a.no_steric,
                     with_pressure=True, device=a.device)
    print(f"[coupled] building cell n_filaments={cfg.n_filaments} steric={cfg.with_steric}", flush=True)
    dump_coupled(cfg, a.out, eps=a.eps, dt_phys=a.dt_phys, activate_myosin=not a.no_activate_myosin, phi=a.phi)


if __name__ == "__main__":
    main()

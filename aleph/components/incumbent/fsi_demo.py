#!/usr/bin/env python3
r"""FSI div(v_s) coupling DEMONSTRATION on the assembled Active Cell — a deforming solid drives the fluid.

Proves the wire built in :mod:`ac.cell.fsi_coupling` is LIVE (no longer stubbed at 0): a controlled
deformation of the composed cortex produces a nonzero solid dilatation rate ``div(v_s)`` on the field grid,
and stepping the conservative Biot substrate with it makes the pore pressure RESPOND. The control (div_vs=0,
the old stub) shows the pressure does NOT respond — so the difference is exactly the FSI back-coupling.

STERIC is turned OFF (``--no-steric``) so the pN-scale poroelastic coupling is isolated from the sigma_EV
interpenetration confound (~63k interpenetrating nodes, up to 1000 pN steric) on the dense builder cortex.

Load: a small uniform radial compression of every node about the actin centroid over one outer step
``dt_phys`` (``pos_new = com + (1-eps)(pos-com)``), giving a nodal solid velocity ``v_s = (pos_new-pos_old)/
dt_phys``. The mass-weighted Peskin spread + Eulerian divergence (fsi_coupling) turn that into ``div(v_s) <
0`` (contraction) inside the cell, so ``-alpha*div(v_s) > 0`` pressurises the pore fluid — the physically
correct sign (squeezing the skeleton raises pore pressure).

Reuses the frozen ``build_cell`` (assemble) + ``BiotSubstrate`` read-only; edits no track / ff module.
Runtime: Warp-CUDA only (I0-A). Runs on the gbook A5000.

    ~/miniconda3/envs/ffn_sim/bin/python -m aleph.components.incumbent.fsi_demo \
        --n-filaments 70686 --no-steric --eps 0.01 --out aleph/outputs/ac/fsi/fsi_demo.npz
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.fsi_coupling import SolidDilatationCoupling
from aleph.components.fluid.field_grid import FLUID

__all__ = ["run_fsi_demo"]


def _sub_cycle_biot(substrate, membrane_bc, dt_phys: float, *, with_membrane: bool, safety: float = 0.9) -> int:
    """Advance the Biot substrate over ``dt_phys`` at its CFL (div_vs held constant); return n_subcycles."""
    dt_cfl = substrate.cfl_dt(safety)
    n_sub = max(1, int(dt_phys / dt_cfl) + (1 if dt_phys % dt_cfl > 0 else 0))
    dt_sub = dt_phys / n_sub
    for _ in range(n_sub):
        if with_membrane and membrane_bc is not None:
            membrane_bc.apply(0.0)
        substrate.step(dt_sub)
    return n_sub


def run_fsi_demo(cfg: CellConfig, *, eps: float = 0.01, dt_phys: float = 0.05,
                 with_membrane: bool = True, out_npz: str | None = None) -> dict:
    """Compress the composed cortex by ``eps`` over one outer step; drive the fluid; report the response."""
    t0 = time.time()
    cell = build_cell(cfg)
    dev = cell.device
    g, sub, mbc = cell.grid, cell.substrate, cell.membrane_bc
    if g is None or sub is None:
        raise RuntimeError("FSI demo requires the Biot substrate (do not pass --no-pressure).")
    n_actin, n_total = cell.n_actin, cell.n_total
    n_myo = n_total - n_actin

    pos_old = cell.pos_d.numpy().copy()
    com = pos_old[:n_actin].mean(axis=0)
    pos_new = com + (1.0 - eps) * (pos_old - com)                 # uniform radial compression (the load)
    vs = (pos_new - pos_old) / dt_phys                            # nodal solid velocity [um/s]
    active = np.concatenate([np.zeros(n_actin, np.int32), np.full(n_myo, -1, np.int32)])

    cell.pos_d.assign(np.ascontiguousarray(pos_new, np.float64))  # move the solid (state.pos aliases pos_d)
    vs_d = wp.array(np.ascontiguousarray(vs, np.float64), dtype=wp.vec3d, device=dev)
    active_d = wp.array(active, dtype=wp.int32, device=dev)

    coupling = SolidDilatationCoupling(g)
    coupling.update(cell.pos_d, vs_d, cell.node_volume_d, active_d)
    div_host = coupling.div_vs_host()
    mask = g.mask.numpy()
    fluid = mask == FLUID
    div_f = div_host[fluid]
    cell_vol = g.dx**3
    div_integral = float(div_f.sum() * cell_vol)                  # ~ -eps * V_enclosed (total dilatation rate)

    # ── FSI-ON: step the Biot substrate with the live div(v_s) over one outer dt_phys ────────────────
    p_before = g.p.numpy().copy()
    n_sub = _sub_cycle_biot(sub, mbc, dt_phys, with_membrane=with_membrane)
    p_after = g.p.numpy().copy()
    dp_on = (p_after - p_before)[fluid]                          # pressure response to the deforming solid

    # ── CONTROL (the OLD stub): reset p, zero div_vs, step again -> pressure must NOT respond ────────
    g.set_pressure(p_before)                                     # restore the pre-step field
    with wp.ScopedDevice(dev):
        g.div_vs.zero_()
    _sub_cycle_biot(sub, mbc, dt_phys, with_membrane=with_membrane)
    dp_stub = (g.p.numpy() - p_before)[fluid]

    report = {
        "n_actin": n_actin, "n_total": n_total, "n_filaments": cfg.n_filaments,
        "with_steric": cfg.with_steric, "with_myosin": cfg.with_myosin,
        "eps_compression": eps, "dt_phys": dt_phys, "n_biot_subcycles": n_sub,
        "c_v_um2_s": sub.c_v, "alpha": sub.alpha, "storage_S": sub.storage_S, "mobility": sub.mobility,
        "div_vs_min_per_s": float(div_f.min()), "div_vs_max_per_s": float(div_f.max()),
        "div_vs_mean_per_s": float(div_f.mean()),
        "div_vs_integral_um3_per_s": div_integral,
        "expected_div_vs_continuum_per_s": -3.0 * eps / dt_phys,   # div of uniform radial compression = -3*eps/dt
        "expected_dp_undrained_Pa": sub.alpha * 3.0 * eps / sub.storage_S,  # alpha*|div|*dt/S = 3*alpha*eps/S
        "dp_FSI_ON_max_Pa": float(np.abs(dp_on).max()), "dp_FSI_ON_mean_Pa": float(dp_on.mean()),
        "dp_FSI_ON_p95_Pa": float(np.percentile(dp_on, 95)),
        "dp_CONTROL_stub_max_Pa": float(np.abs(dp_stub).max()),
        "pressure_responds_to_solid": bool(np.abs(dp_on).max() > 100.0 * max(np.abs(dp_stub).max(), 1e-12)),
        "wall_s": time.time() - t0,
    }
    if out_npz:
        np.savez_compressed(
            out_npz,
            div_vs=div_host.astype(np.float32), dp_on=(p_after - p_before).astype(np.float32),
            mask=mask.astype(np.int32), grid_shape=np.array(g.shape, np.int32),
            grid_dx=np.float64(g.dx), grid_origin=np.array(g.origin, np.float64),
            report_json=np.array(json.dumps(report), dtype=object),
        )
        with open(str(out_npz).rsplit(".npz", 1)[0] + ".json", "w") as fh:
            json.dump(report, fh, indent=2, default=float)
    return report


def _print(r: dict) -> None:
    print("\n=== FSI div(v_s) coupling demo (deforming solid -> pore fluid) ===", flush=True)
    print(f"  cortex n_filaments={r['n_filaments']:,}  n_actin={r['n_actin']:,}  n_total={r['n_total']:,}  "
          f"steric={r['with_steric']}", flush=True)
    print(f"  load: uniform radial compression eps={r['eps_compression']} over dt_phys={r['dt_phys']}s  "
          f"(c_v={r['c_v_um2_s']:.1f} um^2/s, alpha={r['alpha']})", flush=True)
    print(f"  div(v_s):  min={r['div_vs_min_per_s']:.4g}  max={r['div_vs_max_per_s']:.4g}  "
          f"mean={r['div_vs_mean_per_s']:.4g} /s   (0 == the old stub;  continuum expect "
          f"{r['expected_div_vs_continuum_per_s']:.4g} /s)", flush=True)
    print(f"  pore-pressure RESPONSE:  FSI-ON max|dp|={r['dp_FSI_ON_max_Pa']:.4g} Pa  "
          f"mean={r['dp_FSI_ON_mean_Pa']:.4g} Pa  (undrained expect ~{r['expected_dp_undrained_Pa']:.4g} Pa)  "
          f"vs   CONTROL(stub) max|dp|={r['dp_CONTROL_stub_max_Pa']:.4g} Pa", flush=True)
    print(f"  ==> a DEFORMING solid drives the pore fluid: {r['pressure_responds_to_solid']}", flush=True)
    print(f"  wall={r['wall_s']:.1f}s", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description="FSI div(v_s) coupling demo on the assembled cell.")
    p.add_argument("--n-filaments", type=int, default=70686)
    p.add_argument("--eps", type=float, default=0.01, help="radial compression fraction (the load)")
    p.add_argument("--dt-phys", type=float, default=0.05)
    p.add_argument("--no-steric", action="store_true", help="isolate the poroelastic signal (recommended)")
    p.add_argument("--no-myosin", action="store_true")
    p.add_argument("--no-membrane", action="store_true", help="skip the membrane efflux BC")
    p.add_argument("--device", type=str, default=None, help="Warp CUDA device alias (default: current CUDA)")
    p.add_argument("--out", type=str, default="")
    a = p.parse_args()
    cfg = CellConfig(n_filaments=a.n_filaments, with_myosin=not a.no_myosin, with_steric=not a.no_steric,
                     with_pressure=True, device=a.device)
    print(f"[fsi] building composed cell n_filaments={cfg.n_filaments} steric={cfg.with_steric} "
          f"myosin={cfg.with_myosin}", flush=True)
    r = run_fsi_demo(cfg, eps=a.eps, dt_phys=a.dt_phys, with_membrane=not a.no_membrane,
                     out_npz=a.out or None)
    _print(r)
    if a.out:
        print(f"[fsi] wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()

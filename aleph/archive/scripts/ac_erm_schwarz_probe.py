"""Validate the ERM-pair Schwarz block preconditioner against the 2026-07-22d resting-baseline blocker.

The FD probe proved the analytic operator is EXACT; the blocker is the coupled soft-membrane/stiff-cortex
CONDITIONING (the exact projected-Newton step overshoots — native 8× — and a damped Newton stalls at 0.776 pN,
>> the 0.21 gate). This probe checks the fork-2 fix (`ac/cell/erm_schwarz.ERMPairSchwarz`):

  1. SPD / descent: the one-sweep direction ``P M^-1 PF`` is a genuine descent (``PF·dx > 0``) and finite.
  2. NO overshoot: line-search scan ``t·dx → max|PF|`` (Newton overshot 8× at t=1; ERM-Schwarz should not).
  3. Descends PAST the 0.776 stall: iterated damped ERM-Schwarz vs the damped-Newton reference, to the gate.

Diagnostic; launches no committed state change. Warp-CUDA (I0-A) → gbook A5000.
  PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim python scripts/ac_erm_schwarz_probe.py [--native]
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import _preload_erm_resting_balance, _residual_host
from aleph.components.incumbent.erm_schwarz import ERMPairSchwarz
from aleph.components.incumbent.implicit_mechanics import (
    ProjectedAnalyticCG,
    compute_regularization_kernel,
    omitted_regularization_base,
)
from aleph.components.incumbent.inner_mechanics import max_abs_pressure_kernel
from aleph.scripts.ac_fd_operator_probe import project, true_total_force
from aleph.engine.gate_criteria import (
    NATIVE_L500_DT_MU,
    NATIVE_L500_TOLERANCE_UM,
    predicted_force_floor,
)

# D4 (PI 2026-07-28): the gate floor is DERIVED from the run's own convergence predicate,
# never the hand-copied 0.21 literal. These archived figures/probes score data from the
# ell=0.5 um native run, so they use that run's recorded tolerance/step. F_pred scales with
# the mesh spacing, so a finer rung must NOT be scored against this value.
GATE = predicted_force_floor(
    inner_tolerance_um=NATIVE_L500_TOLERANCE_UM, inner_dt_mu=NATIVE_L500_DT_MU
)


def _vec(pos_np, device):
    return wp.array(np.ascontiguousarray(pos_np, np.float64), dtype=wp.vec3d, device=device)


def _pf_max(cell, pos_np):
    return float(np.linalg.norm(project(cell, true_total_force(cell, _vec(pos_np, cell.device))), axis=1).max())


def _regularization(cell):
    d = cell.device
    omitted_base = float(omitted_regularization_base(cell))
    max_pressure_d = wp.zeros(1, dtype=wp.float64, device=d)
    if cell.grid is not None and cell.membrane_pressure is not None:
        wp.launch(max_abs_pressure_kernel, dim=cell.grid.shape,
                  inputs=[cell.grid.p, cell.grid.mask, wp.float64(cell.membrane_pressure.p_ext)],
                  outputs=[max_pressure_d], device=d)
    reg_d = wp.zeros(1, dtype=wp.float64, device=d)
    wp.launch(compute_regularization_kernel, dim=1,
              inputs=[wp.float64(omitted_base), max_pressure_d, wp.float64(cell.pressure_edge_um), reg_d], device=d)
    wp.synchronize_device(d)
    return reg_d


def _damped_descent(cell, pos0, direction_fn, reg_d, n_steps):
    """Iterate a host-side backtracking-damped descent using `direction_fn(pos_d, PF_d) -> dx (N,3)`."""
    d = cell.device
    finite_d = wp.ones(1, dtype=wp.int32, device=d)
    active_d = wp.ones(1, dtype=wp.int32, device=d)
    pos = pos0.copy()
    traj = []
    for step in range(n_steps):
        PF = project(cell, true_total_force(cell, _vec(pos, d)))
        res = float(np.linalg.norm(PF, axis=1).max())
        if res <= GATE:
            traj.append({"step": step, "res": res, "t": 0.0, "converged": True})
            break
        dx = direction_fn(_vec(pos, d), _vec(PF, d), active_d, finite_d)
        best_t, best_res = 0.0, res
        for t in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125):
            m = _pf_max(cell, pos + t * dx)
            if m < best_res:
                best_res, best_t = float(m), t
        if best_t == 0.0:
            traj.append({"step": step, "res": res, "t": 0.0, "converged": False, "stalled": True})
            break
        pos = pos + best_t * dx
        traj.append({"step": step, "res": res, "t": best_t, "converged": False})
    return traj, pos


def run(nf, subdiv, n_steps):
    wp.init()
    cfg = CellConfig(n_filaments=nf, overlap_free_cortex=True, erm_radial_pairing=True,
                     membrane_subdivisions=subdiv)
    cell = build_cell(cfg)
    d = cell.device
    _preload_erm_resting_balance(cell)
    pos0 = cell.pos_d.numpy().copy()
    res0, _ = _residual_host(cell, cell.pos_d, cell.f_d)
    reg_d = _regularization(cell)
    a_val = float(reg_d.numpy()[0])

    erm = ERMPairSchwarz(cell)
    ws = ProjectedAnalyticCG(cell, max_iterations=max(256, 8 * subdiv), coarse_iterations=0, coarse_modes=0)

    def erm_dir(pos_d, PF_d, active_d, finite_d):
        return erm.direction(pos_d, PF_d, reg_d, active_d, finite_d).numpy()

    def newton_dir(pos_d, PF_d, active_d, finite_d):
        return ws.solve(pos_d, PF_d, reg_d, finite_d).numpy()

    report = {"nf": nf, "subdiv": subdiv, "n_total": int(cell.n_total),
              "residual_after_preload": float(res0), "regularizer_a": a_val, "gate": GATE}

    # (1) one ERM-Schwarz direction: descent + SPD/finite + line-search scan (+ region breakdown).
    mem = cell.membrane
    erm_m = np.unique(mem.erm_m_d.numpy().astype(np.int64))
    erm_c = np.unique(mem.erm_c_d.numpy().astype(np.int64))
    PF0 = project(cell, true_total_force(cell, _vec(pos0, d)))
    dx = erm_dir(_vec(pos0, d), _vec(PF0, d), wp.ones(1, dtype=wp.int32, device=d),
                 wp.ones(1, dtype=wp.int32, device=d))
    diag = erm.diagonal.numpy().reshape(-1, 6)   # per node: xx,yy,zz,xy,xz,yz
    trace = diag[:, 0] + diag[:, 1] + diag[:, 2]
    report["diag_trace"] = {
        "a_times3": 3.0 * a_val,
        "cortex_median": float(np.median(trace[erm_c])), "cortex_min": float(trace[erm_c].min()),
        "cortex_max": float(trace[erm_c].max()), "membrane_median": float(np.median(trace[erm_m])),
    }
    dxn = np.linalg.norm(dx, axis=1) * 1e3
    report["descent_dot_PF"] = float(np.einsum("ij,ij->", PF0, dx))   # > 0 ⇒ descent direction
    report["dx_finite"] = bool(np.all(np.isfinite(dx)))
    report["dx_max_nm"] = float(dxn.max())
    report["dx_max_membrane_nm"] = float(dxn[erm_m].max())
    report["dx_max_cortex_nm"] = float(dxn[erm_c].max())
    report["dx_argmax_is_membrane"] = bool(int(dxn.argmax()) in set(erm_m.tolist()))
    scan = []
    for t in (0.0, 1.0, 0.5, 0.25, 0.125, 0.0625):
        m = np.linalg.norm(project(cell, true_total_force(cell, _vec(pos0 + t * dx, d))), axis=1)
        scan.append({"t": t, "PF_max": float(m.max()),
                     "PF_max_membrane": float(m[erm_m].max()), "PF_max_cortex": float(m[erm_c].max())})
    report["erm_line_search_scan"] = scan

    # (2) iterated damped ERM-Schwarz vs damped-Newton reference (side by side).
    erm_traj, _ = _damped_descent(cell, pos0, erm_dir, reg_d, n_steps)
    newton_traj, _ = _damped_descent(cell, pos0, newton_dir, reg_d, min(n_steps, 12))
    report["erm_schwarz_descent"] = {
        "reached_gate": bool(erm_traj and erm_traj[-1].get("converged")),
        "final_res": erm_traj[-1]["res"] if erm_traj else None, "n_steps": len(erm_traj),
        "trajectory": erm_traj,
    }
    report["newton_descent"] = {
        "reached_gate": bool(newton_traj and newton_traj[-1].get("converged")),
        "final_res": newton_traj[-1]["res"] if newton_traj else None, "n_steps": len(newton_traj),
        "trajectory": newton_traj,
    }
    return report


def _fmt(r):
    L = [f"=== ERM-Schwarz probe  nf={r['nf']} subdiv={r['subdiv']}  n_total={r['n_total']} ===",
         f"residual after preload = {r['residual_after_preload']:.4f} pN   a={r['regularizer_a']:.4g}   gate={r['gate']}",
         "",
         f"(1) ONE ERM-Schwarz direction:  PF·dx = {r['descent_dot_PF']:.4g} (>0 ⇒ descent)   "
         f"finite={r['dx_finite']}",
         f"    diag trace: cortex median={r['diag_trace']['cortex_median']:.4g} min={r['diag_trace']['cortex_min']:.4g}"
         f"  membrane median={r['diag_trace']['membrane_median']:.4g}  (3a={r['diag_trace']['a_times3']:.4g}; "
         f"cortex should carry crosslink ~1.7e5)",
         f"    dx_max={r['dx_max_nm']:.4f} nm  (membrane {r['dx_max_membrane_nm']:.4f} / cortex "
         f"{r['dx_max_cortex_nm']:.4f}; argmax is {'MEMBRANE' if r['dx_argmax_is_membrane'] else 'cortex/other'})",
         "    line-search scan  t·dx → max|PF| (all / mem / cortex):"]
    for row in r["erm_line_search_scan"]:
        L.append(f"      t={row['t']:.4f}:  {row['PF_max']:.4f}   "
                 f"mem={row.get('PF_max_membrane', 0):.4f}  cor={row.get('PF_max_cortex', 0):.4f}")
    e = r["erm_schwarz_descent"]; n = r["newton_descent"]
    L += ["", f"(2) DAMPED DESCENT to gate {r['gate']}:",
          f"    ERM-Schwarz : reached_gate={e['reached_gate']}  final_res={e['final_res']:.4f}  steps={e['n_steps']}"]
    for row in e["trajectory"][:24]:
        tag = "  CONVERGED" if row.get("converged") else ("  STALLED" if row.get("stalled") else "")
        L.append(f"        step {row['step']:2d}: res={row['res']:.4f}  t={row['t']:.4f}{tag}")
    L.append(f"    Newton (ref): reached_gate={n['reached_gate']}  final_res={n['final_res']:.4f}  steps={n['n_steps']}  "
             f"(22d: stalls ~0.78 native / 0.92 8k)")
    return "\n".join(L)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--native", action="store_true")
    ap.add_argument("--subdiv", type=int, default=6)
    ap.add_argument("--steps", type=int, default=60)
    ap.add_argument("--json", type=str, default="")
    args = ap.parse_args()
    rep = run(70686 if args.native else 8000, args.subdiv, args.steps)
    print(_fmt(rep))
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rep, fh, indent=2, default=str)
        print(f"\n[json] {args.json}")

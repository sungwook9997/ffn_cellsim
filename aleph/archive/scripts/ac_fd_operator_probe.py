"""FD-vs-analytic-operator probe — P0 critical path for the resting-baseline blocker.

At the fixed subdiv-6 + overlap_free + radial-ERM + ERM-preload state (frozen state / active set),
answer the three PI questions BEFORE any new solver/preconditioner code (AC_DECISION_CARDS Card 4, P0):

  (a) Does the analytic operator action ``K.v`` match the true-force directional derivative
      ``-[F(x+eps v) - F(x-eps v)] / (2 eps)``?  Broken down PER FORCE FAMILY (bending / crosslink /
      ERM) and per node region, to name the FIRST mismatching family.
  (b) ``||P F||`` vs ``||F||`` — is the residual in the projection nullspace?
  (c) The proposed Newton ``dx = (aI + P K P)^-1 P F`` on membrane / loaded-cortex nodes vs the
      ~0.13 nm expected cortex drift and ~turgor/k_erm membrane drift.

This launches NO solver step that changes committed state (it perturbs a scratch copy and restores).
Warp-CUDA only (I0-A); runs on the gbook A5000, not the dev Mac.

Usage (gbook):
    PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
      ~/miniconda3/envs/ffn_sim/bin/python aleph/scripts/ac_fd_operator_probe.py [--native]
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import (
    _accumulate_all,
    _preload_erm_resting_balance,
    _residual_host,
)
from aleph.components.incumbent.erm_tether import erm_tether_force_kernel
from aleph.components.incumbent.implicit_mechanics import (
    ProjectedAnalyticCG,
    add_bending_stiffness_kernel,
    add_erm_stiffness_kernel,
    add_pair_stiffness_kernel,
    compute_regularization_kernel,
    omitted_regularization_base,
)
from aleph.components.incumbent.inner_mechanics import max_abs_pressure_kernel, project_constraint_forces_kernel
from aleph.laws.forces_warp import cytosim_bending_kernel
from aleph.laws.network_warp import link_spring_kernel
from aleph.engine.gate_criteria import (
    NATIVE_L500_DT_MU,
    NATIVE_L500_TOLERANCE_UM,
    predicted_force_floor,
)


# --------------------------------------------------------------------------------------------------
# Force / operator evaluators (each launches ONE kernel into a fresh (N,3) accumulator).
# --------------------------------------------------------------------------------------------------
def _zeros(n: int, device) -> wp.array:
    return wp.zeros(n, dtype=wp.vec3d, device=device)


def true_total_force(cell, pos_d: wp.array) -> np.ndarray:
    f = _zeros(cell.n_total, cell.device)
    _accumulate_all(cell, pos_d, f)
    wp.synchronize_device(cell.device)
    return f.numpy()


def true_family_force(cell, pos_d: wp.array, family: str) -> np.ndarray:
    """Exact runtime force of ONE family, launched in isolation on the given positions."""
    d = cell.device
    f = _zeros(cell.n_total, d)
    if family == "bending" and cell.n_tri:
        wp.launch(cytosim_bending_kernel, dim=cell.n_tri, inputs=[pos_d, cell.tri_d, cell.alpha_d, f], device=d)
    elif family == "crosslink" and cell.n_xl:
        wp.launch(link_spring_kernel, dim=cell.n_xl,
                  inputs=[pos_d, cell.xl_d, cell.kxl_d, cell.r0xl_d, f], device=d)
    elif family == "erm":
        mem = cell.membrane
        if mem is not None and mem.n_erm:
            wp.launch(erm_tether_force_kernel, dim=mem.n_erm,
                      inputs=[pos_d, mem.erm_m_d, mem.erm_c_d, mem.erm_bound_d, wp.float64(mem.k_erm),
                              mem.erm_rest_d, wp.float64(mem.f_rupt), f], device=d)
    wp.synchronize_device(d)
    return f.numpy()


def operator_family_action(cell, pos_d: wp.array, v_d: wp.array, family: str) -> np.ndarray:
    """Analytic operator tangent action K_family . v (no projection, no regularizer)."""
    d = cell.device
    out = _zeros(cell.n_total, d)
    if family == "bending" and cell.n_tri:
        wp.launch(add_bending_stiffness_kernel, dim=cell.n_tri, inputs=[v_d, cell.tri_d, cell.alpha_d, out], device=d)
    elif family == "crosslink" and cell.n_xl:
        wp.launch(add_pair_stiffness_kernel, dim=cell.n_xl,
                  inputs=[pos_d, v_d, cell.xl_d, cell.kxl_d, cell.r0xl_d, out], device=d)
    elif family == "erm":
        mem = cell.membrane
        if mem is not None and mem.n_erm:
            wp.launch(add_erm_stiffness_kernel, dim=mem.n_erm,
                      inputs=[pos_d, v_d, mem.erm_m_d, mem.erm_c_d, mem.erm_bound_d, mem.erm_rest_d,
                              wp.float64(mem.k_erm), wp.float64(mem.f_rupt), out], device=d)
    wp.synchronize_device(d)
    return out.numpy()


def fd_family(cell, pos0: np.ndarray, v: np.ndarray, family: str, eps: float) -> np.ndarray:
    """Central FD of the true family force: -[f(x+eps v) - f(x-eps v)] / (2 eps)  (= K_family . v)."""
    d = cell.device
    pos_d = wp.array(np.ascontiguousarray(pos0 + eps * v, np.float64), dtype=wp.vec3d, device=d)
    f_plus = true_family_force(cell, pos_d, family)
    pos_d.assign(np.ascontiguousarray(pos0 - eps * v, np.float64))
    f_minus = true_family_force(cell, pos_d, family)
    return -(f_plus - f_minus) / (2.0 * eps)


def fd_total(cell, pos0: np.ndarray, v: np.ndarray, eps: float) -> np.ndarray:
    d = cell.device
    pos_d = wp.array(np.ascontiguousarray(pos0 + eps * v, np.float64), dtype=wp.vec3d, device=d)
    f_plus = true_total_force(cell, pos_d)
    pos_d.assign(np.ascontiguousarray(pos0 - eps * v, np.float64))
    f_minus = true_total_force(cell, pos_d)
    return -(f_plus - f_minus) / (2.0 * eps)


def project(cell, vec: np.ndarray) -> np.ndarray:
    """Apply the exact NF2007 inextensibility projector P to a host (N,3) field."""
    d = cell.device
    src = wp.array(np.ascontiguousarray(vec, np.float64), dtype=wp.vec3d, device=d)
    dst = wp.zeros_like(src)
    wp.copy(dst, src)
    if cell.n_fibers:
        diag = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)
        rhs = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)
        finite = wp.ones(1, dtype=wp.int32, device=d)
        wp.launch(project_constraint_forces_kernel, dim=cell.n_fibers,
                  inputs=[cell.pos_d, src, cell.foff_d, cell.soff_d, dst, diag, rhs, finite], device=d)
    wp.synchronize_device(d)
    return dst.numpy()


# --------------------------------------------------------------------------------------------------
def _region_stats(name: str, vec: np.ndarray, idx: np.ndarray) -> dict:
    mag = np.linalg.norm(vec[idx], axis=1) if idx.size else np.zeros(0)
    return {
        "region": name, "n": int(idx.size),
        "max": float(mag.max()) if mag.size else 0.0,
        "l2": float(np.linalg.norm(mag)) if mag.size else 0.0,
        "mean": float(mag.mean()) if mag.size else 0.0,
    }


def _cmp(op: np.ndarray, fd: np.ndarray, idx: np.ndarray) -> dict:
    """Relative agreement of operator action vs FD directional derivative on a node set."""
    o, g = op[idx], fd[idx]
    num = float(np.linalg.norm(o - g))
    den = float(np.linalg.norm(g))
    dmag = np.linalg.norm(o - g, axis=1)
    worst = int(idx[int(dmag.argmax())]) if idx.size else -1
    return {
        "n": int(idx.size),
        "op_l2": float(np.linalg.norm(o)), "fd_l2": den,
        "abs_l2_diff": num, "rel_l2_diff": (num / den if den > 0 else float("nan")),
        "worst_node": worst, "worst_abs": float(dmag.max()) if idx.size else 0.0,
    }


def run(nf: int, subdiv: int, eps_list: list[float]) -> dict:
    wp.init()
    cfg = CellConfig(n_filaments=nf, overlap_free_cortex=True, erm_radial_pairing=True,
                     membrane_subdivisions=subdiv)
    cell = build_cell(cfg)
    d = cell.device
    report: dict = {"nf": nf, "subdiv": subdiv, "n_total": int(cell.n_total),
                    "n_actin": int(getattr(cell, "n_actin", -1)), "n_fibers": int(cell.n_fibers)}

    res_raw_before, _ = _residual_host(cell, cell.pos_d, cell.f_d)
    preload = _preload_erm_resting_balance(cell)
    report["preload"] = preload
    res_raw_after, _ = _residual_host(cell, cell.pos_d, cell.f_d)
    report["residual_raw_before_preload"] = float(res_raw_before)
    report["residual_raw_after_preload"] = float(res_raw_after)

    pos0 = cell.pos_d.numpy().copy()
    mem = cell.membrane
    erm_m = np.unique(mem.erm_m_d.numpy().astype(np.int64))
    erm_c = np.unique(mem.erm_c_d.numpy().astype(np.int64))
    n_actin = int(getattr(cell, "n_actin", cell.n_total))
    all_idx = np.arange(cell.n_total)
    other_actin = np.setdiff1d(np.arange(n_actin), erm_c)
    report["n_membrane_erm"] = int(erm_m.size)
    report["n_loaded_cortex"] = int(erm_c.size)

    # ---- (b) ||P F|| vs ||F|| ---------------------------------------------------------------------
    F = true_total_force(cell, cell.pos_d)
    PF = project(cell, F)
    report["F_regions"] = [_region_stats(n, F, ix) for n, ix in
                           (("all", all_idx), ("membrane", erm_m), ("loaded_cortex", erm_c),
                            ("other_actin", other_actin))]
    report["PF_regions"] = [_region_stats(n, PF, ix) for n, ix in
                            (("all", all_idx), ("membrane", erm_m), ("loaded_cortex", erm_c),
                             ("other_actin", other_actin))]
    report["projection_nullspace"] = {
        "F_l2": float(np.linalg.norm(F)), "PF_l2": float(np.linalg.norm(PF)),
        "ratio_PF_over_F": float(np.linalg.norm(PF) / (np.linalg.norm(F) + 1e-300)),
        "F_max": float(np.linalg.norm(F, axis=1).max()),
        "PF_max": float(np.linalg.norm(PF, axis=1).max()),
    }

    # Descent direction v = normalized P F (the actual RHS direction).
    v = PF / (np.linalg.norm(PF) + 1e-300)

    # ---- (a) FD vs operator, per family, in the descent direction v -------------------------------
    families = ["bending", "crosslink", "erm"]
    v_d = wp.array(np.ascontiguousarray(v, np.float64), dtype=wp.vec3d, device=d)
    fam_report = {}
    for fam in families:
        op_f = operator_family_action(cell, cell.pos_d, v_d, fam)
        per_eps = {}
        for eps in eps_list:
            fd_f = fd_family(cell, pos0, v, fam, eps)
            per_eps[f"{eps:.0e}"] = {
                "all": _cmp(op_f, fd_f, all_idx),
                "loaded_cortex": _cmp(op_f, fd_f, erm_c),
                "membrane": _cmp(op_f, fd_f, erm_m),
            }
        fam_report[fam] = per_eps
    report["family_fd_vs_operator"] = fam_report

    # Total FD vs total operator stiffness K.v (sum of ALL families incl. omitted turgor/wca/etc.),
    # plus Rayleigh stiffness in the residual direction: s = v . (K v). Compares "stiffness the true
    # force delivers" vs "stiffness the analytic K delivers" for the exact descent direction.
    eps_star = eps_list[len(eps_list) // 2]
    fd_tot = fd_total(cell, pos0, v, eps_star)
    # Operator total K.v via the workspace stiffness (regularizer 0, no projection).
    ws = ProjectedAnalyticCG(cell, max_iterations=max(256, 8 * subdiv), coarse_iterations=0, coarse_modes=0)
    zero_reg = wp.zeros(1, dtype=wp.float64, device=d)
    Kv_d = _zeros(cell.n_total, d)
    ws._stiffness(cell.pos_d, v_d, Kv_d, zero_reg)
    wp.synchronize_device(d)
    Kv = Kv_d.numpy()
    s_op = float(np.einsum("ij,ij->", v, Kv))
    s_fd = float(np.einsum("ij,ij->", v, fd_tot))
    report["rayleigh_stiffness_residual_dir"] = {
        "s_operator_vKv": s_op, "s_fd_true_vKv": s_fd,
        "ratio_op_over_fd": (s_op / s_fd if s_fd != 0 else float("nan")),
        "eps": eps_star,
        "Kv_total_all": _cmp(Kv, fd_tot, all_idx),
        "Kv_total_loaded_cortex": _cmp(Kv, fd_tot, erm_c),
        "Kv_total_membrane": _cmp(Kv, fd_tot, erm_m),
    }
    # Sum-of-isolated families vs total FD → how much stiffness lives in OMITTED families (turgor/wca/
    # membrane-area/nucleus/myosin), which K does NOT represent (they sit in the regularizer a).
    op_sum = sum(operator_family_action(cell, cell.pos_d, v_d, fam) for fam in families)
    report["omitted_family_stiffness"] = {
        "loaded_cortex": _cmp(op_sum, fd_tot, erm_c),
        "membrane": _cmp(op_sum, fd_tot, erm_m),
        "all": _cmp(op_sum, fd_tot, all_idx),
    }

    # ---- regularizer a (runtime value) ------------------------------------------------------------
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
    a_val = float(reg_d.numpy()[0])
    report["regularizer"] = {
        "a": a_val, "omitted_base": omitted_base,
        "max_pressure": float(max_pressure_d.numpy()[0]),
        "pressure_edge_um": float(cell.pressure_edge_um),
        "mechanical_kmax": float(cell.mechanical_kmax),
    }

    # ---- (c) Newton dx = (aI + PKP)^-1 P F --------------------------------------------------------
    PF_d = wp.array(np.ascontiguousarray(PF, np.float64), dtype=wp.vec3d, device=d)
    finite_d = wp.ones(1, dtype=wp.int32, device=d)
    dx_d = ws.solve(cell.pos_d, PF_d, reg_d, finite_d)
    wp.synchronize_device(d)
    dx = dx_d.numpy()
    dx_um_to_nm = 1.0e3
    report["newton_dx_nm"] = {
        "membrane": _region_stats("membrane", dx * dx_um_to_nm, erm_m),
        "loaded_cortex": _region_stats("loaded_cortex", dx * dx_um_to_nm, erm_c),
        "all": _region_stats("all", dx * dx_um_to_nm, all_idx),
        "expected_cortex_nm": (float(res_raw_after) / 17653.0) * dx_um_to_nm,
        "expected_membrane_nm": float(np.linalg.norm(F[erm_m], axis=1).mean() / float(mem.k_erm)) * dx_um_to_nm,
    }

    # Actual residual reduction if this Newton dx were applied (perturb scratch, measure, restore).
    pos_try = wp.array(np.ascontiguousarray(pos0 + dx, np.float64), dtype=wp.vec3d, device=d)
    F_try = true_total_force(cell, pos_try)
    PF_try = project(cell, F_try)
    report["newton_step_effect"] = {
        "PF_max_before": float(np.linalg.norm(PF, axis=1).max()),
        "PF_max_after": float(np.linalg.norm(PF_try, axis=1).max()),
        "reduction_frac": 1.0 - float(np.linalg.norm(PF_try, axis=1).max()) /
        (float(np.linalg.norm(PF, axis=1).max()) + 1e-300),
    }

    # Line-search scan: does ANY fraction of dx descend? Distinguishes a bad DIRECTION (ascent ⇒ the SPD
    # operator flipped an indefinite/negative-curvature mode) from a good direction that merely OVERSHOOTS
    # (step-size ⇒ a globalization/scaling fix suffices). Reports region-wise so we see who overshoots.
    pf0_all = float(np.linalg.norm(PF, axis=1).max())
    pf0_mem = float(np.linalg.norm(PF[erm_m], axis=1).max())
    pf0_cor = float(np.linalg.norm(PF[erm_c], axis=1).max())
    tscan = [{"t": 0.0, "PF_max": pf0_all, "PF_max_membrane": pf0_mem, "PF_max_cortex": pf0_cor}]
    for t in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625):
        pos_t = wp.array(np.ascontiguousarray(pos0 + t * dx, np.float64), dtype=wp.vec3d, device=d)
        m = np.linalg.norm(project(cell, true_total_force(cell, pos_t)), axis=1)
        tscan.append({"t": t, "PF_max": float(m.max()),
                      "PF_max_membrane": float(m[erm_m].max()),
                      "PF_max_cortex": float(m[erm_c].max())})
    report["line_search_scan"] = tscan

    # Compression stats: negative-curvature the SPD operator clamps away (transverse_k=max(k·ext/L,0)).
    # A large clamped-away stiffness at the preloaded state is the candidate cause of an indefinite true
    # Hessian that the PSD operator cannot descend.
    xl = cell.xl_d.numpy().astype(np.int64)
    kxl = cell.kxl_d.numpy()
    r0xl = cell.r0xl_d.numpy()
    dvec = pos0[xl[:, 1]] - pos0[xl[:, 0]]
    Lxl = np.linalg.norm(dvec, axis=1)
    ext = Lxl - r0xl
    comp = ext < 0.0
    neg_trans = np.where(comp & (Lxl > 1e-30), kxl * ext / np.maximum(Lxl, 1e-30), 0.0)  # ≤ 0
    report["crosslink_compression"] = {
        "n_xl": int(xl.shape[0]), "frac_compressed": float(comp.mean()),
        "ext_min_nm": float(ext.min() * 1e3), "ext_mean_nm": float(ext.mean() * 1e3),
        "ext_max_nm": float(ext.max() * 1e3),
        "clamped_neg_transverse_k_sum": float(-neg_trans.sum()),
        "clamped_neg_transverse_k_max": float(-neg_trans.min()),
        "axial_k_typical": float(np.median(kxl)),
    }

    # DECISIVE convergence reference: iterate a HOST-side damped Newton using ONLY the existing operator
    # (ws.solve) + a backtracking line search over t in {1,1/2,...}. This is a diagnostic reference (like the
    # NumPy oracles), NOT a production solver: it answers whether the correct operator + proper globalization
    # reaches the strict gate, which decides conditioned-implicit (fork 2) vs targeted preload (fork 1).
    # D4 (PI 2026-07-28): derived from the run's own predicate, not the 0.21 literal.
    # F_pred scales with mesh spacing — do not score a finer rung against this value.
    gate = predicted_force_floor(
        inner_tolerance_um=NATIVE_L500_TOLERANCE_UM, inner_dt_mu=NATIVE_L500_DT_MU
    )
    pos_it = pos0.copy()
    traj = []
    n_newton_steps = 20 if nf >= 70000 else 40
    for step in range(n_newton_steps):
        Fi = true_total_force(cell, wp.array(np.ascontiguousarray(pos_it, np.float64), dtype=wp.vec3d, device=d))
        PFi = project(cell, Fi)
        res = float(np.linalg.norm(PFi, axis=1).max())
        if res <= gate:
            traj.append({"step": step, "res": res, "t": 0.0, "converged": True})
            break
        PFi_d = wp.array(np.ascontiguousarray(PFi, np.float64), dtype=wp.vec3d, device=d)
        pos_it_d = wp.array(np.ascontiguousarray(pos_it, np.float64), dtype=wp.vec3d, device=d)
        dxi = ws.solve(pos_it_d, PFi_d, reg_d, finite_d).numpy()
        best_t, best_res = 0.0, res
        for t in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125):
            m = np.linalg.norm(project(cell, true_total_force(
                cell, wp.array(np.ascontiguousarray(pos_it + t * dxi, np.float64), dtype=wp.vec3d, device=d))), axis=1).max()
            if m < best_res:
                best_res, best_t = float(m), t
        if best_t == 0.0:            # no descending fraction found → stall
            traj.append({"step": step, "res": res, "t": 0.0, "converged": False, "stalled": True})
            break
        pos_it = pos_it + best_t * dxi
        traj.append({"step": step, "res": res, "t": best_t, "converged": False})
    # Characterize the STALLED (or final) state: which region holds the residual max, and did crosslink
    # compression develop under the accumulated motion (SPD clamp candidate)?
    F_fin = true_total_force(cell, wp.array(np.ascontiguousarray(pos_it, np.float64), dtype=wp.vec3d, device=d))
    PF_fin = project(cell, F_fin)
    dfin = pos_it[xl[:, 1]] - pos_it[xl[:, 0]]
    Lfin = np.linalg.norm(dfin, axis=1)
    extfin = Lfin - r0xl
    report["damped_newton_reference"] = {
        "gate": gate, "n_steps": len(traj),
        "final_res": traj[-1]["res"] if traj else None,
        "reached_gate": bool(traj and traj[-1].get("converged")),
        "trajectory": traj,
        "final_state": {
            "PF_max_all": float(np.linalg.norm(PF_fin, axis=1).max()),
            "PF_max_membrane": float(np.linalg.norm(PF_fin[erm_m], axis=1).max()),
            "PF_max_loaded_cortex": float(np.linalg.norm(PF_fin[erm_c], axis=1).max()),
            "PF_max_other_actin": float(np.linalg.norm(PF_fin[other_actin], axis=1).max())
            if other_actin.size else 0.0,
            "xl_frac_compressed": float((extfin < 0.0).mean()),
            "xl_ext_min_nm": float(extfin.min() * 1e3), "xl_ext_max_nm": float(extfin.max() * 1e3),
            "total_drift_nm": float(np.linalg.norm(pos_it - pos0, axis=1).max() * 1e3),
        },
    }
    return report


def _fmt(report: dict) -> str:
    r = report
    lines = [
        f"=== FD-vs-operator probe  nf={r['nf']} subdiv={r['subdiv']}  n_total={r['n_total']} "
        f"n_fibers={r['n_fibers']}  membrane_erm={r['n_membrane_erm']} loaded_cortex={r['n_loaded_cortex']} ===",
        f"residual raw: before_preload={r['residual_raw_before_preload']:.4f}  "
        f"after_preload={r['residual_raw_after_preload']:.4f} pN",
        "",
        "(b) PROJECTION NULLSPACE  ||PF|| vs ||F||:",
    ]
    pn = r["projection_nullspace"]
    lines.append(f"    ||F||={pn['F_l2']:.4g}  ||PF||={pn['PF_l2']:.4g}  ratio={pn['ratio_PF_over_F']:.4f}"
                 f"   maxF={pn['F_max']:.4f}  maxPF={pn['PF_max']:.4f}")
    for row in r["PF_regions"]:
        fr = next(x for x in r["F_regions"] if x["region"] == row["region"])
        lines.append(f"      {row['region']:14s} n={row['n']:6d}  |F|max={fr['max']:.4f} l2={fr['l2']:.4g}"
                     f"  |PF|max={row['max']:.4f} l2={row['l2']:.4g}")
    lines += ["", "(a) FD vs OPERATOR per family (rel_l2_diff should be ~0 if tangent is exact):"]
    for fam, per_eps in r["family_fd_vs_operator"].items():
        lines.append(f"    [{fam}]")
        for eps, sets in per_eps.items():
            c = sets["loaded_cortex"]; m = sets["membrane"]; a = sets["all"]
            lines.append(f"      eps={eps}: all rel={a['rel_l2_diff']:.2e}  "
                         f"loaded_cortex rel={c['rel_l2_diff']:.2e} (fd_l2={c['fd_l2']:.3g})  "
                         f"membrane rel={m['rel_l2_diff']:.2e} (fd_l2={m['fd_l2']:.3g})")
    rq = r["rayleigh_stiffness_residual_dir"]
    lines += ["", "    STIFFNESS in residual direction  s = v.(Kv):",
              f"      s_operator={rq['s_operator_vKv']:.5g}   s_fd_true={rq['s_fd_true_vKv']:.5g}   "
              f"ratio_op/fd={rq['ratio_op_over_fd']:.4f}   (eps={rq['eps']:.0e})",
              f"      total Kv vs FD: all rel={rq['Kv_total_all']['rel_l2_diff']:.2e}  "
              f"loaded_cortex rel={rq['Kv_total_loaded_cortex']['rel_l2_diff']:.2e}  "
              f"membrane rel={rq['Kv_total_membrane']['rel_l2_diff']:.2e}"]
    of = r["omitted_family_stiffness"]
    lines.append(f"      (sum bending+xl+erm) vs total-FD:  loaded_cortex rel={of['loaded_cortex']['rel_l2_diff']:.2e}"
                 f"  membrane rel={of['membrane']['rel_l2_diff']:.2e}  -> nonzero ⇒ omitted turgor/wca/etc carry tangent")
    reg = r["regularizer"]
    lines += ["", f"    regularizer a={reg['a']:.5g}  (omitted_base={reg['omitted_base']:.5g}, "
              f"max_pressure*edge={reg['max_pressure']*reg['pressure_edge_um']:.5g}, kmax={reg['mechanical_kmax']:.4g})"]
    dx = r["newton_dx_nm"]
    lines += ["", "(c) NEWTON dx = (aI+PKP)^-1 PF  [nm]:",
              f"    membrane      max={dx['membrane']['max']:.4f} mean={dx['membrane']['mean']:.4f}  "
              f"(expected≈{dx['expected_membrane_nm']:.4f})",
              f"    loaded_cortex max={dx['loaded_cortex']['max']:.4f} mean={dx['loaded_cortex']['mean']:.4f}  "
              f"(expected≈{dx['expected_cortex_nm']:.4f})",
              f"    all           max={dx['all']['max']:.4f}"]
    ne = r["newton_step_effect"]
    lines.append(f"    ONE full Newton step: PF_max {ne['PF_max_before']:.4f} -> {ne['PF_max_after']:.4f} "
                 f"({100*ne['reduction_frac']:.1f}% reduction)")
    lines += ["", "    LINE-SEARCH SCAN  t*dx -> PF_max (all / membrane / cortex):"]
    for row in r["line_search_scan"]:
        lines.append(f"      t={row['t']:.5f}:  all={row['PF_max']:.4f}   "
                     f"mem={row['PF_max_membrane']:.4f}   cor={row['PF_max_cortex']:.4f}")
    cc = r["crosslink_compression"]
    lines += ["", f"    CROSSLINK COMPRESSION: n_xl={cc['n_xl']}  frac_compressed={cc['frac_compressed']:.3f}"
              f"  ext[min/mean/max]={cc['ext_min_nm']:.3f}/{cc['ext_mean_nm']:.3f}/{cc['ext_max_nm']:.3f} nm",
              f"      clamped-away neg transverse k: sum={cc['clamped_neg_transverse_k_sum']:.4g} "
              f"max={cc['clamped_neg_transverse_k_max']:.4g}  (axial k typ={cc['axial_k_typical']:.4g} pN/µm)"]
    dn = r["damped_newton_reference"]
    lines += ["", f"    DAMPED-NEWTON REFERENCE (existing operator + backtracking LS, gate={dn['gate']}):",
              f"      reached_gate={dn['reached_gate']}  n_steps={dn['n_steps']}  final_res={dn['final_res']:.4f}"]
    for row in dn["trajectory"][:24]:
        tag = "  CONVERGED" if row.get("converged") else ("  STALLED" if row.get("stalled") else "")
        lines.append(f"        step {row['step']:2d}: res={row['res']:.4f}  best_t={row['t']:.4f}{tag}")
    fs = dn.get("final_state", {})
    if fs:
        lines += [f"      STALLED-STATE regions: all={fs['PF_max_all']:.4f}  mem={fs['PF_max_membrane']:.4f}"
                  f"  cortex={fs['PF_max_loaded_cortex']:.4f}  other_actin={fs['PF_max_other_actin']:.4f}",
                  f"      xl @ stall: frac_compressed={fs['xl_frac_compressed']:.3f}  "
                  f"ext[min/max]={fs['xl_ext_min_nm']:.3f}/{fs['xl_ext_max_nm']:.3f} nm  "
                  f"max_drift={fs['total_drift_nm']:.3f} nm"]
    return "\n".join(lines)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--native", action="store_true", help="70,686 cortical filaments (else 8000)")
    ap.add_argument("--subdiv", type=int, default=6)
    ap.add_argument("--json", type=str, default="")
    args = ap.parse_args()
    nf = 70686 if args.native else 8000
    rep = run(nf, args.subdiv, eps_list=[1e-7, 1e-8, 1e-9])
    print(_fmt(rep))
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rep, fh, indent=2, default=str)
        print(f"\n[json] {args.json}")

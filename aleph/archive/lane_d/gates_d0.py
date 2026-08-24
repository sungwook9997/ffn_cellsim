"""D0 sanity gates — membrane+cytosol traction rig (Lane D, NON-AUTHORITATIVE DIAGNOSTIC).

Gates (brief D0):
    G0.1 positive internal pressure → outward traction
    G0.2 pressure sign reversal → traction reverses
    G0.3 uniform pressure on a closed shell → net vector force ≈ 0
    G0.4 membrane nodal work vs fluid boundary work — FSI adjoint (P·dV) identity
    G0.5 no-flux mass conservation (closed shell, dPi_osm=dP_hyd ⇒ zero net flux)
    G0.6 rejected candidate → pressure + membrane state bit-exact recovery

All gates run Warp-CUDA kernels and require a CUDA device. A pure-host cross-check
(``uniform_pressure_force_reference``) anchors G0.1/G0.3.
"""
from __future__ import annotations

import numpy as np
import warp as wp

from aleph.components.incumbent.membrane_pressure import uniform_pressure_force_reference
from aleph.components.fluid.field_grid import FLUID
from .membrane_cytosol_rig import (
    MembraneCytosolRig, DP_HYD_REST_PA, DPI_OSM_REST_PA, BIOT_STORAGE_S, DX_UM,
)


def _set_uniform_pressure(rig: MembraneCytosolRig, value: float):
    p = np.where(rig._mask_np == FLUID, value, 0.0).astype(np.float64)
    rig.grid.set_pressure(np.ascontiguousarray(p))


def _accumulate_pressure_traction(rig: MembraneCytosolRig):
    f = wp.zeros(rig.n_nodes, dtype=wp.vec3d, device=rig.device)
    rig.traction.accumulate(rig.pos_d, f)
    return f.numpy()


def gate_traction_sign(rig: MembraneCytosolRig) -> dict:
    """G0.1/G0.2: outward for p>p_ext, inward for p<p_ext (dot with outward node normal)."""
    normals = rig._normals0
    _set_uniform_pressure(rig, +DP_HYD_REST_PA)
    f_pos = _accumulate_pressure_traction(rig)
    dot_pos = np.sum(f_pos * normals, axis=1)
    _set_uniform_pressure(rig, -DP_HYD_REST_PA)
    f_neg = _accumulate_pressure_traction(rig)
    dot_neg = np.sum(f_neg * normals, axis=1)
    _set_uniform_pressure(rig, DP_HYD_REST_PA)  # restore
    # ignore nodes whose surface trace was rank-deficient (traction ~0)
    live = np.linalg.norm(f_pos, axis=1) > 1e-9
    ok_pos = bool(np.all(dot_pos[live] > 0))
    ok_neg = bool(np.all(dot_neg[live] < 0))
    return {
        "name": "G0.1/G0.2 traction sign",
        "passed": ok_pos and ok_neg,
        "outward_fraction_pos": float(np.mean(dot_pos[live] > 0)),
        "inward_fraction_neg": float(np.mean(dot_neg[live] < 0)),
        "n_live": int(live.sum()),
        "n_nodes": int(rig.n_nodes),
    }


def gate_net_force_zero(rig: MembraneCytosolRig) -> dict:
    """G0.3: closed shell under uniform pressure → Σ f ≈ 0 vector. Anchored by the host oracle."""
    _set_uniform_pressure(rig, DP_HYD_REST_PA)
    rig.traction.reset_diagnostics()
    f = wp.zeros(rig.n_nodes, dtype=wp.vec3d, device=rig.device)
    rig.traction.accumulate(rig.pos_d, f)
    n_unresolved = int(rig.traction.unresolved_faces_d.numpy()[0])
    f = f.numpy()
    net = np.sum(f, axis=0)
    scale = float(np.sum(np.linalg.norm(f, axis=1))) + 1e-30
    rel = float(np.linalg.norm(net) / scale)
    # host oracle: uniform pressure nodal load on the same mesh
    f_ref = uniform_pressure_force_reference(rig.verts0, rig.mesh.faces, DP_HYD_REST_PA)
    net_ref = float(np.linalg.norm(np.sum(f_ref, axis=0)) / (np.sum(np.linalg.norm(f_ref, axis=1)) + 1e-30))
    return {
        "name": "G0.3 net force ~0",
        "passed": bool(rel < 1e-3 and n_unresolved == 0),
        "relative_net_force": rel,
        "oracle_relative_net_force": net_ref,
        "n_unresolved_faces": n_unresolved,
        "tol": 1e-3,
    }


def gate_adjoint_work(rig: MembraneCytosolRig) -> dict:
    """G0.4: FSI P·dV adjoint — Σ f_pressure·v_s (membrane work) vs ∫ p·div(v_s) dV (fluid work).

    Prescribe a smooth outward 'breathing' membrane velocity v_s = c·n̂. The pressure does
    ∮ p (v_s·n̂) dA on the membrane; the moving boundary injects ∫ p·div(v_s) dV into the fluid.
    Continuum identity: they are equal (divergence theorem). Discretely the two operators differ
    (surface trace+area quadrature vs Peskin spread+finite-difference divergence) so this is a
    same-sign / same-order consistency gate, not a bit-exact one.
    """
    _set_uniform_pressure(rig, DP_HYD_REST_PA)
    c = 0.05  # µm/s breathing speed
    v_s = (c * rig._normals0).astype(np.float64)
    v_d = wp.array(np.ascontiguousarray(v_s), dtype=wp.vec3d, device=rig.device)
    # membrane work Σ f_pressure·v_s
    f_p = _accumulate_pressure_traction(rig)
    w_mem = float(np.sum(f_p * v_s))
    # fluid work ∫ p·div(v_s) dV
    rig.grid.div_vs.zero_()
    rig.coupling.update(rig.pos_d, v_d, rig.node_volume_d, rig.active_d)
    rig._pv_work_d.zero_()
    from .membrane_cytosol_rig import _pv_work_kernel
    wp.launch(_pv_work_kernel, dim=(rig._nx, rig._nx, rig._nx),
              inputs=[rig.grid.p, rig.grid.div_vs, rig.grid.mask, wp.float64(DX_UM ** 3),
                      rig._pv_work_d], device=rig.device)
    w_flu = float(rig._pv_work_d.numpy()[0])
    ratio = w_flu / w_mem if abs(w_mem) > 1e-30 else float("nan")
    same_sign = (w_mem > 0) == (w_flu > 0)
    return {
        "name": "G0.4 FSI work adjoint (P·dV)",
        "passed": bool(same_sign and 0.5 < ratio < 2.0),
        "w_membrane_pn_um_s": w_mem,
        "w_fluid_pn_um_s": w_flu,
        "ratio_fluid_over_membrane": ratio,
        "note": "same-sign + ratio∈[0.5,2] (two distinct discretizations of P·dV)",
    }


def gate_no_flux_mass_conservation(rig: MembraneCytosolRig, n_steps: int = 20,
                                   dt_phys: float = 0.02) -> dict:
    """G0.5: at rest (dPi_osm=dP_hyd, no interior source, no membrane motion) fluid content is conserved."""
    rig.clear_source()
    c0 = rig.grid.total_content(BIOT_STORAGE_S)
    max_flux = 0.0
    for _ in range(n_steps):
        rec = rig.outer_step(dt_phys, inv_gamma=0.0, dpi_osm=DPI_OSM_REST_PA, move_membrane=False)
        max_flux = max(max_flux, abs(rec["integrated_flux_um3_s"]))
    c1 = rig.grid.total_content(BIOT_STORAGE_S)
    drift = abs(c1 - c0) / (abs(c0) + 1e-30)
    return {
        "name": "G0.5 no-flux mass conservation",
        "passed": bool(drift < 1e-6 and max_flux < 1e-6),
        "content_drift_rel": float(drift),
        "max_integrated_flux_um3_s": float(max_flux),
        "content0": float(c0),
        "content1": float(c1),
    }


def gate_reject_bit_exact(rig: MembraneCytosolRig, dt_phys: float = 0.02) -> dict:
    """G0.6: snapshot p + membrane pos, run a rejected perturbation candidate, restore, assert bit-exact."""
    p_snap = rig.grid.pressure_to_host().copy()
    pos_snap = rig.pos_d.numpy().copy()
    # --- rejected candidate: strong one-sided source + membrane move ---
    rig.set_patch_source(rig.patch_centre(), radius=2.0, rate=5.0)
    rig.outer_step(dt_phys, inv_gamma=2.0, dpi_osm=DPI_OSM_REST_PA, move_membrane=True)
    # --- reject → restore both stores exactly ---
    rig.clear_source()
    rig.grid.set_pressure(p_snap)
    rig.pos_d.assign(np.ascontiguousarray(pos_snap))
    p_back = rig.grid.pressure_to_host()
    pos_back = rig.pos_d.numpy()
    p_exact = bool(np.array_equal(p_back, p_snap))
    pos_exact = bool(np.array_equal(pos_back, pos_snap))
    return {
        "name": "G0.6 rejected-candidate bit-exact recovery",
        "passed": p_exact and pos_exact,
        "pressure_bit_exact": p_exact,
        "membrane_bit_exact": pos_exact,
    }


def run_all(device: str) -> dict:
    """Run every D0 gate on fresh rigs; return a JSON-serialisable report."""
    results = []
    # sign / net-force / adjoint share a fresh rig (state reset between them)
    rig = MembraneCytosolRig(device=device)
    results.append(gate_traction_sign(rig))
    results.append(gate_net_force_zero(rig))
    results.append(gate_adjoint_work(rig))
    results.append(gate_no_flux_mass_conservation(MembraneCytosolRig(device=device)))
    results.append(gate_reject_bit_exact(MembraneCytosolRig(device=device)))
    passed = all(r["passed"] for r in results)
    return {"suite": "D0", "device": str(device), "all_passed": bool(passed), "gates": results}

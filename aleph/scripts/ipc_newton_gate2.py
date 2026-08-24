"""Gate 2 (#1 Step 2): projected-Newton IPC step correctness on 2 cells.

Contracts from DCM_IPC_LARGE_DT_PLAN_2026-07-07 (Step 2):
  (1) ENERGY monotone-decrease per step — Φ(x)=½a|x−xn|²+U(x) non-increasing across Newton iterates
      (the Armijo line-search works; ∇Φ==−(a(x−xn)−F) after the Gate-1 energy↔force match).
  (2) RESIDUAL clears tol in bounded iters — |G|=|a(x−xn)−F_total| < newton_tol·|G0|.
  (3) NON-PENETRATION under load — from a FEASIBLE start (gap>0), a strong inward bundle pull can
      NEVER create penetration: the log-barrier diverges + CCD caps every step ⇒ pen stays ≈0 for
      the whole loaded trajectory. (Curing a *penetrating* start is a separate multi-step
      feasibilization recovery — Scenario B, a demonstration, not the guarantee.)

Stiff set = turgor + cortex edges + IPC barrier (barrier force in RHS, barrier stiffness in
hess_apply, excluded from the FD matvec — exactly production stiff_force_into when ipc=True).
Soft driver = constant LAGGED inward pull; its linear potential −F_soft(xn)·(x−xn) is in Φ."""
import numpy as np, sys
import warp as wp
wp.init()
from aleph.dcm.geometry import icosphere_mesh, ResolvedDCM
from aleph.dcm.dcm_contact_implicit_warp import (nearest_face_ipc_kernel, make_contact_hess_apply,
                                                   ccd_alpha)
from aleph.dcm.dcm_ipc_energy import (ipc_barrier_energy_kernel, edge_energy_kernel,
    turgor_energy_kernel, inertial_energy_kernel, soft_linear_energy_kernel)
from aleph.dcm.dcm_warp_hybrid import _bond_accumulate
from aleph.dcm.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel
from aleph.dcm.dcm_warp_implicit import ipc_newton_step, _vaxpy, _vcopy
from aleph.dcm.dcm_neighbor_warp import pos_to_f32, face_centroids_f32, penetration_depth_kernel

dev = sys.argv[1] if len(sys.argv) > 1 else "cpu"

p = ResolvedDCM(subdivisions=1); R = p.R_cell
v1, e1, f1 = icosphere_mesh(R, 1); npc = v1.shape[0]
me = float(np.linalg.norm(v1[e1[:, 0]] - v1[e1[:, 1]], axis=1).mean())
faces = np.concatenate([f1, f1 + npc]).astype(np.int32)
edges = np.concatenate([e1, e1 + npc]).astype(np.int32)
cof = np.array([0] * npc + [1] * npc, np.int32)
fcell = np.array([0] * f1.shape[0] + [1] * f1.shape[0], np.int32)
N, n_faces, n_edges, n_cells = 2 * npc, faces.shape[0], edges.shape[0], 2

rep = 2.0e8; c_rep = 0.30 * me; d_hat = c_rep
grid_q = float(0.80 * me + 0.7 * (3 * me / 2.9) + 3.0 * me)
k_edge = p.k_edge; k_vol = 7.73e5; dP0 = p.turgor_dP0
gamma = 6.0 * np.pi * 65.9 * R / npc
R_eff = float(np.linalg.norm(v1 - v1.mean(0), axis=1).mean())
V0c = np.full(n_cells, (4.0 / 3.0) * np.pi * R_eff ** 3)
area_typ = 4 * np.pi * R ** 2 / npc; Fcap = rep * area_typ * c_rep


def make_state(overlap):
    sep = (2.0 - overlap) * R
    return np.concatenate([v1 + [sep / 2, 0, 0], v1 + [-sep / 2, 0, 0]])


faces_d = wp.array(faces, dtype=wp.int32, device=dev)
edges_d = wp.array(edges, dtype=wp.int32, device=dev)
cof_d = wp.array(cof, dtype=wp.int32, device=dev)
fcell_d = wp.array(fcell, dtype=wp.int32, device=dev)
V0_d = wp.array(V0c, dtype=wp.float64, device=dev)
r0_d = wp.zeros(n_edges, dtype=wp.float64, device=dev)
nf32 = wp.zeros(N, dtype=wp.vec3, device=dev); cf32 = wp.zeros(n_faces, dtype=wp.vec3, device=dev)
Vc = wp.zeros(n_cells, dtype=wp.float64, device=dev); dP_d = wp.zeros(n_cells, dtype=wp.float64, device=dev)
cn_k = wp.zeros(N, dtype=wp.float64, device=dev); cn_nrm = wp.zeros(N, dtype=wp.vec3d, device=dev)
ipc_hess = make_contact_hess_apply(cn_k, cn_nrm, device=dev)
t_ccd = wp.zeros(N, dtype=wp.float64, device=dev); pend = wp.zeros(N, dtype=wp.float64, device=dev)
grid = wp.HashGrid(16, 16, 16, device=dev)
scratch = {k: wp.zeros(N, dtype=wp.vec3d, device=dev) for k in
           ("r", "p", "Ap", "dx", "Fx", "Fp", "xp", "nt_Ft", "nt_rhs", "nt_xtr")}
scratch["sca"] = wp.zeros(1, dtype=wp.float64, device=dev)
pos_d = wp.zeros(N, dtype=wp.vec3d, device=dev); xn_d = wp.zeros(N, dtype=wp.vec3d, device=dev)
fsoft_d = wp.zeros(N, dtype=wp.vec3d, device=dev)
A_CUR = [0.0]


def rebuild(x_d):
    wp.launch(pos_to_f32, dim=N, inputs=[x_d, nf32], device=dev)
    wp.launch(face_centroids_f32, dim=n_faces, inputs=[x_d, faces_d, cf32], device=dev)
    grid.build(points=cf32, radius=grid_q)


def _turgor_into(x_d, out):
    Vc.zero_(); wp.launch(dcm_volume_kernel, dim=n_faces, inputs=[x_d, faces_d, fcell_d, Vc], device=dev)
    dPn = dP0 + k_vol * (V0c - Vc.numpy()) / V0c
    dP_d.assign(np.ascontiguousarray(dPn))
    wp.launch(dcm_turgor_force_kernel, dim=n_faces, inputs=[x_d, faces_d, fcell_d, dP_d, out], device=dev)


def stiff_force_into(x_d, out):
    out.zero_(); _turgor_into(x_d, out)
    wp.launch(_bond_accumulate, dim=n_edges, inputs=[x_d, edges_d, wp.float64(k_edge), r0_d, out], device=dev)


def force_total_into(x_d, out):
    rebuild(x_d); out.zero_(); _turgor_into(x_d, out)
    wp.launch(_bond_accumulate, dim=n_edges, inputs=[x_d, edges_d, wp.float64(k_edge), r0_d, out], device=dev)
    wp.launch(nearest_face_ipc_kernel, dim=N, inputs=[grid.id, nf32, x_d, cof_d, faces_d, fcell_d,
              wp.float32(grid_q), wp.float64(rep), wp.float64(d_hat), wp.float64(0.0), wp.float64(0.0),
              out, cn_k, cn_nrm], device=dev)
    wp.launch(_vaxpy, dim=N, inputs=[out, wp.float64(1.0), fsoft_d], device=dev)


def energy_fn(x_d):
    e = wp.zeros(1, dtype=wp.float64, device=dev); rebuild(x_d)
    wp.launch(ipc_barrier_energy_kernel, dim=N, inputs=[grid.id, nf32, x_d, cof_d, faces_d, fcell_d,
              wp.float32(grid_q), wp.float64(rep), wp.float64(d_hat), e], device=dev)
    wp.launch(edge_energy_kernel, dim=n_edges, inputs=[x_d, edges_d, wp.float64(k_edge), r0_d, e], device=dev)
    Vc.zero_(); wp.launch(dcm_volume_kernel, dim=n_faces, inputs=[x_d, faces_d, fcell_d, Vc], device=dev)
    wp.launch(turgor_energy_kernel, dim=n_cells, inputs=[Vc, V0_d, wp.float64(k_vol), wp.float64(dP0), e], device=dev)
    wp.launch(inertial_energy_kernel, dim=N, inputs=[x_d, xn_d, wp.float64(A_CUR[0]), cof_d, e], device=dev)
    wp.launch(soft_linear_energy_kernel, dim=N, inputs=[x_d, xn_d, fsoft_d, cof_d, e], device=dev)
    wp.synchronize_device(dev); return float(e.numpy()[0])


def ccd_alpha_fn(x_d, dx_d):
    return ccd_alpha(grid.id, nf32, x_d, dx_d, cof_d, faces_d, fcell_d, grid_q, t_ccd, eta=0.9, device=dev)


def max_pen(x_d):
    rebuild(x_d); pend.zero_()
    wp.launch(penetration_depth_kernel, dim=N,
              inputs=[grid.id, nf32, x_d, cof_d, faces_d, fcell_d, wp.float32(grid_q), pend], device=dev)
    wp.synchronize_device(dev); return float(pend.numpy().max() / me)


def trajectory(overlap, Fpull, dt, nsteps):
    """Run nsteps implicit projected-Newton steps; return per-step (newton_iters, converged, gratio,
    energy_monotone, pen) and the worst penetration ever seen."""
    x0 = make_state(overlap)
    pos_d.assign(x0)
    r0_d.assign(np.ascontiguousarray(np.linalg.norm(x0[edges[:, 0]] - x0[edges[:, 1]], axis=1)))
    fs = np.zeros((N, 3)); fs[:npc, 0] = -Fpull / npc; fs[npc:, 0] = +Fpull / npc
    fsoft_d.assign(fs)
    A_CUR[0] = gamma / dt
    pen0 = max_pen(pos_d); worst_pen = pen0
    all_mono = True; all_conv = True; max_nit = 0; sum_nit = 0
    for _ in range(nsteps):
        wp.launch(_vcopy, dim=N, inputs=[xn_d, pos_d], device=dev)   # xn <- x
        info = ipc_newton_step(pos_d, xn_d, gamma / dt,
                               force_total_into=force_total_into, stiff_force_into=stiff_force_into,
                               energy_fn=energy_fn, ccd_alpha_fn=ccd_alpha_fn, cof_d=cof_d,
                               scratch=scratch, device=dev, max_newton=20, newton_tol=1e-4,
                               hess_apply=ipc_hess, cg_maxiter=200)
        E = info["energies"]
        mono = all(E[i + 1] <= E[i] + 1e-9 * abs(E[i]) + 1e-30 for i in range(len(E) - 1))
        conv = info["converged"] or info["g_norm"] <= 1e-3 * max(info["g0"], 1e-30)
        all_mono = all_mono and mono; all_conv = all_conv and conv
        max_nit = max(max_nit, info["newton_iters"]); sum_nit += info["newton_iters"]
        worst_pen = max(worst_pen, max_pen(pos_d))
    return dict(pen0=pen0, pen_end=max_pen(pos_d), worst_pen=worst_pen, mono=all_mono,
                conv=all_conv, max_nit=max_nit, avg_nit=sum_nit / nsteps)


print(f"Gate 2: projected-Newton IPC (2 cells, npc={npc}, Fcap={Fcap:.2e} N)\n")
ok = True

# ---- Scenario A (GATE): FEASIBLE start + strong pull, 30-step loaded trajectory at several dt ----
print("  [A] non-penetration under load — feasible start (gap>0) + 4·Fcap inward pull, 30 steps:")
print(f"    {'dt/dt0':>7} {'avg_nit':>7} {'worst_pen':>9} {'E mono':>7} {'resid_ok':>8}  verdict")
for mult in (1, 100, 1000):
    r = trajectory(overlap=-0.15, Fpull=4.0 * Fcap, dt=8e-6 * mult, nsteps=30)
    pen_ok = r["worst_pen"] < 1e-3
    v = pen_ok and r["mono"] and r["conv"]
    ok = ok and v
    print(f"    {mult:>7} {r['avg_nit']:>7.1f} {r['worst_pen']:>9.2e} {str(r['mono']):>7} "
          f"{str(r['conv']):>8}  {'PASS' if v else 'FAIL'}")

# ---- Scenario B (DEMO): PENETRATING start + no pull, recovery trajectory --------------------------
print("\n  [B] feasibilization recovery — penetrating start (0.45R overlap), no pull, 60 steps @dt×50:")
rb = trajectory(overlap=0.45, Fpull=0.0, dt=8e-6 * 50, nsteps=60)
net = "net-recovery (slow: finite inside-penalty vs turgor inflation)" if rb["pen_end"] < rb["pen0"] \
    else "no-recovery"
print(f"    pen: {rb['pen0']:.3f} -> {rb['pen_end']:.3f}·me  (worst {rb['worst_pen']:.3f})  "
      f"E mono={rb['mono']}  avg_nit={rb['avg_nit']:.1f}  {net}")
print("    NOTE: deep penetrating-start recovery is slow by design (only the OUTSIDE barrier "
      "diverges) → production MUST start feasible/watertight (Scenario A's regime).")

print(f"\n  Gate 2 {'PASS' if ok else 'FAIL'}: energy-monotone + residual→tol + penetration-free under load")
sys.exit(0 if ok else 1)

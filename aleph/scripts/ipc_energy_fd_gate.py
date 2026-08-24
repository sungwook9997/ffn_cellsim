"""Gate 1: central-difference ∇E == −F for barrier / edge / turgor (energy↔force consistency)."""
import numpy as np, sys
import warp as wp
wp.init()
from aleph.dcm.geometry import icosphere_mesh
from aleph.dcm.dcm_ipc_energy import (ipc_barrier_energy_kernel, edge_energy_kernel, turgor_energy_kernel)
from aleph.dcm.dcm_contact_implicit_warp import nearest_face_ipc_kernel
from aleph.dcm.dcm_warp_hybrid import _bond_accumulate
from aleph.dcm.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel

dev = sys.argv[1] if len(sys.argv) > 1 else "cuda:0"
R = 7.5e-6
v0, edges0, faces0 = icosphere_mesh(R, subdivisions=1)
npc = v0.shape[0]

# two cells, slightly OVERLAPPING along x so the barrier is active
c0 = v0.copy(); c1 = v0.copy() + np.array([1.85 * R, 0, 0])   # centre spacing 1.85R (<2R → contact)
pos = np.concatenate([c0, c1], 0)
N = pos.shape[0]
faces = np.concatenate([faces0, faces0 + npc], 0).astype(np.int32)
edges = np.concatenate([edges0, edges0 + npc], 0).astype(np.int32)
cof = np.concatenate([np.zeros(npc), np.ones(npc)]).astype(np.int32)
fcell = cof[faces[:, 0]].astype(np.int32)
rng = np.random.default_rng(0)
pos = pos + rng.normal(0, 0.02 * R, pos.shape)          # jitter so no degenerate ties

pos_d = wp.array(pos, dtype=wp.vec3d, device=dev)
faces_d = wp.array(faces, dtype=wp.int32, device=dev)
edges_d = wp.array(edges, dtype=wp.int32, device=dev)
cof_d = wp.array(cof, dtype=wp.int32, device=dev)
fcell_d = wp.array(fcell, dtype=wp.int32, device=dev)
r0_edge = np.linalg.norm(pos[edges[:, 0]] - pos[edges[:, 1]], axis=1) * 0.9   # nonzero tension
r0_d = wp.array(r0_edge, dtype=wp.float64, device=dev)
mean_edge = float(np.linalg.norm(pos[edges[:, 0]] - pos[edges[:, 1]], axis=1).mean())
rep = 2.0e8; d_hat = 0.30 * mean_edge; k_edge = 1.0e-6
K_vol = 7.73e5; dP0 = 40.0
V0c = np.full(2, (4.0 / 3.0) * np.pi * R**3)
V0_d = wp.array(V0c, dtype=wp.float64, device=dev)

def evec():
    return wp.zeros(1, dtype=wp.float64, device=dev)

def barrier_E(P):
    pd = wp.array(P, dtype=wp.vec3d, device=dev)
    nf = wp.array(P.astype(np.float32), dtype=wp.vec3, device=dev)
    cent = ((P[faces[:, 0]] + P[faces[:, 1]] + P[faces[:, 2]]) / 3.0).astype(np.float32)
    cf = wp.array(cent, dtype=wp.vec3, device=dev)
    grid = wp.HashGrid(16, 16, 16, device=dev); grid.build(points=cf, radius=float(d_hat + mean_edge))
    e = evec()
    wp.launch(ipc_barrier_energy_kernel, dim=N, inputs=[grid.id, nf, pd, cof_d, faces_d, fcell_d,
              wp.float32(d_hat + mean_edge), wp.float64(rep), wp.float64(d_hat), e], device=dev)
    wp.synchronize_device(dev); return float(e.numpy()[0]), grid, nf

def barrier_F(P):
    pd = wp.array(P, dtype=wp.vec3d, device=dev)
    _, grid, nf = barrier_E(P)
    f = wp.zeros(N, dtype=wp.vec3d, device=dev)
    cnk = wp.zeros(N, dtype=wp.float64, device=dev); cnn = wp.zeros(N, dtype=wp.vec3d, device=dev)
    wp.launch(nearest_face_ipc_kernel, dim=N, inputs=[grid.id, nf, pd, cof_d, faces_d, fcell_d,
              wp.float32(d_hat + mean_edge), wp.float64(rep), wp.float64(d_hat), wp.float64(0.0),
              wp.float64(0.0), f, cnk, cnn], device=dev)
    wp.synchronize_device(dev); return f.numpy()

def edge_E(P):
    pd = wp.array(P, dtype=wp.vec3d, device=dev); e = evec()
    wp.launch(edge_energy_kernel, dim=edges.shape[0], inputs=[pd, edges_d, wp.float64(k_edge), r0_d, e], device=dev)
    wp.synchronize_device(dev); return float(e.numpy()[0])

def edge_F(P):
    pd = wp.array(P, dtype=wp.vec3d, device=dev); f = wp.zeros(N, dtype=wp.vec3d, device=dev)
    wp.launch(_bond_accumulate, dim=edges.shape[0], inputs=[pd, edges_d, wp.float64(k_edge), r0_d, f], device=dev)
    wp.synchronize_device(dev); return f.numpy()

def turgor_E(P):
    pd = wp.array(P, dtype=wp.vec3d, device=dev)
    Vc = wp.zeros(2, dtype=wp.float64, device=dev)
    wp.launch(dcm_volume_kernel, dim=faces.shape[0], inputs=[pd, faces_d, fcell_d, Vc], device=dev)
    e = evec()
    wp.launch(turgor_energy_kernel, dim=2, inputs=[Vc, V0_d, wp.float64(K_vol), wp.float64(dP0), e], device=dev)
    wp.synchronize_device(dev); return float(e.numpy()[0])

def turgor_F(P):
    pd = wp.array(P, dtype=wp.vec3d, device=dev)
    Vc = wp.zeros(2, dtype=wp.float64, device=dev)
    wp.launch(dcm_volume_kernel, dim=faces.shape[0], inputs=[pd, faces_d, fcell_d, Vc], device=dev)
    Vcn = Vc.numpy(); dP = dP0 + K_vol * (V0c - Vcn) / V0c
    dP_d = wp.array(dP, dtype=wp.float64, device=dev)
    f = wp.zeros(N, dtype=wp.vec3d, device=dev)
    wp.launch(dcm_turgor_force_kernel, dim=faces.shape[0], inputs=[pd, faces_d, fcell_d, dP_d, f], device=dev)
    wp.synchronize_device(dev); return f.numpy()

eps = 1e-11
for name, Efn, Ffn in [("edge", edge_E, edge_F), ("turgor", turgor_E, turgor_F), ("barrier", barrier_E, barrier_F)]:
    E0fn = (lambda P: Efn(P)[0]) if name == "barrier" else Efn
    F = Ffn(pos)
    v = rng.normal(0, 1, pos.shape); v /= np.linalg.norm(v)
    Ep = E0fn(pos + eps * v); Em = E0fn(pos - eps * v)
    fd = (Ep - Em) / (2 * eps)
    analytic = -np.sum(F * v)                      # -F·v = ∇E·v
    rel = abs(fd - analytic) / max(abs(analytic), 1e-30)
    print(f"{name:8s}: FD(∇E·v)={fd:+.4e}  -F·v={analytic:+.4e}  rel_err={rel:.2e}  {'PASS' if rel < 1e-3 else 'FAIL'}")

"""GPU-porting check: exercise every new compartment kernel device-resident on cuda:0 (Warp) + verify a
kernel result matches the CPU closed form (device parity). Confirms the units are GPU-native, not CPU-only."""
import sys
import numpy as np
import warp as wp

DEV = sys.argv[1] if len(sys.argv) > 1 else "cuda:0"
wp.init()
print(f"=== GPU-port check on {DEV} ===")

# 1. Microtubule aster — bending kernel (reused) on cuda
from aleph.laws.microtubule import build_microtubule_aster
from aleph.laws.forces_warp import cytosim_bending_kernel, _per_triple_alpha
from aleph.laws.network_warp import _zero
ast = build_microtubule_aster((0, 0, 0), n_mt=40, L_mt_um=6.0, seg_um=0.5)
tri = np.ascontiguousarray(ast.net.bend_triples, np.int32)
pos = wp.array(np.ascontiguousarray(ast.net.pos, np.float64), dtype=wp.vec3d, device=DEV)
tri_d = wp.array(tri, dtype=wp.int32, ndim=2, device=DEV)
al_d = wp.array(np.ascontiguousarray(_per_triple_alpha(ast.net), np.float64), dtype=wp.float64, device=DEV)
f = wp.zeros(ast.net.n_nodes, dtype=wp.vec3d, device=DEV)
wp.launch(_zero, dim=ast.net.n_nodes, inputs=[f], device=DEV)
wp.launch(cytosim_bending_kernel, dim=tri.shape[0], inputs=[pos, tri_d, al_d, f], device=DEV)
print(f"[MT]        aster {ast.net.n_nodes} nodes, bending kernel on {DEV}: |f|max={np.abs(f.numpy()).max():.2e} OK")

# 2. Substrate anchor relax on cuda
from aleph.laws.substrate import resolve_substrate, substrate_anchor_relax_kernel
s = resolve_substrate(E_pa=5e3)
anc = wp.array(np.zeros((100, 3)), dtype=wp.vec3d, device=DEV); rest = wp.array(np.zeros((100, 3)), dtype=wp.vec3d, device=DEV)
bd = wp.array(np.ones(100, np.int32), dtype=wp.int32, device=DEV); cf = wp.array(np.tile([50., 0, 0], (100, 1)), dtype=wp.vec3d, device=DEV)
for _ in range(2000):
    wp.launch(substrate_anchor_relax_kernel, dim=100, inputs=[anc, rest, bd, cf, wp.float64(s.k_sub), wp.float64(0.1 / s.k_sub)], device=DEV)
u = anc.numpy()[0, 0]; print(f"[substrate] anchor relax on {DEV}: u={u:.4f}µm vs F/k_sub={50/s.k_sub:.4f} {'PARITY OK' if abs(u-50/s.k_sub)<1e-3 else 'FAIL'}")

# 3. FA maturation kernels on cuda
from aleph.laws.fa_maturation import talin_vinculin_kernel, fa_growth_kernel, MaturationParams, talin_unfold_rate
p = MaturationParams()
pu = wp.array(np.zeros(50), dtype=wp.float64, device=DEV); nv = wp.array(np.zeros(50), dtype=wp.float64, device=DEV)
ld = wp.array(np.full(50, 8.0), dtype=wp.float64, device=DEV); bd2 = wp.array(np.ones(50, np.int32), dtype=wp.int32, device=DEV)
ar = wp.array(np.ones(50), dtype=wp.float64, device=DEV)
wp.launch(talin_vinculin_kernel, dim=50, inputs=[ld, pu, nv, bd2, wp.float64(p.ku0), wp.float64(p.dx_um), wp.float64(p.kT), wp.float64(p.k_refold), wp.float64(p.k_rec), wp.float64(p.k_diss), wp.float64(p.n_max), wp.float64(1e-4)], device=DEV)
wp.launch(fa_growth_kernel, dim=50, inputs=[ld, ar, wp.float64(p.kg0), wp.float64(p.kd), wp.float64(p.n_hill), wp.float64(p.fth), wp.float64(0.01)], device=DEV)
kf = pu.numpy()[0] / 1e-4; print(f"[FA-mat]    talin+FA-growth on {DEV}: k_unfold={kf:.3e} vs closed {talin_unfold_rate(8.0):.3e} {'PARITY OK' if abs(kf-talin_unfold_rate(8.0))/talin_unfold_rate(8.0)<0.01 else 'FAIL'}")

# 4. Piezo on cuda
from aleph.laws.piezo import resolve_piezo, piezo_popen_kernel, p_open
pz = resolve_piezo()
t = wp.array(np.array([pz.gamma_half_pn_um]), dtype=wp.float64, device=DEV); a2 = wp.array(np.ones(1), dtype=wp.float64, device=DEV)
po = wp.zeros(1, dtype=wp.float64, device=DEV); oa = wp.zeros(1, dtype=wp.float64, device=DEV)
wp.launch(piezo_popen_kernel, dim=1, inputs=[t, a2, wp.float64(pz.gamma_half_pn_um), wp.float64(pz.gamma_s_pn_um), po, oa], device=DEV)
print(f"[Piezo]     P_open kernel on {DEV}: {po.numpy()[0]:.4f} vs 0.5 {'PARITY OK' if abs(po.numpy()[0]-0.5)<1e-9 else 'FAIL'}")

# 5. Myosin minifilament on cuda
from aleph.laws.myosin_linear import resolve_myosin, minifilament_kernel
m = resolve_myosin()
mp = wp.array(np.array([[0., 0, 0], [1., 0, 0]]), dtype=wp.vec3d, device=DEV); lk = wp.array(np.array([[0, 1]], np.int32), dtype=wp.int32, ndim=2, device=DEV)
vv = wp.array(np.array([0.0]), dtype=wp.float64, device=DEV); fm = wp.zeros(2, dtype=wp.vec3d, device=DEV)
wp.launch(minifilament_kernel, dim=1, inputs=[mp, lk, vv, wp.float64(m.f_stall_pn), wp.float64(m.v0_um_s), fm], device=DEV)
print(f"[myosin]    minifilament on {DEV}: isometric |F|={np.linalg.norm(fm.numpy()[0]):.1f}pN vs stall {m.f_stall_pn:.0f} {'PARITY OK' if abs(np.linalg.norm(fm.numpy()[0])-m.f_stall_pn)<0.1 else 'FAIL'}")

# 6. Nuclear envelope on cuda
from aleph.laws.nucleus_envelope import resolve_nuclear_envelope, nuclear_envelope_kernel
ne = resolve_nuclear_envelope(5.25, subdivisions=3)
np_ = wp.array(ne.mesh.verts, dtype=wp.vec3d, device=DEV); fc = wp.array(ne.mesh.faces.astype(np.int32), dtype=wp.int32, ndim=2, device=DEV)
rupt = wp.zeros(ne.mesh.faces.shape[0], dtype=wp.int32, device=DEV); nf = wp.zeros(ne.mesh.verts.shape[0], dtype=wp.vec3d, device=DEV)
wp.launch(nuclear_envelope_kernel, dim=ne.mesh.faces.shape[0], inputs=[np_, fc, wp.float64(0.3), rupt, nf], device=DEV)
print(f"[nucleus]   envelope {ne.mesh.verts.shape[0]} nodes on {DEV}: |f|max={np.abs(nf.numpy()).max():.2e} OK")
print(f"=== ALL 6 new-unit kernels run device-resident on {DEV} ✓ ===")

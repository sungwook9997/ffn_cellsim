"""Validate cleave_cell: manifold, volume conservation, winding, separation."""
import sys; sys.path.insert(0, '/Users/sw1/ffn_cellsim')
import numpy as np
from ffn_sim.cell.dcm import icosphere_mesh
from ffn_sim.dcm.dcm_remesh import mesh_edges, enclosed_volume
from ffn_sim.dcm.dcm_cleave import cleave_cell

def check_daughter(pos, faces, fc, cid, label):
    f = faces[fc == cid]
    nodes = np.unique(f)
    e, ef, ea, bnd = mesh_edges(f)
    V, E, F = len(nodes), len(e), len(f)
    euler = V - E + F
    vol = enclosed_volume(pos, f)
    closed = (not bnd.any())
    print(f"  {label}: V={V} E={E} F={F}  Euler={euler}(want 2)  closed={closed}  vol={vol:.2f}  "
          f"{'OK' if (euler==2 and closed and vol>0) else 'FAIL'}")
    return euler == 2 and closed and vol > 0, vol, set(nodes.tolist())

verts, _edges, tris = icosphere_mesh(7.5, 2)
npc = verts.shape[0]
print(f"icosphere subdiv-2: {npc} nodes, {len(tris)} faces")

# global arrays: 1 active cell + dormant pool of 100 nodes
POOL = 200
pos = np.vstack([verts, np.tile([1e4, 1e4, 1e4], (POOL, 1))])
faces = tris.astype(np.int64).copy()
face_cell = np.zeros(len(faces), np.int64)
cof = np.concatenate([np.zeros(npc, np.int64), np.full(POOL, -1, np.int64)])

V0 = enclosed_volume(pos, faces)
print(f"mother volume V0 = {V0:.3f}")

# cleave through centroid along a tilted axis (general position)
p0 = verts.mean(0)
n = np.array([0.3, 0.5, 1.0]); n /= np.linalg.norm(n)
pos2, faces2, fc2, cof2, info = cleave_cell(
    pos, faces, face_cell, cof, cell_id=0, daughter_id=1, p0=p0, n=n)

print(f"ring nodes={info['n_ring']}, dormant used={info['nodes_used']}")
ok_m, vm, nodes_m = check_daughter(pos2, faces2, fc2, 0, "mother(-)")
ok_d, vd, nodes_d = check_daughter(pos2, faces2, fc2, 1, "daughter(+)")

print(f"\nvolume conservation: V0={V0:.3f}  vm+vd={vm+vd:.3f}  err={abs(V0-(vm+vd))/V0*100:.4f}%")
print(f"daughters share nodes? {len(nodes_m & nodes_d)} shared (want 0)")
print(f"node-count: pool had {POOL}, used {info['nodes_used']}, remaining dormant={int((cof2<0).sum())}")

passed = (ok_m and ok_d and abs(V0-(vm+vd))/V0 < 1e-3 and len(nodes_m & nodes_d)==0)
print(f"\n{'='*40}\n{'ALL PASS' if passed else 'FAIL'}\n{'='*40}")

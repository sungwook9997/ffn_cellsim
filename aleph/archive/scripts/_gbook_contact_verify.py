"""GPU contact-fix verification: does --conservative-contact close the cortex 'floating above face' gap?
Runs confluent N=8 full-compartment (nucleus + membrane tension + physiological turgor), conservative
tent vs default penalty, measures the min inter-cell surface gap at the last frame. Run on gbook cuda:0.

Physiological (mcf7_baseline): turgor 133 Pa, eta 65.9, E_nuc 4700, R_nuc 0.25R, ratio_lamin 1.4,
cortical/membrane gamma 5e-4 N/m (band-centre that gives Pi0=133 via Young-Laplace).
"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.expanduser("~/ff_scratch"))
from scipy.spatial import cKDTree
from aleph.dcm.dcm_warp_decohesion import run_decohesion

OUT = os.path.expanduser("~/ff_scratch/_contact_out")
os.makedirs(OUT, exist_ok=True)
PHYS = dict(nucleus=True, E_nuc=4700.0, R_nuc_factor=0.25, ratio_lamin=1.4,
            surface_tension=True, gamma_surf=5.0e-4, turgor=None)  # turgor via dP0 default 133 in kernel


def min_gap(npz):
    d = np.load(npz); fr = d["frames"]; cof = d["cof"]; nc = int(cof.max()) + 1
    pos = fr[-1]; gaps = []
    for c in range(nc):
        my = np.where(cof == c)[0]
        oth = np.where(cof != c)[0]
        if my.size < 4 or oth.size < 1:
            continue
        gaps.append(cKDTree(pos[oth]).query(pos[my])[0].min())
    return float(np.median(gaps)) * 1e6, float(np.min(gaps)) * 1e6  # um


for tag, cons in [("penalty_default", False), ("conservative_FIX", True)]:
    t0 = time.time()
    npz = f"{OUT}/contact_{tag}.npz"
    r = run_decohesion(n_cells=8, subdiv=2, steps=8000, frames=6, device="cuda:0", builder="confluent",
                       inset=0.02, v0_from_init=True,
                       nucleus=True, E_nuc=4700.0, R_nuc_factor=0.25, ratio_lamin=1.4,
                       surface_tension=True, gamma_surf=5.0e-4,
                       conservative_contact=cons, warmup=1000, save_frames=npz)
    med, mn = min_gap(npz)
    print(f"[{tag}] {time.time()-t0:.0f}s  median inter-cell surface gap = {med:.3f} um  (min {mn:.3f})  "
          f"GATES compaction_x={r.get('compaction_x'):.4f}", flush=True)
print("CONTACT VERIFY DONE", flush=True)

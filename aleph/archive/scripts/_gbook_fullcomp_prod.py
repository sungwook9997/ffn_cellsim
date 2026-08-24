"""GPU-NATIVE production: N=400 full-compartment DCM spheroid on gbook A5000 (cuda:0).
Full physiological compartment stack + pressing contact:
  - plasma membrane / cortex : icosphere shell + membrane surface tension gamma=5e-4 N/m (band-centre)
  - cytoplasm                : turgor Pi0=133 Pa + K_vol + eta=65.9 Pa.s (MCF7, auto per-node drag)
  - nucleus                  : E_nuc=4700 Pa, R_nuc=0.25R, ratio_lamin=1.4 (MCF7, mcf7_baseline)
  - contact                  : CONSERVATIVE tent (energy min at contact) + tight confluent inset 0.01
                               so the cortex APPOSES the face (no c_rep floating gap)
Saves frames npz + measures the inter-cell surface gap (pressing check). No CPU, no bare shell.
"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.expanduser("~/ff_scratch"))
from scipy.spatial import cKDTree
from aleph.dcm.dcm_warp_decohesion import run_decohesion

OUT = os.path.expanduser("~/ff_scratch/_prod_out")
os.makedirs(OUT, exist_ok=True)
N = int(os.environ.get("NCELLS", "400"))
STEPS = int(os.environ.get("STEPS", "30000"))
INTEG = os.environ.get("INTEG", "baoab")          # 'implicit' = stable stiff contact
TAG = os.environ.get("TAG", "")
# CONTACT: 'penalty' (conservative tent, default) or 'ipc' (Li-2020 log-barrier, penetration-free by CCD).
# --ipc must NOT be combined with conservative_contact (pen=73 broken combo) and needs the implicit integrator.
CONTACT = os.environ.get("CONTACT", "penalty")
USE_IPC = CONTACT == "ipc"
# INIT_NPZ: restart from a pre-built confluent init (frames/faces/cof), e.g. the ghost-seed-corrected
# init from dcm_ghost_seed_init_prototype — lets us run the audit#9 size-gradient fix end-to-end WITHOUT
# touching the shared confluent builder. When set, the builder is bypassed (run_decohesion loads the mesh).
INIT_NPZ = os.environ.get("INIT_NPZ", "") or None
npz = f"{OUT}/fullcomp_n{N}{TAG}.npz"

t0 = time.time()
r = run_decohesion(
    n_cells=N, subdiv=2, steps=STEPS, frames=20, device="cuda:0",
    builder="confluent", inset=0.01, v0_from_init=True, lloyd_iters=6,
    init_npz=INIT_NPZ,                             # ghost-corrected init when set (else fresh confluent build)
    conservative_contact=(not USE_IPC),           # penalty tent (energy min at d~0) UNLESS ipc
    ipc=USE_IPC,                                   # Li-2020 log-barrier contact (CCD, penetration-free)
    integrator=("implicit" if USE_IPC else INTEG), # --ipc needs implicit
    nucleus=True, E_nuc=4700.0, R_nuc_factor=float(os.environ.get("RNUC", "0.25")),  # RNUC=0.7 = lit MCF7 (Moore 2016, N:C 1.9); 0.25 = old ratified
    ratio_lamin=1.4,   # nucleus compartment
    surface_tension=True, gamma_surf=5.0e-4,      # plasma-membrane tension
    diff_tension=(os.environ.get("DIFF", "1") == "1"),  # DAH: contact faces lower tension -> apposition
    warmup=1500, save_frames=npz)

d = np.load(npz); fr = d["frames"]; cof = d["cof"]; nc = int(cof.max()) + 1
pos = fr[-1]; gaps = []
for c in range(nc):
    my = np.where(cof == c)[0]; oth = np.where(cof != c)[0]
    if my.size < 4 or oth.size < 1:
        continue
    gaps.append(cKDTree(pos[oth]).query(pos[my])[0].min())
gaps = np.array(gaps) * 1e6

# per-cell asphericity (faceting) + mean radius
def asph(p):
    dd = p - p.mean(0); ev = np.sort(np.linalg.eigvalsh(dd.T @ dd))[::-1]; s = ev.sum()
    return (ev[0] - 0.5 * (ev[1] + ev[2])) / s if s > 0 else 0.0
A = np.mean([asph(pos[cof == c]) for c in range(nc)])
print(f"[PROD N={N}] {time.time()-t0:.0f}s  cells={nc}  inter-cell gap median={np.median(gaps):.3f}um "
      f"(min {gaps.min():.3f}, p90 {np.percentile(gaps,90):.3f})  asph={A:.4f}  "
      f"vv0={r.get('vv0_final', r.get('aa0_final','?'))}  -> {npz}", flush=True)
print("PROD DONE", flush=True)

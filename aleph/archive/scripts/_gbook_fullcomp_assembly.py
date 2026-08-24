"""GPU-NATIVE production: full-compartment DCM spheroid ASSEMBLED FROM DISPERSED CELLS via
fine-grained CADHERIN catch-bonds (E1, Rakshit KU-4.2 — de-cohesion emergent), on gbook A5000.

Unlike _gbook_fullcomp_prod.py (which starts from a CONFLUENT pre-packed init and holds together
by interfacial-tension adhesion), this run starts cells in a DISPERSED spherical cluster
(builder='fcc', gap*R spacing) and lets EXPLICIT cadherin trans-dimer bonds pull them into a
cohesive spheroid — the actual assembly process, with the full compartment stack:
  - plasma membrane / cortex : surface tension gamma=5e-4 N/m (+ differential tension)
  - cytoplasm                : turgor + K_vol + eta=65.9 Pa.s (MCF7 per-node drag)
  - nucleus                  : E_nuc=4700 Pa, R_nuc=RNUC*R (RNUC=0.7 lit MCF7), bilinear chromatin/lamin
  - cell-cell COHESION       : cadherin=True (catch-bond ensemble; contact runs repulsion-only)
Env: NCELLS, STEPS, RNUC, GAP, TAG, INTEG(baoab default; contact=penalty).
"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.expanduser("~/ff_scratch"))
from scipy.spatial import cKDTree
from aleph.dcm.dcm_warp_decohesion import run_decohesion

OUT = os.path.expanduser("~/ff_scratch/_prod_out")
os.makedirs(OUT, exist_ok=True)
N     = int(os.environ.get("NCELLS", "400"))
STEPS = int(os.environ.get("STEPS", "60000"))          # assembly from dispersed needs more steps than confluent
RNUC  = float(os.environ.get("RNUC", "0.7"))
GAP   = float(os.environ.get("GAP", "2.05"))           # centre spacing gap*R (>2 => cells start with a gap → aggregate)
INTEG = os.environ.get("INTEG", "baoab")
BUILDER  = os.environ.get("BUILDER", "fcc")            # 'fcc' dispersed cluster | 'confluent' Voronoi space-filling spheroid
INIT_NPZ = os.environ.get("INIT_NPZ", "") or None      # for confluent: a pre-built Voronoi init mesh (ghost_TIGHT)
TAG   = os.environ.get("TAG", "")
# --- fidelity-audit fixes (PI-ratified 2026-07-07), all env-overridable ---
ENUC   = float(os.environ.get("ENUC", "399"))      # MCF7 nucleus in-situ ~399 Pa (audit#19); was 4700 (8-12x too stiff)
BUNDLE = float(os.environ.get("BUNDLE", "100"))    # N_cad per contact = rho_cad*A_junction (KB-4.1/4.11/4.17, Buckley2014);
                                                   # bundles the single-molecule catch bond → interface force 1-10nN (was 1)
GAMMA  = float(os.environ.get("GAMMA", "5.0e-4"))  # cortical/surface tension [N/m] — SWEEP variable (band 0.1-1 mN/m ↔ MCF7 1e-2)
INSET  = float(os.environ.get("INSET", "0.01" if os.environ.get("BUILDER", "fcc") == "confluent" else "0.0"))  # >0.01 → cells start non-overlapping (G2 fix)
IPC    = os.environ.get("IPC", "0") == "1"         # log-barrier contact → guarantees non-penetration (G2 fix)
npz   = f"{OUT}/fullcomp_assembly_n{N}{TAG}.npz"

t0 = time.time()
r = run_decohesion(
    n_cells=N, subdiv=2, steps=STEPS, frames=24, device="cuda:0",
    builder=BUILDER, gap=GAP, init_npz=INIT_NPZ,
    inset=INSET,
    v0_from_init=(BUILDER == "confluent"), lloyd_iters=6,   # confluent → Voronoi space-filling spheroid
    substrate_wetting=False, use_substrate_well=False,  # no flat floor
    # U-BOTTOM ULA WELL + GRAVITY = the real spheroid-formation route: cells sediment (Δρ) into a
    # rounded bowl that funnels them to the centre where cadherin catch-bonds cohere them → they
    # ACTUALLY COME TOGETHER from dispersed into a spheroid (env AGGREGATE=1 to enable).
    gravity=(os.environ.get("AGGREGATE", "0") == "1"),
    ubottom=(os.environ.get("AGGREGATE", "0") == "1"),
    ubottom_r_factor=float(os.environ.get("UBR", "1.05")),  # tight bowl → funnels cells inward as they sediment
    ubottom_k=float(os.environ.get("UBK", "5.0")),          # bowl-wall stiffness
    delta_rho=55.0,                                         # MCF7 − medium (SimuCell3D)
    cadherin=True, cad_bundle=BUNDLE,                  # E1 explicit cadherin catch-bonds, N_cad-bundled (de-cohesion emergent)
    conservative_contact=False, ipc=IPC,               # IPC=1 → log-barrier non-penetration; else penalty repulsion-only
    integrator=INTEG,
    accel_dt=(float(os.environ["ACCEL_DT"]) if os.environ.get("ACCEL_DT") else None),  # implicit-only: larger dt → fewer steps
    nucleus=True, E_nuc=ENUC, R_nuc_factor=RNUC, ratio_lamin=1.4,
    surface_tension=True, gamma_surf=GAMMA, diff_tension=True,
    warmup=2000, save_frames=npz)

d = np.load(npz); fr = d["frames"]; cof = d["cof"]; nc = int(cof.max()) + 1
pos = fr[-1]
# aggregation quality: inter-cell nearest gap + cluster compaction (Rg first vs last frame)
gaps = []
for c in range(nc):
    my = np.where(cof == c)[0]; oth = np.where(cof != c)[0]
    if my.size < 4 or oth.size < 1:
        continue
    gaps.append(cKDTree(pos[oth]).query(pos[my])[0].min())
gaps = np.array(gaps) * 1e6
def com_rg(p, cofa, nca):
    cen = np.array([p[cofa == c].mean(0) for c in range(nca)])
    com = cen.mean(0); return float(np.sqrt(((cen - com) ** 2).sum(1).mean())) * 1e6
rg0, rg1 = com_rg(fr[0], cof, nc), com_rg(fr[-1], cof, nc)
print(f"[ASSEMBLY N={N}] {time.time()-t0:.0f}s  cells={nc}  cadherin=ON  "
      f"inter-cell gap median={np.median(gaps):.3f}um (min {gaps.min():.3f})  "
      f"cluster Rg {rg0:.1f}->{rg1:.1f}um (compaction {rg1/max(rg0,1e-9):.3f}x)  "
      f"vv0={r.get('vv0_final','?')}  -> {npz}", flush=True)
print("ASSEMBLY DONE", flush=True)

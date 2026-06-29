import numpy as np, tempfile, os
from ffn_sim.dcm.dcm_warp_decohesion import run_decohesion
print("2-cell contact equilibrium at gap 2.05 (WITHIN adhesion range), rep sweep, adh=5e7:")
print(f"{'rep':>8} {'NN_eq(R)':>9} {'pen_fin':>8} {'V/V0':>6} {'verdict':>12}")
for rep in [4e7, 1e8, 2e8, 5e8, 1e9, 2e9]:
    tmp=os.path.join(tempfile.gettempdir(),f"_2c_{rep:.0e}.npz")
    out=run_decohesion(n_cells=2, subdiv=2, steps=100, frames=2, device="cpu",
        dt=8e-6, warmup=300, settle_steps=6000, settle_frames=0, gap=2.05,
        rep_strength=rep, adh_strength=5e7, surface_tension=True, gamma_surf=1e-4,
        substrate_wetting=False, use_substrate_well=False, lamellipodium=False,
        builder="fcc", integrator="implicit", accel_dt=8e-4, save_frames=tmp)
    d=np.load(tmp,allow_pickle=True); P=d["frames"][-1].astype(float); cof=d["cof"]
    cs=np.unique(cof[cof>=0]); cen=[P[cof==c].mean(0) for c in cs]
    nn=np.linalg.norm(cen[1]-cen[0])/7.5e-6
    v="JUST-TOUCH" if 1.92<=nn<=2.1 else ("OVERLAP" if nn<1.92 else "GAP")
    print(f"{rep:8.0e} {nn:9.2f} {out['pen_frac_final']:8.2f} {out['vv0_final']:6.3f} {v:>12}")

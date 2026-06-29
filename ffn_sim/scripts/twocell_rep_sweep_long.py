"""2-cell LONG-settle rep sweep — does the contact LEAK (real interpenetration that grows
with softer rep) or is the midplane crossing just flat-junction mesh interleaving?

The prior session called it "rigid-sphere overlap, do NOT flatten" from a SHORT run + a 3D
trisurf render (which hides flattening). This reaches equilibrium and uses the cross-section
diagnostic (twocell_overlap_diag) so we can SEE the junction and READ oblate/deep-cross.
"""
from __future__ import annotations
import numpy as np, time
from ffn_sim.dcm.dcm_warp_decohesion import run_decohesion
from ffn_sim.scripts.twocell_overlap_diag import diagnose, render

SETTLE = 10000
print(f"2-cell long settle ({SETTLE}) rep sweep, adh=5e7, gap=2.05, gamma=1e-4")
print(f"{'rep':>8} {'NN/R':>6} {'Psi':>6} {'oblate':>7} {'cross':>6} {'deep/R':>7} {'inmesh':>6}")
rows = []
for rep in [4e7, 2e8, 1e9, 5e9]:
    tmp = f"/tmp/_2c_sweep_{rep:.0e}.npz"
    t0 = time.time()
    out = run_decohesion(n_cells=2, subdiv=2, steps=100, frames=2, device="cpu",
        dt=8e-6, warmup=300, settle_steps=SETTLE, settle_frames=0, gap=2.05,
        rep_strength=rep, adh_strength=5e7, surface_tension=True, gamma_surf=1e-4,
        substrate_wetting=False, use_substrate_well=False, lamellipodium=False,
        builder="fcc", integrator="implicit", accel_dt=8e-4, save_frames=tmp)
    r = diagnose(tmp)
    print(f"{rep:8.0e} {r['NN_over_R']:6.2f} {r['psi0']:6.3f} "
          f"{r['oblate_axial_over_lateral']:7.3f} {r['cross_midplane_nodes']:6d} "
          f"{r['deep_cross_R']:7.3f} {r.get('nodes_inside_other_mesh','NA'):>6} "
          f"  ({time.time()-t0:.0f}s)")
    render(r, f"ffn_sim/outputs/warp_decohesion/figs/twocell_long_rep{rep:.0e}.png")
    rows.append((rep, r))

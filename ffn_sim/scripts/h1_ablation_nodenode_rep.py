"""H1 ablation: does dropping the (double-counted) node-NODE cohesion repulsion materially shift the
2-cell equilibrium? Adversarial verification (wf_ad040f4c) showed node-NODE rep = ~56% of node-FACE
rep (NOT my erroneous 0.2%), so ~36% of junction excluded volume is double-counted. This runs the
gentle 2-cell with node-NODE repulsion ON (baseline) vs OFF (node-FACE owns excluded volume) by
monkeypatching cohesion_grid_kernel to a rep-stripped variant — measuring the real NN/oblate/pen shift.
"""
from __future__ import annotations
import numpy as np
import warp as wp
import ffn_sim.warp_port.dcm_warp_decohesion as drv
from ffn_sim.warp_port.dcm_contact_warp import closest_bary  # noqa (kernel dep parity)
from ffn_sim.scripts.twocell_overlap_diag import diagnose

wp.init()


@wp.kernel
def cohesion_grid_kernel_norep(
    grid: wp.uint64, qpts: wp.array(dtype=wp.vec3), pos: wp.array(dtype=wp.vec3d),
    cof: wp.array(dtype=wp.int32), radius: wp.float32,
    r_contact: wp.float64, c_adh: wp.float64,
    rep: wp.float64, omega: wp.float64, A: wp.float64, force_cap: wp.float64,
    force: wp.array(dtype=wp.vec3d)):
    """cohesion_grid_kernel with the node-NODE REPULSION branch removed (adhesion only).
    node-FACE contact_grid_kernel then solely owns excluded volume → no double-count."""
    i = wp.tid()
    c1 = cof[i]
    if c1 < wp.int32(0):
        return
    z = wp.float64(0.0)
    half = wp.float64(0.5) * c_adh
    pi = pos[i]
    acc = wp.vec3d(z, z, z)
    q = wp.hash_grid_query(grid, qpts[i], radius)
    j = wp.int32(0)
    while wp.hash_grid_query_next(q, j):
        if j != i and cof[j] >= wp.int32(0) and cof[j] != c1:
            r_vec = pi - pos[j]
            d = wp.length(r_vec)
            fmag = z
            if d < r_contact:
                fmag = z                                       # ABLATED: no node-node repulsion
            else:
                if d < c_adh:
                    tent = c_adh - d
                    if d < half:
                        tent = d
                    fmag = -omega * A * tent
            if fmag > force_cap:
                fmag = force_cap
            if fmag < -force_cap:
                fmag = -force_cap
            if d > z and fmag != z:
                acc = acc + r_vec * (fmag / d)
    force[i] = acc


def run(label, save):
    out = drv.run_decohesion(n_cells=2, subdiv=2, steps=200, frames=2, device="cpu",
        dt=8e-6, warmup=300, settle_steps=10000, settle_frames=0, gap=2.05,
        rep_strength=2e8, adh_strength=5e7, surface_tension=True, gamma_surf=1e-4,
        substrate_wetting=False, use_substrate_well=False, lamellipodium=False,
        builder="fcc", integrator="implicit", accel_dt=8e-4, save_frames=save)
    r = diagnose(save)
    print(f"{label:22s} NN={r['NN_over_R']:.3f}R oblate={r['oblate_axial_over_lateral']:.4f} "
          f"Psi={r['psi0']:.4f} cross={r['cross_midplane_nodes']} deep={r['deep_cross_R']:.4f} "
          f"pen={out['pen_frac_final']:.4f} vv0={out['vv0_final']:.4f}")
    return r, out


if __name__ == "__main__":
    orig = drv.cohesion_grid_kernel
    print("=== H1 ablation: node-NODE cohesion repulsion ON (baseline) vs OFF ===")
    r_on, o_on = run("node-node rep ON ", "/tmp/_h1_on.npz")
    drv.cohesion_grid_kernel = cohesion_grid_kernel_norep      # node-FACE owns excluded volume
    r_off, o_off = run("node-node rep OFF", "/tmp/_h1_off.npz")
    drv.cohesion_grid_kernel = orig
    dNN = (r_off["NN_over_R"] - r_on["NN_over_R"]) / r_on["NN_over_R"] * 100
    dpen = o_off["pen_frac_final"] - o_on["pen_frac_final"]
    print(f"\nIMPACT of removing the double-counted node-NODE rep:")
    print(f"  NN: {r_on['NN_over_R']:.3f}R -> {r_off['NN_over_R']:.3f}R ({dNN:+.1f}%)  "
          f"(cells settle {'CLOSER' if dNN<0 else 'FARTHER'})")
    print(f"  oblate: {r_on['oblate_axial_over_lateral']:.4f} -> {r_off['oblate_axial_over_lateral']:.4f}")
    print(f"  pen: {o_on['pen_frac_final']:.4f} -> {o_off['pen_frac_final']:.4f} ({dpen:+.4f})")
    print(f"  deep_cross: {r_on['deep_cross_R']:.4f}R -> {r_off['deep_cross_R']:.4f}R")

"""Overnight proxy-free spreading run — the full mechanistic stack with the cohesion-bundle FORCE
fix (cad×40 / ecm×167, KB-anchored), measured with the peeling-aware spread_eval (so a 7-cell
basal-rim PEELING artifact is NOT misread as collective spread — the retracted A/A0→1.94 lesson).

NO proxy: wetting / substrate-well OFF; spreading must come from the ECM clutch (Pereverzev catch-
slip) + lamellipodium traction. NO outcome-tuning: every value is the committed lit/KB anchor.

Usage (on gbook A5000):
  python -m ffn_sim.scripts.run_spread_overnight --n 100 --steps 20000 --tag n100_proxyfree_spread
"""
from __future__ import annotations
import argparse
import json
import time
from ffn_sim.warp_port.dcm_warp_decohesion import run_decohesion

OUT = "/home/sungwook/ffn_phase_c/ffn_sim/outputs/warp_decohesion"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--settle", type=int, default=1000)
    ap.add_argument("--frames", type=int, default=60)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--tag", default="n100_proxyfree_spread")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--ipc", action="store_true",
                    help="M1: use IPC node-face contact (barrier+CCD) instead of the capped penalty — "
                         "the decisive validation that the strong cadherin/ecm bundle no longer tunnels "
                         "(penalty gave pen 3.1; IPC target pen<0.3, V/V0 stable, A/A0 unconfounded)")
    ap.add_argument("--k-vol", type=float, default=7.73e5, dest="k_vol",
                    help="osmotic bulk modulus (Pa). Default 7.73e5 = the 'fried-egg' near-incompressible "
                         "VOLUME LOCK (a contact band-aid that pins V/V0=1.000). Set to the pre-fried-egg "
                         "soft default 1e3 to REMOVE the lock so turgor (dP0) engages (V/V0 emergent) — "
                         "only valid PAIRED WITH --ipc (which provides the real over-compression guard).")
    ap.add_argument("--polarize", action="store_true",
                    help="apico-basal differential surface tension (Young-Dupre): basal faces wet, "
                         "apical keep cortical gamma. The directional-spread lever. S=w_cs-2*gamma>0 spreads.")
    ap.add_argument("--w-cs-polarize", type=float, default=2.85e-3, dest="w_cs_polarize",
                    help="basal adhesion energy J/m2 (lit MCF7 2.85e-3 -> S<0; raise >2*gamma for S>0 spread test)")
    ap.add_argument("--gamma-surf", type=float, default=1e-4, dest="gamma_surf",
                    help="cortical/apical surface tension N/m (lit MCF7 ~1e-2; default 1e-4)")
    ap.add_argument("--no-bundle", action="store_true",
                    help="turgor-engage / gentle-aggregation test: cadherin + ecm + lamellipodium OFF "
                         "(adh 5e7 node-face cohesion only), so V/V0 engagement under --k-vol is measured "
                         "WITHOUT the strong-bundle contact confound (this is the regime where pen~0.12).")
    ap.add_argument("--well", action="store_true",
                    help="substrate WELL/floor ON (z-anchor + rigid dish): the BOUND that the single-cell "
                         "fried-egg used to keep S>0 spreading finite. Tests if the spheroid spreads with "
                         "volume conservation (stiff k_vol) + substrate floor + polarization.")
    ap.add_argument("--ipc-dhat-factor", type=float, default=1.0, dest="ipc_dhat_factor",
                    help="IPC barrier activation d_hat=factor*c_rep. 1.0=IPC barrier; 0.01=SimuCell3D-style penalty-only (A/B)")
    ap.add_argument("--project", action="store_true",
                    help="M1 lever#3: SimuCell3D-style position-based PROJECTION hard-constraint (post-step "
                         "geometric non-penetration). A few Jacobi sweeps push penetrating nodes out to "
                         "c_rep clearance regardless of force magnitude — clamps pen where penalty/barrier "
                         "plateau at ~1.5-2.6. Compatible with the penalty contact (not paired with --ipc).")
    ap.add_argument("--proj-iter", type=int, default=8, dest="proj_iter",
                    help="projection Jacobi sweeps per step (self-test converged 3.67*c_rep penetration in ~6-12; "
                         "8 = margin since the momentum-conserving 50/50 split moves the node only half/sweep)")
    ap.add_argument("--proj-omega", type=float, default=0.7, dest="proj_omega",
                    help="projection relaxation (0.5-1.0; 0.7 default, lower = gentler/more stable)")
    ap.add_argument("--proj-gap-factor", type=float, default=1.0, dest="proj_gap_factor",
                    help="projection target clearance = factor*c_rep. 1.0=full c_rep shell (fights adhesion -> "
                         "diverged at accel_dt 8e-4); ~0.05 = near-zero gap = PURE non-penetration (only acts on "
                         "actual overlap, does NOT fight the force equilibrium) — the stability diagnostic")
    ap.add_argument("--accel-dt", type=float, default=8e-4, dest="accel_dt",
                    help="implicit accel dt (default 8e-4 = 100x base). Lower (8e-5/8e-6) for the stiff "
                         "polarized+substrate regime that diverges/crawls at 8e-4 (Colab-sweepable across GPUs)")
    ap.add_argument("--cfl-limit", type=float, default=0.0, dest="cfl_limit",
                    help="A3 adaptive substepping: cap per-step node displacement < cfl_limit*c_rep "
                         "(0.3 stabilizes the deformable+polarized regime that diverges at fixed accel_dt)")
    ap.add_argument("--rep", type=float, default=2e8, dest="rep_strength",
                    help="node-face repulsion (default 2e8 stiff). 4e7 = lit-anchored SOFT for physiological "
                         "AGGREGATION (cells deform-pack instead of staying rigid round spheres)")
    ap.add_argument("--adh", type=float, default=5e7, dest="adh_strength",
                    help="node-face adhesion / cohesion (default 5e7 lit-anchored)")
    ap.add_argument("--gap", type=float, default=2.05, dest="gap",
                    help="initial inter-cell center spacing factor (2.05 default; 2.3 = gapped start, "
                         "lets turgor+cohesion compact rather than starting pre-overlapped)")
    ap.add_argument("--coupling", action="store_true",
                    help="node-FACE CONTINUOUS adhesion ON (coh_adh=adh_strength) — the flat-interface "
                         "adhesion that, with surface tension, facets cells into space-filling polyhedra "
                         "(vs sparse node-node point cohesion). Physiological-baseline: should be ON in production.")
    ap.add_argument("--ubottom", action="store_true",
                    help="ULA U-bottom confinement: cells are held in a non-adhesive hemispherical bowl "
                         "(geometric confinement only; the surface never grips). Independent of --well.")
    ap.add_argument("--filopodia", action="store_true",
                    help="explicit filopodia finger protrusions (node-FACE + node-plane tip adhesions). "
                         "ON automatically under --ula (cell-cell junction formation).")
    ap.add_argument("--lamellipodium", action="store_true",
                    help="per-cell advancing-anchor lamellipodium crawl. ON automatically under --ula.")
    ap.add_argument("--division", action="store_true",
                    help="C7 rim-cell proliferation ON. Pair with --div-real-hours for TIME-CONSISTENT "
                         "division (cells divide at the MCF7 cycle rate over the run's represented real time).")
    ap.add_argument("--div-real-hours", type=float, default=0.0, dest="div_real_hours",
                    help="real biological hours this run REPRESENTS (e.g. 24). >0 => time-consistent division: "
                         "p_div = S*dt*div_every/T_cycle so cells divide ~div_real_hours/T_cycle times. 0 = bare div_rate.")
    ap.add_argument("--div-t-cycle", type=float, default=24.0, dest="div_t_cycle_h",
                    help="cell-cycle / doubling time in HOURS (MCF7 ~24h). Used by time-consistent division.")
    ap.add_argument("--accel-real-hours", type=float, default=0.0, dest="accel_real_hours",
                    help="UNIFIED time-acceleration: the real biological hours this ONE feasible run "
                         "REPRESENTS (e.g. 24-48 for spheroid formation). >0 => a single factor "
                         "S=accel_real_hours*3600/(dt*steps) multiplies ALL slow biological RATES "
                         "(filopodia v_poly/p_seed, cadherin k_on/k_off, lamellipodium front, division) "
                         "while the MECHANICAL forces stay at the faithful dt. div-real-hours defaults to "
                         "this so division is not double-accelerated. 0 = native rates (S=1).")
    ap.add_argument("--necrosis", action="store_true",
                    help="C8 3-zone necrosis ON. Depth-from-surface O2 proxy assigns prolif/quiescent/"
                         "necrotic; necrotic core softens turgor and division is gated to the proliferating "
                         "rim. Inert until the spheroid is large enough to develop an anoxic core.")
    ap.add_argument("--active-batch", type=int, default=50, dest="active_batch",
                    help="filopodia+lamellipodium host-update cadence (steps). 50 default; 200 = 4x fewer "
                         "expensive host PROBE updates (advance is batch-consistent so velocity is preserved). "
                         "The single biggest cheap speedup for the active-junction stack.")
    ap.add_argument("--gpu-probe", action="store_true", dest="gpu_probe",
                    help="GPU-resident filopodia PROBE (hash-grid Warp kernel, ~GPU-only). CPU-parity-validated.")
    ap.add_argument("--frozen-neighbors", action="store_true", dest="frozen_neighbors",
                    help="I-opt #2: cache the cohesion/contact neighbour set once per step + reuse it "
                         "in every implicit-CG matvec (no per-iter hash-grid query). Parity-safe (the "
                         "neighbours are frozen across a solve); removes the query-bound matvec cost.")
    ap.add_argument("--precond-diag", action="store_true", dest="precond_diag",
                    help="I-opt #1: analytic-diagonal Jacobi preconditioner for the implicit CG. Same "
                         "converged dx; helps only in heterogeneous-diagonal regimes.")
    ap.add_argument("--ula", action="store_true",
                    help="ULA spheroid formation: U-bottom bowl ON, flat substrate well + wetting + ECM "
                         "clutch OFF (non-adhesive surface), cadherin + filopodia + lamellipodium ON "
                         "(cell-cell adhesion is the ONLY adhesion). Post-centrifuge pellet → spheroid.")
    a = ap.parse_args()
    # --ula is the convenience preset; individual flags OR with it so they also work standalone.
    ula = a.ula
    ubottom = a.ubottom or ula
    filopodia = a.filopodia or ula
    lamellipodium = a.lamellipodium or ula
    if ula:
        # ULA = cell-cell adhesion ONLY (no ECM / no substrate grip); bowl confines geometrically.
        use_substrate_well = False
        substrate_wetting = False
        ecm_clutch = False
        cadherin = True
    else:
        # unchanged legacy behavior (proxy-free mechanistic stack)
        use_substrate_well = a.well
        substrate_wetting = False
        ecm_clutch = not a.no_bundle
        cadherin = not a.no_bundle
        # legacy stack drives lamellipodium off --no-bundle unless explicitly requested above
        lamellipodium = lamellipodium or (not a.no_bundle)
    npz = f"{a.out}/{a.tag}.npz"
    t0 = time.time()
    out = run_decohesion(
        n_cells=a.n, subdiv=2, steps=a.steps, frames=a.frames, device=a.device,
        dt=8e-6, warmup=200, settle_steps=a.settle, gap=a.gap,
        k_vol=a.k_vol,
        rep_strength=a.rep_strength, adh_strength=a.adh_strength,
        surface_tension=True, gamma_surf=a.gamma_surf, bending=True, edge_edge=True, nucleus=True,
        polarize=a.polarize, w_cs_polarize=a.w_cs_polarize, cfl_limit=a.cfl_limit,
        ipc_dhat_factor=a.ipc_dhat_factor,
        cadherin=cadherin, ecm_clutch=ecm_clutch,
        lamellipodium=lamellipodium, filopodia=filopodia, active_batch=a.active_batch,
        gpu_probe=a.gpu_probe,
        coupling=a.coupling,
        cad_bundle=40.0, ecm_bundle=167.0,
        substrate_wetting=substrate_wetting, use_substrate_well=use_substrate_well,
        ubottom=ubottom,                                      # ULA non-adhesive bowl confinement
        division=a.division, div_real_hours=a.div_real_hours, div_t_cycle_h=a.div_t_cycle_h,
        accel_real_hours=a.accel_real_hours,
        necrosis=a.necrosis,
        builder="fcc", integrator="implicit", accel_dt=a.accel_dt,
        frozen_neighbors=a.frozen_neighbors, precond_diag=a.precond_diag,
        ipc=a.ipc, project=a.project, proj_iter=a.proj_iter, proj_omega=a.proj_omega,
        proj_gap_factor=a.proj_gap_factor,
        save_frames=npz)
    dt = time.time() - t0
    rec = {"tag": a.tag, "n": a.n, "steps": a.steps, "elapsed_s": round(dt, 1),
           "aa0_final": out.get("aa0_final"), "aa0_peak": out.get("aa0_peak"),
           "vv0_final": out.get("vv0_final"), "pen_final": out.get("pen_frac_final"),
           "pen_peak": out.get("pen_frac_peak"), "drift_um": out.get("drift_final_um"),
           "npz": npz}
    with open(f"{a.out}/{a.tag}.json", "w") as f:
        json.dump(rec, f, indent=2)
    print("SPREAD_RUN_DONE", json.dumps(rec))
    # peeling-aware verdict inline
    try:
        from ffn_sim.scripts.spread_eval import evaluate, verdict
        m = evaluate(npz)
        print("SPREAD_EVAL", verdict(m),
              f"aa0={m['aa0_final']:.3f} maxZdrop={m['maxZ_drop_pct']:.1f}% "
              f"peel_idx={m['peeling_index']:.1f} Psi={m['psi_init']:.3f}->{m['psi_final']:.3f}")
    except Exception as e:
        print("SPREAD_EVAL_FAILED", repr(e))


if __name__ == "__main__":
    main()

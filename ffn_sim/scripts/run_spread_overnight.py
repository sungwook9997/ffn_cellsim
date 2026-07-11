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
from ffn_sim.dcm.dcm_warp_decohesion import run_decohesion

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
    ap.add_argument("--remesh-period", type=int, default=0, dest="remesh_period",
                    help="remesh cadence (steps): split stretched faces / collapse slivers so the SPREADING "
                         "cortex gets NEW area (nodes) instead of over-stretching a few nodes into ejected "
                         "spikes. The area-source that lets a cell flatten cleanly (0=off).")
    ap.add_argument("--coupling", action="store_true",
                    help="node-FACE CONTINUOUS adhesion ON (coh_adh=adh_strength) — the flat-interface "
                         "adhesion that, with surface tension, facets cells into space-filling polyhedra "
                         "(vs sparse node-node point cohesion). Physiological-baseline: should be ON in production.")
    ap.add_argument("--ubottom", action="store_true",
                    help="ULA U-bottom confinement: cells are held in a non-adhesive hemispherical bowl "
                         "(geometric confinement only; the surface never grips). Independent of --well.")
    ap.add_argument("--no-ubottom", action="store_true",
                    help="Override: force the U-bottom rigid bowl OFF even under --ula (free-floating "
                         "hanging-drop-style aggregate held by cohesion alone, no rigid confinement).")
    ap.add_argument("--osmotic", action="store_true",
                    help="osmotic water-flux volume regulation (KB-3.9): per-cell rest volume relaxes "
                         "toward Vc ∝ media-exposed face fraction + concentration feedback → stable volume "
                         "equilibrium (fixes the faceting over-compression strain that stalled the CG).")
    ap.add_argument("--cleave", action="store_true",
                    help="C7 division by SimuCell3D IN-PLACE MESH CLEAVAGE: carve the mother shell into "
                         "mother+daughter along the Hertwig plane (union=mother → zero neighbour "
                         "displacement, no contact spike). Replaces the parked-icosphere mitotic-round "
                         "insert. Pair with --division. Needs a dormant node pool (div pool sized auto).")
    ap.add_argument("--builder", default="fcc", choices=["cubic", "fcc", "voronoi", "sphere"],
                    help="initial cell-centre packing (D10). fcc=close-pack but cuboctahedral-faceted "
                         "ENVELOPE (angular spheroid); sphere=random-close-pack inside a BALL -> SMOOTH "
                         "spherical envelope (hull-psi 0.99 vs fcc 0.96). Use 'sphere' for a round spheroid.")
    ap.add_argument("--init-npz", default=None,
                    help="Restart from a saved aggregate npz (load its final frame as the initial state, "
                         "dropped onto the substrate) instead of building a fresh ball — for 'spread from "
                         "the aggregate' continuation runs (pair with --well, no --ula/--division).")
    ap.add_argument("--filopodia", action="store_true",
                    help="explicit filopodia finger protrusions (node-FACE + node-plane tip adhesions). "
                         "ON automatically under --ula (cell-cell junction formation).")
    ap.add_argument("--lamellipodium", action="store_true",
                    help="per-cell advancing-anchor lamellipodium crawl. ON automatically under --ula.")
    ap.add_argument("--lamel-clutch", action="store_true", dest="lamel_clutch",
                    help="lamellipodial SUBSTRATE clutch: the advancing actin anchor grips the dish (z0) so "
                         "the protrusion transmits traction to the substrate (node-to-plane), not a floating "
                         "bead. The spreading-traction path (pair with --lamellipodium).")
    ap.add_argument("--lamel-all-cell", action="store_true", dest="lamel_all_cell",
                    help="cryptic-FOLLOWER mode: EVERY active cell (not just the substrate rim) grows a "
                         "basal-outward lamellipodial leading edge at the UNCHANGED Gil-Redondo physiological "
                         "per-node force (Farooqui-Fenteany KB-4.6). Fixes the rim-only peeling; spreading "
                         "diagnostic. Forces the host ratchet path (all-cell geometry is host-only).")
    ap.add_argument("--cad-cluster", action="store_true", dest="cad_cluster",
                    help="fine-grained load-sharing cadherin cluster de-cohesion (S1-S3): junction = bundle_n "
                         "parallel trans-dimers, m→0 death; de-cohesion EMERGES as traction ruptures nascent "
                         "junctions → cells crawl out. The spreading-by-dispersal driver.")
    ap.add_argument("--cad-mature", action="store_true", dest="cad_mature",
                    help="cadherin junction maturation (contact-age): matured junctions resist de-cohesion "
                         "(confined = epithelial/MCF-7); nascent let go (dispersing = mesenchymal).")
    ap.add_argument("--cad-n-nascent", type=int, default=4, dest="cad_n_nascent",
                    help="nascent cluster size (KB-4.3 controlled var).")
    ap.add_argument("--active-motility", action="store_true", dest="active_motility",
                    help="per-cell active self-propulsion (persistent random walk) = the active-matter UNJAMMING "
                         "lever: cells crawl in slowly-reorienting directions, fluidising the jammed aggregate so "
                         "it can rearrange / disperse. Pair with de-cohesion (--cad-cluster). Tests if cell MOTILITY "
                         "is the missing spreading piece.")
    ap.add_argument("--f-active-nn", type=float, default=10.0, dest="f_active_nn",
                    help="whole-cell self-propulsion force [nN] (physiological single-cell traction ~1-100 nN, "
                         "KB-2.12; controlled variable, swept — never tuned to a spreading target).")
    ap.add_argument("--motility-tau", type=float, default=600.0, dest="motility_tau",
                    help="polarity persistence time tau_p [s] (cell directional persistence ~10 min).")
    ap.add_argument("--motility-seed", type=int, default=7, dest="motility_seed",
                    help="seed for the random per-cell polarity (ensemble over this for seed-robustness).")
    ap.add_argument("--v0-um-min", type=float, default=0.0, dest="v0_um_min",
                    help="SPV v0-mode: target physiological migration SPEED [µm/min] (~0.5-2). The per-node force is "
                         "set from gamma_node·v0 so cells self-propel at v0, bounding CFL — the SPV-faithful, "
                         "CFL-safe knob (overrides --f-active-nn when >0).")
    ap.add_argument("--bdf2", action="store_true", dest="bdf2",
                    help="TIMESCALE ATTACK (A-stable, PI 2026-07-11): 2nd-order stiffly-stable BDF2 integrator "
                         "instead of backward-Euler. Fixes BE's over-damping of slow active forcing (the freeze that "
                         "broke naive large-dt) fundamentally — so the FORCE-based motility itself stays accurate at "
                         "large dt (no operator-split needed). Same Newton/CG cost per step.")
    ap.add_argument("--t1-rate", action="store_true", dest="t1_rate",
                    help="TIMESCALE ATTACK biology-time route (PI 2026-07-12): supply the slow tissue rearrangement "
                         "the large-dt BDF2 mechanics correctly freeze, as a physical KMC of T1 events (k_T1(s) rate "
                         "law, grounded k0←k_endo/KB-4.13). Pair with --bdf2 + a large --accel-dt. See "
                         "DCM_T1_RATE_COARSEGRAIN_DESIGN_2026-07-12.")
    ap.add_argument("--t1-k0", type=float, default=0.03, dest="t1_k0",
                    help="T1 attempt/gating frequency [1/s] (k_endo anchor 0.01-0.1, KB-4.13).")
    ap.add_argument("--t1-barrier-b", type=float, default=3.0, dest="t1_barrier_b",
                    help="T1 shape-index barrier stiffness B (calibrated from the small-dt k_T1(s) measurement).")
    ap.add_argument("--t1-cadence", type=int, default=200, dest="t1_cadence",
                    help="steps between KMC T1 passes (host round-trip cadence).")
    ap.add_argument("--t1-seed", type=int, default=13, dest="t1_seed")
    ap.add_argument("--motility-split", action="store_true", dest="motility_split",
                    help="TIMESCALE ATTACK: apply the v0-mode drift as a position translation (x*=xₙ+v0·dt·p̂) "
                         "BEFORE the implicit relax (Lie-Trotter operator split) instead of a self-propulsion force. "
                         "The force path freezes at large dt (implicit equilibrates F_active vs contact); the split "
                         "preserves the drift at ANY dt → enables large-dt physiological-time runs. v0-mode only.")
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
    ubottom = (a.ubottom or ula) and not a.no_ubottom        # --no-ubottom overrides the bowl off
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
        lamel_clutch=a.lamel_clutch, lamel_all_cell=a.lamel_all_cell,
        gpu_probe=a.gpu_probe,
        coupling=a.coupling,
        cad_bundle=(20.0 if a.cad_cluster else 40.0), ecm_bundle=167.0,
        cad_cluster=a.cad_cluster, cad_mature=a.cad_mature, cad_n_nascent=a.cad_n_nascent,
        active_motility=a.active_motility, f_active_N=a.f_active_nn * 1e-9,
        motility_persistence_s=a.motility_tau, motility_seed=a.motility_seed,
        motility_v0_um_s=a.v0_um_min / 60.0, motility_split=a.motility_split,
        t1_rate=a.t1_rate, t1_k0=a.t1_k0, t1_barrier_b=a.t1_barrier_b,
        t1_cadence=a.t1_cadence, t1_seed=a.t1_seed,
        remesh_period=a.remesh_period,
        substrate_wetting=substrate_wetting, use_substrate_well=use_substrate_well,
        ubottom=ubottom,                                      # ULA non-adhesive bowl confinement
        division=a.division, div_real_hours=a.div_real_hours, div_t_cycle_h=a.div_t_cycle_h,
        cleave=a.cleave,
        init_npz=a.init_npz,
        osmotic=a.osmotic,
        accel_real_hours=a.accel_real_hours,
        necrosis=a.necrosis,
        builder=a.builder, integrator=("bdf2" if a.bdf2 else "implicit"), accel_dt=a.accel_dt,
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

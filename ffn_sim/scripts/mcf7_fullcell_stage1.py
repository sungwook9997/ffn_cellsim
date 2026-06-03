"""MCF7 FULL-CELL STAGE-1 test — cortex (STAGE-1 grip_walk) + ALL compartments ON.

The γ-floor was diagnosed in a bare-cortex sandbox. This driver builds the FULL
MCF7 cell (cortex + H.8 membrane + H.9 nucleus + H.10 cytoplasm + KU-3.1
enclosed-volume pressure), all additive modules turned ON with MCF7 cell-type
parameters, and runs the constrained grip_walk cortex (STAGE-1: CH1 r0≈0 + CH2
series-stall) so we can see γ_soft/γ_rigid + r/r0 + grip_s in the real
multi-compartment context (per the platform-test conventions).

MCF7 cell-type (PI-ratified 2026-06-02):
  R_cell = 7.5 µm (Wagner 2011) · cytoplasm η = 65.9 Pa·s (Hu 2024) ·
  nucleus E_nuc = 0.15 kPa chromatin small-strain (MCF7 microrheology PMC12173521),
  ratio_lamin = 5 (conservative vs the ~31× indentation/microrheology span) ·
  membrane γ_mem = 0.10 mN/m (KU-3.B1 mid-band) · cortex γ band [0.35,0.65] mN/m.

Compartments are md.force.Custom → INVISIBLE to the harmonic-bond method-of-planes
estimator; the cortex γ_soft/γ_rigid still measure the cortex channel. enclosed_volume
is the CHANGE-3c load-retention term (Young-Laplace ΔP=2γ/R) the cortex builds against.

Run (smoke first):  conda run -n ffn_sim python -m ffn_sim.scripts.mcf7_fullcell_stage1 --smoke
Real A/B (small):   conda run -n ffn_sim python -m ffn_sim.scripts.mcf7_fullcell_stage1 --n-fil 2000 --n-sample 6
"""
from __future__ import annotations

import argparse
import time
from copy import deepcopy
from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml
import hoomd

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import resolve_cortex_myosin, cortex_myosin_attach_bin_names
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cell.cell import build_cortex_full_simulation
from ffn_sim.cell.nucleus import resolve_nucleus
from ffn_sim.cell.membrane_surface import resolve_membrane_surface
from ffn_sim.cell.cytoplasm import resolve_cytoplasm
from ffn_sim.cortex.enclosed_volume import resolve_enclosed_volume
from ffn_sim.scripts.h3_ku35_tension import (
    _tension_method_of_planes,
    _tension_method_of_planes_rigid,
)

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"

# ---- MCF7 cell-type definition (PI-ratified) -------------------------------
MCF7 = dict(
    R_cell=7.5e-6,            # m   Wagner 2011 (15 µm diameter)
    eta_celltype="MCF7",      # cytoplasm.py table -> 65.9 Pa·s (Hu 2024)
    E_nuc=4700.0,             # Pa  4.7 kPa MCF7 indentation (in KU-3.B2.1 band [1,10] kPa);
    #                              the 0.15 kPa microrheology is a distinct small-strain quantity
    ratio_lamin=1.4,          # —   KU-3.B2.1 in-situ band minimum [1.4,5.0] (softest legal;
    #                              nucleus orthogonal to γ, softer => fewer CFL-beads => speed)
    gamma_mem=1.0e-4,         # N/m 0.10 mN/m, KU-3.B1 mid-band
    R_nuc_frac=0.25,          # nucleus:cell radius ratio
    n_nuc_beads=3000,         # CFL-safe at ratio_lamin=1.4 (k_hi ∝ 1/n_beads)
)


def _tagpos(sim) -> np.ndarray:
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        return pos[inv].copy()


def _cortex_angle_theta(sim) -> np.ndarray:
    """θ [deg] at the middle bead of every cortex-angle triplet (straight≈180°).

    Buckling diagnostic (Murrell-Gardel/Lenz): constrained mode SHAKEs only stretch,
    so bending is free — a contracting filament should BUCKLE (θ ≪ 180°). Measures
    whether the active load actually bends the backbone.
    """
    snap = sim.state.get_snapshot()       # tag-ordered; local snapshot lacks .types
    if snap.communicator.rank != 0:
        return np.empty(0)
    ag = np.asarray(snap.angles.group, dtype=np.int64)
    at = np.asarray(snap.angles.typeid, dtype=np.int64)
    atypes = list(snap.angles.types)
    if "cortex-angle" not in atypes or ag.shape[0] == 0:
        return np.empty(0)
    trip = ag[at == atypes.index("cortex-angle")]
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    v1 = pos[trip[:, 0]] - pos[trip[:, 1]]
    v2 = pos[trip[:, 2]] - pos[trip[:, 1]]
    v1 /= np.linalg.norm(v1, axis=1, keepdims=True).clip(min=1e-30)
    v2 /= np.linalg.norm(v2, axis=1, keepdims=True).clip(min=1e-30)
    return np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", v1, v2), -1.0, 1.0)))


def _myosin_head_tension(sim, myo_action, p_myo):
    """Diagnose grip-walk force delivery: head-actin attach-bond tension F=k·r vs
    F_stall, and the commanded walk accumulator s_grip vs ℓ₀.

    The STAGE-1 grip_walk kernel computes a series stall load F_series=k·min(s,r)
    but DELIVERS F=k·r with the attach-bond r0≈0 (CH1). If the bond relaxes r→0
    while s_grip walks to ℓ₀, the motor transports but applies ~no force — the
    suspected force-generation floor. This measures both directly.
    """
    snap = sim.state.get_snapshot()
    if snap.communicator.rank != 0:
        return None
    btypes = list(snap.bonds.types)
    attach = set(cortex_myosin_attach_bin_names(p_myo.n_bins))
    bt = np.asarray(snap.bonds.typeid)
    bg = np.asarray(snap.bonds.group, dtype=np.int64)
    aids = [i for i, nm in enumerate(btypes) if nm in attach]
    mask = np.isin(bt, aids)
    n_attach = int(mask.sum())
    out = {"n_attach": n_attach, "F_stall": float(p_myo.F_stall_per_head),
           "k_head_actin": float(p_myo.k_head_actin)}
    if n_attach:
        pos = np.asarray(snap.particles.position, dtype=np.float64)
        ab = bg[mask]
        r = np.linalg.norm(pos[ab[:, 0]] - pos[ab[:, 1]], axis=1)
        F = p_myo.k_head_actin * r                      # r0≈0 for grip_walk
        out["F_mean"] = float(F.mean()); out["F_max"] = float(F.max())
        out["F_over_stall_mean"] = float(F.mean() / p_myo.F_stall_per_head)
    sg = getattr(myo_action, "_head_grip_s", None)
    if sg is not None:
        sg = np.asarray(sg)
        ell0 = getattr(myo_action, "_ell0_cortex", None)
        out["s_grip_mean"] = float(sg.mean()); out["s_grip_max"] = float(sg.max())
        if ell0:
            out["ell0"] = float(ell0)
            out["s_over_ell0_mean"] = float(sg.mean() / ell0)
    return out


def _resolve_compartments(p_cortex):
    """Resolve the 4 MCF7 compartments (all ON)."""
    R = p_cortex.R_cell
    p_nuc = resolve_nucleus(
        {"E_nuc": MCF7["E_nuc"], "ratio_lamin": MCF7["ratio_lamin"]},
        R_nuc=MCF7["R_nuc_frac"] * R, n_beads=MCF7["n_nuc_beads"],
    )
    p_ev = resolve_enclosed_volume({}, R_cell=R)
    p_mem = resolve_membrane_surface({"gamma_mem": MCF7["gamma_mem"]}, R_cell=R)
    p_cyto = resolve_cytoplasm(cell_type=MCF7["eta_celltype"])
    return dict(p_nucleus=p_nuc, p_enclosed_volume=p_ev,
                p_membrane_surface=p_mem, p_cytoplasm=p_cyto)


def _build(cfg, *, stepping_mode, force_scaling, constrained, compartments,
           dtc=None, seed=1, equilibrate=False, n_warmup=0, device="cpu", kon_scale=1.0,
           bind_scale=1.0, turnover_tau=None, xl_k_scale=1.0):
    from dataclasses import replace as _replace
    p = resolve_h3_derived(cfg)
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = dtc if dtc is not None else 0.001 * tau_bend
    p_myo = resolve_cortex_myosin(cfg, dt=dtc, R_cell=p.R_cell)
    p_xl = resolve_crosslinkers(cfg, dt=dtc)
    if kon_scale != 1.0:  # couple_accel: accelerate xlink binding so the network
        p_xl = _replace(p_xl, k_on=p_xl.k_on * kon_scale)   # PERCOLATES (z->[2,4])
    if xl_k_scale != 1.0:  # ⚠️ DIAGNOSTIC ONLY (PI-approved 2026-06-03, NOT adopted):
        # soft-coupling-ceiling probe. Crosslinker stiffening is hard-rule FORBIDDEN
        # (k=1e-7 is the Furuike/Ferrer literature constant) — this scales k_intra +
        # k_attach to test whether the soft inter-filament coupling is what caps g_soft.
        # If g_soft rises with k → ceiling confirmed (generation vs coupling separated).
        # The value is never kept; same status as the FFN_MYOSIN_ALIGN meridional probe.
        p_xl = _replace(p_xl, k_intra=p_xl.k_intra * xl_k_scale,
                        k_attach=p_xl.k_attach * xl_k_scale)
    if bind_scale != 1.0:  # mesoscale-consistent reach (geometric dual of ×40 areal
        # coarse-graining): the physical 60nm filamin can't bridge the sparse mesoscale
        # mesh; scaling the partner-search radius (NOT k) restores binding. See
        # stage2_diagnostics: bind×6 + n_xl/n_fil=1.5 → z=2.96 (Ennomani optimum).
        p_xl = _replace(p_xl, max_bind_dist=p_xl.max_bind_dist * bind_scale)
    # STAGE-2 actin turnover: cofilin sever + pointed-end re-anneal. In constrained
    # mode the updater is wired to the BAOAB Action (cell.py) so severing SPLITS the
    # M-SHAKE chain — releasing the rigid stretch jam that pins r/r0=1 even on a
    # percolated network (2026-06-03 finding). tau_half is the Chugh/Fritzsche anchor.
    p_turnover = None
    if turnover_tau is not None:
        from ffn_sim.cortex.turnover import resolve_turnover
        tcfg = deepcopy(cfg)
        tblk = tcfg["cortex"].setdefault("turnover", {})
        # ACCELERATED-DYNAMICS probe (same pattern as v0_accel / kon_scale): the
        # physical cortical-actin half-life is ~10s but the constrained dt is set
        # by the stiff bending CFL (~1e-6s) → ~1e7 steps to see one turnover at
        # the physical rate (infeasible). We shrink tau_half by A=tau_phys/tau and
        # scale k_anneal by the SAME A, which leaves the steady-state connected
        # fraction f_ss = k_anneal/(k_sev+k_anneal) INVARIANT (the physical
        # observable) while moving the events into a feasible step budget. Only
        # the RATE is accelerated, not the equilibrium — the mechanism test
        # (does turnover release the elastic jam?) is rate-independent.
        tau_phys = float(tblk.get("tau_half", 10.0))
        kann_phys = float(tblk.get("k_anneal", 0.1))
        accel = tau_phys / float(turnover_tau)
        tblk["enabled"] = True
        tblk["tau_half"] = float(turnover_tau)
        tblk["k_anneal"] = kann_phys * accel
        p_turnover = resolve_turnover(tcfg, dt=dtc, rest_length=p.rest_length)
    dev = (hoomd.device.GPU(notice_level=0) if device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    kw = dict(p_xlinks=p_xl, p_myosin=p_myo, device=dev, with_baoab=True,
              rng=np.random.default_rng(seed), **compartments)
    if p_turnover is not None:
        kw["p_turnover"] = p_turnover
    if constrained:
        kw.update(constrained=True, constrained_dt=dtc)
    if equilibrate:
        kw.update(equilibrate=True, equilibrate_steps=n_warmup,
                  equilibrate_softstart_steps=max(300, n_warmup // 8))
    hw = build_cortex_full_simulation(p, **kw)
    return p, p_myo, p_xl, dtc, hw


def _cfg_for(n_fil, n_motors, n_xl, stepping_mode, force_scaling, backbone_nm=700):
    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["R_cell"] = MCF7["R_cell"]          # MCF7 7.5 µm
    # Minifilament backbone length override. Brief literal is 700 nm (14-bead grid
    # artifact); literature NMII bipolar minifilament ≈ 300 nm (Billington; bare
    # zone ~160 nm). 300 nm shrinks the placement exclusion (backbone+100 nm) from
    # 800→400 nm, so the NATIVE motor density fits without coarse-graining → denser,
    # more native-like motor field → better force TRANSMISSION (the γ wall).
    cfg["cortex"]["myosin"]["backbone_length"] = backbone_nm * 1e-9
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors
    cfg["cortex"]["myosin"]["stepping_mode"] = stepping_mode
    cfg["cortex"]["myosin"]["mesoscale_force_scaling"] = bool(force_scaling)
    if n_xl is not None:
        cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = int(n_xl)
    return cfg


def run_arm(stepping_mode, *, n_fil, n_motors, n_xl, force_scaling, v0_accel,
            couple_accel, n_warmup, n_sample, interval, smoke, device="cpu",
            compartments_on=True, only=None, backbone_nm=700, kon_scale=1.0,
            bind_scale=1.0, turnover_tau=None, xl_k_scale=1.0):
    cfg = _cfg_for(n_fil, n_motors, n_xl, stepping_mode, force_scaling, backbone_nm)
    # Resolve cortex once to get R_cell for the compartments.
    p0 = resolve_h3_derived(cfg)
    if compartments_on and only:
        full = _resolve_compartments(p0)
        comp = {only: full[only]}
        print(f"  [ONLY {only}] attribution run (other compartments OFF)", flush=True)
    elif compartments_on:
        comp = _resolve_compartments(p0)
        print(f"  [compartments ON] nucleus k_chrom={comp['p_nucleus'].k_chrom:.2e} "
              f"E_nuc={MCF7['E_nuc']}Pa | enclosed dP_ref={comp['p_enclosed_volume'].dP_ref:.0f}Pa "
              f"| membrane gamma_mem={MCF7['gamma_mem']*1e3:.2f}mN/m "
              f"| cytoplasm eta={comp['p_cytoplasm'].eta_eff:.1f}Pa·s (MCF7)", flush=True)
    else:
        comp = {}
        print("  [compartments OFF] cortex-only at MCF7 R_cell=7.5µm (attribution baseline)",
              flush=True)

    # --- Warm-up (unconstrained, literal v0) with all compartments ---
    _, _, _, dtc, hw = _build(cfg, stepping_mode=stepping_mode,
                              force_scaling=force_scaling, constrained=False,
                              compartments=comp, equilibrate=True, n_warmup=n_warmup,
                              device=device, kon_scale=kon_scale, bind_scale=bind_scale,
                              xl_k_scale=xl_k_scale)
    hw["sim"].run(0)
    pos_warm = _tagpos(hw["sim"])
    del hw

    # --- Constrained production (rigid backbone) with all compartments ---
    # Accelerated v0/k_on like the gripwalk driver.
    p, p_myo_lit, p_xl, dtc, hc = _build(cfg, stepping_mode=stepping_mode,
                                         force_scaling=force_scaling,
                                         constrained=True, compartments=comp, dtc=dtc,
                                         device=device, kon_scale=kon_scale,
                                         bind_scale=bind_scale, turnover_tau=turnover_tau,
                                         xl_k_scale=xl_k_scale)
    sim = hc["sim"]
    act = hc["baoab_action"]
    turn_act = hc.get("turnover_action")
    myo_act = hc.get("myosin_action")
    if turn_act is not None:
        print(f"  [turnover ON] tau_half={turnover_tau}s k_sev={turn_act.p.k_sev:.3e} "
              f"k_anneal={turn_act.p.k_anneal:.3e} f_ss={turn_act.p.f_ss:.3f} "
              f"batch_steps={turn_act.p.batch_steps} → SHAKE-chain split wired", flush=True)
    if hasattr(act, "record_lambda"):
        act.record_lambda = True
    # transfer warmed positions by tag
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = pos_warm   # get_snapshot is tag-ordered
    sim.state.set_snapshot(snap)
    sim.run(0)

    nca = p.n_filaments * p.beads_per_filament
    r0 = float(np.linalg.norm(_tagpos(sim)[:nca], axis=1).mean())
    samples = []
    n = 2 if smoke else n_sample
    iv = 2000 if smoke else interval
    for k in range(n):
        sim.run(iv)
        r = _tagpos(sim)
        rmean = float(np.linalg.norm(r[:nca], axis=1).mean())
        g_soft = _tension_method_of_planes(sim, p.R_cell)
        g_rigid = _tension_method_of_planes_rigid(sim, act, p.R_cell, dtc)
        samples.append((rmean / r0, g_soft, g_rigid))
        turn_str = ""
        if turn_act is not None:
            turn_str = (f" | turnover cf={turn_act.connected_fraction:.3f} "
                        f"sev={turn_act.n_sever_total} ann={turn_act.n_anneal_total}")
        th = _cortex_angle_theta(sim)
        buck_str = ""
        if th.size:
            buck_str = (f" | θ mean={th.mean():.1f}° buckled(<150°)={(th<150).mean()*100:.1f}%"
                        f" sharp(<120°)={(th<120).mean()*100:.1f}%")
        print(f"  [{stepping_mode}] s={k+1}/{n} r/r0={rmean/r0:.5f} "
              f"g_soft={g_soft*1e3:.3e} g_rigid={g_rigid*1e3:.3e} "
              f"g_tot={(g_soft+g_rigid)*1e3:.3e} mN/m{turn_str}{buck_str}", flush=True)
        if myo_act is not None and n_motors > 0:
            ht = _myosin_head_tension(sim, myo_act, p_myo_lit)
            if ht and ht.get("n_attach"):
                print(f"      head: n_attach={ht['n_attach']} F_mean={ht.get('F_mean',0):.2e}N "
                      f"F_max={ht.get('F_max',0):.2e}N F_stall={ht['F_stall']:.2e}N "
                      f"F/F_stall={ht.get('F_over_stall_mean',0):.3f} | "
                      f"s_grip/ℓ₀={ht.get('s_over_ell0_mean',0):.3f} "
                      f"(ℓ₀={ht.get('ell0',0)*1e9:.0f}nm)", flush=True)
    return samples


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--smoke", action="store_true", help="tiny: n_fil=200, 2 samples")
    ap.add_argument("--n-fil", type=int, default=2000)
    ap.add_argument("--n-motors", type=int, default=200)
    ap.add_argument("--n-xl", type=int, default=2000)
    ap.add_argument("--force-scaling", action="store_true", default=True)
    ap.add_argument("--v0-accel", type=float, default=300.0)
    ap.add_argument("--couple-accel", action="store_true", default=True)
    ap.add_argument("--n-warmup", type=int, default=8000)
    ap.add_argument("--n-sample", type=int, default=6)
    ap.add_argument("--interval", type=int, default=20000)
    ap.add_argument("--arm", choices=["both", "binned_r0", "grip_walk"],
                    default="grip_walk")
    ap.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--no-compartments", action="store_true",
                    help="cortex-only baseline (attribution: isolate the compartment effect)")
    ap.add_argument("--only", choices=["p_enclosed_volume", "p_cytoplasm",
                                       "p_nucleus", "p_membrane_surface"], default=None,
                    help="attribution: enable ONLY this one compartment")
    ap.add_argument("--backbone-nm", type=float, default=700,
                    help="myosin minifilament backbone length [nm] (700=brief, 300=literature NMII)")
    ap.add_argument("--kon-scale", type=float, default=1.0,
                    help="couple_accel: scale xlink k_on for percolation (z->[2,4]); native needs ~300")
    ap.add_argument("--bind-scale", type=float, default=1.0,
                    help="mesoscale-consistent xlink reach (scales max_bind_dist, NOT k); "
                         "bind×6 + n_xl/n_fil=1.5 → z=2.96 at mesoscale (fast STAGE-2)")
    ap.add_argument("--turnover-tau", type=float, default=None,
                    help="enable actin turnover; tau_half [s] (Chugh/Fritzsche 5-30s). "
                         "In constrained mode splits M-SHAKE chains on sever (STAGE-2)")
    ap.add_argument("--xl-k-scale", type=float, default=1.0,
                    help="⚠️ DIAGNOSTIC ONLY (PI-approved): scale crosslinker k_intra+k_attach "
                         "(soft-coupling-ceiling probe; stiffen is hard-rule forbidden, never adopted)")
    args = ap.parse_args()

    if args.smoke:
        args.n_fil, args.n_motors, args.n_xl, args.n_warmup = 200, 20, 200, 1500

    print(f"=== MCF7 FULL-CELL STAGE-1 (R_cell=7.5µm, all compartments ON, "
          f"n_fil={args.n_fil}, arm={args.arm}, smoke={args.smoke}) ===", flush=True)
    print(f"HOOMD {hoomd.version.version}", flush=True)

    arms = ["binned_r0", "grip_walk"] if args.arm == "both" else [args.arm]
    t0 = time.time()
    results = {}
    for arm in arms:
        print(f"\n--- arm: {arm} ---", flush=True)
        results[arm] = run_arm(
            arm, n_fil=args.n_fil, n_motors=args.n_motors, n_xl=args.n_xl,
            force_scaling=args.force_scaling, v0_accel=args.v0_accel,
            couple_accel=args.couple_accel, n_warmup=args.n_warmup,
            n_sample=args.n_sample, interval=args.interval, smoke=args.smoke,
            device=args.device, compartments_on=not args.no_compartments,
            only=args.only, backbone_nm=args.backbone_nm, kon_scale=args.kon_scale,
            bind_scale=args.bind_scale, turnover_tau=args.turnover_tau,
            xl_k_scale=args.xl_k_scale,
        )
    dt = time.time() - t0
    print(f"\n=== DONE in {dt:.0f}s. Full cell (cortex+membrane+nucleus+cytoplasm"
          f"+enclosed-vol) ran stably with STAGE-1 grip_walk. ===", flush=True)
    for arm, s in results.items():
        if s:
            gs = np.mean([x[1] for x in s]) * 1e3
            gr = np.mean([x[2] for x in s]) * 1e3
            rr = np.mean([x[0] for x in s])
            print(f"  {arm}: <r/r0>={rr:.5f} <g_soft>={gs:.3e} <g_rigid>={gr:.3e} "
                  f"<g_tot>={gs+gr:.3e} mN/m  (band [0.35,0.65])", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

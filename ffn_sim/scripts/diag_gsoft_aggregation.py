#!/usr/bin/env python
"""DIAGNOSTIC — decompose g_soft (KU-3.5 active cortical tension) into
multiplicative factors for the MYOSIN-ATTACH channel and find the cap.

NOT a production fix. Read-only physics decomposition. The established fact
is that the KU-3.5 ACTIVE channel g_soft floors ~1150x under band even though
(a) adding motors does not raise g_soft (force-budget sweep) and (b) per-head
force IS delivered sub-stall. Open question: why don't the delivered per-head
forces AGGREGATE into the measured shell tension g_soft?

This script answers it by decomposing g_soft_attach into:
  N_engaged  · mean per-head F · N_crossing/N_engaged (geometric plane-fraction)
into the method-of-planes plane-sum, and reconciling that against the actual
estimator output run on the attach channel ONLY.

Two channels are instrumented side by side on the SAME isotropic planes the
estimator uses:
  * ATTACH  = harmonic bonds whose type name starts cortex_myosin_attach_b
              (the dynamic head->actin grip-walk bonds; SHORT, head_off~200nm).
  * BACKBONE = the cortex actin backbone bonds (type 'cortex-bond'); the long
               load-bearing network that carries g_rigid.

Recruitment-lever sweeps (discriminator):
  * n_motors sweep        — count knob.
  * capture_perp sweep    — myosin binding RANGE knob (head_actin_capture_perp +
                            head_actin_max_bind_dist). Wider range => more engaged
                            attach bonds at fixed motor count.

Config: matches stage2_forcebudget_sweep (n_fil=300, grip_walk, force_scaling,
backbone_nm=300, bind_scale=6, kon_scale=300) unless overridden. Use --fast for
a quick n_fil=60 tier1 smoke.

Run (worktree, no cd):
  PYTHONPATH=/Users/sw1/ffn_aggdiag \
    /Users/sw1/miniconda3/envs/ffn_sim/bin/python -P \
    /Users/sw1/ffn_aggdiag/ffn_sim/scripts/diag_gsoft_aggregation.py --n-fil 300
"""
from __future__ import annotations

import argparse
import json
from dataclasses import replace as _replace
from pathlib import Path

import numpy as np

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import (
    resolve_cortex_myosin,
    cortex_myosin_attach_bin_names,
)
from ffn_sim.scripts.mcf7_fullcell_stage1 import (
    _build, _cfg_for, _resolve_compartments, _tagpos,
)

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "h3" / "production"

# Fibonacci-isotropic plane normals — IDENTICAL convention to the estimator
# (_tension_method_of_planes in h3_ku35_tension.py).
def _plane_normals(n_planes: int = 12) -> np.ndarray:
    phi = (1 + 5 ** 0.5) / 2
    i = np.arange(n_planes, dtype=np.float64)
    z = 1 - 2 * (i + 0.5) / n_planes
    rxy = np.sqrt(1 - z * z)
    theta = 2 * np.pi * i / phi
    return np.stack([rxy * np.cos(theta), rxy * np.sin(theta), z], axis=1)


def _channel_decomp(sim, R_cell, type_prefixes, n_planes=12):
    """Per-channel method-of-planes decomposition.

    Replays the ESTIMATOR's exact arithmetic but restricted to bonds whose
    type name starts with any of ``type_prefixes`` (e.g. 'cortex_myosin_attach_b'
    for the attach channel, 'cortex-bond' for the backbone). Returns the factor
    decomposition needed to reconcile g_soft for this channel.
    """
    btypes = list(sim.state.bond_types)
    # read the integrator's harmonic bond force params (k, r0 per type)
    bond_force = None
    for f in sim.operations.integrator.forces:
        if hasattr(f, "params") and any(t in btypes for t in getattr(f, "params", {})):
            bond_force = f
            break
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        pos_byTag = pos[inv]
        bg = np.asarray(s.bonds.group)
        bt = np.asarray(s.bonds.typeid)
    if bond_force is None or bg.shape[0] == 0:
        return None

    # channel mask over bonds, by type-name prefix
    sel_typeids = [i for i, t in enumerate(btypes)
                   if any(t.startswith(pref) for pref in type_prefixes)]
    if not sel_typeids:
        return dict(n_bonds=0, n_engaged=0, mean_F=0.0, mean_F_over_stall=0.0,
                    mean_N_crossing=0.0, geo_frac=0.0, ideal_budget_mN_m=0.0,
                    g_channel_mN_m=0.0, obs_over_ideal=float("nan"))
    chan_mask = np.isin(bt, sel_typeids)

    # per-bond k, r0 in the selected channel
    k_arr = np.zeros(bg.shape[0])
    r0_arr = np.zeros(bg.shape[0])
    for i, tname in enumerate(btypes):
        try:
            kp = bond_force.params[tname]
            m = (bt == i)
            k_arr[m] = float(kp["k"])
            r0_arr[m] = float(kp["r0"])
        except Exception:
            pass

    rA = pos_byTag[bg[:, 0]]
    rB = pos_byTag[bg[:, 1]]
    d = rB - rA
    L = np.linalg.norm(d, axis=1)
    L_safe = np.where(L > 0, L, 1.0)
    u = d / L_safe[:, None]
    T = k_arr * (L - r0_arr)          # signed scalar tension (N) per bond

    cm = chan_mask
    n_bonds_chan = int(cm.sum())
    # ENGAGED = channel bonds that actually carry load (extension > 0 within
    # this channel). For attach (grip-walk r0~0) essentially all are engaged
    # springs; report the count and the per-bond |F| = |T| over the channel.
    F_chan = np.abs(T[cm])
    mean_F = float(np.mean(F_chan)) if n_bonds_chan else 0.0

    normals = _plane_normals(n_planes)
    g_chan_list = []
    ncross_list = []
    planesum_list = []
    for n_hat in normals:
        a_side = rA @ n_hat
        b_side = rB @ n_hat
        crossing = ((a_side * b_side) < 0) & cm
        ncross_list.append(int(crossing.sum()))
        if not crossing.any():
            g_chan_list.append(0.0)
            planesum_list.append(0.0)
            continue
        f_cut = T[crossing] * np.abs(u[crossing] @ n_hat)
        planesum = float(np.sum(f_cut))
        planesum_list.append(planesum)
        g_chan_list.append(planesum / (2.0 * np.pi * R_cell))
    g_chan = float(np.mean(np.abs(np.asarray(g_chan_list))))     # mN/m via *1e3 later
    mean_ncross = float(np.mean(ncross_list))
    geo_frac = mean_ncross / n_bonds_chan if n_bonds_chan else 0.0
    # IDEAL budget: every channel bond crosses and projects fully.
    # g_ideal = (N_engaged * mean_F) / (2 pi R)
    ideal_budget = (n_bonds_chan * mean_F) / (2.0 * np.pi * R_cell)
    return dict(
        n_bonds=n_bonds_chan,
        mean_F=mean_F,
        mean_planesum=float(np.mean(np.abs(planesum_list))),
        mean_N_crossing=mean_ncross,
        geo_frac=geo_frac,
        ideal_budget_mN_m=ideal_budget * 1e3,
        g_channel_mN_m=g_chan * 1e3,
        obs_over_ideal=(g_chan / ideal_budget) if ideal_budget > 0 else float("nan"),
    )


def _bound_F_over_stall(ma, sim):
    """Mean F/F_stall on bound heads (mirrors stage2_forcebudget_sweep
    _bound_metrics). The per-head force is k_head_actin * |head - bound bead|."""
    p = ma.p
    bound = ma._head_bound_to_actin >= 0
    n_bound = int(bound.sum())
    n_heads = int(ma._head_bound_to_actin.size)
    if n_bound == 0:
        return 0.0, 0.0, 0, n_heads
    pos = _tagpos(sim)
    bidx = np.flatnonzero(bound)
    htags = np.array([ma._head_global_tag(int(h)) for h in bidx], dtype=np.int64)
    atags = ma._head_bound_to_actin[bidx]
    r = np.linalg.norm(pos[htags] - pos[atags], axis=1)
    F = p.k_head_actin * np.clip(r, 0.0, None)
    return (float(np.mean(F)),
            float(np.mean(F / p.F_stall_per_head)),
            n_bound, n_heads)


def _run_condition(n_fil, n_motors, n_xl, *, n_warmup, n_sample, interval,
                   bind_scale, kon_scale, capture_perp_scale=1.0,
                   device="cpu", seed=1, n_planes=12, backbone_nm=300):
    cfg = _cfg_for(n_fil, n_motors, n_xl, "grip_walk", force_scaling=True,
                   backbone_nm=backbone_nm)
    # binding-RANGE recruitment knob for the MYOSIN attach channel.
    if capture_perp_scale != 1.0:
        myo = cfg["cortex"]["myosin"]
        myo["head_actin_capture_perp"] = (
            float(myo["head_actin_capture_perp"]) * capture_perp_scale)
        myo["head_actin_max_bind_dist"] = (
            float(myo["head_actin_max_bind_dist"]) * capture_perp_scale)
    p0 = resolve_h3_derived(cfg)
    comp = _resolve_compartments(p0)

    # warm-up (unconstrained) then constrained production — same as the sweep.
    _, _, _, dtc, hw = _build(cfg, stepping_mode="grip_walk", force_scaling=True,
                              constrained=False, compartments=comp, equilibrate=True,
                              n_warmup=n_warmup, device=device, kon_scale=kon_scale,
                              bind_scale=bind_scale, seed=seed)
    hw["sim"].run(0)
    pos_warm = _tagpos(hw["sim"])
    del hw
    p, _pm, _px, dtc, hc = _build(cfg, stepping_mode="grip_walk", force_scaling=True,
                                  constrained=True, compartments=comp, dtc=dtc,
                                  device=device, kon_scale=kon_scale,
                                  bind_scale=bind_scale, seed=seed)
    sim = hc["sim"]
    ma = hc.get("myosin_action")
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = pos_warm
    sim.state.set_snapshot(snap)
    sim.run(0)

    R = p.R_cell
    samples = []
    for k in range(n_sample):
        sim.run(interval)
        attach = _channel_decomp(sim, R, ["cortex_myosin_attach_b"], n_planes)
        backbone = _channel_decomp(sim, R, ["cortex-bond"], n_planes)
        if ma is not None and n_motors > 0:
            mean_F, mean_fos, n_bound, n_heads = _bound_F_over_stall(ma, sim)
            n_engaged = int(ma.n_engaged)
        else:
            mean_F = mean_fos = 0.0
            n_bound = n_engaged = 0
        rec = dict(
            sample=k + 1,
            n_engaged=n_engaged,
            n_bound_heads=n_bound,
            mean_perhead_F_pN=mean_F * 1e12,
            mean_F_over_stall=mean_fos,
            attach_n_bonds=attach["n_bonds"] if attach else 0,
            attach_mean_F_pN=(attach["mean_F"] * 1e12) if attach else 0.0,
            attach_N_crossing=attach["mean_N_crossing"] if attach else 0.0,
            attach_geo_frac=attach["geo_frac"] if attach else 0.0,
            attach_ideal_budget_mN_m=attach["ideal_budget_mN_m"] if attach else 0.0,
            attach_g_mN_m=attach["g_channel_mN_m"] if attach else 0.0,
            attach_obs_over_ideal=attach["obs_over_ideal"] if attach else float("nan"),
            backbone_n_bonds=backbone["n_bonds"] if backbone else 0,
            backbone_N_crossing=backbone["mean_N_crossing"] if backbone else 0.0,
            backbone_geo_frac=backbone["geo_frac"] if backbone else 0.0,
            backbone_g_mN_m=backbone["g_channel_mN_m"] if backbone else 0.0,
        )
        samples.append(rec)
        print(f"    s={k+1}/{n_sample} n_eng={n_engaged:4d} "
              f"F/Fstall={mean_fos:.3f} attachN={rec['attach_n_bonds']:4d} "
              f"Ncross={rec['attach_N_crossing']:.2f} "
              f"geo={rec['attach_geo_frac']:.4f} "
              f"g_attach={rec['attach_g_mN_m']:.3e} "
              f"ideal={rec['attach_ideal_budget_mN_m']:.3e} "
              f"obs/ideal={rec['attach_obs_over_ideal']:.3e} | "
              f"bbN={rec['backbone_n_bonds']} bbGeo={rec['backbone_geo_frac']:.4f}",
              flush=True)

    def _m(key):
        vals = [s[key] for s in samples if np.isfinite(s[key])]
        return float(np.mean(vals)) if vals else float("nan")

    agg = {k: _m(k) for k in samples[0] if k != "sample"}
    agg["n_motors"] = n_motors
    agg["capture_perp_scale"] = capture_perp_scale
    agg["samples"] = samples
    return agg


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=300)
    ap.add_argument("--motors", type=int, nargs="+", default=[50, 100, 200, 400])
    ap.add_argument("--capture-scales", type=float, nargs="+", default=[0.5, 1.0, 2.0])
    ap.add_argument("--n-samples", type=int, default=4)
    ap.add_argument("--interval", type=int, default=6000)
    ap.add_argument("--n-warmup", type=int, default=6000)
    ap.add_argument("--bind-scale", type=float, default=6.0)
    ap.add_argument("--kon-scale", type=float, default=300.0)
    ap.add_argument("--n-planes", type=int, default=12)
    ap.add_argument("--backbone-nm", type=int, default=300)
    ap.add_argument("--device", default="cpu", choices=["cpu", "gpu"])
    ap.add_argument("--fast", action="store_true",
                    help="n_fil=60 tier1 smoke (overrides --n-fil, fewer motors)")
    ap.add_argument("--ref-motors", type=int, default=100,
                    help="fixed motor count for the capture-perp range sweep")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    if args.fast:
        args.n_fil = 60
        args.motors = [50, 200]
        args.capture_scales = [1.0, 2.0]
        args.n_samples = 3
    n_xl = int(1.5 * args.n_fil)

    print(f"[gsoft-aggregation diag] n_fil={args.n_fil} n_xl={n_xl} "
          f"motors={args.motors} capture_scales={args.capture_scales} "
          f"n_samples={args.n_samples} interval={args.interval} "
          f"n_planes={args.n_planes} backbone_nm={args.backbone_nm}\n", flush=True)

    # ---- Sweep 1: n_motors (count recruitment lever) at capture_scale=1 ----
    print("=== SWEEP 1: n_motors (count lever), capture_perp_scale=1.0 ===", flush=True)
    motor_sweep = []
    for nm in args.motors:
        print(f"  [n_motors={nm}]", flush=True)
        r = _run_condition(args.n_fil, nm, n_xl, n_warmup=args.n_warmup,
                           n_sample=args.n_samples, interval=args.interval,
                           bind_scale=args.bind_scale, kon_scale=args.kon_scale,
                           capture_perp_scale=1.0, device=args.device,
                           n_planes=args.n_planes, backbone_nm=args.backbone_nm)
        motor_sweep.append(r)
        print(f"  >> n_motors={nm}: n_eng={r['n_engaged']:.0f} "
              f"g_attach={r['attach_g_mN_m']:.3e} mN/m "
              f"F/Fstall={r['mean_F_over_stall']:.3f} "
              f"geo_frac={r['attach_geo_frac']:.4f}\n", flush=True)

    # ---- Sweep 2: capture_perp range lever at fixed ref-motors ----
    print(f"=== SWEEP 2: capture_perp range lever, n_motors={args.ref_motors} ===",
          flush=True)
    range_sweep = []
    for cs in args.capture_scales:
        print(f"  [capture_perp_scale={cs}]", flush=True)
        r = _run_condition(args.n_fil, args.ref_motors, n_xl, n_warmup=args.n_warmup,
                           n_sample=args.n_samples, interval=args.interval,
                           bind_scale=args.bind_scale, kon_scale=args.kon_scale,
                           capture_perp_scale=cs, device=args.device,
                           n_planes=args.n_planes, backbone_nm=args.backbone_nm)
        range_sweep.append(r)
        print(f"  >> capture_scale={cs}: n_eng={r['n_engaged']:.0f} "
              f"g_attach={r['attach_g_mN_m']:.3e} mN/m "
              f"F/Fstall={r['mean_F_over_stall']:.3f} "
              f"geo_frac={r['attach_geo_frac']:.4f}\n", flush=True)

    res = dict(
        config=dict(n_fil=args.n_fil, n_xl=n_xl, motors=args.motors,
                    capture_scales=args.capture_scales, ref_motors=args.ref_motors,
                    n_samples=args.n_samples, interval=args.interval,
                    n_warmup=args.n_warmup, bind_scale=args.bind_scale,
                    kon_scale=args.kon_scale, n_planes=args.n_planes,
                    backbone_nm=args.backbone_nm),
        motor_sweep=motor_sweep,
        range_sweep=range_sweep,
    )
    _OUT.mkdir(parents=True, exist_ok=True)
    out = Path(args.out) if args.out else (_OUT / "diag_gsoft_aggregation.json")
    out.write_text(json.dumps(res, indent=2, default=float))
    print(f"\n[json] {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

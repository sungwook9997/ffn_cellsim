"""STAGE-2 / pre-fix gating measurement — SIGNED cortical contractility (contraction vs extension).

The gated KU-3.5 estimator (`h3_ku35_tension._tension_method_of_planes`) returns
``mean(|γ_plane|)`` — a tension-scalar MAGNITUDE. By construction it **cannot tell net
CONTRACTION (tensile hoop stress, the cell pulling itself together) from net EXTENSION
(compressive, pushing apart) or from locally-cancelling dipoles** (the diagnosed root cause:
each minifilament is a local force dipole whose two ends cancel → net ~0 even when |stress| is
nonzero). Stam 2017 proves the SAME rigid actomyosin network is EXTENSILE at low cross-link
density and contractile only above a threshold — a magnitude-only readout is blind to which
regime the cortex is in.

This module measures the SIGNED method-of-planes stress — the SAME f_cut = T·|û·n̂| per
crossing bond (T = k(L−r0), tension +, compression −; the |û·n̂| abs is the IK-validated
storage-order fix, retained), summed per plane / 2πR, but reported as **mean(γ_plane) keeping
the sign** alongside the gated mean(|γ_plane|). It does NOT modify the gated estimator (no
gate-loosening). Applied motors ON vs OFF on the SAME network:

    Δ_signed = signed(ON) − signed(OFF)
      > 0  → myosin adds NET TENSILE (contractile) hoop stress   [what we want]
      ≈ 0  → local dipoles cancel (force non-aggregation — the diagnosed wall)
      < 0  → myosin is EXTENSILE (pushing apart — Stam low-crosslink rigid regime)

This is the gating diagnostic the critic flagged as prerequisite #1: without a SIGNED
observable, no contraction-fix can be validated (a raw g_soft jump can be a structural relabel
or even an extensile stress reading as positive magnitude).

Run:  conda activate ffn_sim
      python -m ffn_sim.scripts.stage2_signed_contractility --smoke
      python -m ffn_sim.scripts.stage2_signed_contractility --n-fil 600 --n-samples 6
Outputs: ffn_sim/outputs/h3/production/stage2_signed_contractility.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import hoomd

from ffn_sim.common.production_policy import (
    add_production_device_args,
    require_full_cell_physiological_baseline,
    validate_production_device_args,
)
from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.scripts.mcf7_fullcell_stage1 import (
    _build,
    _cfg_for,
    _resolve_compartments,
    _tagpos,
)

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "h3" / "production"


def signed_cortical_stress(sim, R_cell: float, n_planes: int = 12) -> dict:
    """Method-of-planes cortical stress, BOTH magnitude (gated) and SIGNED (net).

    Replicates ``h3_ku35_tension._tension_method_of_planes`` bond extraction + f_cut exactly,
    but returns mean(γ) signed (net contraction/extension) in addition to mean(|γ|) (the gated
    magnitude). Does not modify the gated estimator. Units: mN/m.
    """
    btypes = list(sim.state.bond_types)
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        pos_byTag = pos[inv]
        bg = np.asarray(s.bonds.group)
        bt = np.asarray(s.bonds.typeid)
    bond_force = None
    for f in sim.operations.integrator.forces:
        if hasattr(f, "params") and any(t in btypes for t in getattr(f, "params", {})):
            bond_force = f
            break
    if bond_force is None or bg.shape[0] == 0:
        return {"magnitude_mN_m": float("nan"), "signed_mN_m": float("nan"),
                "n_planes_positive": 0, "n_bonds": int(bg.shape[0])}
    k_arr = np.zeros(bg.shape[0])
    r0_arr = np.zeros(bg.shape[0])
    for i, tname in enumerate(btypes):
        try:
            kp = bond_force.params[tname]
            mask = (bt == i)
            k_arr[mask] = float(kp["k"])
            r0_arr[mask] = float(kp["r0"])
        except Exception:
            pass
    rA = pos_byTag[bg[:, 0]]
    rB = pos_byTag[bg[:, 1]]
    d = rB - rA
    L = np.linalg.norm(d, axis=1)
    L_safe = np.where(L > 0, L, 1.0)
    u = d / L_safe[:, None]
    T = k_arr * (L - r0_arr)            # + tension / − compression
    phi = (1 + 5 ** 0.5) / 2
    i = np.arange(n_planes, dtype=np.float64)
    z = 1 - 2 * (i + 0.5) / n_planes
    rxy = np.sqrt(1 - z * z)
    theta = 2 * np.pi * i / phi
    normals = np.stack([rxy * np.cos(theta), rxy * np.sin(theta), z], axis=1)
    gammas = []
    for n_hat in normals:
        a_side = rA @ n_hat
        b_side = rB @ n_hat
        crossing = (a_side * b_side) < 0
        if not crossing.any():
            gammas.append(0.0)
            continue
        f_cut = (T[crossing] * np.abs(u[crossing] @ n_hat))
        gammas.append(float(np.sum(f_cut)) / (2.0 * np.pi * R_cell))
    arr = np.asarray(gammas)
    return {
        "magnitude_mN_m": float(np.mean(np.abs(arr))) * 1e3,   # = the gated KU-3.5 value
        "signed_mN_m": float(np.mean(arr)) * 1e3,              # net: + contractile / − extensile
        "n_planes_positive": int(np.sum(arr > 0)),
        "n_planes": int(n_planes),
    }


def _run_arm(n_fil, n_motors, n_xl, *, n_warmup, n_sample, interval, bind_scale,
             kon_scale, device="gpu", seed=1):
    """Warm-up (unconstrained) → constrained production; sample SIGNED + magnitude stress."""
    cfg = _cfg_for(n_fil, n_motors, n_xl, "grip_walk", force_scaling=True, backbone_nm=300)
    p0 = resolve_h3_derived(cfg)
    comp = _resolve_compartments(p0)
    require_full_cell_physiological_baseline(comp)
    # warm-up (unconstrained) to relax overlaps
    _, _, _, dtc, hw = _build(cfg, stepping_mode="grip_walk", force_scaling=True,
                              constrained=False, compartments=comp, equilibrate=True,
                              n_warmup=n_warmup, device=device, kon_scale=kon_scale,
                              bind_scale=bind_scale, seed=seed)
    hw["sim"].run(0)
    pos_warm = _tagpos(hw["sim"])
    del hw
    # constrained production (the real KU-3.5 mode)
    p, _pm, _px, dtc, hc = _build(cfg, stepping_mode="grip_walk", force_scaling=True,
                                  constrained=True, compartments=comp, dtc=dtc,
                                  device=device, kon_scale=kon_scale, bind_scale=bind_scale,
                                  seed=seed)
    sim = hc["sim"]
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = pos_warm
    sim.state.set_snapshot(snap)
    sim.run(0)
    nca = p.n_filaments * p.beads_per_filament
    r0 = float(np.linalg.norm(_tagpos(sim)[:nca], axis=1).mean())
    rows = []
    for k in range(n_sample):
        sim.run(interval)
        st = signed_cortical_stress(sim, p.R_cell)
        rmean = float(np.linalg.norm(_tagpos(sim)[:nca], axis=1).mean())
        st["r_over_r0"] = rmean / r0
        rows.append(st)
        print(f"    [n_motors={n_motors}] s={k+1}/{n_sample} r/r0={rmean/r0:.5f} "
              f"|γ|={st['magnitude_mN_m']:.4e} SIGNED γ={st['signed_mN_m']:+.4e} mN/m "
              f"(+planes {st['n_planes_positive']}/{st['n_planes']})", flush=True)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description="Signed cortical contractility (contraction vs extension).")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--n-fil", type=int, default=600)
    ap.add_argument("--n-motors", type=int, default=100)
    ap.add_argument("--n-samples", type=int, default=6)
    ap.add_argument("--interval", type=int, default=8000)
    ap.add_argument("--n-warmup", type=int, default=8000)
    ap.add_argument("--bind-scale", type=float, default=6.0)
    ap.add_argument("--kon-scale", type=float, default=300.0)
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args, hoomd_module=hoomd)

    if args.smoke:
        n_fil, n_motors, n_sample, interval, n_warmup = 300, 60, 3, 3000, 3000
    else:
        n_fil, n_motors, n_sample, interval, n_warmup = (
            args.n_fil, args.n_motors, args.n_samples, args.interval, args.n_warmup)
    n_xl = int(1.5 * n_fil)

    print(f"[signed-contractility] n_fil={n_fil} n_xl={n_xl} | SIGNED method-of-planes "
          f"(contraction vs extension), motors ON vs OFF on the SAME constrained network.\n")
    print("  --- ARM OFF (motors=0, passive baseline) ---")
    off = _run_arm(n_fil, 0, n_xl, n_warmup=n_warmup, n_sample=n_sample, interval=interval,
                   bind_scale=args.bind_scale, kon_scale=args.kon_scale, device=args.device)
    print("  --- ARM ON (motors active) ---")
    on = _run_arm(n_fil, n_motors, n_xl, n_warmup=n_warmup, n_sample=n_sample, interval=interval,
                  bind_scale=args.bind_scale, kon_scale=args.kon_scale, device=args.device)

    def _mean(rows, key):
        return float(np.mean([r[key] for r in rows])) if rows else float("nan")

    signed_off, signed_on = _mean(off, "signed_mN_m"), _mean(on, "signed_mN_m")
    mag_off, mag_on = _mean(off, "magnitude_mN_m"), _mean(on, "magnitude_mN_m")
    delta_signed = signed_on - signed_off
    verdict = (
        "CONTRACTILE (myosin adds net tensile hoop stress)" if delta_signed > 1e-5 else
        "EXTENSILE (myosin net pushes apart — Stam low-crosslink rigid regime)" if delta_signed < -1e-5 else
        "LOCAL-DIPOLE CANCELLATION (net ~0 despite nonzero |stress| — the force-aggregation wall)"
    )
    res = {
        "config": {"n_fil": n_fil, "n_motors_on": n_motors, "n_xl": n_xl,
                   "bind_scale": args.bind_scale, "kon_scale": args.kon_scale,
                   "n_sample": n_sample, "interval": interval},
        "motors_off": {"signed_mN_m": signed_off, "magnitude_mN_m": mag_off,
                       "r_over_r0": _mean(off, "r_over_r0"), "rows": off},
        "motors_on": {"signed_mN_m": signed_on, "magnitude_mN_m": mag_on,
                      "r_over_r0": _mean(on, "r_over_r0"), "rows": on},
        "delta_signed_mN_m": delta_signed,
        "delta_magnitude_mN_m": mag_on - mag_off,
        "verdict": verdict,
        "interpretation": (
            "Δ_signed = signed_stress(motors ON) − signed_stress(OFF). The gated estimator "
            "reports |stress| and cannot see this sign. A positive Δ_signed is the FIRST "
            "evidence of net contraction; ≈0 confirms local-dipole cancellation (the wall); "
            "<0 means the current cortex is EXTENSILE. This gates whether the selector fix is "
            "needed and what it must achieve."
        ),
    }
    print(f"\n=== SIGNED cortical contractility (motors ON vs OFF) ===")
    print(f"  signed γ:   OFF={signed_off:+.4e}  ON={signed_on:+.4e}  Δ={delta_signed:+.4e} mN/m")
    print(f"  |γ| (gated): OFF={mag_off:.4e}  ON={mag_on:.4e}  Δ={mag_on-mag_off:+.4e} mN/m")
    print(f"  → VERDICT: {verdict}")

    _OUT.mkdir(parents=True, exist_ok=True)
    out = _OUT / "stage2_signed_contractility.json"
    out.write_text(json.dumps(res, indent=2))
    print(f"\n[json] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
"""KU-3.5 — myosin-active cortex under constrained-BD (foundational driver).

Builds a cortex + myosin cell, warms it up with the standard BAOAB at the
CFL dt to relax construction overlaps, then hands the relaxed configuration
to the constrained integrator (rigid actin backbone, M-SHAKE + Fixman) for a
fast-dt production run where myosin actively contracts the network.

Saves the cortex-actin trajectory (for visualisation) plus per-snapshot
diagnostics: myosin engagement, total binding/stepping, mean cortex radius
(a contraction proxy), constraint drift.

NOTE — cortical tension γ: the rigorous Young-Laplace / method-of-planes
measurement with RIGID constraints (the actin "tension" is a Lagrange
multiplier, not a harmonic-bond force) is a Sanity-Gate measurement-protocol
decision pending PI ratification (like the L_p band). This driver reports a
CONTRACTION PROXY (mean radius vs time) and myosin engagement, NOT a
ratified γ — KU-3.5 PASS is not declared here.

ERM is OFF by default: at the LJ-limited production dt (~0.001·τ_bend = 0.7 μs)
the ratified k_ERM=1e-4 (τ_ERM=3.9 μs) re-violates CFL; the k_ERM for this dt
is a PI item. Myosin-active contraction is demonstrable without ERM.

Usage:
    python ffn_sim/scripts/h3_ku35_tension.py --n-fil 150 --dt-factor 0.001
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yaml
import hoomd

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import resolve_cortex_myosin
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cortex.erm import resolve_erm, attach_erm_to_simulation
from ffn_sim.cell.cell import build_cortex_full_simulation

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"


def _tension_method_of_planes(sim, R_cell: float, n_planes: int = 12) -> float:
    """Mechanical cortical tension γ by method-of-planes (KU-3.5).

    For each of n_planes diametral planes through the cell center (axes
    spread on the sphere), sum the bond tensions of every harmonic bond
    that crosses the plane, projected onto the cut normal, and divide by
    the cut circumference 2πR_cell. Average across plane orientations.

    Includes ALL harmonic bonds present in the snap (myosin head springs +
    head-actin attach + xlink intra + xlink attach + any soft cortex
    component). Rigid backbone constraints (M-SHAKE Lagrange multipliers)
    are NOT included in this v1 — their shell-tension contribution is a
    separate calculation flagged for PI follow-up. Soft-bond tension is
    the dominant contractile signal once myosin engages.
    """
    btypes = list(sim.state.bond_types)
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position); tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
        pos_byTag = pos[inv]
        bg = np.asarray(s.bonds.group)
        bt = np.asarray(s.bonds.typeid)
    # Harmonic bond params: read from the integrator's bond force.
    bond_force = None
    for f in sim.operations.integrator.forces:
        if hasattr(f, "params") and any(t in btypes for t in getattr(f, "params", {})):
            bond_force = f; break
    if bond_force is None:
        return float("nan")
    # Per-bond k, r0 (in tag-ordered position frame).
    k_arr = np.zeros(bg.shape[0]); r0_arr = np.zeros(bg.shape[0])
    for i, tname in enumerate(btypes):
        try:
            kp = bond_force.params[tname]
            mask = (bt == i)
            k_arr[mask] = float(kp["k"]); r0_arr[mask] = float(kp["r0"])
        except Exception:
            pass
    # Bond endpoints + direction + scalar tension T = k·(|r|−r0).
    rA = pos_byTag[bg[:, 0]]; rB = pos_byTag[bg[:, 1]]
    d = rB - rA
    L = np.linalg.norm(d, axis=1); L_safe = np.where(L > 0, L, 1.0)
    u = d / L_safe[:, None]
    T = k_arr * (L - r0_arr)   # +T = tension (extension), −T = compression
    # Sample plane normals isotropically (Fibonacci-like).
    phi = (1 + 5 ** 0.5) / 2
    i = np.arange(n_planes, dtype=np.float64)
    z = 1 - 2 * (i + 0.5) / n_planes
    rxy = np.sqrt(1 - z * z); theta = 2 * np.pi * i / phi
    normals = np.stack([rxy * np.cos(theta), rxy * np.sin(theta), z], axis=1)
    gammas = []
    for n_hat in normals:
        a_side = rA @ n_hat; b_side = rB @ n_hat
        crossing = (a_side * b_side) < 0
        if not crossing.any():
            gammas.append(0.0); continue
        # Force along cut normal: |T| · |u · n_hat| (signed by tension dir).
        # The total force across the cut (tensile) = sum over crossing bonds.
        f_cut = (T[crossing] * (u[crossing] @ n_hat))
        # The cortex tension γ = force / circumference of the great circle.
        gammas.append(float(np.sum(f_cut)) / (2.0 * np.pi * R_cell))
    arr = np.asarray(gammas)
    # Use magnitude (sign tracks tension/compression); report mean of |γ|
    # — the tension scalar (Codex method-of-planes recommendation).
    return float(np.mean(np.abs(arr)))


def _tagpos(sim) -> np.ndarray:
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position); tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
        return pos[inv].copy()


def run(n_fil: int, seed: int, *, dt_factor: float = 0.001, with_xlinks: bool = True,
        with_erm: bool = True, k_erm_fast: float = 5.6e-5,
        n_warmup: int = 40_000, n_sample: int = 80, interval: int = 5_000,
        out: str | None = None) -> dict:
    cfg = yaml.safe_load(open(CFG))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    N = p.beads_per_filament; F = p.n_filaments; nca = F * N
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = dt_factor * tau_bend
    p_myo = resolve_cortex_myosin(cfg, dt=dtc)
    p_xl = resolve_crosslinkers(cfg, dt=dtc) if with_xlinks else None
    # ERM at fast dt: re-derive k_ERM from CFL (PI 2026-05-29).
    # k_ERM = cfl_safety · γ_b / dt → k=0.1·3.9e-10/dt ≈ 5.6e-5 at dt=0.7μs.
    p_erm = None
    if with_erm:
        p_erm = resolve_erm(cfg, kT=p.kT, R_cell=p.R_cell)
        from dataclasses import replace
        try: p_erm = replace(p_erm, k_ERM=k_erm_fast)
        except Exception: p_erm.k_ERM = k_erm_fast
    dev = hoomd.device.CPU(notice_level=0)

    # Phase 1 — warm-up (standard BAOAB, CFL dt) to relax overlaps.
    hw = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo, device=dev, with_baoab=True,
        constrained=False, rng=np.random.default_rng(seed))
    # Attach ERM at WARM-UP dt (cfl) with config k_ERM (=1e-4); stable here.
    if with_erm:
        from dataclasses import replace as _replace
        p_erm_warm = resolve_erm(cfg, kT=p.kT, R_cell=p.R_cell)
        attach_erm_to_simulation(
            hw["sim"], p_erm_warm, actin_cortex_tag_range=(0, nca),
            gamma_b=p.gamma_b, cfl_safety_factor=p.cfl_safety_factor,
            cfl_strict=True,
        )
    hw["sim"].run(0); hw["sim"].run(n_warmup)
    pos_warm = _tagpos(hw["sim"])
    del hw

    # Phase 2 — constrained production (rigid actin backbone, fast dt).
    hc = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo, device=dev, with_baoab=True,
        constrained=True, constrained_dt=dtc, rng=np.random.default_rng(seed))
    # Attach ERM at fast dt with re-derived k_ERM.
    if with_erm:
        attach_erm_to_simulation(
            hc["sim"], p_erm, actin_cortex_tag_range=(0, nca),
            gamma_b=p.gamma_b, cfl_safety_factor=p.cfl_safety_factor,
            cfl_strict=True,
        )
    sim = hc["sim"]; ma = hc["myosin_action"]; act = hc["baoab_action"]
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = pos_warm
    sim.state.set_snapshot(snap)
    sim.run(0)

    frames = np.empty((n_sample, F, N, 3))
    diag = []
    r0 = float(np.linalg.norm(_tagpos(sim)[:nca], axis=1).mean())
    for k in range(n_sample):
        sim.run(interval)
        r = _tagpos(sim)
        frames[k] = r[:nca].reshape(F, N, 3)
        rmean = float(np.linalg.norm(r[:nca], axis=1).mean())
        # KU-3.5 method-of-planes cortical tension (soft-bond contribution).
        gamma = _tension_method_of_planes(sim, p.R_cell)
        diag.append(dict(
            step=int(sim.timestep), r_cortex_um=rmean * 1e6,
            r_over_r0=rmean / r0, myosin_engaged=int(ma.n_engaged),
            bind_total=int(ma.n_bind_total), step_advances=int(ma.n_step_advances_total),
            tension_mN_per_m=gamma * 1e3,
            max_drift=float(act.max_constraint_drift),
        ))

    outdir = PKG / "outputs" / "h3" / "production" / "ku35"
    outdir.mkdir(parents=True, exist_ok=True)
    if out is None:
        out_json = outdir / f"ku35_n{n_fil}_s{seed}.json"
        frames_npz = outdir / f"ku35_frames_n{n_fil}_s{seed}.npz"
    else:
        out_json = Path(out)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        frames_npz = out_json.with_suffix(".npz")
    np.savez_compressed(frames_npz, frames=frames)
    # Plateau (last third) tension mean.
    last_third = diag[max(0, len(diag) * 2 // 3):]
    g_plateau = float(np.mean([d["tension_mN_per_m"] for d in last_third]))
    result = dict(
        n_fil=int(F), seed=int(seed), dt_s=float(dtc), dt_factor=float(dt_factor),
        dt_speedup_vs_cfl=float(dtc / p.dt_cfl), r0_um=r0 * 1e6,
        with_xlinks=bool(with_xlinks), with_erm=bool(with_erm),
        k_erm_fast=float(k_erm_fast) if with_erm else None,
        r_final_over_r0=diag[-1]["r_over_r0"], myosin_engaged_final=diag[-1]["myosin_engaged"],
        step_advances_final=diag[-1]["step_advances"],
        tension_plateau_mN_per_m=g_plateau, tension_final_mN_per_m=diag[-1]["tension_mN_per_m"],
        ku35_target_mN_per_m=0.5, ku35_band_mN_per_m=[0.35, 0.65],
        n_motors=int(p_myo.n_motors_per_cell), max_drift=diag[-1]["max_drift"],
        note="method-of-planes γ (soft bonds; rigid-constraint Lagrange shell-tension separate; PI ratification pending)",
        diag=diag,
    )
    out_json.write_text(json.dumps(result, indent=2))
    print("RESULT " + json.dumps({k: v for k, v in result.items() if k != "diag"}), flush=True)
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dt-factor", type=float, default=0.001)
    ap.add_argument("--no-xlinks", dest="with_xlinks", action="store_false", default=True)
    ap.add_argument("--no-erm", dest="with_erm", action="store_false", default=True)
    ap.add_argument("--k-erm", type=float, default=5.6e-5)
    ap.add_argument("--n-warmup", type=int, default=40_000)
    ap.add_argument("--n-sample", type=int, default=80)
    ap.add_argument("--interval", type=int, default=5_000)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    run(args.n_fil, args.seed, dt_factor=args.dt_factor,
        with_xlinks=args.with_xlinks, with_erm=args.with_erm, k_erm_fast=args.k_erm,
        n_warmup=args.n_warmup, n_sample=args.n_sample, interval=args.interval, out=args.out)


if __name__ == "__main__":
    main()

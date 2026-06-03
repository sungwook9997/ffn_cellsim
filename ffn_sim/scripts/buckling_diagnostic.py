"""STAGE-2 buckling diagnostic — does the constrained cortex actually buckle?

The constrained M-SHAKE backbone is rigid in STRETCH but FLEXIBLE in BENDING (only
the stretch DOF is SHAKE'd; the cortex-angle harmonic is preserved). That asymmetric
compliance is the whole design rationale for contractility (Murrell-Gardel 2012, Lenz
2012): a semiflexible filament BUCKLES under compression, breaking the tension/
compression symmetry → net contraction. This script CHECKS the assumption is live:
it measures the cortex-angle distribution θ (angle at the middle bead of each
backbone triplet; straight filament → θ≈180°, buckled → θ ≪ 180°) at t=0 vs after
running under motor load, for motors ON vs OFF. If motor activity drives θ away from
180° (filaments bend/buckle), buckling is happening; if θ stays ≈180°, the filaments
stay straight and buckling-mediated contractility is absent.

Run:  PYTHONPATH=. python -m ffn_sim.scripts.buckling_diagnostic --n-motors 100
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml
import hoomd

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cortex.myosin import resolve_cortex_myosin
from ffn_sim.cell.cell import build_cortex_full_simulation

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"
CORTEX_ANGLE = "cortex-angle"


def _angles_theta(sim) -> np.ndarray:
    """Return θ [deg] at the middle bead of every cortex-angle triplet."""
    snap = sim.state.get_snapshot()              # tag-ordered; has .types
    if snap.communicator.rank != 0:
        return np.empty(0)
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    ag = np.asarray(snap.angles.group, dtype=np.int64)
    at = np.asarray(snap.angles.typeid, dtype=np.int64)
    atypes = list(snap.angles.types)
    if CORTEX_ANGLE not in atypes or ag.shape[0] == 0:
        return np.empty(0)
    cid = atypes.index(CORTEX_ANGLE)
    trip = ag[at == cid]
    # get_snapshot is tag-ordered, so group entries index pos directly.
    p0 = pos[trip[:, 0]]; p1 = pos[trip[:, 1]]; p2 = pos[trip[:, 2]]
    v1 = p0 - p1; v2 = p2 - p1
    v1 /= np.linalg.norm(v1, axis=1, keepdims=True).clip(min=1e-30)
    v2 /= np.linalg.norm(v2, axis=1, keepdims=True).clip(min=1e-30)
    cos = np.clip(np.einsum("ij,ij->i", v1, v2), -1.0, 1.0)
    return np.degrees(np.arccos(cos))


def _summary(theta: np.ndarray, label: str) -> None:
    if theta.size == 0:
        print(f"  [{label}] no cortex angles"); return
    buckled = (theta < 150.0).mean() * 100.0
    sharp = (theta < 120.0).mean() * 100.0
    print(f"  [{label}] n={theta.size} mean θ={theta.mean():.1f}° "
          f"median={np.median(theta):.1f}° min={theta.min():.1f}° "
          f"| buckled(<150°)={buckled:.1f}% sharp(<120°)={sharp:.1f}%", flush=True)


def build(n_fil, n_motors, n_xl, kon_scale, bind_scale, seed=1):
    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["R_cell"] = 7.5e-6
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors
    cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    cfg["cortex"]["myosin"]["mesoscale_force_scaling"] = True
    cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = int(n_xl)
    p = resolve_h3_derived(cfg)
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = 0.001 * tau_bend
    p_xl = resolve_crosslinkers(cfg, dt=dtc)
    if kon_scale != 1.0:
        p_xl = replace(p_xl, k_on=p_xl.k_on * kon_scale)
    if bind_scale != 1.0:
        p_xl = replace(p_xl, max_bind_dist=p_xl.max_bind_dist * bind_scale)
    p_myo = resolve_cortex_myosin(cfg, dt=dtc, R_cell=p.R_cell) if n_motors > 0 else None
    dev = hoomd.device.CPU(notice_level=0)
    kw = dict(p_xlinks=p_xl, device=dev, with_baoab=True, constrained=True,
              constrained_dt=dtc, rng=np.random.default_rng(seed))
    if p_myo is not None:
        kw["p_myosin"] = p_myo
    hw = build_cortex_full_simulation(p, **kw)
    return hw["sim"], p


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--n-motors", type=int, default=100)
    ap.add_argument("--n-xl", type=int, default=1500)
    ap.add_argument("--kon-scale", type=float, default=300.0)
    ap.add_argument("--bind-scale", type=float, default=6.0)
    ap.add_argument("--steps", type=int, default=40000)
    args = ap.parse_args()

    print("=== STAGE-2 buckling diagnostic (cortex-angle θ; straight≈180°) ===")
    print(f"HOOMD {hoomd.version.version}", flush=True)
    for n_mot, tag in ((args.n_motors, f"motors ON (n={args.n_motors})"), (0, "motors OFF")):
        sim, p = build(args.n_fil, n_mot, args.n_xl, args.kon_scale, args.bind_scale)
        sim.run(0)
        _summary(_angles_theta(sim), f"{tag} t=0")
        sim.run(args.steps)
        _summary(_angles_theta(sim), f"{tag} after {args.steps} steps")
        del sim
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

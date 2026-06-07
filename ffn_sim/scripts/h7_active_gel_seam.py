"""H.7 active-gel γ-seam M1 — DIAGNOSTIC postprocessor (PI 2026-06-07).

Builds the (connected-mesh) FG cortex, extracts the active-stress drive + the
generation ceiling, relaxes the 0-D active-Maxwell gel over the seconds timescale
the MD can't reach, and reports whether the active cortical-tension floor is a
TIMESCALE GAP or a GENERATION LIMIT. This is a DIAGNOSTIC, not a solution: if the
generation ceiling is below band, the seam cannot close it (→ M2 / re-target).

Usage:
    python -m ffn_sim.scripts.h7_active_gel_seam --n-filaments 200 --warmup 1500 \
        --sample 1000 --device cpu --allow-cpu-dev
"""

from __future__ import annotations

import argparse
import math
from copy import deepcopy
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest  # noqa: E402
from ffn_sim.common.production_policy import (  # noqa: E402
    add_production_device_args, validate_production_device_args,
)
from ffn_sim.cortex.active_gel_seam import (  # noqa: E402
    BAND_N_PER_M, MCF7_INTERPHASE_GAMMA, diagnose_seam, relax_active_maxwell,
)
from ffn_sim.cortex.cortical_tension import measure_cortical_tension  # noqa: E402

_MN = 1.0e3


def _engaged_heads(cell) -> int:
    snap = cell.simulation.state.get_snapshot()
    types = list(snap.bonds.types)
    tid = np.asarray(snap.bonds.typeid)
    import collections
    c = collections.Counter(tid.tolist())
    return sum(c.get(i, 0) for i, t in enumerate(types) if "myosin_attach" in t)


def run(n_filaments, warmup, sample, device, seed):
    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    if n_filaments is not None:
        # MERGE (preserve the production myosin config from mcf7_baseline.yaml).
        co = manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})
        co["n_filaments"] = int(n_filaments)
        co["demo_mode"] = True
    cell = build_baseline_cell(manifest=manifest, device=device, seed=seed,
                               connected_mesh=True, equilibrate=True,
                               equilibrate_steps=warmup,
                               equilibrate_softstart_steps=max(50, warmup // 2))
    if sample > 0:
        cell.simulation.run(sample)
    pm, pc = cell.p_myosin, cell.p_cortex
    g = measure_cortical_tension(cell.simulation, R_cell=pc.R_cell,
                                 p_enclosed_volume=cell.p_enclosed_volume)
    n_eng = _engaged_heads(cell)
    n_total = pm.n_motors_per_cell * 2 * pm.n_heads_per_side
    # turnover relaxation time: tau = tau_half/ln2 (Chugh 2017 tau_half=10 s if turnover off).
    tau_half = getattr(cell.p_turnover, "tau_half", None) or 10.0
    tau = tau_half / math.log(2.0)
    dg = diagnose_seam(
        gamma_soft=g["gamma_soft"], gamma_rigid=g["gamma_rigid"],
        n_eng=n_eng, n_total_heads=n_total, f_stall_per_head=pm.F_stall_per_head,
        ell_dipole=pm.backbone_length, area=pc.cortex_surface_area,
        thickness=pc.cortex_thickness, tau=tau,
        heads_per_minifilament=pm.n_heads_per_side * 2,
    )
    return cell, g, dg


def _figure(dg, out: Path):
    h = dg.thickness
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    # Panel 1: steady γ (realized / capacity / ceiling) vs band + MCF7 datum (log).
    names = ["realized\n(FG soft)", f"capacity\n({dg.n_eng} eng heads)", f"ceiling\n(all {dg.n_total_heads})"]
    vals = [max(dg.gamma_ss_realized, 1e-9) * _MN, max(dg.gamma_ss_capacity, 1e-9) * _MN,
            max(dg.gamma_ss_ceiling, 1e-9) * _MN]
    ax1.bar(names, vals, color=["#bbb", "#e8902a", "#8a1f1f"])
    ax1.axhspan(BAND_N_PER_M[0] * _MN, BAND_N_PER_M[1] * _MN, color="#2f6fb0", alpha=0.18,
                label="band [0.35,0.65]")
    ax1.axhline(MCF7_INTERPHASE_GAMMA * _MN, color="green", ls="--", label="MCF7 datum 0.27 (Hosseini)")
    ax1.set_yscale("log"); ax1.set_ylabel("steady γ (mN/m, log)")
    ax1.set_title("M1 diagnostic — active generation vs band"); ax1.legend(fontsize=8)
    # Panel 2: the active-Maxwell relaxation transient (ceiling drive), seconds-scale.
    sigma_a = dg.gamma_ss_ceiling / h
    sigma_pre = dg.gamma_rigid / h
    t, sig, _ = relax_active_maxwell(sigma_a, sigma_pre, dg.tau)
    ax2.plot(t, sig * h * _MN, "-", color="#8a1f1f")
    ax2.axhspan(BAND_N_PER_M[0] * _MN, BAND_N_PER_M[1] * _MN, color="#2f6fb0", alpha=0.18)
    ax2.set_xlabel("time (s)"); ax2.set_ylabel("γ (mN/m)")
    ax2.set_title(f"active-Maxwell relaxation (τ={dg.tau:.1f}s) toward the ceiling\n"
                  "(MD reaches only ~ms ≪ τ)")
    fig.suptitle("H.7 active-gel seam M1 (DIAGNOSTIC, not a solution)", fontweight="bold")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150); plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=None)
    ap.add_argument("--warmup", type=int, default=2000)
    ap.add_argument("--sample", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1]
                    / "outputs" / "h7" / "figs" / "h7_active_gel_seam_M1.png"))
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)
    import hoomd
    dev = hoomd.device.GPU(notice_level=0) if args.device == "gpu" else hoomd.device.CPU(notice_level=0)
    cell, g, dg = run(args.n_filaments, args.warmup, args.sample, dev, args.seed)
    _figure(dg, Path(args.out))
    print("=" * 70, flush=True)
    print("H.7 ACTIVE-GEL SEAM M1 — DIAGNOSTIC (not a solution)", flush=True)
    print(f"  motor density = {dg.motor_density_per_um2:.3f} /µm² (Salbreux target 3/µm²); "
          f"engaged heads {dg.n_eng}/{dg.n_total_heads}; F_stall {dg.f_stall_per_head*1e12:.2f} pN/head; "
          f"ℓ_dipole {dg.ell_dipole*1e9:.0f} nm", flush=True)
    print(f"  steady γ:  realized {dg.gamma_ss_realized*_MN:.4f} | capacity {dg.gamma_ss_capacity*_MN:.4f} "
          f"| CEILING {dg.gamma_ss_ceiling*_MN:.4f} mN/m   vs band [0.35,0.65], MCF7 0.27", flush=True)
    print(f"  VERDICT: {dg.verdict}", flush=True)
    print(f"  fig: {args.out}", flush=True)
    print("=" * 70, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
"""KU-3.5 assay ladder — GATE 2: PATCH CONTRACTILITY (FREE actin).

GATE 1 (`tests/test_ku35_dipole_gate.py`) established that the *single*
minifilament, hand-built into a complete antiparallel pair on RIGIDLY-ANCHORED
filaments, SUSTAINS a contractile dipole (P<0, coherence −1). GATE 1 is
isometric (the actin cannot move) and the pair is built by hand. GATE 2 asks the
next question up the ladder:

    Does a FREE (unconstrained) actin cortex+myosin PATCH — heads binding
    de-novo, actin free to move under the BAOAB integrator at the physiological
    MCF7 cytoplasm viscosity (η = 65.9 Pa·s) — actually CONTRACT when motors are
    ON, and is there a BIPHASIC connectivity sweet-spot (Ennomani 2016:
    contraction maximal at INTERMEDIATE crosslink connectivity, low at very-low
    and very-high)?

This is a STANDALONE assay fixture built on the production assembler
(`cell.build_cortex_full_simulation`, ``constrained=False`` = FREE actin,
``connected_mesh=True``, all at the MCF7 cytoplasm η). It is NOT a production
γ run — it is the diagnostic that localises the KU-3.5-active floor by asking
whether contraction emerges from the FREE mechanism at patch scale.

Observables (ACTIVE — never g_rigid / Lagrange structural tension):
  (a) CONTRACTION = the cortex-actin radius-of-gyration ``Rg`` and mean radius
      ``r`` over time, motors ON vs OFF (n_motors=0). Contraction ⇒ Rg/Rg0 and
      r/r0 DECREASE more with motors ON than OFF.
  (b) STRESSLET ΣP (`h3_ku35_stresslet.minifilament_stresslets`) over complete
      bipolar pairs: ΣP < 0 = net contractile; coherence = ΣP/Σ|P| ∈ [−1,+1]
      (−1 = perfectly contractile). Plus frac_complete_pairs (the recruitment
      ceiling) and max F_bond/F_stall (Hill validity).

Connectivity sweep: vary the connected-mesh structural coordination ``cm_z``
across under- / intermediate- / over-connected, measure the contraction
magnitude (ΔRg and ΣP) vs connectivity, look for a biphasic peak.

v0 NOTE (HARD constraint compliance): at the literal production NMIIA velocity
(0.2 µm/s) and the bend-CFL dt (~7e-7 s) the grip walks one ℓ₀ bead in ~3.6M
ticks — far beyond a seconds-long assay. As in the canonical STAGE-1 fixture
(`tests/test_stage1_grip_walk.py`), this STANDALONE fixture cranks ``v0`` purely
to set the mechanism RATE; the contractile *mechanism* (bipolar dipole → P<0) is
v0-independent and the Hill/Bell-Evans force ceiling is UNCHANGED (we report
max F_bond/F_stall so any super-stall is visible). NO force_scaling, v0_accel,
or passive-stiffness is tuned to inflate FORCE or make a gate pass.

Run:
  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_patch_contractility --smoke
  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_patch_contractility --sweep
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

from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation
from ffn_sim.archive.hoomd_legacy.cell.cytoplasm import resolve_cytoplasm
from ffn_sim.common.production_policy import (
    add_production_device_args,
    require_production_device,
    validate_production_device_args,
)
from ffn_sim.scripts.h3_ku35_stresslet import stresslet_ledger

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"

# MCF7 patch defaults (small-scale standalone fixture; the cytoplasm η is the
# physiological MCF7 value, the stabilizer — NOT water).
MCF7_R_CELL = 7.5e-6
MCF7_BACKBONE_NM = 300.0          # literature NMII bipolar minifilament length
# v0 multiplier: rate-only crank so the grip-walk mechanism runs in a
# seconds-long assay (mechanism is v0-independent; force ceiling unchanged).
DEFAULT_V0_ACCEL = 3.0e6


def _tagpos(sim) -> np.ndarray:
    """Tag-ordered positions (so cortex-actin tags [0, nca) are the actin cloud)."""
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        return pos[inv].copy()


def _patch_cfg(*, n_fil: int, n_motors: int, n_xl: int) -> dict:
    """MCF7 patch config (small scale, grip_walk, literature backbone)."""
    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["R_cell"] = MCF7_R_CELL
    cfg["cortex"]["n_filaments"] = int(n_fil)
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = int(n_motors)
    cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    cfg["cortex"]["myosin"]["backbone_length"] = MCF7_BACKBONE_NM * 1e-9
    cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = int(n_xl)
    return cfg


def build_patch(
    *,
    motors_off: bool,
    n_fil: int = 120,
    n_motors: int = 60,
    n_xl: int = 400,
    cm_z: float = 3.3,
    cm_bundle_mult: int = 2,
    v0_accel: float = DEFAULT_V0_ACCEL,
    seed: int = 1,
    device: str = "gpu",
    allow_cpu_dev: bool = False,
):
    """Build a FREE-actin MCF7 cortex+myosin patch (connected mesh, cytoplasm ON).

    Returns ``(p_cortex, p_myo, handles)``.  ``motors_off`` builds the passive
    baseline (no myosin) for the ON−OFF delta.  ``constrained=False`` ⇒ the
    actin is FREE to move under the BAOAB integrator at the MCF7 cytoplasm η.
    """
    require_production_device(
        device, allow_cpu_dev=allow_cpu_dev, hoomd_module=hoomd
    )

    cfg = _patch_cfg(
        n_fil=n_fil, n_motors=(0 if motors_off else n_motors), n_xl=n_xl
    )
    p = resolve_h3_derived(cfg)
    dtc = 0.001 * p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    p_myo = resolve_cortex_myosin(cfg, dt=dtc, R_cell=p.R_cell)
    if v0_accel != 1.0:
        # rate-only crank (mechanism RATE, NOT force) — F ceiling unchanged.
        p_myo = replace(p_myo, v0_per_head=p_myo.v0_per_head * v0_accel)
    p_xl = resolve_crosslinkers(cfg, dt=dtc)
    p_cyto = resolve_cytoplasm(cell_type="MCF7")        # η = 65.9 Pa·s
    dev = (hoomd.device.GPU(notice_level=0) if device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    hw = build_cortex_full_simulation(
        p,
        p_xlinks=p_xl,
        p_myosin=(None if motors_off else p_myo),
        device=dev,
        with_baoab=True,
        constrained=False,                 # FREE actin
        connected_mesh=True,
        cm_z_struct=cm_z,
        cm_bundle_mult=cm_bundle_mult,
        p_cytoplasm=p_cyto,
        rng=np.random.default_rng(seed),
    )
    return p, p_myo, hw


def actin_radius_metrics(p_cortex, sim) -> tuple[float, float]:
    """Return ``(Rg, r_mean)`` of the cortex-actin bead cloud [m].

    ``Rg`` = sqrt(mean |x − x_cm|²); ``r_mean`` = mean |x| (radius about the
    box origin, where the cortex sphere is centred at construction).
    """
    nca = p_cortex.n_filaments * p_cortex.beads_per_filament
    pos = _tagpos(sim)[:nca]
    cm = pos.mean(axis=0)
    rg = float(np.sqrt(((pos - cm) ** 2).sum(axis=1).mean()))
    r_mean = float(np.linalg.norm(pos, axis=1).mean())
    return rg, r_mean


def run_arm(
    *,
    motors_off: bool,
    n_steps: int,
    n_blocks: int = 1,
    **build_kw,
) -> dict:
    """Build a patch and run it ``n_blocks`` × ``n_steps``, recording the active
    observables.  Returns a dict with the Rg/r traces and the final stresslet
    ledger summary (motors-ON only)."""
    p, p_myo, hw = build_patch(motors_off=motors_off, **build_kw)
    sim = hw["sim"]
    myo = hw.get("myosin_action")
    sim.run(0)
    rg0, r0 = actin_radius_metrics(p, sim)
    rg_trace, r_trace, P_trace, coh_trace, complete_trace, maxF_trace = (
        [rg0], [r0], [], [], [], [],
    )
    for _ in range(n_blocks):
        sim.run(n_steps)
        rg, r = actin_radius_metrics(p, sim)
        rg_trace.append(rg)
        r_trace.append(r)
        if not motors_off and myo is not None:
            led = stresslet_ledger(
                sim, p_myo=p_myo, myosin_action=myo,
                beads_per_filament=p.beads_per_filament,
            )
            s = led["summary"]
            P_trace.append(s["sumP"])
            coh_trace.append(s["coherence"])
            complete_trace.append(s["frac_complete_pairs"])
            maxF = max(
                (x["max_Fbond_over_Fstall"] for x in led["stresslets"]
                 if np.isfinite(x["max_Fbond_over_Fstall"])),
                default=float("nan"),
            )
            maxF_trace.append(maxF)
    out = dict(
        rg0=rg0, r0=r0,
        rg_final=rg_trace[-1], r_final=r_trace[-1],
        rg_ratio=rg_trace[-1] / rg0, r_ratio=r_trace[-1] / r0,
        rg_trace=rg_trace, r_trace=r_trace,
    )
    if not motors_off:
        out.update(
            sumP=(P_trace[-1] if P_trace else float("nan")),
            coherence=(coh_trace[-1] if coh_trace else float("nan")),
            frac_complete=(complete_trace[-1] if complete_trace else float("nan")),
            max_Fbond_over_Fstall=(maxF_trace[-1] if maxF_trace else float("nan")),
            P_trace=P_trace, coherence_trace=coh_trace,
            complete_trace=complete_trace,
        )
    return out


def contraction_delta(
    *, n_steps: int, n_blocks: int = 1, **build_kw
) -> dict:
    """Motors-ON minus motors-OFF contraction at one connectivity.

    Returns the ON / OFF arm dicts plus the ON−OFF Rg-ratio and r-ratio deltas
    (negative ⇒ motors ON contract MORE than the passive baseline).
    """
    on = run_arm(motors_off=False, n_steps=n_steps, n_blocks=n_blocks, **build_kw)
    off = run_arm(motors_off=True, n_steps=n_steps, n_blocks=n_blocks, **build_kw)
    return dict(
        on=on, off=off,
        d_rg_ratio=on["rg_ratio"] - off["rg_ratio"],   # <0 ⇒ ON contracts more
        d_r_ratio=on["r_ratio"] - off["r_ratio"],
        sumP=on["sumP"], coherence=on["coherence"],
        frac_complete=on["frac_complete"],
        max_Fbond_over_Fstall=on["max_Fbond_over_Fstall"],
    )


def connectivity_sweep(
    *,
    cm_z_values=(1.5, 3.3, 6.0),
    n_steps: int = 20000,
    n_blocks: int = 1,
    **build_kw,
) -> list[dict]:
    """Sweep the connected-mesh structural coordination cm_z (under/inter/over)
    and measure contraction (ΔRg, ΣP) vs connectivity — biphasic check."""
    rows = []
    for cm_z in cm_z_values:
        t0 = time.time()
        d = contraction_delta(
            n_steps=n_steps, n_blocks=n_blocks, cm_z=cm_z, **build_kw
        )
        rows.append(dict(
            cm_z=cm_z,
            on_rg_ratio=d["on"]["rg_ratio"], off_rg_ratio=d["off"]["rg_ratio"],
            on_r_ratio=d["on"]["r_ratio"], off_r_ratio=d["off"]["r_ratio"],
            d_rg_ratio=d["d_rg_ratio"], d_r_ratio=d["d_r_ratio"],
            sumP=d["sumP"], coherence=d["coherence"],
            frac_complete=d["frac_complete"],
            max_Fbond_over_Fstall=d["max_Fbond_over_Fstall"],
            wall_s=time.time() - t0,
        ))
        print(
            f"  cm_z={cm_z:>4}: complete={d['frac_complete']:.3f} "
            f"ΣP={d['sumP']:.3e} coh={d['coherence']} "
            f"Rg/Rg0(ON)={d['on']['rg_ratio']:.7f} (OFF)={d['off']['rg_ratio']:.7f} "
            f"ΔRg(ON−OFF)={d['d_rg_ratio']:+.2e} maxF/Fst={d['max_Fbond_over_Fstall']:.2f} "
            f"[{rows[-1]['wall_s']:.0f}s]",
            flush=True,
        )
    return rows


def _make_figure(rows: list[dict], out_path: Path) -> None:
    """Connectivity curve: ΣP and ΔRg(ON−OFF) vs cm_z, plus complete-pair frac."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # noqa: BLE001
        print(f"  [VIZ] WARN matplotlib unavailable: {exc}", flush=True)
        return
    cm_z = [r["cm_z"] for r in rows]
    sumP = [r["sumP"] for r in rows]
    d_rg = [r["d_rg_ratio"] for r in rows]
    frac = [r["frac_complete"] for r in rows]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(cm_z, sumP, "o-", color="C3")
    ax[0].axhline(0.0, ls="--", lw=0.8, color="0.5")
    ax[0].set_xlabel("connected-mesh coordination z (struct)")
    ax[0].set_ylabel("ΣP over complete pairs  [N·m]  (<0 = contractile)")
    ax[0].set_title("Stresslet contractility vs connectivity")
    ax[1].plot(cm_z, d_rg, "s-", color="C0")
    ax[1].axhline(0.0, ls="--", lw=0.8, color="0.5")
    ax[1].set_xlabel("connected-mesh coordination z (struct)")
    ax[1].set_ylabel("ΔRg/Rg0  (motors ON − OFF)  (<0 = ON contracts more)")
    ax[1].set_title("Bulk contraction ON−OFF vs connectivity")
    ax[2].plot(cm_z, frac, "^-", color="C2")
    ax[2].set_xlabel("connected-mesh coordination z (struct)")
    ax[2].set_ylabel("frac_complete_pairs (bipolar recruitment)")
    ax[2].set_title("Recruitment ceiling vs connectivity")
    ax[2].set_ylim(bottom=0.0)
    fig.suptitle(
        "KU-3.5 GATE 2 — FREE-actin PATCH contractility + connectivity sweep "
        "(MCF7 η=65.9 Pa·s)",
        fontsize=11,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)
    print(f"  [VIZ] wrote {out_path}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--smoke", action="store_true",
                    help="tiny single ON/OFF delta at cm_z=3.3 (fast)")
    ap.add_argument("--sweep", action="store_true",
                    help="full connectivity sweep + figure")
    ap.add_argument("--n-fil", type=int, default=120)
    ap.add_argument("--n-motors", type=int, default=60)
    ap.add_argument("--n-xl", type=int, default=400)
    ap.add_argument("--n-steps", type=int, default=20000)
    ap.add_argument("--n-blocks", type=int, default=1)
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    print(f"=== KU-3.5 GATE 2 — FREE-actin PATCH contractility "
          f"(MCF7 η=65.9 Pa·s, n_fil={args.n_fil}, n_motors={args.n_motors}) ===",
          flush=True)
    print(f"HOOMD {hoomd.version.version}", flush=True)

    if args.smoke:
        d = contraction_delta(
            n_steps=4000, n_blocks=1, n_fil=args.n_fil, n_motors=args.n_motors,
            n_xl=args.n_xl, cm_z=3.3, device=args.device,
            allow_cpu_dev=args.allow_cpu_dev,
        )
        print(f"  [SMOKE] ΣP={d['sumP']:.3e} coh={d['coherence']} "
              f"complete={d['frac_complete']:.3f} "
              f"Rg/Rg0(ON)={d['on']['rg_ratio']:.7f} (OFF)={d['off']['rg_ratio']:.7f} "
              f"ΔRg={d['d_rg_ratio']:+.2e} maxF/Fst={d['max_Fbond_over_Fstall']:.2f}",
              flush=True)
        return 0

    cm_z_values = (1.5, 3.3, 6.0) if args.sweep else (3.3,)
    t0 = time.time()
    rows = connectivity_sweep(
        cm_z_values=cm_z_values, n_steps=args.n_steps, n_blocks=args.n_blocks,
        n_fil=args.n_fil, n_motors=args.n_motors, n_xl=args.n_xl,
        device=args.device, allow_cpu_dev=args.allow_cpu_dev,
    )
    print(f"\n=== sweep DONE in {time.time()-t0:.0f}s ===", flush=True)

    # honest verdict
    sumPs = [r["sumP"] for r in rows]
    cohs = [r["coherence"] for r in rows]
    d_rgs = [r["d_rg_ratio"] for r in rows]
    contractile_stresslet = all((not np.isfinite(c)) or c <= 0 for c in cohs)
    bulk_contracts = any(d < -1e-5 for d in d_rgs)
    print(f"  complete pairs contractile (coherence ≤ 0): {contractile_stresslet}",
          flush=True)
    print(f"  FREE patch bulk-contracts (ON−OFF ΔRg < −1e-5): {bulk_contracts}",
          flush=True)
    # biphasic? peak |ΣP| at the INTERIOR cm_z, not the endpoints.
    if len(rows) >= 3:
        absP = [abs(p) for p in sumPs]
        biphasic = absP[1] > absP[0] and absP[1] > absP[-1]
        print(f"  biphasic |ΣP| (interior > both endpoints): {biphasic} "
              f"(|ΣP|={[f'{a:.2e}' for a in absP]})", flush=True)

    if args.sweep:
        _make_figure(
            rows,
            PKG / "outputs" / "h3" / "figs" / "ku35_patchcontract_sweep.png",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

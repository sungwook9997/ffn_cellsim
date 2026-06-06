"""H.7 GATE-B — emergent cortical tension (gamma) at the settled operating point.

Assembles the FULL physiological MCF7 cell from configs/mcf7_baseline.yaml (with
FA-adhesion ENABLED) via the manifest loader -> A1 Cell.build, settles it onto
the substrate with an equilibration prelude (a raw FA-adhered free run trips the
BAOAB guard), runs the rigid M-SHAKE backbone with Lagrange-lambda recording,
then measures the THREE cortical-tension channels SEPARATELY via the B4 unified
estimator (ffn_sim/cortex/cortical_tension.py).

GATE-B CONTRACT (HARD, CLAUDE.md + PI 2026-06-05/06):
  * gamma is reported as DISTINCT channels — active (soft-MOP), rigid (Lagrange),
    passive (turgor Young-Laplace). The turgor channel is NEVER folded into a
    single "total": an in-band total must NOT be read as active-cortex closure.
  * This is the EMERGENT operating-point measurement, not a dial. No KU-3.5
    conclusion is drawn here; the band [0.35,0.65] mN/m is an overlay only.
  * AUTHORITATIVE gamma requires the FULL production scale on GPU (gbook). A
    small --n-filaments run is a PIPELINE SMOKE, explicitly labelled as such.

Usage (smoke):
    python -m ffn_sim.scripts.h7_gate_b --n-filaments 160 --n-nuc-beads 400 \
        --device cpu --allow-cpu-dev --warmup 300 --sample 200
Usage (full production scale, gbook GPU):
    python -m ffn_sim.scripts.h7_gate_b --device gpu   # manifest defaults (1000 fil)
"""

from __future__ import annotations

import argparse
from copy import deepcopy

import numpy as np

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)
from ffn_sim.cortex.cortical_tension import measure_cortical_tension

_MN_PER_M = 1.0e3  # N/m -> mN/m


def run_gate_b(
    *,
    n_filaments: int | None,
    n_nuc_beads: int | None,
    warmup: int,
    sample: int,
    device,
    seed: int,
    softstart: int | None = None,
    constrained: bool = True,
    dt_safety: float = 1.0,
    turnover: bool = False,
    with_fa: bool = True,
) -> dict:
    """Build the FA-adhered cell, settle it, run constrained + record lambda,
    measure the 3 gamma channels. Returns the channel dict + metadata.

    ``softstart`` = number of force-ramped softstart steps in the equilibration
    prelude (default warmup//4). The constrained (M-SHAKE) FA-adhered build is
    seed-sensitive at full scale (round-2: 4/5 seeds blew the BAOAB int32 guard);
    a longer softstart lets stiff FA-clutch/turgor overlaps relax before full-
    force dynamics, stabilizing more seeds."""
    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    manifest["optional_subsystems"]["fa"]["enabled"] = bool(with_fa)  # adhered operating point
    if turnover:
        # Chugh-2017 cortical-actin turnover (tau_half from phase1_h3.yaml cortex.turnover):
        # the one reviewed mechanism that could un-floor the active channel.
        manifest["optional_subsystems"]["turnover"] = {
            "enabled": True, "base_config": "phase1_h3.yaml",
        }
    if n_filaments is not None:
        manifest["cortex_overrides"] = {
            "cortex": {"n_filaments": int(n_filaments), "demo_mode": True}
        }
    if n_nuc_beads is not None:
        manifest["compartments"]["nucleus"]["n_beads"] = int(n_nuc_beads)

    # Equilibration prelude settles the cell onto the substrate (softstart ramps
    # forces so the FA-adhered construction state does not explode); constrained
    # runs the rigid backbone for the Lagrange channel.
    cell = build_baseline_cell(
        manifest=manifest,
        device=device,
        seed=seed,
        constrained=constrained,
        constrained_dt_safety=dt_safety,
        equilibrate=True,
        equilibrate_steps=warmup,
        equilibrate_softstart_steps=(softstart if softstart is not None
                                     else max(100, warmup // 4)),
    )
    handles = cell.extras["handles"]
    sim = cell.simulation
    R_cell = cell.p_cortex.R_cell
    dt_used = float(handles.get("dt_used") or cell.p_cortex.dt_cfl)

    # Capture the M-SHAKE Lagrange multipliers on the constrained action so the
    # rigid channel can be measured (zero overhead otherwise).
    act = handles.get("baoab_action")
    if act is not None and hasattr(act, "record_lambda"):
        act.record_lambda = True

    # Sample run at the settled operating point.
    if sample > 0:
        sim.run(sample)

    # Soft (active) channel filters by bond TYPE: it counts cortical /
    # actomyosin bonds and EXCLUDES the focal-adhesion load path
    # (integrin_ligand, fa_actin_clutch[_b{i}]) so the adhesion force is not
    # mis-read as cortical tension. cortical_bond_types=None selects the robust
    # adhesion-denylist default (B4 refinement; cortex/cortical_tension.py).
    gamma = measure_cortical_tension(
        sim,
        R_cell=R_cell,
        p_enclosed_volume=cell.p_enclosed_volume,
        lambda_accumulator=act,
        dt=dt_used,
        cortical_bond_types=None,
    )
    gamma["_meta"] = {
        "n_filaments": int(cell.p_cortex.n_filaments),
        "n_clutch_bonds": int(handles.get("n_fa_clutch_bonds") or 0),
        "turgor_dP0_Pa": float(cell.p_enclosed_volume.turgor_dP0),
        "dt_used": dt_used,
        "is_full_scale": n_filaments is None,
    }
    return gamma


def _report(g: dict) -> None:
    m = g["_meta"]
    lo, hi = g["band_N_per_m"]
    scale = ("FULL production x40 scale" if m["is_full_scale"]
             else f"SMOKE scale (n_filaments={m['n_filaments']})")
    print("=" * 66, flush=True)
    print(f"H.7 GATE-B — emergent cortical tension  [{scale}]", flush=True)
    print(f"  operating point: FA-adhered ({m['n_clutch_bonds']} clutches), "
          f"turgor Pi_0={m['turgor_dP0_Pa']:.0f} Pa", flush=True)
    print("-" * 66, flush=True)
    soft_ch = g["channels"]["soft"]
    n_excl = int(soft_ch.get("n_bonds_excluded", 0))
    n_cort = int(soft_ch.get("n_bonds", 0))
    print("  CHANNELS (reported SEPARATELY — turgor NOT folded into a total):", flush=True)
    print(f"    gamma_soft     (active MOP)        = {g['gamma_soft']*_MN_PER_M:+.4f} mN/m"
          f"  [cortical bonds={n_cort}, adhesion excluded={n_excl}]", flush=True)
    print(f"    gamma_rigid    (Lagrange M-SHAKE)  = {g['gamma_rigid']*_MN_PER_M:+.4f} mN/m"
          f"  (available={g['channels']['rigid']['available']})", flush=True)
    print(f"    gamma_passive  (turgor Young-Lap.) = {g['gamma_passive']*_MN_PER_M:+.4f} mN/m", flush=True)
    print(f"    gamma_structural (soft+rigid)      = {g['gamma_structural']*_MN_PER_M:+.4f} mN/m", flush=True)
    print(f"  band overlay (Salbreux/Charras/Paluch): "
          f"[{lo*_MN_PER_M:.2f}, {hi*_MN_PER_M:.2f}] mN/m", flush=True)
    print("-" * 66, flush=True)
    if m["n_clutch_bonds"] > 0:
        print(f"  active channel: cortical-bond-TYPE filter ON — {n_excl} adhesion/", flush=True)
        print("    substrate bonds (integrin_ligand, fa_actin_clutch[_b*]) excluded", flush=True)
        print("    from the soft method-of-planes so the FA load path is NOT counted", flush=True)
        print("    as cortical tension (B4 refinement). gamma_soft is the cortical-", flush=True)
        print("    only active tension.", flush=True)
        if g["gamma_soft"] > hi * 5:
            print("  ⚠ gamma_soft still > 5× band after filtering — investigate the", flush=True)
            print("    cortical bonds (NOT an adhesion-bleed artefact anymore).", flush=True)
        print("-" * 66, flush=True)
    # Honesty caveats (direction-review 2026-06-07): the passive channel is an
    # IDENTITY with the B3 turgor setpoint, and the active channel needs the
    # myosin timescale, not the warmup CFL window.
    print("  CAVEATS (read before interpreting):", flush=True)
    print("   - gamma_passive = turgor_dP0*R/2 is fixed BY the B3 band-implied", flush=True)
    print("     turgor setpoint (133 Pa = 2*0.50 mN/m / R) -> a SETPOINT IDENTITY,", flush=True)
    print("     NOT an independent measurement. The EMERGENT cortical tension is", flush=True)
    print("     the active + rigid channels.", flush=True)
    if g["gamma_soft"] < lo:
        print("   - gamma_soft (active) is below band: myosin contraction has not", flush=True)
        print("     developed at the warmup timescale (us-scale CFL steps << the", flush=True)
        print("     seconds-scale myosin/turnover timescale). This is the active-", flush=True)
        print("     GENERATION floor (KU-3.5), not a settled tension. Report PHYSICAL", flush=True)
        print("     time; turnover OFF here. gamma_rigid may carry turgor pre-tension", flush=True)
        print("     + an M-SHAKE shunt (cross-check vs an unconstrained soft run).", flush=True)
    print("  PROVISIONAL. No KU-3.5 conclusion here. Band [0.35,0.65] mN/m is a", flush=True)
    print("  rounded/de-adhered non-MCF7 proxy (gate-contract review pending PI).", flush=True)
    print("=" * 66, flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=None,
                    help="override for a smoke; omit for full production x40 scale")
    ap.add_argument("--n-nuc-beads", type=int, default=None)
    ap.add_argument("--warmup", type=int, default=300, help="equilibration baoab steps")
    ap.add_argument("--sample", type=int, default=200, help="sample steps at operating point")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--softstart", type=int, default=None,
                    help="force-ramped softstart steps (default warmup//4)")
    ap.add_argument("--dt-safety", type=float, default=1.0,
                    help="scale the constrained dt (<1 stabilizes: smaller steps keep the "
                         "FA integrin catch-bond force in the Pereverzev pN range)")
    ap.add_argument("--no-constrained", dest="constrained", action="store_false",
                    help="unconstrained soft-backbone GATE-B (stable across seeds; the "
                         "M-SHAKE-shunt cross-check — gives gamma_soft, not gamma_rigid)")
    ap.set_defaults(constrained=True)
    ap.add_argument("--turnover", action="store_true",
                    help="enable Chugh-2017 actin turnover (does it un-floor the active channel?)")
    ap.add_argument("--no-fa", dest="with_fa", action="store_false",
                    help="free pressurized cortex (no FA adhesion) — clean gamma_rigid without "
                         "the rigid-backbone integrin overload")
    ap.set_defaults(with_fa=True)
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    g = run_gate_b(
        n_filaments=args.n_filaments, n_nuc_beads=args.n_nuc_beads,
        warmup=args.warmup, sample=args.sample, device=dev, seed=args.seed,
        softstart=args.softstart, constrained=args.constrained, dt_safety=args.dt_safety,
        turnover=args.turnover, with_fa=args.with_fa,
    )
    _report(g)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

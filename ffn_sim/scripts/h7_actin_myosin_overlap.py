"""H.7 actin-myosin OVERLAP probe — quantify the Truong-Quang / Chugh architecture regime.

The active-γ literature reconciliation (H7_MYOSIN_DENSITY_LITERATURE_2026-06-09.md) found the
gap is NOT density but actin-myosin OVERLAP / filament-length architecture (Chugh 2017,
Truong Quang 2021): "myosin motors that incompletely overlap with actin generate stress
inefficiently". This probe measures, on the BUILT physiological cell (Gate-B operating point,
loading phase — no contraction needed because overlap is a STRUCTURAL property), how much
actin each engaged minifilament overlaps, to decide whether the model sits in the low-overlap
regime and whether the deficit is CONSTRUCTION-fixable (longer engagement) or geometry-capped.

Reported per engaged minifilament:
  * n_sides_engaged  — 0/1/2 (a bipolar dipole needs both sides on actin to transmit).
  * realized dipole arm ℓ_eff — contour distance between the +side and −side bound actin beads
    (the lever over which the dipole stress acts in the virial). Compared to the minifilament
    backbone (301 nm, the assumed arm) and the actin filament length (3 µm, the ceiling if the
    contraction propagated network-wide).
  * engagement span per side — contour length of actin a side's heads cover (the overlap).
  * same-filament vs different-filament (a real dipole pulls two DIFFERENT filaments together).

No tuning; constants from config. Loading-phase smoke is sufficient (structural statistic).

Usage:
    python -m ffn_sim.scripts.h7_actin_myosin_overlap --n-filaments 200 --warmup 2500 \
        --device cpu --allow-cpu-dev
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)
from ffn_sim.scripts.h7_active_force_budget import _build_settled_cell, _engaged_attach_bonds

_NM = 1e9
_UM = 1e6


def _tag_to_fil_pos(action, tag):
    """Decompose a cortex actin bead tag → (filament, pos) using the updater's map."""
    if action._cortex_filament_starts is not None:
        fil = int(np.searchsorted(action._cortex_filament_starts, tag, side="right") - 1)
        return fil, int(tag - action._cortex_filament_starts[fil])
    nb = action._cortex_beads_per_filament
    return tag // nb, tag % nb


def probe(*, cell, warmup_extra):
    sim = cell.simulation
    if warmup_extra > 0:
        sim.run(warmup_extra)
    action = cell.myosin_action
    p = cell.p_myosin
    H = p.n_heads_per_side
    n_cortex = int(cell.n_cortex_actin)
    ell0 = float(cell.p_cortex.rest_length)        # cortex bead spacing
    backbone = float(p.backbone_length)            # 301 nm (assumed dipole arm)
    fil_len = float(cell.p_cortex.L_filament) if hasattr(cell.p_cortex, "L_filament") else 3.0e-6

    snap = sim.state.get_snapshot()
    bonds, _ = _engaged_attach_bonds(snap, n_cortex)
    # head_local for each engaged bond → which motor, which side.
    # _head_bound_to_actin holds the bead tag per head_local; reconstruct from updater state.
    bound = action._head_bound_to_actin            # (n_heads_total,) bead tag or -1
    bound_fil = action._head_bound_filament
    bound_pos = action._head_bound_bead_pos

    # Group engaged heads by motor.
    per_motor = defaultdict(lambda: {"plus": [], "minus": []})
    n_heads_total = bound.shape[0]
    for h in range(n_heads_total):
        if bound[h] < 0:
            continue
        motor = h // (2 * H)
        side = "plus" if (h % (2 * H)) < H else "minus"
        fil = int(bound_fil[h]) if bound_fil[h] >= 0 else _tag_to_fil_pos(action, int(bound[h]))[0]
        pos = int(bound_pos[h]) if bound_pos[h] >= 0 else _tag_to_fil_pos(action, int(bound[h]))[1]
        per_motor[motor]["plus" if side == "plus" else "minus"].append((fil, pos))

    rows = []
    for motor, sides in per_motor.items():
        nplus, nminus = len(sides["plus"]), len(sides["minus"])
        n_sides = (nplus > 0) + (nminus > 0)
        # engagement span per side (contour covered) in units of ell0
        def span(side_list):
            if not side_list:
                return 0.0
            byfil = defaultdict(list)
            for fil, pos in side_list:
                byfil[fil].append(pos)
            # max contiguous-ish span = (max-min) over the dominant filament
            return max((max(ps) - min(ps)) for ps in byfil.values()) * ell0
        span_plus = span(sides["plus"])
        span_minus = span(sides["minus"])
        # realized dipole arm: contour distance between the + centroid bead and − centroid bead
        ell_eff = None
        diff_fil = None
        if n_sides == 2:
            pf = [f for f, _ in sides["plus"]]; pp = [p for _, p in sides["plus"]]
            mf = [f for f, _ in sides["minus"]]; mp = [p for _, p in sides["minus"]]
            # dominant filament per side
            from statistics import mode
            try:
                pfil, mfil = mode(pf), mode(mf)
            except Exception:
                pfil, mfil = pf[0], mf[0]
            diff_fil = pfil != mfil
            if pfil == mfil:
                ppos = np.mean([p for f, p in sides["plus"] if f == pfil])
                mpos = np.mean([p for f, p in sides["minus"] if f == mfil])
                ell_eff = abs(ppos - mpos) * ell0   # same filament: contour arm
            else:
                # different filaments: arm ≈ backbone (the minifilament spans between them)
                ell_eff = backbone
        rows.append({
            "motor": int(motor), "n_sides_engaged": int(n_sides),
            "n_plus": nplus, "n_minus": nminus,
            "span_plus_nm": span_plus * _NM, "span_minus_nm": span_minus * _NM,
            "ell_eff_nm": (ell_eff * _NM) if ell_eff is not None else None,
            "diff_filament": diff_fil,
        })

    n_eng_motors = len(rows)
    complete = [r for r in rows if r["n_sides_engaged"] == 2]
    diff = [r for r in complete if r["diff_filament"]]
    ell_effs = [r["ell_eff_nm"] for r in complete if r["ell_eff_nm"] is not None]
    spans = [r["span_plus_nm"] for r in rows] + [r["span_minus_nm"] for r in rows]
    spans = [s for s in spans if s > 0]

    result = {
        "operating_point": "suspended/rounded, FA OFF, turgor ON, connected mesh, grip_walk, LOADING",
        "scale": {"n_filaments": int(cell.p_cortex.n_filaments),
                  "n_cortex_actin": n_cortex,
                  "n_motors": int(p.n_motors_per_cell)},
        "geometry_refs_nm": {
            "cortex_bead_spacing_ell0": ell0 * _NM,
            "minifilament_backbone_assumed_arm": backbone * _NM,
            "actin_filament_length_ceiling": fil_len * _NM,
        },
        "engaged_motors": n_eng_motors,
        "complete_bipolar": len(complete),
        "frac_complete": (len(complete) / n_eng_motors) if n_eng_motors else 0.0,
        "frac_complete_diff_filament": (len(diff) / len(complete)) if complete else 0.0,
        "mean_engagement_span_nm": float(np.mean(spans)) if spans else 0.0,
        "median_engagement_span_nm": float(np.median(spans)) if spans else 0.0,
        "mean_ell_eff_nm": float(np.mean(ell_effs)) if ell_effs else None,
        "median_ell_eff_nm": float(np.median(ell_effs)) if ell_effs else None,
        "overlap_ratio_vs_backbone": (float(np.mean(spans)) / (backbone * _NM)) if spans else 0.0,
        "ell_eff_vs_filament_ceiling": (
            (float(np.mean(ell_effs)) / (fil_len * _NM)) if ell_effs else None),
    }
    return result


def _report(r):
    g = r["geometry_refs_nm"]
    print("=" * 72, flush=True)
    print("H.7 ACTIN-MYOSIN OVERLAP PROBE (architecture-regime, loading phase)", flush=True)
    print(f"  scale: n_fil={r['scale']['n_filaments']}, {r['scale']['n_cortex_actin']} beads, "
          f"{r['scale']['n_motors']} motors", flush=True)
    print(f"  refs: ℓ0={g['cortex_bead_spacing_ell0']:.0f}nm  backbone(arm)="
          f"{g['minifilament_backbone_assumed_arm']:.0f}nm  "
          f"actin-fil(ceiling)={g['actin_filament_length_ceiling']:.0f}nm", flush=True)
    print("-" * 72, flush=True)
    print(f"  engaged motors            : {r['engaged_motors']}", flush=True)
    print(f"  complete bipolar (2 sides): {r['complete_bipolar']} "
          f"({r['frac_complete']*100:.1f}%)  diff-filament {r['frac_complete_diff_filament']*100:.0f}%",
          flush=True)
    print(f"  engagement span / side    : mean {r['mean_engagement_span_nm']:.0f}nm, "
          f"median {r['median_engagement_span_nm']:.0f}nm  "
          f"(= {r['overlap_ratio_vs_backbone']:.2f}× backbone)", flush=True)
    if r["mean_ell_eff_nm"] is not None:
        print(f"  realized dipole arm ℓ_eff : mean {r['mean_ell_eff_nm']:.0f}nm, "
              f"median {r['median_ell_eff_nm']:.0f}nm  "
              f"(= {r['ell_eff_vs_filament_ceiling']*100:.1f}% of actin-fil ceiling)", flush=True)
    print("=" * 72, flush=True)
    print("  READING:", flush=True)
    ratio = r["overlap_ratio_vs_backbone"]            # per-side engagement span / backbone
    arm_ceiling = r["ell_eff_vs_filament_ceiling"]    # dipole arm / actin-filament length
    frac_complete = r["frac_complete"]
    # Two distinct architecture limits: per-side OVERLAP (construction-fixable) vs the
    # transmitted ARM (network-transmission, rigid-backbone-gated). The dominant one is the arm.
    if ratio >= 0.5:
        print(f"   • per-side overlap is ADEQUATE ({ratio:.2f}× backbone, ~{ratio:.0f} beads) — so", flush=True)
        print(f"     lengthening engagement (construction) is NOT the missing factor.", flush=True)
    else:
        print(f"   • per-side overlap is LOW ({ratio:.2f}× backbone) — construction could lengthen it.",
              flush=True)
    if arm_ceiling is not None and arm_ceiling < 0.3:
        print(f"   → DOMINANT LIMIT = SHORT TRANSMITTED ARM: realized dipole arm is only "
              f"{arm_ceiling*100:.0f}% of the", flush=True)
        print(f"     actin-filament length — the dipole acts over the minifilament backbone, NOT the", flush=True)
        print(f"     long actin network. With only {frac_complete*100:.0f}% complete bipolar pairs, the", flush=True)
        print(f"     limits are (a) bipolar completion [binding-throughput, prior-refuted as a γ lever]", flush=True)
        print(f"     and (b) NETWORK TRANSMISSION — myosin tension doesn't load the actin network into a", flush=True)
        print(f"     long-arm prestress (rigid M-SHAKE backbone). ⇒ the lever is TRANSMISSION (relax", flush=True)
        print(f"     backbone, integrator-gated), NOT construction overlap.", flush=True)
    else:
        print(f"   → transmitted arm is a substantial fraction of the filament; transmission is engaging.",
              flush=True)
    print("=" * 72, flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=200)
    ap.add_argument("--n-nuc-beads", type=int, default=400)
    ap.add_argument("--warmup", type=int, default=2500)
    ap.add_argument("--softstart", type=int, default=200)
    ap.add_argument("--warmup-extra", type=int, default=8000,
                    help="extra steps after build so the binder engages more heads")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-json", type=str,
                    default="ffn_sim/outputs/h7/production/h7_actin_myosin_overlap.json")
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    n_fil = None if (args.n_filaments is None or args.n_filaments <= 0) else args.n_filaments
    cell = _build_settled_cell(
        n_filaments=n_fil, n_nuc_beads=args.n_nuc_beads,
        warmup=args.warmup, softstart=args.softstart, device=dev, seed=args.seed,
    )
    r = probe(cell=cell, warmup_extra=args.warmup_extra)
    _report(r)
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(r, indent=2))
    print(f"  json → {args.out_json}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

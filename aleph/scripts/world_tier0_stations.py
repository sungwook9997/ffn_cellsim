#!/usr/bin/env python
r"""TIER 0 — does the built cell admit a legal bipolar NMII station at all, and how many?

**Needs NONE of the eight crossbridge PI-GAPs.** Positions, per-strand polarity and a declared reach.
No stiffness, no rate, no force, so nothing here is a magnitude and nothing here can be tuned. The
family module (`world/families/nmii_cortex_crossbridge.py`) asserts in a test that
:func:`station_census`'s signature stays disjoint from its ``BLOCKED_BY`` set, so this tier cannot
quietly acquire a constant.

**Why it runs before the eight constants are ruled on, not after.** ⚠ *The answer may be zero, and
zero is the result.* `STATE.md` (e) 5 records the precedent — SF geometry returned ``N_stations: 0``
for a MODEL reason, with the standing instruction *do not force it*. **If the built cell admits no
legal station, no value of any of the eight constants makes this family contractile**, and knowing
that is worth more than the ruling it precedes.

⚠ **REACH IS A DECLARED AXIS AND ONE POINT IS NOT A RESULT.** The station count rises monotonically
with ``reach_um``, so a single value produces a number whose axis is invisible — and picking the value
that makes the count acceptable is exactly what `CLAUDE.md` forbids. This driver sweeps
:data:`REACH_SWEEP_UM`, **declared here before the run and never narrowed after**, and reports the
whole curve including the values that answer badly.

⚠ **The two failure reasons are different problems and are reported separately**, because collapsing
them sends the fix to the wrong builder: ``n_no_cortex_in_reach`` is a PLACEMENT problem owned by
``build/nmii.py``; ``n_no_antiparallel_partner`` is a POLARITY problem owned by ``build/cortex.py``'s
``polarity_mix`` (PI decision 4).

Sanity Gate:
    * dimensions — µm throughout. The reach ladder is built from two named length scales in the
      build itself (the NMII head offset and the cortex thickness), never typed as a bare number.
    * boundary cases — a uniform-polarity cortex holds no anti-parallel pair at all and the census
      refuses rather than returning zero; that refusal is recorded as a refusal, not as ``0``.
    * conservation — the census's ``station`` array is ``(n_mf, 2)`` with ``-1`` in BOTH columns for
      an unstationed minifilament, so a missing station cannot be read as filament 0. This driver
      re-derives the unstationed count from that array rather than trusting the scalar beside it.
    * CFL/precision — not applicable: no time integration, no force. Distances are float64 host.
    * sign sense — polarity is the whole question; the census counts anti-parallel partners and this
      driver records the cortex's own ``n_plus``/``n_minus`` beside the result so a degenerate mix is
      visible in the artifact rather than inferred.
    * measurement protocol — build once, read positions to the host BETWEEN steps (there are no
      steps), sweep reach over the built cell. One build, many reaches: the cell is identical across
      the curve by construction rather than by assertion.

⚠ **No magnitude, at any reach.** This counts what could ever be paired at build time. The KINETIC
capture radius that decides what binds on a tick is a different quantity and is a PI-GAP.

engine units: µm. Runtime: NVIDIA Warp on CUDA.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

#: The reach ladder, DECLARED BEFORE THE RUN. Multipliers on the two length scales the family's own
#: docstring names as the physically meaningful bracket — the built NMII head offset (0.200 µm) and
#: the cortex thickness. Written as multipliers rather than as absolute µm so the ladder follows the
#: build if either scale changes, instead of silently becoming a different question.
REACH_MULTIPLIERS = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0)


def main(argv: list[str] | None = None) -> int:
    """Build the cell, then sweep the pairing reach over it."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--nmii-radius-um", type=float, default=None,
                    help="⚠ DECLARED A/B AXIS for PI queue item 14, not a tuning knob. Default: the "
                         "shell geometry.footprint() derives (r_cell - tip_clearance - 0.30). That "
                         "0.30 charges the head arm radially, but build/nmii.py's own ey = c x ex is "
                         "TANGENTIAL and costs offset^2/2r radially, so the two files disagree about "
                         "the same arm and nothing owns the relationship. Pass the cortex shell "
                         "radius to measure what that costs. The ACCEPTANCE QUESTION IS FIXED BEFORE "
                         "EITHER RUN and is not a count: does nearest_cortex_node_um.min fall below "
                         "the built head offset at all. Everything else -- support, count, n_bb, "
                         "head offset -- is held identical across the pair, so a difference in the "
                         "curve has exactly one cause.")
    ap.add_argument("--nmii-heads-per-side", type=int, required=True,
                    help="⚠ UNRATIFIED (PI queue item 2). No default: this driver will not invent "
                         "the number that forks the head population 2.8x.")
    args = ap.parse_args(argv)

    import numpy as np
    import warp as wp

    from aleph.scripts.run_provenance import stamp
    from aleph.world.arena import Kind, WorldArena
    from aleph.world.bond import BondCount, SourceClass
    from aleph.world.build import build_all, cortex_shell
    from aleph.world.build.nmii import build_nmii
    from aleph.world.families.nmii_cortex_crossbridge import station_census
    from aleph.world.geometry import footprint

    wp.init()
    dev = wp.get_device(args.device)
    if not dev.is_cuda:                       # the arena raises too; this says WHY before it does
        raise RuntimeError(f"refused: {dev} is not CUDA. There is no CPU path and must never be one.")

    t0 = time.perf_counter()
    arena = WorldArena(capacity={
        Kind.NODE: 12_000_000, Kind.SEGMENT: 12_000_000, Kind.ANGLE3: 12_000_000,
        Kind.ANGLE4: 2_000_000, Kind.FACE: 1_500_000, Kind.STRAND: 200_000,
        Kind.BOND: 1_000_000, Kind.GRID_CELL: 4_000_000}, device=args.device)

    cell = build_all(args.device, arena=arena)
    fp = footprint(7.5)
    H = int(args.nmii_heads_per_side)
    # ⚠ From `build.cortex_shell()`, the ONE derivation. It used to come from `footprint()`, which
    # carried a second formula for the same shell; the two drifted and the motors ended up clear of
    # the cortex. The override stays, as a declared A/B axis, not as a source of truth.
    cortex_r_um, cortex_t_um = cortex_shell()
    nmii_radius_um = cortex_r_um if args.nmii_radius_um is None else float(args.nmii_radius_um)
    HEAD_OFFSET_UM = 0.200                    # the built offset, passed to build_nmii below
    nmii = build_nmii(
        arena,
        count=BondCount(value=0.625, basis="areal", source_class=SourceClass.PI_GAP,
                        scope="MCF7 cortical shell, resting",
                        provenance="Nie 2015 areal density 0.625 um^-2; the density the frozen "
                                   "incumbent already carried, not a new number"),
        support=706.86, n_bb=14, n_heads_per_side=H,
        backbone_length_um=0.301, head_offset_um=HEAD_OFFSET_UM,
        radius_um=nmii_radius_um, thickness_um=cortex_t_um)
    wp.synchronize_device(dev)
    build_s = time.perf_counter() - t0

    # Host copy BETWEEN steps — there are no steps here, and the census is host-side (cKDTree).
    pos = arena.node_arrays["position"].numpy()
    n_lo, n_hi = nmii["claims"]["node"]
    c_lo, c_hi = cell.cortex.nodes.lo, cell.cortex.nodes.hi

    # ⚠ REACH and CONTAINMENT in ONE record. Item 14 happened because they never were: geometry.py
    # asserts against the MEMBRANE and build/nmii.py's excursion guard against its OWN shell, each
    # correct alone, so `n_outside: 0` was true while the motors sat out of reach. This is the
    # built radius of the outermost NMII node -- a geometric fact about the build, not a gate;
    # world_phase1_native.py's assert_inside_membrane remains the gate that refuses.
    nmii_r = np.linalg.norm(pos[n_lo:n_hi], axis=1)
    containment = {
        "max_nmii_node_radius_um": float(nmii_r.max()),
        "r_cell_um": fp.r_cell_um,
        "n_nmii_nodes_outside_membrane": int((nmii_r > fp.r_cell_um).sum()),
        "n_nmii_nodes": int(nmii_r.size),
        "guard_predicted_shell_excursion_um": float(nmii["shell_excursion_um"]),
        "note": ("MEASURED, with build/nmii.py's own predicted excursion beside it. ⚠ THE GUARD'S "
                 "NUMBER IS RIGHT AND ITS QUESTION IS WRONG. It raises when excursion > thickness/2 "
                 "-- a FIT test, 'does a minifilament fit inside the shell' -- but nmii.py:156 draws "
                 "the per-filament radius over the FULL shell, `radius + thickness*(randf - 0.5)`, "
                 "with no inset for the excursion. So a filament drawn at the outer edge puts its "
                 "beads that excursion PAST the shell no matter what the guard returns, and the "
                 "guard cannot fail for the thing it appears to protect. Same family as balance_ok "
                 "and descent_ratio, which geometry.py's own self-check already names. Measured: "
                 "the outermost node sits a CONSTANT ~3.3 nm past radius + thickness/2 at every "
                 "radius on the ladder, against a predicted excursion of 4.2 nm -- so the magnitude "
                 "was never the problem. ⚠ PI queue 14's arithmetic still had the SIGN backwards: "
                 "it reasoned the heads land 3 nm INSIDE the membrane at radius 7.40 and they land "
                 "3 nm outside, which assert_inside_membrane refuses. This driver REPORTS; "
                 "world_phase1_native.py's assert_inside_membrane stays the gate."),
    }

    cortex_thickness_um = float(cortex_t_um)
    scales = {"nmii_head_offset_um": HEAD_OFFSET_UM, "cortex_thickness_um": cortex_thickness_um}
    # One ladder per named scale, so the curve brackets BOTH and neither is privileged.
    ladder = sorted({round(m * v, 6) for m in REACH_MULTIPLIERS for v in scales.values()})

    rows: list[dict] = []
    for reach_um in ladder:
        t = time.perf_counter()
        try:
            cen = station_census(nmii_inventory=nmii, nmii_pos_um=pos[n_lo:n_hi],
                                 cortex=cell.cortex, cortex_pos_um=pos[c_lo:c_hi],
                                 reach_um=float(reach_um))
        except Exception as exc:              # noqa: BLE001 — a refusal IS the datum at this reach
            rows.append({"reach_um": reach_um, "refused": f"{type(exc).__name__}: {exc}",
                         "wall_s": round(time.perf_counter() - t, 3)})
            print(f"[tier0] reach {reach_um:7.4f} um  REFUSED  {type(exc).__name__}", flush=True)
            continue
        # Re-derive the unstationed count from the station array rather than trusting the scalar,
        # because "-1 in both columns" is the representation the census promises and a scalar that
        # disagreed with it would be the interesting failure.
        st = np.asarray(cen["station"], np.int64)
        derived_unstationed = int(((st[:, 0] < 0) & (st[:, 1] < 0)).sum())
        row = {k: (int(v) if isinstance(v, (int, np.integer)) else v)
               for k, v in cen.items() if k != "station" and not hasattr(v, "shape")}
        row.update({"reach_um": reach_um, "unstationed_rederived_from_station_array":
                    derived_unstationed, "wall_s": round(time.perf_counter() - t, 3)})
        rows.append(row)
        print(f"[tier0] reach {reach_um:7.4f} um  stationed "
              f"{int(cen.get('n_stationed', -1)):>5}  unstationed {derived_unstationed:>5}", flush=True)

    record = {
        "schema": "ffn-world-tier0-stations@1",
        "provenance": stamp(__file__),
        "kind": "diagnostic",
        "quantitative_claim_status": "BLOCKED",
        "not_a_claim": (
            "TIER 0: existence and sign only. No stiffness, no rate, no force, no magnitude at any "
            "reach. This counts what COULD ever be paired at build time; the kinetic capture radius "
            "that decides what binds on a tick is a different quantity and is a PI-GAP."),
        "device": str(dev), "build_s": round(build_s, 3),
        # ⚠ The two shells, side by side, because that RELATIONSHIP is the thing no file owns. The
        # membrane containment check and the shell excursion guard are each correct in their own
        # scope, which is how phase1_inside.json reported n_outside: 0 truthfully while the motors
        # sat out of reach of the only population they exist to pull.
        "containment": containment,
        "shells_um": {
            "nmii_radius_used": nmii_radius_um,
            "nmii_radius_from_cortex_shell": cortex_r_um,
            "nmii_radius_overridden": args.nmii_radius_um is not None,
            "nmii_thickness": cortex_t_um,
            "nmii_span": [nmii_radius_um - cortex_t_um / 2.0,
                          nmii_radius_um + cortex_t_um / 2.0],
            "cortex_radius": getattr(cell.cortex, "radius_um", None),
            "cortex_thickness": getattr(cell.cortex, "thickness_um", None),
            "r_cell": fp.r_cell_um,
            "note": ("PI queue 14. support/count/n_bb/head_offset are IDENTICAL across the A/B pair; "
                     "only nmii_radius_used differs, so the curve has one cause. The count is NOT "
                     "the acceptance question -- PI queue 19 is open on the support area and is "
                     "deliberately NOT changed here."),
        },
        "n_heads_per_side": H,
        "n_heads_per_side_status": "UNRATIFIED — PI queue item 2. A declared test point.",
        "reach_axis": {
            "declared_before_run": True,
            "multipliers": list(REACH_MULTIPLIERS),
            "scales_um": scales,
            "ladder_um": ladder,
            "note": ("REACH IS A DECLARED AXIS. The count rises monotonically with reach, so a single "
                     "value hides its own axis. The whole curve is reported, including the reaches "
                     "that answer badly; no value is dropped after seeing the result."),
        },
        # The census returns its own cortex_polarity block INCLUDING `basis`, and that is what is
        # recorded — a mixed population is a claim about structure, so the basis travels with the
        # counts. Recomputing it here would be a second copy of the same fact, which is how the two
        # placement copies drifted on 2026-08-20.
        "cortex_polarity_note": ("per-reach, inside each row — the census supplies it with its basis"),
        # ⚠ All three of these are PI-GAPs and all three DECIDE this run's answer, so the label
        # travels with the value rather than beside it.
        "nmii": {
            "n_minifilaments": int(nmii["n_minifilaments"]),
            "n_bb": {"value": int(nmii["n_bb"]), "source_class": "PI_GAP"},
            "head_offset_um": {"value": HEAD_OFFSET_UM, "source_class": "PI_GAP"},
            "backbone_length_um": {"value": 0.301, "source_class": "PI_GAP"},
            "why_labelled": ("n_bb, head_offset_um and backbone_length_um set the straddle geometry, "
                             "so they set what a reach can reach. A station count quoted without "
                             "them is a number whose inputs are invisible."),
        },
        "failure_reasons_are_not_summable": (
            "n_no_cortex_in_reach is a PLACEMENT problem owned by build/nmii.py; "
            "n_no_antiparallel_partner is a POLARITY problem owned by build/cortex.py's "
            "polarity_mix (PI decision 4). ⚠ They are never added. A combined 'unstationed' count "
            "sends the fix to the wrong builder, which is why the census returns both and no total."),
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1, default=str))
    ok = [r for r in rows if "refused" not in r]
    print(f"[tier0] {len(ok)}/{len(rows)} reaches answered, wrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

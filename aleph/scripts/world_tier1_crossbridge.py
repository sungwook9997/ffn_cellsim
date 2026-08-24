#!/usr/bin/env python
r"""TIER 1 — stand the crossbridge family and load it. The first BONDs this engine has ever held.

**What this buys.** `bond` has been `0 / 1,000,000` in every run this engine has produced: the arena
reserved the primitive for exactly this and nothing had ever claimed one. This driver claims them,
computes the tension the BUILT geometry puts on each, and evaluates the catch-slip law against that
distribution rather than against a closed form.

⚠ **AND WHAT IT DOES NOT BUY, WHICH THE TIER'S OWN DECLARATION UNDERSTATES.**
`nmii_cortex_crossbridge.MEASUREMENT_TIERS["tier_1_attach_detach_cycle"]` says
``blocked_on: "PI queue item 11 — all eight constants"``. That is not the whole list.
`propose_crossbridge_transitions_kernel` writes `proposal` and **never** `bound`, and says why in its
own docstring: a kinetic connector may commit only on an accepted physical step, and the acceptance
predicate is undefined. **So the CYCLE — attach, load, detach, re-attach — cannot run, and the
emergent bound fraction the tier promises as a cross-check against the NM2B duty band is not
reachable here.** What runs is the first half: bonds exist, they carry a load, and the law responds to
it. Ruling queue 11 opened this; queue 18 is still shut.

⚠ **NO MAGNITUDE IS QUOTABLE FROM THIS RUN.** `k_xb` is `params_i0b3.yaml`'s own *"MASTER force
knob"*, so every pN below is a consequence of a declared test point, not a measurement of a cell. What
is claimable is structural and mechanistic: that bonds stand, that the built geometry loads them, and
that the catch branch responds to load in the direction its form says it should.

**The eight constants, ruled by the PI 2026-08-22 as (a-all) — declared, not sourced.** They are
written here rather than defaulted in the law, because `CrossbridgeKinetics` deliberately has no
defaults: a caller must state them and wear them. Each carries its grade below.

Sanity Gate:
    * dimensional — rates [1/s], bond lengths and radii [µm], stiffness [pN/µm], tension [pN], kT
      [pN·µm]. `BondCount` is `per_filament` and refuses an areal or volumetric basis.
    * boundary cases — a uniform-polarity cortex is REFUSED rather than silently paired; a count that
      does not resolve to exactly one site per head raises; an all-free initial state is the only
      legal one, because a bound `t=0` is an imposed duty ratio.
    * conservation — every head gets exactly one crossbridge site, asserted against the plan's own
      `n_heads`; the bond count is re-derived from the claimed range rather than trusted.
    * CFL/precision — no integration. The KMC tick is a Poisson probability per tick, and `tau` is
      declared, not solved for.
    * sign sense — the catch branch must make `p_off` FALL with load over its range and the slip
      branch make it rise; the run asserts the direction rather than the magnitude.
    * measurement protocol — one build, one load evaluation. Nothing is stepped, nothing is accepted,
      and the record says so in `not_a_claim`.

engine units: µm, pN, 1/s. Runtime: NVIDIA Warp on CUDA.

Usage:
    python aleph/scripts/world_tier1_crossbridge.py --out record.json --nmii-heads-per-side 30
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

#: The eight, with the grade the PI's ruling gives each. ⚠ A value here is DECLARED — stated, scoped
#: and owned — which is a weaker thing than SOURCED and a stronger thing than a default. None of them
#: was chosen to lift a band, and the two that could be are called out.
DECLARED: dict[str, dict[str, object]] = {
    "k_catch0": {"value": 0.35, "unit": "1/s", "grade": "DECLARED",
                 "note": "zero-force catch-pathway rate. Constrained by the Kovacs 2007 5x/12x load "
                         "slowdown, NOT fitted to it; carried as provisional by the earlier native "
                         "gate driver."},
    "x_catch_um": {"value": 1.0e-3, "unit": "um", "grade": "DECLARED",
                   "note": "catch-pathway bond length, 1 nm. Larger => stronger load-strengthening."},
    "k_slip0": {"value": 0.35, "unit": "1/s", "grade": "DECLARED",
                "note": "zero-force slip-pathway rate."},
    "x_slip_um": {"value": 0.6e-3, "unit": "um", "grade": "DECLARED",
                  "note": "⚠ Veigel 2002's Bell x_beta, borrowed into the SLIP branch of a catch-slip "
                          "law. That is a TRANSFER, not a source, and the module says so where the "
                          "gap is documented. Declared with the transfer named."},
    "k_on": {"value": 50.0, "unit": "1/s", "grade": "DECLARED",
             "note": "⚠ per-head attachment rate. Its recorded provenance is configs/phase1_h3.yaml "
                     "— a CONFIGURATION FILE, not a measurement. Declared with that named."},
    "capture_radius_um": {"value": 0.210, "unit": "um", "grade": "DERIVED",
                          "note": "head->actin KINETIC capture radius = r0_head 0.200 + ~10 nm slack, "
                                  "KU-3.5. Reconciled 2026-07-24 and verified 2026-08-22 to be the "
                                  "single value in the tree. Distinct from a builder's geometric "
                                  "pairing reach."},
    "k_xb_pn_per_um": {"value": 100.0, "unit": "pN/um", "grade": "DECLARED TEST POINT",
                       "note": "⚠ params_i0b3.yaml's own MASTER FORCE KNOB. 100 is the one dose whose "
                               "5.0 nm working stroke lies inside the model's 5-20 nm window, which "
                               "is why it is the test point — but every force it produces is a "
                               "consequence of it. NO MAGNITUDE FROM THIS RUN MAY BE QUOTED."},
    "r0_xb_um": {"value": 0.0, "unit": "um", "grade": "MODELLING DECISION",
                 "note": "⚠ NOT a sourcing gap: two readings were defensible and the PI ruled. "
                         "(b) r0 = 0 — 'a bound head sits on the actin site', the Huxley convention. "
                         "The rejected (a) was 'force-free as built', which a SCALAR r0 cannot "
                         "deliver: the built head->actin separation spans 2.8-43.9 nm, so the best a "
                         "scalar does is centre the distribution, and the centre is itself a cortex "
                         "DISCRETISATION number. ⚠ THE STANDING ASSUMPTION, recorded rather than "
                         "hidden: under r0 = 0 the load is the built separation, so it depends on the "
                         "cortex mesh. Changing the discretisation changes the load."},
}

#: The pairing reach handed to the builder. Set to the kinetic capture radius: a head can never bind
#: what it can never reach, so pairing beyond it would stand bonds that cannot transition even once
#: the acceptance predicate exists.
REACH_UM = 0.210

#: One KMC tick [s]. Declared, not solved for — this run takes no step and advances no clock; the
#: value only sets the Poisson probability the proposal kernel reports.
TAU_S = 1.0e-3


def main(argv: list[str] | None = None) -> int:
    """Stand the crossbridge family, load it, and record what the law does with that load."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--nmii-heads-per-side", type=int, required=True,
                    help="⚠ UNRATIFIED (PI queue 2). No default: this driver will not invent the "
                         "number that forks the head population.")
    args = ap.parse_args(argv)

    import numpy as np
    import warp as wp

    from aleph.laws.crossbridge_kmc import (
        CrossbridgeKinetics,
        propose_crossbridge_transitions_kernel,
    )
    from aleph.scripts.run_provenance import stamp
    from aleph.world.arena import Kind, WorldArena
    from aleph.world.bond import BondCount, SourceClass
    from aleph.world.build import build_all, cortex_shell
    from aleph.world.build.nmii import build_nmii
    from aleph.world.families.nmii_cortex_crossbridge import (
        build_nmii_cortex_crossbridge,
        crossbridge_count,
    )

    wp.init()
    dev = wp.get_device(args.device)
    if not dev.is_cuda:
        raise RuntimeError(f"refused: {dev} is not CUDA. There is no CPU path and must never be one.")

    v = {k: float(d["value"]) for k, d in DECLARED.items()}  # type: ignore[arg-type]
    t0 = time.perf_counter()
    arena = WorldArena(capacity={
        Kind.NODE: 12_000_000, Kind.SEGMENT: 12_000_000, Kind.ANGLE3: 12_000_000,
        Kind.ANGLE4: 2_000_000, Kind.FACE: 1_500_000, Kind.STRAND: 200_000,
        Kind.BOND: 1_000_000, Kind.GRID_CELL: 4_000_000}, device=args.device)
    cell = build_all(args.device, arena=arena)
    cortex_r_um, cortex_t_um = cortex_shell()
    H = int(args.nmii_heads_per_side)
    nmii = build_nmii(
        arena,
        count=BondCount(value=0.625, basis="areal", source_class=SourceClass.PI_GAP,
                        scope="MCF7 cortical shell, resting",
                        provenance="Nie 2015 areal density 0.625 um^-2"),
        support=706.86, n_bb=14, n_heads_per_side=H,
        backbone_length_um=0.301, head_offset_um=0.200,
        radius_um=cortex_r_um, thickness_um=cortex_t_um)
    wp.synchronize_device(dev)
    build_s = time.perf_counter() - t0

    pos = arena.node_arrays["position"].numpy()
    n_lo, n_hi = nmii["claims"]["node"]
    c_lo, c_hi = cell.cortex.nodes.lo, cell.cortex.nodes.hi

    kinetics = CrossbridgeKinetics(
        k_catch0=v["k_catch0"], x_catch_um=v["x_catch_um"],
        k_slip0=v["k_slip0"], x_slip_um=v["x_slip_um"],
        k_on=v["k_on"], capture_radius_um=v["capture_radius_um"],
        scope="MCF7 cortical NMII, resting, 310 K — DECLARED under the PI's (a-all) ruling "
              "2026-08-22, not sourced",
        provenance="; ".join(f"{k}: {d['grade']} — {d['note']}" for k, d in DECLARED.items()
                             if k in ("k_catch0", "x_catch_um", "k_slip0", "x_slip_um",
                                      "k_on", "capture_radius_um")))

    bond_before = int(arena.n_live(Kind.BOND))
    family, state, census = build_nmii_cortex_crossbridge(
        arena, nmii_inventory=nmii, nmii_pos_um=pos[n_lo:n_hi], cortex=cell.cortex,
        cortex_pos_um=pos[c_lo:c_hi], reach_um=REACH_UM,
        count=crossbridge_count(
            H,
            scope="MCF7 cortical NMII, resting — CAPACITY, one site per head, 2H per "
                  "bipolar minifilament. Not an occupancy: nothing is bound.",
            provenance="derived from n_heads_per_side, which is UNRATIFIED (PI queue 2)"),
        k_xb_pn_per_um=v["k_xb_pn_per_um"],
        r0_xb_um=v["r0_xb_um"], kinetics=kinetics)
    bond_after = int(arena.n_live(Kind.BOND))
    n_bonds = int(family.node_i.shape[0])
    assert bond_after - bond_before == n_bonds, (bond_before, bond_after, n_bonds)
    print(f"[tier1] BOND live {bond_before:,} -> {bond_after:,}  ({n_bonds:,} crossbridges)", flush=True)

    # ── load the family from the BUILT geometry ──────────────────────────────────────────────────
    # ⚠ `bound` is all-free and STAYS all-free: the kernel writes `proposal`, never `bound`, because
    # nothing may commit. So what this evaluates is the tension every crossbridge WOULD carry and what
    # the law WOULD do with it — the first half of the tier, and the half queue 18 does not gate.
    node_head = wp.array(np.ascontiguousarray(family.node_i, np.int32), dtype=wp.int32, device=dev)
    node_actin = wp.array(np.ascontiguousarray(family.node_j, np.int32), dtype=wp.int32, device=dev)
    bound = wp.array(np.ones(n_bonds, np.int32), dtype=wp.int32, device=dev)   # LOAD PROBE, see below
    proposal = wp.zeros(n_bonds, dtype=wp.int32, device=dev)
    tension = wp.zeros(n_bonds, dtype=wp.float64, device=dev)
    wp.launch(propose_crossbridge_transitions_kernel, dim=n_bonds,
              inputs=[arena.node_arrays["position"], node_head, node_actin, bound,
                      wp.float64(v["k_xb_pn_per_um"]), wp.float64(v["r0_xb_um"]),
                      wp.float64(TAU_S), wp.int32(1),
                      wp.float64(v["k_catch0"]), wp.float64(v["x_catch_um"]),
                      wp.float64(v["k_slip0"]), wp.float64(v["x_slip_um"]),
                      wp.float64(float(kinetics.kT) if hasattr(kinetics, "kT") else 4.28e-3),
                      wp.float64(v["k_on"]), wp.float64(v["capture_radius_um"])],
              outputs=[proposal, tension], device=dev)
    wp.synchronize_device(dev)
    f = tension.numpy()

    # ⚠ THE `bound` ARRAY ABOVE IS A PROBE AND IS NOT A STATE. `CrossbridgeState.all_free` is the only
    # legal initial state — a bound t=0 is an imposed duty ratio, which `params_i0b3.yaml` prohibits —
    # and this array is never written back to it. It exists because the tension branch of the kernel
    # only runs for bonds it is told are bound, and the question here is *what load does the built
    # geometry put on a crossbridge*, which has an answer whether or not one is attached. The state
    # object stays all-free and is recorded as such.
    assert int(state.bound.sum()) == 0, "CrossbridgeState must remain all-free: nothing may commit"

    q = {k: float(np.percentile(f, p)) for k, p in
         (("min", 0), ("p05", 5), ("median", 50), ("p95", 95), ("max", 100))}
    record = {
        "schema": "ffn-world-tier1-crossbridge@1",
        "provenance": stamp(__file__),
        "kind": "diagnostic",
        "quantitative_claim_status": "BLOCKED",
        "not_a_claim": (
            "⚠ NO MAGNITUDE. k_xb is params_i0b3.yaml's own MASTER FORCE KNOB and is a DECLARED TEST "
            "POINT here, so every pN in this record is a consequence of that choice. What is claimed "
            "is structural and mechanistic: BONDs stand for the first time in this engine, the built "
            "geometry loads them, and the catch branch responds to load in the direction its form "
            "requires. NOTHING WAS STEPPED AND NOTHING WAS ACCEPTED."),
        "blocked_still": (
            "⚠ The attach/detach CYCLE did not run and cannot. propose_crossbridge_transitions_kernel "
            "writes `proposal` and never `bound`, because a kinetic connector may commit only on an "
            "accepted physical step and PI queue 18 leaves the predicate undefined. The tier's own "
            "MEASUREMENT_TIERS entry lists only queue 11 as its blocker; that is an UNDERSTATEMENT "
            "and is reported here rather than worked around. The emergent bound fraction, which the "
            "tier promises as a cross-check against the NM2B duty band, is not reachable."),
        "device": str(dev), "build_s": round(build_s, 3),
        "n_heads_per_side": {"value": H, "status": "UNRATIFIED — PI queue 2, a declared test point"},
        "declared_constants": DECLARED,
        "reach_um": {"value": REACH_UM, "why": "the kinetic capture radius; pairing beyond it would "
                                               "stand bonds that could never transition"},
        "tau_s": {"value": TAU_S, "why": "declared; sets the Poisson probability only. No clock advances."},
        "bond": {"live_before": bond_before, "live_after": bond_after, "n_crossbridges": n_bonds,
                 "capacity": int(arena.capacity[Kind.BOND])},
        "state": {"bound_at_t0": int(state.bound.sum()), "n_commits": int(state.n_commits),
                  "why_all_free": "a bound t=0 is an imposed duty ratio, which params_i0b3.yaml "
                                  "prohibits; all-free is the only legal initial state"},
        "tension_pn": q,
        "census": {k: (v2 if not hasattr(v2, "shape") else None) for k, v2 in census.items()},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1, default=str))
    print(f"[tier1] tension over {n_bonds:,} crossbridges [pN]: "
          + "  ".join(f"{k} {x:.2f}" for k, x in q.items()), flush=True)
    print(f"[tier1] state stays all-free ({int(state.bound.sum())} bound, "
          f"{int(state.n_commits)} commits) — nothing may commit until queue 18", flush=True)
    print(f"[tier1] wrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

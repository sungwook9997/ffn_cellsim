"""PHASE 1 — stand all eight populations in one arena, natively, and measure what standing costs.

WHAT THIS IS, AND WHY IT LIVES HERE.  Sessions A1 and A2 own five and three builder modules under
``aleph/world/build/``; neither owns an assembly driver, so this file is where the populations are
stood up together and the standing cost is read.

⚠ **The reason recorded here until 2026-08-22 no longer holds, and it was the reason for the file.**
It said ``world/`` could not import the NVML accounting this measurement needs, so a driver in
``aleph/scripts/`` had to call it instead. The probe moved to :mod:`aleph.world.gpu_memory` on
2026-08-22 and the layer objection is gone; what remains is the ORDERING constraint, which is real —
NVML accounting must be requested before the CUDA context exists, so the request happens here, above
``wp.init()``, and not inside a package that a caller may already have imported.

⚠ THE COUNTS BELOW ARE NOT A CLAIM.  Five of the eight populations have NO sourced per-cell structure
count — `PARALLEL_AUDIT_PARAMS_DETECTORS_2026-08-20.md` §3 — and every one of them is passed here as an
``explicit`` :class:`BondCount` carrying ``source_class=PI_GAP``, so the refusal travels into the
artifact rather than being lost in a default. This driver measures what a world of THAT SIZE costs to
stand. It does not assert that a cell has that many microtubules.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — lengths [µm], memory [B], wall times [s]. No force, no physics.
  * boundary — a population whose count cannot resolve raises rather than building a smaller one; the
    arena refuses a claim exceeding capacity with both numbers.
  * conservation/invariant — ``assert_partitioned`` after every claim: the ranges must tile ``[0,
    n_live)`` exactly, no gap and no overlap. Node arrays are span-checked for disjointness.
  * CFL/precision — nothing integrates. float64 positions throughout.
  * sign sense — not applicable; nothing is accumulated.
  * measurement protocol — every wall time brackets ``wp.synchronize_device`` so an asynchronous launch
    cannot be reported as a fast build. Memory is read from the driver, and ``exact_peak_gpu_bytes``
    comes from NVML per-process accounting enabled BEFORE the CUDA context exists, or is reported null
    with the reason.

engine units: length µm.  Runtime: NVIDIA Warp on CUDA only.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--nmii-radius-um", type=float, default=None,
                    help="⚠ DECLARED A/B AXIS for PI queue item 14; same flag and same meaning as "
                         "world_tier0_stations.py. Default: the shell geometry.footprint() derives. "
                         "TIER 0 answers whether the motors can REACH the actin; this driver answers "
                         "whether the same placement is still inside the MEMBRANE, and the two "
                         "questions have been answered by different files with nobody owning the "
                         "relationship -- which is the whole of item 14. Do not move one without "
                         "re-running the other.")
    ap.add_argument("--nmii-heads-per-side", type=int, default=None,
                    help="NMII motor heads per minifilament side. NO DEFAULT, deliberately: this is "
                         "decision-queue item 9 and it is UNRATIFIED. Both candidates fail the citation "
                         "audit -- 10 (AFINES, from a draft KB row that declares its own citations "
                         "absent) and 30 (from KB-3.18, verified as a CORTEX COMPOSITION claim whose "
                         "citations are reviews, not head-count measurements). Omit the flag and the "
                         "NMII population is skipped with that reason recorded. Pass a value and it is "
                         "recorded as an UNRATIFIED test point, never as a sourced count.")
    ap.add_argument("--cortex-seg-um", type=float, default=0.05,
                    help="sourced band is 50-100 nm; 0.05 is its FLOOR and 0.10 its ceiling, and the "
                         "choice halves or doubles the dominant population. Recorded per run.")
    args = ap.parse_args()

    # NVML accounting must be on before the CUDA context exists — hence before importing anything
    # that calls wp.init(). This is why the driver, not the world package, owns the measurement.
    peak_note = ""
    accounting = None
    try:
        from aleph.world import gpu_memory as gma

        accounting = gma
        peak_note = "NVML per-process accounting requested before context creation"
    except Exception as exc:  # noqa: BLE001 — a missing probe is reported, never silently skipped
        peak_note = f"accounting unavailable: {type(exc).__name__}: {exc}"

    import warp as wp

    from aleph.world.arena import Kind, WorldArena
    from aleph.world.bond import BondCount, SourceClass
    from aleph.world.geometry import assert_inside_membrane, footprint
    from aleph.world.build import build_all, cortex_shell
    from aleph.world.populations import build_remaining_populations
    from aleph.world.build.cytosol import CYTOSOL_AXES, PI_DX_UM, build_cytosol
    from aleph.world.build.nmii import build_nmii

    wp.init()
    dev = wp.get_device(args.device)
    if not dev.is_cuda:
        raise RuntimeError(f"PHASE 1 is native-only; {args.device!r} resolved to {dev}, not CUDA")

    def gap(n: int, what: str) -> BondCount:
        """An explicit count that says, in the artifact, that no source exists for it."""
        return BondCount(value=float(n), basis="explicit", source_class=SourceClass.PI_GAP,
                         scope="NO VALUE EXISTS — per-cell structure count is absent from the "
                               f"Contract-Graph for {what}",
                         provenance="PI-GAP 2026-08-20; placeholder for a size measurement, not a claim")

    def derived(n: int, what: str) -> BondCount:
        """A count that FOLLOWS from geometry, so it is not a gap and must not be labelled one."""
        return BondCount(value=float(n), basis="explicit", source_class=SourceClass.DERIVED,
                         scope="MCF7 basal footprint AS BUILT — R_cell 7.5 um, substrate z -7.0 um",
                         provenance=f"DERIVED from the contact disc: {what} = "
                                    "floor(2 r_footprint / pitch) + 1, the pitch held at its "
                                    "physiological value so the sourced quantity stays sourced")

    t_all = time.perf_counter()
    stages: list[dict] = []

    def stage(name, fn):
        wp.synchronize_device(dev)
        t = time.perf_counter()
        out = fn()
        wp.synchronize_device(dev)
        dt = time.perf_counter() - t
        live = {str(k): v for k, v in arena._live.items()} if "arena" in dir() else {}
        stages.append({"population": name, "wall_s": round(dt, 4), "live_after": live})
        print(f"[phase1] {name:<24} {dt:8.3f} s", flush=True)
        return out

    cap = {Kind.NODE: 12_000_000, Kind.SEGMENT: 12_000_000,
           Kind.ANGLE3: 12_000_000, Kind.ANGLE4: 2_000_000, Kind.FACE: 1_500_000,
           Kind.STRAND: 200_000, Kind.BOND: 1_000_000,
           # GRID_CELL, added 2026-08-20 with the PI's Kind approval. Headroom for the whole declared
           # dx band down to 0.1 um (153^3 = 3,581,577); at the PI's 0.25 test point the claim is
           # 63^3 = 250,047. Capacity is not a prediction — the arena refuses an over-claim with both
           # numbers, which is the check that catches a dx typo before it becomes a device allocation.
           Kind.GRID_CELL: 4_000_000}
    arena = WorldArena(capacity=cap, device=args.device)

    cell = stage("A1: cortex+membrane+envelope",
                 lambda: build_all(args.device, arena=arena))

    # ⚠ EVERY placement below is DERIVED from this. PI rule 2026-08-20: every node lies inside the
    # built membrane. The call sites here previously carried origin/length/pitch/width as unlabelled
    # constants, and rendering then measuring on 2026-08-20 found they put 95.5% of the stress fibres
    # and 92.7% of the lamellipodium OUTSIDE the cell -- while assert_partitioned reported success,
    # because partition is a property of ID ranges and says nothing about where a node is.
    fp = footprint(7.5)
    # ⚠ Every length below is DECLARED here rather than defaulted, because A2's builders refuse a
    # missing one — "reach_R_um has no default and must be declared" is what stopped the first run of
    # this driver, and that refusal is the discipline working. The values are the ones A2's own
    # self-checks use; none of them is sourced, and the record says so in `not_a_claim`.
    seg = args.cortex_seg_um
    # Every population below cortex/membrane/envelope comes from ONE place — `world/populations.py`.
    # It was two places until 2026-08-21 (this driver and the renderer), and the two copies drifted
    # from the cell they were placing: 95.5% of the stress fibres ended up outside the membrane while
    # assert_partitioned reported success. The geometry moved to world/geometry.py then; the calls
    # move now, because PHASE 4 needed the same cell and a third copy is how a third divergence starts.
    built: dict[str, object] = build_remaining_populations(arena, fp, seg_um=seg, stage=stage)

    # ── the cytosol, the first FIELD in this arena ────────────────────────────────────────────────
    # Not a node population: `Kind.GRID_CELL` landed 2026-08-20 because `arena.py`'s own argument for
    # closing a seven-member set already cited it by name. dx is the PI's declared TEST POINT and is
    # passed explicitly — `build_cytosol` has no default, deliberately, because cost goes as dx^-5 and
    # a run must never inherit a resolution it did not state.
    #
    # ⚠ D_p is `KB-3.B3.2` (verified, 40-60), and that row says of itself that it is a cross-cell-type
    # anchor, NOT breast-specific. For MCF7 it is a proxy whose scope IS recorded — which is the whole
    # difference between this and the rows that failed today. dt_phys is not physiological at all: it
    # is this driver's outer step, and it is here only because the reported subcycle count is what
    # makes dx cost anything.
    DT_PHYS_S = 0.05
    D_P_UM2_S = 50.0
    cytosol = stage("E: cytosol field",
          lambda: build_cytosol(arena, dx_um=PI_DX_UM, poroelastic_diffusion_um2_s=D_P_UM2_S,
                                dt_phys_s=DT_PHYS_S, membrane=cell.membrane, envelope=cell.envelope,
                                dx_provenance=str(CYTOSOL_AXES["dx_um"]["source"])))

    # ── NMII minifilaments — head-resolved, and the head COUNT is not ours to choose ───────────────
    # G's builder refuses without n_heads_per_side on purpose. This driver refuses to invent one too:
    # omitted, the population is skipped and the record says why. The 442 minifilaments themselves are
    # NOT the open question -- Nie 2015's 0.625/um^2 x 706.86 um^2 = 442 is what the frozen incumbent
    # already carried. What is open is how many heads hang off each one, and it forks the head
    # population 2.8x.
    nmii = None
    if args.nmii_heads_per_side is not None:
        H = int(args.nmii_heads_per_side)
        nmii = built["nmii"] = stage(f"G: nmii (H={H}, UNRATIFIED)",
              lambda: build_nmii(arena,
                    # areal, not explicit: the density is the datum and the arena does the arithmetic,
                    # so 442 is derived in the record instead of typed into this driver.
                    count=BondCount(value=0.625, basis="areal", source_class=SourceClass.PI_GAP,
                        scope="MCF7 cortical shell, resting",
                        provenance="Nie 2015 areal density 0.625 um^-2; the density the frozen "
                                   "incumbent already carried, not a new number"),
                    support=706.86, n_bb=14, n_heads_per_side=H,
                    backbone_length_um=0.301, head_offset_um=0.200,
                    radius_um=(cortex_shell()[0] if args.nmii_radius_um is None
                               else float(args.nmii_radius_um)),
                    thickness_um=cortex_shell()[1]))

    wp.synchronize_device(dev)
    total_s = time.perf_counter() - t_all
    arena.assert_partitioned()

    # ⚠ The PI rule, checked against BUILT positions rather than argued from the arithmetic above.
    # `assert_partitioned` has already passed at this point and it does NOT cover this: it says the ID
    # ranges are disjoint, not that any node is inside the cell. Both were true on 2026-08-20 while
    # 95.5% of the stress fibres sat outside the membrane.
    placement = assert_inside_membrane(arena.node_arrays["position"].numpy(),
                                       {"cortex": cell.cortex, "membrane": cell.membrane,
                                        "nuclear_envelope": cell.envelope, **built}, fp.r_cell_um)
    census = arena.census()

    peak = None
    if accounting is not None:
        # The probe moved into `world/` on 2026-08-22, so this driver no longer reaches into the
        # frozen port source for it. The comment that stood here recorded the opposite arrangement
        # as a PI decision of 2026-08-20; the ruling that the arena is the engine and the rest is
        # archived supersedes it.
        #
        # ⚠ AND THE TWO LINES BELOW USED TO READ ATTRIBUTES THAT DO NOT EXIST. `peak_bytes` and
        # `max_memory_usage_bytes` are on no version of `NvmlAccountingResult` — the field is
        # `exact_peak_bytes` — so both `getattr` calls returned their default and `peak` was
        # unconditionally None on every run this driver has ever done, including on a host where
        # accounting was enabled. It read as a probe that had been tried and come up empty, because
        # the note beside it printed the repr of a result object that did carry the number.
        try:
            res = accounting.query_process_peak(args.device)
            peak = res.exact_peak_bytes if res.exact else None
            peak_note = f"gpu_memory.query_process_peak -> {res!r}"
        except Exception as exc:  # noqa: BLE001 — a probe that cannot read says so
            peak_note = f"query_process_peak raised {type(exc).__name__}: {exc}"

    from aleph.scripts.run_provenance import stamp

    record = {
        "schema": "ffn-world-phase1-native@1",
        # ⚠ See world_phase4_native.py: the GPU tree is hand-deployed, so the record must stamp the
        # bytes it imported rather than a commit that may be hours stale.
        "provenance": stamp(__file__),
        "kind": "diagnostic",
        "quantitative_claim_status": "BLOCKED",
        "device": str(dev),
        "cortex_seg_um": args.cortex_seg_um,
        "stages": stages,
        "total_build_s": round(total_s, 3),
        "census": census,
        "partitioned": True,
        "placement_envelope": fp.as_record(),
        "placement_check": placement,
        "placement_rule": "PI 2026-08-20 — every node of every population lies inside the built "
                          "membrane. Enforced after the build by assert_inside_membrane, which "
                          "raises and names the population; assert_partitioned does not cover it.",
        "nmii": ({
            "n_heads_per_side": int(args.nmii_heads_per_side),
            "head_count_status": "UNRATIFIED — decision-queue item 9. NEITHER candidate traces to a "
                                 "primary measurement passing the citation audit: Billington2013 and "
                                 "StamAlbertsGardelMunro2015 are both source_audit CHECK. This run "
                                 "stands the geometry at a test point and claims nothing about it.",
            "n_minifilaments": int(nmii["n_minifilaments"]),
            "radius_um": (cortex_shell()[0] if args.nmii_radius_um is None
                          else float(args.nmii_radius_um)),
            "radius_um_overridden": args.nmii_radius_um is not None,
            "radius_um_from_cortex_shell": cortex_shell()[0],
            "minifilament_count_basis": "Nie 2015 areal 0.625 um^-2 x 706.86 um^2 shell — DERIVED by "
                                        "BondCount.resolve, not typed in",
            "physics_claim": "NONE — no bond, no crossbridge, no force. n_bonds is 0 by construction.",
        } if nmii is not None else {
            "skipped": True,
            "why": "--nmii-heads-per-side not given. The head count is unratified (decision-queue item "
                   "9) and inventing one to make the driver run is the defect class this engine spent "
                   "2026-08-20 cataloguing. The population is absent by decision, not by oversight.",
        }),
        "cytosol": {
            "dx_um": PI_DX_UM,
            "dx_basis": "PI 2026-08-20 test point inside the sourced band; NOT a derivation",
            "poroelastic_diffusion_um2_s": D_P_UM2_S,
            "poroelastic_diffusion_scope": "KB-3.B3.2 verified 40-60, self-declared CROSS-CELL-TYPE "
                                           "anchor and NOT breast-specific — a scope-recorded proxy "
                                           "for MCF7, not a sourced MCF7 value",
            "dt_phys_s": DT_PHYS_S,
            "dt_phys_basis": "this driver's outer step, not a physiological quantity",
            "shape": list(cytosol.shape),
            "n_cells": int(cytosol.shape[0] * cytosol.shape[1] * cytosol.shape[2]),
            "origin_um": list(cytosol.origin_um),
            "dx_provenance": cytosol.dx_provenance,
            "field_bytes": int(cytosol.field_bytes),
            # The self-contained check: the masked FLUID volume against V_cell - V_nuc taken from the
            # SAME built surfaces the classifier ran against. Deliberately NOT against KB-3.21's
            # V_cyto = 4200 um^3, which is a R = 10.0 um sphere — 3.47x this cell's sourced R = 7.50.
            # That row is verified and out of scope, which is what UNRATIFIED_PROXY names.
            "fluid_volume_um3": float(cytosol.fluid_volume_um3),
            "reference_volume_um3": float(cytosol.reference_volume_um3),
            "volume_error": float(cytosol.volume_error),
            "physics_claim": "NONE — seven fields allocated and zero; no Biot step, no pressure, no flux",
        },
        "exact_peak_gpu_bytes": peak,
        "exact_peak_note": peak_note,
        "not_a_claim": (
            "five of the eight populations carry PI_GAP counts — no sourced per-cell structure count "
            "exists for microtubule, intermediate filament, filopodium, lamellipodium or stress fibre. "
            "This measures what a world of this SIZE costs to stand; it does not say a cell has these "
            "many of anything."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=2, default=str))
    print(f"\n[phase1] total {total_s:.2f} s   partitioned OK   wrote {args.out}", flush=True)
    print(f"[phase1] exact_peak_gpu_bytes = {peak}  ({peak_note})", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""How long one accepted physical step of the composed native cell costs, and how far that is from real time.

Implements: KU-0.0

WHY THIS EXISTS.  The goal "a full cell that steps at a rate comparable to real time" has no meaning
without a number, and this repository has been burned by timing that was not comparable: a solver
stopped early is arbitrarily fast, so **cost is only cost at a stated acceptance**.  ``run-record@2``
already encodes that rule (``comparable: false`` when the run did not meet its own predicate); this
driver reports the same way and additionally refuses to average over the warm-up, where Warp is
compiling kernels and first-touching device memory.

THE NUMBER.  ``RTF = dt_phys / wall_per_step`` — the **real-time factor**.  RTF = 1 means the cell is
simulated as fast as it lives; RTF = 0.01 means 100 s of wall per second of cell.  It is reported
beside its reciprocal (``wall seconds per second of cell life``) because that is the quantity a run
plan is written in.

SCOPE IS PART OF THE NUMBER.  An RTF measured on 3 of 11 components is not the cell's RTF, and quoting
it as one is exactly the promotion this charter forbids.  Every record therefore carries the bound
component and connector census, and ``fraction_of_target_components``.  The target here is the full
cell EXCLUDING ``ecm``, ``extracellular_medium`` and ``world_boundary`` — 11 components and the 31
connectors internal to them.

Sanity Gate (per CLAUDE.md, before first execution):
  * dimensional — ``dt_phys`` [s] of cell time per step, ``wall`` [s]; RTF is dimensionless.
  * boundary — zero timed steps RAISES rather than reporting an undefined mean; a single step reports
    its own value with no spread and says so.
  * conservation — timing only; no physics is read from this driver and no gate is scored by it.
  * CFL/precision — ``dt_phys`` is the caller's; this driver never chooses one.
  * sign sense — n/a.
  * measurement protocol — CUDA only, RAISES otherwise; the device is synchronized before and after
    every timed region, because an unsynchronized Warp launch returns immediately and would report a
    step cost of nearly zero.
"""

from __future__ import annotations

import argparse
import json
import statistics
import time
from pathlib import Path

import numpy as np
import warp as wp

#: Components deliberately outside this goal's "full cell".
EXCLUDED_FROM_FULL_CELL = ("ecm", "extracellular_medium", "world_boundary")


def _target_census() -> tuple[frozenset[str], tuple[str, ...]]:
    """The 11 in-cell components and the connectors internal to them, read from the contract."""
    from aleph.engine.contracts import reference_cell_architecture

    architecture = reference_cell_architecture()
    inside = frozenset(c.name for c in architecture.components) - frozenset(EXCLUDED_FROM_FULL_CELL)
    edges = tuple(
        sorted(
            c.name
            for c in architecture.connectors
            if c.component_a in inside and c.component_b in inside
        )
    )
    return inside, edges


def _timed_steps(world, *, dt_phys: float, n_steps: int, warmup: int, device: str) -> dict:
    """Run ``warmup + n_steps`` accepted steps, timing only the last ``n_steps``.

    The warm-up is excluded because the first steps pay kernel compilation and first-touch page
    faults; averaging them in would report a cost no later step ever pays again.
    """
    if n_steps <= 0:
        raise ValueError("n_steps must be positive; a mean over zero steps is not a measurement")
    accepted = wp.array(np.ones(1, np.int32), dtype=wp.int32, device=device)
    handle = wp.get_device(device)

    for _ in range(max(int(warmup), 0)):
        world.step_once(dt_phys=dt_phys, accepted_d=accepted)
    wp.synchronize_device(handle)

    per_step: list[float] = []
    for _ in range(int(n_steps)):
        start = time.perf_counter()
        world.step_once(dt_phys=dt_phys, accepted_d=accepted)
        wp.synchronize_device(handle)          # an unsynchronized launch returns almost instantly
        per_step.append(time.perf_counter() - start)

    mean = statistics.fmean(per_step)
    return {
        "n_timed_steps": int(n_steps),
        "n_warmup_steps": int(max(int(warmup), 0)),
        "dt_phys_s": float(dt_phys),
        "wall_per_step_s": mean,
        "wall_per_step_stdev_s": statistics.stdev(per_step) if len(per_step) > 1 else None,
        "wall_per_step_min_s": min(per_step),
        "wall_per_step_max_s": max(per_step),
        "real_time_factor": float(dt_phys / mean) if mean > 0 else None,
        "wall_seconds_per_second_of_cell": float(mean / dt_phys) if dt_phys > 0 else None,
        "per_step_s": per_step,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--k-axial", type=float, required=True, help="sf_arc axial stiffness [pN/um] (PI-GAP)")
    ap.add_argument("--cortex-filaments", type=int, default=70686)
    ap.add_argument("--n-fibers-ecm", type=int, default=200)
    ap.add_argument("--ecm-box-um", type=float, default=10.0)
    ap.add_argument("--dt", type=float, default=0.01, help="physical timestep [s]")
    ap.add_argument("--steps", type=int, default=20, help="TIMED accepted steps")
    ap.add_argument("--warmup", type=int, default=3, help="untimed steps first (kernel compile/first touch)")
    ap.add_argument("--relax-sf", action="store_true", help="sf_arc advances its own positions (subcycled)")
    ap.add_argument("--sf-motor", action="store_true",
                    help="bind nmii_sf_motor as well, with head exclusivity DERIVED by minifilament "
                         "proximity. This is the one power port the interface residual reported UNWIRED, "
                         "and its absence is why sf_arc carries no load.")
    ap.add_argument("--interior-column", action="store_true",
                    help="CARD-5 (PI-approved 2026-08-11): bind membrane / cortex / cytosol / nucleus and "
                         "the three fluid connectors into the SAME world, with the incumbent inner solve "
                         "owning their mechanics under D2 `omit=`. This is where the cost lives.")
    ap.add_argument("--n-inner", type=int, default=40,
                    help="incumbent inner iterations for the interior column (PI verdict 2026-08-11: 40)")
    ap.add_argument("--osmotic-pi0-pa", type=float, default=40.0,
                    help="van 't Hoff pi_0 [Pa]; PI Option-A constant (GAP), stated not defaulted silently")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--profile-kernels", action="store_true",
                    help="wrap ONE extra step in wp.TIMING_KERNEL and report per-kernel device time, "
                         "ranked. Optimising before this is guessing; the repo's own 2026-07-28g profile "
                         "is what established that filament physics is a minority of a step.")
    ap.add_argument("--label", type=str, default="baseline")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    wp.init()
    if not wp.get_device().is_cuda:
        raise SystemExit("I0-A: requires CUDA; a step time measured on CPU is not a result.")
    device = str(wp.get_device())

    from aleph.components.ecm.mikado_topology import MikadoInitConfig, MikadoTopologyBuilder
    from aleph.components.incumbent.assemble import (
        NMII_BACKBONE_LP_DIAGNOSTIC_UM, NMII_K_XB_TEST, CellConfig, build_cell,
    )
    from aleph.engine.composed_native import build_native_composed_cell_world
    from aleph.engine.ecm_mechanics import COLLAGEN_MATERIAL_KEY
    from aleph.engine.ecm_world import BoundaryAnchorMode, ECMWorldSettings
    from aleph.engine.microtubule_rig import SF_COMPONENT
    from aleph.engine.overdamped_relax import build_overdamped_step
    from aleph.engine.sf_mechanics import build_sf_mechanics_topology
    from aleph.engine.sf_population import build_sf_arc_population
    from aleph.laws.ecm_library import get_spec
    from aleph.scripts.ac_gate_b_composed_native import _build_nmii_cortex_connector
    from aleph.scripts.ac_gate_b_cortex_motor_native import _build_actuator, _build_port, _params

    build_start = time.perf_counter()
    k_xb = float(NMII_K_XB_TEST)
    cell = build_cell(CellConfig(
        n_filaments=int(args.cortex_filaments), with_myosin=True, overlap_free_cortex=True,
        with_membrane=True, with_nucleus=True, with_pressure=bool(args.interior_column),
        membrane_subdivisions=6, nucleus_subdivisions=3, erm_radial_pairing=True,
        resting_bound_myosin_fraction=None, nmii_straddle_placement=True,
        nmii_backbone_lp_um=NMII_BACKBONE_LP_DIAGNOSTIC_UM,
        nmii_backbone_lp_source="STEP_TIMING_DIAGNOSTIC_FIXTURE_NOT_PRODUCTION_PI_GAP",
        seed=int(args.seed),
    ))
    actuator = _build_actuator(cell, device, k_xb)
    port = _build_port(cell, device)
    connector = _build_nmii_cortex_connector(
        actuator_state=actuator, port=port, params=_params(False, k_xb), device=device,
        max_segment_length_um=float(cell.srest_d.numpy().max()),
    )
    sf_pop = build_sf_arc_population(n_ventral=8, n_dorsal=4, n_arc=4, n_cap=4, n_per_fiber=9)
    sf_pop.assert_partitioned()
    spec = get_spec(COLLAGEN_MATERIAL_KEY)
    ecm_topology = MikadoTopologyBuilder(MikadoInitConfig(
        box_lo_um=(0.0, 0.0, 0.0), box_hi_um=(args.ecm_box_um,) * 3, n_fibers=int(args.n_fibers_ecm),
        fiber_length_um=spec.fiber_len_um, target_segment_um=spec.seg_um,
        crosslink_capture_um=spec.xl_contact_um, pin_faces=("z_lo",), pin_margin_um=1.0,
        rng_seed=int(args.seed), max_refinement_level=0,
    ), device=device).initialize()

    # CARD-5. The interior column is where the cost lives — CLAUDE.md's own profile says filament physics is
    # a small minority of a step and the expense is the fluid grid and the moving membrane boundary. Binding
    # it is therefore both the composition goal and the thing that makes the RTF number honest.
    interior_column = interior_solve = None
    if args.interior_column:
        from aleph.components.incumbent.driver import make_inner_solve
        from aleph.engine.interior_column_slice import build_native_interior_column_parts

        interior_column = build_native_interior_column_parts(
            cell, device=device, osmotic_pi0_pa=float(args.osmotic_pi0_pa))
        # D2 `omit=`: the engine's NMII actuator now supplies the myosin channel on its own arrays, so the
        # incumbent must NOT launch it as well. Any other channel stays with the incumbent — this is a
        # relocation, not a removal, and the default is bit-identical.
        _inner = make_inner_solve(cell, n_inner=int(args.n_inner), reshape_every=20,
                                  inner_solver="explicit", omit=("nmii_minifilament_force",))
        # The composed `solve()` slot is zero-argument by contract (three instruments call it bare), so the
        # timestep is CLOSED OVER here rather than threaded through a changed signature.
        _dt = float(args.dt)
        interior_solve = lambda: _inner(_dt)       # noqa: E731 — a zero-arg closure is the solve contract

    world = build_native_composed_cell_world(
        device=device, base_seed=int(args.seed), sf_population=sf_pop,
        interior_column=interior_column, interior_solve=interior_solve,
        sf_k_axial_pn_per_um=float(args.k_axial), ecm_topology=ecm_topology,
        ecm_settings=ECMWorldSettings(BoundaryAnchorMode.FAR_FIELD_DIRICHLET), ecm_spec=spec,
        nmii_actuator=actuator, cortex_port=port, nmii_cortex_connector=connector,
        cortex_capacity=int(args.cortex_filaments),
        # nmii_sf_motor: the params are the SAME sourced hand card the cortex edge uses; what differs
        # is only the target port and the head set, so no new constant enters with this edge.
        sf_motor_params=_params(False, k_xb) if args.sf_motor else None,
        cortex_positions_um=cell.pos_d.numpy()[: int(cell.n_actin)] if args.sf_motor else None,
        # The immersed control volume for the two new solid<->cytosol edges, read from the SAME
        # fluid grid the cortex edge uses. A different rule per edge would make three bound
        # transfers physically inconsistent while all three still reported as wired.
        cytosol_grid_dx_um=float(cell.grid.dx) if args.interior_column else None,
    )
    build_seconds = time.perf_counter() - build_start
    sf_motor = dict(getattr(world, "sf_motor_report", {"bound": False}))
    if sf_motor.get("bound"):
        print(f"[sf-motor] partition: {sf_motor['sf_heads']} heads -> sf_arc, "
              f"{sf_motor['cortex_heads']} -> cortex (ties to cortex: "
              f"{sf_motor['tie_broken_to_cortex']})")

    relax = {"attached": False}
    if args.relax_sf:
        topology = build_sf_mechanics_topology(sf_pop, k_axial_pn_per_um=float(args.k_axial))
        step = build_overdamped_step(
            topology, n_filaments=max(int(getattr(sf_pop, "n_filaments", 0) or 1), 1),
            borrowed_pairs=world.borrowed_pairs_for(SF_COMPONENT),
        )
        world.attach_relaxer(SF_COMPONENT, step, world.component_owners[SF_COMPONENT].position_d)
        relax = {"attached": True, "dt_max_s": step.dt_max_s,
                 "subcycles_per_outer_step": world.subcycles_for(float(args.dt))}

    timing = _timed_steps(world, dt_phys=float(args.dt), n_steps=int(args.steps),
                          warmup=int(args.warmup), device=device)

    kernel_profile: list[dict] = []
    if args.profile_kernels:
        accepted = wp.array(np.ones(1, np.int32), dtype=wp.int32, device=device)
        wp.timing_begin(wp.TIMING_KERNEL)
        try:
            world.step_once(dt_phys=float(args.dt), accepted_d=accepted)
        finally:
            launched = wp.timing_end(synchronize=True)
        rolled: dict[str, dict] = {}
        for record in launched:
            entry = rolled.setdefault(record.name, {"kernel": record.name, "launches": 0, "device_ms": 0.0})
            entry["launches"] += 1
            entry["device_ms"] += float(record.elapsed)
        total_ms = sum(e["device_ms"] for e in rolled.values()) or 1.0
        kernel_profile = sorted(rolled.values(), key=lambda e: -e["device_ms"])
        for entry in kernel_profile:
            entry["share"] = entry["device_ms"] / total_ms
        print(f"\n[profile] {len(kernel_profile)} distinct kernels, {sum(e['launches'] for e in kernel_profile)} "
              f"launches, {total_ms:.2f} ms device time in ONE step")
        print(f"[profile] {'kernel':52s} {'launches':>9s} {'ms':>9s} {'share':>7s}")
        cumulative = 0.0
        for entry in kernel_profile[:22]:
            cumulative += entry["share"]
            print(f"[profile] {entry['kernel'][:52]:52s} {entry['launches']:9d} "
                  f"{entry['device_ms']:9.2f} {entry['share']:6.1%}  (cum {cumulative:.0%})")

    target_components, target_edges = _target_census()
    bound_components = set(world.component_owners) | {"cortex"}
    bound_connectors = set(world.connector_runtimes)
    in_scope_bound = sorted(bound_components & target_components)

    record = {
        "driver": "ac_fullcell_step_timing",
        "device": device,
        "label": args.label,
        "evidence_class": "TIMING at a FORCE-ACCEPTED step",
        "quantitative_claim_status": "SCOPED",
        "evidence_basis": (
            "wall-clock per accepted step of the composed native cell, device-synchronized, warm-up "
            "excluded. The steps are FORCE-ACCEPTED: this measures the cost of the step the runtime "
            "currently takes, NOT the cost of a converged step, and the two are not the same number."
        ),
        "timing": {**timing, "build_seconds": build_seconds, "comparable": False,
                   "not_comparable_reason": "force-accepted steps; cost is only cost at a stated "
                                            "acceptance, and no acceptance predicate gated these"},
        "scope": {
            "target_full_cell_components": sorted(target_components),
            "n_target_components": len(target_components),
            "n_target_in_cell_connectors": len(target_edges),
            "excluded_by_goal": list(EXCLUDED_FROM_FULL_CELL),
            "bound_components_in_scope": in_scope_bound,
            "fraction_of_target_components": len(in_scope_bound) / len(target_components),
            "bound_connectors": sorted(bound_connectors),
            "n_bound_connectors": len(bound_connectors),
            "bound_connectors_in_target": sorted(bound_connectors & set(target_edges)),
        },
        "census": {"cortex_filaments": int(args.cortex_filaments),
                   "n_actin_nodes": int(cell.n_actin), "n_total_nodes": int(cell.n_total),
                   "n_heads": int(actuator.n_heads),
                   "fraction_of_native_cortex": float(int(args.cortex_filaments) / 70686.0)},
        "sf_relax": relax,
        "nmii_sf_motor": sf_motor,
        "kernel_profile": kernel_profile,
        "config": {"k_axial_pn_per_um": float(args.k_axial), "k_xb_pn_per_um": k_xb,
                   "n_fibers_ecm": int(args.n_fibers_ecm), "seed": int(args.seed)},
        "parameter_provenance": {"k_axial_pn_per_um": "PI_GAP", "k_xb_pn_per_um": "PI_GAP"},
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2), encoding="utf-8")

    rtf = timing["real_time_factor"]
    print(f"\n[timing] {args.label}  device={device}")
    print(f"[timing] bound in-scope components {len(in_scope_bound)}/{len(target_components)}: {in_scope_bound}")
    print(f"[timing] bound connectors in target: {len(record['scope']['bound_connectors_in_target'])}"
          f"/{len(target_edges)}")
    print(f"[timing] wall/step {timing['wall_per_step_s']:.4f} s  "
          f"(min {timing['wall_per_step_min_s']:.4f}  max {timing['wall_per_step_max_s']:.4f})")
    print(f"[timing] dt={args.dt}s  ->  REAL-TIME FACTOR {rtf:.5f}  "
          f"= {timing['wall_seconds_per_second_of_cell']:.1f} wall-s per second of cell life")
    print(f"[timing] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

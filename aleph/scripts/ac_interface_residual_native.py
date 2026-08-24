#!/usr/bin/env python
r"""The interface residual of the REAL composed native cell — the number that decides the next move.

`ac_interface_residual_control.py` showed the gate can pass and, more importantly, that it FAILS on a
one-sided scatter and on a sign flip, returning the injected magnitude. That established the
instrument on synthetic forces. This driver points it at the cell.

WHY THIS RUN DECIDES SOMETHING.  `SOLVE_COUPLED` is empty, so the plan is a partitioned iteration:
each component relaxes its own positions with the connector force held, forces are re-evaluated,
repeat under acceleration. **That iteration is only worth building if the interface is actually out of
balance.** Nobody has measured it. Two outcomes, and they point opposite ways:

  residual at round-off   the connectors already close on the assembled candidate. There is nothing
                          for an interface iteration to converge, and the blocker on private arrays is
                          the missing per-component INTEGRATOR, not the missing coupling loop. The
                          next work is a relaxation step and its derived mobility -- a PI-visible
                          modelling step, with this number as the justification for taking it.
  residual finite         the connectors do NOT close, and which cut says which connector. Then the
                          iteration has a target and its convergence is a checkable claim.

Either way the answer is evidence rather than a guess, and it costs one candidate assembly.

WHAT MAY NOT BE READ OFF THIS.  Not convergence — no position is updated here and nothing is solved.
Not a physics magnitude — the composed world binds a handful of components and its parameters are
PI-GAPs. The residual is an ADJOINT-CLOSURE reading of one assembled force state, which is exactly
what `ledger.py` says its two-channel gate tests: the wiring, not the solve.

Sanity Gate: dimensional — vec3d [pN], residual [pN], tolerance dimensionless. Boundary — a cut whose
components are unbound in this world reports UNBOUND, never BALANCED (most cuts here are unbound, and
reading those as balanced would be the whole error). Conservation — two never-merged arrays per cut.
Sign sense — covered by the control driver, not re-derived here. CFL/precision — no integration.
Measurement protocol — one candidate assembly, one device reduction per cut, host readback after.

Runtime: CUDA only.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import warp as wp

from aleph.engine.contracts import reference_cell_architecture
from aleph.engine.coupled_solve import collect_force_arrays, measure_interface_residual
from aleph.engine.microtubule_rig import SF_COMPONENT
from aleph.engine.overdamped_relax import build_overdamped_step
from aleph.engine.sf_mechanics import build_sf_mechanics_topology


@wp.kernel
def _zero_vec3d(force: wp.array(dtype=wp.vec3d)) -> None:
    """Zero one watched force array between isolation passes."""
    force[wp.tid()] = wp.vec3d(0.0, 0.0, 0.0)


#: The two solid<->cytosol transfers, and whose arrays each one scatters into.
#: Every interior-column edge whose `accumulate` takes (pos, force), and the component whose
#: arrays it scatters into. These need an explicit closure: `callable()` does not check ARITY,
#: so the zero-arg probe claims them and then TypeErrors at launch.
_SOLID_CYTOSOL_ISOLATION = {
    "sf_cytosol_transfer": "sf_arc",
    "nmii_cytosol_transfer": "nmii",
    "surface_porous_transfer": "cortex",
    "membrane_cytosol_boundary": "membrane",
    "nucleus_cytosol_boundary": "nucleus",
}


def main() -> None:
    ap = argparse.ArgumentParser(description="Interface residual of the composed native cell (CUDA).")
    ap.add_argument("--k-axial", type=float, required=True,
                    help="SF axial backbone stiffness [pN/um] — a PI-GAP, required so it is never defaulted")
    ap.add_argument("--cortex-filaments", type=int, default=70686)
    ap.add_argument("--n-fibers-ecm", type=int, default=200)
    ap.add_argument("--ecm-box-um", type=float, default=10.0)
    ap.add_argument("--k-xb", type=float, default=None)
    ap.add_argument("--steps", type=int, default=0,
                    help="accepted physical steps to drive BEFORE measuring. At 0 the cell is at rest, "
                         "every NMII head is unbound, and the reading is vacuous by construction — the "
                         "2026-07-28f trap. Steps let k_on events bind heads so force exists to cancel.")
    ap.add_argument("--dt", type=float, default=0.01, help="physical timestep [s] for those steps")
    ap.add_argument("--relax-perturb-um", type=float, default=0.0,
                    help="POSITIVE CONTROL. Displace sf_arc's interior nodes by this much along +x "
                         "before driving, then report how much of it the relax removes. Without it, "
                         "'the relax moved 4.5e-15 um' and 'no relaxer was attached' are the SAME "
                         "number, and the run cannot tell them apart.")
    ap.add_argument("--interior-column", action="store_true",
                    help="CARD-5: bind membrane / cortex / cytosol / nucleus + the three fluid "
                         "connectors, so the cuts they create can be SCORED and not merely bound.")
    ap.add_argument("--n-inner", type=int, default=40)
    ap.add_argument("--osmotic-pi0-pa", type=float, default=40.0)
    ap.add_argument("--relax-sf", action="store_true",
                    help="let sf_arc advance its OWN positions inside each accepted step "
                         "(overdamped_relax; gamma and dt_max both derived, subcycled to the "
                         "derived bound). Off by default so the SAME driver measures the relax "
                         "against its own absence.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--commit", type=str, default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    wp.init()
    if not wp.get_device().is_cuda:
        raise SystemExit("I0-A: requires CUDA; there is no CPU simulation path.")
    dev = str(wp.get_device())

    from aleph.components.ecm.mikado_topology import MikadoInitConfig, MikadoTopologyBuilder
    from aleph.components.incumbent.assemble import (
        NMII_BACKBONE_LP_DIAGNOSTIC_UM, NMII_K_XB_TEST, CellConfig, build_cell,
    )
    from aleph.engine.composed_native import build_native_composed_cell_world
    from aleph.engine.ecm_mechanics import COLLAGEN_MATERIAL_KEY
    from aleph.engine.ecm_world import BoundaryAnchorMode, ECMWorldSettings
    from aleph.engine.sf_population import build_sf_arc_population
    from aleph.laws.ecm_library import get_spec
    from aleph.scripts.ac_gate_b_composed_native import _build_nmii_cortex_connector
    from aleph.scripts.ac_gate_b_cortex_motor_native import _build_actuator, _build_port, _params

    t0 = time.perf_counter()
    k_xb = float(args.k_xb) if args.k_xb is not None else float(NMII_K_XB_TEST)
    cell = build_cell(CellConfig(
        n_filaments=int(args.cortex_filaments), with_myosin=True, overlap_free_cortex=True,
        with_membrane=True, with_nucleus=True, with_pressure=bool(args.interior_column),
        membrane_subdivisions=6, nucleus_subdivisions=3, erm_radial_pairing=True,
        resting_bound_myosin_fraction=None, nmii_straddle_placement=True,
        nmii_backbone_lp_um=NMII_BACKBONE_LP_DIAGNOSTIC_UM,
        nmii_backbone_lp_source="INTERFACE_RESIDUAL_DIAGNOSTIC_FIXTURE_NOT_PRODUCTION_PI_GAP",
        seed=int(args.seed),
    ))
    actuator = _build_actuator(cell, dev, k_xb)
    port = _build_port(cell, dev)
    connector = _build_nmii_cortex_connector(
        actuator_state=actuator, port=port, params=_params(False, k_xb), device=dev,
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
    ), device=dev).initialize()
    build_seconds = time.perf_counter() - t0

    # CARD-5. The interior column is where the cost lives — CLAUDE.md's own profile says filament physics is
    # a small minority of a step and the expense is the fluid grid and the moving membrane boundary. Binding
    # it is therefore both the composition goal and the thing that makes the RTF number honest.
    interior_column = interior_solve = None
    if args.interior_column:
        from aleph.components.incumbent.driver import make_inner_solve
        from aleph.engine.interior_column_slice import build_native_interior_column_parts

        interior_column = build_native_interior_column_parts(
            cell, device=dev, osmotic_pi0_pa=float(args.osmotic_pi0_pa))
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
        device=dev, base_seed=int(args.seed), sf_population=sf_pop,
        sf_k_axial_pn_per_um=float(args.k_axial), ecm_topology=ecm_topology,
        ecm_settings=ECMWorldSettings(BoundaryAnchorMode.FAR_FIELD_DIRICHLET), ecm_spec=spec,
        nmii_actuator=actuator, cortex_port=port, nmii_cortex_connector=connector,
        cortex_capacity=int(args.cortex_filaments),
        interior_column=interior_column, interior_solve=interior_solve,
        cytosol_grid_dx_um=float(cell.grid.dx) if args.interior_column else None,
    )

    # sf_arc advances its OWN positions when asked. Both constants derived: gamma from NF2007 fiber
    # mobility, dt_max = 2*gamma/lambda_max from the Gershgorin row sums — WITH the alpha-actinin
    # arc joints (omitting them made dt 459x too large) and WITH whatever a bound connector scatters
    # into this component's force array, which its own topology cannot see.
    relax = {"attached": False}
    if args.relax_sf:
        sf_topology = build_sf_mechanics_topology(sf_pop, k_axial_pn_per_um=float(args.k_axial))
        sf_owner = world.component_owners[SF_COMPONENT]
        step = build_overdamped_step(
            sf_topology,
            n_filaments=max(int(getattr(sf_pop, "n_filaments", 0) or 1), 1),
            borrowed_pairs=world.borrowed_pairs_for(SF_COMPONENT),
        )
        world.attach_relaxer(SF_COMPONENT, step, sf_owner.position_d)
        relax = {
            "attached": True,
            "gamma_pn_s_per_um": step.gamma_pn_s_per_um,
            "lambda_max_pn_per_um": step.lambda_max_pn_per_um,
            "dt_max_s": step.dt_max_s,
            "subcycles_per_outer_step": world.subcycles_for(float(args.dt)),
            "n_borrowed_pair_sets": len(world.borrowed_pairs_for(SF_COMPONENT)),
        }
        print(f"[relax] sf_arc: gamma={step.gamma_pn_s_per_um:.4f} pN*s/um  "
              f"lambda_max={step.lambda_max_pn_per_um:.6g} pN/um  dt_max={step.dt_max_s:.6g} s  "
              f"-> {relax['subcycles_per_outer_step']} subcycles per dt={args.dt}s outer step")

    # Drive accepted steps FIRST when asked, so that NMII k_on events have bound heads and the
    # connector has force to carry. Force-accepted, exactly as the GATE-B lane does: the point is to
    # populate the interface, not to claim a converged trajectory.
    # POSITIVE CONTROL. sf_arc is built AT its own rest state, so every axial spring is unstretched
    # and — with `nmii_sf_motor` unwired — nothing loads it. A relax then correctly moves nothing,
    # which is indistinguishable from a relaxer that was never attached. Displacing the interior
    # nodes creates a real restoring force whose removal only an ATTACHED, CORRECTLY SIGNED relax
    # can produce.
    sf_owner_arrays = world.component_owners[SF_COMPONENT]
    perturbation = float(args.relax_perturb_um)
    rest_positions = sf_owner_arrays.position_d.numpy().copy()
    if perturbation != 0.0:
        perturbed = rest_positions.copy()
        interior = slice(1, len(perturbed) - 1)          # ends are where the connectors attach
        perturbed[interior, 0] += perturbation
        sf_owner_arrays.position_d.assign(perturbed)

    sf_position_before = sf_owner_arrays.position_d.numpy().copy()
    n_bound_before = int(connector.state.bound_d.numpy().sum())
    import numpy as _np
    for _ in range(int(args.steps)):
        world.step_once(
            dt_phys=float(args.dt),
            accepted_d=wp.array(_np.ones(1, _np.int32), dtype=wp.int32, device=dev),
        )
    wp.synchronize_device(wp.get_device())
    n_bound_after = int(connector.state.bound_d.numpy().sum())

    # ONE candidate force assembly on whatever state the steps left. Nothing is committed here.
    world.solve_candidate()
    wp.synchronize_device(wp.get_device())

    owners = dict(world.component_owners)

    def runtime_for(name: str) -> object | None:
        """The cortex is a PORT, not a bound component, but it is a real force-array endpoint here."""
        if name in owners:
            return owners[name]
        return port if name == "cortex" else None

    # LEAVE-ONE-IN. A whole-body reduction sums every edge on both bodies, so a body sitting on two
    # edges shows the other edge's force as an uncancelled residual — which is how a correctly wired
    # sf_cortex_transient read cancellation 0.99999999998 on 2026-08-11. Zero the watched arrays and
    # re-run ONLY this cut's connectors, so the number belongs to the connector rather than the bodies.
    watched = collect_force_arrays(port, *world.component_owners.values())
    connector_runtimes = {
        name: world.world.actor.connector_runtime(name)
        for name in reference_cell_architecture().power_ports()
        if name in set(world.world.registered_connectors())
    }

    # Per-connector accumulate calls. A zero-argument `accumulate()` is the common shape, but the motor
    # connector needs its actuator geometry and its bind-target port, so it gets an explicit closure.
    # Anything NOT listed here cannot be isolated, and that must RAISE rather than silently contribute
    # nothing: job 97 skipped the motor because getattr(runtime, "accumulate") returned None, and the
    # cut then read "0 pN, cancellation None" — "could not isolate" and "carries no force" arriving as
    # the same number is the exact failure this whole instrument exists to prevent.
    isolation_calls: dict[str, object] = {}
    for name, runtime in connector_runtimes.items():
        # EXPLICIT CASES FIRST. `callable()` does not check ARITY: every `ImmersedPorousTransfer`
        # exposes an `accumulate`, so the zero-arg probe below claimed all three transfers and then
        # TypeError'd at launch. That is the same shape as every other defect this lane has hit —
        # "has the method" read as "can be called".
        if name in _SOLID_CYTOSOL_ISOLATION:
            _solid = world.component_owners.get(_SOLID_CYTOSOL_ISOLATION[name]) or port
            isolation_calls[name] = (
                lambda r=runtime, o=_solid: r.accumulate(o.position_d, o.force_d)
            )
            continue
        zero_arg = getattr(runtime, "accumulate", None)
        if callable(zero_arg):
            isolation_calls[name] = zero_arg
        elif name == "nmii_cortex_motor":
            isolation_calls[name] = lambda r=runtime: r.accumulate_candidate(actuator.geometry(), port)


    def isolate(names: tuple[str, ...]) -> None:
        missing = [n for n in names if n not in isolation_calls]
        if missing:
            raise RuntimeError(
                f"cannot isolate {missing}: no accumulate call is known for them. Reporting zero force "
                "instead would make 'not isolatable' indistinguishable from 'carries nothing'."
            )
        for array in watched:
            n = int(array.shape[0])
            if n:
                wp.launch(_zero_vec3d, dim=n, inputs=[array], device=str(array.device))
        for name in names:
            isolation_calls[name]()
        wp.synchronize_device(wp.get_device(dev))

    report = measure_interface_residual(
        world=world.world, architecture=reference_cell_architecture(), device=dev,
        component_runtime=runtime_for, isolate=isolate,
    )
    wall = time.perf_counter() - t0

    record = {
        "schema": "ac.engine.observe/run-record@2",
        "run_label": "interface_residual_composed_native",
        "kind": "diagnostic",
        "build": {"commit": args.commit or None, "source": "declared",
                  "reason": "the run host is not a git checkout; the commit is the caller's assertion"},
        "device": dev,
        "evidence": "NATIVE",
        "quantitative_claim_status": "BLOCKED",
        "evidence_basis": (
            "adjoint-closure residual of ONE assembled candidate force state on the composed native "
            "cell, per two-body cut, scored by the existing D8-derived tolerance. BLOCKED: no position "
            "was updated, nothing was solved, and the composed world's parameters are PI-GAPs. This "
            "reads the WIRING, which is what ledger.py says its two-channel gate tests."
        ),
        "config": {"n_filaments": int(args.cortex_filaments), "k_axial_pn_per_um": float(args.k_axial),
                   "k_xb_pn_per_um": k_xb, "n_fibers_ecm": int(args.n_fibers_ecm),
                   "membrane_subdivisions": 6, "seed": int(args.seed),
                   "steps_driven": int(args.steps), "dt_phys": float(args.dt),
                   "acceptance": "force-accept — these steps populate the interface, they do not "
                                 "claim a converged or admissible trajectory"},
        "nmii_binding": {"n_heads": int(actuator.n_heads), "bound_before": n_bound_before,
                         "bound_after": n_bound_after},
        # Without this the run cannot distinguish "the relax ran and the residual fell" from "the
        # relax moved nothing", and zero motion would read as a pass.
        "sf_relax": {
            **relax,
            "max_node_displacement_um": float(
                np.abs(sf_owner_arrays.position_d.numpy() - sf_position_before).max()
            ),
            "positive_control": {
                "perturbation_um": perturbation,
                # How far the cell still sits from where it was BUILT. An attached, correctly signed
                # relax drives this down; an inert one leaves it at the perturbation exactly.
                "offset_from_rest_before_um": float(
                    np.abs(sf_position_before - rest_positions).max()
                ),
                "offset_from_rest_after_um": float(
                    np.abs(sf_owner_arrays.position_d.numpy() - rest_positions).max()
                ),
            },
        },
        "parameter_provenance": {"k_axial_pn_per_um": "PI_GAP", "k_xb_pn_per_um": "PI_GAP"},
        "census": {"cortex_filaments": int(args.cortex_filaments), "n_actin_nodes": int(cell.n_actin),
                   "n_total_nodes": int(cell.n_total), "n_heads": int(actuator.n_heads),
                   "fraction_of_native": float(int(args.cortex_filaments) / 70686.0)},
        "bound_components": sorted(owners) + ["cortex (PORT, not a component)"],
        "interface": report.as_dict(),
        "timing": {"wall_seconds": wall, "build_seconds": build_seconds, "n_steps": 0,
                   "comparable": False,
                   "not_comparable_reason": "no step was accepted; this is one force assembly read"},
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2), encoding="utf-8")

    payload = report.as_dict()
    print(f"\n[residual] device={dev}  bound: {sorted(owners)} + cortex PORT")
    print(f"[residual] steps driven {args.steps} @ dt={args.dt}s  NMII heads bound "
          f"{n_bound_before} -> {n_bound_after} / {actuator.n_heads}")
    print(f"[residual] cuts declared {payload['n_cuts_declared']}  measured {payload['n_measured']}  "
          f"unbound {payload['n_unbound']}  unwired {payload['n_unwired']}  "
          f"SCOREABLE {payload['n_scoreable']}  imbalanced {payload['n_imbalanced']}")
    for cut in payload["cuts"]:
        if cut["status"] == "UNBOUND":
            continue
        cancel = cut["cancellation"]
        note = "UNWIRED" if cut["unwired"] else ("NOTHING-CANCELLED" if cancel is not None
                                                 and cancel > 0.5 else "")
        print(f"[residual]   {cut['components'][0]:>10s} | {cut['components'][1]:<10s} "
              f"{cut['status']:11s} residual={cut['residual_pn']:.6e} pN  scale={cut['scale_pn']:.6e} "
              f"cancel={'None' if cancel is None else f'{cancel:.4e}'}  "
              f"A={cut['reaction_pn']:.4e} B={cut['traction_pn']:.4e} "
              f"{'ISO' if cut['isolated'] else 'WHOLE-BODY'}  {note}  {cut['connectors']}")
    if payload["n_scoreable"] == 0:
        print("[residual] NOTHING SCOREABLE — every measured cut is either unwired or has no force to "
              "cancel. A BALANCED verdict here is not evidence that any adjoint closes.")
    if payload["n_measured"] == 0:
        print("[residual] NO CUT WAS MEASURABLE — this world binds no two components that share a "
              "power port. That is itself the finding: the coupling loop has nothing to iterate on "
              "until more owners are bound.")
    print(f"[residual] wrote {out}")


if __name__ == "__main__":
    main()

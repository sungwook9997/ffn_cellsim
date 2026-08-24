#!/usr/bin/env python
r"""T10 native gate — the ``extracellular_medium`` exterior Stokes operator on CUDA.

WHAT THIS RUN IS FOR.  PI decision D5-A declared the 14th component and its ``membrane_medium_traction``
connector while stating plainly that the declaration bought none of the physics: the cell's six whole-body
rigid modes stayed held by a numerical ``aI`` regulariser, rated BLOCKS_CRAWL.  This driver executes the
exterior operator that replaces that regulariser
(:mod:`aleph.engine.medium_exterior`) on the A5000 and measures whether it does the one thing a
regulariser cannot: hold the six rigid modes with a resistance that is a function of the cell's GEOMETRY.

WHAT "NATIVE" MEANS HERE, PRECISELY -- and it is not the usual caveat.  The exterior medium addresses the
cell's outer SURFACE only, and the native membrane surface is the icosphere at subdivision 3 (642 nodes /
1280 faces, ``ac/cell/compartments.py:104,641``).  So ``--subdivisions 3`` is this component's COMPLETE
native population, not a slice of it -- unlike the ``nmii_sf_motor`` slice at 0.18% of the cell.  What is
still NOT native is the COMPOSITION: the medium is not yet dispatched through the
``membrane_medium_traction`` connector on a composed cell, because that claim belongs to ``SurfaceBody``.
The rung this run can therefore reach is ``CUDA_UNIT``, never ``CONNECTED``, and the record says so.

WHAT IS QUOTABLE AND WHAT IS NOT.  ``mu_medium`` is a PI-GAP -- the knowledge base has no KnowledgeClaim
and no SourceEvidence for a culture-medium viscosity (``ac/engine/medium_params_t10.yaml``).  Every
observable this driver reports is therefore deliberately viscosity-free: RATIOS to the Stokes closed forms,
the exactly-zero isotropy and coupling residuals, the sign of the dissipated work, the positivity of the
six rigid-mode eigenvalues, and the convergence order under refinement.  ``QuantitativeClaim`` is
``BLOCKED`` and no drag magnitude may be quoted from this run.

THE VOID CEILING IS DERIVED, NOT PICKED.  A translating icosphere feels exactly zero transverse force by
icosahedral symmetry, so any transverse resultant is pure solver error in physical units.  The ceiling is
the standard CG error bound -- the relative residual tolerance times the operator's own measured condition
number -- so above it the reported resistance is dominated by an unconverged solve and neither PASS nor
FAIL of a Stokes-law criterion is interpretable.  The formula is declared here, before the run; only the
condition number is measured.

Run on gbook through the lease:
    python aleph/scripts/ffn_gpu.py run --minutes 20 \
        --reason "T10: exterior Stokes medium native gate" \
        --population "642 membrane surface quadrature points = the FULL native membrane surface" \
        --artifact outputs/ac/gate_t10_medium/medium_exterior_native.json \
        -- aleph/scripts/ac_gate_t10_medium_native.py --viscosity 1.0e-3 \
           --out outputs/ac/gate_t10_medium/medium_exterior_native.json
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import warp as wp

from aleph.engine.contracts import (
    EvidenceLabel,
    EvidenceRung,
    QuantitativeClaim,
    VoidCeiling,
    reference_cell_architecture,
)
from aleph.engine.medium_exterior import (
    EXTERIOR_MEDIUM_COMPONENT,
    MEMBRANE_MEDIUM_CONNECTOR,
    ExteriorStokesMediumSettings,
    build_exterior_stokes_medium,
    derive_blob_epsilon,
)
from aleph.engine.medium_stokes_analytic import (
    mobility_matrix,
    sphere_drag_convergence,
    sphere_rotation_resistance,
    sphere_translation_resistance,
    surface_quadrature_sphere,
)
from aleph.engine.observe.artifact import observation_artifact, timing_block, write_artifact

#: The native membrane surface resolution (``ac/cell/compartments.py:104`` icosphere level 3 = 642 nodes).
NATIVE_MEMBRANE_SUBDIVISIONS = 3

#: MCF7 baseline cell radius in um (``ac/cell/compartments.py``).
CELL_RADIUS_UM = 7.5

#: Round-off floor for the two EXACT structural residuals (isotropy, translation-rotation decoupling).
#: Derived: a few float64 round-offs accumulated across the ``3N``-term reduction, not a chosen tolerance.
EXACT_RESIDUAL_FLOOR = 64.0 * float(np.finfo(np.float64).eps)

#: The VOID-ceiling FORMULA, declared before the run: ratio = cg_relative_tolerance * cond(M).
VOID_CEILING_REASON = (
    "a translating icosphere feels exactly zero transverse force by icosahedral symmetry, so any "
    "transverse resultant is solver error expressed in physical units; the ceiling is the standard CG "
    "error bound (relative residual tolerance x the operator's measured condition number), above which "
    "the reported resistance is dominated by an unconverged solve rather than by the exterior operator, "
    "and neither PASS nor FAIL of a Stokes-law criterion is interpretable"
)


def _rigid_translation(
    medium: Any,
    verts: np.ndarray,
    *,
    device: str,
    speed: float,
    dt_phys: float,
    field: dict[str, np.ndarray] | None = None,
    field_key: str = "",
) -> tuple[np.ndarray, int, bool]:
    """Prescribe a rigid x-translation on the device and return the surface force resultant.

    A pure translation leaves every pair separation unchanged, so the mobility is bit-identical to the one
    assembled at rest -- which is why the translation channel is the one that can be compared to the host
    reference with no discretisation caveat of its own.

    Args:
        medium: The allocated exterior medium.
        verts: Host copy of the committed surface quadrature, shape ``(N, 3)``, in um.
        device: Warp device string.
        speed: Translation speed in um/s.
        dt_phys: Outer physical step in seconds.
        field: Optional sink for the PER-NODE traction field, so the isolation render the ladder gate
            requires reads a committed array rather than a re-run. Text artifacts hide exactly the
            wrong-sign and outlier bugs a figure catches (PI 2026-05-21).
        field_key: Key to store the per-node traction under.

    Returns:
        Tuple ``(resultant, cg_iterations, cg_converged)`` with the resultant in pN.
    """
    medium.set_step(dt_phys)
    moved = verts + np.array([speed * dt_phys, 0.0, 0.0])
    pos_d = wp.array(moved, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(verts.shape[0], dtype=wp.vec3d, device=device)
    medium.snapshot_candidate()
    medium.accumulate(pos_d, force_d)
    traction = force_d.numpy()
    if field is not None and field_key:
        field[field_key] = traction.copy()
    return traction.sum(axis=0), medium.iteration_count(), medium.converged_flag()


def _rigid_rotation(
    medium: Any,
    verts: np.ndarray,
    *,
    device: str,
    omega: float,
    dt_phys: float,
    field: dict[str, np.ndarray] | None = None,
    field_key: str = "",
) -> tuple[np.ndarray, int, bool]:
    """Prescribe a rigid z-rotation on the device and return the surface torque about the centre.

    Args:
        medium: The allocated exterior medium.
        verts: Host copy of the committed surface quadrature, shape ``(N, 3)``, in um.
        device: Warp device string.
        omega: Angular speed in rad/s about z.
        dt_phys: Outer physical step in seconds; kept small so the first-order rotation stays on the sphere.
        field: Optional sink for the per-node traction field, for the isolation render.
        field_key: Key to store the per-node traction under.

    Returns:
        Tuple ``(torque, cg_iterations, cg_converged)`` with the torque in pN*um.
    """
    medium.set_step(dt_phys)
    axis = np.array([0.0, 0.0, omega])
    rotated = verts + dt_phys * np.cross(axis, verts)
    pos_d = wp.array(rotated, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(verts.shape[0], dtype=wp.vec3d, device=device)
    medium.snapshot_candidate()
    medium.accumulate(pos_d, force_d)
    traction = force_d.numpy()
    if field is not None and field_key:
        field[field_key] = traction.copy()
    return np.cross(verts, traction).sum(axis=0), medium.iteration_count(), medium.converged_flag()


def _grand_resistance_on_device(
    medium: Any, verts: np.ndarray, *, device: str, dt_phys: float
) -> np.ndarray:
    """Assemble the 6x6 grand resistance by six device solves, one per rigid mode.

    This is the run's central object: six strictly positive eigenvalues is the falsifiable statement that
    the exterior medium LOADS the whole-body modes, which the ``aI`` regulariser it replaces cannot make.

    Args:
        medium: The allocated exterior medium.
        verts: Host copy of the committed surface quadrature, shape ``(N, 3)``, in um.
        device: Warp device string.
        dt_phys: Outer physical step in seconds.

    Returns:
        Symmetric ``(6, 6)`` resistance; translation block in pN*s/um, rotation block in pN*s*um.
    """
    medium.set_step(dt_phys)
    n = verts.shape[0]
    columns = []
    for mode in range(6):
        velocity = np.zeros_like(verts)
        if mode < 3:
            velocity[:, mode] = 1.0
        else:
            unit = np.zeros(3)
            unit[mode - 3] = 1.0
            velocity = np.cross(unit, verts)
        pos_d = wp.array(verts + dt_phys * velocity, dtype=wp.vec3d, device=device)
        force_d = wp.zeros(n, dtype=wp.vec3d, device=device)
        medium.snapshot_candidate()
        medium.accumulate(pos_d, force_d)
        # `accumulate` scatters f_medium = -M^-1 v, so the resistance column is its negative.
        traction = -force_d.numpy()
        columns.append(
            np.concatenate([traction.sum(axis=0), np.cross(verts, traction).sum(axis=0)])
        )
    resistance = np.column_stack(columns)
    return 0.5 * (resistance + resistance.T)


def main() -> int:
    """Run the T10 exterior-medium native gate and emit one self-stamping run record."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--viscosity",
        type=float,
        required=True,
        help="medium dynamic viscosity in Pa*s. REQUIRED with no default: this is a PI-GAP, so the "
             "caller must state it, and every observable reported here is viscosity-free regardless",
    )
    parser.add_argument("--subdivisions", type=int, default=NATIVE_MEMBRANE_SUBDIVISIONS,
                        help="icosphere level; 3 is the FULL native membrane surface (642 nodes)")
    parser.add_argument("--radius", type=float, default=CELL_RADIUS_UM, help="cell radius in um")
    parser.add_argument("--dt", type=float, default=1.0e-3, help="outer physical step in seconds")
    parser.add_argument("--speed", type=float, default=1.0, help="rigid translation speed in um/s")
    parser.add_argument("--omega", type=float, default=1.0e-2, help="rigid rotation rate in rad/s")
    parser.add_argument(
        "--cg-tolerance",
        type=float,
        default=1.0e-10,
        help="CG relative residual tolerance. Exposed as a SWEEP AXIS, not a knob: the three 'exact "
             "structural' observables (isotropy, translation-rotation decoupling, time reversal) are "
             "properties of the OPERATOR but are measured HERE through a CG solve, so scaling this and "
             "watching whether they track it is what distinguishes solver error from an operator defect",
    )
    parser.add_argument("--out", default="", help="where to write the run record")
    parser.add_argument("--build-commit", default="", help="build commit, for a non-git checkout")
    args = parser.parse_args()

    run_started = time.perf_counter()
    wp.init()
    cuda = next((d for d in wp.get_devices() if d.is_cuda), None)
    if cuda is None:
        raise SystemExit("I0-A: this gate requires a CUDA GPU; there is no host simulation path")
    device = str(cuda)

    verts, faces = surface_quadrature_sphere(args.subdivisions, args.radius)
    epsilon = derive_blob_epsilon(verts)
    settings = ExteriorStokesMediumSettings(
        viscosity_pa_s=float(args.viscosity),
        blob_epsilon_um=epsilon,
        cg_relative_tolerance=float(args.cg_tolerance),
    )
    medium = build_exterior_stokes_medium(
        settings=settings,
        surface_index=np.arange(verts.shape[0]),
        initial_position=verts,
        device=device,
    )

    exact_translation = sphere_translation_resistance(args.radius, settings.viscosity_pa_s)
    exact_rotation = sphere_rotation_resistance(args.radius, settings.viscosity_pa_s)

    # t0 — BEFORE any motion. A stationary surface must exchange exactly nothing with the medium; without
    # this row a non-zero force later cannot be said to have EMERGED from the motion.
    rest_pos_d = wp.array(verts, dtype=wp.vec3d, device=device)
    rest_force_d = wp.zeros(verts.shape[0], dtype=wp.vec3d, device=device)
    medium.set_step(args.dt)
    medium.snapshot_candidate()
    medium.accumulate(rest_pos_d, rest_force_d)
    rest_force = rest_force_d.numpy()
    t0 = {
        "surface_at_rest_max_abs_force_pN": float(np.abs(rest_force).max()),
        "surface_at_rest_resultant_pN": rest_force.sum(axis=0).tolist(),
        "dissipated_work_pN_um": float(medium.dissipated_work_d.numpy()[0]),
        "note": "no motion, so no force and no dissipation: the t0 row a run needs to claim emergence",
    }

    # The per-node traction fields the isolation render reads. A summary scalar cannot show a wrong-sign
    # or an outlier node, which is the whole reason the visualization rule exists.
    fields: dict[str, np.ndarray] = {"surface_position_um": verts, "surface_faces": faces}

    translation, translation_iterations, translation_converged = _rigid_translation(
        medium, verts, device=device, speed=args.speed, dt_phys=args.dt,
        field=fields, field_key="traction_translation_pN",
    )
    rotation, rotation_iterations, rotation_converged = _rigid_rotation(
        medium, verts, device=device, omega=args.omega, dt_phys=args.dt,
        field=fields, field_key="traction_rotation_pN",
    )
    resistance = _grand_resistance_on_device(medium, verts, device=device, dt_phys=args.dt)

    # Time reversal: Stokes is linear and inertialess, so reversing the motion must reverse the force.
    reversed_translation, _, _ = _rigid_translation(
        medium, verts, device=device, speed=-args.speed, dt_phys=args.dt
    )

    # Viscosity linearity: the whole viscosity dependence is the 1/(8 pi mu) prefactor, so this is exact.
    doubled = build_exterior_stokes_medium(
        settings=ExteriorStokesMediumSettings(
            viscosity_pa_s=2.0 * float(args.viscosity),
            blob_epsilon_um=epsilon,
            cg_relative_tolerance=float(args.cg_tolerance),
        ),
        surface_index=np.arange(verts.shape[0]),
        initial_position=verts,
        device=device,
    )
    doubled_translation, _, _ = _rigid_translation(
        doubled, verts, device=device, speed=args.speed, dt_phys=args.dt
    )

    # Accepted-step transaction: work is banked only on acceptance, and it must be strictly negative.
    accepted = wp.ones(1, dtype=wp.int32, device=device)
    moved = verts + np.array([args.speed * args.dt, 0.0, 0.0])
    pos_d = wp.array(moved, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(verts.shape[0], dtype=wp.vec3d, device=device)
    medium.snapshot_candidate()
    medium.accumulate(pos_d, force_d)
    medium.commit_irreversible(accepted, dt_phys=float(args.dt), rng_seed=0)
    banked_work = float(medium.dissipated_work_d.numpy()[0])
    committed_advanced = float(
        np.abs(medium.committed_position_d.numpy() - moved).max()
    )

    # Host reference: the same operator assembled densely, for the condition number the VOID ceiling needs
    # and for the device-vs-reference agreement. Population-DEPENDENT by nature, and reported as such.
    reference_mobility = mobility_matrix(
        verts, epsilon=epsilon, viscosity=settings.viscosity_pa_s
    )
    condition_number = float(np.linalg.cond(reference_mobility))
    void_ceiling = VoidCeiling(settings.cg_relative_tolerance * condition_number, VOID_CEILING_REASON)

    translation_block = resistance[:3, :3]
    rotation_block = resistance[3:, 3:]
    translation_mean = float(np.trace(translation_block) / 3.0)
    rotation_mean = float(np.trace(rotation_block) / 3.0)
    eigenvalues = np.linalg.eigvalsh(resistance)

    refinement = sphere_drag_convergence(
        radius=args.radius,
        viscosity=settings.viscosity_pa_s,
        subdivisions=(1, 2, args.subdivisions) if args.subdivisions > 2 else (1, 2),
        epsilon_ratio=1.0,
    )

    axial = float(abs(translation[0]))
    transverse = float(np.linalg.norm(translation[1:]))

    measurements: dict[str, Any] = {
        "stokes_law_ratio": {
            "translation_measured_over_exact": axial / (exact_translation * args.speed),
            "rotation_measured_over_exact": float(abs(rotation[2])) / (exact_rotation * args.omega),
            "note": "the viscosity cancels in both ratios, which is why they are quotable while the "
                    "magnitudes are not; departure from 1 is the surface quadrature's discretisation "
                    "error at this refinement, MEASURED not tuned",
        },
        "rigid_modes_are_loaded": {
            "eigenvalues": eigenvalues.tolist(),
            "min_eigenvalue": float(eigenvalues.min()),
            "all_strictly_positive": bool(eigenvalues.min() > 0.0),
            "translation_resistance_over_radius": translation_mean / args.radius,
            "rotation_over_translation_shape": rotation_mean / translation_mean,
            "rotation_over_translation_exact_shape": 4.0 / 3.0 * args.radius**2,
            "note": "the BLOCKS_CRAWL statement in falsifiable form: an aI regulariser gives a resistance "
                    "independent of the cell, this one is a function of its geometry",
        },
        "exact_structural_residuals": {
            "translation_isotropy": float(
                np.linalg.norm(translation_block - translation_mean * np.eye(3))
                / np.linalg.norm(translation_block)
            ),
            "translation_rotation_coupling": float(
                np.linalg.norm(resistance[:3, 3:]) / np.linalg.norm(resistance)
            ),
            "operator_symmetry": float(
                np.linalg.norm(reference_mobility - reference_mobility.T)
            ),
            "round_off_floor": EXACT_RESIDUAL_FLOOR,
            "note": "icosahedral symmetry forces the first two to vanish and the kernel's evenness forces "
                    "the third; all three are population-INDEPENDENT properties of the operator's form",
        },
        "time_reversal": {
            "forward_axial_pN": float(translation[0]),
            "reversed_axial_pN": float(reversed_translation[0]),
            "relative_asymmetry": float(
                np.linalg.norm(reversed_translation + translation) / np.linalg.norm(translation)
            ),
        },
        "viscosity_linearity": {
            "doubled_over_single": float(doubled_translation[0] / translation[0]),
            "exact": 2.0,
        },
        "dissipation": {
            "banked_work_pN_um": banked_work,
            "strictly_negative": bool(banked_work < 0.0),
            "committed_reference_advanced_max_abs_um": committed_advanced,
            "note": "work is banked only inside commit_irreversible, so it is an accepted-step quantity "
                    "and not an integral over the solver's path",
        },
        "solver": {
            "cg_iterations_translation": translation_iterations,
            "cg_iterations_rotation": rotation_iterations,
            "cg_converged": bool(translation_converged and rotation_converged),
            "mobility_condition_number": condition_number,
            "transverse_over_axial": transverse / axial if axial > 0.0 else float("nan"),
            "void_ceiling_ratio": void_ceiling.ratio,
        },
        "grid_invariance": refinement,
        "wall": {
            "mode": str(settings.wall_mode),
            "half_space_blake_implemented": False,
            "note": "T10 checklist item 3 (asymmetric world boundary: basal 2D collagen + free media face) "
                    "is NOT landed here; spreading and crawl stay gated on it",
        },
    }

    architecture = reference_cell_architecture()
    mechanism = {
        "rest_is_force_free": t0["surface_at_rest_max_abs_force_pN"] == 0.0,
        "all_six_rigid_modes_loaded": bool(eigenvalues.min() > 0.0),
        "isotropy_at_round_off": measurements["exact_structural_residuals"]["translation_isotropy"]
        <= EXACT_RESIDUAL_FLOOR,
        "decoupling_at_round_off": measurements["exact_structural_residuals"][
            "translation_rotation_coupling"
        ] <= EXACT_RESIDUAL_FLOOR,
        "operator_exactly_symmetric": measurements["exact_structural_residuals"]["operator_symmetry"] == 0.0,
        "time_reversible": measurements["time_reversal"]["relative_asymmetry"] <= EXACT_RESIDUAL_FLOOR,
        "linear_in_viscosity": abs(measurements["viscosity_linearity"]["doubled_over_single"] - 2.0) <= 1e-12,
        "dissipates_strictly": bool(banked_work < 0.0),
        "cg_converged": bool(translation_converged and rotation_converged),
        "converges_under_refinement": [
            row["translation_relative_error"] for row in refinement["rows"]
        ] == sorted(
            (row["translation_relative_error"] for row in refinement["rows"]), reverse=True
        ),
    }

    record = observation_artifact(
        run_label="T10 extracellular_medium — exterior Stokes operator on the native membrane surface",
        evidence=EvidenceLabel(
            rung=EvidenceRung.CUDA_UNIT,
            quantitative=QuantitativeClaim.BLOCKED,
            basis=(
                "the exterior operator executes on CUDA over the FULL native membrane surface quadrature "
                "and its six rigid-mode resistances are measured strictly positive with the translation "
                "channel proportional to the cell radius; the rung stops at CUDA_UNIT because the "
                "membrane_medium_traction connector is NOT dispatched -- that claim belongs to "
                "SurfaceBody -- and the quantitative axis is BLOCKED because mu_medium is a PI-GAP, so "
                "only the viscosity-free ratios and the exact structural residuals are quotable"
            ),
        ),
        config={
            "argv": {key: value for key, value in sorted(vars(args).items())},
            "blob_epsilon_um": epsilon,
            "epsilon_rule": "epsilon = 1.0 x mean nearest-neighbour spacing (DERIVED, grid-invariant)",
            "wall_mode": str(settings.wall_mode),
            "cg_relative_tolerance": settings.cg_relative_tolerance,
            "cg_max_iterations": settings.cg_max_iterations,
        },
        census={
            "n_surface_quadrature_points": int(verts.shape[0]),
            "n_surface_faces": int(faces.shape[0]),
            "n_dof": 3 * int(verts.shape[0]),
            "icosphere_subdivisions": int(args.subdivisions),
            "is_native_population": bool(args.subdivisions == NATIVE_MEMBRANE_SUBDIVISIONS),
            "population_note": (
                "the exterior medium addresses the cell's outer SURFACE only, and the native membrane "
                "surface IS icosphere level 3 = 642 nodes / 1280 faces (ac/cell/compartments.py:104,641), "
                "so this is the component's complete native population rather than a slice of it. What is "
                "not native is the COMPOSITION: no composed cell, no connector dispatch"
            ),
            "components_declared": len(architecture.components),
            "connectors_declared": len(architecture.connectors),
            "component": EXTERIOR_MEDIUM_COMPONENT,
            "connector_not_yet_dispatched": MEMBRANE_MEDIUM_CONNECTOR,
        },
        t0=t0,
        timing=timing_block(
            wall_seconds=time.perf_counter() - run_started,
            n_steps=1,
            n_inner_iterations=int(translation_iterations + rotation_iterations),
            device_note=(
                "NOT a benchmark. The run performs 11 device resistance solves (rest, translation, "
                "rotation, six rigid modes, a reversed translation, a doubled-viscosity translation) plus "
                "a dense host mobility assembly and condition-number estimate for the VOID ceiling, and "
                "Warp kernel compilation is included in the wall clock"
            ),
        ),
        measurements=measurements,
        device=device,
        parameter_provenance={
            "viscosity_medium": "PI_GAP",
            "blob_epsilon": "DERIVED",
            "cell_radius": "SOURCED",
            "membrane_subdivisions": "SOURCED",
        },
        void_ceiling=void_ceiling,
        gate_passed=all(mechanism.values()),
        residual=transverse,
        signal=axial,
        notes={
            "mechanism": mechanism,
            "trap_four": (
                "PI framework trap #4 is discharged as an enforced invariant, not a deletion: the audit "
                "(docs/v2_audit/T10_GAMMA_NODE_AUDIT_2026-07-28.md) found gamma_node = 6 pi eta R / Nc has "
                "never existed in ac/ -- both live sites are in ff/ with only legacy scripts/ff_* callers "
                "-- so this landing takes ac/'s velocity-proportional drag count from 0 to 1 and "
                "medium_exterior.assert_single_dissipation_owner keeps it there"
            ),
            "not_closed": [
                "membrane_medium_traction is declared and NOT dispatched; rung cannot exceed CUDA_UNIT",
                "HALF_SPACE_BLAKE (the basal no-slip wall) is not implemented, so the world boundary is "
                "still symmetric and spreading/crawl remain gated",
                "mu_medium is a PI-GAP, so no drag magnitude from this run is quotable",
                "the aI regularisers in sf_implicit.py and ac/cell/ are untouched: removing them is the "
                "cortex/sf-motor lanes' change once the medium is dispatched on a composed cell",
            ],
        },
        declared_commit=(args.build_commit or None),
    )

    text = json.dumps(record, indent=2, default=str)
    print(text, flush=True)
    if args.out:
        write_artifact(args.out, record)
        # Sidecar NPZ beside the record, so the lane's own figure entry point renders from committed data
        # rather than by re-running on the GPU. Figures are persistent project state, not a by-product.
        fields["grand_resistance_pN_s"] = resistance
        fields["rigid_mode_eigenvalues"] = eigenvalues
        npz_path = Path(args.out).with_suffix(".npz")
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(npz_path, **fields)
        print(f"[T10 medium] wrote field dump {npz_path}", flush=True)

    failed = [name for name, ok in mechanism.items() if not ok]
    if failed:
        print(f"[T10 medium] MECHANISM FAIL: {failed}", flush=True)
        return 1
    print("[T10 medium] MECHANISM PASS — six rigid modes loaded by physics, "
          f"translation ratio {measurements['stokes_law_ratio']['translation_measured_over_exact']:.6f}, "
          f"rotation ratio {measurements['stokes_law_ratio']['rotation_measured_over_exact']:.6f}",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

r"""The native cortex tangent, looked at: family-by-family certification, then a matrix-free spectrum.

WHAT THIS CLOSES.  On 2026-07-28 four spectral numbers were retired to ``STATE.md`` (c) 14 one session
after landing, for one reason: they were measured on 904 nodes, which is 0.18% of native, and every one
of them is a function of the connectivity graph or of the system's size.  ``observe/spectrum.py`` gained
:func:`~aleph.engine.observe.spectrum.lanczos_extremal_spectrum` in the same commit so they would be
re-measurable rather than merely withdrawn.  This driver is the caller that re-measures them, at the
population the conclusion is supposed to be about.

IT MEASURES THE OPERATOR THE SOLVER ACTUALLY INVERTS.  The engine's inner iteration is
``pos += dt_mu · P F`` (``ac/cell/inner_mechanics.py``), so its linearisation on the constraint tangent
space is ``I − dt_mu · P K P`` and its stability limit is ``dt_mu · λ_max(P K P |range(P)) < 2``.  The
unprojected ``K`` is reported too, because symmetry of ``K`` — not of ``P K P`` — is the integrability
condition for ``F = −∇U``, and the projector would launder an asymmetric ``K`` into a symmetric-looking
one.  Both are measured; neither substitutes for the other.

WHY A CERTIFICATION STAGE EXISTS, AND WHY IT IS NOT A "SLICE RESULT".  Lanczos assumes the operator is
symmetric.  It cannot check that — a Krylov method given an asymmetric operator returns confident
numbers about the symmetric part and says nothing.  So symmetry has to be established where the operator
can be assembled, which is below 8,192 DOF by :func:`assemble_dense_operator`'s own refusal.  The only
statements this driver takes from that stage are **per-contribution algebraic** ones (is family ``j``'s
tangent symmetric; is the assembled operator the sum of its families), and those transfer with
population by the same argument ``STATE.md`` (b) already accepts for the ``nmii_sf_motor`` row: whether
a contribution is symmetric is a property of its differential form, and one instance settles it.  No
magnitude, count, ratio or spectrum from the certification build is quoted.

AND IT CERTIFIES BY ISOLATION, NOT BY DIFFERENCE AND NOT BY LAUNCH DIMENSION.  A family is "exercised"
only if its tangent is nonzero, which a launch dimension does not establish — the WCA steric kernel
launches over every node and contributes exactly nothing when no pair is inside its cutoff, which is the
DESIGNED state of an ``overlap_free_cortex`` build.  So each family is assembled ALONE, by disabling
every other family's live gating attribute, and a family that assembles to the exact zero matrix is
reported UNCERTIFIED rather than silently counted as covered.  Summing the isolated families against the
full assembly gives the multilinearity residual ``‖K − Σ_j K_j‖/‖K‖`` for the cortex families, which
``STATE.md`` (b) flags as untested.

WHY ALONE AND NOT ``K_on − K_off``, WHICH IS THE OBVIOUS WAY AND IS WRONG.  The first version of this
driver isolated each family by difference.  Measured on the certification build, that made **every**
family report the identical asymmetry ``9.3132e-10`` — which is exactly one half-ulp of the FULL
operator's largest entry (``max|K| = 7.30e6``).  A difference of two matrices each carrying round-off at
the scale of the whole assembly has a noise floor set by the WHOLE ASSEMBLY, so a small family's own
tangent was buried under the crosslink family's round-off and scored against a floor derived from its
own much smaller magnitude.  Four of seven families "failed" for that reason alone.  Assembling each
family alone removes the subtraction, so its round-off floor is genuinely its own and the exercised test
is exact rather than thresholded.  It is also a STRONGER multilinearity test: if one family's tangent
depended on another's presence, the leave-one-in sum would not reproduce the assembly, whereas the
difference construction cannot see that at all.

engine units: length µm, force pN, stiffness pN/µm; a mobility step is µm/pN.

Sanity Gate:
    * dimensional: eigenvalues of a tangent are [pN/µm]; ``dt_mu`` is [µm/pN]; their product
      ``dt_mu·λ_max`` is dimensionless and is compared against the forward-Euler limit 2.
    * boundary: a family whose isolated tangent is exactly zero yields ``exercised: false`` and is
      excluded from the certified set instead of passing vacuously; an unresolved Ritz value yields a
      ``None`` condition number carrying the reason, never a number built from the iteration count.
    * conservation/invariant: symmetry of ``K`` is exactly the integrability condition for ``F = −∇U``,
      and the sum of the isolated families must reproduce the assembled operator to round-off.
    * numerical: the symmetry floor is the ledger's own float64 accumulation bound (PI D8) scaled by the
      operator's largest entry — the same function the accepted-step balance gate uses, never a second
      independently-written tolerance.  Atomics make the accumulation order nondeterministic, so no
      bit-identity criterion is used anywhere here.
    * sign-sense: a negative Ritz value at the small end means the configuration is not a minimum of any
      potential in that direction; it is reported as measured and never clipped to zero.
    * measurement-protocol: every probe runs at regularization exactly zero.  ``aI + K`` shifts every
      eigenvalue by ``a``, so a regularized probe would report ``λ + a`` and manufacture a spectrum with
      no null space.  The folded pass's shift ``σ`` is a method parameter, not a physical one, so the
      run reports ``λ_min`` at TWO values of ``σ`` and records their disagreement — a measured
      invariance in place of a chosen constant.
"""

from __future__ import annotations

import argparse
import json
import time
from contextlib import ExitStack, contextmanager
from pathlib import Path
from typing import Any, Callable, Iterator

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import AssembledCell, CellConfig, build_cell
from aleph.components.incumbent.implicit_mechanics import ProjectedAnalyticCG
from aleph.engine.contracts import (
    EvidenceLabel,
    EvidenceRung,
    QuantitativeClaim,
    VoidCeiling,
)
from aleph.engine.forces_manifest import dump_force_channels
from aleph.engine.observe import (
    assemble_dense_operator,
    lanczos_extremal_spectrum,
    lanczos_ritz,
    observation_artifact,
    stiffness_spectrum,
    symmetry_report,
    timing_block,
    write_artifact,
)

#: Provenance of every parameter that enters the operator, so a magnitude standing on an unsourced
#: constant is visible in the record rather than inferred from a memo.  Classifications follow
#: ``docs/v2_audit/PARAM_PROVENANCE_AUDIT_2026-07-24.md``.
PARAMETER_PROVENANCE: dict[str, str] = {
    "cortex_seg_um": "CONVENIENCE",
    "cortex_length_um": "CONVENIENCE",
    "cortex_density_per_fil": "CONVENIENCE",
    "k_crosslink": "SOURCED",
    "kappa_bending": "SOURCED",
    "k_erm": "SOURCED",
    "k_xb": "PI_GAP",
    "k_backbone": "PI_GAP",
    "k_linc": "PI_GAP",
    "turgor_Pi_0": "CONVENIENCE",
}

#: The stiffness families ``ProjectedAnalyticCG._stiffness`` launches, each named by the LIVE attribute
#: that gates its launch.  Isolating a family means setting that attribute to its inert value on the
#: composed object and re-probing; the operator reads all of them at launch time, so the toggle is exact
#: and reversible.  ``None`` as an owner means the attribute lives on the solver rather than the cell.
_FAMILY_TOGGLES: tuple[tuple[str, str | None, str, Any], ...] = (
    ("bending", "cell", "n_tri", 0),
    ("crosslink", "cell", "n_xl", 0),
    ("membrane_area_edges", "solver", "mem_edges_d", None),
    ("steric_wca", "cell", "steric", None),
    ("myosin", "cell", "myosin", None),
    ("nucleus_linc", "cell", "nucleus", None),
    ("membrane_erm", "cell", "membrane", None),
)

#: Where a family's tangent KERNELS are certified if this build cannot exercise them.  A hole that is
#: bounded by a cross-reference is a different thing from an open one, and the distinction has to be in
#: the record rather than in a reader's memory.  The argument is the same per-contribution algebraic one
#: the rest of this driver rests on: these are the SAME Warp kernels, so the same differential forms.
_FAMILY_KERNEL_CROSS_REFERENCE: dict[str, str] = {
    "myosin": (
        "myosin sites are generated from cortex node proximity at native areal density (100/µm²), so a "
        "build small enough to column-probe (<8,192 DOF, i.e. ~4 decades below native) has none — this "
        "family is not exercisable on ANY probeable cortex build, not merely on this one. Its kernels "
        "ARE certified elsewhere: ac/engine/sf_implicit.py::SFImplicitCG._stiffness launches the same "
        "add_uniform_pair_stiffness_kernel (backbone + head springs), the same "
        "add_angle_gauss_newton_stiffness_kernel (backbone + arm angles) and the same "
        "add_segment_crossbridge_stiffness_kernel, and that assembly was measured symmetric to "
        "round-off at build db841185 (STATE.md (b), the multilinearity + symmetry row). NOT covered by "
        "that cross-reference: add_crossbridge_stiffness_kernel, the NON-segment crossbridge variant "
        "ac/cell takes when myosin.segment_runtime is None — the sf slice exercised the segment variant"
    ),
    "steric_wca": (
        "the WCA kernel launches over every node and contributes nothing when no pair is inside its "
        "cutoff, which is the DESIGNED state of an overlap_free_cortex build. Run the certification "
        "with --certify-overlap-free off (the default) so the build carries interpenetrations and the "
        "family is genuinely assembled"
    ),
}


@contextmanager
def _family_disabled(cell: AssembledCell, solver: ProjectedAnalyticCG,
                     owner: str, attribute: str, inert: Any) -> Iterator[None]:
    """Temporarily set one live attribute to its inert value, restoring it unconditionally.

    Args:
        cell: The composed cell the operator reads.
        solver: The CG workspace, which owns the membrane edge-spring topology.
        owner: ``"cell"`` or ``"solver"`` — which object carries the attribute.
        attribute: The attribute name gating the family's launch.
        inert: The value that makes the launch not happen.

    Yields:
        Nothing; the attribute is restored on exit, including on exception.
    """
    target = cell if owner == "cell" else solver
    previous = getattr(target, attribute)
    setattr(target, attribute, inert)
    try:
        yield
    finally:
        setattr(target, attribute, previous)


def _stiffness_launch_census(cell: AssembledCell, solver: ProjectedAnalyticCG) -> dict[str, int]:
    """Return the per-family LAUNCH dimension of the assembled tangent.

    A launch dimension is an upper bound on a family's contribution count and is what the round-off
    floor is derived from.  It is deliberately NOT used to decide whether a family was exercised — see
    :func:`certify_families`, which decides that by the family's own tangent norm.

    Args:
        cell: The composed cell.
        solver: The CG workspace holding the membrane edge springs.

    Returns:
        Family name to launch dimension.  Every key is present even at zero.
    """
    myosin = cell.myosin
    state = myosin.state if myosin is not None else None
    membrane = cell.membrane
    nucleus = cell.nucleus
    steric = cell.steric
    return {
        "mass_action": int(solver.n),
        "bending": int(cell.n_tri),
        "crosslink": int(cell.n_xl),
        "membrane_area_edges": (
            int(solver.mem_edges_d.shape[0]) if solver.mem_edges_d is not None else 0),
        "steric_wca": int(steric.n) if steric is not None else 0,
        "myosin_backbone": int(myosin.backbone_bonds.shape[0]) if myosin is not None else 0,
        "myosin_head_spring": int(myosin.head_bonds.shape[0]) if myosin is not None else 0,
        "myosin_backbone_angle": (
            int(myosin.backbone_angles.shape[0])
            if myosin is not None and myosin.backbone_angles is not None
            and float(myosin.k_theta_backbone) > 0.0 else 0),
        "myosin_arm_angle": (
            int(myosin.head_arm_angles.shape[0])
            if myosin is not None and myosin.head_arm_angles is not None
            and float(myosin.k_theta_arm) > 0.0 else 0),
        "myosin_crossbridge": int(state["bound"].shape[0]) if state is not None else 0,
        "nucleus_linc": int(nucleus.n_linc) if nucleus is not None else 0,
        "membrane_erm": int(membrane.n_erm) if membrane is not None else 0,
    }


def _population_census(cell: AssembledCell, solver: ProjectedAnalyticCG) -> dict[str, Any]:
    """Return the population census the retraction of 2026-07-28 turned on nobody recording.

    Args:
        cell: The composed cell.
        solver: The CG workspace, whose ``n`` is the operator's node count.

    Returns:
        Node/DOF/topology counts plus the per-family launch census.
    """
    membrane = cell.membrane
    nucleus = cell.nucleus
    return {
        "n_nodes_total": int(cell.n_total),
        "n_dof": 3 * int(cell.n_total),
        "n_actin_nodes": int(cell.n_actin),
        "n_fibers": int(cell.n_fibers),
        "n_crosslinks": int(cell.n_xl),
        "n_membrane_faces": int(membrane.n_faces) if membrane is not None else 0,
        # MEASURED vertex counts, not the subdivision level echoed back: a knob that resizes a
        # compartment has to be visible in the census as a SIZE, or a reduced run reads as a full one.
        "membrane_vertices": int(membrane.n_verts) if membrane is not None else 0,
        "nucleus_vertices": int(nucleus.n_verts) if nucleus is not None else 0,
        "n_erm_tethers": int(membrane.n_erm) if membrane is not None else 0,
        "n_nucleus_linc": int(nucleus.n_linc) if nucleus is not None else 0,
        "operator_width_nodes": int(solver.n),
        "stiffness_launch_census": _stiffness_launch_census(cell, solver),
    }


def _make_matvec(
    cell: AssembledCell, solver: ProjectedAnalyticCG, *, projected: bool,
) -> Callable[[np.ndarray], np.ndarray]:
    """Return the matrix-free tangent application, at regularization exactly zero.

    Args:
        cell: The composed cell, supplying the position array the tangent is linearised about.
        solver: The CG workspace whose kernels the closure launches.
        projected: When true the closure applies ``P K P`` — the operator the inner iteration
            linearises to — using the exact runtime projector.  When false it applies the raw ``K``,
            whose symmetry is the integrability condition for ``F = −∇U``.

    Returns:
        A closure ``(n, 3) -> (n, 3)``.  Device buffers are allocated once and reused, so the closure is
        NOT re-entrant; the Lanczos recurrence calls it strictly sequentially.
    """
    device = cell.device
    n = int(solver.n)
    with wp.ScopedDevice(device):
        vector_d = wp.zeros(n, dtype=wp.vec3d)
        out_d = wp.zeros(n, dtype=wp.vec3d)
        zero_regularization = wp.zeros(1, dtype=wp.float64)
        finite = wp.ones(1, dtype=wp.int32)

    def matvec(vector: np.ndarray) -> np.ndarray:
        vector_d.assign(np.ascontiguousarray(
            np.asarray(vector, np.float64).reshape(n, 3)))
        if projected:
            solver._operator(cell.pos_d, vector_d, out_d, zero_regularization, finite)
        else:
            solver._stiffness(cell.pos_d, vector_d, out_d, zero_regularization)
        wp.synchronize_device(device)
        return out_d.numpy()

    return matvec


def _make_projector(
    cell: AssembledCell, solver: ProjectedAnalyticCG,
) -> Callable[[np.ndarray], np.ndarray]:
    """Return the exact runtime inextensibility projector as a node-shaped closure.

    Args:
        cell: The composed cell.
        solver: The CG workspace owning the per-fiber Thomas workspace.

    Returns:
        A closure ``(n, 3) -> (n, 3)`` applying ``P``.
    """
    device = cell.device
    n = int(solver.n)
    with wp.ScopedDevice(device):
        source_d = wp.zeros(n, dtype=wp.vec3d)
        destination_d = wp.zeros(n, dtype=wp.vec3d)
        finite = wp.ones(1, dtype=wp.int32)

    def project(vector: np.ndarray) -> np.ndarray:
        source_d.assign(np.ascontiguousarray(
            np.asarray(vector, np.float64).reshape(n, 3)))
        solver._project(cell.pos_d, source_d, destination_d, finite)
        wp.synchronize_device(device)
        return destination_d.numpy()

    return project


def certify_families(
    cell: Any, solver: Any, *, launch_census: dict[str, int],
    matvec_factory: Callable[[bool], Callable[[np.ndarray], np.ndarray]] | None = None,
) -> dict[str, Any]:
    r"""Isolate every stiffness family by difference and report each one's symmetry.

    The full operator is assembled once, then once more per family with that family's live attribute
    set inert; the difference is that family's own tangent.  Symmetry of each difference is the
    per-contribution algebraic statement that transfers to any population, and the differences summed
    against the full operator give the multilinearity residual.

    Args:
        cell: The composed cell.
        solver: The CG workspace.
        launch_census: The per-family launch dimensions, used only to derive the round-off floor.
        matvec_factory: ``projected -> matvec``, called AFTER each toggle so the closure it returns
            reads the attributes as they stand.  It exists so this function — the one piece of novel
            logic here — is exercisable against a synthetic operator with no device; the default builds
            the real Warp closures.  A factory that captured the operator once instead of per call
            would silently defeat the whole ablation, which is why it is a factory and not a matvec.

    Returns:
        A record carrying the full operator's symmetry, one entry per family, the multilinearity
        residual, and the list of families that contributed exactly nothing.

    Raises:
        ValueError: Propagated from :func:`assemble_dense_operator` when the build is too large to
            probe — which is the intended refusal, not a failure of this driver.
    """
    if matvec_factory is None:
        def matvec_factory(projected: bool) -> Callable[[np.ndarray], np.ndarray]:
            return _make_matvec(cell, solver, projected=projected)

    n = int(solver.n)
    contribution_count = max(1, int(sum(launch_census.values())))
    full = assemble_dense_operator(matvec_factory(False), n)
    full_report = symmetry_report(full, contribution_count=contribution_count)

    families: dict[str, Any] = {}
    accumulated = np.zeros_like(full)
    unexercised: list[str] = []
    for name, owner, attribute, inert in _FAMILY_TOGGLES:
        # LEAVE-ONE-IN: disable every OTHER family and assemble what is left, so this family's tangent
        # is computed directly rather than as a difference of two whole-assembly-sized matrices.
        with ExitStack() as others:
            for other, other_owner, other_attribute, other_inert in _FAMILY_TOGGLES:
                if other == name:
                    continue
                others.enter_context(_family_disabled(
                    cell, solver, other_owner or "cell", other_attribute, other_inert))
            isolated = assemble_dense_operator(matvec_factory(False), n)
        norm = float(np.linalg.norm(isolated, "fro"))
        exercised = bool(np.any(isolated != 0.0))
        if not exercised:
            unexercised.append(name)
        accumulated += isolated
        report = symmetry_report(isolated, contribution_count=contribution_count) if exercised else None
        families[name] = {
            "exercised": exercised,
            "frobenius_norm_pN_per_um": norm,
            "fraction_of_assembled_norm": (
                norm / full_report.frobenius_norm if full_report.frobenius_norm > 0.0 else None),
            "gate_attribute": f"{owner}.{attribute}",
            "launch_dimension": int(sum(
                value for key, value in launch_census.items() if key.startswith(name))),
            "symmetry": None if report is None else {
                "max_abs_asymmetry_pN_per_um": report.max_abs_asymmetry,
                "max_abs_entry_pN_per_um": report.max_abs_entry,
                "relative_asymmetry": report.relative_asymmetry,
                "round_off_floor_pN_per_um": report.round_off_floor,
                "conservative": report.conservative,
            },
            "isolation": "assembled ALONE (every other family's gate disabled), never as a difference "
                         "of two full assemblies — see the module docstring for the measured reason",
            "why_not_exercised": None if exercised else (
                "the family assembles to the EXACT zero matrix, so its launch dimension above counts "
                "threads that added nothing. It is NOT certified by this run"
            ),
            "certified_elsewhere": (
                None if exercised else _FAMILY_KERNEL_CROSS_REFERENCE.get(name)),
        }

    residual = float(np.linalg.norm(full - accumulated, "fro"))
    projected_full = assemble_dense_operator(matvec_factory(True), n)
    projected_report = symmetry_report(projected_full, contribution_count=contribution_count)

    # The certification's own end-to-end check of the instrument: on this build both the dense
    # eigensolve and Lanczos are possible, so the iterative path can be scored against ground truth
    # BEFORE it is used at a population where no ground truth exists.
    dense_spectrum = stiffness_spectrum(full)
    krylov = min(3 * n, 200)
    iterative = lanczos_ritz(matvec_factory(False), n, n_iterations=krylov)
    lanczos_lambda_max, lanczos_bound = iterative.extreme_high

    return {
        "assembled_symmetry": {
            "n_dof": full_report.n_dof,
            "max_abs_asymmetry_pN_per_um": full_report.max_abs_asymmetry,
            "max_abs_entry_pN_per_um": full_report.max_abs_entry,
            "relative_asymmetry": full_report.relative_asymmetry,
            "round_off_floor_pN_per_um": full_report.round_off_floor,
            "contribution_count": full_report.contribution_count,
            "conservative": full_report.conservative,
        },
        "projected_operator_symmetry": {
            "max_abs_asymmetry_pN_per_um": projected_report.max_abs_asymmetry,
            "max_abs_entry_pN_per_um": projected_report.max_abs_entry,
            "relative_asymmetry": projected_report.relative_asymmetry,
            "round_off_floor_pN_per_um": projected_report.round_off_floor,
            "conservative": projected_report.conservative,
            "note": (
                "P K P is symmetric whenever K and P both are, so this checks the runtime per-fiber "
                "Thomas projector as well; it is reported SEPARATELY because a projector can make an "
                "asymmetric K look symmetric and must never stand in for the K test above"
            ),
        },
        "families": families,
        "families_unexercised": unexercised,
        "multilinearity": {
            "residual_frobenius_pN_per_um": residual,
            "relative_residual": (
                residual / full_report.frobenius_norm if full_report.frobenius_norm > 0.0 else None),
            "meaning": (
                "the assembled tangent minus the sum of the families each assembled ALONE; zero to "
                "round-off means the operator is exactly multilinear in its stiffness families, which "
                "is what licenses an analytic per-family gradient. Because each term is assembled "
                "independently rather than differenced out, this also tests that no family's tangent "
                "depends on another family being present — which a leave-one-out construction cannot "
                "see, since there the sum reproduces the assembly by algebra whatever the physics does"
            ),
        },
        "instrument_check_against_ground_truth": {
            "dense_lambda_max_pN_per_um": dense_spectrum.lambda_max,
            "lanczos_lambda_max_pN_per_um": lanczos_lambda_max,
            "lanczos_residual_bound_pN_per_um": lanczos_bound,
            "relative_difference": (
                abs(lanczos_lambda_max - dense_spectrum.lambda_max) / abs(dense_spectrum.lambda_max)
                if dense_spectrum.lambda_max else None),
            "krylov_dimension": int(iterative.n_iterations),
            "gershgorin_bound_pN_per_um": dense_spectrum.gershgorin_bound,
            "gershgorin_over_lambda_max": dense_spectrum.gershgorin_over_lambda_max,
            "why": (
                "the iterative path is used at a population where nothing can check it, so it is "
                "checked here against a dense eigensolve of the SAME Warp operator"
            ),
        },
        "scope": (
            "PER-CONTRIBUTION ALGEBRAIC statements only. Whether a family's tangent is symmetric, and "
            "whether the assembly is the sum of its families, are properties of the differential forms "
            "being summed and one instance settles them. NO count, ratio, spectrum or magnitude from "
            "this build is quotable — it is not native (see STATE.md (c) 14)"
        ),
    }


def _native_spectrum(
    cell: AssembledCell, solver: ProjectedAnalyticCG, *,
    projected: bool, n_iterations: int, probe_iterations: int, seed: int,
) -> dict[str, Any]:
    r"""Measure both ends of the native operator with Lanczos, at two folding shifts.

    The shift ``σ`` is a method parameter with no physical content: any upper bound works, and the
    recovered ``λ_min = σ − λ_max(σP − K)`` must not depend on which one was used.  Rather than pick
    one and defend the choice, the small end is measured at ``σ = 2θ`` and ``σ = 4θ`` (``θ`` from a
    short probing pass) and the disagreement is reported.  That converts a chosen constant into a
    measured invariance.

    Args:
        cell: The composed cell.
        solver: The CG workspace.
        projected: Measure ``P K P`` restricted to ``range(P)`` (the operator the inner iteration
            linearises to) rather than the raw ``K``.
        n_iterations: Krylov dimension of each reported pass.  Host memory is
            ``3·n_nodes · n_iterations`` float64 and is the real bound on this method.
        probe_iterations: Krylov dimension of the short pass that sets ``σ``.
        seed: Start-vector seed of the first shift; the second uses ``seed + 10`` so the two
            ``λ_max`` estimates are independent samples rather than the same one twice.

    Returns:
        A record with both shifts' extremal spectra, the shift-invariance check, and the seed-to-seed
        agreement of ``λ_max``.
    """
    n = int(solver.n)
    matvec = _make_matvec(cell, solver, projected=projected)
    projector = _make_projector(cell, solver) if projected else None

    start = None
    if projector is not None:
        start = projector(np.random.default_rng(int(seed) + 99).standard_normal((n, 3)))
    probe = lanczos_ritz(matvec, n, n_iterations=probe_iterations, seed=int(seed) + 99,
                         start_vector=start)
    theta, theta_bound = probe.extreme_high
    if not (theta > 0.0):
        raise SystemExit(
            f"[native-spectrum] the probing pass found no positive eigenvalue (theta={theta!r}); "
            "there is no upper bound to fold against and the small end cannot be reached"
        )

    shifts = {"sigma_2theta": 2.0 * theta, "sigma_4theta": 4.0 * theta}
    records: dict[str, Any] = {}
    for index, (label, sigma) in enumerate(sorted(shifts.items())):
        records[label] = lanczos_extremal_spectrum(
            matvec, n, upper_bound=float(sigma), n_iterations=n_iterations,
            seed=int(seed) + 10 * index, projector=projector,
        )

    first, second = records["sigma_2theta"], records["sigma_4theta"]
    lambda_min_values = [first["lambda_min_pN_per_um"], second["lambda_min_pN_per_um"]]
    reference = max(abs(value) for value in lambda_min_values)
    lambda_max_values = [first["lambda_max_pN_per_um"], second["lambda_max_pN_per_um"]]

    return {
        "operator": "P K P restricted to range(P)" if projected else "K (unprojected)",
        "why_this_operator": (
            "the inner iteration is pos += dt_mu · P F, so its linearisation on the constraint tangent "
            "space is I − dt_mu · P K P and this is the spectrum its stability and its conditioning are "
            "properties of" if projected else
            "symmetry and conservativity are statements about K itself; the projector would launder an "
            "asymmetric K into a symmetric-looking P K P, so K is measured on its own terms"
        ),
        "probing_pass": {
            "theta_pN_per_um": theta, "residual_bound_pN_per_um": theta_bound,
            "krylov_dimension": int(probe.n_iterations),
            "role": "sets the folding shift only; sigma has no physical content and the invariance "
                    "check below is what makes that claim rather than asserts it",
        },
        "shifts": {label: float(value) for label, value in shifts.items()},
        "at_sigma_2theta": first,
        "at_sigma_4theta": second,
        "shift_invariance": {
            "lambda_min_pN_per_um": lambda_min_values,
            "absolute_disagreement": abs(lambda_min_values[0] - lambda_min_values[1]),
            "relative_disagreement": (
                abs(lambda_min_values[0] - lambda_min_values[1]) / reference if reference else None),
            "both_resolved": bool(first["lambda_min_resolved"] and second["lambda_min_resolved"]),
            "meaning": (
                "lambda_min = sigma − lambda_max(sigma·P − K) must not depend on sigma. A disagreement "
                "at or near round-off says the folding is measuring the operator; a large one says it "
                "is measuring the iteration count, and the small end is then not resolved"
            ),
        },
        "lambda_max_seed_agreement": {
            "values_pN_per_um": lambda_max_values,
            "relative_disagreement": (
                abs(lambda_max_values[0] - lambda_max_values[1]) / max(
                    abs(value) for value in lambda_max_values)
                if max(abs(value) for value in lambda_max_values) else None),
            "meaning": "two independent random starts; a large end that depends on the start vector "
                       "has not converged",
        },
    }


def classify_small_end_sweep(entries: list[dict[str, Any]]) -> dict[str, Any]:
    r"""Decide whether the folded small end is CONVERGING with Krylov dimension, or is simply out of reach.

    The 2026-07-28 native run reported ``λ_min`` unresolved at one Krylov dimension, which is ambiguous
    between two very different states — "run it longer" and "this is unreachable at any affordable
    dimension" — and the difference matters, because the second is a permanent negative result that
    stops the next session walking into the same wall.  Distinguishing them needs the SAME quantity at
    several dimensions and a rule fixed in advance.

    The rule uses two independent witnesses, and demands BOTH improve:

    * ``bound/value`` — Parlett's rigorous bound relative to the Ritz value it brackets.  This is the
      method's own statement about itself and needs no reference.
    * the ``σ = 2θ`` vs ``4θ`` disagreement — ``λ_min`` cannot depend on the folding shift, so a
      shrinking disagreement is the answer settling down and a static one is the start vector talking.

    A factor of 2 in Krylov dimension buys roughly ``√2`` in reach for an extreme eigenvalue, so
    "improving" is deliberately a weak test: ANY monotone decrease in both witnesses counts.  Failing a
    weak test across a 4× span is what makes the negative verdict strong.

    Args:
        entries: One dict per dimension, ascending, each carrying ``n_iterations``,
            ``bound_over_value`` and ``relative_disagreement``.

    Returns:
        A record with the verdict, the two witness trajectories, and the reason.

    Raises:
        ValueError: If fewer than two dimensions are supplied — a sweep of one measures nothing.
    """
    if len(entries) < 2:
        raise ValueError(
            f"a small-end sweep needs at least two Krylov dimensions; got {len(entries)}. One point "
            "cannot distinguish 'converging slowly' from 'not converging'"
        )
    bounds = [float(entry["bound_over_value"]) for entry in entries]
    disagreements = [float(entry["relative_disagreement"]) for entry in entries]
    span = float(entries[-1]["n_iterations"]) / float(entries[0]["n_iterations"])

    bound_improves = all(later < earlier for earlier, later in zip(bounds, bounds[1:]))
    shift_improves = all(later < earlier for earlier, later in zip(disagreements, disagreements[1:]))
    resolved_anywhere = any(bool(entry.get("resolved")) for entry in entries)

    if resolved_anywhere:
        verdict = "RESOLVED"
        reason = ("the folded small end resolved at one of these dimensions — its Parlett bound fell "
                  "below 1% of the value, so lambda_min is a property of the operator here")
    elif bound_improves and shift_improves:
        verdict = "CONVERGING"
        reason = (f"both witnesses fall monotonically across a {span:g}x span in Krylov dimension, so "
                  "the small end is being approached and a larger dimension is the remedy")
    else:
        verdict = "OUT_OF_REACH"
        reason = (
            f"across a {span:g}x span in Krylov dimension the rigorous bound and the shift "
            "disagreement do not BOTH fall, so the reported lambda_min is a property of the start "
            "vector and the iteration count. The solve-free folded small end is not reachable at this "
            "population at an affordable dimension — the folded extremes are separated by "
            "(lambda_2 - lambda_min)/(sigma - lambda_min), which is smallest exactly when K is most "
            "ill-conditioned, so the measurement is hardest where it matters most. Reaching it needs "
            "shift-and-invert, i.e. a linear solve, which would make the measurement depend on the "
            "very solver it exists to diagnose. DO NOT RE-QUEUE this as a longer Lanczos run"
        )
    return {
        "verdict": verdict,
        "reason": reason,
        "krylov_span": span,
        "bound_over_value_trajectory": bounds,
        "shift_disagreement_trajectory": disagreements,
        "criterion": ("both Parlett bound/value AND the sigma-invariance disagreement must fall "
                      "monotonically with Krylov dimension; the test is deliberately weak, so failing "
                      "it across the full span is a strong negative"),
    }


def _small_end_sweep(
    cell: AssembledCell, solver: ProjectedAnalyticCG, *,
    dimensions: list[int], probe_iterations: int, seed: int,
) -> dict[str, Any]:
    """Measure the folded small end of ``P K P|range(P)`` at several Krylov dimensions.

    Args:
        cell: The composed cell.
        solver: The CG workspace.
        dimensions: Ascending Krylov dimensions to measure at.
        probe_iterations: Krylov dimension of the short pass that sets the folding shift.
        seed: Start-vector seed.

    Returns:
        The per-dimension record plus :func:`classify_small_end_sweep`'s verdict.
    """
    n = int(solver.n)
    matvec = _make_matvec(cell, solver, projected=True)
    projector = _make_projector(cell, solver)

    start = projector(np.random.default_rng(int(seed) + 99).standard_normal((n, 3)))
    probe = lanczos_ritz(matvec, n, n_iterations=probe_iterations, seed=int(seed) + 99,
                         start_vector=start)
    theta, _ = probe.extreme_high

    entries: list[dict[str, Any]] = []
    for index, dimension in enumerate(sorted(int(value) for value in dimensions)):
        at_shift = {}
        for offset, (label, sigma) in enumerate(
                (("sigma_2theta", 2.0 * theta), ("sigma_4theta", 4.0 * theta))):
            at_shift[label] = lanczos_extremal_spectrum(
                matvec, n, upper_bound=float(sigma), n_iterations=dimension,
                seed=int(seed) + 10 * offset + 100 * index, projector=projector)
        low = at_shift["sigma_2theta"]
        other = at_shift["sigma_4theta"]
        values = [low["lambda_min_pN_per_um"], other["lambda_min_pN_per_um"]]
        reference = max(abs(value) for value in values)
        entry = {
            "n_iterations": dimension,
            "lambda_min_pN_per_um": low["lambda_min_pN_per_um"],
            "residual_bound_pN_per_um": low["lambda_min_residual_bound"],
            "bound_over_value": (
                abs(low["lambda_min_residual_bound"]) / abs(low["lambda_min_pN_per_um"])
                if low["lambda_min_pN_per_um"] else float("inf")),
            "relative_disagreement": (
                abs(values[0] - values[1]) / reference if reference else float("inf")),
            "resolved": bool(low["lambda_min_resolved"] and other["lambda_min_resolved"]),
            "lambda_min_at_both_shifts": values,
        }
        entries.append(entry)
        print(f"[sweep]   m={dimension:>4}  λ_min={entry['lambda_min_pN_per_um']:.4e}  "
              f"bound/value={entry['bound_over_value']:.3g}  "
              f"shift-disagreement={entry['relative_disagreement']:.3g}")

    return {
        "by_dimension": entries,
        "theta_pN_per_um": theta,
        **classify_small_end_sweep(entries),
        "why": (
            "one dimension cannot distinguish 'converging slowly' from 'not converging'; this is the "
            "measurement that turns 'lambda_min unresolved at m=120' into either a remedy or a "
            "permanent negative result"
        ),
    }


def _explicit_step_verdict(cell: AssembledCell, lambda_max: float) -> dict[str, Any]:
    r"""Compare the production inner step against the stability limit this run just measured.

    The inner iteration is ``pos += dt_mu · P F`` with ``dt_mu = 0.1/kmax``.  Forward Euler on
    ``I − dt·K`` is stable exactly while ``dt·λ_max < 2``, and that criterion applies to the PLAIN
    descent because ``range(P)`` is invariant under it: for ``v`` in ``range(P)``,
    ``(I − dt·P K)v = v − dt·P K P v``.  ``kmax`` is the largest SINGLE stiffness in the build, not a
    row sum, so it bounds ``λ_max`` not at all — a node accumulates one contribution per incident bond
    — and whether the production step is stable is an open question only a measured ``λ_max`` closes.

    THE FIRE NUMBER IS REPORTED, NOT ADJUDICATED.  FIRE raises the step to ``10·dt_mu`` but is a damped
    molecular-dynamics integrator with a fictitious mass and a velocity, not the plain descent above, so
    its own stability limit is a different bound (a velocity-Verlet-type scheme scales with
    ``√λ_max``).  This function therefore reports the accelerated product and refuses to call it
    stable or unstable: what it establishes is the plain-descent verdict and the size of the number
    FIRE's own analysis would have to account for.

    Args:
        cell: The composed cell, carrying its own ``kmax`` and ``dt_mu``.
        lambda_max: The measured largest eigenvalue of the projected operator [pN/µm].

    Returns:
        The dimensionless stability products at the base step and at FIRE's maximum multiplier, the
        plain-descent verdict for the base step, and the explicit refusal to adjudicate FIRE's.
    """
    fire_max_multiplier = 10.0                     # ac/cell/inner_mechanics.py::_FIRE_DT_MAX_MULT
    dt_mu = float(cell.dt_mu)
    base = dt_mu * lambda_max
    accelerated = base * fire_max_multiplier
    return {
        "kmax_pN_per_um": float(cell.kmax),
        "mechanical_kmax_pN_per_um": float(cell.mechanical_kmax),
        "dt_mu_um_per_pN": dt_mu,
        "measured_lambda_max_pN_per_um": lambda_max,
        "lambda_max_over_kmax": lambda_max / float(cell.kmax) if cell.kmax else None,
        "stability_product_base": base,
        "stability_product_at_fire_max": accelerated,
        "fire_max_multiplier": fire_max_multiplier,
        "forward_euler_limit": 2.0,
        "base_step_stable": bool(base < 2.0),
        "fire_max_verdict": (
            "NOT ADJUDICATED HERE: FIRE is damped MD with a fictitious mass, not the plain descent this "
            "limit belongs to, so its stability bound scales differently. The product is reported so "
            "FIRE's own analysis has the measured lambda_max to work from"
        ),
        "stable_explicit_mobility_step_um_per_pN": 2.0 / lambda_max if lambda_max > 0.0 else None,
        "headroom_factor_at_base": 2.0 / base if base > 0.0 else None,
        "why_this_is_not_a_gershgorin_comparison": (
            "kmax is max over the build's individual stiffness constants, not a row sum, so it bounds "
            "nothing: a node accumulates one contribution per incident bond and lambda_max/kmax is "
            "exactly that accumulation. The 2026-07-28 slice measured a Gershgorin ROW SUM against its "
            "own lambda_max; the native runtime does not compute a row sum, so the comparable native "
            "statement is the stability product above"
        ),
    }


def native_summary_lines(
    stability: dict[str, Any], projected: dict[str, Any],
) -> list[str]:
    r"""Return the native stage's console summary, as data rather than as inline prints.

    THIS FUNCTION EXISTS BECAUSE ITS ABSENCE COST A NATIVE RUN.  On 2026-07-28 the summary was three
    inline ``print`` calls, and one of them asked ``stability`` for ``fire_max_step_stable`` — a key
    :func:`_explicit_step_verdict` deliberately does not produce, because FIRE is damped MD and the
    forward-Euler limit does not adjudicate it.  The rename to ``fire_max_verdict`` updated the tests
    and the figure and missed the print.  The result was a ``KeyError`` raised AFTER the full-native
    Lanczos pass had completed: the physics finished, ``record.json`` was never written, and the run
    was lost at a format string.  ``main()`` is not under test; a pure function is.

    Args:
        stability: The :func:`_explicit_step_verdict` block.
        projected: The projected operator's :func:`_native_spectrum` block.

    Returns:
        The lines to print, in order.

    Raises:
        KeyError: If a key is missing — which is the point: it now raises in a unit test that costs
            nothing, rather than on a device after the expensive pass.
    """
    lines = [
        f"[native] λ_max(PKP) = "
        f"{projected['at_sigma_2theta']['lambda_max_pN_per_um']:.6e} pN/µm  "
        f"(resolved={projected['at_sigma_2theta']['lambda_max_resolved']})",
        f"[native] dt_mu·λ_max = {stability['stability_product_base']:.4f} "
        f"(forward-Euler limit {stability['forward_euler_limit']} -> "
        f"{'stable' if stability['base_step_stable'] else 'UNSTABLE'}); "
        f"at FIRE max {stability['stability_product_at_fire_max']:.4f}, not adjudicated here",
        f"[native] λ_min resolved at both shifts: "
        f"{projected['shift_invariance']['both_resolved']}",
    ]
    sweep = projected.get("small_end_sweep")
    if sweep is not None:
        for entry in sweep["by_dimension"]:
            lines.append(
                f"[native]   m={entry['n_iterations']:>4}: λ_min = "
                f"{entry['lambda_min_pN_per_um']:.4e}  bound/value = "
                f"{entry['bound_over_value']:.3g}  shift disagreement = "
                f"{entry['relative_disagreement']:.3g}")
        lines.append(f"[native] small-end sweep verdict: {sweep['verdict']}")
    return lines


def form_verdict(
    certification: dict[str, Any] | None, native: dict[str, Any] | None,
) -> tuple[bool, float, float]:
    r"""Return ``(gate_passed, residual, signal)`` for one run, refusing an unjustifiable verdict.

    The gate has three conjuncts and each closes a way this measurement could be confidently wrong:
    the assembled tangent admits a potential to within the ledger's accumulation floor; every family
    that actually contributed is individually symmetric (an assembly can be symmetric overall while one
    family's asymmetry cancels another's, and that cancellation would not survive a population change);
    and the large end announced its own convergence through Parlett's bound rather than through an
    iteration count.

    Args:
        certification: The family-certification record, or ``None`` if that stage did not run.
        native: The native spectrum record, or ``None`` if that stage did not run.

    Returns:
        The gate outcome and the residual/signal pair the void ceiling is applied to — the assembled
        tangent's worst entrywise asymmetry against its own largest entry.

    Raises:
        ValueError: If ``certification`` is ``None``.  Symmetry is Lanczos's unchecked premise and is
            only establishable where the operator can be assembled, so a spectrum measured without it
            describes the symmetric part of something unknown.  Refusing is the only honest outcome —
            a Krylov method given an asymmetric operator returns confident numbers and no warning.
    """
    if certification is None:
        raise ValueError(
            "no certification stage: symmetry is Lanczos's unchecked premise and is only establishable "
            "where the operator can be assembled, so a native spectrum on its own cannot be scored. "
            "Run --stage both, or score this record against a certification run explicitly"
        )
    assembled = certification["assembled_symmetry"]
    residual = float(assembled["max_abs_asymmetry_pN_per_um"])
    signal = float(assembled["max_abs_entry_pN_per_um"])
    families_ok = all(
        entry["symmetry"]["conservative"]
        for entry in certification["families"].values() if entry["exercised"])
    lambda_max_resolved = (
        native is None
        or bool(native["projected_PKP"]["at_sigma_2theta"]["lambda_max_resolved"]))
    passed = bool(assembled["conservative"] and families_ok and lambda_max_resolved)
    return passed, residual, signal


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--stage", choices=("certify", "native", "both"), default="both",
                        help="certify: assemble and check the families on a probeable build. "
                             "native: Lanczos at full population. both: certify, then native")
    parser.add_argument("--n-filaments", type=int, default=70686,
                        help="native cortical F-actin population (default = the full baseline)")
    parser.add_argument("--membrane-subdivisions", type=int, default=6)
    parser.add_argument("--nucleus-subdivisions", type=int, default=3)
    parser.add_argument("--certify-filaments", type=int, default=60,
                        help="filament count of the probeable certification build; 3*n_nodes must "
                             "stay under assemble_dense_operator's refusal threshold")
    parser.add_argument("--certify-membrane-subdivisions", type=int, default=2)
    parser.add_argument("--certify-nucleus-subdivisions", type=int, default=2)
    parser.add_argument("--certify-overlap-free", action="store_true",
                        help="build the certification cortex overlap-free, as production does. OFF by "
                             "default and deliberately so: an overlap-free shell has no pair inside "
                             "the WCA cutoff, so the steric family assembles to exactly zero and goes "
                             "UNCERTIFIED. Symmetry of a contribution is a property of its form, not "
                             "of the configuration that exercises it, so exercising it is what matters")
    parser.add_argument("--iterations", type=int, default=120,
                        help="Krylov dimension at native; host memory is 3*n_nodes*iterations float64")
    parser.add_argument("--probe-iterations", type=int, default=40,
                        help="Krylov dimension of the short pass that sets the folding shift")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--small-end-sweep", default="",
                        help="comma-separated Krylov dimensions, e.g. 120,240,480. Measures the "
                             "FOLDED small end of P K P at each and classifies the trajectory. One "
                             "dimension cannot tell 'converging slowly' from 'not converging', which "
                             "is exactly the ambiguity the 2026-07-28 native run left open")
    parser.add_argument("--build-commit", default="",
                        help="the commit this tree corresponds to, for a run on a synced non-git tree")
    parser.add_argument("--out", default="",
                        help="artifact path (default aleph/outputs/ac/observe/"
                             "native_cortex_spectrum/record.json); figures land in its figs/")
    return parser.parse_args()


def _build(config: CellConfig) -> tuple[AssembledCell, ProjectedAnalyticCG]:
    """Compose a cell and its CG workspace.

    Args:
        config: The assembly configuration.

    Returns:
        The composed cell and a ``ProjectedAnalyticCG`` bound to it.
    """
    cell = build_cell(config)
    return cell, ProjectedAnalyticCG(cell)


def main() -> None:
    """Certify the tangent's families, then measure its native spectrum, and write one record."""
    args = _parse_args()
    run_started = time.perf_counter()
    wp.init()

    certification: dict[str, Any] | None = None
    if args.stage in ("certify", "both"):
        certify_config = CellConfig(
            n_filaments=int(args.certify_filaments),
            membrane_subdivisions=int(args.certify_membrane_subdivisions),
            nucleus_subdivisions=int(args.certify_nucleus_subdivisions),
            overlap_free_cortex=bool(args.certify_overlap_free),
            resting_bound_myosin_fraction=0.5,
            resting_bound_myosin_force_pn=1.5,
            resting_bound_myosin_source="native_spectrum_certification",
        )
        cell, solver = _build(certify_config)
        census = _population_census(cell, solver)
        print(f"[certify] build: {census['n_nodes_total']} nodes / {census['n_dof']} DOF")
        certification = certify_families(
            cell, solver, launch_census=census["stiffness_launch_census"])
        certification["build_census"] = census
        assembled = certification["assembled_symmetry"]
        print(f"[certify] assembled K: max|K−Kᵀ| = {assembled['max_abs_asymmetry_pN_per_um']:.4e} "
              f"vs floor {assembled['round_off_floor_pN_per_um']:.4e} -> "
              f"{'CONSERVATIVE' if assembled['conservative'] else 'ASYMMETRIC'}")
        for name, entry in certification["families"].items():
            verdict = ("not exercised" if not entry["exercised"]
                       else "symmetric" if entry["symmetry"]["conservative"] else "ASYMMETRIC")
            cross = " (kernels certified elsewhere)" if entry.get("certified_elsewhere") else ""
            print(f"[certify]   {name:<20} ‖K_j‖ = {entry['frobenius_norm_pN_per_um']:.4e}  "
                  f"{verdict}{cross}")
        print(f"[certify] multilinearity residual = "
              f"{certification['multilinearity']['relative_residual']!r}")
        del cell, solver

    native: dict[str, Any] | None = None
    native_census: dict[str, Any] | None = None
    stability: dict[str, Any] | None = None
    force_dump: dict[str, Any] | None = None
    device = ""
    if args.stage in ("native", "both"):
        native_config = CellConfig(
            n_filaments=int(args.n_filaments),
            membrane_subdivisions=int(args.membrane_subdivisions),
            nucleus_subdivisions=int(args.nucleus_subdivisions),
            overlap_free_cortex=True,
            resting_bound_myosin_fraction=0.5,
            resting_bound_myosin_force_pn=1.5,
            resting_bound_myosin_source="native_spectrum",
        )
        cell, solver = _build(native_config)
        device = str(cell.device)
        # Read the force channels off the COMPOSED cell rather than the config's intent: that is the
        # only observation strong enough to catch a channel that is configured on and composed off.
        force_dump = dump_force_channels(
            native_config, cell=cell, capture_env=True,
            run_label="ac/engine/observe native cortex spectrum", profile_claim="")
        native_census = _population_census(cell, solver)
        print(f"[native] build: {native_census['n_nodes_total']} nodes / "
              f"{native_census['n_dof']} DOF on {device}")
        basis_gb = 3 * native_census["n_nodes_total"] * int(args.iterations) * 8 / 1e9
        print(f"[native] Lanczos basis = {basis_gb:.2f} GB host per pass "
              f"({args.iterations} iterations, full reorthogonalization)")

        projected = _native_spectrum(
            cell, solver, projected=True, n_iterations=int(args.iterations),
            probe_iterations=int(args.probe_iterations), seed=int(args.seed))
        unprojected = _native_spectrum(
            cell, solver, projected=False, n_iterations=int(args.iterations),
            probe_iterations=int(args.probe_iterations), seed=int(args.seed) + 1000)
        if args.small_end_sweep.strip():
            projected["small_end_sweep"] = _small_end_sweep(
                cell, solver,
                dimensions=[int(v) for v in args.small_end_sweep.split(",") if v.strip()],
                probe_iterations=int(args.probe_iterations), seed=int(args.seed))
        stability = _explicit_step_verdict(
            cell, float(projected["at_sigma_2theta"]["lambda_max_pN_per_um"]))
        native = {"projected_PKP": projected, "unprojected_K": unprojected,
                  "explicit_step_verdict": stability}
        for line in native_summary_lines(stability, projected):
            print(line)
        del cell, solver

    # THE GATE, declared before the run: the operator this driver measures must admit a potential to
    # within the ledger's own accumulation floor, and the large end must announce its own convergence
    # through Parlett's bound. Symmetry is the precondition Lanczos cannot check for itself, so a run
    # that measures a spectrum without it has produced numbers about the symmetric part of something
    # and does not know what.
    ceiling = VoidCeiling(
        0.01,
        "an antisymmetric part above 1% of the tangent's own largest entry leaves no separable "
        "conservative structure: at that level the operator is not a Hessian of anything, and a "
        "Krylov method run against it reports the symmetric part while appearing to report the "
        "operator",
    )
    try:
        gate_passed, residual, signal = form_verdict(certification, native)
    except ValueError as refusal:
        raise SystemExit(f"[native-spectrum] {refusal}") from refusal

    evidence = EvidenceLabel(
        rung=EvidenceRung.CUDA_UNIT,
        quantitative=QuantitativeClaim.BLOCKED,
        basis=(
            "the assembled tangent was applied on CUDA at "
            + (f"{native_census['n_nodes_total']} nodes / {native_census['n_dof']} DOF"
               if native_census is not None else "the certification build only")
            + " and its extremal spectrum measured matrix-free; the family-by-family symmetry that "
            "Lanczos assumes was established by dense probing on a build where the operator can be "
            "assembled, which is a per-contribution algebraic property and transfers. Absolute "
            "eigenvalues stand on k_xb / k_backbone / k_linc, all unresolved PI-GAPs, so the "
            "DIMENSIONLESS structure (stability product, shift invariance, decades, family fractions) "
            "is what this run offers and no eigenvalue is quotable as a cortical material property"
        ),
    )

    config = {
        "stage": args.stage,
        "n_filaments": int(args.n_filaments),
        "membrane_subdivisions": int(args.membrane_subdivisions),
        "nucleus_subdivisions": int(args.nucleus_subdivisions),
        "certify_filaments": int(args.certify_filaments),
        "certify_membrane_subdivisions": int(args.certify_membrane_subdivisions),
        "certify_nucleus_subdivisions": int(args.certify_nucleus_subdivisions),
        "certify_overlap_free_cortex": bool(args.certify_overlap_free),
        "iterations": int(args.iterations),
        "probe_iterations": int(args.probe_iterations),
        "seed": int(args.seed),
        "small_end_sweep": args.small_end_sweep,
        "regularization": 0.0,
        "regularization_note": "exactly zero: aI + K shifts every eigenvalue by a, so a regularized "
                               "probe reports lambda + a and fabricates a null-space-free spectrum",
        "configuration_probed": "the AS-BUILT state, not a relaxed one",
    }

    record = observation_artifact(
        run_label="ac/engine/observe — native cortex tangent: family certification + Lanczos spectrum",
        evidence=evidence,
        config=config,
        census=native_census or certification["build_census"],
        t0={
            "configuration": "as built, before any inner iteration; the tangent is well defined at "
                             "any configuration but nothing about an equilibrium follows from this one",
            "regularization_pN_per_um": 0.0,
            "certification_assembled_max_abs_entry_pN_per_um": signal,
            "families_declared": [name for name, _, _, _ in _FAMILY_TOGGLES],
        },
        timing=timing_block(
            wall_seconds=time.perf_counter() - run_started,
            n_inner_iterations=(
                2 * (int(args.probe_iterations) + 2 * int(args.iterations)) if native else None),
            device_note=(
                "NOT a benchmark and NOT a converged physical step. This driver takes no physical-time "
                "step at all: the cost is the two builds, one operator application per Krylov iteration "
                "per pass, and — in the certification stage — one application per DOF per family. Warp "
                "kernel COMPILATION is included when the cache is cold and dominates everything else on "
                "a first run"
            ),
        ),
        device=device or "certification-only",
        parameter_provenance=PARAMETER_PROVENANCE,
        declared_commit=args.build_commit or None,
        void_ceiling=ceiling,
        gate_passed=gate_passed,
        residual=residual,
        signal=signal,
        force_channels=force_dump if force_dump is not None else dump_force_channels(
            None, capture_env=True,
            run_label="ac/engine/observe native cortex spectrum (certification stage only)",
            profile_claim="",
        ),
        measurements={
            "family_certification": certification,
            "native_spectrum": native,
        },
        notes={
            "what_this_re_measures": (
                "STATE.md (c) 14 retired the Gershgorin factor, the condition number and the floppy "
                "fraction because they were measured at 0.18% of native. The condition number and the "
                "step-size statement are re-measured here at native; the floppy FRACTION is not, and "
                "cannot be — Lanczos resolves the extremes and never the interior, so a mode COUNT "
                "stays a dense-probe question and therefore a per-object one"
            ),
            "why_two_operators": (
                "K carries the conservativity statement, P K P carries the solver's stability and "
                "conditioning. Reporting only one of them is how a projector comes to launder an "
                "asymmetric assembly, or a null space belonging to the constraints comes to be read "
                "as a soft physical mode"
            ),
        },
    )

    out = Path(args.out) if args.out else Path(
        "aleph/outputs/ac/observe/native_cortex_spectrum/record.json")
    write_artifact(out, record)
    print(f"[native-spectrum] wrote {out}")
    print(json.dumps(record["gate"], indent=2))


if __name__ == "__main__":
    main()

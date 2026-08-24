r"""The outer-physical / inner-mechanical driver for the composed Active Cell — one physical clock (P5).

Wires the :class:`ac.cell.assemble.AssembledCell` into the :class:`ac.fluid.scheduler.PhysicalScheduler`:

  * OUTER step (physical seconds, ``dt_phys``): membrane water-flux BC → sub-cycled conservative Biot
    p/mass drainage → the INNER mechanical solve → conservative live-domain remap. The
    physical clock comes from the process, NOT a drag-scaled descent step (the retired 6πηR clock).
  * INNER solve (dimensionless overdamped relaxation, NOT physical time — I0-A): each iteration SUMS the
    §1.4 force primitives into one per-node accumulator and, by default, advances the reused explicit
    overdamped step ``x += dt_mu·P·F`` (``dt_mu = 0.1/kmax``), with the exact NF2007 constraint-force
    projector every iteration and periodic NF2007 reshape for finite-step drift. The optional diagnostic
    accelerated paths propose either a convergence-gated analytic-PCG displacement, a live fiber-block
    preconditioned descent, a depth-one Anderson extrapolation of that descent, an exact vector-block solve on
    every overlapping WCA contact pair, or an RKC super-step. Tournament paths generate their declared
    proposals from the same state. The unchanged state, matched explicit update, and every accelerated trial
    compete under the exact projected nonlinear residual. The default max-force objective therefore remains
    monotone; the opt-in smooth L2-squared diagnostic objective can cross a max-norm plateau, while the C-2
    acceptance predicate still reads the unchanged maximum projected force. The summed
    primitives are: Cytosim bending + Hookean crosslink link-spring (reused ff/ kernels the engine KEEPS)
    + MyosinForce (I3) + StericForce (I2b) + PressureCoupling (I1a). NO lumped mechanism is launched
    (myosin_kernel / turgor_kernel / soft_contact_kernel / nucleus_shell_kernel are retired by non-use).

``--from-resting`` runs ONE outer physical step and checks both numerical integrity and actual projected
mechanical convergence.  Uniform Π₀ has zero bulk gradient but non-zero membrane traction, so an initializer
that merely remains finite is not called stable: if the membrane/cortex is not equilibrated at the
Young-Laplace preload, NG-3 and ``inner_converged`` remain false.

Runtime: Warp-CUDA only (I0-A). Runs on the gbook A5000; the dev Mac cannot launch these kernels.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import PI_0_PA, AssembledCell, CellConfig, build_cell
from aleph.components.incumbent.contact_schwarz import (
    ContactPairSchwarz,
    FiberContactSchwarz,
    RigidFiberContactSchwarz,
)
from aleph.components.incumbent.erm_gauss_seidel import ERMGaussSeidel
from aleph.components.incumbent.erm_schwarz import ERMPairSchwarz
from aleph.components.incumbent.fsi_coupling import SolidDilatationCoupling
from aleph.world.gpu_memory import query_process_peak
from aleph.components.incumbent.implicit_mechanics import (
    ProjectedAnalyticCG,
    commit_better_line_search_trial_kernel,
    compute_regularization_kernel,
    conditional_copy_vec3_kernel,
    conditional_displacement_kernel,
    copy_validity_kernel,
    decide_better_line_search_trial_kernel,
    finalize_line_search_kernel,
    omitted_regularization_base,
)
from aleph.components.incumbent.inner_mechanics import (
    begin_inner_attempt_kernel,
    compute_descent_step_kernel,
    conditional_axpy_kernel,
    conditional_reshape_kernel,
    conditional_rollback_vec3_kernel,
    finite_vec3_kernel,
    fire_adapt_kernel,
    fire_mix_and_step_kernel,
    fire_power_kernel,
    fire_reset_state_kernel,
    fire_sumsq_kernel,
    inner_state_init_kernel,
    invalidate_convergence_kernel,
    max_abs_pressure_kernel,
    max_constraint_error_kernel,
    max_displacement_kernel,
    max_force_kernel,
    project_constraint_forces_kernel,
    record_convergence_history_kernel,
    reset_active_reduction_kernel,
    update_convergence_kernel,
)
from aleph.components.incumbent.nonlinear_acceleration import (
    AndersonDepthOne,
    displacement_between_kernel,
    rkc1_recurrence_kernel,
)
from aleph.components.fluid.scheduler import DeviceInnerSolveReport, PhysicalScheduler
from aleph.components.weave.branch_angle import ARP23_K_THETA, ARP23_THETA0_RAD
from aleph.components.weave.branch_angle_warp import branch_angle_kernel
from aleph.laws.forces_warp import cytosim_bending_kernel
from aleph.laws.network_warp import _zero, link_spring_kernel

__all__ = ["StabilityReport", "run_from_resting", "make_inner_solve", "step_myosin_kinetics"]


def step_myosin_kinetics(
    cell: AssembledCell, nearest_dist: wp.array, nearest_idx: wp.array, tau: float, seed: int,
) -> None:
    """One NMII KMC tick WITH the I4 walk_dir hand-off: attach → fill_walk_dir → step+detach (device).

    Replaces ``MyosinForce.step_kinetics`` at the lead level so ``fill_walk_dir_kernel`` runs right AFTER the
    I3 ``attach_kernel`` (motor INTEGRATION.md §1, HARD): a newly-bound head's ``walk_dir`` is overwritten from
    the ACTUAL actin filament's barbed-end polarity (``unit(pos[barbed] − pos[anchor])``), so the power stroke
    walks along the actin it grabbed, not the frozen minifilament-axis default — the directed contraction is
    physically correct. ``compute_loads(pos)`` MUST be called before this (the step uses the current loads; a
    freshly-attached head is unloaded ⇒ its first stroke is v0, correct). RNG seeds are decorrelated between the
    attach and step/detach streams (matching ``step_kinetics``).
    """
    from aleph.components.motor.hand import attach_kernel, step_detach_kernel
    from aleph.components.weave.presets import fill_walk_dir_kernel

    myo = cell.myosin
    if myo.segment_runtime is not None:
        # Diagnostic/manual callers use an always-accepted scalar; production calls the same runtime only from
        # ``commit_irreversible`` with the scheduler-owned final outer predicate.
        myo.commit_segment_kinetics(
            cell.pos_d,
            myo.segment_runtime.accepted_diagnostic,
            tau,
            seed,
        )
        return
    st = myo.state
    n = int(st["bound"].shape[0])
    seed_attach = wp.int32(seed & 0x7FFFFFFF)
    seed_detach = wp.int32((seed ^ 0x5BD1E995) & 0x7FFFFFFF)
    with wp.ScopedDevice(cell.device):
        wp.launch(attach_kernel, dim=n,
                  inputs=[st["bound"], st["anchor"], st["abscissa"], nearest_dist, nearest_idx,
                          myo.params, wp.float64(tau), seed_attach])
        wp.launch(fill_walk_dir_kernel, dim=n,
                  inputs=[st["bound"], st["anchor"], cell.node_fiber_d, cell.barbed_node_d,
                          cell.pos_d, st["walk_dir"]])
        wp.launch(step_detach_kernel, dim=n,     # F4/F5 split: tangential (Hill) + full |F| (Bell) loads
                  inputs=[st["bound"], st["anchor"], st["abscissa"], myo.loads, myo.loads_full,
                          myo.params, wp.float64(tau), seed_detach])


#: Force-channel names ``_accumulate_all`` will omit on request, in the manifest's vocabulary
#: (`ac/engine/forces_manifest.py::_FORCE_TERMS`).  ONE vocabulary on purpose: a second spelling here
#: would let a typo omit nothing silently, which is precisely the double-count this parameter exists to
#: prevent.
OMITTABLE_CHANNELS: frozenset[str] = frozenset({
    "actin_bending_cytosim",
    "actin_crosslink_link_spring",
    "arp23_branch_angle",
    "nmii_minifilament_force",
    "nucleus_compartment_force",
    "membrane_compartment_force",
    "membrane_pressure_traction",
    "steric_wca",
    "biot_pressure_coupling",
})


@wp.kernel
def _sum_force_squared_kernel(force: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.float64)) -> None:
    """Smooth line-search objective ``sum_i |P F_i|^2`` [pN^2]; never an acceptance predicate."""
    i = wp.tid()
    wp.atomic_add(out, 0, wp.length_sq(force[i]))


def _accumulate_all(
    cell: AssembledCell,
    pos: wp.array,
    f: wp.array,
    *,
    omit: frozenset[str] | set[str] | tuple[str, ...] = (),
    additional_force_accumulators: tuple[Callable[[wp.array, wp.array], None], ...] = (),
) -> None:
    """Zero ``f`` then SUM every composed §1.4 primitive into it (the inner force-assembly contract).

    ``omit`` names channels this call must NOT launch, because an engine component now owns them and
    computes them on its OWN arrays.  Without it there is no way to hand one component to the engine
    while the rest of the cell keeps the incumbent path: the component's force would be added twice, and
    the only alternative — build-time ``--no-myosin``-style flags — removes the physics instead of
    relocating it.  Four validation tracks were blocked behind this parameter (T2; PI decision
    `COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md` §6 D2,
    approved 2026-07-29, option (a): additive, bit-identical default, per-component).

    **The default is bit-identical.** With ``omit`` empty every launch is reached exactly as before, in
    the same order, so a run that does not opt in cannot change.

    **What omitting does NOT do.** It does not disable the physics — it asserts that someone else
    computes it. Nothing here can check that claim, which is why T2's negative control is to leave the
    mask OFF while the engine owner is bound and require parity to fail by an exact factor of two: the
    double count must be *detectable*, not silently absorbed.

    Args:
        cell: The assembled cell.
        pos: Positions to evaluate at.
        f: Force array, zeroed here.
        omit: Channel names from :data:`OMITTABLE_CHANNELS`.
        additional_force_accumulators: Experiment-scoped additive force launchers called with the exact
            candidate ``(pos, f)`` after every incumbent channel. The empty default is launch-order
            identical. This is an apparatus/instrumentation seam, not a new biological force owner.

    Raises:
        ValueError: If a name is not a known channel.  A typo must not be read as "omit nothing" — that
            failure mode is invisible and produces exactly the doubled force this guards against.
    """
    omit = frozenset(omit)
    unknown = omit - OMITTABLE_CHANNELS
    if unknown:
        raise ValueError(
            f"unknown force channel(s) in omit={sorted(unknown)}; known: {sorted(OMITTABLE_CHANNELS)}. "
            "Refusing to treat an unrecognised name as 'omit nothing' — that reads as success and "
            "doubles the force the caller believes it moved to the engine"
        )
    d = cell.device
    wp.launch(_zero, dim=cell.n_total, inputs=[f], device=d)
    if cell.n_tri and "actin_bending_cytosim" not in omit:    # Cytosim bending (reused ff/ kernel; actin)
        wp.launch(cytosim_bending_kernel, dim=cell.n_tri, inputs=[pos, cell.tri_d, cell.alpha_d, f], device=d)
    if cell.n_xl and "actin_crosslink_link_spring" not in omit:  # Hookean crosslink link-spring (ff/ kernel)
        wp.launch(link_spring_kernel, dim=cell.n_xl,
                  inputs=[pos, cell.xl_d, cell.kxl_d, cell.r0xl_d, f], device=d)
    if (cell.branch_triples_d is not None and cell.n_branch
            and "arp23_branch_angle" not in omit):        # Arp2/3 70° angle-harmonic (mixed cortex only)
        # SAME kernel/arg-order/units as LamellipodiumBranchAngleMechanics.accumulate (protrusion.py) — the
        # Faessler-2020 θ₀/k_θ anchors. Gated on branch_triples_d so the formin-only default launches nothing.
        wp.launch(branch_angle_kernel, dim=cell.n_branch,
                  inputs=[pos, cell.branch_triples_d, cell.branch_active_d,
                          ARP23_THETA0_RAD, ARP23_K_THETA, f], device=d)
    if cell.myosin is not None and "nmii_minifilament_force" not in omit:   # MyosinForce (I3)
        cell.myosin.accumulate(pos, f)
    if cell.nucleus is not None and "nucleus_compartment_force" not in omit:  # NucleusCompartment (P4)
        cell.nucleus.accumulate(pos, f)
    if cell.membrane is not None and "membrane_compartment_force" not in omit:  # MembraneCompartment (P4)
        cell.membrane.accumulate(pos, f)
    if cell.membrane_pressure is not None and "membrane_pressure_traction" not in omit:   # NG-3
        cell.membrane_pressure.accumulate(pos, f)
    if cell.steric is not None and "steric_wca" not in omit:      # StericForce (I2b) — all-fiber WCA
        cell.steric.accumulate(cell.state, f)
    if cell.pressure is not None and "biot_pressure_coupling" not in omit:  # PressureCoupling (I1a)
        cell.pressure.accumulate(cell.state, f)
    for accumulator in additional_force_accumulators:
        accumulator(pos, f)


def _residual_host(cell: AssembledCell, pos: wp.array, f: wp.array) -> tuple[float, np.ndarray]:
    """Post-loop whole-cell force readback (never called from the physical/inner hot loop)."""
    _accumulate_all(cell, pos, f)
    wp.synchronize_device(cell.device)
    fh = f.numpy()
    mag = np.linalg.norm(fh, axis=1)
    return float(np.nanmax(mag)) if mag.size else 0.0, mag


def _preload_erm_resting_balance(cell: AssembledCell, sweeps: int = 2) -> dict:
    """Pre-stretch ERM tethers so the resting membrane turgor is held at ``t=0`` (physiological baseline).

    At build each ERM rest length equals its formation length, so every tether is force-free and the
    membrane carries the full turgor pressure (~ΔP·A_node per node, ~40 pN here) unbalanced. Physiologically
    the actomyosin cortex bears cortical tension that holds turgor, and the membrane rides on it through the
    ERM linkers; a relaxed cortex plus a positive turgor is not a valid resting baseline (PI
    physiological-baseline rule) — there is no force-balanced equilibrium there, which is why no inner solver
    converges from it.

    This sets each tether's rest length ``rest = |gap| + f_in/k_erm`` where ``f_in`` is the membrane node's
    net force projected onto the inward (membrane→cortex) direction, so the tether pre-tension exactly cancels
    its membrane node's net outward force. The equal-and-opposite reaction loads the cortex, which — being
    stiff (crosslink k ~ 8e5 pN/µm) — settles in negligible motion, whereas the soft membrane would otherwise
    have to travel ~ f/k_erm ≈ 9 nm through a CFL throttled by that same cortex stiffness. Every quantity is
    derived from the live pressure field and the sourced ``k_erm`` — nothing is tuned to make a gate pass.

    The pairing is one tether per membrane node, so no force-sharing correction is needed. Two sweeps make the
    balance self-consistent after the tether force enters the accumulation.
    """
    mem = getattr(cell, "membrane", None)
    if mem is None or getattr(mem, "n_erm", 0) == 0:
        return {"n_erm": 0, "res_before": None, "res_after": None}
    k = float(mem.k_erm)
    erm_m = mem.erm_m_d.numpy().astype(np.int64)
    erm_c = mem.erm_c_d.numpy().astype(np.int64)
    res_before, _ = _residual_host(cell, cell.pos_d, cell.f_d)
    # PER-SWEEP TRACE. The fixed point claims that adding f_in/k to rest drives f_in to zero, so the
    # sweeps should show |f_in| collapsing. Measured 2026-08-10 it does not: rest moves ~8 nm with no
    # clipping, the force is almost entirely ALONG the tether axis (p50 37.3 of 39.4 pN, so it is the
    # component a tether can cancel), and the residual is unchanged to fifteen digits. Recording each
    # sweep is what distinguishes "the correction is not applied" from "the correction is applied and
    # the force does not respond", and only the second is a statement about the model.
    sweep_trace = []
    for _ in range(max(1, sweeps)):
        cell.f_d.zero_()
        _accumulate_all(cell, cell.pos_d, cell.f_d)
        wp.synchronize_device(cell.device)
        f = cell.f_d.numpy()
        pos = cell.pos_d.numpy()
        d = pos[erm_c] - pos[erm_m]                       # membrane -> cortex (inward)
        length = np.linalg.norm(d, axis=1)
        ok = length > 1e-12
        u = np.zeros_like(d)
        u[ok] = d[ok] / length[ok, None]
        f_in = np.einsum("ij,ij->i", f[erm_m], u)         # membrane net force along inward (neg if outward)
        # Incremental fixed point: adjusting rest by δ changes the tether force by -k·δ, so rest += f_in/k
        # drives the residual membrane force (which already includes the current tether) toward zero. Using
        # the current rest (not the fixed gap) is essential — otherwise a second sweep, seeing a now-balanced
        # membrane, would reset rest back to the force-free formation length and undo the first sweep.
        rest_cur = mem.erm_rest_d.numpy()
        rest_want = rest_cur + f_in / k
        rest_new = np.clip(rest_want, 0.0, length)  # never a compressive (pushing) tether at rest
        mem.erm_rest_d.assign(np.ascontiguousarray(rest_new, dtype=np.float64))
        sweep_trace.append({
            "f_in_p50_pN": float(np.percentile(f_in, 50)),
            "f_in_absmax_pN": float(np.abs(f_in).max()),
            "rest_p50_um": float(np.percentile(rest_cur, 50)),
            "rest_new_p50_um": float(np.percentile(rest_new, 50)),
            "extension_p50_um": float(np.percentile(length - rest_new, 50)),
            "implied_tether_tension_p50_pN": float(k * np.percentile(length - rest_new, 50)),
        })
    res_after, _ = _residual_host(cell, cell.pos_d, cell.f_d)
    # WHY THE PRELOAD DID OR DID NOT MOVE THE RESIDUAL. Measured 2026-08-10: it reports
    # res_before == res_after to fifteen digits on a full native build, i.e. it runs and changes
    # nothing, and the same signature appears in the gate-b driver's log where res_after came out
    # HIGHER. A preload that silently no-ops is worse than one that is absent, because the ledger
    # entry reads as evidence that the resting baseline was established. These diagnostics say which
    # of the three ways it can fail actually happened: the correction was too small to matter, the
    # clip refused it, or it was applied to nodes whose force is not what the residual is made of.
    clipped_lo = int(np.count_nonzero(rest_want < 0.0))
    clipped_hi = int(np.count_nonzero(rest_want > length))
    delta = rest_new - rest_cur
    # THE DECOMPOSITION THAT SETTLES IT. A tether can only pull along its own axis, so it can cancel
    # the component of a membrane node's force ALONG that axis and nothing else. Measured: the preload
    # changes rest by ~8 nm with no clipping and moves the residual by zero, which means the 47.4 pN
    # is not in the direction the tether acts. Splitting the post-preload membrane force into its
    # along-tether and perpendicular parts says so directly instead of by elimination.
    _f_post = cell.f_d.numpy()
    _f_mem = _f_post[erm_m]
    _along = np.einsum("ij,ij->i", _f_mem, u)
    _perp = np.linalg.norm(_f_mem - _along[:, None] * u, axis=1)
    _mag = np.linalg.norm(_f_mem, axis=1)
    residual_direction = {
        "membrane_node_force_pN": {"p50": float(np.percentile(_mag, 50)),
                                   "max": float(_mag.max())},
        "along_tether_abs_pN": {"p50": float(np.percentile(np.abs(_along), 50)),
                                "max": float(np.abs(_along).max())},
        "perpendicular_pN": {"p50": float(np.percentile(_perp, 50)), "max": float(_perp.max())},
        "perp_fraction_of_magnitude_p50": float(np.percentile(_perp / np.maximum(_mag, 1e-300), 50)),
        "why": "a tether acts only along its axis; whatever sits perpendicular cannot be preloaded away",
    }
    return {
        "n_erm": int(erm_m.shape[0]),
        "res_before": float(res_before), "res_after": float(res_after),
        "residual_unchanged": bool(res_before == res_after),
        "f_in_pN": {"min": float(f_in.min()), "p50": float(np.percentile(f_in, 50)),
                    "max": float(f_in.max())},
        "rest_delta_um": {"min": float(delta.min()), "p50": float(np.percentile(delta, 50)),
                          "max": float(delta.max())},
        "n_clipped_at_zero": clipped_lo, "n_clipped_at_length": clipped_hi,
        "k_erm_pN_per_um": k,
        "erm_tether_length_um": {"min": float(length.min()), "p50": float(np.percentile(length, 50)),
                                 "max": float(length.max())},
        "residual_direction": residual_direction,
        "sweep_trace": sweep_trace,
    }


def _preload_cortex_pretension(cell: AssembledCell, prestrain: float) -> dict:
    """Give the cortex shell hoop tension by pre-straining its crosslinks (physiological Laplace pre-stress).

    A resting cell holds its turgor ΔP by cortical tension ``T = ΔP·R/2`` (Young-Laplace); a force-free
    cortex (crosslink rest = formation length) carries no tension and therefore cannot balance the membrane
    turgor — the root cause of the ~40 pN/node membrane residual. Scaling every crosslink rest length by
    ``(1 - prestrain)`` pre-stretches the network so it carries isotropic in-plane (hoop) tension. The
    prestrain is calibrated to the derived Laplace tension, never tuned to make a gate pass; the cortex is
    the physiological turgor-bearing structure (cortical-tension measurements ARE cortex tension).

    Crosslink turnover would relax this pre-stress over physical time; that is the physiological creep the
    active myosin maintains against. During the resting inner solve turnover only commits at accepted outer
    steps, so the pre-stress holds while the baseline is found.
    """
    n_xl = int(getattr(cell, "n_xl", 0))
    if n_xl == 0 or prestrain == 0.0:
        return {"n_xl": n_xl, "prestrain": 0.0}
    r0 = cell.r0xl_d.numpy()
    cell.r0xl_d.assign(np.ascontiguousarray(r0 * (1.0 - float(prestrain)), dtype=np.float64))
    return {"n_xl": n_xl, "prestrain": float(prestrain)}


def make_inner_solve(
    cell: AssembledCell,
    n_inner: int,
    reshape_every: int = 20,
    max_inner_retries: int = 0,
    *,
    capture_candidate: bool = False,
    inner_solver: str = "explicit",
    implicit_cg_max_iterations: int = 32,
    implicit_line_search_steps: int = 4,
    implicit_coarse_iterations: int = 8,
    implicit_coarse_modes: int = 0,
    rkc_stages: int = 4,
    line_search_objective: str = "max_force",
    omit: frozenset[str] | set[str] | tuple[str, ...] = (),
    additional_force_accumulators: tuple[Callable[[wp.array, wp.array], None], ...] = (),
    additional_stiffness_pn_per_um: float = 0.0,
):
    """Return device-converged mechanics with fixed-launch, device-predicated numerical retry chunks.

    ``omit`` and ``additional_force_accumulators`` are forwarded verbatim to every
    :func:`_accumulate_all` this solve performs, so a component
    that owns a channel can take it over INSIDE the inner loop rather than only in an out-of-loop probe.
    Without the passthrough the mask reaches nothing that iterates, and `omit=` would be usable only by
    callers that assemble forces once — which is not where a cortex channel lives.  Default empty:
    bit-identical, and a run that does not opt in cannot change.  (PI `COMPARTMENT_VALIDATION_TRACKS` §6 D2,
    2026-07-29.) Additional force accumulators are restricted to experiment apparatus and own no irreversible
    biology. ``additional_stiffness_pn_per_um`` expands the explicit stability scale by the declared apparatus
    tangent; adding an accumulator without that bound is refused.

    ``n_inner`` is the nonlinear iteration count per chunk and ``max_inner_retries`` is an explicit compute
    budget. ``analytic_implicit`` solves ``(aI + P K_analytic P) dx = P F``; ``block_descent`` proposes the
    source-independent ``P M^-1 P F`` direction; ``anderson`` applies depth-one residual-minimizing Anderson
    acceleration to that block map. ``contact_schwarz`` applies the exact degree-weighted overlapping 6x6
    WCA-pair block solve, and ``rkc1`` composes the already-validated explicit step through the first-order
    Chebyshev recurrence. ``tournament`` retains the prior three-candidate benchmark;
    ``contact_tournament`` adds the contact-block proposal. A fixed device-resident geometric line search tests
    ``1, 1/2, ...`` times every active
    direction against the exact nonlinear projected residual and retains the global best among the unchanged
    state, matched explicit-CFL update, and accelerated trials. Candidate selection defaults to maximum force;
    opt-in ``l2_squared`` changes only that competition objective, never the convergence predicate. Force
    assembly, convergence, reshape,
    transaction, and physical time semantics remain identical. A
    retry continues the same quasi-static candidate inside one outer transaction; it neither re-runs the fluid
    update nor advances physical time.  The host launches every configured chunk, while ``active_d`` makes all
    chunks after convergence or invalidation no-ops without a device-to-host decision. ``capture_candidate`` is
    a diagnostic-only D2D snapshot taken immediately before rollback; it is disabled in every production path
    and therefore adds neither storage nor a copy to the physical-time loop by default.
    """
    if n_inner <= 0:
        raise ValueError("n_inner is a positive maximum-iteration budget")
    if reshape_every <= 0:
        raise ValueError("reshape_every must be positive")
    if max_inner_retries < 0:
        raise ValueError("max_inner_retries must be nonnegative")
    solver_choices = {
        "explicit", "fire", "analytic_implicit", "block_descent", "anderson", "rkc1", "contact_schwarz",
        "fiber_contact_schwarz", "rigid_contact_schwarz", "tournament", "contact_tournament",
        "cluster_tournament", "rigid_cluster_tournament", "erm_schwarz", "erm_tournament",
        "augmented_block", "augmented_tournament",
        "erm_jacobi_block", "erm_jacobi_tournament",
        "erm_jacobi_pure", "erm_jacobi_pure_tournament",
        "erm_gauss_seidel", "erm_gauss_seidel_tournament",
    }
    if inner_solver not in solver_choices:
        raise ValueError(f"inner_solver must be one of {sorted(solver_choices)}")
    if implicit_cg_max_iterations <= 0:
        raise ValueError("implicit_cg_max_iterations must be positive")
    if implicit_line_search_steps <= 0:
        raise ValueError("implicit_line_search_steps must be positive")
    if implicit_coarse_iterations < 0:
        raise ValueError("implicit_coarse_iterations must be nonnegative")
    if not 0 <= implicit_coarse_modes <= 12:
        raise ValueError("implicit_coarse_modes must be between 0 and 12")
    if rkc_stages < 2:
        raise ValueError("rkc_stages must be at least 2")
    if line_search_objective not in {"max_force", "l2_squared"}:
        raise ValueError("line_search_objective must be 'max_force' or 'l2_squared'")
    additional_force_accumulators = tuple(additional_force_accumulators)
    if any(not callable(accumulator) for accumulator in additional_force_accumulators):
        raise TypeError("every additional force accumulator must be callable")
    additional_stiffness_pn_per_um = float(additional_stiffness_pn_per_um)
    if not math.isfinite(additional_stiffness_pn_per_um) or additional_stiffness_pn_per_um < 0.0:
        raise ValueError("additional stiffness must be finite and nonnegative [pN/um]")
    if additional_force_accumulators and additional_stiffness_pn_per_um <= 0.0:
        raise ValueError(
            "experiment force accumulators require a positive declared stiffness bound; refusing an "
            "unscaled force in the explicit stability step"
        )
    mechanical_stiffness_bound = float(cell.mechanical_kmax) + additional_stiffness_pn_per_um

    max_attempts = max_inner_retries + 1
    max_total_inner = n_inner * max_attempts
    # FIRE, like the plain explicit descent, proposes a single self-contained update and never builds a
    # tangent/tournament workspace; both are excluded from the accelerated (line-search) path.
    accelerated_solver = inner_solver not in {"explicit", "fire"}
    candidate_families = {
        "explicit": (),
        "fire": (),
        "analytic_implicit": ("analytic_implicit",),
        "block_descent": ("block_descent",),
        "anderson": ("anderson",),
        "rkc1": ("rkc1",),
        "contact_schwarz": ("contact_schwarz",),
        "fiber_contact_schwarz": ("fiber_contact_schwarz",),
        "rigid_contact_schwarz": ("rigid_contact_schwarz",),
        "tournament": ("block_descent", "anderson", "rkc1"),
        "contact_tournament": ("block_descent", "anderson", "contact_schwarz", "rkc1"),
        "cluster_tournament": ("block_descent", "anderson", "fiber_contact_schwarz", "rkc1"),
        "rigid_cluster_tournament": ("block_descent", "anderson", "rigid_contact_schwarz", "rkc1"),
        "erm_schwarz": ("erm_schwarz",),
        "erm_tournament": ("block_descent", "anderson", "erm_schwarz", "rkc1"),
        # Backbone-aware ERM-augmented fiber block (the fiber Cholesky block extended with its ERM-tethered
        # membrane DOFs). block_descent's ``P M^-1 P F`` direction is now the augmented-block descent.
        "augmented_block": ("block_descent",),
        "augmented_tournament": ("block_descent", "anderson", "rkc1"),
        # Tension-side ERM membrane Jacobi (no block augmentation): plain fiber block for actin +
        # k_erm-scaled scalar-Jacobi membrane step so the coupled membrane->ERM->cortex mode converges.
        "erm_jacobi_block": ("block_descent",),
        "erm_jacobi_tournament": ("block_descent", "anderson", "rkc1"),
        # Pure per-node scalar-Jacobi (fiber block DISABLED so the actin over-step cannot throttle the
        # shared line-search t) + tension-side membrane diagonal.
        "erm_jacobi_pure": ("block_descent",),
        "erm_jacobi_pure_tournament": ("block_descent", "anderson", "rkc1"),
        # ★ Multiplicative (coloured symmetric block Gauss-Seidel) sweep over the crosslink+ERM+membrane graph:
        # coordinates the stiff crosslinks so the soft-membrane step cannot detonate them (the additive
        # erm_schwarz over-steps; this does not). Pure = the GS direction alone; tournament adds the standard
        # anderson/rkc extrapolators of the explicit step.
        "erm_gauss_seidel": ("erm_gauss_seidel",),
        "erm_gauss_seidel_tournament": ("block_descent", "anderson", "erm_gauss_seidel", "rkc1"),
    }[inner_solver]
    acceptance_slots = len(candidate_families) * implicit_line_search_steps

    d = cell.device
    with wp.ScopedDevice(d):
        outer_start_d = wp.empty_like(cell.pos_d)
        pos_prev_d = wp.empty_like(cell.pos_d)
        max_disp_d = wp.zeros(1, dtype=wp.float64)
        max_constraint_d = wp.zeros(1, dtype=wp.float64)
        convergence_force_d = wp.zeros(1, dtype=wp.float64)
        residual_d = wp.zeros(1, dtype=wp.float64)
        active_d = wp.zeros(1, dtype=wp.int32)
        converged_d = wp.zeros(1, dtype=wp.int32)
        finite_d = wp.ones(1, dtype=wp.int32)
        iters_d = wp.zeros(1, dtype=wp.int32)
        attempts_d = wp.zeros(1, dtype=wp.int32)
        max_pressure_d = wp.zeros(1, dtype=wp.float64)
        candidate_pos_d = wp.empty_like(cell.pos_d) if capture_candidate else None
        candidate_preprojection_pos_d = wp.empty_like(cell.pos_d) if capture_candidate else None
        history_slots = max_total_inner // reshape_every + 2 if capture_candidate else 0
        history_iteration_d = wp.zeros(history_slots, dtype=wp.int32) if capture_candidate else None
        history_displacement_d = wp.zeros(history_slots, dtype=wp.float64) if capture_candidate else None
        history_constraint_d = wp.zeros(history_slots, dtype=wp.float64) if capture_candidate else None
        history_projected_force_d = wp.zeros(history_slots, dtype=wp.float64) if capture_candidate else None
        projected_f_d = wp.empty_like(cell.f_d)
        projector_diag_d = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64)
        projector_rhs_d = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64)
        dt_mu_d = wp.array(np.array([cell.dt_mu], dtype=np.float64), dtype=wp.float64, device=d)
        valid_default_d = wp.ones(1, dtype=wp.int32, device=d)
        is_fire = inner_solver == "fire"
        fire_velocity_d = wp.zeros(cell.n_total, dtype=wp.vec3d) if is_fire else None
        fire_dt_mult_d = wp.zeros(1, dtype=wp.float64) if is_fire else None
        fire_alpha_d = wp.zeros(1, dtype=wp.float64) if is_fire else None
        fire_n_positive_d = wp.zeros(1, dtype=wp.int32) if is_fire else None
        fire_reset_d = wp.zeros(1, dtype=wp.int32) if is_fire else None
        fire_power_d = wp.zeros(1, dtype=wp.float64) if is_fire else None
        fire_vel_sumsq_d = wp.zeros(1, dtype=wp.float64) if is_fire else None
        fire_force_sumsq_d = wp.zeros(1, dtype=wp.float64) if is_fire else None
        implicit_regularization_d = wp.zeros(1, dtype=wp.float64) if accelerated_solver else None
        base_projected_f_d = wp.empty_like(cell.f_d) if accelerated_solver else None
        explicit_candidate_d = wp.empty_like(cell.pos_d) if accelerated_solver else None
        best_candidate_d = wp.empty_like(cell.pos_d) if accelerated_solver else None
        implicit_residual_d = wp.zeros(1, dtype=wp.float64) if accelerated_solver else None
        explicit_residual_d = wp.zeros(1, dtype=wp.float64) if accelerated_solver else None
        best_residual_d = wp.zeros(1, dtype=wp.float64) if accelerated_solver else None
        base_finite_d = wp.ones(1, dtype=wp.int32) if accelerated_solver else None
        implicit_finite_d = wp.ones(1, dtype=wp.int32) if accelerated_solver else None
        explicit_finite_d = wp.ones(1, dtype=wp.int32) if accelerated_solver else None
        best_finite_d = wp.ones(1, dtype=wp.int32) if accelerated_solver else None
        take_trial_d = wp.zeros(1, dtype=wp.int32) if accelerated_solver else None
        best_trial_index_d = wp.zeros(1, dtype=wp.int32) if accelerated_solver else None
        implicit_accept_count_d = wp.zeros(1, dtype=wp.int32) if accelerated_solver else None
        explicit_accept_count_d = wp.zeros(1, dtype=wp.int32) if accelerated_solver else None
        stationary_accept_count_d = wp.zeros(1, dtype=wp.int32) if accelerated_solver else None
        implicit_scale_accept_counts_d = (
            wp.zeros(acceptance_slots, dtype=wp.int32)
            if accelerated_solver else None
        )
        implicit_line_search_scales_d = (
            tuple(
                wp.array(np.array([2.0 ** (-index)], dtype=np.float64), dtype=wp.float64, device=d)
                for index in range(implicit_line_search_steps)
            )
            if accelerated_solver else ()
        )
        has_block_candidate = inner_solver in {
            "block_descent", "anderson", "tournament", "contact_tournament", "cluster_tournament",
            "rigid_cluster_tournament", "erm_tournament", "augmented_block", "augmented_tournament",
            "erm_jacobi_block", "erm_jacobi_tournament",
            "erm_jacobi_pure", "erm_jacobi_pure_tournament",
            "erm_gauss_seidel_tournament",
        }
        has_contact_candidate = inner_solver in {"contact_schwarz", "contact_tournament"}
        has_erm_candidate = inner_solver in {"erm_schwarz", "erm_tournament"}
        has_gs_candidate = inner_solver in {"erm_gauss_seidel", "erm_gauss_seidel_tournament"}
        has_fiber_contact_candidate = inner_solver in {"fiber_contact_schwarz", "cluster_tournament"}
        has_rigid_contact_candidate = inner_solver in {
            "rigid_contact_schwarz", "rigid_cluster_tournament",
        }
        has_rkc_candidate = inner_solver in {
            "rkc1", "tournament", "contact_tournament", "cluster_tournament",
            "rigid_cluster_tournament", "erm_tournament", "augmented_tournament",
            "erm_jacobi_tournament", "erm_jacobi_pure_tournament",
            "erm_gauss_seidel_tournament",
        }
        block_authorized_d = wp.ones(1, dtype=wp.int32) if has_block_candidate else None
        anderson_authorized_d = (
            wp.ones(1, dtype=wp.int32)
            if inner_solver in {
                "tournament", "contact_tournament", "cluster_tournament", "rigid_cluster_tournament",
                "erm_tournament", "augmented_tournament", "erm_jacobi_tournament",
                "erm_jacobi_pure_tournament", "erm_gauss_seidel_tournament",
            } else None
        )
        contact_authorized_d = wp.ones(1, dtype=wp.int32) if has_contact_candidate else None
        erm_authorized_d = wp.ones(1, dtype=wp.int32) if has_erm_candidate else None
        gs_authorized_d = wp.ones(1, dtype=wp.int32) if has_gs_candidate else None
        fiber_contact_authorized_d = (
            wp.ones(1, dtype=wp.int32) if has_fiber_contact_candidate else None
        )
        rigid_contact_authorized_d = (
            wp.ones(1, dtype=wp.int32) if has_rigid_contact_candidate else None
        )
        rkc_authorized_d = wp.ones(1, dtype=wp.int32) if has_rkc_candidate else None
        analytic_authorized_d = wp.ones(1, dtype=wp.int32) if inner_solver == "analytic_implicit" else None
        rkc_previous_d = wp.empty_like(cell.pos_d) if has_rkc_candidate else None
        rkc_next_d = wp.empty_like(cell.pos_d) if has_rkc_candidate else None
        rkc_direction_d = wp.empty_like(cell.pos_d) if has_rkc_candidate else None
    tangent_accelerator = inner_solver in {
        "analytic_implicit", "block_descent", "anderson", "tournament", "contact_tournament",
        "cluster_tournament", "rigid_cluster_tournament", "erm_tournament",
        "augmented_block", "augmented_tournament",
        "erm_jacobi_block", "erm_jacobi_tournament",
        "erm_jacobi_pure", "erm_jacobi_pure_tournament",
        "erm_gauss_seidel_tournament",
    }
    augment_fiber_block_erm = inner_solver in {"augmented_block", "augmented_tournament"}
    erm_tension_side_precond = inner_solver in {
        "erm_jacobi_block", "erm_jacobi_tournament", "erm_jacobi_pure", "erm_jacobi_pure_tournament"}
    disable_fiber_block = inner_solver in {"erm_jacobi_pure", "erm_jacobi_pure_tournament"}
    analytic_workspace = (
        ProjectedAnalyticCG(
            cell,
            max_iterations=implicit_cg_max_iterations,
            coarse_iterations=implicit_coarse_iterations,
            coarse_modes=implicit_coarse_modes,
            augment_fiber_block_erm=augment_fiber_block_erm,
            erm_tension_side_precond=erm_tension_side_precond,
            disable_fiber_block=disable_fiber_block,
        )
        if tangent_accelerator else None
    )
    implicit_cg = analytic_workspace if inner_solver == "analytic_implicit" else None
    anderson = (
        AndersonDepthOne(cell.n_total, d)
        if inner_solver in {
            "anderson", "tournament", "contact_tournament", "cluster_tournament",
            "rigid_cluster_tournament", "erm_tournament", "augmented_tournament",
            "erm_jacobi_tournament", "erm_jacobi_pure_tournament",
            "erm_gauss_seidel_tournament",
        } else None
    )
    contact_schwarz = ContactPairSchwarz(cell) if has_contact_candidate else None
    fiber_contact_schwarz = FiberContactSchwarz(cell) if has_fiber_contact_candidate else None
    rigid_contact_schwarz = RigidFiberContactSchwarz(cell) if has_rigid_contact_candidate else None
    erm_schwarz = ERMPairSchwarz(cell) if has_erm_candidate else None
    erm_gs = ERMGaussSeidel(cell) if has_gs_candidate else None
    implicit_omitted_base = (
        omitted_regularization_base(cell)
        if analytic_workspace is not None or contact_schwarz is not None
        or fiber_contact_schwarz is not None or rigid_contact_schwarz is not None
        or erm_schwarz is not None or erm_gs is not None else 0.0
    )
    solid_coupling = SolidDilatationCoupling(cell.grid) if cell.grid is not None else None
    solid_velocity_d = solid_coupling.allocate_velocity(cell.n_total) if solid_coupling is not None else None
    tolerance_um = float(np.sqrt(np.finfo(np.float64).eps) * cell.convergence_length_um)
    current_dt_phys = 0.0

    def _check_after_projected_step(iteration: int) -> None:
        wp.launch(reset_active_reduction_kernel, dim=1, inputs=[max_disp_d, active_d], device=d)
        wp.launch(reset_active_reduction_kernel, dim=1, inputs=[max_constraint_d, active_d], device=d)
        wp.launch(max_displacement_kernel, dim=cell.n_total,
                  inputs=[pos_prev_d, cell.pos_d], outputs=[max_disp_d], device=d)
        wp.launch(finite_vec3_kernel, dim=cell.n_total, inputs=[cell.pos_d, finite_d], device=d)
        if cell.n_fibers:
            wp.launch(max_constraint_error_kernel, dim=cell.n_fibers,
                      inputs=[cell.pos_d, cell.foff_d, cell.soff_d, cell.srest_d],
                      outputs=[max_constraint_d], device=d)
        _evaluate_projected_residual(finite_d, convergence_force_d)
        wp.launch(update_convergence_kernel, dim=1,
                  inputs=[max_disp_d, max_constraint_d, convergence_force_d, dt_mu_d,
                          wp.float64(tolerance_um), wp.int32(iteration), finite_d, active_d,
                          converged_d, iters_d],
                  device=d)
        if history_iteration_d is not None:
            slot = iteration // reshape_every
            wp.launch(
                record_convergence_history_kernel,
                dim=1,
                inputs=[
                    wp.int32(slot), wp.int32(iteration), max_disp_d, max_constraint_d,
                    convergence_force_d, history_iteration_d, history_displacement_d,
                    history_constraint_d, history_projected_force_d,
                ],
                device=d,
            )

    def _reduce_projected_metric(
        candidate_force_d: wp.array, candidate_residual_d: wp.array, *, selection: bool,
    ) -> None:
        candidate_residual_d.zero_()
        if selection and line_search_objective == "l2_squared":
            wp.launch(_sum_force_squared_kernel, dim=cell.n_total,
                      inputs=[candidate_force_d], outputs=[candidate_residual_d], device=d)
        else:
            wp.launch(max_force_kernel, dim=cell.n_total,
                      inputs=[candidate_force_d], outputs=[candidate_residual_d], device=d)

    def _evaluate_projected_residual(
        candidate_finite_d: wp.array, candidate_residual_d: wp.array, *, selection: bool = False,
    ) -> None:
        """Evaluate one trial geometry with the exact runtime force assembly and NF2007 projector."""
        _accumulate_all(
            cell,
            cell.pos_d,
            cell.f_d,
            omit=omit,
            additional_force_accumulators=additional_force_accumulators,
        )
        wp.launch(finite_vec3_kernel, dim=cell.n_total, inputs=[cell.f_d, candidate_finite_d], device=d)
        wp.copy(projected_f_d, cell.f_d)
        if cell.n_fibers:
            wp.launch(project_constraint_forces_kernel, dim=cell.n_fibers,
                      inputs=[cell.pos_d, cell.f_d, cell.foff_d, cell.soff_d, projected_f_d,
                              projector_diag_d, projector_rhs_d, candidate_finite_d], device=d)
        wp.launch(finite_vec3_kernel, dim=cell.n_total,
                  inputs=[projected_f_d, candidate_finite_d], device=d)
        _reduce_projected_metric(projected_f_d, candidate_residual_d, selection=selection)

    def _rkc_direction(candidate_finite_d: wp.array) -> wp.array:
        """Generate one RKC direction from ``pos_prev_d`` without sharing validity with other candidates."""
        wp.copy(cell.pos_d, pos_prev_d)
        wp.copy(rkc_previous_d, pos_prev_d)
        wp.launch(conditional_axpy_kernel, dim=cell.n_total,
                  inputs=[cell.pos_d, dt_mu_d, base_projected_f_d, active_d], device=d)
        wp.launch(finite_vec3_kernel, dim=cell.n_total,
                  inputs=[cell.pos_d, candidate_finite_d], device=d)
        for _ in range(2, rkc_stages + 1):
            _accumulate_all(
                cell,
                cell.pos_d,
                cell.f_d,
                omit=omit,
                additional_force_accumulators=additional_force_accumulators,
            )
            wp.launch(finite_vec3_kernel, dim=cell.n_total,
                      inputs=[cell.f_d, candidate_finite_d], device=d)
            wp.copy(projected_f_d, cell.f_d)
            if cell.n_fibers:
                wp.launch(
                    project_constraint_forces_kernel,
                    dim=cell.n_fibers,
                    inputs=[
                        cell.pos_d, cell.f_d, cell.foff_d, cell.soff_d, projected_f_d,
                        projector_diag_d, projector_rhs_d, candidate_finite_d,
                    ],
                    device=d,
                )
            wp.launch(rkc1_recurrence_kernel, dim=cell.n_total,
                      inputs=[rkc_previous_d, cell.pos_d, projected_f_d, dt_mu_d,
                              active_d, candidate_finite_d, rkc_next_d], device=d)
            wp.copy(rkc_previous_d, cell.pos_d)
            wp.copy(cell.pos_d, rkc_next_d)
            wp.launch(finite_vec3_kernel, dim=cell.n_total,
                      inputs=[cell.pos_d, candidate_finite_d], device=d)
        wp.launch(displacement_between_kernel, dim=cell.n_total,
                  inputs=[pos_prev_d, cell.pos_d, rkc_direction_d], device=d)
        return rkc_direction_d

    def inner_solve(dt_phys: float, fluid_valid_d: wp.array | None = None) -> DeviceInnerSolveReport:
        nonlocal current_dt_phys
        current_dt_phys = float(dt_phys)
        # ``dt_phys`` is used only to convert the accepted outer displacement to solid velocity for Biot;
        # it never scales or reinterprets the quasi-static inner iterations themselves.
        pos, f = cell.pos_d, cell.f_d
        precondition_d = valid_default_d if fluid_valid_d is None else fluid_valid_d
        with wp.ScopedDevice(d):
            wp.copy(outer_start_d, pos)
            wp.launch(inner_state_init_kernel, dim=1,
                      inputs=[active_d, converged_d, iters_d, attempts_d, finite_d,
                              wp.int32(max_total_inner), precondition_d], device=d)
            if implicit_accept_count_d is not None:
                implicit_accept_count_d.zero_()
                explicit_accept_count_d.zero_()
                stationary_accept_count_d.zero_()
                implicit_scale_accept_counts_d.zero_()
            if anderson is not None:
                anderson.reset()
            if is_fire:
                fire_velocity_d.zero_()
                wp.launch(fire_reset_state_kernel, dim=1,
                          inputs=[fire_dt_mult_d, fire_alpha_d, fire_n_positive_d, fire_reset_d], device=d)
            if cell.grid is not None and cell.membrane_pressure is not None:
                max_pressure_d.zero_()
                wp.launch(max_abs_pressure_kernel, dim=cell.grid.shape,
                          inputs=[cell.grid.p, cell.grid.mask, wp.float64(cell.membrane_pressure.p_ext)],
                          outputs=[max_pressure_d], device=d)
                wp.launch(compute_descent_step_kernel, dim=1,
                          inputs=[wp.float64(mechanical_stiffness_bound), max_pressure_d,
                                  wp.float64(cell.pressure_edge_um)],
                          outputs=[dt_mu_d], device=d)
            if implicit_regularization_d is not None:
                wp.launch(compute_regularization_kernel, dim=1,
                          inputs=[wp.float64(implicit_omitted_base), max_pressure_d,
                                  wp.float64(cell.pressure_edge_um), implicit_regularization_d], device=d)
            for attempt in range(max_attempts):
                wp.launch(begin_inner_attempt_kernel, dim=1, inputs=[active_d, attempts_d], device=d)
                for local_it in range(n_inner):
                    iteration = attempt * n_inner + local_it + 1
                    wp.copy(pos_prev_d, pos)
                    _accumulate_all(
                        cell,
                        pos,
                        f,
                        omit=omit,
                        additional_force_accumulators=additional_force_accumulators,
                    )
                    wp.copy(projected_f_d, f)
                    if cell.n_fibers:
                        wp.launch(
                            project_constraint_forces_kernel,
                            dim=cell.n_fibers,
                            inputs=[
                                pos, f, cell.foff_d, cell.soff_d, projected_f_d,
                                projector_diag_d, projector_rhs_d, finite_d,
                            ],
                            device=d,
                        )
                    if inner_solver == "explicit":
                        wp.launch(conditional_axpy_kernel, dim=cell.n_total,
                                  inputs=[pos, dt_mu_d, projected_f_d, active_d], device=d)
                    elif is_fire:
                        # FIRE momentum step on the SAME projected force P F: power sign -> adaptive
                        # dt/damping -> velocity mixing + semi-implicit Euler. No tangent, no line search;
                        # the strict projected-force convergence gate below is unchanged.
                        fire_power_d.zero_()
                        wp.launch(fire_power_kernel, dim=cell.n_total,
                                  inputs=[projected_f_d, fire_velocity_d, fire_power_d], device=d)
                        fire_vel_sumsq_d.zero_()
                        wp.launch(fire_sumsq_kernel, dim=cell.n_total,
                                  inputs=[fire_velocity_d, fire_vel_sumsq_d], device=d)
                        fire_force_sumsq_d.zero_()
                        wp.launch(fire_sumsq_kernel, dim=cell.n_total,
                                  inputs=[projected_f_d, fire_force_sumsq_d], device=d)
                        wp.launch(fire_adapt_kernel, dim=1,
                                  inputs=[fire_power_d, fire_dt_mult_d, fire_alpha_d,
                                          fire_n_positive_d, fire_reset_d, active_d], device=d)
                        wp.launch(fire_mix_and_step_kernel, dim=cell.n_total,
                                  inputs=[pos, fire_velocity_d, projected_f_d, dt_mu_d,
                                          fire_dt_mult_d, fire_alpha_d, fire_vel_sumsq_d,
                                          fire_force_sumsq_d, fire_reset_d, active_d], device=d)
                    else:
                        wp.copy(base_projected_f_d, projected_f_d)
                        wp.launch(copy_validity_kernel, dim=1,
                                  inputs=[finite_d, base_finite_d], device=d)
                        candidate_directions: list[tuple[wp.array, wp.array]] = []
                        if implicit_cg is not None:
                            wp.launch(copy_validity_kernel, dim=1,
                                      inputs=[base_finite_d, analytic_authorized_d], device=d)
                            displacement_d = implicit_cg.solve(
                                pos, base_projected_f_d, implicit_regularization_d,
                                analytic_authorized_d)
                            candidate_directions.append((displacement_d, implicit_cg.converged))
                        if has_block_candidate:
                            wp.copy(pos, pos_prev_d)
                            wp.launch(copy_validity_kernel, dim=1,
                                      inputs=[base_finite_d, block_authorized_d], device=d)
                            block_direction_d = analytic_workspace.preconditioned_direction(
                                pos, base_projected_f_d, implicit_regularization_d, block_authorized_d)
                            if inner_solver in {
                                "block_descent", "tournament", "contact_tournament", "cluster_tournament",
                                "rigid_cluster_tournament", "augmented_block", "augmented_tournament",
                                "erm_jacobi_block", "erm_jacobi_tournament",
                                "erm_jacobi_pure", "erm_jacobi_pure_tournament",
                            }:
                                candidate_directions.append((block_direction_d, block_authorized_d))
                            if anderson is not None:
                                if anderson_authorized_d is not None:
                                    wp.launch(copy_validity_kernel, dim=1,
                                              inputs=[block_authorized_d, anderson_authorized_d], device=d)
                                    anderson_finite_d = anderson_authorized_d
                                else:
                                    anderson_finite_d = block_authorized_d
                                anderson_direction_d = anderson.direction(
                                    pos, block_direction_d, active_d, anderson_finite_d)
                                candidate_directions.append((anderson_direction_d, anderson_finite_d))
                        if contact_schwarz is not None:
                            wp.copy(pos, pos_prev_d)
                            wp.launch(copy_validity_kernel, dim=1,
                                      inputs=[base_finite_d, contact_authorized_d], device=d)
                            contact_direction_d = contact_schwarz.direction(
                                pos, base_projected_f_d, implicit_regularization_d, active_d,
                                contact_authorized_d)
                            candidate_directions.append((contact_direction_d, contact_authorized_d))
                        if erm_schwarz is not None:
                            wp.copy(pos, pos_prev_d)
                            wp.launch(copy_validity_kernel, dim=1,
                                      inputs=[base_finite_d, erm_authorized_d], device=d)
                            erm_direction_d = erm_schwarz.direction(
                                pos, base_projected_f_d, implicit_regularization_d, active_d,
                                erm_authorized_d)
                            candidate_directions.append((erm_direction_d, erm_authorized_d))
                        if erm_gs is not None:
                            wp.copy(pos, pos_prev_d)
                            wp.launch(copy_validity_kernel, dim=1,
                                      inputs=[base_finite_d, gs_authorized_d], device=d)
                            gs_direction_d = erm_gs.direction(
                                pos, base_projected_f_d, implicit_regularization_d, active_d,
                                gs_authorized_d)
                            candidate_directions.append((gs_direction_d, gs_authorized_d))
                        if fiber_contact_schwarz is not None:
                            wp.copy(pos, pos_prev_d)
                            wp.launch(copy_validity_kernel, dim=1,
                                      inputs=[base_finite_d, fiber_contact_authorized_d], device=d)
                            fiber_contact_direction_d = fiber_contact_schwarz.direction(
                                pos, base_projected_f_d, implicit_regularization_d, active_d,
                                fiber_contact_authorized_d)
                            candidate_directions.append(
                                (fiber_contact_direction_d, fiber_contact_authorized_d))
                        if rigid_contact_schwarz is not None:
                            wp.copy(pos, pos_prev_d)
                            wp.launch(copy_validity_kernel, dim=1,
                                      inputs=[base_finite_d, rigid_contact_authorized_d], device=d)
                            rigid_contact_direction_d = rigid_contact_schwarz.direction(
                                pos, base_projected_f_d, implicit_regularization_d, active_d,
                                rigid_contact_authorized_d)
                            candidate_directions.append(
                                (rigid_contact_direction_d, rigid_contact_authorized_d))
                        if has_rkc_candidate:
                            wp.launch(copy_validity_kernel, dim=1,
                                      inputs=[base_finite_d, rkc_authorized_d], device=d)
                            candidate_directions.append(
                                (_rkc_direction(rkc_authorized_d), rkc_authorized_d))
                        # The unchanged state initializes the competition, making the nonlinear residual truly
                        # monotone. The explicit-CFL state is then tested as the ordinary progress candidate.
                        # Accelerated directions are not clipped to its displacement envelope; fixed geometric
                        # backtracking globalizes every proposal against the same exact nonlinear residual.
                        _reduce_projected_metric(base_projected_f_d, best_residual_d, selection=True)
                        wp.copy(best_candidate_d, pos_prev_d)
                        wp.copy(best_finite_d, base_finite_d)
                        best_trial_index_d.fill_(-2)

                        wp.copy(pos, pos_prev_d)
                        wp.launch(copy_validity_kernel, dim=1,
                                  inputs=[base_finite_d, explicit_finite_d], device=d)
                        wp.launch(conditional_axpy_kernel, dim=cell.n_total,
                                  inputs=[pos, dt_mu_d, base_projected_f_d, active_d], device=d)
                        wp.copy(explicit_candidate_d, pos)
                        wp.launch(finite_vec3_kernel, dim=cell.n_total,
                                  inputs=[pos, explicit_finite_d], device=d)
                        _evaluate_projected_residual(
                            explicit_finite_d, explicit_residual_d, selection=True)
                        wp.launch(decide_better_line_search_trial_kernel, dim=1,
                                  inputs=[explicit_residual_d, best_residual_d, explicit_finite_d,
                                          best_finite_d, base_finite_d, active_d, take_trial_d], device=d)
                        wp.launch(conditional_copy_vec3_kernel, dim=cell.n_total,
                                  inputs=[explicit_candidate_d, best_candidate_d, take_trial_d], device=d)
                        wp.launch(commit_better_line_search_trial_kernel, dim=1,
                                  inputs=[explicit_residual_d, explicit_finite_d, wp.int32(-1), take_trial_d,
                                          best_residual_d, best_finite_d, best_trial_index_d], device=d)

                        for family_index, (displacement_d, direction_authorized_d) in enumerate(
                                candidate_directions):
                            for scale_index, trial_scale_d in enumerate(implicit_line_search_scales_d):
                                trial_index = family_index * implicit_line_search_steps + scale_index
                                wp.copy(pos, pos_prev_d)
                                wp.launch(copy_validity_kernel, dim=1,
                                          inputs=[base_finite_d, implicit_finite_d], device=d)
                                wp.launch(conditional_displacement_kernel, dim=cell.n_total,
                                          inputs=[pos, displacement_d, trial_scale_d, active_d,
                                                  implicit_finite_d], device=d)
                                wp.launch(finite_vec3_kernel, dim=cell.n_total,
                                          inputs=[pos, implicit_finite_d], device=d)
                                _evaluate_projected_residual(
                                    implicit_finite_d, implicit_residual_d, selection=True)
                                wp.launch(decide_better_line_search_trial_kernel, dim=1,
                                          inputs=[implicit_residual_d, best_residual_d, implicit_finite_d,
                                                  best_finite_d, direction_authorized_d, active_d,
                                                  take_trial_d], device=d)
                                wp.launch(conditional_copy_vec3_kernel, dim=cell.n_total,
                                          inputs=[pos, best_candidate_d, take_trial_d], device=d)
                                wp.launch(commit_better_line_search_trial_kernel, dim=1,
                                          inputs=[implicit_residual_d, implicit_finite_d,
                                                  wp.int32(trial_index), take_trial_d, best_residual_d,
                                                  best_finite_d, best_trial_index_d], device=d)
                        wp.copy(pos, best_candidate_d)
                        wp.launch(finalize_line_search_kernel, dim=1,
                                  inputs=[best_finite_d, best_trial_index_d, finite_d,
                                          implicit_accept_count_d,
                                          explicit_accept_count_d, stationary_accept_count_d,
                                          implicit_scale_accept_counts_d], device=d)
                    if cell.n_fibers and iteration % reshape_every == 0:
                        wp.launch(
                            conditional_reshape_kernel,
                            dim=cell.n_fibers,
                            inputs=[pos, cell.foff_d, cell.soff_d, cell.srest_d, wp.int32(2), active_d],
                            device=d,
                        )
                        _check_after_projected_step(iteration)
                    elif not cell.n_fibers:
                        _check_after_projected_step(iteration)

                # Retry chunks partition only the compute budget; they must not change the numerical trajectory.
                # The stronger closeout reshape therefore runs once at the FINAL configured boundary, not at
                # every intermediate chunk. Each boundary still checks acceptance, so an early convergence latch
                # device-disables all later launches without a host decision.
                if cell.n_fibers and attempt == max_attempts - 1:
                    if candidate_preprojection_pos_d is not None:
                        wp.copy(candidate_preprojection_pos_d, pos)
                    wp.launch(
                        conditional_reshape_kernel,
                        dim=cell.n_fibers,
                        inputs=[pos, cell.foff_d, cell.soff_d, cell.srest_d, wp.int32(4), active_d],
                        device=d,
                    )
                _check_after_projected_step((attempt + 1) * n_inner)

            _accumulate_all(
                cell,
                pos,
                f,
                omit=omit,
                additional_force_accumulators=additional_force_accumulators,
            )
            wp.launch(finite_vec3_kernel, dim=cell.n_total, inputs=[f, finite_d], device=d)
            wp.launch(invalidate_convergence_kernel, dim=1,
                      inputs=[finite_d, active_d, converged_d], device=d)
            residual_d.zero_()
            wp.launch(max_force_kernel, dim=cell.n_total, inputs=[f], outputs=[residual_d], device=d)
            if candidate_pos_d is not None:
                wp.copy(candidate_pos_d, pos)
            # Candidate diagnostics above are retained, but a rejected solve cannot mutate the authoritative
            # mechanical state seen by domain remap or the next physical step.
            wp.launch(conditional_rollback_vec3_kernel, dim=cell.n_total,
                      inputs=[pos, outer_start_d, converged_d], device=d)
        return DeviceInnerSolveReport(
            residual_d=residual_d, iters_d=iters_d, attempts_d=attempts_d, converged_d=converged_d,
            max_displacement_d=max_disp_d, constraint_residual_d=max_constraint_d,
            finite_d=finite_d, dt_mu_d=dt_mu_d)

    def post_remap(dt_phys: float, accepted_d: wp.array) -> None:
        """Spread accepted solid motion after the scheduler installs the new live-domain mask."""
        if solid_coupling is not None and solid_velocity_d is not None:
            solid_coupling.update_from_displacement(
                outer_start_d, cell.pos_d, dt_phys, cell.node_volume_d, cell.solid_active_d,
                solid_velocity_d, accepted_d)

    def rollback(accepted_d: wp.array) -> None:
        """Restore the outer-start geometry when any scheduler-owned acceptance latch rejects the candidate."""
        wp.launch(conditional_rollback_vec3_kernel, dim=cell.n_total,
                  inputs=[cell.pos_d, outer_start_d, accepted_d], device=d)

    def commit_irreversible(accepted_d: wp.array) -> None:
        """Launch irreversible biology once, using only the device outer-acceptance predicate."""
        if cell.myosin is not None and cell.myosin.segment_runtime is not None:
            # ``0x4E4D4949`` is only an RNG namespace tag (ASCII "NMII"), not a physical parameter.
            seed = (int(cell.cfg.seed) ^ 0x4E4D4949) & 0x7FFFFFFF
            # Final scheduler order: converged mechanics -> every outer validator -> accept/rollback ->
            # load -> live segment query/polarity -> accepted-predicated attach/step/detach.
            cell.myosin.commit_segment_kinetics(cell.pos_d, accepted_d, current_dt_phys, seed)
        if cell.nucleus is not None:
            cell.nucleus.rupture_step(cell.pos_d, accepted_d)
        if cell.membrane is not None:
            # ``0x45524D`` is only an RNG namespace tag (ASCII "ERM"), not a physical parameter.
            seed = (int(cell.cfg.seed) ^ 0x45524D) & 0x7FFFFFFF
            cell.membrane.commit_kinetics(cell.pos_d, accepted_d, current_dt_phys, seed)

    failure_counters: list[wp.array] = []
    if cell.membrane_pressure is not None:
        failure_counters.append(cell.membrane_pressure.unresolved_faces_d)
    if cell.membrane_bc is not None:
        failure_counters.append(cell.membrane_bc.unresolved_faces_d)
    if cell.membrane_mask_provider is not None:
        failure_counters.append(cell.membrane_mask_provider.query_failures_d)
    if cell.nucleus_mask_provider is not None:
        failure_counters.append(cell.nucleus_mask_provider.query_failures_d)

    # Scheduler-owned ordering hook: mechanics -> live remap -> div(v_s) on the NEW mask. The callback carries
    # no host state and performs only prescribed CUDA launches.
    inner_solve.post_remap = post_remap  # type: ignore[attr-defined]
    inner_solve.rollback = rollback  # type: ignore[attr-defined]
    inner_solve.commit_irreversible = commit_irreversible  # type: ignore[attr-defined]
    inner_solve.failure_counters = tuple(failure_counters)  # type: ignore[attr-defined]
    inner_solve.max_attempts = max_attempts  # type: ignore[attr-defined]
    inner_solve.iteration_budget = max_total_inner  # type: ignore[attr-defined]
    inner_solve.candidate_pos_d = candidate_pos_d  # type: ignore[attr-defined]
    inner_solve.candidate_preprojection_pos_d = candidate_preprojection_pos_d  # type: ignore[attr-defined]
    inner_solve.history_iteration_d = history_iteration_d  # type: ignore[attr-defined]
    inner_solve.history_displacement_d = history_displacement_d  # type: ignore[attr-defined]
    inner_solve.history_constraint_d = history_constraint_d  # type: ignore[attr-defined]
    inner_solve.history_projected_force_d = history_projected_force_d  # type: ignore[attr-defined]
    inner_solve.convergence_force_d = convergence_force_d  # type: ignore[attr-defined]
    inner_solve.inner_solver = inner_solver  # type: ignore[attr-defined]
    inner_solve.implicit_cg = implicit_cg  # type: ignore[attr-defined]
    inner_solve.accelerator_workspace = analytic_workspace  # type: ignore[attr-defined]
    inner_solve.anderson = anderson  # type: ignore[attr-defined]
    inner_solve.contact_schwarz = contact_schwarz  # type: ignore[attr-defined]
    inner_solve.fiber_contact_schwarz = fiber_contact_schwarz  # type: ignore[attr-defined]
    inner_solve.rigid_contact_schwarz = rigid_contact_schwarz  # type: ignore[attr-defined]
    inner_solve.erm_schwarz = erm_schwarz  # type: ignore[attr-defined]
    inner_solve.erm_gs = erm_gs  # type: ignore[attr-defined]
    inner_solve.rkc_stages = rkc_stages  # type: ignore[attr-defined]
    inner_solve.implicit_regularization_d = implicit_regularization_d  # type: ignore[attr-defined]
    inner_solve.implicit_residual_d = implicit_residual_d  # type: ignore[attr-defined]
    inner_solve.explicit_residual_d = explicit_residual_d  # type: ignore[attr-defined]
    inner_solve.best_residual_d = best_residual_d  # type: ignore[attr-defined]
    inner_solve.implicit_accept_count_d = implicit_accept_count_d  # type: ignore[attr-defined]
    inner_solve.explicit_accept_count_d = explicit_accept_count_d  # type: ignore[attr-defined]
    inner_solve.stationary_accept_count_d = stationary_accept_count_d  # type: ignore[attr-defined]
    inner_solve.implicit_line_search_scales = tuple(  # type: ignore[attr-defined]
        2.0 ** (-index) for index in range(implicit_line_search_steps))
    inner_solve.accelerator_candidate_families = candidate_families  # type: ignore[attr-defined]
    inner_solve.implicit_coarse_iterations = implicit_coarse_iterations  # type: ignore[attr-defined]
    inner_solve.implicit_coarse_modes = implicit_coarse_modes  # type: ignore[attr-defined]
    inner_solve.implicit_scale_accept_counts_d = implicit_scale_accept_counts_d  # type: ignore[attr-defined]
    inner_solve.line_search_objective = line_search_objective  # type: ignore[attr-defined]
    return inner_solve


@dataclass
class StabilityReport:
    """Resting --from-resting stability diagnostics (all host reads are out-of-hot-loop)."""

    ran: bool
    n_actin: int
    n_total: int
    dt_phys: float
    n_biot_subcycles: int
    inner_iters: int
    inner_attempts: int
    inner_max_attempts: int
    inner_iteration_budget: int
    inner_converged: bool
    inner_tolerance_um: float
    inner_max_displacement_um: float
    inner_constraint_residual_um: float
    inner_finite: bool
    fluid_finite: bool
    inner_dt_mu: float
    residual_start: float
    residual_candidate: float
    residual_end: float
    force_finite: bool
    com_drift_um: float
    r_mean_start: float
    r_mean_end: float
    r_min_end: float
    r_max_end: float
    pos_finite: bool
    content_start: float
    content_end: float
    content_delta: float
    membrane_surface_flux_integral: float
    membrane_grid_source_integral: float
    interior_source_integral: float
    solid_dilatation_integral: float
    moving_face_delta: float
    membrane_deposition_error: float
    solver_conservation_error: float
    physical_conservation_error: float
    membrane_flux_unresolved_faces: int
    membrane_pressure_unresolved_faces: int
    membrane_mask_query_failures: int
    membrane_mask_surface_ties: int
    nucleus_mask_query_failures: int
    nucleus_mask_surface_ties: int
    hot_loop_dtoh_copies: int
    hot_loop_dtoh_elapsed_ms: float
    profiler_positive_control_dtoh_copies: int
    outer_accepted: bool
    outer_rolled_back: bool
    committed_time_s: float
    build_wall_s: float
    precheck_wall_s: float
    outer_step_wall_s: float
    postcheck_wall_s: float
    wall_s: float
    stable: bool
    ledger: dict

    def to_dict(self) -> dict:
        return {k: (v if not isinstance(v, np.floating) else float(v)) for k, v in self.__dict__.items()}


def _radii(pos: np.ndarray) -> tuple[float, float, float, float, np.ndarray]:
    c = pos.mean(axis=0)
    r = np.linalg.norm(pos - c, axis=1)
    return float(r.mean()), float(r.min()), float(r.max()), c, r


def run_from_resting(
    cfg: CellConfig | None = None, *, n_inner: int = 60, dt_phys: float = 0.05,
    reshape_every: int = 20, max_inner_retries: int = 0,
    biot_cfl_safety: float = 0.9, profile_memcpy: bool = False,
    inner_solver: str = "explicit", implicit_cg_max_iterations: int = 32,
    implicit_line_search_steps: int = 4,
    implicit_coarse_iterations: int = 8,
    implicit_coarse_modes: int = 0,
    rkc_stages: int = 4,
    preload_erm_balance: bool = False,
    cortex_prestrain: float = 0.0,
    capture_convergence_history: bool = False,
) -> StabilityReport:
    """Build the composed cell and advance ONE outer physical step; return the resting-stability report.

    Args:
        capture_convergence_history: Record the device-side ``(iteration, max|Δx|, constraint, max|PF|)``
            checkpoint every ``reshape_every`` iterations and read it back into
            ``ledger["convergence_history"]`` AFTER the loop.  Instrumentation only — the capture is the
            already-existing :func:`~aleph.components.incumbent.inner_mechanics.record_convergence_history_kernel`
            writing to device arrays, so the hot loop stays free of host reads and the GPU-residency gate
            is unaffected.  Off by default because it also allocates two candidate position buffers.

            Why this needed plumbing at all: the solve has recorded this trajectory since it was written,
            but nothing surfaced it, so every resting record carries only ``residual_start`` and
            ``residual_end``.  Those two numbers cannot distinguish "converging slowly" from "stalled",
            which is the exact question the cost of a converged step depends on.
    """
    total_t0 = time.perf_counter()
    build_t0 = time.perf_counter()
    cell = build_cell(cfg)
    if cortex_prestrain:
        cell.ledger["cortex_pretension"] = _preload_cortex_pretension(cell, cortex_prestrain)
    if preload_erm_balance:
        cell.ledger["erm_preload"] = _preload_erm_resting_balance(cell)
    cell.ledger["inner_solver"] = inner_solver
    cell.ledger["implicit_cg_max_iterations"] = int(implicit_cg_max_iterations)
    cell.ledger["implicit_line_search_steps"] = int(implicit_line_search_steps)
    cell.ledger["implicit_coarse_iterations"] = int(implicit_coarse_iterations)
    cell.ledger["implicit_coarse_modes"] = int(implicit_coarse_modes)
    cell.ledger["rkc_stages"] = int(rkc_stages)
    build_wall = time.perf_counter() - build_t0
    precheck_t0 = time.perf_counter()
    pos0 = cell.pos_d.numpy()[: cell.n_actin].copy()
    rms0, rmin0, rmax0, com0, _ = _radii(pos0)
    res_start, _ = _residual_host(cell, cell.pos_d, cell.f_d)
    inner_tolerance_um = float(np.sqrt(np.finfo(np.float64).eps) * cell.convergence_length_um)
    precheck_wall = time.perf_counter() - precheck_t0

    n_sub, moving_delta, iters, attempts, inner_converged = 0, 0.0, 0, 0, False
    inner_max_displacement = inner_constraint_residual = 0.0
    inner_finite, fluid_finite, inner_dt_mu = True, True, cell.dt_mu
    hot_dtoh_copies = positive_dtoh_copies = 0
    hot_dtoh_elapsed = 0.0
    residual_candidate = res_start
    outer_accepted = outer_rolled_back = False
    committed_time = 0.0
    content_start = content_end = content_delta = 0.0
    membrane_surface_integral = membrane_grid_integral = 0.0
    interior_integral = dilatation_integral = moving_delta = 0.0
    membrane_deposition_error = solver_conservation_error = physical_conservation_error = 0.0
    outer_wall = 0.0
    history_owner = None
    if cell.substrate is not None and cell.domain is not None:
        # Hoisted out of the PhysicalScheduler(...) call so the convergence-history arrays it hangs on
        # itself are reachable for the post-loop readback below.
        history_owner = make_inner_solve(
            cell,
            n_inner,
            reshape_every,
            max_inner_retries,
            inner_solver=inner_solver,
            implicit_cg_max_iterations=implicit_cg_max_iterations,
            implicit_line_search_steps=implicit_line_search_steps,
            implicit_coarse_iterations=implicit_coarse_iterations,
            implicit_coarse_modes=implicit_coarse_modes,
            rkc_stages=rkc_stages,
            capture_candidate=capture_convergence_history,
        )
        sched = PhysicalScheduler(
            substrate=cell.substrate, domain=cell.domain, membrane_bc=cell.membrane_bc,
            inner_solve=history_owner,
            # Resting osmotic balance ⇒ zero net membrane flux. Read from the LEDGER, not from the module
            # constant: `cfg.turgor_pi0_pa` makes Π₀ a declared axis, and a build seeded at one Π₀ whose
            # membrane balances against another would leak flux the arms would then differ by. The ledger
            # carries the value actually built, so the two cannot drift apart.
            osmotic_difference=lambda t: float(cell.ledger["resting_pressure_Pa"]),
        )
        if profile_memcpy:
            wp.synchronize_device(cell.device)
            wp.timing_begin(wp.TIMING_MEMCPY)
        else:
            wp.synchronize_device(cell.device)
        outer_t0 = time.perf_counter()
        rep_d = sched.outer_step(dt_phys, biot_cfl_safety=biot_cfl_safety)
        if profile_memcpy:
            hot_records = wp.timing_end(synchronize=True)
            dtoh = [record for record in hot_records if record.name == "memcpy DtoH"]
            hot_dtoh_copies = len(dtoh)
            hot_dtoh_elapsed = float(sum(record.elapsed for record in dtoh))
        else:
            wp.synchronize_device(cell.device)
        outer_wall = time.perf_counter() - outer_t0
        if profile_memcpy:
            # Non-vacuity control: the same profiler must see an intentional post-loop scalar download.
            wp.timing_begin(wp.TIMING_MEMCPY)
            rep_d.inner_residual_d.numpy()
            positive_records = wp.timing_end(synchronize=True)
            positive_dtoh_copies = sum(record.name == "memcpy DtoH" for record in positive_records)
        rep = rep_d.readback()
        n_sub, moving_delta, iters = rep.n_biot_subcycles, rep.moving_face_delta, rep.inner_iters
        attempts = rep.inner_attempts
        residual_candidate, inner_converged = rep.inner_residual, rep.inner_converged
        inner_max_displacement = rep.inner_max_displacement
        inner_constraint_residual = rep.inner_constraint_residual
        inner_finite, inner_dt_mu = rep.inner_finite, rep.inner_dt_mu
        fluid_finite = rep.fluid_finite
        outer_accepted, outer_rolled_back = rep.outer_accepted, rep.outer_rolled_back
        committed_time = rep.t
        content_start, content_end, content_delta = rep.content_start, rep.total_content, rep.content_delta
        membrane_surface_integral = rep.membrane_surface_flux_integral
        membrane_grid_integral = rep.membrane_grid_source_integral
        interior_integral = rep.interior_source_integral
        dilatation_integral = rep.solid_dilatation_integral
        membrane_deposition_error = rep.membrane_deposition_error
        solver_conservation_error = rep.solver_conservation_error
        physical_conservation_error = rep.physical_conservation_error
    else:  # no-fluid core loop (still valid): run the inner solve directly
        wp.synchronize_device(cell.device)
        outer_t0 = time.perf_counter()
        inner_solve = make_inner_solve(
            cell,
            n_inner,
            reshape_every,
            max_inner_retries,
            inner_solver=inner_solver,
            implicit_cg_max_iterations=implicit_cg_max_iterations,
            implicit_line_search_steps=implicit_line_search_steps,
            implicit_coarse_iterations=implicit_coarse_iterations,
            implicit_coarse_modes=implicit_coarse_modes,
            rkc_stages=rkc_stages,
            capture_candidate=capture_convergence_history,
        )
        history_owner = inner_solve
        inner_d = inner_solve(dt_phys)
        inner_solve.commit_irreversible(inner_d.converged_d)  # type: ignore[attr-defined]
        wp.synchronize_device(cell.device)
        outer_wall = time.perf_counter() - outer_t0
        residual_candidate = float(inner_d.residual_d.numpy()[0])
        iters = int(inner_d.iters_d.numpy()[0])
        attempts = int(inner_d.attempts_d.numpy()[0])
        inner_converged = bool(inner_d.converged_d.numpy()[0])
        inner_max_displacement = float(inner_d.max_displacement_d.numpy()[0])
        inner_constraint_residual = float(inner_d.constraint_residual_d.numpy()[0])
        inner_finite = bool(inner_d.finite_d.numpy()[0])
        inner_dt_mu = float(inner_d.dt_mu_d.numpy()[0])
        outer_accepted = inner_converged
        outer_rolled_back = not inner_converged
        committed_time = dt_phys if inner_converged else 0.0

    postcheck_t0 = time.perf_counter()
    pos_all1 = cell.pos_d.numpy()
    pos1 = pos_all1[: cell.n_actin]
    rms1, rmin1, rmax1, com1, _ = _radii(pos1)
    pos_finite = bool(np.all(np.isfinite(pos_all1)))       # actin + nucleus + membrane all finite
    res_end, _ = _residual_host(cell, cell.pos_d, cell.f_d)
    # STABLE force check spans the WHOLE cell (actin + myosin + nucleus + membrane): a compartment force
    # blow-up must fail the verdict, not just an actin one (the resting compartments are force-free, so this
    # is a hardening gate that stays green today and catches divergence once GAP magnitudes rise, e.g. I7).
    force_finite = bool(np.all(np.isfinite(cell.f_d.numpy())))
    if cell.membrane is not None:
        cell.ledger["erm_bound_after_step"] = int(np.count_nonzero(cell.membrane.erm_bound_d.numpy()))
        cell.ledger["erm_bell_detach_events"] = int(cell.membrane.erm_detach_events_d.numpy()[0])
        cell.ledger["erm_bell_attach_events"] = int(cell.membrane.erm_attach_events_d.numpy()[0])
        cell.ledger["erm_rng_accepted_epoch"] = int(cell.membrane.erm_rng_epoch_d.numpy()[0])
    if cell.myosin is not None and cell.myosin.segment_runtime is not None:
        runtime = cell.myosin.segment_runtime
        cell.ledger["myosin_bound_after_step"] = int(np.count_nonzero(runtime.state["bound"].numpy()))
        cell.ledger["myosin_rng_accepted_epoch"] = int(runtime.rng_epoch.numpy()[0])
    unresolved_faces = 0
    if cell.membrane_pressure is not None:
        unresolved_faces = int(cell.membrane_pressure.unresolved_faces_d.numpy()[0])
    flux_unresolved_faces = 0
    if cell.membrane_bc is not None:
        flux_unresolved_faces = int(cell.membrane_bc.unresolved_faces_d.numpy()[0])
    membrane_query_failures = membrane_surface_ties = 0
    if cell.membrane_mask_provider is not None:
        membrane_query_failures = int(cell.membrane_mask_provider.query_failures_d.numpy()[0])
        membrane_surface_ties = int(cell.membrane_mask_provider.surface_ties_d.numpy()[0])
    nucleus_query_failures = nucleus_surface_ties = 0
    if cell.nucleus_mask_provider is not None:
        nucleus_query_failures = int(cell.nucleus_mask_provider.query_failures_d.numpy()[0])
        nucleus_surface_ties = int(cell.nucleus_mask_provider.surface_ties_d.numpy()[0])
    com_drift = float(np.linalg.norm(com1 - com0))
    memory_diagnostics_ok = True
    try:
        from warp._src.context import runtime as warp_runtime
        dev = wp.get_device(cell.device)
        cell.ledger["warp_mempool_used_current_bytes"] = int(
            warp_runtime.core.wp_cuda_device_get_mempool_used_mem_current(dev.ordinal))
        cell.ledger["warp_mempool_used_high_bytes"] = int(
            warp_runtime.core.wp_cuda_device_get_mempool_used_mem_high(dev.ordinal))
        cell.ledger["gpu_bytes_observed_after_step"] = int(dev.total_memory - dev.free_memory)
    except Exception as exc:  # noqa: BLE001 - record an explicit open ledger; never silently omit HARD evidence
        memory_diagnostics_ok = False
        cell.ledger["gpu_memory_accounting_status"] = "UNAVAILABLE"
        cell.ledger["gpu_memory_accounting_error"] = f"{type(exc).__name__}: {exc}"
    accounting = query_process_peak(cell.device)
    cell.ledger.update(accounting.ledger_fields())
    if accounting.exact and accounting.exact_peak_bytes is not None:
        cell.ledger["gpu_memory_accounting_status"] = "EXACT_WHOLE_DEVICE_PEAK"
        cell.ledger["gpu_memory_accounting_reason"] = (
            "NVML process-lifetime maxMemoryUsage includes this process's pooled and non-pooled CUDA allocations"
        )
    elif memory_diagnostics_ok:
        cell.ledger["gpu_memory_accounting_status"] = "PARTIAL_MEMPOOL_ONLY"
        cell.ledger["gpu_memory_accounting_reason"] = (
            "Warp mempool high-water excludes non-mempool CUDA allocations; "
            f"NVML exact probe status={accounting.probe_status}"
        )
    if history_owner is not None and getattr(history_owner, "history_iteration_d", None) is not None:
        # Post-loop readback (never inside the hot loop): the kernel wrote these device-side every
        # `reshape_every` iterations. Slot 0 stays zero-initialised, so a positive iteration is the
        # marker of a slot that was actually visited.
        it_h = history_owner.history_iteration_d.numpy()
        keep = it_h > 0
        cell.ledger["convergence_history"] = {
            "record_every": int(reshape_every),
            "iteration": [int(v) for v in it_h[keep]],
            "max_displacement_um": [float(v) for v in history_owner.history_displacement_d.numpy()[keep]],
            "constraint_residual_um": [float(v) for v in history_owner.history_constraint_d.numpy()[keep]],
            "max_projected_force_pn": [
                float(v) for v in history_owner.history_projected_force_d.numpy()[keep]
            ],
        }
    cell.ledger["inner_iteration_chunk_size"] = int(n_inner)
    cell.ledger["inner_max_retries"] = int(max_inner_retries)
    cell.ledger["inner_max_attempts"] = int(max_inner_retries + 1)
    cell.ledger["inner_iteration_budget"] = int(n_inner * (max_inner_retries + 1))
    postcheck_wall = time.perf_counter() - postcheck_t0
    wall = time.perf_counter() - total_t0

    stable = bool(
        pos_finite and force_finite and inner_finite and fluid_finite and inner_converged
        and outer_accepted and not outer_rolled_back
        and unresolved_faces == 0 and flux_unresolved_faces == 0
        and membrane_query_failures == 0 and nucleus_query_failures == 0 and np.isfinite(res_end)
        and (not profile_memcpy or (hot_dtoh_copies == 0 and positive_dtoh_copies > 0))
        and res_end <= max(2.0 * res_start, res_start + 1.0)   # committed forces not diverging
        and com_drift < 0.5                                    # no net internal-force drift
        and 0.5 * rms0 < rms1 < 2.0 * rms0                     # shell not puffing/collapsing (mean)
        and rmax1 < 3.0 * rms0                                 # no local node explosion
    )
    return StabilityReport(
        ran=True, n_actin=cell.n_actin, n_total=cell.n_total, dt_phys=dt_phys,
        n_biot_subcycles=int(n_sub), inner_iters=int(iters), inner_attempts=int(attempts),
        inner_max_attempts=int(max_inner_retries + 1),
        inner_iteration_budget=int(n_inner * (max_inner_retries + 1)), inner_converged=inner_converged,
        inner_tolerance_um=inner_tolerance_um, inner_max_displacement_um=inner_max_displacement,
        inner_constraint_residual_um=inner_constraint_residual,
        inner_finite=inner_finite, fluid_finite=fluid_finite, inner_dt_mu=inner_dt_mu,
        residual_start=res_start, residual_candidate=float(residual_candidate), residual_end=float(res_end),
        force_finite=force_finite,
        com_drift_um=com_drift, r_mean_start=rms0, r_mean_end=rms1, r_min_end=rmin1, r_max_end=rmax1,
        pos_finite=pos_finite,
        content_start=float(content_start), content_end=float(content_end), content_delta=float(content_delta),
        membrane_surface_flux_integral=float(membrane_surface_integral),
        membrane_grid_source_integral=float(membrane_grid_integral),
        interior_source_integral=float(interior_integral),
        solid_dilatation_integral=float(dilatation_integral), moving_face_delta=float(moving_delta),
        membrane_deposition_error=float(membrane_deposition_error),
        solver_conservation_error=float(solver_conservation_error),
        physical_conservation_error=float(physical_conservation_error),
        membrane_flux_unresolved_faces=flux_unresolved_faces,
        membrane_pressure_unresolved_faces=unresolved_faces,
        membrane_mask_query_failures=membrane_query_failures,
        membrane_mask_surface_ties=membrane_surface_ties,
        nucleus_mask_query_failures=nucleus_query_failures,
        nucleus_mask_surface_ties=nucleus_surface_ties,
        hot_loop_dtoh_copies=hot_dtoh_copies,
        hot_loop_dtoh_elapsed_ms=hot_dtoh_elapsed,
        profiler_positive_control_dtoh_copies=positive_dtoh_copies,
        outer_accepted=outer_accepted, outer_rolled_back=outer_rolled_back,
        committed_time_s=committed_time,
        build_wall_s=build_wall, precheck_wall_s=precheck_wall, outer_step_wall_s=outer_wall,
        postcheck_wall_s=postcheck_wall, wall_s=wall, stable=stable,
        ledger=cell.ledger)


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Composed Active Cell — one --from-resting outer step (P5).")
    p.add_argument("--from-resting", action="store_true", help="run the resting-stability outer step")
    p.add_argument("--n-filaments", type=int, default=70686, help="cortical F-actin count (full=70686)")
    p.add_argument("--n-inner", type=int, default=60, help="inner mechanical iterations for the outer step")
    p.add_argument(
        "--max-inner-retries", type=int, default=0,
        help="additional same-time mechanics chunks after the initial n-inner budget (default: 0)",
    )
    p.add_argument("--dt-phys", type=float, default=0.05, help="outer physical step [s]")
    p.add_argument(
        "--inner-solver",
        choices=("explicit", "fire", "analytic_implicit", "block_descent", "anderson", "rkc1",
                 "contact_schwarz", "fiber_contact_schwarz", "tournament", "contact_tournament",
                 "cluster_tournament", "rigid_contact_schwarz", "rigid_cluster_tournament",
                 "erm_schwarz", "erm_tournament", "augmented_block", "augmented_tournament",
                 "erm_jacobi_block", "erm_jacobi_tournament",
                 "erm_jacobi_pure", "erm_jacobi_pure_tournament",
                 "erm_gauss_seidel", "erm_gauss_seidel_tournament"),
        default="explicit",
        help="inner mechanics update (fire = tangent-free FIRE static minimizer; every accelerated "
             "candidate remains exact-residual gated)",
    )
    p.add_argument(
        "--implicit-cg-max-iterations", type=int, default=32,
        help="fixed PCG launch budget when --inner-solver=analytic_implicit",
    )
    p.add_argument(
        "--implicit-line-search-steps", type=int, default=4,
        help="fixed full, half, quarter, ... exact-residual trials for accelerated mechanics",
    )
    p.add_argument(
        "--implicit-coarse-iterations", type=int, default=8,
        help="fixed 1/2-Jacobi sweeps on the per-fiber translation contact graph",
    )
    p.add_argument(
        "--implicit-coarse-modes", type=int, default=0,
        help="global rigid-body+constant-strain (l<=2) coarse deflation modes in the analytic "
             "preconditioner (0=off; 12=full). Preconditioner-only: accelerates the coupled global "
             "smooth mode, never moves the fixed point or relaxes a residual gate",
    )
    p.add_argument(
        "--rkc-stages", type=int, default=4,
        help="first-order Chebyshev recurrence stages (each uses the validated explicit dt_mu)",
    )
    p.add_argument("--no-myosin", action="store_true", help="omit MyosinForce composition")
    p.add_argument("--no-steric", action="store_true", help="omit StericForce composition")
    p.add_argument("--preload-erm-balance", action="store_true",
                   help="pre-stretch ERM tethers so resting membrane turgor is held at t=0 "
                        "(physiological baseline; derived from the live pressure and k_erm)")
    p.add_argument("--overlap-free-cortex", action=argparse.BooleanOptionalAction, default=True,
                   help="build the cortex overlap-free (radial ~0.2µm thickness + WCA relaxation) so it starts "
                        "with zero steric force — removes the ~64k build-interpenetration artifact (2634 pN). "
                        "Default ON (physiological baseline); pass --no-overlap-free-cortex for the legacy "
                        "zero-thickness γ/Gate-1-parity build.")
    p.add_argument("--cortex-overlap-mode", choices=["transverse", "radial_span"], default="transverse",
                   dest="cortex_overlap_mode",
                   help="HOW overlap_free_cortex resolves build interpenetrations: 'transverse' (default, "
                        "byte-identical) shoves crossing nodes IN-plane; 'radial_span' separates them OUT-of-plane "
                        "with a smooth radial bump so fine-mesh fibers stay unkinked.")
    p.add_argument("--cortex-overlap-span", type=int, default=2, dest="cortex_overlap_span",
                   help="radial_span bump half-window in arc-neighbours (default 2 ⇒ 5-node bump).")
    p.add_argument("--membrane-subdivisions", type=int, default=3,
                   help="plasma-membrane icosphere level (3=642 resting, 6=40,962 dynamic baseline, 7=163,842 "
                        "validation; PI 2026-07-22)")
    p.add_argument("--nucleus-subdivisions", type=int, default=3,
                   help="nuclear-envelope icosphere level (same 3/6/7 scaling)")
    p.add_argument("--no-nucleus", action="store_true", help="omit the deformable-mesh NucleusCompartment")
    p.add_argument("--no-membrane", action="store_true", help="omit the plasma-membrane Helfrich sheet")
    p.add_argument("--no-pressure", action="store_true", help="omit the Biot substrate + PressureCoupling")
    p.add_argument(
        "--turgor-pi0", type=float, default=None,
        help="Π₀ axis point [Pa]. Omit for the gated default (byte-identical to before this flag "
             "existed). 0 keeps the Biot substrate and PressureCoupling composed while removing the "
             "turgor LOAD, which --no-pressure cannot do because it removes both at once.")
    p.add_argument(
        "--erm-density", type=float, default=None,
        help="explicit ERM linker density [um^-2]; production use requires --erm-density-source",
    )
    p.add_argument(
        "--erm-density-source", type=str, default="",
        help="literature/provenance label for a non-null --erm-density",
    )
    p.add_argument(
        "--erm-density-mcf7-production", action="store_true",
        help="assert that the supplied density is an MCF7 production contract, not a cross-cell proxy",
    )
    p.add_argument("--erm-k-on", type=float, default=None, help="source-gated ERM attach rate [s^-1]")
    p.add_argument("--erm-k-off0", type=float, default=None, help="source-gated ERM zero-force off-rate [s^-1]")
    p.add_argument("--erm-bell-force", type=float, default=None, help="source-gated ERM Bell force F0 [pN]")
    p.add_argument(
        "--erm-capture-radius", type=float, default=None,
        help="source-gated ERM rebind capture radius [um]",
    )
    p.add_argument(
        "--erm-kinetics-source", type=str, default="",
        help="literature/provenance label covering the complete ERM Bell kinetics contract",
    )
    p.add_argument(
        "--nmii-backbone-lp", type=float, default=None,
        help="source-gated NMII backbone persistence length [um]; no production default",
    )
    p.add_argument(
        "--nmii-backbone-lp-source", type=str, default="",
        help="literature/provenance label for --nmii-backbone-lp",
    )
    p.add_argument(
        "--resting-bound-myosin-fraction", type=float, default=None,
        help="RESTING cortical-tension source: fraction of NMII heads bound at rest (PI-GAP; no default). "
             "Requires --resting-bound-myosin-force and --resting-bound-myosin-source; engages radial ERM.",
    )
    p.add_argument(
        "--resting-bound-myosin-force", type=float, default=None,
        help="resting isometric per-head tangential tension [pN] carried by each bound head (PI-GAP; no default)",
    )
    p.add_argument(
        "--resting-bound-myosin-source", type=str, default="",
        help="literature/provenance label covering the resting bound-myosin fraction + per-head force",
    )
    p.add_argument(
        "--resting-bound-myosin-capture", type=float, default=None,
        help="resting bind reach [um] (default: the motor's params.capture_radius)",
    )
    p.add_argument("--device", type=str, default=None, help="Warp CUDA device alias (default: current CUDA)")
    p.add_argument("--profile-memcpy", action="store_true",
                   help="profile hot-loop copies and verify a post-loop DtoH positive control")
    p.add_argument("--json", type=str, default="", help="write the report JSON to this path")
    return p


def main() -> None:
    args = _build_argparser().parse_args()
    cfg = CellConfig(
        n_filaments=args.n_filaments, with_myosin=not args.no_myosin, with_steric=not args.no_steric,
        overlap_free_cortex=args.overlap_free_cortex,
        cortex_overlap_mode=args.cortex_overlap_mode,
        cortex_overlap_span=args.cortex_overlap_span,
        membrane_subdivisions=args.membrane_subdivisions, nucleus_subdivisions=args.nucleus_subdivisions,
        with_nucleus=not args.no_nucleus, with_membrane=not args.no_membrane,
        with_pressure=not args.no_pressure, turgor_pi0_pa=args.turgor_pi0,
        erm_density_per_um2=args.erm_density,
        erm_density_source=args.erm_density_source,
        erm_density_mcf7_production=args.erm_density_mcf7_production,
        erm_k_on_s=args.erm_k_on, erm_k_off0_s=args.erm_k_off0,
        erm_bell_force_pn=args.erm_bell_force, erm_capture_radius_um=args.erm_capture_radius,
        erm_kinetics_source=args.erm_kinetics_source,
        nmii_backbone_lp_um=args.nmii_backbone_lp,
        nmii_backbone_lp_source=args.nmii_backbone_lp_source,
        resting_bound_myosin_fraction=args.resting_bound_myosin_fraction,
        resting_bound_myosin_force_pn=args.resting_bound_myosin_force,
        resting_bound_myosin_source=args.resting_bound_myosin_source,
        resting_bound_myosin_capture_um=args.resting_bound_myosin_capture,
        device=args.device)
    print(f"[ac.cell] building composed cell: n_filaments={cfg.n_filaments} myosin={cfg.with_myosin} "
          f"steric={cfg.with_steric} nucleus={cfg.with_nucleus} membrane={cfg.with_membrane} "
          f"pressure={cfg.with_pressure} erm_density={cfg.erm_density_per_um2} device={cfg.device}", flush=True)
    rep = run_from_resting(
        cfg, n_inner=args.n_inner, dt_phys=args.dt_phys, max_inner_retries=args.max_inner_retries,
        profile_memcpy=args.profile_memcpy, inner_solver=args.inner_solver,
        implicit_cg_max_iterations=args.implicit_cg_max_iterations,
        implicit_line_search_steps=args.implicit_line_search_steps,
        implicit_coarse_iterations=args.implicit_coarse_iterations,
        implicit_coarse_modes=args.implicit_coarse_modes,
        rkc_stages=args.rkc_stages, preload_erm_balance=args.preload_erm_balance)
    print("\n=== population ledger ===", flush=True)
    for k, v in rep.ledger.items():
        if k == "per_region":
            continue
        print(f"  {k:28s} {v:,}" if isinstance(v, int) else f"  {k:28s} {v}", flush=True)
    print("\n=== --from-resting stability ===", flush=True)
    for k, v in rep.to_dict().items():
        if k == "ledger":
            continue
        print(f"  {k:22s} {v}", flush=True)
    print(f"\nSTABLE: {rep.stable}", flush=True)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"schema": "ffn-ac-native-resting-v2", "config": asdict(cfg),
                       "report": rep.to_dict(), "ledger": rep.ledger}, fh, indent=2, default=float)
        print(f"wrote {args.json}", flush=True)


if __name__ == "__main__":
    main()

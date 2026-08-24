#!/usr/bin/env python
r"""GATE-B NATIVE driver — the interior column membrane → cortex → cytosol → nucleus as ONE transaction.

The THIRD dynamic-runtime slice (after the NMII cortex-motor and the ERM membrane--cortex slices).  The first
two grew ACTIVE tension along the SURFACE; this one drives :class:`aleph.engine.transaction.CellTransaction`
INWARD, tying the four load-bearing interior compartments together under ONE device-resident acceptance
predicate on the real composed cell (CUDA / shared Slurm workstation)::

    resting turgor pressurises the cytosol Biot field  →  the field's pressure traction rides outward onto the
    membrane + cortex and inward onto the nuclear envelope (Peskin −α·V·grad(p), Newton's 3rd law)  →  water
    permeates the SEMIPERMEABLE membrane (Kedem–Katchalsky J = L_p(σΔΠ − ΔP)) while the IMPERMEABLE envelope
    stays no-flux (mask-structural)  →  the nucleus lamina/volume/LINC balance the pressure, and any face whose
    areal strain exceeds ε_rupt tears (accepted-gated)  →  one coupled interior-column steady state.

This EXTENDS the proven surface slice inward: it introduces NO new physics — every kernel it launches is one
already bound by ``ac/engine/fluid_core.py`` + ``ac/engine/cytosol_connected.py`` (Biot p/mass, the Peskin
pressure traction, the KK membrane flux, the conservative moving-domain remap) or by ``ac/nucleus`` (per-face
lamina tension, nucleoplasm volume, LINC, the device-gated envelope-rupture commit).  The slice is the
COMPOSITION: it re-expresses the incumbent single-owner coupled fluid step as a :class:`CellTransaction` whose
four owners + three fluid connectors snapshot / roll back / commit-on-accept under one predicate.

────────────────────────────────────────────────────────────────────────────────────────────────────────
RUN ON THE SHARED WORKSTATION THROUGH SLURM (this needs CUDA; it will NOT run on the dev Mac):

    MEM_GB=48 CPUS_PER_TASK=8 gpu-submit <approved-card> <approved-time> \
      ~/miniforge3/envs/ffn_sim/bin/python \
      aleph/scripts/ac_gate_b_interior_column_native.py --steps 40 --dt 0.01 --outer 40 \
      --membrane-subdiv 8 --nucleus-subdiv 3

The PI must name the card and duration in-session before submission; monitor the Slurm log.
────────────────────────────────────────────────────────────────────────────────────────────────────────

NO DOUBLE COUNT (the ERM-slice precedent).  The built cell already owns the whole coupled interior mechanics:
the driver's ``_accumulate_all`` (through ``make_inner_solve``) assembles the membrane Helfrich/area, the cortex
network, the nucleus lamina/volume/LINC, AND the real Biot pressure coupling (``PressureCoupling``,
``MembranePressureTraction``).  So this driver lets the driver's inner solve OWN the coupled mechanics, and the
SLICE owns the accepted-step TRANSACTION on top of it: it snapshots + reject-restores the Biot field
(``ac/fluid/scheduler`` restore kernels) and the nucleus surface + ``ruptured`` flag, and it commits the
accepted-gated envelope rupture (``conditional_rupture_update_kernel``) + advances the device clock — all four
compartments under ONE predicate.  The slice's coupling connectors are bound as transaction participants and
hold the SAME real primitives (``cell.pressure`` / ``cell.membrane_bc`` / ``cell.domain``); their standalone
adjoint-force scatter (``InteriorColumnSlice.accumulate``) is unit-real (CPU structural gate) but is NOT
re-invoked inside this loop — the driver's ``_accumulate_all`` already scatters that traction, and re-adding it
would double count.  This is exactly how the ERM native driver keeps ONE ERM force owner.

ACCEPTANCE IS PHYSICAL.  The inner solve's device-resident ``converged_d`` is copied into the ONE transaction
predicate during ``solve``.  A rejected candidate rolls back the field, geometry, rupture state and clock and
the driver stops with C-2 BLOCKED.  This is deliberately stricter than the former mechanism demo, which
force-accepted a non-converged candidate and could print ``DONE`` for repeated identical ``maxPF`` values.

PROVISIONAL PARAMS (PI/KB GAP).  The nucleus lamina/volume moduli + the rupture strain ε_rupt are the open
I0-B2 GAP-PI (``ac/nucleus/params_i0b2.yaml``: "do NOT keep the ff proxy 0.50"); the osmotic driver Π₀ / L_p
are the membrane-boundary GAP.  This driver reads them from the built cell (the incumbent's provisional TEST
values) and prints a loud banner — a MECHANISM demonstration, NOT quantitative production.  k_erm-style sourced
constants are read from the cell, never invented here.

MEMBRANE FIDELITY.  The grid-converged C-2 default is subdivision 8, not 6.  The native grid-convergence audit
(``GATE_A_RESTING_CONVERGENCE_2026-07-23.md`` section 4) measured the subdivision-6 membrane residual as
0.7766 pN and the subdivision-8 residual as 0.0768 pN under the same physical setpoint.  Subdivision 6 is
therefore a diagnosed discretization artifact and cannot support this gate.  The CLI remains exposed only for
the already-established grid-convergence control; a production C-2 record uses 8.

Runtime: NVIDIA Warp on CUDA only (I0-A).  Authored on the dev Mac (no CUDA); native runs use Slurm on the
shared workstation.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import warp as wp

from aleph.components.fluid.boundary import NucleusNoFluxBC
from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import (
    _preload_cortex_pretension,
    _preload_erm_resting_balance,
    make_inner_solve,
)
from aleph.engine.cortex_state import assert_component_state_disjoint
from aleph.engine.afm_cortical_tension import (
    AFM_INDENTER_COMPONENT,
    AFM_MEMBRANE_CONTACT,
    afm_experiment_architecture,
)
from aleph.engine.afm_indenter import AFMIndenterForceHook, AFMIndenterRuntime, AFMIndenterRuntimeSpec
from aleph.engine.afm_loading import derive_force_free_initial_centre, load_ready_afm_path
from aleph.engine.cytosol_connected import (
    MembranePressureFluxAdjointBoundary,
    MovingSemipermeableFluidBoundaryFacade,
)
from aleph.engine.erm_cortex_slice import build_native_surface_owner
from aleph.engine.fluid_core import (
    BiotSubstrateFluidSolver,
    MovingImpermeableFluidBoundaryFacade,
    NucleusPressureAdjointBoundary,
    lumped_surface_control_volumes,
)
from aleph.engine.interior_column_slice import (
    RecordingLedgerStub,
    assert_disjoint_blocks,
    build_core_surface_owner,
    build_cortex_cytosol_transfer,
    build_fluid_field_owner,
    build_interior_column_slice,
    component_block_views,
)
from aleph.engine.medium_exterior import (
    EXTERIOR_MEDIUM_COMPONENT,
    MEMBRANE_MEDIUM_CONNECTOR,
    ExteriorStokesMediumSettings,
    MembraneMediumTraction,
    build_exterior_stokes_medium,
    derive_blob_epsilon,
)
from aleph.engine.medium_stokes_analytic import resistance_spectral_upper_bound
from aleph.laws.surface_manifold import n_vert_for_subdivisions

_FLUID_CODE = 1  # FieldGrid mask: OUTSIDE=0 / FLUID=1 / NUCLEUS=2 (ac.fluid.field_grid)


class _NoOpParticipant:
    """A no-op :class:`TransactionParticipant` for a moving-boundary cache the driver already manages.

    The nucleus/membrane fluid-boundary delegates forward snapshot/rollback/commit to an injected participant.
    In this no-double-count driver the driver's inner solve owns the moving-domain remap + pressure caches, so
    the injected participant is a no-op (the field + nucleus reversible state is snapshot/rolled-back by the
    slice's own FluidFieldTransaction + CoreSurfaceTransaction).
    """

    def snapshot_candidate(self) -> None: ...
    def rollback(self, accepted: wp.array) -> None: ...
    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None: ...


@wp.kernel
def _ramp_turgor_fluid_kernel(
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    d_pa: wp.float64,
) -> None:
    """Raise the interior osmotic-turgor setpoint by ``d_pa`` [Pa == pN/µm²] on every FLUID cell (device-only).

    This is the LOAD RAMP.  ``p`` is the same conservative Biot field the membrane traction + the pressure
    coupling read, so bumping the FLUID interior raises the outward pressure the inner solve must balance — the
    membrane genuinely bulges, the cortex feels the porous traction, and the nuclear envelope feels the inward
    pressure.  Nothing is injected into any tether; every coupling force is the real reaction of the moved field.
    OUTSIDE/NUCLEUS cells are untouched (no spurious body gradient across the impermeable envelope).
    """
    i, j, k = wp.tid()
    if mask[i, j, k] == _FLUID_CODE:
        p[i, j, k] = p[i, j, k] + d_pa


def _banner(cell) -> None:
    """Print the loud PROVISIONAL-PARAMS banner listing the interior-column GAP magnitudes + their status."""
    nuc = cell.nucleus
    rows = [
        ("nucleus k_soft [pN/µm]", float(nuc.k_soft), "I0-B2 GAP-PI (chromatin+lamin-B areal modulus)"),
        ("nucleus k_ac   [pN/µm]", float(nuc.k_ac), "I0-B2 GAP-PI (lamin-A/C strain-stiffened tangent)"),
        ("nucleus knee   [-]     ", float(nuc.knee), "I0-B2 provisional (strain crossover)"),
        ("nucleus eps_rupt [-]   ", float(nuc.eps_rupt), "I0-B2 GAP-PI — do NOT keep the ff proxy 0.50"),
        ("nucleus k_vol [pN/µm²] ", float(nuc.k_vol), "I0-B2 (nucleoplasm incompressible penalty)"),
        ("nucleus k_linc[pN/µm]  ", float(nuc.k_linc), "I0-B2 HALT→PI (nesprin nonlinear, no Hookean k)"),
        ("cytosol Π₀    [Pa]     ", 40.0, "GAP (HeLa proxy; osmotic vertical deferred POST-GATE-A)"),
    ]
    print("=" * 104)
    print("  GATE-B PROVISIONAL PARAMETERS  —  interior column membrane→cytosol→nucleus (fluid + Core Body)")
    print("  MECHANISM DEMONSTRATION, NOT quantitative production — the nucleus + osmotic magnitudes are PI GAPs.")
    print("-" * 104)
    for label, value, status in rows:
        print(f"    {label} = {value:>10.5g}   {status}")
    print("=" * 104, flush=True)


def _local_faces_verts(pos_np: np.ndarray, faces_d: wp.array, node_off: int, n_verts: int):
    """Return (verts_local, faces_local) for a global-indexed triangle surface — for the node-volume weights."""
    faces_g = faces_d.numpy().astype(np.int64)
    faces_local = faces_g - int(node_off)
    verts_local = pos_np[int(node_off):int(node_off) + int(n_verts)]
    return verts_local, faces_local


def _column_stats(cell, slc) -> dict:
    """Emergent interior-column telemetry (out-of-hot-loop host readback; I0-A honoured inside the loop).

    Reports the Biot field pressure over the FLUID interior [Pa], the net membrane water permeation
    (``cell.membrane_bc.integrated_flux_d`` [µm³/s·area]), the would-be no-flux teeth the impermeable envelope
    holds at 0, and the ruptured-face count — all emergent from the coupled candidate, no imposed prestress.
    """
    p = cell.grid.p.numpy()
    mask = cell.grid.mask.numpy()
    fluid = mask == _FLUID_CODE
    flux = float(cell.membrane_bc.integrated_flux_d.numpy()[0]) if cell.membrane_bc is not None else 0.0
    teeth = slc.nucleus_boundary.delegate.no_flux_teeth()
    ruptured = slc.ruptured_face_count()
    return {
        "p_mean_pa": float(p[fluid].mean()) if fluid.any() else 0.0,
        "p_max_pa": float(p[fluid].max()) if fluid.any() else 0.0,
        "membrane_flux": flux,
        "no_flux_teeth": teeth,
        "ruptured_faces": ruptured,
    }


def _ownership_gate(cell, blocks: dict, owners: list, out_dir: Path, device: str, build_s: float,
                    build_commit: str | None) -> None:
    """Certify the Card-5 array split at full native, with the positive control, and write its artifact.

    The gate is the pair, never the assertion alone: a check that cannot fail certifies nothing, and this
    tree has now found five guards that returned green for a reason unrelated to what they claimed.  So the
    control rebuilds the owners the way this driver built them BEFORE Card-5 — all three handed the whole
    ``cell.pos_d``/``cell.f_d`` — and requires :func:`assert_component_state_disjoint` to refuse it.

    What the artifact may be quoted for: that membrane / cortex / nucleus address DISJOINT blocks at the
    full native population, and that the same check rejects the arrangement they replaced.  What it may NOT
    be quoted for: private allocation (``arrays_private`` is recorded FALSE — these are slice views of the
    one array the incumbent's inner solve relaxes), any rung, any physics magnitude, or ``STATE.md`` (c) 4
    being closed.  Only the PI closes that.
    """
    from aleph.engine.cortex_state import _byte_span

    control_raised = None
    try:
        control = [
            build_native_surface_owner(name="membrane", position_d=cell.pos_d, force_d=cell.f_d, device=device),
            build_native_surface_owner(name="cortex", position_d=cell.pos_d, force_d=cell.f_d, device=device),
        ]
        assert_component_state_disjoint(control)
    except ValueError as exc:
        control_raised = str(exc)
    if control_raised is None:
        raise RuntimeError(
            "POSITIVE CONTROL DID NOT FAIL: the pre-Card-5 arrangement (every surface owner handed the whole "
            "global array) was accepted by assert_component_state_disjoint. The gate cannot certify the split "
            "when the check it uses accepts what the split replaced."
        )

    spans = {
        owner.name: {
            "position": list(_byte_span(owner.position_d) or ()),
            "force": list(_byte_span(owner.force_d) or ()),
            "n_nodes": int(owner.position_d.shape[0]),
        }
        for owner in owners
    }
    # The run host is not a git checkout, so the build cannot be READ here — it is DECLARED by the caller
    # and its evidence is `run_provenance.remote_mismatches`, which compares the driver's first-party
    # import closure file-by-file against this host before launch. That is the same `source: declared`
    # convention the cortex-ownership row carries, and a record with no commit at all is worse: it cannot
    # be corrected afterwards, because the run is over and the number is in a file.
    if not build_commit:
        raise RuntimeError(
            "--build-commit is REQUIRED for the ownership gate. The run host is not a git checkout, so "
            "nothing here can read the build; an unstamped record cannot be traced to the code that "
            "produced it. Verify with run_provenance.remote_mismatches first, then declare the commit."
        )
    record = {
        "record": "run-record@2",
        "kind": "gate",
        "gate": "card5_interior_column_array_ownership",
        "device": str(device),
        "build": {
            "commit": str(build_commit),
            "source": "declared",
            "why": "the run host is an rsync'd tree, not a git checkout — the caller asserts the commit "
                   "and run_provenance.remote_mismatches is what checks the assertion",
        },
        "population": {
            "n_actin": int(cell.n_actin), "n_total": int(cell.n_total),
            "membrane_verts": int(cell.membrane.n_verts), "nucleus_verts": int(cell.nucleus.n_verts),
            "nucleus_faces": int(cell.nucleus.n_faces),
        },
        "blocks": {name: {"node_off": int(off), "n_nodes": int(n)} for name, (off, n) in sorted(blocks.items())},
        "byte_spans": spans,
        "arrays_private": False,
        "ownership_class": "disjoint_slice_views",
        "what_this_is_not": (
            "Slice views of the incumbent's global arrays, NOT private allocations. The incumbent's inner "
            "solve relaxes one global pos_d; a private allocation without splitting that solve would be a "
            "passenger (the T2 negative control). This certifies exclusive ADDRESSING only."
        ),
        "positive_control": {
            "arrangement": "membrane and cortex both handed the whole cell.pos_d/cell.f_d (pre-Card-5)",
            "refused": True,
            "message": control_raised,
        },
        "verdict": "PASS",
        "timing": {"build_s": float(build_s), "comparable": True},
        "may_not_be_quoted_for": [
            "any CONNECTED rung", "private device allocations", "STATE.md (c) 4 being closed",
            "any physics magnitude",
        ],
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "native_record.json").write_text(json.dumps(record, indent=2))
    _plot_blocks(blocks, int(cell.n_total), out_dir / "array_ownership_blocks.png")
    print(f"[gate-b-column] OWNERSHIP GATE PASS -> {out_dir}/native_record.json (+1 fig). Positive control "
          f"refused the pre-Card-5 arrangement: {control_raised[:90]}...", flush=True)


def _c2_diagnostic_step_path(path: Path, *, step_index: int, total_steps: int) -> Path:
    """Keep every physical-step record instead of overwriting one diagnostic in a multi-step gate."""
    if total_steps <= 1:
        return path
    return path.with_name(f"{path.stem}.step-{step_index:03d}{path.suffix}")


def _write_c2_diagnostic(path: Path, *, args: argparse.Namespace, cell, inner_solve, report,
                         max_projected_force_pn: float, accepted: bool, physical_step_index: int) -> None:
    """Persist the already device-captured C-2 convergence trace after transaction resolution."""
    iterations = int(report.iters_d.numpy()[0])
    iteration_history = inner_solve.history_iteration_d.numpy()
    valid = (iteration_history > 0) & (iteration_history <= iterations)
    dt_mu = float(report.dt_mu_d.numpy()[0])
    tolerance_um = float(np.sqrt(np.finfo(np.float64).eps) * cell.convergence_length_um)
    trial_counts = None
    if inner_solve.implicit_accept_count_d is not None:
        families = tuple(str(name) for name in inner_solve.accelerator_candidate_families)
        scales = tuple(float(scale) for scale in inner_solve.implicit_line_search_scales)
        flat_counts = inner_solve.implicit_scale_accept_counts_d.numpy().astype(np.int64, copy=False)
        by_family = {
            family: {
                f"{scale:.17g}": int(flat_counts[family_index * len(scales) + scale_index])
                for scale_index, scale in enumerate(scales)
            }
            for family_index, family in enumerate(families)
        }
        trial_counts = {
            "accelerated": int(inner_solve.implicit_accept_count_d.numpy()[0]),
            "explicit": int(inner_solve.explicit_accept_count_d.numpy()[0]),
            "stationary": int(inner_solve.stationary_accept_count_d.numpy()[0]),
            "accelerated_by_family_and_scale": by_family,
        }
    payload = {
        "status": "ACCEPTED" if accepted else "C2_REJECTED_ROLLED_BACK",
        "physical_step_index": int(physical_step_index),
        "device": str(cell.device),
        "full_native": {
            "n_filaments": int(cell.cfg.n_filaments),
            "n_total_nodes": int(cell.n_total),
            "membrane_subdivisions": int(args.membrane_subdiv),
            "nucleus_subdivisions": int(args.nucleus_subdiv),
        },
        "solver": {
            "name": str(args.inner_solver),
            "dt_phys_s": float(args.dt),
            "iterations_per_chunk": int(args.outer),
            "max_inner_retries": int(args.max_inner_retries),
            "iteration_budget": int(inner_solve.iteration_budget),
            "line_search_steps": int(args.line_search_steps),
            "line_search_objective": str(args.line_search_objective),
            "coarse_modes": int(args.coarse_modes),
            "coarse_iterations": int(args.coarse_iterations),
            "reported_iterations": iterations,
            "attempts": int(report.attempts_d.numpy()[0]),
            "candidate_residual_pN": float(report.residual_d.numpy()[0]),
            "max_projected_force_pN": float(max_projected_force_pn),
            "max_displacement_um": float(report.max_displacement_d.numpy()[0]),
            "constraint_residual_um": float(report.constraint_residual_d.numpy()[0]),
            "dt_mu_um_per_pN": dt_mu,
            "displacement_tolerance_um": tolerance_um,
            "projected_force_tolerance_pN": tolerance_um / dt_mu,
            "accepted_trial_counts": trial_counts,
            "convergence_history": {
                "iteration": iteration_history[valid].tolist(),
                "max_displacement_um": inner_solve.history_displacement_d.numpy()[valid].tolist(),
                "constraint_residual_um": inner_solve.history_constraint_d.numpy()[valid].tolist(),
                "projected_force_pN": inner_solve.history_projected_force_d.numpy()[valid].tolist(),
            },
        },
        "acceptance": {
            "accepted": bool(accepted),
            "candidate_geometry_captured_before_rollback": inner_solve.candidate_pos_d is not None,
            "rejected_candidate_never_advanced_biological_time": not accepted,
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"[gate-b-column] C-2 diagnostic: {path}", flush=True)


def _plot_blocks(blocks: dict, n_total: int, path: Path) -> None:
    """Draw the global node array with each component's block, at true scale (no axis truncation)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 3.2))
    colours = {"cortex": "#4C72B0", "nucleus": "#C44E52", "membrane": "#55A868"}
    ax.broken_barh([(0, n_total)], (0.55, 0.28), facecolors="#DDDDDD", edgecolor="none")
    ax.text(n_total / 2, 0.69, f"cell.pos_d / cell.f_d — ONE allocation, {n_total:,} nodes",
            ha="center", va="center", fontsize=9)
    # Blocks are drawn at TRUE scale, so a small compartment stays a sliver — that is the fact, not a
    # rendering defect. Labels are staggered by block order because at true scale they would otherwise
    # overlap, which is the only thing about them that may be adjusted.
    for row, (name, (off, n)) in enumerate(sorted(blocks.items(), key=lambda kv: kv[1][0])):
        ax.broken_barh([(off, n)], (0.15, 0.28), facecolors=colours.get(name, "#8172B2"), edgecolor="none")
        ax.annotate(f"{name}  [{off:,}, {off + n:,})  {n:,} nodes",
                    xy=(off + n / 2, 0.15), xytext=(off + n / 2, -0.12 - 0.22 * row), ha="center", fontsize=8,
                    arrowprops=dict(arrowstyle="-", lw=0.6))
    ax.set_xlim(0, n_total)
    ax.set_ylim(-0.16 - 0.22 * max(len(blocks) - 1, 0), 0.95)
    ax.set_yticks([])
    ax.set_xlabel("node index in the global array [-]")
    ax.set_title("Card-5: membrane / cortex / nucleus address DISJOINT blocks — slice views, not private "
                 "allocations", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="GATE-B: the interior column membrane→cytosol→nucleus (native).")
    ap.add_argument("--steps", type=int, default=40, help="number of ACCEPTED physical steps to run")
    ap.add_argument("--dt", type=float, default=0.01, help="physical timestep dt_phys [s]")
    ap.add_argument("--outer", type=int, default=40, help="inner mechanical relax iterations per step (n_inner)")
    ap.add_argument(
        "--inner-solver", type=str, default="explicit", dest="inner_solver",
        choices=(
            "explicit", "fire", "analytic_implicit", "block_descent", "anderson", "rkc1",
            "contact_schwarz", "fiber_contact_schwarz", "rigid_contact_schwarz", "tournament",
            "contact_tournament", "cluster_tournament", "rigid_cluster_tournament", "erm_schwarz",
            "erm_tournament", "augmented_block", "augmented_tournament", "erm_jacobi_block",
            "erm_jacobi_tournament", "erm_jacobi_pure", "erm_jacobi_pure_tournament",
            "erm_gauss_seidel", "erm_gauss_seidel_tournament",
        ),
        help="existing make_inner_solve numerical operator; physics and convergence predicate stay unchanged",
    )
    ap.add_argument("--line-search-steps", type=int, default=4, dest="line_search_steps",
                    help="fixed geometric line-search depth for accelerated solvers")
    ap.add_argument("--line-search-objective", choices=("max_force", "l2_squared"),
                    default="max_force", dest="line_search_objective",
                    help="candidate-selection metric only; C-2 acceptance always retains max projected force")
    ap.add_argument("--max-inner-retries", type=int, default=0, dest="max_inner_retries",
                    help="additional fixed numerical chunks in the same candidate; never advances time")
    ap.add_argument("--coarse-modes", type=int, default=0, dest="coarse_modes",
                    help="existing geometry-derived global rigid/strain preconditioner modes (0-12)")
    ap.add_argument("--coarse-iterations", type=int, default=8, dest="coarse_iterations",
                    help="existing fiber-translation coarse preconditioner iterations")
    ap.add_argument("--filaments", type=int, default=70686, help="cortex F-actin count (native = 70686)")
    ap.add_argument("--membrane-subdiv", type=int, default=8, dest="membrane_subdiv",
                    help="membrane icosphere level (production C-2 default 8=655,362 verts; 6 is a diagnosed "
                         "grid-truncation control and cannot support acceptance)")
    ap.add_argument("--nucleus-subdiv", type=int, default=3, dest="nucleus_subdiv",
                    help="nucleus oblate-mesh subdivision level")
    ap.add_argument("--turgor-ramp", type=float, default=0.0, dest="turgor_ramp",
                    help="LOAD RAMP: raise the interior turgor Π by this many Pa PER ACCEPTED STEP (0.0 = "
                         "constant-load baseline). >0 pressurises the cytosol so the pressure traction rides "
                         "onto membrane/cortex and inward onto the envelope, and high-strain faces may rupture.")
    ap.add_argument("--preload", action=argparse.BooleanOptionalAction, default=True, dest="preload",
                    help="establish the FORCE-BALANCED resting baseline before the loop (cortex pretension + ERM "
                         "resting balance), the SAME preloads the production driver uses. DEFAULT ON "
                         "(PI physiological-baseline rule).")
    ap.add_argument("--cortex-prestrain", type=float, default=0.0, dest="cortex_prestrain",
                    help="cortex crosslink pre-strain fed to _preload_cortex_pretension [dimensionless].")
    ap.add_argument("--seed", type=int, default=0, help="base RNG seed")
    ap.add_argument("--ownership-gate", action="store_true", dest="ownership_gate",
                    help="CARD-5: build the full native cell, certify that membrane/cortex/nucleus address "
                         "DISJOINT blocks (with the positive control that must fail), write the record + "
                         "figure, and exit WITHOUT the coupled loop. The coupled step does not converge "
                         "today (C-2), and this gate is about addressing, not about the step.")
    ap.add_argument("--ownership-out", type=Path,
                    default=Path("aleph/outputs/ac/interior_column_ownership"), dest="ownership_out",
                    help="directory for the ownership gate's record + figure")
    ap.add_argument("--build-commit", type=str, default=None, dest="build_commit",
                    help="REQUIRED with --ownership-gate: the commit whose code is on this host, verified "
                         "with run_provenance.remote_mismatches before launch. Declared, not read — the "
                         "run host is an rsync'd tree and has no git.")
    ap.add_argument("--diagnostic-out", type=Path, default=None, dest="diagnostic_out",
                    help="optional JSON path for the device-captured C-2 convergence history; instrumentation "
                         "only, written after acceptance/rollback and never used by the solve")
    ap.add_argument("--afm-manifest", type=Path, default=None, dest="afm_manifest",
                    help="READY afm-sweep-preflight@1 manifest; enables the spherical loading path")
    ap.add_argument("--afm-speed-um-s", type=float, default=None, dest="afm_speed_um_s",
                    help="one speed present in --afm-manifest; one Slurm job runs one speed/seed path")
    ap.add_argument("--afm-output", type=Path, default=None, dest="afm_output",
                    help="accepted F-delta JSON for one speed/seed path; required with --afm-manifest")
    args = ap.parse_args()

    afm_requested = any(
        value is not None for value in (args.afm_manifest, args.afm_speed_um_s, args.afm_output)
    )
    if afm_requested and not all(
        value is not None for value in (args.afm_manifest, args.afm_speed_um_s, args.afm_output)
    ):
        ap.error("--afm-manifest, --afm-speed-um-s and --afm-output must be supplied together")
    afm_path = (
        load_ready_afm_path(args.afm_manifest, seed=args.seed, speed_um_s=args.afm_speed_um_s)
        if afm_requested else None
    )
    if afm_path is not None and args.inner_solver != "explicit":
        ap.error("the AFM loading path currently requires the exact explicit solver; external tangents are "
                 "not assembled in accelerated candidates")
    if afm_path is not None and args.turgor_ramp != 0.0:
        ap.error("AFM indentation and a turgor ramp are distinct perturbations and cannot share one run")
    if afm_path is not None and not args.build_commit:
        ap.error("--build-commit is required for every AFM path; an unstamped force curve is not evidence")
    if afm_path is not None and (
        args.filaments != 70686 or args.membrane_subdiv != 8 or args.nucleus_subdiv != 3
    ):
        ap.error(
            "AFM paths require the C-2-accepted full-native configuration exactly: "
            "--filaments 70686 --membrane-subdiv 8 --nucleus-subdiv 3"
        )

    wp.init()
    dev = str(wp.get_device())
    if not wp.get_device().is_cuda:
        raise RuntimeError(f"GATE-B interior-column driver needs a CUDA GPU (I0-A). Resolved {dev!r} is not CUDA.")
    print(f"[gate-b-column] device={dev}  steps={args.steps}  dt={args.dt}s  inner={args.outer}  "
          f"solver={args.inner_solver} retries={args.max_inner_retries} "
          f"membrane_subdiv={args.membrane_subdiv}  nucleus_subdiv={args.nucleus_subdiv}", flush=True)

    # 1. Build the FULL native cell — membrane + cortex + cytosol (Biot field) + nucleus + pressure coupling.
    t0 = time.time()
    # Quantitative cortical-tension execution remains fail-closed upstream because no accepted
    # post-induction/dt-converged active-cortex state handoff exists (STATE (c) 3/17).  A READY manifest
    # can therefore currently represent only the explicitly labelled passive apparatus mechanism demo.
    cfg = CellConfig(
        n_filaments=args.filaments, with_myosin=False, overlap_free_cortex=True,
        with_membrane=True, membrane_subdivisions=args.membrane_subdiv,
        with_nucleus=True, nucleus_subdivisions=args.nucleus_subdiv,
        with_pressure=True, erm_radial_pairing=True, seed=args.seed,
    )
    cell = build_cell(cfg)
    for name, obj in (("grid", cell.grid), ("substrate", cell.substrate), ("pressure", cell.pressure),
                      ("membrane_bc", cell.membrane_bc), ("domain", cell.domain), ("nucleus", cell.nucleus),
                      ("membrane", cell.membrane)):
        if obj is None:
            raise RuntimeError(f"interior column needs the {name!r} compartment/primitive "
                               f"(build with with_pressure/with_membrane/with_nucleus=True)")
    print(f"[gate-b-column] built cell in {time.time() - t0:.1f}s: n_actin={cell.n_actin} n_total={cell.n_total} "
          f"nucleus_faces={cell.nucleus.n_faces} membrane_verts={cell.membrane.n_verts} "
          f"field_cells={int(np.prod(cell.grid.shape))}", flush=True)
    _banner(cell)

    # 1b. RESTING BASELINE — the physiological-baseline preloads (PI rule), same as the production driver.
    if args.preload:
        pre_cortex = _preload_cortex_pretension(cell, args.cortex_prestrain)
        pre_erm = _preload_erm_resting_balance(cell) if cell.membrane and cell.membrane.n_erm else None
        print(f"[gate-b-column] RESTING PRELOAD ON: cortex_prestrain={args.cortex_prestrain:g} ({pre_cortex}); "
              f"ERM balance={pre_erm}", flush=True)

    # 2. Build the four interior-column state-owners over the built cell's REAL arrays/primitives.
    #
    #    CARD-5 (PI-approved 2026-08-09).  Until this change all three surface owners were handed the WHOLE
    #    `cell.pos_d`/`cell.f_d`, so each claimed every node in the cell and the charter's line — co-location
    #    in a shared array is NEVER a connection — was violated by construction (`STATE.md` (c) 4).  Each owner
    #    now addresses its OWN contiguous block through `component_block_views`.
    #
    #    THESE ARE SLICE VIEWS, NOT PRIVATE ALLOCATIONS, and this driver may not report them as private.  The
    #    incumbent's inner solve relaxes ONE global `pos_d`; giving a component a private allocation without
    #    splitting that solve would make it a passenger — distinct arrays the solve never updates, which is
    #    precisely the T2 negative control.  What the views buy is exclusive ADDRESSING: a kernel launched over
    #    the membrane's view cannot reach a cortex node, and each owner's snapshot / reject-restore covers its
    #    own block instead of silently restoring the whole cell three times over.
    blocks = {
        "cortex": (0, int(cell.n_actin)),
        "nucleus": (int(cell.nucleus.node_off), int(cell.nucleus.n_verts)),
        "membrane": (int(cell.membrane.node_off), int(cell.membrane.n_verts)),
    }
    assert_disjoint_blocks(blocks, n_total=int(cell.n_total))
    membrane_pos_d, membrane_force_d = component_block_views(
        cell.pos_d, cell.f_d, node_off=blocks["membrane"][0], n_nodes=blocks["membrane"][1])
    cortex_pos_d, cortex_force_d = component_block_views(
        cell.pos_d, cell.f_d, node_off=blocks["cortex"][0], n_nodes=blocks["cortex"][1])
    nucleus_pos_d, nucleus_force_d = component_block_views(
        cell.pos_d, cell.f_d, node_off=blocks["nucleus"][0], n_nodes=blocks["nucleus"][1])

    # The nucleus rupture kernel indexes `position_d[faces[t, k]]`, so a LOCAL position view needs a LOCAL
    # face table.  `cell.nucleus.faces_d` stays untouched — the incumbent still reads it with global indices.
    pos_np = cell.pos_d.numpy()
    dx = float(cell.grid.dx)
    m_verts, m_faces = _local_faces_verts(pos_np, cell.membrane.faces_d, cell.membrane.node_off, cell.membrane.n_verts)
    n_verts, n_faces = _local_faces_verts(pos_np, cell.nucleus.faces_d, cell.nucleus.node_off, cell.nucleus.n_verts)
    nucleus_faces_local_d = wp.array(np.ascontiguousarray(n_faces, np.int32), dtype=wp.int32, device=dev)

    membrane_owner = build_native_surface_owner(
        name="membrane", position_d=membrane_pos_d, force_d=membrane_force_d, device=dev)
    cortex_owner = build_native_surface_owner(
        name="cortex", position_d=cortex_pos_d, force_d=cortex_force_d, device=dev)
    cytosol_owner = build_fluid_field_owner(
        grid=cell.grid, solver=BiotSubstrateFluidSolver(cell.substrate), device=dev)
    nucleus_owner = build_core_surface_owner(
        position_d=nucleus_pos_d, force_d=nucleus_force_d, faces_d=nucleus_faces_local_d,
        a0_d=cell.nucleus.a0_d, ruptured_d=cell.nucleus.ruptured_d,
        eps_rupture=float(cell.nucleus.eps_rupt), device=dev)

    # The definition of done, checked on the built objects rather than declared: the three node-array owners
    # address disjoint bytes.  The cytosol is a FIELD owner and holds no node arrays, so it is not a party to
    # this check.  Overlapping views would pass a pointer-identity test and are refused here by byte span.
    assert_component_state_disjoint([membrane_owner, cortex_owner, nucleus_owner])
    print(f"[gate-b-column] ARRAY OWNERSHIP: membrane/cortex/nucleus address DISJOINT blocks of the global "
          f"arrays — {', '.join(f'{k}=[{o}, {o + n})' for k, (o, n) in sorted(blocks.items()))} of "
          f"{int(cell.n_total)} nodes. These are SLICE VIEWS, not private allocations: the incumbent's inner "
          f"solve still relaxes one global pos_d, so `arrays_private` is FALSE and STATE.md (c) 4 stands "
          f"until the PI says otherwise.", flush=True)

    if args.ownership_gate:
        _ownership_gate(cell, blocks, [membrane_owner, cortex_owner, nucleus_owner],
                        args.ownership_out, dev, time.time() - t0, args.build_commit)
        return

    pi0_pa = afm_path.pi0_pa if afm_path is not None else 40.0
    additional_force_accumulators = ()
    additional_stiffness_pn_per_um = 0.0
    additional_component_bindings = ()
    additional_connector_bindings = ()
    apparatus = None
    exterior_medium = None
    if afm_path is not None:
        if args.steps != 1:
            raise ValueError("AFM mode uses exactly one accepted zero-load baseline step; pass --steps 1")
        contact_spec = AFMIndenterRuntimeSpec(
            radius_um=afm_path.indenter_radius_um,
            surface_contact_gap_um=afm_path.surface_contact_gap_um,
            stiffness_pn_per_um=afm_path.contact_stiffness_pn_per_um,
            provenance=f"{afm_path.apparatus_source}; {afm_path.contact_source}",
        )
        initial_centre = derive_force_free_initial_centre(
            m_verts,
            centre_contact_distance_um=contact_spec.centre_contact_distance_um,
        )
        apparatus = AFMIndenterRuntime(
            membrane_position_d=membrane_pos_d,
            membrane_force_d=membrane_force_d,
            initial_centre_um=initial_centre,
            spec=contact_spec,
            device=dev,
        )
        apparatus_hook = AFMIndenterForceHook(
            runtime=apparatus,
            membrane_node_offset=blocks["membrane"][0],
            total_node_count=int(cell.n_total),
        )

        # T10 owns a 642-point native exterior quadrature. Recursive icosphere construction preserves all
        # previous vertices as a prefix, so subdivision-8 membrane mechanics can carry the already-gated
        # subdivision-3 exterior operator without down-counting the membrane component itself.
        n_medium = n_vert_for_subdivisions(3)
        if int(cell.membrane.n_verts) < n_medium:
            raise ValueError("AFM exterior holder requires at least the T10 subdivision-3 quadrature")
        medium_local_index = np.arange(n_medium, dtype=np.int32)
        medium_global_index = medium_local_index + blocks["membrane"][0]
        medium_position = np.ascontiguousarray(m_verts[medium_local_index], dtype=np.float64)
        medium_epsilon_um = derive_blob_epsilon(medium_position)
        medium_settings = ExteriorStokesMediumSettings(
            viscosity_pa_s=float(afm_path.exterior_viscosity_pa_s),
            blob_epsilon_um=float(medium_epsilon_um),
        )
        exterior_medium = build_exterior_stokes_medium(
            settings=medium_settings,
            surface_index=medium_global_index,
            initial_position=medium_position,
            device=dev,
        )
        medium_connector = MembraneMediumTraction(delegate=exterior_medium)
        positive_depths = afm_path.depths_um[1:]
        dt_values = tuple(
            (depth - previous) / afm_path.speed_um_s
            for previous, depth in zip(afm_path.depths_um[:-1], positive_depths, strict=True)
        )
        resistance_bound = resistance_spectral_upper_bound(
            medium_position,
            epsilon=medium_epsilon_um,
            viscosity=afm_path.exterior_viscosity_pa_s,
        )
        medium_stiffness_bound = resistance_bound / min((float(args.dt), *dt_values))
        additional_force_accumulators = (apparatus_hook, medium_connector.accumulate)
        additional_stiffness_pn_per_um = (
            apparatus_hook.stiffness_bound_pn_per_um + medium_stiffness_bound
        )
        additional_component_bindings = (
            (AFM_INDENTER_COMPONENT, apparatus),
            (EXTERIOR_MEDIUM_COMPONENT, medium_connector),
        )
        additional_connector_bindings = (
            (AFM_MEMBRANE_CONTACT, apparatus),
            (MEMBRANE_MEDIUM_CONNECTOR, medium_connector),
        )
        print(
            "[gate-b-column] AFM APPARATUS READY: zero-depth centre derived force-free; "
            f"T10 exterior quadrature={n_medium}, epsilon={medium_epsilon_um:.6g} um; "
            f"roundoff-enclosed added stiffness={additional_stiffness_pn_per_um:.6g} pN/um",
            flush=True,
        )

    # 3. Build the three fluid coupling connectors over the built cell's REAL primitives (bidirectional+adjoint).
    #    node-volume immersed quadrature weights: triangle-mesh lumped one-ring area × dx for the two surfaces
    #    (membrane, nucleus); a uniform dx³ immersed cell for the cortex network nodes.  Every weight array is
    #    the length of ITS OWN component's block, which is what the sliced views now make consistent:
    #    `PressureCoupling.accumulate` launches `dim=node_pos.shape[0]` and indexes `node_volume[i]`, so the
    #    global position array it used to be handed ran that launch off the end of the weights.
    membrane_node_vol = wp.array(lumped_surface_control_volumes(m_verts, m_faces, dx), dtype=wp.float64, device=dev)
    nucleus_node_vol = wp.array(lumped_surface_control_volumes(n_verts, n_faces, dx), dtype=wp.float64, device=dev)
    cortex_node_vol = wp.full(int(cell.n_actin), dx ** 3, dtype=wp.float64, device=dev)

    cortex_transfer = build_cortex_cytosol_transfer(pressure_coupling=cell.pressure, node_volume_d=cortex_node_vol)
    membrane_boundary = MovingSemipermeableFluidBoundaryFacade(
        delegate=MembranePressureFluxAdjointBoundary(
            cell.membrane_bc, cell.pressure, membrane_node_vol, _NoOpParticipant(), RecordingLedgerStub(),
            osmotic_difference=pi0_pa,
        ))
    no_flux_bc = NucleusNoFluxBC(cell.grid, mobility=float(cell.substrate.mobility))
    nucleus_boundary = MovingImpermeableFluidBoundaryFacade(
        delegate=NucleusPressureAdjointBoundary(
            cell.domain, cell.pressure, no_flux_bc, nucleus_node_vol, _NoOpParticipant(), RecordingLedgerStub()))

    # 4. Compose the interior-column slice — the four owners + three connectors under one accepted-step clock.
    residual_d = wp.zeros(cell.grid.shape, dtype=wp.float64, device=dev)
    cytosol_endpoint = cytosol_owner.make_immersed_transfer_endpoint(residual_d)
    slc = build_interior_column_slice(
        membrane_owner=membrane_owner, cortex_owner=cortex_owner, cytosol_owner=cytosol_owner,
        nucleus_owner=nucleus_owner, cortex_transfer=cortex_transfer, membrane_boundary=membrane_boundary,
        nucleus_boundary=nucleus_boundary, cytosol_endpoint=cytosol_endpoint, base_seed=args.seed, device=dev,
        architecture=afm_experiment_architecture() if afm_path is not None else None,
        additional_component_bindings=additional_component_bindings,
        additional_connector_bindings=additional_connector_bindings)
    print(f"[gate-b-column] interior column composed: {len(slc.transaction.world.participants)} transaction "
          f"participants (membrane+cortex+cytosol+nucleus + surface_porous_transfer + membrane/nucleus "
          f"cytosol boundaries), one accepted-step predicate.", flush=True)

    # 5. Drive the accepted-step transaction; the driver's inner solve OWNS the coupled mechanics (no double
    #    count — see the header).  solve = the inner mechanical relax (Biot subcycles + pressure coupling +
    #    membrane/cortex/nucleus forces).  The slice adds the field snapshot/restore + the accepted-gated
    #    nucleus rupture + the clock advance, all four compartments under one predicate.
    inner_solve = make_inner_solve(
        cell,
        n_inner=args.outer,
        reshape_every=20,
        max_inner_retries=args.max_inner_retries,
        capture_candidate=args.diagnostic_out is not None,
        inner_solver=args.inner_solver,
        implicit_line_search_steps=args.line_search_steps,
        line_search_objective=args.line_search_objective,
        implicit_coarse_modes=args.coarse_modes,
        implicit_coarse_iterations=args.coarse_iterations,
        additional_force_accumulators=additional_force_accumulators,
        additional_stiffness_pn_per_um=additional_stiffness_pn_per_um,
    )
    accepted_d = wp.zeros(1, dtype=wp.int32, device=dev)
    last_inner_report = None
    ramp_pa = float(args.turgor_ramp)
    current_dt_phys = float(args.dt)
    current_afm_depth_um = 0.0

    def solve() -> None:
        nonlocal last_inner_report
        # The load is candidate state: CellTransaction has already snapshotted the field, so rejection
        # restores it.  Applying the ramp before ``slc.step`` made the raised pressure authoritative even
        # when mechanics rejected the step.
        if ramp_pa > 0.0:
            wp.launch(_ramp_turgor_fluid_kernel, dim=cell.grid.shape,
                      inputs=[cell.grid.p, cell.grid.mask, wp.float64(ramp_pa)], device=dev)
        if apparatus is not None and exterior_medium is not None:
            apparatus.set_candidate_depth(current_afm_depth_um)
            exterior_medium.set_step(current_dt_phys)
        last_inner_report = inner_solve(current_dt_phys)
        wp.copy(accepted_d, last_inner_report.converged_d)

    if ramp_pa > 0.0:
        print(f"[gate-b-column] LOAD RAMP ON: +{ramp_pa:g} Pa/step on the FLUID interior "
              f"(Π {pi0_pa:g}→{pi0_pa + ramp_pa * args.steps:g} Pa over {args.steps} steps) — the coupled pressure "
              f"traction rides onto membrane/cortex/envelope; high-strain faces may rupture (accepted-gated).",
              flush=True)

    # 6. Physical-step event loop — only converged candidates become accepted interior-column state.
    print(f"\n{'step':>4} {'p_mean_Pa':>10} {'p_max_Pa':>10} {'memb_flux':>12} {'noflux_teeth':>13} "
          f"{'ruptured':>9} {'maxPF':>10}", flush=True)
    accepted_steps = 0
    for step in range(args.steps):
        accepted_d.zero_()
        slc.step(solve, dt_phys=current_dt_phys, accepted_d=accepted_d)

        st = _column_stats(cell, slc)  # out-of-hot-loop readback (I0-A honoured inside the loop)
        max_pf = float(inner_solve.convergence_force_d.numpy()[0])
        print(f"{step:>4d} {st['p_mean_pa']:>10.4g} {st['p_max_pa']:>10.4g} {st['membrane_flux']:>12.4g} "
              f"{st['no_flux_teeth']:>13.4g} {st['ruptured_faces']:>9d} {max_pf:>10.4g}", flush=True)
        accepted = bool(int(accepted_d.numpy()[0]))
        if args.diagnostic_out is not None:
            diagnostic_path = _c2_diagnostic_step_path(
                args.diagnostic_out, step_index=step, total_steps=args.steps)
            _write_c2_diagnostic(
                diagnostic_path, args=args, cell=cell, inner_solve=inner_solve,
                report=last_inner_report, max_projected_force_pn=max_pf, accepted=accepted,
                physical_step_index=step,
            )
        if not accepted:
            candidate_residual = (
                float(last_inner_report.residual_d.numpy()[0]) if last_inner_report is not None else float("nan")
            )
            raise RuntimeError(
                "C-2 BLOCKED: the coupled inner solve did not converge, so the candidate was rolled back "
                f"and biological time did not advance (step={step}, maxPF={max_pf:.6g} pN, "
                f"candidate residual={candidate_residual:.6g} pN). No F-delta or cortical-tension "
                "magnitude may be read from this trajectory."
            )
        accepted_steps += 1

    if afm_path is not None and apparatus is not None and exterior_medium is not None:
        tolerance_um = float(np.sqrt(np.finfo(np.float64).eps) * cell.convergence_length_um)
        loading_time_s = 0.0

        def afm_readout(depth_um: float, dt_phys_s: float) -> dict[str, object]:
            reaction = apparatus.reaction_d.numpy()
            dt_mu = float(last_inner_report.dt_mu_d.numpy()[0])
            return {
                "depth_um": float(depth_um),
                "loading_time_s": float(loading_time_s),
                "dt_phys_s": float(dt_phys_s),
                "reaction_z_pn": float(reaction[0]),
                "contact_count": int(round(float(reaction[1]))),
                "max_projected_force_pn": float(inner_solve.convergence_force_d.numpy()[0]),
                "projected_force_tolerance_pn": tolerance_um / dt_mu,
                "max_displacement_um": float(last_inner_report.max_displacement_d.numpy()[0]),
                "displacement_tolerance_um": tolerance_um,
                "reported_iterations": int(last_inner_report.iters_d.numpy()[0]),
                "accepted": True,
            }

        afm_records = [afm_readout(0.0, float(args.dt))]
        previous_depth = 0.0
        for point_index, depth_um in enumerate(afm_path.depths_um[1:], start=1):
            current_afm_depth_um = float(depth_um)
            current_dt_phys = (current_afm_depth_um - previous_depth) / afm_path.speed_um_s
            accepted_d.zero_()
            slc.step(solve, dt_phys=current_dt_phys, accepted_d=accepted_d)
            accepted = bool(int(accepted_d.numpy()[0]))
            max_pf = float(inner_solve.convergence_force_d.numpy()[0])
            if not accepted:
                blocked = {
                    "schema": "afm-force-indentation-path@1",
                    "status": "REJECTED_ROLLED_BACK",
                    "quantitative_claim": "BLOCKED",
                    "seed": afm_path.seed,
                    "speed_um_s": afm_path.speed_um_s,
                    "accepted_points": [],
                    "rejected_candidate": {
                        "depth_um": current_afm_depth_um,
                        "dt_phys_s": current_dt_phys,
                        "max_projected_force_pn": max_pf,
                        "candidate_residual_pn": float(last_inner_report.residual_d.numpy()[0]),
                    },
                    "reason": "Rejected candidates are never exposed as a partial material curve.",
                }
                args.afm_output.parent.mkdir(parents=True, exist_ok=True)
                args.afm_output.write_text(json.dumps(blocked, indent=2) + "\n", encoding="utf-8")
                raise RuntimeError(
                    "AFM PATH REJECTED: candidate rolled back and no partial F-delta curve was emitted "
                    f"(depth={current_afm_depth_um:g} um, maxPF={max_pf:.6g} pN)"
                )
            loading_time_s += current_dt_phys
            afm_records.append(afm_readout(current_afm_depth_um, current_dt_phys))
            previous_depth = current_afm_depth_um
            print(
                f"[afm] point={point_index} depth={current_afm_depth_um:.6g} um "
                f"reaction={afm_records[-1]['reaction_z_pn']:.6g} pN "
                f"contacts={afm_records[-1]['contact_count']} maxPF={max_pf:.6g} pN ACCEPTED",
                flush=True,
            )

        afm_payload = {
            "schema": "afm-force-indentation-path@1",
            "status": "ACCEPTED_PATH",
            "quantitative_claim": (
                "ELIGIBLE_FOR_PROTOCOL_AGGREGATION"
                if afm_path.result_class == "quantitative"
                else "MECHANISM_DEMO_NOT_QUANTITATIVE"
            ),
            "device": dev,
            "build_commit": args.build_commit,
            "seed": afm_path.seed,
            "speed_um_s": afm_path.speed_um_s,
            "population": {
                "n_filaments": int(cell.cfg.n_filaments),
                "n_actin_nodes": int(cell.n_actin),
                "n_total_nodes": int(cell.n_total),
                "membrane_subdivisions": int(args.membrane_subdiv),
                "membrane_vertices": int(cell.membrane.n_verts),
                "nucleus_subdivisions": int(args.nucleus_subdiv),
                "full_compartments": True,
            },
            "protocol": {
                "pi0_pa": afm_path.pi0_pa,
                "pi0_source": afm_path.pi0_source,
                "indenter_radius_um": afm_path.indenter_radius_um,
                "surface_contact_gap_um": afm_path.surface_contact_gap_um,
                "contact_stiffness_pn_per_um": afm_path.contact_stiffness_pn_per_um,
                "contact_source": afm_path.contact_source,
                "apparatus_source": afm_path.apparatus_source,
                "exterior_viscosity_pa_s": afm_path.exterior_viscosity_pa_s,
                "exterior_quadrature_points": int(exterior_medium.n_quadrature),
                "reaction_channel": 0,
            },
            "points": afm_records,
            "scope": (
                "One accepted speed/seed loading path. Cortical tension requires >=3-seed aggregation and "
                "a predeclared inversion; this record alone is not gamma."
            ),
        }
        args.afm_output.parent.mkdir(parents=True, exist_ok=True)
        args.afm_output.write_text(json.dumps(afm_payload, indent=2) + "\n", encoding="utf-8")
        print(f"[afm] ACCEPTED PATH -> {args.afm_output}", flush=True)

    final = _column_stats(cell, slc)
    print(f"\n[gate-b-column] DONE: {accepted_steps} accepted steps. Interior column drove membrane→cytosol→nucleus "
          f"as one transaction: field Π_mean={final['p_mean_pa']:.4g} Pa, membrane permeation flux "
          f"{final['membrane_flux']:.4g}, no-flux teeth {final['no_flux_teeth']:.4g} (impermeable envelope held "
          f"at 0 applied flux), {final['ruptured_faces']} ruptured envelope faces. Coupled state EMERGED from "
          f"the pressure field, no imposed prestress. (Nucleus + osmotic magnitudes PROVISIONAL — a PI decision.)",
          flush=True)


if __name__ == "__main__":
    main()

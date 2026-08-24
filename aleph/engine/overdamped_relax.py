r"""One component advances its OWN positions from its OWN forces — the missing keystone.

WHY THIS IS THE KEYSTONE.  Every engine component owns ``position_d`` and ``force_d`` and **nothing
integrates one from the other**; only the feature-frozen incumbent driver has an integrator, over its
own global arrays.  That single absence appeared three times on 2026-08-11 wearing three faces:

  * ``SOLVE_COUPLED`` is an empty phase — because there is nothing to couple, only forces to sum;
  * a component cannot own private arrays — because owning them means being able to move them;
  * there is no multirate — because a component with no clock cannot be given a different one.

It is the same thing each time, and it is why FMI puts ``doStep`` on the first line of its contract: a
part that cannot step is not a component of a composition, it is a subroutine of someone else's solver.

NOTHING HERE IS CHOSEN.  Two numbers are needed and both are derived:

``γ`` — the per-point overdamped drag, from :func:`~aleph.laws.units.fiber_point_drag`, which is
NF2007 §5.2's single fiber mobility ``μ = log(L_h/δ)/(3πηL)`` with ``μ_p = (p+1)μ``.  ``laws/relax.py``
already derives it this way for a ``FiberNetwork``; this reads the same law off a component's own
topology.  **Every filament-graph component shares it** — cortex, sf_arc, microtubule,
intermediate_filament, lamellipodium, filopodium — so six of the ten unbound components need no new
physics to become steppable.  Membrane, cytosol and nucleus are different physics and are NOT covered
here; claiming otherwise would be a lumped drag wearing a filament's name.

``dt_max`` — the EXACT explicit-Euler stability bound for ``γẋ = −Kx``, which is ``2γ/λ_max``, with
``λ_max`` bounded by Gershgorin: the largest per-node row sum of the assembled stiffness.  No safety
factor is applied, because a safety factor is a chosen constant.  ``laws/relax.py`` uses ``0.5``; that
belongs to its caller's judgement, not to a stability statement.

The Gershgorin row sum is the same quantity the 2026-07-28f native spectrum row measured — *"the
production explicit step is stable at full native, and ``kmax`` was never a bound; the ratio to it IS
the per-node bond accumulation."*  This computes that accumulation instead of assuming a bound.

WHAT THIS IS NOT.  Not a coupled solve — one component, its own arrays, forces held fixed.  Not an
implicit step, so a stiff configuration will demand a small ``dt`` and this will say so rather than
quietly going unstable.  Not a claim that any resulting state is an equilibrium.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — γ [pN·s·µm⁻¹], force [pN], dt [s], so ``dt·f/γ`` is µm.  ``λ_max`` [pN/µm]; ``2γ/λ``
    is seconds.
  * boundary — ``dt = 0`` moves nothing; a ``dt`` at or above the derived bound RAISES rather than
    integrating unstably, and the message carries the bound so the caller can subcycle knowingly.
  * conservation — positions only; forces are the caller's, untouched.  This writes no committed
    state: the owner's snapshot/rollback machinery still decides whether the move survives.
  * CFL/precision — the bound is exact for explicit Euler on a linear overdamped system, and
    Gershgorin is an upper bound on ``λ_max``, so the admitted ``dt`` is conservative, never optimistic.
    The system is NOT linear — a stretched bond carries a transverse geometric stiffness
    ``k_t = k(L−r0)/L`` absent from the scalar row sum — so that the row sum still bounds the true
    3N×3N tangent is a claim, and it is MEASURED rather than assumed
    (``test_the_row_sum_bounds_the_true_nonlinear_tangent``: worst ratio 0.9976 over 3,000 stiff
    high-coordination configurations stretched up to 50×).  Nearly tight is the wanted outcome: a
    4×-loose bound would buy nothing and cost 4× the subcycles.
  * sign sense — ``x ← x + (dt/γ)·f``: a node moves ALONG the force on it.  The opposite sign would
    grow every mode and is caught by the descent test.
  * measurement protocol — device-resident; γ and the bound are derived once at build from host
    topology, which is where topology work belongs.

engine units: length µm, force pN, time s, stiffness pN/µm.  Runtime: NVIDIA Warp on CUDA only.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.laws import units as U

__all__ = [
    "OverdampedStep",
    "derive_gershgorin_lambda_max",
    "derive_point_drag",
    "build_overdamped_step",
]


@wp.kernel
def _overdamped_advance_kernel(
    position: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
    dt_over_gamma: wp.float64,
) -> None:
    """``x ← x + (dt/γ)·f`` — one node per thread, positions only."""
    t = wp.tid()
    position[t] = position[t] + dt_over_gamma * force[t]


def derive_point_drag(
    *,
    segment_rest_um: npt.NDArray[np.float64],
    n_nodes: int,
    n_filaments: int,
) -> float:
    """Per-point overdamped drag γ [pN·s·µm⁻¹] from NF2007 §5.2, read off a component's topology.

    The same route ``laws/relax.py::_default_gamma`` takes for a ``FiberNetwork``: the mean fiber
    contour length and the mean segments-per-fiber determine ``μ_p = (p+1)μ``, and ``γ = 1/μ_p``.

    Args:
        segment_rest_um: per-segment rest lengths [µm]; their mean × segments-per-fiber is the mean
            contour length.
        n_nodes: node count of this component.
        n_filaments: filament count, so segments-per-fiber is measured rather than assumed.

    Returns:
        γ [pN·s·µm⁻¹].

    Raises:
        ValueError: on an empty topology — a component with no segments has no fiber to take a
            mobility from, and returning a default would be inventing one.
    """
    if segment_rest_um.size == 0 or n_filaments <= 0 or n_nodes <= 0:
        raise ValueError(
            "cannot derive gamma from an empty topology; a component with no segments has no fiber "
            "mobility, and a default here would be a chosen constant"
        )
    seg_mean = float(np.mean(segment_rest_um))
    segments_per_fiber = max(int(round(segment_rest_um.size / n_filaments)), 1)
    length_um = max(seg_mean * segments_per_fiber, 1e-3)
    return float(U.fiber_point_drag(length_um, segments_per_fiber))


def derive_gershgorin_lambda_max(
    *,
    links: npt.NDArray[np.integer],
    link_k: npt.NDArray[np.float64],
    n_nodes: int,
    bend_triples: npt.NDArray[np.integer] | None = None,
    bend_alpha: npt.NDArray[np.float64] | None = None,
    extra_pairs: tuple[tuple[npt.NDArray[np.integer], npt.NDArray[np.float64]], ...] = (),
) -> float:
    """Largest per-node stiffness row sum — an upper bound on ``λ_max`` [pN/µm].

    Gershgorin: every eigenvalue lies within ``Σ_j |K_ij|`` of a diagonal entry, so the largest row
    sum bounds ``λ_max``.  A Hookean link contributes ``2k`` to each endpoint's row (its own diagonal
    plus the off-diagonal to its partner).  A Cytosim bending triple is a 4th-difference stencil whose
    row sum is ``16α`` at the middle node (``λ_max = 16κ/seg³`` and ``α = κ/seg³``, the same identity
    ``laws/relax.py::_bending_cfl_dt`` uses).

    ``extra_pairs`` carries any FURTHER ``(endpoints, k)`` set acting on the same node array — e.g.
    ``sf_arc``'s α-actinin ``arc_joints``, whose 4.6e5 pN/µm is ~460× the axial k.  Such a term is not
    optional detail: leave it out and the result is not an upper bound at all.

    This is an UPPER bound, so the ``dt`` it admits is conservative.  The 2026-07-28f native spectrum
    row measured exactly this quantity and found the previously-assumed ``kmax`` was never a bound.
    """
    rows = np.zeros(int(n_nodes), dtype=np.float64)
    for endpoints, stiffness in ((links, link_k), *extra_pairs):
        if endpoints is None or not np.asarray(endpoints).size:
            continue
        pairs = np.asarray(endpoints, dtype=np.int64).reshape(-1, 2)
        k = np.asarray(stiffness, dtype=np.float64).reshape(-1)
        np.add.at(rows, pairs[:, 0], 2.0 * k)
        np.add.at(rows, pairs[:, 1], 2.0 * k)
    if bend_triples is not None and bend_triples.size and bend_alpha is not None:
        triples = np.asarray(bend_triples, dtype=np.int64).reshape(-1, 3)
        alpha = np.asarray(bend_alpha, dtype=np.float64).reshape(-1)
        for column in range(3):
            np.add.at(rows, triples[:, column], 16.0 * alpha)
    return float(rows.max()) if rows.size else 0.0


@dataclass(frozen=True, slots=True)
class OverdampedStep:
    """A component's own explicit overdamped integrator, with its derived stability bound.

    Attributes:
        gamma_pn_s_per_um: the NF2007 per-point drag γ.
        lambda_max_pn_per_um: the Gershgorin bound on the assembled stiffness.
        n_nodes: the array length this step is valid for.
        launch: injected launcher (a CUDA-free gate substitutes a recorder).
    """

    gamma_pn_s_per_um: float
    lambda_max_pn_per_um: float
    n_nodes: int
    launch: object = wp.launch

    def __post_init__(self) -> None:
        if not self.gamma_pn_s_per_um > 0.0:
            raise ValueError("gamma must be positive")
        if self.lambda_max_pn_per_um < 0.0:
            raise ValueError("the Gershgorin bound cannot be negative")

    @property
    def dt_max_s(self) -> float:
        """The EXACT explicit-Euler stability bound ``2γ/λ_max`` [s]; ``inf`` for a stiffness-free set.

        No safety factor. A safety factor is a chosen constant, and choosing one here would put it
        inside a stability statement where a reader would take it for physics.
        """
        if self.lambda_max_pn_per_um <= 0.0:
            return float("inf")
        return 2.0 * self.gamma_pn_s_per_um / self.lambda_max_pn_per_um

    def advance(self, position_d: wp.array, force_d: wp.array, dt_phys: float) -> None:
        """Advance positions in place. Refuses an unstable ``dt`` instead of integrating it.

        Raises:
            ValueError: on a negative ``dt``, on mismatched arrays, or on ``dt >= dt_max_s`` — the
                message carries the bound, so a caller that needs a larger step knows exactly how many
                subcycles it owes rather than discovering instability in a trajectory.
        """
        if dt_phys < 0.0:
            raise ValueError("dt_phys must be nonnegative")
        if dt_phys == 0.0:
            return
        if dt_phys >= self.dt_max_s:
            raise ValueError(
                f"dt_phys={dt_phys!r} s is at or above the derived explicit-stability bound "
                f"{self.dt_max_s!r} s (gamma={self.gamma_pn_s_per_um!r} pN·s/µm, "
                f"lambda_max={self.lambda_max_pn_per_um!r} pN/µm). Subcycle "
                f"{int(dt_phys / self.dt_max_s) + 1} times, or use an implicit step — this will not "
                "integrate an unstable configuration quietly."
            )
        if int(position_d.shape[0]) != self.n_nodes or int(force_d.shape[0]) != self.n_nodes:
            raise ValueError(
                f"this step was derived for {self.n_nodes} nodes; got "
                f"{int(position_d.shape[0])} positions and {int(force_d.shape[0])} forces"
            )
        self.launch(
            _overdamped_advance_kernel,
            dim=self.n_nodes,
            inputs=[position_d, force_d, wp.float64(dt_phys / self.gamma_pn_s_per_um)],
            device=str(position_d.device),
        )


def build_overdamped_step(
    topology: object,
    *,
    n_filaments: int,
    borrowed_pairs: tuple[tuple[npt.NDArray[np.integer], npt.NDArray[np.float64]], ...],
    launch: object = wp.launch,
) -> OverdampedStep:
    """Derive a filament-graph component's integrator from its own topology. No constant is supplied.

    Args:
        topology: anything exposing ``n_nodes``, ``links``, ``link_k``, ``link_r0`` and optionally
            ``bend_triples`` / ``bend_alpha`` — the shape ``SFMechanicsTopology`` already has.
        n_filaments: measured filament count, so segments-per-fiber is read rather than assumed.
        borrowed_pairs: stiffness a CONNECTOR scatters into this component's ``force_d`` while owning
            no array here — ``(node_pairs, k)`` per bound connector, in THIS component's node indices.
            **Required, with no default, on purpose.**  ``_refuse_unsummed_stiffness`` scans the
            topology object, so it is structurally blind to a connector: ``sf_cortex_transient``
            borrows ``sf_owner.position_d``/``force_d`` through ``BorrowedSegmentActorView`` and
            scatters α-actinin at 4.6e5 pN/µm, none of which appears on the topology at all.  A
            default of ``()`` would let that be forgotten silently — which is the exact defect that
            put ``arc_k`` outside the bound and admitted a ``dt`` 459× too large.  Pass the connector's
            CANDIDATE pairs, not its engaged ones: engagement is transient, and a bound that changes
            with state is not a bound.  ``()`` is a legal answer for a component with no borrower.
        launch: injected launcher.
    """
    n_nodes = int(topology.n_nodes)
    gamma = derive_point_drag(
        segment_rest_um=np.asarray(topology.link_r0, dtype=np.float64).reshape(-1),
        n_nodes=n_nodes, n_filaments=int(n_filaments),
    )
    extra = tuple(
        (getattr(topology, pairs), getattr(topology, stiffness))
        for pairs, stiffness in _EXTRA_PAIR_ARRAYS
        if getattr(topology, stiffness, None) is not None
    ) + tuple(borrowed_pairs)
    _refuse_unsummed_stiffness(topology)
    lambda_max = derive_gershgorin_lambda_max(
        links=np.asarray(topology.links), link_k=np.asarray(topology.link_k, dtype=np.float64),
        n_nodes=n_nodes,
        bend_triples=getattr(topology, "bend_triples", None),
        bend_alpha=getattr(topology, "bend_alpha", None),
        extra_pairs=extra,
    )
    return OverdampedStep(gamma, lambda_max, n_nodes, launch)


#: Pair-stiffness arrays beyond the axial ``links``/``link_k`` that act on the SAME node array and so
#: enter the same stiffness matrix.  ``arc_joints`` is α-actinin at 4.6e5 pN/µm — ~460× the axial k —
#: so omitting it does not make the bound loose, it makes it not a bound.
_EXTRA_PAIR_ARRAYS: tuple[tuple[str, str], ...] = (("arc_joints", "arc_k"),)

#: Stiffness arrays accounted for above; anything else ending `_k` is unsummed and refused.
_SUMMED_STIFFNESS = frozenset({"link_k", *(s for _, s in _EXTRA_PAIR_ARRAYS)})


def _refuse_unsummed_stiffness(topology: object) -> None:
    """Raise if the topology carries a stiffness array this builder does not sum.

    The bound is only a bound if EVERY term touching a node is in its row sum.  A silently-omitted
    stiffness admits a `dt` above the true stability limit, and the resulting trajectory looks like
    physics for a while — so a new stiffness must break the build, not the run.
    """
    unsummed = sorted(
        name
        for name in dir(topology)
        if name.endswith("_k") and name not in _SUMMED_STIFFNESS and not name.startswith("_")
        and getattr(topology, name, None) is not None
        and np.asarray(getattr(topology, name)).size
    )
    if unsummed:
        raise ValueError(
            f"topology carries stiffness arrays this integrator does not sum into its Gershgorin "
            f"bound: {unsummed}. Add each to `_EXTRA_PAIR_ARRAYS` with its endpoint pairs — an "
            f"omitted stiffness admits a dt above the true stability limit."
        )

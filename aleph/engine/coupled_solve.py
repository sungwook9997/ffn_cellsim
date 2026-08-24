r"""``SOLVE_COUPLED`` — the interface residual across every attributable cut, measured.

WHAT THIS CLOSES, AND WHAT IT DELIBERATELY DOES NOT.  ``dispatch.CandidatePhase.SOLVE_COUPLED`` has
been declared since the phase enum was written and **no production task has ever claimed it**.  The
cortex-motor driver says why in one sentence (``ac_gate_b_cortex_motor_native.py``): *"a genuinely
private force array needs the inner solve to relax TWO arrays together, and that is `SOLVE_COUPLED` —
an empty phase"*.  That is the reason components cannot own private arrays, and therefore the reason
the composed cell binds a handful of connectors instead of its declared graph.

The obvious move is to write the partitioned iteration straight away — each component relaxes its own
positions with the connector force held, forces are re-evaluated, repeat under Aitken acceleration
(the standard partitioned-FSI construction; preCICE's IQN-ILS is the same family).  **That is not what
this module does yet, and the reason is a constant.**  Relaxing a component's own positions needs that
component's overdamped mobility ``γ``, and ``γ`` is a *derived physical quantity per component*, not a
wiring detail.  Choosing one here would change the physics silently, which the charter forbids and
which is exactly how a tuning constant enters a codebase.

So this module does the measurement that must come first and that no run has ever made: **is the
interface out of balance at all, and if so, across which cut and by how much?**  An accelerator built
before that number exists would be an accelerator for an unmeasured residual.

WHAT IS MEASURED.  For each cut returned by :meth:`CellArchitecture.attributable_cuts` — the 26 pairs
whose two components are joined by exactly ONE power port, so a failure names its own connector — the
two sub-bodies' force resultants are reduced into the two INDEPENDENTLY sourced channels of a fresh
:class:`~aleph.engine.ledger.GlobalCellLedger` and scored by the existing D8 gate.  Nothing numerical
is supplied by this module: the tolerance is derived on the device from the run's own ``Σ|f|`` via
:meth:`GlobalCellLedger.derive_balance_tolerance`, and the term count comes from
:func:`~aleph.engine.ledger.assemble_balance_tolerance_ratio`.

WHY THE CUT AND NOT THE CONNECTOR.  ``add_body_force`` reduces a whole force array, so what it scores
is everything crossing the cut.  That resolves a single connector only where that connector is the
pair's sole coupling — 27 of 28 power ports, measured.  The exception (``cortex``/``membrane``, joined
by ``membrane_erm_cortex`` and ``membrane_cortex_contact``) is scored jointly and reported as the pair,
never attributed to one edge.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — every resultant is ``wp.vec3d`` [pN]; the residual is ``|reaction + traction|`` [pN];
    the tolerance is RELATIVE and carries no unit.
  * boundary — a cut whose two components are not both bound in this world is reported ``UNBOUND``,
    never as a balanced cut.  Absence and balance are different facts and a zero could mean either.
  * conservation/invariant — the two channels are separately accumulated from two never-merged arrays,
    so the gate cannot self-satisfy; this is the property ``ledger.py``'s header protects.
  * CFL/precision — no time integration and no position update happens here.  This is a read of one
    already-assembled candidate force state.
  * sign sense — a one-sided scatter, a double count and a sign flip are exactly what the gate fails
    on; the negative control injects one and requires the failure.
  * measurement protocol — one composed candidate assembly, then one device reduction per cut, then a
    single host readback after all cuts complete.  No readback inside a physical-time loop.

engine units: force pN, length µm.  Runtime: NVIDIA Warp on CUDA only.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import warp as wp

from aleph.engine.contracts import CellArchitecture
from aleph.engine.ledger import assemble_balance_tolerance_ratio, make_global_cell_ledger

__all__ = [
    "CutResidual",
    "InterfaceResidualReport",
    "collect_force_arrays",
    "measure_interface_residual",
]


@dataclass(frozen=True, slots=True)
class CutResidual:
    """One two-body cut's force-balance reading.

    Attributes:
        components: the cut, endpoints sorted.
        connectors: the power ports crossing it.  One name means the reading attributes to that
            connector; more than one means it does not, and the artifact must say the pair.
        status: ``BALANCED`` / ``IMBALANCED`` / ``UNBOUND`` (a component of this cut has no runtime
            here, so nothing was measured — NOT a balanced cut).
        residual_pn: ``|reaction + traction|`` [pN], or ``None`` when unbound.
        scale_pn: the ``Σ|f|`` bracket the tolerance was derived against [pN], or ``None``.
        tolerance_ratio: the dimensionless D8 floor ``γ_n`` used, or ``None``.
        reaction_pn / traction_pn: the two halves SEPARATELY [pN].  Recorded because a residual and a
            scale alone cannot say WHICH side carries the force, and on 2026-08-11 that ambiguity cost
            a wrong root cause: a cut reading 47.5 pN with cancellation 1.0 was read as the sf_arc side
            carrying an unpaired load, when the load was on the cortex side all along — the other
            body's half of a different connector.  Two numbers would have shown it immediately.
        isolated: whether this reading was taken with ONLY this cut's connector accumulating.  A
            whole-body reduction cannot attribute force to a connector once a body sits on more than
            one edge, so a non-isolated reading is a property of the BODIES, not of the connector.
    """

    components: tuple[str, str]
    connectors: tuple[str, ...]
    status: str
    residual_pn: float | None
    scale_pn: float | None
    tolerance_ratio: float | None
    connectors_bound: tuple[str, ...] = ()
    reaction_pn: float | None = None
    traction_pn: float | None = None
    isolated: bool = False

    @property
    def attributable(self) -> bool:
        """Whether a failure here names one connector rather than a pair."""
        return len(self.connectors) == 1

    @property
    def cancellation(self) -> float | None:
        """``residual / scale`` — how much of the force present actually cancelled. Dimensionless.

        This is the reading that says whether a ``BALANCED`` verdict means anything, and it needs no
        chosen threshold because it is a ratio the run supplies both halves of.  A genuine Newton pair
        gives ``≈ 0``: the parts are finite and their sum is not.  **A value near 1 means nothing
        cancelled** — the "sum" is as big as the "parts", which is what happens when one channel is
        empty or both are round-off noise.

        The repository has already paid for the absence of this number.  On 2026-07-28f the cortex
        lane's negative control turned out vacuous because at rest every NMII head is unbound, so
        myosin's force on actin is exactly zero — a control that cannot fail and a control that always
        passes are the same artifact, and only the scale distinguishes them.
        """
        if self.residual_pn is None or not self.scale_pn:
            return None
        return self.residual_pn / self.scale_pn

    @property
    def unwired(self) -> bool:
        """Whether this cut has no bound connector runtime crossing it.

        Structural and threshold-free: with nothing crossing the cut the two arrays cannot fail to
        balance, so a ``BALANCED`` here is an absence of wiring rather than a closed adjoint.
        """
        return not self.connectors_bound


@dataclass(frozen=True, slots=True)
class InterfaceResidualReport:
    """Every cut's reading plus the counts a caller should quote instead of re-deriving them."""

    cuts: tuple[CutResidual, ...]

    @property
    def measured(self) -> tuple[CutResidual, ...]:
        """Cuts where both components were bound, so a number exists."""
        return tuple(c for c in self.cuts if c.status != "UNBOUND")

    @property
    def imbalanced(self) -> tuple[CutResidual, ...]:
        """Cuts whose two channels did not cancel within the derived floor."""
        return tuple(c for c in self.cuts if c.status == "IMBALANCED")

    @property
    def scoreable(self) -> tuple[CutResidual, ...]:
        """Measured cuts that COULD say something: a connector is bound and a cancellation exists.

        Structural only. It deliberately does not grade the cancellation — deciding that 0.3 is
        "cancelled" and 0.6 is not would be a threshold nobody derived, and this module's whole point
        is that the reader gets the ratio and judges it. What it does exclude is the two cases where
        no judgement is possible at all: nothing crossing the cut, and no force to cancel.
        """
        return tuple(c for c in self.measured if not c.unwired and c.cancellation is not None)

    @property
    def worst(self) -> CutResidual | None:
        """The largest measured residual, or ``None`` when nothing was measured."""
        scored = [c for c in self.measured if c.residual_pn is not None]
        return max(scored, key=lambda c: c.residual_pn or 0.0) if scored else None

    def as_dict(self) -> dict[str, Any]:
        """A JSON-ready summary for a ``run-record@2`` artifact."""
        return {
            "n_cuts_declared": len(self.cuts),
            "n_measured": len(self.measured),
            "n_unbound": len(self.cuts) - len(self.measured),
            "n_imbalanced": len(self.imbalanced),
            "n_unwired": sum(1 for c in self.measured if c.unwired),
            "n_scoreable": len(self.scoreable),
            "worst_residual_pn": self.worst.residual_pn if self.worst else None,
            "worst_cut": list(self.worst.components) if self.worst else None,
            "cuts": [
                {
                    "components": list(c.components),
                    "connectors": list(c.connectors),
                    "attributable": c.attributable,
                    "status": c.status,
                    "residual_pn": c.residual_pn,
                    "reaction_pn": c.reaction_pn,
                    "traction_pn": c.traction_pn,
                    "isolated": c.isolated,
                    "scale_pn": c.scale_pn,
                    "cancellation": c.cancellation,
                    "tolerance_ratio": c.tolerance_ratio,
                    "connectors_bound": list(c.connectors_bound),
                    "unwired": c.unwired,
                }
                for c in self.cuts
            ],
            "scope": (
                "the interface residual of ONE assembled candidate force state. NOT a convergence "
                "claim, NOT a physics magnitude, and NOT a statement that any component reached "
                "equilibrium — no position was updated here. READ `cancellation` BEFORE `status`: a "
                "BALANCED cut whose cancellation is near 1 cancelled nothing and is evidence of "
                "absent force or absent wiring, not of a closed adjoint."
            ),
        }


def collect_force_arrays(*objects: object) -> tuple[wp.array, ...]:
    """Every distinct ``vec3d`` device force array reachable on the supplied runtimes.

    Distinctness is by device POINTER, not by object identity: two views onto one allocation must be
    collected once, or a leave-one-in pass would zero the same storage twice and count it twice.

    Public because the isolation pass belongs to the caller — only the driver knows which runtimes make
    up its world — while the definition of "a watched force array" belongs here, beside the reduction
    that consumes it.  Keeping them apart is how the driver and the gate came to disagree once already.
    """
    found: list[wp.array] = []
    seen: set[int] = set()
    for obj in objects:
        if obj is None:
            continue
        for attr in ("force_d", "force", "node_force_d"):
            arr = getattr(obj, attr, None)
            if arr is None or getattr(arr, "dtype", None) is not wp.vec3d:
                continue
            ptr = int(getattr(arr, "ptr", 0) or 0)
            key = ptr if ptr else id(arr)
            if key not in seen:
                seen.add(key)
                found.append(arr)
    return tuple(found)


def _force_array(runtime: object) -> wp.array | None:
    """Return a runtime's own vec3d force array, or ``None`` when it owns none."""
    for attr in ("force_d", "force", "node_force_d"):
        arr = getattr(runtime, attr, None)
        if arr is not None and getattr(arr, "dtype", None) is wp.vec3d:
            return arr
    return None


def measure_interface_residual(
    world: Any,
    architecture: CellArchitecture,
    *,
    device: str,
    component_runtime: Callable[[str], object] | None = None,
    connector_bound: Callable[[str], bool] | None = None,
    isolate: Callable[[tuple[str, ...]], None] | None = None,
) -> InterfaceResidualReport:
    """Score every attributable cut's force balance on the current candidate state.

    Call AFTER a candidate force assembly (``ComposedNativeCell.solve_candidate``) and BEFORE any
    commit.  Reads only; updates no position and commits nothing.

    Args:
        world: the composed world, used for its ``actor`` when ``component_runtime`` is omitted.
        architecture: the declaration the cuts come from.
        device: resolved CUDA device string.
        component_runtime: optional lookup ``name -> runtime``; defaults to the actor's binding.
        connector_bound: optional predicate ``connector name -> is a runtime bound for it``; defaults
            to the actor's connector binding. A cut with none is reported ``unwired``, because two
            arrays with nothing crossing between them cannot fail to balance.
        isolate: optional ``(connector names) -> None`` that leaves ONLY those connectors'
            contributions in the watched force arrays before the cut is scored — zero the arrays and
            re-run just that edge.  **Without it a reading cannot be attributed to a connector.** A
            whole-body reduction sums every edge on both bodies, so once a body sits on two edges the
            other edge's force appears as an uncancelled residual; that is exactly how a correctly
            wired ``sf_cortex_transient`` read ``cancellation = 0.99999999998`` on 2026-08-11.
            LEAVE-ONE-IN, not leave-one-out — the repository adopted that distinction on 2026-07-28e
            after leave-one-out inherited the whole assembly's round-off.

    Returns:
        An :class:`InterfaceResidualReport`, one entry per declared cut including the unbound ones —
        a cut that could not be measured is reported as such rather than omitted, because omission
        would make an incomplete world look balanced.
    """
    actor = getattr(world, "actor", None)
    if connector_bound is None:
        def connector_bound(name: str) -> bool:  # type: ignore[misc]
            if actor is None:
                return False
            try:
                actor.connector_runtime(name)
            except KeyError:
                return False
            return True

    if component_runtime is None:

        def component_runtime(name: str) -> object | None:  # type: ignore[misc]
            try:
                return actor.component_runtime(name)  # type: ignore[union-attr]
            except (KeyError, AttributeError):
                return None

    readings: list[CutResidual] = []
    for (comp_a, comp_b), connectors in architecture.balance_cuts().items():
        bound = tuple(n for n in connectors if connector_bound(n))
        isolated = False
        if isolate is not None and bound:
            isolate(bound)
            isolated = True
        runtime_a, runtime_b = component_runtime(comp_a), component_runtime(comp_b)
        force_a = _force_array(runtime_a) if runtime_a is not None else None
        force_b = _force_array(runtime_b) if runtime_b is not None else None
        if force_a is None or force_b is None:
            readings.append(
                CutResidual((comp_a, comp_b), connectors, "UNBOUND", None, None, None, bound)
            )
            continue

        # A FRESH ledger per cut: the two channels must hold this cut's two halves and nothing else.
        ledger = make_global_cell_ledger(device=device)
        ledger.add_body_force(force_a, side="reaction")
        ledger.add_body_force(force_b, side="traction")

        n_terms = int(force_a.shape[0]) + int(force_b.shape[0])
        ratio = assemble_balance_tolerance_ratio(n_terms)
        gamma_n_d = wp.array([ratio], dtype=wp.float64, device=device)
        tol_sq_d = wp.zeros(1, dtype=wp.float64, device=device)
        ledger.derive_balance_tolerance(gamma_n_d, tol_sq_d)
        ledger.assemble_balance(tol_sq_d)

        wp.synchronize_device(wp.get_device(device))
        residual = float(ledger.balance_residual_sq_d.numpy()[0]) ** 0.5
        reaction = ledger.reaction_resultant_d.numpy()[0]
        traction = ledger.traction_resultant_d.numpy()[0]
        reaction_pn = float(sum(v * v for v in reaction) ** 0.5)
        traction_pn = float(sum(v * v for v in traction) ** 0.5)
        scale = reaction_pn + traction_pn
        balanced = bool(ledger.balance_ok_d.numpy()[0])
        readings.append(
            CutResidual(
                (comp_a, comp_b), connectors,
                "BALANCED" if balanced else "IMBALANCED",
                residual, scale, ratio, bound,
                reaction_pn=reaction_pn, traction_pn=traction_pn, isolated=isolated,
            )
        )
    return InterfaceResidualReport(tuple(readings))

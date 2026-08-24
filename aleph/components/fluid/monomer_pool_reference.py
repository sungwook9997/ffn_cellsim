"""Multi-consumer G-actin conserved-pool REFERENCE (pure NumPy, no Warp) — whole-cell slice-3.

Host-side acceptance oracle for the ``CHEMICAL_FLUX`` conserved-pool convention of
``WHOLE_CELL_COMMON_CONTRACTS_SPEC_2026-07-23.md`` §4: the cytosolic **G-actin** monomer field is owned by
ONE component (Cytosol / :class:`~aleph.components.fluid.transport.MonomerField`); cortex, SF, lamellipodium, and
filopodium each **draw and return** monomer through an EXPLICIT per-consumer chemical-flux connector — never
a shared array.  This oracle implements the exact discrete accounting the Warp-CUDA
``monomer_flux_commit_kernel`` (``ac/engine/monomer_flux.py``) ports to the device, so the HARD invariant a
closed form cannot express —

    A_total  =  integral phi c dV  +  sum_over_consumers  bound_c        (conserved to round-off)

— is gated on the dev Mac across an arbitrary sequence of per-consumer draws/returns AND the interior RAD
transport, with the additional **separability** guarantee that a consumer can only move monomer between the
field and ITS OWN bound pool (no consumer can draw against another consumer's pool — the "no shared array"
property, verified structurally by the per-consumer reaction split).

This is the DISCRETE companion to ``transport_reference.RADTransportReference`` (which carries the single
``bound_monomer`` polymer pool): here the single bound scalar is generalised to one bound pool per consumer,
and the field's net reaction is the SUM of the per-consumer reactions, so total actin is invariant regardless
of how the individual consumers exchange.  It is NOT a runtime and never a CPU simulation fallback (I0-A).

Sanity Gate (self-tested in tests/ac/fluid/test_monomer_pool_reference.py):
  * total actin conserved to round-off across N consumers, arbitrary draws/returns + transport;
  * separability: consumer c's bound pool changes by exactly ``-integral R_c dV dt`` and nothing else;
  * positivity: a non-negative field with bounded draw stays non-negative under the CFL;
  * balanced treadmill (barbed draw == pointed return, per consumer) leaves field and every pool unchanged;
  * a draw against an unknown consumer is rejected (no implicit shared pool).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from aleph.components.fluid.transport_reference import RADTransportReference

__all__ = ["MonomerPoolConservationReference"]


@dataclass
class MonomerPoolConservationReference:
    """Conservative multi-consumer G-actin pool oracle over the RAD transport stencil.

    The Cytosol owns ONE monomer density field ``u = phi c`` (transported by :class:`RADTransportReference`);
    each named consumer owns its OWN scalar bound pool.  Every accepted flux moves monomer between the field
    and exactly one consumer's pool via that consumer's private reaction field, so ``A_total`` is invariant.

    Attributes:
        shape: Grid shape.
        dx: Uniform spacing.
        phi: Porosity (uniform; ``u = phi c``).
        d_c: Monomer diffusivity (I0-B1c draft; a GAP for the native magnitude, not for this discrete gate).
        consumers: The disjoint consumer names owning bound pools (cortex/sf/lamellipodium/filopodium/...).
        mask: FLUID/OUTSIDE/NUCLEUS classification (defaults all-fluid); non-fluid faces are no-flux.
    """

    shape: tuple[int, ...]
    dx: float
    phi: float
    d_c: float
    consumers: Sequence[str]
    mask: npt.NDArray[np.int_] | None = None
    _bound: dict[str, float] = field(default_factory=dict, init=False, repr=False)
    _transport: RADTransportReference = field(init=False, repr=False)

    def __post_init__(self) -> None:
        names = tuple(self.consumers)
        if len(set(names)) != len(names):
            raise ValueError("consumer names must be unique (each owns a disjoint bound pool)")
        if not names:
            raise ValueError("need at least one monomer consumer")
        self.consumers = names
        self._bound = {name: 0.0 for name in names}
        # The transport oracle carries the SUMMED reaction; its single bound_monomer mirror is the sum of the
        # per-consumer pools (checked in assert_conserved), so it never becomes an independent shared pool.
        self._transport = RADTransportReference(
            shape=self.shape, dx=self.dx, phi=self.phi, d_c=self.d_c, mask=self.mask
        )

    @property
    def bound(self) -> Mapping[str, float]:
        """Read-only view of each consumer's polymer-bound monomer count."""
        return dict(self._bound)

    def field_content(self, u: npt.NDArray[np.float64]) -> float:
        """Field monomer content ``sum_i u_i V_i`` over the fluid cells."""
        return self._transport.field_content(u)

    def total_actin(self, u: npt.NDArray[np.float64]) -> float:
        """Total actin ``integral phi c dV + sum_c bound_c`` — the conserved quantity."""
        return self.field_content(u) + float(sum(self._bound.values()))

    def _consumer_reaction(
        self,
        u: npt.NDArray[np.float64],
        reactions: Mapping[str, npt.NDArray[np.float64] | float] | None,
    ) -> npt.NDArray[np.float64]:
        """Sum the per-consumer reaction fields into the field's net source (rejecting unknown consumers)."""
        total = np.zeros_like(u, dtype=np.float64)
        if reactions is None:
            return total
        for name, r_c in reactions.items():
            if name not in self._bound:
                raise KeyError(
                    f"unknown monomer consumer {name!r}: a flux must name one of {tuple(self._bound)} "
                    "(no implicit shared pool)"
                )
            total = total + np.broadcast_to(np.asarray(r_c, dtype=np.float64), u.shape)
        return total

    def step(
        self,
        u: npt.NDArray[np.float64],
        dt: float,
        per_consumer_reactions: Mapping[str, npt.NDArray[np.float64] | float] | None = None,
        v_f: npt.NDArray[np.float64] | None = None,
    ) -> npt.NDArray[np.float64]:
        """One conservative RAD step; move each consumer's flux to/from ITS OWN pool; return the new ``u``.

        Args:
            u: Current density field ``phi c``.
            dt: Time step (<= the RAD CFL).
            per_consumer_reactions: ``{consumer -> R_c}`` net field source for that consumer (release positive,
                consumption negative).  Each ``R_c`` integral is transferred to/from ``bound[consumer]`` alone,
                so the pools stay separable.  An unnamed consumer is rejected.
            v_f: Cell-centred pore-fluid velocity ``(*grid, dim)``; ``None`` = quiescent.

        Returns:
            The updated density field.
        """
        net = self._consumer_reaction(u, per_consumer_reactions)
        # Transport the field once with the SUMMED reaction (its bound mirror tracks the aggregate exchange).
        u_new = self._transport.step(u, dt, v_f=v_f, reaction=net)
        # Split the exact exchange back to each consumer's OWN pool: bound_c -= integral R_c dV dt.
        if per_consumer_reactions is not None:
            fluid = self._transport.mask == self._transport._FLUID
            cell_vol = self.dx**self._transport.dim
            for name, r_c in per_consumer_reactions.items():
                r = np.broadcast_to(np.asarray(r_c, dtype=np.float64), u.shape)
                self._bound[name] -= float(np.sum(r[fluid]) * cell_vol * dt)
        return u_new

    def assert_conserved(self, u: npt.NDArray[np.float64], a_total0: float, atol: float = 1e-9) -> None:
        """Assert ``A_total`` matches its reference value and the pool sum mirrors the transport oracle."""
        if abs(self.total_actin(u) - a_total0) > atol:
            raise AssertionError(
                f"A_total drifted: {self.total_actin(u)!r} != {a_total0!r} (round-off exchange violated)"
            )
        pool_sum = float(sum(self._bound.values()))
        if abs(pool_sum - self._transport.bound_monomer) > atol:
            raise AssertionError(
                "per-consumer pool sum diverged from the transport oracle's aggregate bound mirror — a "
                "consumer moved monomer outside its own pool (shared-array leak)"
            )

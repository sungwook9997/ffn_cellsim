"""Membrane area response: an effective rest area with a SWEPT reservoir capacity.

PI decision **D2, option B** (2026-07-28).

**The defect this closes.** The ``ac/`` membrane has no area stiffness at all. Its in-plane law
(``ff/membrane_surface.membrane_area_kernel``) takes a scalar tension, and the production path passes the
constant plateau ``gamma_mem`` because the ``K_A`` upturn is default-OFF ("hard-truth #8"). So
``dsigma/dA = 0`` **everywhere** — the infinite-reservoir limit, and the exact opposite of the cell
framework's "the membrane is very stiff in area".

Note this is stronger than the execution plan's own description, which said the membrane sits at a
constant-tension plateau *below strain 0.60*. In ``ac/`` there is no 0.60: ``RESERVOIR_STRAIN`` lives in
``ff/`` and the ``ac/`` path never reaches it. The unsourced magic number was never the real problem —
the missing stiffness was.

**What option B does.** The membrane holds an excess-area reservoir (folds, microvilli, caveolae). Rather
than resolving that population mechanistically — option A, the rule-consistent choice, which needs new
state and an implicit solve — the reservoir is represented by an effective rest area that follows the
apparent area until the reservoir is exhausted::

    A0_eff(t) = min(A(t), A0 * (1 + capacity))        # clamped below by A0
    sigma     = min(gamma_mem + K_A * (A - A0_eff) / A0_eff, tau_lysis)

Below capacity the sheet pays no elastic cost (``A0_eff = A`` so the strain is zero) and the plateau is
reproduced exactly; above it the reservoir is spent and ``K_A`` engages. This IS a lumped stand-in, and
the PI ratified it as one with eyes open.

**Why lumping is defensible here, and only here.** Under the 2026-07-25 reframe the project's purpose is
to INFER molecular parameters, so ``capacity`` stops being an input to choose and becomes a per-cell-type,
per-cell-state signature to recover — the weakness becomes the product. That argument does not
generalise: it works because capacity is an observable of the cell, not a mechanism being replaced by a
proxy. It would NOT license lumping, say, a motor stall force.

**``capacity`` has no default, deliberately.** Supplying one would recreate ``RESERVOIR_STRAIN = 0.60`` —
an unsourced constant that silently governs the physics — under a new name. It is a sweep axis and an
inference target; the caller states it, the artifact records it, and the gates below never assume it.

Sanity Gate:
    * dimensional: ``capacity`` and strain are dimensionless; ``gamma_mem``, ``K_A``, ``tau_lysis`` and
      the returned tension are all [pN/um]; areas are [um^2].
    * boundary: ``capacity -> inf`` recovers today's pure plateau exactly, so this is a strict
      GENERALISATION of the incumbent, not a different model. ``capacity -> 0`` recovers pure ``K_A``
      elasticity measured from ``A0``. Both limits are gated.
    * continuity: sigma is continuous at ``A = A0*(1 + capacity)``, where it equals ``gamma_mem`` from
      both sides — a jump there would be a spurious force.
    * monotonicity: sigma is non-decreasing in area, and STRICTLY increasing above capacity. That strict
      increase is the area stiffness that is currently identically zero.
    * conservation/sign: below ``A0`` the sheet is slack, not in compression — the law floors at
      ``gamma_mem`` rather than returning a negative tension a bilayer cannot carry.
    * not-tuned: no constant in this module is chosen to make a gate pass; ``capacity`` is required from
      the caller and every gate is a shape or a limit, never a magnitude.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class MembraneAreaCard:
    """Parameters of the effective-rest-area law, with ``capacity`` as a declared sweep axis.

    Attributes:
        gamma_mem_pn_per_um: The plateau tension [pN/um] while reservoir remains (KB-3.B1.1).
        k_area_pn_per_um: Areal stretch modulus ``K_A`` [pN/um] once the reservoir is spent (Rawicz).
        tau_lysis_pn_per_um: Lysis tension [pN/um]; the law saturates here rather than growing without
            bound, because a real bilayer ruptures.
        capacity: Excess-area fraction the reservoir can supply before ``K_A`` engages. **Required, and
            an OUTPUT of inference rather than an input to choose** — see the module docstring.
        capacity_provenance: How this capacity value was obtained (a sweep point, a posterior, a cited
            measurement). Required, so an artifact can never carry a bare number.
    """

    gamma_mem_pn_per_um: float
    k_area_pn_per_um: float
    tau_lysis_pn_per_um: float
    capacity: float
    capacity_provenance: str

    def __post_init__(self) -> None:
        for name in ("gamma_mem_pn_per_um", "k_area_pn_per_um", "tau_lysis_pn_per_um"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive; got {value!r}")
        if self.tau_lysis_pn_per_um <= self.gamma_mem_pn_per_um:
            raise ValueError(
                "tau_lysis must exceed the resting plateau, else the membrane starts already lysed"
            )
        if math.isnan(self.capacity) or self.capacity < 0.0:
            raise ValueError(
                f"capacity is an excess-area fraction and cannot be negative; got {self.capacity!r}. "
                "math.inf is allowed and reproduces the incumbent infinite-reservoir plateau"
            )
        if not self.capacity_provenance.strip():
            raise ValueError(
                "capacity needs a recorded provenance — an unexplained value is RESERVOIR_STRAIN=0.60 "
                "under a new name"
            )


def effective_rest_area(area_um2: float, rest_area_um2: float, capacity: float) -> float:
    """Return ``A0_eff`` [um^2]: the rest area after the reservoir has yielded what it can.

    Args:
        area_um2: Current apparent membrane area ``A`` [um^2].
        rest_area_um2: Unstretched reference area ``A0`` [um^2].
        capacity: Excess-area fraction available. ``math.inf`` means an unlimited reservoir.

    Returns:
        ``min(A, A0*(1 + capacity))``, floored at ``A0`` — the reservoir can pay out area but cannot
        absorb it, so a sub-rest membrane is slack rather than pre-compressed.

    Raises:
        ValueError: If either area is non-finite or non-positive, or ``capacity`` is negative.
    """
    if not math.isfinite(area_um2) or area_um2 <= 0.0:
        raise ValueError(f"area must be finite and positive; got {area_um2!r}")
    if not math.isfinite(rest_area_um2) or rest_area_um2 <= 0.0:
        raise ValueError(f"rest area must be finite and positive; got {rest_area_um2!r}")
    if math.isnan(capacity) or capacity < 0.0:
        raise ValueError(f"capacity cannot be negative; got {capacity!r}")
    ceiling = math.inf if math.isinf(capacity) else rest_area_um2 * (1.0 + capacity)
    return max(rest_area_um2, min(area_um2, ceiling))


def area_tension(area_um2: float, rest_area_um2: float, card: MembraneAreaCard) -> float:
    """Return the in-plane tension ``sigma`` [pN/um] to hand to ``membrane_area_kernel``.

    Args:
        area_um2: Current apparent membrane area ``A`` [um^2].
        rest_area_um2: Unstretched reference area ``A0`` [um^2].
        card: The parameter card, carrying the declared ``capacity``.

    Returns:
        ``gamma_mem`` while reservoir remains, then ``gamma_mem + K_A * (A - A0_eff)/A0_eff``, saturating
        at ``tau_lysis``.

    Raises:
        ValueError: If the areas are not finite and positive.

    Note:
        This returns a tension, not a force. It replaces the constant scalar the production path passes
        today; the force law itself (``ff/membrane_surface.membrane_area_kernel``) is unchanged, which is
        why this is additive and cannot perturb any existing result while the incumbent keeps passing a
        constant.
    """
    a0_eff = effective_rest_area(area_um2, rest_area_um2, card.capacity)
    if math.isinf(a0_eff):  # pragma: no cover - unreachable: a0_eff is bounded by a finite area
        return card.gamma_mem_pn_per_um
    strain = (area_um2 - a0_eff) / a0_eff
    if strain <= 0.0:
        return card.gamma_mem_pn_per_um
    return min(
        card.gamma_mem_pn_per_um + card.k_area_pn_per_um * strain,
        card.tau_lysis_pn_per_um,
    )


def area_stiffness(area_um2: float, rest_area_um2: float, card: MembraneAreaCard) -> float:
    """Return ``dsigma/dA`` [pN/um^3] — the quantity that is identically zero in the incumbent.

    Computed analytically rather than by finite difference so a gate can assert it is EXACTLY zero
    inside the reservoir and strictly positive outside, with no differencing tolerance to argue about.

    Args:
        area_um2: Current apparent membrane area ``A`` [um^2].
        rest_area_um2: Unstretched reference area ``A0`` [um^2].
        card: The parameter card.

    Returns:
        ``0`` while reservoir remains or once lysis tension saturates the law; ``K_A / A0_eff`` in the
        elastic band between them.
    """
    a0_eff = effective_rest_area(area_um2, rest_area_um2, card.capacity)
    strain = (area_um2 - a0_eff) / a0_eff
    if strain <= 0.0:
        return 0.0
    if card.gamma_mem_pn_per_um + card.k_area_pn_per_um * strain >= card.tau_lysis_pn_per_um:
        return 0.0                       # saturated at lysis; the law is flat again
    return card.k_area_pn_per_um / a0_eff


__all__ = ["MembraneAreaCard", "area_stiffness", "area_tension", "effective_rest_area"]

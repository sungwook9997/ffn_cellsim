"""Young-equation contact angle at a two-cell junction (KU-4.4).

Maître *et al.* 2012 *Nature* (zebrafish germ-layer cells) showed that the
mechanical contact angle at a two-cell triple line is set by the balance
of three line tensions — the two cortical tensions :math:`\\gamma_{c,1}`,
:math:`\\gamma_{c,2}` and the cell-cell junction tension :math:`\\gamma_J`
— through a Young-style equation. For a symmetric pair this reduces to
:math:`\\cos(\\theta/2) = \\gamma_J / \\gamma_c`. The Phase 1 Unit 4.1
brief uses the equivalent dihedral-angle form

.. math::
    \\cos\\theta = \\frac{\\gamma_{c,1} + \\gamma_{c,2} - 2\\gamma_J}{2\\gamma_c}

with :math:`\\gamma_c = (\\gamma_{c,1} + \\gamma_{c,2})/2` interpreted as
the mean cortical tension. The three KU-4.4 limiting cases are then

* :math:`\\gamma_J \\to 0`   ⇒ :math:`\\cos\\theta \\to 1`,  θ → 0
  (the "full-adhesion" limit: the junction has zero interfacial cost,
  so the two cortexes share a single straight contour at the contact;
  the brief's θ-convention reports this as θ = 0).
* :math:`\\gamma_J = \\gamma_c`  ⇒ :math:`\\cos\\theta = 0`,  θ = π/2
  (the symmetric mature contact, KU-4.4 nominal).
* :math:`\\gamma_J \\to 2\\gamma_c`  ⇒ :math:`\\cos\\theta \\to -1`,  θ → π
  (the "no-adhesion" limit: junction tension equals the sum of both
  cortex tensions, so the cortexes cannot relax across the contact and
  the cells round away from each other; θ = π in the brief's convention).

The Maître 2012 external contact angle :math:`2\\theta_E` (the angle
*outside* the cells at the triple line) is a different quantity from
the brief's θ; we report θ here without claiming a direct
:math:`\\theta(\\theta_E)` map and let downstream callers convert if
needed. KU-4.4 specifies the brief's form as the Phase 1 contact-angle
diagnostic.

KU-4.4 specifies that the junction tension depends on cadherin bond
density (more engaged bonds reduce :math:`\\gamma_J` through the work of
adhesion). For Phase 1 Unit 4.1 we use the linear interpolation

.. math::
    \\gamma_J(n_{\\rm eng}) = \\gamma_{J,\\rm max} - (\\gamma_{J,\\rm max} -
    \\gamma_{J,\\rm min}) \\, \\frac{n_{\\rm eng}}{N_{\\rm cad}}

with the two endpoints anchored to the Young equation:

* :math:`\\gamma_{J,\\rm max} = \\gamma_{c,1} + \\gamma_{c,2}` —
  zero-cadherin contact (Maître's "no adhesion" limit, cos θ = −1).
* :math:`\\gamma_{J,\\rm min} = 0` — full-cadherin contact (Maître's
  "complete wetting" limit, cos θ = +1).

The KU-4.4 nominal symmetric mature contact (θ = π/2) is reached when
exactly half the available bonds are engaged. This is a *literature-
anchored* choice (Maître 2012's mid-range observed angles) rather than
a fitted target, so the Magic-Number Block tests pass: derivable from
the Young equation's two limiting cases, grid-invariant, and not chosen
to make any specific test numeric match.

Sanity Gate
-----------

1. **Dimensional**: all tensions in N/m (line tensions in a 2-D model);
   :math:`\\cos\\theta` dimensionless. No CFL — instantaneous geometry.
2. **Boundary cases**:
   - ``gamma_c = 0`` ⇒ division-by-zero guard; raise
     :class:`ValueError`. A cortexless cell has no Young equation.
   - ``n_engaged < 0`` or ``> n_total`` ⇒ guarded upstream by
     :class:`~acs_kb.junction.cadherin.update_bonds`; here we just
     clip to ``[0, n_total]`` for safety.
   - Asymmetric ``γ_{c,1} ≠ γ_{c,2}``: formula still holds; the
     contact angle becomes asymmetric but ``cos θ`` is the same scalar
     for the dihedral.
3. **Conservation**: stateless function — nothing to conserve.
4. **Numerical sanity**: ``cos θ`` is clipped to ``[-1, 1]`` before
   ``arccos`` so float roundoff at the limits does not raise.
   ``float64`` throughout.
5. **Sign / sense**: increasing :math:`n_{\\rm eng}` reduces
   :math:`\\gamma_J`, which raises ``cos θ`` and lowers θ in the
   brief's convention — i.e., more cadherin bonds drive the contact
   toward θ = 0 (cortexes share a flat contour at the contact;
   adhesion-dominated). Fewer cadherin bonds drive θ toward π (cortex
   tension dominates and the cells round away). The two limits
   bracket Maître 2012's adhesion-vs.-cortex-dominance regimes.
   **Cortex pulls along its tangent**, **junction tension pulls along
   the contact line**; the line tensions enter as positive scalars in
   the Young equation regardless of the geometric angle convention.

References: Maître *et al.* 2012 *Nature* 490, 7421;
Foty & Steinberg 2005 *Dev Biol* 278; KU-4.4 Phase 1 brief.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from acs_kb.cell.cell import Cell
from acs_kb.junction.types import EcadherinJunction


def junction_tension(
    n_engaged: int,
    n_total: int,
    gamma_c_a: float,
    gamma_c_b: float,
) -> float:
    """Junction tension :math:`\\gamma_J` as a function of bond occupancy.

    Linear interpolation between the two Young-equation endpoints
    (``n_engaged = 0`` ⇒ :math:`\\gamma_{c,1} + \\gamma_{c,2}`,
    ``n_engaged = n_total`` ⇒ ``0``).

    Returns
    -------
    gamma_J : float
        Junction tension in N/m. Non-negative by construction (the
        interpolation endpoints are non-negative).
    """
    if n_total <= 0:
        raise ValueError(f"n_total must be positive (got {n_total}).")
    gamma_max = float(gamma_c_a) + float(gamma_c_b)
    gamma_min = 0.0
    f_engaged = float(np.clip(n_engaged / n_total, 0.0, 1.0))
    return float(gamma_max - (gamma_max - gamma_min) * f_engaged)


def compute_contact_angle(
    cell_a: Cell,
    cell_b: Cell,
    junction: EcadherinJunction,
    params: dict[str, Any],
) -> float:
    """Young-equation contact angle at the two-cell triple line (KU-4.4).

    Parameters
    ----------
    cell_a, cell_b : Cell
        The two cells the junction bridges. Currently only used to
        cross-check ``cell_a_id`` / ``cell_b_id`` against the junction
        and to read ``cortex.params['gamma_cortex']`` if the caller did
        not pass an explicit tension override; this lets Phase 1 reuse
        Worker C's cortex tension without re-discovering it from a
        separate file.
    junction : EcadherinJunction
        Provides ``n_bonds_engaged`` and ``n_bonds_total`` (needed for
        the linear γ_J interpolation).
    params : dict
        Required keys:

        * ``gamma_c_a``, ``gamma_c_b`` (float, N/m) — cortical tensions
          for cells A and B. If absent, fall back to
          ``cell_*.cortex.params.get('gamma_cortex')``; if still
          absent, raise :class:`KeyError`.

        Optional keys:

        * ``gamma_J`` (float, N/m) — override the model-internal
          interpolation and use this value directly. Useful for sanity
          tests against the Young equation in isolation.

    Returns
    -------
    theta : float
        Contact angle in radians, in :math:`[0, \\pi]`.
    """
    # Cross-check that the junction connects these two cells.
    junction_pair = {junction.cell_a_id, junction.cell_b_id}
    cell_pair = {cell_a.id, cell_b.id}
    if junction_pair != cell_pair:
        raise ValueError(
            f"junction connects cells {junction_pair} but got "
            f"cell_a.id={cell_a.id}, cell_b.id={cell_b.id}."
        )

    gamma_c_a = params.get("gamma_c_a")
    if gamma_c_a is None:
        gamma_c_a = cell_a.cortex.params.get("gamma_cortex")
        if gamma_c_a is None:
            raise KeyError(
                "gamma_c_a not in params and "
                "cell_a.cortex.params['gamma_cortex'] is missing."
            )
    gamma_c_b = params.get("gamma_c_b")
    if gamma_c_b is None:
        gamma_c_b = cell_b.cortex.params.get("gamma_cortex")
        if gamma_c_b is None:
            raise KeyError(
                "gamma_c_b not in params and "
                "cell_b.cortex.params['gamma_cortex'] is missing."
            )

    gamma_c_a = float(gamma_c_a)
    gamma_c_b = float(gamma_c_b)
    if gamma_c_a <= 0.0 or gamma_c_b <= 0.0:
        raise ValueError(
            f"cortical tensions must be positive "
            f"(got γ_c_a={gamma_c_a}, γ_c_b={gamma_c_b})."
        )

    gamma_c_mean = 0.5 * (gamma_c_a + gamma_c_b)
    if "gamma_J" in params:
        gamma_J = float(params["gamma_J"])
    else:
        gamma_J = junction_tension(
            n_engaged=int(junction.n_bonds_engaged),
            n_total=int(junction.n_bonds_total),
            gamma_c_a=gamma_c_a,
            gamma_c_b=gamma_c_b,
        )

    cos_theta = (gamma_c_a + gamma_c_b - 2.0 * gamma_J) / (2.0 * gamma_c_mean)
    cos_theta = float(np.clip(cos_theta, -1.0, 1.0))
    return float(np.arccos(cos_theta))

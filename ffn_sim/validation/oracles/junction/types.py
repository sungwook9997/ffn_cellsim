"""Junction value objects (week-9 interface contract).

The :class:`EcadherinJunction` dataclass is the shared data contract
between Worker D's Unit 4.1 (this unit, two-cell pair) and the upcoming
Unit 4.2 (multicellular aggregate) plus any Phase 2 Junction Track worker.
Field names and types here are **frozen at week 9**; any change must be
coordinated with downstream consumers via Notion alert before merge.

All quantities are SI (metres, seconds, Newtons, radians). The 2-D
contact plane is assumed for Phase 1; the 3-D extension is reserved for
Phase 2 along with catch-bond detail (KU-4.2 follow-up).

The junction is identified by the integer ``id`` pair of the two
``<v1 archived Cell duck type>`` instances it bridges. We deliberately store
``id`` rather than a Cell reference so the junction object is hashable,
serialisable, and free of import-time cycles with the cell module.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from typing import Any as Cell  # v1 Cortex-coupled Cell stripped; oracle uses duck typing


@dataclass(slots=True)
class EcadherinJunction:
    """Single E-cadherin junction state (Phase 1 Unit 4.1).

    Attributes
    ----------
    cell_a_id, cell_b_id : int
        :attr:`<v1 archived Cell duck type>.id` of the two cells this junction
        bridges. ``cell_a_id < cell_b_id`` is the canonical ordering used
        by the factory; downstream code should not rely on it but tests
        do so for deterministic iteration.
    contact_position_a, contact_position_b : np.ndarray, shape (2,)
        Contact-point coordinates (m) on cell A's and cell B's cortex
        boundary respectively. Sampled at construction time via
        :meth:`Cell.compute_cortex_boundary_position` and refreshed each
        time the junction is re-paired (Unit 4.2). For a symmetric pair
        on the +x / −x axis, ``contact_position_a`` is at cell A's
        boundary facing cell B (positive x for left-of-pair cell).
    n_bonds_total : int
        Total cadherin pairs available at this contact (engaged + free
        reservoir). KU-4.17 mature-junction default ``N_cad_per_contact
        = 100``.
    n_bonds_engaged : int
        Currently bonded cadherin pairs. Updated stochastically by
        :func:`ffn_sim.validation.oracles.junction.cadherin (ARCHIVED v1).update_bonds`.
    bond_forces : np.ndarray, shape ``(n_bonds_engaged,)``
        Per-bond axial tensile force in Newtons. Length follows
        :attr:`n_bonds_engaged` exactly so iteration is
        ``for f in junction.bond_forces:`` without masking. Phase 1
        loads all engaged bonds equally with ``F_total / n_engaged``
        (mean-field per KU-4.2 simplification); per-bond variance is a
        Phase 2 extension.
    contact_angle : float
        Cortex-cortex contact angle θ at the junction triple line in
        radians, measured by the KU-4.4 Young equation
        (:func:`ffn_sim.validation.oracles.junction.contact_angle.compute_contact_angle`).
        Convention: θ = 0 → cortices parallel, no spreading; θ = π/2 →
        symmetric mature contact (acceptance); θ → π → full wetting.
    age : float
        Seconds since junction instantiation. Used by KU-4.11 junction
        maturation (Phase 1 informational; Phase 2 may gate the
        slip-to-catch switch on it).
    """

    cell_a_id: int
    cell_b_id: int
    contact_position_a: np.ndarray
    contact_position_b: np.ndarray
    n_bonds_total: int
    n_bonds_engaged: int
    bond_forces: np.ndarray
    contact_angle: float
    age: float = 0.0


def make_ecadherin_junction(
    cell_a: Cell,
    cell_b: Cell,
    *,
    n_bonds_total: int = 100,
    contact_angle_init: float = 0.0,
) -> EcadherinJunction:
    """Construct a fresh junction between two cells (KU-4.17 defaults).

    The contact points are sampled from each cell's cortex along the
    centre-to-centre direction: cell A's contact lies on its boundary in
    the +(b−a) direction, cell B's on its boundary in the −(b−a)
    direction, both reported in lab coordinates by
    :meth:`Cell.compute_cortex_boundary_position`.

    All cadherin pairs start disengaged
    (``n_bonds_engaged = 0``, ``bond_forces`` is a length-0 ``float64``
    array); the first call to
    :func:`ffn_sim.validation.oracles.junction.cadherin (ARCHIVED v1).update_bonds` engages bonds according
    to KU-4.17's :math:`k_{\\rm on}`.

    Parameters
    ----------
    cell_a, cell_b : Cell
        Two distinct Worker C cells (``cell_a.id != cell_b.id``). The
        canonical ordering ``cell_a_id < cell_b_id`` is enforced inside
        the junction to keep junctions hashable by an ordered pair.
    n_bonds_total : int, default 100
        KU-4.17 ``N_cad_per_contact`` for a mature junction.
    contact_angle_init : float, default 0.0
        Initial contact angle in radians. Real value is computed by the
        Young-equation routine; 0.0 is the geometric pre-spreading state.
    """
    if cell_a.id == cell_b.id:
        raise ValueError("cell_a and cell_b must have distinct ids.")
    if n_bonds_total <= 0:
        raise ValueError(f"n_bonds_total must be positive (got {n_bonds_total}).")

    # Canonicalise: cell_a_id < cell_b_id for stable iteration.
    if cell_a.id > cell_b.id:
        cell_a, cell_b = cell_b, cell_a

    delta = cell_b.center_position - cell_a.center_position
    angle_ab = float(np.arctan2(delta[1], delta[0]))  # a → b direction
    angle_ba = angle_ab + np.pi                       # b → a direction

    contact_a = cell_a.compute_cortex_boundary_position(angle_ab)
    contact_b = cell_b.compute_cortex_boundary_position(angle_ba)

    return EcadherinJunction(
        cell_a_id=int(cell_a.id),
        cell_b_id=int(cell_b.id),
        contact_position_a=np.asarray(contact_a, dtype=np.float64).reshape(2),
        contact_position_b=np.asarray(contact_b, dtype=np.float64).reshape(2),
        n_bonds_total=int(n_bonds_total),
        n_bonds_engaged=0,
        bond_forces=np.zeros(0, dtype=np.float64),
        contact_angle=float(contact_angle_init),
        age=0.0,
    )

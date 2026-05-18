"""ECM adapter — swap-in replacement for the Phase 1 substrate stub.

Unit 2.1's :class:`acs_kb.bridge.substrate_stub.LinearElasticSubstrate`
returned ``k_sub = π · E_eff · a`` from a Boussinesq half-space. This
Unit replaces it with an adapter that consumes Worker A's discrete
fibre network (``acs_kb.ecm.fiber_network.FiberNetwork``,
``acs_kb.ecm.fiber_mechanics.compute_forces``) and reports a *linearised*
local stiffness at the FA position.

Method (linear response, KU-1.21 linear regime)
-----------------------------------------------
The network sits at its rest configuration (compute_forces ≈ 0 modulo
thermal noise). A small finite-difference probe ``δ`` applied to the
nearest bead in the +x direction gives ``ΔF`` along x; the diagonal
Hessian element is::

    k_local = − (F_x(r + δ x̂) − F_x(r)) / δ                  [N/m]

evaluated by calling Worker A's ``compute_forces`` twice. We
symmetrise over (+δ, −δ) for a centred difference, and average the x-
and y-direction probes for isotropic networks. The probe magnitude
``δ`` is one bond-fraction (default 1 % of ``rest_length``), well
inside the linear-response regime.

This makes the adapter API-identical to the substrate stub: it exposes
``stiffness`` (the cached scalar k_local) and ``compute_displacement
(position, force) -> force / stiffness``. The motor-clutch stepper
treats both interchangeably.

Caveats / deferred to Phase 2
-----------------------------
- Single-bead Hessian, no off-diagonal coupling. Real Green's function
  ``G(r, r')`` away from the FA is non-trivial.
- Periodic-image fibre intersections are not enumerated (matches
  Unit 1.1 cross-link generator's limitation).
- Network drift due to the applied force is not fed back into the
  fibre dynamics — the substrate is held at rest. Unit 2.x simulation
  durations are short (≤ 60 s) so this is consistent with Worker A's
  Phase 1 "energy and forces only" scope.

Sanity Gate
-----------
1. Dimensional analysis: ``δ`` [m], ``ΔF`` [N] ⇒ ``k_local`` [N/m] ✓.
2. Boundary cases: ``compute_forces`` returns 0 at rest ⇒ probe ratio
   is well-defined. If the nearest bead is isolated (no bonds, no
   cross-links) the diagonal Hessian is 0 — guarded by raising
   ``ECMAdapterError`` rather than silently returning ∞ compliance.
3. Conservation: probe is symmetric (+δ, −δ); each compute_forces call
   satisfies Newton 3rd law internally (Worker A's invariant).
4. Numerical sanity: float64; two compute_forces evaluations cached at
   construction.
5. Sign / sense: tensile force at the bead → displacement in the same
   direction (positive k_local).
6. Measurement protocol: ``stiffness`` is computed once at __init__ and
   cached. The consistency integration test compares it against a
   stub's k_sub at a matched E and asserts the *traction* under
   identical motor-clutch input agrees within ±10 % (the Brief Task 4
   acceptance).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from acs_kb.ecm.fiber_mechanics import compute_forces
from acs_kb.ecm.fiber_network import FiberNetwork


class ECMAdapterError(RuntimeError):
    """Raised when the local Hessian cannot be evaluated (isolated bead)."""


class SubstrateProtocol(Protocol):
    """Shared interface between the stub and the ECM adapter (Unit 2.2)."""

    @property
    def stiffness(self) -> float: ...
    def compute_displacement(
        self, position: np.ndarray, force: float | np.ndarray
    ) -> float | np.ndarray: ...


@dataclass(frozen=True, slots=True)
class ECMAdapterParams:
    """Configuration for the ECM adapter (see ``configs/phase1_unit2_2.yaml``)."""

    stretching_modulus: float        # μ, KU-1.2
    bending_modulus: float           # κ, KU-1.1
    local_radius: float = 5.0e-6     # 5 μm probe radius (not yet used in Phase 1)
    probe_fraction: float = 0.01     # δ = probe_fraction · rest_length

    @classmethod
    def from_config(cls, cfg: dict) -> "ECMAdapterParams":
        b = cfg["bridge"] if "bridge" in cfg else cfg
        a = b["ecm_adapter"]
        return cls(
            stretching_modulus=float(a["stretching_modulus"]),
            bending_modulus=float(a["bending_modulus"]),
            local_radius=float(a.get("local_radius", 5.0e-6)),
            probe_fraction=float(a.get("probe_fraction", 0.01)),
        )


class ECMAdapter:
    """ECM-backed substrate (drop-in for :class:`LinearElasticSubstrate`).

    Parameters
    ----------
    network : FiberNetwork
        Worker A's discrete fibre network at rest. Held by reference;
        not mutated by the adapter.
    fa_position : np.ndarray, shape (2,)
        FA centre in the lab frame. The adapter finds the nearest bead
        and uses its diagonal Hessian as the local effective stiffness.
    params : ECMAdapterParams
        Stretching / bending moduli and the FD probe size.
    cross_links : list or None
        Optional cross-link list (Worker A KU-1.28) — if provided, the
        Hessian includes the cross-link contribution at the probed bead.
    """

    __slots__ = (
        "network", "fa_position", "params", "cross_links",
        "_nearest_fiber", "_nearest_bead", "_stiffness",
    )

    def __init__(
        self,
        network: FiberNetwork,
        fa_position: np.ndarray,
        params: ECMAdapterParams,
        cross_links: list | None = None,
    ) -> None:
        self.network = network
        self.fa_position = np.asarray(fa_position, dtype=np.float64).reshape(2)
        self.params = params
        self.cross_links = cross_links if cross_links is not None else []
        self._nearest_fiber, self._nearest_bead = self._find_nearest_bead()
        self._stiffness = self._compute_local_stiffness()

    # ------------------------------------------------------------------ #
    # Construction helpers                                               #
    # ------------------------------------------------------------------ #

    def _find_nearest_bead(self) -> tuple[int, int]:
        beads = self.network.bead_positions
        delta = beads - self.fa_position
        delta -= self.network.box_size * np.round(delta / self.network.box_size)
        d2 = (delta * delta).sum(axis=-1)
        f, b = np.unravel_index(int(np.argmin(d2)), d2.shape)
        return int(f), int(b)

    def _compute_local_stiffness(self) -> float:
        """Centred-difference diagonal Hessian at the nearest bead.

        ``k = − (F(+δ) − F(−δ)) / (2δ)``  averaged over +x and +y probes.
        Worker A's ``compute_forces`` is called four times (one per
        signed probe direction); both calls preserve Newton 3rd law.
        """
        net = self.network
        p = self.params
        delta = p.probe_fraction * net.rest_length
        if delta <= 0.0:
            raise ECMAdapterError("probe_fraction · rest_length must be > 0")

        f_idx, b_idx = self._nearest_fiber, self._nearest_bead
        k_vals: list[float] = []
        for axis in (0, 1):
            pos_plus = net.bead_positions.copy()
            pos_plus[f_idx, b_idx, axis] += delta
            pos_minus = net.bead_positions.copy()
            pos_minus[f_idx, b_idx, axis] -= delta
            F_plus = compute_forces(
                pos_plus, net.rest_length, p.stretching_modulus,
                p.bending_modulus, net.box_size, cross_links=self.cross_links,
            )
            F_minus = compute_forces(
                pos_minus, net.rest_length, p.stretching_modulus,
                p.bending_modulus, net.box_size, cross_links=self.cross_links,
            )
            dF = F_plus[f_idx, b_idx, axis] - F_minus[f_idx, b_idx, axis]
            k_axis = -dF / (2.0 * delta)
            k_vals.append(float(k_axis))
        k_eff = 0.5 * (k_vals[0] + k_vals[1])
        if k_eff <= 0.0:
            raise ECMAdapterError(
                f"Non-positive effective stiffness at bead (f={f_idx}, b={b_idx}): "
                f"k_x={k_vals[0]:.3e}, k_y={k_vals[1]:.3e}. "
                f"Network may be isolated or unstably configured at this bead."
            )
        return k_eff

    # ------------------------------------------------------------------ #
    # SubstrateProtocol surface                                          #
    # ------------------------------------------------------------------ #

    @property
    def stiffness(self) -> float:
        return self._stiffness

    def compute_displacement(
        self, position: np.ndarray, force: float | np.ndarray
    ) -> float | np.ndarray:
        """``u = F / k_local``. ``position`` is API parity only."""
        del position
        return force / self._stiffness

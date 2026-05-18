"""Linear-elastic substrate stub (KU-1.21, Phase 1 Unit 2.1).

This is a deliberately minimal substrate model so Worker B can run the
motor-clutch dynamics before Worker A's full discrete ECM is ready. The
stub treats the substrate as a single-DOF linear spring at every clutch
contact, computed from the Boussinesq point-force solution for an
infinite linear-elastic half-space (KU-1.21 linear regime, Storm-MacKintosh
small-strain) collapsed to a per-clutch compliance::

    u = F / k_sub          [m]
    k_sub = π · E_eff · a  [N/m]
    E_eff = E / (1 − ν²)   [Pa]

with the **substrate-compliance length scale** ``a`` (Phase 1 default
1 μm per Brief Task 2) playing the role of an effective point-contact
radius in the Boussinesq form. This is a *linear*, isotropic,
half-space approximation: no Hertz contact, no finite thickness
correction, no inter-clutch elastic coupling. The Brief calls out that
those are intentionally deferred — Hertz contact validation sits with
Worker A's KU-1.21 module and need not be re-derived here.

``a`` (the substrate-stub compliance length) is **distinct** from the
FocalAdhesion's ``area`` field (the FA's projected footprint on the
substrate, KU-2.18 default 1 μm²). The two are *independent* model
parameters in this stub: ``a`` sets the clutch-to-substrate spring
constant, while ``area`` is an FA state variable used (from Unit 2.2
onward) for Hill-function FA growth (KU-2.17) and traction stress
reporting. A geometric disc-radius derivation ``a = √(A/π) ≈ 0.564 μm``
is also defensible but is *not* what Phase 1 uses — switching to it
would shift ``k_sub`` by 1/0.564 = 1.77× and re-position the biphasic
peak. Unit 2.2's ECM adapter eliminates this ambiguity by computing
``k_sub`` directly from the local fibre network.

Sanity Gate
-----------
1. Dimensional analysis: ``E`` in Pa, ``a`` in m → ``k_sub`` in N/m;
   ``F`` in N → ``u`` in m. No CFL (this module is stateless, called by
   the motor-clutch step).
2. Boundary cases: ``E → 0`` → ``k_sub → 0`` → diverging displacement,
   handled by the caller's substrate-stiffness sweep guard (E > 0). ``F
   = 0`` returns ``u = 0`` exactly. ν = 0.5 → E_eff = ∞·E/3, still
   finite for ν < 0.5 (Phase 1 uses 0.45).
3. Conservation: no energy stored — instantaneous compliance only.
4. Numerical sanity: float64; closed-form algebra, no iteration.
5. Sign / sense: positive ``F`` (tensile, into the substrate from the
   clutch) produces positive ``u`` (substrate displaces toward actin).
6. Measurement protocol: ``stiffness`` and ``compute_displacement`` are
   the only outputs the motor-clutch module reads; both are functions
   of (E, ν, a) only. No state, no hidden defaults.

Replacement plan
----------------
Phase 1 Unit 2.2+ will swap this for ``acs_kb.bridge.ecm_adapter`` that
projects an ECM fibre network into a local effective stiffness sampled
at the FA position. The motor-clutch module reads only ``stiffness`` and
``compute_displacement``, so the swap is API-compatible.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class LinearElasticSubstrate:
    """Linear elastic half-space stub (KU-1.21).

    The four attributes are **independent model parameters**, not
    geometric projections of each other:

    Attributes
    ----------
    young_modulus : float
        Bulk Young's modulus ``E`` in pascals. Phase 1 default 5 kPa
        (Engler PAA mid-range, KU-1.5).
    poisson_ratio : float
        Poisson's ratio ν. Phase 1 default 0.45 (PAA, KU-1.21).
    thickness : float
        Substrate thickness in metres. Carried for downstream
        traceability; not used in the linear half-space stub.
    contact_radius : float
        Phenomenological **substrate-compliance length scale** ``a`` in
        metres, used by the Boussinesq form ``k_sub = π · E_eff · a``.
        This is an *effective* point-contact radius that sets the
        clutch-to-substrate spring constant — it is **not** a geometric
        projection of the FA's projected area.

        For Phase 1 we follow the Brief Task 2 specification of
        ``a = 1 μm``. The FocalAdhesion's ``area`` field (Phase 1
        default 1 μm², KU-2.18) is a separate quantity describing the
        FA's footprint on the substrate; the two are independent model
        parameters in this stub (a geometric disc-radius derivation
        ``a = √(A/π) ≈ 0.564 μm`` is also defensible but is *not* what
        Phase 1 uses, see ``docstring of substrate_stub`` and the
        Unit 2.1 REPORT.md for rationale).

        Unit 2.2's ECM adapter will compute the local stiffness from
        the discrete fibre network and obviate ``contact_radius``
        altogether.
    """

    young_modulus: float = 5.0e3
    poisson_ratio: float = 0.45
    thickness: float = 100.0e-6
    contact_radius: float = 1.0e-6

    @classmethod
    def from_config(cls, cfg: dict) -> "LinearElasticSubstrate":
        """Build from a resolved ``bridge.substrate`` config block."""
        s = cfg["bridge"]["substrate"] if "bridge" in cfg else cfg
        return cls(
            young_modulus=float(s["young_modulus"]),
            poisson_ratio=float(s["poisson_ratio"]),
            thickness=float(s["thickness"]),
            contact_radius=float(s["contact_radius"]),
        )

    @property
    def effective_modulus(self) -> float:
        """``E_eff = E / (1 − ν²)``  [Pa]."""
        return self.young_modulus / (1.0 - self.poisson_ratio * self.poisson_ratio)

    @property
    def stiffness(self) -> float:
        """Per-clutch substrate stiffness ``k_sub = π · E_eff · a``  [N/m]."""
        return math.pi * self.effective_modulus * self.contact_radius

    def compute_displacement(
        self, position: np.ndarray, force: float | np.ndarray
    ) -> float | np.ndarray:
        """Linear point compliance ``u = F / k_sub``.

        Parameters
        ----------
        position : np.ndarray, shape (2,)
            Lab-frame point where ``F`` is applied. Linear stub does not
            use it; carried so the API matches the ECM adapter (Unit 2.2+).
        force : float or np.ndarray
            Axial force at ``position``, in Newtons. Positive = tensile
            (clutch pulls substrate toward actin).

        Returns
        -------
        displacement : same type as ``force``
            Substrate displacement at the loading point, in metres.
        """
        del position  # signature parity with the future ECM adapter
        k = self.stiffness
        return force / k

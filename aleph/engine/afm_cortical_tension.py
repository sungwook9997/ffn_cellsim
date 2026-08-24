"""Pre-declared AFM cortical-tension protocol for a suspended round cell.

This module is the experiment contract, not an indentation solver.  It prevents the old modulus harness from
being reused under a cortical-tension label when the measurement has changed in four load-bearing ways:

* Young--Laplace ``gamma = deltaP * R / 2`` is admitted only for a suspended round cell.  An adherent spread
  cell is an SF/FA apparent-stiffness experiment and is rejected at construction.
* the holder is the T10 exterior Stokes medium that loads all six rigid modes.  A focal-adhesion holder is not
  an alternative value of the same experiment.
* the AFM indenter is an experiment-only 15th component joined to the membrane by a sixth CONTACT edge.  The
  reference cell remains 14 components / 38 connectors; the extension is scoped to the apparatus and the edge
  is serviced by :mod:`aleph.engine.connector_joints`, hence by the same gap-gated unilateral kernel as the
  five biological contacts.
* no force curve is enumerable until Pi_0 provenance, contact distance/stiffness calibration, the C-2 physical
  acceptance predicate, a physical load speed and a three-seed ensemble are all present.

No magnitude is defaulted here.  In particular, the HeLa 40 Pa value is representable only as a proxy unless
the exact target protocol supplies it directly; a proxy can drive a mechanism demonstration but cannot be
promoted to a quantitative cortical-tension sweep.  Contact distance and penalty stiffness are likewise
required inputs with provenance rather than constants hidden in the connector.

Sanity Gate (host structural tests in ``test_afm_cortical_tension.py``):
  * dimensions: Pi_0 [Pa], speed [um/s], indentation depth/radius/contact distance [um], contact stiffness
    [pN/um], exterior viscosity [Pa*s]; every numeric value must be finite and positive (depth may include 0).
  * boundary cases: adherent geometry, an FA holder, fewer than three distinct seeds, a zero/negative apparatus
    value, missing provenance, or an unconverged/rejected C-2 record all refuse sweep enumeration.
  * conservation/sign: the apparatus edge is CONTACT, bidirectional and adjoint, and
    :func:`force_kernel_for` selects the existing unilateral law; the protocol introduces no second drag owner.
  * measurement: reaction channel 0 is declared before a run; every sweep point is a seed replicate and carries
    the same immutable apparatus/provenance block.  Biological-parameter axes are deliberately absent -- stage
    4 owns those after Fisher identifiability and stage-2 time convergence.

CPU-importable: this file allocates and launches nothing.  CUDA execution belongs to a later driver after the
blockers returned by :meth:`AFMSweepDeclaration.blockers` are empty.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum
from itertools import product

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentContract,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
    reference_cell_architecture,
)

__all__ = [
    "AFMContactCalibration",
    "AFMSweepDeclaration",
    "AFMSweepPoint",
    "C2Acceptance",
    "CellGeometry",
    "Holder",
    "OsmoticSetpoint",
    "Pi0Evidence",
    "ResultClass",
    "SphericalIndenterProtocol",
    "afm_experiment_architecture",
]

AFM_INDENTER_COMPONENT = "afm_indenter"
AFM_MEMBRANE_CONTACT = "indenter_membrane_contact"
AFM_REACTION_CHANNEL = 0


class CellGeometry(StrEnum):
    """Geometry/state identity of the measurement."""

    SUSPENDED_ROUND = "suspended_round"
    ADHERENT_SPREAD = "adherent_spread"


class Holder(StrEnum):
    """What removes whole-cell rigid motion during the measurement."""

    EXTERIOR_STOKES = "exterior_stokes_six_rigid_modes"
    FOCAL_ADHESION = "focal_adhesion"


class Pi0Evidence(StrEnum):
    """How the osmotic setpoint is related to the exact target protocol."""

    DIRECT_TARGET_PROTOCOL = "direct_target_protocol"
    CROSS_PROTOCOL_PROXY = "cross_protocol_proxy"
    DERIVED_ESTIMATE = "derived_estimate"


class ResultClass(StrEnum):
    """The strongest use the declared inputs license before execution."""

    QUANTITATIVE = "quantitative"
    MECHANISM_DEMO = "mechanism_demo_not_quantitative"


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive; got {value!r}")
    return value


def _source(value: str, name: str) -> str:
    value = str(value).strip()
    if not value:
        raise ValueError(f"{name} provenance/source is required")
    return value


@dataclass(frozen=True, slots=True)
class OsmoticSetpoint:
    """Explicit Pi_0 plus the evidence that licenses its use."""

    value_pa: float
    evidence: Pi0Evidence
    source: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "value_pa", _positive(self.value_pa, "Pi_0 [Pa]"))
        if not isinstance(self.evidence, Pi0Evidence):
            raise TypeError("Pi_0 evidence must be a Pi0Evidence member")
        object.__setattr__(self, "source", _source(self.source, "Pi_0"))

    @property
    def result_class(self) -> ResultClass:
        """A proxy/derived Pi_0 is a mechanism demonstration, never quantitative."""
        return (
            ResultClass.QUANTITATIVE
            if self.evidence is Pi0Evidence.DIRECT_TARGET_PROTOCOL
            else ResultClass.MECHANISM_DEMO
        )


@dataclass(frozen=True, slots=True)
class AFMContactCalibration:
    """Numerical contact distance/stiffness declaration for the unilateral connector."""

    contact_distance_um: float
    stiffness_pn_per_um: float
    provenance: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "contact_distance_um", _positive(self.contact_distance_um, "contact distance [um]")
        )
        object.__setattr__(
            self, "stiffness_pn_per_um", _positive(self.stiffness_pn_per_um, "contact stiffness [pN/um]")
        )
        object.__setattr__(self, "provenance", _source(self.provenance, "contact calibration"))


@dataclass(frozen=True, slots=True)
class SphericalIndenterProtocol:
    """Apparatus values fixed by the reproduced AFM protocol."""

    target_protocol: str
    target_cell_line: str
    geometry: CellGeometry
    holder: Holder
    indenter_radius_um: float
    exterior_viscosity_pa_s: float
    apparatus_source: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "target_protocol", _source(self.target_protocol, "target protocol"))
        object.__setattr__(self, "target_cell_line", _source(self.target_cell_line, "target cell line"))
        if self.geometry is not CellGeometry.SUSPENDED_ROUND:
            raise ValueError(
                "cortical-tension AFM requires a suspended round cell; adherent spread cells measure "
                "SF/FA-dominated apparent stiffness"
            )
        if self.holder is not Holder.EXTERIOR_STOKES:
            raise ValueError(
                "a cortical-tension protocol holds rigid modes with the exterior Stokes medium, not focal "
                "adhesions"
            )
        object.__setattr__(
            self, "indenter_radius_um", _positive(self.indenter_radius_um, "indenter radius [um]")
        )
        object.__setattr__(
            self,
            "exterior_viscosity_pa_s",
            _positive(self.exterior_viscosity_pa_s, "exterior viscosity [Pa*s]"),
        )
        object.__setattr__(self, "apparatus_source", _source(self.apparatus_source, "apparatus"))


@dataclass(frozen=True, slots=True)
class C2Acceptance:
    """Full-native physical-step evidence required before a force curve can be read."""

    build_commit: str
    inner_converged: bool
    outer_accepted: bool
    committed_time_s: float
    max_projected_force_pn: float
    projected_force_tolerance_pn: float
    membrane_subdivisions: int
    full_native_population: bool
    full_compartments: bool
    accepted_seeds: tuple[int, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "build_commit", _source(self.build_commit, "C-2 build commit"))
        committed = float(self.committed_time_s)
        if not math.isfinite(committed) or committed < 0.0:
            raise ValueError("C-2 committed physical time must be finite and nonnegative")
        object.__setattr__(self, "committed_time_s", committed)
        max_force = float(self.max_projected_force_pn)
        tolerance = _positive(self.projected_force_tolerance_pn, "C-2 projected-force tolerance [pN]")
        if not math.isfinite(max_force) or max_force < 0.0:
            raise ValueError("C-2 max projected force must be finite and nonnegative")
        object.__setattr__(self, "max_projected_force_pn", max_force)
        object.__setattr__(self, "projected_force_tolerance_pn", tolerance)
        subdivisions = int(self.membrane_subdivisions)
        if subdivisions < 0:
            raise ValueError("C-2 membrane subdivisions must be nonnegative")
        object.__setattr__(self, "membrane_subdivisions", subdivisions)
        seeds = tuple(int(seed) for seed in self.accepted_seeds)
        if any(seed < 0 for seed in seeds) or len(seeds) != len(set(seeds)):
            raise ValueError("C-2 accepted seeds must be distinct nonnegative integers")
        object.__setattr__(self, "accepted_seeds", seeds)

    @property
    def passed(self) -> bool:
        """Whether the inner candidate became an accepted physical step."""
        return (
            self.inner_converged
            and self.outer_accepted
            and self.committed_time_s > 0.0
            and self.max_projected_force_pn <= self.projected_force_tolerance_pn
            and self.membrane_subdivisions >= 8
            and self.full_native_population
            and self.full_compartments
            and len(self.accepted_seeds) >= 3
        )


@dataclass(frozen=True, slots=True)
class AFMSweepPoint:
    """One protocol-validation replicate; never a biological-parameter grid point."""

    speed_um_s: float
    depth_um: float
    seed: int


@dataclass(frozen=True, slots=True)
class AFMSweepDeclaration:
    """All prerequisites and apparatus axes for a rate/depth/seed validation sweep."""

    protocol: SphericalIndenterProtocol
    pi0: OsmoticSetpoint | None
    contact: AFMContactCalibration | None
    c2: C2Acceptance | None
    speeds_um_s: tuple[float, ...]
    depths_um: tuple[float, ...]
    seeds: tuple[int, ...]

    def __post_init__(self) -> None:
        speeds = tuple(_positive(v, "loading speed [um/s]") for v in self.speeds_um_s)
        depths = tuple(float(v) for v in self.depths_um)
        if not speeds:
            raise ValueError("at least one physical loading speed is required")
        if not depths or any(not math.isfinite(v) or v < 0.0 for v in depths):
            raise ValueError("indentation depths must be a non-empty tuple of finite nonnegative values")
        seeds = tuple(int(v) for v in self.seeds)
        if len(seeds) < 3 or len(set(seeds)) != len(seeds):
            raise ValueError("cortical tension requires at least three distinct seed replicates")
        if any(v < 0 for v in seeds):
            raise ValueError("RNG seeds must be nonnegative")
        object.__setattr__(self, "speeds_um_s", speeds)
        object.__setattr__(self, "depths_um", depths)
        object.__setattr__(self, "seeds", seeds)

    @property
    def reaction_channel(self) -> int:
        """The pre-declared cantilever readout: vertical indenter reaction ``reaction[0]``."""
        return AFM_REACTION_CHANNEL

    def blockers(self, *, require_quantitative: bool = True) -> tuple[str, ...]:
        """Return why CUDA force-curve enumeration is not yet licensed."""
        blockers: list[str] = []
        if require_quantitative:
            blockers.append(
                "ACTIVE_CORTEX_STATE_MISSING: no post-induction, equilibrated, dt-converged active-cortex "
                "state handoff exists; STATE (c) 3/17 therefore still block cortical-tension magnitudes"
            )
        if self.pi0 is None:
            blockers.append("PI0_MISSING: the target protocol has no osmotic setpoint")
        elif require_quantitative and self.pi0.result_class is not ResultClass.QUANTITATIVE:
            blockers.append(
                "PI0_PROXY: this Pi_0 licenses a mechanism demonstration, not a quantitative result"
            )
        if self.contact is None:
            blockers.append("CONTACT_CALIBRATION_MISSING: rest_um and stiffness are undeclared")
        if self.c2 is None:
            blockers.append("C2_MISSING: no physical-step acceptance record was supplied")
        elif not self.c2.passed:
            blockers.append(
                "C2_FAILED: require accepted full-native/full-compartment mechanics, maxPF<=the unchanged "
                "predicate, membrane subdivision>=8, committed time, and >=3 accepted seeds"
            )
        return tuple(blockers)

    def points(self, *, require_quantitative: bool = True) -> tuple[AFMSweepPoint, ...]:
        """Enumerate rate/depth/seed points only after every preflight blocker is closed."""
        blockers = self.blockers(require_quantitative=require_quantitative)
        if blockers:
            raise RuntimeError("AFM sweep is blocked: " + "; ".join(blockers))
        return tuple(
            AFMSweepPoint(speed_um_s=speed, depth_um=depth, seed=seed)
            for speed, depth, seed in product(self.speeds_um_s, self.depths_um, self.seeds)
        )


def afm_experiment_architecture(
    base: CellArchitecture | None = None,
) -> CellArchitecture:
    """Extend the reference cell with the apparatus component and sixth unilateral CONTACT edge.

    The returned object is experiment-scoped.  :func:`reference_cell_architecture` is not mutated, so its
    canonical census remains 14/38.  The spherical runtime declares one fixed pair per membrane vertex against
    the single sphere centre.  Prescribed centre motion changes gaps but never pair identity, so this edge does
    not request topology generation or accepted-step remapping; the unilateral kernel gates separated pairs.
    """
    architecture = base or reference_cell_architecture()
    if any(component.name == AFM_INDENTER_COMPONENT for component in architecture.components):
        raise ValueError(f"architecture already declares {AFM_INDENTER_COMPONENT!r}")
    if any(connector.name == AFM_MEMBRANE_CONTACT for connector in architecture.connectors):
        raise ValueError(f"architecture already declares {AFM_MEMBRANE_CONTACT!r}")
    indenter = ComponentContract(
        AFM_INDENTER_COMPONENT,
        ComponentRole.ENVIRONMENT,
        "prescribed rigid spherical AFM indenter surface",
        "prescribed-motion rigid-body apparatus",
        owns_geometry=True,
        dynamically_evolving=True,
    )
    contact = ConnectorContract(
        AFM_MEMBRANE_CONTACT,
        ConnectorFamily.CONTACT,
        AFM_INDENTER_COMPONENT,
        "membrane",
        False,
        False,
        endpoint_role_a="rigid indenter surface material point",
        endpoint_role_b="membrane contact quadrature",
        generation_required=False,
        remap_on_accept=False,
    )
    return CellArchitecture(
        components=architecture.components + (indenter,),
        connectors=architecture.connectors + (contact,),
    )

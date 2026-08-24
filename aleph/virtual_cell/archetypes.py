"""Ring A/B/C cell-archetype manifests, priors as distributions, and missing-law refusals.

The programme's destination is molecular parameters *per cell type and cell state*, so the archetype
registry is where the layer either widens past one cell line or quietly re-narrows to it.  This module
is built to prevent the second.  It names no cell line as the target; the two adherent lines the repo
has measured are at most one regression benchmark inside a panel, and nothing here privileges them.

Two rings, two opposite jobs.

**Ring A** are archetypes the canonical actomyosin grammar should already carry: an immune-like
amoeboid cell, an endothelial cell, a fibroblast/myofibroblast, a junctional epithelial cell, a
stem/iPSC-like cell, and a pulsatile contractile cell.  Each emits a real
:class:`~aleph.virtual_cell.contracts.CellStateManifest` whose components and connectors are the
INDUCED SUBGRAPH of the census read from :mod:`aleph.virtual_cell.grammar` -- so an archetype cannot
name a component the engine does not declare, and cannot leave a declared component dangling with no
connector.

A finding worth stating plainly, because it is the honest result of doing this: **at the component
level the canonical 14-component census barely distinguishes these six.**  Four of the six are the
same component set or one component away from it.  Everything that actually separates an endothelial
cell from a fibroblast from an epithelial cell lives in the ``missing_laws`` block, not in the
component census -- junctional force transmission, apicobasal polarity, ECM deposition, shear sensing.
The component census is archetype-coarse, and the applicability block is where the archetype is.

**Ring B/C** are the opposite test: an erythrocyte, a neuron/growth cone, a ciliated cell, and a
cardiomyocyte.  These are NOT forced into the grammar.  Each REFUSES to produce a manifest and emits a
:class:`MissingLawVerdict` naming the specific constitutive laws the actomyosin-centric component set
cannot express, each justified by quoting the representation string the engine itself declares (a
``fluid Helfrich surface FEM`` has no in-plane shear modulus; an ``active rod/cable graph`` has no
sarcomeric register).  A successful Ring B/C entry is a refusal with a named law.

**Priors are distributions, and almost all of them are ``PI_GAP``.**  A prior is only ``SOURCED`` when
a source can be named; this module cites no literature, because citing a paper this lane has not read
would be worse than the gap.  The repo has already been burned by exactly the alternative: the live
resting turgor is a HeLa measurement standing in for a different line, with a ~3.3x spread across the
three values in the tree (``aleph/components/incumbent/assemble.py``:135-147).  A ``PI_GAP`` with a stated reason
is the useful output; an invented lognormal is not.

Nothing here runs physics, imports Warp, or grants an evidence rung.

Sanity Gate:
  * dimensional: every prior carries an explicit ``unit`` string and no prior carries a bare number.
    A ``PI_GAP`` prior carries no parameters at all, so no unit-less magnitude can escape it.
  * boundary: a ``SOURCED`` prior with no source, a ``PI_GAP`` prior with parameters, an unsupported
    Ring A archetype, and a Ring B/C entry that produced a manifest are each rejected at construction.
  * conservation-invariant: the topological invariant is delegated -- every Ring A manifest is built
    from the induced subgraph of the canonical census and is asserted to have no isolated component,
    so no archetype can rely on co-location instead of a connector.
  * measurement-protocol: an archetype declares no measurement.  ``unsupported_claims`` names the
    observables that MAY NOT be quoted from a run of this archetype, and it is populated from the
    census (e.g. the exterior medium is declared-only, so no migration speed is quotable).
  * numerical: prior parameters are finite floats or the prior is a gap; there is no default value,
    no fallback, and no code path that substitutes a number for a missing one.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from types import MappingProxyType

from aleph.virtual_cell.contracts import CellStateManifest
from aleph.virtual_cell.grammar import EngineGrammar, load_engine_grammar

__all__ = [
    "ApplicabilityBlock",
    "ArchetypeBuild",
    "ArchetypeRing",
    "DistributionFamily",
    "MissingLaw",
    "MissingLawVerdict",
    "PriorDistribution",
    "PriorStatus",
    "RING_A_ARCHETYPES",
    "RING_BC_ARCHETYPES",
    "ARCHETYPE_IDS",
    "ARCHETYPE_EXPANSION",
    "build_all_archetypes",
    "build_archetype",
]

#: Where an archetype failure routes (shared lane brief, "failure routing").
ARCHETYPE_EXPANSION = "archetype failure -> missing component / missing law project"

#: Why no archetype-specific prior in this module carries a number.
_PI_GAP_PREAMBLE = (
    "no archetype-specific measurement is registered for this quantity and this lane cites no "
    "literature it has not read; the repo's own worst case is the resting turgor, where a HeLa "
    "measurement stands in for a different line across a ~3.3x spread "
    "(aleph/components/incumbent/assemble.py:135-147)"
)


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _unique(values: Sequence[str], *, what: str) -> tuple[str, ...]:
    normalized = tuple(_nonempty(value, what=f"{what} entry") for value in values)
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{what} must not contain duplicates")
    return normalized


class ArchetypeRing(StrEnum):
    """Which ring an archetype belongs to, and therefore what a success looks like."""

    RING_A = "ring-a"
    RING_B = "ring-b"
    RING_C = "ring-c"


class PriorStatus(StrEnum):
    """Whether a prior carries a sourced distribution or is an open PI gap."""

    SOURCED = "sourced"
    PI_GAP = "pi-gap"


class DistributionFamily(StrEnum):
    """The distribution family a sourced prior uses; ``UNRESOLVED`` belongs to a gap."""

    UNRESOLVED = "unresolved"
    NORMAL = "normal"
    LOGNORMAL = "lognormal"
    UNIFORM = "uniform"
    BETA = "beta"


@dataclass(frozen=True, slots=True)
class PriorDistribution:
    """One population or geometry prior, as a DISTRIBUTION or as an explicit gap.

    A point value is never representable here: a sourced prior must name a family and its parameters,
    and an unsourced prior must be a :attr:`PriorStatus.PI_GAP` carrying no number at all.
    """

    prior_id: str
    quantity: str
    unit: str
    status: PriorStatus
    family: DistributionFamily = DistributionFamily.UNRESOLVED
    parameters: Mapping[str, float] = MappingProxyType({})
    source: str | None = None
    pi_gap_reason: str | None = None

    def __post_init__(self) -> None:
        for name in ("prior_id", "quantity", "unit"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        if not isinstance(self.status, PriorStatus):
            raise TypeError("status must be a PriorStatus")
        if not isinstance(self.family, DistributionFamily):
            raise TypeError("family must be a DistributionFamily")
        if not isinstance(self.parameters, Mapping):
            raise TypeError("parameters must be a mapping")
        frozen: dict[str, float] = {}
        for key, value in self.parameters.items():
            name = _nonempty(key, what="parameter name")
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError(f"parameter {name!r} must be a real number")
            number = float(value)
            if number != number or number in (float("inf"), float("-inf")):
                raise ValueError(f"parameter {name!r} must be finite")
            frozen[name] = number
        object.__setattr__(self, "parameters", MappingProxyType(frozen))

        if self.status is PriorStatus.PI_GAP:
            if frozen:
                raise ValueError(
                    f"prior {self.prior_id!r} is a PI gap and must carry no parameters; a gap with "
                    "a number is an invented value"
                )
            if self.family is not DistributionFamily.UNRESOLVED:
                raise ValueError(f"prior {self.prior_id!r} is a PI gap and has no family")
            if self.source is not None:
                raise ValueError(f"prior {self.prior_id!r} is a PI gap and cannot name a source")
            object.__setattr__(
                self, "pi_gap_reason", _nonempty(self.pi_gap_reason, what="pi_gap_reason")
            )
        else:
            if not frozen:
                raise ValueError(f"sourced prior {self.prior_id!r} needs distribution parameters")
            if self.family is DistributionFamily.UNRESOLVED:
                raise ValueError(f"sourced prior {self.prior_id!r} needs a distribution family")
            object.__setattr__(self, "source", _nonempty(self.source, what="source"))
            if self.pi_gap_reason is not None:
                raise ValueError(f"sourced prior {self.prior_id!r} cannot carry a gap reason")

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-able record of this prior."""
        return {
            "prior_id": self.prior_id,
            "quantity": self.quantity,
            "unit": self.unit,
            "status": self.status.value,
            "family": self.family.value,
            "parameters": dict(sorted(self.parameters.items())),
            "source": self.source,
            "pi_gap_reason": self.pi_gap_reason,
        }


@dataclass(frozen=True, slots=True)
class MissingLaw:
    """One constitutive law the canonical grammar cannot express, and why."""

    law_id: str
    statement: str
    blocked_by: str
    nearest_existing: str | None = None

    def __post_init__(self) -> None:
        for name in ("law_id", "statement", "blocked_by"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        if self.nearest_existing is not None:
            object.__setattr__(
                self, "nearest_existing", _nonempty(self.nearest_existing, what="nearest_existing")
            )

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-able record of this missing law."""
        return {
            "law_id": self.law_id,
            "statement": self.statement,
            "blocked_by": self.blocked_by,
            "nearest_existing": self.nearest_existing,
        }


@dataclass(frozen=True, slots=True)
class ApplicabilityBlock:
    """What an archetype cannot express and what may not be claimed from running it."""

    missing_components: tuple[str, ...]
    missing_laws: tuple[MissingLaw, ...]
    unsupported_claims: tuple[str, ...]
    pi_gap_prior_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "missing_components", _unique(self.missing_components, what="missing_components")
        )
        laws = tuple(self.missing_laws)
        if any(not isinstance(law, MissingLaw) for law in laws):
            raise TypeError("missing_laws must contain MissingLaw values")
        law_ids = tuple(law.law_id for law in laws)
        if len(set(law_ids)) != len(law_ids):
            raise ValueError("missing_laws must not repeat a law_id")
        object.__setattr__(self, "missing_laws", laws)
        object.__setattr__(
            self, "unsupported_claims", _unique(self.unsupported_claims, what="unsupported_claims")
        )
        object.__setattr__(
            self, "pi_gap_prior_ids", _unique(self.pi_gap_prior_ids, what="pi_gap_prior_ids")
        )

    @property
    def law_ids(self) -> tuple[str, ...]:
        """Return the missing-law identifiers in declaration order."""
        return tuple(law.law_id for law in self.missing_laws)

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-able record of this applicability block."""
        return {
            "missing_components": list(self.missing_components),
            "missing_laws": [law.to_dict() for law in self.missing_laws],
            "unsupported_claims": list(self.unsupported_claims),
            "pi_gap_prior_ids": list(self.pi_gap_prior_ids),
        }


@dataclass(frozen=True, slots=True)
class MissingLawVerdict:
    """A Ring B/C refusal: no manifest, and the named laws that would be required to build one."""

    archetype_id: str
    ring: ArchetypeRing
    missing_components: tuple[str, ...]
    missing_laws: tuple[MissingLaw, ...]
    unsupported_claims: tuple[str, ...]
    expansion: str = ARCHETYPE_EXPANSION
    refused: bool = True

    def __post_init__(self) -> None:
        object.__setattr__(self, "archetype_id", _nonempty(self.archetype_id, what="archetype_id"))
        if not isinstance(self.ring, ArchetypeRing):
            raise TypeError("ring must be an ArchetypeRing")
        if self.ring is ArchetypeRing.RING_A:
            raise ValueError("a Ring A archetype does not emit a refusal verdict")
        object.__setattr__(
            self, "missing_components", _unique(self.missing_components, what="missing_components")
        )
        laws = tuple(self.missing_laws)
        if not laws or any(not isinstance(law, MissingLaw) for law in laws):
            raise ValueError(
                f"{self.archetype_id!r} refuses without naming a missing law; a generic failure is "
                "not a verdict"
            )
        law_ids = tuple(law.law_id for law in laws)
        if len(set(law_ids)) != len(law_ids):
            raise ValueError("missing_laws must not repeat a law_id")
        object.__setattr__(self, "missing_laws", laws)
        object.__setattr__(
            self, "unsupported_claims", _unique(self.unsupported_claims, what="unsupported_claims")
        )
        object.__setattr__(self, "expansion", _nonempty(self.expansion, what="expansion"))
        if self.refused is not True:
            raise ValueError("a MissingLawVerdict is a refusal by construction")

    @property
    def law_ids(self) -> tuple[str, ...]:
        """Return the missing-law identifiers in declaration order."""
        return tuple(law.law_id for law in self.missing_laws)

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-able record of this refusal."""
        return {
            "archetype_id": self.archetype_id,
            "ring": self.ring.value,
            "refused": True,
            "missing_components": list(self.missing_components),
            "missing_laws": [law.to_dict() for law in self.missing_laws],
            "unsupported_claims": list(self.unsupported_claims),
            "expansion": self.expansion,
        }


@dataclass(frozen=True, slots=True)
class ArchetypeBuild:
    """The result of asking the canonical grammar for one archetype: a manifest, or a refusal."""

    archetype_id: str
    ring: ArchetypeRing
    supported: bool
    priors: tuple[PriorDistribution, ...]
    applicability: ApplicabilityBlock
    manifest: CellStateManifest | None = None
    verdict: MissingLawVerdict | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "archetype_id", _nonempty(self.archetype_id, what="archetype_id"))
        if not isinstance(self.ring, ArchetypeRing):
            raise TypeError("ring must be an ArchetypeRing")
        if not isinstance(self.supported, bool):
            raise TypeError("supported must be bool")
        priors = tuple(self.priors)
        if any(not isinstance(prior, PriorDistribution) for prior in priors):
            raise TypeError("priors must contain PriorDistribution values")
        prior_ids = tuple(prior.prior_id for prior in priors)
        if len(set(prior_ids)) != len(prior_ids):
            raise ValueError("priors must not repeat a prior_id")
        object.__setattr__(self, "priors", priors)
        if not isinstance(self.applicability, ApplicabilityBlock):
            raise TypeError("applicability must be an ApplicabilityBlock")

        if self.supported:
            if self.manifest is None:
                raise ValueError(f"{self.archetype_id!r} is supported but produced no manifest")
            if self.verdict is not None:
                raise ValueError(f"{self.archetype_id!r} is supported and cannot also refuse")
            if self.ring is not ArchetypeRing.RING_A:
                raise ValueError(
                    f"{self.archetype_id!r} is Ring B/C; it must not be forced into the grammar"
                )
        else:
            if self.manifest is not None:
                raise ValueError(f"{self.archetype_id!r} refused but still produced a manifest")
            if self.verdict is None:
                raise ValueError(
                    f"{self.archetype_id!r} refused without a MissingLawVerdict; a bare failure is "
                    "not an answer"
                )
            if self.verdict.archetype_id != self.archetype_id:
                raise ValueError("verdict archetype_id must match the build")

    @property
    def pi_gap_prior_ids(self) -> tuple[str, ...]:
        """Return the identifiers of every prior that is an open PI gap."""
        return tuple(prior.prior_id for prior in self.priors if prior.status is PriorStatus.PI_GAP)

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-able record of this build."""
        return {
            "archetype_id": self.archetype_id,
            "ring": self.ring.value,
            "supported": self.supported,
            "manifest": self.manifest.to_dict() if self.manifest is not None else None,
            "manifest_hash": self.manifest.manifest_hash if self.manifest is not None else None,
            "priors": [prior.to_dict() for prior in self.priors],
            "applicability": self.applicability.to_dict(),
            "verdict": self.verdict.to_dict() if self.verdict is not None else None,
        }


def _gap(prior_id: str, quantity: str, unit: str, detail: str) -> PriorDistribution:
    """Build one open-gap prior; a gap carries a reason, never a number."""
    return PriorDistribution(
        prior_id=prior_id,
        quantity=quantity,
        unit=unit,
        status=PriorStatus.PI_GAP,
        pi_gap_reason=f"{detail}. {_PI_GAP_PREAMBLE}",
    )


def _population_priors(
    archetype_id: str, extra: Sequence[tuple[str, str, str, str]] = ()
) -> tuple[PriorDistribution, ...]:
    """Return the geometry/population prior slots every archetype must declare.

    Every slot is a DISTRIBUTION slot, and every one of them is currently an open gap.  The slots
    exist so that the gap is enumerable: an unnamed gap cannot be closed, and a manifest with no
    prior slots would look complete.
    """
    base: tuple[tuple[str, str, str, str], ...] = (
        (
            "cell_radius_um",
            "equivalent spherical cell radius",
            "um",
            f"{archetype_id} geometry is not measured in this tree",
        ),
        (
            "cortex_thickness_nm",
            "cortical actin layer thickness",
            "nm",
            f"{archetype_id} cortical thickness is not measured in this tree",
        ),
        (
            "cortical_tension_pn_per_um",
            "resting cortical tension",
            "pN/um",
            f"{archetype_id} cortical tension is not measured in this tree, and it is set by the "
            "resting turgor, which is itself an open PI gap",
        ),
        (
            "cortical_f_actin_density_per_um3",
            "cortical F-actin filament number density",
            "1/um^3",
            f"{archetype_id} cortical filament density is not measured in this tree; the native "
            "population is derived from density x geometry, so this gap gates the population count",
        ),
        (
            "nmii_minifilament_density_per_um3",
            "NMII minifilament number density",
            "1/um^3",
            f"{archetype_id} minifilament density is not measured in this tree",
        ),
        (
            "nucleus_volume_fraction",
            "nuclear volume fraction",
            "1",
            f"{archetype_id} nuclear volume fraction is not measured in this tree",
        ),
        (
            "environment_stiffness_kpa",
            "environment (substrate or matrix) Young modulus",
            "kPa",
            f"the physiological niche stiffness for {archetype_id} is a PI decision, not a "
            "convenience default; starting from a non-physiological baseline invalidates the "
            "comparison rather than weakening it",
        ),
    )
    return tuple(_gap(*slot) for slot in (*base, *extra))


@dataclass(frozen=True, slots=True)
class _RingASpec:
    """Internal declaration of one Ring A archetype."""

    archetype_id: str
    species: str
    lineage: str
    identity: Mapping[str, str]
    biological_state: Mapping[str, str]
    environment: Mapping[str, str]
    excluded_components: frozenset[str]
    missing_components: tuple[str, ...]
    missing_laws: tuple[MissingLaw, ...]
    unsupported_claims: tuple[str, ...]
    extra_prior_slots: tuple[tuple[str, str, str, str], ...] = ()
    observations: tuple[str, ...] = ()


#: Claims no archetype may quote today, because the census says the mechanism is declared only.
_MEDIUM_DECLARED_ONLY = (
    "any whole-cell translocation speed: the exterior medium component is declared with no solve "
    "behind it, so the six rigid-body modes are regularised numerically rather than loaded"
)
_TURGOR_GAP_CLAIM = (
    "any absolute resting cortical tension: it is set by the resting turgor, which is an open gap"
)


_AMOEBOID = _RingASpec(
    archetype_id="amoeboid-immune",
    species="Homo sapiens",
    lineage="haematopoietic (neutrophil / T lymphocyte class)",
    identity={
        "archetype_class": "amoeboid",
        "adhesion_mode": "adhesion-independent / low-affinity integrin",
    },
    biological_state={"activation": "motile", "polarity": "front-back"},
    environment={
        "confinement": "3D confining channel or interstitial space",
        "matrix_engagement": "friction / intercalation rather than focal adhesion",
    },
    excluded_components=frozenset({"focal_adhesion", "ecm", "world_boundary", "sf_arc"}),
    missing_components=("nonspecific_surface_friction",),
    missing_laws=(
        MissingLaw(
            law_id="law.adhesion_independent_friction_traction",
            statement=(
                "a traction law that converts retrograde cortical flow into propulsion through "
                "non-specific surface friction and topographic intercalation, with no ligand-bound "
                "clutch"
            ),
            blocked_by=(
                "every traction path in the census runs through focal_adhesion, whose only "
                "load-bearing edge is the kinetic integrin_collagen_clutch; there is no "
                "ligand-free friction connector at all"
            ),
            nearest_existing="integrin_collagen_clutch",
        ),
        MissingLaw(
            law_id="law.confinement_pressure_boundary",
            statement=(
                "a confining boundary that applies a normal load over the whole cell surface, so a "
                "cell can push against a channel wall"
            ),
            blocked_by=(
                "world_boundary is declared as a 'physiological far-field reference frame' with "
                "owns_geometry=False and dynamically_evolving=False, so it cannot represent a "
                "moving confining wall"
            ),
            nearest_existing="ecm_far_field_anchor",
        ),
        MissingLaw(
            law_id="law.chemotactic_gradient_sensing",
            statement=(
                "a receptor-level gradient sensing law that biases protrusion and cortical "
                "contractility from an external chemoattractant field"
            ),
            blocked_by=(
                "no component owns a signalling field; the census's only fields are the "
                "Biot/Darcy cytosol and the declared-only exterior medium"
            ),
        ),
    ),
    unsupported_claims=(
        _MEDIUM_DECLARED_ONLY,
        _TURGOR_GAP_CLAIM,
        "migration persistence or directionality under a gradient: no gradient-sensing law exists",
    ),
    observations=("cell_shape_silhouette", "cortical_flow_field"),
)

_ENDOTHELIAL = _RingASpec(
    archetype_id="endothelial",
    species="Homo sapiens",
    lineage="vascular endothelium",
    identity={"archetype_class": "endothelial", "junction_type": "VE-cadherin adherens"},
    biological_state={"activation": "quiescent-or-angiogenic", "polarity": "apicobasal"},
    environment={
        "matrix": "basement-membrane-like",
        "mechanical_load": "luminal shear + circumferential stretch",
    },
    excluded_components=frozenset(),
    missing_components=("cell_cell_junction", "glycocalyx", "lumen"),
    missing_laws=(
        MissingLaw(
            law_id="law.ve_cadherin_adherens_junction",
            statement=(
                "a bidirectional cell-cell force-transmitting junction with catch-bond kinetics, "
                "so a monolayer carries tension between cells rather than each cell in isolation"
            ),
            blocked_by=(
                "the canonical census declares NO cell-cell connector of any family; every "
                "connector joins two components of ONE cell, or a cell component to the shared "
                "environment"
            ),
            nearest_existing="integrin_collagen_clutch",
        ),
        MissingLaw(
            law_id="law.glycocalyx_shear_transduction",
            statement=(
                "conversion of luminal wall shear stress into a cortical/junctional mechanical "
                "signal through a compliant surface layer"
            ),
            blocked_by=(
                "the membrane is a 'fluid Helfrich surface FEM' with no surface-attached polymer "
                "layer, and the exterior medium that would carry the shear is declared only"
            ),
            nearest_existing="membrane_medium_traction",
        ),
        MissingLaw(
            law_id="law.apicobasal_polarity_field",
            statement=(
                "an apicobasal polarity field that assigns distinct material and kinetic "
                "properties to apical vs basal membrane and cortex"
            ),
            blocked_by=(
                "membrane and cortex are each ONE homogeneous surface component; the census has no "
                "sub-surface domain partition and no polarity field owner"
            ),
        ),
        MissingLaw(
            law_id="law.barrier_permeability",
            statement="paracellular transport across a junctional barrier as a function of junction tension",
            blocked_by=(
                "the only fluid boundary connectors are membrane_cytosol_boundary and "
                "nucleus_cytosol_boundary, both single-cell; there is no inter-cell flux path"
            ),
        ),
    ),
    unsupported_claims=(
        _MEDIUM_DECLARED_ONLY,
        _TURGOR_GAP_CLAIM,
        "monolayer-level tension, barrier function, or collective alignment under shear: no "
        "cell-cell junction law exists",
    ),
    observations=("cell_shape_silhouette", "junction_intensity_map"),
)

_FIBROBLAST = _RingASpec(
    archetype_id="fibroblast-myofibroblast",
    species="Homo sapiens",
    lineage="mesenchymal fibroblast",
    identity={"archetype_class": "mesenchymal", "transition_axis": "fibroblast->myofibroblast"},
    biological_state={"activation": "quiescent-or-activated", "polarity": "front-back"},
    environment={"matrix": "3D fibrillar collagen", "mechanical_load": "self-generated traction"},
    excluded_components=frozenset(),
    missing_components=("matrix_synthesis_source", "sf_isoform_composition"),
    missing_laws=(
        MissingLaw(
            law_id="law.alpha_sma_incorporation",
            statement=(
                "a mechanosensitive law that changes stress-fiber composition (alpha-SMA "
                "incorporation) as a function of sustained tension, and thereby changes the fiber's "
                "own constitutive response"
            ),
            blocked_by=(
                "sf_arc is an 'active rod/cable graph' with fixed material parameters; the census "
                "has no composition state on a fiber and no law that evolves it"
            ),
            nearest_existing="nmii_sf_motor",
        ),
        MissingLaw(
            law_id="law.matrix_deposition_and_degradation",
            statement=(
                "a source/sink law by which the cell adds and removes ECM material, so matrix "
                "topology is an OUTPUT of cell activity rather than a fixed world"
            ),
            blocked_by=(
                "ecm_crosslink is a crosslink kinetics edge on an existing fiber population; there "
                "is no fiber creation or removal law and no secretion path"
            ),
            nearest_existing="ecm_crosslink",
        ),
        MissingLaw(
            law_id="law.plastic_matrix_remodelling",
            statement="irreversible (plastic) matrix reorganisation under sustained cell traction",
            blocked_by=(
                "ecm_mechanics declares an elastic fiber/crosslink world; permanent set is not "
                "representable without a plasticity state"
            ),
        ),
    ),
    unsupported_claims=(
        _MEDIUM_DECLARED_ONLY,
        _TURGOR_GAP_CLAIM,
        "any myofibroblast-transition threshold or timescale: the composition law does not exist",
    ),
    observations=("traction_field", "cell_shape_silhouette", "fiber_alignment_map"),
)

_EPITHELIAL = _RingASpec(
    archetype_id="junctional-epithelial",
    species="Homo sapiens",
    lineage="simple epithelium",
    identity={"archetype_class": "epithelial", "junction_type": "adherens + tight + desmosome"},
    biological_state={"activation": "confluent", "polarity": "apicobasal"},
    environment={"matrix": "basement membrane", "neighbours": "confluent sheet"},
    excluded_components=frozenset({"lamellipodium"}),
    missing_components=("cell_cell_junction", "tight_junction_barrier", "desmosome"),
    missing_laws=(
        MissingLaw(
            law_id="law.e_cadherin_adherens_junction",
            statement=(
                "an E-cadherin junction that transmits force between two cells' cortices with "
                "catch-bond kinetics and tension-dependent reinforcement"
            ),
            blocked_by=(
                "the census declares no cell-cell connector; wiring one cell's cortex directly to "
                "another's is exactly the cross-cell defect the population validator refuses"
            ),
            nearest_existing="sf_cortex_transient",
        ),
        MissingLaw(
            law_id="law.desmosome_if_coupling",
            statement=(
                "a desmosomal junction that couples the intermediate-filament networks of two "
                "neighbouring cells into one tissue-scale load path"
            ),
            blocked_by=(
                "intermediate_filament has exactly three declared edges (if_nucleus_linc, "
                "if_sf_plectin, if_cytosol_transfer) and none of them leaves the cell"
            ),
            nearest_existing="if_sf_plectin",
        ),
        MissingLaw(
            law_id="law.tight_junction_barrier",
            statement="a selective paracellular barrier whose permeability depends on junction tension",
            blocked_by="no inter-cell flux path exists in the census",
        ),
        MissingLaw(
            law_id="law.apicobasal_polarity_field",
            statement=(
                "an apicobasal polarity field partitioning membrane and cortex into apical, "
                "lateral and basal domains with distinct properties"
            ),
            blocked_by=(
                "membrane and cortex are single homogeneous surface components with no domain "
                "partition"
            ),
        ),
    ),
    unsupported_claims=(
        _MEDIUM_DECLARED_ONLY,
        _TURGOR_GAP_CLAIM,
        "apical constriction, sheet folding, or any tissue-scale mechanics: junction and polarity "
        "laws are both absent",
    ),
    observations=("cell_shape_silhouette", "apical_area"),
)

_STEM = _RingASpec(
    archetype_id="stem-ipsc-like",
    species="Homo sapiens",
    lineage="pluripotent stem / iPSC",
    identity={"archetype_class": "pluripotent", "colony_context": "compact colony"},
    biological_state={"activation": "self-renewing", "polarity": "colony-edge vs interior"},
    environment={"matrix": "defined coating", "neighbours": "compact colony"},
    excluded_components=frozenset({"filopodium", "sf_arc"}),
    missing_components=("cell_cell_junction", "lineage_state"),
    missing_laws=(
        MissingLaw(
            law_id="law.lineage_state_to_mechanics_coupling",
            statement=(
                "a coupling from transcriptional/lineage state to mechanical parameters, so "
                "differentiation changes cortical and nuclear mechanics rather than being a label"
            ),
            blocked_by=(
                "the runtime CellState explicitly owns no physical state; it conditions rates and "
                "inventories, and there is no declared law from a lineage axis to a material "
                "parameter"
            ),
        ),
        MissingLaw(
            law_id="law.lamin_stoichiometry_nuclear_stiffness",
            statement=(
                "a law relating lamin A/C stoichiometry to nuclear envelope stiffness, which is the "
                "mechanical signature that separates a pluripotent from a differentiated nucleus"
            ),
            blocked_by=(
                "nucleus is declared as 'native lamina/chromatin FEM' with fixed material "
                "parameters; lamina composition is not a state variable"
            ),
            nearest_existing="actin_cap_linc",
        ),
        MissingLaw(
            law_id="law.colony_junction_and_edge_tension",
            statement=(
                "a colony-level junctional belt whose tension differs at the colony edge, which is "
                "the dominant mechanical feature of a pluripotent colony"
            ),
            blocked_by="the census declares no cell-cell connector",
        ),
        MissingLaw(
            law_id="law.cell_cycle_coupled_cortical_stiffening",
            statement=(
                "cell-cycle-dependent cortical tension, including mitotic rounding, driven by an "
                "explicit cycle clock"
            ),
            blocked_by=(
                "the cell-cycle axis conditions rates only; no component owns a cycle-driven "
                "mechanical transition and there is no division law"
            ),
        ),
    ),
    unsupported_claims=(
        _MEDIUM_DECLARED_ONLY,
        _TURGOR_GAP_CLAIM,
        "any differentiation-linked mechanical difference: the lineage-to-mechanics coupling does "
        "not exist",
        "colony-scale organisation: no cell-cell junction law exists",
    ),
    observations=("cell_shape_silhouette", "nuclear_shape"),
)

_PULSATILE = _RingASpec(
    archetype_id="pulsatile-contractile",
    species="Homo sapiens",
    lineage="apical-medial pulsatile actomyosin (generic)",
    identity={"archetype_class": "pulsatile", "contractility_mode": "medioapical pulses"},
    biological_state={"activation": "pulsing", "polarity": "apical"},
    environment={"matrix": "compliant", "mechanical_load": "self-generated"},
    excluded_components=frozenset({"lamellipodium", "filopodium"}),
    missing_components=("rho_signalling_field", "actomyosin_advection_field"),
    missing_laws=(
        MissingLaw(
            law_id="law.rho_excitable_medium",
            statement=(
                "an excitable RhoA/ROCK reaction-advection field with delayed negative feedback, "
                "which is what MAKES the pulses; without it, pulsatility must be imposed"
            ),
            blocked_by=(
                "the census's activation is per-head NMII kinetics (Stam-Hocky backbone + "
                "individual heads, Hill/Bell); no component owns a signalling field and no "
                "connector carries a reaction-advection flux"
            ),
            nearest_existing="nmii_cortex_motor",
        ),
        MissingLaw(
            law_id="law.mechanochemical_feedback_to_kinetics",
            statement=(
                "a feedback path from local cortical strain rate back to the signalling field that "
                "sets motor recruitment, closing the oscillator loop"
            ),
            blocked_by=(
                "connector kinetics are conditioned on load per bond (Bell/Hill); there is no "
                "declared path from a mechanical field back to a recruitment rate"
            ),
        ),
        MissingLaw(
            law_id="law.ratchet_stabilisation",
            statement=(
                "the stabilisation step that converts a reversible pulse into a net shape change, "
                "i.e. an irreversible strain accumulator"
            ),
            blocked_by=(
                "cortex is a crosslinked F-actin network with reversible transient crosslinks; no "
                "declared law accumulates irreversible strain"
            ),
            nearest_existing="sf_cortex_transient",
        ),
    ),
    unsupported_claims=(
        _MEDIUM_DECLARED_ONLY,
        _TURGOR_GAP_CLAIM,
        "pulse period, amplitude, or duty cycle: the oscillator that would set them does not exist "
        "in the census, so any pulsatility would be imposed rather than emergent",
    ),
    extra_prior_slots=(
        (
            "pulse_period_s",
            "medioapical contraction pulse period",
            "s",
            "no oscillator law exists to generate a period, so a prior on it would be a prior on "
            "an imposed forcing rather than on the cell",
        ),
    ),
    observations=("apical_area", "cortical_flow_field"),
)

_RING_A_SPECS: tuple[_RingASpec, ...] = (
    _AMOEBOID,
    _ENDOTHELIAL,
    _FIBROBLAST,
    _EPITHELIAL,
    _STEM,
    _PULSATILE,
)

#: Ring A archetype identifiers, in registry order.
RING_A_ARCHETYPES: tuple[str, ...] = tuple(spec.archetype_id for spec in _RING_A_SPECS)


_RING_BC_VERDICTS: tuple[MissingLawVerdict, ...] = (
    MissingLawVerdict(
        archetype_id="erythrocyte",
        ring=ArchetypeRing.RING_B,
        missing_components=("spectrin_membrane_skeleton", "haemoglobin_cytoplasm"),
        missing_laws=(
            MissingLaw(
                law_id="law.membrane_in_plane_shear_elasticity",
                statement=(
                    "an in-plane shear modulus on the cell surface; the RBC's defining mechanic is "
                    "shear elasticity of the membrane-skeleton composite"
                ),
                blocked_by=(
                    "membrane is declared as a 'fluid Helfrich surface FEM' -- a fluid surface has "
                    "ZERO in-plane shear modulus by construction, so the quantity is not merely "
                    "unparameterised, it is unrepresentable"
                ),
                nearest_existing="membrane",
            ),
            MissingLaw(
                law_id="law.fixed_connectivity_entropic_spring_network",
                statement=(
                    "a triangulated network of fixed 6-fold connectivity whose links are entropic "
                    "worm-like-chain springs and whose connectivity does not turn over"
                ),
                blocked_by=(
                    "cortex is declared as an 'explicit crosslinked F-actin network' with "
                    "kinetic transient crosslinks and filament turnover; a spectrin skeleton has "
                    "neither, and its elasticity is entropic rather than filament-bending"
                ),
                nearest_existing="cortex",
            ),
            MissingLaw(
                law_id="law.membrane_skeleton_fixed_anchor",
                statement=(
                    "permanent anchoring of the skeleton to the bilayer at discrete "
                    "band3/ankyrin and glycophorin/4.1 nodes"
                ),
                blocked_by=(
                    "membrane_erm_cortex is a KINETIC tether that binds and unbinds, and the "
                    "engine's own note records that a tether is force-free in compression; a "
                    "skeletal anchor is neither kinetic nor unilateral"
                ),
                nearest_existing="membrane_erm_cortex",
            ),
            MissingLaw(
                law_id="law.enucleate_cytoplasm_viscosity",
                statement=(
                    "a dense enucleate haemoglobin cytoplasm as a simple viscous fluid with no "
                    "cytoskeletal interior"
                ),
                blocked_by=(
                    "cytosol is a 'Biot/Darcy finite-volume field', i.e. flow through a porous "
                    "SOLID network; with no interior network the Darcy permeability is undefined "
                    "rather than large"
                ),
                nearest_existing="cytosol",
            ),
        ),
        unsupported_claims=(
            "membrane shear modulus, area-expansion modulus, or any ektacytometry deformability index",
            "tank-treading, parachute or slipper shapes in flow",
            "haemolysis or fragmentation thresholds",
        ),
    ),
    MissingLawVerdict(
        archetype_id="neuron-growth-cone",
        ring=ArchetypeRing.RING_B,
        missing_components=("neurite_process", "transported_cargo", "mt_polarity_field"),
        missing_laws=(
            MissingLaw(
                law_id="law.processive_motor_cargo_transport",
                statement=(
                    "processive kinesin/dynein translocation of a CARGO along the microtubule "
                    "lattice over hundreds of micrometres, with run length, stall force and "
                    "directional switching"
                ),
                blocked_by=(
                    "the census's only microtubule motor edge is mt_cortex_capture, a cortical "
                    "dynein CAPTURE site; no component owns cargo state and no connector "
                    "translocates anything along a lattice"
                ),
                nearest_existing="mt_cortex_capture",
            ),
            MissingLaw(
                law_id="law.mt_polarity_sorting",
                statement=(
                    "uniform plus-end-out axonal vs mixed dendritic microtubule polarity, which is "
                    "what makes transport directional"
                ),
                blocked_by=(
                    "microtubule is declared as a 'dynamic rod graph' with no per-filament polarity "
                    "field and no sorting law"
                ),
                nearest_existing="microtubule",
            ),
            MissingLaw(
                law_id="law.high_aspect_ratio_process_geometry",
                statement=(
                    "a process whose length exceeds the cell body by three to five orders of "
                    "magnitude, with its own membrane, cortex and transport domain"
                ),
                blocked_by=(
                    "the architecture assumes ONE closed cell surface with a single membrane and "
                    "cortex component; there is no sub-cellular domain decomposition, so an axon "
                    "cannot be declared at all"
                ),
            ),
            MissingLaw(
                law_id="law.growth_cone_substrate_clutch_over_process",
                statement=(
                    "a point-contact clutch at a growth cone that advances a process, distinct from "
                    "a focal adhesion holding a spread cell body"
                ),
                blocked_by=(
                    "focal_adhesion is declared owns_geometry=False and is wired to sf_arc, "
                    "lamellipodium, filopodium and ecm of the SAME single cell body; with no "
                    "process component there is nothing for it to advance"
                ),
                nearest_existing="lamellipodium_nascent_fa",
            ),
        ),
        unsupported_claims=(
            "axonal transport velocity, run length or flux",
            "growth-cone advance rate or turning angle",
            "neurite outgrowth length or branching statistics",
        ),
    ),
    MissingLawVerdict(
        archetype_id="ciliated-cell",
        ring=ArchetypeRing.RING_C,
        missing_components=("axoneme", "basal_body", "ciliary_membrane_compartment"),
        missing_laws=(
            MissingLaw(
                law_id="law.axonemal_dynein_sliding",
                statement=(
                    "dynein motors walking on microtubule DOUBLETS, generating inter-doublet "
                    "sliding rather than actin crossbridging"
                ),
                blocked_by=(
                    "every MOTOR-family connector in the census except mt_cortex_capture is NMII "
                    "acting on a 'live polar actin material coordinate' with Hill/Bell kinetics; "
                    "there is no motor law whose track is a microtubule doublet"
                ),
                nearest_existing="nmii_cortex_motor",
            ),
            MissingLaw(
                law_id="law.sliding_to_bending_constraint",
                statement=(
                    "the nexin/radial-spoke constraint that converts inter-doublet sliding into "
                    "axonemal bending -- without it, sliding produces telescoping, not a beat"
                ),
                blocked_by=(
                    "the census has no constraint family; its intra-component edges "
                    "(dorsal_arc_crosslink, ecm_crosslink) are kinetic crosslinks between free "
                    "material points, not a fixed geometric shear constraint"
                ),
                nearest_existing="dorsal_arc_crosslink",
            ),
            MissingLaw(
                law_id="law.self_sustained_beat_oscillator",
                statement=(
                    "a load- and curvature-dependent dynein detachment law that produces a "
                    "self-sustained periodic beat with a phase, rather than approaching a "
                    "steady state"
                ),
                blocked_by=(
                    "no component owns a periodic active stress or a phase variable; the engine's "
                    "acceptance vocabulary is built around an accepted quasi-static step, and "
                    "an oscillatory attractor has no declared owner"
                ),
            ),
            MissingLaw(
                law_id="law.external_low_reynolds_hydrodynamics",
                statement=(
                    "exterior Stokes flow around a beating appendage, including the drag anisotropy "
                    "that makes a beat produce net flow"
                ),
                blocked_by=(
                    "cytosol is an INTERIOR Biot/Darcy porous field, and extracellular_medium -- "
                    "the only exterior fluid -- is declared with owns_geometry=False and no solve "
                    "behind it"
                ),
                nearest_existing="membrane_medium_traction",
            ),
        ),
        unsupported_claims=(
            "beat frequency, waveform or amplitude",
            "generated fluid flow rate or transport efficiency",
            "ciliary bending stiffness or work output",
        ),
    ),
    MissingLawVerdict(
        archetype_id="cardiomyocyte",
        ring=ArchetypeRing.RING_C,
        missing_components=(
            "sarcomere",
            "calcium_transient_field",
            "sarcoplasmic_reticulum",
            "intercalated_disc",
        ),
        missing_laws=(
            MissingLaw(
                law_id="law.sarcomeric_registration",
                statement=(
                    "a periodic register of thick and thin filaments anchored at Z-discs and "
                    "M-lines, with lattice spacing coupled to sarcomere length"
                ),
                blocked_by=(
                    "sf_arc is an 'active rod/cable graph'; a graph of rods has no periodic "
                    "register, no Z-disc, and no filament-lattice geometry, so length-tension is "
                    "not derivable from it"
                ),
                nearest_existing="sf_arc",
            ),
            MissingLaw(
                law_id="law.calcium_activated_thin_filament_cooperativity",
                statement=(
                    "Ca2+ binding to troponin C and cooperative tropomyosin shifting that gates "
                    "crossbridge availability along the thin filament"
                ),
                blocked_by=(
                    "activation in the census is per-head NMII duty ratio under Hill/Bell; there "
                    "is no Ca2+ field component, no regulatory protein state, and no cooperative "
                    "neighbour coupling along a filament"
                ),
                nearest_existing="nmii_sf_motor",
            ),
            MissingLaw(
                law_id="law.titin_passive_nonlinear_elasticity",
                statement=(
                    "a strongly nonlinear passive element in parallel with the active one, which "
                    "sets diastolic stiffness and restoring force"
                ),
                blocked_by=(
                    "the census's passive elements are cortical network mechanics and cable/rod "
                    "graphs; no parallel sarcomeric passive element is declared"
                ),
            ),
            MissingLaw(
                law_id="law.excitation_contraction_coupling",
                statement=(
                    "membrane excitation propagating to a Ca2+ release event that drives a periodic "
                    "whole-cell active stress with a physiological period"
                ),
                blocked_by=(
                    "the membrane is a mechanical Helfrich surface with no electrical state, and "
                    "no component owns a periodic active stress"
                ),
            ),
            MissingLaw(
                law_id="law.intercalated_disc_force_transmission",
                statement=(
                    "end-to-end mechanical and electrical coupling between adjacent myocytes, so a "
                    "tissue contracts as one syncytium"
                ),
                blocked_by="the census declares no cell-cell connector of any family",
            ),
        ),
        unsupported_claims=(
            "twitch force, peak systolic stress or sarcomere shortening velocity",
            "force-frequency or length-tension relationships",
            "any tissue-scale contraction: the intercalated disc law does not exist",
        ),
    ),
)

#: Ring B/C archetype identifiers, in registry order.
RING_BC_ARCHETYPES: tuple[str, ...] = tuple(verdict.archetype_id for verdict in _RING_BC_VERDICTS)

#: Every archetype the registry knows.
ARCHETYPE_IDS: tuple[str, ...] = (*RING_A_ARCHETYPES, *RING_BC_ARCHETYPES)

_RING_A_BY_ID: Mapping[str, _RingASpec] = MappingProxyType(
    {spec.archetype_id: spec for spec in _RING_A_SPECS}
)
_RING_BC_BY_ID: Mapping[str, MissingLawVerdict] = MappingProxyType(
    {verdict.archetype_id: verdict for verdict in _RING_BC_VERDICTS}
)


def _build_ring_a(spec: _RingASpec, grammar: EngineGrammar) -> ArchetypeBuild:
    """Build one Ring A archetype from the induced subgraph of the canonical census."""
    components = tuple(sorted(grammar.component_names - spec.excluded_components))
    if not components:
        raise ValueError(f"{spec.archetype_id!r} excludes every declared component")
    isolated = grammar.isolated_components(components)
    if isolated:
        raise ValueError(
            f"{spec.archetype_id!r} declares components no connector touches {list(isolated)}; "
            "co-location in a component list is not a connection"
        )
    connectors = tuple(connector.name for connector in grammar.induced_connectors(components))
    priors = _population_priors(spec.archetype_id, spec.extra_prior_slots)
    applicability = ApplicabilityBlock(
        missing_components=spec.missing_components,
        missing_laws=spec.missing_laws,
        unsupported_claims=spec.unsupported_claims,
        pi_gap_prior_ids=tuple(
            prior.prior_id for prior in priors if prior.status is PriorStatus.PI_GAP
        ),
    )
    manifest = CellStateManifest(
        manifest_id=spec.archetype_id,
        species=spec.species,
        lineage=spec.lineage,
        identity=dict(spec.identity),
        biological_state=dict(spec.biological_state),
        environment=dict(spec.environment),
        components=components,
        connectors=connectors,
        observations=spec.observations,
        missing_components=spec.missing_components,
        missing_laws=applicability.law_ids,
        provenance=(
            f"canonical census: {grammar.source}",
            "aleph.virtual_cell.archetypes",
        ),
    )
    return ArchetypeBuild(
        archetype_id=spec.archetype_id,
        ring=ArchetypeRing.RING_A,
        supported=True,
        priors=priors,
        applicability=applicability,
        manifest=manifest,
    )


def _build_ring_bc(verdict: MissingLawVerdict) -> ArchetypeBuild:
    """Return the refusal build for one Ring B/C archetype.

    The refusal is the deliverable.  No manifest is produced and no prior carries a number: a prior
    on a quantity whose governing law does not exist would be a prior on nothing.
    """
    priors = _population_priors(verdict.archetype_id)
    applicability = ApplicabilityBlock(
        missing_components=verdict.missing_components,
        missing_laws=verdict.missing_laws,
        unsupported_claims=verdict.unsupported_claims,
        pi_gap_prior_ids=tuple(prior.prior_id for prior in priors),
    )
    return ArchetypeBuild(
        archetype_id=verdict.archetype_id,
        ring=verdict.ring,
        supported=False,
        priors=priors,
        applicability=applicability,
        verdict=verdict,
    )


def build_archetype(archetype_id: str, *, grammar: EngineGrammar | None = None) -> ArchetypeBuild:
    """Build one archetype against the canonical engine grammar.

    Args:
        archetype_id: One of :data:`ARCHETYPE_IDS`.
        grammar: The canonical grammar; loaded from the engine source when omitted.

    Returns:
        An :class:`ArchetypeBuild`.  Ring A returns ``supported=True`` with a manifest; Ring B/C
        returns ``supported=False`` with a :class:`MissingLawVerdict` naming the specific laws the
        current component set cannot express.

    Raises:
        KeyError: If ``archetype_id`` is not registered.
        ~aleph.virtual_cell.grammar.CensusUnavailableError: If the canonical census cannot be read
            and no grammar was supplied.
    """
    identifier = _nonempty(archetype_id, what="archetype_id")
    if identifier in _RING_BC_BY_ID:
        return _build_ring_bc(_RING_BC_BY_ID[identifier])
    if identifier not in _RING_A_BY_ID:
        raise KeyError(f"unknown archetype {identifier!r}; known: {list(ARCHETYPE_IDS)}")
    resolved = grammar if grammar is not None else load_engine_grammar()
    if not isinstance(resolved, EngineGrammar):
        raise TypeError("grammar must be an EngineGrammar")
    return _build_ring_a(_RING_A_BY_ID[identifier], resolved)


def build_all_archetypes(*, grammar: EngineGrammar | None = None) -> tuple[ArchetypeBuild, ...]:
    """Build every registered archetype, Ring A first.

    Args:
        grammar: The canonical grammar; loaded once from the engine source when omitted.

    Returns:
        One :class:`ArchetypeBuild` per entry in :data:`ARCHETYPE_IDS`, in that order.
    """
    resolved = grammar if grammar is not None else load_engine_grammar()
    return tuple(build_archetype(archetype_id, grammar=resolved) for archetype_id in ARCHETYPE_IDS)

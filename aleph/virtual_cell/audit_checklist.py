"""Machine-checkable register of the adversarial audit A1-A15 and the three §10.3 removals.

§10.2 of ``FFN_PROBABILISTIC_VIRTUAL_CELL_MASTER_PLAN.md`` is prose.  Prose is re-read by whoever
remembers to and silently violated by whoever does not, which is the failure mode §10.1 says the
audit exists to prevent.  This module is the same fifteen items plus the three prohibitions of
§10.3 as data, with one property that prose cannot have: **the status of each item is DERIVED, never
asserted.**

Three statuses, and the distinction is the whole design:

* :data:`AuditStatus.MECHANICALLY_ENFORCED` — code makes the violation impossible or raises.  The
  item must name the module and callable, at least one binding must declare the exception it
  actually raises, and a test in ``tests/ac/virtual_cell/test_audit_checklist.py`` must fire that
  enforcement.  It must also state a ``residual_gap``: an enforcement claim with no stated residual
  is an assertion wearing a status.
* :data:`AuditStatus.EVIDENCED` — a measurement supports the item but nothing prevents regression.
  The number and where it was measured are carried in the item, and
  :attr:`Measurement.reproduced_here` says whether this register's own test re-derives it or quotes
  a lane result.
* :data:`AuditStatus.UNENFORCED` — nothing checks it.  This is a first-class, prominently reported
  status.  An audit register whose unenforced items are invisible is the failure it exists to
  prevent, so :func:`unenforced_gap_queue` and :func:`format_audit_register` list every one of them
  with its expansion route.

Most items are UNENFORCED or EVIDENCED, and that is the useful output: §10.1 states the audit "does
not ask what to discard", so the unenforced set is the queue of what to build next, not a list of
things to drop.

**Binding, not scanning.**  Enforcement targets are bound by importing the actual object and holding
it (:attr:`EnforcementBinding.target`).  No string is grepped out of the repository to decide a
status, because a heuristic scan degrades silently when a symbol is renamed, whereas an import
fails loudly at collection time.  :meth:`EnforcementBinding.resolve` re-resolves each target through
``importlib`` and checks object identity, so a rename or a moved symbol breaks the register instead
of quietly downgrading an item to UNENFORCED.

**This module is ``oracle-only`` bookkeeping.**  It runs no physics, measures no cell, grants no
evidence rung and closes no gate.  A MECHANICALLY_ENFORCED row means one named callable refuses one
named violation; it does not mean the audit item is discharged, and it is never native evidence.
The register is a document about code, and it is wrong exactly when the code moves without it.

Sanity Gate (recorded before first execution):

* **Dimensional** — the register holds no physical quantity of its own.  Every
  :class:`Measurement` carries an explicit ``units`` string (``"1"`` for dimensionless ratios and
  fractions), so a ratio can never be read as a length, an error, or a time.
* **Boundary** — an item with no enforcement and no measurement is UNENFORCED rather than absent;
  an item whose enforcement tuple is non-empty must carry a ``residual_gap`` and at least one
  binding that names a raised exception; ``partial_enforcement`` is admissible only on items that
  are NOT enforced, because machinery that closes the item belongs in ``enforcement`` and machinery
  that does not must not be able to masquerade as it.
* **Conservation / invariant** — the eighteen ids ``A1..A15`` and ``R1..R3`` are exactly the
  registered ids, each appears once, and the three status counts sum to eighteen.  Status is a
  computed property with no setter, so a row cannot claim a status its content does not support.
* **Measurement protocol** — a measurement records ``quantity``, ``value``, ``units``, the
  ``protocol`` that produced it, ``measured_in`` (document section, lane, or module), and
  ``reproduced_here``.  A number quoted from another lane is marked as quoted; a number this
  register's tests re-derive is marked as reproduced, and the test asserts the stored value.
* **Numerical** — every measurement value must be finite.  No tolerance, threshold, or gate value
  is defined by this module; the numbers it stores are outcomes, and none of them is compared
  against a limit chosen here.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from importlib import import_module
from types import MappingProxyType
from typing import Any

from aleph.virtual_cell.archetypes import MissingLawVerdict
from aleph.virtual_cell.cell_state_posterior import (
    CellStatePosterior,
    LossyMergeRefused,
    UncertaintyBudget,
)
from aleph.virtual_cell.contracts import (
    ArtifactEnvelope,
    CellStateManifest,
    RepresentationDescriptor,
    RepresentationKind,
    SandboxExperimentCard,
)
from aleph.virtual_cell.grammar import CensusUnavailableError, load_engine_grammar
from aleph.virtual_cell.observation_operator import ObservationRefused, UnimplementedOperator
from aleph.virtual_cell.sandbox_inverse import UnidentifiableDirectionError
from aleph.virtual_cell.source_registry import EvidenceGraph, HeldOutValidationRefused
from aleph.virtual_cell.surrogate import (
    GateAuthorityError,
    GateEvidenceLedger,
    SurrogateAuthorityError,
    SurrogatePrediction,
    SurrogateResponse,
)
from aleph.virtual_cell.telemetry import AxisOrderSearchReport, search_axis_order

__all__ = [
    "AUDIT_IDS",
    "AUDIT_REGISTER",
    "AuditItem",
    "AuditKind",
    "AuditRegisterError",
    "AuditRegisterSummary",
    "AuditStatus",
    "EnforcementBinding",
    "GapEntry",
    "Measurement",
    "REMOVAL_IDS",
    "audit_item",
    "audit_register_to_dict",
    "enforcement_target_paths",
    "format_audit_register",
    "items_with_status",
    "resolve_enforcement_targets",
    "status_counts",
    "summarize_audit_register",
    "unenforced_gap_queue",
]

_AUDIT_ID = re.compile(r"^(?:A(?:[1-9]|1[0-5])|R[1-3])$")


class AuditRegisterError(RuntimeError):
    """An enforcement binding no longer resolves to the object the register recorded."""


class AuditStatus(StrEnum):
    """How well an audit item is actually held, derived from the item's own content."""

    MECHANICALLY_ENFORCED = "mechanically-enforced"
    EVIDENCED = "evidenced"
    UNENFORCED = "unenforced"


class AuditKind(StrEnum):
    """Whether an item routes a failure to an expansion (§10.2) or forbids a thing (§10.3)."""

    EXPANSION = "expansion"
    PROHIBITION = "prohibition"


def _one_line(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    text = value.strip()
    if "\n" in text:
        raise ValueError(f"{what} must be a single line so the register prints as a table")
    return text


def _lines(values: object, *, what: str) -> tuple[str, ...]:
    if isinstance(values, str) or not isinstance(values, Sequence):
        raise TypeError(f"{what} must be a sequence of strings")
    return tuple(_one_line(value, what=f"{what} entry") for value in values)


@dataclass(frozen=True, slots=True)
class EnforcementBinding:
    """One live object that refuses one named violation, held by reference rather than by name.

    ``target`` is the imported object itself.  ``attribute`` names the method or property on it when
    the enforcement lives there.  Neither is a string looked up lazily against the repository: a
    renamed symbol breaks the import of this module, and :meth:`resolve` additionally checks that
    the dotted path still leads back to the same object.
    """

    target: Any
    guarantee: str
    proof_test: str
    attribute: str | None = None
    raises: type[BaseException] | None = None

    def __post_init__(self) -> None:
        module_name = getattr(self.target, "__module__", None)
        qualname = getattr(self.target, "__qualname__", None)
        if not isinstance(module_name, str) or not isinstance(qualname, str):
            raise TypeError("target must be an imported class or function, not a name or a value")
        if not module_name.startswith("aleph.virtual_cell."):
            raise ValueError(
                f"{module_name}.{qualname} is outside aleph/virtual_cell; this register only "
                "binds enforcement it can import without a GPU or a run record"
            )
        object.__setattr__(self, "guarantee", _one_line(self.guarantee, what="guarantee"))
        proof = _one_line(self.proof_test, what="proof_test")
        if not proof.startswith("test_"):
            raise ValueError("proof_test must name a test function, so the claim is checkable")
        object.__setattr__(self, "proof_test", proof)
        if self.attribute is not None:
            attribute = _one_line(self.attribute, what="attribute")
            if not hasattr(self.target, attribute):
                raise AttributeError(
                    f"{module_name}.{qualname} no longer exposes {attribute!r}; the register is "
                    "stale and must be corrected rather than downgraded"
                )
            object.__setattr__(self, "attribute", attribute)
        if self.raises is not None and not (
            isinstance(self.raises, type) and issubclass(self.raises, BaseException)
        ):
            raise TypeError("raises must be an exception type or None")

    @property
    def module_name(self) -> str:
        """Dotted module the enforcement lives in."""
        return str(self.target.__module__)

    @property
    def dotted_path(self) -> str:
        """Fully qualified path of the enforcing object, attribute included when there is one."""
        base = f"{self.target.__module__}.{self.target.__qualname__}"
        return base if self.attribute is None else f"{base}.{self.attribute}"

    def resolve(self) -> Any:
        """Re-resolve this binding through ``importlib`` and confirm object identity.

        Returns:
            The object the dotted path leads to, which must be the recorded ``target``.

        Raises:
            AuditRegisterError: If the module or any qualname part is gone, or if the path now
                leads to a different object.  Both are register defects, not audit outcomes.
        """
        try:
            found: Any = import_module(self.module_name)
            for part in str(self.target.__qualname__).split("."):
                found = getattr(found, part)
        except (ImportError, AttributeError) as error:
            raise AuditRegisterError(f"{self.dotted_path} no longer resolves: {error}") from error
        if found is not self.target:
            raise AuditRegisterError(
                f"{self.dotted_path} resolves to a different object than the register bound"
            )
        if self.attribute is not None and not hasattr(found, self.attribute):
            raise AuditRegisterError(f"{self.dotted_path} lost attribute {self.attribute!r}")
        return found

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of this binding."""
        return {
            "target": self.dotted_path,
            "guarantee": self.guarantee,
            "raises": None if self.raises is None else self.raises.__name__,
            "proof_test": self.proof_test,
        }


@dataclass(frozen=True, slots=True)
class Measurement:
    """One number that supports an audit item, with its units, protocol, and provenance."""

    quantity: str
    value: float
    units: str
    protocol: str
    measured_in: str
    reproduced_here: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "quantity", _one_line(self.quantity, what="quantity"))
        object.__setattr__(self, "units", _one_line(self.units, what="units"))
        object.__setattr__(self, "protocol", _one_line(self.protocol, what="protocol"))
        object.__setattr__(self, "measured_in", _one_line(self.measured_in, what="measured_in"))
        value = float(self.value)
        if value != value or value in (float("inf"), float("-inf")):
            raise ValueError("a measurement value must be finite")
        object.__setattr__(self, "value", value)
        if not isinstance(self.reproduced_here, bool):
            raise TypeError("reproduced_here must be bool")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of this measurement."""
        return {
            "quantity": self.quantity,
            "value": self.value,
            "units": self.units,
            "protocol": self.protocol,
            "measured_in": self.measured_in,
            "reproduced_here": self.reproduced_here,
        }


@dataclass(frozen=True, slots=True)
class AuditItem:
    """One audit finding whose status is computed from what it actually carries."""

    audit_id: str
    kind: AuditKind
    attack: str
    verdict: str
    expansion_route: str
    enforcement: tuple[EnforcementBinding, ...] = ()
    measurements: tuple[Measurement, ...] = ()
    partial_enforcement: tuple[str, ...] = ()
    residual_gap: str | None = None
    opens_lane: str | None = None

    def __post_init__(self) -> None:
        audit_id = _one_line(self.audit_id, what="audit_id")
        if not _AUDIT_ID.fullmatch(audit_id):
            raise ValueError(f"audit_id {audit_id!r} must be one of A1..A15 or R1..R3")
        object.__setattr__(self, "audit_id", audit_id)
        if not isinstance(self.kind, AuditKind):
            raise TypeError("kind must be an AuditKind")
        expected = AuditKind.EXPANSION if audit_id.startswith("A") else AuditKind.PROHIBITION
        if self.kind is not expected:
            raise ValueError(f"{audit_id} must be registered as {expected.value}")
        for name in ("attack", "verdict", "expansion_route"):
            object.__setattr__(self, name, _one_line(getattr(self, name), what=name))
        bindings = tuple(self.enforcement)
        if any(not isinstance(binding, EnforcementBinding) for binding in bindings):
            raise TypeError("enforcement must contain EnforcementBinding values")
        object.__setattr__(self, "enforcement", bindings)
        measurements = tuple(self.measurements)
        if any(not isinstance(item, Measurement) for item in measurements):
            raise TypeError("measurements must contain Measurement values")
        object.__setattr__(self, "measurements", measurements)
        object.__setattr__(
            self,
            "partial_enforcement",
            _lines(self.partial_enforcement, what="partial_enforcement"),
        )
        if self.residual_gap is not None:
            object.__setattr__(self, "residual_gap", _one_line(self.residual_gap, what="residual"))
        if self.opens_lane is not None:
            object.__setattr__(self, "opens_lane", _one_line(self.opens_lane, what="opens_lane"))

        if bindings:
            if self.residual_gap is None:
                raise ValueError(
                    f"{audit_id} claims mechanical enforcement without naming what is still open; "
                    "an enforcement claim with no residual gap is an assertion"
                )
            if not any(binding.raises is not None for binding in bindings):
                raise ValueError(
                    f"{audit_id} claims mechanical enforcement but no binding names an exception "
                    "it raises; a status that nothing can fail is not enforcement"
                )
            if self.partial_enforcement:
                raise ValueError(
                    f"{audit_id} is enforced, so partial machinery belongs in residual_gap"
                )
        elif self.residual_gap is not None:
            raise ValueError(
                f"{audit_id} has no enforcement, so the whole item is the gap; use "
                "partial_enforcement for machinery that does not close it"
            )

    @property
    def status(self) -> AuditStatus:
        """Derived status: enforcement wins, then measurement, otherwise the gap is reported."""
        if self.enforcement:
            return AuditStatus.MECHANICALLY_ENFORCED
        if self.measurements:
            return AuditStatus.EVIDENCED
        return AuditStatus.UNENFORCED

    @property
    def enforcement_paths(self) -> tuple[str, ...]:
        """Dotted paths of every object this item binds."""
        return tuple(binding.dotted_path for binding in self.enforcement)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of this item, status included."""
        return {
            "audit_id": self.audit_id,
            "kind": self.kind.value,
            "status": self.status.value,
            "attack": self.attack,
            "verdict": self.verdict,
            "expansion_route": self.expansion_route,
            "enforcement": [binding.to_dict() for binding in self.enforcement],
            "measurements": [item.to_dict() for item in self.measurements],
            "partial_enforcement": list(self.partial_enforcement),
            "residual_gap": self.residual_gap,
            "opens_lane": self.opens_lane,
        }


@dataclass(frozen=True, slots=True)
class GapEntry:
    """One UNENFORCED item, reduced to what a session needs to pick it up."""

    audit_id: str
    kind: AuditKind
    attack: str
    expansion_route: str
    partial_enforcement: tuple[str, ...]
    opens_lane: str | None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of this gap."""
        return {
            "audit_id": self.audit_id,
            "kind": self.kind.value,
            "attack": self.attack,
            "expansion_route": self.expansion_route,
            "partial_enforcement": list(self.partial_enforcement),
            "opens_lane": self.opens_lane,
        }


@dataclass(frozen=True, slots=True)
class AuditRegisterSummary:
    """Counts by status plus the full gap queue, so printing the register shows the queue."""

    counts: Mapping[AuditStatus, int] = field(default_factory=dict)
    enforced_ids: tuple[str, ...] = ()
    evidenced_ids: tuple[str, ...] = ()
    unenforced: tuple[GapEntry, ...] = ()
    enforcement_paths: Mapping[str, tuple[str, ...]] = field(default_factory=dict)

    @property
    def total(self) -> int:
        """Number of registered items; the three counts sum to this."""
        return sum(self.counts.values())

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able summary."""
        return {
            "counts": {status.value: count for status, count in self.counts.items()},
            "total": self.total,
            "enforced_ids": list(self.enforced_ids),
            "evidenced_ids": list(self.evidenced_ids),
            "unenforced": [entry.to_dict() for entry in self.unenforced],
            "enforcement_paths": {
                audit_id: list(paths) for audit_id, paths in self.enforcement_paths.items()
            },
        }


# --- the register -------------------------------------------------------------------------------
#
# Every enforcement guarantee below was written after firing the enforcement, not before.  Where an
# item could plausibly be argued into MECHANICALLY_ENFORCED but the machinery does not actually
# refuse the attack, it is left UNENFORCED with the machinery recorded in ``partial_enforcement``.

_A1 = AuditItem(
    audit_id="A1",
    kind=AuditKind.EXPANSION,
    attack="Merging thermal, KMC, cell-to-cell, parameter and measurement uncertainty into one "
    "covariance erases both the cause and the timescale.",
    verdict="A single universal Gaussian state is rejected.",
    expansion_route="Source-factored uncertainty plus an ensemble/mixture/flow/TT representation "
    "selector in W3; the decomposition is itself an uncertainty-decomposition paper.",
    measurements=(
        Measurement(
            quantity="ratio of E[x^4] under a two-source mixture to E[x^4] under the single "
            "Gaussian with the same total variance",
            value=1.9607881580237234,
            units="1",
            protocol="closed form: equal-weight zero-mean Gaussian mixture with sigma 0.1 and 1.0 "
            "has E[x^4] = 1.5*(sigma1^4 + sigma2^4) = 1.50015; the moment-matched single Gaussian "
            "(variance 0.505) has 3*var^2 = 0.765075",
            measured_in="aleph/tests/ac/virtual_cell/test_audit_checklist.py",
            reproduced_here=True,
        ),
        Measurement(
            quantity="scalar time average of a 120-cycle limit cycle as a fraction of its own "
            "cycle amplitude",
            value=0.0010204410973641743,
            units="1",
            protocol="phase_decomposition on an analytic limit cycle: time_average 0.001897 "
            "against cycle_amplitude 1.859",
            measured_in="master plan §15.9 (Lane C), aleph/virtual_cell/sweep.py",
            reproduced_here=False,
        ),
    ),
    partial_enforcement=(
        "contracts.UncertaintySource enumerates eight mechanistically distinct sources, but nothing "
        "requires a descriptor to factor an uncertainty across them.",
        "contracts.RepresentationDescriptor forces a Gaussian representation to carry a measured "
        "error report and a fallback, which is R2's prohibition rather than A1's decomposition.",
    ),
    opens_lane="W3 source-factored uncertainty and representation selector",
)

_A2 = AuditItem(
    audit_id="A2",
    kind=AuditKind.EXPANSION,
    attack="Materialising a density grid over the whole native continuous state is impossible.",
    verdict="Full density materialization is rejected; the probabilistic direction is kept.",
    expansion_route="Path ensembles, local closures, latent collective coordinates and a tensorized "
    "parameter/state subspace held as parallel representations; where each breaks is the W4 question.",
    measurements=(
        Measurement(
            quantity="global relative Frobenius error of the accepted rank-7 SVD of jacobian_fd "
            "(2207, 8) at tolerance 0.032",
            value=0.022928,
            units="1",
            protocol="telemetry.gate_compression over the committed archive sf_operator_observables"
            ".npz; the global norm accepted the reduction",
            measured_in="master plan §15.6 (Lane B)",
            reproduced_here=False,
        ),
        Measurement(
            quantity="local relative error at entry [521, 5] of that same accepted reduction",
            value=0.999846,
            units="1",
            protocol="declared VALUE projection on the reconstructed array: original 0.972998 "
            "against reconstructed 1.502616e-04, max|J| = 1.0007",
            measured_in="master plan §15.6 (Lane B)",
            reproduced_here=False,
        ),
    ),
    partial_enforcement=(
        "contracts.RepresentationKind is a closed vocabulary with no full-density member, so a full "
        "density cannot be DECLARED — but nothing stops a module computing one internally.",
        "reduction.select_array_representation and probability.audit_probability_mass fail closed to "
        "dense, which bounds an approximation rather than the dimension of the state space.",
    ),
    opens_lane="W4 parallel representation laboratory",
)

_A3 = AuditItem(
    audit_id="A3",
    kind=AuditKind.EXPANSION,
    attack="TT rank depends on axis order and tolerance, so reading it as a physical quantity is a "
    "category error.",
    verdict="Direct identification such as rank = correlation length is forbidden.",
    expansion_route="Representation-complexity atlas, axis-order robustness and entanglement-like "
    "diagnostics kept as their own research axis.",
    enforcement=(
        EnforcementBinding(
            target=AxisOrderSearchReport,
            guarantee="a search report cannot claim 'exhaustive' while leaving orderings "
            "unevaluated, and unsearched_fraction is derived on every report",
            proof_test="test_a3_axis_order_report_refuses_a_false_exhaustive_claim",
            raises=ValueError,
        ),
        EnforcementBinding(
            target=search_axis_order,
            guarantee="a capped search reports the fraction of the ordering space it never "
            "measured, so a rank is never quotable as THE rank",
            proof_test="test_a3_capped_axis_order_search_reports_its_unsearched_fraction",
        ),
    ),
    measurements=(
        Measurement(
            quantity="distinct TT storage sizes over the 24 axis orderings of one (2, 3, 4, 5) array",
            value=5.0,
            units="1",
            protocol="telemetry.search_axis_order at relative tolerance 1e-9, exhaustive over 4! "
            "orderings; storage sizes span 38 to 46 float64 slots for the same array",
            measured_in="aleph/tests/ac/virtual_cell/test_audit_checklist.py",
            reproduced_here=True,
        ),
        Measurement(
            quantity="accepted SVD rank of jacobian_fd (2207, 8) at tolerance 0.032 and at 0.1",
            value=7.0,
            units="1",
            protocol="same array, two declared tolerances: rank 7 at 0.032 and rank 6 at 0.1, so "
            "the rank is a property of the request as much as of the array",
            measured_in="master plan §14 rehearsal and §15.6 (Lane B)",
            reproduced_here=False,
        ),
    ),
    residual_gap="Nothing forbids prose that calls a rank a correlation length; the enforcement "
    "constrains the report, not the sentence written about it, and no axis-order robustness atlas "
    "exists yet.",
)

_A4 = AuditItem(
    audit_id="A4",
    kind=AuditKind.EXPANSION,
    attack="A fast surrogate becomes convenient and starts producing more conclusions than the "
    "reference it was fitted to.",
    verdict="A surrogate is confined to search/proposal/initialisation and its evidence source is "
    "forced.",
    expansion_route="Learned preconditioner, trust-region monitor, native-query acquisition and a "
    "surrogate failure atlas as separate deliverables.",
    enforcement=(
        EnforcementBinding(
            target=SurrogatePrediction,
            attribute="evidence_source",
            guarantee="evidence_source is a property over __slots__ and not a field, so "
            "dataclasses.replace(prediction, evidence_source=NATIVE_ACCEPTED) cannot be spelled",
            proof_test="test_a4_a_surrogate_prediction_cannot_be_relabelled_as_native",
            raises=TypeError,
        ),
        EnforcementBinding(
            target=SurrogateResponse,
            attribute="value",
            guarantee="reading a number off a refused query raises instead of degrading the "
            "refusal into a guess",
            proof_test="test_a4_a_refused_surrogate_query_has_no_value",
            raises=SurrogateAuthorityError,
        ),
        EnforcementBinding(
            target=GateEvidenceLedger,
            attribute="assert_gate_closable",
            guarantee="an artifact with surrogate provenance anywhere in its transitive ancestry "
            "is refused as gate evidence, and the laundering path is named",
            proof_test="test_a4_r3_surrogate_provenance_is_refused_transitively",
            raises=SurrogateAuthorityError,
        ),
    ),
    residual_gap="Nothing counts how many conclusions came from the surrogate against the "
    "reference; the failure atlas and the learned preconditioner are unbuilt, and training itself "
    "is blocked until an accepted native trace exists.",
)

_A5 = AuditItem(
    audit_id="A5",
    kind=AuditKind.EXPANSION,
    attack="Different molecular parameters can produce the same shape image.",
    verdict="Image-to-parameter point prediction is rejected.",
    expansion_route="Posterior prediction, temporal images, traction/omics fusion, virtual "
    "interventions and expected-information-gain observation design.",
    partial_enforcement=(
        "optics.render_emitters and optics.OpticsCapabilityVerdict give a forward observation model "
        "and refuse confocal/light-sheet by typed verdict — the forward half only.",
        "No inverse map from image to parameter exists in this layer, so there is nothing yet that "
        "could refuse a point estimate for an unobservable.",
    ),
    opens_lane="S11 image-inverse sandbox; W7 observation design",
)

_A6 = AuditItem(
    audit_id="A6",
    kind=AuditKind.EXPANSION,
    attack="An encoder that learns the renderer's artifacts can score itself perfectly.",
    verdict="A test inside one renderer family is a structural test, not biological validation.",
    expansion_route="Renderer holdout, multiple rendering engines, style randomisation, public "
    "real-image domain shift and adversarial artifact controls in W7.",
    enforcement=(
        EnforcementBinding(
            target=EvidenceGraph,
            attribute="citable_as_biological_evidence",
            guarantee="a statistic derived from a synthetic render is excluded from the biological "
            "corpus at arbitrary ancestry depth, naming the blocking source and the whole path",
            proof_test="test_a6_synthetic_ancestry_is_excluded_two_hops_down",
        ),
        EnforcementBinding(
            target=EvidenceGraph,
            attribute="assert_admissible_as_held_out_validation",
            guarantee="the same closure refuses a synthetic descendant, and refuses a "
            "CALIBRATION_ONLY descendant separately, as held-out validation",
            proof_test="test_a6_synthetic_and_calibration_descendants_are_refused_as_held_out",
            raises=HeldOutValidationRefused,
        ),
    ),
    residual_gap="Only one renderer family exists, so renderer holdout, style randomisation and "
    "real-image domain shift have nothing to hold out against; the exclusion is over provenance "
    "records a caller supplies honestly, and nothing yet detects an undeclared synthetic ancestor.",
)

_A7 = AuditItem(
    audit_id="A7",
    kind=AuditKind.EXPANSION,
    attack="Knocking something out inside the model and getting the result the model implies is not "
    "a biological discovery.",
    verdict="A biological claim from a single-formalism intervention alone is forbidden.",
    expansion_route="Analytic oracle, independent discretization, ablation, blind held-out "
    "intervention and public-data phenotype combined in layers; the sandbox's first product is "
    "causal discrimination and identifiability.",
    enforcement=(
        EnforcementBinding(
            target=SandboxExperimentCard,
            guarantee="an experiment cannot be pre-registered without adversarial controls, "
            "held-out interventions, negative controls, a falsifier and a failure expansion",
            proof_test="test_a7_a_sandbox_card_without_adversarial_controls_is_refused",
            raises=ValueError,
        ),
        EnforcementBinding(
            target=GateEvidenceLedger,
            attribute="assert_gate_closable",
            guarantee="no artifact this layer holds can close a gate at all, because the generic "
            "virtual-cell layer is an unverified observer",
            proof_test="test_a7_a9_no_generic_layer_artifact_can_close_a_gate",
            raises=GateAuthorityError,
        ),
    ),
    residual_gap="The card records that controls were DECLARED; nothing checks that the declared "
    "controls were run, or that the intervention was executed in a second independent formalism.",
)

_A8 = AuditItem(
    audit_id="A8",
    kind=AuditKind.EXPANSION,
    attack="An actomyosin-centred mammalian grammar is insufficient for RBC, neuron, ciliated, "
    "plant and fungal cells.",
    verdict="The current component set is not called universal.",
    expansion_route="Ring A/B/C plus a missing-component-law verdict; each application failure "
    "opens a new component-law project.",
    enforcement=(
        EnforcementBinding(
            target=MissingLawVerdict,
            guarantee="a refusal must name at least one specific missing law, so a generic "
            "'unsupported' verdict cannot be recorded",
            proof_test="test_a8_a_refusal_without_a_named_law_is_not_a_verdict",
            raises=ValueError,
        ),
        EnforcementBinding(
            target=load_engine_grammar,
            guarantee="the component/connector census is read from engine source or not at all; "
            "there is no fallback census that could silently become the universal grammar",
            proof_test="test_a8_the_census_has_no_fallback",
            raises=CensusUnavailableError,
        ),
    ),
    residual_gap="Ring B/C verdicts are authored, not derived: nothing checks that a named missing "
    "law is still missing after the engine grammar changes, and §15.5 measured the 14-component "
    "census to be archetype-coarse, with Ring A archetypes separated only by their missing laws.",
)

_A9 = AuditItem(
    audit_id="A9",
    kind=AuditKind.EXPANSION,
    attack="Building AI on force-accepted, aliased, non-admissible traces teaches solver artifacts "
    "in fine detail.",
    verdict="W1 is the platform P0, but it does not halt the other lanes.",
    expansion_route="W4 on analytic oracles, W7 on public/synthetic toy geometry, W9 on the "
    "literature graph; only native biological claims stay locked to W1.",
    enforcement=(
        EnforcementBinding(
            target=ArtifactEnvelope,
            guarantee="no artifact from this layer can carry NATIVE_ACCEPTED or NATIVE_ATTEMPTED "
            "provenance, so an unverified observer cannot mint native evidence",
            proof_test="test_a9_the_generic_layer_cannot_declare_native_evidence",
            raises=ValueError,
        ),
        EnforcementBinding(
            target=GateEvidenceLedger,
            attribute="assert_gate_closable",
            guarantee="a gate closure attempt from this layer is refused even when no surrogate is "
            "involved, and an unregistered ancestor fails closed",
            proof_test="test_a7_a9_no_generic_layer_artifact_can_close_a_gate",
            raises=GateAuthorityError,
        ),
    ),
    residual_gap="The conjunctive accepted-step predicate lives in admissibility.py as an oracle "
    "only: §15.3 measured that the engine's two commit paths each implement a different half of it, "
    "and adopting the conjunction is an unsigned PI decision this layer may not make.",
)

_A10 = AuditItem(
    audit_id="A10",
    kind=AuditKind.EXPANSION,
    attack="Adding archetypes and interventions grows the native budget faster than linearly.",
    verdict="This is not solved by cutting the number of cell types.",
    expansion_route="Hierarchical sharing, adaptive design, tensor parameter response, learned "
    "preconditioner and multi-fidelity evidence; the compute budget is itself an algorithm-paper axis.",
    partial_enforcement=(
        "sweep.ess_interval and sweep.finite_difference_sensitivity make a sample budget auditable, "
        "but neither is connected to a native cost model.",
        "No native run has ever been launched from this layer, so there is no measured per-archetype "
        "cost for anything to check a budget against.",
    ),
    opens_lane="W5 adaptive design; W10 multi-fidelity evidence",
)

_A11 = AuditItem(
    audit_id="A11",
    kind=AuditKind.EXPANSION,
    attack="Mixing cell lines, passages, cycle phases, substrates and timepoints makes a prior that "
    "is a false average.",
    verdict="Evidence is not combined on a cell-type label alone.",
    expansion_route="CellStateManifest applicability, context-conditioned mixture prior and proxy "
    "distance added to the knowledge graph.",
    enforcement=(
        EnforcementBinding(
            target=ArtifactEnvelope,
            guarantee="an artifact must be keyed by the SHA-256 of a full manifest; a bare cell-type "
            "label such as 'MCF7' is refused as a manifest hash",
            proof_test="test_a11_an_artifact_cannot_be_keyed_by_a_cell_type_label",
            raises=ValueError,
        ),
        EnforcementBinding(
            target=CellStateManifest,
            attribute="manifest_hash",
            guarantee="the hash covers identity, biological state and environment, so records that "
            "differ only in passage or substrate stiffness get different keys",
            proof_test="test_a11_passage_and_substrate_change_the_manifest_hash",
        ),
    ),
    residual_gap="Nothing in this layer combines evidence yet, so there is no pooling operation to "
    "refuse; the context-conditioned mixture prior and the proxy distance between two manifests are "
    "unbuilt, and manifest equality is currently all-or-nothing.",
)

_A12 = AuditItem(
    audit_id="A12",
    kind=AuditKind.EXPANSION,
    attack="DFT does not give a whole-cell density or a migration parameter.",
    verdict="The QM/MD lane's scope is stated as molecular free energy and energetics.",
    expansion_route="A QM to MD/CG to connector-prior scale-bridging study, kept separate from the "
    "cell-level posterior.",
    enforcement=(
        EnforcementBinding(
            target=EvidenceGraph,
            attribute="citable_as_biological_evidence",
            guarantee="a COMPUTED_QM_MD source of computation origin is outside the empirical-origin "
            "whitelist and is excluded from the biological corpus, as is anything derived from it",
            proof_test="test_a12_a_computed_qm_md_source_is_not_biological_evidence",
        ),
        EnforcementBinding(
            target=EvidenceGraph,
            attribute="assert_admissible_as_held_out_validation",
            guarantee="the same source is refused as held-out validation for a cell-level claim",
            proof_test="test_a12_a_computed_qm_md_source_is_not_biological_evidence",
            raises=HeldOutValidationRefused,
        ),
    ),
    residual_gap="This constrains only what a QM number may be CITED as. No QM/MD module exists, so "
    "nothing bounds what such a lane would claim internally, and the scale-bridging path from a "
    "molecular free energy to a connector prior is entirely unbuilt.",
)

_A13 = AuditItem(
    audit_id="A13",
    kind=AuditKind.EXPANSION,
    attack="Segmentation, GNN, TT and renderer can each end as a separate unconnected demo.",
    verdict="Every tool must take or return a CellStatePosterior or an ObservationOperator.",
    expansion_route="A shared artifact schema and the intervention sandbox as the common test bed.",
    # FLIPPED 2026-07-29, and the two lines this replaces are worth remembering: they said the
    # verdict named an interface nobody had written.  That was TRUE and it was the register's
    # sharpest finding — a rule quantified over "every tool" that no tool could satisfy.  Both
    # types now exist and ten lanes compose into them, so the row moves; the residual gap below
    # says which half of the verdict is still unchecked, because the universal half is.
    enforcement=(
        EnforcementBinding(
            target=CellStatePosterior,
            attribute="point_estimate",
            guarantee="a direction measured as unidentifiable yields no point estimate, and a "
            "parameter with no identifiability declaration is refused too — 'nobody looked' is "
            "not 'identifiable'; the credible interval stays available either way",
            proof_test=(
                "test_the_posterior_refuses_a_point_estimate_exactly_where_sandbox_inverse_does"
            ),
            raises=UnidentifiableDirectionError,
        ),
        EnforcementBinding(
            target=UncertaintyBudget,
            attribute="collapse_to_gaussian",
            guarantee="collapsing source-factored uncertainty into one width is refused when the "
            "measured cost exceeds the declared tolerance, and the refusal carries that "
            "measurement rather than a warning",
            proof_test=(
                "test_collapsing_the_a1_budget_is_refused_and_the_refusal_carries_the_measurement"
            ),
            raises=LossyMergeRefused,
        ),
        EnforcementBinding(
            target=CellStatePosterior,
            attribute="from_archetype_build",
            guarantee="a Ring B/C archetype that REFUSED cannot be forced into a posterior; the "
            "missing laws are named instead of being silently absent",
            proof_test="test_a_ring_bc_refusal_cannot_be_forced_into_a_posterior",
            raises=ValueError,
        ),
        EnforcementBinding(
            target=UnimplementedOperator,
            attribute="observe",
            guarantee="an operator that is declared but not implemented satisfies the protocol and "
            "still returns no number, so a missing modality cannot be mistaken for a measured zero",
            proof_test="test_an_unimplemented_operator_refuses_every_state_and_returns_no_number",
            raises=ObservationRefused,
        ),
    ),
    residual_gap="The verdict quantifies over EVERY tool and nothing checks that. Ten lanes now "
    "compose INTO a CellStatePosterior, but no lane's entry point is REQUIRED to accept one, so a "
    "new module can still be written that touches neither type. Enforcing the universal half needs "
    "a contract test over virtual_cell.__all__. Separately, sweep/sandbox_invariance/"
    "sandbox_fluctuation do not compose: their verdicts are about a knob ladder or a mode "
    "spectrum, and converting those into a variance in a parameter's own unit needs a declared "
    "projection nobody has written — the posterior names them as unpropagated rather than "
    "inventing one.",
    opens_lane="W3/W12 shared posterior and observation-operator interface",
)

_A14 = AuditItem(
    audit_id="A14",
    kind=AuditKind.EXPANSION,
    attack="One pipeline can be over-split into many small papers.",
    verdict="An item with no independent falsifier, dataset or method contribution is a module, not "
    "a paper.",
    expansion_route="Fold failures into a shared benchmark/resource paper and grow the biological "
    "question by cross-archetype comparison.",
    partial_enforcement=(
        "contracts.SandboxExperimentCard forces a per-experiment falsifier, which is the nearest "
        "machine-checkable half of this verdict — but a card is not a paper claim.",
        "Nothing enumerates intended publications, so there is no object on which the "
        "falsifier/dataset/method test could be run.",
    ),
    opens_lane="PI decision: publication-unit contract",
)

_A15 = AuditItem(
    audit_id="A15",
    kind=AuditKind.EXPANSION,
    attack="A new master document per direction reproduces exactly the drift the previous audit "
    "already failed on.",
    verdict="After approval this plan must be absorbed into ROADMAP.md or replace it, never kept as "
    "a parallel canonical plan.",
    expansion_route="Move to a machine-readable capability graph and a generated status view.",
    partial_enforcement=(
        "CLAUDE.md records the ROADMAP.md reframing as PROPOSED and NOT ratified, so the master "
        "plan is still a second canonical plan running beside ROADMAP.md.",
        "This register is the first machine-readable slice of that capability graph, and it covers "
        "§10 only; no code reads ROADMAP.md or the master plan.",
    ),
    opens_lane="PI ratification of the ROADMAP.md amendment, then a generated status view",
)

_R1 = AuditItem(
    audit_id="R1",
    kind=AuditKind.PROHIBITION,
    attack="Literal quantum position dynamics — amplitudes, Born-rule readout, Hamiltonian or "
    "open-system evolution of cell state.",
    verdict="REMOVE. This is classical probability plus a quantum-INSPIRED tensor representation.",
    expansion_route="Nothing replaces it: §14 records this as a physics boundary rather than a "
    "missing feature, so it is removed and not routed.",
    partial_enforcement=(
        "contracts.RepresentationKind is a closed nine-member vocabulary containing no amplitude, "
        "wavefunction or density-matrix member, and a test in this lane pins that membership.",
        "Pinning the vocabulary stops a quantum representation being DECLARED; nothing stops a "
        "module implementing a propagator, and no code reads the prose claims made about it.",
    ),
    opens_lane="none — this is a boundary, not a gap",
)

_R2 = AuditItem(
    audit_id="R2",
    kind=AuditKind.PROHIBITION,
    attack="A universal single-Gaussian state, or a materialized full density over the native "
    "state, used as THE representation.",
    verdict="REMOVE both. Neither may be the platform's state representation.",
    expansion_route="Replaced by source-factored uncertainty with an explicit selector and measured "
    "approximation error — the A1 and A2 expansions.",
    enforcement=(
        EnforcementBinding(
            target=RepresentationKind,
            guarantee="the representation vocabulary is closed and has no full-density member, so "
            "RepresentationKind('full-density') raises rather than creating one",
            proof_test="test_r2_a_full_density_representation_cannot_be_declared",
            raises=ValueError,
        ),
        EnforcementBinding(
            target=RepresentationDescriptor,
            guarantee="a Gaussian, mixture, TT, latent-density or surrogate representation cannot "
            "be declared without a measured error report and a named fallback",
            proof_test="test_r2_a_gaussian_representation_needs_measured_error_and_a_fallback",
            raises=ValueError,
        ),
    ),
    residual_gap="Both bindings refuse the DECLARATION. A module may still compute a dense density "
    "grid or a single covariance internally; what is refused is publishing one as an artifact "
    "representation, and no dimension budget is checked anywhere.",
)

_R3 = AuditItem(
    audit_id="R3",
    kind=AuditKind.PROHIBITION,
    attack="A surrogate or synthetic dataset closing its own physics or biology gate.",
    verdict="REMOVE. Self-certification is structurally impossible, not merely discouraged.",
    expansion_route="Replaced by native/oracle evidence closing gates while surrogate and synthetic "
    "artifacts stay search, proposal and calibration inputs.",
    enforcement=(
        EnforcementBinding(
            target=GateEvidenceLedger,
            attribute="assert_gate_closable",
            guarantee="surrogate provenance is followed transitively through source_artifact_ids, "
            "an unregistered ancestor fails closed, and a cycle is refused",
            proof_test="test_a4_r3_surrogate_provenance_is_refused_transitively",
            raises=SurrogateAuthorityError,
        ),
        EnforcementBinding(
            target=EvidenceGraph,
            attribute="assert_admissible_as_held_out_validation",
            guarantee="synthetic-render and CALIBRATION_ONLY ancestry are both refused as held-out "
            "validation, so a synthetic set cannot become its own test set",
            proof_test="test_a6_synthetic_and_calibration_descendants_are_refused_as_held_out",
            raises=HeldOutValidationRefused,
        ),
        EnforcementBinding(
            target=ArtifactEnvelope,
            guarantee="a synthetic-image representation must retain SYNTHETIC provenance and a "
            "surrogate representation must retain SURROGATE provenance",
            proof_test="test_r3_synthetic_and_surrogate_representations_keep_their_provenance",
            raises=ValueError,
        ),
    ),
    residual_gap="Enforcement operates on declared provenance: a caller who omits the synthetic "
    "ancestor from source_artifact_ids is not detected, and no runtime-attested adapter exists that "
    "could supply the ancestry independently.",
)

#: Every registered item, in audit order.
_ITEMS: tuple[AuditItem, ...] = (
    _A1,
    _A2,
    _A3,
    _A4,
    _A5,
    _A6,
    _A7,
    _A8,
    _A9,
    _A10,
    _A11,
    _A12,
    _A13,
    _A14,
    _A15,
    _R1,
    _R2,
    _R3,
)

#: The A1-A15 expansion findings of §10.2, in order.
AUDIT_IDS: tuple[str, ...] = tuple(
    item.audit_id for item in _ITEMS if item.kind is AuditKind.EXPANSION
)

#: The three §10.3 prohibitions, which are removals rather than expansions.
REMOVAL_IDS: tuple[str, ...] = tuple(
    item.audit_id for item in _ITEMS if item.kind is AuditKind.PROHIBITION
)

#: The register itself, keyed by audit id.
AUDIT_REGISTER: Mapping[str, AuditItem] = MappingProxyType({item.audit_id: item for item in _ITEMS})

if len(AUDIT_REGISTER) != len(_ITEMS):  # pragma: no cover - construction invariant
    raise AuditRegisterError("duplicate audit id in the register")


def audit_item(audit_id: str) -> AuditItem:
    """Return one registered item.

    Args:
        audit_id: One of ``A1``..``A15`` or ``R1``..``R3``.

    Returns:
        The :class:`AuditItem`.

    Raises:
        KeyError: If the id is not registered; the message lists the known ids.
    """
    try:
        return AUDIT_REGISTER[audit_id]
    except KeyError:
        raise KeyError(f"unknown audit id {audit_id!r}; known: {list(AUDIT_REGISTER)}") from None


def items_with_status(status: AuditStatus) -> tuple[AuditItem, ...]:
    """Return every item whose derived status is ``status``, in audit order.

    Args:
        status: The status to select.

    Returns:
        The matching items.

    Raises:
        TypeError: If ``status`` is not an :class:`AuditStatus`.
    """
    if not isinstance(status, AuditStatus):
        raise TypeError("status must be an AuditStatus")
    return tuple(item for item in _ITEMS if item.status is status)


def status_counts() -> Mapping[AuditStatus, int]:
    """Return the count of items per status, with every status present even at zero."""
    counts = dict.fromkeys(AuditStatus, 0)
    for item in _ITEMS:
        counts[item.status] += 1
    return MappingProxyType(counts)


def unenforced_gap_queue() -> tuple[GapEntry, ...]:
    """Return every UNENFORCED item as a gap entry, in audit order.

    This is the queue §10.1 asks the audit to produce: what to build next, with the expansion each
    failure routes to.  An empty tuple would mean every item is enforced or evidenced, which is not
    the current state and must never be assumed from a summary that omits this list.
    """
    return tuple(
        GapEntry(
            audit_id=item.audit_id,
            kind=item.kind,
            attack=item.attack,
            expansion_route=item.expansion_route,
            partial_enforcement=item.partial_enforcement,
            opens_lane=item.opens_lane,
        )
        for item in items_with_status(AuditStatus.UNENFORCED)
    )


def enforcement_target_paths() -> Mapping[str, tuple[str, ...]]:
    """Return the dotted paths each enforced item binds, keyed by audit id."""
    return MappingProxyType(
        {
            item.audit_id: item.enforcement_paths
            for item in items_with_status(AuditStatus.MECHANICALLY_ENFORCED)
        }
    )


def resolve_enforcement_targets() -> Mapping[str, tuple[Any, ...]]:
    """Re-resolve every enforcement binding and confirm it still leads to the bound object.

    Returns:
        The resolved objects keyed by audit id.

    Raises:
        AuditRegisterError: If any target no longer resolves or resolves elsewhere.  A rename must
            break this register loudly rather than silently downgrade an item to UNENFORCED.
    """
    return MappingProxyType(
        {
            item.audit_id: tuple(binding.resolve() for binding in item.enforcement)
            for item in _ITEMS
            if item.enforcement
        }
    )


def summarize_audit_register() -> AuditRegisterSummary:
    """Return counts by status, the enforced/evidenced ids, and the full unenforced gap queue."""
    return AuditRegisterSummary(
        counts=status_counts(),
        enforced_ids=tuple(
            item.audit_id for item in items_with_status(AuditStatus.MECHANICALLY_ENFORCED)
        ),
        evidenced_ids=tuple(item.audit_id for item in items_with_status(AuditStatus.EVIDENCED)),
        unenforced=unenforced_gap_queue(),
        enforcement_paths=enforcement_target_paths(),
    )


def audit_register_to_dict() -> dict[str, Any]:
    """Return the whole register as a JSON-able record, suitable for an artifact sidecar."""
    return {
        "register": "master plan §10.2 A1-A15 and §10.3 removals",
        "evidence_class": "oracle-only; grants no evidence rung and closes no gate",
        "items": [item.to_dict() for item in _ITEMS],
        "summary": summarize_audit_register().to_dict(),
    }


def format_audit_register() -> str:
    """Return a printable register: status counts, enforced bindings, and the gap queue.

    The gap queue is printed last and in full, because a summary that shows only what is enforced
    is precisely the reading error this register was built to stop.
    """
    summary = summarize_audit_register()
    lines = [
        "Adversarial audit register — master plan §10.2 (A1-A15) + §10.3 removals (R1-R3)",
        "oracle-only bookkeeping: grants no evidence rung, closes no gate",
        "",
        f"total items                : {summary.total}",
        f"  mechanically enforced    : {summary.counts[AuditStatus.MECHANICALLY_ENFORCED]}  "
        f"({', '.join(summary.enforced_ids) or 'none'})",
        f"  evidenced (no guard)     : {summary.counts[AuditStatus.EVIDENCED]}  "
        f"({', '.join(summary.evidenced_ids) or 'none'})",
        f"  UNENFORCED               : {summary.counts[AuditStatus.UNENFORCED]}  "
        f"({', '.join(entry.audit_id for entry in summary.unenforced) or 'none'})",
        "",
        "MECHANICALLY ENFORCED — each binding is proven by a test that fires it",
    ]
    for item in items_with_status(AuditStatus.MECHANICALLY_ENFORCED):
        lines.append(f"  {item.audit_id}: {item.verdict}")
        for binding in item.enforcement:
            raised = "" if binding.raises is None else f" [raises {binding.raises.__name__}]"
            lines.append(f"      -> {binding.dotted_path}{raised}")
            lines.append(f"         {binding.guarantee}")
        lines.append(f"      residual gap: {item.residual_gap}")
    lines.extend(["", "EVIDENCED — a measurement supports it, nothing prevents regression"])
    for item in items_with_status(AuditStatus.EVIDENCED):
        lines.append(f"  {item.audit_id}: {item.verdict}")
        for measurement in item.measurements:
            origin = "reproduced here" if measurement.reproduced_here else "quoted"
            lines.append(
                f"      {measurement.quantity} = {measurement.value:g} "
                f"[{measurement.units}] ({origin}: {measurement.measured_in})"
            )
    lines.extend(["", "UNENFORCED — the gap queue; nothing checks these"])
    for entry in summary.unenforced:
        lines.append(f"  {entry.audit_id}: {entry.attack}")
        lines.append(f"      expansion: {entry.expansion_route}")
        if entry.opens_lane is not None:
            lines.append(f"      opens: {entry.opens_lane}")
        for note in entry.partial_enforcement:
            lines.append(f"      partial: {note}")
    return "\n".join(lines)

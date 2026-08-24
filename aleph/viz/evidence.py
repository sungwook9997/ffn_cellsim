"""Read a CUDA run artefact as evidence: per **channel**, not per case.

`ALEPH-PORT-3612`.

Lane G7 drove 49 law cases on an RTX A5000 and wrote 228 KB of JSON. The evidence is real and
unreadable without a text editor. This module turns it into a model something can draw.

Why a channel and not a case
----------------------------
A law case yields several fields -- an energy and one or two force arrays -- and **a wrong kernel can
be loud in one and silent in another**. Two of them in this artefact:

``chromatin_network[WRONG_ENERGY_HALF]``
    ``energy_pn_um`` disagrees by 8.39e6 ULP against a budget of 8. ``forces_pn`` is **bit-identical
    to the true kernel's**.

``cortex_membrane_contact[WRONG_FROZEN_NORMAL]``
    ``forces_cortex_pn`` disagrees by 5.59e5 against a budget of 2048. ``energy_pn_um`` *and*
    ``forces_membrane_pn`` are bit-identical to the true kernel's.

A rendering that shows one number per case erases exactly that. So the unit here is the
:class:`FieldRow`, and :class:`ChannelState` distinguishes *over budget* from *bit-identical to the
true kernel* from *moved but inside its gate* -- three situations a single verdict collapses into one.

Bit-identity is read off equal ULP, and §4a of the entry argues why that is sound rather than a
proxy: the reference scale is a property of the family and the state, not of the variant, so a variant
that moved any element of a field would move its ULP. The arrays are not in the artefact; §14.2 marks
that limit rather than hiding it.

Log scale, and the one value that is not on it
----------------------------------------------
The measured ULP range is nine decades wide. A linear axis draws every real result on the first pixel.
So the axis is ``log10`` -- with the exception that matters:

    ``ulp == 0.0`` is an **exact** result, not a small one.

``ALEPH-PORT-3604`` shipped an exact zero as its headline. Clamping it into the bottom decade would
draw the strongest statement in the tree as merely a good one, so :func:`log_position` **raises** on
zero and exact rows go in :data:`EXACT_LANE`, drawn left of the axis with a break.

Two margins
-----------
``pressure = ulp / budget`` is how close a result sits to its gate. ``headroom = roundoff_bound /
budget`` is how far the gate sits under the round-off bound its own ledger derived. A budget with
``headroom < THIN_HEADROOM`` is *thin*: less than a binary order of slack, so a change of mesh moves
the gate through the data. Both are computed here so the reader does not divide two columns by eye.

This module does **no** arithmetic on the physics. Every ULP is read, never recomputed, so it cannot
disagree with the runner about a measurement -- only about how it is drawn.
"""

from __future__ import annotations

import json
import math
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

__all__ = [
    "DECADE_HIGH",
    "DECADE_LOW",
    "DEFAULT_ARTEFACT_DIR",
    "DEFAULT_HOST",
    "EXACT_LANE",
    "CHANNEL_STYLE",
    "PARITY_ARTEFACT",
    "OPTIONAL_ARTEFACTS",
    "THIN_HEADROOM",
    "CaseRow",
    "ChannelState",
    "ChannelStyle",
    "DescentWindow",
    "EvidenceBundle",
    "EvidenceError",
    "FieldRow",
    "ParityEvidence",
    "ScatterEvidence",
    "discover_artefacts",
    "load_bundle",
    "log_position",
    "read_descent_window",
    "read_parity",
    "read_scatter",
]

#: The log axis spans these decades. The floor is above the smallest measured value (0.0159 ULP) on
#: purpose: a mark below it is *clamped and flagged*, which keeps that branch live rather than
#: theoretical, and a clamped mark must never read as a measured one.
DECADE_LOW = -1
DECADE_HIGH = 8

#: Exact agreement is not a point on the log axis. It gets its own lane, drawn left of the axis start
#: with a visible break. See the module docstring and `ALEPH-PORT-3612` §4b.
EXACT_LANE = "exact 0"

#: A gate closer than this factor to the round-off bound its ledger derived has no room to move.
THIN_HEADROOM = 2.0

#: Where the unattended runner writes, on the host it runs on and in the local cache.
DEFAULT_ARTEFACT_DIR = Path("~/.aleph_runner/artefacts").expanduser()

#: The ssh alias the artefacts are fetched from when they are not on this machine.
DEFAULT_HOST = "sungwook@100.110.26.26"

#: The one artefact without which there is no figure.
PARITY_ARTEFACT = "j1_law_parity.json"

#: Artefacts that add a panel. An absent one is drawn **as absent**, never dropped.
OPTIONAL_ARTEFACTS = (
    "j2_scatter_determinism.json",
    "g6_descent_window.json",
    "SUMMARY.json",
)


class EvidenceError(RuntimeError):
    """The artefact could not be read as evidence. Carries the reason a caller should print."""


class ChannelState(Enum):
    """What happened to one field of one case.

    Five states, and the reason there are five rather than two is the module docstring's whole
    subject: ``BLIND_TO_TRUE`` and ``WITHIN_BUDGET`` both look like "passed" to a per-case verdict,
    and they are completely different facts. One says the wrong kernel *did not touch this channel*;
    the other says it touched it and stayed inside the gate.
    """

    OVER_BUDGET = "over_budget"
    BLIND_TO_TRUE = "blind_to_true"
    EXACT = "exact"
    WITHIN_BUDGET = "within_budget"
    ABSENT = "absent"


@dataclass(frozen=True, slots=True)
class ChannelStyle:
    """How one channel state is drawn.

    Every field is a separate visual channel and the five states differ on more than one of them,
    following `render.py`'s rule: colour alone fails a greyscale print and a colour-blind reader.

    Attributes:
        fill: Mark fill colour.
        stroke: Mark outline colour.
        glyph: A one-character marker, so the states survive losing colour entirely.
        label: Legend text.
        is_fault: Whether this state is a defect rather than a result.
    """

    fill: str
    stroke: str
    glyph: str
    label: str
    is_fault: bool = False


#: One style per channel state. ``BLIND_TO_TRUE`` is deliberately the loudest non-fault: it is the
#: port's central finding and a reader must not skim past it.
CHANNEL_STYLE: dict[ChannelState, ChannelStyle] = {
    ChannelState.OVER_BUDGET: ChannelStyle(
        fill="#c62828",
        stroke="#7f0000",
        glyph="X",
        label="over budget — the gate caught it in this channel",
    ),
    ChannelState.BLIND_TO_TRUE: ChannelStyle(
        fill="#6a1b9a",
        stroke="#38006b",
        glyph="=",
        label="bit-identical to the true kernel — this channel cannot see the defect",
    ),
    ChannelState.EXACT: ChannelStyle(
        fill="#0277bd",
        stroke="#004c8c",
        glyph="0",
        label="exact — 0.0 ULP, agreement is bitwise",
    ),
    ChannelState.WITHIN_BUDGET: ChannelStyle(
        fill="#2e7d32",
        stroke="#005005",
        glyph="+",
        label="within budget — moved, and stayed inside its gate",
    ),
    ChannelState.ABSENT: ChannelStyle(
        fill="#ffffff",
        stroke="#8d6e63",
        glyph="?",
        label="absent — the artefact carries no such field for this case",
    ),
}


def log_position(ulp: float, *, low: int = DECADE_LOW, high: int = DECADE_HIGH) -> tuple[float, bool]:
    """Place a ULP figure on the log axis.

    Args:
        ulp: The measured disagreement, strictly positive.
        low: Bottom decade exponent.
        high: Top decade exponent.

    Returns:
        ``(fraction, clamped)`` -- the position in ``[0, 1]`` and whether it was clamped to an edge.
        A clamped mark is flagged so it is never read as a measured position.

    Raises:
        EvidenceError: If ``ulp`` is zero or negative. Zero is **exact**, which is a different kind of
            statement from small, and returning the axis origin for it would draw
            `ALEPH-PORT-3604`'s headline result as merely a very good one.
    """
    if ulp == 0.0:
        raise EvidenceError(
            "0.0 ULP is exact agreement, not a small disagreement; it belongs in EXACT_LANE and has "
            "no position on a log axis"
        )
    if not (ulp > 0.0) or not math.isfinite(ulp):
        raise EvidenceError(f"{ulp!r} is not a positive finite ULP figure")
    span = float(high - low)
    fraction = (math.log10(ulp) - low) / span
    if fraction < 0.0:
        return 0.0, True
    if fraction > 1.0:
        return 1.0, True
    return fraction, False


@dataclass(slots=True)
class FieldRow:
    """One field of one case: the unit this module judges on.

    Attributes:
        case: Full case name, e.g. ``chromatin_network[WRONG_ENERGY_HALF]``.
        family: The law family.
        variant: The kernel variant.
        field_name: The measured field, e.g. ``energy_pn_um``.
        is_wrong_on_purpose: Whether the variant is a planted wrong kernel.
        ulp_cuda: Disagreement on the CUDA device, in float32 ULP of the field scale.
        ulp_cpu: The same figure on the warp CPU device, or ``None`` if the case is absent there.
        budget: Effective ULP budget the field was judged against.
        roundoff_bound: The round-off bound the owning ledger derived for this field.
        pressure: ``ulp_cuda / budget``, or ``None`` when the budget is exactly zero.
        headroom: ``roundoff_bound / budget``, or ``None`` when the budget is exactly zero.
        within_budget: The artefact's own verdict for this field on CUDA.
        bit_identical_to_true: ``True``/``False``, or ``None`` when the family has no true variant to
            compare against -- which is *unknown*, not "not blind".
        state: The rendered channel state on CUDA.
        state_cpu: The same on the CPU device, for the cross-device comparison.
        scatter_assembled: Whether the field was assembled by a scatter, so the reader can see which
            marks the ORDERED/ATOMIC question applies to at all.
        shape: The field's array shape; ``[]`` for a scalar.
        absent_reason: Why the field is absent, when it is.
    """

    case: str
    family: str
    variant: str
    field_name: str
    is_wrong_on_purpose: bool
    ulp_cuda: float | None = None
    ulp_cpu: float | None = None
    budget: float | None = None
    roundoff_bound: float | None = None
    pressure: float | None = None
    headroom: float | None = None
    within_budget: bool | None = None
    bit_identical_to_true: bool | None = None
    state: ChannelState = ChannelState.ABSENT
    state_cpu: ChannelState = ChannelState.ABSENT
    scatter_assembled: bool = False
    shape: tuple[int, ...] = ()
    absent_reason: str = ""

    @property
    def is_exact(self) -> bool:
        """Whether the CUDA disagreement is exactly zero."""
        return self.ulp_cuda == 0.0

    @property
    def thin_headroom(self) -> bool:
        """Whether this field's gate sits less than a binary order under its round-off bound."""
        return self.headroom is not None and self.headroom < THIN_HEADROOM

    @property
    def moved_between_devices(self) -> bool:
        """Whether CPU and CUDA disagree about this field's own disagreement."""
        return (
            self.ulp_cpu is not None
            and self.ulp_cuda is not None
            and self.ulp_cpu != self.ulp_cuda
        )

    def label(self) -> str:
        """Row label: family, variant and field, in the order a reader scans them."""
        return f"{self.family}[{self.variant}]·{self.field_name}"


@dataclass(slots=True)
class CaseRow:
    """One law case, with its fields and both devices' gate verdicts.

    Attributes:
        name: Full case name.
        family: The law family.
        variant: The kernel variant.
        is_wrong_on_purpose: Whether this is a planted wrong kernel.
        caught_cuda: Whether the gate caught it on CUDA.
        caught_cpu: Whether the gate caught it on the warp CPU device.
        declared_blind: Whether the artefact **declares** this variant blind on both devices.
        description: The case's own description of the state it ran on.
        sites: Number of sites the case drove.
        fields: The per-field rows.
    """

    name: str
    family: str
    variant: str
    is_wrong_on_purpose: bool
    caught_cuda: bool
    caught_cpu: bool
    declared_blind: bool
    description: str = ""
    sites: int = 0
    fields: list[FieldRow] = field(default_factory=list)

    @property
    def blind_on_both(self) -> bool:
        """A wrong kernel no gate caught, on either device."""
        return self.is_wrong_on_purpose and not self.caught_cuda and not self.caught_cpu

    @property
    def undeclared_blind(self) -> bool:
        """Blind on both devices **and not declared** — which is a defect, unlike a declared one.

        `ALEPH-PORT-3604` §9a ships ``WRONG_NEWTON_MIRROR`` as a declared blind spot: it is
        bit-identical to the true kernel at every input by an IEEE-754 identity, so no gate can see
        it and none should. Collapsing that with an *undeclared* silent mutant would hide the only
        one of the two that is a finding.
        """
        return self.blind_on_both and not self.declared_blind

    @property
    def is_regression(self) -> bool:
        """A true kernel that failed its own gate."""
        return not self.is_wrong_on_purpose and any(
            row.state is ChannelState.OVER_BUDGET for row in self.fields
        )


@dataclass(slots=True)
class ParityEvidence:
    """The law-parity artefact, resolved into rows plus everything that could not be read.

    Attributes:
        cases: Every case, in artefact order.
        rows: Every field row, flattened.
        absent: ``(case, reason)`` for anything named by the artefact and not turned into a row.
        gate: The artefact's own ``gate_comparison`` block.
        recounted: The same five figures, recounted from ``cases``.
        device: The CUDA device block.
        cpu_key: The device key the CPU column was read from.
        cuda_key: The device key the CUDA column was read from.
        verdict: The artefact's own verdict.
        evidence_rung: The artefact's declared evidence rung.
        not_established: The artefact's own list of what it does not establish.
        finished_utc: When the run finished.
        host: The host the run happened on.
        position_precision: The position mode the run carried.
    """

    cases: list[CaseRow] = field(default_factory=list)
    rows: list[FieldRow] = field(default_factory=list)
    absent: list[tuple[str, str]] = field(default_factory=list)
    gate: dict[str, Any] = field(default_factory=dict)
    recounted: dict[str, Any] = field(default_factory=dict)
    device: dict[str, Any] = field(default_factory=dict)
    cpu_key: str = ""
    cuda_key: str = ""
    verdict: str = ""
    evidence_rung: str = ""
    not_established: list[str] = field(default_factory=list)
    finished_utc: str = ""
    host: str = ""
    position_precision: str = ""

    @property
    def wrong_kernels(self) -> list[CaseRow]:
        """Every planted wrong kernel."""
        return [case for case in self.cases if case.is_wrong_on_purpose]

    @property
    def true_kernels(self) -> list[CaseRow]:
        """Every true kernel."""
        return [case for case in self.cases if not case.is_wrong_on_purpose]

    @property
    def blind_rows(self) -> list[FieldRow]:
        """Every field a wrong kernel left bit-identical to the true kernel's."""
        return [row for row in self.rows if row.state is ChannelState.BLIND_TO_TRUE]

    @property
    def thin_budgets(self) -> list[FieldRow]:
        """Every field whose gate has less than a binary order of headroom."""
        return [row for row in self.rows if row.thin_headroom]

    @property
    def exact_rows(self) -> list[FieldRow]:
        """Every field that agreed bitwise."""
        return [row for row in self.rows if row.is_exact]

    @property
    def undeclared_blind(self) -> list[CaseRow]:
        """Every wrong kernel silent on both devices without being declared so."""
        return [case for case in self.cases if case.undeclared_blind]

    @property
    def regressions(self) -> list[CaseRow]:
        """Every true kernel that failed its own gate."""
        return [case for case in self.cases if case.is_regression]

    def faults(self) -> list[str]:
        """Everything a reader must not miss, as sentences."""
        out = [
            f"{case.name}: silent on BOTH devices and NOT declared blind"
            for case in self.undeclared_blind
        ]
        out += [f"{case.name}: true kernel over its own budget" for case in self.regressions]
        out += [f"{name}: {reason}" for name, reason in self.absent]
        return out


def _as_float(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _pick_device_keys(devices: Mapping[str, Any]) -> tuple[str, str]:
    """Choose the CPU and CUDA device keys, and refuse if they are the same column.

    A reader that pulled both columns from ``cuda:0`` would report "no cross-device difference"
    forever and look perfect -- the same failure shape as `ALEPH-PORT-3609` §14.4's zero-copy view,
    where a control compared a buffer with itself. This is where that is made impossible.
    """
    if not devices:
        raise EvidenceError("the artefact names no devices, so there is nothing to compare")
    cuda_keys = [
        key
        for key, block in devices.items()
        if isinstance(block, Mapping)
        and (bool(block.get("is_cuda")) or "cuda" in str(block.get("warp_device", key)).lower())
    ]
    cpu_keys = [key for key in devices if key not in cuda_keys]
    if not cuda_keys:
        raise EvidenceError(
            f"none of the devices {sorted(devices)} is a CUDA device; this figure renders a "
            "CPU-vs-CUDA comparison and there is no CUDA column to draw"
        )
    if not cpu_keys:
        raise EvidenceError(
            f"the artefact carries only CUDA devices {sorted(devices)}; with no CPU column the "
            "cross-device comparison would compare a column with itself"
        )
    cuda_key, cpu_key = cuda_keys[0], cpu_keys[0]
    if cuda_key == cpu_key:
        raise EvidenceError(
            f"the CPU and CUDA columns resolved to the same key {cuda_key!r}; a comparison of a "
            "column with itself is silent by construction"
        )
    if devices[cuda_key] is devices[cpu_key]:
        raise EvidenceError(
            f"devices {cuda_key!r} and {cpu_key!r} are the same object, so the cross-device "
            "comparison would be vacuous"
        )
    return cpu_key, cuda_key


def _field_map(case: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    out: dict[str, Mapping[str, Any]] = {}
    for entry in case.get("fields") or []:
        if isinstance(entry, Mapping) and entry.get("field") is not None:
            out[str(entry["field"])] = entry
    return out


def _state_for(
    *,
    absent: bool,
    within_budget: bool | None,
    blind: bool | None,
    is_wrong: bool,
    ulp: float | None,
) -> ChannelState:
    """Precedence: a fault first, then the finding, then exactness, then the ordinary pass."""
    if absent:
        return ChannelState.ABSENT
    if within_budget is False:
        return ChannelState.OVER_BUDGET
    if is_wrong and blind is True:
        return ChannelState.BLIND_TO_TRUE
    if ulp == 0.0:
        return ChannelState.EXACT
    return ChannelState.WITHIN_BUDGET


def read_parity(payload: Mapping[str, Any]) -> ParityEvidence:
    """Resolve a ``j1_law_parity`` artefact into per-channel rows.

    Args:
        payload: The parsed artefact.

    Returns:
        The :class:`ParityEvidence`.

    Raises:
        EvidenceError: If the artefact is not this shape, if the two device columns cannot be
            distinguished, or if **recounting the gate comparison from the rows disagrees with the
            artefact's own block**. That last one is the strongest control here: it fails if the
            reader drops a case, double-counts one, or reads one device twice.
    """
    measurements = payload.get("measurements")
    if not isinstance(measurements, Mapping):
        raise EvidenceError("the artefact carries no 'measurements' block")
    devices = measurements.get("devices")
    if not isinstance(devices, Mapping):
        raise EvidenceError("the artefact's measurements carry no 'devices' block")

    cpu_key, cuda_key = _pick_device_keys(devices)
    cuda_cases = devices[cuda_key].get("cases") or {}
    cpu_cases = devices[cpu_key].get("cases") or {}
    if not isinstance(cuda_cases, Mapping) or not cuda_cases:
        raise EvidenceError(
            f"device {cuda_key!r} carries no cases, so there is no evidence to draw. An empty "
            "matrix is a refusal, not an empty figure."
        )

    gate = dict(measurements.get("gate_comparison") or {})
    declared_blind = {str(name) for name in (gate.get("blind_on_both_devices") or [])}

    evidence = ParityEvidence(
        gate=gate,
        device=dict(payload.get("device") or {}),
        cpu_key=cpu_key,
        cuda_key=cuda_key,
        verdict=str(measurements.get("verdict", "")),
        evidence_rung=str(payload.get("evidence_rung", "")),
        not_established=[str(item) for item in (payload.get("not_established") or [])],
        finished_utc=str(payload.get("finished_utc", "")),
        host=str((payload.get("environment") or {}).get("hostname", "")),
        position_precision=str(measurements.get("position_precision", "")),
    )

    # The true kernel of each family is the baseline every blind verdict is measured against. A family
    # with none yields `None` -- unknown -- and never `False`, because "not bit-identical to nothing"
    # is not a statement anybody can check.
    # One baseline per device. Comparing the CPU column against the CUDA baseline would be a
    # cross-device comparison wearing a same-device label: the two devices disagree slightly even for
    # the true kernel, so every CPU field would read "not blind" and the finding would vanish from
    # half the figure.
    def _baselines(source: Mapping[str, Any]) -> dict[str, dict[str, float]]:
        out: dict[str, dict[str, float]] = {}
        for case_name, case in source.items():
            if not isinstance(case, Mapping) or case.get("is_wrong_on_purpose"):
                continue
            out[str(case.get("family", case_name))] = {
                key: value
                for key, entry in _field_map(case).items()
                if (value := _as_float(entry.get("ulp"))) is not None
            }
        return out

    true_ulp = _baselines(cuda_cases)
    true_ulp_cpu = _baselines(cpu_cases) if isinstance(cpu_cases, Mapping) else {}

    # Every field name the family uses anywhere, so a field present for one variant and missing for
    # another is rendered ABSENT rather than quietly dropped from that row.
    family_fields: dict[str, set[str]] = {}
    for name, case in cuda_cases.items():
        if not isinstance(case, Mapping):
            continue
        family_fields.setdefault(str(case.get("family", name)), set()).update(_field_map(case))

    for name, case in cuda_cases.items():
        if not isinstance(case, Mapping):
            evidence.absent.append((str(name), "the case is not an object and could not be read"))
            continue
        family = str(case.get("family", name))
        variant = str(case.get("variant", "?"))
        is_wrong = bool(case.get("is_wrong_on_purpose"))
        cpu_case = cpu_cases.get(name) if isinstance(cpu_cases, Mapping) else None
        cpu_fields = _field_map(cpu_case) if isinstance(cpu_case, Mapping) else {}
        row = CaseRow(
            name=str(name),
            family=family,
            variant=variant,
            is_wrong_on_purpose=is_wrong,
            caught_cuda=bool(case.get("caught")),
            caught_cpu=bool(cpu_case.get("caught")) if isinstance(cpu_case, Mapping) else False,
            declared_blind=str(name) in declared_blind,
            description=str(case.get("description", "")),
            sites=int(case.get("sites") or 0),
        )
        if not isinstance(cpu_case, Mapping):
            evidence.absent.append(
                (str(name), f"present on {cuda_key} and absent from {cpu_key}, so it has no "
                 "cross-device comparison")
            )

        cuda_fields = _field_map(case)
        baseline = true_ulp.get(family)
        baseline_cpu = true_ulp_cpu.get(family)
        for field_name in sorted(family_fields.get(family, set())):
            entry = cuda_fields.get(field_name)
            if entry is None:
                row.fields.append(
                    FieldRow(
                        case=str(name),
                        family=family,
                        variant=variant,
                        field_name=field_name,
                        is_wrong_on_purpose=is_wrong,
                        state=ChannelState.ABSENT,
                        state_cpu=ChannelState.ABSENT,
                        absent_reason=(
                            f"other variants of {family} report {field_name} and this case does not"
                        ),
                    )
                )
                continue
            ulp = _as_float(entry.get("ulp"))
            budget = _as_float(entry.get("effective_ulp_budget"))
            roundoff = _as_float(entry.get("roundoff_ulp_bound"))
            within = entry.get("within_budget")
            within = bool(within) if isinstance(within, bool) else None
            cpu_entry = cpu_fields.get(field_name)
            ulp_cpu = _as_float(cpu_entry.get("ulp")) if isinstance(cpu_entry, Mapping) else None
            cpu_within = cpu_entry.get("within_budget") if isinstance(cpu_entry, Mapping) else None
            cpu_within = bool(cpu_within) if isinstance(cpu_within, bool) else None

            blind: bool | None = None
            blind_cpu: bool | None = None
            if baseline is not None and field_name in baseline and ulp is not None:
                # `==` on purpose: exactness is the evidence, and `isclose` would turn an exact
                # statement into an approximate one. See `ALEPH-PORT-3612` §10.
                blind = ulp == baseline[field_name]
            if baseline_cpu is not None and field_name in baseline_cpu and ulp_cpu is not None:
                blind_cpu = ulp_cpu == baseline_cpu[field_name]

            pressure = (
                ulp / budget if ulp is not None and budget not in (None, 0.0) else None
            )
            headroom = (
                roundoff / budget if roundoff is not None and budget not in (None, 0.0) else None
            )
            # The artefact records its own pressure. If ours disagrees we have paired a ULP with the
            # wrong budget, which is a reader defect and not a rounding difference.
            declared_pressure = _as_float(entry.get("pressure"))
            if (
                pressure is not None
                and declared_pressure is not None
                and not math.isclose(pressure, declared_pressure, rel_tol=1.0e-9, abs_tol=0.0)
            ):
                raise EvidenceError(
                    f"{name}·{field_name}: recomputed pressure {pressure!r} disagrees with the "
                    f"artefact's {declared_pressure!r}; the ULP has been paired with the wrong budget"
                )

            shape = entry.get("shape") or []
            row.fields.append(
                FieldRow(
                    case=str(name),
                    family=family,
                    variant=variant,
                    field_name=field_name,
                    is_wrong_on_purpose=is_wrong,
                    ulp_cuda=ulp,
                    ulp_cpu=ulp_cpu,
                    budget=budget,
                    roundoff_bound=roundoff,
                    pressure=pressure,
                    headroom=headroom,
                    within_budget=within,
                    bit_identical_to_true=blind,
                    state=_state_for(
                        absent=False, within_budget=within, blind=blind, is_wrong=is_wrong, ulp=ulp
                    ),
                    state_cpu=_state_for(
                        absent=not isinstance(cpu_entry, Mapping),
                        within_budget=cpu_within,
                        blind=blind_cpu,
                        is_wrong=is_wrong,
                        ulp=ulp_cpu,
                    ),
                    scatter_assembled=bool(entry.get("scatter_assembled")),
                    shape=tuple(int(value) for value in shape),
                )
            )
        evidence.cases.append(row)
        evidence.rows.extend(row.fields)

    evidence.recounted = _recount(evidence)
    _check_recount(evidence)
    return evidence


def _recount(evidence: ParityEvidence) -> dict[str, Any]:
    """Recount the gate comparison from the rows that were actually built."""
    wrong = evidence.wrong_kernels
    caught_cuda = sorted(case.name for case in wrong if case.caught_cuda)
    caught_cpu = sorted(case.name for case in wrong if case.caught_cpu)
    return {
        "wrong_kernels_total": len(wrong),
        "caught_on_cuda": len(caught_cuda),
        "caught_on_cpu": len(caught_cpu),
        "caught_on_cuda_but_not_cpu": sorted(set(caught_cuda) - set(caught_cpu)),
        "caught_on_cpu_but_not_cuda": sorted(set(caught_cpu) - set(caught_cuda)),
        "blind_on_both_devices": sorted(case.name for case in wrong if case.blind_on_both),
    }


def _check_recount(evidence: ParityEvidence) -> None:
    """Require the recount to equal the artefact's own block, on every key the artefact carries."""
    gate = evidence.gate
    if not gate:
        raise EvidenceError(
            "the artefact carries no gate_comparison block, so the recount has nothing to be "
            "checked against and the figure would be asserting on itself"
        )
    for key in (
        "wrong_kernels_total",
        "caught_on_cuda",
        "caught_on_cpu",
        "caught_on_cuda_but_not_cpu",
        "caught_on_cpu_but_not_cuda",
        "blind_on_both_devices",
    ):
        if key not in gate:
            continue
        declared = gate[key]
        found = evidence.recounted[key]
        if isinstance(declared, list):
            declared = sorted(str(item) for item in declared)
        if declared != found:
            raise EvidenceError(
                f"recounting {key} from the {len(evidence.cases)} cases drawn gives {found!r}, and "
                f"the artefact declares {declared!r}. The figure would not be describing the run."
            )


@dataclass(slots=True)
class ScatterEvidence:
    """The scatter-determinism artefact: whether a repeated reduction is bit-identical.

    Attributes:
        modes: One entry per measured mode, each carrying its spread and verdict.
        verdict: The artefact's own verdict.
    """

    modes: list[dict[str, Any]] = field(default_factory=list)
    verdict: str = ""

    @property
    def ordered(self) -> dict[str, Any] | None:
        """The ORDERED mode's row, if measured."""
        return next((row for row in self.modes if row.get("mode") == "ordered"), None)

    @property
    def atomic(self) -> dict[str, Any] | None:
        """The ATOMIC mode's row, if measured."""
        return next((row for row in self.modes if row.get("mode") == "atomic"), None)


def read_scatter(payload: Mapping[str, Any]) -> ScatterEvidence:
    """Resolve a ``j2_scatter_determinism`` artefact.

    Args:
        payload: The parsed artefact.

    Returns:
        The :class:`ScatterEvidence`.

    Raises:
        EvidenceError: If the artefact carries no measurements.
    """
    measurements = payload.get("measurements")
    if not isinstance(measurements, Mapping):
        raise EvidenceError("the scatter artefact carries no 'measurements' block")
    rows = measurements.get("measurements")
    if not isinstance(rows, Sequence) or not rows:
        raise EvidenceError("the scatter artefact records no repeats, so there is nothing to draw")
    modes: list[dict[str, Any]] = []
    for entry in rows:
        if not isinstance(entry, Mapping):
            continue
        if str(entry.get("device", "")).startswith("cuda") or bool(entry.get("is_cuda")):
            modes.append(dict(entry))
    if not modes:
        modes = [dict(entry) for entry in rows if isinstance(entry, Mapping)]
    return ScatterEvidence(modes=modes, verdict=str(measurements.get("verdict", "")))


@dataclass(slots=True)
class DescentWindow:
    """The float32 descent window, measured on the device and compared to the prior CPU figure.

    `ALEPH-PORT-3609` §14.6 raised the disagreement between this measurement and `-3608`'s ratified
    derivation, named three candidate repairs and **took none**. This carries both windows so a
    figure can show a disagreement rather than pick a side.

    Attributes:
        per_config: Converging slack multiples, per configuration.
        intersection: Multiples that converge on every configuration.
        prior_cpu: The prior warp-CPU-device window.
        prior_source: Where the prior figure came from.
        rows: Every measured row.
        max_steps: The step ceiling the sweep ran under.
    """

    per_config: dict[str, list[float]] = field(default_factory=dict)
    intersection: list[float] = field(default_factory=list)
    prior_cpu: list[float] = field(default_factory=list)
    prior_source: str = ""
    rows: list[dict[str, Any]] = field(default_factory=list)
    max_steps: int = 0

    @property
    def disagrees_with_prior(self) -> bool:
        """Whether the device window and the prior CPU window share no multiple."""
        return bool(self.intersection) and bool(self.prior_cpu) and not (
            set(self.intersection) & set(self.prior_cpu)
        )

    @property
    def multiples(self) -> list[float]:
        """Every slack multiple that was probed, in order."""
        seen = {
            value
            for row in self.rows
            if (value := _as_float(row.get("slack_multiple_of_u32"))) is not None
        }
        return sorted(seen)

    @property
    def has_margin(self) -> bool:
        """Whether the window is wider than a single probed multiple.

        A window one multiple wide has **no margin on either side**: the multiple below stalls and
        the one above wanders, and a lane inheriting the number rather than re-measuring it is
        relying on a coincidence.
        """
        return len(self.intersection) > 1


def read_descent_window(payload: Mapping[str, Any]) -> DescentWindow:
    """Resolve a ``g6_descent_window`` artefact.

    Args:
        payload: The parsed artefact.

    Returns:
        The :class:`DescentWindow`.

    Raises:
        EvidenceError: If the artefact records no rows.
    """
    rows = payload.get("rows")
    if not isinstance(rows, Sequence) or not rows:
        raise EvidenceError("the descent-window artefact records no rows, so there is no window")
    per_config = {
        str(key): [float(value) for value in values]
        for key, values in (payload.get("converging_multiples_per_config") or {}).items()
    }
    return DescentWindow(
        per_config=per_config,
        intersection=[
            float(value) for value in (payload.get("converging_multiples_on_both_configs") or [])
        ],
        prior_cpu=[
            float(value) for value in (payload.get("prior_cpu_device_window_multiples") or [])
        ],
        prior_source=str(payload.get("prior_source", "")),
        rows=[dict(row) for row in rows if isinstance(row, Mapping)],
        max_steps=int(payload.get("max_steps") or 0),
    )


@dataclass(slots=True)
class EvidenceBundle:
    """Everything found, plus everything that was looked for and was not there.

    ``missing`` is the point. A panel whose artefact is absent is drawn **as absent** with the reason,
    following `render.py`'s rule that an owner which cannot be drawn is listed rather than skipped.

    Attributes:
        parity: The law-parity evidence. Required.
        scatter: The scatter-determinism evidence, when present.
        descent: The descent-window evidence, when present.
        summary: The runner's own summary block, when present.
        source: Where the artefacts were read from.
        missing: ``(name, reason)`` for every artefact looked for and not read.
    """

    parity: ParityEvidence
    scatter: ScatterEvidence | None = None
    descent: DescentWindow | None = None
    summary: dict[str, Any] = field(default_factory=dict)
    source: str = ""
    missing: list[tuple[str, str]] = field(default_factory=list)


def _fetch_over_ssh(host: str, destination: Path) -> tuple[bool, str]:
    """Copy the runner's artefact directory off ``host``. Returns ``(ok, detail)``."""
    if shutil.which("scp") is None:
        return False, "scp is not on PATH, so the remote artefacts cannot be fetched"
    destination.mkdir(parents=True, exist_ok=True)
    try:
        result = subprocess.run(  # noqa: S603 - fixed argv, no shell
            ["scp", "-q", "-o", "ConnectTimeout=15", f"{host}:~/.aleph_runner/artefacts/*.json",
             str(destination)],
            capture_output=True,
            text=True,
            timeout=180,
        )
    except (OSError, subprocess.SubprocessError) as error:
        return False, f"scp from {host} failed: {error}"
    if result.returncode != 0:
        return False, f"scp from {host} exited {result.returncode}: {result.stderr.strip()[:200]}"
    return True, f"fetched from {host}:~/.aleph_runner/artefacts/"


def discover_artefacts(
    *,
    directory: Path | str | None = None,
    host: str | None = DEFAULT_HOST,
    refresh: bool = False,
) -> tuple[Path, str]:
    """Find the artefact directory, fetching it from the runner's host if it is not here.

    Args:
        directory: An explicit directory. When given it is used as-is and never fetched into.
        host: ssh alias to fall back to; ``None`` disables the fallback.
        refresh: Re-fetch even when a local copy exists.

    Returns:
        ``(directory, provenance)``.

    Raises:
        EvidenceError: If no artefact directory can be found or fetched. The message names what was
            tried, because "not found" without the search is not actionable.
    """
    if directory is not None:
        target = Path(directory).expanduser()
        if not target.is_dir():
            raise EvidenceError(f"{target} is not a directory")
        if not (target / PARITY_ARTEFACT).is_file():
            raise EvidenceError(
                f"{target} holds no {PARITY_ARTEFACT}; that artefact is what this figure renders. "
                f"Present: {', '.join(sorted(p.name for p in target.glob('*.json'))) or 'nothing'}"
            )
        return target, f"local {target}"

    local = DEFAULT_ARTEFACT_DIR
    if not refresh and (local / PARITY_ARTEFACT).is_file():
        return local, f"local {local}"
    tried = [f"{local} (no {PARITY_ARTEFACT})"]
    if host:
        ok, detail = _fetch_over_ssh(host, local)
        if ok and (local / PARITY_ARTEFACT).is_file():
            return local, detail
        tried.append(detail if not ok else f"{host} carries no {PARITY_ARTEFACT}")
    raise EvidenceError(
        "no CUDA run artefacts found. Tried: " + "; ".join(tried) + ". Point --from at a directory "
        f"holding {PARITY_ARTEFACT}, or make `ssh {host or DEFAULT_HOST}` work."
    )


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise EvidenceError(f"{path} could not be read as JSON: {error}") from error


def load_bundle(
    directory: Path | str | None = None,
    *,
    host: str | None = DEFAULT_HOST,
    refresh: bool = False,
) -> EvidenceBundle:
    """Load every artefact this figure can draw, and name the ones that are not there.

    Args:
        directory: An explicit artefact directory, or ``None`` to discover one.
        host: ssh alias for the fallback fetch.
        refresh: Re-fetch even when a local copy exists.

    Returns:
        The :class:`EvidenceBundle`.

    Raises:
        EvidenceError: If the parity artefact is absent or unreadable. The optional ones never raise;
            they are recorded in ``missing`` and drawn as absent.
    """
    target, provenance = discover_artefacts(directory=directory, host=host, refresh=refresh)
    bundle = EvidenceBundle(
        parity=read_parity(_read_json(target / PARITY_ARTEFACT)), source=provenance
    )
    readers = {
        "j2_scatter_determinism.json": ("scatter", read_scatter),
        "g6_descent_window.json": ("descent", read_descent_window),
    }
    for name in OPTIONAL_ARTEFACTS:
        path = target / name
        if not path.is_file():
            bundle.missing.append((name, f"not present in {target}"))
            continue
        if name == "SUMMARY.json":
            payload = _read_json(path)
            bundle.summary = dict(payload) if isinstance(payload, Mapping) else {}
            continue
        attribute, reader = readers[name]
        try:
            setattr(bundle, attribute, reader(_read_json(path)))
        except EvidenceError as error:
            bundle.missing.append((name, str(error)))
    return bundle

"""Resting osmotic turgor Π₀ — the REQUIRED, no-silent-default parameter gate.

PI decision 2026-07-25: **gate Π₀**; label (but do not gate) cortex ``seg_um`` / ``length_um`` /
``density_per_fil``. Ledger: :data:`LEDGER_PATH` (``aleph/laws/params_turgor.yaml``), whose format
mirrors ``ac/motor/params_i0b3.yaml`` exactly (``value: null`` + ``evidence_status: "GAP — PI"`` +
``owner`` + ``expiry_trigger`` + ``demonstrated_by``).

Why this module exists
----------------------
``PARAM_PROVENANCE_AUDIT_2026-07-24.md`` row 1: Π₀ has **three** competing in-tree values — 40 Pa
(HeLa proxy, ``ff/gamma_floor.py``/``ac/cell/assemble.py``/``CLAUDE.md``), 133 Pa (band-implied ⇒
circular, ``configs/mcf7_baseline.yaml``/``common/production_policy.py``/``dcm/``), and ~72 Pa
(48–107, MCF7-geometry estimate, unwired) — and **none is an MCF7 measurement**. Π₀ sets the entire
resting cortical pre-tension via Young-Laplace γ = ΔP·R/2, so a silent choice silently decides the
baseline every measurement is taken from.

What the gate does and does NOT do
----------------------------------
DOES: refuse to supply Π₀ from thin air. :func:`resolve_pi0_pa` with no value **raises**
:class:`TurgorPi0Unspecified` with a message naming all three competing values and pointing at the
audit. Every accepted value is returned wrapped in a :class:`Pi0Resolution` that classifies its
provenance, so run artifacts record the value *and* that it is CONVENIENCE.

DOES NOT: make any Π₀ correct, select one, or change one. ``resolve_pi0_pa(40.0)`` still runs — and
still stamps ``provenance="CONVENIENCE"`` on the artifact. This is claim hygiene, not physics.

Sanity Gate
-----------
* Dimensional: Π₀ in Pa; engine pressure units are pN/µm² and 1 Pa == 1 pN/µm² exactly, so the
  value is unit-identical in both conventions (no conversion, asserted by test).
* Boundary: ``None`` ⇒ raise (the gate). ``0.0`` ⇒ accepted but classified
  ``UNPHYSICAL_ZERO_BASELINE`` — the physiological-baseline HARD rule forbids a force-free shell as a
  production baseline, so the classification is loud rather than silently permissive.
  Negative Π₀ ⇒ raise (a resting cell does not suck inward).
* Sign-sense: positive Π₀ = net OUTWARD (inflating) pressure, matching
  ``gamma_floor._turgor_force`` (radial outward ΔP·area·r̂).
* Invariant: the ledger's gated value stays ``null``; a test asserts it, so this gate cannot be
  quietly converted back into a default by editing the YAML.
* Measurement-protocol: the returned record is what artifacts must serialize; a run that cannot name
  its Π₀ provenance is a run whose resting baseline is unsourced.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

#: The Π₀ ledger (single source of truth for the competing claims).
LEDGER_PATH: Path = Path(__file__).with_name("params_turgor.yaml")

#: The provenance audit that motivated the gate — quoted in the raise message.
AUDIT_DOC: str = "aleph/docs/v2_audit/PARAM_PROVENANCE_AUDIT_2026-07-24.md"

# ── Named claims (values are RECORDS OF WHAT IS IN THE TREE, not recommendations) ────────────────
#: 40 Pa — Fischer-Friedrich 2014 HeLa interphase. CONVENIENCE: measured, wrong cell line.
PI0_CLAIM_HELA_PROXY_40PA: float = 40.0
#: 133 Pa — Young-Laplace inversion of the MCF7 cortical-tension band. Circular by construction.
PI0_CLAIM_BAND_IMPLIED_133PA: float = 133.0
#: ~72 Pa — MCF7-geometry central estimate (48–107) from the provenance audit. Not a measurement.
PI0_CLAIM_MCF7_GEOMETRY_72PA: float = 72.0

_CLAIM_IDS: dict[float, str] = {
    PI0_CLAIM_HELA_PROXY_40PA: "claim_hela_proxy_40pa",
    PI0_CLAIM_BAND_IMPLIED_133PA: "claim_band_implied_133pa",
    PI0_CLAIM_MCF7_GEOMETRY_72PA: "claim_mcf7_geometry_72pa",
}


class TurgorPi0Unspecified(ValueError):
    """Raised when Π₀ is needed but no value was passed (the gate; there is no default)."""


class TurgorPi0Invalid(ValueError):
    """Raised when a Π₀ value cannot be a resting turgor at all (negative / non-finite)."""


@dataclass(frozen=True, slots=True)
class Pi0Resolution:
    """One explicit, provenance-classified Π₀ choice — the thing artifacts must record."""

    value_pa: float
    provenance: str
    claim_id: str
    source: str
    caller: str
    note: str = ""

    @property
    def is_sourced_for_mcf7(self) -> bool:
        """True only if this Π₀ is an MCF7 measurement. Currently False for every in-tree claim."""
        return self.provenance == "SOURCED_MCF7"

    def as_artifact_record(self) -> dict[str, Any]:
        """Serializable provenance record for a run artifact / REPORT.md table."""
        return {
            "turgor_PI_0_Pa": self.value_pa,
            "turgor_PI_0_provenance": self.provenance,
            "turgor_PI_0_claim_id": self.claim_id,
            "turgor_PI_0_source": self.source,
            "turgor_PI_0_sourced_for_MCF7": self.is_sourced_for_mcf7,
            "turgor_PI_0_resolved_by": self.caller,
            "turgor_PI_0_ledger": "aleph/laws/params_turgor.yaml",
            "turgor_PI_0_audit": AUDIT_DOC,
            "turgor_PI_0_note": self.note,
        }


@functools.lru_cache(maxsize=1)
def load_ledger() -> dict[str, Any]:
    """Parse ``params_turgor.yaml`` (cached)."""
    with LEDGER_PATH.open(encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def ledger_gated_value() -> float | None:
    """The ledger's Π₀ ``value`` — ``None`` while the parameter is a GAP (the whole point)."""
    return load_ledger()["parameters"]["PI_0"]["value"]


def competing_claims() -> list[dict[str, Any]]:
    """The three in-tree Π₀ claims, as recorded in the ledger."""
    row = load_ledger()["parameters"]["PI_0"]
    return [row[k] | {"claim_id": k} for k in sorted(row) if k.startswith("claim_")]


def _claim_summary() -> str:
    lines = []
    for c in competing_claims():
        sites = ", ".join(c.get("sites") or ["(no call sites)"])
        lines.append(
            f"    - {c['value']:g} Pa  [{c['provenance_class']}]  {c['source']}\n"
            f"        sites: {sites}"
        )
    return "\n".join(lines)


def unspecified_message(caller: str) -> str:
    """The informative refusal message (also used by tests, so it cannot rot silently)."""
    return (
        f"Resting osmotic turgor Π₀ is REQUIRED and has NO default (caller: {caller}).\n"
        f"Π₀ sets the ENTIRE resting cortical tension via Young-Laplace gamma = dP*R/2, so it may not be\n"
        f"chosen silently. Pass it explicitly, e.g. resolve_pi0_pa(40.0, caller=...).\n"
        f"THREE competing values are in the tree and NONE is an MCF7 measurement:\n"
        f"{_claim_summary()}\n"
        f"    They span 3.3x — the resting pre-tension is undetermined by ~3x before any physics runs.\n"
        f"  Ledger : aleph/laws/params_turgor.yaml  (value: null, evidence_status 'GAP — PI', owner PI)\n"
        f"  Audit  : {AUDIT_DOC}  (row 1: 'the worst convenience default')\n"
        f"  PI must select Π₀ (assay + cell line) or accept a proxy EXPLICITLY; do not pick the one that\n"
        f"  makes a gate pass."
    )


def resolve_pi0_pa(
    value: float | None = None,
    *,
    claim: str | None = None,
    caller: str = "unspecified",
) -> Pi0Resolution:
    """Resolve an explicit Π₀ [Pa] into a provenance-classified record, or refuse.

    Args:
        value: the Π₀ the caller chose, in Pa (== pN/µm²). ``None`` raises — there is no default and
            the ledger deliberately holds ``value: null``.
        claim: optional ledger claim id (e.g. ``"claim_hela_proxy_40pa"``). Omitted ⇒ inferred by
            exact value match; a value matching no recorded claim is classified
            ``UNREGISTERED`` (accepted, but flagged in the artifact — a sweep point or a new
            selection PI has not yet registered).
        caller: free-text origin (module/function) recorded in the artifact.

    Returns:
        :class:`Pi0Resolution` — value + provenance class + claim id + source + caller.

    Raises:
        TurgorPi0Unspecified: ``value is None`` (and the ledger value is still ``null``).
        TurgorPi0Invalid: ``value`` is negative or non-finite.
    """
    if value is None:
        value = ledger_gated_value()          # still None while Π₀ is a GAP — PI
    if value is None:
        raise TurgorPi0Unspecified(unspecified_message(caller))

    value = float(value)
    if value != value or value in (float("inf"), float("-inf")):
        raise TurgorPi0Invalid(f"Π₀ must be finite, got {value!r} (caller: {caller})")
    if value < 0.0:
        raise TurgorPi0Invalid(
            f"Π₀ must be >= 0 (a resting cell does not suck inward); got {value!r} Pa "
            f"(caller: {caller}). Sign convention: positive Π₀ = net OUTWARD pressure."
        )

    claim_id = claim if claim is not None else _CLAIM_IDS.get(value)
    if value == 0.0:
        return Pi0Resolution(
            value_pa=0.0, provenance="UNPHYSICAL_ZERO_BASELINE", claim_id="none",
            source="no source — a force-free unpressurised shell",
            caller=caller,
            note=("Π₀=0 violates the physiological-baseline HARD rule for a production run (the "
                  "2026-06-04 saturation diagnosis came from exactly this). Accepted only as an "
                  "explicit turgor-OFF diagnostic."),
        )
    if claim_id is None:
        return Pi0Resolution(
            value_pa=value, provenance="UNREGISTERED", claim_id="none",
            source=f"not one of the ledger claims ({', '.join(f'{v:g}' for v in _CLAIM_IDS)} Pa)",
            caller=caller,
            note=("Value is not a registered claim — a sweep point, or a new selection that must be "
                  "registered in aleph/laws/params_turgor.yaml before it is used in a headline."),
        )

    row = load_ledger()["parameters"]["PI_0"].get(claim_id)
    if row is None:
        raise TurgorPi0Invalid(
            f"claim id {claim_id!r} is not in aleph/laws/params_turgor.yaml (caller: {caller})"
        )
    if float(row["value"]) != value:
        raise TurgorPi0Invalid(
            f"Π₀={value:g} Pa does not match ledger claim {claim_id!r} "
            f"(= {float(row['value']):g} Pa) (caller: {caller})"
        )
    return Pi0Resolution(
        value_pa=value,
        provenance=str(row["provenance_class"]).split()[0],
        claim_id=claim_id,
        source=str(row["source"]),
        caller=caller,
        note=str(row.get("note", "")).strip(),
    )


# ── CONVENIENCE labels (LABELED, NOT GATED — PI 2026-07-25) ─────────────────────────────────────

def convenience_label(name: str) -> dict[str, Any]:
    """The machine-readable CONVENIENCE provenance record for a labeled-not-gated parameter.

    Args:
        name: ledger key under ``convenience_labels`` (``cortex_seg_um`` / ``cortex_length_um`` /
            ``cortex_density_per_fil``).

    Returns:
        The ledger row (value, sourced value/range, deviation factor, sites, owner, expiry_trigger).

    Raises:
        KeyError: unknown label name.
    """
    labels = load_ledger()["convenience_labels"]
    if name not in labels:
        raise KeyError(f"{name!r} is not a labeled convenience parameter; have {sorted(labels)}")
    return labels[name] | {"label_id": name}


def convenience_labels_record() -> dict[str, Any]:
    """All CONVENIENCE labels, flattened for a run artifact.

    Keys are ``<param>_provenance`` / ``<param>_deviation`` so a REPORT.md or JSON artifact carries
    the deviation next to the value it used.
    """
    out: dict[str, Any] = {}
    for name, row in load_ledger()["convenience_labels"].items():
        out[f"{name}_value"] = row["value"]
        out[f"{name}_provenance"] = row["provenance_class"]
        out[f"{name}_deviation"] = row["deviation_factor"]
        out[f"{name}_gated"] = bool(row["gated"])
    return out

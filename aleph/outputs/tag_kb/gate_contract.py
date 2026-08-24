#!/usr/bin/env python3
"""Versioned gate contracts — reference implementation of the format + the content hash.

WHY THIS EXISTS (audit 2026-07-25, failure mode M2). The project's worst *recurring*
manoeuvre is a gate that changes between a FAIL and its superseding PASS. Four real
instances:

  1. H.2 rebanding — ``angle_ks_p_min: 1.0e-6  # widened from 0.05`` plus
     ``equipartition_rel_tol_3d: 0.60``, and the *already-failing* trajectory was
     "verified to PASS all 3 rebanded gates without re-running" (commit ``444ef010``).
  2. KU-5.5 ``xi_v`` band 50–200 widened to ``[1, 400]`` µm (commit ``9d40c09``).
  3. A "closes by construction" gate substitution the project itself later labelled,
     in writing, a "gate substitution"
     (``AC_ENGINE_INTEGRATION_PLAN_2026-07-22.md`` §6).
  4. GATE A closed after BOTH the OBSERVABLE changed (raw ``max|F|`` → projected
     ``max|PF| = F − Jᵀλ``) AND the CONFIGURATION changed (the discrete resting-myosin
     seed was removed from the static baseline).

Cases 3 and 4 are the important ones: the half that moved was **not the threshold**, it
was the *observable* and the *configuration*. A guard that hashes only thresholds would
have missed both. So a contract here pins all three, separately hashed:

    gate_id + THRESHOLD + OBSERVABLE (exact metric definition) + CONFIGURATION/population

The only guard that existed before this file was a sentence in CLAUDE.md plus PI
sign-off — i.e. cultural. This module is the machine half.

WHAT IS AND IS NOT HASHED
  hashed:      ``gate_id``, ``contract_version``, ``configuration_fixed``, ``claims``
               (every claim's ``observable`` and ``threshold`` subtree included)
  not hashed:  ``title``, ``status``, ``enforcement``, ``evidence``,
               ``configuration_recorded`` (per-run knobs — these are *diffed* between a
               FAIL and its superseding PASS instead), and any key named ``notes`` or
               ending in ``_note`` / ``_notes`` (prose, so an editorial fix does not
               invalidate every stamp).

  ``configuration_recorded`` is NOT a loophole: a run-knob that moves between a FAIL and
  a PASS is reported as ``CONFIG_MOVED_AT_PASS`` and needs a disclosed ledger row — it
  just does not need a PI *signature*, because "we ran the solver 20× longer" is a
  disclosure, not a contract change. Thresholds, observables, and fixed configuration DO
  need a PI signature.

USED BY
  * ``verify_gate_contracts.py`` — the verifier / CI gate.
  * gate drivers, which stamp the hash into their own run artifact, e.g.
    ``aleph/scripts/ac_gate_b_sf_motor_native.py`` (stamps BEFORE any GPU work, so a
    missing/renamed contract fails in milliseconds instead of after a 46-minute run).

Pure stdlib + PyYAML. No Notion token, no network, no duckdb, no CUDA.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import yaml

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]  # tag_kb -> outputs -> ffn_sim -> repo root
CONTRACT_DIR = REPO_ROOT / "aleph" / "docs" / "v2_audit" / "gate_contracts"
CHANGES_FILE = CONTRACT_DIR / "contract_changes.yaml"

#: hash-format tag. Bump when the canonicalisation changes (that invalidates every
#: recorded stamp on purpose — a silent canonicalisation change would be a way to
#: launder a contract edit).
HASH_PREFIX = "gcv1"
SUPPORTED_CONTRACT_VERSIONS = (1,)

#: the subtrees the content hash covers, in a fixed order.
HASHED_KEYS = ("gate_id", "contract_version", "configuration_fixed", "claims")

#: change kinds that require an explicit PI signature to authorise a transition.
PI_SIGNATURE_REQUIRED = ("THRESHOLD", "OBSERVABLE", "CONFIGURATION_FIXED")
CHANGE_KINDS = PI_SIGNATURE_REQUIRED + ("CONFIGURATION_RECORDED",)
#: values that look like a signature but are not one.
NON_SIGNATURES = {"", "pending", "todo", "tbd", "none", "n/a", "unsigned", "-"}

_EXCLUDED_KEY_NAMES = {"notes"}
_EXCLUDED_KEY_SUFFIXES = ("_note", "_notes")


# ----------------------------------------------------------------------------------
# canonicalisation + hashing
# ----------------------------------------------------------------------------------
def _excluded(key: str) -> bool:
    """Prose keys are excluded from the hash (documented in the module docstring)."""
    return key in _EXCLUDED_KEY_NAMES or key.endswith(_EXCLUDED_KEY_SUFFIXES)


def _strip(obj: Any) -> Any:
    """Recursively drop prose keys, sort dicts, and normalise floats to an exactly
    round-tripping text form so ``0.21`` and ``0.2100`` hash identically."""
    if isinstance(obj, dict):
        return {str(k): _strip(v) for k, v in sorted(obj.items(), key=lambda kv: str(kv[0]))
                if not _excluded(str(k))}
    if isinstance(obj, (list, tuple)):
        return [_strip(v) for v in obj]
    if isinstance(obj, bool) or obj is None or isinstance(obj, (int, str)):
        return obj
    if isinstance(obj, float):
        return repr(float(obj))
    return str(obj)


def _canon(obj: Any) -> str:
    return json.dumps(_strip(obj), sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def hashed_subset(contract: dict) -> dict:
    return {k: contract.get(k) for k in HASHED_KEYS}


def contract_hash(contract: dict) -> str:
    """Content hash of the (threshold + observable + configuration) subset."""
    return f"{HASH_PREFIX}:{_sha(_canon(hashed_subset(contract)))[:16]}"


def _list_key(index: int, item: Any) -> str:
    """Key a list element by its stable id when it has one, so reordering claims does
    not masquerade as N field changes (and renaming a claim_id does show up)."""
    if isinstance(item, dict):
        for id_key in ("claim_id", "id", "field", "name"):
            if item.get(id_key):
                return str(item[id_key])
    return str(index)


def leaf_hashes(contract: dict) -> dict[str, str]:
    """Per-leaf hashes of the hashed subset, keyed by dotted path.

    This is the half that lets the verifier say WHICH field moved
    (``claims[quantitative.…].threshold.value``) instead of only "the hash differs".
    """
    out: dict[str, str] = {}

    def walk(node: Any, path: str) -> None:
        if isinstance(node, dict):
            for k, v in sorted(node.items(), key=lambda kv: str(kv[0])):
                if _excluded(str(k)):
                    continue
                walk(v, f"{path}.{k}" if path else str(k))
        elif isinstance(node, (list, tuple)):
            for i, v in enumerate(node):
                walk(v, f"{path}[{_list_key(i, v)}]")
        else:
            out[path] = _sha(f"{path}={_canon(node)}")[:12]

    walk(hashed_subset(contract), "")
    return out


def diff_leaves(recorded: dict[str, str], current: dict[str, str]) -> list[str]:
    """Name the fields that moved between a recorded stamp and the contract on disk."""
    moved: list[str] = []
    for key in sorted(set(recorded) | set(current)):
        was, now = recorded.get(key), current.get(key)
        if was is None:
            moved.append(f"ADDED   {key}")
        elif now is None:
            moved.append(f"REMOVED {key}")
        elif was != now:
            moved.append(f"CHANGED {key}")
    return moved


# ----------------------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------------------
def load_contract(path: Path | str) -> dict:
    path = Path(path)
    doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    doc["_path"] = str(path)
    try:
        doc["_rel_path"] = str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        doc["_rel_path"] = str(path)
    return doc


def load_all_contracts(directory: Path | str | None = None) -> list[dict]:
    directory = Path(directory) if directory else CONTRACT_DIR
    return [load_contract(p) for p in sorted(directory.glob("*.yaml"))
            if p.name != CHANGES_FILE.name]


def load_changes(path: Path | str | None = None) -> dict:
    path = Path(path) if path else CHANGES_FILE
    if not path.exists():
        return {"changes": [], "historical_unsigned_changes": []}
    doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    doc.setdefault("changes", [])
    doc.setdefault("historical_unsigned_changes", [])
    return doc


def is_signed(row: dict) -> bool:
    sig = str(row.get("pi_signature") or "").strip().lower()
    return bool(sig) and sig not in NON_SIGNATURES


def authorising_row(changes: dict, *, gate_id: str, from_hash: str, to_hash: str,
                    kinds: tuple[str, ...] = PI_SIGNATURE_REQUIRED) -> tuple[dict | None, str]:
    """Find a ledger row that authorises a from_hash -> to_hash transition.

    Only rows in the ``changes:`` list can authorise. ``historical_unsigned_changes:``
    is a record of what already happened and authorises nothing — that separation is
    the point: writing the manoeuvre down is not the same as approving it.
    """
    for row in changes.get("changes", []):
        if str(row.get("gate_id")) != gate_id:
            continue
        if str(row.get("from_hash")) != from_hash or str(row.get("to_hash")) != to_hash:
            continue
        kind = str(row.get("kind") or "").upper()
        if kind not in kinds:
            return None, (f"ledger row {row.get('change_id')} has kind={kind!r}, which cannot "
                          f"authorise this transition (need one of {list(kinds)})")
        if kind in PI_SIGNATURE_REQUIRED and not is_signed(row):
            return None, (f"ledger row {row.get('change_id')} is kind={kind} but carries no PI "
                          f"signature (pi_signature={row.get('pi_signature')!r})")
        return row, f"authorised by {row.get('change_id')} (kind={kind}, PI-signed)"
    return None, "no ledger row covers this transition"


def disclosing_row(changes: dict, *, gate_id: str, fields: list[str]) -> tuple[dict | None, str]:
    """Find a CONFIGURATION_RECORDED row disclosing a per-run knob move (no PI
    signature required — disclosure, not a contract change)."""
    want = {f.split()[-1] for f in fields}
    for row in changes.get("changes", []):
        if str(row.get("gate_id")) != gate_id:
            continue
        if str(row.get("kind") or "").upper() != "CONFIGURATION_RECORDED":
            continue
        covered = {str(f) for f in (row.get("fields_moved") or [])}
        if want <= covered:
            return row, f"disclosed by {row.get('change_id')}"
    return None, "no CONFIGURATION_RECORDED ledger row discloses this move"


# ----------------------------------------------------------------------------------
# stamping (what a gate driver writes into its own artifact)
# ----------------------------------------------------------------------------------
def stamp_for(contract: dict, *, provenance: str = "RUNTIME",
              extra: dict | None = None) -> dict:
    """Build the ``gate_contract`` block a gate artifact must carry.

    ``provenance``:
      ``RUNTIME``    — written by the gate driver during the run it describes (the norm).
      ``RETROACTIVE``— attached to an artifact produced BEFORE the contract existed.
                       Kept explicit so nobody can read a retroactive stamp as evidence
                       that the run was gated at the time.
    """
    version = contract.get("contract_version")
    if version not in SUPPORTED_CONTRACT_VERSIONS:
        raise ValueError(f"unsupported contract_version={version!r} "
                         f"(supported: {list(SUPPORTED_CONTRACT_VERSIONS)})")
    stamp = {
        "gate_id": contract.get("gate_id"),
        "contract_path": contract.get("_rel_path"),
        "contract_version": version,
        "contract_hash": contract_hash(contract),
        "hash_prefix": HASH_PREFIX,
        "hashed_keys": list(HASHED_KEYS),
        "leaf_hashes": leaf_hashes(contract),
        "stamped_at": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "stamp_provenance": provenance,
        "verify_with": "python aleph/outputs/tag_kb/verify_gate_contracts.py --gate",
    }
    if extra:
        stamp.update(extra)
    return stamp


def stamp_from_gate_id(gate_id: str, *, directory: Path | str | None = None,
                       provenance: str = "RUNTIME", extra: dict | None = None) -> dict:
    """Load the contract for ``gate_id`` and build its stamp. Raises if absent —
    call this BEFORE the expensive part of a run so an unstamped gate fails fast."""
    for contract in load_all_contracts(directory):
        if str(contract.get("gate_id")) == gate_id:
            return stamp_for(contract, provenance=provenance, extra=extra)
    known = [c.get("gate_id") for c in load_all_contracts(directory)]
    raise FileNotFoundError(
        f"no gate contract for gate_id={gate_id!r} in {directory or CONTRACT_DIR} "
        f"(known: {known}). A gate must have a versioned contract before it may run.")


# ----------------------------------------------------------------------------------
# artifact paths + claim evaluation
# ----------------------------------------------------------------------------------
_INDEX_RE = re.compile(r"^(?P<name>[^\[\]]+)\[(?P<idx>-?\d+)\]$")

MISSING = object()


def resolve_path(doc: Any, path: str) -> Any:
    """Resolve a dotted artifact path (``final.residual_over_tension``, ``trace[-1].x``).
    Returns ``MISSING`` when absent — never a silent zero, which is how a broken gate
    quietly turns into a passing one."""
    cur = doc
    for part in str(path).split("."):
        if cur is MISSING:
            return MISSING
        m = _INDEX_RE.match(part)
        if m:
            name, idx = m.group("name"), int(m.group("idx"))
            if name:
                if not isinstance(cur, dict) or name not in cur:
                    return MISSING
                cur = cur[name]
            if not isinstance(cur, list) or not (-len(cur) <= idx < len(cur)):
                return MISSING
            cur = cur[idx]
            continue
        if isinstance(cur, dict) and part in cur:
            cur = cur[part]
        else:
            return MISSING
    return cur


def resolve_series(doc: Any, path: str) -> Any:
    """``trace[].n_bound`` -> the list of that field over every row."""
    if "[]" not in path:
        value = resolve_path(doc, path)
        return value if value is MISSING else list(value)
    base, _, rest = path.partition("[]")
    rows = resolve_path(doc, base)
    if rows is MISSING or not isinstance(rows, list):
        return MISSING
    rest = rest.lstrip(".")
    out = []
    for row in rows:
        value = resolve_path(row, rest) if rest else row
        if value is MISSING:
            return MISSING
        out.append(value)
    return out


_COMPARATORS = {
    "<": lambda a, b: a < b,
    "<=": lambda a, b: a <= b,
    ">": lambda a, b: a > b,
    ">=": lambda a, b: a >= b,
    "==": lambda a, b: a == b,
    "!=": lambda a, b: a != b,
}


def _compare(value: Any, comparator: str, target: Any) -> bool:
    op = _COMPARATORS.get(str(comparator))
    if op is None:
        raise ValueError(f"unknown comparator {comparator!r}")
    return bool(op(float(value), float(target)))


class ClaimResult:
    """Outcome of independently re-deriving one claim from the artifact."""

    __slots__ = ("claim_id", "claim_set", "ok", "measured", "status", "detail", "vacuous")

    def __init__(self, claim_id: str, claim_set: str, ok: bool | None, measured: Any,
                 status: str, detail: str, vacuous: bool = False) -> None:
        self.claim_id = claim_id
        self.claim_set = claim_set
        self.ok = ok
        self.measured = measured
        self.status = status      # RECOMPUTED | MEASURE_MISSING | UNSUPPORTED_KIND | ERROR
        self.detail = detail
        self.vacuous = vacuous


def evaluate_claim(claim: dict, artifact: dict) -> ClaimResult:
    """Independently re-derive a claim's PASS/FAIL from the artifact's own numbers.

    This is the second half of the guard: the contract's threshold is applied to the
    recorded measurement HERE, and compared against the boolean the gate script wrote.
    A driver whose inline literal drifts from the contract is then a hard failure
    instead of an invisible one.
    """
    claim_id = str(claim.get("claim_id"))
    claim_set = str(claim.get("claim_set") or "default")
    measure = claim.get("measure") or {}
    threshold = claim.get("threshold") or {}
    kind = str(measure.get("kind") or "")
    comparator = threshold.get("comparator")
    value = threshold.get("value")

    def missing(what: str) -> ClaimResult:
        return ClaimResult(claim_id, claim_set, None, None, "MEASURE_MISSING",
                           f"artifact has no {what}")

    try:
        if kind == "scalar":
            got = resolve_path(artifact, measure["path"])
            if got is MISSING:
                return missing(measure["path"])
            return ClaimResult(claim_id, claim_set, _compare(got, comparator, value), got,
                               "RECOMPUTED", f"{measure['path']}={got!r} {comparator} {value}")

        if kind == "all_of":
            got = {}
            for path in measure["paths"]:
                v = resolve_path(artifact, path)
                if v is MISSING:
                    return missing(path)
                got[path] = v
            ok = all(_compare(v, comparator, value) for v in got.values())
            return ClaimResult(claim_id, claim_set, ok, got, "RECOMPUTED",
                               f"all of {list(got)} {comparator} {value}")

        if kind in ("series_max", "series_min"):
            series = resolve_series(artifact, measure["series"])
            if series is MISSING:
                return missing(measure["series"])
            if not series:
                return ClaimResult(claim_id, claim_set, None, [], "MEASURE_MISSING",
                                   f"{measure['series']} is empty")
            agg = max(series) if kind == "series_max" else min(series)
            return ClaimResult(claim_id, claim_set, _compare(agg, comparator, value), agg,
                               "RECOMPUTED", f"{kind}({measure['series']})={agg!r} "
                                             f"{comparator} {value}")

        if kind == "series_conditional":
            target = resolve_series(artifact, measure["series"])
            where = resolve_series(artifact, measure["where_series"])
            if target is MISSING:
                return missing(measure["series"])
            if where is MISSING:
                return missing(measure["where_series"])
            rows = [t for t, w in zip(target, where)
                    if _compare(w, measure["where_comparator"], measure["where_value"])]
            if not rows:
                # VACUOUS: the predicate had zero applicable samples, so it cannot fail.
                # Reported, never silently counted as evidence.
                return ClaimResult(claim_id, claim_set, True, [], "RECOMPUTED",
                                   f"no row satisfies {measure['where_series']} "
                                   f"{measure['where_comparator']} {measure['where_value']} "
                                   f"— claim is VACUOUSLY true", vacuous=True)
            ok = all(_compare(t, comparator, value) for t in rows)
            return ClaimResult(claim_id, claim_set, ok, rows, "RECOMPUTED",
                               f"{len(rows)} conditional row(s) {comparator} {value}")

        if kind == "ordered_gt":
            got = []
            for path in measure["paths"]:
                v = resolve_path(artifact, path)
                if v is MISSING:
                    return missing(path)
                got.append(float(v))
            chain = all(a > b for a, b in zip(got, got[1:]))
            floor_ok = _compare(got[-1], comparator, value)
            return ClaimResult(claim_id, claim_set, bool(chain and floor_ok), got, "RECOMPUTED",
                               f"{' > '.join(f'{g!r}' for g in got)} and last {comparator} {value}")

        if kind == "all_true":
            got = resolve_path(artifact, measure["path"])
            if got is MISSING:
                return missing(measure["path"])
            if not isinstance(got, dict) or not got:
                return ClaimResult(claim_id, claim_set, None, got, "MEASURE_MISSING",
                                   f"{measure['path']} is not a non-empty bool map")
            return ClaimResult(claim_id, claim_set, all(bool(v) for v in got.values()), got,
                               "RECOMPUTED", f"all of {sorted(got)} are true")

        if kind == "difference_equals":
            a = resolve_path(artifact, measure["minuend"])
            b = resolve_path(artifact, measure["subtrahend"])
            ref = resolve_path(artifact, measure["reference"])
            for name, v in (("minuend", a), ("subtrahend", b), ("reference", ref)):
                if v is MISSING:
                    return missing(f"{measure[name]} ({name})")
            return ClaimResult(claim_id, claim_set, _compare(float(a) - float(b), comparator, ref),
                               float(a) - float(b), "RECOMPUTED",
                               f"{measure['minuend']}-{measure['subtrahend']}="
                               f"{float(a) - float(b)!r} {comparator} {ref!r}")

        return ClaimResult(claim_id, claim_set, None, None, "UNSUPPORTED_KIND",
                           f"measure.kind={kind!r} has no evaluator")
    except Exception as exc:  # a broken contract must be loud, not silently skipped
        return ClaimResult(claim_id, claim_set, None, None, "ERROR", f"{type(exc).__name__}: {exc}")


def recorded_verdict(claim: dict, artifact: dict) -> Any:
    """The boolean the gate script itself wrote for this claim (or MISSING)."""
    path = claim.get("verdict_path")
    if not path:
        return MISSING
    return resolve_path(artifact, path)

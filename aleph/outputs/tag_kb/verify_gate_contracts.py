#!/usr/bin/env python3
"""Gate-contract integrity gate — refuses a PASS whose contract moved since the FAIL.

THE THIRD KB GATE. `verify_sources.py` audits the LITERATURE side, `verify_runs.py` the
RESULTS side, `verify_params.py` the CONSTANTS side. This audits the **GATE** side: did
the threshold, the OBSERVABLE, or the CONFIGURATION change between a FAIL and the PASS
that superseded it?

WHY (35-agent audit 2026-07-25, failure mode M2 — four real recurrences):
  H.2 rebanding (KS p 0.05 -> 1.0e-6, +/-60% tol, and the *already-failing* trajectory
  verified to pass the new bands WITHOUT being re-run) - KU-5.5 xi_v 50-200 -> [1,400] -
  a "closes by construction" substitution the project itself later labelled a gate
  substitution in writing - GATE A closed after BOTH the observable (raw max|F| ->
  projected max|PF|) AND the configuration (myosin seed removed) changed.
The only guard was a CLAUDE.md sentence plus PI sign-off. This is the machine half.

WHAT IT CHECKS (per contract in docs/v2_audit/gate_contracts/)
  1. SELF_HASH_MISMATCH   the contract's declared hash != its recomputed content hash
                          (someone edited the contract and not the hash). Always blocking.
  2. CONTRACT_DRIFT       an artifact's recorded stamp hash != the contract now on disk,
                          with no authorising ledger row. Names the exact leaf fields that
                          moved (`claims[...].threshold.value`, `...observable.id`, ...).
  3. VERDICT_DRIFT        the contract's threshold, re-applied to the artifact's OWN
                          recorded measurement, disagrees with the boolean the gate script
                          wrote. Catches a driver literal drifting off the contract.
  4. SUPERSESSION_REFUSED a PASS supersedes a FAIL for the same claim under a DIFFERENT
                          contract hash and no PI-signed change row covers the transition.
  5. CONFIG_MOVED_AT_PASS a per-run knob (relax iterations, steps, device) moved at the
                          moment a FAIL became a PASS, with no disclosed ledger row.
  Reported, non-blocking: NO_ARTIFACT, UNSTAMPED, CONFIG_UNRECORDED, VACUOUS_CLAIM,
  ORDER_UNKNOWN, AUTHORISED_CHANGE. `--strict` promotes them to blocking.

VACUOUS_CLAIM is deliberately in the list: the audit's own worst-case is "a gate written
so a broken build cannot fail it". A claim whose predicate had zero applicable samples
passed without being tested, and this prints it instead of counting it as evidence.

RUN
  python aleph/outputs/tag_kb/verify_gate_contracts.py            # full report
  python aleph/outputs/tag_kb/verify_gate_contracts.py --gate     # CI gate, exit 1 on drift
  python aleph/outputs/tag_kb/verify_gate_contracts.py --check    # one-line drift summary
  python aleph/outputs/tag_kb/verify_gate_contracts.py \
      --stamp GATE-B.sf_motor <artifact.json> [--provenance RETROACTIVE]

  --artifact-root DIR   resolve artifact paths against DIR instead of the repo root
  --contract-dir DIR    use a different contract directory
  --changes-file FILE   use a different change ledger
  (the three overrides exist so the negative control can be run against a scratch tree
   without writing fabricated records into aleph/outputs/)

Pure stdlib + PyYAML. Reads committed artifacts only — no Notion token, no network, no
duckdb, no CUDA.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import gate_contract as gc
from gate_contract import MISSING

BLOCKING_CODES = {
    "SELF_HASH_MISMATCH", "CONTRACT_DRIFT", "VERDICT_DRIFT",
    "SUPERSESSION_REFUSED", "CONFIG_MOVED_AT_PASS", "CONTRACT_INVALID",
}
REPORT_CODES = {
    "NO_ARTIFACT", "UNSTAMPED", "CONFIG_UNRECORDED", "VACUOUS_CLAIM",
    "ORDER_UNKNOWN", "AUTHORISED_CHANGE", "CLAIM_NOT_RECOMPUTABLE",
    "VERDICT_UNDER_DIFFERENT_CONTRACT",
}


class Finding:
    __slots__ = ("code", "gate_id", "where", "message", "blocking")

    def __init__(self, code: str, gate_id: str, where: str, message: str, blocking: bool) -> None:
        self.code = code
        self.gate_id = gate_id
        self.where = where
        self.message = message
        self.blocking = blocking

    def line(self) -> str:
        mark = "FAIL" if self.blocking else "warn"
        return f"[{mark}] {self.code:<22} {self.gate_id} :: {self.where}\n       {self.message}"


# ----------------------------------------------------------------------------------
# stamps
# ----------------------------------------------------------------------------------
def sidecar_path(artifact_path: Path) -> Path:
    """`<dir>/<stem>.contract_stamp.json` — used only for RETROACTIVE stamps, so the run
    record itself stays byte-identical."""
    return artifact_path.parent / f"{artifact_path.stem}.contract_stamp.json"


def read_artifact(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def read_stamp(artifact_path: Path, artifact: dict) -> tuple[dict | None, str]:
    """Prefer the stamp embedded by the driver at run time; fall back to a sidecar."""
    inline = artifact.get("gate_contract") if isinstance(artifact, dict) else None
    if isinstance(inline, dict) and inline.get("contract_hash"):
        return inline, "inline"
    side = sidecar_path(artifact_path)
    if side.exists():
        doc = read_artifact(side)
        if isinstance(doc, dict) and doc.get("contract_hash"):
            return doc, f"sidecar {side.name}"
    return None, "absent"


# ----------------------------------------------------------------------------------
# per-artifact evaluation
# ----------------------------------------------------------------------------------
class ArtifactState:
    def __init__(self, rel_path: str, stamp: dict | None, stamp_source: str,
                 artifact: dict, order: int | None) -> None:
        self.rel_path = rel_path
        self.stamp = stamp
        self.stamp_source = stamp_source
        self.artifact = artifact
        self.order = order
        self.claim_ok: dict[str, bool | None] = {}       # re-derived from the contract
        self.claim_recorded: dict[str, bool | None] = {}  # what the run itself claimed
        self.set_verdict: dict[str, bool | None] = {}
        self.recorded_config: dict[str, object] = {}

    @property
    def hash(self) -> str | None:
        return self.stamp.get("contract_hash") if self.stamp else None

    @property
    def stamped_at(self) -> str:
        return str(self.stamp.get("stamped_at") or "") if self.stamp else ""

    def sort_key(self) -> tuple:
        return (self.stamped_at or "~", self.order if self.order is not None else 1 << 30,
                self.rel_path)


def evaluate_contract(contract: dict, changes: dict, artifact_root: Path,
                      findings: list[Finding], claim_log: list[str] | None = None) -> list[ArtifactState]:
    gate_id = str(contract.get("gate_id"))
    enforcement = str(contract.get("enforcement") or "BLOCKING").upper()
    blocking_gate = enforcement == "BLOCKING"

    def add(code: str, where: str, message: str, force_block: bool = False) -> None:
        blocking = force_block or (code in BLOCKING_CODES and blocking_gate)
        findings.append(Finding(code, gate_id, where, message, blocking))

    # --- 1. contract self-consistency (always blocking: a contract that lies about its
    #        own hash makes every downstream stamp meaningless) -----------------------
    if contract.get("contract_version") not in gc.SUPPORTED_CONTRACT_VERSIONS:
        add("CONTRACT_INVALID", contract.get("_rel_path", "?"),
            f"contract_version={contract.get('contract_version')!r} unsupported", force_block=True)
        return []
    current = gc.contract_hash(contract)
    declared = str(contract.get("contract_hash_declared") or "")
    if declared and declared != current:
        add("SELF_HASH_MISMATCH", contract.get("_rel_path", "?"),
            f"contract_hash_declared={declared} but the content hashes to {current} — the "
            f"contract was edited without updating its declared hash", force_block=True)
    if not declared:
        add("CONTRACT_INVALID", contract.get("_rel_path", "?"),
            f"no contract_hash_declared; content hashes to {current}", force_block=True)

    current_leaves = gc.leaf_hashes(contract)
    claims = contract.get("claims") or []
    recorded_fields = contract.get("configuration_recorded") or []
    declared_artifacts = ((contract.get("evidence") or {}).get("artifacts")) or []

    if not declared_artifacts:
        add("NO_ARTIFACT", "evidence.artifacts",
            f"contract declares NO run artifact. {str((contract.get('evidence') or {}).get('artifacts_status') or '').strip()[:220]}")
        return []

    states: list[ArtifactState] = []
    for entry in declared_artifacts:
        rel = str(entry.get("path"))
        path = artifact_root / rel
        if not path.exists():
            add("NO_ARTIFACT", rel, "declared artifact is not on disk")
            continue
        artifact = read_artifact(path)
        if artifact is None:
            add("NO_ARTIFACT", rel, "declared artifact is not readable JSON")
            continue
        stamp, source = read_stamp(path, artifact)
        state = ArtifactState(rel, stamp, source, artifact, entry.get("order"))

        # --- 2. contract drift -------------------------------------------------------
        if stamp is None:
            add("UNSTAMPED", rel,
                "artifact carries no gate_contract stamp (no inline block, no sidecar) — "
                "it cannot be tied to any contract version")
        else:
            if str(stamp.get("gate_id")) != gate_id:
                add("CONTRACT_DRIFT", rel,
                    f"stamp gate_id={stamp.get('gate_id')!r} != contract gate_id={gate_id!r}")
            elif state.hash != current:
                moved = gc.diff_leaves(dict(stamp.get("leaf_hashes") or {}), current_leaves)
                row, why = gc.authorising_row(changes, gate_id=gate_id,
                                              from_hash=str(state.hash), to_hash=current)
                detail = ("fields moved:\n         " + "\n         ".join(moved[:12])
                          if moved else "no leaf diff available (stamp recorded no leaf_hashes)")
                if row is None:
                    add("CONTRACT_DRIFT", rel,
                        f"artifact was measured under contract {state.hash}; the contract on "
                        f"disk is {current}. {why}.\n       {detail}")
                else:
                    add("AUTHORISED_CHANGE", rel, f"{why}: {state.hash} -> {current}")
            if str(stamp.get("stamp_provenance") or "").upper() == "RETROACTIVE":
                add("UNSTAMPED", rel,
                    "stamp provenance is RETROACTIVE — the run predates the contract, so it "
                    "is NOT evidence that the run was gated at the time")
            if not state.stamped_at:
                add("ORDER_UNKNOWN", rel, "stamp records no stamped_at; supersession ordering "
                                          "falls back to the contract's declared order")

        # --- 3. verdict drift: re-derive every claim from the artifact ---------------
        for claim in claims:
            claim_id = str(claim.get("claim_id"))
            result = gc.evaluate_claim(claim, artifact)
            recorded = gc.recorded_verdict(claim, artifact)
            state.claim_ok[claim_id] = result.ok
            state.claim_recorded[claim_id] = None if recorded is MISSING else bool(recorded)
            if claim_log is not None:
                rec = "MISSING" if recorded is MISSING else str(bool(recorded))
                agree = ("agree" if (recorded is not MISSING and result.ok is not None
                                     and bool(recorded) == bool(result.ok)) else "DISAGREE"
                         if recorded is not MISSING and result.ok is not None else "-")
                flag = "  [VACUOUS]" if result.vacuous else ""
                claim_log.append(
                    f"    {claim_id:<52} contract-recomputed={str(result.ok):<5} "
                    f"artifact-recorded={rec:<7} {agree}{flag}\n"
                    f"        threshold {claim.get('threshold', {}).get('comparator')} "
                    f"{claim.get('threshold', {}).get('value')}   |   {result.detail}")
            if result.vacuous:
                add("VACUOUS_CLAIM", f"{rel}::{claim_id}",
                    f"claim passed with ZERO applicable samples — {result.detail}")
            if result.status != "RECOMPUTED":
                add("CLAIM_NOT_RECOMPUTABLE", f"{rel}::{claim_id}",
                    f"{result.status}: {result.detail}")
                state.claim_ok[claim_id] = None if recorded is MISSING else bool(recorded)
                continue
            if recorded is MISSING:
                add("CLAIM_NOT_RECOMPUTABLE", f"{rel}::{claim_id}",
                    f"artifact has no recorded verdict at {claim.get('verdict_path')} "
                    f"(contract re-derives {result.ok} from {result.detail})")
                continue
            if bool(recorded) != bool(result.ok):
                # VERDICT_DRIFT means the DRIVER's criterion disagrees with the contract it
                # was measured under. If the artifact was measured under an older contract
                # version, the disagreement is expected and is not this code's business —
                # the moved contract is already reported as CONTRACT_DRIFT /
                # SUPERSESSION_REFUSED. Reporting it here too would make an AUTHORISED
                # change turn every historical FAIL permanently red, which is an incentive
                # to delete FAIL records. That would defeat the whole check.
                same_contract = state.hash == current
                add("VERDICT_DRIFT" if same_contract else "VERDICT_UNDER_DIFFERENT_CONTRACT",
                    f"{rel}::{claim_id}",
                    f"artifact records {claim.get('verdict_path')}={bool(recorded)} but the "
                    f"CONTRACT threshold re-applied to the artifact's own measurement gives "
                    f"{bool(result.ok)} ({result.detail})"
                    + (" — the gate script's criterion has drifted from the contract"
                       if same_contract else
                       f" — expected: this record was measured under contract {state.hash}, "
                       f"not the current {current}"))

        # A claim-set's PASS/FAIL for the SUPERSESSION check is what the RUN ITSELF
        # claimed (the recorded verdict), not the re-derived one: "a PASS that supersedes
        # a FAIL" is about the two claims made, and re-deriving the old FAIL under the NEW
        # contract is exactly the substitution being policed (it would make the FAIL
        # retroactively disappear). Disagreement between recorded and re-derived is
        # reported separately as VERDICT_DRIFT.
        for claim in claims:
            cs = str(claim.get("claim_set") or "default")
            claim_id = str(claim.get("claim_id"))
            ok = state.claim_recorded.get(claim_id)
            if ok is None:
                ok = state.claim_ok.get(claim_id)
            prev = state.set_verdict.get(cs, True)
            if ok is None or prev is None:
                state.set_verdict[cs] = None
            else:
                state.set_verdict[cs] = bool(prev and ok)

        # --- 4. per-run configuration fields must be recorded ------------------------
        for field in recorded_fields:
            name, apath = str(field.get("field")), str(field.get("artifact_path"))
            value = gc.resolve_path(artifact, apath)
            if value is MISSING:
                add("CONFIG_UNRECORDED", f"{rel}::{name}",
                    f"contract declares a per-run configuration field recorded at {apath}, "
                    f"but the artifact does not contain it — the run's configuration is not "
                    f"reconstructible from its own record")
            else:
                state.recorded_config[name] = value

        states.append(state)

    # --- 5. supersession: FAIL -> PASS transitions ----------------------------------
    ordered = sorted(states, key=ArtifactState.sort_key)
    claim_sets = sorted({str(c.get("claim_set") or "default") for c in claims})
    for cs in claim_sets:
        last_fail: ArtifactState | None = None
        for state in ordered:
            verdict = state.set_verdict.get(cs)
            if verdict is False:
                last_fail = state
                continue
            if verdict is True and last_fail is not None:
                if state.hash is None or last_fail.hash is None:
                    add("SUPERSESSION_REFUSED", f"{cs}: {last_fail.rel_path} -> {state.rel_path}",
                        "a PASS supersedes a FAIL but at least one of the two artifacts carries "
                        "no contract stamp, so the transition cannot be checked at all")
                elif state.hash != last_fail.hash:
                    moved = gc.diff_leaves(dict((last_fail.stamp or {}).get("leaf_hashes") or {}),
                                           dict((state.stamp or {}).get("leaf_hashes") or {}))
                    row, why = gc.authorising_row(changes, gate_id=gate_id,
                                                  from_hash=str(last_fail.hash),
                                                  to_hash=str(state.hash))
                    if row is None:
                        add("SUPERSESSION_REFUSED",
                            f"{cs}: {last_fail.rel_path} -> {state.rel_path}",
                            f"PASS was measured under contract {state.hash}; the FAIL it "
                            f"supersedes was measured under {last_fail.hash}. {why}. This is "
                            f"the M2 manoeuvre — a gate that moved between a FAIL and its "
                            f"PASS.\n       fields moved:\n         "
                            + "\n         ".join(moved[:12] or ["(no leaf diff available)"]))
                    else:
                        add("AUTHORISED_CHANGE",
                            f"{cs}: {last_fail.rel_path} -> {state.rel_path}", why)
                else:
                    diffs = [f"{k}: {last_fail.recorded_config.get(k)!r} -> {v!r}"
                             for k, v in state.recorded_config.items()
                             if k in last_fail.recorded_config and last_fail.recorded_config[k] != v]
                    if diffs:
                        row, why = gc.disclosing_row(changes, gate_id=gate_id, fields=diffs)
                        if row is None:
                            add("CONFIG_MOVED_AT_PASS",
                                f"{cs}: {last_fail.rel_path} -> {state.rel_path}",
                                "contract hash is unchanged, but per-run configuration moved at "
                                "the moment the FAIL became a PASS and no CONFIGURATION_RECORDED "
                                f"ledger row discloses it ({why}):\n         "
                                + "\n         ".join(diffs))
                        else:
                            add("AUTHORISED_CHANGE",
                                f"{cs}: {last_fail.rel_path} -> {state.rel_path}", why)
                last_fail = None
    return states


# ----------------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------------
def do_stamp(args: argparse.Namespace) -> int:
    contract_dir = Path(args.contract_dir) if args.contract_dir else gc.CONTRACT_DIR
    stamp = gc.stamp_from_gate_id(args.stamp, directory=contract_dir,
                                  provenance=args.provenance)
    target = Path(args.artifact).resolve()
    if not target.exists():
        print(f"artifact not found: {target}", file=sys.stderr)
        return 2
    if args.provenance.upper() == "RETROACTIVE":
        stamp["retroactive_note"] = (
            "The run this stamp describes was produced BEFORE the contract existed. The "
            "stamp attests which contract version the record was TRANSCRIBED INTO, not that "
            "the run was gated at the time. The run record itself is left byte-identical.")
        stamp["stamped_artifact"] = target.name
        out = target.parent / f"{target.stem}.contract_stamp.json"
        out.write_text(json.dumps(stamp, indent=2) + "\n", encoding="utf-8")
        print(f"wrote RETROACTIVE sidecar stamp -> {out}")
    else:
        artifact = read_artifact(target) or {}
        artifact["gate_contract"] = stamp
        target.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
        print(f"embedded RUNTIME stamp in {target}")
    print(f"  gate_id       {stamp['gate_id']}")
    print(f"  contract_hash {stamp['contract_hash']}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gate", action="store_true", help="CI mode: exit 1 on any blocking finding")
    ap.add_argument("--check", action="store_true", help="one-line summary (for refresh.sh)")
    ap.add_argument("--strict", action="store_true",
                    help="promote report-only findings (UNSTAMPED / NO_ARTIFACT / …) to blocking")
    ap.add_argument("--contract-dir", default=None)
    ap.add_argument("--changes-file", default=None)
    ap.add_argument("--artifact-root", default=None,
                    help="resolve artifact paths against this dir instead of the repo root")
    ap.add_argument("--stamp", default=None, metavar="GATE_ID",
                    help="stamp an artifact with a contract hash instead of verifying")
    ap.add_argument("--artifact", default=None, help="artifact to stamp (with --stamp)")
    ap.add_argument("--provenance", default="RUNTIME", choices=["RUNTIME", "RETROACTIVE"])
    ap.add_argument("--hashes", action="store_true", help="print each contract's content hash and exit")
    ap.add_argument("--show-claims", action="store_true",
                    help="print every claim re-derived from the artifact next to the boolean the "
                         "gate script recorded (this is the VERDICT_DRIFT check, shown)")
    args = ap.parse_args()

    contract_dir = Path(args.contract_dir) if args.contract_dir else gc.CONTRACT_DIR
    contracts = gc.load_all_contracts(contract_dir)

    if args.hashes:
        for contract in contracts:
            print(f"{gc.contract_hash(contract)}  {contract.get('gate_id')}  "
                  f"(declared {contract.get('contract_hash_declared')})")
        return 0
    if args.stamp:
        if not args.artifact:
            ap.error("--stamp requires --artifact")
        return do_stamp(args)

    changes = gc.load_changes(args.changes_file)
    artifact_root = Path(args.artifact_root).resolve() if args.artifact_root else gc.REPO_ROOT

    if not contracts:
        print(f"[FAIL] no gate contracts found in {contract_dir}")
        return 1 if args.gate else 0

    findings: list[Finding] = []
    claim_log: list[str] | None = [] if args.show_claims else None
    for contract in contracts:
        evaluate_contract(contract, changes, artifact_root, findings, claim_log)
    if args.strict:
        for f in findings:
            if f.code != "AUTHORISED_CHANGE":
                f.blocking = True

    blocking = [f for f in findings if f.blocking]
    warns = [f for f in findings if not f.blocking]

    if args.check:
        codes: dict[str, int] = {}
        for f in findings:
            codes[f.code] = codes.get(f.code, 0) + 1
        summary = ", ".join(f"{k} {v}" for k, v in sorted(codes.items())) or "clean"
        print(f"gate-contracts: {len(contracts)} contract(s); "
              f"{len(blocking)} blocking, {len(warns)} report-only  [{summary}]")
        return 0

    print(f"=== gate-contract integrity — {len(contracts)} contract(s) in "
          f"{contract_dir.relative_to(gc.REPO_ROOT) if contract_dir.is_relative_to(gc.REPO_ROOT) else contract_dir} ===")
    for contract in contracts:
        print(f"  {contract.get('gate_id'):<28} {gc.contract_hash(contract)}  "
              f"enforcement={contract.get('enforcement')}")
    if claim_log:
        print("--- claims re-derived from each artifact (contract threshold vs recorded verdict) ---")
        for line in claim_log:
            print(line)
        print()
    for f in blocking:
        print(f.line())
    for f in warns:
        print(f.line())
    print()
    print(f"BLOCKING {len(blocking)}   report-only {len(warns)}")
    if blocking:
        print("\nA blocking finding means a gate's threshold, observable, or configuration "
              "moved without a PI-signed change row, or a gate script's criterion has drifted "
              "from its contract. Fix the contract/ledger — do NOT edit the artifact.")
    if args.gate:
        return 1 if blocking else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

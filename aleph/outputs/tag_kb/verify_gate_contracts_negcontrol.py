#!/usr/bin/env python3
"""Negative control for the gate-contract verifier — proves the guard still BITES.

A gate nobody has watched fail is not a gate. This script manufactures the exact M2
manoeuvre against the REAL GATE-B.sf_motor contract and asserts that
``verify_gate_contracts.py --gate`` refuses it, then asserts that the honest cases pass.

It writes every fixture into a TEMP directory and passes ``--artifact-root`` /
``--contract-dir`` / ``--changes-file``, so no fabricated run record is ever created inside
``aleph/outputs/``. Nothing in the repo is modified.

The scenario, in all cases: the real committed sf_motor record is the FAIL (its quantitative
claim set genuinely fails — residual/tension 0.1553 > 0.01, and the record itself says
BLOCKED). A second artifact claims a PASS with the SAME measured physics. Only the contract
moved.

  1 threshold_moved           0.01 -> 0.20, unsigned ledger row      -> must FAIL
  2 threshold_moved_signed    same, with a PI-signed row             -> must PASS
  3 observable_substituted    residual/tension -> residual/(total head load), unsigned
                                                                     -> must FAIL
  4 observable_signed         same, with a PI-signed row             -> must PASS
  5 nothing_moved             contract identical, genuinely converged PASS (0.004)
                                                                     -> must PASS
  6 runknob_moved_at_pass     contract identical, accepted_steps 6 -> 40 at the PASS
                                                                     -> must FAIL
  7 self_hash_mismatch        contract edited, declared hash not updated
                                                                     -> must FAIL
  8 reformat_invariance       YAML round-trip must not change the hash

RUN
  python aleph/outputs/tag_kb/verify_gate_contracts_negcontrol.py
  exit 0 = the guard behaves; exit 1 = the guard has stopped biting (or has become noisy)
"""
from __future__ import annotations

import copy
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import gate_contract as gc
import yaml

HERE = Path(__file__).resolve().parent
VERIFY = HERE / "verify_gate_contracts.py"
REAL_CONTRACT = gc.CONTRACT_DIR / "GATE_B_sf_motor.yaml"
REAL_ARTIFACT = gc.REPO_ROOT / "aleph/outputs/ac/gate_b_sf_motor/sf_motor_native_gate.json"
ART_REL = "aleph/outputs/ac/gate_b_sf_motor"
QUANT = "quantitative.inner_solve_converged_relative_to_tension"


def _claim(contract: dict, claim_id: str) -> dict:
    for c in contract["claims"]:
        if c["claim_id"] == claim_id:
            return c
    raise KeyError(claim_id)


def _write_contract(contract: dict, out_dir: Path) -> str:
    out_dir.mkdir(parents=True, exist_ok=True)
    contract = copy.deepcopy(contract)
    contract.pop("_path", None)
    contract.pop("_rel_path", None)
    contract["contract_hash_declared"] = "PLACEHOLDER"
    contract["contract_hash_declared"] = gc.contract_hash(contract)
    (out_dir / "GATE_B_sf_motor.yaml").write_text(
        yaml.safe_dump(contract, sort_keys=False, allow_unicode=True), encoding="utf-8")
    return contract["contract_hash_declared"]


def _stamp(artifact: Path, contract_dir: Path) -> None:
    subprocess.run([sys.executable, str(VERIFY), "--stamp", "GATE-B.sf_motor",
                    "--artifact", str(artifact), "--contract-dir", str(contract_dir),
                    "--provenance", "RUNTIME"], check=True, capture_output=True)


def _ledger(kind: str, from_hash: str, to_hash: str, signature) -> dict:
    return {"ledger_version": 1, "changes": [{
        "change_id": "GC-NEGCONTROL", "gate_id": "GATE-B.sf_motor", "kind": kind,
        "from_hash": from_hash, "to_hash": to_hash,
        "fields_moved": ["accepted_steps"],
        "rationale": "negative-control fixture", "disclosed_by": "negative control",
        "pi_signature": signature, "date": "2026-07-25"}]}


def _run(contract_dir: Path, root: Path, changes: Path | None) -> tuple[int, str]:
    cmd = [sys.executable, str(VERIFY), "--gate", "--contract-dir", str(contract_dir),
           "--artifact-root", str(root)]
    if changes is not None:
        cmd += ["--changes-file", str(changes)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def build(base: Path, mutate, patch_pass, *, kind: str) -> tuple[Path, Path, Path, Path]:
    """Build {before-contract, after-contract, artifact root, ledgers} for one scenario."""
    orig = gc.load_contract(REAL_CONTRACT)
    before_dir = base / "contract_before"
    from_hash = _write_contract(orig, before_dir)

    after = copy.deepcopy(orig)
    mutate(after)
    after["evidence"]["artifacts"] = [
        {"path": f"{ART_REL}/sf_motor_native_gate.json", "role": "native_run_record", "order": 1},
        {"path": f"{ART_REL}/sf_motor_native_gate_v2.json", "role": "native_run_record", "order": 2},
    ]
    after_dir = base / "contracts"
    to_hash = _write_contract(after, after_dir)

    root = base / "root"
    art_dir = root / ART_REL
    art_dir.mkdir(parents=True, exist_ok=True)
    record = json.loads(REAL_ARTIFACT.read_text(encoding="utf-8"))
    # run-record@2 wraps the original gate payload under ``measurements``.
    # The negative-control fixture deliberately mutates that payload, not the
    # observation envelope; keep accepting pre-wrapper artifacts for history.
    real = (
        copy.deepcopy(record["measurements"])
        if str(record.get("schema", "")).endswith("run-record@2")
        else record
    )
    real.pop("gate_contract", None)

    fail_art = art_dir / "sf_motor_native_gate.json"
    fail_art.write_text(json.dumps(real, indent=2), encoding="utf-8")
    _stamp(fail_art, before_dir)

    passing = copy.deepcopy(real)
    passing["verdict_quantitative"][
        "inner_solve_converged_relative_to_tension"] = True
    passing["verdict"]["inner_solve_converged_relative_to_tension"] = True
    passing["quantitative_claim_status"] = "CONVERGED (negative-control fixture)"
    patch_pass(passing)
    pass_art = art_dir / "sf_motor_native_gate_v2.json"
    pass_art.write_text(json.dumps(passing, indent=2), encoding="utf-8")
    _stamp(pass_art, after_dir)

    pending = base / "changes_pending.yaml"
    signed = base / "changes_signed.yaml"
    pending.write_text(yaml.safe_dump(_ledger(kind, from_hash, to_hash, "PENDING")),
                       encoding="utf-8")
    signed.write_text(yaml.safe_dump(_ledger(kind, from_hash, to_hash, "Sungwook Yoon")),
                      encoding="utf-8")
    return after_dir, root, pending, signed


def main() -> int:
    if not REAL_ARTIFACT.exists() or not REAL_CONTRACT.exists():
        print("negative control cannot run: the real GATE-B contract/artifact is missing")
        return 1
    results: list[tuple[str, bool, str]] = []

    def check(name: str, ok: bool, detail: str) -> None:
        results.append((name, ok, detail))
        print(f"  [{'ok  ' if ok else 'FAIL'}] {name}: {detail}")

    tmp = Path(tempfile.mkdtemp(prefix="gate_contract_negcontrol_"))
    print(f"fixtures -> {tmp}\n")

    # 8. reformat invariance (cheap, and it underpins every stamp)
    orig = gc.load_contract(REAL_CONTRACT)
    rt = _write_contract(orig, tmp / "roundtrip")
    check("reformat_invariance", rt == gc.contract_hash(orig),
          f"yaml round-trip hash {rt} vs {gc.contract_hash(orig)}")

    # 1/2. threshold moved
    def mut_threshold(c: dict) -> None:
        # The current failing artifact measures 0.2286. Move the fixture
        # threshold past it so the synthetic successor is internally a PASS;
        # the control is whether that move is authorised, not VERDICT_DRIFT.
        _claim(c, QUANT)["threshold"]["value"] = 0.25

    cdir, root, pending, signed = build(tmp / "threshold", mut_threshold, lambda a: None,
                                        kind="THRESHOLD")
    code, out = _run(cdir, root, pending)
    named = "threshold.value" in out and "SUPERSESSION_REFUSED" in out
    check("threshold_moved_refused", code == 1 and named,
          f"exit={code}, names threshold.value + SUPERSESSION_REFUSED={named}")
    code, out = _run(cdir, root, signed)
    check("threshold_moved_signed_allowed", code == 0,
          f"exit={code} with a PI-signed ledger row")

    # 3/4. observable substituted
    def mut_observable(c: dict) -> None:
        q = _claim(c, QUANT)
        q["observable"]["id"] = "residual_over_total_head_load"
        q["observable"]["definition"] = ("free-node residual / SUM of all bound-head loads "
                                         "at the final accepted step")
        q["measure"]["path"] = "final.residual_over_total_head_load"

    def patch_metric(a: dict) -> None:
        a["final"]["residual_over_total_head_load"] = 0.0022546

    cdir, root, pending, signed = build(tmp / "observable", mut_observable, patch_metric,
                                        kind="OBSERVABLE")
    code, out = _run(cdir, root, pending)
    named = "observable.id" in out and "SUPERSESSION_REFUSED" in out
    check("observable_substituted_refused", code == 1 and named,
          f"exit={code}, names observable.id + SUPERSESSION_REFUSED={named}")
    code, out = _run(cdir, root, signed)
    check("observable_substituted_signed_allowed", code == 0,
          f"exit={code} with a PI-signed ledger row")

    # 5. nothing moved, and the PASS is a genuinely converged measurement
    cdir, root, pending, _ = build(tmp / "unchanged", lambda c: None,
                                   lambda a: a["final"].__setitem__("residual_over_tension", 0.004),
                                   kind="THRESHOLD")
    code, out = _run(cdir, root, pending)
    check("nothing_moved_allowed", code == 0, f"exit={code} (contract identical, 0.004 < 0.01)")

    # 6. contract identical, but a per-run knob moved at the moment FAIL -> PASS
    def patch_knob(a: dict) -> None:
        a["final"]["residual_over_tension"] = 0.004
        a["clock"]["steps"] = 40
        a["clock"]["after"] = 40

    cdir, root, pending, _ = build(tmp / "runknob", lambda c: None, patch_knob,
                                   kind="THRESHOLD")
    code, out = _run(cdir, root, pending)
    named_knob = "CONFIG_MOVED_AT_PASS" in out and "accepted_steps:" in out and "-> 40" in out
    check("runknob_moved_at_pass_refused",
          code == 1 and named_knob,
          f"exit={code}, names accepted_steps baseline -> 40={named_knob}")

    # 7. contract edited without updating its declared hash
    edited = tmp / "selfhash"
    cdir, root, pending, _ = build(edited, lambda c: None, lambda a: None, kind="THRESHOLD")
    target = cdir / "GATE_B_sf_motor.yaml"
    target.write_text(target.read_text(encoding="utf-8").replace("value: 0.01", "value: 0.2"),
                      encoding="utf-8")
    code, out = _run(cdir, root, pending)
    check("self_hash_mismatch_refused", code == 1 and "SELF_HASH_MISMATCH" in out,
          f"exit={code}, SELF_HASH_MISMATCH present")

    failed = [n for n, ok, _ in results if not ok]
    print(f"\n{len(results) - len(failed)}/{len(results)} controls behaved as required")
    if failed:
        print("THE GUARD IS NOT BITING (or has become noisy): " + ", ".join(failed))
        return 1
    print("the gate-contract guard refuses a moved threshold, a substituted observable, and "
          "a run-knob that moves at the PASS; it allows PI-signed changes and unchanged runs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

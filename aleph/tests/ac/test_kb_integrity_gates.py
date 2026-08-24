"""Regression tests for the KB integrity gates (verify_params / verify_runs).

These guard the *gate logic itself* — the disk-grounded auditors that enforce the
"no empirical magic numbers" + "no over-claimed results" HARD rules. They run in
the CPU CI suite on every push, so a change that silently breaks a verdict (e.g.
stops catching a value drift or a forbidden metric) fails the build.

Pure-Python: stdlib + pyyaml, no HOOMD, no Notion token, no kb.duckdb — the
verdict lookups read the committed source_audit_report.md and the real config
YAMLs. The gate scripts live outside the package (aleph/outputs/tag_kb), so we
add that dir to sys.path the same way the scripts are invoked in CI.

Location: this file lives under ``aleph/tests/ac/`` — NOT because it is an
Active-Cell test, but because ``tests/ac`` is the only path pytest collects
(``pyproject.toml`` ``testpaths``) and the only path CI runs. Sitting in the flat
``aleph/tests/`` directory it was collected ONLY by CI's explicit
``pytest aleph/tests`` argument, which also dragged in the quarantined HOOMD
legacy suite; narrowing CI to the quarantine without moving this file first would
have silently stopped running the KB gate-logic tests (2026-07-25 audit,
sequencing trap). Keep it inside the collected scope.

Run::

    conda activate ffn_sim
    pytest aleph/tests/ac/test_kb_integrity_gates.py -v
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

pytest.importorskip("yaml")

# parents[2] == aleph/ (this file is aleph/tests/ac/<name>.py). Asserted rather
# than assumed: if this file moves again, fail loudly here instead of degrading into
# a confusing ModuleNotFoundError on the verify_* import below.
TKB = Path(__file__).parents[2] / "outputs" / "tag_kb"
assert TKB.is_dir(), f"KB gate scripts not found at {TKB} — fix the parents[] depth after moving this file"
sys.path.insert(0, str(TKB))

import verify_params as VP  # noqa: E402
import verify_runs as VR    # noqa: E402

# a real, committed constant to ground the value check on
_CFG = "aleph/validation/oracles/configs/phase1_unit3.yaml"
_KEY = "cell.gamma_cortex"
_VAL = 0.5e-3


def _claim(**over):
    c = {"id": "T", "config": _CFG, "key": _KEY, "value": _VAL,
         "unit": "N/m", "ku": "KU-3.5", "citation_key": "X", "declared": "verified"}
    c.update(over)
    return c


# --------------------------------------------------------------------------- #
# verify_params — the parameter-provenance gate
# --------------------------------------------------------------------------- #
def test_params_verified_when_value_matches_and_citation_ok():
    v, _ = VP.audit_claim(_claim(), {"X": "OK"})
    assert v == "VERIFIED"


def test_params_value_drift_detected():
    # declared value disagrees with the real on-disk config value
    v, note = VP.audit_claim(_claim(value=9.9e-3), {"X": "OK"})
    assert v == "VALUE_DRIFT"
    assert "drift" in note.lower()


def test_params_absent_constant_is_value_drift():
    v, _ = VP.audit_claim(_claim(key="cell.does_not_exist"), {"X": "OK"})
    assert v == "VALUE_DRIFT"


def test_params_unsourced_when_no_citation():
    v, _ = VP.audit_claim(_claim(citation_key=None), {})
    assert v == "UNSOURCED"


def test_params_source_unverified_on_check_citation():
    v, _ = VP.audit_claim(_claim(), {"X": "CHECK"})
    assert v == "SOURCE_UNVERIFIED"


def test_params_source_suspect_on_fabrication_risk_citation():
    for bad in ("DOI_DEAD", "DOI_MISMATCH", "NO_DOI_NOMATCH"):
        v, _ = VP.audit_claim(_claim(), {"X": bad})
        assert v == "SOURCE_SUSPECT", bad


def test_params_overclaim_is_a_drift():
    # declaring 'verified' when the disk verdict is worse must register as drift
    disk = "SOURCE_UNVERIFIED"
    declared = VP.DECLARED_ALIASES["verified"]
    assert VP.SUSPICION[disk] < VP.SUSPICION[declared]


def test_params_real_manifest_has_no_drift():
    # regression guard: editing a config value without updating the manifest
    # (or over-declaring a citation) must light this up.
    rows = VP.run_audit()
    drifts = [r for r in rows if r[6]]  # r[6] == drift flag
    assert rows, "params_manifest.yaml produced no rows"
    assert not drifts, f"params manifest drifted: {drifts}"


def test_params_verdicts_loaded_from_committed_md():
    # the gate must be checkable without kb.duckdb / Notion (committed md path)
    verds = VP._load_verdicts()
    assert verds.get("Bell1978_Science") == "OK"
    assert verds.get("Chugh2017_NatCellBiol") == "CHECK"


# --------------------------------------------------------------------------- #
# verify_runs — the results-integrity gate
# --------------------------------------------------------------------------- #
def test_runs_forbidden_metric_is_retract():
    for m in ("basal-contact-hull", "basal-footprint", "convex-hull"):
        v, _ = VR.audit_claim({"metric": m, "artifact": "whatever.json"})
        assert v == "RETRACT", m


def test_runs_absent_artifact_is_needs_regen():
    v, _ = VR.audit_claim(
        {"metric": "top-down-silhouette", "artifact": "aleph/outputs/nope_missing.json"})
    assert v == "NEEDS_REGEN"


def test_runs_no_artifact_path_is_needs_regen():
    v, _ = VR.audit_claim({"metric": "n/a", "artifact": None})
    assert v == "NEEDS_REGEN"


def test_runs_real_manifest_has_no_drift():
    rows = VR.run_audit()
    drifts = [r for r in rows if r[5]]  # r[5] == drift flag
    assert rows, "results_manifest.yaml produced no rows"
    assert not drifts, f"results manifest drifted: {drifts}"

#!/usr/bin/env python3
"""Force-manifest integrity gate — the FORCES-side sibling of verify_sources.py (citations),
verify_runs.py (results) and verify_params.py (constants).

WHAT IT CHECKS. It diffs an OBSERVED force-channel dump (produced by
``aleph/ac/engine/forces_manifest.py``) against the checked-in production DECLARATION
(``forces_manifest.yaml``) and fails on any mismatch.

WHY IT EXISTS (2026-07-25). A lumped ``SettlingForce`` proxy at 19,000-47,000x real buoyant weight
was ON in every committed spreading run while the artifact recorded it as ``proxy_off``. It did the
flattening credited to the lamellipodium, the project's own regression gates locked the confound in
(so a CORRECT fix would have registered as a regression), and it invalidated every pre-06-21 DCM
spreading run including the capstone A/A0 law. It was found only by PI pushback. Four weeks later
the class recurred in a different engine (Biot FSI default-off in every landed native run). The
2026-06-20 closeout credits a permanent ``forces_manifest`` recurrence gate; the 2026-07-25 audit
checked and found it existed ONLY in that document and its Obsidian mirror. This is the gate.

THE TWO HISTORICAL SHAPES IT IS BUILT TO FAIL ON
  (i)  a channel declared off that is actually on   -> STATE_MISMATCH
  (ii) a channel declared on at a physiological value that is actually running at a
       convenience/default value                    -> VALUE_DRIFT / PROVENANCE_REGRESSION

Verdicts, suspicion-ranked worst first:
  STATE_MISMATCH        declared on/off disagrees with the observed state           <- worst
  INTENT_MISMATCH       the run's own config intent disagrees with its runtime state
  VALUE_DRIFT           channel on, but its magnitude != the declared production value
  PROVENANCE_REGRESSION channel on with a provenance the declaration forbids
                        (e.g. a CONVENIENCE default where a SOURCED value is required)
  DISQUALIFIED_SOURCE   a provenance label containing TEST / "NOT physiological" / proxy_off /
                        placeholder / provisional / diagnostic on a run claiming production
  MISSING_SOURCE        an enabled source-gated channel with an empty provenance label
  UNDECLARED_CHANNEL    the dump has a channel the declaration does not know about
                        (a force channel added without declaring its production state)
  MISSING_CHANNEL       the declaration requires a channel the dump does not contain
  UNOBSERVED_REQUIRED   the declaration requires this be observed and the dump could not observe it
  INCOMPLETE            explicitly pi_undeclared: no ratified production value exists yet
  OK                    matches the declaration

BLOCKING vs INCOMPLETE. Everything above INCOMPLETE is BLOCKING. INCOMPLETE is reported and
counted but does not block by default, because it means the DECLARATION is knowingly open (a PI
queue item), not that the run is wrong. Crucially, INCOMPLETE can only be reached through an
EXPLICIT ``pi_undeclared: true`` key: a channel simply left out of the declaration lands on
UNDECLARED_CHANNEL, which blocks. Silence cannot be used to dodge a requirement. ``--strict``
blocks on INCOMPLETE too.

RUN:
  conda activate ffn_sim
  # observe the current code defaults, then check them
  python -m aleph.ac.engine.forces_manifest --defaults --profile-claim resting_production \
      --out /tmp/dump.json
  python verify_forces.py --dump /tmp/dump.json --gate

  python verify_forces.py --artifact aleph/outputs/ac/.../run.json --profile legacy_unaudited
  python verify_forces.py --selftest        # negative + positive control; proves the gate can fail

Pure stdlib + PyYAML; reads only committed artifacts. No Notion token, no network, no CUDA.
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path
from typing import Any

import yaml

HERE = Path(__file__).parent
DECL_PATH = HERE / "forces_manifest.yaml"
REPO_ROOT = HERE.parents[2]  # tag_kb -> outputs -> ffn_sim -> repo root
DUMP_SCHEMA = "ffn-ac-forces-manifest-v1"

# Running this file by path (the documented invocation, and how kb-check calls its siblings) does
# not put the repo root on sys.path, so the observer import below would fail. Prepend it.
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# LOWER = worse (mirrors verify_sources / verify_runs / verify_params).
SUSPICION = {
    "STATE_MISMATCH": 0,
    "INTENT_MISMATCH": 1,
    "VALUE_DRIFT": 2,
    "PROVENANCE_REGRESSION": 3,
    "DISQUALIFIED_SOURCE": 4,
    "MISSING_SOURCE": 5,
    "UNDECLARED_CHANNEL": 6,
    "MISSING_CHANNEL": 7,
    "UNOBSERVED_REQUIRED": 8,
    "INCOMPLETE": 9,
    "OK": 10,
}
BLOCKING = {v for v, rank in SUSPICION.items() if rank < SUSPICION["INCOMPLETE"]}

_TRUE_STATES = {"on": True, "off": False}


class DeclarationError(RuntimeError):
    """The declaration itself is malformed — a checker that trusts a broken declaration is worse
    than no checker, so this is fatal rather than a verdict."""


def load_declaration(path: Path = DECL_PATH) -> dict[str, Any]:
    """Load and structurally validate ``forces_manifest.yaml``."""
    doc = yaml.safe_load(path.read_text())
    meta = doc.get("meta") or {}
    if meta.get("schema") != "ffn-ac-forces-manifest-decl-v1":
        raise DeclarationError(f"{path}: unexpected declaration schema {meta.get('schema')!r}")
    profiles = doc.get("profiles") or {}
    if not profiles:
        raise DeclarationError(f"{path}: no profiles declared")
    for pname, prof in profiles.items():
        if "certifies_production" not in prof:
            raise DeclarationError(f"{path}: profile {pname!r} must state certifies_production")
        for cname, req in (prof.get("channels") or {}).items():
            if not isinstance(req, dict):
                raise DeclarationError(f"{path}: {pname}.{cname} requirement must be a mapping")
            if not str(req.get("rationale", "")).strip():
                # A requirement with no ground is a magic number in gate clothing.
                raise DeclarationError(
                    f"{path}: {pname}.{cname} has no rationale — every requirement must cite the "
                    "PI decision, CLAUDE.md rule, or in-code declaration it comes from")
            if req.get("pi_undeclared") and not str(req.get("blocked_on", "")).strip():
                raise DeclarationError(
                    f"{path}: {pname}.{cname} is pi_undeclared but states no blocked_on")
            st = req.get("state")
            if st is not None and st not in ("on", "off", "any", True, False):
                raise DeclarationError(f"{path}: {pname}.{cname} bad state {st!r}")
    return doc


def _norm_state(raw: Any) -> str | None:
    """Normalize a declared state to ``'on'`` / ``'off'`` / ``'any'`` / ``None``.

    YAML 1.1 renders bare ``on``/``off`` as booleans, so both spellings must be accepted.
    """
    if raw is None:
        return None
    if isinstance(raw, bool):
        return "on" if raw else "off"
    s = str(raw).strip().lower()
    return s if s in ("on", "off", "any") else None


def _values_equal(observed: Any, declared: Any) -> bool:
    """Compare an observed magnitude to a declared one, tolerant of int/float spelling only."""
    if isinstance(observed, bool) or isinstance(declared, bool):
        return bool(observed) is bool(declared)
    if isinstance(observed, (int, float)) and isinstance(declared, (int, float)):
        if declared == 0:
            return float(observed) == 0.0
        return abs(float(observed) - float(declared)) <= 1e-12 * max(1.0, abs(float(declared)))
    return str(observed) == str(declared)


def check_dump(dump: dict[str, Any], decl: dict[str, Any], *,
               profile: str | None = None) -> dict[str, Any]:
    """Diff an observed dump against the declaration.  Returns a result record."""
    if dump.get("schema") != DUMP_SCHEMA:
        raise DeclarationError(
            f"dump schema {dump.get('schema')!r} != {DUMP_SCHEMA!r} — refusing to check an "
            "unrecognised dump rather than silently pass it")
    meta = decl.get("meta") or {}
    profiles = decl["profiles"]
    pname = profile or dump.get("profile_claim") or ""
    findings: list[dict[str, Any]] = []

    if not pname:
        findings.append({
            "channel": "<profile_claim>", "verdict": "STATE_MISMATCH",
            "detail": "the dump claims no profile; an unlabelled run cannot be certified against "
                      "any declaration (this is the proxy_off mislabel at the run level)",
        })
        pname = ""
    elif pname not in profiles:
        findings.append({
            "channel": "<profile_claim>", "verdict": "STATE_MISMATCH",
            "detail": f"profile_claim {pname!r} is not declared in forces_manifest.yaml "
                      f"(declared: {sorted(profiles)})",
        })

    prof = profiles.get(pname, {})
    reqs: dict[str, Any] = prof.get("channels") or {}
    certifies = bool(prof.get("certifies_production", False))
    forbidden = [s for s in (meta.get("forbidden_production_source_substrings") or [])]

    observed = {c["name"]: c for c in dump.get("channels", [])}

    # ── channels the declaration knows about ──────────────────────────────────────────────
    for cname, req in reqs.items():
        chan = observed.get(cname)
        if chan is None:
            findings.append({
                "channel": cname, "verdict": "MISSING_CHANNEL",
                "detail": "declared in the manifest but absent from the dump — the observer is "
                          "stale relative to the declaration, or the channel was removed",
            })
            continue
        obs = chan.get("observation")
        req_obs = req.get("require_observation_in")
        if req_obs and obs not in req_obs:
            findings.append({
                "channel": cname, "verdict": "UNOBSERVED_REQUIRED",
                "detail": f"observation={obs!r} but the declaration requires one of {req_obs}; "
                          f"{chan.get('notes', '')}".strip(),
            })
            continue

        if chan.get("intent_mismatch"):
            findings.append({
                "channel": cname, "verdict": "INTENT_MISMATCH",
                "detail": f"the run's config intent was {chan.get('config_intent')!r} but the "
                          f"observed runtime state is enabled={chan.get('enabled')!r} / "
                          f"value={chan.get('value')!r} — this is the SettlingForce shape "
                          f"(declared off, actually on). {chan.get('notes', '')}".strip(),
            })
            continue

        if req.get("pi_undeclared"):
            findings.append({
                "channel": cname, "verdict": "INCOMPLETE",
                "detail": f"no ratified production value; blocked_on: {req.get('blocked_on')}. "
                          f"observed enabled={chan.get('enabled')!r} value={chan.get('value')!r} "
                          f"provenance={chan.get('provenance')!r}",
            })
            # A pi_undeclared channel still owes its source label when it is switched ON.
            src_name = req.get("require_nonempty_source")
            if src_name and chan.get("enabled") is True:
                src = observed.get(src_name, {})
                if not str(src.get("value") or "").strip():
                    findings.append({
                        "channel": cname, "verdict": "MISSING_SOURCE",
                        "detail": f"enabled but {src_name} is empty — a source-gated magnitude is "
                                  "running with no provenance label",
                    })
            continue

        want = _norm_state(req.get("state"))
        got = chan.get("enabled")
        if want in ("on", "off"):
            if got is None:
                findings.append({
                    "channel": cname, "verdict": "UNOBSERVED_REQUIRED",
                    "detail": f"declaration requires state={want} but the dump could not observe "
                              f"this channel ({chan.get('notes', '')})".strip(),
                })
                continue
            if bool(got) is not _TRUE_STATES[want]:
                findings.append({
                    "channel": cname, "verdict": "STATE_MISMATCH",
                    "detail": f"declared {want.upper()} but observed enabled={got!r} "
                              f"(observation={obs}). rationale: {req.get('rationale', '').strip()}",
                })
                continue

        if "value" in req or "value_any_of" in req or "value_min" in req or "value_max" in req:
            v = chan.get("value")
            if v is None and chan.get("enabled") is not False:
                findings.append({
                    "channel": cname, "verdict": "UNOBSERVED_REQUIRED",
                    "detail": "declaration pins a value but the dump observed none",
                })
            elif "value" in req and not _values_equal(v, req["value"]):
                findings.append({
                    "channel": cname, "verdict": "VALUE_DRIFT",
                    "detail": f"observed {v!r} != declared production value {req['value']!r} "
                              f"(provenance={chan.get('provenance')}, observation={obs}). "
                              f"rationale: {req.get('rationale', '').strip()}",
                })
            elif "value_any_of" in req and not any(
                    _values_equal(v, a) for a in req["value_any_of"]):
                findings.append({
                    "channel": cname, "verdict": "VALUE_DRIFT",
                    "detail": f"observed {v!r} not in declared {req['value_any_of']!r}",
                })
            else:
                lo, hi = req.get("value_min"), req.get("value_max")
                if isinstance(v, (int, float)) and not isinstance(v, bool):
                    if lo is not None and float(v) < float(lo):
                        findings.append({"channel": cname, "verdict": "VALUE_DRIFT",
                                         "detail": f"observed {v!r} below declared minimum {lo!r}"})
                    if hi is not None and float(v) > float(hi):
                        findings.append({"channel": cname, "verdict": "VALUE_DRIFT",
                                         "detail": f"observed {v!r} above declared maximum {hi!r}"})

        allow = req.get("provenance_allow")
        if allow and chan.get("provenance") not in allow:
            findings.append({
                "channel": cname, "verdict": "PROVENANCE_REGRESSION",
                "detail": f"provenance {chan.get('provenance')!r} not in declared {allow!r} — an "
                          f"enabled production channel is running on a magnitude whose grounding "
                          f"the declaration does not accept. rationale: "
                          f"{req.get('rationale', '').strip()}",
            })

        src_name = req.get("require_nonempty_source")
        if src_name and chan.get("enabled") is True:
            src = observed.get(src_name, {})
            if not str(src.get("value") or "").strip():
                findings.append({
                    "channel": cname, "verdict": "MISSING_SOURCE",
                    "detail": f"enabled but {src_name} is empty",
                })

        if certifies and req.get("forbid_production_substrings") and forbidden:
            text = str(chan.get("value") or "")
            hit = [s for s in forbidden if s and s in text]
            if hit:
                findings.append({
                    "channel": cname, "verdict": "DISQUALIFIED_SOURCE",
                    "detail": f"label {text!r} contains {hit!r}; this run labels itself "
                              "non-physiological yet claims a production profile",
                })

    # ── channels present in the dump that the declaration does not know ───────────────────
    if prof.get("require_channel_coverage"):
        for cname, chan in observed.items():
            if cname in reqs:
                continue
            findings.append({
                "channel": cname, "verdict": "UNDECLARED_CHANNEL",
                "detail": f"observed (kind={chan.get('kind')}, enabled={chan.get('enabled')!r}, "
                          f"value={chan.get('value')!r}, provenance={chan.get('provenance')}) but "
                          f"NOT declared for profile {pname!r}. A force channel or flag-gated path "
                          f"must have a declared production state before it can run in production. "
                          f"code_site={chan.get('code_site')}",
            })

    worst = min((SUSPICION[f["verdict"]] for f in findings), default=SUSPICION["OK"])
    blocking = [f for f in findings if f["verdict"] in BLOCKING]
    incomplete = [f for f in findings if f["verdict"] == "INCOMPLETE"]
    return {
        "profile": pname,
        "certifies_production": certifies,
        "n_channels_observed": len(observed),
        "n_channels_declared": len(reqs),
        "findings": sorted(findings, key=lambda f: (SUSPICION[f["verdict"]], f["channel"])),
        "blocking": blocking,
        "incomplete": incomplete,
        "worst_rank": worst,
        "coverage": dump.get("coverage", {}),
        "run_label": dump.get("run_label", ""),
    }


def _report(res: dict[str, Any], *, verbose: bool) -> None:
    print(f"profile              : {res['profile'] or '<none>'}")
    print(f"certifies_production : {res['certifies_production']}")
    print(f"run_label            : {res['run_label']}")
    cov = res.get("coverage") or {}
    print(f"channels observed    : {res['n_channels_observed']} "
          f"(declared for this profile: {res['n_channels_declared']})")
    if cov:
        print(f"observation strength : {cov.get('by_observation')}")
        print(f"runtime_observed={cov.get('runtime_observed')}  env_observed={cov.get('env_observed')}")
    counts: dict[str, int] = {}
    for f in res["findings"]:
        counts[f["verdict"]] = counts.get(f["verdict"], 0) + 1
    print(f"verdicts             : {counts or {'OK': 1}}")
    print()
    shown = res["findings"] if verbose else [f for f in res["findings"]
                                            if f["verdict"] in BLOCKING][:60]
    for f in shown:
        print(f"[{f['verdict']}] {f['channel']}")
        print(f"    {f['detail']}")
    if not verbose and res["incomplete"]:
        print(f"\n({len(res['incomplete'])} INCOMPLETE / pi_undeclared channels suppressed; "
              f"use --verbose to list the PI queue)")
    print()
    if res["blocking"]:
        print(f"FORCES-MANIFEST DRIFT: {len(res['blocking'])} blocking finding(s).")
    elif not res["certifies_production"]:
        print("NO BLOCKING DRIFT — but this profile does NOT certify production "
              "(certifies_production: false). This run must not be reported as a production result.")
    else:
        print(f"FORCES MANIFEST OK — {len(res['incomplete'])} pi_undeclared channel(s) remain "
              "open as PI queue items.")


# ─────────────────────────────────────────────────────────────────────────────────────────────
# Self-test: the gate must be watched failing.  Injects the two real historical shapes into a
# TRUE dump and asserts the checker names the exact channels, then asserts the true dump passes.
# ─────────────────────────────────────────────────────────────────────────────────────────────
def _inject_defects(dump: dict[str, Any]) -> tuple[dict[str, Any], dict[str, str]]:
    """Return (corrupted dump, {channel: expected verdict}) for the two historical shapes."""
    bad = copy.deepcopy(dump)
    by_name = {c["name"]: c for c in bad["channels"]}
    expect: dict[str, str] = {}

    # Shape (i): a channel DECLARED OFF that is actually ON. `multigrid` is declared
    # `state: off` in resting_production; flip the observation to on, as an AC_MG=1 env
    # override would.  This is the SettlingForce shape.
    ch = by_name["multigrid"]
    ch["enabled"], ch["value"], ch["observation"] = True, True, "env"
    ch["intent_mismatch"] = False  # force the STATE_MISMATCH path, not the intent path
    expect["multigrid"] = "STATE_MISMATCH"

    # Shape (ii): a channel DECLARED ON AT A PHYSIOLOGICAL VALUE that is actually running at
    # the convenience/default value.  membrane_mesh_subdivisions is declared 8 (PI-ratified
    # 2026-07-24 resting production); 3 is the backward-compat default.
    ch = by_name["membrane_mesh_subdivisions"]
    ch["value"] = 3
    expect["membrane_mesh_subdivisions"] = "VALUE_DRIFT"

    # Shape (ii'), the provenance form: an enabled production channel whose magnitude is a
    # CONVENIENCE default where the declaration requires SOURCED.
    ch = by_name["membrane_erm_stiffness"]
    ch["provenance"] = "CONVENIENCE"
    ch["value"] = 100.0  # the retired ~46x-too-soft bleb-PDE continuum value
    expect["membrane_erm_stiffness"] = "VALUE_DRIFT"
    return bad, expect


class _ComposedStandIn:
    """A duck-typed stand-in for the composed ``AssembledCell``, for the CONTROL PAIR only.

    WHAT IT IS. The observer's runtime reads are ``cell.steric is not None``,
    ``int(cell.n_tri) > 0`` and friends — documented attribute reads with no device access. This
    object carries those exact attribute names so the control pair exercises the REAL predicate
    functions in ``forces_manifest._FORCE_TERMS`` on a CUDA-less host.

    WHAT IT IS NOT. It is not a cell and it computes no force. It certifies that the OBSERVER and
    the CHECKER work end to end; it does not certify any physics. A production certification
    requires a dump taken from a real ``AssembledCell`` on the CUDA lane (the gbook), where these
    same attributes are the actual composed primitives.
    """

    def __init__(self) -> None:
        self.n_tri = 1_413_720          # cortex bending triples present
        self.n_xl = 1_413_720           # crosslinks present
        self.branch_triples_d = None    # formin-only cortex: no Arp2/3 branch junctions
        self.n_branch = 0
        self.myosin = object()
        self.nucleus = object()
        self.membrane = object()
        self.membrane_pressure = object()
        self.steric = object()
        self.pressure = object()


def _selftest(decl: dict[str, Any], *, verbose: bool) -> int:
    from aleph.ac.cell.assemble import CellConfig
    from aleph.ac.engine.forces_manifest import dump_force_channels

    # A dump configured to the declared resting_production profile, as a compliant run emits.
    cfg = CellConfig(
        membrane_subdivisions=8, erm_radial_pairing=True, overlap_free_cortex=True,
        nmii_straddle_placement=True)
    solver_kwargs = {
        "fiber_quotient_coarse": False, "multigrid": False, "mg_assembled_coarse": False,
        "cg_check_every": 0, "coarse_modes": 0, "fq_coarse_mode": "pathA",
        "fq_coarse_iterations": 40, "mg_fiber_quotient": True, "mg_smooth": 2, "mg_omega": 0.8,
    }
    good = dump_force_channels(
        cfg, cell=_ComposedStandIn(), solver_kwargs=solver_kwargs, env={},
        profile_claim="resting_production",
        run_label="selftest control pair (declared resting_production config; FORCE state read "
                  "off a duck-typed composed stand-in, NOT a real cell)")

    bad, expect = _inject_defects(good)

    print("=" * 78)
    print("NEGATIVE CONTROL — deliberately-wrong dump; the gate MUST fail and name the channels")
    print("=" * 78)
    res_bad = check_dump(bad, decl)
    _report(res_bad, verbose=verbose)
    # A channel may legitimately trip more than one independent detector (a value that is both off
    # its declared magnitude AND of a forbidden provenance), so collect the verdict SET per channel
    # and require the expected verdict to be among them.
    named: dict[str, set[str]] = {}
    for f in res_bad["blocking"]:
        named.setdefault(f["channel"], set()).add(f["verdict"])
    ok_neg = bool(res_bad["blocking"])
    for chan, want in expect.items():
        got = named.get(chan, set())
        caught = want in got
        if not caught:
            ok_neg = False
        print(f"  expect {chan:32s} -> {want:22s} got {sorted(got) or None}  "
              f"[{'CAUGHT' if caught else 'MISSED'}]")

    print()
    print("=" * 78)
    print("POSITIVE CONTROL — the true dump; the gate MUST pass")
    print("=" * 78)
    res_good = check_dump(good, decl)
    _report(res_good, verbose=verbose)
    ok_pos = not res_good["blocking"]

    print()
    print(f"negative control failed as expected : {ok_neg}")
    print(f"positive control passed as expected : {ok_pos}")
    if ok_neg and ok_pos:
        print("SELFTEST PASS — the gate can fail, and does not fail the true configuration.")
        return 0
    print("SELFTEST FAIL — this gate is not trustworthy.")
    return 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--dump", type=str, default="", help="observed force-channel dump JSON")
    src.add_argument("--artifact", type=str, default="",
                     help="run report JSON; reconstruct a config-strength dump from it")
    src.add_argument("--defaults", action="store_true",
                     help="observe the CURRENT CellConfig code defaults and check them")
    src.add_argument("--selftest", action="store_true",
                     help="run the negative + positive control (proves the gate can fail)")
    ap.add_argument("--profile", type=str, default="",
                    help="declared profile to check against (default: the dump's profile_claim)")
    ap.add_argument("--gate", action="store_true", help="exit 1 on any BLOCKING finding")
    ap.add_argument("--strict", action="store_true",
                    help="with --gate, also exit 1 on INCOMPLETE (pi_undeclared) channels")
    ap.add_argument("--verbose", action="store_true", help="list every finding incl. INCOMPLETE")
    ap.add_argument("--json-out", type=str, default="", help="write the result record here")
    args = ap.parse_args()

    decl = load_declaration()

    if args.selftest:
        return _selftest(decl, verbose=args.verbose)

    if args.dump:
        dump = json.loads(Path(args.dump).read_text())
    elif args.artifact:
        from aleph.ac.engine.forces_manifest import dump_from_artifact
        dump = dump_from_artifact(args.artifact, profile_claim=args.profile or "legacy_unaudited")
    elif args.defaults:
        from aleph.ac.cell.assemble import CellConfig
        from aleph.ac.engine.forces_manifest import dump_force_channels
        dump = dump_force_channels(
            CellConfig(), profile_claim=args.profile or "resting_production",
            run_label="CellConfig() code defaults")
    else:
        ap.error("one of --dump / --artifact / --defaults / --selftest is required")

    res = check_dump(dump, decl, profile=args.profile or None)
    _report(res, verbose=args.verbose)
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(res, indent=2, default=str) + "\n")
        print(f"wrote {args.json_out}")

    if args.gate:
        if res["blocking"]:
            return 1
        if args.strict and res["incomplete"]:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

r"""The one run-record shape every engine driver writes — self-stamping, two-axis, and gate-classified.

WHAT THIS FIXES.  The 2026-07-25 whole-repo audit's highest-severity finding is that **no artifact in
this repo stamps the build it was measured on**, which is why ``STATE.md`` has to reconstruct a "build
commit" for every tier-(a) row from the implementing commit and why one result sits in the
not-quotable list purely because its build cannot be established.  Two further defects have the same
root: drivers hand-typed their evidence rung as a free string (one stamped ``"evidence": "CONNECTED"``
as a constant, recording an intention rather than an observation), and gate verdicts were bare
``"PASS"``/``"FAIL"`` strings with no single definition, so a run whose residual swamped its signal
could and did report a PASS.

The vocabulary that fixes all three landed on 2026-07-28 (PI decisions D1-B, D4, D7, D8) and had **no
callers**.  This module is the caller: one function that takes a driver's measurements and returns the
record, with the stamping, the two-axis label, and the verdict classification applied on the way
through rather than left to each driver to remember.

WHAT IT REFUSES TO DO.  It does not decide a rung — the caller passes an
:class:`~aleph.engine.contracts.EvidenceLabel` whose ``basis`` names the measurement that earned
it, and that class already rejects an empty basis.  It does not choose a void ceiling — the caller
declares a :class:`~aleph.engine.contracts.VoidCeiling` BEFORE the run, and passing one built from
the run's own residual would be gate-loosening that no type can catch.  It computes no physics.

Sanity Gate:
    * dimensional: this module holds no physical quantity of its own; every measurement passes through
      verbatim under the caller's key, units included in the key name by convention.
    * boundary: an unavailable git commit records ``"unknown"`` with the reason, never a fabricated
      hash; a driver with no gate to classify records ``verdict: null`` rather than a default PASS.
    * conservation/invariant: the config hash covers the canonical JSON of the config block, so a
      record whose config changed cannot keep the old hash.
    * numerical: hashing is over sorted-key canonical JSON so it is invariant to dict ordering.
    * sign-sense: not applicable.
    * measurement-protocol: ``t0`` is a REQUIRED field, because a run that reports only its end state
      cannot show that anything emerged; the plan's step-0 telemetry rule names it explicitly.
      ``timing`` is REQUIRED for the mirror reason — a run that records only its physics cannot say
      what it cost — and it is stamped with the run's verdict so a cost can never be lifted out of a
      record that did not establish the convergence the cost belongs to.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any, Mapping

from aleph.engine.contracts import EvidenceLabel, GateVerdict, VoidCeiling, classify_gate

__all__ = [
    "ARTIFACT_SCHEMA",
    "build_stamp",
    "config_hash",
    "observation_artifact",
    "timing_block",
    "write_artifact",
]

#: Schema identifier, bumped when the record shape changes in a way a reader must notice.
#: ``@2`` (2026-07-28) made ``timing`` a required field — see :func:`timing_block`.
ARTIFACT_SCHEMA: str = "ac.engine.observe/run-record@2"


def timing_block(
    *,
    wall_seconds: float,
    physical_time_s: float | None = None,
    n_steps: int | None = None,
    n_inner_iterations: int | None = None,
    device_note: str | None = None,
) -> dict[str, Any]:
    """Return the timing record, with the rate derived and the comparability question left open.

    WHY THIS IS A REQUIRED FIELD AND NOT A CONVENIENCE.  As of 2026-07-28 **no artifact in this repo
    recorded a wall-clock at all**, so answering "is the new engine faster than the old one" meant
    reconstructing timings from log file mtimes.  Worse, the one place a native timing does exist
    (``outputs/ac/foundation-hardening/native_resting_70686.json``) records ``outer_step_wall_s = 1.13``
    and ``inner_converged: false`` in the same breath — and the two were read apart, which is how a cost
    that belongs to a NON-converged step could be quoted as the engine's speed.

    A timing is not a property of an engine.  It is a property of an engine **at a stated convergence**,
    because any solver is arbitrarily fast if it is allowed to stop early.  So this block is deliberately
    incapable of standing alone: :func:`observation_artifact` stamps it with the run's own verdict and
    marks it non-comparable whenever that verdict is absent, failing, or VOID.

    Args:
        wall_seconds: Measured wall-clock of the run [s].  Positive and finite.
        physical_time_s: Biological time the run covered [s], if it stepped a clock.  Supplying it is
            what yields the real-time factor, which is the number anyone actually asks for.
        n_steps: Accepted outer physical steps, if any.
        n_inner_iterations: Total inner solve iterations, if the driver counts them.  This is usually the
            term that dominates the cost, and omitting it makes two runs incomparable even at equal
            convergence.
        device_note: Anything about the device or configuration a later reader would need — an
            instrumented run is not a benchmark and should say so here.

    Returns:
        The JSON-able timing block, with ``real_time_factor`` and per-step / per-iteration costs derived
        where the inputs allow it.

    Raises:
        ValueError: If ``wall_seconds`` is not positive-finite, or a supplied count is not positive.
    """
    if not (isinstance(wall_seconds, (int, float)) and wall_seconds > 0.0
            and wall_seconds == wall_seconds and wall_seconds != float("inf")):
        raise ValueError(f"wall_seconds must be positive-finite; got {wall_seconds!r}")
    for label, value in (("physical_time_s", physical_time_s), ("n_steps", n_steps),
                         ("n_inner_iterations", n_inner_iterations)):
        if value is not None and not value > 0:
            raise ValueError(f"{label} must be positive when supplied; got {value!r}")

    block: dict[str, Any] = {"wall_seconds": float(wall_seconds)}
    if physical_time_s is not None:
        block["physical_time_s"] = float(physical_time_s)
        block["real_time_factor"] = float(wall_seconds) / float(physical_time_s)
    if n_steps is not None:
        block["n_steps"] = int(n_steps)
        block["wall_seconds_per_step"] = float(wall_seconds) / float(n_steps)
    if n_inner_iterations is not None:
        block["n_inner_iterations"] = int(n_inner_iterations)
        block["wall_seconds_per_inner_iteration"] = float(wall_seconds) / float(n_inner_iterations)
    if device_note:
        block["device_note"] = str(device_note)
    return block


def build_stamp(
    repo_root: Path | str | None = None, *, declared_commit: str | None = None
) -> dict[str, Any]:
    """Return the build the run is being measured on: commit, dirty flag, branch, and how it was got.

    A dirty worktree is recorded rather than rejected — development runs are legitimate — but the flag
    means the commit alone does not reproduce the run, and a reader must treat the result accordingly.

    ``declared_commit`` exists because the machine that RUNS is not always the machine that holds the
    repository: this project's native runs execute on a synced tree that is not a git checkout, where a
    git-only stamp degrades silently to ``"unknown"`` — reintroducing precisely the missing-build defect
    this function exists to close.  A declared commit is recorded with ``source: "declared"`` so it is
    never mistaken for one git verified, and when git IS available the two are cross-checked and any
    disagreement is recorded rather than resolved.

    Args:
        repo_root: Directory to resolve the repository from; defaults to this file's repository.
        declared_commit: The commit the caller asserts the run was made on, for runs on a non-git tree.

    Returns:
        A dict with ``commit``, ``dirty``, ``branch``, ``source``, and — when relevant — ``reason`` or
        ``declared_vs_git_mismatch``.
    """
    root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[4]

    def git(*args: str) -> str | None:
        try:
            done = subprocess.run(
                ("git", "-C", str(root), *args),
                capture_output=True, text=True, timeout=15, check=False,
            )
        except (OSError, subprocess.SubprocessError):
            return None
        return done.stdout.strip() if done.returncode == 0 else None

    commit = git("rev-parse", "HEAD")
    if commit is None:
        if declared_commit:
            return {"commit": str(declared_commit), "dirty": None, "branch": None,
                    "source": "declared",
                    "reason": f"{root} is not a git repository; the commit is the caller's assertion "
                              "and was not verified here"}
        return {"commit": "unknown", "dirty": None, "branch": None, "source": "unavailable",
                "reason": f"git unavailable or {root} is not a repository, and no commit was declared"}
    status = git("status", "--porcelain")
    stamp: dict[str, Any] = {
        "commit": commit,
        "dirty": bool(status) if status is not None else None,
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "source": "git",
    }
    if declared_commit and not commit.startswith(str(declared_commit)):
        stamp["declared_vs_git_mismatch"] = {"declared": str(declared_commit), "git": commit}
    return stamp


def config_hash(config: Mapping[str, Any]) -> str:
    """Return the sha256 of a config block's canonical JSON.

    Args:
        config: The run's configuration mapping.

    Returns:
        The hex digest, prefixed ``sha256:``.
    """
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":"), default=str)
    return "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()


#: Configuration knobs that CHANGE THE POPULATION, and the census key that must record what each one
#: actually produced.  A knob here without its census key is the defect this whole guard exists for: on
#: 2026-07-29 two runs 1.73x apart in node count and 1.95x apart in step cost both stamped
#: ``fraction_of_native: 1.0``, because the fraction saw the cortex and the difference was the membrane.
#: Extend this dict in the SAME change that adds a knob which resizes any compartment — and note that
#: TWO SPELLINGS of the membrane knob exist across drivers (`membrane_subdiv` in the cortex-motor driver,
#: `membrane_subdivisions` everywhere else). The first version of this registry knew one of them, which
#: made the new guard an instance of the very defect it was written to catch: narrower than its name.
POPULATION_KNOBS: dict[str, str] = {
    "membrane_subdiv": "membrane_vertices",
    "membrane_subdivisions": "membrane_vertices",     # the SAME knob, spelled differently across drivers
    "nucleus_subdivisions": "nucleus_vertices",
}

#: Census keys that promise a whole-cell quantity while no driver can currently compute one.  Banned by
#: name rather than corrected, because the honest replacement depends on a PI answer that does not exist:
#: the membrane's PHYSIOLOGICAL subdivision is sourced nowhere, so there is no denominator for a true
#: whole-cell fraction.  Until there is, a per-compartment fraction that says what it measures is the
#: most a record may claim.
BANNED_CENSUS_KEYS: dict[str, str] = {
    "fraction_of_native": (
        "it was computed as cortex_filaments/70686 — the CORTEX only — so it reads 1.0 for a run whose "
        "membrane was reduced to a quarter resolution. Record `fraction_of_native_cortex` (which names "
        "what it measures) plus the realized size of every other compartment, and let the reader see "
        "the whole cell rather than a scalar that cannot"
    ),
}


def _reject_overclaiming_census(config: Mapping[str, Any], census: Mapping[str, Any]) -> None:
    """Raise if the census claims more than it measures, or hides a knob that resized the population.

    Two failures, both observed on real records rather than imagined:

    1. A key whose NAME promises a whole-cell quantity it does not compute (:data:`BANNED_CENSUS_KEYS`).
    2. A population knob set in ``config`` whose realized size is absent from ``census``
       (:data:`POPULATION_KNOBS`), so a reader cannot tell that a compartment was shrunk.

    Refusing is the point.  Both defects PASS every arithmetic check — the numbers are internally
    consistent — so nothing short of a guard on the vocabulary catches them.

    Args:
        config: The run's configuration mapping.
        census: The population census the driver assembled.

    Raises:
        ValueError: On either failure, naming the field and what to record instead.
    """
    for key, why in BANNED_CENSUS_KEYS.items():
        if key in census:
            raise ValueError(f"census key {key!r} may not be recorded: {why}")
    missing = [(knob, need) for knob, need in POPULATION_KNOBS.items()
               if config.get(knob) is not None and need not in census]
    if missing:
        detail = "; ".join(f"config sets {k!r} but census has no {v!r}" for k, v in missing)
        raise ValueError(
            f"a knob that resizes a compartment must have its realized size in the census — {detail}. "
            "Without it the record cannot distinguish a full-population run from a reduced one, which "
            "is exactly how a quarter-resolution membrane was recorded as full native"
        )


def observation_artifact(
    *,
    run_label: str,
    evidence: EvidenceLabel,
    config: Mapping[str, Any],
    census: Mapping[str, Any],
    t0: Mapping[str, Any],
    timing: Mapping[str, Any],
    measurements: Mapping[str, Any],
    device: str,
    parameter_provenance: Mapping[str, str] | None = None,
    void_ceiling: VoidCeiling | None = None,
    gate_passed: bool | None = None,
    residual: float | None = None,
    signal: float | None = None,
    force_channels: Mapping[str, Any] | None = None,
    notes: Mapping[str, Any] | None = None,
    declared_commit: str | None = None,
) -> dict[str, Any]:
    """Assemble one driver's run record, stamped and classified.

    Args:
        run_label: Free-text identity of this run.
        evidence: The two-axis label, whose ``basis`` must name what was measured.
        config: Everything needed to reproduce the run; it is hashed into ``config_sha256``.
        census: Population census — unique active IDs per component, node/head/state totals.
        t0: The state measured BEFORE the first step.  Required: a run reporting only its end state
            cannot demonstrate that anything emerged.
        timing: The :func:`timing_block`.  Required, for the mirror-image reason: a run that records
            only its physics cannot say what it cost, and this repo spent 2026-07 unable to answer
            "is the new engine faster" because no artifact carried a wall-clock.  It is stamped here
            with the run's own verdict and marked non-comparable when that verdict does not support a
            speed comparison — see the ``comparable`` key.
        measurements: The observables this run produced, keyed by name.
        device: The device string the run executed on, recorded verbatim.
        parameter_provenance: Per-parameter ``SOURCED`` / ``DERIVED`` / ``CONVENIENCE`` / ``PI_GAP``
            classification, so a magnitude standing on an unsourced parameter is visible in the record.
        void_ceiling: The residual/signal ceiling declared before the run.  Required to classify a
            verdict.
        gate_passed: Whether the gate's own criterion was met, evaluated by the caller.
        residual: The solver residual, same units as ``signal``.
        signal: What the gate measures, same units as ``residual``.
        force_channels: A :func:`~aleph.engine.forces_manifest.dump_force_channels` dump, so the
            record carries which force channels were actually on.
        notes: Anything else worth recording.
        declared_commit: The build commit, for runs on a machine whose tree is not a git checkout —
            recorded as a declaration, never as a verified stamp.

    Returns:
        The JSON-able record.

    Raises:
        TypeError: If ``evidence`` is not an :class:`EvidenceLabel`.
        ValueError: If ``t0`` is empty, or if a verdict is partially specified — all four of
            ``void_ceiling`` / ``gate_passed`` / ``residual`` / ``signal`` must be given together or
            none of them, because a verdict built from a missing ceiling is exactly the unclassified
            ``"PASS"`` string this record replaces.
    """
    if not isinstance(evidence, EvidenceLabel):
        raise TypeError("evidence must be an EvidenceLabel (two axes), not a free string")
    if not t0:
        raise ValueError(
            "t0 is required: a run that records only its end state cannot show that anything emerged"
        )
    _reject_overclaiming_census(config, census)
    if not timing or "wall_seconds" not in timing:
        raise ValueError(
            "timing is required and must carry wall_seconds: build it with timing_block(). A record "
            "without a cost cannot answer whether anything got faster, and reconstructing one from log "
            "file timestamps is how this repo lost that answer for a month"
        )

    verdict_inputs = (void_ceiling, gate_passed, residual, signal)
    supplied = [item is not None for item in verdict_inputs]
    if any(supplied) and not all(supplied):
        raise ValueError(
            "a gate verdict needs all of void_ceiling / gate_passed / residual / signal, or none; "
            f"got {dict(zip(('void_ceiling', 'gate_passed', 'residual', 'signal'), supplied))}"
        )

    verdict: GateVerdict | None = None
    verdict_block: dict[str, Any] | None = None
    if all(supplied):
        assert void_ceiling is not None and gate_passed is not None  # narrowed by the check above
        verdict = classify_gate(
            passed=bool(gate_passed), residual=float(residual), signal=float(signal),
            ceiling=void_ceiling,
        )
        verdict_block = {
            "verdict": verdict.value,
            "criterion_met": bool(gate_passed),
            "residual": float(residual),
            "signal": float(signal),
            "residual_over_signal": (float(residual) / float(signal)) if signal else None,
            "void_ceiling_ratio": void_ceiling.ratio,
            "void_ceiling_rationale": void_ceiling.rationale,
        }

    # A timing is a property of an engine AT A STATED CONVERGENCE, so it is stamped with this run's own
    # verdict rather than left as a bare number a later reader can lift. Any solver is arbitrarily fast
    # if allowed to stop early, which is exactly how `outer_step_wall_s = 1.13` came to be readable next
    # to `inner_converged: false` without the two being connected.
    timing_record = dict(timing)
    if verdict is GateVerdict.PASS:
        timing_record["comparable"] = True
        timing_record["comparable_against"] = (
            "another run whose gate also PASSED under the same contract and observable; a cost is only "
            "a cost at a stated convergence"
        )
    else:
        timing_record["comparable"] = False
        timing_record["not_comparable_reason"] = (
            f"the gate verdict is {verdict.value if verdict is not None else 'absent'}, so this run did "
            "not establish the convergence its cost would have to be quoted at. Quoting it as the "
            "engine's speed compares the cost of NOT converging"
        )

    record: dict[str, Any] = {
        "schema": ARTIFACT_SCHEMA,
        "run_label": run_label,
        "build": build_stamp(declared_commit=declared_commit),
        "device": device,
        **evidence.as_artifact_fields(),
        "config": dict(config),
        "config_sha256": config_hash(config),
        "census": dict(census),
        "parameter_provenance": dict(parameter_provenance or {}),
        "t0": dict(t0),
        "timing": timing_record,
        "measurements": dict(measurements),
        "gate": verdict_block,
    }
    if force_channels is not None:
        record["force_channels"] = dict(force_channels)
    if notes:
        record["notes"] = dict(notes)
    return record


def write_artifact(path: Path | str, record: Mapping[str, Any]) -> Path:
    """Write a run record as indented JSON, creating the parent directory.

    Args:
        path: Destination file.
        record: The record from :func:`observation_artifact`.

    Returns:
        The written path.
    """
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(record, indent=2, default=str) + "\n")
    return destination

#!/usr/bin/env python
r"""The blebbistatin dose response — computed only from points that actually SETTLED.

WHY THIS IS A SCRIPT AND NOT A NOTEBOOK CELL.  The 2026-07-28 version of this sweep was read off six runs of a
FIXED 60 steps, came out non-monotonic in ``k_xb``, and was retired whole (``STATE.md`` (c) 16).  The cause was
not physics and not the engine: a softer crossbridge relaxes more slowly, so a fixed step count compares the
points at different places on different transients.  The protocol fix — run each point until its own series is
judged STATIONARY — only helps if something REFUSES to compute the ratio when a point did not settle.  That
refusal is this file's whole reason to exist, and it is why the check is a hard error rather than a warning.

WHAT IT REFUSES, and each of these retired a number here at least once:

* **Any point whose verdict is not STATIONARY.**  A mean that the detector declined to return is not a small
  mean; there is no mean.
* **A ratio taken across fewer than two decades of ``k_xb``.**  "gamma at k_xb -> 0" is a LIMIT claim, and two
  adjacent points do not establish one.
* **A sweep whose points ran DIFFERENT CODE**, unless the caller says so — the same physics on two builds is a
  claim to be measured, not assumed, which is what the cross-build parity point is for.  "Different code" means
  a different closure digest, NOT a different commit: ``build_commit`` is the repo HEAD at launch, so two runs
  of byte-identical code carry different commits whenever an unrelated file was committed between them.  This
  check refused its first two real points on exactly that false difference, which is why the run index now
  records ``closure_sha256`` and this reads it when it is there.

WHAT R IS, and what it is not.  ``R = gamma(k_xb -> 0) / gamma(k_xb = reference)`` is the fraction of cortical
tension that SURVIVES when the crossbridge is made compliant — the simulation's analogue of the blebbistatin
experiment, in which myosin II is inhibited and the residual tension is measured.  Experiment puts it at roughly
0.1-0.5 (i.e. tension is 50-90% myosin-dependent).  R is a RATIO of one observable to itself under one parameter
change, so it survives the per-parameter PI-GAPs that block quoting either gamma alone — but it does NOT survive
the force-accept caveat: these are stationary trajectories of a force-accepted integration, and that is what the
record says.

Usage:
    python aleph/scripts/ac_bleb_dose_response.py --dir outputs/ac/bleb_settled [--reference 1000]

Sanity Gate:
    * dimensional: gamma is [pN/µm] throughout; R is dimensionless by construction (a ratio of the same
      observable), which is exactly why it can be compared to an experimental band that this engine's absolute
      magnitudes cannot.
    * boundary: a sweep with one point, or with the reference missing, exits non-zero with the reason — it does
      not fall back to a partial ratio.
    * conservation/invariant: the bound-head fraction is reported alongside every point, because a dose response
      that moved the BINDING rather than the force transmission is a different claim, and the 2026-07-28 reading
      turned on exactly that control.
    * numerical: the uncertainty on R is propagated from the two sems in quadrature on the log, not eyeballed.
    * sign-sense: R < 1 means tension FELL as the crossbridge softened; R > 1 is reported as-is rather than
      clipped, because it would mean the sweep is not measuring what its name says.
    * measurement-protocol: each point's stationarity contract is echoed into the output, so the ratio is never
      seen apart from the criterion its inputs were admitted under.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

#: Reported blebbistatin residual cortical tension, as a FRACTION of untreated.  Used only to say whether the
#: computed R lands inside it — never to tune anything, and never as a pass condition without a signed gate.
EXPERIMENTAL_BAND = (0.1, 0.5)

#: A limit claim needs the swept parameter to move by at least this many decades.
MIN_DECADES = 2.0


def _closure_digest_for(record_path: Path) -> str:
    """Return the run index's recorded closure digest for this artifact, or ``""`` if it has none.

    The digest lives in the launcher's index rather than the run record because the driver that writes the
    record does not know its own import closure — the launcher resolves it in order to verify the remote
    tree, and recording what it already computed costs nothing.
    """
    index = Path.home() / ".ffn" / "runs.jsonl"
    if not index.is_file():
        return ""
    tail = str(record_path).split("outputs/", 1)[-1]
    best = ""
    for line in index.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if tail and tail in str(row.get("artifact", "")):
            best = str(row.get("closure_sha256", "")) or best
    return best


def _recorded_f_stall(record: dict) -> float | None:
    """Return the `f_stall` THIS run used, or ``None`` if the record only carries a provenance label.

    Records written before 2026-07-29 stored ``parameter_provenance = {"f_stall": "PI_GAP"}`` — the label
    without the number — so their working-stroke validity can only be judged against whatever the constant
    happens to be today.  That is a real difference, not a formality: `f_stall` is a PI-GAP, so it is exactly
    the kind of value that moves.  The driver now records ``{provenance, value, units}``; this reads the
    value when it is there and says so when it is not.
    """
    entry = (record.get("parameter_provenance") or {}).get("f_stall")
    if isinstance(entry, dict) and entry.get("value") is not None:
        return float(entry["value"])
    return None


def load_points(directory: Path) -> list[dict]:
    """Return one entry per sweep point, sorted by ``k_xb`` descending.

    Args:
        directory: Directory of ``kxb_*/record.json`` run records.

    Returns:
        Dicts with ``k_xb``, the stationarity report for each observable, the build, and the config.

    Raises:
        FileNotFoundError: If no record is found — an empty sweep must not read as a passing one.
    """
    points: list[dict] = []
    for record_path in sorted(directory.glob("*/record.json")):
        record = json.loads(record_path.read_text())
        config = record["config"]
        points.append({
            "name": record_path.parent.name,
            "k_xb": float(config["k_xb"]),
            "steps_run": int(config.get("steps_run", 0)),
            "build": record["build"].get("commit", "unknown"),
            "engine_cortex": bool(config.get("engine_cortex", False)),
            "closure_sha256": _closure_digest_for(record_path),
            "stationarity": record["measurements"]["stationarity"],
            "driving_tau_s": float(record["measurements"].get("bound_head_lifetime_s", float("nan"))),
            "f_stall_pn": _recorded_f_stall(record),
            "path": str(record_path),
        })
    if not points:
        raise FileNotFoundError(f"no sweep records under {directory}")
    return sorted(points, key=lambda p: -p["k_xb"])


#: Physiological myosin working stroke [nm].  NOT a threshold chosen here: it is the range
#: `ac/motor/minifilament_topology.working_stroke_strain` names in its own docstring, in the function that
#: calls itself "the k_xb MASTER-knob arbiter".  A swept `k_xb` sets the strain a head must take on to bear
#: its stall force (`f_stall / k_xb`); once that exceeds the head's own geometry the stall collapses and the
#: run is no longer simulating a myosin head.
WORKING_STROKE_NM = (5.0, 20.0)


def _outside_working_stroke(points: list[dict]) -> list[tuple[dict, float]]:
    """Return points whose implied working stroke leaves the physiological window, with the strain in nm.

    Measured 2026-07-29 on the five-point sweep, with ``f_stall`` = 0.5 pN/head:

        k_xb  1000 -> 0.5 nm · 100 -> 5.0 nm · 10 -> 50 nm · 3 -> 167 nm · 1 -> 500 nm

    Only ``k_xb=100`` lands inside 5-20 nm, and gamma's exponent flips sign (-0.08 -> +1.20) between the
    50 nm and 167 nm points.  At ``k_xb=1`` the strain (500 nm) exceeds the head offset itself (200 nm) —
    literally the collapse the arbiter's docstring predicts.  The two-regime shape is therefore a property
    of running the model outside its stated domain, not a dose response the model is predicting.

    Args:
        points: De-duplicated sweep points.

    Returns:
        ``(point, strain_nm)`` for each point outside the window, softest first.
    """
    from aleph.components.incumbent.assemble import NMII_F_STALL_TEST
    from aleph.components.motor.minifilament_topology import working_stroke_strain

    # LIMITATION, stated because it is not fixable from here: `f_stall` is read from the CURRENT build, not
    # from the run.  The record names it in `parameter_provenance` as "PI_GAP" and stores no number, so a run
    # made under a different `f_stall` would be judged against today's.  Recording the VALUE beside the
    # provenance label is the fix, and it belongs in the driver's artifact block, not in this analyser.

    low, high = WORKING_STROKE_NM
    out: list[tuple[dict, float]] = []
    for point in sorted(points, key=lambda p: p["k_xb"]):
        f_stall = point.get("f_stall_pn")
        if f_stall is None:                       # pre-2026-07-29 record: label without a number
            f_stall = NMII_F_STALL_TEST
            print(f"[bleb] {point['name']}: record carries no f_stall VALUE, only its provenance label — "
                  f"judging its working stroke against the current build's {f_stall:g} pN. If the run used "
                  f"a different one, this verdict is about today's constant, not about that run.",
                  file=sys.stderr)
        nm = working_stroke_strain(f_stall, point["k_xb"]) * 1e3
        if not low <= nm <= high:
            out.append((point, nm))
    return out


#: How far gamma must RISE toward softer `k_xb`, in combined sems, before it counts as a turning point
#: rather than noise.  Declared 2026-07-29 alongside PARITY_SIGMA and PLATEAU_SIGMA.
MONOTONIC_SIGMA = 3.0


def _monotonicity_breaks(points: list[dict]) -> list[tuple[dict, dict, float]]:
    """Return adjacent pairs where gamma ROSE as ``k_xb`` fell — the premise R needs, violated.

    ``R = gamma(k_xb -> 0) / gamma(reference)`` is a LIMIT claim, and every reading of it ("the fraction of
    tension that survives a compliant crossbridge") presumes gamma descends toward that limit.  A turning
    point inside the swept range breaks the premise, not merely the number: the softest point stops being an
    extreme, so dividing by it answers no question anyone asked.

    Measured 2026-07-29 and the reason this exists: gamma went 4.3256 -> 3.4807 -> 2.8839 -> **12.2425** at
    ``k_xb`` 1000 -> 100 -> 10 -> 3.  Before this check the analyser would have reported R = 2.83 with a
    plateau note and no objection, because every other guard it had was about whether the points SETTLED.

    Args:
        points: De-duplicated sweep points.

    Returns:
        ``(softer, stiffer, ratio)`` for each adjacent pair whose rise clears ``MONOTONIC_SIGMA``.
    """
    ordered = sorted(points, key=lambda p: -p["k_xb"])
    breaks: list[tuple[dict, dict, float]] = []
    for stiffer, softer in zip(ordered, ordered[1:]):
        a = stiffer["stationarity"]["gamma_total_pn_per_um"]
        b = softer["stationarity"]["gamma_total_pn_per_um"]
        if a["mean"] is None or b["mean"] is None or not a["mean"]:
            continue
        rise = b["mean"] - a["mean"]
        combined = math.hypot(a["sem"] or 0.0, b["sem"] or 0.0)
        if rise > 0 and combined and rise / combined > MONOTONIC_SIGMA:
            breaks.append((softer, stiffer, b["mean"] / a["mean"]))
    return breaks


def refuse_unless_admissible(points: list[dict], reference: float) -> list[str]:
    """Return the reasons this sweep may NOT yield a ratio; empty means it may.

    Args:
        points: Loaded sweep points.
        reference: The ``k_xb`` treated as untreated myosin.

    Returns:
        One string per blocking reason, in the order a reader should act on them.
    """
    reasons: list[str] = []
    unsettled = [(p["name"], p["stationarity"]["gamma_total_pn_per_um"]["verdict"])
                 for p in points if p["stationarity"]["gamma_total_pn_per_um"]["verdict"] != "STATIONARY"]
    for name, verdict in unsettled:
        reasons.append(
            f"{name} is {verdict}, not STATIONARY — the detector declined to return a mean, and a mean it "
            f"declined to return is not a small mean. Run it longer; do not average the transient")
    if not any(math.isclose(p["k_xb"], reference) for p in points):
        reasons.append(f"the reference point k_xb={reference:g} is absent, so there is nothing to divide by")
    if len(points) >= 2:
        span = math.log10(max(p["k_xb"] for p in points) / min(p["k_xb"] for p in points))
        if span < MIN_DECADES:
            reasons.append(
                f"the sweep spans {span:.2f} decades of k_xb, under the {MIN_DECADES:g} a LIMIT claim needs — "
                f"'gamma as k_xb -> 0' cannot be read off adjacent points")
    else:
        reasons.append("a single point is not a dose response")
    outside = _outside_working_stroke(points)
    if outside:
        listed = ", ".join(
            f"k_xb={p['k_xb']:g} ({nm:.3g} nm, {'too STIFF' if nm < WORKING_STROKE_NM[0] else 'too SOFT'})"
            for p, nm in outside)
        reasons.append(
            f"{len(outside)} point(s) sit OUTSIDE the model's own validity window for k_xb: {listed}, against "
            f"the physiological working stroke {WORKING_STROKE_NM[0]:g}-{WORKING_STROKE_NM[1]:g} nm. This is not "
            f"a threshold invented here — `ac/motor/minifilament_topology.working_stroke_strain` calls itself "
            f"'the k_xb MASTER-knob arbiter' and states that k_xb=1 pN/um gives a strain larger than the whole "
            f"minifilament, so the stall geometry collapses. A dose response may not be read across points "
            f"where the mechanism has stopped being the mechanism")
    turns = _monotonicity_breaks(points)
    for softer, stiffer, ratio in turns:
        reasons.append(
            f"the dose response is NOT MONOTONIC: gamma at k_xb={softer['k_xb']:g} is {ratio:.2f}x the value "
            f"at k_xb={stiffer['k_xb']:g}, i.e. tension ROSE as the crossbridge softened. R is defined as a "
            f"LIMIT at k_xb -> 0 and that framing presumes gamma descends toward it; with a turning point in "
            f"the swept range there is no limit to quote and the softest point is not an extreme. Report the "
            f"turning point — it is the finding — and do not divide across it")
    digests = {p["closure_sha256"] for p in points if p["closure_sha256"]}
    if digests and len(digests) == len(  # every point has a digest and they agree
            [p for p in points if p["closure_sha256"]]) and len(digests) > 1:
        reasons.append(
            f"points ran {len(digests)} DIFFERENT import closures — the code itself differed, not just the "
            f"commit. Same physics on two closures is a claim to be MEASURED (that is what the cross-build "
            f"parity point is for). Pass --allow-mixed-builds once that parity is established")
    elif not digests:
        builds = {p["build"][:8] for p in points}
        if len(builds) > 1:
            reasons.append(
                f"points carry {len(builds)} different build commits ({sorted(builds)}) and NO closure digest, "
                f"so whether the code differed cannot be established. A commit differs whenever any unrelated "
                f"file was committed, so this may be a false difference — but with no digest it cannot be "
                f"ruled out. Submit through `gpu-submit`, stamp the record with `run_provenance`, or pass "
                f"--allow-mixed-builds if the equality is known by other means")
    return reasons


#: Agreement required of the two softest points before the softest may be read as the `k_xb -> 0` LIMIT.
#: Declared 2026-07-29 alongside PARITY_SIGMA, before R was first computed.
PLATEAU_SIGMA = 3.0

#: Agreement required of a cross-build parity pair, in combined sems. Declared BEFORE the parity run landed
#: (2026-07-29) so the criterion is not chosen after seeing the numbers — the failure this file exists to stop.
PARITY_SIGMA = 3.0


def split_parity_controls(points: list[dict]) -> tuple[list[dict], list[dict]]:
    """Separate genuine sweep points from cross-build parity re-runs of a point already present.

    A second record at the SAME ``k_xb`` is not a second dose, it is a control: the same physics asked of a
    different build.  Left in the list it would (a) make the reference ambiguous — two candidates to divide
    by — and (b) inflate the build-diversity count with a difference that is the very thing being measured.

    Args:
        points: All loaded records, sorted by descending ``k_xb``.

    Returns:
        ``(sweep, controls)``.  The FIRST record at each ``k_xb`` stays in ``sweep``; later ones become
        controls.  "First" is by the sort `load_points` already applied, which is stable on name, so the
        original ``kxb_1000`` precedes ``kxb_1000_rebuild`` and the sweep keeps the point it always had.
    """
    seen: set[float] = set()
    sweep: list[dict] = []
    controls: list[dict] = []
    for point in points:
        key = point["k_xb"]
        (controls if key in seen else sweep).append(point)
        seen.add(key)
    return sweep, controls


def assess_parity(sweep: list[dict], controls: list[dict]) -> list[dict]:
    """Judge each control against the sweep point it duplicates, under the pre-declared criterion.

    Args:
        sweep: The de-duplicated sweep points.
        controls: Re-runs at a ``k_xb`` already in ``sweep``.

    Returns:
        One verdict dict per control: the two means, their difference in combined sems, whether it agrees,
        and — the part that limits what may be concluded — which build commits the pair actually covers.
    """
    verdicts: list[dict] = []
    for control in controls:
        twin = next(p for p in sweep if math.isclose(p["k_xb"], control["k_xb"]))
        a = twin["stationarity"]["gamma_total_pn_per_um"]
        b = control["stationarity"]["gamma_total_pn_per_um"]
        if a["mean"] is None or b["mean"] is None:
            verdicts.append({"control": control["name"], "agrees": False,
                             "why": "one side has no mean — the detector declined it, so there is nothing to compare"})
            continue
        combined = math.hypot(a["sem"], b["sem"])
        delta = abs(a["mean"] - b["mean"])
        verdicts.append({
            "control": control["name"], "twin": twin["name"], "k_xb": control["k_xb"],
            "mean_twin": a["mean"], "mean_control": b["mean"], "delta": delta,
            "combined_sem": combined, "sigma": delta / combined if combined else float("inf"),
            "agrees": bool(combined) and delta <= PARITY_SIGMA * combined,
            "builds_covered": sorted({twin["build"][:8], control["build"][:8]}),
        })
    return verdicts


def dose_response(points: list[dict], reference: float) -> dict:
    """Return R and its propagated uncertainty from the settled points.

    Args:
        points: Loaded, admissible sweep points.
        reference: The ``k_xb`` treated as untreated myosin.

    Returns:
        The ratio block, including the softest point used as the ``k_xb -> 0`` limit and the band verdict.
    """
    by_k = {p["k_xb"]: p for p in points}
    ref = by_k[min(by_k, key=lambda k: abs(k - reference))]
    soft = by_k[min(by_k)]
    g_ref = ref["stationarity"]["gamma_total_pn_per_um"]
    g_soft = soft["stationarity"]["gamma_total_pn_per_um"]
    ratio = g_soft["mean"] / g_ref["mean"]
    # Relative errors add in quadrature for a quotient; the sems are already the correlated-series ones.
    rel = math.hypot(g_soft["sem"] / g_soft["mean"], g_ref["sem"] / g_ref["mean"])
    low, high = EXPERIMENTAL_BAND
    # Is the softest point a LIMIT, or just the softest we could afford?  R is defined at k_xb -> 0, so the
    # band comparison is only meaningful once gamma has stopped moving.  Criterion declared here rather than
    # judged by eye: the two softest points must agree within PLATEAU_SIGMA of their combined sem.
    ordered = sorted(points, key=lambda p: p["k_xb"])
    plateaued, plateau_note = False, "only one point — a plateau cannot be tested"
    if len(ordered) >= 2:
        g_next = ordered[1]["stationarity"]["gamma_total_pn_per_um"]
        step = abs(g_soft["mean"] - g_next["mean"])
        step_sem = math.hypot(g_soft["sem"], g_next["sem"])
        sigmas = step / step_sem if step_sem else float("inf")
        plateaued = sigmas <= PLATEAU_SIGMA
        plateau_note = (
            f"gamma moved {step:.4f} pN/um between the two softest points ({ordered[1]['k_xb']:g} -> "
            f"{soft['k_xb']:g}) = {sigmas:.0f} sigma, against the {PLATEAU_SIGMA:g} declared for a plateau")
    return {
        "reference_k_xb": ref["k_xb"], "softest_k_xb": soft["k_xb"],
        "gamma_reference": g_ref["mean"], "gamma_reference_sem": g_ref["sem"],
        "gamma_softest": g_soft["mean"], "gamma_softest_sem": g_soft["sem"],
        "R": ratio, "R_sem": ratio * rel,
        "experimental_band": list(EXPERIMENTAL_BAND),
        "softest_point_is_a_limit": plateaued,
        "plateau_note": plateau_note,
        "inside_band": bool(low <= ratio <= high) if plateaued else None,
        "band_verdict": (
            ("INSIDE" if low <= ratio <= high else "OUTSIDE") if plateaued
            else "INCONCLUSIVE — R is an UPPER BOUND, not a limit"),
        "reading": (
            "R is the fraction of cortical tension that SURVIVES a compliant crossbridge. Inside the band means "
            "the model's myosin dependence resembles the experiment; ABOVE it means the model's tension is "
            "carried by something other than the motor, which is a reportable finding and not a failure to fix "
            "by tuning. BUT all of that presumes the softest point IS the k_xb -> 0 limit. If gamma is still "
            "falling there, the measured ratio only BOUNDS R from above and the band is neither confirmed nor "
            "excluded — quoting INSIDE/OUTSIDE off a truncated sweep is the same error as reading a mean off a "
            "transient"),
        "not_quotable": (
            "as an absolute tension: every NMII parameter is a PI-GAP. And these are stationary trajectories of "
            "a FORCE-ACCEPTED integration, so the ratio is a property of that trajectory, not of a mechanically "
            "balanced cortex"),
    }


#: How close to the pass threshold counts as "the verdict could have gone either way".
MARGINAL_DRIFT = 0.9


def _disclose_marginal_and_drifting(points: list[dict]) -> None:
    """Print what the STATIONARY/not-STATIONARY label hides — never blocking, always visible.

    Two things the table alone cannot say, both measured on the 2026-07-29 sweep:

    * **A pass can be marginal.**  At ``k_xb=10`` the drift was 0.98 sigma against a 1.0 threshold, so the
      label was decided by 2%.  Three identical "STATIONARY" cells invite the reader to treat the points as
      equally settled when one of them nearly was not.
    * **A sibling observable can have failed.**  The ratio is a ``gamma_total`` ratio, so a drifting
      ``gamma_source`` does not make R wrong and is deliberately NOT a blocking reason — but at ``k_xb=10``
      the detector DECLINED TO RETURN A MEAN for it, and a refusal that never reaches the reader is the
      same defect this whole file exists to prevent, one level down.

    Reporting drift/sigma for every point also shows the more useful fact: at ``k_xb=10`` the three
    observables sat at 0.92 / 0.98 / 1.18 sigma, i.e. the threshold sorted a set of physically similar
    states into pass and fail.  That is a property of the criterion, not of the physics.
    """
    for point in points:
        stat = point["stationarity"]
        gamma = stat["gamma_total_pn_per_um"]
        drift = gamma.get("drift_over_std")
        if drift is not None and gamma["verdict"] == "STATIONARY" and drift >= MARGINAL_DRIFT:
            print(f"[bleb] {point['name']}: gamma_total passed at {drift:.2f} sigma against a 1.0 threshold "
                  f"— a marginal pass, not a comfortable one")
        others = [(name, obs.get("drift_over_std")) for name, obs in stat.items()
                  if name != "gamma_total_pn_per_um" and obs.get("verdict") != "STATIONARY"]
        for name, other_drift in others:
            shown = f"{other_drift:.2f} sigma" if other_drift is not None else "no drift recorded"
            print(f"[bleb] {point['name']}: {name} did NOT settle ({shown}); the detector returned no mean "
                  f"for it. R is a gamma_total ratio so this does not block it — but the point is not "
                  f"'settled' without qualification")


def main() -> int:
    """Load, refuse or compute, and print."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--dir", default="outputs/ac/bleb_settled")
    parser.add_argument("--reference", type=float, default=1000.0,
                        help="the k_xb standing for untreated myosin")
    parser.add_argument("--allow-mixed-builds", action="store_true",
                        help="permit points from different builds; only after cross-build parity is measured")
    parser.add_argument("--out", default=None, help="write the ratio block here as JSON")
    args = parser.parse_args()

    all_records = load_points(Path(args.dir))
    points, controls = split_parity_controls(all_records)
    parity = assess_parity(points, controls)
    print(f"{'point':>12}{'k_xb':>9}{'steps':>7}{'verdict':>12}{'drift/sig':>10}"
          f"{'gamma_total':>13}{'sem':>9}{'bound':>8}")
    for point in points:
        gamma = point["stationarity"]["gamma_total_pn_per_um"]
        bound = point["stationarity"]["bound_fraction"]
        mean = f"{gamma['mean']:.4f}" if gamma["mean"] is not None else "—"
        sem = f"{gamma['sem']:.5f}" if gamma["sem"] is not None else "—"
        bfr = f"{bound['mean']:.4f}" if bound.get("mean") is not None else "—"
        drift = gamma.get("drift_over_std")
        dstr = f"{drift:.2f}" if drift is not None else "—"
        print(f"{point['name']:>12}{point['k_xb']:>9g}{point['steps_run']:>7}"
              f"{gamma['verdict']:>12}{dstr:>10}{mean:>13}{sem:>9}{bfr:>8}")
    _disclose_marginal_and_drifting(points)

    build_parity_established = False
    for verdict in parity:
        if "why" in verdict:
            print(f"[bleb] PARITY {verdict['control']}: {verdict['why']}")
            continue
        mark = "AGREE" if verdict["agrees"] else "DISAGREE"
        print(f"[bleb] PARITY {mark} — {verdict['twin']} vs {verdict['control']} at k_xb="
              f"{verdict['k_xb']:g}: {verdict['mean_twin']:.4f} vs {verdict['mean_control']:.4f}, "
              f"delta {verdict['delta']:.4f} = {verdict['sigma']:.2f} sigma against the {PARITY_SIGMA:g} "
              f"declared before the run. Builds compared: {', '.join(verdict['builds_covered'])}")
        build_parity_established |= verdict["agrees"]

    reasons = refuse_unless_admissible(points, args.reference)
    if build_parity_established and not args.allow_mixed_builds:
        covered = {b for v in parity if v.get("agrees") for b in v["builds_covered"]}
        used = {p["build"][:8] for p in points}
        uncovered = sorted(used - covered)
        before = len(reasons)
        reasons = [r for r in reasons if "different build commits" not in r and "DIFFERENT import closures" not in r]
        if len(reasons) < before:
            print("[bleb] the build-difference block is LIFTED BY MEASUREMENT, not by assertion — a re-run of "
                  "the same point on another build reproduced it inside the declared tolerance.")
            if uncovered:
                print(f"[bleb] LIMIT OF THAT EVIDENCE: the parity pair covers {', '.join(sorted(covered))}. "
                      f"Builds {', '.join(uncovered)} carry sweep points and were NOT re-run, so their "
                      f"equivalence is INFERRED from one pair spanning the range, not measured. Say so "
                      f"wherever R is quoted.")
    if args.allow_mixed_builds:
        reasons = [r for r in reasons if "different builds" not in r and "different build commits" not in r]
    if reasons:
        print(f"\n[bleb] NO RATIO — {len(reasons)} blocking reason(s):", file=sys.stderr)
        for reason in reasons:
            print(f"  {reason}", file=sys.stderr)
        return 1

    block = dose_response(points, args.reference)
    low, high = block["experimental_band"]
    print(f"\n[bleb] R = gamma(k_xb={block['softest_k_xb']:g}) / gamma(k_xb={block['reference_k_xb']:g}) "
          f"= {block['R']:.4f} +/- {block['R_sem']:.4f}")
    print(f"[bleb] plateau test: {block['plateau_note']}")
    print(f"[bleb] experimental band {low}-{high} -> {block['band_verdict']}")
    if block["softest_point_is_a_limit"]:
        if not block["inside_band"]:
            print("[bleb] OUTSIDE is a RESULT, not a failure: it says the model's cortical tension is carried "
                  "by something other than the motor. Report it; do not tune toward the band.")
    else:
        print(f"[bleb] The softest point is NOT a limit, so R = {block['R']:.4f} BOUNDS R(k_xb->0) FROM ABOVE "
              f"and nothing follows about the band in either direction. Reaching the limit needs softer k_xb, "
              f"which this sweep measured to be unaffordable — that is the finding to report, not a band "
              f"verdict the data cannot support.")
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps({"points": points, "dose_response": block}, indent=2) + "\n")
        print(f"[bleb] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

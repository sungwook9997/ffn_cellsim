"""What does a converged native step COST — and if it does not converge, how far off is it?

PI decision 2026-07-28 ("둘 다"): measure first so the preconditioner has a numeric target, then build it.

WHY THIS RUN EXISTS.  The ratified roadmap's stage 2 says "the cost of a converged native step is
MEASURED via ``timing_block``", and the whole stage-4 sweep budget (390,000 runs -> ~100-300) rests on
that number.  But there is no converged native run to time: ``RESTING_BASELINE_DIAGNOSIS_2026-07-22e``
records the native resting solve plateauing at 0.620 pN against a derived gate floor near 0.21 pN, with
the cause refined to backbone-inextensibility dominance rather than sparse-crosslink conditioning.

So this driver deliberately does NOT claim to measure the cost of convergence.  It measures the two
things that are actually available and that together bound it:

1. **The shape of the failure.**  ``max|PF|`` versus inner iteration, captured device-side.  Every
   resting record in the repo carries ``residual_start`` and ``residual_end`` and nothing between, and
   those two numbers cannot tell "converging slowly" from "stalled" — which is exactly the distinction
   the cost of a converged step depends on.  A run that is merely slow has a budget; a stalled one has
   no finite cost at all, and quoting a per-iteration price for it is quoting the cost of not converging.
2. **The per-iteration price**, so that once the preconditioner restores a decay rate, the converged
   cost is arithmetic instead of a new GPU-day.

The extrapolation is self-invalidating BY CONSTRUCTION.  A decay rate is fitted on the tail, and if the
tail is flat (the stall signature) the projected iteration count is withheld and the record says
``STALLED`` rather than emitting a large finite number that reads like a budget.  This matters because
the failure mode being guarded against is precisely a plausible number surviving into a plan.

WHAT THIS RUN CANNOT SAY.  Nothing about tension, and no magnitude beyond the residual trajectory
itself.  ``QuantitativeClaim`` is BLOCKED and the gate verdict will not be PASS, so ``timing_block``
stamps the cost ``comparable: false`` on its own.  That is the correct outcome, not a shortfall.

CUDA-gated -> gbook A5000, native population only.  Clear the Warp kernel cache before running if any
``@wp.kernel`` changed (see reference-gbook-warp-kernel-cache-stale: a stale cache produced a
byte-identical false negative once already).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

# NOTE: `aleph.ac.*`, not the bare `ac.*` the older drivers use. STATE.md (a) records 22 scripts
# entering the engine under the bare identity, which double-loads every module outside pytest; this is
# a new file, so it does not add a 23rd.
from aleph.components.incumbent.assemble import CellConfig
from aleph.components.incumbent.driver import run_from_resting
from aleph.engine.contracts import (
    EvidenceLabel,
    EvidenceRung,
    QuantitativeClaim,
    VoidCeiling,
)
from aleph.engine.gate_criteria import predicted_force_floor
from aleph.engine.observe.artifact import (
    observation_artifact,
    timing_block,
    write_artifact,
)

#: Population of the first baseline (100 um^-2 x 4*pi*(7.5 um)^2). Never lowered to fit memory.
NATIVE_FILAMENTS = 70_686

#: Recording cadence for the device-side convergence checkpoint. Also the solve's reshape cadence, so
#: it cannot be set independently — a denser curve costs remeshes, which would change the trajectory
#: being measured. 20 is the incumbent default and is left alone for exactly that reason.
RECORD_EVERY = 20

#: Resting myosin seed. Both values are PI-GAP (recorded, not sourced) and are carried over verbatim
#: from `ac_resting_converge.py` so this run measures the SAME configuration whose stall is on record.
#: Changing them would measure a different cell and silently invalidate the comparison.
FRACTION = 0.5
FORCE_PN = 1.5
CAPTURE_UM = 0.6
SOLVER = "erm_jacobi_pure"

#: Declared BEFORE the run (contracts.VoidCeiling exists to make that ordering explicit in code).
VOID_CEILING = VoidCeiling(
    ratio=0.5,
    rationale=(
        "residual here is the plateau max|PF| and signal is the initial projected load the solve was "
        "given. A solve that leaves more than half of that load standing has not measured a converged "
        "state at all, so its per-iteration cost describes iterations that are not making progress. "
        "Declared before the run; the recorded plateau of the same configuration (0.620 pN against a "
        "residual_start near 7.85 pN) was NOT used to choose it, which would have been ceiling-fitting."
    ),
)


def _tail_decay(iters: list[int], pf: list[float], tail_fraction: float = 0.5) -> dict[str, Any]:
    """Fit ``log(max|PF|)`` linearly against iteration over the trajectory tail.

    The tail rather than the whole curve because the early transient is fast and would flatter the fit:
    an initial drop followed by a plateau has a healthy-looking global slope while being, in the only
    sense that matters here, stalled.

    Args:
        iters: Iteration index at each recorded checkpoint.
        pf: ``max|PF|`` [pN] at each checkpoint.
        tail_fraction: Fraction of the trajectory, measured from the end, to fit.

    Returns:
        The fit, plus ``regime`` in ``{"DECAYING", "STALLED", "DIVERGING", "INSUFFICIENT"}``.  A
        per-iteration decay is only reported when the regime is ``DECAYING``; the point of the other
        three labels is that a projected iteration count must not be computed from them.
    """
    n = len(iters)
    if n < 6:
        return {"regime": "INSUFFICIENT", "reason": f"only {n} checkpoint(s); need >= 6 to fit a tail"}

    start = int(n * (1.0 - tail_fraction))
    x = np.asarray(iters[start:], dtype=np.float64)
    y = np.asarray(pf[start:], dtype=np.float64)
    if not np.all(np.isfinite(y)) or np.any(y <= 0.0):
        return {"regime": "INSUFFICIENT", "reason": "tail contains a non-positive or non-finite max|PF|"}

    slope, intercept = np.polyfit(x, np.log(y), 1)
    # Relative change the fitted line predicts across the whole fitted window. A window that cannot
    # even halve the residual is a plateau however cleanly it fits; the threshold is stated, not tuned
    # to an outcome, and both branches are recorded so the reader can re-judge.
    window = float(x[-1] - x[0])
    predicted_ratio = float(math.exp(slope * window))

    fit: dict[str, Any] = {
        "tail_fraction": tail_fraction,
        "n_points_fitted": int(len(x)),
        "iteration_window": [int(x[0]), int(x[-1])],
        "log_slope_per_iteration": float(slope),
        "predicted_ratio_across_window": predicted_ratio,
    }
    if slope > 0.0:
        fit["regime"] = "DIVERGING"
        fit["reason"] = "max|PF| grows across the fitted tail; no convergence budget exists"
    elif predicted_ratio > 0.5:
        fit["regime"] = "STALLED"
        fit["reason"] = (
            f"the fitted tail removes only {(1.0 - predicted_ratio) * 100:.2f}% of max|PF| across "
            f"{int(window)} iterations. Extrapolating an iteration count from a slope this shallow "
            "would produce a large finite number that reads like a budget; it is withheld instead"
        )
    else:
        fit["regime"] = "DECAYING"
        fit["half_life_iterations"] = float(math.log(0.5) / slope)
    return fit


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--native", action="store_true",
                    help="full population (70,686 filaments). Without it the run is a slice and its "
                         "numbers are non-authoritative by the develop-on-a-slice/conclude-at-native rule")
    ap.add_argument("--n-inner", type=int, default=15000,
                    help="inner iteration budget; the curve is recorded every RECORD_EVERY of these")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--build-commit", default=None,
                    help="the build this run is measured on, supplied by the enforced `ffn_gpu.py run` "
                         "path. Recorded as a DECLARATION, never a verified stamp: the execution host's "
                         "tree is a Syncthing mirror rather than a git checkout, so nothing on that "
                         "machine can verify it")
    args = ap.parse_args()

    n = NATIVE_FILAMENTS if args.native else 8000
    out = args.out or (
        Path(__file__).resolve().parents[1] / "outputs" / "ac" / "resting_native" / "resting_cost_curve.json"
    )

    cfg = CellConfig(
        n_filaments=n,
        overlap_free_cortex=True,
        erm_radial_pairing=True,
        membrane_subdivisions=6,
        resting_bound_myosin_fraction=FRACTION,
        resting_bound_myosin_force_pn=FORCE_PN,
        resting_bound_myosin_source="resting_cost_curve TEST (PI-GAP, carried from ac_resting_converge)",
        resting_bound_myosin_capture_um=CAPTURE_UM,
    )

    print(f"=== resting cost curve  n={n}  solver={SOLVER}  n_inner={args.n_inner} ===")
    print(f"    recording max|PF| every {RECORD_EVERY} iterations (device-side, no hot-loop host read)")
    t0 = time.perf_counter()
    r = run_from_resting(
        cfg,
        n_inner=args.n_inner,
        reshape_every=RECORD_EVERY,
        inner_solver=SOLVER,
        preload_erm_balance=True,
        capture_convergence_history=True,
    )
    wall = time.perf_counter() - t0

    history = r.ledger.get("convergence_history")
    if not history or not history.get("iteration"):
        raise SystemExit(
            "no convergence history was captured — run_from_resting did not plumb capture_candidate "
            "through to the solve, so this run measured nothing it was written to measure"
        )

    iters: list[int] = history["iteration"]
    pf: list[float] = history["max_projected_force_pn"]

    # D4: the floor comes from the run's OWN executing predicate, never the frozen 0.21 literal, which
    # is proportional to mesh spacing and so silently loosens on every finer rung.
    gate_floor = predicted_force_floor(
        inner_tolerance_um=r.inner_tolerance_um, inner_dt_mu=r.inner_dt_mu
    )

    # Plateau = median of the last decade of checkpoints, which is robust to the oscillation the
    # projected-descent step shows near its floor in a way the single final sample is not.
    tail_n = max(1, len(pf) // 10)
    plateau = float(np.median(np.asarray(pf[-tail_n:], dtype=np.float64)))
    fit = _tail_decay(iters, pf)

    total_inner = int(r.inner_iters) or int(args.n_inner)
    per_iteration_s = wall / total_inner

    cost: dict[str, Any] = {
        "regime": fit["regime"],
        "wall_seconds_per_inner_iteration": per_iteration_s,
        "gate_floor_pn": gate_floor,
        "gate_floor_basis": "predicted_force_floor(inner_tolerance_um, inner_dt_mu) — D4, run's own predicate",
        "plateau_max_pf_pn": plateau,
        "factor_above_gate_floor": plateau / gate_floor if gate_floor else None,
    }
    if fit["regime"] == "DECAYING":
        remaining = math.log(gate_floor / plateau) / fit["log_slope_per_iteration"]
        cost["projected_iterations_to_gate"] = float(remaining)
        cost["projected_converged_step_seconds"] = float(remaining) * per_iteration_s
        cost["projection_basis"] = "fitted tail decay rate held constant to the derived gate floor"
    else:
        cost["projected_iterations_to_gate"] = None
        cost["projected_converged_step_seconds"] = None
        cost["projection_withheld_because"] = fit.get("reason", fit["regime"])

    signal = float(r.residual_start)
    record = observation_artifact(
        run_label="resting_cost_curve",
        evidence=EvidenceLabel(
            rung=EvidenceRung.NATIVE if args.native else EvidenceRung.CUDA_UNIT,
            quantitative=QuantitativeClaim.BLOCKED,
            basis=(
                "recorded max|PF| against inner iteration on the full native resting configuration and "
                "the wall-clock of the same solve. BLOCKED because the solve does not reach its own "
                "derived force floor, so no state it produced is an equilibrium and no magnitude read "
                "from that state is a physical one — including, especially, its cost"
            ),
        ),
        config={
            "n_filaments": n,
            "solver": SOLVER,
            "n_inner": args.n_inner,
            "record_every": RECORD_EVERY,
            "resting_bound_myosin_fraction": FRACTION,
            "resting_bound_myosin_force_pn": FORCE_PN,
            "resting_bound_myosin_capture_um": CAPTURE_UM,
            "overlap_free_cortex": True,
            "erm_radial_pairing": True,
            "membrane_subdivisions": 6,
            "preload_erm_balance": True,
        },
        census={
            "cortex_filaments": n,
            "n_actin_nodes": int(r.n_actin),
            "n_total_nodes": int(r.n_total),
            # The membrane this run built, recorded because a subdivision knob that is invisible in the
            # census lets a resized compartment read as a full-population one (icosphere: 10·4^k + 2).
            "membrane_subdivisions": 6,
            "membrane_vertices": 10 * 4 ** 6 + 2,
            # NAMED for what it measures: this is the cortex fraction, and it moves with nothing else.
            "fraction_of_native_cortex": float(n) / float(NATIVE_FILAMENTS),
        },
        t0={
            "residual_start_pn": signal,
            "r_mean_start_um": float(r.r_mean_start),
            "max_pf_first_checkpoint_pn": float(pf[0]),
            "first_checkpoint_iteration": int(iters[0]),
        },
        timing=timing_block(
            wall_seconds=wall,
            n_steps=1,
            n_inner_iterations=total_inner,
            device_note=(
                "instrumented run, not a benchmark: capture_convergence_history allocates two extra "
                "position buffers and the device-side checkpoint kernel launches every "
                f"{RECORD_EVERY} iterations. Treat the per-iteration cost as an upper bound"
            ),
        ),
        measurements={
            "convergence_history": history,
            "tail_fit": fit,
            "cost": cost,
            "residual_end_pn": float(r.residual_end),
            "inner_converged": bool(r.inner_converged),
            "inner_iters": int(r.inner_iters),
        },
        device=str(r.ledger.get("device", "cuda")),
        parameter_provenance={
            "resting_bound_myosin_fraction": "PI_GAP",
            "resting_bound_myosin_force_pn": "PI_GAP",
            "resting_bound_myosin_capture_um": "PI_GAP",
            "gate_floor_pn": "DERIVED",
        },
        void_ceiling=VOID_CEILING,
        gate_passed=bool(plateau < gate_floor),
        residual=plateau,
        signal=signal,
        notes={
            "purpose": (
                "stage-2 prerequisite measurement. The roadmap asks for the cost of a CONVERGED native "
                "step; this run establishes whether one is reachable from the current solver at all, "
                "and supplies the per-iteration price the answer will be built from either way"
            ),
            "not_a_claim": (
                "no tension, no equilibrium magnitude, and no engine speed. If the regime is STALLED, "
                "the only transferable outputs are the factor above the gate floor — the target the "
                "backbone-aware preconditioner has to close — and the per-iteration cost"
            ),
        },
        declared_commit=args.build_commit,
    )

    write_artifact(out, record)
    print(f"  residual_start = {signal:8.3f} pN   plateau max|PF| = {plateau:8.3f} pN")
    print(f"  gate floor     = {gate_floor:8.3f} pN (derived)   -> factor {cost['factor_above_gate_floor']:.2f}x")
    print(f"  regime         = {fit['regime']}")
    print(f"  wall           = {wall:.0f} s over {total_inner} inner iterations "
          f"({per_iteration_s * 1e3:.3f} ms/iteration)")
    if cost["projected_converged_step_seconds"] is not None:
        print(f"  projected converged step = {cost['projected_converged_step_seconds']:.0f} s")
    else:
        print(f"  projection WITHHELD: {cost['projection_withheld_because']}")
    print(f"  gate verdict   = {record['gate']['verdict']}  (timing comparable={record['timing']['comparable']})")
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())

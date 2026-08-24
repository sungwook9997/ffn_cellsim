#!/usr/bin/env python
r"""Is the resting native inner solve STALLED, or converging slowly?  Nobody has ever looked.

WHY THIS EXISTS.  Every resting record in this repository carries ``residual_start`` and
``residual_end`` and **nothing between them**, so the one question that decides what to do next has
never been askable from an artifact: does the residual keep descending past the 60-iteration budget, or
does it flatten?  Those two answers point at opposite work — a budget/preconditioner problem versus an
operator that cannot close this configuration — and `STATE.md` (b) already carries one negative result
on the second reading (the fiber-quotient coarse operator does not close the static resting residual).

The 2026-08-09 full-native run measured **47.405 pN start -> 47.157 candidate in 60/60 iterations**, a
0.52% reduction, and the per-node force field was **uniform** (max/p50 = 1.50, flat across radial
shells, top 0.4% of nodes carrying 0.80% of |F|^2).  Uniformity rules out a construction cause.  It does
not distinguish stall from slow descent, and no scalar can: only the CURVE can.

WHAT IT DOES.  One outer step, ``capture_convergence_history=True``, at a budget large enough for the
shape to be legible, on the FULL NATIVE population.  The driver already records
``(iteration, max|dx|, constraint, max|PF|)`` on device and hands it back between accepted steps; this
script is the entry point that turns that on and writes the curve beside its record.

WHAT IT DOES NOT DO.  It does not accept a step, tune anything, or produce a physics magnitude.  The
step will almost certainly be REJECTED — that is expected and is not the measurement.  The measurement
is the shape of the descent, and it is a DIAGNOSTIC, not a gate: no threshold is declared here, because
declaring one after seeing a curve is precisely what the charter forbids.  What the curve licenses is
the choice of the next piece of work, and that choice is the PI's.

Read the result as: a straight line on log(residual) vs iteration is slow geometric convergence and the
budget is the lever; a knee that flattens is the operator, and the budget will not save it.

Sanity Gate: no constant is introduced; the tolerance is the driver's own
``sqrt(eps_f64) * convergence_length_um`` and is not touched.  CUDA-only (I0-A) — a non-CUDA device
RAISES rather than producing a number that may not be quoted.  Every host read is out-of-loop, between
the accepted-step boundary and the write.

    MEM_GB=8 CPUS_PER_TASK=8 gpu-submit 4090-1 12:00:00 \
        ~/miniforge3/envs/ffn_sim/bin/python -m aleph.scripts.ac_resting_residual_curve \
            --n-inner 4000 --record-every 10

``MEM_GB`` is a HOST cgroup limit (``--mem``), unrelated to the card's 24 GiB of VRAM.  It was ``48``
here, never measured; a full-native composed step peaks at **733 MB** host RSS (measured 2026-08-15,
``/usr/bin/time -v``, warm Warp kernel cache), so 8 is >10x the observed peak with room for a cold-cache
JIT.  This is not cosmetic: ``gpu-check`` demands ``MEM_GB + 8`` GB FREE and the box caps at 53 GB, so
the old ask required 56 GB — effectively the whole machine — and was rejected on any half-loaded box.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import warp as wp

_ROOT = Path(__file__).resolve().parents[1]
_OUT = _ROOT / "outputs" / "ac" / "resting_residual_curve"


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Resting native inner-solve residual curve (diagnostic).")
    p.add_argument("--n-filaments", type=int, default=70686, help="full native cortex count")
    p.add_argument("--n-inner", type=int, default=4000, help="inner iteration budget for the curve")
    p.add_argument("--record-every", type=int, default=10, help="history stride (driver reshape_every)")
    p.add_argument("--dt-phys", type=float, default=0.05)
    p.add_argument("--inner-solver", type=str, default="explicit")
    p.add_argument("--preload-erm-balance", action=argparse.BooleanOptionalAction, default=False,
                   dest="preload_erm_balance",
                   help="pre-stretch the ERM tethers so the resting turgor is held at t=0. DEFAULT OFF, "
                        "which is what every earlier run of this script used — including the one behind "
                        "STATE.md's C-2 row. `driver._preload_erm_resting_balance`'s own docstring says "
                        "that without it 'the membrane carries the full turgor pressure unbalanced' and "
                        "'there is no force-balanced equilibrium there, which is why no inner solver "
                        "converges from it'. Measured 2026-08-10 on a full native dump: 99.98%% of the "
                        "resting sum|F|^2 sits on the 642 membrane nodes at p50 39.4 pN, against the "
                        "~40 pN per node that docstring predicts.")
    p.add_argument("--cortex-prestrain", type=float, default=0.0, dest="cortex_prestrain",
                   help="cortex crosslink pre-strain fed to _preload_cortex_pretension [dimensionless]")
    p.add_argument("--out", type=Path, default=_OUT / "residual_curve.json")
    p.add_argument("--declared-commit", type=str, default=None,
                   help="build commit, for a run host that is not a git checkout (recorded source=declared)")
    p.add_argument("--closure-verified", type=int, default=None,
                   help="first-party import-closure files hash-matched against this host before launch")
    return p


def main() -> None:
    a = _build_argparser().parse_args()
    wp.init()
    device = wp.get_device()
    if not device.is_cuda:
        raise RuntimeError(
            f"I0-A: the residual curve is a native measurement and needs a CUDA GPU. "
            f"Resolved {str(device)!r} is not CUDA — refusing rather than producing an unquotable number."
        )

    from aleph.components.incumbent.assemble import CellConfig
    from aleph.components.incumbent.driver import run_from_resting
    from aleph.engine.contracts import EvidenceLabel, EvidenceRung, QuantitativeClaim
    from aleph.engine.observe.artifact import (
        observation_artifact, timing_block, write_artifact,
    )

    cfg = CellConfig(
        n_filaments=a.n_filaments, with_myosin=True, with_membrane=True,
        with_nucleus=True, with_pressure=True,
    )
    print(f"[curve] device={device}  N={a.n_filaments}  budget={a.n_inner}  stride={a.record_every}",
          flush=True)

    t0 = time.time()
    report = run_from_resting(
        cfg, n_inner=a.n_inner, dt_phys=a.dt_phys, reshape_every=a.record_every,
        inner_solver=a.inner_solver, capture_convergence_history=True,
        preload_erm_balance=a.preload_erm_balance, cortex_prestrain=a.cortex_prestrain,
    )
    wall = time.time() - t0

    history = getattr(report, "convergence_history", None) or {}
    iters = list(history.get("iteration", []))
    pf = list(history.get("max_projected_force_pn", history.get("projected_force_pn", [])))
    dx = list(history.get("max_displacement_um", []))

    verdict = {
        "inner_converged": bool(report.inner_converged),
        "outer_accepted": bool(report.outer_accepted),
        "inner_iters": int(report.inner_iters),
        "inner_tolerance_um": float(report.inner_tolerance_um),
        "residual_candidate": float(report.residual_candidate),
        "residual_end": float(report.residual_end),
    }
    shape = None
    if len(pf) >= 4:
        # Descriptive only — reported so a reader sees the shape without replotting, never thresholded.
        y = np.log10(np.maximum(np.asarray(pf, dtype=float), 1e-300))
        x = np.asarray(iters, dtype=float)
        half = len(y) // 2
        shape = {
            "log10_slope_first_half": float(np.polyfit(x[:half], y[:half], 1)[0]),
            "log10_slope_second_half": float(np.polyfit(x[half:], y[half:], 1)[0]),
        }

    record = observation_artifact(
        run_label="resting inner-solve residual vs iteration budget (diagnostic)",
        evidence=EvidenceLabel(
            rung=EvidenceRung.CUDA_UNIT,
            quantitative=QuantitativeClaim.BLOCKED,
            basis="one --from-resting outer step at full native under a raised inner budget; the step is "
                  "REJECTED and no magnitude is claimed. What is measured is whether the residual "
                  "descends as the budget grows.",
        ),
        config={
            "n_filaments": int(a.n_filaments), "n_inner": int(a.n_inner),
            "record_every": int(a.record_every), "dt_phys": float(a.dt_phys),
            "inner_solver": a.inner_solver,
            "with_myosin": True, "with_membrane": True, "with_nucleus": True, "with_pressure": True,
        },
        census={
            "n_filaments": int(a.n_filaments), "n_actin": int(report.n_actin),
            "n_total": int(report.n_total),
        },
        # REQUIRED: the state BEFORE the step. Without it the two endpoints below are one number.
        t0={
            "residual_start_pn": float(report.residual_start),
            "inner_tolerance_um": float(report.inner_tolerance_um),
            "r_mean_start_um": float(report.r_mean_start),
        },
        timing=timing_block(
            wall_seconds=wall, n_inner_iterations=int(report.inner_iters),
            device_note="build + one outer step; the step was rejected, so this cost belongs to a "
                        "NON-converged attempt and is not an engine speed.",
        ),
        measurements={
            "residual_start_pn": float(report.residual_start),
            "residual_candidate_pn": float(report.residual_candidate),
            "residual_end_pn": float(report.residual_end),
            "budget_used": int(report.inner_iters),
            "curve": {"iteration": iters, "max_projected_force_pn": pf, "max_displacement_um": dx},
            "curve_samples": len(pf),
            "shape": shape,
            # THE BASELINE THIS WAS MEASURED FROM, recorded because it decides what the curve means.
            # Every earlier run of this script left both OFF — the flags did not exist here — including
            # the one behind STATE.md's C-2 row. With the ERM preload off, the membrane carries the full
            # turgor unbalanced and the driver's own docstring says there is no force-balanced
            # equilibrium to converge to, so a non-descending residual is the expected outcome of the
            # baseline rather than evidence about the operator.
            "preload_erm_balance": bool(a.preload_erm_balance),
            "cortex_prestrain": float(a.cortex_prestrain),
            "erm_preload_ledger": report.ledger.get("erm_preload") if hasattr(report, "ledger") else None,
        },
        device=str(device),
        # The run host is not a git checkout, so the commit is the caller's assertion. This is
        # observation_artifact's OWN parameter, not build_stamp's: re-emitting through the schema
        # made the direct build_stamp call below dead, and the parameter that replaced it is one
        # the caller did not know to fill. The record then reported "no commit was declared",
        # which blames the input for a wire that was never connected.
        declared_commit=a.declared_commit,
        # No gate verdict on purpose. observation_artifact requires all of
        # void_ceiling/gate_passed/residual/signal or NONE, and a diagnostic that declares no threshold
        # has no ceiling to be scored against — so it takes none, and the API enforces exactly the
        # discipline this run's notes claim.
        notes={
            "kind": "diagnostic",
            "no_threshold_declared": "This run scores nothing. Declaring a threshold after seeing the "
                                     "number is what the charter forbids; what the result licenses is "
                                     "the choice of next work, and that choice is the PI's.",
            "curve_not_captured": (
                "capture_convergence_history populates only on the accelerated-solver branch; this ran "
                f"{a.inner_solver!r}, so curve_samples may be 0 and these are ENDPOINTS, not a shape."
            ),
            "verdict_of_the_step": verdict,
            "import_closure_verified_files": a.closure_verified,
            "import_closure_note": (
                "run_provenance.remote_mismatches hash-compared the driver's first-party import closure "
                "against THIS host before launch: 0 mismatched, 0 unverifiable. The build is therefore "
                "established by CONTENT, not by the caller's assertion, even though the host is not a "
                "git checkout."
            ),
        },
    )

    a.out.parent.mkdir(parents=True, exist_ok=True)
    write_artifact(a.out, record)
    print(f"[curve] wrote {a.out}  ({len(pf)} samples, {wall:.0f}s)", flush=True)
    print(f"[curve] step verdict: converged={report.inner_converged} accepted={report.outer_accepted} "
          f"residual {report.residual_start:.4f} -> {report.residual_candidate:.4f}", flush=True)


if __name__ == "__main__":
    main()

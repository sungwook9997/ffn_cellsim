#!/usr/bin/env python
r"""C-1 method 3 — creep: hold the load, watch the flow. The textbook counterpart of relaxation.

WHY A THIRD METHOD.  Methods 1 and 2 are step-relaxation and oscillatory.  The third canonical way a
laboratory characterises a medium is CREEP: hold the load constant and watch the material give way.
In linear viscoelasticity the creep and relaxation experiments do NOT return the same time constant —
``viscoelastic_relaxation.py``'s own Sanity Gate states it as a property to check, "Retardation >
relaxation: tau_eps > tau_sigma for a genuine SLS".  Two experiments, one material, two numbers, and
the textbook says so.  This asks whether the same split appears for a poroelastic medium, which is the
same question C-1 asks and one nobody can answer with a single instrument.

THE MEASUREMENT.  The contact patch is held at fixed pore pressure (a Dirichlet load, re-imposed every
step) while the surrounding medium fills.  What is recorded is the FLUX the patch must supply to hold
that pressure — the rate of change of total pore content — which is what the instrument measuring a
held load actually reads.  As the medium around the patch equilibrates the required flux decays, and
the reported time constant is its ``1/e`` crossing:

    ``tau_creep = 1/e time of dQ/dt``,  where ``Q(t)`` is total pore content.

DECLARED IN ADVANCE: the EXPONENT only, ``tau_creep ~ a^2/D`` by dimensional analysis.  The prefactor
is not declared and must not be — a held-pressure patch in a box has no clean closed form, and the
whole point is to see how far its prefactor sits from the relaxation method's 1.396834.

NOT DECLARED, AND DELIBERATELY: whether ``tau_creep`` is larger or smaller than the relaxation time.
The SLS result says retardation exceeds relaxation, but that is a viscoelastic solid, not a
poroelastic medium, and importing its inequality as a prediction would be assuming the answer to the
question this run exists to ask.

Runtime: NVIDIA Warp on CUDA only (I0-A).

Sanity Gate:
    * dimensions: a, dx [µm]; D [µm²/s]; Q [Pa·µm³]; dQ/dt [Pa·µm³/s]; tau [s].
    * boundary: mass is NOT conserved (the patch is a source), so the step method's drift guard does
      not apply and is not reused; the guard here is the wall margin on the filled region.
    * sign sense: the patch is held ABOVE the surroundings, so the required flux is positive and
      decays; a negative or rising flux raises through the scorer.
    * CFL/precision: ``CFL_SAFETY * dx^2/(6D)``, float64.
    * measurement protocol: the projection (dQ/dt, its 1/e crossing) is fixed before the run; the
      held-patch mask is built from the rest configuration and never follows the field.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import warp as wp

from aleph.laws.biot_fluid_warp import BiotField
from aleph.scripts.c1_psr_scoring import one_over_e_time, power_law_exponent

CFL_SAFETY: float = 0.2
N_SAMPLES: int = 120

#: ── AMENDED 2026-08-10, after the first creep run FAILED its own grid check by a factor of 21. ──
#:
#: The first observable was the ``1/e`` crossing of ``dQ/dt``, and it is grid-dependent for a reason
#: that is physics, not numerics: **the flux out of a held patch does not decay to zero.** For a sphere
#: held at fixed potential in an infinite diffusive medium the flux tends to a finite steady value
#: ``4 pi D a``, approached as
#:
#:     ``flux(t) = 4 pi D a (1 + a / sqrt(pi D t))``
#:
#: so the transient excess falls as ``t^(-1/2)`` — no exponential, and therefore no intrinsic ``1/e``
#: time. What the old observable actually measured was "one part in e of the INITIAL flux", and the
#: initial flux is the gradient at the patch edge, which is ``~1/dx``. Refining the mesh raised the
#: reference and moved the crossing: 0.03405 / 0.01830 / 0.01444 s at dx = 0.30 / 0.20 / 0.15.
#:
#: The replacement has a closed form and no reference to the initial value at all: the time at which
#: the TRANSIENT EXCESS FALLS TO THE STEADY FLUX, ``a / sqrt(pi D t) = 1``, giving
#:
#:     ``tau_creep = a^2 / (pi D)``,  i.e.  C = 1/pi = 0.318310  in a^2/D units.
#:
#: The steady flux is MEASURED from the late plateau rather than taken from ``4 pi D a``, because the
#: patch is a discretised sphere and its effective radius is not exactly ``a`` — using the analytic
#: value would import a geometry assumption the mesh does not honour.
#:
#: Both observables are recorded. The old one is what demonstrates that a creep protocol can hand you
#: the mesh, which is a result about the method and is not deleted because a better one exists.
CREEP_ORACLE_C: float = 1.0 / math.pi

#: Fraction of the early curve dropped before the fit: those samples carry the discretised gradient at
#: the patch edge, which is the mesh, and including them tilts the line.
FIT_SKIP_FRACTION: float = 0.05

#: Relative spread of tau across the dx sweep above which this method is NOT measuring the medium.
#: Part of the top-level verdict — see the comment beside it for why that needed saying.
GRID_INVARIANCE_MAX: float = 0.05
DECLARED_EXPONENT_A: float = 2.0
DECLARED_EXPONENT_D: float = -1.0
EXPONENT_TOLERANCE: float = 0.08
R2_FLOOR: float = 0.98

#: Run length in units of the relaxation-method oracle time, so the two methods are integrated over a
#: comparable physical window rather than each over its own convenient one.
RUN_LENGTH_IN_ORACLE_TAU: float = 6.0

#: The filled region must stay clear of the frozen faces. Same reason as method 1, different measure:
#: here the front advances as sqrt(2Dt) from the patch edge.
WALL_MARGIN_MIN: float = 4.0


@wp.kernel
def _hold_patch_kernel(
    p: wp.array3d(dtype=wp.float64),
    patch: wp.array3d(dtype=wp.int32),
    level: wp.float64,
) -> None:
    """Re-impose the held pore pressure on the loaded patch (a Dirichlet load, every step)."""
    i, j, k = wp.tid()
    if patch[i, j, k] == 1:
        p[i, j, k] = level


@wp.kernel
def _total_kernel(p: wp.array3d(dtype=wp.float64), out: wp.array(dtype=wp.float64)) -> None:
    """Total pore content, reduced on device."""
    i, j, k = wp.tid()
    wp.atomic_add(out, 0, p[i, j, k])


def measure_tau_creep(*, a_um: float, d_um2_s: float, n: int, dx: float, device: str,
                      oracle_tau_s: float) -> dict:
    """Hold the patch, record the flux needed to hold it, and return the flux decay's 1/e time."""
    c = 0.5 * (n - 1) * dx
    ax = (np.arange(n, dtype=np.float64) * dx) - c
    r2 = (ax[:, None, None] ** 2) + (ax[None, :, None] ** 2) + (ax[None, None, :] ** 2)
    patch_np = (r2 <= a_um * a_um).astype(np.int32)
    patch = wp.array(np.ascontiguousarray(patch_np), dtype=wp.int32, device=device)

    field = BiotField(n=n, dx=dx, D=d_um2_s, device=device)
    field.set_field(patch_np.astype(np.float64))          # start with only the patch loaded
    t_end = RUN_LENGTH_IN_ORACLE_TAU * oracle_tau_s
    dt = min(CFL_SAFETY * field.dt_max, t_end / N_SAMPLES)
    steps_per_sample = max(1, int(round((t_end / N_SAMPLES) / dt)))
    dt = (t_end / N_SAMPLES) / steps_per_sample

    acc = wp.zeros(1, dtype=wp.float64, device=device)

    def total() -> float:
        acc.zero_()
        wp.launch(_total_kernel, dim=(n, n, n), inputs=[field._p, acc], device=device)
        return float(acc.numpy()[0]) * dx ** 3

    times, totals = [0.0], [total()]
    for s in range(N_SAMPLES):
        for _ in range(steps_per_sample):
            field.step(dt)
            wp.launch(_hold_patch_kernel, dim=(n, n, n),
                      inputs=[field._p, patch, wp.float64(1.0)], device=device)
        times.append((s + 1) * steps_per_sample * dt)
        totals.append(total())

    t = np.array(times)
    q = np.array(totals)
    # The instrument reads the flux it must supply: dQ/dt, on the sample clock (central where possible).
    flux = np.gradient(q, t)
    tm = 0.5 * (t[:-1] + t[1:])
    flux_mid = np.diff(q) / np.diff(t)
    monotone = bool(np.all(np.diff(flux_mid) <= 1e-12))
    tau = one_over_e_time(tm, flux_mid) if monotone else None
    if tau is None:
        # DISCRETISATION NOISE FALLBACK, and it CHANGES THE DATA, so it is recorded rather than done
        # quietly: the flux is replaced by its non-increasing envelope before the crossing is read.
        # A run that needed this is a run whose sampling was too fine for its own step size, and the
        # record says so, because a hull silently applied is a smoothing nobody voted for.
        keep = flux_mid > 0
        if keep.sum() >= 3:
            fm = np.maximum.accumulate(flux_mid[keep][::-1])[::-1]
            tau = one_over_e_time(tm[keep], fm)

    # The grid-independent observable, by MODEL FIT rather than by a reference level.
    #
    # There is no plateau to reference. The excess decays as t^(-1/2), so at 40 tau it is still 16 % of
    # the steady value and a late-window median overestimates the steady flux by ~25 % — the estimate
    # depends on how long you ran, which is the same disease as depending on how fine your mesh is.
    #
    # So fit the KNOWN FORM instead: flux = A + B * t^(-1/2), a straight line in t^(-1/2). The steady
    # flux is the intercept, the transient is the slope, and the crossover — where B/sqrt(t) = A — is
    # tau = (B/A)^2. This is a model-based inversion and its assumption is explicit, which is precisely
    # what the method inventory says a creep analysis always does.
    #
    # The early samples are dropped: the first few carry the discretised gradient at the patch edge,
    # which is the mesh, and including them tilts the line.
    tau_excess = None
    flux_steady = float("nan")
    skip = max(2, int(FIT_SKIP_FRACTION * flux_mid.size))
    inv_sqrt_t = 1.0 / np.sqrt(tm[skip:])
    if flux_mid.size - skip >= 8 and np.all(flux_mid[skip:] > 0):
        b_slope, a_int = np.polyfit(inv_sqrt_t, flux_mid[skip:], 1)
        flux_steady = float(a_int)
        if a_int > 0.0 and b_slope > 0.0:
            tau_excess = float((b_slope / a_int) ** 2)

    front = a_um + math.sqrt(2.0 * d_um2_s * t_end)
    return {
        "a_um": float(a_um), "D_um2_s": float(d_um2_s), "n": int(n), "dx_um": float(dx),
        "L_um": float((n - 1) * dx),
        "tau_creep_s": None if tau_excess is None else float(tau_excess),
        "tau_creep_oracle_s": float(CREEP_ORACLE_C * a_um * a_um / d_um2_s),
        "flux_steady": flux_steady,
        "tau_one_over_e_of_initial_s": None if tau is None else float(tau),
        "retired_observable_note": ("tau_one_over_e_of_initial_s is the FIRST observable, kept because "
                                    "its grid dependence is the result that motivated the replacement; "
                                    "it references the initial flux, which is ~1/dx"),
        "flux_monotone": monotone,
        "monotone_hull_applied": bool(not monotone and tau is not None),
        "oracle_tau_relaxation_s": float(oracle_tau_s),
        "front_um": float(front),
        "wall_margin": float(((n - 1) * dx / 2.0) / front),
        "dt_s": float(dt), "cfl_ratio": float(dt / field.dt_max),
        "curve": {"t_s": [float(v) for v in tm], "flux": [float(v) for v in flux_mid]},
        "content": {"t_s": [float(v) for v in t], "Q": [float(v) for v in q],
                    "dQdt_endpoints": [float(flux[0]), float(flux[-1])]},
    }


def main() -> None:
    from aleph.scripts.c1_psr_scoring import gaussian_spreading_tau

    ap = argparse.ArgumentParser(description="C-1 method 3: creep under a held pore-pressure load.")
    ap.add_argument("--n", type=int, default=321)
    ap.add_argument("--dx", type=float, default=0.15)
    ap.add_argument("--D", type=float, default=50.0)
    ap.add_argument("--a-list", type=float, nargs="+", default=[0.5, 0.7, 1.0, 1.4], dest="a_list")
    ap.add_argument("--a-boxes", type=str, nargs="+", default=None, dest="a_boxes",
                    help="`a:L:dx` triples [µm]. Overrides --a-list. CREEP NEEDS THIS MOST of the three "
                         "methods: at a single dx its prefactor drifts 0.2766 -> 0.3257 across the "
                         "sweep, tracking a/dx from 3.33 to 12.67, which biases the exponent to ~2.12. "
                         "Holding a/dx fixed removes the confound the same way it did for the step "
                         "method, where the residual then proved to be a pure function of a/dx.")
    ap.add_argument("--d-list", type=float, nargs="+", default=[12.5, 25.0, 50.0, 100.0], dest="d_list")
    ap.add_argument("--a-fixed", type=float, default=1.0, dest="a_fixed")
    ap.add_argument("--dx-list", type=float, nargs="+", default=[0.3, 0.2, 0.15], dest="dx_list",
                    help="grid spacings for the invariance check. CREEP NEEDS THIS MORE THAN THE OTHER "
                         "METHODS: the patch is loaded instantaneously into an empty medium, so the "
                         "initial flux is set by the discretisation, and a 1/e crossing read from a "
                         "curve whose start is grid-dependent may be grid-dependent too.")
    ap.add_argument("--out", type=Path, default=Path("aleph/outputs/ac/c1_biot_creep"))
    ap.add_argument("--build-commit", type=str, default=None, dest="build_commit")
    args = ap.parse_args()

    wp.init()
    dev = str(wp.get_device())
    if not wp.get_device().is_cuda:
        raise RuntimeError(f"C-1 creep needs a CUDA GPU (I0-A). Resolved {dev!r} is not CUDA.")
    if not args.build_commit:
        raise RuntimeError("--build-commit is REQUIRED (the run host is an rsync'd tree with no git).")
    print(f"[c1-creep] device={dev} n={args.n} dx={args.dx} D={args.D}", flush=True)

    t0 = time.time()
    a_pts, d_pts = [], []
    if args.a_boxes:
        plan = []
        for spec in args.a_boxes:
            parts = spec.split(":")
            a_v, l_v = float(parts[0]), float(parts[1])
            dx_v = float(parts[2]) if len(parts) > 2 else args.dx
            plan.append((a_v, int(round(l_v / dx_v)) + 1, dx_v))
    else:
        plan = [(a, args.n, args.dx) for a in args.a_list]

    for a, n_a, dx_a in plan:
        pt = measure_tau_creep(a_um=a, d_um2_s=args.D, n=n_a, dx=dx_a, device=dev,
                               oracle_tau_s=gaussian_spreading_tau(a, args.D))
        a_pts.append(pt)
        c_meas = (pt["tau_creep_s"] * args.D / (a * a)) if pt["tau_creep_s"] else float("nan")
        print(f"[c1-creep] a={a:5.2f} L={pt['L_um']:5.1f} dx={dx_a:6.4f} a/dx={a / dx_a:5.2f}  "
              f"tau_creep={pt['tau_creep_s']}  C={c_meas:.5f} (oracle {CREEP_ORACLE_C:.5f})  "
              f"wall_margin={pt['wall_margin']:.2f}", flush=True)
    for d in args.d_list:
        pt = measure_tau_creep(a_um=args.a_fixed, d_um2_s=d, n=args.n, dx=args.dx, device=dev,
                               oracle_tau_s=gaussian_spreading_tau(args.a_fixed, d))
        d_pts.append(pt)
        print(f"[c1-creep] D={d:6.1f}  tau_creep={pt['tau_creep_s']}", flush=True)

    def score(points: list[dict], key: str, declared: float) -> dict:
        ok = [p for p in points
              if p["tau_creep_s"] is not None and p["wall_margin"] >= WALL_MARGIN_MIN]
        dropped = [{key: p[key], "tau_creep_s": p["tau_creep_s"], "wall_margin": p["wall_margin"],
                    "why": ("flux never crossed 1/e" if p["tau_creep_s"] is None
                            else f"wall margin {p['wall_margin']:.2f} < {WALL_MARGIN_MIN}")}
                   for p in points if p not in ok]
        out = {"sweep": key, "declared_exponent": declared, "n_kept": len(ok), "rejected": dropped}
        if len(ok) < 3:
            out["verdict"] = "INSUFFICIENT"
            return out
        p_exp, cc, r2 = power_law_exponent([q[key] for q in ok], [q["tau_creep_s"] for q in ok])
        out.update({"exponent": p_exp, "prefactor": cc, "r2": r2,
                    "verdict": "PASS" if (r2 >= R2_FLOOR and abs(p_exp - declared) <= EXPONENT_TOLERANCE)
                               else "FAIL"})
        return out

    dx_pts = []
    L_fixed = (args.n - 1) * args.dx
    for dxv in args.dx_list:
        n_dx = int(round(L_fixed / dxv)) + 1
        pt = measure_tau_creep(a_um=args.a_fixed, d_um2_s=args.D, n=n_dx, dx=dxv, device=dev,
                               oracle_tau_s=gaussian_spreading_tau(args.a_fixed, args.D))
        dx_pts.append(pt)
        print(f"[c1-creep] dx={dxv:5.3f} (n={n_dx})  tau_creep={pt['tau_creep_s']}  "
              f"hull={pt['monotone_hull_applied']}", flush=True)
    dx_taus = [p["tau_creep_s"] for p in dx_pts if p["tau_creep_s"] is not None]
    dx_spread = (float((max(dx_taus) - min(dx_taus)) / np.median(dx_taus))
                 if len(dx_taus) >= 2 else None)
    grid_invariance_verdict = ("INSUFFICIENT" if dx_spread is None
                               else ("PASS" if dx_spread <= GRID_INVARIANCE_MAX else "FAIL"))

    scored = {"a_um": score(a_pts, "a_um", DECLARED_EXPONENT_A),
              "D_um2_s": score(d_pts, "D_um2_s", DECLARED_EXPONENT_D)}
    ratios = [p["tau_creep_s"] / p["oracle_tau_relaxation_s"] for p in a_pts if p["tau_creep_s"]]
    record = {
        "record": "run-record@2", "kind": "gate", "gate": "c1_biot_creep_timescale",
        "device": str(dev), "build": {"commit": args.build_commit, "source": "declared"},
        "declared_before_run": {
            "method": "held pore-pressure patch (Dirichlet load); the instrument reads the flux dQ/dt",
            "observable": "the time at which the transient flux excess falls to the MEASURED steady "
                          "flux — grid-independent because it references no initial value",
            "exponents": {"a_um": DECLARED_EXPONENT_A, "D_um2_s": DECLARED_EXPONENT_D},
            "prefactor": f"C = 1/pi = {CREEP_ORACLE_C:.6f} in a^2/D units, from "
                         "flux(t) = 4 pi D a (1 + a/sqrt(pi D t)) for a sphere held at fixed potential",
            "direction_vs_relaxation": "NOT declared — importing the SLS inequality tau_eps > tau_sigma "
                                       "would assume the answer this run exists to measure",
            "r2_floor": R2_FLOOR, "wall_margin_min": WALL_MARGIN_MIN,
        },
        "scored": scored,
        "creep_over_relaxation": {
            "per_point": ratios,
            "median": float(np.median(ratios)) if ratios else None,
            "why": "the textbook split between creep and relaxation, measured for a poroelastic medium",
        },
        "grid_invariance": {
            "relative_spread": dx_spread,
            "bound": GRID_INVARIANCE_MAX,
            "verdict": grid_invariance_verdict,
            "why": "the patch is loaded instantaneously, so the initial flux is grid-set; if tau moves "
                   "with dx then this method's number is partly the mesh, not the medium",
            "points": [{"dx_um": p["dx_um"], "n": p["n"], "tau_creep_s": p["tau_creep_s"],
                        "monotone_hull_applied": p["monotone_hull_applied"]} for p in dx_pts],
        },
        "sweeps": {"a_um": a_pts, "D_um2_s": d_pts, "dx_um": dx_pts},
        # THE GRID CHECK IS PART OF THE VERDICT, and it was not until run 59 exposed the gap: that run
        # reported `grid_invariance: FAIL` at 7.8 % spread against a 5 % bound while the top-level
        # verdict said PASS, because the verdict was assembled from `scored` alone. A gate that
        # computes a check and then does not consult it is worse than one that never ran it — the
        # check's presence in the record reads as evidence that it was applied.
        "verdict": ("PASS" if (all(s.get("verdict") == "PASS" for s in scored.values())
                               and grid_invariance_verdict == "PASS") else "FAIL"),
        "timing": {"wall_s": time.time() - t0, "comparable": True},
        "may_not_be_quoted_for": ["any statement about a cell", "a cytoplasm viscosity", "any rung"],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "record.json").write_text(json.dumps(record, indent=2))
    print(f"\n[c1-creep] VERDICT {record['verdict']} — a:{scored['a_um'].get('exponent')} "
          f"prefactor={scored['a_um'].get('prefactor')} "
          f"creep/relax median={record['creep_over_relaxation']['median']} -> {args.out}", flush=True)


if __name__ == "__main__":
    main()

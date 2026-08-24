#!/usr/bin/env python
r"""C-1 step 1 — does the poroelastic clock carry a LENGTH SCALE, and the one the theory says?

WHY THIS RUNS FIRST.  The whole C-1 discriminator rests on one asymmetry: of the four relaxation clocks
in this engine, exactly one scales with the probe size (``tau_p ~ L^2/D``,
``aleph/laws/network_warp.py:708``).  That is a property of the **Biot pore-pressure field alone** — it
needs no mechanical convergence, no contact law, and no cell — so it is measurable today, while C-2 is
open.  If it does not hold, the C-1 premise is wrong and the right time to find out is now.

WHAT IS MEASURED, DECLARED BEFORE THE RUN.  A Gaussian pore-pressure blob ``p(r,0) = P0 exp(-r^2/2a^2)``
is placed in the device field and relaxed by :class:`aleph.laws.biot_fluid_warp.BiotField`.  The recorded
projection is the normalised L2 amplitude of the deviation from the (conserved) box mean::

    A(t) = || p(t) - mean(p) ||_2  /  || p(0) - mean(p) ||_2

and the scored time is its ``1/e`` crossing.  In free space this has a closed form: the blob stays
Gaussian with ``sigma^2 = a^2 + 2Dt``, so ``A(t) = (a^2/(a^2+2Dt))^(3/4)`` and

    tau = a^2 (e^(4/3) - 1) / (2 D)  ~=  1.39683 a^2 / D

so this gate declares an EXPONENT **and a PREFACTOR** in advance, not a shape to be interpreted after.

Three sweeps, three pre-declared exponents (``aleph/scripts/c1_psr_scoring.power_law_exponent``):

    tau vs a   ->  p = +2     (the length scale: the whole point)
    tau vs D   ->  p = -1     (it is a diffusivity, not a rate constant)
    tau vs dx  ->  p =  0     (grid invariance; a constant that moves with dx is not a constant)

THE TRAP THIS GATE IS BUILT AROUND.  A no-flux box has its own diffusive time, and once ``tau``
approaches it the BOX sets the answer and the ``a`` dependence flattens to zero — the measurement would
then report "no length scale" for a reason that has nothing to do with the physics under test.  Every
point therefore carries ``tau / tau_box`` and is REJECTED above ``BOX_CONTAMINATION_MAX``, declared here
and not adjustable from the command line.  Rejected points are recorded, never dropped.

Runtime: NVIDIA Warp on CUDA only (I0-A).  Authored on the dev Mac; measured on the GPU host.

Sanity Gate:
    * dimensions: a, dx [µm]; D [µm²/s]; t, tau [s]; p [pN/µm² == Pa]; A dimensionless.
    * boundary: the no-flux box conserves total pore mass; the run asserts it to 1e-10 relative, which
      is the control that the solver is diffusing rather than leaking.
    * CFL/precision: every step uses ``CFL_SAFETY * dx^2/(6D)``; float64 throughout.
    * sign sense: A(t) is non-increasing (enforced by the scorer, which raises on a rise).
    * measurement protocol: the projection above is fixed BEFORE the run and is the only scalar read;
      the field is stored so nothing is re-chosen afterwards.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import warp as wp

from aleph.laws.biot_fluid_warp import BiotField
from aleph.scripts.c1_psr_scoring import gaussian_spreading_tau, one_over_e_time, power_law_exponent

#: Fraction of the CFL bound each step takes. Below 1 by construction; the field module raises above it.
CFL_SAFETY: float = 0.2

#: A point whose measured tau exceeds this fraction of the box's own diffusive time is BOX-CONTAMINATED:
#: the finite domain, not the blob, is setting the answer. Declared here, before the run, and not a CLI
#: argument — a validity bound the run could widen after seeing the data is not a validity bound.
BOX_CONTAMINATION_MAX: float = 0.05

#: The SHARPER validity bound, and the one that actually binds: relative drift of total pore mass.
#:
#: FOUND BY THE SETUP RENDER, BEFORE ANY SCORING RUN (2026-08-10).  ``biot_diffusion_kernel`` implements
#: its boundary as ``p_new = p`` on the six faces.  That is a FROZEN face, not a mirrored ghost cell, and
#: it behaves like no-flux only while the field beside the face is unperturbed — the module's docstring
#: says the box "conserves total pore mass", and it does, right up until the blob arrives.  The
#: verification pass measured the two failing together: drift 0.0 / 6.4e-9 / 1.3e-2 at a = 0.75 / 1.5 /
#: 3.0 µm, against tau/oracle 1.0021 / 0.9974 / 0.9746.  The mass invariant moves FIRST and is measured
#: by the run itself, so it is the better guard; ``BOX_CONTAMINATION_MAX`` alone accepted the a=3 point.
#:
#: Tightening a bound is not re-thresholding to pass: it rejects more, the verification pass scores
#: nothing by construction, and this is exactly what that pass exists to find.
MASS_DRIFT_MAX: float = 1.0e-6

#: How far past the expected tau to integrate.  Lowered 4.0 -> 2.0 after the verification pass: the
#: blob's width grows as ``sqrt(a^2 + 2Dt)``, so integrating to 4 tau spreads it to 3.5a and walks it
#: into the frozen face, which is what produced the mass drift above.  At 2 tau it reaches 2.6a, which
#: keeps the wall unperturbed for every point the mass guard admits.  The 1/e crossing is at tau, so
#: nothing the projection reads is lost — only tail the scorer never looks at.
RUN_LENGTH_IN_TAU: float = 2.0

#: Samples recorded per run — enough that the 1/e crossing is interpolated over a small interval.
N_SAMPLES: int = 120

#: Pre-declared fit-quality floor for every log-log exponent below (C-1 section 6.3).
R2_FLOOR: float = 0.99

#: How far a measured exponent may sit from its declared value and still PASS.
EXPONENT_TOLERANCE: float = 0.05

#: Grid invariance is scored on SPREAD, not on a slope — see :func:`_score_sweep`. A tau that moves by
#: more than this across a 2.7x range of dx is not a grid-invariant constant.
INVARIANCE_SPREAD_MAX: float = 0.02

#: The observation operators, and the closed-form exponent q of each in the blob width.
#:
#: This is the multi-method comparison, run on ONE simulated state so that nothing but the observer
#: differs. Every real instrument reads a functional of the field, not the field, and the functional it
#: reads is a design decision its inventor made — an AFM cantilever integrates a traction over its
#: contact, a tracer bead reports a local displacement, a whole-field imaging method sees something
#: closer to a norm. They are all correct and they do not agree.
#:
#: ``q = None`` marks an observer with no closed form. That is not a defect: the patch observer is the
#: one that mimics an instrument with a finite contact, and its lack of a clean inversion is exactly
#: why a real AFM needs a model to turn F(t) into a time constant.
OBSERVERS: dict[str, float | None] = {
    "l2_norm": 0.75,        # whole-field: ||p - mean p||_2 ~ sigma^-3/2
    "peak": 1.5,            # single point at the blob centre: p(0,t) ~ sigma^-3
    "patch_mean": None,     # mean over a finite contact patch — the instrument-like observer
}

#: Radius of the ``patch_mean`` observer's contact, as a multiple of the blob width ``a``. Fixed to the
#: blob so the observer is self-similar under the ``a`` sweep; otherwise the sweep would change the
#: instrument and the probe at the same time and could not attribute the result to either.
PATCH_RADIUS_IN_A: float = 1.0

#: Which observer the top-level ``tau_measured_s`` and the box-contamination guard refer to. Named
#: rather than defaulted so a reader of the record knows which of the three the headline number is.
PRIMARY_OBSERVER: str = "l2_norm"


@wp.kernel
def _sum_kernel(p: wp.array3d(dtype=wp.float64), out: wp.array(dtype=wp.float64)) -> None:
    """Total pore pressure, reduced on device. ``out[0] += p``."""
    i, j, k = wp.tid()
    wp.atomic_add(out, 0, p[i, j, k])


@wp.kernel
def _observe_kernel(
    p: wp.array3d(dtype=wp.float64),
    patch: wp.array3d(dtype=wp.int32),
    mean: wp.array(dtype=wp.float64),
    ctr: wp.int32,
    out: wp.array(dtype=wp.float64),
) -> None:
    """Reduce the observation operators on device: ``out = [sum d^2, sum d over patch, n_patch, centre]``.

    Deviations are taken against the mean supplied in ``mean[0]``, which the caller reduced in the pass
    before.  Two passes rather than one because the mean is not known until the field has been summed,
    and estimating it would put a tolerance into a measurement that does not need one.

    Why on device at all: the real sweep runs a 513^3 grid, so a host copy per sample is 1.08 GB and 120
    samples per point would move 130 GB across PCIe to compute three scalars.  The point observer still
    reads one cell on the host between samples, which is a scalar, not the field.
    """
    i, j, k = wp.tid()
    d = p[i, j, k] - mean[0]
    wp.atomic_add(out, 0, d * d)
    if patch[i, j, k] == 1:
        wp.atomic_add(out, 1, d)
        wp.atomic_add(out, 2, wp.float64(1.0))
    # One thread writes the point observer's cell, so the single-point reading never costs a field copy.
    if i == ctr and j == ctr and k == ctr:
        out[3] = d


def _gaussian_blob(n: int, dx: float, a_um: float, p0: float = 1.0) -> np.ndarray:
    """A Gaussian pore-pressure blob centred in an ``n^3`` grid of spacing ``dx`` [Pa]."""
    c = 0.5 * (n - 1) * dx
    ax = (np.arange(n, dtype=np.float64) * dx) - c
    r2 = (ax[:, None, None] ** 2) + (ax[None, :, None] ** 2) + (ax[None, None, :] ** 2)
    return p0 * np.exp(-r2 / (2.0 * a_um * a_um))


def _box_tau(n: int, dx: float, d_um2_s: float) -> float:
    """The no-flux box's own slowest diffusive time [s] — the contamination yardstick.

    The slowest non-uniform Neumann mode of a cube of side ``L`` decays as ``exp(-D (pi/L)^2 t)``, so
    the box time is ``(L/pi)^2 / D``. Compared against, never fitted to.
    """
    L = (n - 1) * dx
    return (L / np.pi) ** 2 / d_um2_s


def measure_tau(
    *,
    a_um: float,
    d_um2_s: float,
    n: int,
    dx: float,
    device: str,
    dump_slices: int = 0,
) -> dict:
    """Relax one Gaussian blob and return its measured ``1/e`` time with its validity evidence.

    ``dump_slices`` stores that many evenly spaced mid-plane ``z`` slices for the renderer.  The render
    is not decoration: it is how a human confirms the blob is centred, round, inside the box and
    actually spreading — an experiment wired to the wrong array produces a perfectly clean exponent
    from the wrong field, and no scalar in this record would say so.
    """
    field = BiotField(n=n, dx=dx, D=d_um2_s, device=device)
    p0 = _gaussian_blob(n, dx, a_um)
    field.set_field(p0)
    mass0 = field.pore_mass()

    tau_expected = gaussian_spreading_tau(a_um, d_um2_s)
    t_end = RUN_LENGTH_IN_TAU * tau_expected
    dt = min(CFL_SAFETY * field.dt_max, t_end / N_SAMPLES)
    steps_per_sample = max(1, int(round((t_end / N_SAMPLES) / dt)))
    dt = (t_end / N_SAMPLES) / steps_per_sample

    # The patch mask is built once from the rest configuration: the observer is an instrument, and an
    # instrument whose footprint follows the field it is measuring is not an instrument.
    c = 0.5 * (n - 1) * dx
    ax = (np.arange(n, dtype=np.float64) * dx) - c
    r2_grid = (ax[:, None, None] ** 2) + (ax[None, :, None] ** 2) + (ax[None, None, :] ** 2)
    patch_np = (r2_grid <= (PATCH_RADIUS_IN_A * a_um) ** 2).astype(np.int32)
    n_patch_cells = int(patch_np.sum())
    patch_d = wp.array(np.ascontiguousarray(patch_np), dtype=wp.int32, device=device)
    ctr = n // 2
    cells = float(n) ** 3
    acc = wp.zeros(4, dtype=wp.float64, device=device)
    sum_d = wp.zeros(1, dtype=wp.float64, device=device)
    mean_d = wp.zeros(1, dtype=wp.float64, device=device)

    def read() -> dict[str, float]:
        """Every observation operator applied to ONE state, reduced on device (see :func:`_observe_kernel`)."""
        sum_d.zero_()
        wp.launch(_sum_kernel, dim=(n, n, n), inputs=[field._p, sum_d], device=device)
        total = float(sum_d.numpy()[0])
        mean_d.assign(np.array([total / cells], np.float64))
        acc.zero_()
        wp.launch(_observe_kernel, dim=(n, n, n),
                  inputs=[field._p, patch_d, mean_d, wp.int32(ctr), acc], device=device)
        s2, patch_sum, patch_n, centre = (float(v) for v in acc.numpy())
        return {
            "l2_norm": float(np.sqrt(s2)),
            "peak": centre,
            "patch_mean": (patch_sum / patch_n) if patch_n else 0.0,
            "_total": total,
        }

    slice_at = (set(np.linspace(0, N_SAMPLES, dump_slices).astype(int).tolist())
                if dump_slices > 0 else set())
    slices: list[dict] = []
    if 0 in slice_at:
        slices.append({"t_s": 0.0, "z_midplane": p0[:, :, ctr].tolist()})

    first = read()
    mass0 = first.pop("_total") * dx ** 3
    times = [0.0]
    series: dict[str, list[float]] = {k: [v] for k, v in first.items()}
    total_last = mass0
    for s in range(N_SAMPLES):
        for _ in range(steps_per_sample):
            field.step(dt)
        times.append((s + 1) * steps_per_sample * dt)
        r = read()
        total_last = r.pop("_total") * dx ** 3
        for k, v in r.items():
            series[k].append(v)
        if (s + 1) in slice_at:
            slices.append({"t_s": times[-1], "z_midplane": field.get_field()[:, :, ctr].tolist()})

    mass1 = total_last
    tau_box = _box_tau(n, dx, d_um2_s)
    t_arr = np.array(times)
    observers: dict[str, dict] = {}
    for name, q in OBSERVERS.items():
        tau = one_over_e_time(t_arr, np.array(series[name]))
        observers[name] = {
            "q": q,
            "tau_measured_s": None if tau is None else float(tau),
            "tau_oracle_s": None if q is None else float(gaussian_spreading_tau(a_um, d_um2_s, q)),
            "box_contamination": None if tau is None else float(tau / tau_box),
        }
    primary = observers[PRIMARY_OBSERVER]
    return {
        "a_um": float(a_um),
        "D_um2_s": float(d_um2_s),
        "n": int(n),
        "dx_um": float(dx),
        "L_um": float((n - 1) * dx),
        "a_over_L": float(a_um / ((n - 1) * dx)),
        "patch_radius_um": float(PATCH_RADIUS_IN_A * a_um),
        "patch_cells": int(n_patch_cells),
        "dt_s": float(dt),
        "cfl_ratio": float(dt / field.dt_max),
        "tau_measured_s": primary["tau_measured_s"],
        "tau_oracle_s": float(tau_expected),
        "tau_box_s": float(tau_box),
        "box_contamination": primary["box_contamination"],
        "observers": observers,
        "pore_mass_rel_drift": float(abs(mass1 - mass0) / abs(mass0)) if mass0 else 0.0,
        "curve": {"t_s": [float(v) for v in times],
                  **{k: [float(v) for v in s] for k, s in series.items()}},
        "slices": slices,
    }


def _score_sweep(points: list[dict], key: str, declared: float) -> dict:
    """Score one sweep, after removing box-contaminated points (recorded, never dropped).

    A DECLARED EXPONENT OF ZERO IS NOT SCORED BY r2, and this is not a convenience.  Grid invariance
    means ``tau`` does not move, so ``log tau`` is constant, its total variance is ~0, and ``r2`` — a
    ratio to that variance — is degenerate: it reports a number that says nothing about invariance and
    can be arbitrarily negative for a perfectly flat sweep.  The right statistic for "does not move" is
    the SPREAD, so an invariance sweep is scored on relative spread against a bound declared here.
    Applying the fit-quality floor to it would have been a gate that cannot pass for the right reason.
    """
    kept, rejected = [], []
    for p in points:
        row = {k: p[k] for k in
               (key, "tau_measured_s", "box_contamination", "a_over_L", "pore_mass_rel_drift")}
        if p["tau_measured_s"] is None:
            row["why"] = "never reached 1/e within the integrated window"
            rejected.append(row)
        elif p["pore_mass_rel_drift"] > MASS_DRIFT_MAX:
            row["why"] = (f"pore mass drifted by {p['pore_mass_rel_drift']:.3e} > {MASS_DRIFT_MAX} — the "
                          "frozen face has started participating, so this is not free-space spreading")
            rejected.append(row)
        elif p["box_contamination"] > BOX_CONTAMINATION_MAX:
            row["why"] = f"box-contaminated: tau/tau_box > {BOX_CONTAMINATION_MAX}"
            rejected.append(row)
        else:
            kept.append(p)

    out = {
        "sweep": key,
        "declared_exponent": float(declared),
        "n_points": len(points),
        "n_kept": len(kept),
        "rejected": rejected,
    }
    if len(kept) < 3:
        out["verdict"] = "INSUFFICIENT"
        out["why"] = ("fewer than three uncontaminated points; a power law was not fitted rather than "
                      "fitted through whatever survived")
        return out

    taus = np.array([q["tau_measured_s"] for q in kept], dtype=float)
    if declared == 0.0:
        spread = float((taus.max() - taus.min()) / np.median(taus))
        out.update({
            "scored_as": "relative_spread",
            "relative_spread": spread,
            "spread_max": INVARIANCE_SPREAD_MAX,
            "verdict": "PASS" if spread <= INVARIANCE_SPREAD_MAX else "FAIL",
        })
        return out

    p, c, r2 = power_law_exponent([q[key] for q in kept], taus)
    out.update({
        "scored_as": "log_log_slope",
        "exponent": p, "prefactor": c, "r2": r2,
        "r2_floor": R2_FLOOR, "exponent_tolerance": EXPONENT_TOLERANCE,
        "verdict": "PASS" if (r2 >= R2_FLOOR and abs(p - declared) <= EXPONENT_TOLERANCE) else "FAIL",
    })
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="C-1 step 1: the Biot field's poroelastic timescale.")
    ap.add_argument("--n", type=int, default=257, help="grid points per side")
    ap.add_argument("--dump-slices", type=int, default=8, dest="dump_slices",
                    help="mid-plane z slices stored per `a` point for the renderer (0 = none)")
    ap.add_argument("--verify-only", action="store_true", dest="verify_only",
                    help="run only the `a` sweep and stop — the short setup-verification pass whose "
                         "render is checked in a browser BEFORE the full sweep is launched")
    ap.add_argument("--dx", type=float, default=0.25, help="grid spacing [µm]")
    ap.add_argument("--D", type=float, default=50.0, help="poroelastic diffusivity [µm²/s] (Moeendarbary)")
    ap.add_argument("--a-list", type=float, nargs="+", default=[0.5, 0.75, 1.0, 1.5, 2.0, 3.0],
                    dest="a_list", help="Gaussian blob widths to sweep [µm] (uses --n for every point)")
    ap.add_argument("--a-boxes", type=str, nargs="+", default=None, dest="a_boxes",
                    help="`a:L` or `a:L:dx` [µm], e.g. 0.4:24:0.075 0.75:24:0.075 0.75:48:0.15. "
                         "Overrides --a-list; dx defaults to --dx. Several boxes extend the lever arm "
                         "past what one can hold without the frozen face participating, and REPEATING "
                         "ONE `a` ACROSS TWO BOXES is the control that changing the domain did not move "
                         "tau. Giving each box its own dx lets a/dx be held CONSTANT across the whole "
                         "sweep, so an a-dependence cannot be a discretisation-ratio dependence in "
                         "disguise — which a single dx over a decade of `a` could not rule out.")
    ap.add_argument("--d-list", type=float, nargs="+", default=[12.5, 25.0, 50.0, 100.0],
                    dest="d_list", help="diffusivities to sweep [µm²/s]")
    ap.add_argument("--dx-list", type=float, nargs="+", default=[0.5, 0.375, 0.25, 0.1875],
                    dest="dx_list", help="grid spacings for the invariance sweep [µm]")
    ap.add_argument("--a-fixed", type=float, default=1.5, dest="a_fixed",
                    help="blob width held fixed in the D and dx sweeps [µm]")
    ap.add_argument("--out", type=Path, default=Path("aleph/outputs/ac/c1_biot_timescale"))
    ap.add_argument("--build-commit", type=str, default=None, dest="build_commit",
                    help="REQUIRED: the commit whose code is on this host, verified with "
                         "run_provenance.remote_mismatches before launch (the host has no git).")
    args = ap.parse_args()

    wp.init()
    dev = str(wp.get_device())
    if not wp.get_device().is_cuda:
        raise RuntimeError(f"C-1 Biot timescale needs a CUDA GPU (I0-A). Resolved {dev!r} is not CUDA.")
    if not args.build_commit:
        raise RuntimeError("--build-commit is REQUIRED: the run host is an rsync'd tree with no git, so "
                           "nothing here can read the build and an unstamped record cannot be traced.")
    print(f"[c1-biot] device={dev} n={args.n} dx={args.dx} D={args.D}", flush=True)

    t0 = time.time()
    sweeps: dict[str, list[dict]] = {"a_um": [], "D_um2_s": [], "dx_um": []}

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
        pt = measure_tau(a_um=a, d_um2_s=args.D, n=n_a, dx=dx_a, device=dev,
                         dump_slices=args.dump_slices)
        sweeps["a_um"].append(pt)
        obs = pt["observers"]
        print(f"[c1-biot] a={a:6.3f} L={pt['L_um']:5.1f} dx={dx_a:6.4f} a/dx={a / dx_a:5.2f} n={n_a:4d}  "
              f"tau_l2={pt['tau_measured_s']}  "
              f"oracle={pt['tau_oracle_s']:.6g}  tau_peak={obs['peak']['tau_measured_s']}  "
              f"tau_patch={obs['patch_mean']['tau_measured_s']}  "
              f"drift={pt['pore_mass_rel_drift']:.2e}  box={pt['box_contamination']:.4f}", flush=True)

    if args.verify_only:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "verify_record.json").write_text(json.dumps(
            {"record": "run-record@2", "kind": "diagnostic", "verify_only": True, "device": str(dev),
             "build": {"commit": args.build_commit, "source": "declared"},
             "sweeps": {"a_um": sweeps["a_um"]},
             "why": "setup-verification pass: the render of these slices is checked in a real browser "
                    "before the full sweep is launched. Scores nothing."},
            indent=2))
        print(f"[c1-biot] VERIFY-ONLY written -> {args.out}/verify_record.json", flush=True)
        return

    for d in args.d_list:
        pt = measure_tau(a_um=args.a_fixed, d_um2_s=d, n=args.n, dx=args.dx, device=dev)
        sweeps["D_um2_s"].append(pt)
        print(f"[c1-biot] D={d:7.2f}  tau={pt['tau_measured_s']}  oracle={pt['tau_oracle_s']:.6g}",
              flush=True)

    # Grid invariance holds the PHYSICAL box fixed while dx changes, so n moves with it — otherwise the
    # sweep changes the domain and the box time, and would measure that instead of the discretisation.
    L_fixed = (args.n - 1) * args.dx
    for dx in args.dx_list:
        n_dx = int(round(L_fixed / dx)) + 1
        pt = measure_tau(a_um=args.a_fixed, d_um2_s=args.D, n=n_dx, dx=dx, device=dev)
        sweeps["dx_um"].append(pt)
        print(f"[c1-biot] dx={dx:6.4f} µm (n={n_dx})  tau={pt['tau_measured_s']}", flush=True)

    scored = {
        "a_um": _score_sweep(sweeps["a_um"], "a_um", 2.0),
        "D_um2_s": _score_sweep(sweeps["D_um2_s"], "D_um2_s", -1.0),
        "dx_um": _score_sweep(sweeps["dx_um"], "dx_um", 0.0),
    }

    # THE MULTI-METHOD COMPARISON. Every observer is scored on the SAME `a` sweep, so the exponent and
    # the prefactor are separated: if the exponents agree and the prefactors do not, then the scaling is
    # a property of the physics and the magnitude is a property of the instrument — which is the whole
    # reason C-1 scores an exponent. Reported whatever the answer is; a disagreement in the EXPONENT
    # would be the more interesting result and would refute the design, so it is not filtered out.
    per_observer = {}
    for name in OBSERVERS:
        pts = [{**p, "tau_measured_s": p["observers"][name]["tau_measured_s"],
                "box_contamination": p["observers"][name]["box_contamination"]}
               for p in sweeps["a_um"]]
        per_observer[name] = _score_sweep(pts, "a_um", 2.0)
    # The box-overlap control: any `a` measured in two different boxes must give the same tau. If it
    # does not, the box is in the answer and the exponent is not a property of the field.
    by_a: dict[float, list[dict]] = {}
    for p in sweeps["a_um"]:
        by_a.setdefault(round(p["a_um"], 9), []).append(p)
    overlap = [{"a_um": a, "L_um": [q["L_um"] for q in qs],
                "tau_s": [q["tau_measured_s"] for q in qs],
                "rel_spread": (max(t := [q["tau_measured_s"] for q in qs]) - min(t)) / np.median(t)}
               for a, qs in by_a.items() if len(qs) > 1 and all(q["tau_measured_s"] for q in qs)]

    exps = [s["exponent"] for s in per_observer.values() if "exponent" in s]
    pres = [s["prefactor"] for s in per_observer.values() if "prefactor" in s]
    method_comparison = {
        "per_observer": per_observer,
        "exponent_spread": float(max(exps) - min(exps)) if len(exps) > 1 else None,
        "prefactor_ratio_max_over_min": float(max(pres) / min(pres)) if len(pres) > 1 else None,
        "closed_form_ratio_l2_over_peak": float(
            gaussian_spreading_tau(1.0, 1.0, OBSERVERS["l2_norm"])
            / gaussian_spreading_tau(1.0, 1.0, OBSERVERS["peak"])),
        "reading": ("if the exponents agree while the prefactors do not, the SCALING belongs to the "
                    "physics and the MAGNITUDE belongs to the instrument"),
    }
    # The prefactor check is the part that makes this more than a shape: compare the measured tau to the
    # closed form point by point, on the uncontaminated points of the `a` sweep.
    ok = [p for p in sweeps["a_um"]
          if p["tau_measured_s"] is not None and p["box_contamination"] <= BOX_CONTAMINATION_MAX]
    ratios = [p["tau_measured_s"] / p["tau_oracle_s"] for p in ok]
    verdicts = [s["verdict"] for s in scored.values()]

    record = {
        "record": "run-record@2",
        "kind": "gate",
        "gate": "c1_biot_poroelastic_timescale",
        "device": str(dev),
        "build": {"commit": args.build_commit, "source": "declared",
                  "why": "the run host is an rsync'd tree, not a git checkout"},
        "declared_before_run": {
            "projection": "normalised L2 norm of (p - mean p); scored time is its 1/e crossing",
            "oracle": "tau = a^2 (e^(4/3) - 1) / (2D), free-space Gaussian spreading",
            "exponents": {"a_um": 2.0, "D_um2_s": -1.0, "dx_um": 0.0},
            "r2_floor": R2_FLOOR,
            "box_contamination_max": BOX_CONTAMINATION_MAX,
        },
        "scored": scored,
        "box_overlap_control": {
            "points": overlap,
            "max_rel_spread": max((o["rel_spread"] for o in overlap), default=None),
            "bound": INVARIANCE_SPREAD_MAX,
            "verdict": ("PASS" if overlap and max(o["rel_spread"] for o in overlap)
                        <= INVARIANCE_SPREAD_MAX else ("ABSENT" if not overlap else "FAIL")),
            "why": "one `a` measured in two boxes must give one tau, or the box is in the answer",
        },
        "method_comparison": method_comparison,
        "oracle_ratio": {
            "n": len(ratios),
            "median": float(np.median(ratios)) if ratios else None,
            "min": float(np.min(ratios)) if ratios else None,
            "max": float(np.max(ratios)) if ratios else None,
        },
        "max_pore_mass_rel_drift": max(
            (p["pore_mass_rel_drift"] for v in sweeps.values() for p in v), default=0.0),
        "sweeps": sweeps,
        "verdict": "PASS" if all(v == "PASS" for v in verdicts) else "FAIL",
        "timing": {"wall_s": time.time() - t0, "comparable": True},
        "may_not_be_quoted_for": [
            "any statement about the CELL — this is the field alone, with no mechanics and no contact",
            "the C-1 gate itself, which is scored on tau(a) of an indented cell",
            "any rung",
        ],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "record.json").write_text(json.dumps(record, indent=2))
    _plot(sweeps, scored, args.out / "c1_biot_timescale.png")
    print(f"\n[c1-biot] VERDICT {record['verdict']} — "
          f"a:{scored['a_um'].get('exponent')} D:{scored['D_um2_s'].get('exponent')} "
          f"dx:{scored['dx_um'].get('exponent')}  oracle ratio median "
          f"{record['oracle_ratio']['median']}  -> {args.out}/record.json", flush=True)


def _plot(sweeps: dict, scored: dict, path: Path) -> None:
    """Decay curves and the three exponent fits, with the oracle overlaid. No axis truncation."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 4, figsize=(19, 4.2))
    ax = axes[0]
    styles = {"l2_norm": "-", "peak": "--", "patch_mean": ":"}
    mid = sweeps["a_um"][len(sweeps["a_um"]) // 2]
    for name, ls in styles.items():
        c = mid["curve"]
        y = np.array(c[name])
        ax.plot(np.array(c["t_s"]) * 1e3, y / y[0], ls, lw=1.4, label=name)
    ax.axhline(1.0 / np.e, ls="-.", lw=0.8, color="0.4")
    for name in styles:
        tau = mid["observers"][name]["tau_measured_s"]
        if tau is not None:
            ax.axvline(tau * 1e3, lw=0.6, color="0.6")
    ax.set_yscale("log")
    ax.set_xlabel("t [ms]")
    ax.set_ylabel("observer reading (normalised) [-]")
    ax.set_title(f"THREE observers, ONE field (a={mid['a_um']:.2f} µm)\n"
                 "same physics, different 1/e times", fontsize=9)
    ax.legend(fontsize=7)

    for ax, (key, xlabel) in zip(axes[1:], [("a_um", "blob width a [µm]"),
                                            ("D_um2_s", r"D [µm$^2$/s]"),
                                            ("dx_um", "grid spacing dx [µm]")]):
        pts = [p for p in sweeps[key] if p["tau_measured_s"] is not None]
        x = np.array([p[key] for p in pts])
        y = np.array([p["tau_measured_s"] for p in pts])
        keep = np.array([p["box_contamination"] <= BOX_CONTAMINATION_MAX for p in pts])
        ax.loglog(x[keep], y[keep], "o", ms=6, label="measured")
        if (~keep).any():
            ax.loglog(x[~keep], y[~keep], "x", ms=7, color="crimson", label="box-contaminated (excluded)")
        ax.loglog(x, [p["tau_oracle_s"] for p in pts], "-", lw=1.0, color="0.5", label="oracle")
        s = scored[key]
        if "exponent" in s:
            ax.set_title(f"p = {s['exponent']:.3f} (declared {s['declared_exponent']:+.0f}), "
                         f"$r^2$={s['r2']:.4f}", fontsize=9)
        else:
            ax.set_title(f"{s['verdict']}: {s.get('why', '')}", fontsize=8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$\tau_{1/e}$ [s]")
        ax.legend(fontsize=7)

    fig.suptitle("C-1 step 1 — the poroelastic clock carries a length scale (Biot field alone, no cell)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()

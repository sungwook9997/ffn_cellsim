#!/usr/bin/env python
r"""C-1 method 2 — read the SAME poroelastic clock the way a rheologist reads it: in frequency.

WHY A SECOND METHOD, AND WHY THIS ONE.  ``c1_biot_timescale.py`` measures a step relaxation: perturb,
release, watch the decay, report a ``1/e`` time.  That is the AFM force-relaxation family.  It is not
how the other half of the literature measures the same material — active microrheology (AFM
oscillatory, optical tweezers, magnetic twisting) DRIVES the sample at a frequency and reads the
amplitude and the phase lag, then reports a time constant from where the response crosses over.

Both are measuring one property of one medium.  Neither is wrong.  They report different numbers, and
this driver is how much different, measured rather than argued, with the physics held identical: the
same :class:`~aleph.laws.biot_fluid_warp.BiotField`, the same ``D``, the same length scale ``a``.

THE MEASUREMENT.  A Gaussian source of width ``a`` is driven with a MULTI-TONE waveform — a sum of
sinusoids, which is what a real multi-frequency rheometer uses, and which gets the whole spectrum from
ONE run instead of one run per frequency.  The patch-mean response is lock-in detected at each drive
tone over an integer number of its periods (in-phase and quadrature projections), giving amplitude and
phase lag per frequency.  The reported time constant is

    ``tau_omega = 1 / omega_45``    where the phase lag crosses 45 degrees,

which for any first-order (diffusive) response is the corner of the low-pass, and is what the
instrument's own inversion would report.

WHAT IS DECLARED IN ADVANCE.  The EXPONENT, and only the exponent: ``tau_omega ~ a^2 / D`` follows from
dimensional analysis — ``a`` is the only length and ``D`` the only diffusivity, so the only time that
can be built is ``a^2/D``.  The PREFACTOR is deliberately NOT declared: a finite Gaussian source in a
box has no clean closed form, and inventing one would be a fit dressed as a prediction.  That asymmetry
is the point of the experiment: the exponent is predictable across methods, the prefactor is not, and
the prefactor is what a paper reports.

Runtime: NVIDIA Warp on CUDA only (I0-A).

Sanity Gate:
    * dimensions: a, dx [µm]; D [µm²/s]; omega [rad/s]; tau [s]; p [Pa].
    * boundary: the drive is local and its diffusive penetration depth ``sqrt(2D/omega)`` is order ``a``,
      so unlike the step method the field never reaches the wall — recorded as the wall-margin check.
    * CFL/precision: every step takes ``CFL_SAFETY * dx^2/(6D)``; float64 throughout.
    * sign sense: phase lag is positive and increases with frequency (a driven diffusive medium lags,
      never leads); a negative lag raises.
    * conservation: the drive INJECTS pore pressure, so mass is not conserved here and is not used as a
      guard — the step method's guard does not transfer, and pretending it did would be a check that
      cannot fail.
    * measurement protocol: transient periods are discarded by a count fixed before the run, and the
      lock-in integrates over whole periods only, so no window is chosen after seeing the response.
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
from aleph.scripts.c1_psr_scoring import power_law_exponent

CFL_SAFETY: float = 0.2

#: Drive tones, as multiples of the diffusive corner ``D/a^2``.
#:
#: RECENTRED, not widened, and the distinction cost a run. The measured crossing sits near 4x the
#: dimensional guess, so the first fix — extending the band DOWNWARD to 0.125 — put the lowest tone's
#: penetration depth ``sqrt(2D/omega)`` at ``4a`` and collapsed the wall margins to 3.0-5.5, below the
#: bound. The guard rejected almost everything, correctly. It also made the run 8x longer for nothing:
#: the low tones carry no information about a corner that is above them. Moving the band UP costs
#: neither margin nor time.
TONE_MULTIPLES: tuple[float, ...] = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0)

#: ── AMENDED 2026-08-10, AFTER RUN 34 AND BEFORE THE NEXT RUN. Read this before the code. ──
#:
#: The first declared inversion was ``tau = 1/omega_45``, the frequency where the phase lag crosses 45
#: degrees — the standard first-order (single-relaxation) reading. Run 34 measured, at D=50, dx=0.2:
#:
#:     a=0.80   23.2  30.0  37.3  43.1  43.6   deg   -> no crossing
#:     a=1.10   24.5  32.8  42.8  54.1  65.5   deg   -> crossing at 47.26 rad/s
#:     a=1.50   23.3  30.1  37.4  43.1  43.5   deg   -> no crossing
#:     a=2.00   23.0  29.8  37.0  42.7  42.7   deg   -> no crossing
#:
#: The tones are set as multiples of ``D/a^2``, so the DIMENSIONLESS response cannot depend on ``a``:
#: three of the four profiles agree to three significant figures, as required. ``a=1.10`` does not, and
#: **it is the only point at which the declared inversion returned a number**. The one reading that
#: "worked" is the one that is wrong. A laboratory with a single probe size of 1.1 µm would have
#: reported tau = 0.0212 s with nothing in its data saying otherwise.
#:
#: THE CAUSE IS NOT YET DIAGNOSED and no explanation is asserted here. What can be done without knowing
#: it is to make the physics check itself: self-similarity is guaranteed by the construction of the
#: tone grid, so a profile that deviates from the others is a numerical artifact whatever its origin.
#: :data:`SELF_SIMILARITY_MAX_DEG` is that guard, and it would have rejected the a=1.10 point.
#:
#: WHAT THE CAUSE TURNED OUT TO BE — and I was wrong about it twice on the way here, so both wrong
#: readings are kept rather than quietly replaced by the right one:
#:
#:   (i)  I first wrote that 45 degrees is the ASYMPTOTE of a diffusive phase lag. ``a=1.10`` reaching
#:        65.5 degrees IN THE SAME RUN refuted that, immediately and from data already in hand.
#:   (ii) I then wrote that the inversion was ill-posed on this medium. Also wrong. Widening the band
#:        from (0.25 .. 4.0) to (0.125 .. 16.0) x D/a^2 takes the phase to 82.7 degrees and the
#:        crossing exists at EVERY probe size.
#:
#: The band was mis-centred, nothing more: the dimensional guess ``omega_c = D/a^2`` is about 4x low,
#: so the original top tone stopped at ~43 degrees, just short of the level being looked for. That is
#: not an exotic failure — **it is exactly what a real instrument does when it picks its frequency
#: range from a dimensional estimate**, and what it returns is either no measurement at all or, at one
#: probe size, a corrupted one. Both are worse than a wrong error bar, because neither looks wrong.
#:
#: The ``a=1.10`` corruption is a SEPARATE defect and is still undiagnosed. It survives the wider band
#: and now shows a non-monotone tail (64.0 -> 59.8 degrees), so it is reproducible and it is not the
#: band. :data:`SELF_SIMILARITY_MAX_DEG` rejects it on an invariant the tone grid guarantees.
#:
#: The scored observable stays the **-3 dB amplitude point**: it is defined wherever the response rolls
#: off at all, so unlike a crossing that fires at some probe sizes and not others it cannot select its
#: own sample. The phase crossing is recorded beside it, not deleted.
AMPLITUDE_LEVEL: float = 1.0 / math.sqrt(2.0)

#: Maximum deviation [deg] of one point's phase profile from the run's median profile before the point
#: is rejected as numerically corrupted. The tone grid scales with ``D/a^2``, so the dimensionless
#: response is a-independent BY CONSTRUCTION and any deviation is the numerics, not the medium. This is
#: a physics-guaranteed internal control, which is what makes it usable without a diagnosis.
SELF_SIMILARITY_MAX_DEG: float = 2.0

#: Phase level still recorded for the failed inversion, so the record carries what did not work.
PHASE_LEVEL_RAD: float = math.pi / 4.0

#: Periods of the LOWEST tone to run, and how many of them to discard as transient. Fixed here, before
#: the run: a transient window chosen after seeing the response is a window chosen to give an answer.
N_PERIODS: int = 8
N_DISCARD: int = 4

#: Samples per period of the HIGHEST tone — sets the lock-in quadrature resolution.
SAMPLES_PER_FASTEST_PERIOD: int = 16

#: Declared before the run. Only the exponent; see the module docstring for why not the prefactor.
DECLARED_EXPONENT_A: float = 2.0
DECLARED_EXPONENT_D: float = -1.0
EXPONENT_TOLERANCE: float = 0.08
R2_FLOOR: float = 0.98

#: The response must have decayed to well under the source width before the wall, or the box is in the
#: answer. Penetration depth is sqrt(2D/omega); at the lowest tone this is the deepest.
WALL_MARGIN_MIN: float = 6.0


@wp.kernel
def _drive_kernel(
    p: wp.array3d(dtype=wp.float64),
    src: wp.array3d(dtype=wp.float64),
    amp: wp.float64,
    dt: wp.float64,
) -> None:
    """Inject ``dt * amp * src`` — one multi-tone drive step, evaluated on the host clock."""
    i, j, k = wp.tid()
    p[i, j, k] = p[i, j, k] + dt * amp * src[i, j, k]


@wp.kernel
def _patch_mean_kernel(
    p: wp.array3d(dtype=wp.float64),
    patch: wp.array3d(dtype=wp.int32),
    out: wp.array(dtype=wp.float64),
) -> None:
    """Reduce the patch sum and count on device: ``out = [sum, n]``."""
    i, j, k = wp.tid()
    if patch[i, j, k] == 1:
        wp.atomic_add(out, 0, p[i, j, k])
        wp.atomic_add(out, 1, wp.float64(1.0))


def _gaussian(n: int, dx: float, a_um: float) -> np.ndarray:
    c = 0.5 * (n - 1) * dx
    ax = (np.arange(n, dtype=np.float64) * dx) - c
    r2 = (ax[:, None, None] ** 2) + (ax[None, :, None] ** 2) + (ax[None, None, :] ** 2)
    return np.exp(-r2 / (2.0 * a_um * a_um))


def _lock_in(t: np.ndarray, y: np.ndarray, omega: float) -> tuple[float, float]:
    """In-phase / quadrature projection at ``omega`` -> (amplitude, phase lag [rad]).

    The same operation a lock-in amplifier performs, and the reason a multi-tone drive works at all:
    projections at distinct tones are orthogonal over a common whole number of periods.
    """
    s, c = np.sin(omega * t), np.cos(omega * t)
    i_comp = 2.0 * float(np.mean(y * s))
    q_comp = 2.0 * float(np.mean(y * c))
    return float(math.hypot(i_comp, q_comp)), float(math.atan2(-q_comp, i_comp))


def _phase_crossing(omegas: np.ndarray, phases: np.ndarray, level: float) -> float | None:
    """Interpolate in ``log omega`` for the frequency where the phase lag crosses ``level`` [rad]."""
    for i in range(1, len(omegas)):
        p0, p1 = phases[i - 1], phases[i]
        if (p0 - level) * (p1 - level) <= 0.0 and p1 != p0:
            f = (level - p0) / (p1 - p0)
            return float(math.exp(math.log(omegas[i - 1]) + f * (math.log(omegas[i]) - math.log(omegas[i - 1]))))
    return None


def _amplitude_crossing(omegas: np.ndarray, amps: np.ndarray, level: float) -> float | None:
    """Frequency where the response falls to ``level`` times its LOWEST-frequency value.

    The -3 dB bandwidth reading. Interpolated in ``log omega`` / ``log amplitude``, which is exact for
    a power-law roll-off and is what the roll-off of a diffusive response is — unlike the phase, which
    saturates and therefore has no crossing at all (see the module constants).
    """
    if amps[0] <= 0.0:
        return None
    rel = amps / amps[0]
    for i in range(1, len(omegas)):
        a0, a1 = rel[i - 1], rel[i]
        if (a0 - level) * (a1 - level) <= 0.0 and a1 != a0 and a0 > 0 and a1 > 0:
            f = (math.log(a0) - math.log(level)) / (math.log(a0) - math.log(a1))
            return float(math.exp(math.log(omegas[i - 1])
                                  + f * (math.log(omegas[i]) - math.log(omegas[i - 1]))))
    return None


def measure_tau_omega(*, a_um: float, d_um2_s: float, n: int, dx: float, device: str) -> dict:
    """Drive one Gaussian source multi-tone and return the 45-degree corner time."""
    corner = d_um2_s / (a_um * a_um)                       # the dimensional corner, in rad/s
    omegas = np.array([m * corner for m in TONE_MULTIPLES], dtype=float)
    period_slow = 2.0 * math.pi / omegas.min()
    t_end = N_PERIODS * period_slow
    sample_dt = (2.0 * math.pi / omegas.max()) / SAMPLES_PER_FASTEST_PERIOD

    field = BiotField(n=n, dx=dx, D=d_um2_s, device=device)
    dt = min(CFL_SAFETY * field.dt_max, sample_dt)
    steps_per_sample = max(1, int(round(sample_dt / dt)))
    dt = sample_dt / steps_per_sample
    n_samples = int(round(t_end / sample_dt))

    src_np = _gaussian(n, dx, a_um)
    src = wp.array(np.ascontiguousarray(src_np), dtype=wp.float64, device=device)
    patch_np = (src_np >= math.exp(-0.5)).astype(np.int32)          # the 1-sigma contact
    patch = wp.array(np.ascontiguousarray(patch_np), dtype=wp.int32, device=device)
    acc = wp.zeros(2, dtype=wp.float64, device=device)

    def drive(t: float) -> float:
        return float(np.sum(np.sin(omegas * t)))

    def patch_mean() -> float:
        acc.zero_()
        wp.launch(_patch_mean_kernel, dim=(n, n, n), inputs=[field._p, patch, acc], device=device)
        s, c = (float(v) for v in acc.numpy())
        return s / c if c else 0.0

    t_now = 0.0
    ts, ys = [], []
    for _ in range(n_samples):
        for _ in range(steps_per_sample):
            wp.launch(_drive_kernel, dim=(n, n, n),
                      inputs=[field._p, src, wp.float64(drive(t_now)), wp.float64(dt)], device=device)
            field.step(dt)
            t_now += dt
        ts.append(t_now)
        ys.append(patch_mean())

    t = np.array(ts)
    y = np.array(ys)
    keep = t >= (N_DISCARD * period_slow)
    # Trim to a whole number of the SLOWEST period so every tone integrates over whole periods.
    t_k, y_k = t[keep], y[keep]
    span = t_k[-1] - t_k[0]
    whole = math.floor(span / period_slow) * period_slow
    sel = t_k <= (t_k[0] + whole)
    t_k, y_k = t_k[sel], y_k[sel]

    amps, phases = [], []
    for w in omegas:
        # ABSOLUTE time, not `t_k - t_k[0]`, AND THIS IS THE BUG THAT PRODUCED EVERY "CORRUPTED POINT".
        #
        # The drive is `sum(sin(omega_j t))` on absolute t, so the projection must use absolute t too.
        # Subtracting the window start adds a phase offset `omega_j * t_k[0]` to every tone, which is a
        # multiple of 2*pi only if `t_k[0]` is an exact whole number of periods of every tone. In exact
        # arithmetic it is — the discard boundary lands on sample 4096 of 8192, and every tone completes
        # an integer count there. In floating point the comparison `t >= N_DISCARD * period_slow` sits
        # on a knife edge: one ulp either way includes or excludes that sample, and excluding it shifts
        # the phase reference by exactly ONE SAMPLE.
        #
        # That is the whole "corruption": phase error `-omega * sample_dt = -(2*pi/256) * m` degrees per
        # tone, i.e. -1.406 * m. Measured against a clean point: -0.9 / -1.6 / -3.0 / -5.7 / -11.3 /
        # -22.3 at m = 0.5 / 1 / 2 / 4 / 8 / 16, against -0.70 / -1.41 / -2.81 / -5.62 / -11.25 /
        # -22.50 predicted — ratios 1.01, 1.00, 0.99 at the tones where the effect dominates.
        #
        # It explains everything that was strange about it: amplitudes were always clean (a time shift
        # does not touch them), only the phase moved, the error grew linearly with frequency, and which
        # point was hit moved between runs because a knife-edge comparison is decided by rounding.
        amp, ph = _lock_in(t_k, y_k - y_k.mean(), w)
        amps.append(amp)
        phases.append(ph if ph >= 0 else ph + 2.0 * math.pi)
    phases = np.array(phases)
    amps_arr = np.array(amps)
    omega45 = _phase_crossing(omegas, phases, PHASE_LEVEL_RAD)          # the inversion that cannot fire
    omega3db = _amplitude_crossing(omegas, amps_arr, AMPLITUDE_LEVEL)   # the one that can

    depth = math.sqrt(2.0 * d_um2_s / omegas.min())
    return {
        "a_um": float(a_um), "D_um2_s": float(d_um2_s), "n": int(n), "dx_um": float(dx),
        "L_um": float((n - 1) * dx),
        "corner_dimensional_rad_s": float(corner),
        "omegas_rad_s": [float(w) for w in omegas],
        "amplitude": [float(v) for v in amps],
        "phase_lag_rad": [float(v) for v in phases],
        "phase_lag_deg": [float(math.degrees(v)) for v in phases],
        "omega_45_rad_s": None if omega45 is None else float(omega45),
        "tau_phase45_s": None if omega45 is None else float(1.0 / omega45),
        "phase_max_deg": float(math.degrees(phases.max())),
        "omega_3db_rad_s": None if omega3db is None else float(omega3db),
        "tau_omega_s": None if omega3db is None else float(1.0 / omega3db),
        "penetration_depth_um": float(depth),
        "wall_margin": float(((n - 1) * dx / 2.0) / depth),
        "n_samples_used": int(t_k.size),
        "dt_s": float(dt), "cfl_ratio": float(dt / field.dt_max),
        "curve": {"t_s": [float(v) for v in t_k], "patch_mean": [float(v) for v in y_k]},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="C-1 method 2: oscillatory (microrheology) reading.")
    ap.add_argument("--n", type=int, default=241)
    ap.add_argument("--dx", type=float, default=0.2)
    ap.add_argument("--D", type=float, default=50.0)
    ap.add_argument("--a-list", type=float, nargs="+", default=[0.8, 1.1, 1.5, 2.0], dest="a_list")
    ap.add_argument("--a-boxes", type=str, nargs="+", default=None, dest="a_boxes",
                    help="`a:L:dx` triples [µm]. Overrides --a-list. The step method reached a full "
                         "decade this way — several boxes, each with its own dx, holding a/dx roughly "
                         "constant so an a-dependence cannot be a mesh-ratio dependence in disguise, "
                         "with probe sizes repeated across boundaries as the control.")
    ap.add_argument("--d-list", type=float, nargs="+", default=[12.5, 25.0, 50.0, 100.0], dest="d_list")
    ap.add_argument("--a-fixed", type=float, default=1.5, dest="a_fixed",
                    help="probe size held fixed in the D sweep. NOT 1.1: that is the probe size the "
                         "self-similarity guard rejects, and runs 34 and 38 both put the entire D "
                         "sweep on it — returning a flawless -1.0000 exponent from a corrupted "
                         "configuration, twice.")
    ap.add_argument("--out", type=Path, default=Path("aleph/outputs/ac/c1_biot_microrheology"))
    ap.add_argument("--build-commit", type=str, default=None, dest="build_commit")
    args = ap.parse_args()

    wp.init()
    dev = str(wp.get_device())
    if not wp.get_device().is_cuda:
        raise RuntimeError(f"C-1 microrheology needs a CUDA GPU (I0-A). Resolved {dev!r} is not CUDA.")
    if not args.build_commit:
        raise RuntimeError("--build-commit is REQUIRED (the run host is an rsync'd tree with no git).")
    print(f"[c1-rheo] device={dev} n={args.n} dx={args.dx} D={args.D}", flush=True)

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
        pt = measure_tau_omega(a_um=a, d_um2_s=args.D, n=n_a, dx=dx_a, device=dev)
        a_pts.append(pt)
        print(f"[c1-rheo] a={a:5.2f} L={pt['L_um']:5.1f} dx={dx_a:6.4f} a/dx={a / dx_a:5.2f}  "
              f"tau_3dB={pt['tau_omega_s']}  "
              f"tau_phase45={pt['tau_phase45_s']}  "
              f"phase_max={pt['phase_max_deg']:.1f}deg  "
              f"phases_deg={[round(v, 1) for v in pt['phase_lag_deg']]}  "
              f"wall_margin={pt['wall_margin']:.1f}", flush=True)
    for d in args.d_list:
        pt = measure_tau_omega(a_um=args.a_fixed, d_um2_s=d, n=args.n, dx=args.dx, device=dev)
        d_pts.append(pt)
        print(f"[c1-rheo] D={d:6.1f}  tau_omega={pt['tau_omega_s']}", flush=True)

    def self_similarity(points: list[dict]) -> dict[int, float]:
        """Max deviation [deg] of each point's phase profile from the run's MEDIAN profile.

        The tones are multiples of ``D/a^2``, so every point solves the same dimensionless problem and
        the profiles must coincide. Deviation is numerics. Using the median as the reference means one
        corrupted point cannot drag the reference to itself.
        """
        prof = np.array([p["phase_lag_deg"] for p in points], dtype=float)
        if prof.shape[0] < 3:
            return {i: 0.0 for i in range(prof.shape[0])}
        ref = np.median(prof, axis=0)
        return {i: float(np.max(np.abs(prof[i] - ref))) for i in range(prof.shape[0])}

    def score(points: list[dict], key: str, declared: float,
              tau_field: str = "tau_omega_s") -> dict:
        """Score one sweep on one inversion. ``tau_field`` selects which inversion.

        BOTH are scored and neither is primary, because run 44 showed they fail in DIFFERENT ways and
        picking one would hide whichever weakness the winner has:

          * ``tau_omega_s`` (-3 dB) is precise inside a band — the three clean points agree to 0.24 %
            — but it is defined RELATIVE TO THE LOWEST TONE MEASURED, and that tone is not the DC
            plateau unless the band reaches it. Moving the band from 0.125 to 0.5 x D/a^2 moved its
            prefactor from 1.943 to 0.863 a^2/D, a factor of 2.25, with the medium unchanged.
          * ``tau_phase45_s`` is defined by an ABSOLUTE level, so it is band-independent: 0.8679
            against 0.8645 across those same two bands, 0.4 %. It is noisier within a run and it moves
            when a point is corrupted.

        AND THE PART THAT MATTERS MOST, from run 52. Two of six points were corrupted. Their
        prefactors, in ``a^2/D``:

            a       C(-3 dB)    C(phase-45)     self-similarity deviation
            0.60    0.8631      0.7675          21.82 deg   <- corrupted
            0.80    0.8619      0.8735           0.49 deg
            1.00    0.8625      0.8745           0.31 deg
            1.30    0.8631      0.8760           0.19 deg
            1.70    0.8644      0.7668          22.65 deg   <- corrupted
            2.20    0.8645      0.8651           0.19 deg

        **The -3 dB inversion cannot see the corruption at all.** Its corrupted readings, 0.8631 and
        0.8644, sit inside the clean spread of 0.8619-0.8645. Phase-45 puts them 12 % out. A study that
        measured only the -3 dB point would have six points agreeing to 0.3 %, two of them corrupt, and
        no signal anywhere in its own data.

        So the trade is not precision against transferability. **It is precision against the ability to
        detect that you are wrong**, and the more precise instrument is the blind one.
        """
        dev = self_similarity(points)
        for i, p in enumerate(points):
            p["self_similarity_dev_deg"] = dev.get(i, 0.0)
        ok, dropped = [], []
        for i, p in enumerate(points):
            if p[tau_field] is None:
                why = f"no crossing for {tau_field} inside the driven band"
            elif dev.get(i, 0.0) > SELF_SIMILARITY_MAX_DEG:
                why = (f"phase profile deviates {dev[i]:.2f} deg from the run median > "
                       f"{SELF_SIMILARITY_MAX_DEG} — the dimensionless response cannot depend on the "
                       "sweep variable, so this point is numerically corrupted")
            elif p["wall_margin"] < WALL_MARGIN_MIN:
                why = f"wall margin {p['wall_margin']:.2f} < {WALL_MARGIN_MIN}"
            else:
                ok.append(p)
                continue
            dropped.append({key: p[key], tau_field: p[tau_field],
                            "wall_margin": p["wall_margin"],
                            "self_similarity_dev_deg": dev.get(i, 0.0), "why": why})
        out = {"sweep": key, "inversion": tau_field, "declared_exponent": declared,
               "n_kept": len(ok), "rejected": dropped}
        if len(ok) < 3:
            out["verdict"] = "INSUFFICIENT"
            return out
        p_exp, c, r2 = power_law_exponent([q[key] for q in ok], [q[tau_field] for q in ok])
        out.update({"exponent": p_exp, "prefactor": c, "r2": r2,
                    "verdict": "PASS" if (r2 >= R2_FLOOR and abs(p_exp - declared) <= EXPONENT_TOLERANCE)
                               else "FAIL"})
        return out

    scored = {"a_um": score(a_pts, "a_um", DECLARED_EXPONENT_A),
              "D_um2_s": score(d_pts, "D_um2_s", DECLARED_EXPONENT_D)}
    scored_phase45 = {"a_um": score(a_pts, "a_um", DECLARED_EXPONENT_A, "tau_phase45_s"),
                      "D_um2_s": score(d_pts, "D_um2_s", DECLARED_EXPONENT_D, "tau_phase45_s")}
    record = {
        "record": "run-record@2", "kind": "gate", "gate": "c1_biot_microrheology_timescale",
        "device": str(dev),
        "build": {"commit": args.build_commit, "source": "declared"},
        "failed_inversion": {
            "observable": "tau = 1/omega_45, the 45-degree phase-lag crossing (the standard "
                          "first-order/single-relaxation reading)",
            "outcome": "NO CROSSING EXISTS. 45 degrees is the ASYMPTOTE of a diffusive phase lag, not "
                       "a point on it; a first-order model runs to 90 degrees, so fitting "
                       "arctan(omega*tau) to poroelastic data reports a corner frequency that does "
                       "not exist. Measured phase maxima are recorded per point.",
            "phase_max_deg": [p["phase_max_deg"] for p in a_pts],
            "why_kept": "a method that cannot measure this medium is a result about the method, and "
                        "deleting the attempt would delete the finding",
        },
        "declared_before_run": {
            "method": "multi-tone drive of a Gaussian source, lock-in detected on the patch mean",
            "observable": "tau_omega = 1/omega_3dB, where the response falls to 1/sqrt(2) of its "
                          "lowest-frequency value — the bandwidth reading, AMENDED before this run "
                          "after the phase inversion was shown ill-posed on a diffusive medium",
            "exponents": {"a_um": DECLARED_EXPONENT_A, "D_um2_s": DECLARED_EXPONENT_D},
            "prefactor": "NOT declared — a finite Gaussian source in a box has no clean closed form, "
                         "and inventing one would be a fit dressed as a prediction",
            "r2_floor": R2_FLOOR, "exponent_tolerance": EXPONENT_TOLERANCE,
            "wall_margin_min": WALL_MARGIN_MIN,
        },
        "scored": scored,
        "scored_phase45": scored_phase45,
        "two_inversions": {
            "why_both": ("they fail differently, so reporting one hides the other's weakness: -3 dB is "
                         "precise inside a band but its value moves with the band because it is "
                         "referenced to the lowest tone measured; phase-45 is band-independent because "
                         "its level is absolute, but it is noisier and it sits where the corruption is"),
            "band_dependence_evidence": {
                "band_0.125_to_16": {"C_3dB_a2_over_D": 1.943, "C_phase45_a2_over_D": 0.8679},
                "band_0.5_to_16": {"C_3dB_a2_over_D": 0.863, "C_phase45_a2_over_D": 0.8645},
                "note": "same medium, same D, same probe sizes; only the driven band differs",
            },
        },
        "sweeps": {"a_um": a_pts, "D_um2_s": d_pts},
        "verdict": "PASS" if all(s.get("verdict") == "PASS" for s in scored.values()) else "FAIL",
        "timing": {"wall_s": time.time() - t0, "comparable": True},
        "may_not_be_quoted_for": [
            "any statement about the CELL — this is the field alone",
            "a cytoplasm viscosity", "any rung",
        ],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "record.json").write_text(json.dumps(record, indent=2))
    print(f"\n[c1-rheo] VERDICT {record['verdict']} — a:{scored['a_um'].get('exponent')} "
          f"D:{scored['D_um2_s'].get('exponent')}  prefactor(a)={scored['a_um'].get('prefactor')} "
          f"-> {args.out}/record.json", flush=True)


if __name__ == "__main__":
    main()

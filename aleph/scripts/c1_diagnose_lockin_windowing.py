"""Reproduce run 34's lock-in windowing arithmetic on the host and look for what singles out a=1.1.

No physics: the field is replaced by an ANALYTIC first-order response per tone, so any deviation the
windowing produces is the windowing, isolated.
"""
import math
import sys

import numpy as np

sys.path.insert(0, "aleph/scripts")
from c1_biot_microrheology import _lock_in  # noqa: E402

D, DX, N = 50.0, 0.2, 241
TONES = (0.25, 0.5, 1.0, 2.0, 4.0)          # run 34 used these
N_PERIODS, N_DISCARD, SPP = 8, 4, 16
CFL = 0.2
TRUE_TAU_IN_A2D = 0.2                        # arbitrary, identical for every a (self-similar)

print(f"{'a':>5} {'n_samp':>8} {'T_slow/sdt':>12} {'steps/samp':>10} "
      f"{'periods_in_win':>15} {'phase err (deg)':>32}")
for a in (0.8, 1.1, 1.5, 2.0):
    om = np.array([m * D / (a * a) for m in TONES])
    T_slow = 2 * math.pi / om.min()
    t_end = N_PERIODS * T_slow
    sample_dt = (2 * math.pi / om.max()) / SPP
    dt_max = DX * DX / (6.0 * D)
    dt = min(CFL * dt_max, sample_dt)
    steps = max(1, int(round(sample_dt / dt)))
    dt = sample_dt / steps
    n_samples = int(round(t_end / sample_dt))

    t = (np.arange(n_samples) + 1) * sample_dt
    tau = TRUE_TAU_IN_A2D * a * a / D
    y = sum(np.sin(w * t - math.atan(w * tau)) / math.sqrt(1 + (w * tau) ** 2) for w in om)

    keep = t >= N_DISCARD * T_slow
    t_k, y_k = t[keep], y[keep]
    span = t_k[-1] - t_k[0]
    whole = math.floor(span / T_slow) * T_slow
    sel = t_k <= (t_k[0] + whole)
    t_k, y_k = t_k[sel], y_k[sel]

    errs = []
    for w in om:
        _, ph = _lock_in(t_k - t_k[0], y_k - y_k.mean(), w)
        ph = ph if ph >= 0 else ph + 2 * math.pi
        errs.append(math.degrees(ph) - math.degrees(math.atan(w * tau)))
    print(f"{a:>5} {n_samples:>8} {T_slow / sample_dt:>12.3f} {steps:>10} "
          f"{span / T_slow:>15.4f} {str([round(e, 2) for e in errs]):>32}")

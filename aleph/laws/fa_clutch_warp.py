"""FF active movement piece-3/5 — GPU-resident focal-adhesion CLUTCH (integrin catch-slip traction).

PI /goal 2026-07-02: traction = how a cell grips the substrate/ECM and pulls. The molecular-clutch picture
(Mitchison-Kirschner 1988; Chan-Odde 2008): retrograde actin flow / cortical contraction stretches integrin
bonds to a fixed substrate → the bond builds force → it ruptures at a LOAD-DEPENDENT rate (catch-slip) → de-
adhesion is EMERGENT, never a binary latch ([[feedback-junction-switch-fine-grained]]). This is the GPU-resident
port of the host-side `fa_anchor.py` scaffolding (INTEGRIN_A5B1, Kong 2009): a Hookean clutch spring from an
actin END node to a FIXED substrate anchor + a Pereverzev two-pathway catch-slip KMC turnover, both as Warp
kernels so the whole cell stays GPU-native.

Constants (µm·pN·s) — INTEGRIN_A5B1 (hand_kmc, Kong 2009, KB-2.5, PI-gated): k_int=1000 pN/µm, k_on=1/s,
Pereverzev k_catch0=0.4/s x_catch=0.61 nm, k_slip0=0.5/s x_slip=0.14 nm.
⚠️ These params give a catch-slip peak **F*=7 pN** (analytic F* = kT/(x_c+x_s)·ln[(k_c0 x_c)/(k_s0 x_s)]); the
KB claim / Kong 2009 report ~**30 pN**. NOT reconciled — surfaced to PI (do NOT retune x_c/x_s to hit 30 pN, that
would be fitting to a target). The kernel is validated against the analytic off-rate AS RECORDED; the 7-vs-30 pN
KB question is separate (KB content, PI-gated).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

from aleph.laws import units as U
from aleph.laws.hand_kmc import INTEGRIN_A5B1, pereverzev_off_rate


@wp.func
def _pereverzev(f: wp.float64, kc0: wp.float64, xc: wp.float64, ks0: wp.float64, xs: wp.float64,
                kT: wp.float64) -> wp.float64:
    """Pereverzev two-pathway catch-slip off-rate k_catch0·exp(−f·x_c/kT) + k_slip0·exp(+f·x_s/kT)."""
    return kc0 * wp.exp(-f * xc / kT) + ks0 * wp.exp(f * xs / kT)


@wp.kernel
def clutch_spring_kernel(
    pos: wp.array(dtype=wp.vec3d),
    actin: wp.array(dtype=wp.int32),            # (M,) actin end-node index per clutch
    anchor: wp.array(dtype=wp.vec3d),           # (M,) FIXED substrate anchor position [µm]
    bound: wp.array(dtype=wp.int32),            # (M,) 1=engaged, 0=detached
    k_int: wp.float64,                          # clutch stiffness [pN/µm]
    rest: wp.float64,                           # clutch rest length [µm]
    force: wp.array(dtype=wp.vec3d),
):
    """Engaged clutch = Hookean spring pulling the actin node toward its fixed substrate anchor:
    f = k_int·(L−rest)/L · (anchor−actin). Detached clutches exert nothing. The anchor is FIXED (substrate),
    so this is a traction force on the cell (reaction goes to the immovable ECM)."""
    t = wp.tid()
    if bound[t] == 0:
        return
    a = actin[t]
    d = anchor[t] - pos[a]
    L = wp.length(d)
    if L > wp.float64(1e-12):
        wp.atomic_add(force, a, (k_int * (L - rest) / L) * d)


@wp.kernel
def clutch_catchslip_kmc_kernel(
    pos: wp.array(dtype=wp.vec3d),
    actin: wp.array(dtype=wp.int32),
    anchor: wp.array(dtype=wp.vec3d),
    bound: wp.array(dtype=wp.int32),            # (M,) in/out engaged state
    k_int: wp.float64, rest: wp.float64,
    kc0: wp.float64, xc: wp.float64, ks0: wp.float64, xs: wp.float64, kT: wp.float64,
    cap: wp.float64,                            # re-attach capture radius [µm]
    k_on: wp.float64, tau: wp.float64, seed: wp.int32,
):
    """One KMC tick: a BOUND clutch on load F=k_int·|L−rest| detaches with 1−exp(−τ·off_rate(F)) (Pereverzev
    catch-slip → de-adhesion EMERGENT from load); a DETACHED clutch within ``cap`` re-attaches with
    1−exp(−τ·k_on). One thread per clutch (independent RNG stream)."""
    t = wp.tid()
    rs = wp.rand_init(seed, t)
    a = actin[t]
    L = wp.length(anchor[t] - pos[a])
    F = k_int * wp.abs(L - rest)
    if bound[t] == 1:
        p_det = wp.float64(1.0) - wp.exp(-tau * _pereverzev(F, kc0, xc, ks0, xs, kT))
        if wp.float64(wp.randf(rs)) < p_det:
            bound[t] = wp.int32(0)
    else:
        if L < cap:
            if wp.float64(wp.randf(rs)) < wp.float64(1.0) - wp.exp(-tau * k_on):
                bound[t] = wp.int32(1)


@dataclass(frozen=True, slots=True)
class ClutchParams:
    """Resolved integrin-clutch constants (µm·pN·s)."""

    k_int: float; rest_um: float; kc0: float; xc_um: float; ks0: float; xs_um: float
    kT: float; cap_um: float; k_on: float
    F_star_pN: float                            # analytic catch-slip peak (min off-rate)


def resolve_clutch(integrin=INTEGRIN_A5B1, *, rest_um=0.03):
    """Build ClutchParams from a HandParams integrin preset; compute the analytic catch-slip peak F*."""
    kc0, xc, ks0, xs, kT = (integrin.cs_k_catch0, integrin.cs_x_catch_um, integrin.cs_k_slip0,
                            integrin.cs_x_slip_um, U.KBT)
    F_star = kT / (xc + xs) * np.log((kc0 * xc) / (ks0 * xs))    # d(off_rate)/dF = 0
    return ClutchParams(k_int=integrin.link_k, rest_um=rest_um, kc0=kc0, xc_um=xc, ks0=ks0, xs_um=xs,
                        kT=kT, cap_um=integrin.capture_radius_um, k_on=integrin.k_on, F_star_pN=float(F_star))


def clutch_off_rate_np(f, cp: ClutchParams):
    """Analytic off-rate at load f [pN] (host reference for the kernel)."""
    return pereverzev_off_rate(f, cp.kc0, cp.xc_um, cp.ks0, cp.xs_um, cp.kT)

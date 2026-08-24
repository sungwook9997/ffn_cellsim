r"""Per-head NMII hand/KMC device API (I3 OWNS the §1.4 frozen interface) — Warp CUDA source.

This is the device-resident kinetic building block of the active-force root: one *hand* is one myosin head,
carrying its own bound/free state, actin anchor, and abscissa entirely in GPU arrays with NO authoritative
per-step host state (I0-A). Ported from ``ff/hand_kmc.py`` (host NumPy) to device-resident per-head state;
the closed-form device functions here MATCH the host oracles in ``bell_kinetics_analytic`` /
``hill_fv_analytic`` bit-for-formula so the native gate reproduces the CPU-green analytic gates.

Frozen interface (§1.4): downstream tracks (I4 weave, I6 adhesion, I8 mt-coupling) CONSTRUCT hands via
:class:`NMIIHandParams` + the kernels here in their OWN ``presets.py``; they never edit this KMC core.

Device functions (``@wp.func``) — the mechanistic laws, identical to the host oracles:
  * :func:`bell_off_rate`     p_off(f) = k_off0 * exp(|f| / f0)              (Bell 1978 slip)
  * :func:`catch_slip_off_rate` p_off(f) = k_c0*exp(-|f|*x_c/kT) + k_s0*exp(+|f|*x_s/kT)  (Pereverzev catch-slip)
  * :func:`hill_velocity`     v(F) = v0 (1 - F/F_s)/(1 + (F/F_s)/kappa)      (Hill 1938, kappa->inf = linear)
  * :func:`attach_prob`       1 - exp(-tau k_on)                            (Poisson per tick)
  * :func:`detach_prob`       1 - exp(-tau p_off)

⚠ FIDELITY CORRECTION (2026-07-23, sourcing §b): the NMII head-actin bond is CATCH-SLIP, not pure Bell slip.
Kovacs, Thirumurugan, Knight & Sellers 2007, *PNAS* 104(24):9994 (10.1073/pnas.0701181104): resistive load
RAISES the bound-head duty (ADP release, the rate-limiting detachment step, is slowed 5x NM2A / 12x NM2B under
load) — the OPPOSITE sign of a slip bond. :func:`step_detach_kernel` (pure Bell slip) is retained UNCHANGED
so the ERM / alpha-actinin / crosslink consumers of ``bell_off_rate`` are untouched; the corrected NMII head
turnover uses :func:`step_detach_catch_slip_kernel` with :class:`NMIICatchSlipParams`.

KMC sub-step kernels (device RNG, no host draw):
  * :func:`attach_kernel`      free hand within capture radius binds with prob attach_prob
  * :func:`step_detach_kernel` bound hand advances abscissa by tau*v(F) and detaches by Bell SLIP off-rate
  * :func:`step_detach_catch_slip_kernel` same Hill step, detaches by the CATCH-SLIP off-rate (NMII head, Kovacs
    2007). ``attach_kernel`` + this kernel = the full corrected per-head NMII attach/detach/turnover cycle.

⚠ Runtime: NVIDIA Warp on CUDA GPU only (I0-A). Authored on the dev Mac (no CUDA) but NOT launched here — the
native KMC gate (GPU-residency / ensemble-stall emergence) runs on the lead's gbook A5000. HOOMD is never
imported (the archived cortex/myosin.py is a read-only port spec).

Sanity Gate (native, lead): per-head Hill FV == host Hill-1938 oracle; Bell slip monotone; NMII catch-slip
off-rate biphasic (falls then rises, min at F*); engaged fraction emerges from k_on/k_off (NOT an imposed
duty) and RISES under load for the catch-slip head; zero authoritative GPU->CPU roundtrip inside the loop.
"""

from __future__ import annotations

import warp as wp

# engine units: length um, force pN, time s, stiffness pN/um, rate 1/s (matches ff/ and params_i0b3.yaml)


@wp.struct
class NMIIHandParams:
    """Kinetic + mechanical constants of one NMII head type (device-resident scalars, engine units).

    All values are I0-B3 evidence decisions (params_i0b3.yaml); several are GAP — PI and must be supplied by
    the lead at the native gate (never chosen to lift the gamma-floor).
    """

    k_on: wp.float64          # per-head attach rate [1/s]
    k_off0: wp.float64        # zero-force detachment rate [1/s]
    f0: wp.float64            # Bell characteristic force [pN] = kBT / x_beta
    v0: wp.float64            # unloaded stepping velocity [um/s]
    f_stall: wp.float64       # per-head isometric stall force [pN]
    kappa: wp.float64         # Hill curvature a/F0 [-] (kappa -> inf recovers the linear law)
    k_xb: wp.float64          # crossbridge (head->actin) stiffness [pN/um] (MASTER force knob)
    r0_head: wp.float64       # head<->BACKBONE arm rest length [um] (= topology head_offset_um; NOT the crossbridge rest)
    r0_xb: wp.float64         # head<->ACTIN crossbridge rest length [um] (~0: a bound head sits on the actin site)
    capture_radius: wp.float64  # actin capture radius [um]


@wp.struct
class NMIICatchSlipParams:
    """Pereverzev two-pathway catch-slip constants for the NMII head-actin bond (device-resident, engine units).

    The corrected NMII detachment law (Kovacs 2007): the catch pathway is the load-slowed ADP-release
    detachment (off-rate falls with load up to F*), the slip pathway is forced unbinding (off-rate rises past
    F*). Every value is an I0-B3 GAP — PI (magnitudes constrained by Kovacs 5x/12x slowdown + NM2B duty 0.2-0.3,
    but no direct single-molecule fit); NEVER chosen to lift a force/gamma band.
    """

    k_catch0: wp.float64      # zero-force catch-pathway rate [1/s]
    x_catch: wp.float64       # catch-pathway bond length [um] (larger => stronger load-strengthening)
    k_slip0: wp.float64       # zero-force slip-pathway rate [1/s]
    x_slip: wp.float64        # slip-pathway (Bell) bond length [um]
    kT: wp.float64            # thermal energy [pN*um]


@wp.func
def bell_off_rate(f: wp.float64, k_off0: wp.float64, f0: wp.float64) -> wp.float64:
    """Bell slip off-rate ``p_off = k_off0 * exp(|f| / f0)`` [1/s] (matches host ``bell_off_rate``)."""
    return k_off0 * wp.exp(wp.abs(f) / f0)


@wp.func
def catch_slip_off_rate(
    f: wp.float64,
    k_catch0: wp.float64,
    x_catch: wp.float64,
    k_slip0: wp.float64,
    x_slip: wp.float64,
    kT: wp.float64,
) -> wp.float64:
    """Pereverzev two-pathway catch-slip off-rate [1/s] (matches host ``bell_kinetics_analytic.catch_slip_off_rate``).

    ``p_off = k_catch0*exp(-|f|*x_catch/kT) + k_slip0*exp(+|f|*x_slip/kT)`` — falls with load (catch) up to the
    peak-lifetime force F*, then rises (slip). The NMII head-actin correction (Kovacs 2007); same form as the
    device ``erm_spectrum._erm_catch_slip_off_rate`` and the FF ``pereverzev_off_rate`` SoT.
    """
    fa = wp.abs(f)
    return k_catch0 * wp.exp(-fa * x_catch / kT) + k_slip0 * wp.exp(fa * x_slip / kT)


@wp.func
def hill_velocity(f: wp.float64, v0: wp.float64, f_stall: wp.float64, kappa: wp.float64) -> wp.float64:
    """Per-head Hill shortening velocity [um/s], clamped >= 0 past stall (matches host ``hill_velocity``).

    ``v = v0 (1 - F/F_s) / (1 + (F/F_s)/kappa)``. kappa -> inf is the PI-2026-07-07 linear law.
    """
    fr = f / f_stall
    v = v0 * (wp.float64(1.0) - fr) / (wp.float64(1.0) + fr / kappa)
    return wp.max(v, wp.float64(0.0))


@wp.func
def attach_prob(tau: wp.float64, k_on: wp.float64) -> wp.float64:
    """Per-tick attach probability ``1 - exp(-tau k_on)`` (Poisson)."""
    return wp.float64(1.0) - wp.exp(-tau * k_on)


@wp.func
def detach_prob(tau: wp.float64, p_off: wp.float64) -> wp.float64:
    """Per-tick detach probability ``1 - exp(-tau p_off)``."""
    return wp.float64(1.0) - wp.exp(-tau * p_off)


@wp.kernel
def attach_kernel(
    bound: wp.array(dtype=wp.int32),
    anchor: wp.array(dtype=wp.int32),
    abscissa: wp.array(dtype=wp.float64),
    nearest_dist: wp.array(dtype=wp.float64),
    nearest_idx: wp.array(dtype=wp.int32),
    params: NMIIHandParams,
    tau: wp.float64,
    rng_seed: wp.int32,
):
    """One attach sub-step (NF2007 attach(m), device): free hands within the capture radius bind stochastically.

    Device RNG per head (``wp.rand_init(rng_seed, tid)``) — no host random draw, no host state written.
    """
    h = wp.tid()
    if bound[h] == 1:
        return
    if nearest_dist[h] > params.capture_radius:
        return
    state = wp.rand_init(rng_seed, h)
    p = attach_prob(tau, params.k_on)
    if wp.float64(wp.randf(state)) < p:
        bound[h] = 1
        anchor[h] = nearest_idx[h]
        abscissa[h] = wp.float64(0.0)


@wp.kernel
def step_detach_kernel(
    bound: wp.array(dtype=wp.int32),
    anchor: wp.array(dtype=wp.int32),
    abscissa: wp.array(dtype=wp.float64),
    loads_hill: wp.array(dtype=wp.float64),
    loads_bell: wp.array(dtype=wp.float64),
    params: NMIIHandParams,
    tau: wp.float64,
    rng_seed: wp.int32,
):
    r"""One step+detach sub-step (NF2007 step(f), device) — the F4/F5-split step + Bell detachment.

    Two DISTINCT per-head loads from the inner mechanical solve
    (:func:`~aleph.components.motor.minifilament_warp.compute_head_loads_kernel`):

      * ``loads_hill[h]`` — the TANGENTIAL (walk-direction) crossbridge tension that resists the WALK. This
        drives the Hill stepping velocity. The abscissa advanced here is the SAME walked displacement that
        advances the crossbridge attachment point along ``walk_dir`` in the force kernels — so the power stroke
        enters the force, the directed contraction emerges, and force-velocity self-limits: as the walked
        strain lifts ``loads_hill`` toward ``f_stall``, ``hill_velocity`` → 0 and the abscissa (hence the
        force) stalls at ``f_stall``.
      * ``loads_bell[h]`` — the FULL crossbridge tension magnitude |F|. A crossbridge slips under its TOTAL
        load (Bell 1978 reads the bond force magnitude), NOT just the along-actin component (F4). At collinear
        isometric stall the transverse part ~0 so the two coincide; off-axis they differ and the slip must see
        the whole force. Bell slip uses the SAME device off-rate as the host oracle, so the engaged fraction
        EMERGES from k_on/k_off (P2), never an imposed duty ratio.
    """
    h = wp.tid()
    if bound[h] == 0:
        return
    f_hill = loads_hill[h]
    f_bell = loads_bell[h]
    # DERIVED stall-crossing guard (Codex/Claude overshoot fix): a bare Δa = τ·v(F_old) can walk PAST isometric
    # stall in one tick (F_old just below f_stall ⇒ v>0, but the settled load overshoots f_stall). Cap Δa by the
    # Hill-stall bound Δa ≤ (f_stall − F)/k_xb so the along-actin tension cannot cross f_stall within a tick —
    # an implicit stall boundary, NOT an empirical clamp. (v is already ≥0 past stall; this bounds the approach.)
    da = tau * hill_velocity(f_hill, params.v0, params.f_stall, params.kappa)
    da_max = wp.max(wp.float64(0.0), (params.f_stall - f_hill) / params.k_xb)
    abscissa[h] = abscissa[h] + wp.min(da, da_max)
    p_off = bell_off_rate(f_bell, params.k_off0, params.f0)
    pdet = detach_prob(tau, p_off)
    state = wp.rand_init(rng_seed, h)
    if wp.float64(wp.randf(state)) < pdet:
        bound[h] = 0
        anchor[h] = -1
        abscissa[h] = wp.float64(0.0)


@wp.kernel
def step_detach_catch_slip_kernel(
    bound: wp.array(dtype=wp.int32),
    anchor: wp.array(dtype=wp.int32),
    abscissa: wp.array(dtype=wp.float64),
    loads_hill: wp.array(dtype=wp.float64),
    loads_bell: wp.array(dtype=wp.float64),
    params: NMIIHandParams,
    cs: NMIICatchSlipParams,
    tau: wp.float64,
    rng_seed: wp.int32,
):
    r"""One step+detach sub-step for the CORRECTED NMII head — Hill step + Pereverzev CATCH-SLIP detachment.

    Identical to :func:`step_detach_kernel` (same F4/F5 load split, same Hill stepping, same derived
    stall-crossing guard) EXCEPT the detachment off-rate is :func:`catch_slip_off_rate` (from ``cs``), not the
    pure Bell slip. This is the Kovacs 2007 fix: the ``loads_bell`` full-|F| load now SLOWS detachment (catch)
    up to the peak-lifetime force F*, so the engaged fraction RISES under resistive load for tension
    maintenance, then falls past F* (slip). ``params`` supplies only the Hill mechanics (v0, f_stall, kappa,
    k_xb); its ``k_off0``/``f0`` are unused here — the off-rate comes entirely from ``cs``. The engaged
    fraction still EMERGES from k_on/k_off (P2), never an imposed duty.
    """
    h = wp.tid()
    if bound[h] == 0:
        return
    f_hill = loads_hill[h]
    f_bell = loads_bell[h]
    # DERIVED stall-crossing guard — identical to step_detach_kernel: cap Δa by the Hill-stall bound
    # Δa ≤ (f_stall − F)/k_xb so the along-actin tension cannot cross f_stall within a tick.
    da = tau * hill_velocity(f_hill, params.v0, params.f_stall, params.kappa)
    da_max = wp.max(wp.float64(0.0), (params.f_stall - f_hill) / params.k_xb)
    abscissa[h] = abscissa[h] + wp.min(da, da_max)
    p_off = catch_slip_off_rate(f_bell, cs.k_catch0, cs.x_catch, cs.k_slip0, cs.x_slip, cs.kT)
    pdet = detach_prob(tau, p_off)
    state = wp.rand_init(rng_seed, h)
    if wp.float64(wp.randf(state)) < pdet:
        bound[h] = 0
        anchor[h] = -1
        abscissa[h] = wp.float64(0.0)


def allocate_hand_state(n_hands: int, device: str | None = None):
    """Allocate an all-free device hand population (bound=0, anchor=-1, abscissa=0, walk_dir=0).

    Args:
        n_hands: total heads across all minifilaments (e.g. ``n_minifilaments * 2 * N_side``).
        device: Warp CUDA device alias; ``None`` uses the current device without fixing an ordinal.

    Returns:
        dict of device arrays: ``bound`` (int32), ``anchor`` (int32), ``abscissa`` (float64),
        ``walk_dir`` (vec3d).

    ``walk_dir[h]`` is the unit direction along the bound actin filament the head walks toward
    (its barbed-end polarity). The crossbridge attachment point advances along ``walk_dir`` by the
    walked ``abscissa`` (:func:`~aleph.components.motor.minifilament_warp.crossbridge_kernel`), so the
    power stroke becomes a directed force. A ZERO ``walk_dir`` makes the head passive (no directed
    force, only the compliance spring) — the OFF/regression default; the production path MUST fill it
    from the head's actin barbed-end polarity (I4-weave hand-off, INTEGRATION.md §1).
    """
    return {
        "bound": wp.zeros(n_hands, dtype=wp.int32, device=device),
        "anchor": wp.full(n_hands, -1, dtype=wp.int32, device=device),
        "abscissa": wp.zeros(n_hands, dtype=wp.float64, device=device),
        "walk_dir": wp.zeros(n_hands, dtype=wp.vec3d, device=device),
    }


__all__ = [
    "NMIIHandParams",
    "NMIICatchSlipParams",
    "bell_off_rate",
    "catch_slip_off_rate",
    "hill_velocity",
    "attach_prob",
    "detach_prob",
    "attach_kernel",
    "step_detach_kernel",
    "step_detach_catch_slip_kernel",
    "allocate_hand_state",
]
